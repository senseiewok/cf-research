#!/usr/bin/env python3
"""Independent verification audit for benchmarks/cftr_dms_canx_selectivity/.

Separate pure-stdlib code importing NOTHING from `cftr_dms_canx_selectivity.py` or `tools/`:
1. Validates `done_sha256`, locked `split_manifest_sha256`, codon disjointness, and zero
   `PRED-0010` codons in holdout.
2. Checks arithmetic parity and Bonferroni verdict rules for both Arm A and Arm B.
3. When local Tedman et al. 2025 `.xlsx` files exist under `sources/downloads/tedman-2025/`,
   independently re-parses `All_variants` & `All_variants CANX KO` and re-solves Arm A (2x2 &
   3x3) and Arm B (2x2 & 10x10) via Cholesky L L^T decomposition (`tol <= 1e-6`).
"""
from __future__ import annotations

import argparse
from decimal import Decimal
import hashlib
import json
import math
from pathlib import Path
import re
import sys
from typing import Any, Dict, List, Optional, Tuple
import xml.etree.ElementTree as ET
import zipfile

REPO_ROOT = Path(__file__).resolve().parents[2]
MASTER_XLSX = "ALT 9-6-24 CFTR DMS Expanded Library Master Tables.xlsx"
NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
DOMAINS = ("TMD1", "NBD1", "R Domain", "TMD2", "NBD2")
ALL_DOMAINS = ("Lasso",) + DOMAINS
CLASSES = ("nonsense", "inframe_indel", "other")
EPS = 0.05


def _locate_data_dir(explicit: Optional[Path], repo_root: Path) -> Optional[Path]:
    if explicit is not None:
        return Path(explicit)
    direct = repo_root / "sources" / "downloads" / "tedman-2025"
    if (direct / MASTER_XLSX).is_file():
        return direct
    if len(repo_root.parents) >= 3:
        cand = repo_root.parents[2] / "cf-research" / "sources" / "downloads" / "tedman-2025"
        if (cand / MASTER_XLSX).is_file():
            return cand
    return None


def _col_to_idx(col: str) -> int:
    n = 0
    for ch in col:
        n = n * 26 + (ord(ch) - ord("A") + 1)
    return n


def _idx_to_col(idx: int) -> str:
    out: List[str] = []
    n = idx
    while n > 0:
        n, rem = divmod(n - 1, 26)
        out.append(chr(ord("A") + rem))
    return "".join(reversed(out))


def _read_sheet(xlsx_path: Path, sheet_name: str, *, expand_merged: bool = True) -> Dict[Tuple[int, str], str]:
    with zipfile.ZipFile(xlsx_path, "r") as zf:
        strings: List[str] = []
        if "xl/sharedStrings.xml" in zf.namelist():
            sst = ET.fromstring(zf.read("xl/sharedStrings.xml"))
            for si in sst.findall(f"{NS}si"):
                strings.append("".join(t.text or "" for t in si.iter(f"{NS}t")))
        wb = ET.fromstring(zf.read("xl/workbook.xml"))
        rels = ET.fromstring(zf.read("xl/_rels/workbook.xml.rels"))
        rid_map = {r.attrib["Id"]: r.attrib["Target"] for r in rels}
        sheet_xml: Optional[str] = None
        for sh in wb.findall(f".//{NS}sheet"):
            if sh.attrib.get("name") == sheet_name:
                rid = sh.attrib["{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"]
                t = rid_map[rid].lstrip("/")
                sheet_xml = t if t.startswith("xl/") else f"xl/{t}"
                break
        if sheet_xml is None:
            raise ValueError(f"Missing sheet {sheet_name}")

        sroot = ET.fromstring(zf.read(sheet_xml))
        grid: Dict[Tuple[int, str], str] = {}
        for row in sroot.findall(f".//{NS}sheetData/{NS}row"):
            for cell in row.findall(f"{NS}c"):
                coord = cell.attrib.get("r", "")
                m = re.fullmatch(r"([A-Z]+)(\d+)", coord)
                if not m:
                    raise ValueError(f"Invalid coordinate {coord!r}")
                col, rnum = m.group(1), int(m.group(2))
                ctype = cell.attrib.get("t", "")
                val: Optional[str] = None
                if ctype == "inlineStr":
                    val = "".join(t.text or "" for t in cell.iter(f"{NS}t"))
                else:
                    ve = cell.find(f"{NS}v")
                    if ve is not None and ve.text is not None:
                        val = strings[int(ve.text)] if ctype == "s" else ve.text
                if val is not None:
                    grid[(rnum, col)] = val

        if expand_merged:
            for mc in sroot.findall(f".//{NS}mergeCells/{NS}mergeCell"):
                m = re.fullmatch(r"([A-Z]+)(\d+):([A-Z]+)(\d+)", mc.attrib.get("ref", ""))
                if not m:
                    continue
                c1, r1, c2, r2 = (
                    _col_to_idx(m.group(1)),
                    int(m.group(2)),
                    _col_to_idx(m.group(3)),
                    int(m.group(4)),
                )
                top = grid.get((r1, _idx_to_col(c1)))
                if top is not None:
                    for rr in range(r1, r2 + 1):
                        for cc in range(c1, c2 + 1):
                            grid.setdefault((rr, _idx_to_col(cc)), top)
        return grid


def _cholesky_solve_spd(mat: List[List[float]], rhs: List[float]) -> List[float]:
    n = len(rhs)
    chol = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1):
            acc = mat[i][j] - sum(chol[i][k] * chol[j][k] for k in range(j))
            if i == j:
                if acc <= 0.0:
                    raise ValueError("Matrix not positive definite.")
                chol[i][j] = math.sqrt(acc)
            else:
                chol[i][j] = acc / chol[j][j]
    y = [0.0] * n
    for i in range(n):
        y[i] = (rhs[i] - sum(chol[i][k] * y[k] for k in range(i))) / chol[i][i]
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        x[i] = (y[i] - sum(chol[k][i] * x[k] for k in range(i + 1, n))) / chol[i][i]
    return x


def _decimal_pearson(xs: List[float], ys: List[float]) -> float:
    dx_list = [Decimal(str(v)) for v in xs]
    dy_list = [Decimal(str(v)) for v in ys]
    n = Decimal(len(dx_list))
    mx = sum(dx_list) / n
    my = sum(dy_list) / n
    num = sum((x - mx) * (y - my) for x, y in zip(dx_list, dy_list))
    vx = sum((x - mx) ** 2 for x in dx_list)
    vy = sum((y - my) ** 2 for y in dy_list)
    return float(num) / math.sqrt(float(vx) * float(vy))


def _check_verdict_math(he: Dict[str, Any], base_key: str = "baseline_dmso_ols_mae_log2", cand_key: str = "candidate_mae_log2") -> None:
    base = float(he[base_key])
    cand = float(he[cand_key])
    recomputed = round((base - cand) / base, 6)
    # Propagate 0.5e-6 rounding uncertainty of the printed 6-decimal numerator/denominator values.
    max_round_err = (1e-6 * (1.0 + abs(base - cand) / base)) / base + 1e-6
    if abs(recomputed - float(he["relative_mae_reduction"])) > max_round_err:
        raise AssertionError(f"Relative reduction mismatch: {recomputed} vs {he['relative_mae_reduction']}")
    ci_ok = bool(he["residue_clustered_bootstrap"]["ci975_low_gt_zero"])
    if recomputed >= 0.10 and ci_ok:
        if he["ledger_standing"] != "open":
            raise AssertionError("Expected open standing when both gate and CI97.5 pass")
    else:
        if he["ledger_standing"] != "refuted":
            raise AssertionError("Expected refuted standing when threshold or CI97.5 fails")


def run_independent_audit(repo_root: Path = REPO_ROOT, data_dir: Optional[Path] = None) -> Dict[str, Any]:
    bench_dir = repo_root / "benchmarks" / "cftr_dms_canx_selectivity"
    results = json.loads((bench_dir / "results_2026-11-23.json").read_text(encoding="utf-8"))
    manifest = json.loads((repo_root / "benchmarks/cftr_dms_expression/split_manifest.json").read_text(encoding="utf-8"))
    done_txt = (bench_dir / "done.md").read_text(encoding="utf-8", errors="replace").replace("\r\n", "\n")

    if results["done_sha256"] != hashlib.sha256(done_txt.encode("utf-8")).hexdigest():
        raise AssertionError("done_sha256 mismatch.")
    if results["split_manifest_sha256"] != hashlib.sha256(json.dumps(manifest, sort_keys=True).encode("utf-8")).hexdigest():
        raise AssertionError("split_manifest_sha256 mismatch.")

    th = json.loads((repo_root / "benchmarks/theratype/theratype_open_annotations_2026-11-21.json").read_text(encoding="utf-8"))
    exposed = {int(v["codon_position"]) for v in th["variants"]}
    assigns = manifest["assignments"]
    dev_codons = {int(a["codon_position"]) for a in assigns if a["split"] == "dev"}
    hold_codons = {int(a["codon_position"]) for a in assigns if a["split"] == "holdout"}
    if dev_codons & hold_codons or hold_codons & exposed:
        raise AssertionError("Leakage check failed in locked split manifest.")

    ha = results["arm_a_canx_chaperone_coupling"]["holdout_evaluation"]
    hb = results["arm_b_canx_ko_single_corrector_selectivity"]["holdout_evaluation"]
    _check_verdict_math(ha, base_key="baseline_dmso_ols_mae_log2", cand_key="candidate_mae_log2")
    _check_verdict_math(hb, base_key="locked_baseline_mae_s_eps", cand_key="candidate_partial_ridge_mae_s_eps")

    ddir = _locate_data_dir(data_dir, repo_root)
    raw_checked = False
    recomputed: Dict[str, Any] = {}
    if ddir is not None:
        par = _read_sheet(ddir / MASTER_XLSX, "All_variants", expand_merged=True)
        ko = _read_sheet(ddir / MASTER_XLSX, "All_variants CANX KO", expand_merged=True)
        by_row = {int(a["row_num"]): a for a in assigns}
        max_r = max(r for r, _ in par)
        par_sub = ""
        ko_sub = ""
        dev_rows: List[Tuple[str, str, float, float, float, float]] = []
        hold_rows: List[Tuple[str, str, float, float, float, float]] = []
        for rnum in range(3, max_r + 1):
            if (rnum, "A") in par and par[(rnum, "A")].strip():
                par_sub = par[(rnum, "A")].strip()
            if (rnum, "A") in ko and ko[(rnum, "A")].strip():
                ko_sub = ko[(rnum, "A")].strip()
            p_nuc = par.get((rnum, "B"), "").strip()
            p_prot = par.get((rnum, "C"), "").strip()
            k_nuc = ko.get((rnum, "B"), "").strip()
            k_prot = ko.get((rnum, "C"), "").strip()
            if not any((p_nuc, p_prot, k_nuc, k_prot)):
                continue
            if (par_sub, p_nuc, p_prot) != (ko_sub, k_nuc, k_prot) or par_sub not in ALL_DOMAINS:
                raise AssertionError(f"Dual-sheet mismatch at row {rnum}")
            info = by_row[rnum]
            vals = [
                par.get((rnum, "F")),
                par.get((rnum, "N")),
                par.get((rnum, "R")),
                ko.get((rnum, "F")),
                ko.get((rnum, "J")),
                ko.get((rnum, "N")),
                ko.get((rnum, "R")),
            ]
            if any(v is None or float(v) <= 0.0 for v in vals):
                continue
            p_dmso, _, p_combo, k_dmso, k_661, k_445, _ = (float(v) for v in vals)
            rec = (
                info["sub_domain"],
                info["variant_class"],
                math.log2(p_dmso),
                math.log2(k_dmso),
                math.log2(p_combo),
                math.log2((k_445 + EPS) / (k_661 + EPS)),
            )
            if info["split"] == "dev":
                dev_rows.append(rec)
            else:
                hold_rows.append(rec)

        r_canx = _decimal_pearson([x[2] for x in dev_rows], [x[3] for x in dev_rows])
        rep_r = float(results["pre_holdout_controls"]["dev_canx_ko_dmso_stability"]["pearson_r_log2"])
        if abs(r_canx - rep_r) > 1e-6:
            raise AssertionError(f"CANX stability Pearson mismatch: {r_canx} != {rep_r}")

        # Arm A Cholesky 2x2 & 3x3
        n_d = float(len(dev_rows))
        x_a2 = [[1.0, r[2]] for r in dev_rows]
        xtx_a2 = [[sum(x_a2[i][a] * x_a2[i][b] for i in range(len(dev_rows))) for b in range(2)] for a in range(2)]
        xty_a2 = [sum(x_a2[i][a] * dev_rows[i][4] for i in range(len(dev_rows))) for a in range(2)]
        b_a2 = _cholesky_solve_spd(xtx_a2, xty_a2)
        mae_a2 = sum(abs(r[4] - (b_a2[0] + b_a2[1] * r[2])) for r in hold_rows) / len(hold_rows)
        if abs(mae_a2 - float(ha["baseline_dmso_ols_mae_log2"])) > 1e-6:
            raise AssertionError(f"Arm A baseline MAE mismatch: {mae_a2}")

        x_a3 = [[1.0, r[2], r[3]] for r in dev_rows]
        xtx_a3 = [[sum(x_a3[i][a] * x_a3[i][b] for i in range(len(dev_rows))) for b in range(3)] for a in range(3)]
        xty_a3 = [sum(x_a3[i][a] * dev_rows[i][4] for i in range(len(dev_rows))) for a in range(3)]
        b_a3 = _cholesky_solve_spd(xtx_a3, xty_a3)
        mae_a3 = sum(abs(r[4] - (b_a3[0] + b_a3[1] * r[2] + b_a3[2] * r[3])) for r in hold_rows) / len(hold_rows)
        if abs(mae_a3 - float(ha["candidate_mae_log2"])) > 1e-6:
            raise AssertionError(f"Arm A candidate MAE mismatch: {mae_a3}")

        # Arm B Cholesky 2x2 (CANX_KO_DMSO) & 10x10 partial Ridge
        x_b2 = [[1.0, r[3]] for r in dev_rows]
        xtx_b2 = [[sum(x_b2[i][a] * x_b2[i][b] for i in range(len(dev_rows))) for b in range(2)] for a in range(2)]
        xty_b2 = [sum(x_b2[i][a] * dev_rows[i][5] for i in range(len(dev_rows))) for a in range(2)]
        b_b2 = _cholesky_solve_spd(xtx_b2, xty_b2)
        mae_b2 = sum(abs(r[5] - (b_b2[0] + b_b2[1] * r[3])) for r in hold_rows) / len(hold_rows)
        if abs(mae_b2 - float(hb["canx_ko_dmso_ols_mae_s_eps"])) > 1e-6:
            raise AssertionError(f"Arm B OLS MAE mismatch: {mae_b2}")

        alpha = float(
            results["arm_b_canx_ko_single_corrector_selectivity"]["dev_models"]["candidate_partial_ridge"]["selected_alpha"]
        )
        raw_d = [[1.0 if r[0] == s else 0.0 for s in DOMAINS] + [1.0 if r[1] == v else 0.0 for v in CLASSES] for r in dev_rows]
        m = len(raw_d[0])
        means = [sum(row[j] for row in raw_d) / n_d for j in range(m)]
        stds = [math.sqrt(sum((row[j] - means[j]) ** 2 for row in raw_d) / n_d) for j in range(m)]
        xmat = [
            [1.0, dev_rows[i][3]] + [(raw_d[i][j] - means[j]) / stds[j] if stds[j] > 1e-12 else 0.0 for j in range(m)]
            for i in range(len(dev_rows))
        ]
        p = 2 + m
        xtx_b = [[sum(xmat[i][a] * xmat[i][b] for i in range(len(dev_rows))) for b in range(p)] for a in range(p)]
        for j in range(2, p):
            xtx_b[j][j] += alpha if stds[j - 2] > 1e-12 else 1.0
        xty_b = [sum(xmat[i][a] * dev_rows[i][5] for i in range(len(dev_rows))) for a in range(p)]
        beta_b = _cholesky_solve_spd(xtx_b, xty_b)
        cand_errs_b: List[float] = []
        for sd, vc, _, xko, _, sy in hold_rows:
            rh = [1.0 if sd == s else 0.0 for s in DOMAINS] + [1.0 if vc == v else 0.0 for v in CLASSES]
            pred = beta_b[0] + beta_b[1] * xko + sum(
                beta_b[2 + j] * ((rh[j] - means[j]) / stds[j] if stds[j] > 1e-12 else 0.0) for j in range(m)
            )
            cand_errs_b.append(abs(sy - pred))
        mae_b_cand = sum(cand_errs_b) / len(hold_rows)
        if abs(mae_b_cand - float(hb["candidate_partial_ridge_mae_s_eps"])) > 1e-6:
            raise AssertionError(f"Arm B candidate MAE mismatch: {mae_b_cand}")

        raw_checked = True
        recomputed = {
            "dev_canx_stability_r": round(r_canx, 6),
            "arm_a_baseline_mae_log2": round(mae_a2, 6),
            "arm_a_candidate_mae_log2": round(mae_a3, 6),
            "arm_b_ols_mae_s_eps": round(mae_b2, 6),
            "arm_b_candidate_mae_s_eps": round(mae_b_cand, 6),
        }

    return {
        "status": "PASS",
        "raw_local_data_audited": raw_checked,
        "dev_codons": len(dev_codons),
        "holdout_codons": len(hold_codons),
        "arm_a_verdict": ha["gate_verdict"],
        "arm_b_verdict": hb["gate_verdict"],
        "recomputed": recomputed,
    }


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--data-dir", type=Path, default=None)
    parser.add_argument(
        "--check",
        action="store_true",
        default=True,
        help="Run verification checks (default; accepted for CLI parity).",
    )
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args(argv)

    rep = run_independent_audit(repo_root=REPO_ROOT, data_dir=args.data_dir) if args.check else {"status": "SKIP"}
    if args.json:
        print(json.dumps(rep, indent=2, allow_nan=False))
    else:
        print(
            f"PASS independent_audit cftr_dms_canx_selectivity: raw_audited={rep['raw_local_data_audited']} "
            f"dev_codons={rep['dev_codons']} holdout_codons={rep['holdout_codons']} "
            f"arm_a={rep['arm_a_verdict']} arm_b={rep['arm_b_verdict']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
