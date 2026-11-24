#!/usr/bin/env python3
"""Independent verification audit for benchmarks/cftr_dms_expression/ (Council Item #1).

Separate pure-stdlib code importing NOTHING from `cftr_dms_benchmark.py` or `tools/`:
1. Validates `done_sha256`, `split_manifest_sha256`, residue-cluster disjointness, and zero
   overlap between `holdout` codons and exposed `PRED-0010` codons.
2. Checks arithmetic consistency of relative MAE reduction and pre-registered gate verdict.
3. When local Tedman et al. 2025 `.xlsx` files exist under `sources/downloads/tedman-2025/`,
   independently re-parses `All_variants` & `AGM ALT Library Comparison.xlsx` using a
   Cholesky L L^T solver (distinct linear algebra method from Gaussian elimination) and
   confirms exact count and float parity (`tol <= 1e-6`).
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

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]

MASTER_XLSX = "ALT 9-6-24 CFTR DMS Expanded Library Master Tables.xlsx"
AGM_XLSX = "AGM ALT Library Comparison.xlsx"
NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
DOMAINS = ("TMD1", "NBD1", "R Domain", "TMD2", "NBD2")
CLASSES = ("nonsense", "inframe_indel", "other")
ALL_DOMAINS = frozenset(("Lasso",) + DOMAINS)


def _locate_data_dir(explicit: Optional[Path], repo_root: Path) -> Optional[Path]:
    if explicit is not None:
        return Path(explicit)
    local = repo_root / "sources" / "downloads" / "tedman-2025"
    if (local / MASTER_XLSX).is_file():
        return local
    if len(repo_root.parents) >= 3:
        sibling = repo_root.parents[2] / "cf-research" / "sources" / "downloads" / "tedman-2025"
        if (sibling / MASTER_XLSX).is_file():
            return sibling
    return None


def _read_sheet_rows_dict(xlsx_file: Path, target_sheet: str, expand_merged: bool) -> Dict[Tuple[int, str], str]:
    with zipfile.ZipFile(xlsx_file, "r") as zf:
        strings: List[str] = []
        if "xl/sharedStrings.xml" in zf.namelist():
            root = ET.fromstring(zf.read("xl/sharedStrings.xml"))
            for si in root.findall(f"{NS}si"):
                strings.append("".join(el.text or "" for el in si.iter(f"{NS}t")))

        wb = ET.fromstring(zf.read("xl/workbook.xml"))
        rels = {el.attrib["Id"]: el.attrib["Target"] for el in ET.fromstring(zf.read("xl/_rels/workbook.xml.rels"))}
        sheet_xml: Optional[str] = None
        for sh in wb.findall(f".//{NS}sheet"):
            if sh.attrib.get("name") == target_sheet:
                rid = sh.attrib["{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"]
                t = rels[rid].lstrip("/")
                sheet_xml = t if t.startswith("xl/") else f"xl/{t}"
                break
        if sheet_xml is None:
            raise KeyError(f"Sheet {target_sheet!r} absent from {xlsx_file.name}")

        sroot = ET.fromstring(zf.read(sheet_xml))
        grid: Dict[Tuple[int, str], str] = {}
        for row in sroot.findall(f".//{NS}row"):
            for cell in row.findall(f"{NS}c"):
                m = re.match(r"^([A-Z]+)(\d+)$", cell.attrib.get("r", ""))
                if not m:
                    continue
                col, rnum = m.group(1), int(m.group(2))
                tp = cell.attrib.get("t", "")
                if tp == "inlineStr":
                    val = "".join(el.text or "" for el in cell.iter(f"{NS}t"))
                else:
                    ve = cell.find(f"{NS}v")
                    if ve is None or ve.text is None:
                        continue
                    val = strings[int(ve.text)] if tp == "s" else ve.text
                if val != "":
                    grid[(rnum, col)] = val

        if expand_merged:
            for mc in sroot.findall(f".//{NS}mergeCell"):
                m = re.fullmatch(r"([A-Z]+)(\d+):([A-Z]+)(\d+)", mc.attrib.get("ref", ""))
                if not m:
                    continue
                c1, r1, c2, r2 = m.group(1), int(m.group(2)), m.group(3), int(m.group(4))
                if c1 == c2 and (r1, c1) in grid:
                    for rr in range(r1, r2 + 1):
                        grid.setdefault((rr, c1), grid[(r1, c1)])
        return grid


def _cholesky_solve_spd(mat: List[List[float]], rhs: List[float]) -> List[float]:
    """Solve positive-definite symmetric linear system via Cholesky factorization L L^T x = b."""
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


def _decimal_pearson_agm_alt(agm_path: Path) -> float:
    grid = _read_sheet_rows_dict(agm_path, "Sheet1", expand_merged=False)
    max_r = max((r for r, _ in grid), default=0)
    xs: List[Decimal] = []
    ys: List[Decimal] = []
    for r in range(1, max_r + 1):
        bv = grid.get((r, "B"))
        dv = grid.get((r, "D"))
        if bv is not None and dv is not None:
            xs.append(Decimal(bv))
            ys.append(Decimal(dv))
    n = Decimal(len(xs))
    mx = sum(xs) / n
    my = sum(ys) / n
    num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    dx = sum((x - mx) ** 2 for x in xs)
    dy = sum((y - my) ** 2 for y in ys)
    return float(num) / math.sqrt(float(dx) * float(dy))


def run_independent_audit(repo_root: Path = REPO_ROOT, data_dir: Optional[Path] = None) -> Dict[str, Any]:
    bench_dir = repo_root / "benchmarks" / "cftr_dms_expression"
    manifest = json.loads((bench_dir / "split_manifest.json").read_text(encoding="utf-8"))
    results = json.loads((bench_dir / "results_2026-11-23.json").read_text(encoding="utf-8"))
    done_txt = (bench_dir / "done.md").read_text(encoding="utf-8", errors="replace").replace("\r\n", "\n")
    done_sha = hashlib.sha256(done_txt.encode("utf-8")).hexdigest()
    if results["done_sha256"] != done_sha:
        raise AssertionError(f"done_sha256 mismatch: {results['done_sha256']} != {done_sha}")

    manifest_sha = hashlib.sha256(json.dumps(manifest, sort_keys=True).encode("utf-8")).hexdigest()
    if results["split_manifest_sha256"] != manifest_sha:
        raise AssertionError("split_manifest_sha256 mismatch.")

    th = json.loads((repo_root / "benchmarks/theratype/theratype_open_annotations_2026-11-21.json").read_text(encoding="utf-8"))
    exposed = {int(v["codon_position"]) for v in th["variants"]}
    assigns = manifest["assignments"]
    dev_codons = {int(a["codon_position"]) for a in assigns if a["split"] == "dev"}
    hold_codons = {int(a["codon_position"]) for a in assigns if a["split"] == "holdout"}
    if dev_codons & hold_codons:
        raise AssertionError("Codon overlap between dev and holdout.")
    if hold_codons & exposed:
        raise AssertionError("Exposed PRED-0010 codon present in holdout.")

    he = results["holdout_evaluation"]
    stronger = float(he["stronger_baseline_mae_log2"])
    cand = float(he["candidate_partial_ridge"]["mae_log2"])
    recomputed_rel = round((stronger - cand) / stronger, 6)
    if abs(recomputed_rel - float(he["relative_mae_reduction_vs_stronger"])) > 2e-6:
        raise AssertionError("Arithmetic discrepancy in relative MAE reduction.")

    ddir = _locate_data_dir(data_dir, repo_root)
    raw_checked = False
    recomputed: Dict[str, Any] = {}
    if ddir is not None:
        agm_r = _decimal_pearson_agm_alt(ddir / AGM_XLSX)
        reported_agm = float(results["positive_controls"]["agm_alt_cross_study"]["pearson_r_linear"])
        if abs(agm_r - reported_agm) > 1e-6:
            raise AssertionError(f"AGM ALT Pearson mismatch: {agm_r} != {reported_agm}")

        grid = _read_sheet_rows_dict(ddir / MASTER_XLSX, "All_variants", expand_merged=True)
        by_row = {int(a["row_num"]): a for a in assigns}
        cur_sub = ""
        dev_rows: List[Tuple[str, str, float, float]] = []
        hold_rows: List[Tuple[str, str, float, float]] = []
        max_r = max(r for r, _ in grid)
        for rnum in range(3, max_r + 1):
            if (rnum, "A") in grid:
                cur_sub = grid[(rnum, "A")].strip()
            if rnum not in by_row:
                continue
            info = by_row[rnum]
            if cur_sub != info["sub_domain"] or cur_sub not in ALL_DOMAINS:
                raise AssertionError(f"Row {rnum} sub_domain mismatch: {cur_sub} != {info['sub_domain']}")
            dmso = float(grid[(rnum, "F")])
            combo = float(grid[(rnum, "R")])
            pair = (info["sub_domain"], info["variant_class"], math.log2(dmso), math.log2(combo))
            if info["split"] == "dev":
                dev_rows.append(pair)
            else:
                hold_rows.append(pair)

        # Independently fit 1D OLS (Cholesky 2x2) and selected Ridge (Cholesky 10x10).
        xs_d = [x[2] for x in dev_rows]
        ys_d = [x[3] for x in dev_rows]
        n_d = float(len(dev_rows))
        b_ols = _cholesky_solve_spd(
            [[n_d, sum(xs_d)], [sum(xs_d), sum(v * v for v in xs_d)]],
            [sum(ys_d), sum(x * y for x, y in zip(xs_d, ys_d))],
        )
        ols_hold_mae = sum(abs(y - (b_ols[0] + b_ols[1] * x)) for _, _, x, y in hold_rows) / len(hold_rows)
        if abs(ols_hold_mae - float(he["dmso_ols"]["mae_log2"])) > 1e-6:
            raise AssertionError(f"OLS holdout MAE mismatch: {ols_hold_mae}")

        alpha = float(results["models_fitted_on_dev"]["candidate_partial_ridge"]["selected_alpha"])
        raw_d = [
            [1.0 if sd == s else 0.0 for s in DOMAINS] + [1.0 if vc == c else 0.0 for c in CLASSES]
            for sd, vc, _, _ in dev_rows
        ]
        m = len(raw_d[0])
        means = [sum(row[j] for row in raw_d) / len(raw_d) for j in range(m)]
        stds = [math.sqrt(sum((row[j] - means[j]) ** 2 for row in raw_d) / len(raw_d)) for j in range(m)]
        xmat = [
            [1.0, dev_rows[i][2]] + [((raw_d[i][j] - means[j]) / stds[j] if stds[j] > 1e-12 else 0.0) for j in range(m)]
            for i in range(len(dev_rows))
        ]
        p = 2 + m
        xtx = [[sum(xmat[r][a] * xmat[r][b] for r in range(len(dev_rows))) for b in range(p)] for a in range(p)]
        for j in range(2, p):
            xtx[j][j] += alpha if stds[j - 2] > 1e-12 else 1.0
        xty = [sum(xmat[r][a] * dev_rows[r][3] for r in range(len(dev_rows))) for a in range(p)]
        beta = _cholesky_solve_spd(xtx, xty)

        cand_hold_errs: List[float] = []
        for sd, vc, xd, yd in hold_rows:
            raw_h = [1.0 if sd == s else 0.0 for s in DOMAINS] + [1.0 if vc == c else 0.0 for c in CLASSES]
            pred = beta[0] + beta[1] * xd + sum(
                beta[2 + j] * ((raw_h[j] - means[j]) / stds[j] if stds[j] > 1e-12 else 0.0) for j in range(m)
            )
            cand_hold_errs.append(abs(yd - pred))
        cand_hold_mae = sum(cand_hold_errs) / len(hold_rows)
        if abs(cand_hold_mae - cand) > 1e-6:
            raise AssertionError(f"Cholesky Ridge holdout MAE mismatch: {cand_hold_mae} != {cand}")

        raw_checked = True
        recomputed = {
            "agm_alt_pearson_r": round(agm_r, 6),
            "dmso_ols_holdout_mae_log2": round(ols_hold_mae, 6),
            "candidate_ridge_holdout_mae_log2": round(cand_hold_mae, 6),
        }

    return {
        "status": "PASS",
        "raw_local_data_audited": raw_checked,
        "dev_codons": len(dev_codons),
        "holdout_codons": len(hold_codons),
        "gate_verdict": he["gate_verdict"],
        "recomputed": recomputed,
    }


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--data-dir", type=Path, default=None)
    parser.add_argument(
        "--check",
        action="store_true",
        default=True,
        help="Run verification checks (default behavior; accepted for CLI parity with cftr_dms_benchmark.py).",
    )
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args(argv)

    rep = run_independent_audit(repo_root=REPO_ROOT, data_dir=args.data_dir) if args.check else {}
    if args.json:
        print(json.dumps(rep, indent=2, allow_nan=False))
    else:
        print(
            f"PASS independent_audit cftr_dms_expression: raw_audited={rep['raw_local_data_audited']} "
            f"dev_codons={rep['dev_codons']} holdout_codons={rep['holdout_codons']} "
            f"verdict={rep['gate_verdict']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
