#!/usr/bin/env python3
"""Independent verification audit for benchmarks/cftr_dms_substitution_biophysics/.

Pure-stdlib verification code importing NOTHING from `cftr_dms_substitution_biophysics.py` or `tools/`:
1. Validates `done_sha256`, `split_manifest_sha256`, `p13569_sha256`, `ATTRIBUTION.md` digest entry,
   residue-cluster disjointness, `0` `PRED-0010` holdout overlap, and strict WT amino-acid parity
   across all parsed components of all `234` variants (including compound `Q359K/T360K`).
2. Independently reconstructs the 20x20 Henikoff & Henikoff (1992) BLOSUM62 matrix from its
   lower-triangular representation, verifies symmetry and published anchors, verifies Kyte & Doolittle
   (1982) hydropathy and charge rules, and tests feature extraction across all `234` manifest rows (CI + local).
3. Checks rounded-ratio arithmetic bounds (`+/-0.5e-6` input uncertainty propagation) on
   `relative_mae_reduction_vs_stronger` and gate verdict.
4. When local Tedman et al. 2025 `.xlsx` exists under `sources/downloads/tedman-2025/`,
   independently re-parses `All_variants` and solves `dmso_ols` (2x2 Cholesky L L^T),
   `dmso_subst_biophysics_ridge` (9x9 Cholesky L L^T with `dev`-only mean imputation), and
   `dev` VIFs (Cholesky correlation inversion), confirming `<= 1e-6` float parity across holdout MAEs,
   dev weights, 5-fold dev-CV MAE, and VIFs.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import sys
from typing import Any, Dict, List, Mapping, Optional, Sequence, Set, Tuple
import xml.etree.ElementTree as ET
import zipfile

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]

MASTER_XLSX = "ALT 9-6-24 CFTR DMS Expanded Library Master Tables.xlsx"
MASTER_XLSX_SHA256 = "2a702ee3744b8039a9c7c4c74aec9bb30a91301b138d4948e0e2c81a1931e51b"
NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
ALLOWLISTED_WT_EXCEPTIONS = frozenset({(89, 470, "M470V"), (217, 966, "Q966=")})

# Lower-triangular BLOSUM62 (Henikoff & Henikoff 1992, PNAS 89:10915) along A R N D C Q E G H I L K M F P S T W Y V
_AA_AXIS = "ARNDCQEGHILKMFPSTWYV"
_BLOSUM62_LOWER_TRI = (
    (4,),
    (-1, 5),
    (-2, 0, 6),
    (-2, -2, 1, 6),
    (0, -3, -3, -3, 9),
    (-1, 1, 0, 0, -3, 5),
    (-1, 0, 0, 2, -4, 2, 5),
    (0, -2, 0, -1, -3, -2, -2, 6),
    (-2, 0, 1, -1, -3, 0, 0, -2, 8),
    (-1, -3, -3, -3, -1, -3, -3, -4, -3, 4),
    (-1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4),
    (-1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5),
    (-1, -1, -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5),
    (-2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3, 0, 6),
    (-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7),
    (1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -2, 0, -1, -2, -1, 4),
    (0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5),
    (-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3, -2, 11),
    (-2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2, 2, 7),
    (0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0, -3, -1, 4),
)


def _build_blosum62() -> Dict[Tuple[str, str], int]:
    mat: Dict[Tuple[str, str], int] = {}
    for i, row in enumerate(_BLOSUM62_LOWER_TRI):
        if len(row) != i + 1:
            raise AssertionError("Malformed lower-triangular BLOSUM62 row.")
        aa_i = _AA_AXIS[i]
        for j, score in enumerate(row):
            aa_j = _AA_AXIS[j]
            mat[(aa_i, aa_j)] = score
            mat[(aa_j, aa_i)] = score
    return mat


_B62 = _build_blosum62()
_KD = {
    "I": 4.5, "V": 4.2, "L": 3.8, "F": 2.8, "C": 2.5,
    "M": 1.9, "A": 1.8, "G": -0.4, "T": -0.7, "S": -0.8,
    "W": -0.9, "Y": -1.3, "P": -1.6, "H": -3.2, "Q": -3.5,
    "E": -3.5, "D": -3.5, "N": -3.5, "K": -3.9, "R": -4.5,
}
_CHG = {aa: (1.0 if aa in "KR" else (-1.0 if aa in "DE" else 0.0)) for aa in _AA_AXIS}


def _locate_data_dir(explicit: Optional[Path], repo_root: Path) -> Optional[Path]:
    if explicit is not None:
        p = Path(explicit)
        return p if (p / MASTER_XLSX).is_file() else None
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


def cholesky_solve_spd(mat: Sequence[Sequence[float]], rhs: Sequence[float]) -> List[float]:
    """Solve symmetric positive-definite linear system via Cholesky factorization L L^T x = b."""
    n = len(rhs)
    chol = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1):
            acc = mat[i][j] - sum(chol[i][k] * chol[j][k] for k in range(j))
            if i == j:
                if acc <= 0.0:
                    raise ValueError("Matrix is not positive definite.")
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


def _decompose_variant(prot: str) -> List[Tuple[str, str, int, Optional[str]]]:
    items: List[Tuple[str, str, int, Optional[str]]] = []
    for seg in [p.strip() for p in prot.strip().split("/") if p.strip()]:
        m_ind = re.fullmatch(r"([A-Z])(\d+)(?:del|DEL|ins|INS)", seg)
        if m_ind:
            items.append(("indel", m_ind.group(1), int(m_ind.group(2)), None))
            continue
        m_stp = re.fullmatch(r"([A-Z])(\d+)[X\*]", seg)
        if m_stp:
            items.append(("stop", m_stp.group(1), int(m_stp.group(2)), None))
            continue
        m_sub = re.fullmatch(r"([A-Z])(\d+)([A-Z=])", seg)
        if m_sub:
            wt, pos, mut = m_sub.group(1), int(m_sub.group(2)), m_sub.group(3)
            if mut == "=" or mut == wt:
                items.append(("syn", wt, pos, wt))
            else:
                items.append(("sub", wt, pos, mut))
            continue
        raise AssertionError(f"Unparseable variant segment {seg!r} in {prot!r}")
    if not items:
        raise AssertionError(f"Empty variant {prot!r}")
    return items


def _compute_7_features(
    prot: str,
    tm_pos: Set[int],
    helix_or_tm_pos: Set[int],
) -> Tuple[Optional[float], ...]:
    parts = _decompose_variant(prot)
    if all(p[0] == "syn" for p in parts):
        return (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    if len(parts) == 1 and parts[0][0] == "stop":
        return (None, None, None, 0.0, 0.0, 1.0, 0.0)
    if len(parts) == 1 and parts[0][0] == "indel":
        return (None, None, None, 0.0, 0.0, 0.0, 1.0)
    if not all(p[0] == "sub" for p in parts):
        raise AssertionError(f"Unsupported mixed variant components in {prot!r}")

    d_blo = 0.0
    d_kd = 0.0
    d_chg = 0.0
    hb = 0.0
    tm_shk = 0.0
    for _, wt, pos, mut in parts:
        assert mut is not None
        d_blo += float(_B62[(wt, mut)] - _B62[(wt, wt)])
        dk = float(_KD[mut] - _KD[wt])
        dc = abs(_CHG[mut] - _CHG[wt])
        d_kd += dk
        d_chg += dc
        if pos in helix_or_tm_pos and (wt in {"P", "G"} or mut in {"P", "G"}):
            hb = 1.0
        if pos in tm_pos and (abs(dk) >= 3.0 or dc >= 1.0):
            tm_shk = 1.0
    return (d_blo, d_kd, d_chg, hb, tm_shk, 0.0, 0.0)


def _fit_partial_ridge_cholesky(
    xs_dmso: Sequence[float],
    raw_feats: Sequence[Sequence[Optional[float]]],
    ys: Sequence[float],
    alpha: float,
) -> Tuple[List[float], List[float], List[float], List[float]]:
    n = len(ys)
    m = len(raw_feats[0])
    imp_means = [
        sum(float(r[j]) for r in raw_feats if r[j] is not None) / sum(1 for r in raw_feats if r[j] is not None)
        for j in range(m)
    ]
    filled = [[imp_means[j] if r[j] is None else float(r[j]) for j in range(m)] for r in raw_feats]
    means = [sum(r[j] for r in filled) / n for j in range(m)]
    stds: List[float] = []
    for j in range(m):
        var = sum((r[j] - means[j]) ** 2 for r in filled) / n
        stds.append(math.sqrt(var) if var > 1e-12 else 0.0)
    p = 2 + m
    xmat = [
        [1.0, xs_dmso[i]] + [((filled[i][j] - means[j]) / stds[j] if stds[j] > 0.0 else 0.0) for j in range(m)]
        for i in range(n)
    ]
    xtx = [[sum(xmat[r][a] * xmat[r][b] for r in range(n)) for b in range(p)] for a in range(p)]
    for j in range(2, p):
        xtx[j][j] += alpha if stds[j - 2] > 0.0 else 1.0
    xty = [sum(xmat[r][a] * ys[r] for r in range(n)) for a in range(p)]
    beta = cholesky_solve_spd(xtx, xty)
    return beta, imp_means, means, stds


def _predict_cholesky_ridge(
    x_dmso: float,
    raw_f: Sequence[Optional[float]],
    beta: Sequence[float],
    imp_means: Sequence[float],
    means: Sequence[float],
    stds: Sequence[float],
) -> float:
    pred = beta[0] + beta[1] * x_dmso
    for j, val in enumerate(raw_f):
        filled = imp_means[j] if val is None else float(val)
        z = (filled - means[j]) / stds[j] if stds[j] > 0.0 else 0.0
        pred += beta[2 + j] * z
    return pred


def run_independent_audit(repo_root: Path = REPO_ROOT, data_dir: Optional[Path] = None) -> Dict[str, Any]:
    bench_dir = repo_root / "benchmarks" / "cftr_dms_substitution_biophysics"
    results = json.loads((bench_dir / "results_2026-11-24.json").read_text(encoding="utf-8"))
    done_txt = (bench_dir / "done.md").read_text(encoding="utf-8", errors="replace").replace("\r\n", "\n")
    done_sha = hashlib.sha256(done_txt.encode("utf-8")).hexdigest()
    if results["done_sha256"] != done_sha:
        raise AssertionError(f"done_sha256 mismatch: {results['done_sha256']} != {done_sha}")

    manifest_path = repo_root / "benchmarks" / "cftr_dms_expression" / "split_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
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
    if len(hold_codons) < 40:
        raise AssertionError("Insufficient holdout codon clusters.")

    # Verify BLOSUM62 symmetry & published anchors
    for a in _AA_AXIS:
        for b in _AA_AXIS:
            if _B62[(a, b)] != _B62[(b, a)]:
                raise AssertionError(f"BLOSUM62 asymmetry at {a},{b}")
    if _B62[("W", "W")] != 11 or _B62[("C", "C")] != 9 or _B62[("H", "G")] != -2:
        raise AssertionError("BLOSUM62 anchor mismatch.")
    if abs(_KD["I"] - 4.5) > 1e-12 or abs(_KD["R"] - (-4.5)) > 1e-12:
        raise AssertionError("Kyte-Doolittle anchor mismatch.")

    # Verify p13569.json digest against ATTRIBUTION.md and parse sequence + topology
    p135_path = repo_root / "tools/cftr-model/data/p13569.json"
    p135_bytes = p135_path.read_bytes()
    p135_sha = hashlib.sha256(p135_bytes).hexdigest()
    if results["p13569_sha256"] != p135_sha:
        raise AssertionError(f"p13569_sha256 mismatch: {results['p13569_sha256']} != {p135_sha}")
    attr_txt = (repo_root / "tools/cftr-model/data/ATTRIBUTION.md").read_text(encoding="utf-8")
    if "p13569.json" not in attr_txt or p135_sha not in attr_txt:
        raise AssertionError("ATTRIBUTION.md missing expected entry for p13569.json")

    p135 = json.loads(p135_bytes.decode("utf-8"))
    seq = str(p135["sequence"]["value"])
    if len(seq) != 1480:
        raise AssertionError(f"Unexpected P13569 sequence length: {len(seq)}")
    tm_pos: Set[int] = set()
    helix_pos: Set[int] = set()
    for feat in p135.get("features", []):
        ft = feat.get("type")
        if ft in {"Transmembrane", "Helix"}:
            s = int(feat["location"]["start"]["value"])
            e = int(feat["location"]["end"]["value"])
            target = tm_pos if ft == "Transmembrane" else helix_pos
            target.update(range(s, e + 1))
    helix_or_tm = helix_pos | tm_pos

    # Independent P13569 WT sequence parity check across all components & CI 7-feature verification
    seen_exc = set()
    feat_by_row: Dict[int, Tuple[Optional[float], ...]] = {}
    for a in assigns:
        rnum = int(a["row_num"])
        prot = str(a["protein_distinction"])
        for _, wt_aa, pos, _ in _decompose_variant(prot):
            obs = seq[pos - 1]
            if obs != wt_aa:
                key = (rnum, pos, prot)
                if key in ALLOWLISTED_WT_EXCEPTIONS and a["split"] == "dev":
                    seen_exc.add(key)
                else:
                    raise AssertionError(f"Unexpected P13569 WT mismatch: {key} obs={obs} exp={wt_aa}")
        feats = _compute_7_features(prot, tm_pos, helix_or_tm)
        if len(feats) != 7:
            raise AssertionError(f"Expected 7 features for {prot}, got {len(feats)}")
        feat_by_row[rnum] = feats
    if seen_exc != ALLOWLISTED_WT_EXCEPTIONS:
        raise AssertionError("Expected allowlisted WT exceptions not matched.")

    # Verify M470V experimental-WT rule explicitly
    m470v_f = feat_by_row[89]
    if m470v_f[0] != -4.0 or abs((m470v_f[1] or 0.0) - 2.3) > 1e-12:
        raise AssertionError(f"M470V experimental-WT biophysics mismatch: {m470v_f}")

    # Arithmetic check with input rounding uncertainty propagation (+/-0.5e-6)
    he = results["holdout_evaluation"]
    stronger = float(he["stronger_baseline_mae_log2"])
    cand = float(he["candidate_dmso_subst_biophysics_ridge"]["mae_log2"])
    rel_rep = float(he["relative_mae_reduction_vs_stronger"])
    rel_lo = ((stronger - 0.5e-6) - (cand + 0.5e-6)) / (stronger - 0.5e-6) - 1e-6
    rel_hi = ((stronger + 0.5e-6) - (cand - 0.5e-6)) / (stronger + 0.5e-6) + 1e-6
    if not (rel_lo <= rel_rep <= rel_hi):
        raise AssertionError(f"Relative MAE reduction arithmetic check failed: {rel_rep} not in [{rel_lo}, {rel_hi}]")

    ci_low_gt_zero = bool(he["residue_clustered_bootstrap"]["ci95_low_gt_zero"])
    if rel_rep >= 0.10 and ci_low_gt_zero:
        expected_standing = "open"
    else:
        expected_standing = "refuted"
    if he["ledger_standing"] != expected_standing:
        raise AssertionError(f"Gate verdict mismatch: {he['ledger_standing']} != {expected_standing}")

    ddir = _locate_data_dir(data_dir, repo_root)
    raw_checked = False
    recomputed: Dict[str, Any] = {}
    if ddir is not None:
        xlsx_path = ddir / MASTER_XLSX
        if hashlib.sha256(xlsx_path.read_bytes()).hexdigest() != MASTER_XLSX_SHA256:
            raise AssertionError("Master xlsx SHA256 mismatch.")
        grid = _read_sheet_rows_dict(xlsx_path, "All_variants", expand_merged=True)

        dev_rows: List[Tuple[int, float, float, Tuple[Optional[float], ...]]] = []
        hold_rows: List[Tuple[int, float, float, Tuple[Optional[float], ...]]] = []
        for a in assigns:
            rnum = int(a["row_num"])
            codon = int(a["codon_position"])
            x_dmso = math.log2(float(grid[(rnum, "F")]))
            y_combo = math.log2(float(grid[(rnum, "R")]))
            item = (codon, x_dmso, y_combo, feat_by_row[rnum])
            if a["split"] == "dev":
                dev_rows.append(item)
            else:
                hold_rows.append(item)

        # Cholesky 2x2 OLS on dev -> check holdout parity
        xs_d = [x[1] for x in dev_rows]
        ys_d = [x[2] for x in dev_rows]
        n_d = float(len(dev_rows))
        b_ols = cholesky_solve_spd(
            [[n_d, sum(xs_d)], [sum(xs_d), sum(v * v for v in xs_d)]],
            [sum(ys_d), sum(x * y for x, y in zip(xs_d, ys_d))],
        )
        ols_hold_mae = sum(abs(y - (b_ols[0] + b_ols[1] * x)) for _, x, y, _ in hold_rows) / len(hold_rows)
        if abs(ols_hold_mae - float(he["dmso_ols"]["mae_log2"])) > 1e-6:
            raise AssertionError(f"Cholesky OLS holdout MAE mismatch: {ols_hold_mae}")

        # Cholesky 9x9 partial Ridge on dev & 5-fold dev CV
        cand_rep = results["models_fitted_on_dev"]["candidate_dmso_subst_biophysics_ridge"]
        alpha = float(cand_rep["selected_alpha"])
        fold_map = {c: idx % 5 for idx, c in enumerate(sorted(dev_codons))}
        cv_errs: List[float] = []
        for fold in range(5):
            tr = [x for x in dev_rows if fold_map[x[0]] != fold]
            va = [x for x in dev_rows if fold_map[x[0]] == fold]
            fb, fimp, fm, fs = _fit_partial_ridge_cholesky(
                [x[1] for x in tr], [x[3] for x in tr], [x[2] for x in tr], alpha
            )
            for _, xd, yd, ft in va:
                cv_errs.append(abs(yd - _predict_cholesky_ridge(xd, ft, fb, fimp, fm, fs)))
        cv_mae = sum(cv_errs) / len(cv_errs)
        rep_cv_mae = float(cand_rep["dev_cv_mae_log2"])
        if abs(cv_mae - rep_cv_mae) > 1e-6:
            raise AssertionError(f"Cholesky dev CV MAE mismatch: {cv_mae} != {rep_cv_mae}")

        beta, imp_m, fm, fs = _fit_partial_ridge_cholesky(
            [x[1] for x in dev_rows], [x[3] for x in dev_rows], [x[2] for x in dev_rows], alpha
        )
        if abs(beta[0] - float(cand_rep["b0_intercept"])) > 1e-6 or abs(beta[1] - float(cand_rep["b1_dmso_slope"])) > 1e-6:
            raise AssertionError("Cholesky Ridge b0/b1 mismatch vs results JSON.")
        for w_chol, w_rep in zip(beta[2:], cand_rep["weights"]):
            if abs(w_chol - float(w_rep)) > 1e-6:
                raise AssertionError(f"Cholesky Ridge weight mismatch: {w_chol} != {w_rep}")

        cand_hold_mae = sum(
            abs(yd - _predict_cholesky_ridge(xd, ft, beta, imp_m, fm, fs))
            for _, xd, yd, ft in hold_rows
        ) / len(hold_rows)
        if abs(cand_hold_mae - cand) > 1e-6:
            raise AssertionError(f"Cholesky Ridge holdout MAE mismatch: {cand_hold_mae} != {cand}")

        # Cholesky correlation matrix VIF parity check on dev
        z_rows = [
            [((imp_m[j] if ft[j] is None else float(ft[j])) - fm[j]) / fs[j] for j in range(7)]
            for _, _, _, ft in dev_rows
        ]
        n_dev_i = len(z_rows)
        corr7 = [[sum(z_rows[r][i] * z_rows[r][j] for r in range(n_dev_i)) / n_dev_i for j in range(7)] for i in range(7)]
        rep_vifs = cand_rep["dev_vif_context"]["vif_features_only"]
        for j, fname in enumerate(cand_rep["feature_names"]):
            rhs = [1.0 if i == j else 0.0 for i in range(7)]
            inv_col = cholesky_solve_spd(corr7, rhs)
            if abs(inv_col[j] - float(rep_vifs[fname])) > 1e-6:
                raise AssertionError(f"Cholesky VIF mismatch for {fname}: {inv_col[j]} != {rep_vifs[fname]}")

        raw_checked = True
        recomputed = {
            "cholesky_dmso_ols_holdout_mae_log2": round(ols_hold_mae, 6),
            "cholesky_ridge_dev_cv_mae_log2": round(cv_mae, 6),
            "cholesky_ridge_holdout_mae_log2": round(cand_hold_mae, 6),
            "cholesky_is_inframe_indel_weight": round(beta[8], 6),
        }

    if not isinstance(he.get("gate_verdict"), str) or not he["gate_verdict"]:
        raise AssertionError("Missing holdout_evaluation.gate_verdict.")
    return {
        "status": "PASS",
        "raw_local_data_audited": raw_checked,
        "p13569_sequence_length": len(seq),
        "dev_codons": len(dev_codons),
        "holdout_codons": len(hold_codons),
        "gate_verdict": he["gate_verdict"],
        "recomputed": recomputed,
    }


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--data-dir", type=Path, default=None)
    parser.add_argument("--check", action="store_true", help="Run the independent verification audit.")
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args(argv)

    rep = run_independent_audit(repo_root=REPO_ROOT, data_dir=args.data_dir)
    if args.json:
        print(json.dumps(rep, indent=2, allow_nan=False))
    else:
        print(
            f"PASS independent_audit cftr_dms_substitution_biophysics: "
            f"raw_audited={rep['raw_local_data_audited']} "
            f"p13569_len={rep['p13569_sequence_length']} verdict={rep['gate_verdict']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
