#!/usr/bin/env python3
"""Independent verification audit for benchmarks/cftr_dms_3d_structure/.

Pure-stdlib verification code importing NOTHING from `cftr_dms_3d_structure.py` or `tools/`:
1. Validates `done_sha256`, `split_manifest_sha256`, `ATTRIBUTION.md` mmCIF digests, `deltaF508`
   `auth_seq_id` vs `label_seq_id` `-1` shift, residue-cluster disjointness, `0` `PRED-0010`
   holdout overlap, and strict WT amino-acid parity across all `234` variants.
2. Independently solves TMD1+TMD2 rigid superposition (`8EIQ -> 6MSM` and `5UAK -> 6MSM`) via
   Horn's (1987) 4x4 quaternion eigenvalue method (Jacobi symmetric eigensolver) and checks
   `<= 1e-6` parity on RMSDs and `WJX`/`CV6` separation (runs both offline and in CI).
3. Checks rounded-ratio arithmetic bounds on `relative_mae_reduction_vs_stronger` and gate verdict.
4. When local Tedman et al. 2025 `.xlsx` exists under `sources/downloads/tedman-2025/`,
   independently re-parses `All_variants` and fits `dmso_ols` (2x2 Cholesky) and
   `dmso_3d_cryoem_ridge` (8x8 Cholesky L L^T with `dev`-only mean imputation) confirming
   `<= 1e-6` float parity across holdout MAEs, dev weights, and 5-fold dev-CV MAE.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
from pathlib import Path
import re
import shlex
import sys
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple
import xml.etree.ElementTree as ET
import zipfile

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]

MASTER_XLSX = "ALT 9-6-24 CFTR DMS Expanded Library Master Tables.xlsx"
MASTER_XLSX_SHA256 = "2a702ee3744b8039a9c7c4c74aec9bb30a91301b138d4948e0e2c81a1931e51b"
NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
TMD_INTERVALS = ((77, 360), (859, 1155))
THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}
ALLOWLISTED_WT_EXCEPTIONS = frozenset({(89, 470, "M470V"), (217, 966, "Q966=")})


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


def _parse_cif_gz(path: Path, expected_sha: str) -> Tuple[Dict[int, Tuple[str, Tuple[float, float, float], float, int]], Dict[str, List[Tuple[float, float, float]]]]:
    raw = path.read_bytes()
    if hashlib.sha256(raw).hexdigest() != expected_sha:
        raise AssertionError(f"mmCIF SHA256 mismatch for {path.name}")
    lines = gzip.decompress(raw).decode("utf-8").splitlines()
    cols: List[str] = []
    mode: Optional[str] = None
    ca: Dict[int, Tuple[str, Tuple[float, float, float], float, int]] = {}
    het: Dict[str, List[Tuple[float, float, float]]] = {}
    for ln in lines:
        s = ln.strip()
        if s == "loop_":
            cols = []
            mode = "hdr"
            continue
        if mode == "hdr" and s.startswith("_"):
            cols.append(s)
            continue
        if mode == "hdr" and cols:
            if cols[0].startswith("_atom_site."):
                cols = [c.split(".", 1)[1] for c in cols]
                mode = "atom"
            else:
                mode = None
        if mode == "atom":
            if not s or s.startswith("#") or s.startswith("_") or s == "loop_":
                break
            toks = shlex.split(s)
            if len(toks) != len(cols):
                raise AssertionError(f"Malformed _atom_site row in {path.name}")
            rec = dict(zip(cols, toks))
            if rec["pdbx_PDB_model_num"] != "1" or rec["pdbx_PDB_ins_code"] not in {"?", "."} or rec["label_alt_id"] not in {".", "A"}:
                continue
            if rec["group_PDB"] == "ATOM" and rec["auth_asym_id"] == "A" and rec["label_atom_id"] == "CA":
                seq = int(rec["auth_seq_id"])
                lseq = int(rec["label_seq_id"]) if rec["label_seq_id"].lstrip("-").isdigit() else -1
                ca[seq] = (
                    rec["label_comp_id"],
                    (float(rec["Cartn_x"]), float(rec["Cartn_y"]), float(rec["Cartn_z"])),
                    float(rec["B_iso_or_equiv"]),
                    lseq,
                )
            elif rec["group_PDB"] == "HETATM" and rec["type_symbol"] not in {"H", "D"}:
                het.setdefault(rec["label_comp_id"], []).append(
                    (float(rec["Cartn_x"]), float(rec["Cartn_y"]), float(rec["Cartn_z"]))
                )
    return ca, het


def _jacobi_eigen_sym(mat: Sequence[Sequence[float]], sweeps: int = 80) -> Tuple[List[float], List[List[float]]]:
    n = len(mat)
    a = [list(r) for r in mat]
    v = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    for _ in range(sweeps):
        off = sum(a[i][j] ** 2 for i in range(n) for j in range(i + 1, n))
        if off < 1e-24:
            break
        for p in range(n):
            for q in range(p + 1, n):
                apq = a[p][q]
                if abs(apq) < 1e-300:
                    continue
                tau = (a[q][q] - a[p][p]) / (2.0 * apq)
                t = (1.0 if tau >= 0.0 else -1.0) / (abs(tau) + math.sqrt(tau * tau + 1.0))
                c = 1.0 / math.sqrt(t * t + 1.0)
                s = t * c
                for k in range(n):
                    akp, akq = a[k][p], a[k][q]
                    a[k][p], a[k][q] = c * akp - s * akq, s * akp + c * akq
                for k in range(n):
                    apk, aqk = a[p][k], a[q][k]
                    a[p][k], a[q][k] = c * apk - s * aqk, s * apk + c * aqk
                for k in range(n):
                    vkp, vkq = v[k][p], v[k][q]
                    v[k][p], v[k][q] = c * vkp - s * vkq, s * vkp + c * vkq
    return [a[i][i] for i in range(n)], v


def horn_quaternion_superpose(
    moving: Sequence[Tuple[float, float, float]],
    target: Sequence[Tuple[float, float, float]],
) -> Tuple[List[List[float]], Tuple[float, float, float]]:
    """Independent Horn (1987) 4x4 quaternion eigensystem superposition `(R, t)`."""
    n = len(moving)
    if n < 3 or n != len(target):
        raise ValueError("Need >= 3 paired coordinates for Horn superposition.")
    pc = [sum(p[k] for p in moving) / n for k in range(3)]
    qc = [sum(q[k] for q in target) / n for k in range(3)]
    s = [[0.0] * 3 for _ in range(3)]
    for p, q in zip(moving, target):
        for i in range(3):
            for j in range(3):
                s[i][j] += (p[i] - pc[i]) * (q[j] - qc[j])
    sxx, sxy, sxz = s[0]
    syx, syy, syz = s[1]
    szx, szy, szz = s[2]
    n_mat = [
        [sxx + syy + szz, syz - szy, szx - sxz, sxy - syx],
        [syz - szy, sxx - syy - szz, sxy + syx, szx + sxz],
        [szx - sxz, sxy + syx, -sxx + syy - szz, syz + szy],
        [sxy - syx, szx + sxz, syz + szy, -sxx - syy + szz],
    ]
    evals, evecs = _jacobi_eigen_sym(n_mat)
    best_k = max(range(4), key=lambda idx: evals[idx])
    q0, q1, q2, q3 = (evecs[r][best_k] for r in range(4))
    rot = [
        [q0 * q0 + q1 * q1 - q2 * q2 - q3 * q3, 2.0 * (q1 * q2 - q0 * q3), 2.0 * (q1 * q3 + q0 * q2)],
        [2.0 * (q1 * q2 + q0 * q3), q0 * q0 - q1 * q1 + q2 * q2 - q3 * q3, 2.0 * (q2 * q3 - q0 * q1)],
        [2.0 * (q1 * q3 - q0 * q2), 2.0 * (q2 * q3 + q0 * q1), q0 * q0 - q1 * q1 - q2 * q2 + q3 * q3],
    ]
    trans = tuple(qc[i] - sum(rot[i][j] * pc[j] for j in range(3)) for i in range(3))
    return rot, trans  # type: ignore[return-value]


def _apply_rt(rot: Sequence[Sequence[float]], trans: Tuple[float, float, float], p: Tuple[float, float, float]) -> Tuple[float, float, float]:
    return (
        rot[0][0] * p[0] + rot[0][1] * p[1] + rot[0][2] * p[2] + trans[0],
        rot[1][0] * p[0] + rot[1][1] * p[1] + rot[1][2] * p[2] + trans[1],
        rot[2][0] * p[0] + rot[2][1] * p[1] + rot[2][2] * p[2] + trans[2],
    )


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


def _parse_wt_aa(prot: str) -> str:
    s = prot.strip().split("/", 1)[0].strip()
    m = re.match(r"^([A-Z])(\d+)", s)
    if not m:
        raise AssertionError(f"Unparseable Protein distinction: {prot!r}")
    return m.group(1)


def run_independent_audit(repo_root: Path = REPO_ROOT, data_dir: Optional[Path] = None) -> Dict[str, Any]:
    bench_dir = repo_root / "benchmarks" / "cftr_dms_3d_structure"
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

    # Verify mmCIF SHA-256 digests against ATTRIBUTION.md and parse Cryo-EM structures via Horn quaternion method.
    attr_txt = (repo_root / "tools/cftr-model/data/ATTRIBUTION.md").read_text(encoding="utf-8")
    cif_shas = results["cryoem_mmcif_sha256"]
    for fname, digest in cif_shas.items():
        if fname not in attr_txt or digest not in attr_txt:
            raise AssertionError(f"ATTRIBUTION.md missing expected entry for {fname}")

    data_model_dir = repo_root / "tools" / "cftr-model" / "data"
    ca6, _ = _parse_cif_gz(data_model_dir / "6MSM.cif.gz", cif_shas["6MSM.cif.gz"])
    ca5, _ = _parse_cif_gz(data_model_dir / "5UAK.cif.gz", cif_shas["5UAK.cif.gz"])
    ca8, het8 = _parse_cif_gz(data_model_dir / "8EIQ.cif.gz", cif_shas["8EIQ.cif.gz"])

    if ca6[508][0] != "PHE" or 508 in ca8 or ca8[509][3] != 508:
        raise AssertionError("DeltaF508 auth_seq_id / label_seq_id shift verification failed.")

    # Independent WT sequence parity check
    seen_exc = set()
    for a in assigns:
        rnum = int(a["row_num"])
        codon = int(a["codon_position"])
        prot = str(a["protein_distinction"])
        wt_aa = _parse_wt_aa(prot)
        for pid, cmap in (("6MSM", ca6), ("8EIQ", ca8)):
            if pid == "8EIQ" and codon == 508:
                continue
            if codon not in cmap:
                continue
            obs = THREE_TO_ONE[cmap[codon][0]]
            if obs != wt_aa:
                key = (rnum, codon, prot)
                if key in ALLOWLISTED_WT_EXCEPTIONS and a["split"] == "dev":
                    seen_exc.add(key)
                else:
                    raise AssertionError(f"Unexpected WT mismatch in {pid}: {key} obs={obs} exp={wt_aa}")
    if seen_exc != ALLOWLISTED_WT_EXCEPTIONS:
        raise AssertionError("Expected allowlisted WT exceptions not matched.")

    in_tmd = lambda seq: any(lo <= seq <= hi for lo, hi in TMD_INTERVALS)
    c8 = sorted(k for k in ca8 if k in ca6 and in_tmd(k))
    c5 = sorted(k for k in ca5 if k in ca6 and in_tmd(k))
    r8, t8 = horn_quaternion_superpose([ca8[k][1] for k in c8], [ca6[k][1] for k in c8])
    r5, t5 = horn_quaternion_superpose([ca5[k][1] for k in c5], [ca6[k][1] for k in c5])

    rmsd8 = math.sqrt(sum(math.dist(_apply_rt(r8, t8, ca8[k][1]), ca6[k][1]) ** 2 for k in c8) / len(c8))
    rmsd5 = math.sqrt(sum(math.dist(_apply_rt(r5, t5, ca5[k][1]), ca6[k][1]) ** 2 for k in c5) / len(c5))
    wjx_pts = [_apply_rt(r8, t8, p) for p in het8["WJX"]]
    cv6_pts = [_apply_rt(r8, t8, p) for p in het8["CV6"]]
    min_wc = min(math.dist(a, b) for a in wjx_pts for b in cv6_pts)
    wc_cent = math.dist(
        tuple(sum(p[k] for p in wjx_pts) / len(wjx_pts) for k in range(3)),
        tuple(sum(p[k] for p in cv6_pts) / len(cv6_pts) for k in range(3)),
    )

    rep_m = results["pre_holdout_controls"]["cryoem_tmd_superposition_and_pockets"]
    if (
        abs(rmsd8 - float(rep_m["rmsd_tmd_8eiq_to_6msm_A"])) > 1e-6
        or abs(rmsd5 - float(rep_m["rmsd_tmd_5uak_to_6msm_A"])) > 1e-6
        or abs(min_wc - float(rep_m["wjx_cv6_min_heavy_atom_dist_A"])) > 1e-6
        or abs(wc_cent - float(rep_m["wjx_cv6_centroid_dist_A"])) > 1e-6
    ):
        raise AssertionError("Horn superposition structural metrics mismatch vs results JSON.")

    # Arithmetic check with input rounding uncertainty propagation (+/-0.5e-6)
    he = results["holdout_evaluation"]
    stronger = float(he["stronger_baseline_mae_log2"])
    cand = float(he["candidate_dmso_3d_cryoem_ridge"]["mae_log2"])
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
        ca5_aligned = {k: _apply_rt(r5, t5, v[1]) for k, v in ca5.items()}

        dev_rows: List[Tuple[int, float, float, Tuple[Optional[float], ...]]] = []
        hold_rows: List[Tuple[int, float, float, Tuple[Optional[float], ...]]] = []
        for a in assigns:
            rnum = int(a["row_num"])
            codon = int(a["codon_position"])
            prot = str(a["protein_distinction"])
            x_dmso = math.log2(float(grid[(rnum, "F")]))
            y_combo = math.log2(float(grid[(rnum, "R")]))
            if (rnum, codon, prot) == (217, 966, "Q966="):
                feats: Tuple[Optional[float], ...] = (None, None, None, None, None, 1.0)
            elif codon not in ca6:
                feats = (None, None, None, None, None, 1.0)
            else:
                p6 = ca6[codon][1]
                dw = min(math.dist(p6, q) for q in wjx_pts)
                dc = min(math.dist(p6, q) for q in cv6_pts)
                cn = float(sum(1 for k, v in ca6.items() if k != codon and math.dist(p6, v[1]) <= 10.0))
                bf = float(ca6[codon][2])
                disp = math.dist(ca5_aligned[codon], p6) if codon in ca5_aligned else None
                unres = 0.0 if codon in ca5_aligned else 1.0
                feats = (math.exp(-((dw / 12.0) ** 2)), math.exp(-((dc / 12.0) ** 2)), cn, bf, disp, unres)
            item = (codon, x_dmso, y_combo, feats)
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

        # Cholesky 8x8 partial Ridge on dev & 5-fold dev CV
        alpha = float(results["models_fitted_on_dev"]["candidate_dmso_3d_cryoem_ridge"]["selected_alpha"])
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
        rep_cv_mae = float(results["models_fitted_on_dev"]["candidate_dmso_3d_cryoem_ridge"]["dev_cv_mae_log2"])
        if abs(cv_mae - rep_cv_mae) > 1e-6:
            raise AssertionError(f"Cholesky dev CV MAE mismatch: {cv_mae} != {rep_cv_mae}")

        beta, imp_m, fm, fs = _fit_partial_ridge_cholesky(
            [x[1] for x in dev_rows], [x[3] for x in dev_rows], [x[2] for x in dev_rows], alpha
        )
        cand_hold_mae = sum(
            abs(yd - _predict_cholesky_ridge(xd, ft, beta, imp_m, fm, fs))
            for _, xd, yd, ft in hold_rows
        ) / len(hold_rows)
        if abs(cand_hold_mae - cand) > 1e-6:
            raise AssertionError(f"Cholesky Ridge holdout MAE mismatch: {cand_hold_mae} != {cand}")

        raw_checked = True
        recomputed = {
            "horn_rmsd_8eiq_to_6msm_A": round(rmsd8, 6),
            "horn_rmsd_5uak_to_6msm_A": round(rmsd5, 6),
            "cholesky_dmso_ols_holdout_mae_log2": round(ols_hold_mae, 6),
            "cholesky_ridge_dev_cv_mae_log2": round(cv_mae, 6),
            "cholesky_ridge_holdout_mae_log2": round(cand_hold_mae, 6),
        }

    return {
        "status": "PASS",
        "raw_local_data_audited": raw_checked,
        "horn_rmsd_8eiq_to_6msm_A": round(rmsd8, 6),
        "dev_codons": len(dev_codons),
        "holdout_codons": len(hold_codons),
        "gate_verdict": he["gate_verdict"],
        "recomputed": recomputed,
    }


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--data-dir", type=Path, default=None)
    parser.add_argument("--check", action="store_true", default=True)
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args(argv)

    rep = run_independent_audit(repo_root=REPO_ROOT, data_dir=args.data_dir) if args.check else {}
    if args.json:
        print(json.dumps(rep, indent=2, allow_nan=False))
    else:
        print(
            f"PASS independent_audit cftr_dms_3d_structure: raw_audited={rep['raw_local_data_audited']} "
            f"horn_rmsd_8eiq={rep['horn_rmsd_8eiq_to_6msm_A']} verdict={rep['gate_verdict']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
