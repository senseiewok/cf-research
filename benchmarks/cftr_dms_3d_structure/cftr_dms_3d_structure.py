#!/usr/bin/env python3
"""Pure 3D Cryo-EM structure-to-DMS modulator-rescue benchmark on Tedman et al. 2025.

Evaluates whether 6 Cryo-EM structural features (VX-445/WJX & VX-661/CV6 pocket RBFs from 8EIQ,
6MSM local C-alpha packing CN_10A, 6MSM C-alpha B_iso_or_equiv map-resolution proxy, 5UAK->6MSM
conformational displacement, and unresolved-residue flag) improve log2(VX-661+VX-445 Mut/WT)
prediction over pure baseline untreated expression OLS on the locked residue-clustered holdout,
with sub-domain and variant-class indicators excluded from both baseline and candidate.

Pure Python 3.10+ standard library only; zero network calls; raw `.xlsx` stays git-ignored.
Research benchmark only: plasma-membrane expression is not channel conductance or clinical benefit.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
import datetime as dt
import gzip
import hashlib
import json
import math
from pathlib import Path
import random
import re
import shlex
import statistics
import sys
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple
import xml.etree.ElementTree as ET
import zipfile

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]

MASTER_XLSX_NAME = "ALT 9-6-24 CFTR DMS Expanded Library Master Tables.xlsx"
MASTER_XLSX_SHA256 = "2a702ee3744b8039a9c7c4c74aec9bb30a91301b138d4948e0e2c81a1931e51b"
DONE_SHA256_EXPECTED = "a997290fda86f3f5ee69c839df7bf1d2c78074c50497283e564f69b71f60fb62"

SPLIT_MANIFEST_REL = Path("benchmarks/cftr_dms_expression/split_manifest.json")
SPLIT_MANIFEST_SHA256_EXPECTED = "e6d60aa56d2ac87a9b82aaca8b3df5b560131943eb0ffb809f0834e8a9db3543"
EXPR_RESULTS_REL = Path("benchmarks/cftr_dms_expression/results_2026-11-23.json")
THERATYPE_FIXTURE_REL = Path("benchmarks/theratype/theratype_open_annotations_2026-11-21.json")
RESULTS_FILE_NAME = "results_2026-11-24.json"

CIF_EXPECTED_SHA256: Dict[str, str] = {
    "6MSM.cif.gz": "82933988abc3ed24a343819d410b39dc0d86913eaba4426dd945e3273e34fa59",
    "5UAK.cif.gz": "8bf00f2612dd961d3b1ccd2f2193781a534ad44f32666a6de8da948afb8d84cf",
    "8EIQ.cif.gz": "d13cbfc8ba6575ab6228f351578c6f9cc5b9dafe13e5913d9b27b3cacebc18d6",
}

TMD_RANGES: Tuple[Tuple[int, int], ...] = ((77, 360), (859, 1155))
RBF_SCALE_A = 12.0
CONTACT_RADIUS_A = 10.0
RMSD_MAX_8EIQ_6MSM_A = 1.5
PARITY_TOL = 1e-6

FEATURE_NAMES: Tuple[str, ...] = (
    "rbf_wjx",
    "rbf_cv6",
    "cn_10A_6msm",
    "bfactor_6msm",
    "disp_5uak_6msm",
    "struct_unresolved",
)

SUBDOMAIN_ORDER: Tuple[str, ...] = ("Lasso", "TMD1", "NBD1", "R Domain", "TMD2", "NBD2")
ALPHA_GRID: Tuple[float, ...] = (0.01, 0.1, 1.0, 10.0, 100.0, 1000.0)
CV_FOLDS = 5
BOOTSTRAP_SEED = 20261124
BOOTSTRAP_RESAMPLES = 10000
SHUFFLE_SEED = 20261124
GATE_REL_MAE_REDUCTION = 0.10
MIN_HOLDOUT_CLUSTERS = 40

WT_PARITY_ALLOWLIST: Dict[Tuple[int, int, str], str] = {
    (89, 470, "M470V"): (
        "Common CFTR rs213950 polymorphism: Tedman background M470 vs canonical "
        "P13569/6MSM/5UAK/8EIQ VAL470; mapped to auth_seq_id=470."
    ),
    (217, 966, "Q966="): (
        "Source-table discrepancy: Q966= vs THR966 in 6MSM/5UAK/8EIQ/P13569 "
        "(2988G>A corresponds to Q996=); kept in dev without exclusion and treated as "
        "unmapped (struct_unresolved=1.0, 5 continuous 3D features missing)."
    ),
}

THREE_TO_ONE: Dict[str, str] = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

XLSX_NS = {"s": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
CELL_REF_RE = re.compile(r"^([A-Z]+)(\d+)$")
RANGE_REF_RE = re.compile(r"^([A-Z]+)(\d+):([A-Z]+)(\d+)$")


@dataclass(frozen=True)
class CaAtom:
    """Single Chain A C-alpha record indexed by `auth_seq_id`."""

    auth_seq_id: int
    label_seq_id: int
    comp_id: str
    xyz: Tuple[float, float, float]
    bfactor: float


@dataclass(frozen=True)
class Variant3DRecord:
    """Single DMS row paired with pure 3D Cryo-EM structural features (no sub-domain/class predictors)."""

    row_num: int
    sub_domain: str
    nucleic_acid_mutation: str
    protein_distinction: str
    codon_position: int
    variant_class: str
    split: str
    exposed_pred0010_forced: bool
    dmso_rel_mut_wt: float
    combo_rel_mut_wt: float
    x_dmso_log2: float
    y_combo_log2: float
    features_raw: Tuple[Optional[float], ...]


def resolve_tedman_data_dir(
    data_dir: Optional[Path] = None, repo_root: Path = REPO_ROOT
) -> Optional[Path]:
    """Locate git-ignored `tedman-2025/` data folder across main checkout or `.worktrees/cf-research/...`."""
    if data_dir is not None:
        p = Path(data_dir)
        return p if (p / MASTER_XLSX_NAME).is_file() else None
    direct = repo_root / "sources" / "downloads" / "tedman-2025"
    if (direct / MASTER_XLSX_NAME).is_file():
        return direct
    if len(repo_root.parents) >= 3:
        main_sibling = repo_root.parents[2] / "cf-research" / "sources" / "downloads" / "tedman-2025"
        if (main_sibling / MASTER_XLSX_NAME).is_file():
            return main_sibling
    return None


def _col_to_index(col: str) -> int:
    val = 0
    for ch in col:
        val = val * 26 + (ord(ch) - ord("A") + 1)
    return val


def _index_to_col(idx: int) -> str:
    chars: List[str] = []
    n = idx
    while n > 0:
        n, rem = divmod(n - 1, 26)
        chars.append(chr(ord("A") + rem))
    return "".join(reversed(chars))


def _extract_cell_text(cell_elem: ET.Element, shared_strings: List[str]) -> Optional[str]:
    ctype = cell_elem.attrib.get("t", "")
    if ctype == "inlineStr":
        is_elem = cell_elem.find("s:is", XLSX_NS)
        if is_elem is None:
            return None
        parts = [t.text or "" for t in is_elem.iter("{http://schemas.openxmlformats.org/spreadsheetml/2006/main}t")]
        joined = "".join(parts)
        return joined if joined != "" else None
    v_elem = cell_elem.find("s:v", XLSX_NS)
    if v_elem is None or v_elem.text is None:
        return None
    raw = v_elem.text
    if ctype == "s":
        idx = int(raw)
        if 0 <= idx < len(shared_strings):
            return shared_strings[idx]
        raise ValueError(f"Shared string index out of range: {idx}")
    if ctype == "b":
        return "TRUE" if raw == "1" else "FALSE"
    return raw


def parse_xlsx_sheet(
    xlsx_path: Path,
    sheet_name: str,
    *,
    expand_merged: bool = False,
) -> Tuple[Dict[Tuple[int, str], str], int]:
    """Parse sparse `.xlsx` worksheet XML by exact `(row, col)` coordinates."""
    with zipfile.ZipFile(xlsx_path, "r") as xzf:
        shared_strings: List[str] = []
        if "xl/sharedStrings.xml" in xzf.namelist():
            sroot = ET.fromstring(xzf.read("xl/sharedStrings.xml"))
            for si in sroot.findall("s:si", XLSX_NS):
                text = "".join(
                    t.text or ""
                    for t in si.iter("{http://schemas.openxmlformats.org/spreadsheetml/2006/main}t")
                )
                shared_strings.append(text)

        wb_root = ET.fromstring(xzf.read("xl/workbook.xml"))
        rels_root = ET.fromstring(xzf.read("xl/_rels/workbook.xml.rels"))
        rmap = {r.attrib["Id"]: r.attrib["Target"] for r in rels_root}

        target_xml: Optional[str] = None
        for sheet_elem in wb_root.findall(".//s:sheet", XLSX_NS):
            if sheet_elem.attrib.get("name") == sheet_name:
                rid = sheet_elem.attrib["{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"]
                raw_target = rmap[rid].lstrip("/")
                target_xml = raw_target if raw_target.startswith("xl/") else f"xl/{raw_target}"
                break
        if target_xml is None:
            raise ValueError(f"Sheet {sheet_name!r} not found in {xlsx_path.name}")

        sroot = ET.fromstring(xzf.read(target_xml))
        cells: Dict[Tuple[int, str], str] = {}
        max_row = 0
        for row_elem in sroot.findall(".//s:sheetData/s:row", XLSX_NS):
            row_fallback = int(row_elem.attrib.get("r", 0) or 0)
            for c_elem in row_elem.findall("s:c", XLSX_NS):
                coord = c_elem.attrib.get("r", "")
                m = CELL_REF_RE.fullmatch(coord)
                if not m:
                    raise ValueError(f"Missing/invalid sparse coordinate r={coord!r} in {sheet_name!r}")
                col = m.group(1)
                row_idx = int(m.group(2))
                if row_idx > max_row:
                    max_row = row_idx
                if row_fallback and row_idx != row_fallback:
                    raise ValueError(f"Mismatched row coord {row_idx} != {row_fallback}")
                val = _extract_cell_text(c_elem, shared_strings)
                if val is not None:
                    cells[(row_idx, col)] = val

        if expand_merged:
            for mc in sroot.findall(".//s:mergeCells/s:mergeCell", XLSX_NS):
                ref = mc.attrib.get("ref", "")
                m = RANGE_REF_RE.fullmatch(ref)
                if not m:
                    continue
                c1, r1, c2, r2 = (
                    _col_to_index(m.group(1)),
                    int(m.group(2)),
                    _col_to_index(m.group(3)),
                    int(m.group(4)),
                )
                top_left = cells.get((r1, _index_to_col(c1)))
                if top_left is not None:
                    for rr in range(r1, r2 + 1):
                        for cc in range(c1, c2 + 1):
                            cells.setdefault((rr, _index_to_col(cc)), top_left)

        return cells, max_row


def parse_protein_distinction(raw: str) -> Tuple[Optional[int], str, Optional[str]]:
    """Return `(codon_position, variant_class, wt_one_letter)` from `Protein distinction`."""
    s = str(raw or "").strip()
    if not s:
        return None, "other", None
    if "/" in s:
        first = s.split("/", 1)[0].strip()
        codon, _, wt = parse_protein_distinction(first)
        return codon, "other", wt
    m_indel = re.fullmatch(r"([A-Z])(\d+)(del|DEL|ins|INS)", s)
    if m_indel:
        return int(m_indel.group(2)), "inframe_indel", m_indel.group(1)
    m_stop = re.fullmatch(r"([A-Z])(\d+)[X\*]", s)
    if m_stop:
        return int(m_stop.group(2)), "nonsense", m_stop.group(1)
    m_sub = re.fullmatch(r"([A-Z])(\d+)([A-Z=])", s)
    if m_sub:
        ref, codon_s, alt = m_sub.group(1), m_sub.group(2), m_sub.group(3)
        codon = int(codon_s)
        vclass = "other" if alt in ("=", ref) else "missense"
        return codon, vclass, ref
    return None, "other", None


def _parse_pos_float(raw: Optional[str]) -> Optional[float]:
    if raw is None or str(raw).strip() == "":
        return None
    try:
        val = float(raw)
    except ValueError:
        return None
    if not math.isfinite(val) or val <= 0.0:
        return None
    return val


def parse_mmcif_text(
    cif_text: str,
    pdb_id: str,
    *,
    seq_id_field: str = "auth_seq_id",
) -> Tuple[Dict[int, CaAtom], Dict[str, List[Tuple[float, float, float]]]]:
    """Dynamically parse `_atom_site.*` columns into Chain A C-alpha atoms and HETATM ligand heavy atoms."""
    cols: List[str] = []
    mode: Optional[str] = None
    ca_atoms: Dict[int, CaAtom] = {}
    ligands: Dict[str, List[Tuple[float, float, float]]] = {}
    required = {
        "group_PDB", "pdbx_PDB_model_num", "auth_asym_id", "label_atom_id",
        "label_comp_id", "auth_seq_id", "label_seq_id", "Cartn_x", "Cartn_y",
        "Cartn_z", "B_iso_or_equiv", "type_symbol", "label_alt_id", "pdbx_PDB_ins_code",
    }

    for raw in cif_text.splitlines():
        s = raw.strip()
        if s == "loop_":
            cols = []
            mode = "loop_hdr"
            continue
        if mode == "loop_hdr" and s.startswith("_"):
            cols.append(s)
            continue
        if mode == "loop_hdr" and cols:
            if cols[0].startswith("_atom_site."):
                col_names = [c.split(".", 1)[1] for c in cols]
                missing = required - set(col_names)
                if missing:
                    raise ValueError(f"{pdb_id}: missing _atom_site columns {sorted(missing)}")
                cols = col_names
                mode = "atom"
            else:
                mode = None
        if mode == "atom":
            if not s or s.startswith("#") or s.startswith("_") or s == "loop_":
                break
            toks = shlex.split(s)
            if len(toks) != len(cols):
                raise ValueError(f"{pdb_id}: malformed _atom_site row length {len(toks)} != {len(cols)}")
            row = dict(zip(cols, toks))
            if row["pdbx_PDB_model_num"] != "1":
                continue
            if row["pdbx_PDB_ins_code"] not in {"?", "."}:
                continue
            if row["label_alt_id"] not in {".", "A"}:
                continue

            grp = row["group_PDB"]
            if grp == "ATOM":
                if row["auth_asym_id"] != "A" or row["label_atom_id"] != "CA":
                    continue
                raw_seq = row[seq_id_field]
                if not raw_seq.lstrip("-").isdigit():
                    raise ValueError(f"{pdb_id}: non-integer {seq_id_field}={raw_seq!r}")
                seq_num = int(raw_seq)
                if seq_num in ca_atoms:
                    raise ValueError(f"{pdb_id}: duplicate CA at {seq_id_field}={seq_num}")
                raw_label = row["label_seq_id"]
                label_seq = int(raw_label) if raw_label.lstrip("-").isdigit() else -1
                xyz = (float(row["Cartn_x"]), float(row["Cartn_y"]), float(row["Cartn_z"]))
                ca_atoms[seq_num] = CaAtom(
                    auth_seq_id=seq_num,
                    label_seq_id=label_seq,
                    comp_id=row["label_comp_id"],
                    xyz=xyz,
                    bfactor=float(row["B_iso_or_equiv"]),
                )
            elif grp == "HETATM":
                if row["type_symbol"] in {"H", "D"}:
                    continue
                comp = row["label_comp_id"]
                xyz = (float(row["Cartn_x"]), float(row["Cartn_y"]), float(row["Cartn_z"]))
                ligands.setdefault(comp, []).append(xyz)

    if not ca_atoms:
        raise ValueError(f"{pdb_id}: no Chain A CA atoms parsed from _atom_site")
    return ca_atoms, ligands


def load_mmcif_gz(path: Path, pdb_id: str, expected_sha256: str) -> Tuple[Dict[int, CaAtom], Dict[str, List[Tuple[float, float, float]]]]:
    raw_bytes = path.read_bytes()
    digest = hashlib.sha256(raw_bytes).hexdigest()
    if digest != expected_sha256:
        raise ValueError(f"{pdb_id} SHA256 mismatch: {digest} != {expected_sha256}")
    text = gzip.decompress(raw_bytes).decode("utf-8")
    return parse_mmcif_text(text, pdb_id, seq_id_field="auth_seq_id")


def _dot3(a: Sequence[float], b: Sequence[float]) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _cross3(a: Sequence[float], b: Sequence[float]) -> Tuple[float, float, float]:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def _det3(m: Sequence[Sequence[float]]) -> float:
    return (
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    )


def _jacobi_symmetric_3x3(mat: Sequence[Sequence[float]], sweeps: int = 60) -> Tuple[List[float], List[List[float]]]:
    a = [list(row) for row in mat]
    v = [[1.0 if i == j else 0.0 for j in range(3)] for i in range(3)]
    for _ in range(sweeps):
        off = a[0][1] ** 2 + a[0][2] ** 2 + a[1][2] ** 2
        if off < 1e-24:
            break
        for p, q in ((0, 1), (0, 2), (1, 2)):
            apq = a[p][q]
            if abs(apq) < 1e-300:
                continue
            tau = (a[q][q] - a[p][p]) / (2.0 * apq)
            t = (1.0 if tau >= 0.0 else -1.0) / (abs(tau) + math.sqrt(tau * tau + 1.0))
            c = 1.0 / math.sqrt(t * t + 1.0)
            s = t * c
            for k in range(3):
                akp, akq = a[k][p], a[k][q]
                a[k][p], a[k][q] = c * akp - s * akq, s * akp + c * akq
            for k in range(3):
                apk, aqk = a[p][k], a[q][k]
                a[p][k], a[q][k] = c * apk - s * aqk, s * apk + c * aqk
            for k in range(3):
                vkp, vkq = v[k][p], v[k][q]
                v[k][p], v[k][q] = c * vkp - s * vkq, s * vkp + c * vkq
    return [a[0][0], a[1][1], a[2][2]], v


def kabsch_svd_superpose(
    moving_pts: Sequence[Tuple[float, float, float]],
    ref_pts: Sequence[Tuple[float, float, float]],
) -> Tuple[List[List[float]], Tuple[float, float, float]]:
    """Compute proper SO(3) Kabsch rigid motion `(R, t)` minimizing `||R * p + t - q||^2` via 3x3 SVD."""
    n = len(moving_pts)
    if n < 3 or n != len(ref_pts):
        raise ValueError(f"Need >= 3 paired points for Kabsch superposition, got {n} vs {len(ref_pts)}")
    pc = tuple(sum(p[k] for p in moving_pts) / n for k in range(3))
    qc = tuple(sum(q[k] for q in ref_pts) / n for k in range(3))

    # Cross-covariance H = sum (q_i - qc)(p_i - pc)^T so q' ~ R p'.
    h = [[0.0] * 3 for _ in range(3)]
    for p, q in zip(moving_pts, ref_pts):
        dp = (p[0] - pc[0], p[1] - pc[1], p[2] - pc[2])
        dq = (q[0] - qc[0], q[1] - qc[1], q[2] - qc[2])
        for i in range(3):
            for j in range(3):
                h[i][j] += dq[i] * dp[j]

    hth = [[sum(h[k][i] * h[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
    evals, evecs = _jacobi_symmetric_3x3(hth)
    order = sorted(range(3), key=lambda idx: evals[idx], reverse=True)
    sigmas = [math.sqrt(max(0.0, evals[idx])) for idx in order]
    v_cols = [[evecs[r][idx] for r in range(3)] for idx in order]

    # Enforce right-handed V basis.
    v_mat = [[v_cols[j][i] for j in range(3)] for i in range(3)]
    if _det3(v_mat) < 0.0:
        v_cols[2] = [-x for x in v_cols[2]]

    u_cols: List[List[float]] = []
    for col_idx in range(3):
        if sigmas[col_idx] > 1e-10:
            hv = [sum(h[r][c] * v_cols[col_idx][c] for c in range(3)) / sigmas[col_idx] for r in range(3)]
            u_cols.append(hv)
    if len(u_cols) == 2:
        u_cols.append(list(_cross3(u_cols[0], u_cols[1])))
    elif len(u_cols) < 2:
        raise ValueError("Degenerate point cloud in Kabsch SVD superposition (rank < 2)")

    # Modified Gram-Schmidt orthonormalization on U columns for numerical precision.
    u0_norm = math.sqrt(_dot3(u_cols[0], u_cols[0]))
    u0 = [x / u0_norm for x in u_cols[0]]
    proj10 = _dot3(u_cols[1], u0)
    u1_raw = [u_cols[1][i] - proj10 * u0[i] for i in range(3)]
    u1_norm = math.sqrt(_dot3(u1_raw, u1_raw))
    u1 = [x / u1_norm for x in u1_raw]
    u2 = list(_cross3(u0, u1))
    # Preserve original orientation of u_cols[2] before reflection sign guard.
    if _dot3(u2, u_cols[2]) < 0.0:
        u2 = [-x for x in u2]

    u_mat = [[u0[i], u1[i], u2[i]] for i in range(3)]
    v_mat = [[v_cols[0][i], v_cols[1][i], v_cols[2][i]] for i in range(3)]
    uvt = [[sum(u_mat[i][k] * v_mat[j][k] for k in range(3)) for j in range(3)] for i in range(3)]
    if _det3(uvt) < 0.0:
        u_mat[0][2], u_mat[1][2], u_mat[2][2] = -u_mat[0][2], -u_mat[1][2], -u_mat[2][2]

    rot = [[sum(u_mat[i][k] * v_mat[j][k] for k in range(3)) for j in range(3)] for i in range(3)]
    if abs(_det3(rot) - 1.0) > 1e-6:
        raise ValueError(f"Expected proper rotation det=+1, got {_det3(rot)}")
    trans = (
        qc[0] - (rot[0][0] * pc[0] + rot[0][1] * pc[1] + rot[0][2] * pc[2]),
        qc[1] - (rot[1][0] * pc[0] + rot[1][1] * pc[1] + rot[1][2] * pc[2]),
        qc[2] - (rot[2][0] * pc[0] + rot[2][1] * pc[1] + rot[2][2] * pc[2]),
    )
    return rot, trans


def apply_rigid_transform(
    rot: Sequence[Sequence[float]],
    trans: Tuple[float, float, float],
    pt: Tuple[float, float, float],
) -> Tuple[float, float, float]:
    return (
        rot[0][0] * pt[0] + rot[0][1] * pt[1] + rot[0][2] * pt[2] + trans[0],
        rot[1][0] * pt[0] + rot[1][1] * pt[1] + rot[1][2] * pt[2] + trans[1],
        rot[2][0] * pt[0] + rot[2][1] * pt[1] + rot[2][2] * pt[2] + trans[2],
    )


def _in_tmd_core(seq_num: int) -> bool:
    return any(lo <= seq_num <= hi for lo, hi in TMD_RANGES)


def build_structural_context(repo_root: Path = REPO_ROOT) -> Dict[str, Any]:
    """Load offline Cryo-EM mmCIFs, verify 508 & TMD alignment, and precompute ligand/C-alpha frames."""
    data_dir = repo_root / "tools" / "cftr-model" / "data"
    ca_6msm, _ = load_mmcif_gz(data_dir / "6MSM.cif.gz", "6MSM", CIF_EXPECTED_SHA256["6MSM.cif.gz"])
    ca_5uak, _ = load_mmcif_gz(data_dir / "5UAK.cif.gz", "5UAK", CIF_EXPECTED_SHA256["5UAK.cif.gz"])
    ca_8eiq, lig_8eiq = load_mmcif_gz(data_dir / "8EIQ.cif.gz", "8EIQ", CIF_EXPECTED_SHA256["8EIQ.cif.gz"])

    if 508 not in ca_6msm or ca_6msm[508].comp_id != "PHE":
        raise AssertionError("Expected 6MSM auth_seq_id 508 to be PHE.")
    if 508 in ca_8eiq:
        raise AssertionError("Expected 8EIQ auth_seq_id 508 to be absent (deltaF508).")
    if 509 not in ca_8eiq or ca_8eiq[509].label_seq_id != 508 or ca_8eiq[509].comp_id != "GLY":
        raise AssertionError("Expected 8EIQ auth_seq_id 509 GLY with label_seq_id=508 (-1 shift check).")

    common_8_6 = sorted(k for k in ca_8eiq if k in ca_6msm and _in_tmd_core(k))
    common_5_6 = sorted(k for k in ca_5uak if k in ca_6msm and _in_tmd_core(k))
    rot_8, trans_8 = kabsch_svd_superpose(
        [ca_8eiq[k].xyz for k in common_8_6],
        [ca_6msm[k].xyz for k in common_8_6],
    )
    rot_5, trans_5 = kabsch_svd_superpose(
        [ca_5uak[k].xyz for k in common_5_6],
        [ca_6msm[k].xyz for k in common_5_6],
    )

    rmsd_8_6 = math.sqrt(
        sum(math.dist(apply_rigid_transform(rot_8, trans_8, ca_8eiq[k].xyz), ca_6msm[k].xyz) ** 2 for k in common_8_6)
        / len(common_8_6)
    )
    rmsd_5_6 = math.sqrt(
        sum(math.dist(apply_rigid_transform(rot_5, trans_5, ca_5uak[k].xyz), ca_6msm[k].xyz) ** 2 for k in common_5_6)
        / len(common_5_6)
    )
    if rmsd_8_6 >= RMSD_MAX_8EIQ_6MSM_A:
        raise AssertionError(f"TMD core RMSD(8EIQ->6MSM)={rmsd_8_6:.6f} >= {RMSD_MAX_8EIQ_6MSM_A}")

    wjx_pts = [apply_rigid_transform(rot_8, trans_8, p) for p in lig_8eiq.get("WJX", ())]
    cv6_pts = [apply_rigid_transform(rot_8, trans_8, p) for p in lig_8eiq.get("CV6", ())]
    if not wjx_pts or not cv6_pts:
        raise AssertionError("8EIQ missing required WJX or CV6 HETATM heavy atoms.")

    min_wjx_cv6 = min(math.dist(a, b) for a in wjx_pts for b in cv6_pts)
    wjx_cent = tuple(sum(p[k] for p in wjx_pts) / len(wjx_pts) for k in range(3))
    cv6_cent = tuple(sum(p[k] for p in cv6_pts) / len(cv6_pts) for k in range(3))
    cent_wjx_cv6 = math.dist(wjx_cent, cv6_cent)

    ca_5uak_aligned = {k: apply_rigid_transform(rot_5, trans_5, atom.xyz) for k, atom in ca_5uak.items()}
    return {
        "ca_6msm": ca_6msm,
        "ca_5uak": ca_5uak,
        "ca_8eiq": ca_8eiq,
        "ca_5uak_aligned": ca_5uak_aligned,
        "wjx_pts": wjx_pts,
        "cv6_pts": cv6_pts,
        "metrics": {
            "tmd_common_ca_8eiq_6msm": len(common_8_6),
            "tmd_common_ca_5uak_6msm": len(common_5_6),
            "rmsd_tmd_8eiq_to_6msm_A": round(rmsd_8_6, 6),
            "rmsd_tmd_5uak_to_6msm_A": round(rmsd_5_6, 6),
            "wjx_heavy_atoms": len(wjx_pts),
            "cv6_heavy_atoms": len(cv6_pts),
            "wjx_cv6_min_heavy_atom_dist_A": round(min_wjx_cv6, 6),
            "wjx_cv6_centroid_dist_A": round(cent_wjx_cv6, 6),
        },
    }


def verify_wt_sequence_parity(
    assignments: Sequence[Mapping[str, Any]],
    struct_ctx: Mapping[str, Any],
) -> Dict[str, Any]:
    """Fail closed if any resolved codon differs from WT `Protein distinction` outside allowlisted rules."""
    ca_6msm: Mapping[int, CaAtom] = struct_ctx["ca_6msm"]
    ca_8eiq: Mapping[int, CaAtom] = struct_ctx["ca_8eiq"]

    unexpected: List[Dict[str, Any]] = []
    documented_observed: List[Dict[str, Any]] = []
    seen_allowlist_keys = set()

    for row in assignments:
        row_num = int(row["row_num"])
        codon = int(row["codon_position"])
        prot = str(row["protein_distinction"])
        split = str(row["split"])
        _, _, wt_aa = parse_protein_distinction(prot)
        if wt_aa is None:
            unexpected.append({"row_num": row_num, "codon": codon, "prot": prot, "reason": "unparseable_wt_aa"})
            continue

        for pdb_id, ca_map in (("6MSM", ca_6msm), ("8EIQ", ca_8eiq)):
            if pdb_id == "8EIQ" and codon == 508:
                continue
            if codon not in ca_map:
                continue
            obs_aa = THREE_TO_ONE.get(ca_map[codon].comp_id, "?")
            if obs_aa != wt_aa:
                key = (row_num, codon, prot)
                if key in WT_PARITY_ALLOWLIST and split == "dev":
                    seen_allowlist_keys.add(key)
                    documented_observed.append({
                        "pdb_id": pdb_id,
                        "row_num": row_num,
                        "codon_position": codon,
                        "protein_distinction": prot,
                        "split": split,
                        "expected_wt_aa": wt_aa,
                        "observed_mmcif_aa": obs_aa,
                        "ruling": WT_PARITY_ALLOWLIST[key],
                    })
                else:
                    unexpected.append({
                        "pdb_id": pdb_id,
                        "row_num": row_num,
                        "codon_position": codon,
                        "protein_distinction": prot,
                        "split": split,
                        "expected_wt_aa": wt_aa,
                        "observed_mmcif_aa": obs_aa,
                    })

    if unexpected:
        raise AssertionError(f"WT amino-acid sequence parity failed closed: {unexpected}")
    if seen_allowlist_keys != set(WT_PARITY_ALLOWLIST):
        raise AssertionError(
            f"Expected documented WT parity exceptions {sorted(WT_PARITY_ALLOWLIST)}, saw {sorted(seen_allowlist_keys)}"
        )
    return {
        "passed": True,
        "unexpected_mismatches": 0,
        "holdout_mismatches": 0,
        "documented_dev_exceptions": documented_observed,
    }


def extract_6_structural_features(
    row_num: int,
    codon: int,
    prot: str,
    struct_ctx: Mapping[str, Any],
) -> Tuple[Optional[float], ...]:
    """Compute `(rbf_wjx, rbf_cv6, cn_10A_6msm, bfactor_6msm, disp_5uak_6msm, struct_unresolved)`."""
    if (row_num, codon, prot) == (217, 966, "Q966="):
        return (None, None, None, None, None, 1.0)

    ca_6msm: Mapping[int, CaAtom] = struct_ctx["ca_6msm"]
    ca_5uak_aligned: Mapping[int, Tuple[float, float, float]] = struct_ctx["ca_5uak_aligned"]
    wjx_pts: Sequence[Tuple[float, float, float]] = struct_ctx["wjx_pts"]
    cv6_pts: Sequence[Tuple[float, float, float]] = struct_ctx["cv6_pts"]

    in_6msm = codon in ca_6msm
    in_5uak = codon in ca_5uak_aligned
    struct_unresolved = 0.0 if (in_6msm and in_5uak) else 1.0

    if not in_6msm:
        return (None, None, None, None, None, struct_unresolved)

    atom6 = ca_6msm[codon]
    p6 = atom6.xyz
    d_wjx = min(math.dist(p6, q) for q in wjx_pts)
    d_cv6 = min(math.dist(p6, q) for q in cv6_pts)
    rbf_wjx = math.exp(-((d_wjx / RBF_SCALE_A) ** 2))
    rbf_cv6 = math.exp(-((d_cv6 / RBF_SCALE_A) ** 2))
    cn_10a = float(
        sum(1 for seq_o, atom_o in ca_6msm.items() if seq_o != codon and math.dist(p6, atom_o.xyz) <= CONTACT_RADIUS_A)
    )
    bfactor_6msm = float(atom6.bfactor)
    disp_5uak = math.dist(ca_5uak_aligned[codon], p6) if in_5uak else None
    return (rbf_wjx, rbf_cv6, cn_10a, bfactor_6msm, disp_5uak, struct_unresolved)


def load_locked_split_manifest(repo_root: Path = REPO_ROOT) -> Tuple[Dict[str, Any], str]:
    manifest_path = repo_root / SPLIT_MANIFEST_REL
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest_sha = hashlib.sha256(json.dumps(manifest, sort_keys=True).encode("utf-8")).hexdigest()
    if manifest_sha != SPLIT_MANIFEST_SHA256_EXPECTED:
        raise AssertionError(
            f"split_manifest.json SHA256 mismatch: {manifest_sha} != {SPLIT_MANIFEST_SHA256_EXPECTED}"
        )
    th_path = repo_root / THERATYPE_FIXTURE_REL
    exposed = {int(v["codon_position"]) for v in json.loads(th_path.read_text(encoding="utf-8"))["variants"]}
    assigns = manifest["assignments"]
    dev_codons = {int(a["codon_position"]) for a in assigns if a["split"] == "dev"}
    hold_codons = {int(a["codon_position"]) for a in assigns if a["split"] == "holdout"}
    if dev_codons & hold_codons:
        raise AssertionError("Leakage: dev and holdout share codons.")
    if hold_codons & exposed:
        raise AssertionError("Leakage: exposed PRED-0010 codons present in holdout.")
    if len(hold_codons) < MIN_HOLDOUT_CLUSTERS:
        raise AssertionError(f"Insufficient holdout clusters: {len(hold_codons)} < {MIN_HOLDOUT_CLUSTERS}")
    return manifest, manifest_sha


def load_tedman_with_3d_features(
    master_xlsx: Path,
    manifest: Mapping[str, Any],
    struct_ctx: Mapping[str, Any],
) -> List[Variant3DRecord]:
    """Parse `All_variants`, assert 234/234 parity with `split_manifest.json`, and attach 3D features."""
    master_sha = hashlib.sha256(master_xlsx.read_bytes()).hexdigest()
    if master_sha != MASTER_XLSX_SHA256:
        raise ValueError(f"Master xlsx SHA256 mismatch: {master_sha} != {MASTER_XLSX_SHA256}")

    cells, max_row = parse_xlsx_sheet(master_xlsx, "All_variants", expand_merged=True)
    by_row = {int(a["row_num"]): a for a in manifest["assignments"]}
    records: List[Variant3DRecord] = []
    current_sub = ""

    for row_num in range(3, max_row + 1):
        sub_raw = cells.get((row_num, "A"), "").strip()
        if sub_raw:
            current_sub = sub_raw
        nuc = cells.get((row_num, "B"), "").strip()
        prot = cells.get((row_num, "C"), "").strip()
        if not nuc and not prot:
            continue
        if row_num not in by_row:
            raise AssertionError(f"Unexpected extra variant row {row_num} absent from split_manifest.json")
        info = by_row[row_num]
        if current_sub != info["sub_domain"]:
            raise AssertionError(f"Row {row_num} sub_domain mismatch: {current_sub!r} != {info['sub_domain']!r}")
        if prot != info["protein_distinction"] or nuc != info["nucleic_acid_mutation"]:
            raise AssertionError(f"Row {row_num} variant label mismatch vs split_manifest.json")

        codon, vclass, _ = parse_protein_distinction(prot)
        if codon != int(info["codon_position"]) or vclass != info["variant_class"]:
            raise AssertionError(f"Row {row_num} parsed codon/class mismatch vs split_manifest.json")

        dmso_rel = _parse_pos_float(cells.get((row_num, "F")))
        combo_rel = _parse_pos_float(cells.get((row_num, "R")))
        if dmso_rel is None or combo_rel is None:
            raise ValueError(f"Row {row_num} has missing or non-positive DMSO/VX-661+VX-445 ratio (zero imputation).")

        feats = extract_6_structural_features(row_num, codon, prot, struct_ctx)
        records.append(
            Variant3DRecord(
                row_num=row_num,
                sub_domain=current_sub,
                nucleic_acid_mutation=nuc,
                protein_distinction=prot,
                codon_position=codon,
                variant_class=vclass,
                split=str(info["split"]),
                exposed_pred0010_forced=bool(info["exposed_pred0010_forced"]),
                dmso_rel_mut_wt=dmso_rel,
                combo_rel_mut_wt=combo_rel,
                x_dmso_log2=math.log2(dmso_rel),
                y_combo_log2=math.log2(combo_rel),
                features_raw=feats,
            )
        )

    if len(records) != len(by_row):
        raise AssertionError(f"Expected {len(by_row)} variants, parsed {len(records)}")
    return records


def _solve_linear_system(a: Sequence[Sequence[float]], b: Sequence[float]) -> List[float]:
    n = len(b)
    mat = [list(row) + [b[i]] for i, row in enumerate(a)]
    for col in range(n):
        pivot = max(range(col, n), key=lambda r: abs(mat[r][col]))
        if abs(mat[pivot][col]) < 1e-14:
            raise ValueError("Singular linear system.")
        mat[col], mat[pivot] = mat[pivot], mat[col]
        diag = mat[col][col]
        for j in range(col, n + 1):
            mat[col][j] /= diag
        for r in range(n):
            if r == col:
                continue
            factor = mat[r][col]
            if factor != 0.0:
                for j in range(col, n + 1):
                    mat[r][j] -= factor * mat[col][j]
    return [mat[i][n] for i in range(n)]


def fit_dmso_ols(xs_dmso: Sequence[float], ys: Sequence[float]) -> Tuple[float, float]:
    mx = sum(xs_dmso) / len(xs_dmso)
    my = sum(ys) / len(ys)
    denom = sum((x - mx) ** 2 for x in xs_dmso)
    if denom <= 1e-14:
        return my, 0.0
    b1 = sum((x - mx) * (y - my) for x, y in zip(xs_dmso, ys)) / denom
    return my - b1 * mx, b1


@dataclass(frozen=True)
class FittedPartial3DRidgeModel:
    """Partial-penalized Ridge on `x_dmso_log2` + 6 fold-imputed/standardized 3D features."""

    alpha: float
    b0_intercept: float
    b1_dmso_slope: float
    impute_means: Tuple[float, ...]
    feat_means: Tuple[float, ...]
    feat_stds: Tuple[float, ...]
    weights: Tuple[float, ...]

    def transform_features(self, raw: Sequence[Optional[float]]) -> Tuple[float, ...]:
        out: List[float] = []
        for j, val in enumerate(raw):
            filled = self.impute_means[j] if val is None else float(val)
            std = self.feat_stds[j]
            out.append((filled - self.feat_means[j]) / std if std > 0.0 else 0.0)
        return tuple(out)

    def predict_one(self, x_dmso_log2: float, raw_features: Sequence[Optional[float]]) -> float:
        z = self.transform_features(raw_features)
        return self.b0_intercept + self.b1_dmso_slope * x_dmso_log2 + sum(w * zv for w, zv in zip(self.weights, z))


def fit_partial_penalized_3d_ridge(
    xs_dmso: Sequence[float],
    raw_features_rows: Sequence[Sequence[Optional[float]]],
    ys: Sequence[float],
    alpha: float,
) -> FittedPartial3DRidgeModel:
    """Fit partial Ridge (`0` penalty on intercept & `x_dmso`; in-fold mean imputation + standardization)."""
    n = len(ys)
    if n < 5:
        raise ValueError("Need at least 5 observations for partial-penalized Ridge.")
    if alpha < 0.0:
        raise ValueError(f"alpha must be >= 0, got {alpha}")
    p_extra = len(raw_features_rows[0])

    impute_means: List[float] = []
    for j in range(p_extra):
        obs = [float(row[j]) for row in raw_features_rows if row[j] is not None]
        if not obs:
            raise ValueError(f"Feature column {j} is entirely missing in training fold.")
        impute_means.append(sum(obs) / len(obs))

    filled = [
        [impute_means[j] if row[j] is None else float(row[j]) for j in range(p_extra)]
        for row in raw_features_rows
    ]
    means = [sum(row[j] for row in filled) / n for j in range(p_extra)]
    stds: List[float] = []
    for j in range(p_extra):
        var = sum((row[j] - means[j]) ** 2 for row in filled) / n
        stds.append(math.sqrt(var) if var > 1e-12 else 0.0)

    p = 2 + p_extra
    x_mat = [
        [1.0, xs_dmso[i]]
        + [(filled[i][j] - means[j]) / stds[j] if stds[j] > 0.0 else 0.0 for j in range(p_extra)]
        for i in range(n)
    ]
    xtx = [[sum(x_mat[i][a] * x_mat[i][b] for i in range(n)) for b in range(p)] for a in range(p)]
    xty = [sum(x_mat[i][a] * ys[i] for i in range(n)) for a in range(p)]
    for j in range(2, p):
        xtx[j][j] += alpha if stds[j - 2] > 0.0 else 1.0

    beta = _solve_linear_system(xtx, xty)
    return FittedPartial3DRidgeModel(
        alpha=alpha,
        b0_intercept=beta[0],
        b1_dmso_slope=beta[1],
        impute_means=tuple(impute_means),
        feat_means=tuple(means),
        feat_stds=tuple(stds),
        weights=tuple(beta[2:]),
    )


def _dev_fold_map(dev_rows: Sequence[Variant3DRecord]) -> Dict[int, int]:
    sorted_codons = sorted({r.codon_position for r in dev_rows})
    return {c: idx % CV_FOLDS for idx, c in enumerate(sorted_codons)}


def _cv_select_alpha(
    dev_rows: Sequence[Variant3DRecord],
    feature_rows: Sequence[Sequence[Optional[float]]],
    fold_map: Mapping[int, int],
) -> Tuple[float, float, List[Dict[str, Any]]]:
    cv_path: List[Dict[str, Any]] = []
    best_alpha = ALPHA_GRID[0]
    best_mae = float("inf")
    for alpha in ALPHA_GRID:
        errs: List[float] = []
        for fold in range(CV_FOLDS):
            tr_idx = [i for i, r in enumerate(dev_rows) if fold_map[r.codon_position] != fold]
            va_idx = [i for i, r in enumerate(dev_rows) if fold_map[r.codon_position] == fold]
            mod = fit_partial_penalized_3d_ridge(
                [dev_rows[i].x_dmso_log2 for i in tr_idx],
                [feature_rows[i] for i in tr_idx],
                [dev_rows[i].y_combo_log2 for i in tr_idx],
                alpha,
            )
            for i in va_idx:
                pred = mod.predict_one(dev_rows[i].x_dmso_log2, feature_rows[i])
                errs.append(abs(dev_rows[i].y_combo_log2 - pred))
        mae = sum(errs) / len(errs)
        cv_path.append({"alpha": alpha, "dev_cv_mae_log2": round(mae, 6)})
        if mae < best_mae:
            best_mae = mae
            best_alpha = alpha
    return best_alpha, best_mae, cv_path


def compute_vif_context(dev_rows: Sequence[Variant3DRecord]) -> Dict[str, Any]:
    """Compute non-gated `dev` variance inflation factors `(R^-1)_jj` for 3D features."""
    raw = [r.features_raw for r in dev_rows]
    mod = fit_partial_penalized_3d_ridge(
        [r.x_dmso_log2 for r in dev_rows], raw, [r.y_combo_log2 for r in dev_rows], alpha=1.0
    )
    z_3d = [mod.transform_features(row) for row in raw]
    n = len(z_3d)

    def _corr_vifs(cols_data: Sequence[Sequence[float]], names: Sequence[str]) -> Dict[str, float]:
        k = len(names)
        corr = [[sum(cols_data[r][i] * cols_data[r][j] for r in range(n)) / n for j in range(k)] for i in range(k)]
        vifs: Dict[str, float] = {}
        for j, name in enumerate(names):
            e_vec = [1.0 if i == j else 0.0 for i in range(k)]
            inv_col = _solve_linear_system(corr, e_vec)
            vifs[name] = round(inv_col[j], 6)
        return vifs

    xs = [r.x_dmso_log2 for r in dev_rows]
    mx = sum(xs) / n
    sx = math.sqrt(sum((x - mx) ** 2 for x in xs) / n)
    z_all = [[(xs[i] - mx) / sx] + list(z_3d[i]) for i in range(n)]
    return {
        "vif_3d_features_only": _corr_vifs(z_3d, FEATURE_NAMES),
        "vif_with_x_dmso_log2": _corr_vifs(z_all, ("x_dmso_log2",) + FEATURE_NAMES),
    }


def residue_clustered_bootstrap_delta_ci(
    holdout_rows: Sequence[Variant3DRecord],
    stronger_preds: Sequence[float],
    cand_preds: Sequence[float],
    *,
    resamples: int = BOOTSTRAP_RESAMPLES,
    seed: int = BOOTSTRAP_SEED,
) -> Dict[str, Any]:
    """Cluster-bootstrap 95% CI on `delta_MAE = MAE_stronger - MAE_candidate` over holdout codons."""
    by_codon: Dict[int, List[Tuple[float, float]]] = {}
    for rec, sp, cp in zip(holdout_rows, stronger_preds, cand_preds):
        by_codon.setdefault(rec.codon_position, []).append(
            (abs(rec.y_combo_log2 - sp), abs(rec.y_combo_log2 - cp))
        )
    clusters = list(by_codon.values())
    n_clusters = len(clusters)
    if n_clusters < 2:
        raise ValueError("Need at least 2 clusters for residue-clustered bootstrap.")

    rng = random.Random(seed)
    deltas: List[float] = []
    for _ in range(resamples):
        sum_s = 0.0
        sum_c = 0.0
        cnt = 0
        for _k in range(n_clusters):
            for es, ec in clusters[rng.randrange(n_clusters)]:
                sum_s += es
                sum_c += ec
                cnt += 1
        deltas.append((sum_s - sum_c) / cnt)

    deltas.sort()
    lo = deltas[max(0, int(0.025 * resamples))]
    hi = deltas[min(resamples - 1, int(0.975 * resamples))]
    return {
        "n_holdout_clusters": n_clusters,
        "resamples": resamples,
        "seed": seed,
        "delta_mae_ci95_log2": [round(lo, 6), round(hi, 6)],
        "ci95_low_gt_zero": lo > 0.0,
    }


def _mae_log2_and_linear(rows: Sequence[Variant3DRecord], preds_log2: Sequence[float]) -> Tuple[float, float]:
    mae_log2 = sum(abs(r.y_combo_log2 - p) for r, p in zip(rows, preds_log2)) / len(rows)
    mae_lin = sum(abs(r.combo_rel_mut_wt - (2.0 ** p)) for r, p in zip(rows, preds_log2)) / len(rows)
    return mae_log2, mae_lin


def evaluate_3d_benchmark(data_dir: Path, repo_root: Path = REPO_ROOT) -> Dict[str, Any]:
    """Run all pre-holdout controls (fail closed), then single-shot holdout evaluation."""
    manifest, manifest_sha = load_locked_split_manifest(repo_root)
    struct_ctx = build_structural_context(repo_root)
    wt_parity = verify_wt_sequence_parity(manifest["assignments"], struct_ctx)

    records = load_tedman_with_3d_features(data_dir / MASTER_XLSX_NAME, manifest, struct_ctx)
    dev_rows = [r for r in records if r.split == "dev"]
    holdout_rows = [r for r in records if r.split == "holdout"]
    dev_codons = sorted({r.codon_position for r in dev_rows})
    holdout_codons = sorted({r.codon_position for r in holdout_rows})

    # Reference #396 metrics (dev parameters & holdout parity check target)
    expr_ref = json.loads((repo_root / EXPR_RESULTS_REL).read_text(encoding="utf-8"))

    # Pre-holdout `dev` baseline fitting & 5-fold CV
    fold_map = _dev_fold_map(dev_rows)
    dev_median_log2 = statistics.median([r.y_combo_log2 for r in dev_rows])
    ols_b0, ols_b1 = fit_dmso_ols([r.x_dmso_log2 for r in dev_rows], [r.y_combo_log2 for r in dev_rows])

    cv_b1_errs: List[float] = []
    cv_b2_errs: List[float] = []
    for fold in range(CV_FOLDS):
        tr = [r for r in dev_rows if fold_map[r.codon_position] != fold]
        va = [r for r in dev_rows if fold_map[r.codon_position] == fold]
        tr_med = statistics.median([r.y_combo_log2 for r in tr])
        fb0, fb1 = fit_dmso_ols([r.x_dmso_log2 for r in tr], [r.y_combo_log2 for r in tr])
        for r in va:
            cv_b1_errs.append(abs(r.y_combo_log2 - tr_med))
            cv_b2_errs.append(abs(r.y_combo_log2 - (fb0 + fb1 * r.x_dmso_log2)))

    dev_cv_b1_mae = sum(cv_b1_errs) / len(cv_b1_errs)
    dev_cv_b2_mae = sum(cv_b2_errs) / len(cv_b2_errs)
    locked_baseline = "dmso_ols" if dev_cv_b2_mae <= dev_cv_b1_mae else "constant_dev_median"
    locked_dev_cv_mae = min(dev_cv_b1_mae, dev_cv_b2_mae)

    ref_models = expr_ref["models_fitted_on_dev"]
    if (
        abs(dev_median_log2 - float(ref_models["constant_dev_median_log2"])) > PARITY_TOL
        or abs(ols_b0 - float(ref_models["dmso_ols"]["b0_intercept"])) > PARITY_TOL
        or abs(ols_b1 - float(ref_models["dmso_ols"]["b1_dmso_slope"])) > PARITY_TOL
        or abs(dev_cv_b2_mae - float(ref_models["dev_cv_5fold_baselines"]["dmso_ols_mae_log2"])) > PARITY_TOL
    ):
        raise AssertionError("Pre-holdout #396 dev baseline fit parity check failed.")

    # Pre-holdout Control 4: Dev-only codon-permuted 3D null control (seed=20261124)
    codon_to_feats: Dict[int, Tuple[Optional[float], ...]] = {}
    for r in dev_rows:
        if r.codon_position not in codon_to_feats:
            codon_to_feats[r.codon_position] = r.features_raw
    perm_codons = list(dev_codons)
    random.Random(SHUFFLE_SEED).shuffle(perm_codons)
    codon_remap = dict(zip(dev_codons, perm_codons))
    permuted_dev_feats = [codon_to_feats[codon_remap[r.codon_position]] for r in dev_rows]
    shuf_alpha, shuf_cv_mae, _ = _cv_select_alpha(dev_rows, permuted_dev_feats, fold_map)
    shuf_rel_red = (locked_dev_cv_mae - shuf_cv_mae) / locked_dev_cv_mae
    if shuf_rel_red >= GATE_REL_MAE_REDUCTION:
        raise AssertionError(
            f"Dev codon-permuted 3D null control exceeded 10% threshold: {shuf_rel_red:.6f}"
        )

    # Dev CV selection for true 3D candidate (zero holdout contact)
    best_alpha, best_cv_mae, cv_path = _cv_select_alpha(
        dev_rows, [r.features_raw for r in dev_rows], fold_map
    )
    final_cand = fit_partial_penalized_3d_ridge(
        [r.x_dmso_log2 for r in dev_rows],
        [r.features_raw for r in dev_rows],
        [r.y_combo_log2 for r in dev_rows],
        best_alpha,
    )
    vifs = compute_vif_context(dev_rows)

    # ---------------- Single-Shot Holdout Evaluation ----------------
    b1_hold_preds = [dev_median_log2] * len(holdout_rows)
    b2_hold_preds = [ols_b0 + ols_b1 * r.x_dmso_log2 for r in holdout_rows]
    cand_hold_preds = [final_cand.predict_one(r.x_dmso_log2, r.features_raw) for r in holdout_rows]

    b1_mae_log2, b1_mae_lin = _mae_log2_and_linear(holdout_rows, b1_hold_preds)
    b2_mae_log2, b2_mae_lin = _mae_log2_and_linear(holdout_rows, b2_hold_preds)
    cand_mae_log2, cand_mae_lin = _mae_log2_and_linear(holdout_rows, cand_hold_preds)

    ref_hold = expr_ref["holdout_evaluation"]
    if abs(b2_mae_log2 - float(ref_hold["dmso_ols"]["mae_log2"])) > PARITY_TOL:
        raise AssertionError(
            f"#396 holdout dmso_ols MAE parity failed: {b2_mae_log2} != {ref_hold['dmso_ols']['mae_log2']}"
        )
    if abs(b1_mae_log2 - float(ref_hold["constant_dev_median"]["mae_log2"])) > PARITY_TOL:
        raise AssertionError(
            f"#396 holdout constant_dev_median MAE parity failed: {b1_mae_log2} != {ref_hold['constant_dev_median']['mae_log2']}"
        )

    stronger_mae_log2 = b2_mae_log2 if locked_baseline == "dmso_ols" else b1_mae_log2
    stronger_preds = b2_hold_preds if locked_baseline == "dmso_ols" else b1_hold_preds
    rel_improvement = (stronger_mae_log2 - cand_mae_log2) / stronger_mae_log2
    boot = residue_clustered_bootstrap_delta_ci(holdout_rows, stronger_preds, cand_hold_preds)

    if rel_improvement >= GATE_REL_MAE_REDUCTION and boot["ci95_low_gt_zero"]:
        ledger_standing = "open"
        gate_verdict = "PASS_3D_STRUCTURE_BEATS_DMSO_OLS"
    elif rel_improvement >= GATE_REL_MAE_REDUCTION:
        ledger_standing = "refuted"
        gate_verdict = "INCONCLUSIVE_CI_STRADDLES_ZERO"
    else:
        ledger_standing = "refuted"
        gate_verdict = "REFUTED_BELOW_10PCT_MAE_REDUCTION_THRESHOLD"

    subdomain_breakdown: List[Dict[str, Any]] = []
    for sd in SUBDOMAIN_ORDER:
        idxs = [i for i, r in enumerate(holdout_rows) if r.sub_domain == sd]
        if not idxs:
            continue
        sd_codons = len({holdout_rows[i].codon_position for i in idxs})
        sd_base_mae = sum(abs(holdout_rows[i].y_combo_log2 - stronger_preds[i]) for i in idxs) / len(idxs)
        sd_cand_mae = sum(abs(holdout_rows[i].y_combo_log2 - cand_hold_preds[i]) for i in idxs) / len(idxs)
        subdomain_breakdown.append({
            "sub_domain": sd,
            "holdout_variants": len(idxs),
            "holdout_codons": sd_codons,
            "stronger_baseline_mae_log2": round(sd_base_mae, 6),
            "candidate_mae_log2": round(sd_cand_mae, 6),
            "delta_mae_log2": round(sd_base_mae - sd_cand_mae, 6),
        })

    done_sha = hashlib.sha256(
        (repo_root / "benchmarks" / "cftr_dms_3d_structure" / "done.md")
        .read_text(encoding="utf-8", errors="replace")
        .replace("\r\n", "\n")
        .encode("utf-8")
    ).hexdigest()

    return {
        "schema_version": "tedman2025-3d-structure-results-v1",
        "benchmark_id": "cftr-dms-3d-structure",
        "evaluated_on": dt.date.today().isoformat(),
        "evaluated_on_local": dt.date.today().isoformat(),
        "evaluated_on_utc": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%d"),
        "done_sha256": done_sha,
        "split_manifest_sha256": manifest_sha,
        "cryoem_mmcif_sha256": CIF_EXPECTED_SHA256,
        "dataset": {
            "source_id": "cftr-dms-tedman-2025",
            "doi": "10.17632/ksfmxv8zts.1",
            "paper_doi": "10.7554/eLife.107180",
            "license": "CC BY 4.0",
            "assay": "HEK293T parental CFTR DMS plasma-membrane intensity (VX-661+VX-445 vs DMSO, Mut/WT)",
            "scale_boundary": "Plasma-membrane expression in HEK293T cells is not CFTR channel conductance or clinical lung benefit.",
            "bfactor_boundary": "6MSM B_iso_or_equiv is a local Cryo-EM map-resolution/flexibility proxy, not a physical atomic displacement parameter.",
            "total_valid_variants": len(records),
            "excluded_count": 0,
            "dev_variants": len(dev_rows),
            "holdout_variants": len(holdout_rows),
            "dev_unique_codons": len(dev_codons),
            "holdout_unique_codons": len(holdout_codons),
            "holdout_pred0010_overlap_codons": 0,
            "dev_struct_unresolved_rows": sum(1 for r in dev_rows if r.features_raw[5] > 0.5),
            "holdout_struct_unresolved_rows": sum(1 for r in holdout_rows if r.features_raw[5] > 0.5),
            "excluded_predictors_for_3d_isolation": ["sub_domain", "variant_class"],
        },
        "pre_holdout_controls": {
            "wt_amino_acid_sequence_parity": wt_parity,
            "cryoem_tmd_superposition_and_pockets": struct_ctx["metrics"],
            "baseline_396_dev_parity_passed": True,
            "dev_codon_permuted_3d_null_control": {
                "seed": SHUFFLE_SEED,
                "selected_alpha_permuted": shuf_alpha,
                "permuted_dev_cv_mae_log2": round(shuf_cv_mae, 6),
                "relative_cv_mae_reduction": round(shuf_rel_red, 6),
                "threshold": GATE_REL_MAE_REDUCTION,
                "passed": True,
            },
        },
        "models_fitted_on_dev": {
            "constant_dev_median_log2": round(dev_median_log2, 6),
            "dmso_ols": {"b0_intercept": round(ols_b0, 6), "b1_dmso_slope": round(ols_b1, 6)},
            "dev_cv_5fold_baselines": {
                "constant_dev_median_mae_log2": round(dev_cv_b1_mae, 6),
                "dmso_ols_mae_log2": round(dev_cv_b2_mae, 6),
                "locked_stronger_baseline": locked_baseline,
            },
            "candidate_dmso_3d_cryoem_ridge": {
                "selected_alpha": best_alpha,
                "dev_cv_mae_log2": round(best_cv_mae, 6),
                "cv_path": cv_path,
                "b0_intercept": round(final_cand.b0_intercept, 6),
                "b1_dmso_slope": round(final_cand.b1_dmso_slope, 6),
                "feature_names": list(FEATURE_NAMES),
                "dev_impute_means": [round(v, 6) for v in final_cand.impute_means],
                "dev_feat_means": [round(v, 6) for v in final_cand.feat_means],
                "dev_feat_stds": [round(v, 6) for v in final_cand.feat_stds],
                "weights": [round(w, 6) for w in final_cand.weights],
                "dev_vif_context": vifs,
            },
        },
        "holdout_evaluation": {
            "primary_unit": "log2(VX-661+VX-445 Mut/WT)",
            "constant_dev_median": {
                "mae_log2": round(b1_mae_log2, 6),
                "mae_linear_context_only": round(b1_mae_lin, 6),
            },
            "dmso_ols": {
                "mae_log2": round(b2_mae_log2, 6),
                "mae_linear_context_only": round(b2_mae_lin, 6),
            },
            "candidate_dmso_3d_cryoem_ridge": {
                "mae_log2": round(cand_mae_log2, 6),
                "mae_linear_context_only": round(cand_mae_lin, 6),
            },
            "stronger_baseline_on_holdout": locked_baseline,
            "stronger_baseline_mae_log2": round(stronger_mae_log2, 6),
            "relative_mae_reduction_vs_stronger": round(rel_improvement, 6),
            "residue_clustered_bootstrap": boot,
            "pre_registered_gate_threshold": GATE_REL_MAE_REDUCTION,
            "gate_verdict": gate_verdict,
            "ledger_standing": ledger_standing,
            "non_gated_subdomain_breakdown": subdomain_breakdown,
        },
    }


def verify_checked_in_artifacts(repo_root: Path = REPO_ROOT, data_dir: Optional[Path] = None) -> Dict[str, Any]:
    """Verify offline Cryo-EM structural invariants and checked-in results (plus raw `.xlsx` parity if present)."""
    bench_dir = repo_root / "benchmarks" / "cftr_dms_3d_structure"
    done_sha = hashlib.sha256(
        (bench_dir / "done.md").read_text(encoding="utf-8", errors="replace").replace("\r\n", "\n").encode("utf-8")
    ).hexdigest()
    if done_sha != DONE_SHA256_EXPECTED:
        raise AssertionError(f"done.md SHA256 mismatch: {done_sha} != {DONE_SHA256_EXPECTED}")

    results = json.loads((bench_dir / RESULTS_FILE_NAME).read_text(encoding="utf-8"))
    if results["done_sha256"] != done_sha:
        raise AssertionError("results JSON done_sha256 mismatch.")

    manifest, manifest_sha = load_locked_split_manifest(repo_root)
    if results["split_manifest_sha256"] != manifest_sha:
        raise AssertionError("results JSON split_manifest_sha256 mismatch.")

    struct_ctx = build_structural_context(repo_root)
    wt_parity = verify_wt_sequence_parity(manifest["assignments"], struct_ctx)
    if not wt_parity["passed"]:
        raise AssertionError("WT amino-acid parity check failed.")

    reported_m = results["pre_holdout_controls"]["cryoem_tmd_superposition_and_pockets"]
    for key, val in struct_ctx["metrics"].items():
        if isinstance(val, float):
            if abs(val - float(reported_m[key])) > PARITY_TOL:
                raise AssertionError(f"Cryo-EM structural metric mismatch on {key}: {val} != {reported_m[key]}")
        elif int(val) != int(reported_m[key]):
            raise AssertionError(f"Cryo-EM structural count mismatch on {key}: {val} != {reported_m[key]}")

    resolved_data = resolve_tedman_data_dir(data_dir, repo_root=repo_root)
    raw_parity_checked = False
    if resolved_data is not None:
        fresh = evaluate_3d_benchmark(resolved_data, repo_root=repo_root)
        f_he = fresh["holdout_evaluation"]
        r_he = results["holdout_evaluation"]
        if (
            abs(f_he["candidate_dmso_3d_cryoem_ridge"]["mae_log2"] - r_he["candidate_dmso_3d_cryoem_ridge"]["mae_log2"]) > PARITY_TOL
            or abs(f_he["dmso_ols"]["mae_log2"] - r_he["dmso_ols"]["mae_log2"]) > PARITY_TOL
            or f_he["gate_verdict"] != r_he["gate_verdict"]
        ):
            raise AssertionError("Live raw .xlsx evaluation differs from checked-in results_2026-11-24.json.")
        raw_parity_checked = True

    return {
        "status": "PASS",
        "raw_local_data_parity_checked": raw_parity_checked,
        "cryoem_tmd_rmsd_8eiq_6msm_A": struct_ctx["metrics"]["rmsd_tmd_8eiq_to_6msm_A"],
        "dev_variants": manifest["dev_variants"],
        "holdout_variants": manifest["holdout_variants"],
        "holdout_unique_codons": manifest["holdout_unique_codons"],
        "gate_verdict": results["holdout_evaluation"]["gate_verdict"],
    }


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--data-dir", type=Path, default=None, help="Optional path to local tedman-2025 folder.")
    parser.add_argument("--write-artifacts", action="store_true", help="Run evaluation and write results_2026-11-24.json.")
    parser.add_argument("--check", action="store_true", help="Verify checked-in structural invariants and local raw parity.")
    parser.add_argument("--json", action="store_true", help="Print strict JSON (allow_nan=False).")
    args = parser.parse_args(argv)

    if args.check:
        summary = verify_checked_in_artifacts(repo_root=REPO_ROOT, data_dir=args.data_dir)
        if args.json:
            print(json.dumps(summary, indent=2, allow_nan=False))
        else:
            print(
                f"PASS cftr_dms_3d_structure: raw_parity={summary['raw_local_data_parity_checked']} "
                f"rmsd_8eiq_6msm={summary['cryoem_tmd_rmsd_8eiq_6msm_A']} "
                f"verdict={summary['gate_verdict']}"
            )
        return 0

    resolved_data = resolve_tedman_data_dir(args.data_dir, repo_root=REPO_ROOT)
    if resolved_data is None:
        print("Local Tedman 2025 workbook not present under sources/downloads/tedman-2025/.", file=sys.stderr)
        return 1

    results = evaluate_3d_benchmark(resolved_data, repo_root=REPO_ROOT)
    if args.write_artifacts:
        (HERE / RESULTS_FILE_NAME).write_text(
            json.dumps(results, indent=2, allow_nan=False) + "\n", encoding="utf-8"
        )

    if args.json:
        print(json.dumps(results, indent=2, allow_nan=False))
    else:
        he = results["holdout_evaluation"]
        print(
            f"Tedman 2025 3D Cryo-EM: verdict={he['gate_verdict']} | "
            f"stronger({he['stronger_baseline_on_holdout']})={he['stronger_baseline_mae_log2']} | "
            f"candidate={he['candidate_dmso_3d_cryoem_ridge']['mae_log2']} | "
            f"rel_reduction={he['relative_mae_reduction_vs_stronger']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
