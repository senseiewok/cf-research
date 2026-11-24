#!/usr/bin/env python3
"""Quantitative CFTR calnexin dependency (`CANX KO`) & single-corrector (`VX-661` vs. `VX-445`) selectivity benchmark.

Evaluates two independently gated pre-registered arms on the locked CC BY 4.0 Tedman et al. 2025
residue-clustered split (`benchmarks/cftr_dms_expression/split_manifest.json`):
- Arm A (`canx_chaperone_coupling`): does adding `log2(CANX_KO_DMSO)` improve prediction of
  parental `log2(VX-661+VX-445)` over parental `1D dmso_ols`?
- Arm B (`canx_ko_single_corrector_selectivity`): note Parental `All_variants` Col J (`VX-661 Mut/WT`)
  is blank (`0/234`) in the published workbook, whereas `All_variants CANX KO` Cols F/J/N/R are
  `234/234` populated; tests whether `log2(CANX_KO_DMSO)`, CFTR sub-domain, and variant
  consequence class predict dampened `CANX KO` single-corrector selectivity
  `S_eps_canx_ko = log2((CANX_KO_VX445_rel + 0.05) / (CANX_KO_VX661_rel + 0.05))` over a
  dev-CV-locked baseline.

Pure Python 3.10+ standard library only; zero network calls; raw `.xlsx` files remain git-ignored
under `sources/downloads/tedman-2025/`. Research benchmark only: HEK293T plasma-membrane
expression is not CFTR channel conductance or clinical benefit.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
import datetime as dt
import hashlib
import json
import math
from pathlib import Path
import random
import re
import statistics
import sys
from typing import Any, Dict, List, Optional, Sequence, Tuple
import xml.etree.ElementTree as ET
import zipfile

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]

MASTER_XLSX_NAME = "ALT 9-6-24 CFTR DMS Expanded Library Master Tables.xlsx"
MASTER_XLSX_SHA256 = "2a702ee3744b8039a9c7c4c74aec9bb30a91301b138d4948e0e2c81a1931e51b"
DONE_SHA256_EXPECTED = "9d701976cea245957e18b3feb532d6ad7b01338928f40397bb5472116aa81187"

SPLIT_MANIFEST_REL = Path("benchmarks/cftr_dms_expression/split_manifest.json")
SPLIT_MANIFEST_SHA256_EXPECTED = "e6d60aa56d2ac87a9b82aaca8b3df5b560131943eb0ffb809f0834e8a9db3543"
THERATYPE_FIXTURE_REL = Path("benchmarks/theratype/theratype_open_annotations_2026-11-21.json")
RESULTS_FILE_NAME = "results_2026-11-23.json"

SUBDOMAIN_CATEGORIES: Tuple[str, ...] = ("Lasso", "TMD1", "NBD1", "R Domain", "TMD2", "NBD2")
EXTRA_SUBDOMAINS: Tuple[str, ...] = ("TMD1", "NBD1", "R Domain", "TMD2", "NBD2")
VARIANT_CLASSES: Tuple[str, ...] = ("missense", "nonsense", "inframe_indel", "other")
EXTRA_CLASSES: Tuple[str, ...] = ("nonsense", "inframe_indel", "other")

ALPHA_GRID: Tuple[float, ...] = (0.01, 0.1, 1.0, 10.0, 100.0, 1000.0)
EPS_DAMPENING = 0.05
CV_FOLDS = 5
BOOTSTRAP_SEED = 20261123
SHUFFLE_SEED = 20261124
BOOTSTRAP_RESAMPLES = 2000
GATE_REL_MAE_REDUCTION = 0.10
MAX_EXCLUSIONS = 5
MIN_HOLDOUT_CLUSTERS = 40
CANX_STABILITY_MIN_R = 0.70

XLSX_NS = {"s": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
CELL_REF_RE = re.compile(r"^([A-Z]+)(\d+)$")
RANGE_REF_RE = re.compile(r"^([A-Z]+)(\d+):([A-Z]+)(\d+)$")


@dataclass(frozen=True)
class DualSheetVariant:
    """Aligned parental + CANX KO DMS measurement row."""

    row_num: int
    sub_domain: str
    nucleic_acid_mutation: str
    protein_distinction: str
    codon_position: int
    variant_class: str
    parental_dmso_rel: float
    parental_vx445_rel: float
    parental_combo_rel: float
    canx_ko_dmso_rel: float
    canx_ko_vx661_rel: float
    canx_ko_vx445_rel: float
    canx_ko_combo_rel: float
    x_dmso_log2: float
    x_canx_ko_dmso_log2: float
    y_combo_log2: float
    s_eps_canx_ko: float
    s_raw_canx_ko: float
    split: str
    exposed_pred0010_forced: bool


def resolve_tedman_data_dir(
    data_dir: Optional[Path] = None, repo_root: Path = REPO_ROOT
) -> Optional[Path]:
    """Locate local git-ignored `tedman-2025/` data folder across main checkout or worktree."""
    if data_dir is not None:
        return Path(data_dir)
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


def parse_protein_distinction(raw: str) -> Tuple[Optional[int], str]:
    """Classify `Protein distinction` into `(codon_position, variant_class)`."""
    s = str(raw or "").strip()
    if not s:
        return None, "other"
    if "/" in s:
        first = s.split("/", 1)[0].strip()
        codon, _ = parse_protein_distinction(first)
        return codon, "other"
    m_indel = re.fullmatch(r"([A-Z])(\d+)(del|DEL|ins|INS)", s)
    if m_indel:
        return int(m_indel.group(2)), "inframe_indel"
    m_stop = re.fullmatch(r"([A-Z])(\d+)[X\*]", s)
    if m_stop:
        return int(m_stop.group(2)), "nonsense"
    m_sub = re.fullmatch(r"([A-Z])(\d+)([A-Z=])", s)
    if m_sub:
        ref, codon_s, alt = m_sub.group(1), m_sub.group(2), m_sub.group(3)
        codon = int(codon_s)
        return (codon, "other") if alt in ("=", ref) else (codon, "missense")
    return None, "other"


def _parse_pos_float(raw: Optional[str]) -> Optional[float]:
    if raw is None or str(raw).strip() == "":
        return None
    try:
        val = float(raw)
    except ValueError:
        return None
    return val if (math.isfinite(val) and val > 0.0) else None


def compute_selectivity_s_eps(vx445_rel: float, vx661_rel: float, *, eps: float = EPS_DAMPENING) -> float:
    """Compute dampened log2 selectivity S_eps = log2((VX445 + eps) / (VX661 + eps))."""
    if vx445_rel <= 0.0 or vx661_rel <= 0.0 or eps < 0.0:
        raise ValueError(f"Invalid positive inputs for selectivity: vx445={vx445_rel}, vx661={vx661_rel}, eps={eps}")
    return math.log2((vx445_rel + eps) / (vx661_rel + eps))


def load_locked_split_manifest(repo_root: Path = REPO_ROOT) -> Tuple[Dict[str, Any], Dict[int, Dict[str, Any]], str]:
    """Load and SHA-256 verify the locked residue-clustered split manifest from Item #1."""
    manifest_path = repo_root / SPLIT_MANIFEST_REL
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest_sha = hashlib.sha256(json.dumps(manifest, sort_keys=True).encode("utf-8")).hexdigest()
    if manifest_sha != SPLIT_MANIFEST_SHA256_EXPECTED:
        raise AssertionError(
            f"split_manifest_sha256 mismatch: {manifest_sha} != {SPLIT_MANIFEST_SHA256_EXPECTED}"
        )

    th = json.loads((repo_root / THERATYPE_FIXTURE_REL).read_text(encoding="utf-8"))
    exposed = {int(v["codon_position"]) for v in th.get("variants") or []}
    assigns = manifest["assignments"]
    dev_codons = {int(a["codon_position"]) for a in assigns if a["split"] == "dev"}
    hold_codons = {int(a["codon_position"]) for a in assigns if a["split"] == "holdout"}
    if dev_codons & hold_codons:
        raise AssertionError("Locked split manifest has codon overlap between dev and holdout.")
    if hold_codons & exposed:
        raise AssertionError("Locked split manifest places exposed PRED-0010 codon in holdout.")
    if len(hold_codons) < MIN_HOLDOUT_CLUSTERS:
        raise AssertionError(f"Insufficient holdout clusters: {len(hold_codons)}")

    by_row = {int(a["row_num"]): a for a in assigns}
    return manifest, by_row, manifest_sha


def load_dual_sheet_tedman_variants(
    master_xlsx: Path,
    by_row_split: Dict[int, Dict[str, Any]],
) -> Tuple[List[DualSheetVariant], List[Dict[str, Any]]]:
    """Load aligned rows across `All_variants` & `All_variants CANX KO` under pre-registered gates."""
    par_cells, par_max = parse_xlsx_sheet(master_xlsx, "All_variants", expand_merged=True)
    ko_cells, ko_max = parse_xlsx_sheet(master_xlsx, "All_variants CANX KO", expand_merged=True)

    max_row = max(par_max, ko_max)
    par_sub = ""
    ko_sub = ""
    aligned_rows: List[DualSheetVariant] = []
    excluded: List[Dict[str, Any]] = []
    aligned_identity_count = 0

    for row_num in range(3, max_row + 1):
        ps_raw = par_cells.get((row_num, "A"), "").strip()
        ks_raw = ko_cells.get((row_num, "A"), "").strip()
        if ps_raw:
            par_sub = ps_raw
        if ks_raw:
            ko_sub = ks_raw

        p_nuc = par_cells.get((row_num, "B"), "").strip()
        p_prot = par_cells.get((row_num, "C"), "").strip()
        k_nuc = ko_cells.get((row_num, "B"), "").strip()
        k_prot = ko_cells.get((row_num, "C"), "").strip()

        if not any((p_nuc, p_prot, k_nuc, k_prot)):
            continue

        if (par_sub, p_nuc, p_prot) != (ko_sub, k_nuc, k_prot):
            raise ValueError(
                f"Dual-sheet alignment mismatch at row {row_num}: "
                f"parental={(par_sub, p_nuc, p_prot)} vs canx_ko={(ko_sub, k_nuc, k_prot)}"
            )
        if par_sub not in SUBDOMAIN_CATEGORIES:
            raise ValueError(f"Unexpected sub_domain {par_sub!r} at row {row_num}")
        if row_num not in by_row_split:
            raise ValueError(f"Row {row_num} ({p_prot}) missing from locked split_manifest.json")

        split_info = by_row_split[row_num]
        if split_info["protein_distinction"] != p_prot or split_info["sub_domain"] != par_sub:
            raise ValueError(f"Manifest identity mismatch at row {row_num}: {split_info} vs {(par_sub, p_prot)}")

        aligned_identity_count += 1
        codon = int(split_info["codon_position"])
        vclass = str(split_info["variant_class"])

        p_dmso = _parse_pos_float(par_cells.get((row_num, "F")))
        p_vx445 = _parse_pos_float(par_cells.get((row_num, "N")))
        p_combo = _parse_pos_float(par_cells.get((row_num, "R")))
        k_dmso = _parse_pos_float(ko_cells.get((row_num, "F")))
        k_vx661 = _parse_pos_float(ko_cells.get((row_num, "J")))
        k_vx445 = _parse_pos_float(ko_cells.get((row_num, "N")))
        k_combo = _parse_pos_float(ko_cells.get((row_num, "R")))

        missing_cols: List[str] = []
        if p_dmso is None:
            missing_cols.append("Parental_DMSO(F)")
        if p_vx445 is None:
            missing_cols.append("Parental_VX445(N)")
        if p_combo is None:
            missing_cols.append("Parental_Combo(R)")
        if k_dmso is None:
            missing_cols.append("CANX_KO_DMSO(F)")
        if k_vx661 is None:
            missing_cols.append("CANX_KO_VX661(J)")
        if k_vx445 is None:
            missing_cols.append("CANX_KO_VX445(N)")
        if k_combo is None:
            missing_cols.append("CANX_KO_Combo(R)")

        if missing_cols:
            excluded.append({
                "row_num": row_num,
                "protein_distinction": p_prot,
                "nucleic_acid_mutation": p_nuc,
                "split": split_info["split"],
                "missing_or_nonpositive": missing_cols,
            })
            continue

        assert p_dmso is not None and p_vx445 is not None and p_combo is not None
        assert k_dmso is not None and k_vx661 is not None and k_vx445 is not None and k_combo is not None
        aligned_rows.append(
            DualSheetVariant(
                row_num=row_num,
                sub_domain=par_sub,
                nucleic_acid_mutation=p_nuc,
                protein_distinction=p_prot,
                codon_position=codon,
                variant_class=vclass,
                parental_dmso_rel=p_dmso,
                parental_vx445_rel=p_vx445,
                parental_combo_rel=p_combo,
                canx_ko_dmso_rel=k_dmso,
                canx_ko_vx661_rel=k_vx661,
                canx_ko_vx445_rel=k_vx445,
                canx_ko_combo_rel=k_combo,
                x_dmso_log2=math.log2(p_dmso),
                x_canx_ko_dmso_log2=math.log2(k_dmso),
                y_combo_log2=math.log2(p_combo),
                s_eps_canx_ko=compute_selectivity_s_eps(k_vx445, k_vx661, eps=EPS_DAMPENING),
                s_raw_canx_ko=math.log2(k_vx445 / k_vx661),
                split=str(split_info["split"]),
                exposed_pred0010_forced=bool(split_info["exposed_pred0010_forced"]),
            )
        )

    if aligned_identity_count != len(by_row_split):
        raise ValueError(
            f"Dual-sheet aligned variant count {aligned_identity_count} != manifest count {len(by_row_split)}"
        )
    if len(excluded) > MAX_EXCLUSIONS:
        raise RuntimeError(
            f"Non-positive/missing count {len(excluded)} exceeds pre-registered limit MAX_EXCLUSIONS={MAX_EXCLUSIONS}"
        )
    return aligned_rows, excluded


def _pearson_r(xs: Sequence[float], ys: Sequence[float]) -> float:
    n = len(xs)
    if n < 2:
        raise ValueError("Pearson r requires at least 2 paired observations")
    mx = sum(xs) / n
    my = sum(ys) / n
    num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    dx = sum((x - mx) ** 2 for x in xs)
    dy = sum((y - my) ** 2 for y in ys)
    if dx <= 0.0 or dy <= 0.0:
        raise ValueError("Pearson r is undefined for zero-variance input")
    return num / math.sqrt(dx * dy)


def _solve_linear_system(a: List[List[float]], b: List[float]) -> List[float]:
    """Gaussian elimination with partial pivoting."""
    n = len(b)
    mat = [row[:] + [b[i]] for i, row in enumerate(a)]
    for col in range(n):
        pivot = max(range(col, n), key=lambda r: abs(mat[r][col]))
        if abs(mat[pivot][col]) < 1e-14:
            raise ValueError("Singular linear system")
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


def _invert_matrix(a: List[List[float]]) -> List[List[float]]:
    n = len(a)
    inv_cols = []
    for k in range(n):
        e = [1.0 if i == k else 0.0 for i in range(n)]
        inv_cols.append(_solve_linear_system(a, e))
    return [[inv_cols[j][i] for j in range(n)] for i in range(n)]


def fit_unpenalized_ols(x_rows: Sequence[Sequence[float]], y_vec: Sequence[float]) -> Tuple[List[float], List[float]]:
    """Fit unpenalized OLS and return `(coefficients, conventional_ols_standard_errors)`."""
    n = len(x_rows)
    p = len(x_rows[0])
    if n <= p:
        raise ValueError(f"Need n > p for OLS, got n={n}, p={p}")
    xtx = [[sum(x_rows[i][a] * x_rows[i][b] for i in range(n)) for b in range(p)] for a in range(p)]
    xty = [sum(x_rows[i][a] * y_vec[i] for i in range(n)) for a in range(p)]
    beta = _solve_linear_system(xtx, xty)
    rss = sum((y_vec[i] - sum(beta[j] * x_rows[i][j] for j in range(p))) ** 2 for i in range(n))
    sigma2 = rss / (n - p)
    xtx_inv = _invert_matrix(xtx)
    ses = [math.sqrt(max(0.0, sigma2 * xtx_inv[j][j])) for j in range(p)]
    return beta, ses


def _extra_onehot(sub_domain: str, variant_class: str) -> List[float]:
    return [1.0 if sub_domain == s else 0.0 for s in EXTRA_SUBDOMAINS] + [
        1.0 if variant_class == v else 0.0 for v in EXTRA_CLASSES
    ]


@dataclass(frozen=True)
class FittedPartialRidgeModel:
    alpha: float
    b0_intercept: float
    b1_dmso_slope: float
    extra_means: Tuple[float, ...]
    extra_stds: Tuple[float, ...]
    extra_weights: Tuple[float, ...]

    def predict(self, x_dmso_log2: float, sub_domain: str, variant_class: str) -> float:
        raw = _extra_onehot(sub_domain, variant_class)
        total = self.b0_intercept + self.b1_dmso_slope * x_dmso_log2
        for w, val, mean, std in zip(self.extra_weights, raw, self.extra_means, self.extra_stds):
            if std > 0.0:
                total += w * ((val - mean) / std)
        return total


def fit_partial_ridge_on_arrays(
    xs_dmso: Sequence[float],
    extra_raw: Sequence[Sequence[float]],
    ys: Sequence[float],
    alpha: float,
) -> FittedPartialRidgeModel:
    """Fit partial-penalized Ridge with 0 diagonal penalty on intercept and `x_dmso`."""
    n = len(ys)
    if n < 5:
        raise ValueError("Need at least 5 rows for partial-penalized Ridge")
    if alpha < 0.0:
        raise ValueError(f"alpha must be >= 0, got {alpha}")
    p_extra = len(extra_raw[0])
    means = [sum(row[j] for row in extra_raw) / n for j in range(p_extra)]
    stds: List[float] = []
    for j in range(p_extra):
        var = sum((row[j] - means[j]) ** 2 for row in extra_raw) / n
        stds.append(math.sqrt(var) if var > 1e-12 else 0.0)

    p = 2 + p_extra
    x_mat = [
        [1.0, xs_dmso[i]]
        + [(extra_raw[i][j] - means[j]) / stds[j] if stds[j] > 0.0 else 0.0 for j in range(p_extra)]
        for i in range(n)
    ]
    xtx = [[sum(x_mat[i][a] * x_mat[i][b] for i in range(n)) for b in range(p)] for a in range(p)]
    xty = [sum(x_mat[i][a] * ys[i] for i in range(n)) for a in range(p)]
    for j in range(2, p):
        xtx[j][j] += alpha if stds[j - 2] > 0.0 else 1.0

    beta = _solve_linear_system(xtx, xty)
    return FittedPartialRidgeModel(
        alpha=alpha,
        b0_intercept=beta[0],
        b1_dmso_slope=beta[1],
        extra_means=tuple(means),
        extra_stds=tuple(stds),
        extra_weights=tuple(beta[2:]),
    )


def _safe_rel_mae_reduction(baseline_mae: float, candidate_mae: float) -> float:
    if baseline_mae <= 0.0:
        raise ValueError(f"Baseline MAE must be strictly > 0, got {baseline_mae}")
    return (baseline_mae - candidate_mae) / baseline_mae


def residue_clustered_bootstrap_delta_cis(
    codons: Sequence[int],
    targets: Sequence[float],
    baseline_preds: Sequence[float],
    cand_preds: Sequence[float],
    *,
    resamples: int = BOOTSTRAP_RESAMPLES,
    seed: int = BOOTSTRAP_SEED,
) -> Dict[str, Any]:
    """Compute residue-clustered 95% and Bonferroni 97.5% two-sided CIs on MAE_base - MAE_cand."""
    by_codon: Dict[int, List[Tuple[float, float]]] = {}
    for c, y, bp, cp in zip(codons, targets, baseline_preds, cand_preds):
        by_codon.setdefault(c, []).append((abs(y - bp), abs(y - cp)))
    clusters = list(by_codon.values())
    n_clusters = len(clusters)
    if n_clusters < 2:
        raise ValueError("Need at least 2 clusters for bootstrap")

    rng = random.Random(seed)
    deltas: List[float] = []
    for _ in range(resamples):
        sum_b = 0.0
        sum_c = 0.0
        cnt = 0
        for _k in range(n_clusters):
            cl = clusters[rng.randrange(n_clusters)]
            for eb, ec in cl:
                sum_b += eb
                sum_c += ec
                cnt += 1
        deltas.append((sum_b - sum_c) / cnt)

    deltas.sort()
    lo95 = deltas[max(0, int(0.025 * resamples))]
    hi95 = deltas[min(resamples - 1, int(0.975 * resamples))]
    lo975 = deltas[max(0, int(0.0125 * resamples))]
    hi975 = deltas[min(resamples - 1, int(0.9875 * resamples))]
    return {
        "n_holdout_clusters": n_clusters,
        "resamples": resamples,
        "seed": seed,
        "delta_mae_ci95": [round(lo95, 6), round(hi95, 6)],
        "ci95_low_gt_zero": lo95 > 0.0,
        "delta_mae_ci975_bonferroni": [round(lo975, 6), round(hi975, 6)],
        "ci975_low_gt_zero": lo975 > 0.0,
    }


def _gate_outcome(rel_reduction: float, ci975_low_gt_zero: bool, pass_label: str) -> Tuple[str, str]:
    if rel_reduction >= GATE_REL_MAE_REDUCTION and ci975_low_gt_zero:
        return "open", pass_label
    if rel_reduction >= GATE_REL_MAE_REDUCTION:
        return "refuted", "INCONCLUSIVE_CI975_STRADDLES_ZERO"
    return "refuted", "REFUTED_BELOW_10PCT_MAE_REDUCTION_THRESHOLD"


def _dev_fold_map(dev_rows: Sequence[DualSheetVariant]) -> Dict[int, int]:
    sorted_codons = sorted({r.codon_position for r in dev_rows})
    return {c: idx % CV_FOLDS for idx, c in enumerate(sorted_codons)}


def run_pre_holdout_controls(
    dev_rows: Sequence[DualSheetVariant],
    aligned_total: int,
    excluded_count: int,
) -> Dict[str, Any]:
    """Run all pre-holdout positive & shuffled-predictor negative controls strictly on `dev`."""
    canx_r = _pearson_r([r.x_dmso_log2 for r in dev_rows], [r.x_canx_ko_dmso_log2 for r in dev_rows])
    if canx_r <= CANX_STABILITY_MIN_R:
        raise RuntimeError(f"CANX KO DMSO stability control failed on dev: r={canx_r} <= {CANX_STABILITY_MIN_R}")

    # Monotonic pharmacological directionality in dev (Parental + CANX KO)
    p_med_dmso = statistics.median([r.parental_dmso_rel for r in dev_rows])
    p_med_445 = statistics.median([r.parental_vx445_rel for r in dev_rows])
    p_med_combo = statistics.median([r.parental_combo_rel for r in dev_rows])
    k_med_dmso = statistics.median([r.canx_ko_dmso_rel for r in dev_rows])
    k_med_661 = statistics.median([r.canx_ko_vx661_rel for r in dev_rows])
    k_med_445 = statistics.median([r.canx_ko_vx445_rel for r in dev_rows])
    k_med_combo = statistics.median([r.canx_ko_combo_rel for r in dev_rows])

    par_dir_ok = p_med_combo > p_med_445 > p_med_dmso
    ko_dir_ok = (k_med_combo > k_med_445 > k_med_dmso) and (k_med_combo > k_med_661 > k_med_dmso)
    if not (par_dir_ok and ko_dir_ok):
        raise RuntimeError(
            f"Pharmacological directionality control failed on dev: "
            f"parental=({p_med_dmso},{p_med_445},{p_med_combo}) "
            f"canx_ko=({k_med_dmso},{k_med_661},{k_med_445},{k_med_combo})"
        )

    # Shuffled-predictor negative controls strictly inside dev 5-fold CV
    fold_map = _dev_fold_map(dev_rows)
    rng = random.Random(SHUFFLE_SEED)
    shuf_canx = [r.x_canx_ko_dmso_log2 for r in dev_rows]
    rng.shuffle(shuf_canx)
    shuf_feats = [_extra_onehot(r.sub_domain, r.variant_class) for r in dev_rows]
    rng.shuffle(shuf_feats)

    arm_a_base_errs: List[float] = []
    arm_a_shuf_errs: List[float] = []
    arm_b_med_errs: List[float] = []
    arm_b_ols_errs: List[float] = []
    arm_b_shuf_errs: List[float] = []

    for fold in range(CV_FOLDS):
        tr_idx = [i for i, r in enumerate(dev_rows) if fold_map[r.codon_position] != fold]
        va_idx = [i for i, r in enumerate(dev_rows) if fold_map[r.codon_position] == fold]

        # Arm A shuffled CANX_KO_DMSO vs 1D Parental dmso_ols
        beta_a_base, _ = fit_unpenalized_ols(
            [[1.0, dev_rows[i].x_dmso_log2] for i in tr_idx],
            [dev_rows[i].y_combo_log2 for i in tr_idx],
        )
        beta_a_shuf, _ = fit_unpenalized_ols(
            [[1.0, dev_rows[i].x_dmso_log2, shuf_canx[i]] for i in tr_idx],
            [dev_rows[i].y_combo_log2 for i in tr_idx],
        )
        for i in va_idx:
            r = dev_rows[i]
            arm_a_base_errs.append(abs(r.y_combo_log2 - (beta_a_base[0] + beta_a_base[1] * r.x_dmso_log2)))
            arm_a_shuf_errs.append(
                abs(r.y_combo_log2 - (beta_a_shuf[0] + beta_a_shuf[1] * r.x_dmso_log2 + beta_a_shuf[2] * shuf_canx[i]))
            )

        # Arm B shuffled domain/class on CANX_KO_DMSO vs true CANX KO baselines
        tr_med_s = statistics.median([dev_rows[i].s_eps_canx_ko for i in tr_idx])
        beta_b_ols, _ = fit_unpenalized_ols(
            [[1.0, dev_rows[i].x_canx_ko_dmso_log2] for i in tr_idx],
            [dev_rows[i].s_eps_canx_ko for i in tr_idx],
        )
        mod_b_shuf = fit_partial_ridge_on_arrays(
            [dev_rows[i].x_canx_ko_dmso_log2 for i in tr_idx],
            [shuf_feats[i] for i in tr_idx],
            [dev_rows[i].s_eps_canx_ko for i in tr_idx],
            alpha=10.0,
        )
        for i in va_idx:
            r = dev_rows[i]
            arm_b_med_errs.append(abs(r.s_eps_canx_ko - tr_med_s))
            arm_b_ols_errs.append(abs(r.s_eps_canx_ko - (beta_b_ols[0] + beta_b_ols[1] * r.x_canx_ko_dmso_log2)))
            raw_sh = shuf_feats[i]
            pred_sh = mod_b_shuf.b0_intercept + mod_b_shuf.b1_dmso_slope * r.x_canx_ko_dmso_log2 + sum(
                w * ((val - m) / sd) if sd > 0.0 else 0.0
                for w, val, m, sd in zip(mod_b_shuf.extra_weights, raw_sh, mod_b_shuf.extra_means, mod_b_shuf.extra_stds)
            )
            arm_b_shuf_errs.append(abs(r.s_eps_canx_ko - pred_sh))

    a_base_mae = sum(arm_a_base_errs) / len(arm_a_base_errs)
    a_shuf_mae = sum(arm_a_shuf_errs) / len(arm_a_shuf_errs)
    a_shuf_rel = _safe_rel_mae_reduction(a_base_mae, a_shuf_mae)

    b_base_mae = min(sum(arm_b_med_errs) / len(arm_b_med_errs), sum(arm_b_ols_errs) / len(arm_b_ols_errs))
    b_shuf_mae = sum(arm_b_shuf_errs) / len(arm_b_shuf_errs)
    b_shuf_rel = _safe_rel_mae_reduction(b_base_mae, b_shuf_mae)

    if a_shuf_rel >= GATE_REL_MAE_REDUCTION or b_shuf_rel >= GATE_REL_MAE_REDUCTION:
        raise RuntimeError(
            f"Shuffled-predictor negative control exceeded 10% threshold: armA={a_shuf_rel}, armB={b_shuf_rel}"
        )

    return {
        "dual_sheet_alignment": {
            "aligned_variants": aligned_total,
            "excluded_count": excluded_count,
            "passed": aligned_total >= 230 and excluded_count <= MAX_EXCLUSIONS,
        },
        "dev_canx_ko_dmso_stability": {
            "pearson_r_log2": round(canx_r, 6),
            "threshold": CANX_STABILITY_MIN_R,
            "passed": canx_r > CANX_STABILITY_MIN_R,
        },
        "dev_pharmacological_directionality": {
            "parental_medians_rel": {
                "dmso": round(p_med_dmso, 6),
                "vx445": round(p_med_445, 6),
                "combo": round(p_med_combo, 6),
            },
            "canx_ko_medians_rel": {
                "dmso": round(k_med_dmso, 6),
                "vx661": round(k_med_661, 6),
                "vx445": round(k_med_445, 6),
                "combo": round(k_med_combo, 6),
            },
            "passed": par_dir_ok and ko_dir_ok,
        },
        "dev_shuffled_negative_controls": {
            "seed": SHUFFLE_SEED,
            "arm_a_shuffled_canx_dev_cv_rel_reduction": round(a_shuf_rel, 6),
            "arm_b_shuffled_features_dev_cv_rel_reduction": round(b_shuf_rel, 6),
            "passed": True,
        },
    }


def evaluate_canx_selectivity_benchmark(
    data_dir: Path,
    repo_root: Path = REPO_ROOT,
) -> Dict[str, Any]:
    """Run dual-sheet loading, pre-holdout controls, and one-shot evaluation of Arm A & Arm B."""
    master_xlsx = data_dir / MASTER_XLSX_NAME
    master_sha = hashlib.sha256(master_xlsx.read_bytes()).hexdigest()
    if master_sha != MASTER_XLSX_SHA256:
        raise ValueError(f"Master xlsx SHA256 mismatch: {master_sha} != {MASTER_XLSX_SHA256}")

    manifest, by_row_split, manifest_sha = load_locked_split_manifest(repo_root)
    rows, excluded = load_dual_sheet_tedman_variants(master_xlsx, by_row_split)
    dev_rows = [r for r in rows if r.split == "dev"]
    holdout_rows = [r for r in rows if r.split == "holdout"]
    dev_codons = {r.codon_position for r in dev_rows}
    holdout_codons = {r.codon_position for r in holdout_rows}

    if dev_codons & holdout_codons:
        raise RuntimeError("Codon overlap between dev and holdout")
    if len(holdout_codons) < MIN_HOLDOUT_CLUSTERS:
        raise RuntimeError(f"Insufficient holdout clusters: {len(holdout_codons)}")

    controls = run_pre_holdout_controls(dev_rows, len(rows), len(excluded))
    fold_map = _dev_fold_map(dev_rows)

    # ---------------- Arm A: Calnexin chaperone coupling on y = log2(Parental VX-661+VX-445)
    beta_a_base, se_a_base = fit_unpenalized_ols(
        [[1.0, r.x_dmso_log2] for r in dev_rows],
        [r.y_combo_log2 for r in dev_rows],
    )
    beta_a_cand, se_a_cand = fit_unpenalized_ols(
        [[1.0, r.x_dmso_log2, r.x_canx_ko_dmso_log2] for r in dev_rows],
        [r.y_combo_log2 for r in dev_rows],
    )

    a_base_hold_preds = [beta_a_base[0] + beta_a_base[1] * r.x_dmso_log2 for r in holdout_rows]
    a_cand_hold_preds = [
        beta_a_cand[0] + beta_a_cand[1] * r.x_dmso_log2 + beta_a_cand[2] * r.x_canx_ko_dmso_log2
        for r in holdout_rows
    ]
    a_base_mae = sum(abs(r.y_combo_log2 - p) for r, p in zip(holdout_rows, a_base_hold_preds)) / len(holdout_rows)
    a_cand_mae = sum(abs(r.y_combo_log2 - p) for r, p in zip(holdout_rows, a_cand_hold_preds)) / len(holdout_rows)
    a_base_mae_lin = sum(abs(r.parental_combo_rel - (2.0 ** p)) for r, p in zip(holdout_rows, a_base_hold_preds)) / len(holdout_rows)
    a_cand_mae_lin = sum(abs(r.parental_combo_rel - (2.0 ** p)) for r, p in zip(holdout_rows, a_cand_hold_preds)) / len(holdout_rows)
    a_rel_red = _safe_rel_mae_reduction(a_base_mae, a_cand_mae)
    a_boot = residue_clustered_bootstrap_delta_cis(
        [r.codon_position for r in holdout_rows],
        [r.y_combo_log2 for r in holdout_rows],
        a_base_hold_preds,
        a_cand_hold_preds,
    )
    a_standing, a_verdict = _gate_outcome(
        a_rel_red, bool(a_boot["ci975_low_gt_zero"]), "PASS_CANX_COUPLING_BEATS_DMSO_OLS"
    )

    # ---------------- Arm B: Dampened CANX KO Single-Corrector Selectivity S_eps_canx_ko
    cv_med_errs: List[float] = []
    cv_ols_errs: List[float] = []
    for fold in range(CV_FOLDS):
        tr = [r for r in dev_rows if fold_map[r.codon_position] != fold]
        va = [r for r in dev_rows if fold_map[r.codon_position] == fold]
        tr_med = statistics.median([r.s_eps_canx_ko for r in tr])
        fb, _ = fit_unpenalized_ols([[1.0, r.x_canx_ko_dmso_log2] for r in tr], [r.s_eps_canx_ko for r in tr])
        for r in va:
            cv_med_errs.append(abs(r.s_eps_canx_ko - tr_med))
            cv_ols_errs.append(abs(r.s_eps_canx_ko - (fb[0] + fb[1] * r.x_canx_ko_dmso_log2)))

    cv_med_mae = sum(cv_med_errs) / len(cv_med_errs)
    cv_ols_mae = sum(cv_ols_errs) / len(cv_ols_errs)
    locked_b_baseline = "constant_dev_median" if cv_med_mae <= cv_ols_mae else "canx_ko_dmso_ols"

    cv_path_b: List[Dict[str, Any]] = []
    best_alpha_b = ALPHA_GRID[0]
    best_cv_mae_b = float("inf")
    for alpha in ALPHA_GRID:
        errs: List[float] = []
        for fold in range(CV_FOLDS):
            tr = [r for r in dev_rows if fold_map[r.codon_position] != fold]
            va = [r for r in dev_rows if fold_map[r.codon_position] == fold]
            mod = fit_partial_ridge_on_arrays(
                [r.x_canx_ko_dmso_log2 for r in tr],
                [_extra_onehot(r.sub_domain, r.variant_class) for r in tr],
                [r.s_eps_canx_ko for r in tr],
                alpha=alpha,
            )
            for r in va:
                errs.append(abs(r.s_eps_canx_ko - mod.predict(r.x_canx_ko_dmso_log2, r.sub_domain, r.variant_class)))
        mae_cv = sum(errs) / len(errs)
        cv_path_b.append({"alpha": alpha, "dev_cv_mae_s_eps": round(mae_cv, 6)})
        if mae_cv < best_cv_mae_b:
            best_cv_mae_b = mae_cv
            best_alpha_b = alpha

    dev_med_s_eps = statistics.median([r.s_eps_canx_ko for r in dev_rows])
    beta_b_ols, se_b_ols = fit_unpenalized_ols(
        [[1.0, r.x_canx_ko_dmso_log2] for r in dev_rows],
        [r.s_eps_canx_ko for r in dev_rows],
    )
    final_cand_b = fit_partial_ridge_on_arrays(
        [r.x_canx_ko_dmso_log2 for r in dev_rows],
        [_extra_onehot(r.sub_domain, r.variant_class) for r in dev_rows],
        [r.s_eps_canx_ko for r in dev_rows],
        alpha=best_alpha_b,
    )

    b_med_hold_preds = [dev_med_s_eps] * len(holdout_rows)
    b_ols_hold_preds = [beta_b_ols[0] + beta_b_ols[1] * r.x_canx_ko_dmso_log2 for r in holdout_rows]
    b_cand_hold_preds = [final_cand_b.predict(r.x_canx_ko_dmso_log2, r.sub_domain, r.variant_class) for r in holdout_rows]

    b_med_hold_mae = sum(abs(r.s_eps_canx_ko - p) for r, p in zip(holdout_rows, b_med_hold_preds)) / len(holdout_rows)
    b_ols_hold_mae = sum(abs(r.s_eps_canx_ko - p) for r, p in zip(holdout_rows, b_ols_hold_preds)) / len(holdout_rows)
    b_cand_hold_mae = sum(abs(r.s_eps_canx_ko - p) for r, p in zip(holdout_rows, b_cand_hold_preds)) / len(holdout_rows)

    locked_preds = b_med_hold_preds if locked_b_baseline == "constant_dev_median" else b_ols_hold_preds
    locked_mae = b_med_hold_mae if locked_b_baseline == "constant_dev_median" else b_ols_hold_mae
    b_rel_red = _safe_rel_mae_reduction(locked_mae, b_cand_hold_mae)
    b_boot = residue_clustered_bootstrap_delta_cis(
        [r.codon_position for r in holdout_rows],
        [r.s_eps_canx_ko for r in holdout_rows],
        locked_preds,
        b_cand_hold_preds,
    )
    b_standing, b_verdict = _gate_outcome(
        b_rel_red, bool(b_boot["ci975_low_gt_zero"]), "PASS_SELECTIVITY_RIDGE_BEATS_LOCKED_BASELINE"
    )

    # Context only (never gated): subdomain S_eps medians and raw S MAE
    subdomain_medians_dev = {
        sd: round(statistics.median([r.s_eps_canx_ko for r in dev_rows if r.sub_domain == sd]), 6)
        for sd in SUBDOMAIN_CATEGORIES
        if any(r.sub_domain == sd for r in dev_rows)
    }
    raw_s_med_dev = statistics.median([r.s_raw_canx_ko for r in dev_rows])
    raw_s_hold_med_mae = sum(abs(r.s_raw_canx_ko - raw_s_med_dev) for r in holdout_rows) / len(holdout_rows)

    return {
        "schema_version": "tedman2025-canx-selectivity-results-v1",
        "benchmark_id": "tedman2025-canx-dependency-and-corrector-selectivity",
        "evaluated_on_local": "2026-11-23",
        "evaluated_on_utc": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "done_sha256": DONE_SHA256_EXPECTED,
        "split_manifest_path": str(SPLIT_MANIFEST_REL).replace("\\", "/"),
        "split_manifest_sha256": manifest_sha,
        "multiplicity_policy": "bonferroni_2_arms_975pct_two_sided_ci_no_or_gate",
        "source_workbook_note": (
            "Parental All_variants Col J (VX-661 Mut/WT) is blank (0/234) in the Tedman et al. 2025 "
            "workbook; All_variants CANX KO Cols F/J/N/R are 234/234 populated."
        ),
        "dataset": {
            "source_id": "cftr-dms-tedman-2025",
            "source_master_xlsx": MASTER_XLSX_NAME,
            "source_master_sha256": master_sha,
            "scale_boundary": "HEK293T plasma-membrane expression is not CFTR chloride conductance or clinical response.",
            "aligned_variants_total": len(rows),
            "excluded_count": len(excluded),
            "excluded_rows": excluded,
            "dev_variants": len(dev_rows),
            "holdout_variants": len(holdout_rows),
            "dev_unique_codons": len(dev_codons),
            "holdout_unique_codons": len(holdout_codons),
            "dev_pred0010_forced_variants": sum(1 for r in dev_rows if r.exposed_pred0010_forced),
            "holdout_pred0010_overlap_codons": 0,
        },
        "pre_holdout_controls": controls,
        "arm_a_canx_chaperone_coupling": {
            "target": "y_combo_log2 = log2(Parental VX-661+VX-445)",
            "dev_models": {
                "baseline_dmso_ols": {
                    "b0_intercept": round(beta_a_base[0], 6),
                    "b1_parental_dmso_slope": round(beta_a_base[1], 6),
                    "se": [round(s, 6) for s in se_a_base],
                },
                "candidate_dmso_plus_canx_ko_ols": {
                    "b0_intercept": round(beta_a_cand[0], 6),
                    "b1_parental_dmso_slope": round(beta_a_cand[1], 6),
                    "b2_canx_ko_dmso_slope": round(beta_a_cand[2], 6),
                    "se": [round(s, 6) for s in se_a_cand],
                },
            },
            "holdout_evaluation": {
                "baseline_dmso_ols_mae_log2": round(a_base_mae, 6),
                "baseline_dmso_ols_mae_linear_context_only": round(a_base_mae_lin, 6),
                "candidate_mae_log2": round(a_cand_mae, 6),
                "candidate_mae_linear_context_only": round(a_cand_mae_lin, 6),
                "relative_mae_reduction": round(a_rel_red, 6),
                "pre_registered_gate_threshold": GATE_REL_MAE_REDUCTION,
                "residue_clustered_bootstrap": a_boot,
                "gate_verdict": a_verdict,
                "ledger_standing": a_standing,
            },
        },
        "arm_b_canx_ko_single_corrector_selectivity": {
            "target": "S_eps_canx_ko = log2((CANX_KO_VX445_rel + 0.05) / (CANX_KO_VX661_rel + 0.05))",
            "eps_dampening": EPS_DAMPENING,
            "dev_models": {
                "constant_dev_median_s_eps": round(dev_med_s_eps, 6),
                "canx_ko_dmso_ols": {
                    "b0_intercept": round(beta_b_ols[0], 6),
                    "b1_canx_ko_dmso_slope": round(beta_b_ols[1], 6),
                    "se": [round(s, 6) for s in se_b_ols],
                },
                "dev_cv_5fold_baselines": {
                    "constant_dev_median_mae_s_eps": round(cv_med_mae, 6),
                    "canx_ko_dmso_ols_mae_s_eps": round(cv_ols_mae, 6),
                    "locked_baseline_before_holdout": locked_b_baseline,
                },
                "candidate_partial_ridge": {
                    "selected_alpha": best_alpha_b,
                    "cv_path": cv_path_b,
                    "b0_intercept": round(final_cand_b.b0_intercept, 6),
                    "b1_canx_ko_dmso_slope": round(final_cand_b.b1_dmso_slope, 6),
                    "extra_feature_names": list(EXTRA_SUBDOMAINS) + list(EXTRA_CLASSES),
                    "extra_weights": [round(w, 6) for w in final_cand_b.extra_weights],
                },
                "dev_subdomain_s_eps_medians_context_only": subdomain_medians_dev,
            },
            "holdout_evaluation": {
                "locked_baseline_name": locked_b_baseline,
                "locked_baseline_mae_s_eps": round(locked_mae, 6),
                "constant_dev_median_mae_s_eps": round(b_med_hold_mae, 6),
                "canx_ko_dmso_ols_mae_s_eps": round(b_ols_hold_mae, 6),
                "candidate_partial_ridge_mae_s_eps": round(b_cand_hold_mae, 6),
                "s_raw_constant_median_mae_context_only": round(raw_s_hold_med_mae, 6),
                "relative_mae_reduction": round(b_rel_red, 6),
                "pre_registered_gate_threshold": GATE_REL_MAE_REDUCTION,
                "residue_clustered_bootstrap": b_boot,
                "gate_verdict": b_verdict,
                "ledger_standing": b_standing,
            },
        },
    }


def verify_checked_in_artifacts(repo_root: Path = REPO_ROOT, data_dir: Optional[Path] = None) -> Dict[str, Any]:
    """Verify checked-in results JSON, done/manifest hashes, zero leakage, and raw parity if present."""
    bench_dir = repo_root / "benchmarks" / "cftr_dms_canx_selectivity"
    done_path = bench_dir / "done.md"
    results_path = bench_dir / RESULTS_FILE_NAME

    done_sha = hashlib.sha256(
        done_path.read_text(encoding="utf-8", errors="replace").replace("\r\n", "\n").encode("utf-8")
    ).hexdigest()
    results = json.loads(results_path.read_text(encoding="utf-8"))
    manifest, _, manifest_sha = load_locked_split_manifest(repo_root)

    if results["done_sha256"] != done_sha:
        raise AssertionError(f"done_sha256 mismatch: {results['done_sha256']} != {done_sha}")
    if results["split_manifest_sha256"] != manifest_sha:
        raise AssertionError("split_manifest_sha256 mismatch.")

    resolved_data = resolve_tedman_data_dir(data_dir, repo_root=repo_root)
    raw_parity_checked = False
    if resolved_data is not None:
        fresh = evaluate_canx_selectivity_benchmark(resolved_data, repo_root=repo_root)
        arm_a_diff = abs(
            fresh["arm_a_canx_chaperone_coupling"]["holdout_evaluation"]["candidate_mae_log2"]
            - results["arm_a_canx_chaperone_coupling"]["holdout_evaluation"]["candidate_mae_log2"]
        )
        arm_b_diff = abs(
            fresh["arm_b_canx_ko_single_corrector_selectivity"]["holdout_evaluation"]["candidate_partial_ridge_mae_s_eps"]
            - results["arm_b_canx_ko_single_corrector_selectivity"]["holdout_evaluation"]["candidate_partial_ridge_mae_s_eps"]
        )
        if arm_a_diff > 1e-6 or arm_b_diff > 1e-6:
            raise AssertionError(f"Raw .xlsx parity drift: arm_a_diff={arm_a_diff}, arm_b_diff={arm_b_diff}")
        raw_parity_checked = True

    return {
        "status": "PASS",
        "raw_local_data_parity_checked": raw_parity_checked,
        "dev_variants": manifest["dev_variants"],
        "holdout_variants": manifest["holdout_variants"],
        "holdout_unique_codons": manifest["holdout_unique_codons"],
        "arm_a_verdict": results["arm_a_canx_chaperone_coupling"]["holdout_evaluation"]["gate_verdict"],
        "arm_b_verdict": results["arm_b_canx_ko_single_corrector_selectivity"]["holdout_evaluation"]["gate_verdict"],
    }


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--data-dir", type=Path, default=None, help="Optional path to local tedman-2025 folder.")
    parser.add_argument("--write-artifacts", action="store_true", help="Write results JSON once.")
    parser.add_argument("--check", action="store_true", help="Verify checked-in artifacts and local raw parity.")
    parser.add_argument("--json", action="store_true", help="Print strict JSON output.")
    args = parser.parse_args(argv)

    if args.check:
        summary = verify_checked_in_artifacts(repo_root=REPO_ROOT, data_dir=args.data_dir)
        if args.json:
            print(json.dumps(summary, indent=2, allow_nan=False))
        else:
            print(
                f"PASS cftr_dms_canx_selectivity: dev={summary['dev_variants']}, "
                f"holdout={summary['holdout_variants']} (codons={summary['holdout_unique_codons']}), "
                f"raw_parity={summary['raw_local_data_parity_checked']}, "
                f"arm_a={summary['arm_a_verdict']}, arm_b={summary['arm_b_verdict']}"
            )
        return 0

    resolved_data = resolve_tedman_data_dir(args.data_dir, repo_root=REPO_ROOT)
    if resolved_data is None:
        print("Local Tedman 2025 files not present under sources/downloads/tedman-2025/.", file=sys.stderr)
        return 1

    results = evaluate_canx_selectivity_benchmark(resolved_data, repo_root=REPO_ROOT)
    if args.write_artifacts:
        (HERE / RESULTS_FILE_NAME).write_text(
            json.dumps(results, indent=2, allow_nan=False) + "\n", encoding="utf-8"
        )

    if args.json:
        print(json.dumps(results, indent=2, allow_nan=False))
    else:
        ha = results["arm_a_canx_chaperone_coupling"]["holdout_evaluation"]
        hb = results["arm_b_canx_ko_single_corrector_selectivity"]["holdout_evaluation"]
        print(
            f"Arm A: verdict={ha['gate_verdict']} rel_red={ha['relative_mae_reduction']} | "
            f"Arm B: verdict={hb['gate_verdict']} rel_red={hb['relative_mae_reduction']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
