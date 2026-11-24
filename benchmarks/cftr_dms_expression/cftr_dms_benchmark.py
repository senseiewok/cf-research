#!/usr/bin/env python3
"""Quantitative CFTR plasma-membrane expression benchmark on Tedman et al. 2025 (Council Item #1).

Evaluates whether baseline untreated CFTR expression (DMSO), sub-domain membership, and variant
consequence class predict modulator-rescued plasma-membrane expression (VX-661+VX-445 Mut/WT on a
multiplicative log2 scale) on an unexposed, residue-clustered held-out variant panel.

Pure Python 3.10+ standard library only; zero network calls; raw downloaded `.xlsx`/`.tsv` files
stay git-ignored under `sources/downloads/tedman-2025/`.
Research benchmark only: plasma-membrane expression is not channel conductance or clinical benefit.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
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
AGM_ALT_XLSX_NAME = "AGM ALT Library Comparison.xlsx"
AGM_ALT_SHA256 = "6673ab2a480200e7d36133b98dfbeb4864d1a41aa22cf4fdddea59eda8dbe178"
DONE_SHA256_EXPECTED = "d4b8e94679660bf21f7aa55c023802049c16b2397edb1bb31a191a51fe57d39d"

THERATYPE_FIXTURE_REL = Path("benchmarks/theratype/theratype_open_annotations_2026-11-21.json")
RESULTS_FILE_NAME = "results_2026-11-23.json"
MANIFEST_FILE_NAME = "split_manifest.json"

SUBDOMAIN_CATEGORIES: Tuple[str, ...] = ("Lasso", "TMD1", "NBD1", "R Domain", "TMD2", "NBD2")
EXTRA_SUBDOMAINS: Tuple[str, ...] = ("TMD1", "NBD1", "R Domain", "TMD2", "NBD2")
VARIANT_CLASSES: Tuple[str, ...] = ("missense", "nonsense", "inframe_indel", "other")
EXTRA_CLASSES: Tuple[str, ...] = ("nonsense", "inframe_indel", "other")

ALPHA_GRID: Tuple[float, ...] = (0.01, 0.1, 1.0, 10.0, 100.0, 1000.0)
CV_FOLDS = 5
BOOTSTRAP_SEED = 20261123
BOOTSTRAP_RESAMPLES = 2000
GATE_REL_MAE_REDUCTION = 0.10
MIN_HOLDOUT_CLUSTERS = 40

XLSX_NS = {"s": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
CELL_REF_RE = re.compile(r"^([A-Z]+)(\d+)$")
RANGE_REF_RE = re.compile(r"^([A-Z]+)(\d+):([A-Z]+)(\d+)$")


@dataclass(frozen=True)
class VariantRecord:
    """Single parsed variant row from `All_variants` (log2-transformed ratio scale)."""

    row_num: int
    sub_domain: str
    nucleic_acid_mutation: str
    protein_distinction: str
    codon_position: int
    variant_class: str
    dmso_rel_mut_wt: float
    combo_rel_mut_wt: float
    x_dmso_log2: float
    y_combo_log2: float
    split: str = ""
    exposed_pred0010_forced: bool = False


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
    """Parse sparse `.xlsx` worksheet XML by exact `(row, col)` coordinates without external packages."""
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
                if m:
                    col = m.group(1)
                    row_idx = int(m.group(2))
                else:
                    raise ValueError(f"Missing/invalid sparse coordinate r={coord!r} in {sheet_name!r}")
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
    """Classify `Protein distinction` string into `(codon_position, variant_class)`."""
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
        if alt == "=" or alt == ref:
            return codon, "other"
        return codon, "missense"

    return None, "other"


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


def load_exposed_pred0010_codons(repo_root: Path = REPO_ROOT) -> Tuple[int, ...]:
    """Load the 30 exposed `PRED-0010` codons from `benchmarks/theratype/` to force them out of holdout."""
    path = repo_root / THERATYPE_FIXTURE_REL
    payload = json.loads(path.read_text(encoding="utf-8"))
    codons = sorted({int(v["codon_position"]) for v in payload.get("variants") or []})
    return tuple(codons)


def load_tedman_variants(
    master_xlsx: Path,
    exposed_codons: Sequence[int],
) -> Tuple[List[VariantRecord], List[Dict[str, Any]]]:
    """Load valid `All_variants` parental DMS rows (`log2` scale) plus explicit exclusions."""
    cells, max_row = parse_xlsx_sheet(master_xlsx, "All_variants", expand_merged=True)
    exposed_set = set(int(c) for c in exposed_codons)

    records: List[VariantRecord] = []
    excluded: List[Dict[str, Any]] = []
    current_subdomain = ""

    for row_num in range(3, max_row + 1):
        sub_raw = cells.get((row_num, "A"), "").strip()
        if sub_raw:
            current_subdomain = sub_raw
        nuc = cells.get((row_num, "B"), "").strip()
        prot = cells.get((row_num, "C"), "").strip()
        if not nuc and not prot:
            continue

        if current_subdomain not in SUBDOMAIN_CATEGORIES:
            excluded.append({
                "row_num": row_num,
                "protein_distinction": prot,
                "nucleic_acid_mutation": nuc,
                "reason": f"unknown_subdomain:{current_subdomain!r}",
            })
            continue

        codon, vclass = parse_protein_distinction(prot)
        if codon is None:
            excluded.append({
                "row_num": row_num,
                "protein_distinction": prot,
                "nucleic_acid_mutation": nuc,
                "reason": "unparseable_codon",
            })
            continue

        dmso_rel = _parse_pos_float(cells.get((row_num, "F")))
        combo_rel = _parse_pos_float(cells.get((row_num, "R")))
        if dmso_rel is None or combo_rel is None:
            excluded.append({
                "row_num": row_num,
                "protein_distinction": prot,
                "nucleic_acid_mutation": nuc,
                "reason": "missing_or_nonpositive_dmso_or_combo_rel",
            })
            continue

        forced = codon in exposed_set
        if forced:
            split = "dev"
        else:
            digest = hashlib.sha256(f"tedman2025-item1:{codon}".encode("utf-8")).hexdigest()
            bucket = int(digest[:8], 16) % 100
            split = "dev" if bucket < 60 else "holdout"

        records.append(
            VariantRecord(
                row_num=row_num,
                sub_domain=current_subdomain,
                nucleic_acid_mutation=nuc,
                protein_distinction=prot,
                codon_position=codon,
                variant_class=vclass,
                dmso_rel_mut_wt=dmso_rel,
                combo_rel_mut_wt=combo_rel,
                x_dmso_log2=math.log2(dmso_rel),
                y_combo_log2=math.log2(combo_rel),
                split=split,
                exposed_pred0010_forced=forced,
            )
        )

    return records, excluded


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
        raise ValueError("Pearson r is undefined when either input vector has zero variance")
    return num / math.sqrt(dx * dy)


def _average_ranks(vals: Sequence[float]) -> List[float]:
    n = len(vals)
    indexed = sorted(range(n), key=lambda i: vals[i])
    ranks = [0.0] * n
    i = 0
    while i < n:
        j = i
        while j + 1 < n and vals[indexed[j + 1]] == vals[indexed[i]]:
            j += 1
        avg_rank = 0.5 * (i + j) + 1.0
        for k in range(i, j + 1):
            ranks[indexed[k]] = avg_rank
        i = j + 1
    return ranks


def _spearman_rho(xs: Sequence[float], ys: Sequence[float]) -> float:
    return _pearson_r(_average_ranks(xs), _average_ranks(ys))


def check_agm_alt_positive_control(agm_xlsx: Path) -> Dict[str, Any]:
    """Reproduce cross-study correlation between McKee 2023 (AGM) and Tedman 2025 (ALT)."""
    cells, max_row = parse_xlsx_sheet(agm_xlsx, "Sheet1", expand_merged=False)
    agm_vals: List[float] = []
    alt_vals: List[float] = []
    for row_num in range(1, max_row + 1):
        agm = _parse_pos_float(cells.get((row_num, "B")))
        alt = _parse_pos_float(cells.get((row_num, "D")))
        if agm is not None and alt is not None:
            agm_vals.append(agm)
            alt_vals.append(alt)

    if len(agm_vals) < 2:
        raise ValueError("AGM ALT sheet contains fewer than 2 valid numeric pairs")

    r_lin = _pearson_r(agm_vals, alt_vals)
    rho = _spearman_rho(agm_vals, alt_vals)
    r_log2 = _pearson_r([math.log2(v) for v in agm_vals], [math.log2(v) for v in alt_vals])
    passed = len(agm_vals) >= 90 and r_lin > 0.85 and rho > 0.85 and r_log2 > 0.85
    return {
        "file": AGM_ALT_XLSX_NAME,
        "sha256": hashlib.sha256(agm_xlsx.read_bytes()).hexdigest(),
        "pairs_evaluated": len(agm_vals),
        "pearson_r_linear": round(r_lin, 6),
        "pearson_r_log2": round(r_log2, 6),
        "spearman_rho": round(rho, 6),
        "passed": passed,
    }


def _solve_linear_system(a: List[List[float]], b: List[float]) -> List[float]:
    """Pure-stdlib Gaussian elimination with partial pivoting for square positive-definite system."""
    n = len(b)
    mat = [row[:] + [b[i]] for i, row in enumerate(a)]
    for col in range(n):
        pivot = max(range(col, n), key=lambda r: abs(mat[r][col]))
        if abs(mat[pivot][col]) < 1e-14:
            raise ValueError("Singular linear system in Ridge solver.")
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


def _raw_extra_features(rec: VariantRecord) -> List[float]:
    sub_feats = [1.0 if rec.sub_domain == sd else 0.0 for sd in EXTRA_SUBDOMAINS]
    cls_feats = [1.0 if rec.variant_class == vc else 0.0 for vc in EXTRA_CLASSES]
    return sub_feats + cls_feats


@dataclass(frozen=True)
class FittedPartialRidgeModel:
    """Linear model y = b0 + b1 * x_dmso + sum(w_j * ((raw_j - mean_j) / std_j))."""

    alpha: float
    b0_intercept: float
    b1_dmso_slope: float
    extra_means: Tuple[float, ...]
    extra_stds: Tuple[float, ...]
    extra_weights: Tuple[float, ...]

    def predict_one(self, rec: VariantRecord) -> float:
        raw = _raw_extra_features(rec)
        total = self.b0_intercept + self.b1_dmso_slope * rec.x_dmso_log2
        for w, val, mean, std in zip(self.extra_weights, raw, self.extra_means, self.extra_stds):
            if std > 0.0:
                total += w * ((val - mean) / std)
        return total


def fit_partial_penalized_ridge(
    train_rows: Sequence[VariantRecord],
    alpha: float,
) -> FittedPartialRidgeModel:
    """Fit partial-penalized Ridge leaving intercept & x_dmso unpenalized (0 diagonal penalty)."""
    if len(train_rows) < 5:
        raise ValueError("Need at least 5 rows to fit partial-penalized Ridge.")
    if alpha < 0.0:
        raise ValueError(f"alpha must be >= 0, got {alpha}")

    raw_mat = [_raw_extra_features(r) for r in train_rows]
    p_extra = len(raw_mat[0])
    n = len(train_rows)

    means = [sum(row[j] for row in raw_mat) / n for j in range(p_extra)]
    stds: List[float] = []
    for j in range(p_extra):
        var = sum((row[j] - means[j]) ** 2 for row in raw_mat) / n
        stds.append(math.sqrt(var) if var > 1e-12 else 0.0)

    # Full design matrix rows: [1.0, x_dmso_log2, z_1, ..., z_p_extra]
    p = 2 + p_extra
    x_mat: List[List[float]] = []
    y_vec: List[float] = []
    for i, rec in enumerate(train_rows):
        z = [
            (raw_mat[i][j] - means[j]) / stds[j] if stds[j] > 0.0 else 0.0
            for j in range(p_extra)
        ]
        x_mat.append([1.0, rec.x_dmso_log2] + z)
        y_vec.append(rec.y_combo_log2)

    xtx = [[0.0] * p for _ in range(p)]
    xty = [0.0] * p
    for row, y in zip(x_mat, y_vec):
        for a_i in range(p):
            xty[a_i] += row[a_i] * y
            for b_i in range(a_i, p):
                xtx[a_i][b_i] += row[a_i] * row[b_i]
    for a_i in range(p):
        for b_i in range(a_i):
            xtx[a_i][b_i] = xtx[b_i][a_i]

    # Unpenalized intercept (index 0) and x_dmso (index 1); zero-std features receive unit diagonal ridge (1.0) so XtX stays nonsingular.
    for j in range(2, p):
        pen = alpha if stds[j - 2] > 0.0 else 1.0
        xtx[j][j] += pen

    beta = _solve_linear_system(xtx, xty)
    return FittedPartialRidgeModel(
        alpha=alpha,
        b0_intercept=beta[0],
        b1_dmso_slope=beta[1],
        extra_means=tuple(means),
        extra_stds=tuple(stds),
        extra_weights=tuple(beta[2:]),
    )


def fit_dmso_ols(train_rows: Sequence[VariantRecord]) -> Tuple[float, float]:
    """Fit unpenalized 1D OLS: y_combo_log2 = b0 + b1 * x_dmso_log2."""
    xs = [r.x_dmso_log2 for r in train_rows]
    ys = [r.y_combo_log2 for r in train_rows]
    mx = sum(xs) / len(xs)
    my = sum(ys) / len(ys)
    denom = sum((x - mx) ** 2 for x in xs)
    if denom <= 1e-14:
        return my, 0.0
    b1 = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / denom
    b0 = my - b1 * mx
    return b0, b1


def _mae_log2_and_linear(
    rows: Sequence[VariantRecord], preds_log2: Sequence[float]
) -> Tuple[float, float]:
    mae_log2 = sum(abs(r.y_combo_log2 - p) for r, p in zip(rows, preds_log2)) / len(rows)
    mae_lin = sum(abs(r.combo_rel_mut_wt - (2.0 ** p)) for r, p in zip(rows, preds_log2)) / len(rows)
    return mae_log2, mae_lin


def residue_clustered_bootstrap_delta_ci(
    holdout_rows: Sequence[VariantRecord],
    stronger_preds: Sequence[float],
    cand_preds: Sequence[float],
    *,
    resamples: int = BOOTSTRAP_RESAMPLES,
    seed: int = BOOTSTRAP_SEED,
) -> Dict[str, Any]:
    """Cluster-bootstrap 95% CI on delta_MAE = MAE_stronger - MAE_cand over holdout codons."""
    by_codon: Dict[int, List[Tuple[float, float]]] = {}
    for rec, sp, cp in zip(holdout_rows, stronger_preds, cand_preds):
        abs_s = abs(rec.y_combo_log2 - sp)
        abs_c = abs(rec.y_combo_log2 - cp)
        by_codon.setdefault(rec.codon_position, []).append((abs_s, abs_c))

    clusters = list(by_codon.values())
    n_clusters = len(clusters)
    rng = random.Random(seed)
    deltas: List[float] = []
    for _ in range(resamples):
        sum_s = 0.0
        sum_c = 0.0
        cnt = 0
        for _c in range(n_clusters):
            cl = clusters[rng.randrange(n_clusters)]
            for abs_s, abs_c in cl:
                sum_s += abs_s
                sum_c += abs_c
                cnt += 1
        deltas.append((sum_s - sum_c) / cnt)

    deltas.sort()
    lo_idx = max(0, int(0.025 * resamples))
    hi_idx = min(resamples - 1, int(0.975 * resamples))
    return {
        "n_holdout_clusters": n_clusters,
        "resamples": resamples,
        "seed": seed,
        "delta_mae_ci95_log2": [round(deltas[lo_idx], 6), round(deltas[hi_idx], 6)],
        "ci95_low_gt_zero": deltas[lo_idx] > 0.0,
    }


def evaluate_tedman_benchmark(
    data_dir: Path,
    repo_root: Path = REPO_ROOT,
) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    """Run full positive controls, residue-clustered split, dev CV, and one-shot holdout gate."""
    master_xlsx = data_dir / MASTER_XLSX_NAME
    agm_xlsx = data_dir / AGM_ALT_XLSX_NAME

    master_sha = hashlib.sha256(master_xlsx.read_bytes()).hexdigest()
    agm_sha = hashlib.sha256(agm_xlsx.read_bytes()).hexdigest()
    if master_sha != MASTER_XLSX_SHA256:
        raise ValueError(f"Master xlsx SHA256 mismatch: {master_sha} != {MASTER_XLSX_SHA256}")
    if agm_sha != AGM_ALT_SHA256:
        raise ValueError(f"AGM ALT xlsx SHA256 mismatch: {agm_sha} != {AGM_ALT_SHA256}")

    agm_ctrl = check_agm_alt_positive_control(agm_xlsx)
    if not agm_ctrl["passed"]:
        raise RuntimeError(f"AGM ALT positive control failed: {agm_ctrl}")

    exposed_codons = load_exposed_pred0010_codons(repo_root)
    records, excluded = load_tedman_variants(master_xlsx, exposed_codons)

    dev_rows = [r for r in records if r.split == "dev"]
    holdout_rows = [r for r in records if r.split == "holdout"]
    dev_codons = {r.codon_position for r in dev_rows}
    holdout_codons = {r.codon_position for r in holdout_rows}

    if dev_codons & holdout_codons:
        raise RuntimeError("Leakage defect: residue overlap between dev and holdout.")
    if holdout_codons & set(exposed_codons):
        raise RuntimeError("Leakage defect: PRED-0010 exposed codon in holdout.")
    if len(holdout_codons) < MIN_HOLDOUT_CLUSTERS:
        raise RuntimeError(f"Insufficient holdout clusters: {len(holdout_codons)} < {MIN_HOLDOUT_CLUSTERS}")

    # Positive control 2: biological directionality on dev set only.
    dev_dmso_med = statistics.median([r.dmso_rel_mut_wt for r in dev_rows])
    dev_combo_med = statistics.median([r.combo_rel_mut_wt for r in dev_rows])
    dev_direction_ok = dev_dmso_med < 1.0 and dev_combo_med > dev_dmso_med
    if not dev_direction_ok:
        raise RuntimeError("Dev directionality check failed.")

    # Baseline 1: constant dev median in log2 space.
    dev_median_log2 = statistics.median([r.y_combo_log2 for r in dev_rows])

    # Baseline 2: dev-fitted unpenalized OLS on x_dmso_log2.
    ols_b0, ols_b1 = fit_dmso_ols(dev_rows)

    # Candidate: 5-fold residue-clustered CV strictly inside dev to pick alpha.
    sorted_dev_codons = sorted(dev_codons)
    codon_to_fold = {codon: idx % CV_FOLDS for idx, codon in enumerate(sorted_dev_codons)}

    cv_b1_preds_err: List[float] = []
    cv_b2_preds_err: List[float] = []
    for fold in range(CV_FOLDS):
        tr = [r for r in dev_rows if codon_to_fold[r.codon_position] != fold]
        va = [r for r in dev_rows if codon_to_fold[r.codon_position] == fold]
        tr_med = statistics.median([r.y_combo_log2 for r in tr])
        fb0, fb1 = fit_dmso_ols(tr)
        for r in va:
            cv_b1_preds_err.append(abs(r.y_combo_log2 - tr_med))
            cv_b2_preds_err.append(abs(r.y_combo_log2 - (fb0 + fb1 * r.x_dmso_log2)))

    cv_path: List[Dict[str, Any]] = []
    best_alpha = ALPHA_GRID[0]
    best_cv_mae = float("inf")
    for alpha in ALPHA_GRID:
        errs: List[float] = []
        for fold in range(CV_FOLDS):
            tr = [r for r in dev_rows if codon_to_fold[r.codon_position] != fold]
            va = [r for r in dev_rows if codon_to_fold[r.codon_position] == fold]
            mod = fit_partial_penalized_ridge(tr, alpha)
            for r in va:
                errs.append(abs(r.y_combo_log2 - mod.predict_one(r)))
        cv_mae = sum(errs) / len(errs)
        cv_path.append({"alpha": alpha, "dev_cv_mae_log2": round(cv_mae, 6)})
        if cv_mae < best_cv_mae:
            best_cv_mae = cv_mae
            best_alpha = alpha

    final_candidate = fit_partial_penalized_ridge(dev_rows, best_alpha)

    # Single-shot holdout evaluation (no retuning).
    b1_holdout_preds = [dev_median_log2] * len(holdout_rows)
    b2_holdout_preds = [ols_b0 + ols_b1 * r.x_dmso_log2 for r in holdout_rows]
    cand_holdout_preds = [final_candidate.predict_one(r) for r in holdout_rows]

    b1_mae_log2, b1_mae_lin = _mae_log2_and_linear(holdout_rows, b1_holdout_preds)
    b2_mae_log2, b2_mae_lin = _mae_log2_and_linear(holdout_rows, b2_holdout_preds)
    cand_mae_log2, cand_mae_lin = _mae_log2_and_linear(holdout_rows, cand_holdout_preds)

    dev_cv_b1_mae = sum(cv_b1_preds_err) / len(cv_b1_preds_err)
    dev_cv_b2_mae = sum(cv_b2_preds_err) / len(cv_b2_preds_err)
    dev_cv_stronger = "dmso_ols" if dev_cv_b2_mae <= dev_cv_b1_mae else "constant_dev_median"

    if b2_mae_log2 <= b1_mae_log2:
        stronger_name = "dmso_ols"
        stronger_mae_log2 = b2_mae_log2
        stronger_preds = b2_holdout_preds
    else:
        stronger_name = "constant_dev_median"
        stronger_mae_log2 = b1_mae_log2
        stronger_preds = b1_holdout_preds

    rel_improvement = (stronger_mae_log2 - cand_mae_log2) / stronger_mae_log2
    boot = residue_clustered_bootstrap_delta_ci(holdout_rows, stronger_preds, cand_holdout_preds)

    if rel_improvement >= GATE_REL_MAE_REDUCTION and boot["ci95_low_gt_zero"]:
        gate_standing = "open"
        gate_verdict = "PASS_CANDIDATE_BEATS_STRONGER_BASELINE"
    elif rel_improvement >= GATE_REL_MAE_REDUCTION:
        gate_standing = "refuted"
        gate_verdict = "INCONCLUSIVE_CI_STRADDLES_ZERO"
    else:
        gate_standing = "refuted"
        gate_verdict = "REFUTED_BELOW_10PCT_MAE_REDUCTION_THRESHOLD"

    manifest = {
        "schema_version": "tedman2025-split-v1",
        "source_master_xlsx": MASTER_XLSX_NAME,
        "source_master_sha256": master_sha,
        "exposed_pred0010_codons": list(exposed_codons),
        "total_valid_variants": len(records),
        "dev_variants": len(dev_rows),
        "holdout_variants": len(holdout_rows),
        "dev_unique_codons": len(dev_codons),
        "holdout_unique_codons": len(holdout_codons),
        "excluded_count": len(excluded),
        "excluded_rows": excluded,
        "assignments": [
            {
                "row_num": r.row_num,
                "protein_distinction": r.protein_distinction,
                "nucleic_acid_mutation": r.nucleic_acid_mutation,
                "sub_domain": r.sub_domain,
                "codon_position": r.codon_position,
                "variant_class": r.variant_class,
                "split": r.split,
                "exposed_pred0010_forced": r.exposed_pred0010_forced,
            }
            for r in records
        ],
    }

    manifest_sha = hashlib.sha256(
        json.dumps(manifest, sort_keys=True).encode("utf-8")
    ).hexdigest()

    results = {
        "schema_version": "tedman2025-item1-results-v1",
        "benchmark_id": "council-item-1-cftr-dms-expression",
        "evaluated_on": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%d"),
        "done_sha256": DONE_SHA256_EXPECTED,
        "split_manifest_sha256": manifest_sha,
        "dataset": {
            "source_id": "cftr-dms-tedman-2025",
            "doi": "10.17632/ksfmxv8zts.1",
            "paper_doi": "10.7554/eLife.107180",
            "license": "CC BY 4.0",
            "assay": "HEK293T parental CFTR DMS plasma-membrane intensity (VX-661+VX-445 vs DMSO, Mut/WT)",
            "scale_boundary": "Plasma-membrane expression in HEK293T cells is not CFTR channel conductance or clinical lung benefit.",
            "total_valid_variants": len(records),
            "excluded_count": len(excluded),
            "dev_variants": len(dev_rows),
            "holdout_variants": len(holdout_rows),
            "dev_unique_codons": len(dev_codons),
            "holdout_unique_codons": len(holdout_codons),
            "dev_pred0010_forced_variants": sum(1 for r in dev_rows if r.exposed_pred0010_forced),
            "holdout_pred0010_overlap_codons": 0,
        },
        "positive_controls": {
            "agm_alt_cross_study": agm_ctrl,
            "dev_directionality": {
                "dev_dmso_median_rel": round(dev_dmso_med, 6),
                "dev_combo_median_rel": round(dev_combo_med, 6),
                "passed": dev_direction_ok,
            },
        },
        "models_fitted_on_dev": {
            "constant_dev_median_log2": round(dev_median_log2, 6),
            "dmso_ols": {"b0_intercept": round(ols_b0, 6), "b1_dmso_slope": round(ols_b1, 6)},
            "dev_cv_5fold_baselines": {
                "constant_dev_median_mae_log2": round(dev_cv_b1_mae, 6),
                "dmso_ols_mae_log2": round(dev_cv_b2_mae, 6),
                "stronger_baseline_in_dev_cv": dev_cv_stronger,
            },
            "candidate_partial_ridge": {
                "selected_alpha": best_alpha,
                "cv_path": cv_path,
                "b0_intercept": round(final_candidate.b0_intercept, 6),
                "b1_dmso_slope": round(final_candidate.b1_dmso_slope, 6),
                "extra_feature_names": list(EXTRA_SUBDOMAINS) + list(EXTRA_CLASSES),
                "extra_weights": [round(w, 6) for w in final_candidate.extra_weights],
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
            "candidate_partial_ridge": {
                "mae_log2": round(cand_mae_log2, 6),
                "mae_linear_context_only": round(cand_mae_lin, 6),
            },
            "stronger_baseline_on_holdout": stronger_name,
            "stronger_baseline_mae_log2": round(stronger_mae_log2, 6),
            "relative_mae_reduction_vs_stronger": round(rel_improvement, 6),
            "residue_clustered_bootstrap": boot,
            "pre_registered_gate_threshold": GATE_REL_MAE_REDUCTION,
            "gate_verdict": gate_verdict,
            "ledger_standing": gate_standing,
        },
    }
    return manifest, results


def verify_checked_in_artifacts(repo_root: Path = REPO_ROOT, data_dir: Optional[Path] = None) -> Dict[str, Any]:
    """Verify structural invariants of checked-in artifacts (plus raw-data parity when local files exist)."""
    manifest_path = repo_root / "benchmarks" / "cftr_dms_expression" / MANIFEST_FILE_NAME
    results_path = repo_root / "benchmarks" / "cftr_dms_expression" / RESULTS_FILE_NAME
    done_path = repo_root / "benchmarks" / "cftr_dms_expression" / "done.md"

    done_sha = hashlib.sha256(
        done_path.read_text(encoding="utf-8", errors="replace").replace("\r\n", "\n").encode("utf-8")
    ).hexdigest()
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    results = json.loads(results_path.read_text(encoding="utf-8"))

    exposed = set(load_exposed_pred0010_codons(repo_root))
    assigns = manifest["assignments"]
    dev_codons = {int(a["codon_position"]) for a in assigns if a["split"] == "dev"}
    holdout_codons = {int(a["codon_position"]) for a in assigns if a["split"] == "holdout"}

    if dev_codons & holdout_codons:
        raise AssertionError("Split leakage: dev and holdout share codons.")
    if holdout_codons & exposed:
        raise AssertionError("Leakage: PRED-0010 exposed codons in holdout.")
    if results["done_sha256"] != done_sha:
        raise AssertionError("done_sha256 mismatch.")

    manifest_sha = hashlib.sha256(json.dumps(manifest, sort_keys=True).encode("utf-8")).hexdigest()
    if results["split_manifest_sha256"] != manifest_sha:
        raise AssertionError("split_manifest_sha256 mismatch.")

    resolved_data = resolve_tedman_data_dir(data_dir, repo_root=repo_root)
    raw_parity_checked = False
    if resolved_data is not None:
        fresh_manifest, fresh_results = evaluate_tedman_benchmark(resolved_data, repo_root=repo_root)
        if fresh_manifest["assignments"] != manifest["assignments"]:
            raise AssertionError("Live raw .xlsx assignments do not match split_manifest.json.")
        if abs(
            fresh_results["holdout_evaluation"]["candidate_partial_ridge"]["mae_log2"]
            - results["holdout_evaluation"]["candidate_partial_ridge"]["mae_log2"]
        ) > 1e-6:
            raise AssertionError("Live raw .xlsx candidate MAE differs from results JSON.")
        raw_parity_checked = True

    return {
        "status": "PASS",
        "raw_local_data_parity_checked": raw_parity_checked,
        "dev_variants": manifest["dev_variants"],
        "holdout_variants": manifest["holdout_variants"],
        "holdout_unique_codons": manifest["holdout_unique_codons"],
        "gate_verdict": results["holdout_evaluation"]["gate_verdict"],
    }


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--data-dir", type=Path, default=None, help="Optional path to local tedman-2025 folder.")
    parser.add_argument("--write-artifacts", action="store_true", help="Write split_manifest.json and results JSON.")
    parser.add_argument("--check", action="store_true", help="Verify checked-in split manifest, results, and local raw parity if present.")
    parser.add_argument("--json", action="store_true", help="Print strict JSON (allow_nan=False).")
    args = parser.parse_args(argv)

    if args.check:
        summary = verify_checked_in_artifacts(repo_root=REPO_ROOT, data_dir=args.data_dir)
        if args.json:
            print(json.dumps(summary, indent=2, allow_nan=False))
        else:
            print(
                f"PASS cftr_dms_benchmark: dev={summary['dev_variants']}, holdout={summary['holdout_variants']} "
                f"(codons={summary['holdout_unique_codons']}), raw_parity={summary['raw_local_data_parity_checked']}, "
                f"verdict={summary['gate_verdict']}"
            )
        return 0

    resolved_data = resolve_tedman_data_dir(args.data_dir, repo_root=REPO_ROOT)
    if resolved_data is None:
        print("Local Tedman 2025 files not present under sources/downloads/tedman-2025/.", file=sys.stderr)
        return 1

    manifest, results = evaluate_tedman_benchmark(resolved_data, repo_root=REPO_ROOT)
    if args.write_artifacts:
        (HERE / MANIFEST_FILE_NAME).write_text(
            json.dumps(manifest, indent=2, allow_nan=False) + "\n", encoding="utf-8"
        )
        (HERE / RESULTS_FILE_NAME).write_text(
            json.dumps(results, indent=2, allow_nan=False) + "\n", encoding="utf-8"
        )

    if args.json:
        print(json.dumps(results, indent=2, allow_nan=False))
    else:
        he = results["holdout_evaluation"]
        print(
            f"Tedman 2025 Item #1: verdict={he['gate_verdict']} | "
            f"stronger_baseline({he['stronger_baseline_on_holdout']}) MAE_log2={he['stronger_baseline_mae_log2']} | "
            f"candidate MAE_log2={he['candidate_partial_ridge']['mae_log2']} | "
            f"rel_reduction={he['relative_mae_reduction_vs_stronger']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
