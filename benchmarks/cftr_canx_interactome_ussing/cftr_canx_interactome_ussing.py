#!/usr/bin/env python3
"""CFTR calnexin interactome & FRT siCANX Ussing asymmetry benchmark (Tedman et al. 2025).

Evaluates two independently gated pre-registered arms on the untapped CC BY 4.0
Tedman et al. 2025 (DOI 10.7554/eLife.107180; Mendeley Data 10.17632/ksfmxv8zts.1)
workbooks under git-ignored `sources/downloads/tedman-2025/`:
- Arm 1 (`frt_ussing_sicanx_asymmetry`): FRT-F508del vs. FRT-WT Ussing Delta Isc under
  siCANX vs. siNS with strict replicate deduplication (Bonferroni alpha = 0.025).
- Arm 2 (`canx_ko_proteostasis_remodeling`): HEK293T TMT-18plex CFTR interactome
  remodeling across 4 canonical proteostasis modules (Drug 4-test BH-FDR q < 0.05
  AND Variants 20-test BH-FDR q < 0.05).

Pure Python 3.10+ standard library only; zero network calls; zero `.xlsx` files
committed to git. Research benchmark only: macroscopic FRT Delta Isc and HEK293T
co-IP fold changes are not single-channel burst kinetics or clinical patient data.
"""
from __future__ import annotations

import argparse
import bisect
import collections
from dataclasses import dataclass
import datetime as dt
import hashlib
import itertools
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

DONE_SHA256_EXPECTED = "d6eeee7cf61baf2cee4576f6094c278a623b214d3667effe42ff05c14360ca82"
RESULTS_FILE_NAME = "results_2026-11-25.json"

WORKBOOK_SHA256: Dict[str, str] = {
    "Copy of TMT_Labels.xlsx": "5b16d08dfefce36878149259d823604d2153df8b027a48bd8cf69c3f1af5b41e",
    "Copy of CFTR_Variants_LogFoldChanges.xlsx": "492df3ef266f27a5e18e1d9ce61b13356a8101b3e85d89c38ec87f5ca7208f7a",
    "Copy of CFTR_F508del_Drug_LogFoldChanges.xlsx": "d42ae782fe85163909396ea5ab4ee0915eec346e8ae420e419e80b5dd5e270c0",
    "Copy of FRT-F508del_siCANX_Ussing-summary-for-paper.xlsx": "a1ad7e325bab60ba6fe15e001fdddfd4fd14918d571e53ca45844c2f58f1ac1a",
    "Copy of FRT-WT_siCANX_Ussing-summary-for-paper.xlsx": "bddaffb3b68449184b518a923c3c0d18ce51ed409662915d7b16fe7f8fedc628",
}
MASTER_DMS_XLSX = "ALT 9-6-24 CFTR DMS Expanded Library Master Tables.xlsx"
MASTER_DMS_SHA256 = "2a702ee3744b8039a9c7c4c74aec9bb30a91301b138d4948e0e2c81a1931e51b"

SEED = 20261125
BOOTSTRAP_RESAMPLES = 10000
PERM_NULL_DRAWS = 10000
ALPHA_USSING = 0.025
ALPHA_PROTEOMICS_FDR = 0.05
MAX_NON_NUMERIC_PER_WORKBOOK = 5
MIN_USSING_MASS_BALANCE_R = 0.95
MIN_WT_F508DEL_RANGE_MULT = 5.0
MIN_DISPERSION_STD = 0.05
MAX_NULL_FWER = 0.08

CANONICAL_MODULES: Dict[str, frozenset[str]] = {
    "folding": frozenset({"folding", "protein folding"}),
    "ups": frozenset({"ubiquitin proteasome system"}),
    "autophagy_lysosome": frozenset({"autophagy-lysosome pathway"}),
    "trafficking": frozenset({"trafficking"}),
}
CANONICAL_MODULE_ORDER: Tuple[str, ...] = ("folding", "ups", "autophagy_lysosome", "trafficking")

BROAD_PROTEOSTASIS_SET: frozenset[str] = frozenset(
    {
        "folding",
        "protein folding",
        "ubiquitin proteasome system",
        "autophagy-lysosome pathway",
        "trafficking",
        "er proteostasis",
        "cytonuclear proteostasis",
        "mitochondrial proteostasis",
        "nuclear proteostasis",
        "extracellular proteostasis",
        "er",
        "er, trafficking",
        "vesicle trafficking",
        "quality control of membrane proteins",
    }
)

EXPECTED_USSING_GROUPS: Tuple[str, ...] = (
    "NT",
    "siNS, DMSO",
    "siNS, VX-661",
    "siNS, VX-445",
    "siNS, 661+445",
    "siNS, ETI",
    "siCANX, DMSO",
    "siCANX, VX-661",
    "siCANX, VX-445",
    "siCANX, 661+445",
    "siCANX, ETI",
)
EXPECTED_DRUG_ORDER: Tuple[str, str, str] = ("Fsk", "VX-770", "Inh172")

VARIANT_COLS: Tuple[str, ...] = ("WT(-CX)", "V232D (-CX)", "P5L (-CX)", "T1036N (-CX)", "I1366N (-CX)")
DRUG_COLS: Tuple[str, ...] = ("ΔF (-CX)", "DelF (-CX) +445", "DelF (-CX) +661", "DelF (-CX) + combo")

XLSX_NS = {"s": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
CELL_REF_RE = re.compile(r"^([A-Z]+)(\d+)$")
RANGE_REF_RE = re.compile(r"^([A-Z]+)(\d+):([A-Z]+)(\d+)$")


@dataclass(frozen=True)
class UssingReplicate:
    """Single Ussing replicate (Fsk, VX-770, Inh172) from Raw-data-analyzed."""

    col: str
    fsk: float
    vx770: float
    inh172: float

    @property
    def total_act(self) -> float:
        return self.fsk + self.vx770

    @property
    def balance_sum(self) -> float:
        return self.fsk + self.vx770 + self.inh172


@dataclass(frozen=True)
class UssingGroupData:
    """As-deposited and strictly deduplicated replicates for one Ussing Sample ID."""

    workbook: str
    group: str
    row_start: int
    as_deposited: Tuple[UssingReplicate, ...]
    deduplicated: Tuple[UssingReplicate, ...]
    duplicates: Tuple[Dict[str, Any], ...]


@dataclass(frozen=True)
class ProteomicsRow:
    """One non-empty protein row from a LogFoldChanges sheet."""

    row_num: int
    gene_symbol: str
    pathway_raw: str
    pathway_norm: str
    values: Dict[str, Optional[float]]


def resolve_tedman_data_dir(
    data_dir: Optional[Path] = None, repo_root: Path = REPO_ROOT
) -> Optional[Path]:
    """Locate git-ignored `tedman-2025/` data folder across main checkout or worktree."""
    anchor = "Copy of FRT-F508del_siCANX_Ussing-summary-for-paper.xlsx"
    if data_dir is not None:
        p = Path(data_dir)
        return p if (p / anchor).is_file() else None
    direct = repo_root / "sources" / "downloads" / "tedman-2025"
    if (direct / anchor).is_file():
        return direct
    if len(repo_root.parents) >= 3:
        main_sibling = repo_root.parents[2] / "cf-research" / "sources" / "downloads" / "tedman-2025"
        if (main_sibling / anchor).is_file():
            return main_sibling
    return None


def verify_required_workbooks(data_dir: Path) -> Dict[str, str]:
    """Verify exact SHA-256 digests for all five required Tedman et al. 2025 workbooks."""
    digests: Dict[str, str] = {}
    for fname, expected in WORKBOOK_SHA256.items():
        fpath = data_dir / fname
        if not fpath.is_file():
            raise FileNotFoundError(f"Missing required workbook: {fpath}")
        actual = hashlib.sha256(fpath.read_bytes()).hexdigest()
        if actual != expected:
            raise ValueError(f"SHA-256 mismatch for {fname}: {actual} != {expected}")
        digests[fname] = actual
    return digests


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
        joined = "".join(
            t.text or ""
            for t in is_elem.iter("{http://schemas.openxmlformats.org/spreadsheetml/2006/main}t")
        )
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
                rid = sheet_elem.attrib[
                    "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"
                ]
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


def _parse_finite_float(raw: Optional[str], context: str) -> float:
    if raw is None or str(raw).strip() == "":
        raise ValueError(f"Missing numeric value at {context}")
    val = float(raw)
    if not math.isfinite(val):
        raise ValueError(f"Non-finite numeric value {raw!r} at {context}")
    return val


def parse_ussing_workbook(
    xlsx_path: Path,
    workbook_label: str,
    expected_avg_col: str,
) -> Tuple[Dict[str, UssingGroupData], List[Dict[str, Any]]]:
    """Parse raw replicate tuples from `Raw-data-analyzed`, detecting exact duplicates."""
    cells, _ = parse_xlsx_sheet(xlsx_path, "Raw-data-analyzed", expand_merged=False)
    if cells.get((1, "A"), "").strip() != "Sample ID" or cells.get((1, "B"), "").strip() != "Acute drugs":
        raise ValueError(f"Unexpected header A1/B1 in {xlsx_path.name}")
    if cells.get((1, expected_avg_col), "").strip() != "AVG":
        raise ValueError(
            f"Expected AVG header at {expected_avg_col}1 in {xlsx_path.name}, got {cells.get((1, expected_avg_col))!r}"
        )

    avg_idx = _col_to_index(expected_avg_col)
    rep_cols = [_index_to_col(c_idx) for c_idx in range(_col_to_index("C"), avg_idx)]
    for pos, col in enumerate(rep_cols, start=1):
        lbl = cells.get((2, col), "").strip()
        if lbl != str(pos):
            raise ValueError(f"Expected replicate index {pos} at {col}2 in {xlsx_path.name}, got {lbl!r}")

    groups: Dict[str, UssingGroupData] = {}
    dup_records: List[Dict[str, Any]] = []

    for idx, expected_grp in enumerate(EXPECTED_USSING_GROUPS):
        r0 = 3 + 3 * idx
        actual_grp = cells.get((r0, "A"), "").strip()
        if actual_grp != expected_grp:
            raise ValueError(
                f"Expected group {expected_grp!r} at A{r0} in {xlsx_path.name}, got {actual_grp!r}"
            )
        drugs = tuple(cells.get((r0 + k, "B"), "").strip() for k in range(3))
        if drugs != EXPECTED_DRUG_ORDER:
            raise ValueError(
                f"Unexpected drug order {drugs!r} at B{r0}:B{r0+2} in {xlsx_path.name}"
            )

        as_dep: List[UssingReplicate] = []
        dedup: List[UssingReplicate] = []
        seen_tuples: Dict[Tuple[float, float, float], str] = {}
        grp_dups: List[Dict[str, Any]] = []

        for col in rep_cols:
            raw_vals = [cells.get((r0 + k, col)) for k in range(3)]
            non_empty = [v is not None and str(v).strip() != "" for v in raw_vals]
            if not any(non_empty):
                continue
            if not all(non_empty):
                raise ValueError(
                    f"Partial replicate column {col} for {expected_grp!r} in {xlsx_path.name}: {raw_vals!r}"
                )
            fsk = _parse_finite_float(raw_vals[0], f"{xlsx_path.name}:{col}{r0}")
            vx770 = _parse_finite_float(raw_vals[1], f"{xlsx_path.name}:{col}{r0+1}")
            inh172 = _parse_finite_float(raw_vals[2], f"{xlsx_path.name}:{col}{r0+2}")
            rep = UssingReplicate(col=col, fsk=fsk, vx770=vx770, inh172=inh172)
            as_dep.append(rep)

            key = (fsk, vx770, inh172)
            if key in seen_tuples:
                rec = {
                    "workbook": workbook_label,
                    "group": expected_grp,
                    "kept_col": seen_tuples[key],
                    "dropped_col": col,
                    "tuple_fsk_vx770_inh172": [fsk, vx770, inh172],
                }
                grp_dups.append(rec)
                dup_records.append(rec)
            else:
                seen_tuples[key] = col
                dedup.append(rep)

        if not dedup:
            raise ValueError(f"Zero replicates parsed for {expected_grp!r} in {xlsx_path.name}")

        groups[expected_grp] = UssingGroupData(
            workbook=workbook_label,
            group=expected_grp,
            row_start=r0,
            as_deposited=tuple(as_dep),
            deduplicated=tuple(dedup),
            duplicates=tuple(grp_dups),
        )

    return groups, dup_records


def _pearson_r(xs: Sequence[float], ys: Sequence[float]) -> float:
    n = len(xs)
    if n < 2 or len(ys) != n:
        raise ValueError("Pearson r requires >= 2 paired values")
    mx = sum(xs) / n
    my = sum(ys) / n
    num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    dx = sum((x - mx) ** 2 for x in xs)
    dy = sum((y - my) ** 2 for y in ys)
    if dx <= 0.0 or dy <= 0.0:
        raise ValueError("Pearson r undefined for zero-variance input")
    return num / math.sqrt(dx * dy)


def _rank_average_ties(values: Sequence[float]) -> Tuple[List[float], List[int]]:
    """Return 1-based average midranks and tie group sizes."""
    n = len(values)
    indexed = sorted(enumerate(values), key=lambda kv: kv[1])
    ranks = [0.0] * n
    tie_sizes: List[int] = []
    i = 0
    while i < n:
        j = i + 1
        while j < n and indexed[j][1] == indexed[i][1]:
            j += 1
        k = j - i
        if k > 1:
            tie_sizes.append(k)
        avg_rank = 0.5 * ((i + 1) + j)
        for pos in range(i, j):
            ranks[indexed[pos][0]] = avg_rank
        i = j
    return ranks, tie_sizes


def spearman_rank_correlation(xs: Sequence[float], ys: Sequence[float]) -> float:
    """Compute Spearman rank correlation with average midranks for ties."""
    if len(xs) != len(ys) or len(xs) < 2:
        raise ValueError("Spearman correlation requires >= 2 paired items")
    rx, _ = _rank_average_ties(xs)
    ry, _ = _rank_average_ties(ys)
    return _pearson_r(rx, ry)


def _betacf(a: float, b: float, x: float) -> float:
    """Evaluate continued fraction for regularized incomplete beta I_x(a, b)."""
    max_iter = 200
    eps = 3.0e-14
    fpmin = 1.0e-300
    qab = a + b
    qap = a + 1.0
    qam = a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < fpmin:
        d = fpmin
    d = 1.0 / d
    h = d
    for m in range(1, max_iter + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < fpmin:
            d = fpmin
        c = 1.0 + aa / c
        if abs(c) < fpmin:
            c = fpmin
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < fpmin:
            d = fpmin
        c = 1.0 + aa / c
        if abs(c) < fpmin:
            c = fpmin
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < eps:
            break
    return h


def regularized_incomplete_beta(a: float, b: float, x: float) -> float:
    """Compute regularized incomplete beta function I_x(a, b) for x in [0, 1]."""
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    bt = math.exp(
        math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b) + a * math.log(x) + b * math.log(1.0 - x)
    )
    if x < (a + 1.0) / (a + b + 2.0):
        return bt * _betacf(a, b, x) / a
    return 1.0 - bt * _betacf(b, a, 1.0 - x) / b


def welch_t_test(xs: Sequence[float], ys: Sequence[float]) -> Dict[str, float]:
    """Two-sample Welch's t-test with Satterthwaite df and exact Student-t p-value."""
    n1, n2 = len(xs), len(ys)
    if n1 < 2 or n2 < 2:
        raise ValueError("Welch t-test requires at least 2 observations per group")
    m1 = sum(xs) / n1
    m2 = sum(ys) / n2
    v1 = sum((x - m1) ** 2 for x in xs) / (n1 - 1)
    v2 = sum((y - m2) ** 2 for y in ys) / (n2 - 1)
    se2 = v1 / n1 + v2 / n2
    if se2 <= 0.0:
        raise ValueError("Zero pooled variance in Welch t-test")
    t_stat = (m1 - m2) / math.sqrt(se2)
    num_df = se2 ** 2
    den_df = ((v1 / n1) ** 2) / (n1 - 1) + ((v2 / n2) ** 2) / (n2 - 1)
    df = num_df / den_df
    x_beta = df / (df + t_stat * t_stat)
    p_two_sided = regularized_incomplete_beta(0.5 * df, 0.5, x_beta)
    return {
        "t_stat": t_stat,
        "df": df,
        "p_two_sided": min(1.0, max(0.0, p_two_sided)),
    }


def exact_two_sample_permutation_p(xs: Sequence[float], ys: Sequence[float]) -> Tuple[float, int, int]:
    """Exact two-sided permutation enumeration over C(n1+n2, n1) mean differences."""
    n1, n2 = len(xs), len(ys)
    pooled = list(xs) + list(ys)
    total = sum(pooled)
    obs_abs = abs((sum(xs) / n1) - (sum(ys) / n2))
    n_tot = n1 + n2
    ge_count = 0
    total_splits = 0
    for comb in itertools.combinations(range(n_tot), n1):
        s1 = sum(pooled[i] for i in comb)
        s2 = total - s1
        diff_abs = abs((s1 / n1) - (s2 / n2))
        if diff_abs >= obs_abs - 1e-12:
            ge_count += 1
        total_splits += 1
    return ge_count / total_splits, ge_count, total_splits


def bootstrap_log_ratio_ci(
    xs_num: Sequence[float],
    xs_den: Sequence[float],
    *,
    resamples: int = BOOTSTRAP_RESAMPLES,
    seed: int = SEED,
) -> Dict[str, Any]:
    """Stratified bootstrap on ln(mean(xs_num) / mean(xs_den))."""
    n1, n2 = len(xs_num), len(xs_den)
    obs = math.log((sum(xs_num) / n1) / (sum(xs_den) / n2))
    rng = random.Random(seed)
    draws: List[float] = []
    for _ in range(resamples):
        m1 = sum(xs_num[rng.randrange(n1)] for _ in range(n1)) / n1
        m2 = sum(xs_den[rng.randrange(n2)] for _ in range(n2)) / n2
        if m1 <= 0.0 or m2 <= 0.0:
            raise ValueError("Non-positive bootstrap group mean for log-ratio")
        draws.append(math.log(m1 / m2))
    draws.sort()
    return {
        "point_estimate": round(obs, 6),
        "ci95": [round(draws[int(0.025 * resamples)], 6), round(draws[int(0.975 * resamples)], 6)],
        "ci975": [round(draws[int(0.0125 * resamples)], 6), round(draws[int(0.9875 * resamples)], 6)],
        "ci975_high_lt_zero": draws[int(0.9875 * resamples)] < 0.0,
    }


def exact_stratified_did_permutation_p(
    f_canx: Sequence[float],
    f_ns: Sequence[float],
    w_canx: Sequence[float],
    w_ns: Sequence[float],
) -> Tuple[float, float, int]:
    """Exact stratified DiD permutation over all within-cell-line label swaps."""
    f_pool = list(f_canx) + list(f_ns)
    w_pool = list(w_canx) + list(w_ns)
    nf1, nf2 = len(f_canx), len(f_ns)
    nw1, nw2 = len(w_canx), len(w_ns)
    f_tot, w_tot = sum(f_pool), sum(w_pool)

    obs_lf = math.log((sum(f_canx) / nf1) / (sum(f_ns) / nf2))
    obs_lw = math.log((sum(w_canx) / nw1) / (sum(w_ns) / nw2))
    did_obs = obs_lf - obs_lw
    abs_obs = abs(did_obs)

    f_lrs: List[float] = []
    for comb in itertools.combinations(range(nf1 + nf2), nf1):
        s1 = sum(f_pool[i] for i in comb)
        s2 = f_tot - s1
        f_lrs.append(math.log((s1 / nf1) / (s2 / nf2)))

    w_lrs: List[float] = []
    for comb in itertools.combinations(range(nw1 + nw2), nw1):
        s1 = sum(w_pool[i] for i in comb)
        s2 = w_tot - s1
        w_lrs.append(math.log((s1 / nw1) / (s2 / nw2)))
    w_lrs.sort()

    nw_splits = len(w_lrs)
    total_splits = len(f_lrs) * nw_splits
    one_sided_count = 0
    two_sided_count = 0
    for lf in f_lrs:
        # lf - lw <= did_obs  <=>  lw >= lf - did_obs
        idx_le = bisect.bisect_left(w_lrs, lf - did_obs - 1e-12)
        one_sided_count += nw_splits - idx_le
        # |lf - lw| >= abs_obs <=> lw <= lf - abs_obs or lw >= lf + abs_obs
        idx_low = bisect.bisect_right(w_lrs, lf - abs_obs + 1e-12)
        idx_high = bisect.bisect_left(w_lrs, lf + abs_obs - 1e-12)
        if idx_low >= idx_high:
            two_sided_count += nw_splits
        else:
            two_sided_count += idx_low + (nw_splits - idx_high)

    return one_sided_count / total_splits, two_sided_count / total_splits, total_splits


def bootstrap_did_log_ratio_ci(
    f_canx: Sequence[float],
    f_ns: Sequence[float],
    w_canx: Sequence[float],
    w_ns: Sequence[float],
    *,
    resamples: int = BOOTSTRAP_RESAMPLES,
    seed: int = SEED,
) -> Dict[str, Any]:
    """Stratified 4-group bootstrap on ln(F_c/F_s) - ln(W_c/W_s)."""
    nf1, nf2 = len(f_canx), len(f_ns)
    nw1, nw2 = len(w_canx), len(w_ns)
    obs = math.log((sum(f_canx) / nf1) / (sum(f_ns) / nf2)) - math.log(
        (sum(w_canx) / nw1) / (sum(w_ns) / nw2)
    )
    rng = random.Random(seed)
    draws: List[float] = []
    for _ in range(resamples):
        m_fc = sum(f_canx[rng.randrange(nf1)] for _ in range(nf1)) / nf1
        m_fs = sum(f_ns[rng.randrange(nf2)] for _ in range(nf2)) / nf2
        m_wc = sum(w_canx[rng.randrange(nw1)] for _ in range(nw1)) / nw1
        m_ws = sum(w_ns[rng.randrange(nw2)] for _ in range(nw2)) / nw2
        if m_fc <= 0.0 or m_fs <= 0.0 or m_wc <= 0.0 or m_ws <= 0.0:
            raise ValueError("Non-positive bootstrap mean in DiD log-ratio")
        draws.append(math.log(m_fc / m_fs) - math.log(m_wc / m_ws))
    draws.sort()
    q_one_sided_975 = draws[int(0.975 * resamples)]
    return {
        "point_estimate": round(obs, 6),
        "one_sided_975_upper": round(q_one_sided_975, 6),
        "one_sided_975_high_lt_zero": q_one_sided_975 < 0.0,
        "two_sided_ci975": [round(draws[int(0.0125 * resamples)], 6), round(draws[int(0.9875 * resamples)], 6)],
        "two_sided_ci95": [round(draws[int(0.025 * resamples)], 6), round(draws[int(0.975 * resamples)], 6)],
    }


def _evaluate_ussing_pair(
    f508: Dict[str, UssingGroupData],
    wt: Dict[str, UssingGroupData],
    *,
    use_dedup: bool,
) -> Dict[str, Any]:
    def _vals(grp: UssingGroupData) -> List[float]:
        reps = grp.deduplicated if use_dedup else grp.as_deposited
        return [r.total_act for r in reps]

    f_canx = _vals(f508["siCANX, 661+445"])
    f_ns = _vals(f508["siNS, 661+445"])
    w_canx = _vals(wt["siCANX, DMSO"])
    w_ns = _vals(wt["siNS, DMSO"])

    # Arm 1a: F508del 661+445 siCANX vs siNS
    welch_1a = welch_t_test(f_canx, f_ns)
    p_perm_1a, ge_1a, splits_1a = exact_two_sample_permutation_p(f_canx, f_ns)
    boot_1a = bootstrap_log_ratio_ci(f_canx, f_ns, resamples=BOOTSTRAP_RESAMPLES, seed=SEED)
    pass_1a = (p_perm_1a <= ALPHA_USSING) and bool(boot_1a["ci975_high_lt_zero"])

    # Arm 1b: F508del 661+445 vs WT DMSO proportional siCANX drop
    p_did_1s, p_did_2s, splits_1b = exact_stratified_did_permutation_p(f_canx, f_ns, w_canx, w_ns)
    boot_1b = bootstrap_did_log_ratio_ci(
        f_canx, f_ns, w_canx, w_ns, resamples=BOOTSTRAP_RESAMPLES, seed=SEED
    )
    pass_1b = (p_did_1s <= ALPHA_USSING) and bool(boot_1b["one_sided_975_high_lt_zero"])

    supported = pass_1a and pass_1b
    return {
        "arm_1a_f508del_combo_sicanx_drop": {
            "n_sicanx": len(f_canx),
            "n_sins": len(f_ns),
            "mean_sicanx_total_act": round(sum(f_canx) / len(f_canx), 6),
            "mean_sins_total_act": round(sum(f_ns) / len(f_ns), 6),
            "diff_mean": round((sum(f_canx) / len(f_canx)) - (sum(f_ns) / len(f_ns)), 6),
            "pct_change_vs_sins": round(
                100.0 * (((sum(f_canx) / len(f_canx)) / (sum(f_ns) / len(f_ns))) - 1.0), 6
            ),
            "exact_perm_p_two_sided": round(p_perm_1a, 6),
            "exact_perm_ge_count": ge_1a,
            "exact_perm_total_splits": splits_1a,
            "welch_t_context": {
                "t_stat": round(welch_1a["t_stat"], 6),
                "df": round(welch_1a["df"], 6),
                "p_two_sided": round(welch_1a["p_two_sided"], 6),
            },
            "bootstrap_log_ratio": boot_1a,
            "alpha_bonferroni": ALPHA_USSING,
            "passed_gate": pass_1a,
        },
        "arm_1b_sicanx_proportion_f508del_vs_wt": {
            "n_f508del_sicanx": len(f_canx),
            "n_f508del_sins": len(f_ns),
            "n_wt_sicanx_dmso": len(w_canx),
            "n_wt_sins_dmso": len(w_ns),
            "f508del_combo_log_ratio": round(math.log((sum(f_canx) / len(f_canx)) / (sum(f_ns) / len(f_ns))), 6),
            "wt_dmso_log_ratio": round(math.log((sum(w_canx) / len(w_canx)) / (sum(w_ns) / len(w_ns))), 6),
            "wt_dmso_pct_change": round(
                100.0 * (((sum(w_canx) / len(w_canx)) / (sum(w_ns) / len(w_ns))) - 1.0), 6
            ),
            "did_log_ratio": boot_1b["point_estimate"],
            "exact_stratified_perm_p_one_sided": round(p_did_1s, 6),
            "exact_stratified_perm_p_two_sided_context": round(p_did_2s, 6),
            "exact_stratified_perm_total_splits": splits_1b,
            "bootstrap_did_log_ratio": boot_1b,
            "alpha_bonferroni": ALPHA_USSING,
            "passed_gate": pass_1b,
        },
        "overall_verdict": "supported" if supported else "refuted",
        "ledger_standing": "open" if supported else "refuted",
    }


def _summarize_ussing_groups(groups: Dict[str, UssingGroupData], use_dedup: bool = True) -> Dict[str, Any]:
    summary: Dict[str, Any] = {}
    for gname in EXPECTED_USSING_GROUPS:
        g = groups[gname]
        reps = g.deduplicated if use_dedup else g.as_deposited
        n = len(reps)
        totals = [r.total_act for r in reps]
        inhs = [r.inh172 for r in reps]
        bals = [r.balance_sum for r in reps]
        m_tot = sum(totals) / n
        m_inh = sum(inhs) / n
        sem_tot = math.sqrt(sum((x - m_tot) ** 2 for x in totals) / (n * (n - 1))) if n > 1 else 0.0
        summary[gname] = {
            "n": n,
            "mean_fsk": round(sum(r.fsk for r in reps) / n, 6),
            "mean_vx770": round(sum(r.vx770 for r in reps) / n, 6),
            "mean_total_act": round(m_tot, 6),
            "sem_total_act": round(sem_tot, 6),
            "mean_inh172": round(m_inh, 6),
            "mean_mass_balance_residual": round(sum(bals) / n, 6),
        }
    return summary


def run_arm1_positive_controls(
    f508: Dict[str, UssingGroupData],
    wt: Dict[str, UssingGroupData],
    wt_xlsx: Path,
) -> Dict[str, Any]:
    """Evaluate all Arm 1 positive controls before interpreting Arm 1 hypotheses."""
    # 1. Opposite-sign Inh172 reversal & mass-balance correlation > 0.95
    wb_checks: Dict[str, Any] = {}
    pooled_pos: List[float] = []
    pooled_neg_inh: List[float] = []
    for label, gmap in (("FRT-F508del", f508), ("FRT-WT", wt)):
        xs_pos: List[float] = []
        ys_neg_inh: List[float] = []
        xs_non_eti: List[float] = []
        ys_non_eti: List[float] = []
        all_signs_ok = True
        max_abs_resid = 0.0
        for gname in EXPECTED_USSING_GROUPS:
            reps = gmap[gname].deduplicated
            m_act = sum(r.total_act for r in reps) / len(reps)
            m_inh = sum(r.inh172 for r in reps) / len(reps)
            if not (m_act > 0.0 and m_inh < 0.0):
                all_signs_ok = False
            xs_pos.append(m_act)
            ys_neg_inh.append(-m_inh)
            pooled_pos.append(m_act)
            pooled_neg_inh.append(-m_inh)
            if "ETI" not in gname:
                xs_non_eti.append(m_act)
                ys_non_eti.append(-m_inh)
            max_abs_resid = max(max_abs_resid, abs(m_act + m_inh))
        r_all = _pearson_r(xs_pos, ys_neg_inh)
        r_non_eti = _pearson_r(xs_non_eti, ys_non_eti)
        if not all_signs_ok:
            raise RuntimeError(f"Arm 1 Inh172 sign control failed for {label}")
        if label == "FRT-F508del" and r_all <= MIN_USSING_MASS_BALANCE_R:
            raise RuntimeError(f"Arm 1 Inh172 reversal control failed for {label}: r={r_all}")
        if label == "FRT-WT" and r_non_eti <= MIN_USSING_MASS_BALANCE_R:
            raise RuntimeError(f"Arm 1 Inh172 non-ETI reversal control failed for {label}: r={r_non_eti}")
        wb_checks[label] = {
            "all_11_groups_pos_act_neg_inh": all_signs_ok,
            "pearson_r_all_11_groups": round(r_all, 6),
            "pearson_r_9_non_eti_groups": round(r_non_eti, 6),
            "max_abs_group_mean_balance_residual": round(max_abs_resid, 6),
            "passed": True,
        }
    r_pooled_22 = _pearson_r(pooled_pos, pooled_neg_inh)
    if r_pooled_22 <= MIN_USSING_MASS_BALANCE_R:
        raise RuntimeError(f"Arm 1 pooled 22-group Inh172 reversal control failed: r={r_pooled_22}")
    wb_checks["pooled_22_groups_pearson_r"] = round(r_pooled_22, 6)

    # 2. WT NT dynamic range over F508del NT (> 5x) on both C..R (N=15) and C..L (N=10)
    wt_nt_all = [r.total_act for r in wt["NT"].deduplicated]
    wt_nt_c_to_l = [r.total_act for r in wt["NT"].deduplicated if _col_to_index(r.col) <= _col_to_index("L")]
    f508_nt = [r.total_act for r in f508["NT"].deduplicated]
    m_wt_all = sum(wt_nt_all) / len(wt_nt_all)
    m_wt_cl = sum(wt_nt_c_to_l) / len(wt_nt_c_to_l)
    m_f508_nt = sum(f508_nt) / len(f508_nt)
    ratio_all = m_wt_all / m_f508_nt
    ratio_cl = m_wt_cl / m_f508_nt
    if ratio_all <= MIN_WT_F508DEL_RANGE_MULT or ratio_cl <= MIN_WT_F508DEL_RANGE_MULT:
        raise RuntimeError(f"WT NT vs F508del NT dynamic range failed: {ratio_all}, {ratio_cl}")

    # 3. F508del siNS 661+445 > F508del siNS DMSO
    m_sins_combo = sum(r.total_act for r in f508["siNS, 661+445"].deduplicated) / len(f508["siNS, 661+445"].deduplicated)
    m_sins_dmso = sum(r.total_act for r in f508["siNS, DMSO"].deduplicated) / len(f508["siNS, DMSO"].deduplicated)
    if not (m_sins_combo > m_sins_dmso):
        raise RuntimeError(f"F508del combo rescue control failed: {m_sins_combo} <= {m_sins_dmso}")

    return {
        "inh172_reversal_parity": wb_checks,
        "wt_vs_f508del_nt_dynamic_range": {
            "frt_wt_nt_mean_c_to_r_n15": round(m_wt_all, 6),
            "frt_wt_nt_mean_c_to_l_n10": round(m_wt_cl, 6),
            "frt_f508del_nt_mean_n8": round(m_f508_nt, 6),
            "ratio_c_to_r": round(ratio_all, 6),
            "ratio_c_to_l": round(ratio_cl, 6),
            "threshold": MIN_WT_F508DEL_RANGE_MULT,
            "passed": True,
        },
        "f508del_sins_combo_gt_dmso": {
            "sins_661_445_mean": round(m_sins_combo, 6),
            "sins_dmso_mean": round(m_sins_dmso, 6),
            "fold_rescue": round(m_sins_combo / m_sins_dmso, 6),
            "passed": True,
        },
    }


def compute_dms_f508del_concordance(
    data_dir: Path,
    f508_ussing: Dict[str, UssingGroupData],
) -> Dict[str, Any]:
    """Non-gated sensitivity context: Spearman rank vs HEK293T DMS F508del (Row 213)."""
    dms_path = data_dir / MASTER_DMS_XLSX
    if not dms_path.is_file():
        return {"available": False}
    dms_sha = hashlib.sha256(dms_path.read_bytes()).hexdigest()
    if dms_sha != MASTER_DMS_SHA256:
        raise ValueError(f"DMS master xlsx SHA-256 mismatch: {dms_sha}")

    par_cells, _ = parse_xlsx_sheet(dms_path, "All_variants", expand_merged=True)
    ko_cells, _ = parse_xlsx_sheet(dms_path, "All_variants CANX KO", expand_merged=True)
    if par_cells.get((213, "C"), "").strip() != "F508del" or ko_cells.get((213, "C"), "").strip() != "F508del":
        raise ValueError("Expected F508del at Row 213 of DMS Master Tables")

    # Note: Parental All_variants Col J (VX-661) is blank (0/234) in published workbook.
    par_3 = {
        "DMSO": float(par_cells[(213, "F")]),
        "VX-445": float(par_cells[(213, "N")]),
        "661+445": float(par_cells[(213, "R")]),
    }
    ko_4 = {
        "DMSO": float(ko_cells[(213, "F")]),
        "VX-661": float(ko_cells[(213, "J")]),
        "VX-445": float(ko_cells[(213, "N")]),
        "661+445": float(ko_cells[(213, "R")]),
    }
    ussing_sins = {
        t: sum(r.total_act for r in f508_ussing[f"siNS, {t}"].deduplicated)
        / len(f508_ussing[f"siNS, {t}"].deduplicated)
        for t in ("DMSO", "VX-661", "VX-445", "661+445")
    }
    ussing_sicanx = {
        t: sum(r.total_act for r in f508_ussing[f"siCANX, {t}"].deduplicated)
        / len(f508_ussing[f"siCANX, {t}"].deduplicated)
        for t in ("DMSO", "VX-661", "VX-445", "661+445")
    }

    t3 = ("DMSO", "VX-445", "661+445")
    t4 = ("DMSO", "VX-661", "VX-445", "661+445")
    return {
        "available": True,
        "master_xlsx_sha256": dms_sha,
        "parental_col_j_vx661_status": "blank_0_of_234_in_published_workbook",
        "parental_dms_3_treatments": {k: round(v, 6) for k, v in par_3.items()},
        "canx_ko_dms_4_treatments": {k: round(v, 6) for k, v in ko_4.items()},
        "frt_f508del_sins_ussing_means": {k: round(v, 6) for k, v in ussing_sins.items()},
        "frt_f508del_sicanx_ussing_means": {k: round(v, 6) for k, v in ussing_sicanx.items()},
        "spearman_parental_dms_vs_sins_3tr": round(
            spearman_rank_correlation([par_3[t] for t in t3], [ussing_sins[t] for t in t3]), 6
        ),
        "spearman_canx_ko_dms_vs_sicanx_4tr": round(
            spearman_rank_correlation([ko_4[t] for t in t4], [ussing_sicanx[t] for t in t4]), 6
        ),
        "spearman_canx_ko_dms_vs_sins_4tr": round(
            spearman_rank_correlation([ko_4[t] for t in t4], [ussing_sins[t] for t in t4]), 6
        ),
    }


def parse_proteomics_workbook(
    xlsx_path: Path,
    workbook_label: str,
    expected_cols: Sequence[str],
) -> Tuple[List[ProteomicsRow], ProteomicsRow, List[Dict[str, Any]]]:
    """Parse LogFoldChanges Sheet1, normalizing header NBSP and auditing non-numeric cells."""
    cells, max_row = parse_xlsx_sheet(xlsx_path, "Sheet1", expand_merged=False)
    col_letters = [_index_to_col(i) for i in range(3, 3 + len(expected_cols))]
    hdr_gene = cells.get((2, "A"), "").replace("\xa0", " ").strip()
    hdr_pw = cells.get((2, "B"), "").replace("\xa0", " ").strip()
    if hdr_gene != "Gene Symbol" or hdr_pw != "Pathway":
        raise ValueError(f"Unexpected header A2/B2 in {xlsx_path.name}: {(hdr_gene, hdr_pw)!r}")
    for c_let, exp_col in zip(col_letters, expected_cols):
        got = cells.get((2, c_let), "").replace("\xa0", " ").strip()
        if got != exp_col:
            raise ValueError(f"Unexpected column {c_let}2 in {xlsx_path.name}: {got!r} != {exp_col!r}")

    interactors: List[ProteomicsRow] = []
    cftr_anchor: Optional[ProteomicsRow] = None
    non_numeric: List[Dict[str, Any]] = []
    seen_genes: set[str] = set()

    for rnum in range(3, max_row + 1):
        gene = cells.get((rnum, "A"), "").replace("\xa0", " ").strip()
        pw_raw = cells.get((rnum, "B"), "").replace("\xa0", " ").strip()
        raw_vals = [cells.get((rnum, c_let)) for c_let in col_letters]
        if not gene and not pw_raw and all(v is None or str(v).strip() == "" for v in raw_vals):
            continue
        if not gene:
            raise ValueError(f"Missing Gene Symbol at row {rnum} of {xlsx_path.name}")
        if gene in seen_genes:
            raise ValueError(f"Duplicate Gene Symbol {gene!r} in {xlsx_path.name}")
        seen_genes.add(gene)

        parsed_vals: Dict[str, Optional[float]] = {}
        for col_name, c_let, raw in zip(expected_cols, col_letters, raw_vals):
            if raw is None or str(raw).strip() == "":
                raise ValueError(f"Blank numeric cell {c_let}{rnum} ({gene}) in {xlsx_path.name}")
            s = str(raw).strip()
            try:
                num = float(s)
                if not math.isfinite(num):
                    raise ValueError(f"Non-finite value {s!r} at {c_let}{rnum}")
                parsed_vals[col_name] = num
            except ValueError:
                if s != "#DIV/0!":
                    raise ValueError(f"Unexpected non-numeric value {s!r} at {c_let}{rnum} ({gene})")
                non_numeric.append(
                    {
                        "workbook": workbook_label,
                        "row_num": rnum,
                        "gene_symbol": gene,
                        "column": col_name,
                        "cell": f"{c_let}{rnum}",
                        "raw_value": s,
                    }
                )
                parsed_vals[col_name] = None

        if len(non_numeric) > MAX_NON_NUMERIC_PER_WORKBOOK:
            raise RuntimeError(
                f"Non-numeric count in {xlsx_path.name} exceeded limit {MAX_NON_NUMERIC_PER_WORKBOOK}"
            )

        rec = ProteomicsRow(
            row_num=rnum,
            gene_symbol=gene,
            pathway_raw=pw_raw,
            pathway_norm=pw_raw.strip().lower(),
            values=parsed_vals,
        )
        if gene == "CFTR":
            cftr_anchor = rec
        else:
            interactors.append(rec)

    if cftr_anchor is None:
        raise ValueError(f"CFTR anchor row missing from {xlsx_path.name}")
    return interactors, cftr_anchor, non_numeric


def parse_tmt_labels_metadata(xlsx_path: Path) -> Dict[str, Any]:
    """Audit Copy of TMT_Labels.xlsx 18-plex x 6-run design table."""
    cells, max_row = parse_xlsx_sheet(xlsx_path, "Sheet1", expand_merged=False)
    if cells.get((1, "A"), "").strip() != "TMT 18-plex Labels":
        raise ValueError("Unexpected A1 header in Copy of TMT_Labels.xlsx")
    runs = [cells.get((1, _index_to_col(c)), "").strip() for c in range(2, 8)]
    if runs != [f"Run {i}" for i in range(1, 7)]:
        raise ValueError(f"Unexpected TMT run headers: {runs}")
    channels = [cells.get((r, "A"), "").strip() for r in range(2, 20)]
    if any(ch == "" for ch in channels) or max_row != 19:
        raise ValueError(f"Expected 18 non-empty TMT channels (rows 2..19), got {channels}, max_row={max_row}")
    populated_cells = sum(
        1 for r in range(2, 20) for c in range(2, 8) if cells.get((r, _index_to_col(c)), "").strip() != ""
    )
    if populated_cells != 108:
        raise ValueError(f"Expected 108 TMT design entries (18 x 6), got {populated_cells}")
    return {
        "runs": len(runs),
        "channels_per_run": len(channels),
        "total_channel_entries": populated_cells,
        "channels": channels,
    }


def wilcoxon_signed_rank_two_sided(diffs: Sequence[float]) -> Dict[str, Any]:
    """Two-sided Wilcoxon signed-rank test (exact doubled-rank DP for n<=25, tie-corrected normal >25)."""
    nz = [d for d in diffs if d != 0.0]
    n = len(nz)
    if n == 0:
        return {"n_total": len(diffs), "n_nonzero": 0, "w_plus": 0.0, "p_two_sided": 1.0, "method": "degenerate"}
    abs_vals = [abs(d) for d in nz]
    ranks, tie_sizes = _rank_average_ties(abs_vals)
    w_plus = sum(r for d, r in zip(nz, ranks) if d > 0.0)
    total_rank_sum = n * (n + 1) / 2.0
    w_min = min(w_plus, total_rank_sum - w_plus)

    if n <= 25:
        # Exact symmetric signed-rank DP using doubled integer ranks (all midranks * 2 are integers)
        doubled = [int(round(2.0 * r)) for r in ranks]
        s_double = sum(doubled)
        w_obs_double = int(round(2.0 * w_min))
        counts = [1] + [0] * s_double
        running_max = 0
        for dr in doubled:
            for s in range(running_max, -1, -1):
                c = counts[s]
                if c:
                    counts[s + dr] += c
            running_max += dr
        total_ways = 1 << n
        half_s = s_double / 2.0
        obs_dev = abs(w_obs_double - half_s)
        extreme = sum(c for s, c in enumerate(counts) if abs(s - half_s) >= obs_dev - 1e-9)
        p_val = min(1.0, extreme / total_ways)
        method = "exact_dp"
    else:
        mu = total_rank_sum / 2.0
        tie_term = sum(t * (t * t - 1) for t in tie_sizes)
        var = (n * (n + 1) * (2 * n + 1) / 24.0) - (tie_term / 48.0)
        if var <= 0.0:
            return {"n_total": len(diffs), "n_nonzero": n, "w_plus": w_plus, "p_two_sided": 1.0, "method": "zero_var"}
        diff_mu = abs(w_plus - mu)
        z = max(0.0, diff_mu - 0.5) / math.sqrt(var)
        p_val = min(1.0, max(0.0, math.erfc(z / math.sqrt(2.0))))
        method = "normal_tie_corrected"

    return {
        "n_total": len(diffs),
        "n_nonzero": n,
        "w_plus": round(w_plus, 6),
        "w_min": round(w_min, 6),
        "p_two_sided": p_val,
        "method": method,
    }


def mann_whitney_u_two_sided(group_vals: Sequence[float], bg_vals: Sequence[float]) -> Dict[str, Any]:
    """Two-sided Mann-Whitney U test using tie-corrected normal approximation."""
    n1, n2 = len(group_vals), len(bg_vals)
    if n1 < 1 or n2 < 1:
        raise ValueError("Mann-Whitney U requires non-empty groups")
    pooled = list(group_vals) + list(bg_vals)
    ranks, tie_sizes = _rank_average_ties(pooled)
    return _mann_whitney_from_ranks(ranks[:n1], n1, n2, tie_sizes)


def _mann_whitney_from_ranks(
    group_ranks: Sequence[float],
    n1: int,
    n2: int,
    tie_sizes: Sequence[int],
) -> Dict[str, Any]:
    r1 = sum(group_ranks)
    u1 = r1 - (n1 * (n1 + 1)) / 2.0
    n_tot = n1 + n2
    mu = (n1 * n2) / 2.0
    tie_sum = sum(t * (t * t - 1) for t in tie_sizes)
    var = (n1 * n2 / 12.0) * ((n_tot + 1.0) - (tie_sum / (n_tot * (n_tot - 1.0))))
    if var <= 0.0:
        return {"n_group": n1, "n_bg": n2, "u1": u1, "z_signed": 0.0, "p_two_sided": 1.0}
    z_signed = (u1 - mu) / math.sqrt(var)
    z_cc = max(0.0, abs(u1 - mu) - 0.5) / math.sqrt(var)
    p_val = min(1.0, max(0.0, math.erfc(z_cc / math.sqrt(2.0))))
    return {
        "n_group": n1,
        "n_bg": n2,
        "u1": round(u1, 6),
        "z_signed": round(z_signed, 6),
        "p_two_sided": p_val,
    }


def benjamini_hochberg(pvals: Sequence[float]) -> List[float]:
    """Benjamini-Hochberg step-up FDR q-values preserving input order."""
    m = len(pvals)
    if m == 0:
        return []
    indexed = sorted(enumerate(pvals), key=lambda kv: kv[1])
    q_sorted = [0.0] * m
    running_min = 1.0
    for rank_0 in range(m - 1, -1, -1):
        rank_1 = rank_0 + 1
        val = min(1.0, indexed[rank_0][1] * m / rank_1)
        if val < running_min:
            running_min = val
        q_sorted[rank_0] = running_min
    out = [0.0] * m
    for pos, (orig_idx, _) in enumerate(indexed):
        out[orig_idx] = q_sorted[pos]
    return out


def _canonical_module_name(pw_norm: str) -> Optional[str]:
    for mod_name, terms in CANONICAL_MODULES.items():
        if pw_norm in terms:
            return mod_name
    return None


def _fast_wilcoxon_p_for_null(diffs: Sequence[float]) -> float:
    """Fast calibrated Wilcoxon p-value used in 10,000-draw pathway permutation null."""
    nz = [d for d in diffs if d != 0.0]
    n = len(nz)
    if n == 0:
        return 1.0
    if n <= 25:
        return wilcoxon_signed_rank_two_sided(nz)["p_two_sided"]
    abs_vals = [abs(d) for d in nz]
    ranks, tie_sizes = _rank_average_ties(abs_vals)
    w_plus = sum(r for d, r in zip(nz, ranks) if d > 0.0)
    mu = n * (n + 1) / 4.0
    tie_term = sum(t * (t * t - 1) for t in tie_sizes)
    var = (n * (n + 1) * (2 * n + 1) / 24.0) - (tie_term / 48.0)
    if var <= 0.0:
        return 1.0
    z = max(0.0, abs(w_plus - mu) - 0.5) / math.sqrt(var)
    return min(1.0, max(0.0, math.erfc(z / math.sqrt(2.0))))


def run_arm2_controls(
    var_rows: Sequence[ProteomicsRow],
    var_anchor: ProteomicsRow,
    drug_rows: Sequence[ProteomicsRow],
    drug_anchor: ProteomicsRow,
    tmt_meta: Dict[str, Any],
) -> Dict[str, Any]:
    """Run Arm 2 CFTR anchor, dispersion, TMT parity, and 10,000-draw permutation null controls."""
    # 1. CFTR anchor == 0.0 and non-zero dispersion across all 9 columns
    anchor_checks: Dict[str, Any] = {}
    dispersion_checks: Dict[str, Any] = {}
    for col in VARIANT_COLS:
        av = var_anchor.values[col]
        if av is None or abs(av) > 1e-12:
            raise RuntimeError(f"Variants CFTR anchor non-zero in {col}: {av}")
        vals = [r.values[col] for r in var_rows if r.values[col] is not None]
        sd = statistics.stdev(vals)
        if sd <= MIN_DISPERSION_STD:
            raise RuntimeError(f"Variants dispersion control failed in {col}: std={sd}")
        anchor_checks[col] = av
        dispersion_checks[col] = round(sd, 6)

    for col in DRUG_COLS:
        av = drug_anchor.values[col]
        if av is None or abs(av) > 1e-12:
            raise RuntimeError(f"Drug CFTR anchor non-zero in {col}: {av}")
        vals = [r.values[col] for r in drug_rows if r.values[col] is not None]
        sd = statistics.stdev(vals)
        if sd <= MIN_DISPERSION_STD:
            raise RuntimeError(f"Drug dispersion control failed in {col}: std={sd}")
        anchor_checks[col] = av
        dispersion_checks[col] = round(sd, 6)

    # 2. 10,000-draw pathway-label permutation null (seed = SEED)
    rng = random.Random(SEED)

    # Drug paired attenuation null (4 tests BH-FDR q < 0.05)
    drug_pairs_idx: List[int] = []
    drug_diffs_all: List[float] = []
    drug_mod_labels: List[Optional[str]] = []
    for idx, r in enumerate(drug_rows):
        v_dmso = r.values["ΔF (-CX)"]
        v_combo = r.values["DelF (-CX) + combo"]
        if v_dmso is not None and v_combo is not None:
            drug_pairs_idx.append(idx)
            drug_diffs_all.append(abs(v_combo) - abs(v_dmso))
            drug_mod_labels.append(_canonical_module_name(r.pathway_norm))

    # Variants precomputed midranks per variant column
    var_precomputed: List[Tuple[List[int], List[float], List[int], float]] = []
    for col in VARIANT_COLS:
        finite_idx = [i for i, r in enumerate(var_rows) if r.values[col] is not None]
        vals = [float(var_rows[i].values[col]) for i in finite_idx]
        ranks, tie_sizes = _rank_average_ties(vals)
        n_tot = len(vals)
        tie_sum = sum(t * (t * t - 1) for t in tie_sizes)
        # Map row index -> rank
        row_to_rank = {ri: rk for ri, rk in zip(finite_idx, ranks)}
        full_ranks = [row_to_rank.get(i, 0.0) for i in range(len(var_rows))]
        var_precomputed.append((finite_idx, full_ranks, tie_sizes, float(tie_sum)))

    var_labels = [_canonical_module_name(r.pathway_norm) for r in var_rows]
    # Permute only non-None module labels or all row labels (fixing exact module counts)
    drug_fwer_hits = 0
    var_fwer_hits = 0
    for _ in range(PERM_NULL_DRAWS):
        shuffled_d = drug_mod_labels[:]
        rng.shuffle(shuffled_d)
        pvals_d: List[float] = []
        for mod in CANONICAL_MODULE_ORDER:
            sub = [drug_diffs_all[k] for k, lbl in enumerate(shuffled_d) if lbl == mod]
            pvals_d.append(_fast_wilcoxon_p_for_null(sub))
        if any(q < ALPHA_PROTEOMICS_FDR for q in benjamini_hochberg(pvals_d)):
            drug_fwer_hits += 1

        shuffled_v = var_labels[:]
        rng.shuffle(shuffled_v)
        pvals_v: List[float] = []
        for finite_idx, full_ranks, _, tie_sum in var_precomputed:
            n_tot = len(finite_idx)
            mod_sums: Dict[str, float] = {m: 0.0 for m in CANONICAL_MODULE_ORDER}
            mod_cnts: Dict[str, int] = {m: 0 for m in CANONICAL_MODULE_ORDER}
            for ri in finite_idx:
                lbl = shuffled_v[ri]
                if lbl is not None:
                    mod_sums[lbl] += full_ranks[ri]
                    mod_cnts[lbl] += 1
            for mod in CANONICAL_MODULE_ORDER:
                n1 = mod_cnts[mod]
                n2 = n_tot - n1
                if n1 <= 0 or n2 <= 0 or n_tot <= 1:
                    pvals_v.append(1.0)
                    continue
                r1 = mod_sums[mod]
                u1 = r1 - (n1 * (n1 + 1)) / 2.0
                mu = (n1 * n2) / 2.0
                var = (n1 * n2 / 12.0) * ((n_tot + 1.0) - (tie_sum / (n_tot * (n_tot - 1.0))))
                if var <= 0.0:
                    pvals_v.append(1.0)
                    continue
                z_cc = max(0.0, abs(u1 - mu) - 0.5) / math.sqrt(var)
                pvals_v.append(min(1.0, max(0.0, math.erfc(z_cc / math.sqrt(2.0)))))
        if any(q < ALPHA_PROTEOMICS_FDR for q in benjamini_hochberg(pvals_v)):
            var_fwer_hits += 1

    drug_fwer = drug_fwer_hits / PERM_NULL_DRAWS
    var_fwer = var_fwer_hits / PERM_NULL_DRAWS
    if drug_fwer > MAX_NULL_FWER or var_fwer > MAX_NULL_FWER:
        raise RuntimeError(
            f"Arm 2 permutation negative control FWER exceeded {MAX_NULL_FWER}: drug={drug_fwer}, var={var_fwer}"
        )

    return {
        "cftr_anchor_all_zero": anchor_checks,
        "non_cftr_dispersion_std": dispersion_checks,
        "tmt_18plex_design_parity": tmt_meta,
        "pathway_permutation_null": {
            "draws": PERM_NULL_DRAWS,
            "seed": SEED,
            "drug_family_empirical_fwer_q05": round(drug_fwer, 6),
            "variants_family_empirical_fwer_q05": round(var_fwer, 6),
            "max_allowed_fwer": MAX_NULL_FWER,
            "passed": True,
        },
    }


def evaluate_arm2_proteomics(
    var_rows: Sequence[ProteomicsRow],
    drug_rows: Sequence[ProteomicsRow],
) -> Dict[str, Any]:
    """Evaluate Family 2a (Drug attenuation), Family 2b (Variants enrichment), and broad sensitivity."""
    # Family 2a: 4 canonical modules on paired |LFC(combo)| - |LFC(DMSO)|
    drug_tests: List[Dict[str, Any]] = []
    drug_pvals: List[float] = []
    for mod in CANONICAL_MODULE_ORDER:
        abs_diffs: List[float] = []
        signed_diffs: List[float] = []
        dmso_vals: List[float] = []
        combo_vals: List[float] = []
        for r in drug_rows:
            if _canonical_module_name(r.pathway_norm) != mod:
                continue
            vd = r.values["ΔF (-CX)"]
            vc = r.values["DelF (-CX) + combo"]
            if vd is None or vc is None:
                continue
            dmso_vals.append(vd)
            combo_vals.append(vc)
            abs_diffs.append(abs(vc) - abs(vd))
            signed_diffs.append(vc - vd)
        w_abs = wilcoxon_signed_rank_two_sided(abs_diffs)
        w_signed = wilcoxon_signed_rank_two_sided(signed_diffs)
        drug_pvals.append(w_abs["p_two_sided"])
        drug_tests.append(
            {
                "module": mod,
                "n_paired": len(abs_diffs),
                "median_dmso_lfc": round(statistics.median(dmso_vals), 6),
                "median_combo_lfc": round(statistics.median(combo_vals), 6),
                "median_abs_diff_combo_minus_dmso": round(statistics.median(abs_diffs), 6),
                "mean_abs_diff_combo_minus_dmso": round(sum(abs_diffs) / len(abs_diffs), 6),
                "wilcoxon_abs_method": w_abs["method"],
                "wilcoxon_abs_w_plus": w_abs["w_plus"],
                "raw_p_two_sided": w_abs["p_two_sided"],
                "signed_shift_context": {
                    "median_signed_diff": round(statistics.median(signed_diffs), 6),
                    "p_two_sided": round(w_signed["p_two_sided"], 6),
                },
            }
        )

    drug_qs = benjamini_hochberg(drug_pvals)
    for row, q in zip(drug_tests, drug_qs):
        row["raw_p_two_sided"] = round(float(row["raw_p_two_sided"]), 6)
        row["bh_fdr_q"] = round(q, 6)
        row["significant_q05"] = q < ALPHA_PROTEOMICS_FDR

    # Family 2b: 20 tests (4 canonical modules x 5 variants) on signed LFC vs background
    variant_tests: List[Dict[str, Any]] = []
    var_pvals: List[float] = []
    for col in VARIANT_COLS:
        finite_rows = [r for r in var_rows if r.values[col] is not None]
        for mod in CANONICAL_MODULE_ORDER:
            grp = [float(r.values[col]) for r in finite_rows if _canonical_module_name(r.pathway_norm) == mod]
            bg = [float(r.values[col]) for r in finite_rows if _canonical_module_name(r.pathway_norm) != mod]
            mwu = mann_whitney_u_two_sided(grp, bg)
            mwu_abs = mann_whitney_u_two_sided([abs(v) for v in grp], [abs(v) for v in bg])
            var_pvals.append(mwu["p_two_sided"])
            variant_tests.append(
                {
                    "variant": col,
                    "module": mod,
                    "n_module": len(grp),
                    "n_background": len(bg),
                    "median_module_lfc": round(statistics.median(grp), 6),
                    "median_background_lfc": round(statistics.median(bg), 6),
                    "mean_module_lfc": round(sum(grp) / len(grp), 6),
                    "mean_background_lfc": round(sum(bg) / len(bg), 6),
                    "u1": mwu["u1"],
                    "z_signed": mwu["z_signed"],
                    "raw_p_two_sided": mwu["p_two_sided"],
                    "abs_lfc_context": {
                        "median_module_abs_lfc": round(statistics.median([abs(v) for v in grp]), 6),
                        "median_bg_abs_lfc": round(statistics.median([abs(v) for v in bg]), 6),
                        "z_signed": mwu_abs["z_signed"],
                        "p_two_sided": round(mwu_abs["p_two_sided"], 6),
                    },
                }
            )

    var_qs = benjamini_hochberg(var_pvals)
    for row, q in zip(variant_tests, var_qs):
        row["raw_p_two_sided"] = round(float(row["raw_p_two_sided"]), 6)
        row["bh_fdr_q"] = round(q, 6)
        row["significant_q05"] = q < ALPHA_PROTEOMICS_FDR

    drug_sig_count = sum(1 for row in drug_tests if row["significant_q05"])
    var_sig_count = sum(1 for row in variant_tests if row["significant_q05"])
    arm2_supported = (drug_sig_count >= 1) and (var_sig_count >= 1)

    # Non-gated sensitivity context: broad proteostasis union
    broad_drug_diffs = [
        abs(float(r.values["DelF (-CX) + combo"])) - abs(float(r.values["ΔF (-CX)"]))
        for r in drug_rows
        if r.pathway_norm in BROAD_PROTEOSTASIS_SET
        and r.values["ΔF (-CX)"] is not None
        and r.values["DelF (-CX) + combo"] is not None
    ]
    broad_drug_w = wilcoxon_signed_rank_two_sided(broad_drug_diffs)
    broad_var_tests: List[Dict[str, Any]] = []
    broad_var_pvals: List[float] = []
    for col in VARIANT_COLS:
        finite_rows = [r for r in var_rows if r.values[col] is not None]
        grp = [float(r.values[col]) for r in finite_rows if r.pathway_norm in BROAD_PROTEOSTASIS_SET]
        bg = [float(r.values[col]) for r in finite_rows if r.pathway_norm not in BROAD_PROTEOSTASIS_SET]
        mwu = mann_whitney_u_two_sided(grp, bg)
        broad_var_pvals.append(mwu["p_two_sided"])
        broad_var_tests.append(
            {
                "variant": col,
                "n_broad": len(grp),
                "n_background": len(bg),
                "median_broad_lfc": round(statistics.median(grp), 6),
                "median_bg_lfc": round(statistics.median(bg), 6),
                "z_signed": mwu["z_signed"],
                "raw_p_two_sided": mwu["p_two_sided"],
            }
        )
    for row, q in zip(broad_var_tests, benjamini_hochberg(broad_var_pvals)):
        row["raw_p_two_sided"] = round(float(row["raw_p_two_sided"]), 6)
        row["bh_fdr_q"] = round(q, 6)

    return {
        "family_2a_drug_attenuation": {
            "num_tests": len(drug_tests),
            "fdr_threshold": ALPHA_PROTEOMICS_FDR,
            "significant_count": drug_sig_count,
            "tests": drug_tests,
        },
        "family_2b_variant_enrichment": {
            "num_tests": len(variant_tests),
            "fdr_threshold": ALPHA_PROTEOMICS_FDR,
            "significant_count": var_sig_count,
            "tests": variant_tests,
        },
        "overall_verdict": "supported" if arm2_supported else "refuted",
        "ledger_standing": "open" if arm2_supported else "refuted",
        "broad_proteostasis_sensitivity_context": {
            "included_normalized_pathways": sorted(BROAD_PROTEOSTASIS_SET),
            "drug_paired_abs_attenuation": {
                "n_paired": len(broad_drug_diffs),
                "median_abs_diff": round(statistics.median(broad_drug_diffs), 6),
                "w_plus": broad_drug_w["w_plus"],
                "p_two_sided": round(broad_drug_w["p_two_sided"], 6),
                "method": broad_drug_w["method"],
            },
            "variant_enrichment_5_tests": broad_var_tests,
        },
    }


def evaluate_benchmark(data_dir: Path, repo_root: Path = REPO_ROOT) -> Dict[str, Any]:
    """Run full pre-registered Arm 1 + Arm 2 evaluation on Tedman et al. 2025 workbooks."""
    done_path = repo_root / "benchmarks" / "cftr_canx_interactome_ussing" / "done.md"
    done_sha = hashlib.sha256(
        done_path.read_text(encoding="utf-8", errors="replace").replace("\r\n", "\n").encode("utf-8")
    ).hexdigest()
    if done_sha != DONE_SHA256_EXPECTED:
        raise AssertionError(f"done_sha256 mismatch: {done_sha} != {DONE_SHA256_EXPECTED}")

    digests = verify_required_workbooks(data_dir)

    # Parse Arm 1 raw Ussing sheets
    f508_path = data_dir / "Copy of FRT-F508del_siCANX_Ussing-summary-for-paper.xlsx"
    wt_path = data_dir / "Copy of FRT-WT_siCANX_Ussing-summary-for-paper.xlsx"
    f508_groups, f508_dups = parse_ussing_workbook(f508_path, "FRT-F508del", "M")
    wt_groups, wt_dups = parse_ussing_workbook(wt_path, "FRT-WT", "S")
    all_dups = f508_dups + wt_dups

    arm1_controls = run_arm1_positive_controls(f508_groups, wt_groups, wt_path)
    arm1_primary = _evaluate_ussing_pair(f508_groups, wt_groups, use_dedup=True)
    arm1_n4_sens = _evaluate_ussing_pair(f508_groups, wt_groups, use_dedup=False)
    dms_concordance = compute_dms_f508del_concordance(data_dir, f508_groups)

    # Parse Arm 2 raw proteomics sheets + TMT design metadata
    var_path = data_dir / "Copy of CFTR_Variants_LogFoldChanges.xlsx"
    drug_path = data_dir / "Copy of CFTR_F508del_Drug_LogFoldChanges.xlsx"
    tmt_path = data_dir / "Copy of TMT_Labels.xlsx"

    var_rows, var_anchor, var_non_num = parse_proteomics_workbook(var_path, "CFTR_Variants", VARIANT_COLS)
    drug_rows, drug_anchor, drug_non_num = parse_proteomics_workbook(drug_path, "CFTR_F508del_Drug", DRUG_COLS)
    tmt_meta = parse_tmt_labels_metadata(tmt_path)

    arm2_controls = run_arm2_controls(var_rows, var_anchor, drug_rows, drug_anchor, tmt_meta)
    arm2_primary = evaluate_arm2_proteomics(var_rows, drug_rows)

    return {
        "schema_version": "tedman2025-canx-interactome-ussing-results-v1",
        "benchmark_id": "tedman2025-cftr-canx-interactome-ussing",
        "evaluated_on_local": str(dt.date.today()),
        "evaluated_on_utc": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "done_sha256": done_sha,
        "dataset": {
            "source_id": "cftr-dms-tedman-2025",
            "paper_doi": "10.7554/eLife.107180",
            "mendeley_doi": "10.17632/ksfmxv8zts.1",
            "license": "CC BY 4.0",
            "research_boundary": (
                "Macroscopic FRT Ussing Delta Isc and HEK293T TMT-18plex co-IP fold changes are not "
                "single-channel patch-clamp burst durations, clinical endpoints, or medical advice."
            ),
            "workbook_sha256": digests,
            "variants_non_empty_rows_incl_cftr": len(var_rows) + 1,
            "variants_interactors_excl_cftr": len(var_rows),
            "drug_non_empty_rows_incl_cftr": len(drug_rows) + 1,
            "drug_interactors_excl_cftr": len(drug_rows),
        },
        "data_integrity_audit": {
            "duplicate_replicate_columns_count": len(all_dups),
            "duplicate_replicate_columns": all_dups,
            "proteomics_non_numeric_cells": var_non_num + drug_non_num,
        },
        "arm_1_frt_ussing_sicanx_asymmetry": {
            "positive_controls": arm1_controls,
            "primary_strict_dedup": arm1_primary,
            "group_summaries_dedup": {
                "frt_f508del": _summarize_ussing_groups(f508_groups, use_dedup=True),
                "frt_wt": _summarize_ussing_groups(wt_groups, use_dedup=True),
            },
            "sensitivity_as_deposited_n4": arm1_n4_sens,
            "dms_f508del_concordance_context": dms_concordance,
        },
        "arm_2_canx_ko_proteostasis_remodeling": {
            "controls": arm2_controls,
            "primary_evaluation": arm2_primary,
        },
    }


def verify_checked_in_artifacts(repo_root: Path = REPO_ROOT, data_dir: Optional[Path] = None) -> Dict[str, Any]:
    """Verify done_sha256, JSON invariants, BH-FDR consistency, and raw .xlsx parity when present."""
    bench_dir = repo_root / "benchmarks" / "cftr_canx_interactome_ussing"
    done_path = bench_dir / "done.md"
    results_path = bench_dir / RESULTS_FILE_NAME

    done_sha = hashlib.sha256(
        done_path.read_text(encoding="utf-8", errors="replace").replace("\r\n", "\n").encode("utf-8")
    ).hexdigest()
    if done_sha != DONE_SHA256_EXPECTED:
        raise AssertionError(f"done_sha256 mismatch: {done_sha} != {DONE_SHA256_EXPECTED}")
    results = json.loads(results_path.read_text(encoding="utf-8"))
    if results["done_sha256"] != done_sha:
        raise AssertionError(f"results done_sha256 mismatch: {results['done_sha256']} != {done_sha}")

    # Verify internal BH-FDR, ledger, & verdict logic even in offline CI without .xlsx
    a1p = results["arm_1_frt_ussing_sicanx_asymmetry"]["primary_strict_dedup"]
    exp_a1a = (
        a1p["arm_1a_f508del_combo_sicanx_drop"]["exact_perm_p_two_sided"] <= ALPHA_USSING
        and a1p["arm_1a_f508del_combo_sicanx_drop"]["bootstrap_log_ratio"]["ci975_high_lt_zero"] is True
    )
    exp_a1b = (
        a1p["arm_1b_sicanx_proportion_f508del_vs_wt"]["exact_stratified_perm_p_one_sided"] <= ALPHA_USSING
        and a1p["arm_1b_sicanx_proportion_f508del_vs_wt"]["bootstrap_did_log_ratio"]["one_sided_975_high_lt_zero"] is True
    )
    if a1p["arm_1a_f508del_combo_sicanx_drop"]["passed_gate"] != exp_a1a:
        raise AssertionError("Arm 1a gate state mismatch")
    if a1p["arm_1b_sicanx_proportion_f508del_vs_wt"]["passed_gate"] != exp_a1b:
        raise AssertionError("Arm 1b gate state mismatch")
    exp_a1_verdict = "supported" if (exp_a1a and exp_a1b) else "refuted"
    exp_a1_standing = "open" if (exp_a1a and exp_a1b) else "refuted"
    if a1p["overall_verdict"] != exp_a1_verdict or a1p["ledger_standing"] != exp_a1_standing:
        raise AssertionError("Arm 1 overall verdict or standing mismatch")

    a2p = results["arm_2_canx_ko_proteostasis_remodeling"]["primary_evaluation"]
    d_fam = a2p["family_2a_drug_attenuation"]
    v_fam = a2p["family_2b_variant_enrichment"]
    d_tests = d_fam["tests"]
    v_tests = v_fam["tests"]
    if len(d_tests) != 4 or len(v_tests) != 20:
        raise AssertionError("Unexpected test counts in Arm 2 results")

    d_q_chk = benjamini_hochberg([float(t["raw_p_two_sided"]) for t in d_tests])
    for t, qc in zip(d_tests, d_q_chk):
        if abs(float(t["bh_fdr_q"]) - qc) > len(d_tests) * 0.5e-6 + 1e-6:
            raise AssertionError(f"Offline Drug BH-FDR mismatch for {t['module']}: {t['bh_fdr_q']} vs {qc}")
        if t["significant_q05"] != (float(t["bh_fdr_q"]) < ALPHA_PROTEOMICS_FDR):
            raise AssertionError(f"Drug significant_q05 flag mismatch for {t['module']}")

    v_q_chk = benjamini_hochberg([float(t["raw_p_two_sided"]) for t in v_tests])
    for t, qc in zip(v_tests, v_q_chk):
        if abs(float(t["bh_fdr_q"]) - qc) > len(v_tests) * 0.5e-6 + 1e-6:
            raise AssertionError(f"Offline Variants BH-FDR mismatch for {t['variant']}:{t['module']}")
        if t["significant_q05"] != (float(t["bh_fdr_q"]) < ALPHA_PROTEOMICS_FDR):
            raise AssertionError(f"Variants significant_q05 flag mismatch for {t['variant']}:{t['module']}")

    d_sig = sum(1 for t in d_tests if t["significant_q05"])
    v_sig = sum(1 for t in v_tests if t["significant_q05"])
    if d_fam["significant_count"] != d_sig or v_fam["significant_count"] != v_sig:
        raise AssertionError("Arm 2 significant_count mismatch")
    exp_a2_verdict = "supported" if (d_sig >= 1 and v_sig >= 1) else "refuted"
    exp_a2_standing = "open" if (d_sig >= 1 and v_sig >= 1) else "refuted"
    if a2p["overall_verdict"] != exp_a2_verdict or a2p["ledger_standing"] != exp_a2_standing:
        raise AssertionError("Arm 2 overall verdict or standing mismatch")

    ledger = json.loads((bench_dir / "ledger.json").read_text(encoding="utf-8"))
    if ledger.get("done_sha256") != done_sha:
        raise AssertionError("ledger.json done_sha256 mismatch")
    by_id = {c["id"]: c for c in ledger.get("claims", [])}
    if by_id["c4-arm-1-ussing-sicanx-asymmetry"]["standing"] != exp_a1_standing:
        raise AssertionError("ledger Arm 1 standing mismatch")
    if by_id["c5-arm-2-proteostasis-remodeling"]["standing"] != exp_a2_standing:
        raise AssertionError("ledger Arm 2 standing mismatch")

    resolved_data = resolve_tedman_data_dir(data_dir, repo_root=repo_root)
    if data_dir is not None and resolved_data is None:
        raise FileNotFoundError(f"Explicit --data-dir invalid or missing Tedman workbooks: {data_dir}")
    raw_parity_checked = False
    if resolved_data is not None:
        fresh = evaluate_benchmark(resolved_data, repo_root=repo_root)
        if fresh["data_integrity_audit"]["duplicate_replicate_columns_count"] != results["data_integrity_audit"]["duplicate_replicate_columns_count"]:
            raise AssertionError("Duplicate replicate column count parity drift")
        f_p = fresh["arm_1_frt_ussing_sicanx_asymmetry"]["primary_strict_dedup"]
        if f_p["overall_verdict"] != a1p["overall_verdict"] or f_p["ledger_standing"] != a1p["ledger_standing"]:
            raise AssertionError("Arm 1 verdict parity drift")
        f_a1a = f_p["arm_1a_f508del_combo_sicanx_drop"]
        r_a1a = a1p["arm_1a_f508del_combo_sicanx_drop"]
        if (
            abs(f_a1a["exact_perm_p_two_sided"] - r_a1a["exact_perm_p_two_sided"]) > 1e-6
            or abs(f_a1a["mean_sicanx_total_act"] - r_a1a["mean_sicanx_total_act"]) > 1e-6
            or f_a1a["passed_gate"] != r_a1a["passed_gate"]
            or f_a1a["bootstrap_log_ratio"]["ci975"] != r_a1a["bootstrap_log_ratio"]["ci975"]
        ):
            raise AssertionError("Arm 1a raw parity drift")
        f_a1b = f_p["arm_1b_sicanx_proportion_f508del_vs_wt"]
        r_a1b = a1p["arm_1b_sicanx_proportion_f508del_vs_wt"]
        if (
            abs(f_a1b["did_log_ratio"] - r_a1b["did_log_ratio"]) > 1e-6
            or abs(f_a1b["exact_stratified_perm_p_one_sided"] - r_a1b["exact_stratified_perm_p_one_sided"]) > 1e-6
            or f_a1b["passed_gate"] != r_a1b["passed_gate"]
            or f_a1b["bootstrap_did_log_ratio"]["one_sided_975_upper"] != r_a1b["bootstrap_did_log_ratio"]["one_sided_975_upper"]
        ):
            raise AssertionError("Arm 1b raw parity drift")
        f_a2 = fresh["arm_2_canx_ko_proteostasis_remodeling"]["primary_evaluation"]
        if f_a2["overall_verdict"] != a2p["overall_verdict"] or f_a2["ledger_standing"] != a2p["ledger_standing"]:
            raise AssertionError("Arm 2 verdict parity drift")
        f_d = f_a2["family_2a_drug_attenuation"]["tests"]
        for ft, rt in zip(f_d, d_tests):
            if abs(ft["bh_fdr_q"] - rt["bh_fdr_q"]) > 1e-6 or abs(ft["raw_p_two_sided"] - rt["raw_p_two_sided"]) > 1e-6:
                raise AssertionError(f"Arm 2a BH-FDR parity drift for {rt['module']}")
        f_v = f_a2["family_2b_variant_enrichment"]["tests"]
        for ft, rt in zip(f_v, v_tests):
            if abs(ft["bh_fdr_q"] - rt["bh_fdr_q"]) > 1e-6 or abs(ft["raw_p_two_sided"] - rt["raw_p_two_sided"]) > 1e-6:
                raise AssertionError(f"Arm 2b BH-FDR parity drift for {rt['variant']}:{rt['module']}")
        raw_parity_checked = True

    return {
        "status": "PASS",
        "done_sha256": done_sha,
        "raw_local_data_parity_checked": raw_parity_checked,
        "duplicate_replicate_columns_count": results["data_integrity_audit"]["duplicate_replicate_columns_count"],
        "arm_1_verdict": a1p["overall_verdict"],
        "arm_2_verdict": a2p["overall_verdict"],
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
                f"PASS cftr_canx_interactome_ussing: raw_parity={summary['raw_local_data_parity_checked']} "
                f"dup_cols={summary['duplicate_replicate_columns_count']} "
                f"arm1={summary['arm_1_verdict']} arm2={summary['arm_2_verdict']}"
            )
        return 0

    resolved_data = resolve_tedman_data_dir(args.data_dir, repo_root=REPO_ROOT)
    if resolved_data is None:
        print("Local Tedman 2025 files not present under sources/downloads/tedman-2025/.", file=sys.stderr)
        return 1

    results = evaluate_benchmark(resolved_data, repo_root=REPO_ROOT)
    if args.write_artifacts:
        (HERE / RESULTS_FILE_NAME).write_text(
            json.dumps(results, indent=2, allow_nan=False) + "\n", encoding="utf-8"
        )

    if args.json:
        print(json.dumps(results, indent=2, allow_nan=False))
    else:
        a1 = results["arm_1_frt_ussing_sicanx_asymmetry"]["primary_strict_dedup"]
        a2 = results["arm_2_canx_ko_proteostasis_remodeling"]["primary_evaluation"]
        print(
            f"Arm 1: verdict={a1['overall_verdict']} (1a_pass={a1['arm_1a_f508del_combo_sicanx_drop']['passed_gate']}, "
            f"1b_pass={a1['arm_1b_sicanx_proportion_f508del_vs_wt']['passed_gate']}) | "
            f"Arm 2: verdict={a2['overall_verdict']} (drug_sig={a2['family_2a_drug_attenuation']['significant_count']}/4, "
            f"var_sig={a2['family_2b_variant_enrichment']['significant_count']}/20)"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
