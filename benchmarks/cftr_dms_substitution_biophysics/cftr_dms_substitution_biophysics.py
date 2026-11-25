#!/usr/bin/env python3
"""Substitution biophysics + UniProt P13569 topology DMS rescue benchmark on Tedman et al. 2025.

Evaluates whether 7 sequence-substitution biophysics and UniProt P13569 topology features
(delta_blosum62, delta_kd_hydropathy, abs_delta_charge, helix_breaker_in_helix,
tm_polarity_shock, is_nonsense, is_inframe_indel) computed from each variant's experimental WT
residue improve log2(VX-661+VX-445 Mut/WT) prediction over baseline untreated expression OLS on
the locked Tedman et al. 2025 residue-clustered holdout.

Pure Python 3.10+ standard library only; zero network calls; raw `.xlsx` stays git-ignored.
Research benchmark only: plasma-membrane expression is not channel conductance or clinical benefit.
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
from typing import Any, Dict, FrozenSet, List, Mapping, Optional, Sequence, Tuple
import xml.etree.ElementTree as ET
import zipfile

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]

MASTER_XLSX_NAME = "ALT 9-6-24 CFTR DMS Expanded Library Master Tables.xlsx"
MASTER_XLSX_SHA256 = "2a702ee3744b8039a9c7c4c74aec9bb30a91301b138d4948e0e2c81a1931e51b"
DONE_SHA256_EXPECTED = "f956d9440f9a46e5de7569ab91a491ebd0e9793ebc94469a7b1885b2efea68c1"

SPLIT_MANIFEST_REL = Path("benchmarks/cftr_dms_expression/split_manifest.json")
SPLIT_MANIFEST_SHA256_EXPECTED = "e6d60aa56d2ac87a9b82aaca8b3df5b560131943eb0ffb809f0834e8a9db3543"
P13569_REL = Path("tools/cftr-model/data/p13569.json")
P13569_SHA256_EXPECTED = "c8feabfabe8690ef4a06af7de5b5f8dc3c8da554246c5ab3d403edf265d4a893"
ATTRIBUTION_REL = Path("tools/cftr-model/data/ATTRIBUTION.md")
EXPR_RESULTS_REL = Path("benchmarks/cftr_dms_expression/results_2026-11-23.json")
CRYOEM_3D_RESULTS_REL = Path("benchmarks/cftr_dms_3d_structure/results_2026-11-24.json")
THERATYPE_FIXTURE_REL = Path("benchmarks/theratype/theratype_open_annotations_2026-11-21.json")
RESULTS_FILE_NAME = "results_2026-11-24.json"

FEATURE_NAMES: Tuple[str, ...] = (
    "delta_blosum62",
    "delta_kd_hydropathy",
    "abs_delta_charge",
    "helix_breaker_in_helix",
    "tm_polarity_shock",
    "is_nonsense",
    "is_inframe_indel",
)

SUBDOMAIN_ORDER: Tuple[str, ...] = ("Lasso", "TMD1", "NBD1", "R Domain", "TMD2", "NBD2")
ALPHA_GRID: Tuple[float, ...] = (0.01, 0.1, 1.0, 10.0, 100.0, 1000.0)
CV_FOLDS = 5
BOOTSTRAP_SEED = 20261124
BOOTSTRAP_RESAMPLES = 10000
SHUFFLE_SEED = 20261124
GATE_REL_MAE_REDUCTION = 0.10
MIN_HOLDOUT_CLUSTERS = 40
PARITY_TOL = 1e-6
VIF_FLAG_THRESHOLD = 10.0

TM_POLARITY_KD_THRESHOLD = 3.0
TM_POLARITY_CHARGE_THRESHOLD = 1.0
HELIX_BREAKER_AAS: FrozenSet[str] = frozenset({"P", "G"})
P13569_EXPECTED_LENGTH = 1480

WT_PARITY_ALLOWLIST: Dict[Tuple[int, int, str], str] = {
    (89, 470, "M470V"): (
        "Common CFTR rs213950 polymorphism: Tedman background M470 vs canonical "
        "P13569 V470; substitution deltas computed from experimental WT M -> V."
    ),
    (217, 966, "Q966="): (
        "Source-table discrepancy: Q966= vs THR966 in P13569 "
        "(2988G>A corresponds to Q996=); kept in dev as synonymous (all 7 features = 0.0)."
    ),
}

# Published Henikoff & Henikoff (1992, PNAS 89:10915) 20x20 BLOSUM62 integer matrix.
BLOSUM62_ORDER = "ARNDCQEGHILKMFPSTWYV"
_BLOSUM62_ROWS: Dict[str, Tuple[int, ...]] = {
    "A": ( 4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0),
    "R": (-1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3),
    "N": (-2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3),
    "D": (-2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3),
    "C": ( 0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1),
    "Q": (-1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2),
    "E": (-1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2),
    "G": ( 0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3),
    "H": (-2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3),
    "I": (-1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3),
    "L": (-1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1),
    "K": (-1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2),
    "M": (-1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1),
    "F": (-2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1),
    "P": (-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2),
    "S": ( 1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2),
    "T": ( 0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0),
    "W": (-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3),
    "Y": (-2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1),
    "V": ( 0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4),
}
BLOSUM62: Dict[str, Dict[str, int]] = {
    aa_i: {aa_j: _BLOSUM62_ROWS[aa_i][idx_j] for idx_j, aa_j in enumerate(BLOSUM62_ORDER)}
    for aa_i in BLOSUM62_ORDER
}

# Canonical signed Kyte & Doolittle (1982, J. Mol. Biol. 157:105) hydropathy index.
KYTE_DOOLITTLE: Dict[str, float] = {
    "I": 4.5, "V": 4.2, "L": 3.8, "F": 2.8, "C": 2.5,
    "M": 1.9, "A": 1.8, "G": -0.4, "T": -0.7, "S": -0.8,
    "W": -0.9, "Y": -1.3, "P": -1.6, "H": -3.2, "Q": -3.5,
    "E": -3.5, "D": -3.5, "N": -3.5, "K": -3.9, "R": -4.5,
}

# Formal side-chain charge class (`K,R=+1`, `D,E=-1`, others `0`).
RESIDUE_CHARGE: Dict[str, float] = {
    aa: (1.0 if aa in {"K", "R"} else (-1.0 if aa in {"D", "E"} else 0.0))
    for aa in BLOSUM62_ORDER
}

STANDARD_AAS: FrozenSet[str] = frozenset(BLOSUM62_ORDER)

XLSX_NS = {"s": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
CELL_REF_RE = re.compile(r"^([A-Z]+)(\d+)$")
RANGE_REF_RE = re.compile(r"^([A-Z]+)(\d+):([A-Z]+)(\d+)$")


@dataclass(frozen=True)
class VariantComponent:
    """Single parsed variant component from `Protein distinction`."""

    kind: str  # "missense" | "synonymous" | "nonsense" | "inframe_indel"
    wt_aa: str
    codon: int
    mut_aa: Optional[str]
    suffix: str


@dataclass(frozen=True)
class VariantSubstRecord:
    """Single DMS row paired with 7 substitution biophysics & P13569 topology features."""

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
    """Parse sparse `.xlsx` worksheet XML by `(row, col)` coordinates."""
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


def parse_variant_components(raw: str) -> Tuple[VariantComponent, ...]:
    """Parse single or `/`-separated compound `Protein distinction` entries into components."""
    s = str(raw or "").strip()
    if not s:
        raise ValueError("Empty Protein distinction value.")
    parts = [part.strip() for part in s.split("/") if part.strip()]
    if not parts:
        raise ValueError(f"Unparseable Protein distinction: {raw!r}")

    components: List[VariantComponent] = []
    for part in parts:
        m_indel = re.fullmatch(r"([A-Z])(\d+)(del|DEL|ins|INS)", part)
        if m_indel:
            wt_aa = m_indel.group(1)
            if wt_aa not in STANDARD_AAS:
                raise ValueError(f"Unsupported amino acid {wt_aa!r} in {raw!r}")
            components.append(
                VariantComponent(
                    kind="inframe_indel",
                    wt_aa=wt_aa,
                    codon=int(m_indel.group(2)),
                    mut_aa=None,
                    suffix=m_indel.group(3),
                )
            )
            continue

        m_stop = re.fullmatch(r"([A-Z])(\d+)([X\*])", part)
        if m_stop:
            wt_aa = m_stop.group(1)
            if wt_aa not in STANDARD_AAS:
                raise ValueError(f"Unsupported amino acid {wt_aa!r} in {raw!r}")
            components.append(
                VariantComponent(
                    kind="nonsense",
                    wt_aa=wt_aa,
                    codon=int(m_stop.group(2)),
                    mut_aa=None,
                    suffix=m_stop.group(3),
                )
            )
            continue

        m_sub = re.fullmatch(r"([A-Z])(\d+)([A-Z=])", part)
        if m_sub:
            wt_aa, codon_s, alt = m_sub.group(1), m_sub.group(2), m_sub.group(3)
            if wt_aa not in STANDARD_AAS:
                raise ValueError(f"Unsupported WT amino acid {wt_aa!r} in {raw!r}")
            if alt == "=" or alt == wt_aa:
                components.append(
                    VariantComponent(
                        kind="synonymous",
                        wt_aa=wt_aa,
                        codon=int(codon_s),
                        mut_aa=wt_aa,
                        suffix=alt,
                    )
                )
            else:
                if alt not in STANDARD_AAS:
                    raise ValueError(f"Unsupported mutant amino acid {alt!r} in {raw!r}")
                components.append(
                    VariantComponent(
                        kind="missense",
                        wt_aa=wt_aa,
                        codon=int(codon_s),
                        mut_aa=alt,
                        suffix=alt,
                    )
                )
            continue

        raise ValueError(f"Unparseable variant component {part!r} in {raw!r}")

    return tuple(components)


def parse_protein_distinction(raw: str) -> Tuple[Optional[int], str, Optional[str]]:
    """Return `(codon_position, variant_class, wt_one_letter)` matching `split_manifest.json`."""
    s = str(raw or "").strip()
    if not s:
        return None, "other", None
    try:
        comps = parse_variant_components(s)
    except ValueError:
        return None, "other", None
    first = comps[0]
    if len(comps) > 1 or first.kind == "synonymous":
        return first.codon, "other", first.wt_aa
    return first.codon, first.kind, first.wt_aa


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


def load_p13569_context(repo_root: Path = REPO_ROOT) -> Dict[str, Any]:
    """Load tracked UniProt P13569 JSON, verify SHA-256 and ATTRIBUTION.md, and parse topology."""
    p13569_path = repo_root / P13569_REL
    raw_bytes = p13569_path.read_bytes()
    digest = hashlib.sha256(raw_bytes).hexdigest()
    if digest != P13569_SHA256_EXPECTED:
        raise AssertionError(f"p13569.json SHA256 mismatch: {digest} != {P13569_SHA256_EXPECTED}")

    attr_path = repo_root / ATTRIBUTION_REL
    attr_txt = attr_path.read_text(encoding="utf-8")
    if "p13569.json" not in attr_txt or P13569_SHA256_EXPECTED not in attr_txt:
        raise AssertionError("ATTRIBUTION.md missing expected p13569.json SHA-256 entry.")

    payload = json.loads(raw_bytes.decode("utf-8"))
    if payload.get("primaryAccession") != "P13569":
        raise AssertionError(f"Unexpected UniProt primaryAccession: {payload.get('primaryAccession')!r}")
    seq = str(payload.get("sequence", {}).get("value", ""))
    if len(seq) != P13569_EXPECTED_LENGTH:
        raise AssertionError(f"Unexpected P13569 sequence length: {len(seq)} != {P13569_EXPECTED_LENGTH}")

    intervals: Dict[str, List[Tuple[int, int]]] = {
        "Transmembrane": [],
        "Helix": [],
        "Beta strand": [],
        "Binding site": [],
    }
    pos_sets: Dict[str, set[int]] = {k: set() for k in intervals}

    for feat in payload.get("features", []):
        ftype = feat.get("type")
        if ftype not in intervals:
            continue
        loc = feat.get("location", {})
        start = int(loc["start"]["value"])
        end = int(loc["end"]["value"])
        if start < 1 or end > len(seq) or start > end:
            raise AssertionError(f"Invalid P13569 {ftype} interval [{start}, {end}]")
        intervals[ftype].append((start, end))
        pos_sets[ftype].update(range(start, end + 1))

    if not intervals["Transmembrane"] or not intervals["Helix"]:
        raise AssertionError("P13569 JSON missing Transmembrane or Helix annotations.")

    return {
        "sha256": digest,
        "sequence": seq,
        "sequence_length": len(seq),
        "intervals": intervals,
        "tm_positions": frozenset(pos_sets["Transmembrane"]),
        "helix_positions": frozenset(pos_sets["Helix"]),
        "helix_or_tm_positions": frozenset(pos_sets["Helix"] | pos_sets["Transmembrane"]),
        "beta_positions": frozenset(pos_sets["Beta strand"]),
        "binding_site_positions": frozenset(pos_sets["Binding site"]),
        "feature_counts": {k: len(v) for k, v in intervals.items()},
    }


def verify_p13569_wt_parity(
    assignments: Sequence[Mapping[str, Any]],
    p13569_ctx: Mapping[str, Any],
) -> Dict[str, Any]:
    """Fail closed if any parsed component's WT residue mismatches P13569 outside allowlisted dev rows."""
    seq = str(p13569_ctx["sequence"])
    unexpected: List[Dict[str, Any]] = []
    documented_observed: List[Dict[str, Any]] = []
    seen_allowlist_keys: set[Tuple[int, int, str]] = set()
    components_checked = 0

    for row in assignments:
        row_num = int(row["row_num"])
        prot = str(row["protein_distinction"])
        split = str(row["split"])
        try:
            comps = parse_variant_components(prot)
        except ValueError as exc:
            unexpected.append({"row_num": row_num, "prot": prot, "reason": str(exc)})
            continue

        for comp in comps:
            components_checked += 1
            if comp.codon < 1 or comp.codon > len(seq):
                unexpected.append(
                    {"row_num": row_num, "codon_position": comp.codon, "prot": prot, "reason": "codon_out_of_range"}
                )
                continue
            obs_aa = seq[comp.codon - 1]
            if obs_aa != comp.wt_aa:
                key = (row_num, comp.codon, prot)
                if key in WT_PARITY_ALLOWLIST and split == "dev":
                    seen_allowlist_keys.add(key)
                    documented_observed.append({
                        "row_num": row_num,
                        "codon_position": comp.codon,
                        "protein_distinction": prot,
                        "split": split,
                        "expected_experimental_wt_aa": comp.wt_aa,
                        "observed_p13569_aa": obs_aa,
                        "ruling": WT_PARITY_ALLOWLIST[key],
                    })
                else:
                    unexpected.append({
                        "row_num": row_num,
                        "codon_position": comp.codon,
                        "protein_distinction": prot,
                        "split": split,
                        "expected_experimental_wt_aa": comp.wt_aa,
                        "observed_p13569_aa": obs_aa,
                    })

    if unexpected:
        raise AssertionError(f"P13569 WT sequence parity failed closed: {unexpected}")
    if seen_allowlist_keys != set(WT_PARITY_ALLOWLIST):
        raise AssertionError(
            f"Expected documented WT parity exceptions {sorted(WT_PARITY_ALLOWLIST)}, saw {sorted(seen_allowlist_keys)}"
        )

    return {
        "passed": True,
        "variants_checked": len(assignments),
        "components_checked": components_checked,
        "unexpected_mismatches": 0,
        "holdout_mismatches": 0,
        "documented_dev_exceptions": documented_observed,
    }


def extract_7_biophysics_features(
    prot: str,
    p13569_ctx: Mapping[str, Any],
) -> Tuple[Optional[float], ...]:
    """Compute the fixed 7-feature vector from experimental-WT biophysics and P13569 topology."""
    comps = parse_variant_components(prot)

    if all(c.kind == "synonymous" for c in comps):
        return (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

    if len(comps) == 1 and comps[0].kind == "nonsense":
        return (None, None, None, 0.0, 0.0, 1.0, 0.0)

    if len(comps) == 1 and comps[0].kind == "inframe_indel":
        return (None, None, None, 0.0, 0.0, 0.0, 1.0)

    if not all(c.kind == "missense" for c in comps):
        raise ValueError(f"Unsupported mixed variant components in {prot!r}")

    helix_or_tm: FrozenSet[int] = p13569_ctx["helix_or_tm_positions"]
    tm_pos: FrozenSet[int] = p13569_ctx["tm_positions"]

    delta_blosum = 0.0
    delta_kd = 0.0
    abs_dcharge = 0.0
    helix_breaker = 0.0
    tm_shock = 0.0

    for comp in comps:
        assert comp.mut_aa is not None
        wt = comp.wt_aa
        mut = comp.mut_aa
        pos = comp.codon

        d_blo = float(BLOSUM62[wt][mut] - BLOSUM62[wt][wt])
        d_kd = float(KYTE_DOOLITTLE[mut] - KYTE_DOOLITTLE[wt])
        d_chg = abs(RESIDUE_CHARGE[mut] - RESIDUE_CHARGE[wt])

        delta_blosum += d_blo
        delta_kd += d_kd
        abs_dcharge += d_chg

        if pos in helix_or_tm and (wt in HELIX_BREAKER_AAS or mut in HELIX_BREAKER_AAS):
            helix_breaker = 1.0
        if pos in tm_pos and (abs(d_kd) >= TM_POLARITY_KD_THRESHOLD or d_chg >= TM_POLARITY_CHARGE_THRESHOLD):
            tm_shock = 1.0

    return (delta_blosum, delta_kd, abs_dcharge, helix_breaker, tm_shock, 0.0, 0.0)


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


def load_tedman_with_biophysics_features(
    master_xlsx: Path,
    manifest: Mapping[str, Any],
    p13569_ctx: Mapping[str, Any],
) -> List[VariantSubstRecord]:
    """Parse `All_variants`, verify 234/234 parity with `split_manifest.json`, and extract 7 features."""
    master_sha = hashlib.sha256(master_xlsx.read_bytes()).hexdigest()
    if master_sha != MASTER_XLSX_SHA256:
        raise ValueError(f"Master xlsx SHA256 mismatch: {master_sha} != {MASTER_XLSX_SHA256}")

    cells, max_row = parse_xlsx_sheet(master_xlsx, "All_variants", expand_merged=True)
    by_row = {int(a["row_num"]): a for a in manifest["assignments"]}
    records: List[VariantSubstRecord] = []
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

        feats = extract_7_biophysics_features(prot, p13569_ctx)
        records.append(
            VariantSubstRecord(
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
class FittedPartialSubstRidgeModel:
    """Partial-penalized Ridge on `x_dmso_log2` + 7 fold-imputed/standardized biophysics features."""

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


def fit_partial_penalized_subst_ridge(
    xs_dmso: Sequence[float],
    raw_features_rows: Sequence[Sequence[Optional[float]]],
    ys: Sequence[float],
    alpha: float,
) -> FittedPartialSubstRidgeModel:
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
    return FittedPartialSubstRidgeModel(
        alpha=alpha,
        b0_intercept=beta[0],
        b1_dmso_slope=beta[1],
        impute_means=tuple(impute_means),
        feat_means=tuple(means),
        feat_stds=tuple(stds),
        weights=tuple(beta[2:]),
    )


def _dev_fold_map(dev_rows: Sequence[VariantSubstRecord]) -> Dict[int, int]:
    sorted_codons = sorted({r.codon_position for r in dev_rows})
    return {c: idx % CV_FOLDS for idx, c in enumerate(sorted_codons)}


def _cv_select_alpha(
    dev_rows: Sequence[VariantSubstRecord],
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
            mod = fit_partial_penalized_subst_ridge(
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


def compute_vif_context(dev_rows: Sequence[VariantSubstRecord]) -> Dict[str, Any]:
    """Compute non-gated `dev` variance inflation factors `(R^-1)_jj` and flag any `VIF > 10`."""
    raw = [r.features_raw for r in dev_rows]
    mod = fit_partial_penalized_subst_ridge(
        [r.x_dmso_log2 for r in dev_rows], raw, [r.y_combo_log2 for r in dev_rows], alpha=1.0
    )
    z_7 = [mod.transform_features(row) for row in raw]
    n = len(z_7)

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
    z_all = [[(xs[i] - mx) / sx] + list(z_7[i]) for i in range(n)]

    vif_7 = _corr_vifs(z_7, FEATURE_NAMES)
    vif_with_dmso = _corr_vifs(z_all, ("x_dmso_log2",) + FEATURE_NAMES)
    flagged_7 = sorted(k for k, v in vif_7.items() if v > VIF_FLAG_THRESHOLD)
    flagged_all = sorted(k for k, v in vif_with_dmso.items() if v > VIF_FLAG_THRESHOLD)
    return {
        "vif_threshold": VIF_FLAG_THRESHOLD,
        "vif_features_only": vif_7,
        "vif_with_x_dmso_log2": vif_with_dmso,
        "max_vif_features_only": round(max(vif_7.values()), 6),
        "max_vif_with_x_dmso_log2": round(max(vif_with_dmso.values()), 6),
        "flagged_gt_10_features_only": flagged_7,
        "flagged_gt_10_with_x_dmso": flagged_all,
        "any_vif_gt_10": bool(flagged_7 or flagged_all),
    }


def residue_clustered_bootstrap_delta_ci(
    holdout_rows: Sequence[VariantSubstRecord],
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


def _mae_log2_and_linear(rows: Sequence[VariantSubstRecord], preds_log2: Sequence[float]) -> Tuple[float, float]:
    mae_log2 = sum(abs(r.y_combo_log2 - p) for r, p in zip(rows, preds_log2)) / len(rows)
    mae_lin = sum(abs(r.combo_rel_mut_wt - (2.0 ** p)) for r, p in zip(rows, preds_log2)) / len(rows)
    return mae_log2, mae_lin


def evaluate_substitution_biophysics_benchmark(
    data_dir: Path, repo_root: Path = REPO_ROOT
) -> Dict[str, Any]:
    """Run all pre-holdout controls (fail closed), then single-shot holdout evaluation."""
    manifest, manifest_sha = load_locked_split_manifest(repo_root)
    p13569_ctx = load_p13569_context(repo_root)
    wt_parity = verify_p13569_wt_parity(manifest["assignments"], p13569_ctx)

    records = load_tedman_with_biophysics_features(data_dir / MASTER_XLSX_NAME, manifest, p13569_ctx)
    dev_rows = [r for r in records if r.split == "dev"]
    holdout_rows = [r for r in records if r.split == "holdout"]
    dev_codons = sorted({r.codon_position for r in dev_rows})
    holdout_codons = sorted({r.codon_position for r in holdout_rows})

    # Reference #396 & #401 JSON payloads for parity & non-gated context comparison
    expr_ref = json.loads((repo_root / EXPR_RESULTS_REL).read_text(encoding="utf-8"))
    cryoem_ref = json.loads((repo_root / CRYOEM_3D_RESULTS_REL).read_text(encoding="utf-8"))

    # Pre-holdout Control 2a: `dev` baseline fitting & 5-fold CV parity with #396
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

    # Pre-holdout Control 3: Null biophysics control (dev-only variant-row permutation, seed=20261124)
    perm_idx = list(range(len(dev_rows)))
    random.Random(SHUFFLE_SEED).shuffle(perm_idx)
    permuted_dev_feats = [dev_rows[idx].features_raw for idx in perm_idx]
    shuf_alpha, shuf_cv_mae, _ = _cv_select_alpha(dev_rows, permuted_dev_feats, fold_map)
    shuf_rel_red = (locked_dev_cv_mae - shuf_cv_mae) / locked_dev_cv_mae
    if shuf_rel_red >= GATE_REL_MAE_REDUCTION:
        raise AssertionError(
            f"Dev variant-permuted null biophysics control exceeded 10% threshold: {shuf_rel_red:.6f}"
        )

    # Dev CV selection for true substitution biophysics candidate (zero holdout contact)
    best_alpha, best_cv_mae, cv_path = _cv_select_alpha(
        dev_rows, [r.features_raw for r in dev_rows], fold_map
    )
    final_cand = fit_partial_penalized_subst_ridge(
        [r.x_dmso_log2 for r in dev_rows],
        [r.features_raw for r in dev_rows],
        [r.y_combo_log2 for r in dev_rows],
        best_alpha,
    )
    vifs = compute_vif_context(dev_rows)
    weights_by_name = {
        name: round(w, 6) for name, w in zip(FEATURE_NAMES, final_cand.weights)
    }

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
        gate_verdict = "PASS_SUBST_BIOPHYSICS_BEATS_DMSO_OLS"
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

    ref_396_rel = float(ref_hold["relative_mae_reduction_vs_stronger"])
    ref_401_rel = float(cryoem_ref["holdout_evaluation"]["relative_mae_reduction_vs_stronger"])

    done_sha = hashlib.sha256(
        (repo_root / "benchmarks" / "cftr_dms_substitution_biophysics" / "done.md")
        .read_text(encoding="utf-8", errors="replace")
        .replace("\r\n", "\n")
        .encode("utf-8")
    ).hexdigest()

    return {
        "schema_version": "tedman2025-substitution-biophysics-results-v1",
        "benchmark_id": "cftr-dms-substitution-biophysics",
        "evaluated_on": dt.date.today().isoformat(),
        "evaluated_on_local": dt.date.today().isoformat(),
        "evaluated_on_utc": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%d"),
        "done_sha256": done_sha,
        "split_manifest_sha256": manifest_sha,
        "p13569_sha256": p13569_ctx["sha256"],
        "dataset": {
            "source_id": "cftr-dms-tedman-2025",
            "doi": "10.17632/ksfmxv8zts.1",
            "paper_doi": "10.7554/eLife.107180",
            "license": "CC BY 4.0",
            "uniprot_accession": "P13569",
            "uniprot_license": "CC BY 4.0",
            "uniprot_sequence_length": p13569_ctx["sequence_length"],
            "uniprot_feature_counts": p13569_ctx["feature_counts"],
            "assay": "HEK293T parental CFTR DMS plasma-membrane intensity (VX-661+VX-445 vs DMSO, Mut/WT)",
            "scale_boundary": "Plasma-membrane expression in HEK293T cells is not CFTR channel conductance or clinical lung benefit.",
            "total_valid_variants": len(records),
            "excluded_count": 0,
            "dev_variants": len(dev_rows),
            "holdout_variants": len(holdout_rows),
            "dev_unique_codons": len(dev_codons),
            "holdout_unique_codons": len(holdout_codons),
            "holdout_pred0010_overlap_codons": 0,
            "dev_missing_continuous_scale_rows": sum(1 for r in dev_rows if r.features_raw[0] is None),
            "holdout_missing_continuous_scale_rows": sum(1 for r in holdout_rows if r.features_raw[0] is None),
        },
        "pre_holdout_controls": {
            "p13569_wt_sequence_parity": wt_parity,
            "baseline_396_dev_parity_passed": True,
            "dev_variant_permuted_null_control": {
                "seed": SHUFFLE_SEED,
                "permutation_unit": "dev_variants",
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
            "candidate_dmso_subst_biophysics_ridge": {
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
                "weights_by_feature": weights_by_name,
                "is_inframe_indel_dev_weight": weights_by_name["is_inframe_indel"],
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
            "candidate_dmso_subst_biophysics_ridge": {
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
            "non_gated_context_comparison": {
                "pr_396_domain_class_ridge_rel_reduction": ref_396_rel,
                "pr_396_domain_class_ridge_pct": round(ref_396_rel * 100.0, 4),
                "pr_401_3d_cryoem_ridge_rel_reduction": ref_401_rel,
                "pr_401_3d_cryoem_ridge_pct": round(ref_401_rel * 100.0, 4),
                "this_subst_biophysics_rel_reduction": round(rel_improvement, 6),
                "this_subst_biophysics_pct": round(rel_improvement * 100.0, 4),
            },
            "non_gated_subdomain_breakdown": subdomain_breakdown,
        },
    }


def verify_checked_in_artifacts(repo_root: Path = REPO_ROOT, data_dir: Optional[Path] = None) -> Dict[str, Any]:
    """Verify offline P13569, BLOSUM62/KD, and checked-in results (plus raw `.xlsx` parity if present)."""
    bench_dir = repo_root / "benchmarks" / "cftr_dms_substitution_biophysics"
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

    p13569_ctx = load_p13569_context(repo_root)
    if results["p13569_sha256"] != p13569_ctx["sha256"]:
        raise AssertionError("results JSON p13569_sha256 mismatch.")

    wt_parity = verify_p13569_wt_parity(manifest["assignments"], p13569_ctx)
    if not wt_parity["passed"]:
        raise AssertionError("P13569 WT sequence parity check failed.")

    # Verify all 234 manifest rows parse cleanly into 7-feature vectors even in CI without .xlsx
    for row in manifest["assignments"]:
        feats = extract_7_biophysics_features(str(row["protein_distinction"]), p13569_ctx)
        if len(feats) != len(FEATURE_NAMES):
            raise AssertionError(f"Expected {len(FEATURE_NAMES)} features, got {len(feats)}")

    resolved_data = resolve_tedman_data_dir(data_dir, repo_root=repo_root)
    raw_parity_checked = False
    if resolved_data is not None:
        fresh = evaluate_substitution_biophysics_benchmark(resolved_data, repo_root=repo_root)
        f_he = fresh["holdout_evaluation"]
        r_he = results["holdout_evaluation"]
        if (
            abs(
                f_he["candidate_dmso_subst_biophysics_ridge"]["mae_log2"]
                - r_he["candidate_dmso_subst_biophysics_ridge"]["mae_log2"]
            )
            > PARITY_TOL
            or abs(f_he["dmso_ols"]["mae_log2"] - r_he["dmso_ols"]["mae_log2"]) > PARITY_TOL
            or f_he["gate_verdict"] != r_he["gate_verdict"]
        ):
            raise AssertionError("Live raw .xlsx evaluation differs from checked-in results_2026-11-24.json.")
        raw_parity_checked = True

    return {
        "status": "PASS",
        "raw_local_data_parity_checked": raw_parity_checked,
        "p13569_sequence_length": p13569_ctx["sequence_length"],
        "wt_components_checked": wt_parity["components_checked"],
        "dev_variants": manifest["dev_variants"],
        "holdout_variants": manifest["holdout_variants"],
        "holdout_unique_codons": manifest["holdout_unique_codons"],
        "gate_verdict": results["holdout_evaluation"]["gate_verdict"],
    }


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--data-dir", type=Path, default=None, help="Optional path to local tedman-2025 folder.")
    parser.add_argument("--write-artifacts", action="store_true", help="Run evaluation and write results_2026-11-24.json.")
    parser.add_argument("--check", action="store_true", help="Verify checked-in P13569 invariants and local raw parity.")
    parser.add_argument("--json", action="store_true", help="Print strict JSON (allow_nan=False).")
    args = parser.parse_args(argv)

    if args.check:
        summary = verify_checked_in_artifacts(repo_root=REPO_ROOT, data_dir=args.data_dir)
        if args.json:
            print(json.dumps(summary, indent=2, allow_nan=False))
        else:
            print(
                f"PASS cftr_dms_substitution_biophysics: raw_parity={summary['raw_local_data_parity_checked']} "
                f"p13569_len={summary['p13569_sequence_length']} "
                f"verdict={summary['gate_verdict']}"
            )
        return 0

    resolved_data = resolve_tedman_data_dir(args.data_dir, repo_root=REPO_ROOT)
    if resolved_data is None:
        print("Local Tedman 2025 workbook not present under sources/downloads/tedman-2025/.", file=sys.stderr)
        return 1

    results = evaluate_substitution_biophysics_benchmark(resolved_data, repo_root=REPO_ROOT)
    if args.write_artifacts:
        (HERE / RESULTS_FILE_NAME).write_text(
            json.dumps(results, indent=2, allow_nan=False) + "\n", encoding="utf-8"
        )

    if args.json:
        print(json.dumps(results, indent=2, allow_nan=False))
    else:
        he = results["holdout_evaluation"]
        print(
            f"Tedman 2025 Substitution Biophysics: verdict={he['gate_verdict']} | "
            f"stronger({he['stronger_baseline_on_holdout']})={he['stronger_baseline_mae_log2']} | "
            f"candidate={he['candidate_dmso_subst_biophysics_ridge']['mae_log2']} | "
            f"rel_reduction={he['relative_mae_reduction_vs_stronger']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
