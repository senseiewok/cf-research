#!/usr/bin/env python3
"""Modulator-Free vs Trikafta-Bound deltaF508 Cryo-EM Rescue Benchmark (PRED-0014).

100% Pure Python 3.10+ Standard Library Runtime. Zero network calls. Zero patient/PHI data.
Evaluates whether SE(3)xSE(3)-invariant, swap-symmetric, zero-free-parameter local backbone
substructure conservation s_stab(i) = lDDT_sym(i; 8EJ1, 8EIQ) (Mariani et al. 2013 published
defaults R_0 = 15.0 A, thresholds {0.5, 1.0, 2.0, 4.0} A) between dephosphorylated modulator-free
deltaF508 human CFTR (PDB 8EJ1, 6.9 A) and phosphorylated Trikafta+ATP/Mg-bound deltaF508 human
CFTR (PDB 8EIQ, 3.0 A; Fiedorczuk & Chen Science 2022 DOI 10.1126/science.ade2216) predicts
US FDA TRIKAFTA v23 Table 6 Structured Product Label (SPL) listing concordance across the
AST-blinded 24-variant untouched holdout panel (holdout_open_annotations_2026-11-23.json).

Research benchmark only; not clinical or diagnostic guidance.
"""
from __future__ import annotations

import argparse
import copy
from dataclasses import dataclass
import gzip
import hashlib
import json
import math
from pathlib import Path
import re
import shlex
import sys
from typing import Any, Mapping, Sequence

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]

DEFAULT_CIF_8EJ1 = REPO_ROOT / "tools" / "cftr-model" / "data" / "8EJ1.cif.gz"
DEFAULT_CIF_8EIQ = REPO_ROOT / "tools" / "cftr-model" / "data" / "8EIQ.cif.gz"
DEFAULT_UNIPROT_JSON = REPO_ROOT / "tools" / "cftr-model" / "data" / "p13569.json"
DEFAULT_DEV_JSON = REPO_ROOT / "benchmarks" / "theratype" / "theratype_open_annotations_2026-11-21.json"
DEFAULT_HOLDOUT_JSON = HERE / "holdout_open_annotations_2026-11-23.json"
DEFAULT_CATALOG_YAML = REPO_ROOT / "sources" / "catalog.yaml"

RESEARCH_BOUNDARY = "Research benchmark only; not clinical or diagnostic guidance. Zero patient data."

EXPECTED_SHA256 = {
    "8EJ1": "9582203b74c0add00a7917bed765cbc8a943fcc3a88ffbac8ae8dde7de2bc159",
    "8EIQ": "d13cbfc8ba6575ab6228f351578c6f9cc5b9dafe13e5913d9b27b3cacebc18d6",
    "P13569": "c8feabfabe8690ef4a06af7de5b5f8dc3c8da554246c5ab3d403edf265d4a893",
    "HOLDOUT_VARIANTS": "8d410145e5219bfa405dbf7a24f68cedb7a5f88c5035b0da6936b5d402024bb3",
    "PRE_REGISTERED_HOLDOUT_VARIANTS": "9b815b5bc581363fc8cc7c7c4b5c31fc346d60ee2047256bc264d544062eb703",
}

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

ALLOWED_FAMILIES = frozenset({"missense", "inframe_indel", "splice_regulatory", "class_i_ptc_frameshift"})
ALLOWED_SPL_STATES = frozenset({"listed_on_fda_spl", "not_listed_on_fda_spl"})
ALLOWED_EXPOSURE_STATES = frozenset({"untouched_holdout", "exposed_development_reference"})
MANUAL_UNADMITTED_SOURCES = frozenset({"mavedb-dump-zenodo", "cftr-dms-mckee-2023", "cftr-dms-tedman-2025"})
CONSTRUCT_MODIFIED_CODONS = frozenset({508, 1371})
RESERVED_CONTROL_CODONS = frozenset(
    {38, 195, 267, 311, 361, 450, 464, 507, 508, 509, 582, 1030, 1102, 1246, 1293, 1347, 1348, 1371}
)
SCAFFOLD_CONTROL_CODONS = (38, 311, 361, 450, 1030, 1102)
ATP_SITE_CONTROL_CODONS = (464, 582, 1246, 1347, 1348)

R0_INCLUSION_ANGSTROM = 15.0
LDDT_THRESHOLDS_ANGSTROM = (0.5, 1.0, 2.0, 4.0)
MIN_CUTOFF_MARGIN_ANGSTROM = 1e-9
MIN_CODON_BUFFER = 5

CANONICAL_SPLICE_RE = re.compile(r"^c\.\d+[+-][12][ACGT]>[ACGT]$")
INDEL_SPAN_RE = re.compile(r"^p\.[A-Z][a-z]{2}(\d+)_[A-Z][a-z]{2}(\d+)del$")


class CryoEMF508delRescueError(ValueError):
    """Raised when construct provenance, symmetry, split hygiene, or falsification checks fail."""


@dataclass(frozen=True)
class DeltaF508ChainAStructure:
    pdb_id: str
    ca_by_seq: dict[int, tuple[str, tuple[float, float, float]]]
    deleted_codons: frozenset[int]
    engineered_mutations: dict[int, str]


def verify_sha256(path: Path, expected_hex: str, label: str) -> str:
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != expected_hex:
        raise CryoEMF508delRescueError(f"{label} SHA-256 mismatch: expected {expected_hex}, got {digest}")
    return digest


def load_uniprot_p13569_sequence(uniprot_path: Path = DEFAULT_UNIPROT_JSON) -> str:
    verify_sha256(uniprot_path, EXPECTED_SHA256["P13569"], "P13569")
    payload = json.loads(uniprot_path.read_text(encoding="utf-8"))
    seq = str(payload.get("sequence", {}).get("value", ""))
    if len(seq) != 1480:
        raise CryoEMF508delRescueError(f"Expected 1480-aa P13569 sequence, got {len(seq)}")
    return seq


def is_covalent_deltaf508_pair(seq_i: int, seq_j: int, exclude_507_509: bool = True) -> bool:
    if abs(seq_i - seq_j) <= 1:
        return True
    if exclude_507_509 and ((seq_i == 507 and seq_j == 509) or (seq_i == 509 and seq_j == 507)):
        return True
    return False


def euclidean_distance(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return math.hypot(a[0] - b[0], a[1] - b[1], a[2] - b[2])


def verify_proper_rotation_matrix(rot: Sequence[Sequence[float]], tol: float = 1e-12) -> None:
    if len(rot) != 3 or any(len(r) != 3 for r in rot):
        raise CryoEMF508delRescueError("Rotation matrix must be 3x3")
    for r in range(3):
        for c in range(3):
            dot = sum(rot[k][r] * rot[k][c] for k in range(3))
            target = 1.0 if r == c else 0.0
            if abs(dot - target) > tol:
                raise CryoEMF508delRescueError(f"Non-orthogonal matrix entry ({r},{c}): {dot}")
    det = (
        rot[0][0] * (rot[1][1] * rot[2][2] - rot[1][2] * rot[2][1])
        - rot[0][1] * (rot[1][0] * rot[2][2] - rot[1][2] * rot[2][0])
        + rot[0][2] * (rot[1][0] * rot[2][1] - rot[1][1] * rot[2][0])
    )
    if abs(det - 1.0) > tol:
        raise CryoEMF508delRescueError(f"Improper rotation det={det:.12f} (expected +1.0 for SO(3))")


def apply_rigid_motion(
    coords: Mapping[int, tuple[float, float, float]],
    rot: Sequence[Sequence[float]],
    trans: tuple[float, float, float],
) -> dict[int, tuple[float, float, float]]:
    verify_proper_rotation_matrix(rot)
    out: dict[int, tuple[float, float, float]] = {}
    for k, (x, y, z) in coords.items():
        nx = rot[0][0] * x + rot[0][1] * y + rot[0][2] * z + trans[0]
        ny = rot[1][0] * x + rot[1][1] * y + rot[1][2] * z + trans[1]
        nz = rot[2][0] * x + rot[2][1] * y + rot[2][2] * z + trans[2]
        out[k] = (nx, ny, nz)
    return out


def parse_deltaf508_chain_a_mmcif(
    cif_text: str,
    pdb_id: str,
    uniprot_seq: str,
    *,
    allow_hetatm: bool = False,
    allow_non_chain_a: bool = False,
    use_label_seq_id: bool = False,
) -> DeltaF508ChainAStructure:
    lines = cif_text.splitlines()
    cols: list[str] = []
    mode: str | None = None
    seq_difs: list[dict[str, str]] = []
    ca_by_seq: dict[int, tuple[str, tuple[float, float, float]]] = {}

    for raw in lines:
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
                mode = "atom"
            elif cols[0].startswith("_struct_ref_seq_dif."):
                mode = "dif"
            else:
                mode = None
        if mode in ("atom", "dif"):
            if not s or s.startswith("#") or s.startswith("_") or s == "loop_":
                mode = None
                cols = []
                continue
            toks = shlex.split(s)
            if len(toks) != len(cols):
                continue
            row = {c.split(".", 1)[1]: v for c, v in zip(cols, toks)}
            if mode == "dif":
                seq_difs.append(row)
            elif mode == "atom":
                grp = row.get("group_PDB", "")
                allowed_groups = {"ATOM", "HETATM"} if allow_hetatm else {"ATOM"}
                if grp not in allowed_groups:
                    continue
                if row.get("pdbx_PDB_model_num") != "1":
                    continue
                if not allow_non_chain_a:
                    if row.get("auth_asym_id") != "A":
                        continue
                    if grp == "ATOM" and row.get("label_asym_id") != "A":
                        continue
                if row.get("pdbx_PDB_ins_code") not in {"?", "."}:
                    continue
                if row.get("label_alt_id") not in {".", "A"}:
                    continue
                if grp == "ATOM" and row.get("label_atom_id") != "CA":
                    continue
                seq_field = "label_seq_id" if use_label_seq_id else "auth_seq_id"
                raw_seq = row.get(seq_field, "?")
                if raw_seq in {"?", "."}:
                    continue
                seq_no = int(raw_seq)
                if not (1 <= seq_no <= 1480):
                    if grp == "HETATM":
                        raise CryoEMF508delRescueError(
                            f"{pdb_id}: non-polypeptide HETATM ligand {row.get('label_comp_id')} at seq {seq_no} outside 1..1480"
                        )
                    continue
                if seq_no in ca_by_seq:
                    raise CryoEMF508delRescueError(f"{pdb_id}: duplicate CA atom at residue {seq_no}")
                comp = row["label_comp_id"]
                pt = (float(row["Cartn_x"]), float(row["Cartn_y"]), float(row["Cartn_z"]))
                ca_by_seq[seq_no] = (comp, pt)

    deleted_codons: set[int] = set()
    engineered_mutations: dict[int, str] = {}
    for d in seq_difs:
        if d.get("pdbx_pdb_strand_id") != "A" or d.get("pdbx_seq_db_accession_code") != "P13569":
            continue
        det = str(d.get("details", "")).lower()
        if "expression tag" in det or "cloning" in det:
            continue
        db_num = str(d.get("pdbx_seq_db_seq_num", ""))
        if not db_num.isdigit():
            continue
        pos = int(db_num)
        if not (1 <= pos <= 1480):
            continue
        mon = d.get("mon_id", "?")
        if mon in {"?", "."}:
            deleted_codons.add(pos)
        else:
            engineered_mutations[pos] = mon

    if 508 not in deleted_codons:
        raise CryoEMF508delRescueError(f"{pdb_id}: missing required Phe508 deletion in _struct_ref_seq_dif")
    if engineered_mutations.get(1371) != "GLN":
        raise CryoEMF508delRescueError(f"{pdb_id}: missing required E1371Q GLN declaration in _struct_ref_seq_dif")
    if 508 in ca_by_seq:
        raise CryoEMF508delRescueError(f"{pdb_id}: deleted codon 508 present in Chain A coordinates")

    for pos, (comp, _) in sorted(ca_by_seq.items()):
        if pos in engineered_mutations:
            if comp != engineered_mutations[pos]:
                raise CryoEMF508delRescueError(
                    f"{pdb_id}: residue {pos} {comp} disagrees with declared {engineered_mutations[pos]}"
                )
        else:
            aa = THREE_TO_ONE.get(comp)
            if aa != uniprot_seq[pos - 1]:
                raise CryoEMF508delRescueError(
                    f"{pdb_id}: residue {pos} ({comp}={aa}) disagrees with P13569 ({uniprot_seq[pos - 1]})"
                )

    if 507 in ca_by_seq and 509 in ca_by_seq:
        d_junc = euclidean_distance(ca_by_seq[507][1], ca_by_seq[509][1])
        if not (3.70 <= d_junc <= 3.90):
            raise CryoEMF508delRescueError(
                f"{pdb_id}: deltaF508 (507, 509) peptide bond distance {d_junc:.4f} A outside [3.70, 3.90] A"
            )

    return DeltaF508ChainAStructure(
        pdb_id=pdb_id,
        ca_by_seq=ca_by_seq,
        deleted_codons=frozenset(deleted_codons),
        engineered_mutations=engineered_mutations,
    )


def load_verified_structures(
    cif_8ej1_path: Path = DEFAULT_CIF_8EJ1,
    cif_8eiq_path: Path = DEFAULT_CIF_8EIQ,
    uniprot_path: Path = DEFAULT_UNIPROT_JSON,
) -> tuple[DeltaF508ChainAStructure, DeltaF508ChainAStructure, str]:
    verify_sha256(cif_8ej1_path, EXPECTED_SHA256["8EJ1"], "8EJ1")
    verify_sha256(cif_8eiq_path, EXPECTED_SHA256["8EIQ"], "8EIQ")
    uniprot_seq = load_uniprot_p13569_sequence(uniprot_path)
    with gzip.open(cif_8ej1_path, "rt", encoding="utf-8") as fh_8ej1:
        text_8ej1 = fh_8ej1.read()
    with gzip.open(cif_8eiq_path, "rt", encoding="utf-8") as fh_8eiq:
        text_8eiq = fh_8eiq.read()
    s_8ej1 = parse_deltaf508_chain_a_mmcif(text_8ej1, "8EJ1", uniprot_seq)
    s_8eiq = parse_deltaf508_chain_a_mmcif(text_8eiq, "8EIQ", uniprot_seq)
    return s_8ej1, s_8eiq, uniprot_seq


def compute_lddt_sym_conservation_map(
    coords_a: Mapping[int, tuple[float, float, float]],
    coords_b: Mapping[int, tuple[float, float, float]],
    *,
    r0: float = R0_INCLUSION_ANGSTROM,
    thresholds: Sequence[float] = LDDT_THRESHOLDS_ANGSTROM,
    exclude_507_509: bool = True,
    symmetric_inclusion: bool = True,
    unmodeled_fallback: float | None = None,
) -> tuple[dict[int, float | None], dict[int, int], float]:
    if unmodeled_fallback is not None:
        raise CryoEMF508delRescueError("Numeric fallback on unmodeled or empty-neighborhood residues is forbidden")

    common = sorted(set(coords_a.keys()) & set(coords_b.keys()))
    score_map: dict[int, float | None] = {}
    neighbor_counts: dict[int, int] = {}
    min_margin = float("inf")

    for i in common:
        pt_ai = coords_a[i]
        pt_bi = coords_b[i]
        diffs: list[float] = []
        for j in common:
            if is_covalent_deltaf508_pair(i, j, exclude_507_509=exclude_507_509):
                continue
            da = euclidean_distance(pt_ai, coords_a[j])
            db = euclidean_distance(pt_bi, coords_b[j])
            if da == 0.0 or db == 0.0:
                raise CryoEMF508delRescueError(f"Zero non-covalent distance between {i} and {j}")
            inc_d = min(da, db) if symmetric_inclusion else da
            margin_r0 = abs(inc_d - r0)
            if margin_r0 < min_margin:
                min_margin = margin_r0
            if inc_d <= r0:
                dd = abs(da - db)
                for t in thresholds:
                    mt = abs(dd - t)
                    if mt < min_margin:
                        min_margin = mt
                diffs.append(dd)
        neighbor_counts[i] = len(diffs)
        if not diffs:
            score_map[i] = None
        else:
            preserved = sum(sum(1 for t in thresholds if dd < t) for dd in diffs)
            score_map[i] = preserved / (float(len(thresholds)) * float(len(diffs)))

    if min_margin < MIN_CUTOFF_MARGIN_ANGSTROM and common:
        raise CryoEMF508delRescueError(
            f"Pairwise distance margin {min_margin:.3e} A violates {MIN_CUTOFF_MARGIN_ANGSTROM:.1e} A floor"
        )
    return score_map, neighbor_counts, (min_margin if common else float("inf"))


def load_catalog_permissions(catalog_path: Path = DEFAULT_CATALOG_YAML) -> dict[str, dict[str, str]]:
    entries: dict[str, dict[str, str]] = {}
    current_id: str | None = None
    for raw_line in catalog_path.read_text(encoding="utf-8").splitlines():
        s = raw_line.strip()
        if raw_line.startswith("  - id:"):
            current_id = raw_line.split(":", 1)[1].strip().strip("'\"")
            entries[current_id] = {}
        elif current_id is not None and raw_line.startswith("    ") and not raw_line.startswith("      ") and ":" in s:
            k, v = s.split(":", 1)
            entries[current_id][k.strip()] = v.strip().strip("'\"")
    return entries


def check_manual_catalog_gate(
    source_id: str,
    catalog_path: Path = DEFAULT_CATALOG_YAML,
    catalog_entries: Mapping[str, Mapping[str, str]] | None = None,
) -> None:
    if source_id in MANUAL_UNADMITTED_SOURCES:
        raise PermissionError(f"Source {source_id!r} is manual/unverified in catalog.yaml and cannot be auto-ingested")
    entries = catalog_entries if catalog_entries is not None else load_catalog_permissions(catalog_path)
    meta = entries.get(source_id)
    if meta is None:
        raise PermissionError(f"Unknown source {source_id!r} absent from catalog.yaml")
    access = meta.get("access", "")
    claim_label = meta.get("claim_label", "")
    if access not in ("fetch", "api") or claim_label != "verified":
        raise PermissionError(
            f"Source {source_id!r} not admitted for automated benchmark ingestion (access={access!r}, claim_label={claim_label!r})"
        )


def validate_holdout_annotations(
    holdout_payload: Mapping[str, Any],
    dev_payload: Mapping[str, Any],
    *,
    require_exposed_status: bool = False,
    expected_variants_sha256: str | None = EXPECTED_SHA256["HOLDOUT_VARIANTS"],
    catalog_entries: Mapping[str, Mapping[str, str]] | None = None,
) -> list[dict[str, Any]]:
    cat_map = catalog_entries if catalog_entries is not None else load_catalog_permissions()
    prov = holdout_payload.get("provenance", {})
    cat_sources = prov.get("catalog_sources")
    if not isinstance(cat_sources, list) or not cat_sources:
        raise CryoEMF508delRescueError("Holdout provenance.catalog_sources must be a non-empty list")
    for src_id in cat_sources:
        check_manual_catalog_gate(str(src_id), catalog_entries=cat_map)
    for dev_src_id in dev_payload.get("provenance", {}).get("catalog_sources", []):
        check_manual_catalog_gate(str(dev_src_id), catalog_entries=cat_map)

    exp_status = prov.get("panel_exposure_status")
    if exp_status not in ALLOWED_EXPOSURE_STATES:
        raise CryoEMF508delRescueError(f"Invalid panel_exposure_status: {exp_status!r}")
    if require_exposed_status and exp_status != "exposed_development_reference":
        raise CryoEMF508delRescueError(
            f"Evaluated holdout fixture must record panel_exposure_status='exposed_development_reference', got {exp_status!r}"
        )

    variants = holdout_payload.get("variants")
    if not isinstance(variants, list) or len(variants) != 24:
        raise CryoEMF508delRescueError(f"Expected 24 holdout variants, got {len(variants) if isinstance(variants, list) else None}")

    exposed_codons: set[int] = set(RESERVED_CONTROL_CODONS)
    for r in dev_payload.get("variants", []):
        exposed_codons.add(int(r["codon_position"]))
    for c in dev_payload.get("mcdonald_2024_table1_positive_controls", []):
        exposed_codons.add(int(c["codon_position"]))

    seen_ids: set[str] = set()
    seen_spans: list[tuple[int, int, str]] = []
    out: list[dict[str, Any]] = []

    for row in variants:
        vid = str(row["variant_id"])
        if vid in seen_ids:
            raise CryoEMF508delRescueError(f"Duplicate holdout variant_id: {vid}")
        seen_ids.add(vid)

        fam = str(row["variant_family"])
        if fam not in ALLOWED_FAMILIES:
            raise CryoEMF508delRescueError(f"Invalid variant_family {fam!r} on {vid}")

        raw_codon = row.get("codon_position")
        if isinstance(raw_codon, bool) or not isinstance(raw_codon, int):
            raise CryoEMF508delRescueError(f"codon_position must be int on {vid}, got {raw_codon!r}")
        codon = raw_codon
        hgvs = str(row.get("hgvs_p_or_c", ""))
        span_m = INDEL_SPAN_RE.match(hgvs)
        span_lo = int(span_m.group(1)) if span_m else codon
        span_hi = int(span_m.group(2)) if span_m else codon
        if not (1 <= span_lo <= codon <= span_hi <= 1480):
            raise CryoEMF508delRescueError(f"Codon/span out of range on {vid}: codon={codon}, span={span_lo}..{span_hi}")
        for span_pos in range(span_lo, span_hi + 1):
            for exc in exposed_codons:
                if abs(span_pos - exc) <= MIN_CODON_BUFFER:
                    raise CryoEMF508delRescueError(
                        f"Holdout variant {vid} (codon {span_pos}) within <={MIN_CODON_BUFFER} aa of forbidden codon {exc}"
                    )
            for prev_lo, prev_hi, prev_id in seen_spans:
                if min(abs(span_pos - p) for p in range(prev_lo, prev_hi + 1)) <= MIN_CODON_BUFFER:
                    raise CryoEMF508delRescueError(
                        f"Holdout codon collision <={MIN_CODON_BUFFER} aa: {vid} ({span_lo}..{span_hi}) vs {prev_id} ({prev_lo}..{prev_hi})"
                    )
        seen_spans.append((span_lo, span_hi, vid))

        ls = row.get("label_status")
        if not isinstance(ls, dict):
            raise CryoEMF508delRescueError(f"Missing label_status on {vid}")
        if ls.get("experimental_rescue_status") != "unknown":
            raise CryoEMF508delRescueError(
                f"Invalid experimental_rescue_status={ls.get('experimental_rescue_status')!r} on {vid} (SPL 0 is not experimental nonresponse)"
            )
        for flag_k, spl_k in (
            ("trikafta_v23_responsive", "trikafta_v23_spl"),
            ("alyftrek_v9_responsive", "alyftrek_v9_spl"),
            ("kalydeco_v36_responsive", "kalydeco_v36_spl"),
            ("symdeko_v22_responsive", "symdeko_v22_spl"),
        ):
            raw_bit = row[flag_k]
            if isinstance(raw_bit, bool) or not isinstance(raw_bit, int) or raw_bit not in (0, 1):
                raise CryoEMF508delRescueError(f"Invalid binary SPL flag {flag_k}={raw_bit!r} on {vid}")
            state = ls.get(spl_k)
            expected_state = "listed_on_fda_spl" if raw_bit == 1 else "not_listed_on_fda_spl"
            if state not in ALLOWED_SPL_STATES or state != expected_state:
                raise CryoEMF508delRescueError(f"SPL provenance mismatch on {vid}: {flag_k}={raw_bit}, {spl_k}={state!r}")

        if fam != "missense" and (row.get("alphamissense_score") is not None or row.get("revel_score") is not None):
            raise CryoEMF508delRescueError(f"Non-missense variant {vid} must carry None for missense-only scores")
        out.append(dict(row))

    canonical = json.dumps(variants, sort_keys=True, separators=(",", ":"))
    v_sha = hashlib.sha256(canonical.encode("utf-8")).hexdigest()
    if expected_variants_sha256 is not None:
        if v_sha != expected_variants_sha256 or prov.get("variants_sha256") != expected_variants_sha256:
            raise CryoEMF508delRescueError(
                f"Holdout variants_sha256 mismatch: computed {v_sha}, provenance {prov.get('variants_sha256')}"
            )
        if expected_variants_sha256 == EXPECTED_SHA256["HOLDOUT_VARIANTS"]:
            if prov.get("pre_registered_variants_sha256") != EXPECTED_SHA256["PRE_REGISTERED_HOLDOUT_VARIANTS"]:
                raise CryoEMF508delRescueError("Missing or mismatched pre_registered_variants_sha256 from Commit 9a23e16")
            allowed_errata_fields = {
                "kalydeco_v36_responsive",
                "symdeko_v22_responsive",
                "label_status.kalydeco_v36_spl",
                "label_status.symdeko_v22_spl",
            }
            reverted = copy.deepcopy(variants)
            by_id = {str(r["variant_id"]): r for r in reverted}
            for err in prov.get("post_audit_spl_errata", []):
                target_row = by_id.get(str(err.get("variant_id")))
                if target_row is None:
                    raise CryoEMF508delRescueError(f"Unknown errata variant_id: {err.get('variant_id')}")
                for fk, pair in err.get("fields_corrected", {}).items():
                    if fk not in allowed_errata_fields or not isinstance(pair, list) or len(pair) != 2:
                        raise CryoEMF508delRescueError(f"Disallowed errata field {fk!r}; evaluated target columns are immutable")
                    old_v, new_v = pair
                    if fk.startswith("label_status."):
                        sub_k = fk.split(".", 1)[1]
                        if target_row["label_status"].get(sub_k) != new_v:
                            raise CryoEMF508delRescueError(f"Errata current value mismatch on {fk}")
                        target_row["label_status"][sub_k] = old_v
                    else:
                        if target_row.get(fk) != new_v:
                            raise CryoEMF508delRescueError(f"Errata current value mismatch on {fk}")
                        target_row[fk] = old_v
            rev_sha = hashlib.sha256(json.dumps(reverted, sort_keys=True, separators=(",", ":")).encode("utf-8")).hexdigest()
            if rev_sha != EXPECTED_SHA256["PRE_REGISTERED_HOLDOUT_VARIANTS"]:
                raise CryoEMF508delRescueError(
                    f"Errata reversal digest {rev_sha} failed to reproduce pre-registered 9a23e16 digest"
                )

    return out


def score_variant_s_stab(
    row: Mapping[str, Any],
    s_stab_by_seq: Mapping[int, float | None],
    *,
    valid_s_common: frozenset[int] | None = None,
) -> tuple[float | None, float | None]:
    fam = str(row["variant_family"])
    codon = int(row["codon_position"])
    hgvs = str(row["hgvs_p_or_c"])
    if fam == "missense":
        if codon in CONSTRUCT_MODIFIED_CODONS:
            return None, None
        s_val = s_stab_by_seq.get(codon)
        if s_val is None:
            return None, None
        if isinstance(s_val, bool) or not isinstance(s_val, (int, float)) or not math.isfinite(float(s_val)):
            raise CryoEMF508delRescueError(f"Non-finite s_stab score on codon {codon}: {s_val!r}")
        s_f = float(s_val)
        if not (0.0 <= s_f <= 1.0):
            raise CryoEMF508delRescueError(f"s_stab score out of [0, 1] on codon {codon}: {s_f}")
        if valid_s_common is not None and codon not in valid_s_common:
            raise CryoEMF508delRescueError(
                f"Numeric fallback forbidden: codon {codon} is absent from S_common yet received s_stab={s_f}"
            )
        return s_f, 0.50 + 0.50 * s_f
    if fam == "class_i_ptc_frameshift" or CANONICAL_SPLICE_RE.match(hgvs):
        return None, 0.0
    if fam in ("inframe_indel", "splice_regulatory"):
        return None, 0.50
    raise CryoEMF508delRescueError(f"Unsupported variant_family: {fam!r}")


def score_variant_baselines(row: Mapping[str, Any]) -> tuple[float | None, float | None, float]:
    fam = str(row["variant_family"])
    hgvs = str(row["hgvs_p_or_c"])
    am = row.get("alphamissense_score")
    if fam == "missense":
        s_am = (1.0 - float(am)) if am is not None else None
        g_am = (0.50 + 0.50 * s_am) if s_am is not None else None
        return s_am, g_am, 0.75
    if fam == "class_i_ptc_frameshift" or CANONICAL_SPLICE_RE.match(hgvs):
        return None, 0.0, 0.0
    if fam in ("inframe_indel", "splice_regulatory"):
        return None, 0.50, 0.50
    raise CryoEMF508delRescueError(f"Unsupported variant_family: {fam!r}")


def compute_mann_whitney_auroc(
    pairs: Sequence[tuple[float | None, int]],
    *,
    abstention_aware: bool,
) -> float:
    for s, y in pairs:
        if isinstance(y, bool) or not isinstance(y, int) or y not in (0, 1):
            raise CryoEMF508delRescueError(f"Binary outcome y must be int in {{0, 1}}, got {y!r}")
        if s is not None and (isinstance(s, bool) or not isinstance(s, (int, float)) or not math.isfinite(float(s))):
            raise CryoEMF508delRescueError(f"Score must be None or finite float, got {s!r}")
    filtered = list(pairs) if abstention_aware else [(s, y) for s, y in pairs if s is not None]
    pos = [s for s, y in filtered if y == 1]
    neg = [s for s, y in filtered if y == 0]
    if not pos or not neg:
        raise CryoEMF508delRescueError("ROC-AUC requires at least one positive and one negative case")
    wins = 0.0
    for sp in pos:
        for sn in neg:
            if sp is None or sn is None:
                wins += 0.5
            elif sp > sn:
                wins += 1.0
            elif sp == sn:
                wins += 0.5
    return wins / float(len(pos) * len(neg))


def run_positive_controls(
    s_8ej1: DeltaF508ChainAStructure,
    s_8eiq: DeltaF508ChainAStructure,
    s_stab_map: Mapping[int, float | None],
    min_margin: float,
    dev_payload: Mapping[str, Any],
) -> dict[str, Any]:
    c_8ej1 = {k: pt for k, (_, pt) in s_8ej1.ca_by_seq.items()}
    c_8eiq = {k: pt for k, (_, pt) in s_8eiq.ca_by_seq.items()}
    common = sorted(set(c_8ej1) & set(c_8eiq))
    if len(c_8ej1) != 1127 or len(c_8eiq) != 1162 or len(common) != 1121:
        raise CryoEMF508delRescueError(
            f"Unexpected Chain A counts: 8EJ1={len(c_8ej1)}, 8EIQ={len(c_8eiq)}, S_common={len(common)}"
        )

    # Suite 2: Published deltaF508 Trikafta NBD1:NBD2 dimer closure contacts
    d_site1_apo = euclidean_distance(c_8ej1[464], c_8ej1[1348])
    d_site1_eti = euclidean_distance(c_8eiq[464], c_8eiq[1348])
    d_site2_apo = euclidean_distance(c_8ej1[582], c_8ej1[1246])
    d_site2_eti = euclidean_distance(c_8eiq[582], c_8eiq[1246])
    if (d_site1_apo - d_site1_eti) <= 15.0 or d_site1_eti >= 12.0:
        raise CryoEMF508delRescueError("ATP Site 1 (K464, H1348) dimer-closure control failed")
    if (d_site2_apo - d_site2_eti) <= 15.0 or d_site2_eti >= 12.0:
        raise CryoEMF508delRescueError("ATP Site 2 (T582, G1246) dimer-closure control failed")

    # Suite 3: Rigid corrector/helical scaffolding vs ATP-dimerization interface ordering
    scaf_vals = [float(s_stab_map[k]) for k in SCAFFOLD_CONTROL_CODONS if s_stab_map.get(k) is not None]
    atp_vals = [float(s_stab_map[k]) for k in ATP_SITE_CONTROL_CODONS if s_stab_map.get(k) is not None]
    if len(scaf_vals) != len(SCAFFOLD_CONTROL_CODONS) or len(atp_vals) != len(ATP_SITE_CONTROL_CODONS):
        raise CryoEMF508delRescueError("Missing reserved positive-control residue in S_common")
    min_scaf = min(scaf_vals)
    max_atp = max(atp_vals)
    if min_scaf <= 0.70 or max_atp >= 0.56 or min_scaf <= max_atp:
        raise CryoEMF508delRescueError(
            f"Scaffolding vs ATP-site ordering failed: min_scaf={min_scaf:.6f}, max_atp={max_atp:.6f}"
        )

    # Suite 4: Independent SE(3) x SE(3) invariance, swap symmetry, and dev_exposed calibration lock
    rot1 = ((3.0 / 5.0, -4.0 / 5.0, 0.0), (4.0 / 5.0, 3.0 / 5.0, 0.0), (0.0, 0.0, 1.0))
    rot2 = (
        (2.0 / 3.0, -2.0 / 3.0, 1.0 / 3.0),
        (2.0 / 3.0, 1.0 / 3.0, -2.0 / 3.0),
        (1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0),
    )
    t_8ej1 = apply_rigid_motion(c_8ej1, rot1, (17.25, -9.5, 42.125))
    t_8eiq = apply_rigid_motion(c_8eiq, rot2, (-31.0, 14.75, -8.25))
    moved_map, _, _ = compute_lddt_sym_conservation_map(t_8ej1, t_8eiq)
    max_se3_drift = max(abs(float(moved_map[k]) - float(s_stab_map[k])) for k in common)
    if max_se3_drift >= 1e-12:
        raise CryoEMF508delRescueError(f"SE(3)xSE(3) drift {max_se3_drift:.3e} >= 1e-12")

    swapped_map, _, _ = compute_lddt_sym_conservation_map(c_8eiq, c_8ej1)
    max_swap_diff = max(abs(float(swapped_map[k]) - float(s_stab_map[k])) for k in common)
    if max_swap_diff != 0.0:
        raise CryoEMF508delRescueError(f"Swap symmetry violated: {max_swap_diff:.3e}")

    s_a455 = float(s_stab_map[455] or 0.0)
    s_l467 = float(s_stab_map[467] or 0.0)
    if not (s_a455 > s_l467):
        raise CryoEMF508delRescueError(f"dev_exposed A455E ({s_a455}) failed to exceed L467P ({s_l467})")

    # Confirm dev_exposed calibration directional lock (s_stab = 0.575758 > 1 - s_stab = 0.424242)
    dev_miss = [
        (score_variant_s_stab(r, s_stab_map)[0], int(r["trikafta_v23_responsive"]))
        for r in dev_payload["variants"]
        if r["variant_family"] == "missense"
    ]
    dev_cond_auroc = compute_mann_whitney_auroc(dev_miss, abstention_aware=False)
    if not (dev_cond_auroc > 0.50):
        raise CryoEMF508delRescueError(f"dev_exposed directional calibration failed: {dev_cond_auroc:.6f}")

    # Suite 5: Analytic 4-point delta508 chain control & ROC-AUC unit controls
    chain_d508 = {506: (0.0, 0.0, 0.0), 507: (3.8, 0.0, 0.0), 509: (7.6, 0.0, 0.0), 510: (11.4, 0.0, 0.0)}
    _, deg_excl, _ = compute_lddt_sym_conservation_map(chain_d508, chain_d508, exclude_507_509=True)
    _, deg_leak, _ = compute_lddt_sym_conservation_map(chain_d508, chain_d508, exclude_507_509=False)
    if deg_excl[507] != 1 or deg_excl[509] != 1 or deg_leak[507] != 2:
        raise CryoEMF508delRescueError(f"Analytic delta508 junction neighbor control failed: {deg_excl} vs {deg_leak}")

    auc_none = compute_mann_whitney_auroc([(None, 1), (None, 0)], abstention_aware=True)
    auc_pos = compute_mann_whitney_auroc([(0.9, 1), (0.1, 0)], abstention_aware=True)
    auc_inv = compute_mann_whitney_auroc([(0.1, 1), (0.9, 0)], abstention_aware=True)
    if auc_none != 0.5 or auc_pos != 1.0 or auc_inv != 0.0:
        raise CryoEMF508delRescueError(f"ROC-AUC unit controls failed: {auc_none}, {auc_pos}, {auc_inv}")

    return {
        "placed_ca_8ej1": len(c_8ej1),
        "placed_ca_8eiq": len(c_8eiq),
        "shared_s_common": len(common),
        "only_8eiq_count": len(set(c_8eiq) - set(c_8ej1)),
        "only_8ej1_count": len(set(c_8ej1) - set(c_8eiq)),
        "min_cutoff_margin_angstrom": min_margin,
        "atp_site1_contraction_angstrom": round(d_site1_apo - d_site1_eti, 6),
        "atp_site2_contraction_angstrom": round(d_site2_apo - d_site2_eti, 6),
        "min_scaffold_s_stab": round(min_scaf, 6),
        "max_atp_site_s_stab": round(max_atp, 6),
        "max_se3_x_se3_drift": max_se3_drift,
        "max_swap_symmetry_diff": max_swap_diff,
        "dev_exposed_a455_s_stab": round(s_a455, 6),
        "dev_exposed_l467_s_stab": round(s_l467, 6),
        "dev_exposed_cond_missense_auroc": round(dev_cond_auroc, 6),
    }


def evaluate_holdout_panel(
    variants: Sequence[Mapping[str, Any]],
    s_stab_map: Mapping[int, float | None],
    target_key: str = "trikafta_v23_responsive",
    *,
    valid_s_common: frozenset[int] | None = None,
) -> dict[str, Any]:
    if valid_s_common is not None:
        for k, v in s_stab_map.items():
            if v is not None and k not in valid_s_common:
                raise CryoEMF508delRescueError(
                    f"Numeric fallback forbidden in s_stab_map: residue {k} outside S_common has value {v!r}"
                )

    evaluated_rows: list[dict[str, Any]] = []
    missense_placed_pairs: list[tuple[float | None, int]] = []
    missense_all14_stab_pairs: list[tuple[float | None, int]] = []
    missense_am_pairs: list[tuple[float | None, int]] = []
    all24_g_stab_pairs: list[tuple[float | None, int]] = []
    all24_g_am_pairs: list[tuple[float | None, int]] = []
    all24_cls_pairs: list[tuple[float | None, int]] = []

    unmodeled_missense: list[str] = []

    for row in variants:
        if target_key not in row:
            raise CryoEMF508delRescueError(f"Missing target_key {target_key!r} in variant row {row.get('variant_id')}")
        raw_y = row[target_key]
        if isinstance(raw_y, bool) or not isinstance(raw_y, int) or raw_y not in (0, 1):
            raise CryoEMF508delRescueError(f"Invalid binary target {target_key}={raw_y!r}")
        y = raw_y
        s_stab, g_stab = score_variant_s_stab(row, s_stab_map, valid_s_common=valid_s_common)
        s_am, g_am, cls_b = score_variant_baselines(row)
        if row["variant_family"] == "missense":
            missense_all14_stab_pairs.append((s_stab, y))
            missense_am_pairs.append((s_am, y))
            if s_stab is not None:
                missense_placed_pairs.append((s_stab, y))
            else:
                unmodeled_missense.append(str(row["legacy_spl_name"]))
        all24_g_stab_pairs.append((g_stab, y))
        all24_g_am_pairs.append((g_am, y))
        all24_cls_pairs.append((cls_b, y))

        evaluated_rows.append({
            "variant_id": row["variant_id"],
            "legacy_spl_name": row["legacy_spl_name"],
            "codon_position": int(row["codon_position"]),
            "variant_family": row["variant_family"],
            "target_spl_listed": y,
            "s_stab": round(s_stab, 6) if s_stab is not None else None,
            "g_stab": round(g_stab, 6) if g_stab is not None else None,
            "s_am": round(s_am, 6) if s_am is not None else None,
            "g_am": round(g_am, 6) if g_am is not None else None,
            "class_baseline": cls_b,
        })

    cond_placed_auroc = compute_mann_whitney_auroc(missense_placed_pairs, abstention_aware=False)
    missense_all14_abst_auroc = compute_mann_whitney_auroc(missense_all14_stab_pairs, abstention_aware=True)
    missense_am_auroc = compute_mann_whitney_auroc(missense_am_pairs, abstention_aware=True)
    all24_g_stab_auroc = compute_mann_whitney_auroc(all24_g_stab_pairs, abstention_aware=True)
    all24_g_am_auroc = compute_mann_whitney_auroc(all24_g_am_pairs, abstention_aware=True)
    all24_cls_auroc = compute_mann_whitney_auroc(all24_cls_pairs, abstention_aware=True)

    # Pre-registered PRED-0014 decision rule (done.md line 18):
    supported = (
        missense_all14_abst_auroc >= 0.75
        and missense_all14_abst_auroc > missense_am_auroc
        and all24_g_stab_auroc >= 0.75
        and all24_g_stab_auroc > all24_cls_auroc
        and all24_g_stab_auroc > all24_g_am_auroc
    )

    return {
        "target_key": target_key,
        "total_variants": len(variants),
        "positive_count": sum(int(r[target_key]) == 1 for r in variants),
        "negative_count": sum(int(r[target_key]) == 0 for r in variants),
        "missense_total": len(missense_all14_stab_pairs),
        "missense_placed_in_s_common": len(missense_placed_pairs),
        "unmodeled_missense_variants": unmodeled_missense,
        "conditional_placed_missense_auroc_s_stab": round(cond_placed_auroc, 6),
        "missense_all14_abstention_auroc_s_stab": round(missense_all14_abst_auroc, 6),
        "standalone_missense_auroc_s_am": round(missense_am_auroc, 6),
        "gated_all24_abstention_auroc_g_stab": round(all24_g_stab_auroc, 6),
        "gated_all24_abstention_auroc_g_am": round(all24_g_am_auroc, 6),
        "sequence_consequence_baseline_all24_auroc": round(all24_cls_auroc, 6),
        "prediction_outcome": "supported" if supported else "refuted",
        "variants": evaluated_rows,
    }


def run_falsification_suite(
    s_8ej1: DeltaF508ChainAStructure,
    s_8eiq: DeltaF508ChainAStructure,
    uniprot_seq: str,
    s_stab_map: Mapping[int, float | None],
    holdout_payload: Mapping[str, Any],
    dev_payload: Mapping[str, Any],
    holdout_variants: Sequence[Mapping[str, Any]],
) -> dict[str, str]:
    verdicts: dict[str, str] = {}
    with gzip.open(DEFAULT_CIF_8EIQ, "rt", encoding="utf-8") as fh_8eiq:
        raw_8eiq_text = fh_8eiq.read()
    c_8ej1 = {k: pt for k, (_, pt) in s_8ej1.ca_by_seq.items()}
    c_8eiq = {k: pt for k, (_, pt) in s_8eiq.ca_by_seq.items()}
    common_set = frozenset(set(c_8ej1) & set(c_8eiq))

    # 1. f508del_507_509_peptide_bond_or_shift_rejected
    c1_ok = True
    try:
        parse_deltaf508_chain_a_mmcif(raw_8eiq_text, "8EIQ", uniprot_seq, use_label_seq_id=True)
        c1_ok = False
    except CryoEMF508delRescueError:
        pass
    first_atom_nl = raw_8eiq_text.find("\nATOM ")
    next_nl = raw_8eiq_text.find("\n", first_atom_nl + 1) if first_atom_nl >= 0 else -1
    if first_atom_nl < 0 or next_nl < 0:
        raise CryoEMF508delRescueError("Failed to locate _atom_site ATOM record in 8EIQ")
    sample_atom = raw_8eiq_text[first_atom_nl + 1 : next_nl]
    stoks = shlex.split(sample_atom)
    a_cols = [ln.strip().split(".", 1)[1] for ln in raw_8eiq_text.splitlines() if ln.strip().startswith("_atom_site.")]
    fake_map = dict(zip(a_cols, stoks))
    fake_map["id"] = "999999"
    fake_map["label_atom_id"] = "CA"
    fake_map["auth_atom_id"] = "CA"
    fake_map["label_comp_id"] = "PHE"
    fake_map["auth_comp_id"] = "PHE"
    fake_map["auth_asym_id"] = "A"
    fake_map["label_asym_id"] = "A"
    fake_map["auth_seq_id"] = "508"
    fake_map["label_seq_id"] = "508"
    fake_line = " ".join(fake_map[c] for c in a_cols)
    injected_508 = raw_8eiq_text[:next_nl] + "\n" + fake_line + raw_8eiq_text[next_nl:]
    try:
        parse_deltaf508_chain_a_mmcif(injected_508, "8EIQ", uniprot_seq)
        c1_ok = False
    except CryoEMF508delRescueError:
        pass
    chain_d508 = {506: (0.0, 0.0, 0.0), 507: (3.8, 0.0, 0.0), 509: (7.6, 0.0, 0.0), 510: (11.4, 0.0, 0.0)}
    _, deg_excl, _ = compute_lddt_sym_conservation_map(chain_d508, chain_d508, exclude_507_509=True)
    _, deg_leak, _ = compute_lddt_sym_conservation_map(chain_d508, chain_d508, exclude_507_509=False)
    if deg_excl[507] == deg_leak[507]:
        c1_ok = False
    if not c1_ok:
        raise CryoEMF508delRescueError("Control 1 failed")
    verdicts["f508del_507_509_peptide_bond_or_shift_rejected"] = "PASS"

    # 2. hetatm_ligand_or_chain_b_leak_rejected (exercising downstream sequence/range checks without flag short-circuits)
    c2_ok = True
    try:
        parse_deltaf508_chain_a_mmcif(raw_8eiq_text, "8EIQ", uniprot_seq, allow_hetatm=True)
        c2_ok = False
    except CryoEMF508delRescueError:
        pass
    hetatm_ca_map = dict(fake_map)
    hetatm_ca_map["group_PDB"] = "HETATM"
    hetatm_ca_map["label_comp_id"] = "VX7"
    hetatm_ca_map["auth_comp_id"] = "VX7"
    hetatm_ca_map["auth_seq_id"] = "900"
    hetatm_ca_map["label_seq_id"] = "900"
    injected_hetatm_ca = raw_8eiq_text[:next_nl] + "\n" + " ".join(hetatm_ca_map[c] for c in a_cols) + raw_8eiq_text[next_nl:]
    try:
        parse_deltaf508_chain_a_mmcif(injected_hetatm_ca, "8EIQ", uniprot_seq, allow_hetatm=True)
        c2_ok = False
    except CryoEMF508delRescueError:
        pass
    try:
        parse_deltaf508_chain_a_mmcif(raw_8eiq_text, "8EIQ", uniprot_seq, allow_non_chain_a=True)
        c2_ok = False
    except CryoEMF508delRescueError:
        pass
    if not c2_ok:
        raise CryoEMF508delRescueError("Control 2 failed")
    verdicts["hetatm_ligand_or_chain_b_leak_rejected"] = "PASS"

    # 3. unmodeled_loop_numeric_fallback_rejected (exercising downstream evaluate_holdout_panel & score_variant_s_stab)
    c3_ok = True
    for imputed_val in (0.0, 1.0):
        poisoned_map = dict(s_stab_map)
        poisoned_map[110] = imputed_val  # ETI-ordered ECL1 residue 110 in C_8EIQ \ C_8EJ1
        try:
            evaluate_holdout_panel(holdout_variants, poisoned_map, valid_s_common=common_set)
            c3_ok = False
        except CryoEMF508delRescueError:
            pass
    s_110, g_110 = score_variant_s_stab(
        {"variant_family": "missense", "codon_position": 110, "hgvs_p_or_c": "p.Asp110His"},
        s_stab_map,
        valid_s_common=common_set,
    )
    if s_110 is not None or g_110 is not None:
        c3_ok = False
    isolated = {10: (0.0, 0.0, 0.0), 30: (50.0, 0.0, 0.0)}
    iso_map, _, _ = compute_lddt_sym_conservation_map(isolated, isolated)
    if iso_map[10] is not None:
        c3_ok = False
    if not c3_ok:
        raise CryoEMF508delRescueError("Control 3 failed")
    verdicts["unmodeled_loop_numeric_fallback_rejected"] = "PASS"

    # 4. asymmetric_neighborhood_swap_rejected
    asym_ab, _, _ = compute_lddt_sym_conservation_map(c_8ej1, c_8eiq, symmetric_inclusion=False)
    asym_ba, _, _ = compute_lddt_sym_conservation_map(c_8eiq, c_8ej1, symmetric_inclusion=False)
    max_asym_diff = max(abs(float(asym_ab[k]) - float(asym_ba[k])) for k in asym_ab)
    if max_asym_diff <= 0.0:
        raise CryoEMF508delRescueError("Control 4 failed: asymmetric neighborhood unexpectedly swap-symmetric")
    verdicts["asymmetric_neighborhood_swap_rejected"] = "PASS"

    # 5. non_missense_or_construct_codon_leak_rejected
    for synthetic_row in (
        {"variant_family": "inframe_indel", "codon_position": 1140, "hgvs_p_or_c": "p.Met1140del"},
        {"variant_family": "splice_regulatory", "codon_position": 1047, "hgvs_p_or_c": "c.3140-26A>G"},
        {"variant_family": "class_i_ptc_frameshift", "codon_position": 831, "hgvs_p_or_c": "p.Glu831Ter"},
        {"variant_family": "missense", "codon_position": 1371, "hgvs_p_or_c": "p.Glu1371Gln"},
        {"variant_family": "missense", "codon_position": 508, "hgvs_p_or_c": "p.Phe508Cys"},
    ):
        sv, _ = score_variant_s_stab(synthetic_row, s_stab_map)
        if sv is not None:
            raise CryoEMF508delRescueError(f"Control 5 failed on {synthetic_row}")
    verdicts["non_missense_or_construct_codon_leak_rejected"] = "PASS"

    # 6. post_hoc_sign_flip_rejected
    inv_auc = compute_mann_whitney_auroc([(0.1, 1), (0.2, 1), (0.8, 0), (0.9, 0)], abstention_aware=False)
    if inv_auc != 0.0:
        raise CryoEMF508delRescueError(f"Control 6 failed: inverted separator gave {inv_auc}")
    verdicts["post_hoc_sign_flip_rejected"] = "PASS"

    # 7. codon_buffer_exposure_or_manual_source_rejected
    c7_ok = True
    bad_buf = copy.deepcopy(dict(holdout_payload))
    bad_buf["variants"] = copy.deepcopy(list(holdout_payload["variants"]))
    bad_buf["variants"][0]["codon_position"] = 69  # within 2 aa of exposed P67L (67)
    try:
        validate_holdout_annotations(bad_buf, dev_payload, expected_variants_sha256=None)
        c7_ok = False
    except CryoEMF508delRescueError as exc:
        if "within <=" not in str(exc):
            c7_ok = False

    bad_exp_rescue = copy.deepcopy(dict(holdout_payload))
    bad_exp_rescue["variants"] = copy.deepcopy(list(holdout_payload["variants"]))
    bad_exp_rescue["variants"][1]["label_status"]["experimental_rescue_status"] = "measured_experimental_nonresponse"
    try:
        validate_holdout_annotations(bad_exp_rescue, dev_payload, expected_variants_sha256=None)
        c7_ok = False
    except CryoEMF508delRescueError as exc:
        if "experimental_rescue_status" not in str(exc):
            c7_ok = False

    bad_sha = copy.deepcopy(dict(holdout_payload))
    bad_sha["provenance"] = dict(holdout_payload["provenance"])
    bad_sha["provenance"]["variants_sha256"] = "0" * 64
    try:
        validate_holdout_annotations(bad_sha, dev_payload)
        c7_ok = False
    except CryoEMF508delRescueError as exc:
        if "variants_sha256 mismatch" not in str(exc):
            c7_ok = False

    bad_status = copy.deepcopy(dict(holdout_payload))
    bad_status["provenance"] = dict(holdout_payload["provenance"])
    bad_status["provenance"]["panel_exposure_status"] = "untouched_holdout"
    try:
        validate_holdout_annotations(bad_status, dev_payload, require_exposed_status=True)
        c7_ok = False
    except CryoEMF508delRescueError as exc:
        if "panel_exposure_status" not in str(exc):
            c7_ok = False

    cat_entries = load_catalog_permissions()
    bad_prov_src = copy.deepcopy(dict(holdout_payload))
    bad_prov_src["provenance"] = dict(holdout_payload["provenance"])
    bad_prov_src["provenance"]["catalog_sources"] = ["myvariant", "mavedb-dump-zenodo"]
    try:
        validate_holdout_annotations(bad_prov_src, dev_payload, catalog_entries=cat_entries)
        c7_ok = False
    except PermissionError:
        pass

    for src_id in MANUAL_UNADMITTED_SOURCES:
        try:
            check_manual_catalog_gate(src_id, catalog_entries=cat_entries)
            c7_ok = False
        except PermissionError:
            pass
    if not c7_ok:
        raise CryoEMF508delRescueError("Control 7 failed")
    verdicts["codon_buffer_exposure_or_manual_source_rejected"] = "PASS"

    return verdicts


def run_benchmark(
    *,
    run_falsify: bool = False,
    require_exposed_status: bool = False,
) -> dict[str, Any]:
    cat_entries = load_catalog_permissions()
    for admitted_id in ("rcsb-pdb-8ej1", "rcsb-pdb-8eiq", "uniprot-p13569", "openfda-label", "myvariant", "clinvar"):
        check_manual_catalog_gate(admitted_id, catalog_entries=cat_entries)

    s_8ej1, s_8eiq, uniprot_seq = load_verified_structures()
    c_8ej1 = {k: pt for k, (_, pt) in s_8ej1.ca_by_seq.items()}
    c_8eiq = {k: pt for k, (_, pt) in s_8eiq.ca_by_seq.items()}
    common_set = frozenset(set(c_8ej1) & set(c_8eiq))
    s_stab_map, _, min_margin = compute_lddt_sym_conservation_map(c_8ej1, c_8eiq)

    dev_payload = json.loads(DEFAULT_DEV_JSON.read_text(encoding="utf-8"))
    holdout_payload = json.loads(DEFAULT_HOLDOUT_JSON.read_text(encoding="utf-8"))

    pos_controls = run_positive_controls(s_8ej1, s_8eiq, s_stab_map, min_margin, dev_payload)
    holdout_variants = validate_holdout_annotations(
        holdout_payload,
        dev_payload,
        require_exposed_status=require_exposed_status,
        catalog_entries=cat_entries,
    )

    trikafta_eval = evaluate_holdout_panel(
        holdout_variants,
        s_stab_map,
        "trikafta_v23_responsive",
        valid_s_common=common_set,
    )
    alyftrek_eval = evaluate_holdout_panel(
        holdout_variants,
        s_stab_map,
        "alyftrek_v9_responsive",
        valid_s_common=common_set,
    )

    falsify_verdicts = (
        run_falsification_suite(
            s_8ej1,
            s_8eiq,
            uniprot_seq,
            s_stab_map,
            holdout_payload,
            dev_payload,
            holdout_variants,
        )
        if run_falsify
        else {}
    )

    return {
        "benchmark_id": "PRED-0014",
        "research_boundary": RESEARCH_BOUNDARY,
        "positive_controls": pos_controls,
        "holdout_trikafta_v23": trikafta_eval,
        "holdout_alyftrek_v9_secondary": alyftrek_eval,
        "falsification_controls": falsify_verdicts,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="Verify positive controls, provenance, and holdout evaluation")
    parser.add_argument("--falsify", action="store_true", help="Run all 7 planted negative controls")
    parser.add_argument("--json", action="store_true", help="Print full JSON benchmark report")
    args = parser.parse_args(argv)

    rep = run_benchmark(
        run_falsify=args.falsify or args.check,
        require_exposed_status=args.check,
    )
    if args.json:
        print(json.dumps(rep, indent=2))
        return 0

    t_ev = rep["holdout_trikafta_v23"]
    print(
        f"PASS cryoem_f508del_rescue_benchmark (PRED-0014): outcome={t_ev['prediction_outcome']} | "
        f"S_common={rep['positive_controls']['shared_s_common']} | "
        f"placed_missense_AUC={t_ev['conditional_placed_missense_auroc_s_stab']:.6f} "
        f"(all14_abst={t_ev['missense_all14_abstention_auroc_s_stab']:.6f} vs s_AM={t_ev['standalone_missense_auroc_s_am']:.6f}) | "
        f"all24_g_stab_AUC={t_ev['gated_all24_abstention_auroc_g_stab']:.6f} "
        f"(vs class={t_ev['sequence_consequence_baseline_all24_auroc']:.6f}, g_AM={t_ev['gated_all24_abstention_auroc_g_am']:.6f}) | "
        f"falsify={len(rep['falsification_controls'])}/7"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
