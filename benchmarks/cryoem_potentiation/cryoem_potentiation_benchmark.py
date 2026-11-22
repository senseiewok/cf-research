#!/usr/bin/env python3
"""Cryo-EM 3D Potentiation-vs-Folding Pocket Proximity Benchmark (PRED-0011).

Tests whether E(3)-invariant intra-8EIQ 3D pocket proximity
    s_3D(i) = d_fold(i) / (d_pot(i) + d_fold(i))
predicts US FDA KALYDECO v36 Table 3 ivacaftor monotherapy responsiveness across the
fixed 30-variant domain-disjoint split (2026-11-21.v2: 10 train in Lasso/TMD1/R_domain,
20 held-out test in NBD1/TMD2/NBD2_C_terminal) better than the deterministic
sequence-consequence class baseline and standalone/gated AlphaMissense.

Method 1 implementation details (disjoint from independent_audit.py):
- Parses 8EIQ.cif.gz via shlex-tokenized mmCIF loops and typed dataclasses.
- Indexes 3D heavy-atom anchor sets in an 8.0 A spatial voxel hash grid with exact
  nearest-neighbor bounding-sphere pruning.
- Evaluates conditional and abstention-aware ROC-AUC via O(N1*N0) pairwise Mann-Whitney
  concordance enumeration (no post-hoc sign flipping).

100% pure Python 3.10+ standard library, offline, zero patient/PHI data.
"""
from __future__ import annotations

import argparse
import dataclasses
import gzip
import hashlib
import json
import math
import re
import shlex
import sys
from pathlib import Path
from typing import Any, Callable, Iterable

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
DEFAULT_CIF_8EIQ = REPO_ROOT / "tools" / "cftr-model" / "data" / "8EIQ.cif.gz"
DEFAULT_UNIPROT_JSON = REPO_ROOT / "tools" / "cftr-model" / "data" / "p13569.json"
DEFAULT_CATALOG_YAML = REPO_ROOT / "sources" / "catalog.yaml"
DEFAULT_ANNOTATIONS_JSON = (
    REPO_ROOT / "benchmarks" / "theratype" / "theratype_open_annotations_2026-11-21.json"
)

EXPECTED_8EIQ_SHA256 = "d13cbfc8ba6575ab6228f351578c6f9cc5b9dafe13e5913d9b27b3cacebc18d6"
EXPECTED_P13569_SHA256 = "c8feabfabe8690ef4a06af7de5b5f8dc3c8da554246c5ab3d403edf265d4a893"

THREE_TO_ONE: dict[str, str] = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

POTENTIATION_HETATM_IDS: frozenset[str] = frozenset({"ATP", "MG", "VX7"})
FOLDING_HETATM_IDS: frozenset[str] = frozenset({"CV6", "WJX"})
FOLDING_SOCKET_RESIDUES: frozenset[int] = frozenset({1063, 1070})

TRAIN_DOMAINS: frozenset[str] = frozenset({"Lasso_N_terminal", "TMD1", "R_domain"})
TEST_DOMAINS: frozenset[str] = frozenset({"NBD1", "TMD2", "NBD2_C_terminal"})
ALLOWED_FAMILIES: frozenset[str] = frozenset(
    {"missense", "inframe_indel", "splice_regulatory", "class_i_ptc_frameshift"}
)
CANONICAL_SPLICE_RE = re.compile(r"^c\.\d+[+-][12](?:[^\d]|$)")
CATALOG_ID_RE = re.compile(r"^\s*-\s*id:\s*([^\s#]+)")
CATALOG_ACCESS_RE = re.compile(r"^\s*(?:access|mode):\s*([^\s#]+)")
CATALOG_LABEL_RE = re.compile(r"^\s*claim_label:\s*([^\s#]+)")


class CryoEMPotentiationError(ValueError):
    """Raised when structural provenance, positive controls, or negative controls fail."""


@dataclasses.dataclass(frozen=True)
class HeavyAtom:
    group: str
    comp_id: str
    atom_id: str
    seq_id: int | None
    xyz: tuple[float, float, float]


@dataclasses.dataclass(frozen=True)
class Structure8EIQ:
    ca_by_seq: dict[int, tuple[str, tuple[float, float, float]]]
    potentiation_atoms: tuple[HeavyAtom, ...]
    folding_atoms: tuple[HeavyAtom, ...]
    comp_counts: dict[str, int]
    socket_counts: dict[int, int]
    declared_mutations: dict[int, str]
    declared_deletions: frozenset[int]


class VoxelGrid3D:
    """Spatial 8.0 A cubic voxel hash grid for exact nearest-neighbor queries."""

    def __init__(self, points: Iterable[tuple[float, float, float]], cell_size: float = 8.0) -> None:
        self.points: tuple[tuple[float, float, float], ...] = tuple(points)
        if not self.points:
            raise CryoEMPotentiationError("VoxelGrid3D initialized with empty anchor set")
        self.cell_size = float(cell_size)
        buckets: dict[tuple[int, int, int], list[tuple[float, float, float]]] = {}
        for pt in self.points:
            key = self._cell(pt)
            buckets.setdefault(key, []).append(pt)
        self.buckets = buckets

    def _cell(self, pt: tuple[float, float, float]) -> tuple[int, int, int]:
        cs = self.cell_size
        return (math.floor(pt[0] / cs), math.floor(pt[1] / cs), math.floor(pt[2] / cs))

    def min_distance(self, query: tuple[float, float, float]) -> float:
        qx, qy, qz = self._cell(query)
        cs = self.cell_size
        best = math.inf
        for radius in range(0, 16):
            for dx in range(-radius, radius + 1):
                for dy in range(-radius, radius + 1):
                    for dz in range(-radius, radius + 1):
                        if max(abs(dx), abs(dy), abs(dz)) != radius:
                            continue
                        for pt in self.buckets.get((qx + dx, qy + dy, qz + dz), ()):
                            d = math.dist(query, pt)
                            if d < best:
                                best = d
            if math.isfinite(best) and (radius * cs) > best + cs:
                return best
        return min(math.dist(query, pt) for pt in self.points)


def load_admitted_automated_catalog_ids(catalog_path: Path = DEFAULT_CATALOG_YAML) -> frozenset[str]:
    """Allowlist strictly sources in sources/catalog.yaml with access in {'fetch', 'api'} and claim_label != 'unverified'."""
    if not catalog_path.is_file():
        raise CryoEMPotentiationError(f"Catalog YAML not found: {catalog_path}")
    admitted: set[str] = set()
    cur_id: str | None = None
    cur_access: str | None = None
    cur_label: str | None = None

    def flush() -> None:
        if cur_id and cur_access in {"fetch", "api"} and cur_label != "unverified":
            admitted.add(cur_id)

    for ln in catalog_path.read_text(encoding="utf-8").splitlines():
        m_id = CATALOG_ID_RE.match(ln)
        if m_id:
            flush()
            cur_id = m_id.group(1).strip("'\"")
            cur_access = None
            cur_label = None
            continue
        m_acc = CATALOG_ACCESS_RE.match(ln)
        if m_acc and cur_access is None:
            val = m_acc.group(1).strip("'\"")
            if val in {"fetch", "api", "manual", "forbidden", "request"}:
                cur_access = val
        m_lbl = CATALOG_LABEL_RE.match(ln)
        if m_lbl:
            cur_label = m_lbl.group(1).strip("'\"")
    flush()
    return frozenset(admitted)


def map_p13569_domain(codon: int) -> str:
    if 1 <= codon <= 80:
        return "Lasso_N_terminal"
    if 81 <= codon <= 380:
        return "TMD1"
    if 381 <= codon <= 630:
        return "NBD1"
    if 631 <= codon <= 845:
        return "R_domain"
    if 846 <= codon <= 1190:
        return "TMD2"
    if 1191 <= codon <= 1480:
        return "NBD2_C_terminal"
    raise CryoEMPotentiationError(f"Codon {codon} outside human CFTR 1..1480")


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_uniprot_sequence(uniprot_path: Path = DEFAULT_UNIPROT_JSON, verify_hash: bool = True) -> str:
    if verify_hash:
        actual = sha256_file(uniprot_path)
        if actual != EXPECTED_P13569_SHA256:
            raise CryoEMPotentiationError(
                f"p13569.json SHA-256 mismatch: expected {EXPECTED_P13569_SHA256}, got {actual}"
            )
    payload = json.loads(uniprot_path.read_text(encoding="utf-8"))
    seq = str(payload["sequence"]["value"])
    if len(seq) != 1480:
        raise CryoEMPotentiationError(f"Expected 1480-aa P13569 sequence, got {len(seq)}")
    return seq


def verify_ca_against_uniprot(
    ca_by_seq: dict[int, tuple[str, tuple[float, float, float]]],
    uniprot_seq: str,
    declared_mutations: dict[int, str],
    declared_deletions: frozenset[int],
    seq_id_column: str = "auth_seq_id",
) -> None:
    for del_pos in declared_deletions:
        if del_pos in ca_by_seq:
            raise CryoEMPotentiationError(
                f"Deleted construct residue {del_pos} unexpectedly placed in Chain A ({seq_id_column})"
            )
    mismatches: list[tuple[int, str, str]] = []
    for pos, (comp_id, _) in sorted(ca_by_seq.items()):
        if pos < 1 or pos > len(uniprot_seq):
            mismatches.append((pos, comp_id, "OUT_OF_RANGE"))
            continue
        if pos in declared_mutations:
            if comp_id != declared_mutations[pos]:
                mismatches.append((pos, comp_id, f"DECLARED:{declared_mutations[pos]}"))
            continue
        aa = THREE_TO_ONE.get(comp_id)
        if aa != uniprot_seq[pos - 1]:
            mismatches.append((pos, comp_id, uniprot_seq[pos - 1]))
    if mismatches:
        raise CryoEMPotentiationError(
            f"Residue sequence mismatch against UniProt P13569 ({seq_id_column}): {mismatches[:5]}"
        )


def parse_8eiq_structure(
    cif_path: Path = DEFAULT_CIF_8EIQ,
    uniprot_seq: str | None = None,
    *,
    seq_id_column: str = "auth_seq_id",
    pot_hetatm_ids: frozenset[str] = POTENTIATION_HETATM_IDS,
    fold_hetatm_ids: frozenset[str] = FOLDING_HETATM_IDS,
    fold_socket_residues: frozenset[int] = FOLDING_SOCKET_RESIDUES,
    verify_hash: bool = True,
    rigid_transform: Callable[[tuple[float, float, float]], tuple[float, float, float]] | None = None,
) -> Structure8EIQ:
    if verify_hash:
        actual = sha256_file(cif_path)
        if actual != EXPECTED_8EIQ_SHA256:
            raise CryoEMPotentiationError(
                f"8EIQ.cif.gz SHA-256 mismatch: expected {EXPECTED_8EIQ_SHA256}, got {actual}"
            )
    if pot_hetatm_ids & fold_hetatm_ids:
        raise CryoEMPotentiationError(
            f"Overlapping HETATM anchor sets: {sorted(pot_hetatm_ids & fold_hetatm_ids)}"
        )

    lines = gzip.open(cif_path, "rt", encoding="utf-8").read().splitlines()

    diff_cols: list[str] = []
    i = 0
    while i < len(lines) and not lines[i].startswith("_struct_ref_seq_dif."):
        i += 1
    while i < len(lines) and lines[i].startswith("_struct_ref_seq_dif."):
        diff_cols.append(lines[i].strip()[len("_struct_ref_seq_dif."):])
        i += 1
    declared_mutations: dict[int, str] = {}
    declared_deletions: set[int] = set()
    while diff_cols and i < len(lines) and not lines[i].startswith(("#", "_", "loop_")):
        lexemes = shlex.split(lines[i])
        i += 1
        if len(lexemes) == len(diff_cols):
            row = dict(zip(diff_cols, lexemes))
            if row.get("pdbx_pdb_strand_id") != "A" or row.get("pdbx_seq_db_accession_code") != "P13569":
                continue
            if row.get("details") == "expression tag":
                continue
            pos_str = row.get("pdbx_seq_db_seq_num", "")
            if pos_str.isdigit():
                p_int = int(pos_str)
                mon = row.get("mon_id", "?")
                if mon in ("?", "."):
                    declared_deletions.add(p_int)
                else:
                    declared_mutations[p_int] = mon

    atom_cols: list[str] = []
    atom_rows: list[list[str]] = []
    started = False
    for ln in lines:
        if ln.startswith("_atom_site."):
            atom_cols.append(ln.strip()[len("_atom_site."):])
            started = True
        elif started and ln.startswith(("ATOM ", "HETATM ")):
            atom_rows.append(shlex.split(ln))
        elif started and atom_rows and not ln.startswith(("ATOM ", "HETATM ")):
            break

    ix = {c: idx for idx, c in enumerate(atom_cols)}
    if seq_id_column not in ix:
        raise CryoEMPotentiationError(f"Missing sequence column {seq_id_column!r} in 8EIQ mmCIF")

    ca_by_seq: dict[int, tuple[str, tuple[float, float, float]]] = {}
    pot_atoms: list[HeavyAtom] = []
    fold_atoms: list[HeavyAtom] = []
    comp_counts: dict[str, int] = {k: 0 for k in (pot_hetatm_ids | fold_hetatm_ids)}
    socket_counts: dict[int, int] = {r: 0 for r in fold_socket_residues}

    for r in atom_rows:
        if r[ix["pdbx_PDB_model_num"]] != "1":
            continue
        if r[ix["auth_asym_id"]] != "A":
            continue
        if r[ix["label_alt_id"]] not in {".", "A"}:
            continue
        if r[ix["type_symbol"]].upper() == "H":
            continue

        xyz = (
            float(r[ix["Cartn_x"]]),
            float(r[ix["Cartn_y"]]),
            float(r[ix["Cartn_z"]]),
        )
        if rigid_transform is not None:
            xyz = rigid_transform(xyz)

        group = r[ix["group_PDB"]]
        comp_id = r[ix["label_comp_id"]]
        atom_id = r[ix["label_atom_id"]]
        raw_seq = r[ix[seq_id_column]]
        seq_id = int(raw_seq) if raw_seq.lstrip("-").isdigit() else None

        if group == "HETATM":
            ha = HeavyAtom(group=group, comp_id=comp_id, atom_id=atom_id, seq_id=seq_id, xyz=xyz)
            if comp_id in pot_hetatm_ids:
                pot_atoms.append(ha)
                comp_counts[comp_id] = comp_counts.get(comp_id, 0) + 1
            if comp_id in fold_hetatm_ids:
                fold_atoms.append(ha)
                comp_counts[comp_id] = comp_counts.get(comp_id, 0) + 1
        elif group == "ATOM":
            if seq_id is None:
                continue
            if atom_id == "CA" and seq_id not in ca_by_seq:
                ca_by_seq[seq_id] = (comp_id, xyz)
            if seq_id in fold_socket_residues:
                fold_atoms.append(
                    HeavyAtom(group=group, comp_id=comp_id, atom_id=atom_id, seq_id=seq_id, xyz=xyz)
                )
                socket_counts[seq_id] = socket_counts.get(seq_id, 0) + 1

    pot_coords = {a.xyz for a in pot_atoms}
    fold_coords = {a.xyz for a in fold_atoms}
    if pot_coords & fold_coords:
        raise CryoEMPotentiationError("Anchor coordinate sets A_pot and A_fold intersect")

    if uniprot_seq is not None:
        verify_ca_against_uniprot(
            ca_by_seq,
            uniprot_seq,
            declared_mutations,
            frozenset(declared_deletions),
            seq_id_column=seq_id_column,
        )

    return Structure8EIQ(
        ca_by_seq=ca_by_seq,
        potentiation_atoms=tuple(pot_atoms),
        folding_atoms=tuple(fold_atoms),
        comp_counts=comp_counts,
        socket_counts=socket_counts,
        declared_mutations=declared_mutations,
        declared_deletions=frozenset(declared_deletions),
    )


def compute_residue_s3d(
    codon: int,
    struct: Structure8EIQ,
    pot_grid: VoxelGrid3D,
    fold_grid: VoxelGrid3D,
) -> dict[str, float] | None:
    if codon not in struct.ca_by_seq:
        return None
    _, ca_xyz = struct.ca_by_seq[codon]
    d_pot = pot_grid.min_distance(ca_xyz)
    d_fold = fold_grid.min_distance(ca_xyz)
    denom = d_pot + d_fold
    if denom <= 0.0 or not math.isfinite(denom):
        raise CryoEMPotentiationError(f"Zero or non-finite d_pot + d_fold at codon {codon}")
    return {
        "d_pot_A": d_pot,
        "d_fold_A": d_fold,
        "s_3d": d_fold / denom,
    }


def sequence_consequence_score(variant: dict[str, Any]) -> float:
    fam = variant["variant_family"]
    hgvs = str(variant["hgvs_p_or_c"])
    if fam == "class_i_ptc_frameshift":
        return 0.0
    if fam == "splice_regulatory":
        return 0.0 if CANONICAL_SPLICE_RE.match(hgvs) else 0.50
    if fam == "inframe_indel":
        return 0.50
    if fam == "missense":
        return 0.75
    raise CryoEMPotentiationError(f"Unknown variant_family: {fam!r}")


def gated_score_from_missense_unit(variant: dict[str, Any], missense_unit_score: float | None) -> float | None:
    fam = variant["variant_family"]
    if fam != "missense":
        if missense_unit_score is not None:
            raise CryoEMPotentiationError(
                f"Non-missense variant {variant.get('variant_id', fam)} cannot accept a missense unit score"
            )
        return sequence_consequence_score(variant)
    if missense_unit_score is None:
        return None
    if isinstance(missense_unit_score, bool) or not math.isfinite(missense_unit_score) or not (0.0 <= missense_unit_score <= 1.0):
        raise CryoEMPotentiationError(f"Missense unit score {missense_unit_score} outside finite [0, 1]")
    return 0.50 + 0.50 * missense_unit_score


def compute_auroc_metrics(
    labels: list[int],
    scores: list[float | None],
) -> dict[str, Any]:
    """Direct pairwise Mann-Whitney ROC-AUC with 0.5 tie credit on unscored pairs (zero sign-flipping)."""
    if len(labels) != len(scores) or not labels:
        raise CryoEMPotentiationError("Labels and scores must be non-empty and equal length")
    for y in labels:
        if isinstance(y, bool) or y not in (0, 1):
            raise CryoEMPotentiationError(f"Labels must be binary integers in {{0, 1}}, got {y!r}")
    for s in scores:
        if s is not None and (isinstance(s, bool) or not math.isfinite(s)):
            raise CryoEMPotentiationError(f"Scores must be None or finite real numbers, got {s!r}")

    positives = [i for i, y in enumerate(labels) if y == 1]
    negatives = [i for i, y in enumerate(labels) if y == 0]
    if not positives or not negatives:
        raise CryoEMPotentiationError("Both positive and negative labels are required for ROC-AUC")

    scored_pos = [scores[i] for i in positives if scores[i] is not None]
    scored_neg = [scores[i] for i in negatives if scores[i] is not None]

    cond_auc: float | None = None
    if scored_pos and scored_neg:
        cond_sum = 0.0
        for sp in scored_pos:
            assert sp is not None
            for sn in scored_neg:
                assert sn is not None
                if sp > sn:
                    cond_sum += 1.0
                elif sp == sn:
                    cond_sum += 0.5
        cond_auc = cond_sum / (len(scored_pos) * len(scored_neg))

    abst_sum = 0.0
    for i in positives:
        si = scores[i]
        for j in negatives:
            sj = scores[j]
            if si is None or sj is None:
                abst_sum += 0.5
            elif si > sj:
                abst_sum += 1.0
            elif si == sj:
                abst_sum += 0.5
    abst_auc = abst_sum / (len(positives) * len(negatives))
    scored_count = sum(1 for s in scores if s is not None)
    return {
        "n_total": len(labels),
        "n_scored": scored_count,
        "coverage": scored_count / len(labels),
        "n_pos": len(positives),
        "n_neg": len(negatives),
        "n_scored_pos": len(scored_pos),
        "n_scored_neg": len(scored_neg),
        "conditional_auroc": cond_auc,
        "abstention_aware_auroc": abst_auc,
    }


def validate_variants_and_sources(
    variants: list[dict[str, Any]],
    catalog_sources: Iterable[str] = (),
    catalog_path: Path = DEFAULT_CATALOG_YAML,
) -> None:
    src_list = list(catalog_sources)
    if src_list:
        admitted = load_admitted_automated_catalog_ids(catalog_path)
        disallowed = sorted(set(src_list) - admitted)
        if disallowed:
            raise PermissionError(
                f"Source(s) not admitted for automated fetch/api in sources/catalog.yaml: {disallowed}"
            )

    train_codons: set[int] = set()
    test_codons: set[int] = set()
    train_domains: set[str] = set()
    test_domains: set[str] = set()

    for v in variants:
        fam = str(v.get("variant_family", ""))
        if fam not in ALLOWED_FAMILIES:
            raise CryoEMPotentiationError(f"Invalid variant_family {fam!r} on {v.get('variant_id')}")
        if v.get("evo2_zero_shot_score") is not None or v.get("alphagenome_avi_score") is not None:
            raise CryoEMPotentiationError(
                f"Uncataloged zero-shot model score present on {v.get('variant_id')}"
            )
        if fam != "missense":
            for miss_only_key in ("s_3d", "s_3d_override", "alphamissense_score", "revel_score"):
                if v.get(miss_only_key) is not None:
                    raise CryoEMPotentiationError(
                        f"Non-missense variant {v.get('variant_id')} cannot carry non-None {miss_only_key}"
                    )

        codon = int(v["codon_position"])
        dom = map_p13569_domain(codon)
        split = v.get("split")
        if split == "train":
            if dom not in TRAIN_DOMAINS:
                raise CryoEMPotentiationError(f"Train variant {v['variant_id']} in non-train domain {dom}")
            if codon in train_codons:
                raise CryoEMPotentiationError(f"Duplicate train codon {codon}")
            train_codons.add(codon)
            train_domains.add(dom)
        elif split == "test":
            if dom not in TEST_DOMAINS:
                raise CryoEMPotentiationError(f"Test variant {v['variant_id']} in non-test domain {dom}")
            if codon in test_codons:
                raise CryoEMPotentiationError(f"Duplicate test codon {codon}")
            test_codons.add(codon)
            test_domains.add(dom)
        else:
            raise CryoEMPotentiationError(f"Unknown split {split!r}")

    if train_codons & test_codons:
        raise CryoEMPotentiationError(f"Train/test codon overlap: {sorted(train_codons & test_codons)}")
    if train_domains & test_domains:
        raise CryoEMPotentiationError(f"Train/test domain overlap: {sorted(train_domains & test_domains)}")
    min_gap = min(abs(tc - xc) for tc in train_codons for xc in test_codons)
    if min_gap <= 5:
        raise CryoEMPotentiationError(f"Cross-split minimum codon gap {min_gap} <= 5")


def run_positive_controls(
    cif_path: Path = DEFAULT_CIF_8EIQ,
    uniprot_path: Path = DEFAULT_UNIPROT_JSON,
    benchmark_codons: frozenset[int] = frozenset(),
) -> dict[str, Any]:
    uniprot_seq = load_uniprot_sequence(uniprot_path, verify_hash=True)
    struct = parse_8eiq_structure(cif_path, uniprot_seq=uniprot_seq, verify_hash=True)

    expected_comp = {"ATP": 62, "MG": 2, "VX7": 29, "CV6": 37, "WJX": 41}
    if struct.comp_counts != expected_comp:
        raise CryoEMPotentiationError(
            f"8EIQ ligand heavy-atom count mismatch: expected {expected_comp}, got {struct.comp_counts}"
        )
    expected_socket = {1063: 14, 1070: 11}
    if struct.socket_counts != expected_socket:
        raise CryoEMPotentiationError(
            f"8EIQ socket heavy-atom count mismatch: expected {expected_socket}, got {struct.socket_counts}"
        )
    if len(struct.potentiation_atoms) != 93 or len(struct.folding_atoms) != 103:
        raise CryoEMPotentiationError("Expected |A_pot|=93 and |A_fold|=103")
    if 508 in struct.ca_by_seq:
        raise CryoEMPotentiationError("Residue 508 must be absent in deltaF508 8EIQ construct")

    pot_grid = VoxelGrid3D(a.xyz for a in struct.potentiation_atoms)
    fold_grid = VoxelGrid3D(a.xyz for a in struct.folding_atoms)

    vx7_pts = [a.xyz for a in struct.potentiation_atoms if a.comp_id == "VX7"]
    atp_pts = [a.xyz for a in struct.potentiation_atoms if a.comp_id == "ATP"]
    cv6_pts = [a.xyz for a in struct.folding_atoms if a.comp_id == "CV6"]
    wjx_pts = [a.xyz for a in struct.folding_atoms if a.comp_id == "WJX"]

    pocket_residues = (304, 307, 311, 931, 464, 1250, 361, 1102)
    if set(pocket_residues) & set(benchmark_codons):
        raise CryoEMPotentiationError("Positive-control pocket residue overlaps benchmark split")

    for r_vx7 in (304, 307, 311, 931):
        ca_xyz = struct.ca_by_seq[r_vx7][1]
        d_vx7 = min(math.dist(ca_xyz, q) for q in vx7_pts)
        m = compute_residue_s3d(r_vx7, struct, pot_grid, fold_grid)
        if m is None or not (d_vx7 < 8.0 and m["s_3d"] > 0.70):
            raise CryoEMPotentiationError(f"VX7 pocket positive control failed at {r_vx7}: d={d_vx7}, m={m}")

    for r_atp in (464, 1250):
        ca_xyz = struct.ca_by_seq[r_atp][1]
        d_atp = min(math.dist(ca_xyz, q) for q in atp_pts)
        m = compute_residue_s3d(r_atp, struct, pot_grid, fold_grid)
        if m is None or not (d_atp < 4.0 and m["s_3d"] > 0.75):
            raise CryoEMPotentiationError(f"ATP Walker-A positive control failed at {r_atp}: d={d_atp}, m={m}")

    d_w361 = min(math.dist(struct.ca_by_seq[361][1], q) for q in cv6_pts)
    m_w361 = compute_residue_s3d(361, struct, pot_grid, fold_grid)
    if m_w361 is None or not (d_w361 < 5.0 and m_w361["s_3d"] < 0.20):
        raise CryoEMPotentiationError(f"CV6 pocket positive control failed at W361: d={d_w361}, m={m_w361}")

    d_r1102 = min(math.dist(struct.ca_by_seq[1102][1], q) for q in wjx_pts)
    m_r1102 = compute_residue_s3d(1102, struct, pot_grid, fold_grid)
    if m_r1102 is None or not (d_r1102 < 6.0 and m_r1102["s_3d"] < 0.20):
        raise CryoEMPotentiationError(f"WJX pocket positive control failed at R1102: d={d_r1102}, m={m_r1102}")

    ax, ay, az = 0.61, -0.43, 0.89
    cx, sx = math.cos(ax), math.sin(ax)
    cy, sy = math.cos(ay), math.sin(ay)
    cz, sz = math.cos(az), math.sin(az)
    rot = (
        (cy * cz, cz * sx * sy - cx * sz, cx * cz * sy + sx * sz),
        (cy * sz, cx * cz + sx * sy * sz, -cz * sx + cx * sy * sz),
        (-sy, cy * sx, cx * cy),
    )
    shift = (137.25, -84.5, 42.125)

    def apply_se3(p: tuple[float, float, float]) -> tuple[float, float, float]:
        return (
            rot[0][0] * p[0] + rot[0][1] * p[1] + rot[0][2] * p[2] + shift[0],
            rot[1][0] * p[0] + rot[1][1] * p[1] + rot[1][2] * p[2] + shift[1],
            rot[2][0] * p[0] + rot[2][1] * p[1] + rot[2][2] * p[2] + shift[2],
        )

    moved = parse_8eiq_structure(cif_path, uniprot_seq=uniprot_seq, verify_hash=False, rigid_transform=apply_se3)
    moved_pot = VoxelGrid3D(a.xyz for a in moved.potentiation_atoms)
    moved_fold = VoxelGrid3D(a.xyz for a in moved.folding_atoms)
    max_se3_drift = 0.0
    for seq_id in struct.ca_by_seq:
        r0 = compute_residue_s3d(seq_id, struct, pot_grid, fold_grid)
        r1 = compute_residue_s3d(seq_id, moved, moved_pot, moved_fold)
        assert r0 is not None and r1 is not None
        drift = abs(r0["s_3d"] - r1["s_3d"])
        if drift > max_se3_drift:
            max_se3_drift = drift
    if max_se3_drift >= 1e-12:
        raise CryoEMPotentiationError(f"SE(3) invariance drift {max_se3_drift:.3e} >= 1e-12")

    m_g178 = compute_residue_s3d(178, struct, pot_grid, fold_grid)
    m_g85 = compute_residue_s3d(85, struct, pot_grid, fold_grid)
    if m_g178 is None or m_g85 is None or not (m_g178["s_3d"] > m_g85["s_3d"]):
        raise CryoEMPotentiationError(f"Train control G178 > G85 failed: G178={m_g178}, G85={m_g85}")

    auc_none = compute_auroc_metrics([1, 1, 0, 0], [None, None, None, None])["abstention_aware_auroc"]
    auc_pos = compute_auroc_metrics([1, 1, 0, 0], [0.9, 0.8, 0.2, 0.1])["abstention_aware_auroc"]
    auc_inv = compute_auroc_metrics([1, 1, 0, 0], [0.1, 0.2, 0.8, 0.9])["abstention_aware_auroc"]
    if auc_none != 0.5 or auc_pos != 1.0 or auc_inv != 0.0:
        raise CryoEMPotentiationError(
            f"ROC-AUC unit controls failed: none={auc_none}, pos={auc_pos}, inv={auc_inv}"
        )

    return {
        "passed": True,
        "placed_chain_a_residues": len(struct.ca_by_seq),
        "potentiation_heavy_atoms": len(struct.potentiation_atoms),
        "folding_heavy_atoms": len(struct.folding_atoms),
        "pocket_controls_verified": len(pocket_residues),
        "max_se3_drift": max_se3_drift,
        "train_g178_s3d": round(m_g178["s_3d"], 6),
        "train_g85_s3d": round(m_g85["s_3d"], 6),
    }


def run_falsification_suite(
    cif_path: Path = DEFAULT_CIF_8EIQ,
    uniprot_path: Path = DEFAULT_UNIPROT_JSON,
    annotations_path: Path = DEFAULT_ANNOTATIONS_JSON,
) -> dict[str, bool]:
    uniprot_seq = load_uniprot_sequence(uniprot_path, verify_hash=False)
    annot = json.loads(annotations_path.read_text(encoding="utf-8"))
    base_variants = annot["variants"]
    results: dict[str, bool] = {}

    # 1. label_seq_id_shift_rejected
    try:
        parse_8eiq_structure(cif_path, uniprot_seq=uniprot_seq, seq_id_column="label_seq_id", verify_hash=False)
        results["label_seq_id_shift_rejected"] = False
    except CryoEMPotentiationError:
        results["label_seq_id_shift_rejected"] = True

    # 2. overlapping_anchor_sets_rejected
    try:
        parse_8eiq_structure(
            cif_path,
            uniprot_seq=uniprot_seq,
            pot_hetatm_ids=frozenset({"ATP", "VX7"}),
            fold_hetatm_ids=frozenset({"VX7", "CV6"}),
            verify_hash=False,
        )
        results["overlapping_anchor_sets_rejected"] = False
    except CryoEMPotentiationError:
        results["overlapping_anchor_sets_rejected"] = True

    # 3. non_missense_3d_score_rejected
    poisoned_non_missense = [dict(v) for v in base_variants]
    for pv in poisoned_non_missense:
        if pv["legacy_spl_name"] == "F508del":
            pv["s_3d"] = 0.42
            break
    try:
        validate_variants_and_sources(poisoned_non_missense)
        results["non_missense_3d_score_rejected"] = False
    except CryoEMPotentiationError:
        results["non_missense_3d_score_rejected"] = True

    # 4. post_hoc_sign_flip_rejected
    inv_metrics = compute_auroc_metrics([1, 1, 0, 0], [0.1, 0.2, 0.8, 0.9])
    results["post_hoc_sign_flip_rejected"] = (
        inv_metrics["conditional_auroc"] == 0.0 and inv_metrics["abstention_aware_auroc"] == 0.0
    )

    # 5. codon_or_domain_leakage_rejected
    poisoned_gap = [dict(v) for v in base_variants]
    poisoned_gap[0] = dict(poisoned_gap[0], codon_position=452)
    try:
        validate_variants_and_sources(poisoned_gap)
        results["codon_or_domain_leakage_rejected"] = False
    except CryoEMPotentiationError:
        results["codon_or_domain_leakage_rejected"] = True

    # 6. unadmitted_manual_source_rejected
    try:
        validate_variants_and_sources(base_variants, catalog_sources=["openfda-label", "mavedb-dump-zenodo"])
        results["unadmitted_manual_source_rejected"] = False
    except PermissionError:
        results["unadmitted_manual_source_rejected"] = True

    if not all(results.values()):
        failed = [k for k, ok in results.items() if not ok]
        raise CryoEMPotentiationError(f"Falsification negative control(s) failed to trigger: {failed}")
    return results


def evaluate_benchmark(
    cif_path: Path = DEFAULT_CIF_8EIQ,
    uniprot_path: Path = DEFAULT_UNIPROT_JSON,
    annotations_path: Path = DEFAULT_ANNOTATIONS_JSON,
) -> dict[str, Any]:
    annot = json.loads(annotations_path.read_text(encoding="utf-8"))
    variants = annot["variants"]
    sources = annot["provenance"]["catalog_sources"]
    validate_variants_and_sources(variants, catalog_sources=sources)

    all_codons = frozenset(int(v["codon_position"]) for v in variants)
    pos_controls = run_positive_controls(cif_path, uniprot_path, benchmark_codons=all_codons)

    uniprot_seq = load_uniprot_sequence(uniprot_path, verify_hash=False)
    struct = parse_8eiq_structure(cif_path, uniprot_seq=uniprot_seq, verify_hash=False)
    pot_grid = VoxelGrid3D(a.xyz for a in struct.potentiation_atoms)
    fold_grid = VoxelGrid3D(a.xyz for a in struct.folding_atoms)

    split_reports: dict[str, Any] = {}
    variant_rows: list[dict[str, Any]] = []

    for split_name in ("train", "test"):
        sub = [v for v in variants if v["split"] == split_name]
        labels = [int(v["kalydeco_v36_responsive"]) for v in sub]

        s3d_scores: list[float | None] = []
        sam_scores: list[float | None] = []
        seq_scores: list[float | None] = []
        g3d_scores: list[float | None] = []
        gam_scores: list[float | None] = []

        for v in sub:
            codon = int(v["codon_position"])
            fam = str(v["variant_family"])
            am = v.get("alphamissense_score")

            s_am = (1.0 - float(am)) if (fam == "missense" and am is not None) else None
            geom = compute_residue_s3d(codon, struct, pot_grid, fold_grid) if fam == "missense" else None
            s_3d = geom["s_3d"] if geom is not None else None

            seq_s = sequence_consequence_score(v)
            g_3d = gated_score_from_missense_unit(v, s_3d)
            g_am = gated_score_from_missense_unit(v, s_am)

            s3d_scores.append(s_3d)
            sam_scores.append(s_am)
            seq_scores.append(seq_s)
            g3d_scores.append(g_3d)
            gam_scores.append(g_am)

            variant_rows.append(
                {
                    "variant_id": v["variant_id"],
                    "legacy_spl_name": v["legacy_spl_name"],
                    "split": split_name,
                    "domain": map_p13569_domain(codon),
                    "codon_position": codon,
                    "variant_family": fam,
                    "kalydeco_v36_responsive": int(v["kalydeco_v36_responsive"]),
                    "trikafta_v23_responsive": int(v["trikafta_v23_responsive"]),
                    "alyftrek_v9_responsive": int(v["alyftrek_v9_responsive"]),
                    "d_pot_A": round(geom["d_pot_A"], 6) if geom else None,
                    "d_fold_A": round(geom["d_fold_A"], 6) if geom else None,
                    "s_3d": round(s_3d, 6) if s_3d is not None else None,
                    "s_am": round(s_am, 6) if s_am is not None else None,
                    "seq_baseline": seq_s,
                    "g_3d": round(g_3d, 6) if g_3d is not None else None,
                    "g_am": round(g_am, 6) if g_am is not None else None,
                }
            )

        missense_idx = [idx for idx, v in enumerate(sub) if v["variant_family"] == "missense"]
        m_labels = [labels[idx] for idx in missense_idx]
        m_s3d = [s3d_scores[idx] for idx in missense_idx]
        m_sam = [sam_scores[idx] for idx in missense_idx]
        m_seq = [seq_scores[idx] for idx in missense_idx]

        split_reports[split_name] = {
            "n_variants": len(sub),
            "n_kalydeco_responsive": sum(labels),
            "n_kalydeco_nonresponsive": len(labels) - sum(labels),
            "missense_only": {
                "s_3d": compute_auroc_metrics(m_labels, m_s3d),
                "s_am": compute_auroc_metrics(m_labels, m_sam),
                "seq_baseline": compute_auroc_metrics(m_labels, m_seq),
            },
            "all_variants": {
                "s_3d_ungated": compute_auroc_metrics(labels, s3d_scores),
                "s_am_ungated": compute_auroc_metrics(labels, sam_scores),
                "seq_baseline": compute_auroc_metrics(labels, seq_scores),
                "g_3d_gated": compute_auroc_metrics(labels, g3d_scores),
                "g_am_gated": compute_auroc_metrics(labels, gam_scores),
            },
        }

    test_m_s3d = split_reports["test"]["missense_only"]["s_3d"]["conditional_auroc"]
    test_m_sam = split_reports["test"]["missense_only"]["s_am"]["conditional_auroc"]
    test_g3d = split_reports["test"]["all_variants"]["g_3d_gated"]["abstention_aware_auroc"]
    test_gam = split_reports["test"]["all_variants"]["g_am_gated"]["abstention_aware_auroc"]
    test_seq = split_reports["test"]["all_variants"]["seq_baseline"]["abstention_aware_auroc"]

    cond_ok = (
        test_m_s3d is not None
        and test_m_sam is not None
        and test_m_s3d >= 0.75
        and test_m_s3d > test_m_sam
    )
    gated_ok = (
        test_g3d >= 0.75
        and test_g3d > test_seq
        and test_g3d > test_gam
    )
    outcome = "supported" if (cond_ok and gated_ok) else "refuted"

    return {
        "prediction_id": "PRED-0011",
        "split_version": annot["schema_version"],
        "target_label": "kalydeco_v36_responsive (US FDA KALYDECO v36 Table 3, NDA 203188)",
        "structure_source": "PDB 8EIQ (tools/cftr-model/data/8EIQ.cif.gz, CC0 1.0)",
        "positive_controls": pos_controls,
        "splits": split_reports,
        "pred_0011_evaluation": {
            "test_missense_conditional_auroc_s3d": round(float(test_m_s3d), 6),
            "test_missense_conditional_auroc_sam": round(float(test_m_sam), 6),
            "test_all_abstention_auroc_g3d": round(float(test_g3d), 6),
            "test_all_abstention_auroc_gam": round(float(test_gam), 6),
            "test_all_abstention_auroc_seq_baseline": round(float(test_seq), 6),
            "held_out_error_1_minus_g3d_auroc": round(1.0 - float(test_g3d), 6),
            "held_out_missense_error_1_minus_s3d_auroc": round(1.0 - float(test_m_s3d), 6),
            "criterion_a_missense_passed": cond_ok,
            "criterion_b_gated_passed": gated_ok,
            "outcome": outcome,
        },
        "variants": variant_rows,
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Cryo-EM 3D potentiation-vs-folding pocket benchmark (PRED-0011)")
    ap.add_argument("--check", action="store_true", help="Run positive controls and full benchmark evaluation")
    ap.add_argument("--falsify", action="store_true", help="Run all 6 planted negative controls")
    ap.add_argument("--json", action="store_true", help="Emit full JSON evaluation payload")
    args = ap.parse_args(argv)

    falsify_res = run_falsification_suite() if args.falsify else None
    report = evaluate_benchmark()
    if falsify_res is not None:
        report["falsification_controls"] = falsify_res

    if args.json:
        print(json.dumps(report, indent=2))
        return 0

    ev = report["pred_0011_evaluation"]
    print(
        f"PASS cryoem_potentiation_benchmark: outcome={ev['outcome']} "
        f"(test missense AUC_cond s_3D={ev['test_missense_conditional_auroc_s3d']:.4f} vs s_AM={ev['test_missense_conditional_auroc_sam']:.4f}; "
        f"all-20 AUC_abst g_3D={ev['test_all_abstention_auroc_g3d']:.4f} vs seq={ev['test_all_abstention_auroc_seq_baseline']:.4f}, g_AM={ev['test_all_abstention_auroc_gam']:.4f})"
    )
    if falsify_res is not None:
        print(f"PASS falsify controls ({sum(falsify_res.values())}/{len(falsify_res)} rejected)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
