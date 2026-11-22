#!/usr/bin/env python3
"""Modulator-Free Multi-State Cryo-EM Gating Kinematics Benchmark (PRED-0012).

Tests whether SE(3) x SE(3)-invariant modulator-free apo-to-ATP conformational
reorganization
    s_kin(i) = 1 - lDDT_sym(i; 5UAK, 6MSM)
computed strictly from polypeptide Chain A C-alpha coordinates in dephosphorylated
ATP-free WT human CFTR (5UAK, CC0 1.0) and phosphorylated ATP/Mg-bound E1371Q human
CFTR (6MSM, CC0 1.0) using Mariani et al. (2013) published defaults (R_0 = 15.0 A,
thresholds {0.5, 1.0, 2.0, 4.0} A, |i - j| > 1) predicts US FDA KALYDECO v36 Table 3
ivacaftor monotherapy responsiveness on the fixed 20-variant domain-disjoint held-out
test split (2026-11-21.v2: NBD1, TMD2, NBD2_C_terminal).

Method 1 implementation details (disjoint from independent_audit.py):
- Parses 5UAK.cif.gz and 6MSM.cif.gz via shlex-tokenized mmCIF loops/tags and typed dataclasses.
- Indexes shared Chain A C-alpha coordinates in 15.0 A 3D cubic spatial voxel hash grids
  to assemble symmetrized local inclusion neighborhoods N_sym(i).
- Evaluates conditional and abstention-aware ROC-AUC via O(N1*N0) pairwise Mann-Whitney
  concordance enumeration (forbidding post-hoc sign flipping).

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
from typing import Any, Callable

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
DEFAULT_CIF_5UAK = REPO_ROOT / "tools" / "cftr-model" / "data" / "5UAK.cif.gz"
DEFAULT_CIF_6MSM = REPO_ROOT / "tools" / "cftr-model" / "data" / "6MSM.cif.gz"
DEFAULT_UNIPROT_JSON = REPO_ROOT / "tools" / "cftr-model" / "data" / "p13569.json"
DEFAULT_CATALOG_YAML = REPO_ROOT / "sources" / "catalog.yaml"
DEFAULT_ANNOTATIONS_JSON = (
    REPO_ROOT / "benchmarks" / "theratype" / "theratype_open_annotations_2026-11-21.json"
)

EXPECTED_5UAK_SHA256 = "8bf00f2612dd961d3b1ccd2f2193781a534ad44f32666a6de8da948afb8d84cf"
EXPECTED_6MSM_SHA256 = "82933988abc3ed24a343819d410b39dc0d86913eaba4426dd945e3273e34fa59"
EXPECTED_P13569_SHA256 = "c8feabfabe8690ef4a06af7de5b5f8dc3c8da554246c5ab3d403edf265d4a893"

LDDT_INCLUSION_RADIUS_A: float = 15.0
LDDT_THRESHOLDS_A: tuple[float, float, float, float] = (0.5, 1.0, 2.0, 4.0)
BOUNDARY_EPSILON_A: float = 1e-9

THREE_TO_ONE: dict[str, str] = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

TRAIN_DOMAINS: frozenset[str] = frozenset({"Lasso_N_terminal", "TMD1", "R_domain"})
TEST_DOMAINS: frozenset[str] = frozenset({"NBD1", "TMD2", "NBD2_C_terminal"})
ALLOWED_FAMILIES: frozenset[str] = frozenset(
    {"missense", "inframe_indel", "splice_regulatory", "class_i_ptc_frameshift"}
)
CANONICAL_SPLICE_RE = re.compile(r"^c\.\d+[+-][12](?:[^\d]|$)")
CATALOG_ID_RE = re.compile(r"^\s*-\s*id:\s*([^\s#]+)")
CATALOG_ACCESS_RE = re.compile(r"^\s*(?:access|mode):\s*([^\s#]+)")
CATALOG_LABEL_RE = re.compile(r"^\s*claim_label:\s*([^\s#]+)")

DIMER_CLOSURE_PAIRS: tuple[tuple[int, int], ...] = ((464, 1348), (582, 1246))
GATING_HINGE_CONTROLS: tuple[int, ...] = (267, 1347, 1348, 1371)
RIGID_BETA_CONTROLS: tuple[int, ...] = (361, 450, 1102)


class CryoEMGatingKinematicsError(ValueError):
    """Raised when structural provenance, positive controls, or negative controls fail."""


@dataclasses.dataclass(frozen=True)
class ChainAStructure:
    pdb_id: str
    ca_by_seq: dict[int, tuple[str, tuple[float, float, float]]]
    declared_mutations: dict[int, str]
    declared_deletions: frozenset[int]


class SpatialGrid15A:
    """3D cubic voxel hash grid for exact C-alpha neighborhood queries within cell_size A."""

    def __init__(
        self,
        coords_by_seq: dict[int, tuple[float, float, float]],
        cell_size: float = LDDT_INCLUSION_RADIUS_A,
    ) -> None:
        if not coords_by_seq:
            raise CryoEMGatingKinematicsError("SpatialGrid15A initialized with empty residue map")
        self.coords_by_seq = dict(coords_by_seq)
        self.cell_size = float(cell_size)
        buckets: dict[tuple[int, int, int], list[tuple[int, tuple[float, float, float]]]] = {}
        for seq_id, pt in self.coords_by_seq.items():
            buckets.setdefault(self._cell(pt), []).append((seq_id, pt))
        self.buckets = buckets

    def _cell(self, pt: tuple[float, float, float]) -> tuple[int, int, int]:
        cs = self.cell_size
        return (math.floor(pt[0] / cs), math.floor(pt[1] / cs), math.floor(pt[2] / cs))

    def neighbors_within(self, query_seq: int, radius: float = LDDT_INCLUSION_RADIUS_A) -> set[int]:
        pt_q = self.coords_by_seq[query_seq]
        qx, qy, qz = self._cell(pt_q)
        span = max(1, math.ceil(radius / self.cell_size))
        hits: set[int] = set()
        for dx in range(-span, span + 1):
            for dy in range(-span, span + 1):
                for dz in range(-span, span + 1):
                    for seq_j, pt_j in self.buckets.get((qx + dx, qy + dy, qz + dz), ()):
                        if abs(seq_j - query_seq) <= 1:
                            continue
                        d = math.dist(pt_q, pt_j)
                        if abs(d - radius) <= BOUNDARY_EPSILON_A:
                            raise CryoEMGatingKinematicsError(
                                f"Boundary-knife-edge pair ({query_seq}, {seq_j}): d={d:.12f} vs radius={radius}"
                            )
                        if d <= radius:
                            hits.add(seq_j)
        return hits


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


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
    raise CryoEMGatingKinematicsError(f"Codon {codon} outside human CFTR 1..1480")


def load_admitted_automated_catalog_ids(catalog_path: Path = DEFAULT_CATALOG_YAML) -> frozenset[str]:
    if not catalog_path.is_file():
        raise CryoEMGatingKinematicsError(f"Catalog YAML not found: {catalog_path}")
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


def assert_catalog_sources_admitted(
    source_ids: list[str] | tuple[str, ...],
    admitted_ids: frozenset[str],
) -> None:
    for sid in source_ids:
        if sid not in admitted_ids:
            raise PermissionError(
                f"Source {sid!r} is not admitted as automated fetch/api in sources/catalog.yaml"
            )


def load_uniprot_sequence(uniprot_path: Path = DEFAULT_UNIPROT_JSON, verify_hash: bool = True) -> str:
    if verify_hash:
        actual = sha256_file(uniprot_path)
        if actual != EXPECTED_P13569_SHA256:
            raise CryoEMGatingKinematicsError(
                f"p13569.json SHA-256 mismatch: expected {EXPECTED_P13569_SHA256}, got {actual}"
            )
    payload = json.loads(uniprot_path.read_text(encoding="utf-8"))
    seq = str(payload["sequence"]["value"])
    if len(seq) != 1480:
        raise CryoEMGatingKinematicsError(f"Expected 1480-aa P13569 sequence, got {len(seq)}")
    return seq


def verify_ca_against_uniprot(
    ca_by_seq: dict[int, tuple[str, tuple[float, float, float]]],
    uniprot_seq: str,
    declared_mutations: dict[int, str],
    declared_deletions: frozenset[int],
    pdb_id: str,
) -> None:
    for del_pos in declared_deletions:
        if del_pos in ca_by_seq:
            raise CryoEMGatingKinematicsError(
                f"{pdb_id}: deleted construct residue {del_pos} unexpectedly placed in Chain A"
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
        raise CryoEMGatingKinematicsError(
            f"{pdb_id}: sequence mismatch against UniProt P13569: {mismatches[:5]}"
        )


def _parse_seq_dif_rows(lines: list[str]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    i = 0
    n = len(lines)
    while i < n and not lines[i].startswith("_struct_ref_seq_dif."):
        i += 1
    if i >= n:
        return rows
    # Check whether preceded by loop_
    prev = lines[i - 1].strip() if i > 0 else ""
    if prev == "loop_":
        cols: list[str] = []
        while i < n and lines[i].startswith("_struct_ref_seq_dif."):
            cols.append(lines[i].strip()[len("_struct_ref_seq_dif."):])
            i += 1
        while i < n and not lines[i].startswith(("#", "_", "loop_")):
            if lines[i].strip():
                toks = shlex.split(lines[i])
                if len(toks) == len(cols):
                    rows.append(dict(zip(cols, toks)))
            i += 1
    else:
        single: dict[str, str] = {}
        while i < n and lines[i].startswith("_struct_ref_seq_dif."):
            toks = shlex.split(lines[i])
            key = toks[0][len("_struct_ref_seq_dif."):]
            if len(toks) >= 2:
                single[key] = toks[1]
            elif i + 1 < n and not lines[i + 1].startswith(("#", "_", "loop_")):
                val_toks = shlex.split(lines[i + 1])
                if val_toks:
                    single[key] = val_toks[0]
                i += 1
            i += 1
        if single:
            rows.append(single)
    return rows


def parse_chain_a_ca(
    cif_path: Path,
    pdb_id: str,
    expected_sha256: str,
    uniprot_seq: str | None = None,
    *,
    verify_hash: bool = True,
    enforce_chain_a: bool = True,
    honor_seq_dif: bool = True,
    seq_offset: int = 0,
    rigid_transform: Callable[[tuple[float, float, float]], tuple[float, float, float]] | None = None,
) -> ChainAStructure:
    if verify_hash:
        actual = sha256_file(cif_path)
        if actual != expected_sha256:
            raise CryoEMGatingKinematicsError(
                f"{pdb_id} SHA-256 mismatch: expected {expected_sha256}, got {actual}"
            )
    lines = gzip.open(cif_path, "rt", encoding="utf-8").read().splitlines()

    declared_mutations: dict[int, str] = {}
    declared_deletions: set[int] = set()
    if honor_seq_dif:
        for row in _parse_seq_dif_rows(lines):
            if row.get("pdbx_pdb_strand_id") != "A" or row.get("pdbx_seq_db_accession_code") != "P13569":
                continue
            details = (row.get("details") or "").lower()
            if "expression tag" in details or "cloning" in details:
                continue
            pos_str = row.get("pdbx_seq_db_seq_num", "")
            if pos_str.isdigit():
                p_int = int(pos_str)
                if 1 <= p_int <= 1480:
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
    ca_by_seq: dict[int, tuple[str, tuple[float, float, float]]] = {}

    for r in atom_rows:
        if r[ix["group_PDB"]] != "ATOM":
            continue
        if r[ix["pdbx_PDB_model_num"]] != "1":
            continue
        if enforce_chain_a:
            if r[ix["auth_asym_id"]] != "A" or r[ix["label_asym_id"]] != "A":
                continue
        if r[ix["pdbx_PDB_ins_code"]] not in {"?", "."}:
            continue
        if r[ix["label_alt_id"]] not in {".", "A"}:
            continue
        if r[ix["label_atom_id"]] != "CA":
            continue

        raw_seq = r[ix["auth_seq_id"]]
        if not raw_seq.lstrip("-").isdigit():
            continue
        seq_id = int(raw_seq) + seq_offset
        if seq_id < 1 or seq_id > 1480:
            continue
        comp_id = r[ix["label_comp_id"]]
        if seq_id in ca_by_seq:
            raise CryoEMGatingKinematicsError(
                f"{pdb_id}: duplicate C-alpha record at codon {seq_id} ({comp_id} vs {ca_by_seq[seq_id][0]})"
            )

        xyz = (
            float(r[ix["Cartn_x"]]),
            float(r[ix["Cartn_y"]]),
            float(r[ix["Cartn_z"]]),
        )
        if rigid_transform is not None:
            xyz = rigid_transform(xyz)
        ca_by_seq[seq_id] = (comp_id, xyz)

    if uniprot_seq is not None:
        verify_ca_against_uniprot(
            ca_by_seq,
            uniprot_seq,
            declared_mutations,
            frozenset(declared_deletions),
            pdb_id=pdb_id,
        )

    return ChainAStructure(
        pdb_id=pdb_id,
        ca_by_seq=ca_by_seq,
        declared_mutations=declared_mutations,
        declared_deletions=frozenset(declared_deletions),
    )


def compute_symmetric_lddt_kinematics(
    struct_a: ChainAStructure,
    struct_b: ChainAStructure,
    *,
    symmetric_neighborhood: bool = True,
    empty_fallback: float | None = None,
) -> tuple[dict[int, float | None], dict[int, int]]:
    """Compute s_kin(i) = 1 - lDDT_sym(i; struct_a, struct_b) over S_common = C_a intersect C_b."""
    common_residues = sorted(set(struct_a.ca_by_seq.keys()) & set(struct_b.ca_by_seq.keys()))
    coords_a = {i: struct_a.ca_by_seq[i][1] for i in common_residues}
    coords_b = {i: struct_b.ca_by_seq[i][1] for i in common_residues}

    grid_a = SpatialGrid15A(coords_a, cell_size=LDDT_INCLUSION_RADIUS_A)
    grid_b = SpatialGrid15A(coords_b, cell_size=LDDT_INCLUSION_RADIUS_A)

    s_kin_by_seq: dict[int, float | None] = {}
    neighbor_counts: dict[int, int] = {}

    for i in common_residues:
        na = grid_a.neighbors_within(i, radius=LDDT_INCLUSION_RADIUS_A)
        if symmetric_neighborhood:
            nb = grid_b.neighbors_within(i, radius=LDDT_INCLUSION_RADIUS_A)
            n_sym = sorted(na | nb)
        else:
            n_sym = sorted(na)
        neighbor_counts[i] = len(n_sym)
        if not n_sym:
            if empty_fallback is not None:
                raise CryoEMGatingKinematicsError(
                    f"Numeric fallback {empty_fallback} forbidden on empty N_sym({i}); must return None"
                )
            s_kin_by_seq[i] = None
            continue

        preserved = 0
        pa = coords_a[i]
        pb = coords_b[i]
        for j in n_sym:
            da = math.dist(pa, coords_a[j])
            db = math.dist(pb, coords_b[j])
            diff = abs(da - db)
            for t in LDDT_THRESHOLDS_A:
                if abs(diff - t) <= BOUNDARY_EPSILON_A:
                    raise CryoEMGatingKinematicsError(
                        f"Threshold knife-edge pair ({i}, {j}): |da-db|={diff:.12f} vs t={t}"
                    )
                if diff < t:
                    preserved += 1
        denom = 4 * len(n_sym)
        s_kin_by_seq[i] = (denom - preserved) / denom

    return s_kin_by_seq, neighbor_counts


def score_variant_kinematics(
    variant: dict[str, Any],
    s_kin_by_seq: dict[int, float | None],
    *,
    empty_fallback: float | None = None,
) -> float | None:
    fam = variant["variant_family"]
    if fam not in ALLOWED_FAMILIES:
        raise CryoEMGatingKinematicsError(f"Invalid variant_family {fam!r}")
    override = variant.get("_override_s_kin")
    if fam != "missense":
        if override is not None:
            raise CryoEMGatingKinematicsError(
                f"Non-missense variant {variant['variant_id']} ({fam}) cannot carry non-None s_kin={override!r}"
            )
        return None
    codon = int(variant["codon_position"])
    if codon not in s_kin_by_seq:
        if empty_fallback is not None:
            raise CryoEMGatingKinematicsError(
                f"Numeric fallback {empty_fallback} forbidden for unplaced codon {codon}; must return None"
            )
        return None
    val = s_kin_by_seq[codon]
    if val is None and empty_fallback is not None:
        raise CryoEMGatingKinematicsError(
            f"Numeric fallback {empty_fallback} forbidden on empty N_sym({codon}); must return None"
        )
    return val


def sequence_consequence_baseline_score(variant: dict[str, Any]) -> float:
    fam = variant["variant_family"]
    if fam not in ALLOWED_FAMILIES:
        raise CryoEMGatingKinematicsError(f"Invalid variant_family {fam!r}")
    hgvs = str(variant["hgvs_p_or_c"])
    if fam == "class_i_ptc_frameshift":
        return 0.0
    if fam == "splice_regulatory":
        return 0.0 if CANONICAL_SPLICE_RE.match(hgvs) else 0.50
    if fam == "inframe_indel":
        return 0.50
    return 0.75


def class_gated_score(variant: dict[str, Any], missense_score: float | None) -> float | None:
    fam = variant["variant_family"]
    if fam not in ALLOWED_FAMILIES:
        raise CryoEMGatingKinematicsError(f"Invalid variant_family {fam!r}")
    if fam != "missense":
        if missense_score is not None:
            raise CryoEMGatingKinematicsError(
                f"Non-missense variant {variant['variant_id']} ({fam}) received non-None missense_score"
            )
        return sequence_consequence_baseline_score(variant)
    if missense_score is None:
        return None
    return 0.50 + 0.50 * float(missense_score)


def evaluate_mann_whitney_auroc(
    labels: list[int],
    scores: list[float | None],
) -> dict[str, float | int | None]:
    if len(labels) != len(scores) or not labels:
        raise CryoEMGatingKinematicsError("Labels and scores must be non-empty and equal length")
    for y in labels:
        if isinstance(y, bool) or y not in (0, 1):
            raise CryoEMGatingKinematicsError(f"Labels must be binary integers in {{0, 1}}, got {y!r}")
    for s in scores:
        if s is not None and (isinstance(s, bool) or not math.isfinite(float(s))):
            raise CryoEMGatingKinematicsError(f"Score must be finite float or None, got {s!r}")

    pos_idx = [idx for idx, y in enumerate(labels) if y == 1]
    neg_idx = [idx for idx, y in enumerate(labels) if y == 0]
    if not pos_idx or not neg_idx:
        raise CryoEMGatingKinematicsError("Both positive and negative labels are required for ROC-AUC")

    total_pairs = len(pos_idx) * len(neg_idx)
    abst_u = 0.0
    cond_u = 0.0
    cond_pairs = 0
    scored_count = sum(1 for s in scores if s is not None)

    for pi in pos_idx:
        sp = scores[pi]
        for ni in neg_idx:
            sn = scores[ni]
            if sp is None or sn is None:
                abst_u += 0.5
            else:
                cond_pairs += 1
                if sp > sn:
                    abst_u += 1.0
                    cond_u += 1.0
                elif sp == sn:
                    abst_u += 0.5
                    cond_u += 0.5

    cond_auc = (cond_u / cond_pairs) if cond_pairs > 0 else None
    abst_auc = abst_u / total_pairs
    return {
        "n_total": len(labels),
        "n_scored": scored_count,
        "coverage": scored_count / len(labels),
        "conditional_pairs": cond_pairs,
        "conditional_u": cond_u,
        "conditional_auroc": cond_auc,
        "total_pairs": total_pairs,
        "abstention_u": abst_u,
        "abstention_aware_auroc": abst_auc,
    }


def validate_split_integrity(variants: list[dict[str, Any]]) -> int:
    train_codons: list[int] = []
    test_codons: list[int] = []
    seen_codons: set[int] = set()

    for v in variants:
        fam = v.get("variant_family")
        if fam not in ALLOWED_FAMILIES:
            raise CryoEMGatingKinematicsError(f"Disallowed variant_family {fam!r}")
        if fam != "missense" and (
            v.get("alphamissense_score") is not None or v.get("revel_score") is not None
        ):
            raise CryoEMGatingKinematicsError(
                f"Non-missense {v['variant_id']} carries non-None missense pathogenicity score"
            )
        codon = int(v["codon_position"])
        if codon in seen_codons:
            raise CryoEMGatingKinematicsError(f"Duplicate codon {codon} across benchmark")
        seen_codons.add(codon)
        dom = map_p13569_domain(codon)
        split = v.get("split")
        if split == "train":
            if dom not in TRAIN_DOMAINS:
                raise CryoEMGatingKinematicsError(f"Train variant {v['variant_id']} in non-train domain {dom}")
            train_codons.append(codon)
        elif split == "test":
            if dom not in TEST_DOMAINS:
                raise CryoEMGatingKinematicsError(f"Test variant {v['variant_id']} in non-test domain {dom}")
            test_codons.append(codon)
        else:
            raise CryoEMGatingKinematicsError(f"Unknown split {split!r}")

    min_gap = min(abs(tc - hc) for tc in train_codons for hc in test_codons)
    if min_gap <= 5:
        raise CryoEMGatingKinematicsError(f"Cross-split codon gap {min_gap} <= 5 codons")
    return min_gap


def _make_proper_rigid_motion(
    theta_xy: float,
    phi_yz: float,
    shift: tuple[float, float, float],
) -> Callable[[tuple[float, float, float]], tuple[float, float, float]]:
    c1, s1 = math.cos(theta_xy), math.sin(theta_xy)
    c2, s2 = math.cos(phi_yz), math.sin(phi_yz)
    tx, ty, tz = shift

    def apply(pt: tuple[float, float, float]) -> tuple[float, float, float]:
        x, y, z = pt
        x1 = c1 * x - s1 * y
        y1 = s1 * x + c1 * y
        z1 = z
        x2 = x1
        y2 = c2 * y1 - s2 * z1
        z2 = s2 * y1 + c2 * z1
        return (x2 + tx, y2 + ty, z2 + tz)

    return apply


def verify_positive_controls(
    struct_5uak: ChainAStructure,
    struct_6msm: ChainAStructure,
    s_kin_by_seq: dict[int, float | None],
    benchmark_codons: frozenset[int],
    uniprot_seq: str,
) -> dict[str, Any]:
    if len(struct_5uak.ca_by_seq) != 1139:
        raise CryoEMGatingKinematicsError(f"Expected 1,139 C-alpha in 5UAK, got {len(struct_5uak.ca_by_seq)}")
    if len(struct_6msm.ca_by_seq) != 1181:
        raise CryoEMGatingKinematicsError(f"Expected 1,181 C-alpha in 6MSM, got {len(struct_6msm.ca_by_seq)}")
    if struct_6msm.declared_mutations.get(1371) != "GLN":
        raise CryoEMGatingKinematicsError(
            f"Expected E1371Q (GLN at 1371) in 6MSM _struct_ref_seq_dif, got {struct_6msm.declared_mutations}"
        )
    if len(s_kin_by_seq) != 1130:
        raise CryoEMGatingKinematicsError(f"Expected 1,130 shared C-alpha in S_common, got {len(s_kin_by_seq)}")

    control_residues = (
        {r for pair in DIMER_CLOSURE_PAIRS for r in pair}
        | set(GATING_HINGE_CONTROLS)
        | set(RIGID_BETA_CONTROLS)
    )
    min_control_gap = min(abs(c - b) for c in control_residues for b in benchmark_codons)
    if min_control_gap == 0:
        raise CryoEMGatingKinematicsError("Non-benchmark positive control overlaps a benchmark codon")

    closure_deltas: dict[str, dict[str, float]] = {}
    for r1, r2 in DIMER_CLOSURE_PAIRS:
        d_apo = math.dist(struct_5uak.ca_by_seq[r1][1], struct_5uak.ca_by_seq[r2][1])
        d_atp = math.dist(struct_6msm.ca_by_seq[r1][1], struct_6msm.ca_by_seq[r2][1])
        contraction = d_apo - d_atp
        if contraction <= 15.0 or d_atp >= 12.0:
            raise CryoEMGatingKinematicsError(
                f"NBD dimer closure control failed on ({r1}, {r2}): d_5UAK={d_apo:.4f}, d_6MSM={d_atp:.4f}"
            )
        closure_deltas[f"{r1}_{r2}"] = {
            "d_5uak": round(d_apo, 4),
            "d_6msm": round(d_atp, 4),
            "contraction_a": round(contraction, 4),
        }

    for g_res in GATING_HINGE_CONTROLS:
        sg = s_kin_by_seq.get(g_res)
        if sg is None or sg <= 0.45:
            raise CryoEMGatingKinematicsError(f"Gating control {g_res} failed s_kin > 0.45: {sg}")
        for b_res in RIGID_BETA_CONTROLS:
            sb = s_kin_by_seq.get(b_res)
            if sb is None or sb >= 0.25 or not (sg > sb):
                raise CryoEMGatingKinematicsError(
                    f"Gating vs scaffolding ordering failed: s_kin({g_res})={sg} vs s_kin({b_res})={sb}"
                )

    # Independent SE(3) x SE(3) rigid-motion invariance across all 1,130 shared residues
    rot_5uak = parse_chain_a_ca(
        DEFAULT_CIF_5UAK,
        "5UAK",
        EXPECTED_5UAK_SHA256,
        uniprot_seq,
        verify_hash=False,
        rigid_transform=_make_proper_rigid_motion(0.73, -0.41, (19.5, -42.0, 8.25)),
    )
    rot_6msm = parse_chain_a_ca(
        DEFAULT_CIF_6MSM,
        "6MSM",
        EXPECTED_6MSM_SHA256,
        uniprot_seq,
        verify_hash=False,
        rigid_transform=_make_proper_rigid_motion(-1.12, 0.89, (-33.0, 14.75, -61.5)),
    )
    s_rot, _ = compute_symmetric_lddt_kinematics(rot_5uak, rot_6msm)
    max_se3_drift = max(
        abs(float(s_kin_by_seq[k]) - float(s_rot[k]))
        for k in s_kin_by_seq
        if s_kin_by_seq[k] is not None and s_rot[k] is not None
    )
    if max_se3_drift >= 1e-12:
        raise CryoEMGatingKinematicsError(f"SE(3)xSE(3) drift {max_se3_drift:.3e} >= 1e-12")

    # Exact swap symmetry s_kin(5UAK, 6MSM) == s_kin(6MSM, 5UAK)
    s_swapped, _ = compute_symmetric_lddt_kinematics(struct_6msm, struct_5uak)
    max_swap_diff = max(
        abs(float(s_kin_by_seq[k]) - float(s_swapped[k]))
        for k in s_kin_by_seq
        if s_kin_by_seq[k] is not None and s_swapped[k] is not None
    )
    if max_swap_diff != 0.0:
        raise CryoEMGatingKinematicsError(f"State swap asymmetry detected: {max_swap_diff}")

    # Train-split mechanistic control: ICL1 G178 (KALYDECO=1) > TM1 G85 (KALYDECO=0)
    s178 = s_kin_by_seq[178]
    s85 = s_kin_by_seq[85]
    if s178 is None or s85 is None or not (s178 > s85):
        raise CryoEMGatingKinematicsError(
            f"Train-split mechanistic control failed: s_kin(178)={s178} <= s_kin(85)={s85}"
        )

    # ROC-AUC unit controls
    auc_none = evaluate_mann_whitney_auroc([1, 1, 0, 0], [None, None, None, None])["abstention_aware_auroc"]
    auc_pos = evaluate_mann_whitney_auroc([1, 1, 0, 0], [0.9, 0.8, 0.2, 0.1])["abstention_aware_auroc"]
    auc_inv = evaluate_mann_whitney_auroc([1, 1, 0, 0], [0.1, 0.2, 0.8, 0.9])["abstention_aware_auroc"]
    if auc_none != 0.5 or auc_pos != 1.0 or auc_inv != 0.0:
        raise CryoEMGatingKinematicsError(
            f"ROC-AUC unit control failed: none={auc_none}, pos={auc_pos}, inv={auc_inv}"
        )

    return {
        "n_5uak_ca": len(struct_5uak.ca_by_seq),
        "n_6msm_ca": len(struct_6msm.ca_by_seq),
        "n_common_ca": len(s_kin_by_seq),
        "min_control_codon_gap": min_control_gap,
        "dimer_closure_pairs": closure_deltas,
        "gating_hinge_controls": {str(r): round(float(s_kin_by_seq[r]), 6) for r in GATING_HINGE_CONTROLS},
        "rigid_beta_controls": {str(r): round(float(s_kin_by_seq[r]), 6) for r in RIGID_BETA_CONTROLS},
        "train_g178_vs_g85": {"G178": round(float(s178), 6), "G85": round(float(s85), 6)},
        "max_se3_x_se3_drift": max_se3_drift,
        "max_swap_symmetry_diff": max_swap_diff,
    }


def run_falsification_suite(
    struct_5uak: ChainAStructure,
    struct_6msm: ChainAStructure,
    s_kin_by_seq: dict[int, float | None],
    variants: list[dict[str, Any]],
    uniprot_seq: str,
    admitted_ids: frozenset[str],
) -> dict[str, bool]:
    results: dict[str, bool] = {}

    # 1. Stripping _struct_ref_seq_dif (E1371Q) in 6MSM fails P13569 verification
    try:
        parse_chain_a_ca(
            DEFAULT_CIF_6MSM,
            "6MSM",
            EXPECTED_6MSM_SHA256,
            uniprot_seq,
            verify_hash=False,
            honor_seq_dif=False,
        )
        raise AssertionError("Expected CryoEMGatingKinematicsError when E1371Q _struct_ref_seq_dif is ignored")
    except CryoEMGatingKinematicsError:
        results["missing_e1371q_seq_dif_rejected"] = True

    # 2. Omitting Chain A filter admits unassigned R-domain UNK chains (1..19) and fails duplicate/sequence check
    try:
        parse_chain_a_ca(
            DEFAULT_CIF_5UAK,
            "5UAK",
            EXPECTED_5UAK_SHA256,
            uniprot_seq,
            verify_hash=False,
            enforce_chain_a=False,
        )
        raise AssertionError("Expected CryoEMGatingKinematicsError when unassigned UNK chains are admitted")
    except CryoEMGatingKinematicsError:
        results["unk_chain_collision_rejected"] = True

    # 3. Numeric fallback (0.0 or 1.0) on unplaced residue or empty N_sym is rejected
    isolated_a = ChainAStructure("A", {100: ("ALA", (0.0, 0.0, 0.0)), 200: ("GLY", (90.0, 0.0, 0.0))}, {}, frozenset())
    isolated_b = ChainAStructure("B", {100: ("ALA", (0.0, 0.0, 0.0)), 200: ("GLY", (90.0, 0.0, 0.0))}, {}, frozenset())
    s_iso, _ = compute_symmetric_lddt_kinematics(isolated_a, isolated_b)
    if s_iso[100] is not None:
        raise AssertionError("Isolated residue with empty N_sym must return None")
    try:
        compute_symmetric_lddt_kinematics(isolated_a, isolated_b, empty_fallback=0.0)
        raise AssertionError("Expected CryoEMGatingKinematicsError on numeric empty_fallback=0.0")
    except CryoEMGatingKinematicsError:
        results["empty_neighborhood_numeric_fallback_rejected"] = True

    # 4. Asymmetric one-sided neighborhood N_5UAK(i) violates state-swap symmetry
    asym_forward, _ = compute_symmetric_lddt_kinematics(struct_5uak, struct_6msm, symmetric_neighborhood=False)
    asym_reverse, _ = compute_symmetric_lddt_kinematics(struct_6msm, struct_5uak, symmetric_neighborhood=False)
    max_asym = max(
        abs(float(asym_forward[k]) - float(asym_reverse[k]))
        for k in asym_forward
        if asym_forward[k] is not None and asym_reverse[k] is not None
    )
    if max_asym <= 1e-6:
        raise AssertionError("Expected one-sided neighborhood to violate swap symmetry")
    results["asymmetric_state_swap_rejected"] = True

    # 5. Non-None s_kin on non-missense variant (F508del) is rejected
    f508_row = dict(next(v for v in variants if v["legacy_spl_name"] == "F508del"))
    f508_row["_override_s_kin"] = 0.5
    try:
        score_variant_kinematics(f508_row, s_kin_by_seq)
        raise AssertionError("Expected CryoEMGatingKinematicsError on non-missense _override_s_kin")
    except CryoEMGatingKinematicsError:
        results["non_missense_kinematic_score_rejected"] = True

    # 6. Post-hoc sign flip on strict separator yields ROC-AUC 0.0 (never max(AUC, 1-AUC))
    inv_auc = evaluate_mann_whitney_auroc([1, 1, 1, 0, 0, 0], [0.05, 0.10, 0.15, 0.85, 0.90, 0.95])[
        "abstention_aware_auroc"
    ]
    if inv_auc != 0.0:
        raise AssertionError(f"Expected 0.0 on inverted separator, got {inv_auc}")
    results["post_hoc_sign_flip_rejected"] = True

    # 7. Codon gap <= 5 leakage AND manual source auto-ingest rejected
    poisoned_split = [dict(v) for v in variants]
    poisoned_split.append(
        {
            "variant_id": "CFTR-LEAK-383",
            "legacy_spl_name": "L383P",
            "hgvs_p_or_c": "p.Leu383Pro",
            "codon_position": 383,
            "variant_family": "missense",
            "split": "test",
            "alphamissense_score": 0.5,
            "revel_score": 0.5,
            "kalydeco_v36_responsive": 0,
        }
    )
    poisoned_split.append(
        {
            "variant_id": "CFTR-LEAK-380",
            "legacy_spl_name": "L380P",
            "hgvs_p_or_c": "p.Leu380Pro",
            "codon_position": 380,
            "variant_family": "missense",
            "split": "train",
            "alphamissense_score": 0.5,
            "revel_score": 0.5,
            "kalydeco_v36_responsive": 0,
        }
    )
    codon_caught = False
    try:
        validate_split_integrity(poisoned_split)
    except CryoEMGatingKinematicsError:
        codon_caught = True
    manual_caught = False
    try:
        assert_catalog_sources_admitted(["mavedb-dump-zenodo"], admitted_ids)
    except PermissionError:
        manual_caught = True
    if not (codon_caught and manual_caught):
        raise AssertionError("Expected both codon leakage and manual source auto-ingest to be rejected")
    results["codon_leakage_or_manual_source_rejected"] = True

    return results


def run_benchmark(run_falsify: bool = False) -> dict[str, Any]:
    admitted_ids = load_admitted_automated_catalog_ids()
    assert_catalog_sources_admitted(
        ["rcsb-pdb-5uak", "rcsb-pdb-6msm", "uniprot-p13569", "openfda-label", "myvariant"],
        admitted_ids,
    )
    ann_payload = json.loads(DEFAULT_ANNOTATIONS_JSON.read_text(encoding="utf-8"))
    assert_catalog_sources_admitted(ann_payload["provenance"]["catalog_sources"], admitted_ids)

    variants: list[dict[str, Any]] = list(ann_payload["variants"])
    min_cross_split_gap = validate_split_integrity(variants)
    benchmark_codons = frozenset(int(v["codon_position"]) for v in variants)

    uniprot_seq = load_uniprot_sequence(DEFAULT_UNIPROT_JSON, verify_hash=True)
    struct_5uak = parse_chain_a_ca(DEFAULT_CIF_5UAK, "5UAK", EXPECTED_5UAK_SHA256, uniprot_seq, verify_hash=True)
    struct_6msm = parse_chain_a_ca(DEFAULT_CIF_6MSM, "6MSM", EXPECTED_6MSM_SHA256, uniprot_seq, verify_hash=True)
    s_kin_by_seq, neighbor_counts = compute_symmetric_lddt_kinematics(struct_5uak, struct_6msm)

    # Positive controls MUST succeed before any held-out test evaluation
    pos_controls = verify_positive_controls(
        struct_5uak,
        struct_6msm,
        s_kin_by_seq,
        benchmark_codons,
        uniprot_seq,
    )

    def evaluate_split(split_name: str) -> dict[str, Any]:
        rows = [v for v in variants if v["split"] == split_name]
        labels = [int(v["kalydeco_v36_responsive"]) for v in rows]
        s_kin_scores = [score_variant_kinematics(v, s_kin_by_seq) for v in rows]
        s_am_scores = [
            (1.0 - float(v["alphamissense_score"])) if v["alphamissense_score"] is not None else None
            for v in rows
        ]
        g_kin_scores = [
            class_gated_score(v, sk if v["variant_family"] == "missense" else None)
            for v, sk in zip(rows, s_kin_scores)
        ]
        g_am_scores = [
            class_gated_score(v, sa if v["variant_family"] == "missense" else None)
            for v, sa in zip(rows, s_am_scores)
        ]
        base_scores = [sequence_consequence_baseline_score(v) for v in rows]

        missense_idx = [idx for idx, v in enumerate(rows) if v["variant_family"] == "missense"]
        m_labels = [labels[idx] for idx in missense_idx]
        m_skin = [s_kin_scores[idx] for idx in missense_idx]
        m_sam = [s_am_scores[idx] for idx in missense_idx]

        per_variant = []
        for v, y, sk, sa, gk, ga, bs in zip(
            rows, labels, s_kin_scores, s_am_scores, g_kin_scores, g_am_scores, base_scores
        ):
            c_pos = int(v["codon_position"])
            per_variant.append(
                {
                    "variant_id": v["variant_id"],
                    "legacy_spl_name": v["legacy_spl_name"],
                    "codon_position": c_pos,
                    "domain": map_p13569_domain(c_pos),
                    "variant_family": v["variant_family"],
                    "kalydeco_v36_responsive": y,
                    "n_sym_neighbors": neighbor_counts.get(c_pos) if v["variant_family"] == "missense" else None,
                    "s_kin": round(sk, 6) if sk is not None else None,
                    "s_am": round(sa, 6) if sa is not None else None,
                    "g_kin": round(gk, 6) if gk is not None else None,
                    "g_am": round(ga, 6) if ga is not None else None,
                    "class_baseline": bs,
                }
            )

        return {
            "n_variants": len(rows),
            "n_responsive": sum(labels),
            "n_non_responsive": len(labels) - sum(labels),
            "missense_s_kin": evaluate_mann_whitney_auroc(m_labels, m_skin),
            "missense_s_am": evaluate_mann_whitney_auroc(m_labels, m_sam),
            "all_standalone_s_kin": evaluate_mann_whitney_auroc(labels, s_kin_scores),
            "all_gated_g_kin": evaluate_mann_whitney_auroc(labels, g_kin_scores),
            "all_gated_g_am": evaluate_mann_whitney_auroc(labels, g_am_scores),
            "all_sequence_consequence_baseline": evaluate_mann_whitney_auroc(labels, base_scores),
            "variants": per_variant,
        }

    train_eval = evaluate_split("train")
    test_eval = evaluate_split("test")

    test_cond_skin = test_eval["missense_s_kin"]["conditional_auroc"]
    test_cond_sam = test_eval["missense_s_am"]["conditional_auroc"]
    test_abst_gkin = test_eval["all_gated_g_kin"]["abstention_aware_auroc"]
    test_abst_gam = test_eval["all_gated_g_am"]["abstention_aware_auroc"]
    test_abst_base = test_eval["all_sequence_consequence_baseline"]["abstention_aware_auroc"]

    pred_supported = bool(
        test_cond_skin is not None
        and test_cond_sam is not None
        and test_cond_skin >= 0.75
        and test_cond_skin > test_cond_sam
        and test_abst_gkin >= 0.75
        and test_abst_gkin > test_abst_base
        and test_abst_gkin > test_abst_gam
    )
    outcome = "supported" if pred_supported else "refuted"

    falsify_results: dict[str, bool] = {}
    if run_falsify:
        falsify_results = run_falsification_suite(
            struct_5uak, struct_6msm, s_kin_by_seq, variants, uniprot_seq, admitted_ids
        )

    return {
        "prediction_id": "PRED-0012",
        "split_version": ann_payload["schema_version"],
        "min_cross_split_codon_gap": min_cross_split_gap,
        "positive_controls": pos_controls,
        "train": train_eval,
        "test": test_eval,
        "outcome": outcome,
        "held_out_error_1_minus_g_kin_auroc": round(1.0 - float(test_abst_gkin), 6),
        "falsification_controls": falsify_results,
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Modulator-free 5UAK <-> 6MSM gating kinematics benchmark (PRED-0012).")
    ap.add_argument("--check", action="store_true", help="Run positive controls and pre-registered PRED-0012 test.")
    ap.add_argument("--falsify", action="store_true", help="Run all 7 planted negative controls.")
    ap.add_argument("--json", action="store_true", help="Print full JSON report.")
    args = ap.parse_args(argv)

    rep = run_benchmark(run_falsify=args.falsify)
    if args.json:
        print(json.dumps(rep, indent=2))
        return 0

    t_skin = rep["test"]["missense_s_kin"]
    t_sam = rep["test"]["missense_s_am"]
    t_gkin = rep["test"]["all_gated_g_kin"]
    t_gam = rep["test"]["all_gated_g_am"]
    t_base = rep["test"]["all_sequence_consequence_baseline"]
    print(
        f"PASS cryoem_gating_kinematics ({rep['prediction_id']} -> outcome={rep['outcome']}): "
        f"S_common={rep['positive_controls']['n_common_ca']} (5UAK={rep['positive_controls']['n_5uak_ca']}, "
        f"6MSM={rep['positive_controls']['n_6msm_ca']}); "
        f"test missense s_kin AUC_cond={t_skin['conditional_auroc']:.6f} ({t_skin['conditional_u']}/{t_skin['conditional_pairs']}) "
        f"vs s_AM={t_sam['conditional_auroc']:.6f}; "
        f"test all-20 g_kin AUC_abst={t_gkin['abstention_aware_auroc']:.6f} ({t_gkin['abstention_u']}/{t_gkin['total_pairs']}) "
        f"vs baseline={t_base['abstention_aware_auroc']:.6f} and g_AM={t_gam['abstention_aware_auroc']:.6f}"
        + (f"; falsify={len(rep['falsification_controls'])}/7 PASS" if args.falsify else "")
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
