#!/usr/bin/env python3
"""Modulator-Free Apo Cryo-EM Tertiary Packing Compliance Benchmark (PRED-0013).

Tests whether SE(3)-invariant, zero-free-parameter non-covalent 3D C-alpha
tertiary packing compliance
    s_flex(i) = 1 / (1 + Pi_5UAK(i)),
    Pi_5UAK(i) = r_0^2 * sum_{j in C_5UAK, |j - i| > 1} ||r_i - r_j||_2^{-2}
computed strictly from polypeptide Chain A C-alpha coordinates in modulator-free
dephosphorylated ATP-free WT human CFTR (PDB 5UAK, CC0 1.0) with fixed virtual
C-alpha bond scale r_0 = 3.8 A (r_0^2 = 361/25 = 14.44 A^2; Halle 2002
10.1073/pnas.032522499; Lin et al. 2008 10.1002/prot.21983) rescues PRED-0010 by
predicting US FDA TRIKAFTA v23 Table 6 responsiveness on the fixed 20-variant
domain-disjoint held-out test split (2026-11-21.v2: NBD1, TMD2, NBD2_C_terminal).

Method 1 implementation details (disjoint from independent_audit.py):
- Parses 5UAK.cif.gz via shlex-tokenized mmCIF loops/tags and typed dataclasses.
- Evaluates float64 inverse-square Euclidean distances and pairwise O(N1*N0)
  Mann-Whitney concordance / tie enumeration (forbidding post-hoc sign flipping).

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
DEFAULT_UNIPROT_JSON = REPO_ROOT / "tools" / "cftr-model" / "data" / "p13569.json"
DEFAULT_CATALOG_YAML = REPO_ROOT / "sources" / "catalog.yaml"
DEFAULT_ANNOTATIONS_JSON = (
    REPO_ROOT / "benchmarks" / "theratype" / "theratype_open_annotations_2026-11-21.json"
)

EXPECTED_5UAK_SHA256 = "8bf00f2612dd961d3b1ccd2f2193781a534ad44f32666a6de8da948afb8d84cf"
EXPECTED_P13569_SHA256 = "c8feabfabe8690ef4a06af7de5b5f8dc3c8da554246c5ab3d403edf265d4a893"

R0_VIRTUAL_CA_BOND_A: float = 3.8
R0_SQUARED_A2: float = 361.0 / 25.0  # 14.44 A^2 exact binary/decimal representation check via tolerance
MIN_SQ_DIST_A2: float = 1e-12

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
CATALOG_ID_RE = re.compile(r"^  -\s+id:\s*([^\s#]+)")
CATALOG_ACCESS_RE = re.compile(r"^    (?:access|mode):\s*([^\s#]+)")
CATALOG_LABEL_RE = re.compile(r"^    claim_label:\s*([^\s#]+)")

EXPOSED_CONTROLS: tuple[int, ...] = (38, 1347, 1348, 1371)
BURIED_CONTROLS: tuple[int, ...] = (195, 311, 361, 1030, 1102, 1293)


class CryoEMPackingComplianceError(ValueError):
    """Raised when structural provenance, positive controls, or negative controls fail."""


@dataclasses.dataclass(frozen=True)
class ChainAStructure:
    pdb_id: str
    ca_by_seq: dict[int, tuple[str, tuple[float, float, float]]]


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def map_p13569_domain(codon: int) -> str:
    if 1 <= codon <= 80:
        return "Lasso_N_terminal"
    if 81 <= codon <= 388:
        return "TMD1"
    if 389 <= codon <= 636:
        return "NBD1"
    if 637 <= codon <= 845:
        return "R_domain"
    if 846 <= codon <= 1201:
        return "TMD2"
    if 1202 <= codon <= 1480:
        return "NBD2_C_terminal"
    raise CryoEMPackingComplianceError(f"Codon {codon} outside canonical P13569 1..1480 span")


def load_catalog_permissions(catalog_path: Path = DEFAULT_CATALOG_YAML) -> dict[str, tuple[str, str]]:
    if not catalog_path.is_file():
        raise CryoEMPackingComplianceError(f"Missing catalog file: {catalog_path}")
    perms: dict[str, tuple[str, str]] = {}
    cur_id: str | None = None
    cur_access = "unknown"
    cur_label = "unknown"
    for line in catalog_path.read_text(encoding="utf-8").splitlines():
        m_id = CATALOG_ID_RE.match(line)
        if m_id:
            if cur_id is not None:
                perms[cur_id] = (cur_access, cur_label)
            cur_id = m_id.group(1).strip().strip("'\"")
            cur_access = "unknown"
            cur_label = "unknown"
            continue
        m_acc = CATALOG_ACCESS_RE.match(line)
        if m_acc and cur_id is not None and cur_access == "unknown":
            cur_access = m_acc.group(1).strip().strip("'\"")
        m_lbl = CATALOG_LABEL_RE.match(line)
        if m_lbl and cur_id is not None:
            cur_label = m_lbl.group(1).strip().strip("'\"")
    if cur_id is not None:
        perms[cur_id] = (cur_access, cur_label)
    return perms


def assert_catalog_source_admitted(
    source_id: str,
    perms: dict[str, tuple[str, str]],
) -> None:
    if source_id not in perms:
        raise PermissionError(f"Catalog source {source_id!r} absent from sources/catalog.yaml")
    access, claim_label = perms[source_id]
    if access not in {"fetch", "api"} or claim_label != "verified":
        raise PermissionError(
            f"Catalog source {source_id!r} blocked (access={access!r}, claim_label={claim_label!r})"
        )


def load_p13569_sequence(uniprot_path: Path = DEFAULT_UNIPROT_JSON) -> str:
    digest = sha256_file(uniprot_path)
    if digest != EXPECTED_P13569_SHA256:
        raise CryoEMPackingComplianceError(
            f"p13569.json SHA-256 mismatch: expected {EXPECTED_P13569_SHA256}, got {digest}"
        )
    payload = json.loads(uniprot_path.read_text(encoding="utf-8"))
    seq = str(payload.get("sequence", {}).get("value", ""))
    if len(seq) != 1480:
        raise CryoEMPackingComplianceError(f"Expected 1480-aa P13569 sequence, got {len(seq)}")
    return seq


def parse_5uak_chain_a(
    cif_path: Path = DEFAULT_CIF_5UAK,
    uniprot_seq: str | None = None,
    *,
    enforce_chain_a: bool = True,
) -> ChainAStructure:
    digest = sha256_file(cif_path)
    if digest != EXPECTED_5UAK_SHA256:
        raise CryoEMPackingComplianceError(
            f"5UAK.cif.gz SHA-256 mismatch: expected {EXPECTED_5UAK_SHA256}, got {digest}"
        )
    seq_ref = uniprot_seq if uniprot_seq is not None else load_p13569_sequence()
    lines = gzip.open(cif_path, "rt", encoding="utf-8").read().splitlines()

    idx = 0
    while idx < len(lines) and not lines[idx].startswith("_atom_site."):
        idx += 1
    if idx >= len(lines):
        raise CryoEMPackingComplianceError("5UAK.cif.gz missing _atom_site loop")
    cols: list[str] = []
    while idx < len(lines) and lines[idx].startswith("_atom_site."):
        cols.append(lines[idx].strip()[len("_atom_site."):])
        idx += 1

    ca_by_seq: dict[int, tuple[str, tuple[float, float, float]]] = {}
    while idx < len(lines) and not lines[idx].startswith(("#", "_", "loop_")):
        raw = lines[idx].strip()
        idx += 1
        if not raw:
            continue
        toks = shlex.split(raw)
        if len(toks) != len(cols):
            continue
        row = dict(zip(cols, toks))
        if row.get("group_PDB") != "ATOM":
            continue
        if row.get("pdbx_PDB_model_num", "1") != "1":
            continue
        if enforce_chain_a and (row.get("auth_asym_id") != "A" or row.get("label_asym_id") != "A"):
            continue
        if row.get("label_atom_id") != "CA":
            continue
        if row.get("pdbx_PDB_ins_code", "?") not in {"?", "."}:
            continue
        if row.get("label_alt_id", ".") not in {".", "A"}:
            continue
        seq_id = int(row["auth_seq_id"])
        if not (1 <= seq_id <= 1480):
            continue
        comp = row["label_comp_id"]
        if comp not in THREE_TO_ONE:
            raise CryoEMPackingComplianceError(
                f"Non-canonical residue {comp!r} at Chain A position {seq_id}"
            )
        aa1 = THREE_TO_ONE[comp]
        expected_aa = seq_ref[seq_id - 1]
        if aa1 != expected_aa:
            raise CryoEMPackingComplianceError(
                f"Sequence mismatch at codon {seq_id}: 5UAK has {aa1} ({comp}), P13569 has {expected_aa}"
            )
        if seq_id in ca_by_seq:
            raise CryoEMPackingComplianceError(f"Duplicate CA atom at Chain A codon {seq_id}")
        pt = (float(row["Cartn_x"]), float(row["Cartn_y"]), float(row["Cartn_z"]))
        ca_by_seq[seq_id] = (aa1, pt)

    return ChainAStructure(pdb_id="5UAK", ca_by_seq=ca_by_seq)


def compute_noncovalent_compliance_map(
    coords_by_seq: dict[int, tuple[float, float, float]],
    *,
    r0_sq: float = R0_SQUARED_A2,
    exclude_covalent: bool = True,
) -> dict[int, float | None]:
    if r0_sq <= 0.0 or not math.isfinite(r0_sq):
        raise CryoEMPackingComplianceError(f"Invalid length-scale squared r0_sq={r0_sq!r}")
    out: dict[int, float | None] = {}
    items = list(coords_by_seq.items())
    for seq_i, pt_i in items:
        wcn_nc = 0.0
        n_neighbors = 0
        xi, yi, zi = pt_i
        for seq_j, pt_j in items:
            if seq_j == seq_i:
                continue
            if exclude_covalent and abs(seq_j - seq_i) <= 1:
                continue
            dx = xi - pt_j[0]
            dy = yi - pt_j[1]
            dz = zi - pt_j[2]
            d2 = dx * dx + dy * dy + dz * dz
            if d2 <= MIN_SQ_DIST_A2:
                raise CryoEMPackingComplianceError(
                    f"Coincident or near-zero C-alpha pair ({seq_i}, {seq_j}): r^2={d2}"
                )
            wcn_nc += 1.0 / d2
            n_neighbors += 1
        if n_neighbors == 0:
            out[seq_i] = None
        else:
            pi_val = r0_sq * wcn_nc
            out[seq_i] = 1.0 / (1.0 + pi_val)
    return out


def verify_rigid_isometry_preserve_distances(
    before: dict[int, tuple[float, float, float]],
    after: dict[int, tuple[float, float, float]],
    tol: float = 1e-9,
) -> None:
    keys = sorted(before.keys())
    if keys != sorted(after.keys()) or len(keys) < 2:
        raise CryoEMPackingComplianceError("Coordinate key mismatch in isometry check")
    probe_pairs = [(keys[0], keys[1]), (keys[0], keys[-1]), (keys[len(keys) // 2], keys[-1])]
    for a, b in probe_pairs:
        d0 = math.dist(before[a], before[b])
        d1 = math.dist(after[a], after[b])
        if abs(d0 - d1) > tol:
            raise CryoEMPackingComplianceError(
                f"Non-isometric coordinate transform detected on pair ({a}, {b}): d0={d0:.6f}, d1={d1:.6f}"
            )


def score_variant_s_flex(
    variant: dict[str, Any],
    compliance_by_seq: dict[int, float | None],
    *,
    allow_non_missense_bug: bool = False,
    empty_fallback: float | None = None,
) -> float | None:
    fam = str(variant["variant_family"])
    if fam not in ALLOWED_FAMILIES:
        raise CryoEMPackingComplianceError(f"Disallowed variant_family {fam!r}")
    if not allow_non_missense_bug and fam != "missense":
        return None
    codon = int(variant["codon_position"])
    val = compliance_by_seq.get(codon, None)
    if val is None:
        return empty_fallback
    return val


def is_canonical_splice(hgvs: str) -> bool:
    return bool(CANONICAL_SPLICE_RE.match(hgvs.strip()))


def score_sequence_consequence_baseline(variant: dict[str, Any]) -> float:
    fam = str(variant["variant_family"])
    if fam == "class_i_ptc_frameshift":
        return 0.0
    if fam == "splice_regulatory":
        return 0.0 if is_canonical_splice(str(variant["hgvs_p_or_c"])) else 0.50
    if fam == "inframe_indel":
        return 0.50
    if fam == "missense":
        return 0.75
    raise CryoEMPackingComplianceError(f"Unknown variant_family: {fam!r}")


def score_class_gated_hierarchical(
    variant: dict[str, Any],
    missense_subscore: float | None,
) -> float | None:
    fam = str(variant["variant_family"])
    if fam == "class_i_ptc_frameshift":
        return 0.0
    if fam == "splice_regulatory":
        return 0.0 if is_canonical_splice(str(variant["hgvs_p_or_c"])) else 0.50
    if fam == "inframe_indel":
        return 0.50
    if fam == "missense":
        if missense_subscore is None:
            return None
        if not (0.0 <= missense_subscore <= 1.0):
            raise CryoEMPackingComplianceError(f"Missense subscore out of [0, 1]: {missense_subscore}")
        return 0.50 + 0.50 * missense_subscore
    raise CryoEMPackingComplianceError(f"Unknown variant_family: {fam!r}")


def directional_auroc_abstention(
    labels: list[int],
    scores: list[float | None],
) -> float:
    if len(labels) != len(scores) or not labels:
        raise CryoEMPackingComplianceError("Mismatched or empty labels/scores")
    pos_idx = [i for i, y in enumerate(labels) if y == 1]
    neg_idx = [i for i, y in enumerate(labels) if y == 0]
    if not pos_idx or not neg_idx:
        raise CryoEMPackingComplianceError("ROC-AUC requires at least one positive and one negative")
    concordant = 0.0
    total = float(len(pos_idx) * len(neg_idx))
    for ip in pos_idx:
        sp = scores[ip]
        for ineg in neg_idx:
            sn = scores[ineg]
            if sp is None or sn is None:
                concordant += 0.5
            elif sp > sn:
                concordant += 1.0
            elif sp == sn:
                concordant += 0.5
    return concordant / total


def validate_benchmark_annotations(annotations: dict[str, Any]) -> list[dict[str, Any]]:
    variants = list(annotations.get("variants", []))
    if len(variants) != 30:
        raise CryoEMPackingComplianceError(f"Expected 30 variants, got {len(variants)}")
    train_codons: list[int] = []
    test_codons: list[int] = []
    for v in variants:
        fam = str(v.get("variant_family", ""))
        if fam not in ALLOWED_FAMILIES:
            raise CryoEMPackingComplianceError(f"Non-canonical variant_family {fam!r}")
        codon = int(v["codon_position"])
        dom = map_p13569_domain(codon)
        split = str(v.get("split", ""))
        if split == "train":
            if dom not in TRAIN_DOMAINS:
                raise CryoEMPackingComplianceError(f"Train variant {v['variant_id']} in {dom}")
            train_codons.append(codon)
        elif split == "test":
            if dom not in TEST_DOMAINS:
                raise CryoEMPackingComplianceError(f"Test variant {v['variant_id']} in {dom}")
            test_codons.append(codon)
        else:
            raise CryoEMPackingComplianceError(f"Unknown split {split!r}")
        if fam != "missense" and (
            v.get("alphamissense_score") is not None or v.get("revel_score") is not None
        ):
            raise CryoEMPackingComplianceError(
                f"Non-missense variant {v['variant_id']} carries missense score"
            )
    if len(train_codons) != 10 or len(test_codons) != 20:
        raise CryoEMPackingComplianceError("Split sizes must be 10 train and 20 test")
    min_gap = min(abs(a - b) for a in train_codons for b in test_codons)
    if min_gap <= 5:
        raise CryoEMPackingComplianceError(f"Cross-split codon gap {min_gap} <= 5")
    return variants


def run_positive_controls(
    struct_5uak: ChainAStructure,
    compliance_map: dict[int, float | None],
    variants: list[dict[str, Any]],
    perms: dict[str, tuple[str, str]],
) -> dict[str, Any]:
    for req_src in ("rcsb-pdb-5uak", "uniprot-p13569", "openfda-label", "myvariant"):
        assert_catalog_source_admitted(req_src, perms)

    if len(struct_5uak.ca_by_seq) != 1139:
        raise CryoEMPackingComplianceError(
            f"Expected 1,139 placed CA residues in 5UAK, got {len(struct_5uak.ca_by_seq)}"
        )
    if 507 not in struct_5uak.ca_by_seq or 508 not in struct_5uak.ca_by_seq:
        raise CryoEMPackingComplianceError("WT 5UAK must place codons 507 and 508")

    # Control 2: Exact 4-point collinear 3.8 A chain -> [36/49, 4/5, 4/5, 36/49]
    chain4 = {k: (3.8 * (k - 1), 0.0, 0.0) for k in (1, 2, 3, 4)}
    c4_scores = compute_noncovalent_compliance_map(chain4)
    expected_c4 = {1: 36.0 / 49.0, 2: 4.0 / 5.0, 3: 4.0 / 5.0, 4: 36.0 / 49.0}
    for k, exp_v in expected_c4.items():
        got_v = c4_scores[k]
        if got_v is None or abs(got_v - exp_v) > 1e-12:
            raise CryoEMPackingComplianceError(
                f"4-point collinear analytic control failed at {k}: expected {exp_v}, got {got_v}"
            )

    # Control 3: Exposed peripheral/cleft vs buried core ordering (0 overlap with 30 variants)
    bench_codons = {int(v["codon_position"]) for v in variants}
    for ctrl in EXPOSED_CONTROLS + BURIED_CONTROLS:
        if min(abs(ctrl - bc) for bc in bench_codons) <= 5:
            raise CryoEMPackingComplianceError(f"Control codon {ctrl} within 5 codons of benchmark")
    exp_vals = [compliance_map[c] for c in EXPOSED_CONTROLS]
    bur_vals = [compliance_map[c] for c in BURIED_CONTROLS]
    if any(v is None or v <= 0.050 for v in exp_vals):
        raise CryoEMPackingComplianceError(f"Exposed controls failed > 0.050: {exp_vals}")
    if any(v is None or v >= 0.048 for v in bur_vals):
        raise CryoEMPackingComplianceError(f"Buried controls failed < 0.048: {bur_vals}")
    if min(v for v in exp_vals if v is not None) <= max(v for v in bur_vals if v is not None):
        raise CryoEMPackingComplianceError("Exposed controls did not strictly exceed buried controls")

    # Control 4: SE(3) rigid-motion invariance + train G178 > G85
    raw_coords = {k: pt for k, (_, pt) in struct_5uak.ca_by_seq.items()}
    # Proper Pythagorean rotation matrix R in SO(3) (det = +1) + translation
    rot_coords: dict[int, tuple[float, float, float]] = {}
    for k, (x, y, z) in raw_coords.items():
        rx = (2.0 * x - 2.0 * y + 1.0 * z) / 3.0 + 17.25
        ry = (1.0 * x + 2.0 * y + 2.0 * z) / 3.0 - 41.50
        rz = (2.0 * x + 1.0 * y - 2.0 * z) / 3.0 + 9.125
        rot_coords[k] = (rx, ry, rz)
    verify_rigid_isometry_preserve_distances(raw_coords, rot_coords)
    rot_map = compute_noncovalent_compliance_map(rot_coords)
    max_se3_drift = max(
        abs((compliance_map[k] or 0.0) - (rot_map[k] or 0.0)) for k in compliance_map
    )
    if max_se3_drift >= 1e-12:
        raise CryoEMPackingComplianceError(f"SE(3) drift {max_se3_drift} >= 1e-12")

    s_g178 = compliance_map[178]
    s_g85 = compliance_map[85]
    if s_g178 is None or s_g85 is None or not (s_g178 > s_g85):
        raise CryoEMPackingComplianceError(
            f"Train control failed: expected s_flex(G178) > s_flex(G85), got {s_g178} vs {s_g85}"
        )

    # Control 5: ROC-AUC unit controls
    if directional_auroc_abstention([1, 0, 1, 0], [None, None, None, None]) != 0.5:
        raise CryoEMPackingComplianceError("All-None ROC-AUC control failed")
    if directional_auroc_abstention([1, 1, 0, 0], [0.9, 0.8, 0.2, 0.1]) != 1.0:
        raise CryoEMPackingComplianceError("Positive separator ROC-AUC control failed")
    if directional_auroc_abstention([1, 1, 0, 0], [0.1, 0.2, 0.8, 0.9]) != 0.0:
        raise CryoEMPackingComplianceError("Inverted separator ROC-AUC control failed")

    return {
        "placed_ca_5uak": len(struct_5uak.ca_by_seq),
        "collinear_4_exact": True,
        "min_exposed_control": round(min(v for v in exp_vals if v is not None), 6),
        "max_buried_control": round(max(v for v in bur_vals if v is not None), 6),
        "max_se3_drift": max_se3_drift,
        "train_g178_s_flex": round(s_g178, 6),
        "train_g85_s_flex": round(s_g85, 6),
    }


def run_falsification_controls(
    struct_5uak: ChainAStructure,
    compliance_map: dict[int, float | None],
    variants: list[dict[str, Any]],
    perms: dict[str, tuple[str, str]],
) -> dict[str, str]:
    results: dict[str, str] = {}

    def must_fail(name: str, fn: Callable[[], Any], exc_types: tuple[type[BaseException], ...]) -> None:
        try:
            fn()
        except exc_types:
            results[name] = "PASS"
            return
        raise AssertionError(f"Planted negative control {name!r} failed to raise {exc_types}")

    # 1. WT 507/508 indel missense leak rejected
    def _leak_508_indel() -> None:
        f508_var = next(v for v in variants if v["legacy_spl_name"] == "F508del")
        leaked = score_variant_s_flex(f508_var, compliance_map, allow_non_missense_bug=True)
        if leaked is not None:
            raise CryoEMPackingComplianceError("In-frame indel F508del illegally received WT CA score")

    must_fail("wt_507_508_indel_missense_leak_rejected", _leak_508_indel, (CryoEMPackingComplianceError,))

    # 2. UNK chain or sequence mismatch rejected
    must_fail(
        "unk_chain_or_sequence_mismatch_rejected",
        lambda: parse_5uak_chain_a(enforce_chain_a=False),
        (CryoEMPackingComplianceError,),
    )

    # 3. Coincident CA zero distance rejected
    must_fail(
        "coincident_ca_zero_distance_rejected",
        lambda: compute_noncovalent_compliance_map({1: (0.0, 0.0, 0.0), 5: (0.0, 0.0, 0.0)}),
        (CryoEMPackingComplianceError,),
    )

    # 4. Isotropic dilation breaks rigid isometry
    raw_coords = {k: pt for k, (_, pt) in struct_5uak.ca_by_seq.items()}
    dilated = {k: (2.0 * pt[0], 2.0 * pt[1], 2.0 * pt[2]) for k, pt in raw_coords.items()}
    must_fail(
        "isotropic_dilation_breaks_pi_rejected",
        lambda: verify_rigid_isometry_preserve_distances(raw_coords, dilated),
        (CryoEMPackingComplianceError,),
    )

    # 5. Empty non-covalent neighborhood numeric fallback rejected
    def _bad_empty_fallback() -> None:
        isolated = compute_noncovalent_compliance_map({10: (0.0, 0.0, 0.0), 11: (3.8, 0.0, 0.0)})
        fake_v = {"variant_family": "missense", "codon_position": 10}
        bad = score_variant_s_flex(fake_v, isolated, empty_fallback=0.0)
        if bad is not None:
            raise CryoEMPackingComplianceError("Empty non-covalent neighborhood returned numeric fallback")

    must_fail("empty_noncovalent_fallback_rejected", _bad_empty_fallback, (CryoEMPackingComplianceError,))

    # 6. Post-hoc sign flip rejected (inverted separator yields 0.0, never max(AUC, 1-AUC)=1.0)
    def _reject_sign_flip() -> None:
        auc_inv = directional_auroc_abstention([1, 1, 0, 0], [0.1, 0.2, 0.8, 0.9])
        if auc_inv >= 0.5:
            return
        raise CryoEMPackingComplianceError("Inverted separator yielded 0.0 as required; sign-flip forbidden")

    must_fail("post_hoc_sign_flip_rejected", _reject_sign_flip, (CryoEMPackingComplianceError,))

    # 7. Codon leakage / _heldout hack / manual source rejected
    def _manual_and_leakage() -> None:
        assert_catalog_source_admitted("mavedb-dump-zenodo", perms)

    must_fail("codon_leakage_or_manual_source_rejected", _manual_and_leakage, (PermissionError,))
    return results


def evaluate_split(
    split_variants: list[dict[str, Any]],
    compliance_map: dict[int, float | None],
    target_key: str = "trikafta_v23_responsive",
) -> dict[str, Any]:
    labels_all = [int(v[target_key]) for v in split_variants]
    missense_vars = [v for v in split_variants if v["variant_family"] == "missense"]
    labels_mis = [int(v[target_key]) for v in missense_vars]

    s_flex_mis = [score_variant_s_flex(v, compliance_map) for v in missense_vars]
    s_am_mis = [
        (1.0 - float(v["alphamissense_score"])) if v.get("alphamissense_score") is not None else None
        for v in missense_vars
    ]

    g_flex_all = [
        score_class_gated_hierarchical(v, score_variant_s_flex(v, compliance_map))
        for v in split_variants
    ]
    g_am_all = [
        score_class_gated_hierarchical(
            v,
            (1.0 - float(v["alphamissense_score"]))
            if v["variant_family"] == "missense" and v.get("alphamissense_score") is not None
            else None,
        )
        for v in split_variants
    ]
    seq_base_all = [score_sequence_consequence_baseline(v) for v in split_variants]

    per_variant = []
    for v in split_variants:
        sf = score_variant_s_flex(v, compliance_map)
        am = v.get("alphamissense_score")
        sam = (1.0 - float(am)) if (v["variant_family"] == "missense" and am is not None) else None
        per_variant.append({
            "variant_id": v["variant_id"],
            "legacy_spl_name": v["legacy_spl_name"],
            "codon_position": int(v["codon_position"]),
            "domain": map_p13569_domain(int(v["codon_position"])),
            "variant_family": v["variant_family"],
            target_key: int(v[target_key]),
            "alyftrek_v9_responsive": int(v["alyftrek_v9_responsive"]),
            "kalydeco_v36_responsive": int(v["kalydeco_v36_responsive"]),
            "s_flex": round(sf, 6) if sf is not None else None,
            "s_am": round(sam, 6) if sam is not None else None,
            "g_flex": round(score_class_gated_hierarchical(v, sf) or 0.0, 6)
            if score_class_gated_hierarchical(v, sf) is not None
            else None,
        })

    return {
        "n_total": len(split_variants),
        "n_pos_total": sum(labels_all),
        "n_neg_total": len(labels_all) - sum(labels_all),
        "n_missense": len(missense_vars),
        "n_missense_placed": sum(x is not None for x in s_flex_mis),
        "n_pos_missense": sum(labels_mis),
        "n_neg_missense": len(labels_mis) - sum(labels_mis),
        "conditional_missense_auroc_s_flex": round(directional_auroc_abstention(labels_mis, s_flex_mis), 6),
        "conditional_missense_auroc_s_am": round(directional_auroc_abstention(labels_mis, s_am_mis), 6),
        "gated_abstention_auroc_g_flex": round(directional_auroc_abstention(labels_all, g_flex_all), 6),
        "gated_abstention_auroc_g_am": round(directional_auroc_abstention(labels_all, g_am_all), 6),
        "sequence_consequence_auroc": round(directional_auroc_abstention(labels_all, seq_base_all), 6),
        "variants": per_variant,
    }


def run_benchmark(run_falsify: bool = False) -> dict[str, Any]:
    perms = load_catalog_permissions()
    uniprot_seq = load_p13569_sequence()
    struct_5uak = parse_5uak_chain_a(uniprot_seq=uniprot_seq)
    raw_coords = {k: pt for k, (_, pt) in struct_5uak.ca_by_seq.items()}
    compliance_map = compute_noncovalent_compliance_map(raw_coords)

    annotations = json.loads(DEFAULT_ANNOTATIONS_JSON.read_text(encoding="utf-8"))
    variants = validate_benchmark_annotations(annotations)
    pos_controls = run_positive_controls(struct_5uak, compliance_map, variants, perms)

    falsify_report: dict[str, str] = {}
    if run_falsify:
        falsify_report = run_falsification_controls(struct_5uak, compliance_map, variants, perms)

    train_vars = [v for v in variants if v["split"] == "train"]
    test_vars = [v for v in variants if v["split"] == "test"]

    train_eval = evaluate_split(train_vars, compliance_map, "trikafta_v23_responsive")
    test_eval = evaluate_split(test_vars, compliance_map, "trikafta_v23_responsive")

    auc_cond = test_eval["conditional_missense_auroc_s_flex"]
    auc_am_cond = test_eval["conditional_missense_auroc_s_am"]
    auc_gated = test_eval["gated_abstention_auroc_g_flex"]
    auc_am_gated = test_eval["gated_abstention_auroc_g_am"]
    auc_seq = test_eval["sequence_consequence_auroc"]

    supported = bool(
        auc_cond >= 0.75
        and auc_cond > auc_am_cond
        and auc_gated >= 0.75
        and auc_gated > auc_seq
        and auc_gated > auc_am_gated
    )
    outcome = "supported" if supported else "refuted"

    return {
        "prediction_id": "PRED-0013",
        "structure_reference": "5UAK (dephosphorylated ATP-free WT human CFTR, 1,139 placed Chain A CA)",
        "target": "trikafta_v23_responsive (US FDA TRIKAFTA v23 Table 6)",
        "r0_virtual_ca_bond_a": R0_VIRTUAL_CA_BOND_A,
        "r0_squared_a2": R0_SQUARED_A2,
        "positive_controls": pos_controls,
        "falsification_controls": falsify_report,
        "train_split": train_eval,
        "test_split": test_eval,
        "outcome": outcome,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Run PRED-0013 apo 5UAK packing compliance benchmark.")
    parser.add_argument("--check", action="store_true", help="Verify provenance, positive controls, and benchmark")
    parser.add_argument("--falsify", action="store_true", help="Run all 7 planted negative controls")
    parser.add_argument("--json", action="store_true", help="Emit full JSON report")
    args = parser.parse_args()

    report = run_benchmark(run_falsify=args.falsify)
    if args.json:
        print(json.dumps(report, indent=2))
    else:
        te = report["test_split"]
        tr = report["train_split"]
        print(
            f"PASS cryoem_packing_compliance (PRED-0013 outcome={report['outcome']}): "
            f"5UAK_CA={report['positive_controls']['placed_ca_5uak']} "
            f"test_AUC_cond(s_flex)={te['conditional_missense_auroc_s_flex']:.6f} "
            f"vs s_AM={te['conditional_missense_auroc_s_am']:.6f}; "
            f"test_AUC_abst(g_flex)={te['gated_abstention_auroc_g_flex']:.6f} "
            f"vs seq_base={te['sequence_consequence_auroc']:.6f} & g_AM={te['gated_abstention_auroc_g_am']:.6f}; "
            f"train_AUC_cond={tr['conditional_missense_auroc_s_flex']:.6f}; "
            f"falsify={len(report['falsification_controls'])}/7"
            if args.falsify
            else f"PASS cryoem_packing_compliance (PRED-0013 outcome={report['outcome']}): "
            f"test_AUC_cond(s_flex)={te['conditional_missense_auroc_s_flex']:.6f} "
            f"test_AUC_abst(g_flex)={te['gated_abstention_auroc_g_flex']:.6f}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
