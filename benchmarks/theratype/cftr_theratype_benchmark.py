#!/usr/bin/env python3
"""Open CFTR Theratype Benchmark (`T-0925` / `PRED-0010`) — Fixed-Split, Real-Data Evaluator.

100% Pure Python Standard Library Runtime. Zero network calls at test time. Zero patient/PHI data.
Joins recorded open variant-level annotations from admitted `sources/catalog.yaml` entries:
  - `myvariant` (`dbnsfp.alphamissense` & `dbnsfp.revel` on canonical `ENST00000003084` / `P13569`, `cadd.phred`)
  - `clinvar` / `ncbi-eutils` (`ClinVar Variation IDs`)
  - `uniprot-p13569` (`CFTR` 6-region structural domain architecture on 1..1480 aa)
  - `europe-pmc` (`PMC10810519`, McDonald et al., PLOS ONE 2024, DOI `10.1371/journal.pone.0297560`, Table 1 VUS controls)
  - `openfda-label` (US FDA Structured Product Labels: `TRIKAFTA` v23 Table 6, `ALYFTREK` v9 Table 5,
    `KALYDECO` v36 Table 3, `SYMDEKO` v22 Table 6)

Tests `PRED-0010` under the lab's Discovery Protocol (`docs/lab-focus-2026-11-14.md`):
whether AlphaMissense (`s_AM = 1 - AM_score`) predicts CFTR modulator responsiveness across
missense, in-frame indel, splicing, and Class I nonsense/frameshift variants on a domain- and
position-disjoint held-out benchmark without post-hoc sign flipping or synthetic scores.

Research benchmark only; not clinical or diagnostic guidance.
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import json
import math
from pathlib import Path
import re
import sys
from typing import Any, Dict, List, Sequence, Tuple

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
DEFAULT_FIXTURE_PATH = HERE / "theratype_open_annotations_2026-11-21.json"
DEFAULT_CATALOG_PATH = REPO_ROOT / "sources" / "catalog.yaml"

SPLIT_VERSION = "2026-11-21.v2"
RESEARCH_BOUNDARY = "Research benchmark only; not clinical or diagnostic guidance."

ALLOWED_FAMILIES: frozenset[str] = frozenset(
    {
        "missense",
        "inframe_indel",
        "splice_regulatory",
        "class_i_ptc_frameshift",
    }
)

P13569_DOMAIN_SPANS: Tuple[Tuple[str, int, int], ...] = (
    ("Lasso_N_terminal", 1, 80),
    ("TMD1", 81, 380),
    ("NBD1", 381, 630),
    ("R_domain", 631, 845),
    ("TMD2", 846, 1190),
    ("NBD2_C_terminal", 1191, 1480),
)

MIN_CROSS_SPLIT_CODON_BUFFER = 5
CANONICAL_SPLICE_DINUCLEOTIDE_RE = re.compile(r"^c\.\d+[+-][12][ACGT]>[ACGT]$")


class TheratypeBenchmarkError(ValueError):
    """Raised when a positive control, provenance invariant, or catalog gate fails."""


@dataclass(frozen=True)
class CftrVariantCase:
    variant_id: str
    hgvs_p_or_c: str
    codon_position: int
    variant_family: str
    split: str
    modulator_responsive: int  # 1 = responsive on US FDA TRIKAFTA v23 Table 6, 0 = non-responsive
    alphamissense_score: float | None
    evo2_zero_shot_score: float | None = None
    alphagenome_avi_score: float | None = None
    revel_score: float | None = None
    cadd_phred: float | None = None
    clinvar_variation_id: int | None = None
    alyftrek_v9_responsive: int = 0
    kalydeco_v36_responsive: int = 0
    symdeko_v22_responsive: int = 0
    structural_domain: str = ""


@dataclass(frozen=True)
class McDonaldVusControl:
    short_name: str
    hgvs_p: str
    codon_position: int
    myvariant_am_score: float
    myvariant_am_pred: str
    expected_table1_score: float
    expected_table1_class: str


def map_codon_to_p13569_domain(codon_position: int) -> str:
    """Deterministically map a 1-based CFTR codon (1..1480) to its UniProt P13569 structural domain."""
    if not isinstance(codon_position, int) or isinstance(codon_position, bool) or not (1 <= codon_position <= 1480):
        raise TheratypeBenchmarkError(f"CFTR codon_position out of bounds [1, 1480]: {codon_position!r}")
    for domain_id, start_aa, end_aa in P13569_DOMAIN_SPANS:
        if start_aa <= codon_position <= end_aa:
            return domain_id
    raise TheratypeBenchmarkError(f"Unmapped CFTR codon_position: {codon_position!r}")


def classify_alphamissense_threshold(score: float) -> str:
    """Apply AlphaMissense CFTR thresholds reported in Cheng et al. 2023 / McDonald et al. 2024."""
    if score < 0.34:
        return "benign"
    if score <= 0.56:
        return "ambiguous"
    return "pathogenic"


def load_benchmark_fixture(
    fixture_path: Path = DEFAULT_FIXTURE_PATH,
) -> Tuple[Tuple[McDonaldVusControl, ...], Tuple[CftrVariantCase, ...], Dict[str, Any]]:
    """Load recorded open annotations and construct typed VUS controls and benchmark cases."""
    payload = json.loads(fixture_path.read_text(encoding="utf-8"))
    controls: List[McDonaldVusControl] = []
    for row in payload.get("mcdonald_2024_table1_positive_controls", []):
        controls.append(
            McDonaldVusControl(
                short_name=str(row["short_name"]),
                hgvs_p=str(row["hgvs_p"]),
                codon_position=int(row["codon_position"]),
                myvariant_am_score=float(row["myvariant_am_score"]),
                myvariant_am_pred=str(row["myvariant_am_pred"]),
                expected_table1_score=float(row["mcdonald_2024_table1_score"]),
                expected_table1_class=str(row["mcdonald_2024_table1_class"]),
            )
        )

    cases: List[CftrVariantCase] = []
    for row in payload.get("variants", []):
        codon = int(row["codon_position"])
        dom = map_codon_to_p13569_domain(codon)
        am_raw = row.get("alphamissense_score")
        rev_raw = row.get("revel_score")
        cadd_raw = row.get("cadd_phred")
        evo_raw = row.get("evo2_zero_shot_score")
        ag_raw = row.get("alphagenome_avi_score")
        cv_raw = row.get("clinvar_variation_id")
        cases.append(
            CftrVariantCase(
                variant_id=str(row["variant_id"]),
                hgvs_p_or_c=str(row["hgvs_p_or_c"]),
                codon_position=codon,
                variant_family=str(row["variant_family"]),
                split=str(row["split"]),
                modulator_responsive=int(row["trikafta_v23_responsive"]),
                alphamissense_score=float(am_raw) if am_raw is not None else None,
                evo2_zero_shot_score=float(evo_raw) if evo_raw is not None else None,
                alphagenome_avi_score=float(ag_raw) if ag_raw is not None else None,
                revel_score=float(rev_raw) if rev_raw is not None else None,
                cadd_phred=float(cadd_raw) if cadd_raw is not None else None,
                clinvar_variation_id=int(cv_raw) if cv_raw is not None else None,
                alyftrek_v9_responsive=int(row["alyftrek_v9_responsive"]),
                kalydeco_v36_responsive=int(row["kalydeco_v36_responsive"]),
                symdeko_v22_responsive=int(row["symdeko_v22_responsive"]),
                structural_domain=dom,
            )
        )
    return tuple(controls), tuple(cases), payload


MCDONALD_2024_TABLE1_POSITIVE_CONTROLS, OPEN_CFTR_THERATYPE_BENCHMARK, _FIXTURE_META = load_benchmark_fixture()


def _compute_pairwise_auroc(
    scored_labels: Sequence[Tuple[float | None, int]],
) -> Dict[str, Any]:
    """Compute exact Mann-Whitney U conditional ROC-AUC and abstention-aware ROC-AUC (no sign flip)."""
    pos_all = [s for s, y in scored_labels if y == 1]
    neg_all = [s for s, y in scored_labels if y == 0]
    n_total = len(scored_labels)
    n_scored = sum(1 for s, _ in scored_labels if s is not None)
    coverage = n_scored / max(n_total, 1)

    if not pos_all or not neg_all:
        return {
            "total_cases": n_total,
            "scored_cases": n_scored,
            "coverage": round(coverage, 4),
            "coverage_raw": coverage,
            "conditional_auroc": None,
            "conditional_auroc_raw": None,
            "abstention_aware_auroc": None,
            "abstention_aware_auroc_raw": None,
            "youden_j_eff": None,
        }

    pos_scored = [s for s in pos_all if s is not None]
    neg_scored = [s for s in neg_all if s is not None]

    if pos_scored and neg_scored:
        cond_wins = 0.0
        for sp in pos_scored:
            for sn in neg_scored:
                if sp > sn:
                    cond_wins += 1.0
                elif sp == sn:
                    cond_wins += 0.5
        cond_auroc: float | None = cond_wins / (len(pos_scored) * len(neg_scored))
    else:
        cond_auroc = None

    total_wins = 0.0
    total_pairs = len(pos_all) * len(neg_all)
    for sp in pos_all:
        for sn in neg_all:
            if sp is None or sn is None:
                total_wins += 0.5
            elif sp > sn:
                total_wins += 1.0
            elif sp == sn:
                total_wins += 0.5
    abst_auroc = total_wins / total_pairs
    j_eff = 2.0 * abst_auroc - 1.0

    return {
        "total_cases": n_total,
        "scored_cases": n_scored,
        "coverage": round(coverage, 4),
        "coverage_raw": coverage,
        "conditional_auroc": round(cond_auroc, 4) if cond_auroc is not None else None,
        "conditional_auroc_raw": cond_auroc,
        "abstention_aware_auroc": round(abst_auroc, 4),
        "abstention_aware_auroc_raw": abst_auroc,
        "youden_j_eff": round(j_eff, 4),
    }


def verify_positive_controls(
    controls: Sequence[McDonaldVusControl] = MCDONALD_2024_TABLE1_POSITIVE_CONTROLS,
) -> Dict[str, Any]:
    """Verify 7 McDonald et al. 2024 Table 1 VUS scores and exact ROC-AUC directional sanity checks."""
    if len(controls) != 7:
        raise TheratypeBenchmarkError(f"Expected 7 McDonald 2024 Table 1 controls, got {len(controls)}")

    pred_code_map = {"B": "benign", "A": "ambiguous", "P": "pathogenic"}
    verified: List[Dict[str, Any]] = []
    for ctrl in controls:
        map_codon_to_p13569_domain(ctrl.codon_position)
        rounded = round(ctrl.myvariant_am_score, 2)
        threshold_class = classify_alphamissense_threshold(ctrl.myvariant_am_score)
        dbnsfp_class = pred_code_map.get(ctrl.myvariant_am_pred, "unknown")
        if rounded != ctrl.expected_table1_score:
            raise TheratypeBenchmarkError(
                f"Positive control {ctrl.short_name} score mismatch: rounded={rounded} != {ctrl.expected_table1_score}"
            )
        if threshold_class != ctrl.expected_table1_class or dbnsfp_class != ctrl.expected_table1_class:
            raise TheratypeBenchmarkError(
                f"Positive control {ctrl.short_name} class mismatch: {threshold_class}/{dbnsfp_class} != {ctrl.expected_table1_class}"
            )
        verified.append(
            {
                "short_name": ctrl.short_name,
                "codon_position": ctrl.codon_position,
                "rounded_score": rounded,
                "class": threshold_class,
            }
        )

    # Exact ROC-AUC directional & abstention sanity checks (forbids post-hoc sign flipping)
    all_none = _compute_pairwise_auroc([(None, 1), (None, 0), (None, 1), (None, 0)])
    perfect = _compute_pairwise_auroc([(0.9, 1), (0.8, 1), (0.2, 0), (0.1, 0)])
    inverted = _compute_pairwise_auroc([(0.1, 1), (0.2, 1), (0.8, 0), (0.9, 0)])
    if all_none["abstention_aware_auroc"] != 0.5:
        raise TheratypeBenchmarkError("Abstention ROC-AUC failed all-None == 0.5 invariant")
    if perfect["conditional_auroc"] != 1.0 or inverted["conditional_auroc"] != 0.0:
        raise TheratypeBenchmarkError("Directional ROC-AUC failed strict no-sign-flip invariant")

    return {
        "status": "PASS",
        "controls_verified": len(verified),
        "controls": verified,
    }


def _is_non_missense_case(case: CftrVariantCase) -> bool:
    hgvs = case.hgvs_p_or_c.strip()
    return (
        case.variant_family != "missense"
        or "Ter" in hgvs
        or "fs" in hgvs
        or "*" in hgvs
        or "del" in hgvs
        or "dup" in hgvs
        or "ins" in hgvs
        or hgvs.startswith("c.")
    )


def audit_split_leakage(
    cases: Sequence[CftrVariantCase],
    controls: Sequence[McDonaldVusControl] = MCDONALD_2024_TABLE1_POSITIVE_CONTROLS,
    min_codon_buffer: int = MIN_CROSS_SPLIT_CODON_BUFFER,
) -> Dict[str, Any]:
    """Verify codon, 5-aa boundary buffer, P13569 structural-domain, family allowlist, and score-provenance hygiene."""
    invalid_families: List[str] = []
    domain_mismatches: List[str] = []
    provenance_violations: List[str] = []
    out_of_bounds_codons: List[str] = []

    recomputed_domains: Dict[str, set[str]] = {"train": set(), "test": set()}
    for c in cases:
        if c.split not in ("train", "test"):
            provenance_violations.append(f"{c.variant_id}:invalid_split:{c.split}")
            continue
        if c.variant_family not in ALLOWED_FAMILIES or c.variant_family.endswith("_heldout"):
            invalid_families.append(f"{c.variant_id}:{c.variant_family}")
        try:
            real_dom = map_codon_to_p13569_domain(c.codon_position)
            recomputed_domains[c.split].add(real_dom)
            if c.structural_domain and c.structural_domain != real_dom:
                domain_mismatches.append(f"{c.variant_id}:{c.structural_domain}!={real_dom}")
        except TheratypeBenchmarkError:
            out_of_bounds_codons.append(f"{c.variant_id}:{c.codon_position}")

        for reg_name, reg_val in (
            ("trikafta_v23", c.modulator_responsive),
            ("alyftrek_v9", c.alyftrek_v9_responsive),
            ("kalydeco_v36", c.kalydeco_v36_responsive),
            ("symdeko_v22", c.symdeko_v22_responsive),
        ):
            if isinstance(reg_val, bool) or reg_val not in (0, 1):
                provenance_violations.append(f"{c.variant_id}:invalid_{reg_name}_label:{reg_val!r}")

        if _is_non_missense_case(c):
            if c.alphamissense_score is not None:
                provenance_violations.append(f"{c.variant_id}:non_missense_has_alphamissense")
            if c.revel_score is not None:
                provenance_violations.append(f"{c.variant_id}:non_missense_has_revel")
        else:
            for sc_name, sc_val in (("alphamissense", c.alphamissense_score), ("revel", c.revel_score)):
                if sc_val is not None and (not math.isfinite(sc_val) or not (0.0 <= sc_val <= 1.0)):
                    provenance_violations.append(f"{c.variant_id}:out_of_range_{sc_name}:{sc_val!r}")
        if c.cadd_phred is not None and (not math.isfinite(c.cadd_phred) or c.cadd_phred < 0.0):
            provenance_violations.append(f"{c.variant_id}:out_of_range_cadd:{c.cadd_phred!r}")
        if c.evo2_zero_shot_score is not None or c.alphagenome_avi_score is not None:
            provenance_violations.append(f"{c.variant_id}:uncataloged_foundation_score_present")

    train_pos = {c.codon_position for c in cases if c.split == "train"}
    test_pos = {c.codon_position for c in cases if c.split == "test"}
    pos_overlap = sorted(train_pos & test_pos)

    control_pos = {ctrl.codon_position for ctrl in controls}
    control_overlap = sorted((train_pos | test_pos) & control_pos)

    if train_pos and test_pos and not out_of_bounds_codons:
        min_gap = min(abs(tp - qp) for tp in train_pos for qp in test_pos)
    else:
        min_gap = 0
    buffer_violations = [
        f"train:{tp}_test:{qp}"
        for tp in sorted(train_pos)
        for qp in sorted(test_pos)
        if abs(tp - qp) <= min_codon_buffer
    ]

    domain_overlap = sorted(recomputed_domains["train"] & recomputed_domains["test"])

    train_fam = {c.variant_family for c in cases if c.split == "train"}
    test_fam = {c.variant_family for c in cases if c.split == "test"}
    heldout_only_families = sorted(test_fam - train_fam)

    passed = (
        len(pos_overlap) == 0
        and len(control_overlap) == 0
        and len(buffer_violations) == 0
        and len(domain_overlap) == 0
        and len(invalid_families) == 0
        and len(domain_mismatches) == 0
        and len(provenance_violations) == 0
        and len(out_of_bounds_codons) == 0
    )
    return {
        "leakage_audit_passed": passed,
        "train_count": sum(1 for c in cases if c.split == "train"),
        "test_count": sum(1 for c in cases if c.split == "test"),
        "codon_position_overlap": pos_overlap,
        "positive_control_codon_overlap": control_overlap,
        "min_cross_split_codon_gap": min_gap,
        "codon_buffer_violations": buffer_violations,
        "structural_domain_overlap": domain_overlap,
        "train_domains": sorted(recomputed_domains["train"]),
        "test_domains": sorted(recomputed_domains["test"]),
        "heldout_novel_families_in_test": heldout_only_families,
        "invalid_families": invalid_families,
        "domain_mismatches": domain_mismatches,
        "provenance_violations": provenance_violations,
        "out_of_bounds_codons": out_of_bounds_codons,
    }


def _score_case_for_responsiveness(case: CftrVariantCase, model_name: str) -> float | None:
    """Pre-registered directional responsiveness score s(x) in [0, 1] (higher = predicted responsive)."""
    if model_name == "alphamissense":
        if case.alphamissense_score is None or not math.isfinite(case.alphamissense_score):
            return None
        return round(1.0 - float(case.alphamissense_score), 6)
    if model_name == "revel":
        if case.revel_score is None or not math.isfinite(case.revel_score):
            return None
        return round(1.0 - float(case.revel_score), 6)
    if model_name == "cadd_phred":
        if case.cadd_phred is None or not math.isfinite(case.cadd_phred):
            return None
        clipped = min(max(float(case.cadd_phred) / 50.0, 0.0), 1.0)
        return round(1.0 - clipped, 6)
    if model_name == "sequence_consequence_baseline":
        if case.variant_family == "class_i_ptc_frameshift":
            return 0.0
        if CANONICAL_SPLICE_DINUCLEOTIDE_RE.match(case.hgvs_p_or_c.strip()):
            return 0.0
        if case.variant_family in ("splice_regulatory", "inframe_indel"):
            return 0.5
        if case.variant_family == "missense":
            return 0.75
        return 0.0
    raise ValueError(f"Unknown model_name: {model_name}")


def evaluate_model_coverage_and_auroc(
    cases: Sequence[CftrVariantCase],
    model_name: str,
    regimen: str = "trikafta_v23",
    split_filter: str | None = "test",
) -> Dict[str, Any]:
    """Evaluate pre-registered responsiveness coverage, conditional ROC-AUC, and abstention-aware ROC-AUC."""
    eval_cases = [c for c in cases if split_filter is None or c.split == split_filter]
    scored_labels: List[Tuple[float | None, int]] = []
    for c in eval_cases:
        if regimen == "trikafta_v23":
            label = int(c.modulator_responsive)
        elif regimen == "alyftrek_v9":
            label = int(c.alyftrek_v9_responsive)
        elif regimen == "kalydeco_v36":
            label = int(c.kalydeco_v36_responsive)
        elif regimen == "symdeko_v22":
            label = int(c.symdeko_v22_responsive)
        else:
            raise ValueError(f"Unknown regimen: {regimen}")
        score = _score_case_for_responsiveness(c, model_name)
        scored_labels.append((score, label))

    metrics = _compute_pairwise_auroc(scored_labels)
    return {
        "model": model_name,
        "regimen": regimen,
        "split": split_filter or "all",
        "test_total": metrics["total_cases"],
        "test_scored": metrics["scored_cases"],
        "test_coverage": metrics["coverage"],
        "test_coverage_raw": metrics["coverage_raw"],
        "conditional_auroc": metrics["conditional_auroc"],
        "conditional_auroc_raw": metrics["conditional_auroc_raw"],
        "abstention_aware_auroc": metrics["abstention_aware_auroc"],
        "abstention_aware_auroc_raw": metrics["abstention_aware_auroc_raw"],
        "youden_j_eff": metrics["youden_j_eff"],
    }


def evaluate_lofo_splits(
    cases: Sequence[CftrVariantCase],
    model_name: str = "alphamissense",
) -> Dict[str, Dict[str, Any]]:
    """Run 4-fold Leave-One-Family-Out (LOFO) evaluation with zero family overlap per fold."""
    out: Dict[str, Dict[str, Any]] = {}
    for fam in sorted(ALLOWED_FAMILIES):
        heldout = [c for c in cases if c.variant_family == fam]
        complement_fams = {c.variant_family for c in cases if c.variant_family != fam}
        scored_labels = [(_score_case_for_responsiveness(c, model_name), int(c.modulator_responsive)) for c in heldout]
        metrics = _compute_pairwise_auroc(scored_labels)
        out[fam] = {
            "heldout_count": len(heldout),
            "family_overlap_with_train": sorted({fam} & complement_fams),
            "coverage": metrics["coverage"],
            "conditional_auroc": metrics["conditional_auroc"],
            "abstention_aware_auroc": metrics["abstention_aware_auroc"],
        }
    return out


def _strip_yaml_scalar(raw_val: str) -> str:
    val = raw_val.split("#", 1)[0].strip()
    if len(val) >= 2 and val[0] == val[-1] and val[0] in ("'", '"'):
        return val[1:-1].strip()
    return val


def _read_catalog_source_status(catalog_path: Path, source_id: str) -> Tuple[str | None, str | None]:
    """Pure-stdlib indentation-aware parser reading `access` and `claim_label` in `sources/catalog.yaml`."""
    if not catalog_path.is_file():
        raise FileNotFoundError(f"Catalog file not found: {catalog_path}")
    current_id: str | None = None
    access: str | None = None
    claim_label: str | None = None
    in_access_block = False
    for raw_line in catalog_path.read_text(encoding="utf-8").splitlines():
        if not raw_line.strip() or raw_line.lstrip().startswith("#"):
            continue
        indent = len(raw_line) - len(raw_line.lstrip(" "))
        stripped = raw_line.strip()
        if indent in (0, 2) and stripped.startswith("- id:"):
            if current_id == source_id:
                return access, claim_label
            current_id = _strip_yaml_scalar(stripped.split(":", 1)[1])
            access = None
            claim_label = None
            in_access_block = False
        elif current_id == source_id:
            if indent == 4:
                in_access_block = False
                if stripped.startswith("access:"):
                    val = _strip_yaml_scalar(stripped.split(":", 1)[1])
                    if val:
                        access = val
                    else:
                        in_access_block = True
                elif stripped.startswith("claim_label:"):
                    claim_label = _strip_yaml_scalar(stripped.split(":", 1)[1])
            elif indent == 6 and in_access_block and stripped.startswith("mode:"):
                access = _strip_yaml_scalar(stripped.split(":", 1)[1])
    if current_id == source_id:
        return access, claim_label
    return None, None


def load_mavedb_cc0_csv(
    csv_path: Path,
    catalog_path: Path = DEFAULT_CATALOG_PATH,
    source_id: str = "mavedb-dump-zenodo",
) -> List[Dict[str, str]]:
    """Load MaveDB CC0 CSV records ONLY when catalog explicitly admits access in {'fetch', 'api'} and claim_label == 'verified'."""
    access, claim_label = _read_catalog_source_status(catalog_path, source_id)
    if access is None:
        raise PermissionError(f"Catalog entry {source_id!r} not found in {catalog_path}")
    if access not in ("fetch", "api") or claim_label != "verified":
        raise PermissionError(
            f"Refusing {source_id!r}: catalog status is access={access!r}, claim_label={claim_label!r} "
            "(strict allowlist requires access in {'fetch', 'api'} and claim_label == 'verified')"
        )
    with csv_path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def evaluate_pred_0010(cases: Sequence[CftrVariantCase] = OPEN_CFTR_THERATYPE_BENCHMARK) -> Dict[str, Any]:
    """Evaluate pre-registered prediction PRED-0010 (AlphaMissense as a CFTR theratype predictor)."""
    verify_positive_controls()
    leakage = audit_split_leakage(cases)
    if not leakage["leakage_audit_passed"]:
        raise TheratypeBenchmarkError(f"Split leakage audit failed: {leakage}")

    am_eti = evaluate_model_coverage_and_auroc(cases, "alphamissense", regimen="trikafta_v23")
    am_vtd = evaluate_model_coverage_and_auroc(cases, "alphamissense", regimen="alyftrek_v9")
    revel_eti = evaluate_model_coverage_and_auroc(cases, "revel", regimen="trikafta_v23")
    cadd_eti = evaluate_model_coverage_and_auroc(cases, "cadd_phred", regimen="trikafta_v23")
    baseline_eti = evaluate_model_coverage_and_auroc(cases, "sequence_consequence_baseline", regimen="trikafta_v23")
    baseline_vtd = evaluate_model_coverage_and_auroc(cases, "sequence_consequence_baseline", regimen="alyftrek_v9")

    # Pre-registered thresholds in done.md:
    # test_coverage >= 0.80 AND conditional_auroc >= 0.75 AND abstention_aware_auroc >= 0.75
    # AND abstention_aware_auroc >= baseline abstention_aware_auroc
    cond_auc = am_eti["conditional_auroc"]
    abst_auc = am_eti["abstention_aware_auroc"]
    base_auc = baseline_eti["abstention_aware_auroc"]
    supported = (
        am_eti["test_coverage"] >= 0.80
        and cond_auc is not None
        and cond_auc >= 0.75
        and abst_auc is not None
        and abst_auc >= 0.75
        and base_auc is not None
        and abst_auc >= base_auc
    )
    return {
        "prediction_id": "PRED-0010",
        "outcome": "supported" if supported else "refuted",
        "alphamissense_trikafta_v23": am_eti,
        "alphamissense_alyftrek_v9": am_vtd,
        "revel_trikafta_v23": revel_eti,
        "cadd_phred_trikafta_v23": cadd_eti,
        "sequence_consequence_baseline_trikafta_v23": baseline_eti,
        "sequence_consequence_baseline_alyftrek_v9": baseline_vtd,
        "lofo_alphamissense": evaluate_lofo_splits(cases, "alphamissense"),
    }


def evaluate_falsification_controls() -> Dict[str, Any]:
    """Run all 6 deterministic poisoned negative controls required by done.md."""
    controls: List[Dict[str, Any]] = []

    # 1. Poisoned codon leakage between train and test must fail audit_split_leakage
    leaked_codon = list(OPEN_CFTR_THERATYPE_BENCHMARK) + [
        CftrVariantCase("POISON-CODON-117", "p.Arg117Cys", 117, "missense", "test", 1, 0.25, structural_domain="TMD1")
    ]
    controls.append({"id": "codon_leakage_detected", "falsified": not audit_split_leakage(leaked_codon)["leakage_audit_passed"]})

    # 2. Poisoned P13569 domain overlap or <= 5-codon boundary buffer breach must fail
    leaked_domain_buffer = list(OPEN_CFTR_THERATYPE_BENCHMARK) + [
        CftrVariantCase("POISON-BUFFER-350", "p.Arg352Gln", 350, "missense", "test", 1, 0.60, structural_domain="TMD1")
    ]
    controls.append(
        {
            "id": "domain_or_buffer_leakage_detected",
            "falsified": not audit_split_leakage(leaked_domain_buffer)["leakage_audit_passed"],
        }
    )

    # 3. Non-canonical `_heldout` family-renaming hack must be rejected
    hack_family = list(OPEN_CFTR_THERATYPE_BENCHMARK) + [
        CftrVariantCase(
            "POISON-HELDOUT-HACK",
            "p.Gly1349Asp",
            1349,
            "gating_missense_heldout",
            "test",
            1,
            0.9801,
            structural_domain="NBD2_C_terminal",
        )
    ]
    controls.append({"id": "heldout_suffix_hack_rejected", "falsified": not audit_split_leakage(hack_family)["leakage_audit_passed"]})

    # 4. Non-None AlphaMissense score on an in-frame deletion (F508del) or PTC must be rejected
    synth_indel = [
        CftrVariantCase(
            c.variant_id,
            c.hgvs_p_or_c,
            c.codon_position,
            c.variant_family,
            c.split,
            c.modulator_responsive,
            0.94 if c.variant_id == "CFTR-F508del" else c.alphamissense_score,
            structural_domain=c.structural_domain,
        )
        for c in OPEN_CFTR_THERATYPE_BENCHMARK
    ]
    controls.append(
        {
            "id": "non_missense_alphamissense_rejected",
            "falsified": not audit_split_leakage(synth_indel)["leakage_audit_passed"],
        }
    )

    # 5. Unadmitted `access: manual` catalog sources (T-0928) must fail closed with PermissionError
    manual_blocked = True
    for manual_id in ("mavedb-dump-zenodo", "cftr-dms-mckee-2023", "cftr-dms-tedman-2025"):
        try:
            load_mavedb_cc0_csv(DEFAULT_FIXTURE_PATH, DEFAULT_CATALOG_PATH, source_id=manual_id)
            manual_blocked = False
        except PermissionError:
            pass
    controls.append({"id": "unadmitted_manual_source_rejected", "falsified": manual_blocked})

    # 6. Unknown model name must raise ValueError
    try:
        evaluate_model_coverage_and_auroc(OPEN_CFTR_THERATYPE_BENCHMARK, "nonexistent_model")
        controls.append({"id": "unknown_model_rejected", "falsified": False})
    except ValueError:
        controls.append({"id": "unknown_model_rejected", "falsified": True})

    all_passed = all(c["falsified"] for c in controls)
    return {
        "status": "PASS" if all_passed else "FAIL",
        "controls_passed": sum(1 for c in controls if c["falsified"]),
        "total_controls": len(controls),
        "controls": controls,
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Open CFTR Theratype Benchmark (T-0925 / PRED-0010)")
    parser.add_argument("--check", action="store_true", help="Run positive controls, leakage audit, and PRED-0010 evaluation")
    parser.add_argument("--falsify", action="store_true", help="Run 6 poisoned negative controls")
    parser.add_argument("--json", action="store_true", help="Emit JSON report")
    args = parser.parse_args(list(argv) if argv is not None else None)

    pos_ctrl = verify_positive_controls()
    leakage = audit_split_leakage(OPEN_CFTR_THERATYPE_BENCHMARK)
    pred = evaluate_pred_0010(OPEN_CFTR_THERATYPE_BENCHMARK)
    falsify = (
        evaluate_falsification_controls()
        if args.falsify
        else {"status": "SKIPPED", "controls_passed": 0, "total_controls": 0, "controls": []}
    )

    am_eti = pred["alphamissense_trikafta_v23"]
    base_eti = pred["sequence_consequence_baseline_trikafta_v23"]
    ok = (
        pos_ctrl["status"] == "PASS"
        and leakage["leakage_audit_passed"]
        and base_eti["test_coverage"] == 1.0
        and am_eti["test_coverage"] < base_eti["test_coverage"]
        and (not args.falsify or falsify["status"] == "PASS")
    )
    payload = {
        "status": "PASS" if ok else "FAIL",
        "split_version": SPLIT_VERSION,
        "boundary": RESEARCH_BOUNDARY,
        "positive_controls": pos_ctrl,
        "leakage_audit": leakage,
        "pred_0010": pred,
        "falsification": falsify,
    }
    if args.json:
        print(json.dumps(payload, indent=2))
    else:
        print(
            f"CFTR Theratype Benchmark ({SPLIT_VERSION}): {payload['status']} "
            f"(train={leakage['train_count']}, test={leakage['test_count']}, "
            f"PRED-0010={pred['outcome']}, "
            f"AM_cov={am_eti['test_coverage']}, AM_cond_AUROC={am_eti['conditional_auroc']}, "
            f"AM_abst_AUROC={am_eti['abstention_aware_auroc']}, "
            f"baseline_AUROC={base_eti['abstention_aware_auroc']})"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
