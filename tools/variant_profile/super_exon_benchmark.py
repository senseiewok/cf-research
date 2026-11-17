#!/usr/bin/env python3
"""RefSeq NM_000492.4 Exon-Span Coverage Benchmark for Partial & Full CFTR cDNA Strategies.

Grounded in authoritative NCBI RefSeq `NM_000492.4` 27-exon CDS nucleotide boundaries
(`CDS 71..4510`, 4,440 coding nucleotides, 1,480 codons) and `CFTR_EXON_SPANS` in
`tools/variant_profile/benchmark_predictors.py`, the 92 canonical Class I CFTR nonsense
alleles (`CANONICAL_92_NONSENSE_ALLELES` in `tools/discovery/round122_...py`), and
Barilla et al. (Dec 2026, DOI 10.1016/j.omtn.2026.103123, PMC13717110, gold open access),
which targeted insertion of partial CFTR cDNA spanning exons 2-27 into CFTR intron 1
(`TI-1`) so native exon 1 splices onto the therapeutic exons 2-27 cassette.

Also evaluates a distal C-terminal exons 23-27 comparison window (`c.3718..4440`, codons
`1240..1480`) and full exons 1-27 coding replacement (`NCT06237335` RCT2100 mRNA,
`NCT05504837` KB407 HSV-1; noted alongside `NCT05248230` 4D-710 `CFTRdeltaR`) with strict
upstream negative controls (`Q2X` in exon 1; `R1162X` in exon 22; `G542X` in exon 12;
`F508del` in exon 11).

Preclinical Research Boundary:
- Sequence-coordinate span coverage measures only whether a variant codon falls within
  an inserted cDNA exon cassette under RefSeq NM_000492.4. It does not model CRISPR
  editing efficiency, airway vector tropism, clinical efficacy, or personal eligibility.

Usage:
    python tools/variant_profile/super_exon_benchmark.py --check
    python tools/variant_profile/super_exon_benchmark.py --falsify
    python tools/variant_profile/super_exon_benchmark.py --json
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

try:
    from tools.variant_profile.benchmark_predictors import CFTR_EXON_SPANS, CFTR_STRUCTURAL_DOMAINS
except ImportError:
    from benchmark_predictors import CFTR_EXON_SPANS, CFTR_STRUCTURAL_DOMAINS


# Authoritative NCBI RefSeq NM_000492.4 CDS nucleotide spans (1-based within CDS 71..4510; stop TAG 4441..4443 excluded)
NM_000492_4_EXON_CDS_NT_SPANS: tuple[tuple[int, int, int], ...] = (
    (1, 1, 53),
    (2, 54, 164),
    (3, 165, 273),
    (4, 274, 489),
    (5, 490, 579),
    (6, 580, 743),
    (7, 744, 869),
    (8, 870, 1116),
    (9, 1117, 1209),
    (10, 1210, 1392),
    (11, 1393, 1584),
    (12, 1585, 1679),
    (13, 1680, 1766),
    (14, 1767, 2490),
    (15, 2491, 2619),
    (16, 2620, 2657),
    (17, 2658, 2908),
    (18, 2909, 2988),
    (19, 2989, 3139),
    (20, 3140, 3367),
    (21, 3368, 3468),
    (22, 3469, 3717),
    (23, 3718, 3873),
    (24, 3874, 3963),
    (25, 3964, 4136),
    (26, 4137, 4242),
    (27, 4243, 4440),
)


def _load_92_nonsense_alleles() -> tuple[tuple[str, int], ...]:
    mod_path = (
        REPO_ROOT
        / "tools"
        / "discovery"
        / "round122_doob_h_transform_92_nonsense_alleles_mrna_survival_bridge.py"
    )
    if mod_path.is_file():
        spec = importlib.util.spec_from_file_location("round122_nonsense", mod_path)
        if spec and spec.loader:
            mod = importlib.util.module_from_spec(spec)
            sys.modules[spec.name] = mod
            spec.loader.exec_module(mod)
            alleles = getattr(mod, "CANONICAL_92_NONSENSE_ALLELES", None)
            if isinstance(alleles, tuple) and len(alleles) == 92:
                return alleles
    raise RuntimeError(f"Could not load CANONICAL_92_NONSENSE_ALLELES from {mod_path}")


CANONICAL_92_NONSENSE_ALLELES = _load_92_nonsense_alleles()

CLINICAL_BOUNDARY = (
    "Preclinical sequence-span arithmetic on RefSeq NM_000492.4 only; "
    "independent open-source research, not medical advice and never an indicator of clinical readiness."
)

# Exon-span regimes tested
REGIME_DEFINITIONS: dict[str, dict[str, Any]] = {
    "intron1_partial_cdna_exons_2_27": {
        "label": "TI-1 Partial cDNA Exons 2-27 (Intron 1 Insertion, Barilla et al. 2026, PMC13717110)",
        "start_exon": 2,
        "end_exon": 27,
        "start_codon": 19,
        "end_codon": 1480,
        "prior_art_doi": "10.1016/j.omtn.2026.103123",
    },
    "distal_super_exon_23_27": {
        "label": "Distal C-Terminal Partial cDNA Exons 23-27 (NBD2/C-tail Comparison Window)",
        "start_exon": 23,
        "end_exon": 27,
        "start_codon": 1240,
        "end_codon": 1480,
        "comparison_basis": "RefSeq NM_000492.4 Exons 23-27 (c.3718..4440, codons 1240..1480) counterfactual comparison window",
    },
    "full_coding_exons_1_27": {
        "label": "Full-Length CFTR Coding Replacement Exons 1-27 (mRNA / Full-Coding Vector Reference)",
        "start_exon": 1,
        "end_exon": 27,
        "start_codon": 1,
        "end_codon": 1480,
        "prior_art_ncts": ["NCT06237335", "NCT05504837", "NCT05248230"],
    },
}

# Sentinel variant test panel (positive & negative locus controls; W1282X, R1162X, G542X, Q2X, F508del, N1303K, S1455X)
SENTINEL_VARIANT_CONTROLS: tuple[dict[str, Any], ...] = (
    {"variant": "Q2X", "codon": 2, "expected_exon": 1, "expected_in_exons_2_27": False, "expected_in_exons_23_27": False},
    {"variant": "F508del", "codon": 508, "expected_exon": 11, "expected_in_exons_2_27": True, "expected_in_exons_23_27": False},
    {"variant": "G542X", "codon": 542, "expected_exon": 12, "expected_in_exons_2_27": True, "expected_in_exons_23_27": False},
    {"variant": "R1162X", "codon": 1162, "expected_exon": 22, "expected_in_exons_2_27": True, "expected_in_exons_23_27": False},
    {"variant": "W1204X", "codon": 1204, "expected_exon": 22, "expected_in_exons_2_27": True, "expected_in_exons_23_27": False},
    {"variant": "W1282X", "codon": 1282, "expected_exon": 23, "expected_in_exons_2_27": True, "expected_in_exons_23_27": True},
    {"variant": "N1303K", "codon": 1303, "expected_exon": 24, "expected_in_exons_2_27": True, "expected_in_exons_23_27": True},
    {"variant": "S1455X", "codon": 1455, "expected_exon": 27, "expected_in_exons_2_27": True, "expected_in_exons_23_27": True},
)


class SuperExonBenchmarkError(ValueError):
    """Raised when RefSeq NM_000492.4 exon boundaries or negative controls fail."""


def codon_to_exon(codon: int) -> int:
    """Map 1-based CFTR codon (1..1480) to RefSeq NM_000492.4 exon number (1..27)."""
    c = int(codon)
    for exon_no, start_aa, end_aa in CFTR_EXON_SPANS:
        if start_aa <= c <= end_aa:
            return exon_no
    raise SuperExonBenchmarkError(f"Codon {codon} lies outside RefSeq NM_000492.4 coding span 1..1480")


def codon_covered_by_exon_range(codon: int, start_exon: int, end_exon: int) -> bool:
    """Return True strictly when codon maps to an exon inside [start_exon, end_exon]."""
    exon_no = codon_to_exon(codon)
    return int(start_exon) <= exon_no <= int(end_exon)


def codon_to_domain(codon: int) -> str:
    """Return CFTR structural domain_id covering codon."""
    c = int(codon)
    for dom in CFTR_STRUCTURAL_DOMAINS:
        if int(dom["start_aa"]) <= c <= int(dom["end_aa"]):
            return str(dom["domain_id"])
    return "unknown"


def verify_nm_000492_4_exon_table() -> None:
    """Verify CFTR_EXON_SPANS matches authoritative NCBI RefSeq NM_000492.4 CDS nucleotide spans."""
    if len(CFTR_EXON_SPANS) != 27 or len(NM_000492_4_EXON_CDS_NT_SPANS) != 27:
        raise SuperExonBenchmarkError("Expected 27 exons in NM_000492.4")

    expected_aa = 1
    expected_nt = 1
    for (ex_aa, start_aa, end_aa), (ex_nt, start_nt, end_nt) in zip(
        CFTR_EXON_SPANS, NM_000492_4_EXON_CDS_NT_SPANS
    ):
        if ex_aa != ex_nt:
            raise SuperExonBenchmarkError(f"Exon index mismatch: {ex_aa} != {ex_nt}")
        if start_nt != expected_nt:
            raise SuperExonBenchmarkError(f"Exon {ex_nt} CDS start {start_nt} != {expected_nt}")
        expected_nt = end_nt + 1

        if start_aa != expected_aa:
            raise SuperExonBenchmarkError(f"Exon {ex_aa} codon start {start_aa} != {expected_aa}")
        expected_aa = end_aa + 1

        # Majority-nucleotide rule: central nucleotide (3*c - 1) of codon c falls in [start_nt, end_nt]
        derived_codons = [c for c in range(1, 1481) if start_nt <= 3 * c - 1 <= end_nt]
        if not derived_codons or (derived_codons[0], derived_codons[-1]) != (start_aa, end_aa):
            raise SuperExonBenchmarkError(
                f"Exon {ex_aa} span {start_aa}..{end_aa} diverges from NM_000492.4 CDS "
                f"c.{start_nt}..{end_nt} ({derived_codons[0]}..{derived_codons[-1]})"
            )

    if expected_aa != 1481 or expected_nt != 4441:
        raise SuperExonBenchmarkError(f"NM_000492.4 end mismatch: aa={expected_aa - 1}, nt={expected_nt - 1}")

    for rkey, cfg in REGIME_DEFINITIONS.items():
        s_ex, e_ex = int(cfg["start_exon"]), int(cfg["end_exon"])
        exp_start = next(sa for ex, sa, _ in CFTR_EXON_SPANS if ex == s_ex)
        exp_end = next(ea for ex, _, ea in CFTR_EXON_SPANS if ex == e_ex)
        if int(cfg["start_codon"]) != exp_start or int(cfg["end_codon"]) != exp_end:
            raise SuperExonBenchmarkError(
                f"Regime {rkey} codon bounds {cfg['start_codon']}..{cfg['end_codon']} "
                f"do not match exon {s_ex}..{e_ex} ({exp_start}..{exp_end})"
            )


def verify_control_panel(panel: tuple[dict[str, Any], ...] = SENTINEL_VARIANT_CONTROLS) -> list[dict[str, Any]]:
    """Verify exact RefSeq exon mapping and negative controls across curated CFTR loci."""
    verified_rows: list[dict[str, Any]] = []
    for item in panel:
        codon = int(item["codon"])
        obs_exon = codon_to_exon(codon)
        cov_2_27 = codon_covered_by_exon_range(codon, 2, 27)
        cov_23_27 = codon_covered_by_exon_range(codon, 23, 27)
        if obs_exon != int(item["expected_exon"]):
            raise SuperExonBenchmarkError(
                f"{item['variant']} codon {codon}: observed exon {obs_exon} != expected exon {item['expected_exon']}"
            )
        if cov_2_27 != bool(item["expected_in_exons_2_27"]):
            raise SuperExonBenchmarkError(
                f"{item['variant']} exon {obs_exon}: exons 2-27 coverage mismatch ({cov_2_27} != {item['expected_in_exons_2_27']})"
            )
        if cov_23_27 != bool(item["expected_in_exons_23_27"]):
            raise SuperExonBenchmarkError(
                f"{item['variant']} exon {obs_exon}: exons 23-27 coverage mismatch ({cov_23_27} != {item['expected_in_exons_23_27']})"
            )
        verified_rows.append(
            {
                "variant": item["variant"],
                "codon": codon,
                "exon": obs_exon,
                "domain": codon_to_domain(codon),
                "covered_exons_2_27": cov_2_27,
                "covered_exons_23_27": cov_23_27,
            }
        )
    return verified_rows


def evaluate_super_exon_coverage() -> dict[str, Any]:
    """Compute deterministic coverage across all 92 canonical Class I nonsense alleles and sentinel loci."""
    verify_nm_000492_4_exon_table()
    controls_verified = verify_control_panel(SENTINEL_VARIANT_CONTROLS)

    regime_summaries: dict[str, Any] = {}
    for rkey, cfg in REGIME_DEFINITIONS.items():
        s_ex = int(cfg["start_exon"])
        e_ex = int(cfg["end_exon"])
        covered: list[str] = []
        excluded: list[dict[str, Any]] = []
        for allele, codon in CANONICAL_92_NONSENSE_ALLELES:
            ex = codon_to_exon(codon)
            if s_ex <= ex <= e_ex:
                covered.append(allele)
            else:
                excluded.append({"allele": allele, "codon": codon, "exon": ex})
        total = len(CANONICAL_92_NONSENSE_ALLELES)
        pct = round(100.0 * len(covered) / total, 2)
        regime_summaries[rkey] = {
            "label": cfg["label"],
            "exon_range": f"{s_ex}..{e_ex}",
            "codon_range": f"{cfg['start_codon']}..{cfg['end_codon']}",
            "covered_nonsense_count": len(covered),
            "excluded_nonsense_count": len(excluded),
            "total_nonsense_alleles": total,
            "covered_nonsense_pct": pct,
            "excluded_examples": excluded[:5],
        }

    return {
        "schema_version": 1,
        "tool": "super_exon_benchmark",
        "prediction_id": "PRED-0007",
        "reference_transcript": "RefSeq NM_000492.4 (CDS 71..4510, 27 exons, 1480 codons)",
        "epistemic_tier": "T2",
        "clinical_boundary": CLINICAL_BOUNDARY,
        "prior_art": [
            {
                "citation": "Barilla et al. (2026) Mol Ther Nucleic Acids, PMC13717110 (gold open access)",
                "catalog_source_id": "barilla-2026-omtn",
                "doi": "10.1016/j.omtn.2026.103123",
                "scope": "Partial CFTR cDNA exons 2-27 targeted insertion at CFTR intron 1 (TI-1)",
            },
            {
                "citation": "Bihler et al. (2026) J Cyst Fibros, PMID 42728175 (abstract/metadata only; publisher terms pending for bulk table ingestion)",
                "doi": "10.1016/j.jcf.2026.08.010",
                "scope": "In vitro FRT responsiveness profile across 417 CFTR variants to vanzacaftor/tezacaftor/ivacaftor",
            },
        ],
        "sentinel_variant_controls": controls_verified,
        "regime_coverage": regime_summaries,
    }


def run_falsification_suite() -> dict[str, Any]:
    """Run deterministic negative controls (--falsify) catching exon-boundary misattributions."""
    controls: list[dict[str, Any]] = []

    # Control 1: Exon 1 early PTC Q2X mislabeled as covered by Intron-1 TI-1 Exon 2-27 cassette must fail
    poisoned_q2x = (
        {"variant": "Q2X", "codon": 2, "expected_exon": 1, "expected_in_exons_2_27": True, "expected_in_exons_23_27": False},
    )
    try:
        verify_control_panel(poisoned_q2x)
        controls.append({"control": "exon1_q2x_ti1_exclusion_guard", "passed": False, "detail": "Q2X misattribution passed"})
    except SuperExonBenchmarkError as exc:
        controls.append({"control": "exon1_q2x_ti1_exclusion_guard", "passed": True, "detail": str(exc)})

    # Control 2: Exon 22 R1162X (c.3484C>T, ClinVar VCV000007137) mislabeled as covered by distal Exons 23-27 cassette must fail
    poisoned_r1162x = (
        {"variant": "R1162X", "codon": 1162, "expected_exon": 22, "expected_in_exons_2_27": True, "expected_in_exons_23_27": True},
    )
    try:
        verify_control_panel(poisoned_r1162x)
        controls.append({"control": "exon22_r1162x_distal_exclusion_guard", "passed": False, "detail": "R1162X misattribution passed"})
    except SuperExonBenchmarkError as exc:
        controls.append({"control": "exon22_r1162x_distal_exclusion_guard", "passed": True, "detail": str(exc)})

    # Control 3: Phase-0 Exon 22/23 splice boundary (c.3717 codon 1239 end vs c.3718 codon 1240 start)
    boundary_1239 = codon_to_exon(1239)
    boundary_1240 = codon_to_exon(1240)
    boundary_ok = boundary_1239 == 22 and boundary_1240 == 23 and not codon_covered_by_exon_range(1239, 23, 27)
    controls.append(
        {
            "control": "exon22_23_splice_boundary_1239_vs_1240",
            "passed": boundary_ok,
            "detail": f"codon_1239=exon_{boundary_1239}, codon_1240=exon_{boundary_1240}",
        }
    )

    all_passed = all(c["passed"] for c in controls)
    return {
        "tool": "super_exon_benchmark",
        "falsify_passed": all_passed,
        "negative_controls": controls,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", action="store_true", help="Emit full JSON benchmark report")
    parser.add_argument("--check", action="store_true", help="Run deterministic self-check")
    parser.add_argument("--falsify", action="store_true", help="Run negative-control falsification suite")
    args = parser.parse_args(argv)

    if args.falsify:
        f_rep = run_falsification_suite()
        if args.json:
            print(json.dumps(f_rep, indent=2))
        else:
            print(f"PASS falsify super_exon_benchmark ({len(f_rep['negative_controls'])}/3 negative controls)" if f_rep["falsify_passed"] else "FAIL falsify super_exon_benchmark")
        return 0 if f_rep["falsify_passed"] else 1

    rep = evaluate_super_exon_coverage()
    if args.check:
        f_rep = run_falsification_suite()
        ti1 = rep["regime_coverage"]["intron1_partial_cdna_exons_2_27"]
        dist = rep["regime_coverage"]["distal_super_exon_23_27"]
        full = rep["regime_coverage"]["full_coding_exons_1_27"]
        if (
            not f_rep["falsify_passed"]
            or ti1["covered_nonsense_count"] != 89
            or dist["covered_nonsense_count"] != 17
            or full["covered_nonsense_count"] != 92
        ):
            print("FAIL super_exon_benchmark self-check", file=sys.stderr)
            return 1
        print(
            f"PASS super_exon_benchmark: exon2_27={ti1['covered_nonsense_count']}/92 ({ti1['covered_nonsense_pct']}%) "
            f"exon23_27={dist['covered_nonsense_count']}/92 ({dist['covered_nonsense_pct']}%) "
            f"exon1_27={full['covered_nonsense_count']}/92 falsify_controls=3/3"
        )
        return 0

    if args.json:
        print(json.dumps(rep, indent=2))
    else:
        for rkey, summary in rep["regime_coverage"].items():
            print(
                f"{rkey}: {summary['covered_nonsense_count']}/{summary['total_nonsense_alleles']} "
                f"({summary['covered_nonsense_pct']}%, codons {summary['codon_range']})"
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
