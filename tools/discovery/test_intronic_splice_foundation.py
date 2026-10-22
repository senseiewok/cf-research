#!/usr/bin/env python3
"""Unit tests for Deep Intronic Cryptic Exon & Splice Foundation Model Benchmark (T-0269)."""
from __future__ import annotations

import pytest

from intronic_splice_foundation import (
    CANONICAL_INTRONIC_MUTATIONS,
    NLM_ATTRIBUTION,
    RESEARCH_BOUNDARY,
    calculate_nearest_neighbor_thermodynamics,
    design_aso_for_mutation,
    evaluate_splice_foundation_benchmark,
    simulate_dose_response,
)


def test_canonical_mutations_inventory() -> None:
    """Verify all 4 canonical Class V deep-intronic / consensus splice mutations are registered."""
    expected_alleles = {
        "3849+10kbC->T",
        "2789+5G->A",
        "1811+1.6kbA->G",
        "3272-26A->G",
    }
    assert set(CANONICAL_INTRONIC_MUTATIONS.keys()) == expected_alleles

    for allele_id, prof in CANONICAL_INTRONIC_MUTATIONS.items():
        assert prof.mutation_id == allele_id
        assert prof.hgvs_cdna.startswith("c.")
        assert len(prof.pdb_accessions) >= 3
        assert len(prof.primary_dois) >= 1
        assert prof.clinical_trial_nct.startswith("NCT")


def test_spliceai_and_evo2_scores() -> None:
    """Test that SpliceAI and Evo-2 metrics correctly identify cryptic donor/acceptor mechanics."""
    prof_3849 = CANONICAL_INTRONIC_MUTATIONS["3849+10kbC->T"]
    assert prof_3849.spliceai_donor_gain >= 0.85
    assert prof_3849.spliceai_acceptor_gain >= 0.60
    assert prof_3849.evo2_log_odds < -2.5

    prof_2789 = CANONICAL_INTRONIC_MUTATIONS["2789+5G->A"]
    assert prof_2789.spliceai_donor_loss >= 0.75
    assert prof_2789.spliceai_donor_gain < 0.10
    assert prof_2789.evo2_log_odds < -2.0

    prof_1811 = CANONICAL_INTRONIC_MUTATIONS["1811+1.6kbA->G"]
    assert prof_1811.spliceai_donor_gain >= 0.80

    prof_3272 = CANONICAL_INTRONIC_MUTATIONS["3272-26A->G"]
    assert prof_3272.spliceai_acceptor_gain >= 0.75


def test_nearest_neighbor_thermodynamics() -> None:
    """Test SantaLucia nearest-neighbor thermodynamic calculations and chemistry adjustments."""
    test_seq = "TCCTTACTTACTTACCTC"
    gc_moe, tm_moe, raw_moe, net_moe = calculate_nearest_neighbor_thermodynamics(test_seq, "2_MOE_PS")
    gc_pmo, tm_pmo, raw_pmo, net_pmo = calculate_nearest_neighbor_thermodynamics(test_seq, "PMO")

    assert gc_moe == gc_pmo
    assert 30.0 <= gc_moe <= 60.0
    # 2'-MOE should have higher Tm boost than PMO (+1.85 C/nt vs +1.20 C/nt)
    assert tm_moe > tm_pmo
    assert tm_moe > 55.0
    assert net_moe < -15.0  # Spontaneous hybridization (negative dG)


def test_aso_design_for_all_mutations() -> None:
    """Test automated steric-blocking ASO design across all alleles and chemistries."""
    for mut_id in CANONICAL_INTRONIC_MUTATIONS:
        for chem in ("2_MOE_PS", "PMO"):
            aso = design_aso_for_mutation(mut_id, chem)
            assert aso.target_mutation == mut_id
            assert aso.chemistry == chem
            assert 15 <= aso.length_nt <= 25
            assert aso.net_effective_affinity_kcal_mol < -10.0
            assert 5.0 <= aso.ec50_nanomolar <= 150.0
            assert 0.40 <= aso.max_rescue_fraction <= 0.85


def test_dose_response_monotonicity_and_curative_thresholds() -> None:
    """Test concentration-dependent transcript restoration and crossing of 10% and 15% thresholds."""
    concs = [0.0, 10.0, 30.0, 60.0, 120.0, 250.0]
    points = simulate_dose_response("3849+10kbC->T", "2_MOE_PS", concs)

    assert len(points) == len(concs)
    # Monotonicity check
    for i in range(len(points) - 1):
        assert points[i].fraction_aso_bound <= points[i + 1].fraction_aso_bound
        assert points[i].canonical_transcript_fraction_f_wt <= points[i + 1].canonical_transcript_fraction_f_wt
        assert points[i].chloride_transport_pct_wt <= points[i + 1].chloride_transport_pct_wt

    # Baseline is below 10%
    assert not points[0].crosses_10pct_threshold
    assert not points[0].crosses_15pct_threshold

    # High concentration reaches therapeutic and curative rescue
    top = points[-1]
    assert top.crosses_10pct_threshold
    assert top.crosses_15pct_threshold
    assert top.chloride_transport_pct_wt >= 15.0
    assert "CURATIVE" in top.curative_therapeutic_status


def test_full_benchmark_execution() -> None:
    """Verify full benchmark evaluation across all 4 intronic mutations."""
    bench = evaluate_splice_foundation_benchmark()
    assert bench["answer_tier"] == "T2"
    assert "mutations_evaluated" in bench
    assert len(bench["mutations_evaluated"]) == 4

    for mut_id, data in bench["mutations_evaluated"].items():
        assert "profile" in data
        assert "aso_designs" in data
        assert "curative_rescue_achieved" in data
        for chem in ("2_MOE_PS", "PMO"):
            rescue = data["curative_rescue_achieved"][chem]
            assert rescue["reaches_10pct_curative_threshold"] is True
            assert rescue["max_chloride_pct_wt_at_250nm"] >= 10.0


def test_nlm_attribution_and_boundaries() -> None:
    """Verify strict clinical trial metadata provenance and disclaimer boundaries."""
    assert "ClinicalTrials" in NLM_ATTRIBUTION
    assert "National Library of Medicine" in NLM_ATTRIBUTION
    assert "Research, not medical advice" in RESEARCH_BOUNDARY
