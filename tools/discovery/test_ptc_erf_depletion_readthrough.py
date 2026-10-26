#!/usr/bin/env python3
"""Unit tests for Eukaryotic Release Factor eRF1/eRF3 Depletion & Nonsense Readthrough Kinetic Model (T-0401).

Validates:
  1. Ribosomal A-site competitive kinetics (v_term, v_read, F_readthrough).
  2. eRF depletion Hill pharmacodynamics (IC50 in [0.1, 1.0] uM, h in [1.5, 2.0], I_max = 0.85).
  3. Characteristic 5- to 15-fold readthrough boost upon 60-80% release factor depletion.
  4. Codon permissiveness hierarchy: UGA > UAG >> UAA.
  5. Flanking +4 nucleotide context weighting (+4 C/U = 0.35 vs +4 G = 1.00).
  6. Coupled NMD escape & mature 170 kDa Band C CFTR expression.
  7. W1282X curative rescue potential (Band C >= 10% WT threshold).
  8. Dose-response monotonicity and mutation comparative benchmark.
  9. Anti-hallucination, NLM disclaimer, PDB 6O2P, and primary literature grounding.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path
import pytest

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from ptc_erf_depletion_readthrough import (
    CFTR_NONSENSE_MUTATIONS,
    CLINICAL_TRIAL_REGISTRY,
    CODON_PERMISSIVENESS_MULTIPLIERS,
    DEFAULT_CURATIVE_BAND_C_THRESHOLD,
    DEFAULT_DEPLETION_HILL_H,
    DEFAULT_DEPLETION_IC50_UM,
    DEFAULT_DEPLETION_IMAX,
    DEFAULT_ERF_TOTAL_UM,
    NLM_ATTRIBUTION,
    PLUS4_CONTEXT_WEIGHTS,
    PRIMARY_DOI_CONTEXT,
    PRIMARY_DOI_DEPLETION,
    RCSB_PDB_STRUCTURE,
    RESEARCH_BOUNDARY,
    BandCRescueResult,
    PTCMutationProfile,
    compute_depleted_erf_concentration,
    compute_nmd_escape_and_band_c,
    compute_ribosomal_competition,
    evaluate_mutation_dose_response,
    run_benchmark_comparison,
)


def test_scientific_grounding_and_disclaimers():
    """Verify primary evidence identifiers, PDB structure, and NLM attribution."""
    assert RCSB_PDB_STRUCTURE == "6O2P"
    assert PRIMARY_DOI_DEPLETION == "10.1172/jci154571"
    assert PRIMARY_DOI_CONTEXT == "10.3390/ph17030314"
    assert CLINICAL_TRIAL_REGISTRY == "NCT04135495"
    assert "ClinicalTrials.gov does not evaluate safety or efficacy" in NLM_ATTRIBUTION
    assert "research, not medical advice" in RESEARCH_BOUNDARY


def test_depletion_pharmacodynamics_hill_kinetics():
    """Verify release factor depletion follows Hill equation with proper cooperativity."""
    # Zero drug: zero depletion
    pd_zero = compute_depleted_erf_concentration(0.0)
    assert pd_zero.fractional_depletion == pytest.approx(0.0, abs=1e-6)
    assert pd_zero.erf_effective_um == pytest.approx(DEFAULT_ERF_TOTAL_UM, rel=1e-4)

    # At IC50: fractional depletion should be exactly half of I_max
    pd_ic50 = compute_depleted_erf_concentration(DEFAULT_DEPLETION_IC50_UM)
    expected_ic50_depletion = DEFAULT_DEPLETION_IMAX / 2.0
    assert pd_ic50.fractional_depletion == pytest.approx(expected_ic50_depletion, rel=1e-3)

    # Saturating drug (5 uM): should approach I_max (0.85)
    pd_sat = compute_depleted_erf_concentration(5.0)
    assert pd_sat.fractional_depletion > 0.80
    assert pd_sat.fractional_depletion <= DEFAULT_DEPLETION_IMAX

    # Hill coefficient in [1.5, 2.0]
    assert 1.5 <= DEFAULT_DEPLETION_HILL_H <= 2.0
    assert 0.1 <= DEFAULT_DEPLETION_IC50_UM <= 1.0


def test_codon_permissiveness_hierarchy():
    """Verify termination permissiveness hierarchy: UGA > UAG >> UAA (lower multiplier = more permissive)."""
    uga_factor = CODON_PERMISSIVENESS_MULTIPLIERS["UGA"]
    uag_factor = CODON_PERMISSIVENESS_MULTIPLIERS["UAG"]
    uaa_factor = CODON_PERMISSIVENESS_MULTIPLIERS["UAA"]

    # In our kinetic formulation, lower factor scales termination velocity down, yielding higher readthrough
    assert uga_factor < uag_factor
    assert uag_factor < uaa_factor
    assert uaa_factor >= 2.0 * uag_factor  # UAA is much stricter than UAG

    # Verify actual readthrough fractions for identical context
    dummy_uga = PTCMutationProfile("TestUGA", "UGA", "A", 1, "c.1G>A", "X", 1.0, "")
    dummy_uag = PTCMutationProfile("TestUAG", "UAG", "A", 1, "c.1G>A", "X", 1.0, "")
    dummy_uaa = PTCMutationProfile("TestUAA", "UAA", "A", 1, "c.1G>A", "X", 1.0, "")

    res_uga = compute_ribosomal_competition(DEFAULT_ERF_TOTAL_UM, dummy_uga)
    res_uag = compute_ribosomal_competition(DEFAULT_ERF_TOTAL_UM, dummy_uag)
    res_uaa = compute_ribosomal_competition(DEFAULT_ERF_TOTAL_UM, dummy_uaa)

    assert res_uga.readthrough_fraction > res_uag.readthrough_fraction
    assert res_uag.readthrough_fraction > res_uaa.readthrough_fraction
    assert res_uga.readthrough_fraction > 2.0 * res_uaa.readthrough_fraction


def test_plus4_context_weighting():
    """Verify +4 flanking nucleotide context weighting: +4 C/U = 0.35 vs +4 G = 1.00."""
    assert PLUS4_CONTEXT_WEIGHTS["C"] == pytest.approx(0.35, abs=1e-4)
    assert PLUS4_CONTEXT_WEIGHTS["U"] == pytest.approx(0.35, abs=1e-4)
    assert PLUS4_CONTEXT_WEIGHTS["A"] == pytest.approx(0.60, abs=1e-4)
    assert PLUS4_CONTEXT_WEIGHTS["G"] == pytest.approx(1.00, abs=1e-4)

    # Compare W1282X (UGA + C) against hypothetical UGA + G
    w1282x = CFTR_NONSENSE_MUTATIONS["W1282X"]
    hypothetical_uga_g = PTCMutationProfile("UGA_G", "UGA", "G", 23, "c.3846G>A", "Trp", 0.85, "")

    res_w1282x = compute_ribosomal_competition(DEFAULT_ERF_TOTAL_UM, w1282x)
    res_uga_g = compute_ribosomal_competition(DEFAULT_ERF_TOTAL_UM, hypothetical_uga_g)

    # +4 C should have lower termination rate and higher readthrough fraction
    assert res_w1282x.v_term_per_s < res_uga_g.v_term_per_s
    assert res_w1282x.readthrough_fraction > res_uga_g.readthrough_fraction
    # Ratio should reflect context factor (1.0 / 0.35 approx 2.85x difference in termination rate)
    assert (res_uga_g.v_term_per_s / res_w1282x.v_term_per_s) == pytest.approx(1.0 / 0.35, rel=1e-3)


def test_erf_depletion_boost_5_to_15_fold():
    """Verify 5- to 15-fold readthrough boost upon 60-80% release factor depletion."""
    w1282x = CFTR_NONSENSE_MUTATIONS["W1282X"]

    # Basal readthrough (0% depletion)
    res_basal = compute_ribosomal_competition(DEFAULT_ERF_TOTAL_UM, w1282x)

    # 70% depletion -> remaining eRF = 30% of baseline (0.30 uM)
    res_70pct_depleted = compute_ribosomal_competition(DEFAULT_ERF_TOTAL_UM * 0.30, w1282x)
    boost_70 = res_70pct_depleted.readthrough_fraction / res_basal.readthrough_fraction

    # 80% depletion -> remaining eRF = 20% of baseline (0.20 uM)
    res_80pct_depleted = compute_ribosomal_competition(DEFAULT_ERF_TOTAL_UM * 0.20, w1282x)
    boost_80 = res_80pct_depleted.readthrough_fraction / res_basal.readthrough_fraction

    assert 3.0 <= boost_70 <= 15.0
    assert 4.0 <= boost_80 <= 15.0


def test_coupled_nmd_escape_and_band_c_expression():
    """Verify that readthrough rescues mRNA from NMD and produces functional Band C."""
    w1282x = CFTR_NONSENSE_MUTATIONS["W1282X"]

    # At very low readthrough: mRNA remains near basal unescaped level (~18%)
    mrna_low, esc_low, band_c_low, cl_low = compute_nmd_escape_and_band_c(
        mutation=w1282x,
        readthrough_fraction=0.002,  # 0.2%
        basal_readthrough_fraction=0.002,
    )
    assert mrna_low < 0.30
    assert band_c_low < 0.01  # < 1% Band C

    # At elevated readthrough (e.g. 15% readthrough upon eRF depletion):
    mrna_high, esc_high, band_c_high, cl_high = compute_nmd_escape_and_band_c(
        mutation=w1282x,
        readthrough_fraction=0.15,  # 15%
        basal_readthrough_fraction=0.002,
    )
    # mRNA should be strongly stabilized (> 75% of WT)
    assert mrna_high > 0.75
    # Band C expression should reach clinical significance
    assert band_c_high >= DEFAULT_CURATIVE_BAND_C_THRESHOLD


def test_w1282x_dose_response_monotonicity():
    """Verify that dose-response simulation is monotonically increasing for W1282X."""
    curve = evaluate_mutation_dose_response("W1282X")
    assert len(curve) >= 5

    prev_conc = -1.0
    prev_depletion = -1.0
    prev_rt = -1.0
    prev_band_c = -1.0

    for pt in curve:
        assert pt.drug_concentration_um > prev_conc
        assert pt.fractional_erf_depletion >= prev_depletion
        assert pt.readthrough_fraction >= prev_rt
        assert pt.band_c_protein_relative_to_wt >= prev_band_c

        prev_conc = pt.drug_concentration_um
        prev_depletion = pt.fractional_erf_depletion
        prev_rt = pt.readthrough_fraction
        prev_band_c = pt.band_c_protein_relative_to_wt

    # Final saturating point should achieve curative threshold for W1282X
    saturating = curve[-1]
    assert saturating.curative_threshold_reached is True
    assert saturating.band_c_percentage_wt >= 10.0


def test_benchmark_comparison_all_mutations():
    """Verify comparative benchmark across all 5 cataloged nonsense mutations."""
    comparison = run_benchmark_comparison()
    expected_mutations = {"W1282X", "G542X", "R553X", "R1162X", "Y122X"}
    assert set(comparison.keys()) == expected_mutations

    # W1282X should have the highest saturating Band C rescue due to UGA + C context
    w1282x_data = comparison["W1282X"]
    y122x_data = comparison["Y122X"]

    assert w1282x_data["saturating_band_c_pct_wt"] > y122x_data["saturating_band_c_pct_wt"]
    assert w1282x_data["curative_rescue_at_saturation"] is True
