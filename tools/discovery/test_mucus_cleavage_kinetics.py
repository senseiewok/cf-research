#!/usr/bin/env python3
"""Unit tests for CF Airway Mucus Micro-Rheology & Polymer Cleavage Kinetics (T-0253)."""
from __future__ import annotations

import pytest

from mucus_cleavage_kinetics import (
    DEFAULT_DISULFIDE_CONC_MM,
    DEFAULT_EDNA_INITIAL_MW_BP,
    NLM_ATTRIBUTION,
    RESEARCH_BOUNDARY,
    calculate_ogston_diffusivity,
    evaluate_mucolytic_regimens,
    simulate_mucus_cleavage,
)


def test_ogston_diffusivity_boundary_and_scaling() -> None:
    """Verify Ogston hydrodynamic/steric obstruction formulas and boundary bounds."""
    # Zero pore radius gives zero diffusion
    assert calculate_ogston_diffusivity(13.0, 0.0) == 0.0
    # Small pore (10 nm) strongly hinders AAV (13 nm)
    d_tight = calculate_ogston_diffusivity(13.0, 10.0)
    assert d_tight < 0.10
    # Large pore (200 nm) allows substantial diffusion
    d_open = calculate_ogston_diffusivity(13.0, 200.0)
    assert d_open > 0.80
    assert d_open > d_tight
    # Smaller particle (AAV r=13nm) diffuses faster than larger particle (LNP r=40nm)
    d_lnp = calculate_ogston_diffusivity(40.0, 200.0)
    assert d_open > d_lnp


def test_baseline_untreated_hydrogel_state() -> None:
    """Verify untreated CF sputum exhibits high elastic stiffness and impaired clearance."""
    states = simulate_mucus_cleavage(
        rh_dnase_conc_ug_ml=0.0,
        nac_conc_mm=0.0,
        duration_minutes=30.0,
    )
    assert len(states) > 1
    init_state = states[0]
    final_state = states[-1]

    # Baseline properties
    assert init_state.g_prime_elastic_pa >= 7.0
    assert init_state.loss_tangent_tan_delta < 0.30
    assert init_state.rheological_regime == "IMPAIRED_ELASTIC_GEL"
    assert init_state.mucociliary_clearability_index < 0.10

    # Untreated controls remain virtually unchanged
    assert abs(final_state.g_prime_elastic_pa - init_state.g_prime_elastic_pa) < 0.1


def test_rh_dnase_kinetics_and_actin_resistance() -> None:
    """Verify rhDNase I endonucleolytic cleavage and enhanced rate of actin-resistant variant."""
    # Standard DNase (inhibited by G-actin)
    standard_states = simulate_mucus_cleavage(
        rh_dnase_conc_ug_ml=5.0,
        nac_conc_mm=0.0,
        actin_resistant_dnase=False,
        duration_minutes=30.0,
    )
    # Actin-resistant DNase (Ulmer et al. 1993)
    resistant_states = simulate_mucus_cleavage(
        rh_dnase_conc_ug_ml=5.0,
        nac_conc_mm=0.0,
        actin_resistant_dnase=True,
        duration_minutes=30.0,
    )

    std_final = standard_states[-1]
    res_final = resistant_states[-1]

    # Both degrade eDNA
    assert std_final.edna_mw_bp < DEFAULT_EDNA_INITIAL_MW_BP
    assert res_final.edna_mw_bp < DEFAULT_EDNA_INITIAL_MW_BP
    # Actin-resistant variant degrades eDNA faster / more completely
    assert res_final.edna_mw_bp <= std_final.edna_mw_bp
    assert res_final.g_prime_elastic_pa <= std_final.g_prime_elastic_pa


def test_nac_thiol_reduction_kinetics() -> None:
    """Verify NAC concentration-dependent reduction of intermolecular mucin disulfide crosslinks."""
    nac_states = simulate_mucus_cleavage(
        rh_dnase_conc_ug_ml=0.0,
        nac_conc_mm=20.0,
        duration_minutes=30.0,
    )
    init_state = nac_states[0]
    final_state = nac_states[-1]

    assert final_state.disulfide_crosslinks_mm < init_state.disulfide_crosslinks_mm
    assert final_state.fraction_disulfides_intact < 1.0
    # Mucin G' reduction
    assert final_state.g_prime_elastic_pa < init_state.g_prime_elastic_pa


def test_dual_rh_dnase_and_nac_synergy() -> None:
    """Verify combination therapy produces synergistic viscoelasticity reduction and vector penetration."""
    dnase_only = simulate_mucus_cleavage(rh_dnase_conc_ug_ml=5.0, nac_conc_mm=0.0)
    nac_only = simulate_mucus_cleavage(rh_dnase_conc_ug_ml=0.0, nac_conc_mm=20.0)
    combo = simulate_mucus_cleavage(rh_dnase_conc_ug_ml=5.0, nac_conc_mm=20.0)

    # Complex modulus reduction in combo should surpass monotherapies
    assert combo[-1].complex_modulus_g_star_pa < dnase_only[-1].complex_modulus_g_star_pa
    assert combo[-1].complex_modulus_g_star_pa < nac_only[-1].complex_modulus_g_star_pa

    # Dual therapy liquefies sputum and expands mesh pores
    assert combo[-1].mesh_pore_radius_nm > dnase_only[-1].mesh_pore_radius_nm
    assert combo[-1].aav_relative_diffusivity > dnase_only[-1].aav_relative_diffusivity
    assert combo[-1].lnp_relative_diffusivity > dnase_only[-1].lnp_relative_diffusivity
    assert combo[-1].loss_tangent_tan_delta > 0.50


def test_evaluate_mucolytic_regimens_benchmark() -> None:
    """Verify multi-regimen benchmark evaluation structure and clinical trial grounding."""
    bench = evaluate_mucolytic_regimens()
    assert bench["answer_tier"] == "T2"
    assert len(bench["regimen_summaries"]) == 5

    combo_summary = next(s for s in bench["regimen_summaries"] if "Synergy" in s["regimen_name"])
    assert combo_summary["percent_reduction_viscoelasticity"] >= 65.0
    assert combo_summary["aav_diffusivity_fold_increase"] >= 1.2
    assert combo_summary["reaches_physiologic_clearability"] is True

    # Check trial grounding
    for s in bench["regimen_summaries"]:
        assert "NCT06350474" in s["clinical_trial_grounding"]  # SIMPLIFY trial
        assert "NCT00996424" in s["clinical_trial_grounding"]  # Inhaled NAC trial
        assert "1ATN" in s["structural_pdb_evidence"]


def test_nlm_attribution_and_boundaries() -> None:
    """Verify strict clinical trial metadata provenance and disclaimer boundaries."""
    assert "ClinicalTrials" in NLM_ATTRIBUTION
    assert "National Library of Medicine" in NLM_ATTRIBUTION
    assert "Research, not medical advice" in RESEARCH_BOUNDARY
