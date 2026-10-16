"""Unit tests for the Curative Discovery Council Engine and Biophysical Models."""
from __future__ import annotations

import pytest
from biophysical_models import (
    stokes_einstein_diffusion,
    cf_mucus_steric_retardation,
    electrostatic_mucin_adhesion_penalty,
    calculate_mucus_penetration_index,
    nonsense_readthrough_kinetics,
)
from discovery_council import (
    phase1_generative_synthesis,
    phase2_orthogonal_falsification,
    phase3_epistemic_verification,
    execute_discovery_council_round,
)


def test_stokes_einstein_diffusion():
    """Verify physical Stokes-Einstein free diffusion coefficient calculations."""
    d0_30nm = stokes_einstein_diffusion(30.0)  # Rh = 30 nm (diameter 60 nm)
    assert 1e-12 < d0_30nm < 2e-11, f"Expected D0 in 10^-12 to 10^-11 range, got {d0_30nm}"
    with pytest.raises(ValueError):
        stokes_einstein_diffusion(-10.0)


def test_cf_mucus_steric_retardation():
    """Particles under 70 nm must diffuse through 110 nm pores; >150 nm particles entrapped."""
    retard_65nm = cf_mucus_steric_retardation(65.0, mesh_pore_diameter_nm=110.0)
    retard_160nm = cf_mucus_steric_retardation(160.0, mesh_pore_diameter_nm=110.0)
    
    assert retard_65nm > 0.25, f"Expected high transmission for 65 nm particle, got {retard_65nm}"
    assert retard_160nm <= 0.01, f"Expected entrapment for 160 nm particle, got {retard_160nm}"


def test_electrostatic_mucin_adhesion_penalty():
    """Near-neutral zeta (-2 mV) transmits freely; cationic (+10 mV) trapped by polyanionic mucins."""
    neutral_pen = electrostatic_mucin_adhesion_penalty(-2.0)
    cationic_pen = electrostatic_mucin_adhesion_penalty(+10.0)
    
    assert neutral_pen > 0.85, f"Expected >85% transmission for near-neutral particle, got {neutral_pen}"
    assert cationic_pen < 0.05, f"Expected severe entrapment for cationic particle, got {cationic_pen}"


def test_mucus_penetration_index_optimal_bounds():
    """Candidate 65 nm, -1.2 mV formulation must achieve OPTIMAL mucus-penetrating verdict."""
    idx_optimal = calculate_mucus_penetration_index(65.0, -1.2)
    assert idx_optimal["mucus_penetrating_verdict"] == "OPTIMAL"
    assert idx_optimal["effective_diffusivity_ratio"] > 0.25

    idx_failed = calculate_mucus_penetration_index(130.0, +8.0)
    assert idx_failed["mucus_penetrating_verdict"] == "ENTRAPPED"


def test_nonsense_readthrough_kinetics_g542x():
    """G542X requires tri-modal rescue (NMD inhibition + suppressor tRNA + potentiator) to hit 10% WT current."""
    # Baseline alone
    base = nonsense_readthrough_kinetics("G542X", nmd_inhibition_pct=0.0, suppressor_trna_efficiency_pct=0.0)
    assert base["clinical_therapeutic_threshold_met"] is False
    assert base["restored_isc_pct_wt"] < 1.0

    # Tri-modal synergy
    synergy = nonsense_readthrough_kinetics(
        "G542X",
        nmd_inhibition_pct=60.0,
        suppressor_trna_efficiency_pct=9.0,
        potentiator_fold_increase=2.0
    )
    assert synergy["clinical_therapeutic_threshold_met"] is True
    assert synergy["restored_isc_pct_wt"] >= 10.0
    assert synergy["salvaged_mrna_pct_wt"] > 60.0


def test_nonsense_readthrough_kinetics_w1282x():
    """W1282X UGA A context provides higher readthrough permissiveness than G542X UGA G."""
    g542x = nonsense_readthrough_kinetics("G542X", nmd_inhibition_pct=40.0, suppressor_trna_efficiency_pct=5.0)
    w1282x = nonsense_readthrough_kinetics("W1282X", nmd_inhibition_pct=40.0, suppressor_trna_efficiency_pct=5.0)

    assert w1282x["effective_readthrough_pct"] > g542x["effective_readthrough_pct"]
    assert w1282x["salvaged_mrna_pct_wt"] > g542x["salvaged_mrna_pct_wt"]


def test_discovery_council_round_execution():
    """Verify end-to-end Discovery Council deliberation, falsification, and Tier T2 verification."""
    round_lnp = execute_discovery_council_round("aerosol_sort_lnp")
    assert round_lnp["council_status"] == "APPROVED"
    assert round_lnp["phase_2_falsification"]["falsification_status"] == "STRESS_TESTED_AND_RESOLVED"
    assert round_lnp["phase_3_verification"]["answer_tier"] == "T2"
    assert len(round_lnp["phase_3_verification"]["verified_claims"]) >= 2

    round_ptc = execute_discovery_council_round("nonsense_trimodal")
    assert round_ptc["council_status"] == "APPROVED"
    assert round_ptc["phase_2_falsification"]["falsification_status"] == "STRESS_TESTED_AND_RESOLVED"
    assert round_ptc["phase_3_verification"]["answer_tier"] == "T2"
