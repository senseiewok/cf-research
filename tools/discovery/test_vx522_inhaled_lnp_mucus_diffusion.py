#!/usr/bin/env python3
"""Unit test suite for Inhaled mRNA LNP Mucus Mesh Diffusion & Surface Charge Neutralization Model (T-0400)."""

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

from vx522_inhaled_lnp_mucus_diffusion import (
    CANONICAL_FORMULATIONS,
    CLINICAL_TRIAL_REGISTRY,
    NLM_ATTRIBUTION,
    PRIMARY_DOI_KIM,
    PRIMARY_DOI_SUK,
    RESEARCH_BOUNDARY,
    calculate_mobile_fraction,
    dlvo_electrostatic_trap_energy,
    endosomal_escape_probability,
    evaluate_lnp_diffusion_profile,
    ionizable_lipid_protonation,
    ogston_sieving_factor,
    peg_shedding_kinetics,
    renkin_hindrance_polynomial,
    run_benchmark_comparison,
    stokes_einstein_diffusivity,
)


def test_stokes_einstein_diffusivity():
    """Verify Stokes-Einstein hydrodynamic diffusivity calculation across particle sizes."""
    # 80 nm LNP at 37 deg C in water (exact value approx 8.215 um^2/s)
    d_0_80nm = stokes_einstein_diffusivity(80.0)
    d_0_80nm_um2_s = d_0_80nm * 1e12
    assert pytest.approx(d_0_80nm_um2_s, rel=0.01) == 8.215

    # Diffusivity must scale inversely with particle diameter
    d_0_160nm = stokes_einstein_diffusivity(160.0)
    ratio = d_0_80nm / d_0_160nm
    assert pytest.approx(ratio, rel=1e-3) == 2.0


def test_ogston_sieving_and_renkin_hindrance():
    """Verify Ogston steric sieving and Renkin hydrodynamic hindrance scaling in mucus mesh."""
    pore_nm = 200.0

    # Small vs large particles
    s_small = ogston_sieving_factor(40.0, pore_diameter_nm=pore_nm)
    s_med = ogston_sieving_factor(80.0, pore_diameter_nm=pore_nm)
    s_large = ogston_sieving_factor(220.0, pore_diameter_nm=pore_nm)

    assert 0.0 < s_large < s_med < s_small <= 1.0

    f_small = renkin_hindrance_polynomial(40.0, pore_diameter_nm=pore_nm)
    f_med = renkin_hindrance_polynomial(80.0, pore_diameter_nm=pore_nm)
    f_large = renkin_hindrance_polynomial(220.0, pore_diameter_nm=pore_nm)

    assert 0.0 < f_med < f_small <= 1.0
    # Over-pore particle (220 nm > 200 nm pore) must have zero Renkin penetration
    assert f_large == 0.0


def test_dlvo_surface_charge_neutralization():
    """Verify DLVO electrostatic trap energy and near-neutral charge neutralization."""
    # Near-neutral LNP (zeta = -0.5 mV) against anionic mucus (psi = -15 mV)
    u_neutral = dlvo_electrostatic_trap_energy(80.0, -0.5)
    assert abs(u_neutral) < 1.0  # Must be strictly under thermal trap threshold
    f_mobile_neutral = calculate_mobile_fraction(u_neutral)
    assert f_mobile_neutral > 0.60

    # Cationic unshielded LNP (zeta = +18.5 mV)
    # Opposite charges produce strong negative attractive energy: U_elec << -1.0 k_B*T
    u_cationic = dlvo_electrostatic_trap_energy(80.0, 18.5)
    assert u_cationic < -5.0
    f_mobile_cationic = calculate_mobile_fraction(u_cationic)
    assert f_mobile_cationic < 0.01  # Strongly immobilized in mucus mesh


def test_anionic_repulsive_charge_mobility():
    """Verify strongly anionic LNP experiences electrostatic repulsion without adhesive trapping."""
    u_anionic = dlvo_electrostatic_trap_energy(80.0, -16.0)
    assert u_anionic > 0.0  # Repulsive interaction (like charges)
    f_mobile_anionic = calculate_mobile_fraction(u_anionic)
    assert f_mobile_anionic == 1.0  # Zero adhesive entrapment


def test_peg_shedding_kinetics():
    """Verify PEG-lipid shedding dynamics and uptake competence acquisition."""
    k_shed = 0.12  # h^-1 (typical DMG-PEG2000)

    # At t = 0, fully shielded
    bound_0, comp_0 = peg_shedding_kinetics(k_shed, 0.0)
    assert bound_0 == 1.0
    assert comp_0 == 0.0

    # At t = 4 hours
    bound_4, comp_4 = peg_shedding_kinetics(k_shed, 4.0)
    assert bound_4 < 0.65
    assert comp_4 > 0.35  # Competent for cellular uptake

    # Non-shedding lipid (DSPE-PEG2000, k_shed = 0.005 h^-1)
    k_nonshed = 0.005
    bound_ns_4, comp_ns_4 = peg_shedding_kinetics(k_nonshed, 4.0)
    assert bound_ns_4 > 0.95
    assert comp_ns_4 < 0.05


def test_ionizable_lipid_protonation_and_endosomal_escape():
    """Verify pH-dependent ionizable lipid protonation switch and endosomal membrane disruption."""
    pka = 6.38

    # Extracellular airway surface liquid pH (7.2)
    alpha_extra = ionizable_lipid_protonation(7.2, pka)
    assert alpha_extra < 0.15  # Uncharged in transit to avoid toxicity

    # Late endosomal pH (5.5)
    alpha_endo = ionizable_lipid_protonation(5.5, pka)
    assert alpha_endo > 0.85  # Heavily protonated to disrupt endosomal membrane

    # Endosomal escape probability
    p_esc_extra = endosomal_escape_probability(alpha_extra)
    p_esc_endo = endosomal_escape_probability(alpha_endo)
    assert p_esc_extra < 0.02
    assert p_esc_endo >= 0.20


def test_vx522_benchmark_curative_competence():
    """Verify that the clinical VX-522 aerosol benchmark achieves curative delivery competence."""
    form = CANONICAL_FORMULATIONS["VX-522_Aerosol_Benchmark"]
    result = evaluate_lnp_diffusion_profile(form)

    assert result.mucus_penetrating_particle_status is True
    assert result.curative_delivery_competent is True
    assert result.mucus_penetration_efficiency_pct >= 10.0
    assert result.uptake_competence_4h_pct >= 25.0
    assert result.endosomal_escape_probability_pct >= 10.0
    assert result.d_eff_over_d0 >= 0.02


def test_failed_formulations_sub_threshold():
    """Verify that cationic and oversized LNPs fail curative delivery competence."""
    cationic = CANONICAL_FORMULATIONS["Cationic_Unshielded_LNP"]
    res_cat = evaluate_lnp_diffusion_profile(cationic)
    assert res_cat.mucus_penetrating_particle_status is False
    assert res_cat.curative_delivery_competent is False
    assert res_cat.mobile_fraction < 0.05

    oversized = CANONICAL_FORMULATIONS["Oversized_Aerosol_LNP"]
    res_over = evaluate_lnp_diffusion_profile(oversized)
    assert res_over.renkin_hindrance_factor == 0.0
    assert res_over.curative_delivery_competent is False


def test_nlm_attribution_and_epistemic_boundary():
    """Verify NLM attribution string, registry grounding, and clinical boundary compliance."""
    assert "ClinicalTrials.gov does not evaluate safety or efficacy" in NLM_ATTRIBUTION
    assert "Preclinical in-silico" in RESEARCH_BOUNDARY
    assert CLINICAL_TRIAL_REGISTRY == "NCT06237335"
    assert PRIMARY_DOI_SUK == "10.1073/pnas.1901667116"
    assert PRIMARY_DOI_KIM == "10.1038/s41551-021-00724-5"


def test_run_benchmark_comparison():
    """Verify multi-formulation comparative benchmark output structure."""
    bench = run_benchmark_comparison()
    assert len(bench) >= 4
    assert "VX-522_Aerosol_Benchmark" in bench
    assert "Cationic_Unshielded_LNP" in bench
    assert "Oversized_Aerosol_LNP" in bench
    for name, data in bench.items():
        assert "d_eff_over_d0" in data
        assert "curative_delivery_competent" in data
        assert "endosomal_escape_probability_pct" in data
