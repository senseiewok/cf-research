#!/usr/bin/env python3
"""Unit tests for aerosol nebulizer mesh pore shear stress & LNP integrity dynamics (T-0413).

Validates physical boundary invariants, Hagen-Poiseuille shear scaling, capillary deformation,
critical bilayer strain, mRNA cargo leakage, colloidal aggregation, and platform benchmarks.
"""

from __future__ import annotations

import math
import pytest

from nebulizer_shear_stress_integrity import (
    CRITICAL_CAPILLARY_NUMBER,
    CRITICAL_LYSIS_STRAIN,
    DEFAULT_AREA_MODULUS_N_M,
    DEFAULT_BILAYER_TENSION_N_M,
    DEFAULT_BENDING_RIGIDITY_KBT,
    NLM_ATTRIBUTION,
    PRESET_MESH_PLATFORMS,
    PRIMARY_DOI_BIOENG,
    PRIMARY_DOI_HEYDER,
    PRIMARY_DOI_KIM,
    CLINICAL_TRIAL_VX522,
    RESEARCH_BOUNDARY,
    ColloidAggregationProfile,
    LNPBiophysicalParameters,
    LNPIntegrityResult,
    MeshApertureParameters,
    PoreShearProfile,
    benchmark_nebulizer_platforms,
    calculate_bilayer_areal_strain,
    calculate_capillary_number,
    calculate_hagen_poiseuille_shear_stress,
    calculate_taylor_deformation_parameter,
    evaluate_lnp_vesicle_deformation_and_leakage,
    evaluate_pore_radius_sweep,
    simulate_colloid_aggregation_and_mucus_clearance,
    simulate_pore_extrusion_shear_profile,
)


def test_grounding_identifiers_and_disclaimers():
    """Verify primary DOIs, clinical trial IDs, and mandatory NLM disclaimer."""
    assert PRIMARY_DOI_KIM == "10.1038/s41551-021-00724-5"
    assert PRIMARY_DOI_BIOENG == "10.3389/fbioe.2026.1855907"
    assert PRIMARY_DOI_HEYDER == "10.1513/pats.200409-046ta"
    assert CLINICAL_TRIAL_VX522 == "NCT05668741"
    assert "ClinicalTrials.gov does not evaluate safety or efficacy" in NLM_ATTRIBUTION
    assert "research, not medical advice" in RESEARCH_BOUNDARY


def test_hagen_poiseuille_shear_scaling():
    """Verify Hagen-Poiseuille formula tau_wall = 4 * mu * Q / (pi * r^3)."""
    r_pore = 2.0e-6             # 2 um
    q_flow = 1.0e-11            # 10 pL / s
    mu = 1.002e-3               # Pa * s

    v_mean, gamma_wall, tau_wall = calculate_hagen_poiseuille_shear_stress(r_pore, q_flow, mu)

    expected_area = math.pi * (r_pore ** 2)
    expected_v = q_flow / expected_area
    expected_gamma = (4.0 * q_flow) / (math.pi * (r_pore ** 3))
    expected_tau = mu * expected_gamma

    assert math.isclose(v_mean, expected_v, rel_tol=1e-5)
    assert math.isclose(gamma_wall, expected_gamma, rel_tol=1e-5)
    assert math.isclose(tau_wall, expected_tau, rel_tol=1e-5)

    # Invariance: halving pore radius at fixed Q increases shear stress 8-fold (r^3 scaling)
    _, _, tau_half = calculate_hagen_poiseuille_shear_stress(r_pore / 2.0, q_flow, mu)
    assert math.isclose(tau_half, tau_wall * 8.0, rel_tol=1e-4)


def test_mesh_extrusion_shear_profile_physics():
    """Verify flow simulation across micro-mesh nozzles maintains physical constraints."""
    mesh = PRESET_MESH_PLATFORMS["vx522_clinical_mesh"]
    profile = simulate_pore_extrusion_shear_profile(mesh, total_flow_rate_ml_min=0.40)

    # Reynolds number in micro-pores must be strictly laminar (Re << 100)
    assert 0.0 < profile.reynolds_number < 50.0
    # Mean and peak velocities must be positive and follow pulsatility
    assert profile.mean_pore_velocity_m_s > 0.0
    assert profile.peak_pore_velocity_m_s > profile.mean_pore_velocity_m_s
    assert math.isclose(
        profile.peak_pore_velocity_m_s,
        profile.mean_pore_velocity_m_s * mesh.pulsatility_factor,
        rel_tol=1e-4
    )
    # Peak shear stress must scale directly with peak flow
    assert profile.peak_wall_shear_stress_pa > profile.wall_shear_stress_pa > 0.0
    # Extrusion residence time through ~45 um plate should be on sub-millisecond timescale
    assert 0.01 < profile.pore_residence_time_ms < 5.0


def test_capillary_number_and_deformation():
    """Verify dimensionless Capillary number and Taylor deformation scaling."""
    r_lnp = 50.0e-9             # 50 nm radius (100 nm diameter)
    sigma = DEFAULT_BILAYER_TENSION_N_M  # 0.015 N / m

    tau_low = 100.0             # 100 Pa
    tau_high = 2000.0           # 2,000 Pa

    ca_low = calculate_capillary_number(tau_low, r_lnp, sigma)
    ca_high = calculate_capillary_number(tau_high, r_lnp, sigma)

    assert ca_low > 0.0
    assert math.isclose(ca_high, ca_low * 20.0, rel_tol=1e-5)

    # Taylor deformation parameter
    d_low = calculate_taylor_deformation_parameter(ca_low)
    d_high = calculate_taylor_deformation_parameter(ca_high)
    assert 0.0 < d_low < d_high < 1.0


def test_bilayer_areal_strain_and_lysis():
    """Verify membrane areal strain alpha and critical rupture threshold."""
    lnp = LNPBiophysicalParameters()
    tau_mild = 300.0
    tau_severe = 100000.0

    strain_mild = calculate_bilayer_areal_strain(tau_mild, lnp)
    strain_severe = calculate_bilayer_areal_strain(tau_severe, lnp)

    # Mild shear strain should be well below 4% critical lysis threshold
    assert strain_mild < CRITICAL_LYSIS_STRAIN
    # Severe shear strain exceeds critical lysis threshold
    assert strain_severe > CRITICAL_LYSIS_STRAIN


def test_evaluate_lnp_vesicle_leakage_conservation():
    """Verify that cargo leakage and intact payload conserve 100% total mass."""
    mesh = PRESET_MESH_PLATFORMS["vx522_clinical_mesh"]
    lnp = LNPBiophysicalParameters()

    result = evaluate_lnp_vesicle_deformation_and_leakage(mesh, lnp, flow_rate_ml_min=0.40)

    total_pct = result.fractional_mrna_leakage_pct + result.intact_payload_fraction_pct
    assert math.isclose(total_pct, 100.0, abs_tol=0.01)

    # In clinical vibrating mesh, intact payload should be high (>85%)
    assert result.intact_payload_fraction_pct >= 85.0
    assert result.integrity_status in ("HIGH_INTEGRITY", "MODERATE_SHEAR")
    assert 0.0 <= result.pore_opening_probability <= 1.0


def test_colloid_aggregation_and_mucus_penetration():
    """Verify secondary aggregation kinetics and CF mucus steric cutoff."""
    lnp = LNPBiophysicalParameters(initial_hydrodynamic_diameter_nm=90.0)
    mesh = PRESET_MESH_PLATFORMS["vx522_clinical_mesh"]
    shear = simulate_pore_extrusion_shear_profile(mesh, total_flow_rate_ml_min=0.40)

    colloid = simulate_colloid_aggregation_and_mucus_clearance(
        lnp, shear, residence_time_s=0.30, mucus_mesh_pore_nm=120.0
    )

    # Post-aerosol diameter must be >= initial diameter
    assert colloid.post_aerosol_diameter_nm >= colloid.initial_diameter_nm
    assert colloid.diameter_shift_nm >= 0.0
    assert 0.0 <= colloid.peg_loss_fraction_pct <= 100.0
    assert 0.0 <= colloid.mucus_mesh_penetration_fraction_pct <= 100.0

    # Penetration through 120 nm mucus pores should be higher for 90 nm LNP than for 250 nm aggregates
    colloid_large = simulate_colloid_aggregation_and_mucus_clearance(
        LNPBiophysicalParameters(initial_hydrodynamic_diameter_nm=250.0), shear, mucus_mesh_pore_nm=120.0
    )
    assert colloid.mucus_mesh_penetration_fraction_pct > colloid_large.mucus_mesh_penetration_fraction_pct
    assert colloid_large.mcc_evasion_status == "STERICALLY_BLOCKED"


def test_platform_benchmark_comparative_efficacy():
    """Verify nebulizer benchmark reflects established clinical and biophysical hierarchy."""
    bench = benchmark_nebulizer_platforms()

    assert "vibrating_mesh_vx522" in bench
    assert "optimized_soft_mesh" in bench
    assert "pneumatic_jet_nebulizer" in bench
    assert "ultrasonic_nebulizer" in bench

    mesh_res = bench["vibrating_mesh_vx522"]
    jet_res = bench["pneumatic_jet_nebulizer"]
    ultra_res = bench["ultrasonic_nebulizer"]
    soft_res = bench["optimized_soft_mesh"]

    # Intact mRNA hierarchy: Soft Mesh >= Vibrating Mesh > Pneumatic Jet > Ultrasonic
    assert soft_res.intact_mrna_pct >= mesh_res.intact_mrna_pct
    assert mesh_res.intact_mrna_pct > jet_res.intact_mrna_pct
    assert jet_res.intact_mrna_pct > ultra_res.intact_mrna_pct

    # Ultrasonic cavitation has highest temperature rise and severe aggregation
    assert ultra_res.temperature_rise_c > 10.0
    assert ultra_res.hydrodynamic_diameter_shift_nm > mesh_res.hydrodynamic_diameter_shift_nm

    # Overall delivery score ranking
    assert mesh_res.effective_cftr_delivery_score > jet_res.effective_cftr_delivery_score
    assert jet_res.effective_cftr_delivery_score > ultra_res.effective_cftr_delivery_score


def test_pore_radius_sweep_monotonicity():
    """Verify that expanding pore radius monotonically decreases peak shear stress."""
    sweep = evaluate_pore_radius_sweep(pore_radii_um=[1.5, 2.0, 2.5, 3.0, 3.5, 4.0])

    assert len(sweep) == 6
    for i in range(len(sweep) - 1):
        assert sweep[i]["pore_radius_um"] < sweep[i + 1]["pore_radius_um"]
        # Peak shear stress decreases monotonically with increasing pore aperture
        assert sweep[i]["peak_wall_shear_stress_pa"] > sweep[i + 1]["peak_wall_shear_stress_pa"]
        # Intact mRNA payload increases or stays high as shear drops
        assert sweep[i]["intact_mrna_pct"] <= sweep[i + 1]["intact_mrna_pct"]
