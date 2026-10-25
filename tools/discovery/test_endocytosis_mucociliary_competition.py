#!/usr/bin/env python3
"""Unit tests for tools/discovery/endocytosis_mucociliary_competition.py."""
import math
import pytest

from endocytosis_mucociliary_competition import (
    AirwayTransportParameters,
    CompetitionOutcome,
    calculate_effective_endocytosis_rate,
    calculate_pcl_hindered_diffusivity,
    calculate_stokes_einstein_diffusion,
    compare_healthy_vs_cf_phenotype,
    simulate_endocytosis_vs_clearance,
)


def test_stokes_einstein_diffusivity_scaling():
    """Verify Stokes-Einstein inverse diameter scaling and temperature dependence."""
    d_small = calculate_stokes_einstein_diffusion(diameter_nm=30.0)
    d_large = calculate_stokes_einstein_diffusion(diameter_nm=60.0)
    assert d_small > 0
    assert d_large > 0
    assert math.isclose(d_small / d_large, 2.0, rel_tol=1e-3)

    # Temperature scaling: higher temp -> higher diffusivity
    d_cold = calculate_stokes_einstein_diffusion(diameter_nm=50.0, temperature_k=293.15)
    d_warm = calculate_stokes_einstein_diffusion(diameter_nm=50.0, temperature_k=310.15)
    assert d_warm > d_cold


def test_pcl_hindered_diffusivity_steric_partitioning():
    """Verify that steric mesh hindrance retards large particles relative to small."""
    d_eff_small = calculate_pcl_hindered_diffusivity(
        diameter_nm=25.0, mesh_size_nm=40.0
    )
    d_eff_large = calculate_pcl_hindered_diffusivity(
        diameter_nm=80.0, mesh_size_nm=40.0
    )
    assert d_eff_small > d_eff_large
    # 80 nm particle should be dramatically hindered in 40 nm mesh
    assert d_eff_small / d_eff_large > 10.0


def test_peg_dilemma_attenuation():
    """Verify PEG dilemma: dense surface PEG suppresses apical cell endocytosis."""
    k_unshielded = calculate_effective_endocytosis_rate(
        base_rate_min_inv=0.04, diameter_nm=50.0, peg_density_chains_nm2=0.01
    )
    k_shielded = calculate_effective_endocytosis_rate(
        base_rate_min_inv=0.04, diameter_nm=50.0, peg_density_chains_nm2=0.20
    )
    assert k_unshielded > k_shielded
    assert k_shielded > 0


def test_conservation_of_mass():
    """Verify that all fraction outcomes sum to 1.0 (internalized + cleared + trapped)."""
    outcome = simulate_endocytosis_vs_clearance()
    total_fraction = (
        outcome.fraction_internalized
        + outcome.fraction_cleared
        + outcome.fraction_trapped_mucus
    )
    assert math.isclose(total_fraction, 1.0, abs_tol=1e-6)
    assert 0.0 <= outcome.fraction_internalized <= 1.0
    assert 0.0 <= outcome.fraction_cleared <= 1.0
    assert 0.0 <= outcome.fraction_trapped_mucus <= 1.0


def test_damkohler_monotonicity_with_mcc_velocity():
    """Verify that faster mucociliary sweep increases Damköhler number and reduces uptake."""
    slow_p = AirwayTransportParameters(mucociliary_velocity_mm_per_min=0.5)
    fast_p = AirwayTransportParameters(mucociliary_velocity_mm_per_min=5.0)

    res_slow = simulate_endocytosis_vs_clearance(slow_p)
    res_fast = simulate_endocytosis_vs_clearance(fast_p)

    assert res_fast.damkohler_number > res_slow.damkohler_number
    assert res_slow.fraction_internalized > res_fast.fraction_internalized
    assert res_fast.fraction_cleared > res_slow.fraction_cleared


def test_healthy_vs_cf_phenotype_comparison():
    """Verify CF pathology exhibits impaired diffusion and lower cellular internalization."""
    comp = compare_healthy_vs_cf_phenotype(diameter_nm=65.0, peg_density=0.05)
    healthy = comp["healthy"]
    cf = comp["cystic_fibrosis"]

    assert healthy.diffusion_coefficient_pcl_m2_s > cf.diffusion_coefficient_pcl_m2_s
    assert cf.fraction_trapped_mucus > healthy.fraction_trapped_mucus
    assert healthy.fraction_internalized > cf.fraction_internalized


def test_asymptotic_diameter_behavior():
    """Verify delivery dynamics across small AAV (~25 nm) vs large LNP (~100 nm)."""
    aav_p = AirwayTransportParameters(nanoparticle_diameter_nm=25.0)
    lnp_p = AirwayTransportParameters(nanoparticle_diameter_nm=100.0)

    res_aav = simulate_endocytosis_vs_clearance(aav_p)
    res_lnp = simulate_endocytosis_vs_clearance(lnp_p)

    assert res_aav.diffusion_coefficient_pcl_m2_s > res_lnp.diffusion_coefficient_pcl_m2_s
    assert res_aav.pcl_transit_time_seconds < res_lnp.pcl_transit_time_seconds


def test_invalid_inputs_raise_errors():
    """Verify robust parameter validation and physical boundary assertions."""
    with pytest.raises(ValueError, match="Nanoparticle diameter must be positive"):
        calculate_stokes_einstein_diffusion(diameter_nm=-10.0)

    with pytest.raises(ValueError, match="Temperature must be positive"):
        calculate_stokes_einstein_diffusion(diameter_nm=50.0, temperature_k=-5.0)

    with pytest.raises(ValueError, match="PCL mesh size must be positive"):
        calculate_pcl_hindered_diffusivity(diameter_nm=50.0, mesh_size_nm=-20.0)

    with pytest.raises(ValueError, match="Base endocytosis rate cannot be negative"):
        calculate_effective_endocytosis_rate(base_rate_min_inv=-0.05, diameter_nm=50.0, peg_density_chains_nm2=0.05)

    with pytest.raises(ValueError, match="PCL height must be positive"):
        bad_p = AirwayTransportParameters(pcl_height_um=-2.0)
        simulate_endocytosis_vs_clearance(bad_p)
