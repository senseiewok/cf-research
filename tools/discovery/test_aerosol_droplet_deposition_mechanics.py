#!/usr/bin/env python3
"""Unit tests for aerosol droplet aerodynamic size distribution & mucociliary impaction mechanics (T-0408).

Tests physical boundary invariants, Weibel Model A morphological scaling, Stokes impaction,
gravitational settling, Brownian diffusion, CF airway narrowing shifts, and platform benchmarks.
"""

from __future__ import annotations

import math
import pytest

from aerosol_droplet_deposition_mechanics import (
    AirwayGeneration,
    DepositionProfile,
    DropletDistribution,
    GenerationDeposition,
    NLM_ATTRIBUTION,
    PRESET_PLATFORMS,
    PRIMARY_DOI_BIOENG,
    PRIMARY_DOI_HEYDER,
    PRIMARY_DOI_KIM,
    CLINICAL_TRIAL_4D710,
    CLINICAL_TRIAL_VX522,
    RESEARCH_BOUNDARY,
    WEIBEL_CONDUCTING_AIRWAYS,
    benchmark_delivery_platforms,
    build_weibel_tree,
    cunningham_slip_correction,
    droplet_brownian_diffusivity,
    droplet_relaxation_time,
    droplet_settling_velocity,
    evaluate_cf_obstruction_sweep,
    extrathoracic_deposition_efficiency,
    log_normal_size_bins,
    simulate_deposition,
    single_size_generation_mechanics,
)


def test_grounding_identifiers_and_disclaimers():
    """Verify primary DOIs, clinical trial IDs, and mandatory NLM disclaimer."""
    assert PRIMARY_DOI_KIM == "10.1002/adma.201201800"
    assert PRIMARY_DOI_BIOENG == "10.3389/fbioe.2026.1855907"
    assert PRIMARY_DOI_HEYDER == "10.1513/pats.200409-046ta"
    assert CLINICAL_TRIAL_VX522 == "NCT05668741"
    assert CLINICAL_TRIAL_4D710 == "NCT05248230"
    assert "ClinicalTrials.gov does not evaluate safety or efficacy" in NLM_ATTRIBUTION
    assert "research, not medical advice" in RESEARCH_BOUNDARY


def test_cunningham_slip_correction():
    """Verify Cunningham slip correction factor Cc behaviour across submicron and micron ranges."""
    cc_submicron = cunningham_slip_correction(0.5e-6)
    cc_micron = cunningham_slip_correction(3.0e-6)
    cc_large = cunningham_slip_correction(10.0e-6)

    assert cc_submicron > cc_micron > cc_large > 1.0
    assert 1.3 < cc_submicron < 1.5
    assert 1.03 < cc_micron < 1.10
    assert 1.01 < cc_large < 1.03
    assert cunningham_slip_correction(0.0) == 1.0


def test_droplet_relaxation_and_settling_velocity():
    """Verify aerodynamic particle relaxation time and Stokes terminal settling velocity."""
    d_small = 1.0e-6   # 1 um
    d_large = 4.0e-6   # 4 um

    tau_small = droplet_relaxation_time(d_small)
    tau_large = droplet_relaxation_time(d_large)

    assert tau_small > 0.0
    assert tau_large > tau_small * 10.0   # scales roughly as d^2

    v_sed_small = droplet_settling_velocity(d_small)
    v_sed_large = droplet_settling_velocity(d_large)

    assert math.isclose(v_sed_small, tau_small * 9.80665, rel_tol=1e-5)
    assert v_sed_large > v_sed_small


def test_droplet_brownian_diffusivity():
    """Verify Brownian diffusion coefficient increases for smaller particles."""
    d_small = 0.2e-6   # 200 nm
    d_large = 2.0e-6   # 2 um

    d_diff_small = droplet_brownian_diffusivity(d_small)
    d_diff_large = droplet_brownian_diffusivity(d_large)

    assert d_diff_small > d_diff_large > 0.0


def test_weibel_tree_morphology():
    """Verify Weibel Model A conducting airway morphology instantiation and CF narrowing."""
    tree_healthy = build_weibel_tree(cf_narrowing_factor=1.0)
    assert len(tree_healthy) == 17  # G0 to G16
    assert tree_healthy[0].generation == 0
    assert tree_healthy[0].num_airways == 1
    assert math.isclose(tree_healthy[0].diameter_m, 0.018, rel_tol=1e-3)
    assert tree_healthy[16].generation == 16
    assert tree_healthy[16].num_airways == 65536

    # Test narrowing factor
    tree_obstructed = build_weibel_tree(cf_narrowing_factor=0.60)
    assert math.isclose(tree_obstructed[0].diameter_m, 0.018 * 0.60, rel_tol=1e-3)
    assert math.isclose(tree_obstructed[0].cross_sectional_area_m2, tree_healthy[0].cross_sectional_area_m2 * 0.36, rel_tol=1e-3)


def test_extrathoracic_deposition():
    """Verify oropharyngeal impaction increases with particle size and flow rate."""
    flow = 30.0  # L/min
    dep_1um = extrathoracic_deposition_efficiency(1.0, flow)
    dep_3um = extrathoracic_deposition_efficiency(3.0, flow)
    dep_6um = extrathoracic_deposition_efficiency(6.0, flow)

    assert 0.0 <= dep_1um < dep_3um < dep_6um <= 1.0

    # Flow rate scaling
    dep_3um_low = extrathoracic_deposition_efficiency(3.0, 15.0)
    dep_3um_high = extrathoracic_deposition_efficiency(3.0, 60.0)
    assert dep_3um_low < dep_3um_high


def test_log_normal_distribution_normalization():
    """Verify discretized log-normal size distribution integrates to 1.0."""
    bins = log_normal_size_bins(mmad_um=3.0, gsd=1.8, num_bins=50)
    total_mass = sum(weight for _, weight in bins)
    assert math.isclose(total_mass, 1.0, rel_tol=1e-5)
    assert len(bins) == 50
    assert all(d > 0.0 and w > 0.0 for d, w in bins)


def test_simulate_deposition_mass_conservation():
    """Verify that deposition simulation conserves mass across all compartments."""
    dist = PRESET_PLATFORMS["vx522_mesh_lnp"]
    profile = simulate_deposition(dist, inspiratory_flow_l_min=30.0, cf_narrowing_factor=1.0)

    total_sum = (
        profile.extrathoracic_fraction
        + profile.tracheobronchial_total_fraction
        + profile.alveolar_fraction
        + profile.exhaled_fraction
    )
    assert math.isclose(total_sum, 1.0, abs_tol=0.01)

    assert profile.tracheobronchial_total_fraction == pytest.approx(
        profile.upper_conducting_fraction + profile.small_conducting_fraction + profile.terminal_bronchiolar_fraction,
        abs=0.001
    )
    assert 0.0 <= profile.bronchial_delivery_index <= 1.0
    assert len(profile.generation_breakdown) == 17


def test_cf_narrowing_impaction_shift():
    """Verify that CF airway narrowing shifts deposition proximally due to increased local velocity."""
    dist = PRESET_PLATFORMS["vx522_mesh_lnp"]
    profile_healthy = simulate_deposition(dist, inspiratory_flow_l_min=30.0, cf_narrowing_factor=1.0)
    profile_cf = simulate_deposition(dist, inspiratory_flow_l_min=30.0, cf_narrowing_factor=0.60)

    # In narrowed airways, velocity increases inversely with diameter squared,
    # increasing Stokes number and impaction in upper conducting airways
    assert profile_cf.upper_conducting_fraction > profile_healthy.upper_conducting_fraction
    # Penetration into terminal bronchioles is consequently reduced
    assert profile_cf.terminal_bronchiolar_fraction < profile_healthy.terminal_bronchiolar_fraction


def test_platform_benchmark():
    """Verify benchmarking across preset delivery platforms."""
    bench = benchmark_delivery_platforms(flow_l_min=30.0, cf_narrowing=0.75)
    assert "vx522_mesh_lnp" in bench
    assert "4d710_mesh_aav" in bench
    assert "standard_jet_nebulizer" in bench
    assert "optimized_monodisperse_mesh" in bench

    # Monodisperse mesh should achieve higher bronchial delivery index than polydisperse jet nebulizer
    assert bench["optimized_monodisperse_mesh"]["bdi"] > bench["standard_jet_nebulizer"]["bdi"]


def test_obstruction_sweep_monotonicity():
    """Verify that obstruction sweep accurately captures progressive airway narrowing."""
    dist = PRESET_PLATFORMS["4d710_mesh_aav"]
    sweep = evaluate_cf_obstruction_sweep(dist, flow_l_min=30.0, narrowing_steps=[1.0, 0.8, 0.6, 0.4])

    assert len(sweep) == 4
    for i in range(len(sweep) - 1):
        assert sweep[i]["luminal_occlusion_pct"] < sweep[i + 1]["luminal_occlusion_pct"]
        # Upper conducting impaction steadily increases as lumen narrows
        assert sweep[i]["upper_conducting_pct"] <= sweep[i + 1]["upper_conducting_pct"]
