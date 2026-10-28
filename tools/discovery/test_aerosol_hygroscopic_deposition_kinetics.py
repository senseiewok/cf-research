#!/usr/bin/env python3
"""Unit tests for aerosol droplet aerodynamic size distribution, hygroscopic kinetics, and airway bifurcation impaction (T-0473).

Tests Weibel Model A morphometry, Stokes-Cunningham particle aerodynamics,
Köhler dynamic hygroscopic equilibration, multi-mechanism deposition,
CF pathological airway narrowing shifts, delivery platform benchmarks,
and AH-10 negative controls.
"""

from __future__ import annotations

import math
import pytest

from aerosol_hygroscopic_deposition_kinetics import (
    AirwayGeneration,
    DepositionProfile,
    DropletDistribution,
    GenerationDeposition,
    NLM_ATTRIBUTION,
    PDB_STRUCTURE_CFTR,
    PRESET_PLATFORMS,
    PRIMARY_DOI_BUTTON,
    PRIMARY_DOI_DARQUENNE,
    PRIMARY_DOI_MIDDLETON,
    CLINICAL_TRIAL_VANZACAFTOR,
    CLINICAL_TRIAL_VX522,
    RESEARCH_BOUNDARY,
    WEIBEL_CONDUCTING_AIRWAYS,
    benchmark_delivery_platforms,
    build_weibel_tree,
    calculate_equilibrium_growth_ratio,
    calculate_extrathoracic_loss,
    calculate_generation_deposition_efficiencies,
    compare_healthy_vs_cf_obstruction,
    cunningham_slip_correction,
    droplet_brownian_diffusivity,
    droplet_relaxation_time,
    droplet_settling_velocity,
    evaluate_falsification_negative_control,
    generate_log_normal_bins,
    simulate_deposition,
    simulate_droplet_hygroscopic_evolution,
)


def test_grounding_identifiers_and_disclaimers():
    """Verify primary DOIs, clinical trial IDs, PDB structure, and mandatory disclaimers."""
    assert PRIMARY_DOI_BUTTON == "10.1152/physrev.00004.2021"
    assert PRIMARY_DOI_DARQUENNE == "10.1016/j.jaerosci.2018.06.002"
    assert PRIMARY_DOI_MIDDLETON == "10.1056/nejmoa1908639"
    assert CLINICAL_TRIAL_VX522 == "NCT05668741"
    assert CLINICAL_TRIAL_VANZACAFTOR == "NCT05033080"
    assert PDB_STRUCTURE_CFTR == "6MSM"
    assert "ClinicalTrials.gov does not evaluate safety or efficacy" in NLM_ATTRIBUTION
    assert "research, not medical advice" in RESEARCH_BOUNDARY


def test_cunningham_slip_correction():
    """Verify Cunningham slip correction factor Cc behaviour across submicron and micron ranges."""
    cc_submicron = cunningham_slip_correction(0.3e-6)
    cc_fine = cunningham_slip_correction(1.0e-6)
    cc_coarse = cunningham_slip_correction(5.0e-6)

    assert cc_submicron > cc_fine > cc_coarse > 1.0
    assert 1.4 < cc_submicron < 1.7
    assert 1.10 < cc_fine < 1.25
    assert 1.02 < cc_coarse < 1.06
    assert cunningham_slip_correction(0.0) == 1.0


def test_droplet_relaxation_and_settling_velocity():
    """Verify particle aerodynamic relaxation time and Stokes terminal settling velocity."""
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
    d_submicron = 0.2e-6   # 200 nm
    d_micron = 2.0e-6      # 2 um

    d_diff_sub = droplet_brownian_diffusivity(d_submicron)
    d_diff_mic = droplet_brownian_diffusivity(d_micron)

    assert d_diff_sub > d_diff_mic > 0.0
    assert droplet_brownian_diffusivity(0.0) == 0.0


def test_kohler_hygroscopic_growth():
    """Verify equilibrium growth factor G_eq and time-dependent hygroscopic evolution."""
    g_low = calculate_equilibrium_growth_ratio(0.50, kappa=0.35)
    g_high = calculate_equilibrium_growth_ratio(0.995, kappa=0.35)
    assert g_high > g_low > 1.0

    # Droplet of 3.5 um entering 99.5% RH grows towards equilibrium
    d_init = 3.5e-6
    d_evolved_short = simulate_droplet_hygroscopic_evolution(d_init, transit_time_s=0.05, airway_rh=0.995)
    d_evolved_long = simulate_droplet_hygroscopic_evolution(d_init, transit_time_s=1.5, airway_rh=0.995)

    assert d_evolved_long > d_evolved_short >= d_init


def test_weibel_tree_morphology():
    """Verify Weibel Model A conducting airway morphology instantiation and CF narrowing."""
    tree_healthy = build_weibel_tree(inspiratory_flow_l_min=30.0, cf_occlusion_factor=0.0)
    assert len(tree_healthy) == 17  # G0 to G16
    assert tree_healthy[0].generation == 0
    assert tree_healthy[0].num_airways == 1
    assert math.isclose(tree_healthy[0].diameter_m, 0.018, rel_tol=1e-3)
    assert tree_healthy[16].generation == 16
    assert tree_healthy[16].num_airways == 65536

    # Verify flow velocity decelerates in peripheral generations due to total cross-sectional area expansion
    v_trachea = tree_healthy[0].air_velocity_m_s
    v_terminal = tree_healthy[16].air_velocity_m_s
    assert v_trachea > v_terminal > 0.0

    # Verify CF occlusion narrows lumen and elevates velocity in conducting bronchioles
    tree_cf = build_weibel_tree(inspiratory_flow_l_min=30.0, cf_occlusion_factor=0.35)
    assert tree_cf[8].diameter_m < tree_healthy[8].diameter_m
    assert tree_cf[8].air_velocity_m_s > tree_healthy[8].air_velocity_m_s


def test_generation_deposition_mechanisms():
    """Verify individual impaction, sedimentation, and diffusion efficiencies are bounded in [0, 1]."""
    tree = build_weibel_tree(inspiratory_flow_l_min=30.0)
    g3 = tree[3]  # Segmental bronchi

    eta_imp, eta_sed, eta_diff, eta_tot, stk = calculate_generation_deposition_efficiencies(
        droplet_diameter_m=3.5e-6,
        airway=g3,
    )

    assert 0.0 <= eta_imp <= 1.0
    assert 0.0 <= eta_sed <= 1.0
    assert 0.0 <= eta_diff <= 1.0
    assert 0.0 <= eta_tot <= 1.0
    assert stk > 0.0
    # Combined total efficiency must be >= each individual mechanism
    assert eta_tot >= eta_imp
    assert eta_tot >= eta_sed


def test_extrathoracic_loss():
    """Verify extrathoracic throat loss increases with MMAD and flow rate."""
    loss_small = calculate_extrathoracic_loss(mmad_um=2.5, inspiratory_flow_l_min=30.0)
    loss_large = calculate_extrathoracic_loss(mmad_um=6.0, inspiratory_flow_l_min=30.0)
    loss_high_flow = calculate_extrathoracic_loss(mmad_um=2.5, inspiratory_flow_l_min=60.0)

    assert 0.0 <= loss_small < loss_large <= 1.0
    assert loss_high_flow > loss_small


def test_log_normal_distribution_generation():
    """Verify log-normal size bin generation mass balance."""
    bins = generate_log_normal_bins(mmad_um=3.5, gsd=1.6, num_bins=30)
    assert len(bins) == 30
    total_mass = sum(w for _, w in bins)
    assert math.isclose(total_mass, 1.0, rel_tol=1e-4)


def test_full_deposition_simulation_mass_balance():
    """Verify comprehensive deposition simulation satisfies conservation of aerosol mass."""
    platform = PRESET_PLATFORMS["vx522_vibrating_mesh"]
    prof = simulate_deposition(platform, inspiratory_flow_l_min=30.0, cf_occlusion_factor=0.0)

    total_accounted = (
        prof.extrathoracic_fraction
        + prof.tracheobronchial_total_fraction
        + prof.alveolar_fraction
        + prof.exhaled_fraction
    )
    assert math.isclose(total_accounted, 1.0, abs_tol=1e-3)

    assert prof.small_conducting_fraction > 0.0
    assert prof.terminal_bronchiolar_fraction > 0.0
    assert 0.0 < prof.bronchial_delivery_index <= 1.0


def test_cf_pathological_airway_occlusion_shift():
    """Verify CF luminal occlusion increases upper airway impaction and decreases BDI."""
    platform = PRESET_PLATFORMS["vx522_vibrating_mesh"]
    healthy = simulate_deposition(platform, inspiratory_flow_l_min=30.0, cf_occlusion_factor=0.0)
    cf_severe = simulate_deposition(platform, inspiratory_flow_l_min=30.0, cf_occlusion_factor=0.40)

    # Upper conducting airway impaction is heightened in CF due to elevated local velocities
    assert cf_severe.upper_conducting_fraction > healthy.upper_conducting_fraction
    # Bronchial delivery index drops under severe obstruction
    assert cf_severe.bronchial_delivery_index < healthy.bronchial_delivery_index


def test_delivery_platform_benchmark():
    """Verify vibrating mesh nebulizer achieves superior bronchial delivery compared to coarse spray."""
    bench = benchmark_delivery_platforms(inspiratory_flow_l_min=30.0, cf_occlusion_factor=0.20)
    assert "vx522_vibrating_mesh" in bench
    assert "jet_nebulizer_standard" in bench
    assert "coarse_spray_control" in bench

    bdi_mesh = bench["vx522_vibrating_mesh"]["bdi"]
    bdi_coarse = bench["coarse_spray_control"]["bdi"]

    assert bdi_mesh > bdi_coarse
    assert bench["coarse_spray_control"]["small_conducting_pct"] < 1.0


def test_falsification_negative_controls():
    """Verify AH-10 negative control suite passes deterministically."""
    fals_res = evaluate_falsification_negative_control()
    assert fals_res["status"] == "PASS"
    assert fals_res["all_negative_controls_satisfied"] is True
    assert fals_res["coarse_extrathoracic_upper_pct"] > 90.0
    assert fals_res["coarse_small_conducting_pct"] < 2.0
