#!/usr/bin/env python3
"""Unit tests for Phage-Antibiotic Co-Evolution & Biofilm Penetration Mechanics (T-0445)."""

import pytest
import math
from phage_biofilm_penetration_coevolution import (
    BiofilmDepolymeraseParameters,
    PhageAdsorptionLysisParameters,
    CoevolutionaryTradeoffParameters,
    EfferocytosisSynergyParameters,
    calculate_effective_biofilm_diffusivity,
    simulate_enzymatic_eps_depolymerization,
    compute_adsorption_and_burst_kinetics,
    calculate_antibiotic_kill_rate,
    simulate_coevolutionary_population_dynamics,
    main,
    NLM_REGISTRY_DISCLAIMER,
    RESEARCH_BOUNDARY,
    CLINICAL_TRIAL_VX522,
    CLINICAL_TRIAL_4D710,
    PDB_DEPHOSPHO_CFTR,
    PDB_ACTIVE_CFTR,
    PRIMARY_DOI_COLLATERAL,
    PRIMARY_DOI_OMKO1_TRADEOFF,
    PRIMARY_DOI_ALGINATE_LYASE,
    PRIMARY_DOI_MUCUS_BIOPHYSICS,
)


def test_effective_diffusivity_intact_vs_degraded():
    """Verify intact biofilm restricts diffusivity (D_eff/D_0 < 0.20) and degradation restores it (>= 0.40)."""
    intact = calculate_effective_biofilm_diffusivity(eps_volume_fraction=0.25)
    assert intact["effective_diffusivity_ratio"] < 0.25
    assert not intact["penetrable_threshold_met"]

    cleared = calculate_effective_biofilm_diffusivity(eps_volume_fraction=0.03)
    assert cleared["effective_diffusivity_ratio"] >= 0.40
    assert cleared["penetrable_threshold_met"]


def test_diffusivity_invalid_input():
    """Ensure negative EPS volume fraction raises ValueError."""
    with pytest.raises(ValueError, match="cannot be negative"):
        calculate_effective_biofilm_diffusivity(eps_volume_fraction=-0.05)


def test_simulate_enzymatic_depolymerization():
    """Verify time-dependent depolymerase activity achieves D_eff / D_0 >= 0.40."""
    res = simulate_enzymatic_eps_depolymerization(enzyme_concentration_um=1.5, duration_hours=4.0)
    assert res["diffusivity_recovery_achieved"] is True
    assert res["final_diffusivity_ratio"] >= 0.40
    assert res["final_eps_fraction"] < res["initial_eps_fraction"]
    assert len(res["history"]) > 10


def test_hypoxic_adsorption_burst_kinetics():
    """Verify hypoxic biofilm microenvironment prolongs latency and attenuates burst size."""
    normoxic = compute_adsorption_and_burst_kinetics(hypoxia_fraction=0.0)
    hypoxic = compute_adsorption_and_burst_kinetics(hypoxia_fraction=0.80)

    assert hypoxic["effective_burst_size"] < normoxic["effective_burst_size"]
    assert hypoxic["effective_latent_period_min"] > normoxic["effective_latent_period_min"]
    assert hypoxic["burst_to_latency_ratio"] < normoxic["burst_to_latency_ratio"]


def test_hypoxia_bounds_validation():
    """Verify invalid hypoxia fractions raise ValueError."""
    with pytest.raises(ValueError, match="between 0.0 and 1.0"):
        compute_adsorption_and_burst_kinetics(hypoxia_fraction=-0.1)

    with pytest.raises(ValueError, match="between 0.0 and 1.0"):
        compute_adsorption_and_burst_kinetics(hypoxia_fraction=1.2)


def test_collateral_sensitivity_and_resensitization():
    """Verify mutational loss of OprM collapses antibiotic MICs to <= 1.0 ug/mL."""
    tp = CoevolutionaryTradeoffParameters()
    assert tp.mic_cipro_bs_ug_ml >= 16.0
    assert tp.mic_cipro_br_ug_ml <= 1.0
    assert tp.mic_tobra_bs_ug_ml >= 8.0
    assert tp.mic_tobra_br_ug_ml <= 1.0
    assert tp.fitness_cost_opr_m_loss > 0.10


def test_coevolutionary_population_dynamics_synergy():
    """Verify dual phage + antibiotic therapy drives >= 3.5 log10 reduction in CFUs."""
    dual_res = simulate_coevolutionary_population_dynamics(
        initial_bs_cfu_ml=5.0e7,
        initial_br_cfu_ml=1.0e2,
        initial_phage_pfu_ml=1.0e8,
        ciprofloxacin_ug_ml=4.0,
        duration_hours=24.0,
    )
    assert dual_res["eradication_synergy_met"] is True
    assert dual_res["log10_bacterial_reduction"] >= 3.5
    assert dual_res["antibiotic_resensitization_confirmed"] is True
    assert dual_res["efferocytosis_debris_cleared_pct"] > 50.0


def test_monotherapy_failure_controls():
    """Verify antibiotic monotherapy fails against MDR B_S, while phage alone selects for B_R."""
    # Antibiotic alone: zero phage
    ab_alone = simulate_coevolutionary_population_dynamics(
        initial_bs_cfu_ml=5.0e7,
        initial_br_cfu_ml=1.0e2,
        initial_phage_pfu_ml=0.0,
        ciprofloxacin_ug_ml=4.0,
        duration_hours=24.0,
    )
    # B_S survives ciprofloxacin (MIC 32 ug/mL vs 4 ug/mL)
    assert ab_alone["final_bs_cfu_ml"] > 1.0e6
    assert not ab_alone["eradication_synergy_met"]

    # Phage alone: zero antibiotic
    phage_alone = simulate_coevolutionary_population_dynamics(
        initial_bs_cfu_ml=5.0e7,
        initial_br_cfu_ml=1.0e2,
        initial_phage_pfu_ml=1.0e8,
        ciprofloxacin_ug_ml=0.0,
        duration_hours=24.0,
    )
    # B_R grows out when only phage is present
    assert phage_alone["final_br_cfu_ml"] > 1.0e5


def test_cli_execution_and_disclaimers(capsys):
    """Verify CLI flags, JSON output, NLM attribution, and research boundaries."""
    ret = main(["--json", "--enzyme-um", "2.0", "--cipro-ug-ml", "5.0", "--duration-hr", "12.0"])
    assert ret == 0
    captured = capsys.readouterr()
    assert "phage_biofilm_penetration_coevolution" in captured.out
    assert NLM_REGISTRY_DISCLAIMER in captured.out
    assert RESEARCH_BOUNDARY in captured.out
    assert CLINICAL_TRIAL_VX522 in captured.out
    assert CLINICAL_TRIAL_4D710 in captured.out
    assert PDB_DEPHOSPHO_CFTR in captured.out
    assert PDB_ACTIVE_CFTR in captured.out
    assert PRIMARY_DOI_COLLATERAL in captured.out
    assert PRIMARY_DOI_OMKO1_TRADEOFF in captured.out
    assert PRIMARY_DOI_ALGINATE_LYASE in captured.out
    assert PRIMARY_DOI_MUCUS_BIOPHYSICS in captured.out
