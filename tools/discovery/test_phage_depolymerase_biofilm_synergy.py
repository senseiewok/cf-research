#!/usr/bin/env python3
"""Unit tests for Phage Alginate Depolymerase & Tobramycin Penetration Synergy Model (T-0402).

Validates:
  1. Primary literature, RCSB PDB 6QW8 structure, and clinical registry grounding.
  2. Alginate depolymerase Michaelis-Menten cleavage kinetics (k_cat / K_m in [10^4, 10^6] M^-1 s^-1).
  3. Diffusivity hindrance (D_eff,intact / D0 < 0.05) and recovery (D_eff / D0 >= 0.70).
  4. 1D spatiotemporal reaction-diffusion across biofilm thickness (z >= 100 um).
  5. Bactericidal synergy quantification (>= 3.0 log10 CFU synergistic reduction over monotherapy).
  6. Parameter sensitivity and enzyme dose-response monotonicity.
  7. Verification benchmark suite execution and CLI integration.
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

from phage_depolymerase_biofilm_synergy import (
    CLINICAL_TRIAL_REGISTRY,
    DEFAULT_BIOFILM_THICKNESS_UM,
    DEFAULT_C_BULK_TOBRA_UM,
    DEFAULT_D0_TOBRAMYCIN_UM2_S,
    DEFAULT_EC50_TOBRA_UM,
    DEFAULT_ENZYME_CONC_UM,
    DEFAULT_EPS_VOLUME_FRACTION,
    DEFAULT_INITIAL_ALGINATE_SITES_UM,
    DEFAULT_KD_BINDING_UM,
    NLM_ATTRIBUTION,
    PRIMARY_DOI_LYASE,
    PRIMARY_DOI_SYNERGY,
    RCSB_PDB_STRUCTURE,
    RESEARCH_BOUNDARY,
    compute_alginate_lyase_kinetics,
    compute_diffusivity_recovery,
    evaluate_synergy,
    run_benchmark_suite,
    simulate_biofilm_reaction_diffusion,
)


def test_scientific_grounding_and_disclaimers():
    """Verify primary literature DOIs, PDB accession, and NLM attribution."""
    assert RCSB_PDB_STRUCTURE == "6QW8"
    assert PRIMARY_DOI_SYNERGY == "10.1038/s41591-025-03678-8"
    assert PRIMARY_DOI_LYASE == "10.1111/j.1600-0463.2006.apm_356.x"
    assert CLINICAL_TRIAL_REGISTRY == "NCT04596319"
    assert "ClinicalTrials.gov does not evaluate safety or efficacy" in NLM_ATTRIBUTION
    assert "research, not medical advice" in RESEARCH_BOUNDARY


def test_alginate_lyase_kinetics():
    """Verify Michaelis-Menten depolymerase cleavage and catalytic efficiency in [10^4, 10^6] M^-1 s^-1."""
    res = compute_alginate_lyase_kinetics(
        k_cat_per_s=25.0,
        k_m_um=50.0,
        enzyme_conc_um=1.0,
        substrate_conc_um=600.0,
    )
    # k_cat / K_m = 25 / (50e-6) = 5.0e5 M^-1 s^-1
    assert 1.0e4 <= res.catalytic_efficiency_m_s <= 1.0e6
    assert math.isclose(res.catalytic_efficiency_m_s, 5.0e5, rel_tol=1e-3)
    assert res.initial_rate_um_per_s > 0.0
    assert res.depolymerization_half_life_s > 0.0

    # Rate scales with enzyme concentration
    res_double_e = compute_alginate_lyase_kinetics(enzyme_conc_um=2.0)
    assert res_double_e.initial_rate_um_per_s > res.initial_rate_um_per_s

    with pytest.raises(ValueError):
        compute_alginate_lyase_kinetics(k_cat_per_s=-1.0)


def test_diffusivity_recovery():
    """Verify intact diffusivity hindrance < 0.05 and recovery >= 0.70 upon stripping."""
    assert 15.0 <= DEFAULT_KD_BINDING_UM <= 30.0

    res = compute_diffusivity_recovery(
        d0_um2_s=DEFAULT_D0_TOBRAMYCIN_UM2_S,
        eps_volume_fraction=DEFAULT_EPS_VOLUME_FRACTION,
        initial_alginate_sites_um=DEFAULT_INITIAL_ALGINATE_SITES_UM,
        kd_um=DEFAULT_KD_BINDING_UM,
        residual_alginate_fraction=0.03,
    )

    # Intact biofilm diffusivity severely hindered: D_eff / D0 < 0.05
    assert res.d_eff_intact_ratio < 0.05
    assert res.d_eff_intact_um2_s < 0.05 * DEFAULT_D0_TOBRAMYCIN_UM2_S

    # De-densified biofilm restores diffusivity: D_eff / D0 >= 0.70
    assert res.d_eff_recovered_ratio >= 0.70
    assert res.d_eff_recovered_um2_s >= 0.70 * DEFAULT_D0_TOBRAMYCIN_UM2_S

    # Monotonic relationship between intact sites and diffusivity
    res_more_intact = compute_diffusivity_recovery(residual_alginate_fraction=0.50)
    assert res_more_intact.d_eff_recovered_ratio < res.d_eff_recovered_ratio


def test_spatiotemporal_reaction_diffusion():
    """Verify 1D reaction-diffusion across biofilm depth z >= 100 um."""
    thickness = 100.0
    assert thickness >= DEFAULT_BIOFILM_THICKNESS_UM

    prof_ctrl = simulate_biofilm_reaction_diffusion("control", duration_s=1800.0, thickness_um=thickness, n_nodes=21)
    prof_tobra = simulate_biofilm_reaction_diffusion("tobramycin", duration_s=1800.0, thickness_um=thickness, n_nodes=21)
    prof_combo = simulate_biofilm_reaction_diffusion("combination", duration_s=1800.0, thickness_um=thickness, n_nodes=21)

    assert len(prof_ctrl.z_nodes_um) == 21
    assert prof_ctrl.z_nodes_um[0] == 0.0
    assert math.isclose(prof_ctrl.z_nodes_um[-1], thickness, rel_tol=1e-3)

    # Apical boundary condition maintained
    assert prof_tobra.c_tobra_free_um[0] == DEFAULT_C_BULK_TOBRA_UM
    assert prof_combo.c_tobra_free_um[0] == DEFAULT_C_BULK_TOBRA_UM

    # Tobramycin penetrates deeper in combination than monotherapy due to alginate stripping
    deep_node = -1  # at z = L
    assert prof_combo.c_tobra_free_um[deep_node] > prof_tobra.c_tobra_free_um[deep_node]

    # Depolymerase strips alginate sites across the profile
    assert prof_combo.alginate_sites_um[0] < prof_tobra.alginate_sites_um[0]


def test_bactericidal_synergy_quantification():
    """Verify >= 3.0 log10 CFU synergistic reduction for depolymerase + tobramycin over monotherapy."""
    synergy = evaluate_synergy(
        duration_s=3600.0,
        thickness_um=DEFAULT_BIOFILM_THICKNESS_UM,
        c_bulk_tobra_um=DEFAULT_C_BULK_TOBRA_UM,
        enzyme_conc_um=DEFAULT_ENZYME_CONC_UM,
    )

    # Monotherapies have limited kill across thick biofilm
    assert synergy.log10_kill_tobra_mono < 2.5
    assert synergy.log10_kill_depol_mono < 1.0

    # Combination therapy achieves profound bactericidal kill
    assert synergy.log10_kill_combo >= 3.5

    # Synergistic reduction over best monotherapy >= 3.0 log10 CFU
    assert synergy.synergistic_reduction_log10 >= 3.0

    # Combination yields higher fraction of biofilm reaching effective therapeutic concentration
    assert synergy.penetration_depth_ratio_combo > synergy.penetration_depth_ratio_tobra


def test_parameter_sensitivity_and_enzyme_dose():
    """Verify dose-response monotonicity of depolymerase and antibiotic."""
    # Higher enzyme concentration produces faster stripping and greater synergy
    synergy_low_e = evaluate_synergy(duration_s=1800.0, enzyme_conc_um=0.2)
    synergy_high_e = evaluate_synergy(duration_s=1800.0, enzyme_conc_um=2.0)

    assert synergy_high_e.log10_kill_combo > synergy_low_e.log10_kill_combo

    # Higher bulk tobramycin concentration increases killing
    prof_low = simulate_biofilm_reaction_diffusion("tobramycin", duration_s=1800.0, c_bulk_tobra_um=100.0)
    prof_high = simulate_biofilm_reaction_diffusion("tobramycin", duration_s=1800.0, c_bulk_tobra_um=800.0)
    avg_cfu_low = sum(prof_low.bacterial_density_cfu_ml) / len(prof_low.bacterial_density_cfu_ml)
    avg_cfu_high = sum(prof_high.bacterial_density_cfu_ml) / len(prof_high.bacterial_density_cfu_ml)
    assert avg_cfu_high < avg_cfu_low


def test_benchmark_suite_execution():
    """Verify run_benchmark_suite returns full structured result satisfying all acceptance checks."""
    bench = run_benchmark_suite()
    assert bench["task_id"] == "T-0402"
    assert bench["structural_model"] == "6QW8"
    assert bench["verification_summary"]["status"] == "PASS"
    assert bench["verification_summary"]["catalytic_efficiency_in_range"] is True
    assert bench["verification_summary"]["intact_diffusivity_hindered_under_0_05"] is True
    assert bench["verification_summary"]["recovered_diffusivity_over_0_70"] is True
    assert bench["verification_summary"]["synergistic_reduction_over_3_0_log10"] is True
