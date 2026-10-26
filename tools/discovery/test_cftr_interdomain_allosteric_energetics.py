#!/usr/bin/env python3
"""Unit test suite for CFTR Inter-Domain Allosteric Energetics Engine (T-0447)."""
from __future__ import annotations

import sys
from pathlib import Path
import pytest

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from cftr_interdomain_allosteric_energetics import (
    CFTRAllostericEnergeticsEngine,
    CLINICAL_TRIAL_VX522,
    CLINICAL_TRIAL_4D710,
    DOI_AIRWAY_BIOPHYSICS,
    DOI_ETI_PHASE3_TRIAL,
    DOI_POTENTIATOR_HOTSPOT,
    DOI_W1282X_RESCUE_BARRIER,
    NLM_ATTRIBUTION,
    PDB_ACTIVE_CFTR,
    PDB_DEPHOSPHO_CFTR,
    RESEARCH_BOUNDARY,
)


@pytest.fixture
def engine() -> CFTRAllostericEnergeticsEngine:
    """Fixture providing initialized allosteric energetics engine."""
    return CFTRAllostericEnergeticsEngine()


def test_interdomain_packing_conformations(engine: CFTRAllostericEnergeticsEngine) -> None:
    """Verify inter-domain packing metrics between 5UAK (active) and 6MSM (inactive)."""
    state_5uak = engine.compute_interdomain_coupling("5UAK")
    state_6msm = engine.compute_interdomain_coupling("6MSM")

    # 5UAK active ATP-bound conformation
    assert state_5uak["conformation_id"] == "5UAK"
    assert state_5uak["pdb_id"] == "5UAK"
    assert state_5uak["total_buried_surface_area_a2"] > 5000.0
    assert state_5uak["total_coupling_delta_g_kcal_mol"] <= -20.0
    assert state_5uak["tmd1_nbd1"]["delta_g_coupling_kcal_mol"] <= -9.0
    assert state_5uak["tmd2_nbd2"]["delta_g_coupling_kcal_mol"] <= -8.0

    # 6MSM inactive uncoupled conformation
    assert state_6msm["conformation_id"] == "6MSM"
    assert state_6msm["pdb_id"] == "6MSM"
    assert state_6msm["total_buried_surface_area_a2"] < 4000.0
    assert state_6msm["total_coupling_delta_g_kcal_mol"] >= -15.0

    # Differential coupling
    comp = engine.compare_conformations()
    assert comp["buried_surface_area_difference_a2"] > 2000.0
    assert comp["coupling_free_energy_stabilization_kcal_mol"] <= -12.0


def test_corrector_stabilization_threshold(engine: CFTRAllostericEnergeticsEngine) -> None:
    """Verify corrector stabilization Delta-Delta G <= -3.5 kcal/mol across Class II variants."""
    class2_variants = ["F508del", "N1303K", "L206W", "V520F"]

    for var in class2_variants:
        resp_trikafta = engine.simulate_corrector_stabilization(var, "trikafta")
        assert resp_trikafta.stabilization_threshold_met is True
        assert resp_trikafta.delta_delta_g_stabilization_kcal_mol <= -3.50
        assert resp_trikafta.restored_mature_folding_fraction >= 0.25

        # Check monotherapy vs combination synergy
        resp_tez = engine.simulate_corrector_stabilization(var, "tezacaftor")
        resp_elx = engine.simulate_corrector_stabilization(var, "elexacaftor")
        assert resp_tez.delta_delta_g_stabilization_kcal_mol > resp_trikafta.delta_delta_g_stabilization_kcal_mol
        assert resp_elx.delta_delta_g_stabilization_kcal_mol > resp_trikafta.delta_delta_g_stabilization_kcal_mol
        assert resp_trikafta.restored_mature_folding_fraction > resp_tez.restored_mature_folding_fraction
        assert resp_trikafta.restored_mature_folding_fraction > resp_elx.restored_mature_folding_fraction


def test_truncation_erad_kinetics(engine: CFTRAllostericEnergeticsEngine) -> None:
    """Verify ERAD clearance and surface half-life dynamics for truncated W1282X CFTR."""
    wt_kinetics = engine.model_truncation_erad_kinetics(is_truncated=False, treated_with_correctors=False)
    w1282x_untreated = engine.model_truncation_erad_kinetics(is_truncated=True, treated_with_correctors=False)
    w1282x_treated = engine.model_truncation_erad_kinetics(is_truncated=True, treated_with_correctors=True)

    # Untreated truncated CFTR has severe ERAD entrapment and short surface persistence
    assert w1282x_untreated.is_truncated is True
    assert w1282x_untreated.k_erad_per_min > wt_kinetics.k_erad_per_min * 5.0
    assert w1282x_untreated.erad_escape_probability < 0.05
    assert w1282x_untreated.surface_half_life_hours < 3.0

    # Corrector stabilization significantly boosts ERAD escape
    assert w1282x_treated.erad_escape_probability > w1282x_untreated.erad_escape_probability * 3.0
    assert w1282x_treated.surface_half_life_hours > w1282x_untreated.surface_half_life_hours


def test_potentiator_gating_coupling(engine: CFTRAllostericEnergeticsEngine) -> None:
    """Verify allosteric coupling between TM8 hotspot binding and channel gating."""
    g551d_resp = engine.calculate_potentiator_gating_coupling("G551D", potentiator_conc_nm=1000.0)
    wt_resp = engine.calculate_potentiator_gating_coupling("WT", potentiator_conc_nm=1000.0)

    # TM8 hotspot thermodynamic constants
    assert g551d_resp.kd_open_nm == 8.5
    assert g551d_resp.kd_closed_nm == 120.0
    assert g551d_resp.allosteric_ratio_alpha > 14.0
    assert g551d_resp.delta_delta_g_gating_kcal_mol <= -1.50

    # Potentiator stimulates G551D Po by over 15-fold
    assert g551d_resp.basal_po < 0.01
    assert g551d_resp.potentiated_po >= 0.08
    assert g551d_resp.fold_stimulation >= 15.0

    # WT Po reaches high open probability
    assert wt_resp.potentiated_po > 0.60


def test_comprehensive_suite_grounding(engine: CFTRAllostericEnergeticsEngine) -> None:
    """Verify structural, literature, and regulatory grounding compliance."""
    suite = engine.run_comprehensive_thermodynamic_suite()

    assert suite["task_id"] == "T-0447"
    assert PDB_ACTIVE_CFTR in suite["structural_grounding"]
    assert PDB_DEPHOSPHO_CFTR in suite["structural_grounding"]
    assert DOI_POTENTIATOR_HOTSPOT in suite["primary_literature"]
    assert DOI_W1282X_RESCUE_BARRIER in suite["primary_literature"]
    assert DOI_ETI_PHASE3_TRIAL in suite["primary_literature"]
    assert DOI_AIRWAY_BIOPHYSICS in suite["primary_literature"]
    assert CLINICAL_TRIAL_VX522 in suite["clinical_trials"]
    assert CLINICAL_TRIAL_4D710 in suite["clinical_trials"]
    assert "ClinicalTrials.gov does not evaluate safety or efficacy" in suite["nlms_attribution"]
    assert "research, not medical advice" in suite["research_boundary"]
