#!/usr/bin/env python3
"""Unit tests for Coupled ENaC-TMEM16A Transepithelial Equilibrium & Purinergic ASL Model (T-0252).

Verifies:
1. Trial provenance: 100% of cited clinical trials (NCT00846781, NCT00425165, NCT00625612,
   NCT03647228, NCT04375514) verified against tools/trial_atlas/snapshot.json.
2. P2Y2 purinergic receptor occupancy and calcium transient kinetics.
3. TMEM16A voltage-dependent and allosteric Ca2+ gating (PDB: 5OYB, 6QP6).
4. ENaC proteolytic activation and purinergic feedback suppression (PDB: 6BQN, 6WTH).
5. Self-consistent apical Kirchhoff circuit balance and membrane potential V_m.
6. 24-hour longitudinal ASL dynamics: Denufosol Phase 3 inter-dose collapse vs dual-target synergy.
7. Disclaimers and anti-hallucination compliance.
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path
import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
DISCOVERY_DIR = Path(__file__).resolve().parent

if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(DISCOVERY_DIR) not in sys.path:
    sys.path.insert(0, str(DISCOVERY_DIR))

from enac_tmem16a_equilibrium import (
    DEFAULT_ASL_COLLAPSED_UM,
    DEFAULT_ASL_PHYSIOLOGICAL_UM,
    DEFAULT_ASL_TARGET_UM,
    DEFAULT_BASAL_CA_NM,
    DEFAULT_DENUFOSOL_KD_UM,
    DEFAULT_DENUFOSOL_T_HALF_MIN,
    DEFAULT_ENAC_KINH_CA_NM,
    DEFAULT_ENAC_PO_PROTEOLYTIC_CF,
    DEFAULT_PEAK_PURINERGIC_CA_UM,
    DEFAULT_TMEM_HILL_NH,
    DEFAULT_TMEM_K0_NM,
    NLM_REGISTRY_DISCLAIMER,
    RESEARCH_BOUNDARY,
    TRIAL_SURVEILLANCE_CATALOG,
    compute_cytosolic_calcium_nm,
    compute_enac_po,
    compute_p2y2_occupancy,
    compute_tmem16a_apparent_khalf_nm,
    compute_tmem16a_po,
    run_full_comparative_evaluation,
    simulate_longitudinal_asl_regimen,
    solve_apical_circuit_equilibrium,
)


class TestENaCTMEM16AEquilibrium:
    """Test suite for coupled apical electrophysiology, purinergic signaling, and ASL homeostasis."""

    def test_01_trial_provenance_and_snapshot_grounding(self) -> None:
        """Every cited NCT ID must exist in tools/trial_atlas/snapshot.json."""
        snapshot_path = REPO_ROOT / "tools" / "trial_atlas" / "snapshot.json"
        assert snapshot_path.exists(), f"Snapshot file not found: {snapshot_path}"

        with open(snapshot_path, "r", encoding="utf-8") as f:
            snap = json.load(f)

        study_ncts = {s["nct_id"] for s in snap.get("studies", []) if "nct_id" in s}

        expected_ncts = [
            "NCT00846781",  # Denufosol Phase 3 extension (P08642)
            "NCT00425165",  # Denufosol Phase 3 mild/moderate CF
            "NCT00625612",  # Denufosol Phase 3 TIGER-1
            "NCT03647228",  # IONIS-ENaCRx ASO Phase 1/2
            "NCT04375514",  # ARO-ENaC siRNA Phase 1/2
        ]

        for nct in expected_ncts:
            assert nct in study_ncts, f"Cited trial {nct} missing from snapshot.json"
            assert nct in TRIAL_SURVEILLANCE_CATALOG, f"Trial {nct} missing from model catalog"

    def test_02_p2y2_purinergic_occupancy_and_calcium_kinetics(self) -> None:
        """P2Y2 receptor occupancy follows Hill saturation and drives cytosolic [Ca2+]_i."""
        # Baseline: zero agonist -> zero occupancy, basal resting Ca2+ (80 nM)
        assert compute_p2y2_occupancy(0.0) == 0.0
        ca_basal = compute_cytosolic_calcium_nm(0.0)
        assert pytest.approx(ca_basal, rel=1e-3) == DEFAULT_BASAL_CA_NM

        # At K_d (0.85 uM), occupancy is exactly 50%
        theta_kd = compute_p2y2_occupancy(DEFAULT_DENUFOSOL_KD_UM)
        assert pytest.approx(theta_kd, rel=1e-3) == 0.50

        # At K_d, Ca2+ is halfway between basal (80 nM) and peak (1400 nM)
        expected_mid_ca = DEFAULT_BASAL_CA_NM + (DEFAULT_PEAK_PURINERGIC_CA_UM * 1000.0 - DEFAULT_BASAL_CA_NM) * 0.5
        ca_mid = compute_cytosolic_calcium_nm(DEFAULT_DENUFOSOL_KD_UM)
        assert pytest.approx(ca_mid, rel=1e-3) == expected_mid_ca

        # Saturating agonist -> approaches 100% occupancy and peak Ca2+ (1400 nM)
        theta_sat = compute_p2y2_occupancy(100.0)
        assert pytest.approx(theta_sat, abs=0.01) == 1.0
        ca_sat = compute_cytosolic_calcium_nm(100.0)
        assert pytest.approx(ca_sat, abs=15.0) == DEFAULT_PEAK_PURINERGIC_CA_UM * 1000.0

    def test_03_tmem16a_voltage_and_calcium_gating(self) -> None:
        """TMEM16A open probability reflects dual Ca2+ cooperativity and voltage left-shift."""
        # At resting calcium (80 nM) and neutral voltage, P_o is low (< 0.02)
        po_basal = compute_tmem16a_po(DEFAULT_BASAL_CA_NM, vm_mv=0.0)
        assert po_basal < 0.02

        # P_o increases monotonically with [Ca2+]_i
        po_500nm = compute_tmem16a_po(500.0, vm_mv=0.0)
        po_1400nm = compute_tmem16a_po(1400.0, vm_mv=0.0)
        assert po_1400nm > po_500nm > po_basal
        assert po_1400nm > 0.60  # Substantial activation at peak purinergic calcium

        # Membrane depolarization (e.g. from -40 mV to 0 mV to +40 mV) left-shifts K_1/2
        khalf_neg40 = compute_tmem16a_apparent_khalf_nm(vm_mv=-40.0)
        khalf_zero = compute_tmem16a_apparent_khalf_nm(vm_mv=0.0)
        khalf_pos40 = compute_tmem16a_apparent_khalf_nm(vm_mv=40.0)
        assert khalf_neg40 > khalf_zero > khalf_pos40

        # Allosteric potentiator (ETD001-like) left-shifts K_1/2 at resting calcium
        khalf_potentiated = compute_tmem16a_apparent_khalf_nm(vm_mv=-40.0, potentiator_conc_nm=100.0)
        assert khalf_potentiated < khalf_neg40
        po_potentiated_resting = compute_tmem16a_po(DEFAULT_BASAL_CA_NM, vm_mv=-40.0, potentiator_conc_nm=100.0)
        assert po_potentiated_resting > compute_tmem16a_po(DEFAULT_BASAL_CA_NM, vm_mv=-40.0, potentiator_conc_nm=0.0)

    def test_04_enac_proteolytic_and_purinergic_inhibition(self) -> None:
        """ENaC P_o is hyperactive in CF proteolytic state and suppressed by Ca2+ and DAG/PKC."""
        # Uninhibited CF baseline: cleaved P_o = 0.85
        po_cf_uninhibited = compute_enac_po(ca_conc_nm=DEFAULT_BASAL_CA_NM, dag_occupancy=0.0)
        assert po_cf_uninhibited > 0.80

        # Purinergic surge: high Ca2+ (1400 nM) downregulates ENaC via negative allostery
        po_ca_inhibited = compute_enac_po(ca_conc_nm=1400.0, dag_occupancy=0.0)
        assert po_ca_inhibited < 0.15  # Marked Ca2+-dependent suppression

        # DAG/PKC phosphorylation adds further inhibition
        po_dual_inhibited = compute_enac_po(ca_conc_nm=1400.0, dag_occupancy=0.8)
        assert po_dual_inhibited < po_ca_inhibited

        # RNAi knockdown (70%) reduces baseline open probability proportionally
        po_rnai = compute_enac_po(ca_conc_nm=DEFAULT_BASAL_CA_NM, rna_knockdown_fraction=0.70)
        assert pytest.approx(po_rnai, rel=1e-2) == po_cf_uninhibited * 0.30

    def test_05_apical_circuit_kirchhoff_balance(self) -> None:
        """Kirchhoff current law converges to steady-state apical V_m."""
        # Solve for CF baseline state when ASL is at physiological height (7.0 um)
        state_cf = solve_apical_circuit_equilibrium(
            ca_conc_nm=DEFAULT_BASAL_CA_NM,
            current_asl_um=DEFAULT_ASL_PHYSIOLOGICAL_UM,
        )
        # Apical V_m must lie within physiological bounds (-80 mV to 0 mV)
        assert -80.0 < state_cf.apical_vm_mv < 0.0
        # In CF baseline, ENaC sodium absorption dominates
        assert abs(state_cf.i_na_enac_pa) > abs(state_cf.i_cl_tmem_pa)
        # Secretory net flux is negative (absorptive), driving ASL dehydration
        assert state_cf.net_secretory_ion_flux_pa < 0.0
        assert state_cf.dh_asl_dt_um_per_hour < 0.0

        # Solve for purinergic stimulated state (Ca2+ = 1400 nM)
        state_purinergic = solve_apical_circuit_equilibrium(
            ca_conc_nm=1400.0,
            dag_occupancy=0.9,
            current_asl_um=DEFAULT_ASL_COLLAPSED_UM,
        )
        # TMEM16A chloride current becomes prominent
        assert state_purinergic.g_tmem16a_ps > state_cf.g_tmem16a_ps
        # Net secretory driving force reverses to positive (secretory)
        assert state_purinergic.dh_asl_dt_um_per_hour > 0.0

    def test_06_longitudinal_regimen_denufosol_failure_vs_dual_synergy(self) -> None:
        """24-hr longitudinal simulation captures Denufosol clinical trial failure vs dual synergy."""
        # 1. CF untreated: stays collapsed (~3.2 um), 0% ciliary competence
        sum_cf, pts_cf = simulate_longitudinal_asl_regimen("cf_untreated")
        assert sum_cf.mean_asl_height_um < 4.0
        assert sum_cf.ciliary_competence_pct_time == 0.0

        # 2. WT healthy: maintains >= 7.0 um throughout
        sum_wt, pts_wt = simulate_longitudinal_asl_regimen("wt_healthy")
        assert sum_wt.mean_asl_height_um >= 7.0
        assert sum_wt.ciliary_competence_pct_time == 100.0

        # 3. Denufosol TID (Phase 3 NCT00846781 / NCT00625612):
        # Generates acute transient peaks, but rapid degradation (t_1/2 = 20 min)
        # leads to inter-dose collapse; ciliary competence < 30%
        sum_den, pts_den = simulate_longitudinal_asl_regimen("denufosol_tid")
        assert sum_den.max_asl_height_um >= 7.0  # Acute rehydration peak
        assert sum_den.min_asl_height_um <= 3.5  # Inter-dose collapse
        assert sum_den.ciliary_competence_pct_time < 30.0  # Inadequate 24-hr coverage explains trial failure

        # 4. Synergistic dual-target (ENaC RNAi + TMEM16A potentiator):
        # Sustained ASL >= 7.0 um without calcium toxicity
        sum_syn, pts_syn = simulate_longitudinal_asl_regimen("synergistic_dual_target")
        assert sum_syn.mean_asl_height_um >= 7.0
        assert sum_syn.ciliary_competence_pct_time >= 90.0
        assert sum_syn.peak_cytosolic_ca_nm == DEFAULT_BASAL_CA_NM  # No toxic calcium spike

    def test_07_disclaimers_and_anti_hallucination_metadata(self) -> None:
        """Full comparative output asserts NLM disclaimer, research boundary, and structural anchors."""
        data = run_full_comparative_evaluation()
        assert NLM_REGISTRY_DISCLAIMER in data["nlm_disclaimer"]
        assert "ClinicalTrials.gov does not evaluate safety or efficacy" in data["nlm_disclaimer"]
        assert RESEARCH_BOUNDARY in data["research_boundary"]

        # Check PDB structural grounding
        pdbs = data["structural_anchors"]
        assert "PDB_5OYB" in pdbs
        assert "PDB_6BQN" in pdbs
        assert "PDB_6WTH" in pdbs
        assert "PDB_7L1V" in pdbs

        # Check all 6 standard regimens are simulated
        summaries = data["regimen_summaries"]
        assert len(summaries) == 6
        assert "denufosol_tid" in summaries
        assert "synergistic_dual_target" in summaries
