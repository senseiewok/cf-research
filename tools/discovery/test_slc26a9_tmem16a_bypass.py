#!/usr/bin/env python3
"""Unit tests for Constitutive SLC26A9 & TMEM16A Dual-Channel Bypass Model (T-0256).

Verifies:
1. Trial provenance: 100% of cited clinical trials (NCT00846781, NCT04375514) verified against
   tools/trial_atlas/snapshot.json.
2. SLC26A9 constitutive, uncoupled anion conductance and permeability ratio (PDB: 7CH1).
3. TMEM16A allosteric voltage- and Ca2+-dependent cooperative gating (PDB: 5OYB, 6QP6).
4. Self-consistent Kirchhoff apical membrane potential V_m and current conservation.
5. 24-hour longitudinal ASL kinetics: Denufosol desensitization collapse vs dual-channel durability.
6. Henderson-Hasselbalch ASL pH buffering and bicarbonate restoration (pH >= 7.1).
7. Dual-channel and triple-combination synergy achieving physiological ASL height (>= 7.0 um).
8. Anti-hallucination compliance, NLM disclaimer, and research boundary adherence.
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

from slc26a9_tmem16a_bypass import (
    DEFAULT_ASL_COLLAPSED_UM,
    DEFAULT_ASL_THRESHOLD_UM,
    DEFAULT_BASAL_CA_NM,
    DEFAULT_DENUFOSOL_KD_UM,
    DEFAULT_DENUFOSOL_T_HALF_MIN,
    DEFAULT_PEAK_DENUFOSOL_CA_UM,
    DEFAULT_SLC_DENSITY_CFTR_NULL,
    DEFAULT_SLC_DENSITY_RESCUED,
    DEFAULT_SLC_GAMMA_PS,
    DEFAULT_SLC_PERM_RATIO_HCO3_CL,
    DEFAULT_SLC_PO_CONSTITUTIVE,
    DEFAULT_TARGET_PH_THRESHOLD,
    DEFAULT_TMEM_ALLOSTERIC_ALPHA,
    DEFAULT_TMEM_HILL_NH,
    DEFAULT_TMEM_K0_NM,
    NLM_REGISTRY_DISCLAIMER,
    RESEARCH_BOUNDARY,
    TRIAL_SURVEILLANCE_CATALOG,
    compute_asl_ph,
    compute_slc26a9_conductance,
    compute_tmem16a_apparent_khalf_nm,
    compute_tmem16a_conductance,
    compute_tmem16a_po,
    evaluate_regimen_suite,
    simulate_longitudinal_24h_asl,
    solve_apical_circuit_equilibrium,
)


class TestSLC26A9TMEM16ABypass:
    """Test suite for coupled SLC26A9 constitutive and TMEM16A allosteric anion conductance in CFTR-null airways."""

    def test_01_trial_provenance_and_snapshot_grounding(self) -> None:
        """Every cited NCT ID must exist in tools/trial_atlas/snapshot.json."""
        snapshot_path = REPO_ROOT / "tools" / "trial_atlas" / "snapshot.json"
        assert snapshot_path.exists(), f"Snapshot file not found: {snapshot_path}"

        with open(snapshot_path, "r", encoding="utf-8") as f:
            snap = json.load(f)

        study_ncts = {s["nct_id"] for s in snap.get("studies", []) if "nct_id" in s}

        expected_ncts = [
            "NCT00846781",  # Denufosol Phase 3 extension (P08642)
            "NCT04375514",  # ARO-ENaC siRNA Phase 1/2
        ]

        for nct in expected_ncts:
            assert nct in study_ncts, f"Cited trial {nct} missing from snapshot.json"
            assert nct in TRIAL_SURVEILLANCE_CATALOG, f"Trial {nct} missing from model catalog"

    def test_02_slc26a9_constitutive_conductance_and_selectivity(self) -> None:
        """SLC26A9 exhibits constitutive open probability and fixed anion permeability ratio."""
        g_cl, g_hco3 = compute_slc26a9_conductance(
            channel_density=DEFAULT_SLC_DENSITY_RESCUED,
            gamma_ps=DEFAULT_SLC_GAMMA_PS,
            po=DEFAULT_SLC_PO_CONSTITUTIVE,
            perm_ratio_hco3=DEFAULT_SLC_PERM_RATIO_HCO3_CL,
        )
        expected_g_cl = DEFAULT_SLC_DENSITY_RESCUED * DEFAULT_SLC_GAMMA_PS * DEFAULT_SLC_PO_CONSTITUTIVE
        assert pytest.approx(g_cl, rel=1e-3) == expected_g_cl
        assert pytest.approx(g_hco3 / g_cl, rel=1e-3) == DEFAULT_SLC_PERM_RATIO_HCO3_CL

        # Verify basal CFTR-null density (unrescued) is lower than rescued
        g_cl_null, _ = compute_slc26a9_conductance(channel_density=DEFAULT_SLC_DENSITY_CFTR_NULL)
        assert g_cl_null < g_cl

    def test_03_tmem16a_allosteric_potentiation_and_cooperative_gating(self) -> None:
        """TMEM16A exhibits allosteric potentiator sensitization and voltage-dependent K_1/2."""
        # Baseline at 0 mV without potentiator
        khalf_base = compute_tmem16a_apparent_khalf_nm(vm_mv=0.0, potentiator_conc_nm=0.0)
        assert pytest.approx(khalf_base, rel=1e-3) == DEFAULT_TMEM_K0_NM

        # Depolarization shifts K_1/2 left (increases Ca2+ affinity)
        khalf_depol = compute_tmem16a_apparent_khalf_nm(vm_mv=+20.0, potentiator_conc_nm=0.0)
        assert khalf_depol < khalf_base

        # Allosteric potentiator dramatically sensitizes channel to resting Ca2+
        khalf_pot = compute_tmem16a_apparent_khalf_nm(vm_mv=-40.0, potentiator_conc_nm=100.0)
        khalf_unpot = compute_tmem16a_apparent_khalf_nm(vm_mv=-40.0, potentiator_conc_nm=0.0)
        assert khalf_pot < khalf_unpot * 0.40

        # Resting calcium (80 nM): unpotentiated P_o is near 0, potentiated is significantly elevated
        po_unpot = compute_tmem16a_po(ca_conc_nm=DEFAULT_BASAL_CA_NM, potentiator_conc_nm=0.0)
        po_pot = compute_tmem16a_po(ca_conc_nm=DEFAULT_BASAL_CA_NM, potentiator_conc_nm=100.0)
        assert po_unpot < 0.01
        assert po_pot > 0.015

    def test_04_circuit_equilibrium_and_membrane_potential(self) -> None:
        """Solves self-consistent Kirchhoff apical circuit with zero net current divergence."""
        state = solve_apical_circuit_equilibrium(
            slc26a9_density=DEFAULT_SLC_DENSITY_RESCUED,
            slc26a9_rescued=True,
            ca_conc_nm=DEFAULT_BASAL_CA_NM,
            tmem_potentiator_nm=100.0,
        )
        assert abs(state.net_apical_current_pa) < 0.05
        assert -75.0 <= state.apical_vm_mv <= -15.0

    def test_05_henderson_hasselbalch_asl_ph_buffering(self) -> None:
        """Verifies ASL pH calculation and bicarbonate buffering restoration."""
        # Baseline CF: acidic pH ~6.84
        ph_cf = compute_asl_ph(6.8)
        assert ph_cf < 6.90

        # Normal airway: physiological alkaline pH ~7.35
        ph_normal = compute_asl_ph(22.0)
        assert ph_normal > 7.30

        # Dual bypass elevates ASL bicarbonate and restores pH >= 7.10
        state = solve_apical_circuit_equilibrium(
            slc26a9_rescued=True,
            tmem_potentiator_nm=100.0,
        )
        assert state.asl_ph >= DEFAULT_TARGET_PH_THRESHOLD

    def test_06_regimen_suite_comparative_synergy(self) -> None:
        """Compares baseline, monotherapies, and dual/triple combination bypass regimens."""
        regimens = evaluate_regimen_suite()

        base = regimens["baseline_cftr_null"]
        assert base.steady_state_asl_h_um < 4.0
        assert not base.restores_clearance
        assert not base.restores_buffering

        denuf = regimens["denufosol_tid"]
        # Denufosol has low 24h durability due to rapid degradation
        assert denuf.durability_24h_pct < 25.0
        assert not denuf.restores_clearance

        dual = regimens["dual_channel_bypass"]
        assert dual.steady_state_asl_h_um >= DEFAULT_ASL_THRESHOLD_UM
        assert dual.asl_ph >= DEFAULT_TARGET_PH_THRESHOLD
        assert dual.restores_clearance
        assert dual.restores_buffering
        assert dual.durability_24h_pct == 100.0

        triple = regimens["triple_combination_bypass"]
        assert triple.steady_state_asl_h_um > dual.steady_state_asl_h_um
        assert triple.asl_ph >= dual.asl_ph
        assert triple.restores_clearance
        assert triple.restores_buffering

    def test_07_longitudinal_24h_kinetics_denufosol_vs_dual_bypass(self) -> None:
        """Denufosol exhibits periodic inter-dose collapses below 7.0 um, whereas dual bypass stays >= 7.0 um."""
        kinetics_denuf = simulate_longitudinal_24h_asl("denufosol_tid", time_steps_hours=24)
        kinetics_dual = simulate_longitudinal_24h_asl("dual_channel_bypass", time_steps_hours=24)

        # Denufosol troughs fall below 7.0 um (e.g., at 7h just before next dose)
        trough_denuf = min(pt["asl_height_um"] for pt in kinetics_denuf[16:])
        assert trough_denuf < 6.0

        # Dual bypass trajectory expands to >= 7.0 um and remains sustained throughout steady-state
        late_dual_heights = [pt["asl_height_um"] for pt in kinetics_dual if pt["time_hours"] >= 6.0]
        assert all(h >= DEFAULT_ASL_THRESHOLD_UM for h in late_dual_heights)

    def test_08_disclaimers_and_anti_hallucination_compliance(self) -> None:
        """Validates NLM disclaimer and scientific research boundary declarations."""
        assert "ClinicalTrials.gov does not evaluate safety or efficacy" in NLM_REGISTRY_DISCLAIMER
        assert "listing on the registry does not denote government endorsement" in NLM_REGISTRY_DISCLAIMER
        assert "research, not medical advice" in RESEARCH_BOUNDARY
