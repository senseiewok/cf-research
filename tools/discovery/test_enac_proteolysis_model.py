#!/usr/bin/env python3
"""Unit tests for enac_proteolysis_model.py (T-0245)."""

from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))

from tools.discovery.enac_proteolysis_model import (
    DEFAULT_PO_0,
    DEFAULT_PO_ALPHA,
    DEFAULT_PO_ALPHAGAMMA,
    DEFAULT_NE_CF_NM,
    DEFAULT_NE_WT_NM,
    DEFAULT_PRSS8_NM,
    DEFAULT_FURIN_NM,
    DEFAULT_TARGET_ASL_UM,
    PHYSIOLOGICAL_ASL_THRESHOLD_UM,
    COLLAPSED_ASL_CF_UM,
    NLM_REGISTRY_DISCLAIMER,
    RESEARCH_BOUNDARY,
    calculate_cleavage_kinetics,
    simulate_cleavage_ode_dynamics,
    calculate_electrophysiology,
    calculate_asl_homeostasis,
    simulate_enac_proteolysis_scenario,
    get_standard_scenarios,
    evaluate_dual_target_synergy_matrix,
)


class TestENaCProteolysisModel(unittest.TestCase):
    """Test suite for ENaC multi-site proteolytic activation kinetics and ASL homeostasis."""

    def test_01_trial_provenance_and_snapshot_grounding(self):
        """Clinical trial IDs cited must exist in tools/trial_atlas/snapshot.json."""
        snapshot_file = REPO_ROOT / "tools" / "trial_atlas" / "snapshot.json"
        self.assertTrue(snapshot_file.exists(), "snapshot.json must exist in tools/trial_atlas/")

        with open(snapshot_file, "r", encoding="utf-8") as f:
            snapshot_data = json.load(f)

        studies = snapshot_data.get("studies", [])
        study_ncts = {s["nct_id"] for s in studies if "nct_id" in s}

        # Check IONIS-ENaCRx (ASO) and ARO-ENaC (siRNA) and BI 1265162 (pore blocker)
        self.assertIn("NCT03647228", study_ncts, "IONIS-ENaCRx NCT03647228 must be in snapshot.json")
        self.assertIn("NCT04375514", study_ncts, "ARO-ENaC NCT04375514 must be in snapshot.json")
        self.assertIn("NCT02343445", study_ncts, "BI 1265162 NCT02343445 must be in snapshot.json")

    def test_02_cleavage_kinetics_steady_state_fractions(self):
        """Steady-state cleavage fractions must conserve mass and reflect protease concentrations."""
        # Baseline / WT: NE = 0
        wt_cleavage = calculate_cleavage_kinetics(
            ne_conc_nm=DEFAULT_NE_WT_NM,
            prss8_conc_nm=DEFAULT_PRSS8_NM,
            furin_conc_nm=DEFAULT_FURIN_NM,
        )
        # Sum of fractions must equal 1.0
        self.assertAlmostEqual(
            wt_cleavage.f_0 + wt_cleavage.f_alpha + wt_cleavage.f_alphagamma,
            1.0,
            places=6,
        )
        # Mean Po must be bounded within [P_o,0, P_o,alphagamma]
        self.assertGreaterEqual(wt_cleavage.mean_po, DEFAULT_PO_0)
        self.assertLessEqual(wt_cleavage.mean_po, DEFAULT_PO_ALPHAGAMMA)

        # CF sputum with severe neutrophil elastase (180 nM)
        cf_cleavage = calculate_cleavage_kinetics(
            ne_conc_nm=DEFAULT_NE_CF_NM,
            prss8_conc_nm=DEFAULT_PRSS8_NM,
            furin_conc_nm=DEFAULT_FURIN_NM,
        )
        self.assertAlmostEqual(
            cf_cleavage.f_0 + cf_cleavage.f_alpha + cf_cleavage.f_alphagamma,
            1.0,
            places=6,
        )
        # In CF with high NE, gamma-cleaved channel fraction should dominate
        self.assertGreater(cf_cleavage.f_alphagamma, 0.70)
        self.assertGreater(cf_cleavage.mean_po, 0.65)
        self.assertGreater(cf_cleavage.mean_po, wt_cleavage.mean_po)

    def test_03_zero_protease_boundary_conditions(self):
        """When furin and elastase/prostasin are zero, channels stay 100% uncleaved."""
        cleavage_zero = calculate_cleavage_kinetics(
            ne_conc_nm=0.0,
            prss8_conc_nm=0.0,
            furin_conc_nm=0.0,
        )
        self.assertAlmostEqual(cleavage_zero.f_0, 1.0, places=6)
        self.assertAlmostEqual(cleavage_zero.f_alpha, 0.0, places=6)
        self.assertAlmostEqual(cleavage_zero.f_alphagamma, 0.0, places=6)
        self.assertAlmostEqual(cleavage_zero.mean_po, DEFAULT_PO_0, places=6)

    def test_04_ode_dynamics_convergence_to_analytical_steady_state(self):
        """Dynamic ODE solution must converge to the analytical steady state over time."""
        ode_res = simulate_cleavage_ode_dynamics(
            time_span_sec=7200.0,
            time_step_sec=5.0,
            initial_fractions=(1.0, 0.0, 0.0),
            ne_conc_nm=DEFAULT_NE_CF_NM,
            prss8_conc_nm=DEFAULT_PRSS8_NM,
            furin_conc_nm=DEFAULT_FURIN_NM,
        )
        analytical = calculate_cleavage_kinetics(
            ne_conc_nm=DEFAULT_NE_CF_NM,
            prss8_conc_nm=DEFAULT_PRSS8_NM,
            furin_conc_nm=DEFAULT_FURIN_NM,
        )

        final_e0 = ode_res["e0"][-1]
        final_e_alpha = ode_res["e_alpha"][-1]
        final_e_alphagamma = ode_res["e_alphagamma"][-1]
        final_po = ode_res["mean_po"][-1]

        # Check convergence within 0.02 absolute difference
        self.assertAlmostEqual(final_e0, analytical.f_0, delta=0.02)
        self.assertAlmostEqual(final_e_alpha, analytical.f_alpha, delta=0.02)
        self.assertAlmostEqual(final_e_alphagamma, analytical.f_alphagamma, delta=0.02)
        self.assertAlmostEqual(final_po, analytical.mean_po, delta=0.02)

    def test_05_electrophysiology_and_knockdown_scaling(self):
        """Single-channel current and macroscopic current must obey Ohm's law and knockdown scaling."""
        cf_cleavage = calculate_cleavage_kinetics(ne_conc_nm=DEFAULT_NE_CF_NM)
        electro_cf = calculate_electrophysiology(
            cleavage=cf_cleavage,
            enac_knockdown_fraction=0.0,
            pore_block_fraction=0.0,
        )

        # Expected single channel current: gamma_Na (4.5 pS) * |-50 - 60| mV = 4.5 * 110 * 1e-3 pA = 0.495 pA
        self.assertAlmostEqual(electro_cf.i_na_single_pa, 0.495, places=3)
        self.assertGreater(electro_cf.i_na_density_ua_cm2, 30.0)

        # 70% knockdown (ASO)
        electro_aso = calculate_electrophysiology(
            cleavage=cf_cleavage,
            enac_knockdown_fraction=0.70,
            pore_block_fraction=0.0,
        )
        self.assertAlmostEqual(
            electro_aso.i_na_density_ua_cm2,
            electro_cf.i_na_density_ua_cm2 * 0.30,
            delta=0.01,
        )

        # 65% pore block (BI 1265162)
        electro_block = calculate_electrophysiology(
            cleavage=cf_cleavage,
            enac_knockdown_fraction=0.0,
            pore_block_fraction=0.65,
        )
        self.assertAlmostEqual(
            electro_block.i_na_density_ua_cm2,
            electro_cf.i_na_density_ua_cm2 * 0.35,
            delta=0.01,
        )

    def test_06_asl_homeostasis_and_rescue_thresholds(self):
        """ASL height must reflect physiological hydration boundaries."""
        scenarios = get_standard_scenarios()

        wt = scenarios["WT_HEALTHY"]
        self.assertGreaterEqual(wt.asl.h_asl_um, PHYSIOLOGICAL_ASL_THRESHOLD_UM)
        self.assertTrue(wt.asl.is_rescued)
        self.assertEqual(wt.asl.ciliary_beat_competence, "UNHINDERED_CILIARY_CLEARANCE")

        cf = scenarios["CF_UNTREATED"]
        self.assertLess(cf.asl.h_asl_um, 3.5)
        self.assertFalse(cf.asl.is_rescued)
        self.assertEqual(cf.asl.ciliary_beat_competence, "COLLAPSED_CILIARY_ADHESION")

        # Inhaled siRNA (ARO-ENaC / 85% knockdown)
        aro = scenarios["CF_ENAC_SIRNA_ARO"]
        self.assertGreater(aro.asl.h_asl_um, 4.5)
        self.assertGreater(aro.asl.rehydration_percentage, 40.0)

        # Dual-target rescue (70% knockdown + 25% CFTR rescue)
        dual = scenarios["CF_DUAL_TARGET_RESCUE"]
        self.assertGreaterEqual(dual.asl.h_asl_um, PHYSIOLOGICAL_ASL_THRESHOLD_UM)
        self.assertTrue(dual.asl.is_rescued)
        self.assertEqual(dual.asl.ciliary_beat_competence, "UNHINDERED_CILIARY_CLEARANCE")

    def test_07_dual_target_synergy_matrix(self):
        """Dual-target parameter sweep matrix must produce valid monotonic response surfaces."""
        matrix = evaluate_dual_target_synergy_matrix(
            knockdown_steps=(0.0, 0.50, 0.85),
            cftr_rescue_steps=(0.0, 0.20, 0.40),
            ne_conc_nm=DEFAULT_NE_CF_NM,
        )
        self.assertEqual(len(matrix), 9)

        # Untreated corner (0, 0)
        p_0_0 = next(p for p in matrix if p["enac_knockdown_fraction"] == 0.0 and p["cftr_rescue_fraction"] == 0.0)
        self.assertFalse(p_0_0["is_rescued"])
        self.assertLess(p_0_0["h_asl_um"], 3.5)

        # Maximal dual-target corner (0.85, 0.40)
        p_max = next(p for p in matrix if p["enac_knockdown_fraction"] == 0.85 and p["cftr_rescue_fraction"] == 0.40)
        self.assertTrue(p_max["is_rescued"])
        self.assertGreaterEqual(p_max["h_asl_um"], PHYSIOLOGICAL_ASL_THRESHOLD_UM)

    def test_08_cli_execution_and_disclaimers(self):
        """CLI execution must return exit code 0 and include all mandatory disclaimers."""
        script = REPO_ROOT / "tools" / "discovery" / "enac_proteolysis_model.py"
        self.assertTrue(script.exists())

        # Test --json run
        cmd = [sys.executable, str(script), "--json"]
        proc = subprocess.run(cmd, capture_output=True, text=True, check=True)
        self.assertEqual(proc.returncode, 0)

        data = json.loads(proc.stdout)
        self.assertIn("WT_HEALTHY", data)
        self.assertIn("CF_UNTREATED", data)
        self.assertIn("CF_DUAL_TARGET_RESCUE", data)

        cf_entry = data["CF_UNTREATED"]
        self.assertIn("nlm_disclaimer", cf_entry)
        self.assertIn("research_boundary", cf_entry)
        self.assertIn("ClinicalTrials.gov", cf_entry["nlm_disclaimer"])
        self.assertIn("not medical advice", cf_entry["research_boundary"])

        # Test --matrix flag
        cmd_matrix = [sys.executable, str(script), "--matrix", "--json"]
        proc_matrix = subprocess.run(cmd_matrix, capture_output=True, text=True, check=True)
        self.assertEqual(proc_matrix.returncode, 0)
        matrix_data = json.loads(proc_matrix.stdout)
        self.assertIn("matrix", matrix_data)
        self.assertGreater(len(matrix_data["matrix"]), 0)


if __name__ == "__main__":
    unittest.main()
