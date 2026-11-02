"""Unit tests for Round 39 (T-0583): Non-Equilibrium Steady-State Transition Bounds & Sasa-Tasaki Equality in CFTR Catalytic Assemblies in cf-research.

Primary Biophysical Grounding & Evidence:
  - DOI:10.1103/PhysRevLett.95.130602 (Sasa & Tasaki: Equality Satisfied in Steps of a Non-equilibrium Steady State, Phys. Rev. Lett. 95, 130602, 2005)
  - DOI:10.1085/jgp.200609558 (Csanády et al.: Strict coupling between CFTR channel opening and ATP hydrolysis)
  - DOI:10.1126/science.ade2216 (Fiedorczuk & Chen: Cryo-EM structures of rescued human CFTR, PDB 8EIQ)
  - PDB 6MSM (Open phosphorylated human CFTR)
  - PDB 5UAK (Dephosphorylated human CFTR)
  - PDB 8EIQ (Rescued human CFTR bound to ETI)
  - Lean 4 Formal Math Proof Benchmarks: fermats_last_theorem, formal_math_zeta23, formal_math_3sum_apsp
"""

from __future__ import annotations

import unittest
from tools.discovery.formal_proof_thermodynamic_landscape import (
    evaluate_sasa_tasaki_equality,
)


class SasaTasakiEqualityTests(unittest.TestCase):
    """Test suite verifying Round 39 Sasa-Tasaki non-equilibrium transition bounds in CFTR catalytic assemblies."""

    def setUp(self) -> None:
        self.report = evaluate_sasa_tasaki_equality()

    def test_task_id_and_overall_bridge_verification(self) -> None:
        self.assertEqual(self.report["task_id"], "T-0583")
        self.assertTrue(self.report["sasa_tasaki_bridge_verified"])
        self.assertTrue(self.report["ness_transition_stratification_verified"])
        self.assertTrue(self.report["proof_dag_sasa_tasaki_parity_verified"])
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_sasa_tasaki_equality_and_ness_stratification(self) -> None:
        profiles = self.report["gating_profiles"]
        self.assertIn("wt_cftr", profiles)
        self.assertIn("f508del_uncorrected", profiles)
        self.assertIn("g551d_uncoupled_cftr", profiles)
        self.assertIn("rescued_cftr_ivacaftor", profiles)
        self.assertIn("w1282x_nonsense_baseline", profiles)

        wt = profiles["wt_cftr"]
        f508del = profiles["f508del_uncorrected"]
        g551d = profiles["g551d_uncoupled_cftr"]
        rescued = profiles["rescued_cftr_ivacaftor"]
        w1282x = profiles["w1282x_nonsense_baseline"]

        # WT CFTR: Sasa-Tasaki equality satisfied (<exp(-beta[W_ex - Delta F])> ~ 1.0, W_ex <= 3.5 kcal/mol, eta >= 70%, rate >= 1.50 s^-1)
        self.assertTrue(wt["sasa_tasaki_bound_satisfied"])
        self.assertEqual(wt["gating_regime"], "EFFICIENT_STEADY_STATE_TRANSITION")
        self.assertLessEqual(wt["excess_work_w_ex_kcal_mol"], 3.50)
        self.assertGreaterEqual(wt["excess_work_w_ex_kcal_mol"], wt["delta_f_ness_kcal_mol"])
        self.assertAlmostEqual(wt["sasa_tasaki_exp_ratio"], 1.00, delta=0.05)
        self.assertGreaterEqual(wt["active_power_stroke_efficiency_pct"], 70.0)
        self.assertGreaterEqual(wt["opening_transition_rate_k_s_inv"], 1.50)
        self.assertAlmostEqual(wt["excess_dissipation_kcal_mol"], wt["excess_work_w_ex_kcal_mol"] - wt["delta_f_ness_kcal_mol"], places=2)

        # F508del: Degenerate work dissipation trap (W_ex >= 13.0 kcal/mol, eta < 10%, rate < 0.01 s^-1, <exp> << 1.0)
        self.assertFalse(f508del["sasa_tasaki_bound_satisfied"])
        self.assertEqual(f508del["gating_regime"], "DEGENERATE_WORK_DISSIPATION_TRAP")
        self.assertGreaterEqual(f508del["excess_work_w_ex_kcal_mol"], 13.0)
        self.assertLess(f508del["active_power_stroke_efficiency_pct"], 10.0)
        self.assertLess(f508del["opening_transition_rate_k_s_inv"], 0.01)
        self.assertLess(f508del["sasa_tasaki_exp_ratio"], 0.001)

        # G551D: Degenerate work dissipation trap (W_ex >= 13.0 kcal/mol, eta < 10%, rate < 0.01 s^-1, <exp> << 1.0)
        self.assertFalse(g551d["sasa_tasaki_bound_satisfied"])
        self.assertEqual(g551d["gating_regime"], "DEGENERATE_WORK_DISSIPATION_TRAP")
        self.assertGreaterEqual(g551d["excess_work_w_ex_kcal_mol"], 13.0)
        self.assertLess(g551d["active_power_stroke_efficiency_pct"], 10.0)
        self.assertLess(g551d["opening_transition_rate_k_s_inv"], 0.01)
        self.assertLess(g551d["sasa_tasaki_exp_ratio"], 0.001)

        # Rescued CFTR (Ivacaftor / ETI): Restores Sasa-Tasaki transition bounds matching WT kinetics
        self.assertTrue(rescued["sasa_tasaki_bound_satisfied"])
        self.assertEqual(rescued["gating_regime"], "EFFICIENT_STEADY_STATE_TRANSITION")
        self.assertLessEqual(rescued["excess_work_w_ex_kcal_mol"], 3.50)
        self.assertGreaterEqual(rescued["excess_work_w_ex_kcal_mol"], rescued["delta_f_ness_kcal_mol"])
        self.assertAlmostEqual(rescued["sasa_tasaki_exp_ratio"], 1.00, delta=0.05)
        self.assertGreaterEqual(rescued["active_power_stroke_efficiency_pct"], 60.0)
        self.assertGreaterEqual(rescued["opening_transition_rate_k_s_inv"], 1.50)

        # W1282X: Passive equilibrium baseline
        self.assertFalse(w1282x["sasa_tasaki_bound_satisfied"])
        self.assertEqual(w1282x["gating_regime"], "PASSIVE_EQUILIBRIUM_APO")
        self.assertEqual(w1282x["excess_work_w_ex_kcal_mol"], 0.0)
        self.assertEqual(w1282x["delta_f_ness_kcal_mol"], 0.0)
        self.assertEqual(w1282x["sasa_tasaki_exp_ratio"], 0.0)
        self.assertEqual(w1282x["opening_transition_rate_k_s_inv"], 0.0)

    def test_formal_proof_dag_sasa_tasaki_parity(self) -> None:
        challenges = self.report["challenge_profiles"]
        expected_benchmarks = [
            "anthropics/fermats-last-theorem",
            "anthropics/formal-math/zeta23",
            "anthropics/formal-math/3sum-apsp",
        ]
        for c_id in expected_benchmarks:
            self.assertIn(c_id, challenges)
            c_prof = challenges[c_id]
            self.assertTrue(c_prof["sasa_tasaki_proof_parity_satisfied"])
            self.assertTrue(c_prof["excess_work_bound_satisfied"])
            self.assertGreaterEqual(c_prof["excess_derivation_work_dag"], c_prof["delta_f_proof_dag"])
            self.assertAlmostEqual(c_prof["sasa_tasaki_equality_exp_ratio"], 1.00, delta=0.05)
            self.assertGreaterEqual(c_prof["inferential_work_reduction_pct"], 50.0)

    def test_ah10_negative_controls(self) -> None:
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_invalid_regimes_count_raises_value_error(self) -> None:
        with self.assertRaises(ValueError):
            evaluate_sasa_tasaki_equality(n_gating_regimes=2)


if __name__ == "__main__":
    unittest.main()
