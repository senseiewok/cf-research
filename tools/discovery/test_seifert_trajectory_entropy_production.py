"""Unit tests for Round 40 (T-0585): Seifert Trajectory Stochastic Entropy Production & Fluctuation Theorem in CFTR Catalytic Cycles in cf-research.

Primary Biophysical Grounding & Evidence:
  - DOI:10.1088/0034-4885/75/12/126001 (Seifert U: Stochastic thermodynamics, fluctuation theorems and molecular machines, Rep. Prog. Phys. 75, 126001, 2012)
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
    evaluate_seifert_trajectory_entropy_production,
)


class SeifertTrajectoryEntropyProductionTests(unittest.TestCase):
    """Test suite verifying Round 40 Seifert trajectory stochastic entropy production & fluctuation theorem in CFTR."""

    def setUp(self) -> None:
        self.report = evaluate_seifert_trajectory_entropy_production()

    def test_task_id_and_overall_bridge_verification(self) -> None:
        self.assertEqual(self.report["task_id"], "T-0585")
        self.assertTrue(self.report["seifert_trajectory_bridge_verified"])
        self.assertTrue(self.report["trajectory_entropy_stratification_verified"])
        self.assertTrue(self.report["proof_dag_seifert_parity_verified"])
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_seifert_trajectory_entropy_and_gating_stratification(self) -> None:
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

        # WT CFTR: Active directional cycling (<exp(-Delta s_tot/k_B)> ~ 1.0, <Delta s_tot> >= 4.5 k_B, P/P_rev >= 100, rate >= 1.50 s^-1)
        self.assertTrue(wt["trajectory_fluctuation_bound_satisfied"])
        self.assertEqual(wt["gating_regime"], "ACTIVE_IRREVERSIBLE_CATALYTIC_CYCLE")
        self.assertGreaterEqual(wt["mean_trajectory_entropy_production_kb"], 4.50)
        self.assertAlmostEqual(wt["integral_fluctuation_exp_ratio"], 1.00, delta=0.05)
        self.assertGreaterEqual(wt["forward_backward_path_prob_ratio"], 100.0)
        self.assertGreaterEqual(wt["active_cycle_fraction_pct"], 75.0)
        self.assertGreaterEqual(wt["opening_transition_rate_k_s_inv"], 1.50)

        # F508del: Degenerate futile dissipation trap (<Delta s_tot> < 0.2 k_B, active cycle < 5%, rate < 0.01 s^-1, IFT fails)
        self.assertFalse(f508del["trajectory_fluctuation_bound_satisfied"])
        self.assertEqual(f508del["gating_regime"], "DEGENERATE_FUTILE_DISSIPATION_TRAP")
        self.assertLess(f508del["mean_trajectory_entropy_production_kb"], 0.20)
        self.assertLess(f508del["active_cycle_fraction_pct"], 5.0)
        self.assertLess(f508del["opening_transition_rate_k_s_inv"], 0.01)
        self.assertLess(f508del["integral_fluctuation_exp_ratio"], 0.01)

        # G551D: Degenerate futile dissipation trap (<Delta s_tot> < 0.2 k_B, active cycle < 5%, rate < 0.01 s^-1, IFT fails)
        self.assertFalse(g551d["trajectory_fluctuation_bound_satisfied"])
        self.assertEqual(g551d["gating_regime"], "DEGENERATE_FUTILE_DISSIPATION_TRAP")
        self.assertLess(g551d["mean_trajectory_entropy_production_kb"], 0.20)
        self.assertLess(g551d["active_cycle_fraction_pct"], 5.0)
        self.assertLess(g551d["opening_transition_rate_k_s_inv"], 0.01)
        self.assertLess(g551d["integral_fluctuation_exp_ratio"], 0.01)

        # Rescued CFTR (Ivacaftor / ETI): Restores trajectory entropy production matching WT kinetics
        self.assertTrue(rescued["trajectory_fluctuation_bound_satisfied"])
        self.assertEqual(rescued["gating_regime"], "ACTIVE_IRREVERSIBLE_CATALYTIC_CYCLE")
        self.assertGreaterEqual(rescued["mean_trajectory_entropy_production_kb"], 4.00)
        self.assertAlmostEqual(rescued["integral_fluctuation_exp_ratio"], 1.00, delta=0.05)
        self.assertGreaterEqual(rescued["forward_backward_path_prob_ratio"], 50.0)
        self.assertGreaterEqual(rescued["active_cycle_fraction_pct"], 65.0)
        self.assertGreaterEqual(rescued["opening_transition_rate_k_s_inv"], 1.50)

        # W1282X: Passive equilibrium APO baseline
        self.assertFalse(w1282x["trajectory_fluctuation_bound_satisfied"])
        self.assertEqual(w1282x["gating_regime"], "PASSIVE_EQUILIBRIUM_APO")
        self.assertEqual(w1282x["mean_trajectory_entropy_production_kb"], 0.0)
        self.assertEqual(w1282x["integral_fluctuation_exp_ratio"], 0.0)
        self.assertEqual(w1282x["forward_backward_path_prob_ratio"], 1.0)
        self.assertEqual(w1282x["active_cycle_fraction_pct"], 0.0)
        self.assertEqual(w1282x["opening_transition_rate_k_s_inv"], 0.0)

    def test_formal_proof_dag_seifert_parity(self) -> None:
        challenges = self.report["challenge_profiles"]
        expected_benchmarks = [
            "anthropics/fermats-last-theorem",
            "anthropics/formal-math/zeta23",
            "anthropics/formal-math/3sum-apsp",
        ]
        for c_id in expected_benchmarks:
            self.assertIn(c_id, challenges)
            c_prof = challenges[c_id]
            self.assertTrue(c_prof["seifert_proof_parity_satisfied"])
            self.assertTrue(c_prof["trajectory_entropy_bound_satisfied"])
            self.assertGreaterEqual(c_prof["proof_trajectory_entropy_reduction_dag"], 4.0)
            self.assertAlmostEqual(c_prof["integral_fluctuation_exp_ratio"], 1.00, delta=0.05)
            self.assertGreaterEqual(c_prof["inferential_path_asymmetry_ratio"], 100.0)
            self.assertGreaterEqual(c_prof["inferential_irreversibility_pct"], 70.0)

    def test_ah10_negative_controls(self) -> None:
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_invalid_regimes_count_raises_value_error(self) -> None:
        with self.assertRaises(ValueError):
            evaluate_seifert_trajectory_entropy_production(n_gating_regimes=2)


if __name__ == "__main__":
    unittest.main()
