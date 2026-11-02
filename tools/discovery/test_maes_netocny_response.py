"""Unit tests for Round 33 (T-0577): Maes-Netočný Non-Equilibrium Response Formula & Excess Entropy Production for Rescued CFTR Kinetics.

Primary Biophysical Grounding & Evidence:
  - DOI:10.1063/1.2738753 (Maes & Netočný: Minimum entropy production principle from a dynamical fluctuation law, J. Math. Phys. 48, 053306, 2007)
  - DOI:10.1085/jgp.200609558 (Csanády et al.: Strict coupling between CFTR channel opening and ATP hydrolysis)
  - DOI:10.1126/science.ade2216 (Fiedorczuk & Chen: Cryo-EM structures of rescued human CFTR, PDB 8EIQ)
  - PDB 6MSM (Open phosphorylated human CFTR)
  - PDB 5UAK (Dephosphorylated human CFTR)
  - Lean 4 Formal Math Proof Benchmarks: fermats_last_theorem, formal_math_zeta23, formal_math_3sum_apsp
"""

from __future__ import annotations

import unittest
from tools.discovery.formal_proof_thermodynamic_landscape import (
    evaluate_maes_netocny_response,
)


class MaesNetocnyResponseTests(unittest.TestCase):
    """Test suite verifying Round 33 Maes-Netočný non-equilibrium response formula & excess entropy production in CFTR kinetics."""

    def setUp(self) -> None:
        self.report = evaluate_maes_netocny_response()

    def test_task_id_and_overall_bridge_verification(self) -> None:
        self.assertEqual(self.report["task_id"], "T-0577")
        self.assertTrue(self.report["maes_netocny_bridge_verified"])
        self.assertTrue(self.report["frenetic_response_stratification_verified"])
        self.assertTrue(self.report["proof_dag_frenetic_acceleration_verified"])
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_maes_netocny_response_and_frenetic_stratification(self) -> None:
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

        # WT CFTR: Balanced entropic-frenetic response (R >= 1.0, K >= 1.20, excess S >= 2.0, k >= 1.50 s^-1)
        self.assertTrue(wt["maes_netocny_criterion_satisfied"])
        self.assertEqual(wt["gating_regime"], "BALANCED_ENTROPIC_FRENETIC_RESPONSE")
        self.assertGreaterEqual(wt["net_response_coefficient_r"], 1.0)
        self.assertGreaterEqual(wt["frenetic_activity_k"], 1.20)
        self.assertGreaterEqual(wt["excess_entropy_production_s"], 2.0)
        self.assertGreaterEqual(wt["opening_transition_rate_k_s_inv"], 1.50)
        self.assertGreater(wt["frenetic_acceleration_ratio"], 1.0)

        # F508del: Frenetic freeze kinetic trap (R < 0.05, K < 0.10, k < 0.01 s^-1)
        self.assertFalse(f508del["maes_netocny_criterion_satisfied"])
        self.assertEqual(f508del["gating_regime"], "FRENETIC_FREEZE_KINETIC_TRAP")
        self.assertLess(f508del["net_response_coefficient_r"], 0.05)
        self.assertLess(f508del["frenetic_activity_k"], 0.10)
        self.assertLess(f508del["opening_transition_rate_k_s_inv"], 0.01)

        # G551D: Frenetic freeze kinetic trap (R < 0.05, K < 0.10, k < 0.01 s^-1)
        self.assertFalse(g551d["maes_netocny_criterion_satisfied"])
        self.assertEqual(g551d["gating_regime"], "FRENETIC_FREEZE_KINETIC_TRAP")
        self.assertLess(g551d["net_response_coefficient_r"], 0.05)
        self.assertLess(g551d["frenetic_activity_k"], 0.10)
        self.assertLess(g551d["opening_transition_rate_k_s_inv"], 0.01)

        # Rescued CFTR (Ivacaftor / ETI): Restores balanced entropic-frenetic response matching WT bounds
        self.assertTrue(rescued["maes_netocny_criterion_satisfied"])
        self.assertEqual(rescued["gating_regime"], "BALANCED_ENTROPIC_FRENETIC_RESPONSE")
        self.assertGreaterEqual(rescued["net_response_coefficient_r"], 1.0)
        self.assertGreaterEqual(rescued["frenetic_activity_k"], 1.20)
        self.assertGreaterEqual(rescued["opening_transition_rate_k_s_inv"], 1.50)

        # W1282X: Passive equilibrium baseline
        self.assertFalse(w1282x["maes_netocny_criterion_satisfied"])
        self.assertEqual(w1282x["gating_regime"], "PASSIVE_EQUILIBRIUM_APO")
        self.assertEqual(w1282x["net_response_coefficient_r"], 0.0)
        self.assertEqual(w1282x["frenetic_activity_k"], 0.0)
        self.assertEqual(w1282x["opening_transition_rate_k_s_inv"], 0.0)

    def test_formal_proof_dag_frenetic_acceleration_parity(self) -> None:
        challenges = self.report["challenge_profiles"]
        expected_benchmarks = [
            "anthropics/fermats-last-theorem",
            "anthropics/formal-math/zeta23",
            "anthropics/formal-math/3sum-apsp",
        ]
        for c_id in expected_benchmarks:
            self.assertIn(c_id, challenges)
            c_prof = challenges[c_id]
            self.assertTrue(c_prof["maes_netocny_proof_parity_satisfied"])
            self.assertTrue(c_prof["frenetic_bound_satisfied"])
            self.assertGreaterEqual(c_prof["net_response_dag_r"], 1.0)
            self.assertGreaterEqual(c_prof["frenetic_throughput_k"], 2.0)
            self.assertGreaterEqual(c_prof["inferential_traffic_acceleration_pct"], 50.0)

    def test_ah10_negative_controls(self) -> None:
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_invalid_regimes_count_raises_value_error(self) -> None:
        with self.assertRaises(ValueError):
            evaluate_maes_netocny_response(n_gating_regimes=2)


if __name__ == "__main__":
    unittest.main()
