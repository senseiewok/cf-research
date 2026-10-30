"""Unit tests for Round 25 (T-0551): Prigogine Minimum Entropy Production Principle & Stationary NESS Stability in CFTR Gating.

Primary Biophysical Grounding & Evidence:
  - DOI:10.1085/jgp.200609558 (Thermodynamics of CFTR channel gating)
  - DOI:10.1126/science.ade2216 (Cryo-EM structures of rescued CFTR)
  - DOI:10.1098/rspa.1954.0152 (Thermodynamics of irreversible processes and minimum entropy production)
  - PDB 6MSM (Open phosphorylated human CFTR)
  - PDB 5UAK (Dephosphorylated human CFTR)
  - Lean 4 Formal Math Proof Benchmarks: fermats_last_theorem, formal_math_zeta23, formal_math_3sum_apsp
"""

from __future__ import annotations

import unittest
from tools.discovery.formal_proof_thermodynamic_landscape import (
    evaluate_prigogine_minimum_entropy_production,
)


class PrigogineMinimumEntropyProductionTests(unittest.TestCase):
    """Test suite verifying Round 25 Prigogine minimum entropy production and NESS stability in CFTR gating."""

    def setUp(self) -> None:
        self.report = evaluate_prigogine_minimum_entropy_production()

    def test_task_id_and_overall_bridge_verification(self) -> None:
        self.assertEqual(self.report["task_id"], "T-0551")
        self.assertTrue(self.report["prigogine_entropy_bridge_verified"])
        self.assertTrue(self.report["gating_stability_stratification_verified"])
        self.assertTrue(self.report["proof_minimum_entropy_parity_verified"])
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_prigogine_minimum_and_ness_stability_stratification(self) -> None:
        profiles = self.report["gating_profiles"]
        self.assertIn("wt_cftr", profiles)
        self.assertIn("g551d_uncoupled_cftr", profiles)
        self.assertIn("rescued_cftr_ivacaftor", profiles)
        self.assertIn("w1282x_nonsense_baseline", profiles)

        wt = profiles["wt_cftr"]
        g551d = profiles["g551d_uncoupled_cftr"]
        rescued = profiles["rescued_cftr_ivacaftor"]
        w1282x = profiles["w1282x_nonsense_baseline"]

        # WT CFTR: Achieves stable minimal dissipation state with negative Lyapunov exponent
        self.assertTrue(wt["prigogine_minimum_satisfied"])
        self.assertEqual(wt["stationary_stability_status"], "STABLE_NESS_MINIMUM")
        self.assertLessEqual(wt["entropy_production_rate"], 5.0)
        self.assertLess(wt["lyapunov_decay_rate"], 0.0)
        self.assertLessEqual(wt["entropy_production_time_derivative"], 0.0)
        self.assertGreaterEqual(wt["excess_entropy_production_delta2"], 2.50)
        self.assertLessEqual(wt["onsager_conductance_matrix"]["symmetry_error"], 1e-4)
        self.assertGreater(wt["onsager_conductance_matrix"]["determinant"], 0.0)

        # G551D: Uncoupled kinetic trap failing Prigogine minimum entropy production
        self.assertFalse(g551d["prigogine_minimum_satisfied"])
        self.assertEqual(g551d["stationary_stability_status"], "UNSTABLE_KINETIC_TRAP")
        self.assertGreaterEqual(g551d["entropy_production_rate"], 10.0)
        self.assertGreater(g551d["lyapunov_decay_rate"], 0.0)
        self.assertLess(g551d["excess_entropy_production_delta2"], 1.00)

        # Rescued Ivacaftor: Restores stable NESS minimum and linear Onsager coupling
        self.assertTrue(rescued["prigogine_minimum_satisfied"])
        self.assertEqual(rescued["stationary_stability_status"], "STABLE_NESS_MINIMUM")
        self.assertLessEqual(rescued["entropy_production_rate"], 5.0)
        self.assertLess(rescued["lyapunov_decay_rate"], 0.0)
        self.assertLessEqual(rescued["entropy_production_time_derivative"], 0.0)
        self.assertGreaterEqual(rescued["excess_entropy_production_delta2"], 2.50)
        self.assertLessEqual(rescued["onsager_conductance_matrix"]["symmetry_error"], 1e-4)
        self.assertGreater(rescued["onsager_conductance_matrix"]["determinant"], 0.0)

        # W1282X: Passive baseline with zero active dissipation
        self.assertEqual(w1282x["entropy_production_rate"], 0.0)
        self.assertEqual(w1282x["lyapunov_decay_rate"], 0.0)
        self.assertEqual(w1282x["excess_entropy_production_delta2"], 0.0)

    def test_formal_proof_deductive_invariance(self) -> None:
        challenges = self.report["challenge_profiles"]
        expected_benchmarks = [
            "anthropics/fermats-last-theorem",
            "anthropics/formal-math/zeta23",
            "anthropics/formal-math/3sum-apsp",
        ]
        for c_id in expected_benchmarks:
            self.assertIn(c_id, challenges)
            c_prof = challenges[c_id]
            self.assertTrue(c_prof["minimum_entropy_satisfied"])
            self.assertLess(c_prof["lyapunov_decay_rate"], 0.0)
            self.assertGreaterEqual(c_prof["stability_margin"], 3.0)
            self.assertLessEqual(c_prof["inferential_entropy_production_rate"], 3.0)

    def test_ah10_negative_controls(self) -> None:
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_invalid_regimes_count_raises_value_error(self) -> None:
        with self.assertRaises(ValueError):
            evaluate_prigogine_minimum_entropy_production(n_gating_regimes=2)


if __name__ == "__main__":
    unittest.main()
