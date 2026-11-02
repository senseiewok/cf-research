"""Unit tests for Round 26 (T-0553): Glansdorff-Prigogine Universal Evolution Criterion & Non-Linear Gating Dynamics in CFTR.

Primary Biophysical Grounding & Evidence:
  - DOI:10.1016/0031-8914(64)90154-1 (Glansdorff & Prigogine: General evolution criterion in macroscopic physics)
  - DOI:10.1085/jgp.200609558 (Thermodynamics of CFTR channel gating)
  - DOI:10.1126/science.ade2216 (Cryo-EM structures of rescued CFTR)
  - PDB 6MSM (Open phosphorylated human CFTR)
  - PDB 5UAK (Dephosphorylated human CFTR)
  - Lean 4 Formal Math Proof Benchmarks: fermats_last_theorem, formal_math_zeta23, formal_math_3sum_apsp
"""

from __future__ import annotations

import unittest
from tools.discovery.formal_proof_thermodynamic_landscape import (
    evaluate_glansdorff_prigogine_evolution,
)


class GlansdorffPrigogineEvolutionTests(unittest.TestCase):
    """Test suite verifying Round 26 Glansdorff-Prigogine universal evolution criterion & non-linear gating dynamics in CFTR."""

    def setUp(self) -> None:
        self.report = evaluate_glansdorff_prigogine_evolution()

    def test_task_id_and_overall_bridge_verification(self) -> None:
        self.assertEqual(self.report["task_id"], "T-0553")
        self.assertTrue(self.report["glansdorff_prigogine_bridge_verified"])
        self.assertTrue(self.report["gating_stability_stratification_verified"])
        self.assertTrue(self.report["proof_universal_evolution_parity_verified"])
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_glansdorff_prigogine_evolution_and_dynamic_stability_stratification(self) -> None:
        profiles = self.report["gating_profiles"]
        self.assertIn("wt_cftr", profiles)
        self.assertIn("g551d_uncoupled_cftr", profiles)
        self.assertIn("rescued_cftr_ivacaftor", profiles)
        self.assertIn("w1282x_nonsense_baseline", profiles)

        wt = profiles["wt_cftr"]
        g551d = profiles["g551d_uncoupled_cftr"]
        rescued = profiles["rescued_cftr_ivacaftor"]
        w1282x = profiles["w1282x_nonsense_baseline"]

        # WT CFTR: Satisfies Glansdorff-Prigogine evolution with negative force variation rate and positive excess entropy production
        self.assertTrue(wt["glansdorff_criterion_satisfied"])
        self.assertEqual(wt["stationary_stability_status"], "STABLE_NONLINEAR_CYCLE")
        self.assertLessEqual(wt["d_x_sigma_rate"], 0.0)
        self.assertGreaterEqual(wt["excess_entropy_production_rate"], 0.0)
        self.assertLessEqual(wt["entropy_production_rate"], 5.0)

        # G551D: Uncoupled kinetic trap failing Glansdorff-Prigogine universal evolution
        self.assertFalse(g551d["glansdorff_criterion_satisfied"])
        self.assertEqual(g551d["stationary_stability_status"], "UNSTABLE_KINETIC_TRAP")
        self.assertGreater(g551d["d_x_sigma_rate"], 0.0)
        self.assertLess(g551d["excess_entropy_production_rate"], 0.0)
        self.assertGreaterEqual(g551d["entropy_production_rate"], 10.0)

        # Rescued Ivacaftor: Restores non-linear stability and non-positive force evolution rate
        self.assertTrue(rescued["glansdorff_criterion_satisfied"])
        self.assertEqual(rescued["stationary_stability_status"], "STABLE_NONLINEAR_CYCLE")
        self.assertLessEqual(rescued["d_x_sigma_rate"], 0.0)
        self.assertGreaterEqual(rescued["excess_entropy_production_rate"], 0.0)
        self.assertLessEqual(rescued["entropy_production_rate"], 5.0)

        # W1282X: Passive equilibrium baseline
        self.assertEqual(w1282x["stationary_stability_status"], "PASSIVE_EQUILIBRIUM_APO")
        self.assertEqual(w1282x["d_x_sigma_rate"], 0.0)
        self.assertEqual(w1282x["excess_entropy_production_rate"], 0.0)
        self.assertEqual(w1282x["entropy_production_rate"], 0.0)

    def test_formal_proof_universal_evolution_deductive_invariance(self) -> None:
        challenges = self.report["challenge_profiles"]
        expected_benchmarks = [
            "anthropics/fermats-last-theorem",
            "anthropics/formal-math/zeta23",
            "anthropics/formal-math/3sum-apsp",
        ]
        for c_id in expected_benchmarks:
            self.assertIn(c_id, challenges)
            c_prof = challenges[c_id]
            self.assertTrue(c_prof["universal_evolution_satisfied"])
            self.assertLessEqual(c_prof["inferential_d_x_sigma_rate"], 0.0)
            self.assertGreaterEqual(c_prof["inferential_excess_entropy_rate"], 0.0)
            self.assertGreaterEqual(c_prof["stability_margin"], 3.0)

    def test_ah10_negative_controls(self) -> None:
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_invalid_regimes_count_raises_value_error(self) -> None:
        with self.assertRaises(ValueError):
            evaluate_glansdorff_prigogine_evolution(n_gating_regimes=2)


if __name__ == "__main__":
    unittest.main()
