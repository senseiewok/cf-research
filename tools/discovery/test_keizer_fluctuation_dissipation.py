"""Unit tests for Round 27 (T-0559): Keizer's Fluctuation-Dissipation Theory for Non-Linear Reaction-Diffusion in CFTR Gating & Airway ASL Transport.

Primary Biophysical Grounding & Evidence:
  - DOI:10.1063/1.432650 (Keizer: A theory of non-equilibrium thermodynamics for chemical reactions and diffusion)
  - DOI:10.1085/jgp.200609558 (Thermodynamics of CFTR channel gating)
  - DOI:10.1126/science.ade2216 (Cryo-EM structures of rescued CFTR)
  - PDB 6MSM (Open phosphorylated human CFTR)
  - PDB 5UAK (Dephosphorylated human CFTR)
  - Lean 4 Formal Math Proof Benchmarks: fermats_last_theorem, formal_math_zeta23, formal_math_3sum_apsp
"""

from __future__ import annotations

import unittest
from tools.discovery.formal_proof_thermodynamic_landscape import (
    evaluate_keizer_fluctuation_dissipation,
)


class KeizerFluctuationDissipationTests(unittest.TestCase):
    """Test suite verifying Round 27 Keizer fluctuation-dissipation theory and reaction-diffusion coupling in CFTR gating."""

    def setUp(self) -> None:
        self.report = evaluate_keizer_fluctuation_dissipation()

    def test_task_id_and_overall_bridge_verification(self) -> None:
        self.assertEqual(self.report["task_id"], "T-0559")
        self.assertTrue(self.report["keizer_fluctuation_dissipation_bridge_verified"])
        self.assertTrue(self.report["gating_stratification_verified"])
        self.assertTrue(self.report["proof_variance_bounds_parity_verified"])
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_keizer_fluctuation_dissipation_and_genotype_stratification(self) -> None:
        profiles = self.report["gating_profiles"]
        self.assertIn("wt_cftr", profiles)
        self.assertIn("g551d_uncoupled_cftr", profiles)
        self.assertIn("rescued_cftr_ivacaftor", profiles)
        self.assertIn("w1282x_nonsense_baseline", profiles)

        wt = profiles["wt_cftr"]
        g551d = profiles["g551d_uncoupled_cftr"]
        rescued = profiles["rescued_cftr_ivacaftor"]
        w1282x = profiles["w1282x_nonsense_baseline"]

        # WT CFTR: High fluctuation-dissipation coupling efficiency with bounded covariance spread
        self.assertTrue(wt["keizer_criterion_satisfied"])
        self.assertEqual(wt["transport_regime"], "STABLE_COUPLED_CYCLE")
        self.assertGreater(wt["atp_cl_covariance"], 0.0)
        self.assertGreaterEqual(wt["coupling_efficiency"], 0.75)
        self.assertGreater(wt["atp_affinity_kj_mol"], 15.0)
        self.assertLessEqual(wt["dissipation_fluctuation_ratio"], 2.0)

        # G551D: Disconnected cross-correlation covariance representing uncoupled gating dissipation
        self.assertFalse(g551d["keizer_criterion_satisfied"])
        self.assertEqual(g551d["transport_regime"], "UNCOUPLED_DISSIPATIVE_TRAP")
        self.assertLess(g551d["atp_cl_covariance"], 0.5)
        self.assertLess(g551d["coupling_efficiency"], 0.20)
        self.assertGreaterEqual(g551d["dissipation_fluctuation_ratio"], 5.0)

        # Rescued Ivacaftor: Restores non-linear reaction-diffusion covariance matching wild-type bounds
        self.assertTrue(rescued["keizer_criterion_satisfied"])
        self.assertEqual(rescued["transport_regime"], "STABLE_COUPLED_CYCLE")
        self.assertGreater(rescued["atp_cl_covariance"], 0.0)
        self.assertGreaterEqual(rescued["coupling_efficiency"], 0.75)
        self.assertGreater(rescued["atp_affinity_kj_mol"], 15.0)
        self.assertLessEqual(rescued["dissipation_fluctuation_ratio"], 2.0)

        # W1282X: Passive equilibrium baseline
        self.assertEqual(w1282x["transport_regime"], "PASSIVE_EQUILIBRIUM_APO")
        self.assertEqual(w1282x["atp_cl_covariance"], 0.0)
        self.assertEqual(w1282x["coupling_efficiency"], 0.0)
        self.assertFalse(w1282x["keizer_criterion_satisfied"])

    def test_formal_proof_variance_bounds_parity(self) -> None:
        challenges = self.report["challenge_profiles"]
        expected_benchmarks = [
            "anthropics/fermats-last-theorem",
            "anthropics/formal-math/zeta23",
            "anthropics/formal-math/3sum-apsp",
        ]
        for c_id in expected_benchmarks:
            self.assertIn(c_id, challenges)
            c_prof = challenges[c_id]
            self.assertTrue(c_prof["inferential_variance_bound_satisfied"])
            self.assertGreater(c_prof["branching_covariance"], 0.0)
            self.assertGreaterEqual(c_prof["coupling_efficiency"], 0.70)
            self.assertGreater(c_prof["dissipation_margin"], 2.0)

    def test_ah10_negative_controls(self) -> None:
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_invalid_regimes_count_raises_value_error(self) -> None:
        with self.assertRaises(ValueError):
            evaluate_keizer_fluctuation_dissipation(n_reaction_diffusion_regimes=2)


if __name__ == "__main__":
    unittest.main()
