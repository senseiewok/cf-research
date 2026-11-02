"""Unit tests for Round 29 (T-0567): Stochastic Path Integral & Onsager-Machlup Action for Active CFTR Dimerization.

Primary Biophysical Grounding & Evidence:
  - DOI:10.1103/PhysRev.91.1505 (Machlup & Onsager: Fluctuations and irreversible processes. II. Systems with kinetic energy)
  - DOI:10.1085/jgp.200609558 (Strict coupling between CFTR channel opening and ATP hydrolysis)
  - DOI:10.1126/science.ade2216 (Cryo-EM structures of rescued human CFTR, PDB 8EIQ)
  - PDB 6MSM (Open phosphorylated human CFTR)
  - PDB 5UAK (Dephosphorylated human CFTR)
  - Lean 4 Formal Math Proof Benchmarks: fermats_last_theorem, formal_math_zeta23, formal_math_3sum_apsp
"""

from __future__ import annotations

import unittest
from tools.discovery.formal_proof_thermodynamic_landscape import (
    evaluate_onsager_machlup_action,
)


class OnsagerMachlupActionTests(unittest.TestCase):
    """Test suite verifying Round 29 Onsager-Machlup stochastic path integral & action minimization for active CFTR dimerization."""

    def setUp(self) -> None:
        self.report = evaluate_onsager_machlup_action()

    def test_task_id_and_overall_bridge_verification(self) -> None:
        self.assertEqual(self.report["task_id"], "T-0567")
        self.assertTrue(self.report["onsager_machlup_action_bridge_verified"])
        self.assertTrue(self.report["instanton_path_minimization_verified"])
        self.assertTrue(self.report["proof_action_geodesic_parity_verified"])
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_onsager_machlup_action_and_instanton_path_stratification(self) -> None:
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

        # WT CFTR: Low action barrier (S_min <= 6.5 kcal/mol) along the ATP-coupled instanton pathway
        self.assertTrue(wt["onsager_machlup_criterion_satisfied"])
        self.assertEqual(wt["gating_regime"], "ACTIVE_INSTANTON_CROSSING")
        self.assertLessEqual(wt["instanton_action_s_min_kcal_mol"], 6.5)
        self.assertGreaterEqual(wt["instanton_traversal_probability"], 0.85)
        self.assertGreaterEqual(wt["drift_velocity_v_mu_m_s"], 1.0)
        self.assertGreaterEqual(wt["diffusion_tensor_d"], 0.35)

        # F508del: Distorted potential energy landscape increasing transition action (S >= 18.0 kcal/mol)
        self.assertFalse(f508del["onsager_machlup_criterion_satisfied"])
        self.assertEqual(f508del["gating_regime"], "DISTORTED_ENERGY_LANDSCAPE_TRAP")
        self.assertGreaterEqual(f508del["instanton_action_s_min_kcal_mol"], 18.0)
        self.assertLess(f508del["instanton_traversal_probability"], 0.05)

        # G551D: Distorted energy landscape trap preventing instanton crossing (S >= 18.0 kcal/mol)
        self.assertFalse(g551d["onsager_machlup_criterion_satisfied"])
        self.assertEqual(g551d["gating_regime"], "DISTORTED_ENERGY_LANDSCAPE_TRAP")
        self.assertGreaterEqual(g551d["instanton_action_s_min_kcal_mol"], 18.0)
        self.assertLess(g551d["instanton_traversal_probability"], 0.05)

        # Rescued CFTR (Ivacaftor / ETI): Restores minimal action trajectory matching WT bounds (S_min <= 6.5 kcal/mol)
        self.assertTrue(rescued["onsager_machlup_criterion_satisfied"])
        self.assertEqual(rescued["gating_regime"], "ACTIVE_INSTANTON_CROSSING")
        self.assertLessEqual(rescued["instanton_action_s_min_kcal_mol"], 6.5)
        self.assertGreaterEqual(rescued["instanton_traversal_probability"], 0.85)
        self.assertGreaterEqual(rescued["drift_velocity_v_mu_m_s"], 1.0)

        # W1282X: Passive equilibrium baseline
        self.assertFalse(w1282x["onsager_machlup_criterion_satisfied"])
        self.assertEqual(w1282x["gating_regime"], "PASSIVE_EQUILIBRIUM_APO")
        self.assertEqual(w1282x["instanton_action_s_min_kcal_mol"], 0.0)
        self.assertEqual(w1282x["instanton_traversal_probability"], 0.0)

    def test_formal_proof_dag_action_geodesic_parity(self) -> None:
        challenges = self.report["challenge_profiles"]
        expected_benchmarks = [
            "anthropics/fermats-last-theorem",
            "anthropics/formal-math/zeta23",
            "anthropics/formal-math/3sum-apsp",
        ]
        for c_id in expected_benchmarks:
            self.assertIn(c_id, challenges)
            c_prof = challenges[c_id]
            self.assertTrue(c_prof["instanton_proof_geodesic_satisfied"])
            self.assertTrue(c_prof["euler_lagrange_stationarity_verified"])
            self.assertGreaterEqual(c_prof["action_barrier_reduction_pct"], 50.0)
            self.assertLess(c_prof["minimal_tactic_action_s_om"], 5.0)

    def test_ah10_negative_controls(self) -> None:
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_invalid_regimes_count_raises_value_error(self) -> None:
        with self.assertRaises(ValueError):
            evaluate_onsager_machlup_action(n_gating_regimes=2)


if __name__ == "__main__":
    unittest.main()
