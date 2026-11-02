"""Unit tests for Round 31 (T-0575): Graham-Haken Non-Equilibrium Generalized Potential for Multi-Stable CFTR Gating.

Primary Biophysical Grounding & Evidence:
  - DOI:10.1007/BF01394858 (Graham & Haken: Generalized thermodynamic potential for markov systems in detailed balance and far from thermal equilibrium)
  - DOI:10.1085/jgp.200609558 (Csanady et al.: Strict coupling between CFTR channel opening and ATP hydrolysis)
  - DOI:10.1126/science.ade2216 (Fiedorczuk & Chen: Cryo-EM structures of rescued human CFTR, PDB 8EIQ)
  - PDB 6MSM (Open phosphorylated human CFTR)
  - PDB 5UAK (Dephosphorylated human CFTR)
  - Lean 4 Formal Math Proof Benchmarks: fermats_last_theorem, formal_math_zeta23, formal_math_3sum_apsp
"""

from __future__ import annotations

import unittest
from tools.discovery.formal_proof_thermodynamic_landscape import (
    evaluate_graham_haken_potential,
)


class GrahamHakenPotentialTests(unittest.TestCase):
    """Test suite verifying Round 31 Graham-Haken non-equilibrium generalized potential & multi-stable CFTR gating."""

    def setUp(self) -> None:
        self.report = evaluate_graham_haken_potential()

    def test_task_id_and_overall_bridge_verification(self) -> None:
        self.assertEqual(self.report["task_id"], "T-0575")
        self.assertTrue(self.report["graham_haken_bridge_verified"])
        self.assertTrue(self.report["generalized_potential_stratification_verified"])
        self.assertTrue(self.report["proof_dag_potential_reduction_verified"])
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_graham_haken_potential_and_circulation_stratification(self) -> None:
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

        # WT CFTR: Low generalized barrier (Phi <= 4.5 k_B T), active cyclic circulation (r >= 0.50)
        self.assertTrue(wt["graham_haken_criterion_satisfied"])
        self.assertEqual(wt["gating_regime"], "ACTIVE_NON_EQUILIBRIUM_CIRCULATION")
        self.assertLessEqual(wt["phi_generalized_potential_kbt"], 4.5)
        self.assertGreaterEqual(wt["circulation_drift_norm_r"], 0.50)
        self.assertGreaterEqual(wt["opening_transition_rate_k_s_inv"], 1.50)
        self.assertGreater(wt["reversible_drift_norm_d"], 0.0)

        # F508del: High generalized barrier (Phi >= 15.0 k_B T), degenerate reversible collapse (r < 0.10)
        self.assertFalse(f508del["graham_haken_criterion_satisfied"])
        self.assertEqual(f508del["gating_regime"], "DEGENERATE_REVERSIBLE_COLLAPSE")
        self.assertGreaterEqual(f508del["phi_generalized_potential_kbt"], 15.0)
        self.assertLess(f508del["circulation_drift_norm_r"], 0.10)
        self.assertLess(f508del["opening_transition_rate_k_s_inv"], 0.01)

        # G551D: High generalized barrier (Phi >= 15.0 k_B T), degenerate reversible collapse (r < 0.10)
        self.assertFalse(g551d["graham_haken_criterion_satisfied"])
        self.assertEqual(g551d["gating_regime"], "DEGENERATE_REVERSIBLE_COLLAPSE")
        self.assertGreaterEqual(g551d["phi_generalized_potential_kbt"], 15.0)
        self.assertLess(g551d["circulation_drift_norm_r"], 0.10)
        self.assertLess(g551d["opening_transition_rate_k_s_inv"], 0.01)

        # Rescued CFTR (Ivacaftor / ETI): Restores active circulation matching WT (Phi <= 4.5, r >= 0.50)
        self.assertTrue(rescued["graham_haken_criterion_satisfied"])
        self.assertEqual(rescued["gating_regime"], "ACTIVE_NON_EQUILIBRIUM_CIRCULATION")
        self.assertLessEqual(rescued["phi_generalized_potential_kbt"], 4.5)
        self.assertGreaterEqual(rescued["circulation_drift_norm_r"], 0.50)
        self.assertGreaterEqual(rescued["opening_transition_rate_k_s_inv"], 1.50)

        # W1282X: Passive equilibrium baseline
        self.assertFalse(w1282x["graham_haken_criterion_satisfied"])
        self.assertEqual(w1282x["gating_regime"], "PASSIVE_EQUILIBRIUM_APO")
        self.assertEqual(w1282x["phi_generalized_potential_kbt"], 0.0)
        self.assertEqual(w1282x["circulation_drift_norm_r"], 0.0)
        self.assertEqual(w1282x["opening_transition_rate_k_s_inv"], 0.0)

    def test_formal_proof_dag_potential_reduction_parity(self) -> None:
        challenges = self.report["challenge_profiles"]
        expected_benchmarks = [
            "anthropics/fermats-last-theorem",
            "anthropics/formal-math/zeta23",
            "anthropics/formal-math/3sum-apsp",
        ]
        for c_id in expected_benchmarks:
            self.assertIn(c_id, challenges)
            c_prof = challenges[c_id]
            self.assertTrue(c_prof["graham_haken_proof_geodesic_satisfied"])
            self.assertTrue(c_prof["generalized_potential_bound_satisfied"])
            self.assertGreaterEqual(c_prof["action_barrier_reduction_pct"], 50.0)
            self.assertLess(c_prof["phi_generalized_potential_dag"], 4.0)
            self.assertGreaterEqual(c_prof["circulation_flux_r"], 1.0)

    def test_ah10_negative_controls(self) -> None:
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_invalid_regimes_count_raises_value_error(self) -> None:
        with self.assertRaises(ValueError):
            evaluate_graham_haken_potential(n_gating_regimes=2)


if __name__ == "__main__":
    unittest.main()
