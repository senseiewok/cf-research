"""Unit tests for Round 44 (T-0593): Dechant-Sasa Wasserstein Fluctuation Dissipation Bounds for CFTR Gating in cf-research.

Primary Biophysical Grounding & Evidence:
  - DOI:10.1088/1742-5468/ac0b4e (Dechant A, Sasa S: Improving thermodynamic uncertainty relations using optimal transport and Wasserstein distance bounds, J. Stat. Mech. 2021, ac0b4e)
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
    evaluate_dechant_sasa_bounds,
)


class DechantSasaBoundsTests(unittest.TestCase):
    """Test suite verifying Round 44 Dechant-Sasa Wasserstein distance & dissipation bounds for CFTR gating transitions."""

    def setUp(self) -> None:
        self.report = evaluate_dechant_sasa_bounds()

    def test_task_id_and_overall_bridge_verification(self) -> None:
        self.assertEqual(self.report["task_id"], "T-0593")
        self.assertTrue(self.report["dechant_sasa_bridge_verified"])
        self.assertTrue(self.report["wasserstein_dissipation_stratification_verified"])
        self.assertTrue(self.report["proof_dag_wasserstein_parity_verified"])
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_wasserstein_dissipation_and_gating_stratification(self) -> None:
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

        # WT CFTR: Active bounded transport (W_2 >= 1.8 states, sigma >= 4.0 kB/s, tight bound satisfied)
        self.assertTrue(wt["dechant_sasa_bound_satisfied"])
        self.assertEqual(wt["gating_regime"], "ACTIVE_BOUNDED_TRANSPORT")
        self.assertGreaterEqual(wt["wasserstein_distance_w2_states"], 1.8)
        self.assertGreaterEqual(wt["dissipation_rate_sigma_kb_per_s"], 4.0)
        self.assertGreater(wt["optimal_transport_cost_kbt"], 7.0)
        self.assertGreater(wt["catalytic_turnover_frequency_hz"], 1.0)

        # F508del: Passive equilibrium collapse (W_2 < 0.1 states, dissipation uncoupled)
        self.assertFalse(f508del["dechant_sasa_bound_satisfied"])
        self.assertEqual(f508del["gating_regime"], "PASSIVE_EQUILIBRIUM_COLLAPSE")
        self.assertLess(f508del["wasserstein_distance_w2_states"], 0.1)
        self.assertLess(f508del["dissipation_rate_sigma_kb_per_s"], 0.5)
        self.assertLess(f508del["optimal_transport_cost_kbt"], 0.2)
        self.assertLess(f508del["catalytic_turnover_frequency_hz"], 0.1)

        # G551D: Degenerate futile trap (W_2 < 0.1 states, dissipation uncoupled)
        self.assertFalse(g551d["dechant_sasa_bound_satisfied"])
        self.assertEqual(g551d["gating_regime"], "DEGENERATE_FUTILE_TRAP")
        self.assertLess(g551d["wasserstein_distance_w2_states"], 0.1)
        self.assertLess(g551d["dissipation_rate_sigma_kb_per_s"], 0.5)
        self.assertLess(g551d["optimal_transport_cost_kbt"], 0.2)
        self.assertLess(g551d["catalytic_turnover_frequency_hz"], 0.1)

        # Rescued CFTR (Ivacaftor / ETI): Restores optimal transport bounds matching WT
        self.assertTrue(rescued["dechant_sasa_bound_satisfied"])
        self.assertEqual(rescued["gating_regime"], "ACTIVE_BOUNDED_TRANSPORT")
        self.assertGreaterEqual(rescued["wasserstein_distance_w2_states"], 1.8)
        self.assertGreaterEqual(rescued["dissipation_rate_sigma_kb_per_s"], 4.0)
        self.assertGreater(rescued["optimal_transport_cost_kbt"], 7.0)
        self.assertGreater(rescued["catalytic_turnover_frequency_hz"], 1.0)

        # W1282X: Passive equilibrium baseline (W_2 == 0, detailed balance)
        self.assertFalse(w1282x["dechant_sasa_bound_satisfied"])
        self.assertEqual(w1282x["gating_regime"], "PASSIVE_EQUILIBRIUM_APO")
        self.assertEqual(w1282x["wasserstein_distance_w2_states"], 0.0)
        self.assertEqual(w1282x["dissipation_rate_sigma_kb_per_s"], 0.0)
        self.assertEqual(w1282x["optimal_transport_cost_kbt"], 0.0)
        self.assertEqual(w1282x["catalytic_turnover_frequency_hz"], 0.0)

    def test_formal_proof_dag_wasserstein_parity(self) -> None:
        challenges = self.report["challenge_profiles"]
        expected_benchmarks = [
            "anthropics/fermats-last-theorem",
            "anthropics/formal-math/zeta23",
            "anthropics/formal-math/3sum-apsp",
        ]
        for c_id in expected_benchmarks:
            self.assertIn(c_id, challenges)
            c_prof = challenges[c_id]
            self.assertTrue(c_prof["wasserstein_proof_parity_satisfied"])
            self.assertTrue(c_prof["dechant_sasa_proof_bound_satisfied"])
            self.assertGreaterEqual(c_prof["proof_wasserstein_graph_edit_distance"], 2.0)
            self.assertGreaterEqual(c_prof["deductive_optimal_transport_cost"], 4.0)
            self.assertGreaterEqual(c_prof["inferential_dissipation_rate"], 4.0)
            self.assertGreaterEqual(c_prof["topological_flux_preservation_pct"], 75.0)

    def test_ah10_negative_controls(self) -> None:
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_invalid_regimes_count_raises_value_error(self) -> None:
        with self.assertRaises(ValueError):
            evaluate_dechant_sasa_bounds(n_gating_regimes=3)


if __name__ == "__main__":
    unittest.main()
