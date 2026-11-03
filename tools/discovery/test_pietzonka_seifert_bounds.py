"""Unit tests for Round 46 (T-0597): Pietzonka-Seifert Finite-Time Dissipation Bounds & Molecular Machine Efficiency in cf-research.

Primary Biophysical Grounding & Evidence:
  - DOI:10.1088/0034-4885/75/12/126001 (Seifert U: Stochastic thermodynamics, fluctuation theorems and molecular machines, Rep. Prog. Phys. 75, 126001, 2012; Sec. 5 Efficiency and power of molecular machines)
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
    evaluate_pietzonka_seifert_bounds,
)


class PietzonkaSeifertBoundsTests(unittest.TestCase):
    """Test suite verifying Round 46 Pietzonka-Seifert finite-time dissipation bounds & molecular machine efficiency for CFTR channel transitions."""

    def setUp(self) -> None:
        self.report = evaluate_pietzonka_seifert_bounds()

    def test_task_id_and_overall_bridge_verification(self) -> None:
        self.assertEqual(self.report["task_id"], "T-0597")
        self.assertTrue(self.report["pietzonka_seifert_bounds_bridge_verified"])
        self.assertTrue(self.report["efficiency_stratification_verified"])
        self.assertTrue(self.report["proof_dag_efficiency_parity_verified"])
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_pietzonka_seifert_bounds_and_gating_stratification(self) -> None:
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

        # WT CFTR: Efficient active conversion with bounded finite-time dissipation (eta >= 0.45, sigma >= 4.0 kB/s)
        self.assertTrue(wt["finite_time_dissipation_bound_satisfied"])
        self.assertEqual(wt["gating_regime"], "ACTIVE_EFFICIENT_TRANSDUCTION")
        self.assertGreaterEqual(wt["thermodynamic_efficiency_eta"], 0.45)
        self.assertGreaterEqual(wt["maximum_finite_time_efficiency_eta_max"], 0.50)
        self.assertGreaterEqual(wt["dissipation_rate_sigma_kb_per_s"], 4.0)
        self.assertGreaterEqual(wt["atp_coupling_ratio"], 0.90)
        self.assertGreaterEqual(wt["power_output_pw"], 1.5)

        # Rescued CFTR (Ivacaftor / ETI): Restores coherent power output and energetic efficiency matching WT
        self.assertTrue(rescued["finite_time_dissipation_bound_satisfied"])
        self.assertEqual(rescued["gating_regime"], "ACTIVE_EFFICIENT_TRANSDUCTION")
        self.assertGreaterEqual(rescued["thermodynamic_efficiency_eta"], 0.40)
        self.assertGreaterEqual(rescued["dissipation_rate_sigma_kb_per_s"], 4.0)
        self.assertGreaterEqual(rescued["atp_coupling_ratio"], 0.90)
        self.assertGreaterEqual(rescued["power_output_pw"], 1.5)

        # F508del: Futile energetic dissipation with vanishing power output (eta < 0.05, uncoupled ATP turnover)
        self.assertFalse(f508del["finite_time_dissipation_bound_satisfied"])
        self.assertEqual(f508del["gating_regime"], "PASSIVE_EQUILIBRIUM_COLLAPSE")
        self.assertLess(f508del["thermodynamic_efficiency_eta"], 0.05)
        self.assertLess(f508del["power_output_pw"], 0.1)
        self.assertLess(f508del["dissipation_rate_sigma_kb_per_s"], 0.5)
        self.assertLess(f508del["atp_coupling_ratio"], 0.1)

        # G551D: Futile energetic dissipation with vanishing power output (eta < 0.05, uncoupled ATP turnover)
        self.assertFalse(g551d["finite_time_dissipation_bound_satisfied"])
        self.assertEqual(g551d["gating_regime"], "DEGENERATE_FUTILE_TRAP")
        self.assertLess(g551d["thermodynamic_efficiency_eta"], 0.05)
        self.assertLess(g551d["power_output_pw"], 0.1)
        self.assertLess(g551d["dissipation_rate_sigma_kb_per_s"], 0.5)
        self.assertLess(g551d["atp_coupling_ratio"], 0.1)

        # W1282X: Passive equilibrium baseline (detailed balance, zero power, zero dissipation)
        self.assertFalse(w1282x["finite_time_dissipation_bound_satisfied"])
        self.assertEqual(w1282x["gating_regime"], "PASSIVE_EQUILIBRIUM_APO")
        self.assertEqual(w1282x["power_output_pw"], 0.0)
        self.assertEqual(w1282x["power_input_atp_pw"], 0.0)
        self.assertEqual(w1282x["thermodynamic_efficiency_eta"], 0.0)
        self.assertEqual(w1282x["maximum_finite_time_efficiency_eta_max"], 0.0)
        self.assertEqual(w1282x["dissipation_rate_sigma_kb_per_s"], 0.0)
        self.assertEqual(w1282x["atp_coupling_ratio"], 0.0)

    def test_formal_proof_dag_efficiency_parity(self) -> None:
        challenges = self.report["challenge_profiles"]
        expected_benchmarks = [
            "anthropics/fermats-last-theorem",
            "anthropics/formal-math/zeta23",
            "anthropics/formal-math/3sum-apsp",
        ]
        for c_id in expected_benchmarks:
            self.assertIn(c_id, challenges)
            c_prof = challenges[c_id]
            self.assertTrue(c_prof["finite_time_proof_bound_satisfied"])
            self.assertTrue(c_prof["pietzonka_seifert_efficiency_parity_satisfied"])
            self.assertGreaterEqual(c_prof["deductive_power_efficiency"], 0.80)
            self.assertGreaterEqual(c_prof["proof_step_compression_ratio"], 2.0)
            self.assertGreaterEqual(c_prof["inferential_dissipation_rate"], 4.0)

    def test_ah10_negative_controls(self) -> None:
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_invalid_regimes_count_raises_value_error(self) -> None:
        with self.assertRaises(ValueError):
            evaluate_pietzonka_seifert_bounds(n_gating_regimes=3)


if __name__ == "__main__":
    unittest.main()
