"""Unit tests for Round 45 (T-0595): Seifert-Gingrich Current Fluctuation Theorems & Generalized TUR for CFTR Channel Kinetics in cf-research.

Primary Biophysical Grounding & Evidence:
  - DOI:10.1088/0034-4885/75/12/126001 (Seifert U: Stochastic thermodynamics, fluctuation theorems and molecular machines, Rep. Prog. Phys. 75, 126001, 2012; Sec. 3 Fluctuation theorems for currents & Sec. 5.4 Machine efficiency)
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
    evaluate_current_fluctuation_tur,
)


class CurrentFluctuationTURTests(unittest.TestCase):
    """Test suite verifying Round 45 Seifert-Gingrich current fluctuation theorems & generalized TUR for CFTR channel kinetics."""

    def setUp(self) -> None:
        self.report = evaluate_current_fluctuation_tur()

    def test_task_id_and_overall_bridge_verification(self) -> None:
        self.assertEqual(self.report["task_id"], "T-0595")
        self.assertTrue(self.report["current_fluctuation_tur_bridge_verified"])
        self.assertTrue(self.report["tur_stratification_verified"])
        self.assertTrue(self.report["proof_dag_tur_parity_verified"])
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_current_fluctuation_tur_and_gating_stratification(self) -> None:
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

        # WT CFTR: Active driven transport with bounded current variance and tight TUR bound (sigma >= 4.0 kB/s, Q >= 1.0)
        self.assertTrue(wt["current_tur_bound_satisfied"])
        self.assertEqual(wt["gating_regime"], "ACTIVE_TUR_BOUNDED_TRANSPORT")
        self.assertGreaterEqual(wt["dissipation_rate_sigma_kb_per_s"], 4.0)
        self.assertGreaterEqual(wt["tur_ratio_q"], 1.0)
        self.assertGreaterEqual(wt["mean_chloride_current_pa"], 1.5)
        self.assertGreaterEqual(wt["gallavotti_cohen_slope"], 0.95)
        self.assertGreaterEqual(wt["large_deviation_rate_function_i"], 0.5)

        # F508del: Degenerate uncoupled gating with divergent fluctuations (sigma < 0.5 kB/s, Q < 0.1)
        self.assertFalse(f508del["current_tur_bound_satisfied"])
        self.assertEqual(f508del["gating_regime"], "PASSIVE_EQUILIBRIUM_COLLAPSE")
        self.assertLess(f508del["mean_chloride_current_pa"], 0.1)
        self.assertLess(f508del["dissipation_rate_sigma_kb_per_s"], 0.5)
        self.assertLess(f508del["tur_ratio_q"], 0.1)
        self.assertLess(f508del["gallavotti_cohen_slope"], 0.1)

        # G551D: Degenerate futile trap with divergent fluctuations (sigma < 0.5 kB/s, Q < 0.1)
        self.assertFalse(g551d["current_tur_bound_satisfied"])
        self.assertEqual(g551d["gating_regime"], "DEGENERATE_FUTILE_TRAP")
        self.assertLess(g551d["mean_chloride_current_pa"], 0.1)
        self.assertLess(g551d["dissipation_rate_sigma_kb_per_s"], 0.5)
        self.assertLess(g551d["tur_ratio_q"], 0.1)
        self.assertLess(g551d["gallavotti_cohen_slope"], 0.1)

        # Rescued CFTR (Ivacaftor / ETI): Restores coherent current statistics and bounded TUR ratio matching WT
        self.assertTrue(rescued["current_tur_bound_satisfied"])
        self.assertEqual(rescued["gating_regime"], "ACTIVE_TUR_BOUNDED_TRANSPORT")
        self.assertGreaterEqual(rescued["dissipation_rate_sigma_kb_per_s"], 4.0)
        self.assertGreaterEqual(rescued["tur_ratio_q"], 1.0)
        self.assertGreaterEqual(rescued["mean_chloride_current_pa"], 1.5)
        self.assertGreaterEqual(rescued["gallavotti_cohen_slope"], 0.95)

        # W1282X: Passive equilibrium baseline (detailed balance, zero current, zero dissipation)
        self.assertFalse(w1282x["current_tur_bound_satisfied"])
        self.assertEqual(w1282x["gating_regime"], "PASSIVE_EQUILIBRIUM_APO")
        self.assertEqual(w1282x["mean_chloride_current_pa"], 0.0)
        self.assertEqual(w1282x["current_variance_pa2"], 0.0)
        self.assertEqual(w1282x["dissipation_rate_sigma_kb_per_s"], 0.0)
        self.assertEqual(w1282x["tur_ratio_q"], 0.0)
        self.assertEqual(w1282x["large_deviation_rate_function_i"], 0.0)
        self.assertEqual(w1282x["gallavotti_cohen_slope"], 0.0)

    def test_formal_proof_dag_tur_parity(self) -> None:
        challenges = self.report["challenge_profiles"]
        expected_benchmarks = [
            "anthropics/fermats-last-theorem",
            "anthropics/formal-math/zeta23",
            "anthropics/formal-math/3sum-apsp",
        ]
        for c_id in expected_benchmarks:
            self.assertIn(c_id, challenges)
            c_prof = challenges[c_id]
            self.assertTrue(c_prof["tur_proof_parity_satisfied"])
            self.assertTrue(c_prof["seifert_gingrich_bound_satisfied"])
            self.assertGreaterEqual(c_prof["deductive_current_precision"], 2.0)
            self.assertGreaterEqual(c_prof["proof_tur_ratio_q"], 1.0)
            self.assertGreaterEqual(c_prof["inferential_dissipation_rate"], 4.0)
            self.assertGreaterEqual(c_prof["topological_flux_preservation_pct"], 75.0)

    def test_ah10_negative_controls(self) -> None:
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_invalid_regimes_count_raises_value_error(self) -> None:
        with self.assertRaises(ValueError):
            evaluate_current_fluctuation_tur(n_gating_regimes=3)


if __name__ == "__main__":
    unittest.main()
