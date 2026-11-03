"""Unit tests for Round 53 (T-0611): Shiraishi-Funo-Saito Universal Thermodynamic Speed Limits & Operational Cost Bounds in CFTR Channel Gating in cf-research.

Primary Biophysical Grounding & Evidence:
  - DOI:10.1103/PhysRevLett.121.070601 (Shiraishi N, Funo K, Saito K: Universal Trade-Off Relation between Speed and Operational Cost for Dissipative Systems, Phys. Rev. Lett. 121, 070601, 2018)
  - DOI:10.1088/0034-4885/75/12/126001 (Seifert U: Stochastic thermodynamics, fluctuation theorems and molecular machines, Rep. Prog. Phys. 75, 126001, 2012)
  - DOI:10.1085/jgp.200609558 (Csanády et al.: Strict coupling between CFTR channel opening and ATP hydrolysis)
  - DOI:10.1126/science.ade2216 (Fiedorczuk & Chen: Molecular structures reveal synergistic rescue of Δ508 CFTR by Trikafta modulators, Science 2022)
  - PDB 6MSM (Open phosphorylated human CFTR)
  - PDB 5UAK (Dephosphorylated human CFTR)
  - PDB 8EIQ (Rescued human CFTR bound to ETI)
  - Lean 4 Formal Math Proof Benchmarks: fermats_last_theorem, formal_math_zeta23, formal_math_3sum_apsp
"""

from __future__ import annotations

from pathlib import Path
import sys
import unittest

_REPO_ROOT = Path(__file__).resolve().parent.parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from tools.discovery.formal_proof_thermodynamic_landscape import (
    evaluate_shiraishi_funo_saito_speed_limits,
)


class ShiraishiFunoSaitoSpeedLimitsTests(unittest.TestCase):
    """Test suite verifying Round 53 Shiraishi-Funo-Saito universal thermodynamic speed limits and operational cost bounds for CFTR gating."""

    def setUp(self) -> None:
        self.report = evaluate_shiraishi_funo_saito_speed_limits()

    def test_task_id_and_overall_bridge_verification(self) -> None:
        self.assertEqual(self.report["task_id"], "T-0611")
        self.assertTrue(self.report["shiraishi_funo_saito_speed_limits_bridge_verified"])
        self.assertTrue(self.report["speed_limit_bounds_verified"])
        self.assertTrue(self.report["operational_cost_stratification_verified"])
        self.assertTrue(self.report["proof_dag_speed_limit_parity_verified"])
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_speed_limit_bounds_and_operational_cost_stratification(self) -> None:
        profiles = self.report["gating_profiles"]
        self.assertIn("wt_cftr", profiles)
        self.assertIn("rescued_cftr_ivacaftor", profiles)
        self.assertIn("f508del_uncorrected", profiles)
        self.assertIn("g551d_uncoupled_cftr", profiles)
        self.assertIn("w1282x_nonsense_baseline", profiles)

        wt = profiles["wt_cftr"]
        rescued = profiles["rescued_cftr_ivacaftor"]
        f508del = profiles["f508del_uncorrected"]
        g551d = profiles["g551d_uncoupled_cftr"]
        w1282x = profiles["w1282x_nonsense_baseline"]

        # WT CFTR: Shiraishi-Funo-Saito speed limit satisfied (tau >= 10.0 ms, operational cost ratio >= 1.0)
        self.assertTrue(wt["speed_limit_bound_satisfied"])
        self.assertEqual(wt["gating_regime"], "DISCRETE_SPEED_LIMIT_OPTIMAL")
        self.assertGreaterEqual(wt["transition_duration_tau_ms"], 10.0)
        self.assertGreaterEqual(
            wt["transition_duration_tau_ms"],
            wt["shiraishi_funo_saito_lower_bound_ms"],
        )
        self.assertGreaterEqual(wt["operational_cost_ratio"], 1.0)
        self.assertGreaterEqual(wt["allosteric_coupling_fidelity"], 0.90)

        # Rescued CFTR (Ivacaftor / ETI): Restores operational bounds matching WT
        self.assertTrue(rescued["speed_limit_bound_satisfied"])
        self.assertEqual(rescued["gating_regime"], "DISCRETE_SPEED_LIMIT_OPTIMAL")
        self.assertGreaterEqual(rescued["transition_duration_tau_ms"], 10.0)
        self.assertGreaterEqual(
            rescued["transition_duration_tau_ms"],
            rescued["shiraishi_funo_saito_lower_bound_ms"],
        )
        self.assertGreaterEqual(rescued["operational_cost_ratio"], 1.0)
        self.assertGreaterEqual(rescued["allosteric_coupling_fidelity"], 0.90)

        # F508del: Kinetic trapping and diverging transition time (tau > 500 ms)
        self.assertFalse(f508del["speed_limit_bound_satisfied"])
        self.assertEqual(f508del["gating_regime"], "KINETIC_TRAPPING_DIVERGENCE")
        self.assertGreater(f508del["transition_duration_tau_ms"], 500.0)
        self.assertLess(f508del["total_variation_distance_dtv"], 0.15)
        self.assertLess(f508del["allosteric_coupling_fidelity"], 0.15)

        # G551D: Kinetic trapping and diverging transition time (tau > 500 ms)
        self.assertFalse(g551d["speed_limit_bound_satisfied"])
        self.assertEqual(g551d["gating_regime"], "KINETIC_TRAPPING_DIVERGENCE")
        self.assertGreater(g551d["transition_duration_tau_ms"], 500.0)
        self.assertLess(g551d["total_variation_distance_dtv"], 0.15)
        self.assertLess(g551d["allosteric_coupling_fidelity"], 0.15)

        # W1282X: Passive equilibrium baseline (tau = dtv = a = sigma = 0)
        self.assertFalse(w1282x["speed_limit_bound_satisfied"])
        self.assertEqual(w1282x["gating_regime"], "PASSIVE_EQUILIBRIUM_APO")
        self.assertEqual(w1282x["transition_duration_tau_ms"], 0.0)
        self.assertEqual(w1282x["total_variation_distance_dtv"], 0.0)
        self.assertEqual(w1282x["dynamical_activity_a_per_ms"], 0.0)
        self.assertEqual(w1282x["entropy_production_rate_sigma_kb_per_ms"], 0.0)
        self.assertEqual(w1282x["allosteric_coupling_fidelity"], 0.0)

    def test_formal_proof_dag_speed_limit_parity(self) -> None:
        challenges = self.report["challenge_profiles"]
        expected_benchmarks = [
            "anthropics/fermats-last-theorem",
            "anthropics/formal-math/zeta23",
            "anthropics/formal-math/3sum-apsp",
        ]
        for c_id in expected_benchmarks:
            self.assertIn(c_id, challenges)
            c_prof = challenges[c_id]
            self.assertTrue(c_prof["proof_speed_limit_parity_satisfied"])
            self.assertGreaterEqual(c_prof["proof_operational_duration_step"], 5.0)
            self.assertLessEqual(c_prof["proof_total_variation_bound"], 8.0)

    def test_ah10_negative_controls(self) -> None:
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_invalid_regimes_count_raises_value_error(self) -> None:
        with self.assertRaises(ValueError):
            evaluate_shiraishi_funo_saito_speed_limits(n_gating_regimes=3)


if __name__ == "__main__":
    unittest.main()
