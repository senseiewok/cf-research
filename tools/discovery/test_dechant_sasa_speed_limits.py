"""Unit tests for Round 52 (T-0607): Dechant-Sasa Thermodynamic Speed Limits & Geometric Bounds in CFTR Gating in cf-research.

Primary Biophysical Grounding & Evidence:
  - DOI:10.1088/1742-5468/ac0b4e (Dechant A, Sasa S: Continuous time-reversal and thermodynamic uncertainty relations, J. Stat. Mech. 2021)
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
    evaluate_dechant_sasa_speed_limits,
)


class DechantSasaSpeedLimitsTests(unittest.TestCase):
    """Test suite verifying Round 52 Dechant-Sasa thermodynamic speed limits and geometric bounds for CFTR gating."""

    def setUp(self) -> None:
        self.report = evaluate_dechant_sasa_speed_limits()

    def test_task_id_and_overall_bridge_verification(self) -> None:
        self.assertEqual(self.report["task_id"], "T-0607")
        self.assertTrue(self.report["dechant_sasa_speed_limits_bridge_verified"])
        self.assertTrue(self.report["speed_limit_bounds_verified"])
        self.assertTrue(self.report["geometric_transport_stratification_verified"])
        self.assertTrue(self.report["proof_dag_speed_limit_parity_verified"])
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_speed_limit_bounds_and_transport_stratification(self) -> None:
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

        # WT CFTR: Dechant-Sasa speed limit satisfied (tau_gate >= W_1 / sqrt(2 D sigma), tau_gate >= 12 ms)
        self.assertTrue(wt["speed_limit_bound_satisfied"])
        self.assertEqual(wt["gating_regime"], "GEOMETRIC_SPEED_LIMIT_OPTIMAL")
        self.assertGreaterEqual(wt["conformational_transition_time_tau_ms"], 12.0)
        self.assertGreaterEqual(
            wt["conformational_transition_time_tau_ms"],
            wt["dechant_sasa_lower_bound_ms"],
        )
        self.assertGreaterEqual(wt["allosteric_coupling_fidelity"], 0.90)

        # Rescued CFTR (Ivacaftor / ETI): Restores transport bounds matching WT
        self.assertTrue(rescued["speed_limit_bound_satisfied"])
        self.assertEqual(rescued["gating_regime"], "GEOMETRIC_SPEED_LIMIT_OPTIMAL")
        self.assertGreaterEqual(rescued["conformational_transition_time_tau_ms"], 12.0)
        self.assertGreaterEqual(
            rescued["conformational_transition_time_tau_ms"],
            rescued["dechant_sasa_lower_bound_ms"],
        )
        self.assertGreaterEqual(rescued["allosteric_coupling_fidelity"], 0.90)

        # F508del: Kinetic trapping and diverging transition time (tau_gate > 500 ms)
        self.assertFalse(f508del["speed_limit_bound_satisfied"])
        self.assertEqual(f508del["gating_regime"], "KINETIC_TRAPPING_DIVERGENCE")
        self.assertGreater(f508del["conformational_transition_time_tau_ms"], 500.0)
        self.assertLess(f508del["wasserstein_distance_w1_nm"], 0.15)
        self.assertLess(f508del["allosteric_coupling_fidelity"], 0.15)

        # G551D: Kinetic trapping and diverging transition time (tau_gate > 500 ms)
        self.assertFalse(g551d["speed_limit_bound_satisfied"])
        self.assertEqual(g551d["gating_regime"], "KINETIC_TRAPPING_DIVERGENCE")
        self.assertGreater(g551d["conformational_transition_time_tau_ms"], 500.0)
        self.assertLess(g551d["wasserstein_distance_w1_nm"], 0.15)
        self.assertLess(g551d["allosteric_coupling_fidelity"], 0.15)

        # W1282X: Passive equilibrium baseline (tau = W_1 = D = sigma = 0)
        self.assertFalse(w1282x["speed_limit_bound_satisfied"])
        self.assertEqual(w1282x["gating_regime"], "PASSIVE_EQUILIBRIUM_APO")
        self.assertEqual(w1282x["conformational_transition_time_tau_ms"], 0.0)
        self.assertEqual(w1282x["wasserstein_distance_w1_nm"], 0.0)
        self.assertEqual(w1282x["diffusion_coefficient_d_nm2_per_ms"], 0.0)
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
            self.assertGreaterEqual(c_prof["proof_transition_time_step"], 5.0)
            self.assertLessEqual(c_prof["proof_wasserstein_distance_w1"], 8.0)

    def test_ah10_negative_controls(self) -> None:
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_invalid_regimes_count_raises_value_error(self) -> None:
        with self.assertRaises(ValueError):
            evaluate_dechant_sasa_speed_limits(n_gating_regimes=3)


if __name__ == "__main__":
    unittest.main()
