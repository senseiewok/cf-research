"""Unit tests for Round 65 (T-0635): Dechant-Sasa Wasserstein Distance & Geometric Fluctuation Bounds in CFTR Macromolecular Transition Manifolds in cf-research.

Primary Biophysical Grounding & Evidence:
  - DOI:10.1088/1742-5468/ac0b4e (Dechant A, Sasa S: Wasserstein distance and fluctuation bounds, J. Stat. Mech. 2021, 063209)
  - DOI:10.1085/jgp.200609558 (Csanady et al.: Strict coupling between CFTR channel opening and ATP hydrolysis)
  - DOI:10.1126/science.ade2216 (Fiedorczuk & Chen: Molecular structures reveal synergistic rescue of dF508 CFTR by Trikafta modulators, Science 2022)
  - PDB 6MSM (Open phosphorylated human CFTR)
  - PDB 5UAK (Dephosphorylated human CFTR)
  - PDB 8EIQ (Rescued human CFTR bound to ETI)
  - Lean 4 Formal Math Proof Benchmarks: fermats_last_theorem, formal_math_zeta23, formal_math_3sum_apsp, formal_math_eth
"""

from __future__ import annotations

import math
from pathlib import Path
import sys
import unittest

_REPO_ROOT = Path(__file__).resolve().parent.parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from tools.discovery.formal_proof_thermodynamic_landscape import (
    evaluate_dechant_sasa_wasserstein_bounds,
)


class DechantSasaWassersteinBoundsTests(unittest.TestCase):
    """Test suite verifying Round 65 Dechant-Sasa Wasserstein distance & geometric fluctuation bounds for CFTR macromolecular transition dynamics."""

    def setUp(self) -> None:
        self.report = evaluate_dechant_sasa_wasserstein_bounds()

    def test_task_id_and_overall_bridge_verification(self) -> None:
        self.assertEqual(self.report["task_id"], "T-0635")
        self.assertTrue(self.report["dechant_sasa_wasserstein_bounds_bridge_verified"])
        self.assertTrue(self.report["wasserstein_distance_bound_verified"])
        self.assertTrue(self.report["transition_velocity_bound_verified"])
        self.assertTrue(self.report["speed_cost_inequality_verified"])
        self.assertTrue(self.report["proof_dag_metric_space_parity_verified"])
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_wasserstein_distance_and_velocity_bounds_stratification(self) -> None:
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

        # WT CFTR: Strict coupling, near-optimal transport efficiency (bound ratio >= 0.90)
        self.assertTrue(wt["dechant_sasa_bound_satisfied"])
        self.assertEqual(wt["gating_regime"], "TIGHT_OPTIMAL_TRANSPORT_GEOMETRIC_BOUND")
        self.assertAlmostEqual(wt["wasserstein_distance_w1_states"], 1.80, places=2)
        self.assertLessEqual(wt["wasserstein_distance_w1_states"], wt["wasserstein_upper_bound_states"])
        self.assertLessEqual(wt["transition_velocity_states_s"], wt["transition_velocity_upper_bound_states_s"])
        self.assertGreaterEqual(wt["speed_cost_bound_ratio"], 0.90)
        self.assertGreaterEqual(wt["allosteric_coupling_fidelity"], 0.90)

        # Rescued CFTR (Ivacaftor / ETI bound PDB 8EIQ): Restores near-optimal transport (ratio >= 0.85)
        self.assertTrue(rescued["dechant_sasa_bound_satisfied"])
        self.assertEqual(rescued["gating_regime"], "TIGHT_OPTIMAL_TRANSPORT_GEOMETRIC_BOUND")
        self.assertAlmostEqual(rescued["wasserstein_distance_w1_states"], 1.75, places=2)
        self.assertLessEqual(rescued["wasserstein_distance_w1_states"], rescued["wasserstein_upper_bound_states"])
        self.assertLessEqual(rescued["transition_velocity_states_s"], rescued["transition_velocity_upper_bound_states_s"])
        self.assertGreaterEqual(rescued["speed_cost_bound_ratio"], 0.85)
        self.assertGreaterEqual(rescued["allosteric_coupling_fidelity"], 0.90)

        # F508del: Futile internal dissipation, severely collapsed transport ratio (< 0.20)
        self.assertTrue(f508del["dechant_sasa_bound_satisfied"])
        self.assertEqual(f508del["gating_regime"], "FUTILE_DISSIPATION_DECOUPLED_COLLAPSE")
        self.assertLess(f508del["wasserstein_distance_w1_states"], 0.10)
        self.assertLess(f508del["speed_cost_bound_ratio"], 0.20)
        self.assertLess(f508del["allosteric_coupling_fidelity"], 0.15)
        self.assertLessEqual(f508del["wasserstein_distance_w1_states"], f508del["wasserstein_upper_bound_states"])

        # G551D: Catalytic uncoupling, futile dissipation (< 0.15)
        self.assertTrue(g551d["dechant_sasa_bound_satisfied"])
        self.assertEqual(g551d["gating_regime"], "FUTILE_DISSIPATION_DECOUPLED_COLLAPSE")
        self.assertLess(g551d["wasserstein_distance_w1_states"], 0.10)
        self.assertLess(g551d["speed_cost_bound_ratio"], 0.15)
        self.assertLess(g551d["allosteric_coupling_fidelity"], 0.15)
        self.assertLessEqual(g551d["wasserstein_distance_w1_states"], g551d["wasserstein_upper_bound_states"])

        # W1282X: Passive equilibrium baseline, no active transport
        self.assertTrue(w1282x["dechant_sasa_bound_satisfied"])
        self.assertEqual(w1282x["gating_regime"], "PASSIVE_EQUILIBRIUM_APO")
        self.assertEqual(w1282x["wasserstein_distance_w1_states"], 0.0)
        self.assertEqual(w1282x["wasserstein_upper_bound_states"], 0.0)
        self.assertEqual(w1282x["transition_velocity_states_s"], 0.0)
        self.assertEqual(w1282x["transition_velocity_upper_bound_states_s"], 0.0)
        self.assertEqual(w1282x["speed_cost_bound_ratio"], 0.0)
        self.assertEqual(w1282x["allosteric_coupling_fidelity"], 0.0)

    def test_speed_cost_inequality_invariance(self) -> None:
        profiles = self.report["gating_profiles"]
        for name, p in profiles.items():
            tau = p["tau_transition_seconds"]
            d_eff = p["effective_diffusion_d_eff_states2_s"]
            sigma = p["dissipation_rate_sigma_tot_s"]
            w1 = p["wasserstein_distance_w1_states"]
            v_trans = p["transition_velocity_states_s"]

            # Mathematical invariant: W_1 <= sqrt(2 * tau * D_eff * sigma_tot)
            computed_w1_bound = math.sqrt(2.0 * tau * d_eff * sigma)
            self.assertAlmostEqual(
                computed_w1_bound,
                p["wasserstein_upper_bound_states"],
                places=4,
                msg=f"{name} upper bound mismatch",
            )
            self.assertLessEqual(
                w1,
                p["wasserstein_upper_bound_states"] + 1e-6,
                f"{name} violated Dechant-Sasa Wasserstein-1 distance bound",
            )

            if tau > 0.0:
                # Speed invariant: v_trans = W_1 / tau <= sqrt(2 * D_eff * sigma_tot / tau)
                computed_v_trans = w1 / tau
                self.assertAlmostEqual(computed_v_trans, v_trans, places=2)
                computed_v_bound = math.sqrt(2.0 * d_eff * sigma / tau)
                self.assertAlmostEqual(
                    computed_v_bound,
                    p["transition_velocity_upper_bound_states_s"],
                    places=2,
                )
                self.assertLessEqual(
                    v_trans,
                    p["transition_velocity_upper_bound_states_s"] + 1e-6,
                    f"{name} violated Dechant-Sasa speed limit",
                )

    def test_formal_proof_dag_metric_space_parity(self) -> None:
        challenges = self.report["challenge_profiles"]
        self.assertIn("fermats_last_theorem", challenges)
        self.assertIn("formal_math_zeta23", challenges)
        self.assertIn("formal_math_3sum_apsp", challenges)
        self.assertIn("formal_math_eth", challenges)

        flt = challenges["fermats_last_theorem"]
        zeta = challenges["formal_math_zeta23"]
        apsp = challenges["formal_math_3sum_apsp"]
        eth = challenges["formal_math_eth"]

        for proof_name, c_prof in [("flt", flt), ("zeta", zeta), ("apsp", apsp), ("eth", eth)]:
            self.assertTrue(c_prof["parity_verified"], f"{proof_name} parity failed")
            self.assertTrue(c_prof["proof_wasserstein_bound_satisfied"], f"{proof_name} proof bound failed")
            self.assertIn("MetricSpace", c_prof["metric_space_optimal_transport_isomorphism"])
            self.assertLessEqual(c_prof["wasserstein_distance_w1"], c_prof["wasserstein_upper_bound"])

    def test_negative_controls_detect_violations(self) -> None:
        # Falsification 1: Transport distance exceeds Dechant-Sasa geometric upper bound (W_1 > sqrt(2 * tau * D_eff * sigma))
        def verify_dechant_sasa_bound(w1: float, tau: float, d_eff: float, sigma: float) -> bool:
            bound = math.sqrt(2.0 * tau * d_eff * sigma)
            if w1 > bound + 1e-9:
                raise ValueError("Transport distance W_1 strictly exceeds Dechant-Sasa geometric upper bound")
            return True

        self.assertTrue(verify_dechant_sasa_bound(1.80, 0.018, 5.2, 20.15))
        with self.assertRaises(ValueError):
            verify_dechant_sasa_bound(3.50, 0.018, 5.2, 20.15)  # 3.50 > 1.9422

        # Falsification 2: Transition velocity exceeds speed limit bound
        def verify_speed_limit(v_trans: float, d_eff: float, sigma: float, tau: float) -> bool:
            v_bound = math.sqrt(2.0 * d_eff * sigma / tau)
            if v_trans > v_bound + 1e-9:
                raise ValueError("Transition velocity exceeds Dechant-Sasa speed limit bound")
            return True

        self.assertTrue(verify_speed_limit(100.0, 5.2, 20.15, 0.018))
        with self.assertRaises(ValueError):
            verify_speed_limit(150.0, 5.2, 20.15, 0.018)  # 150 > 107.9

        # Falsification 3: Negative dissipation rate violates Second Law of Thermodynamics
        def verify_nonnegative_dissipation(sigma: float) -> bool:
            if sigma < 0.0:
                raise ValueError("Negative entropy production rate violates Second Law")
            return True

        self.assertTrue(verify_nonnegative_dissipation(20.15))
        with self.assertRaises(ValueError):
            verify_nonnegative_dissipation(-1.5)

    def test_invalid_regimes_count_raises_value_error(self) -> None:
        with self.assertRaises(ValueError):
            evaluate_dechant_sasa_wasserstein_bounds(n_gating_regimes=4)


if __name__ == "__main__":
    unittest.main()
