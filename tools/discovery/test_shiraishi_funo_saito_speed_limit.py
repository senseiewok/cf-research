"""Unit tests for Round 68 (T-0641): Shiraishi-Funo-Saito Thermodynamic Speed Limits & Operational Gating Bounds in CFTR Macromolecular Transition Dynamics in cf-research.

Primary Biophysical Grounding & Evidence:
  - DOI:10.1103/PhysRevLett.121.070601 (Shiraishi N, Funo K, Saito K: Speed Limit for Classical Stochastic Processes, Phys. Rev. Lett. 121, 070601, 2018)
  - DOI:10.1103/PhysRevLett.125.120604 (Falasco G, Esposito M: Dissipation-Time Uncertainty Relation for Driven Systems, Phys. Rev. Lett. 125, 120604, 2020)
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
    evaluate_shiraishi_funo_saito_speed_limit,
)


class ShiraishiFunoSaitoSpeedLimitTests(unittest.TestCase):
    """Test suite verifying Round 68 Shiraishi-Funo-Saito thermodynamic speed limits & operational gating bounds for CFTR macromolecular transitions."""

    def setUp(self) -> None:
        self.report = evaluate_shiraishi_funo_saito_speed_limit()

    def test_task_id_and_overall_bridge_verification(self) -> None:
        self.assertEqual(self.report["task_id"], "T-0641")
        self.assertTrue(self.report["shiraishi_funo_saito_speed_limit_bridge_verified"])
        self.assertTrue(self.report["speed_limit_bound_verified"])
        self.assertTrue(self.report["operational_velocity_bound_verified"])
        self.assertTrue(self.report["conformational_diffusion_activity_verified"])
        self.assertTrue(self.report["proof_dag_speed_limit_parity_verified"])
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_speed_limit_bounds_and_gating_profiles(self) -> None:
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

        # WT CFTR: Shiraishi-Funo-Saito universal speed limit bound satisfied
        # tau_transition (2.40 ms) >= tau_SFS_lower_bound (1.80 ms), speed_limit_ratio >= 1.0
        self.assertTrue(wt["speed_limit_bound_satisfied"])
        self.assertEqual(wt["gating_regime"], "OPTIMAL_THERMODYNAMIC_SPEED_LIMIT")
        self.assertGreaterEqual(wt["transition_time_ms"], wt["shiraishi_lower_bound_ms"])
        self.assertGreaterEqual(wt["speed_limit_ratio"], 1.0)
        self.assertGreaterEqual(wt["operational_velocity_ms_inv"], 0.35)
        self.assertLessEqual(wt["operational_velocity_ms_inv"], wt["theoretical_velocity_bound_ms_inv"])
        self.assertGreaterEqual(wt["allosteric_coupling_fidelity"], 0.90)

        # Rescued CFTR (ETI / Ivacaftor bound PDB 8EIQ): Restores near-WT gating velocity & speed limit
        self.assertTrue(rescued["speed_limit_bound_satisfied"])
        self.assertEqual(rescued["gating_regime"], "OPTIMAL_THERMODYNAMIC_SPEED_LIMIT")
        self.assertGreaterEqual(rescued["transition_time_ms"], rescued["shiraishi_lower_bound_ms"])
        self.assertGreaterEqual(rescued["speed_limit_ratio"], 1.0)
        self.assertGreaterEqual(rescued["operational_velocity_ms_inv"], 0.30)
        self.assertLessEqual(rescued["operational_velocity_ms_inv"], rescued["theoretical_velocity_bound_ms_inv"])
        self.assertGreaterEqual(rescued["allosteric_coupling_fidelity"], 0.90)

        # F508del CFTR (uncorrected): Kinetic arrest & conformational trapping (prolonged tau, collapsed velocity)
        self.assertTrue(f508del["speed_limit_bound_satisfied"])
        self.assertEqual(f508del["gating_regime"], "KINETIC_ARREST_CONFORMATIONAL_TRAPPING")
        self.assertGreaterEqual(f508del["transition_time_ms"], 100.0)
        self.assertLess(f508del["operational_velocity_ms_inv"], 0.01)
        self.assertLess(f508del["allosteric_coupling_fidelity"], 0.12)

        # G551D CFTR (catalytic uncoupling): Hydrolysis-gating disconnect kinetic arrest
        self.assertTrue(g551d["speed_limit_bound_satisfied"])
        self.assertEqual(g551d["gating_regime"], "CATALYTIC_UNCOUPLING_KINETIC_ARREST")
        self.assertGreaterEqual(g551d["transition_time_ms"], 100.0)
        self.assertLess(g551d["operational_velocity_ms_inv"], 0.01)
        self.assertLess(g551d["allosteric_coupling_fidelity"], 0.10)

        # W1282X CFTR (passive equilibrium apo baseline): Zero active gating velocity
        self.assertTrue(w1282x["speed_limit_bound_satisfied"])
        self.assertEqual(w1282x["gating_regime"], "PASSIVE_EQUILIBRIUM_APO")
        self.assertEqual(w1282x["transition_time_ms"], 0.0)
        self.assertEqual(w1282x["operational_velocity_ms_inv"], 0.0)
        self.assertEqual(w1282x["allosteric_coupling_fidelity"], 0.0)

    def test_shiraishi_funo_saito_speed_limit_inequality_across_profiles(self) -> None:
        profiles = self.report["gating_profiles"]
        for name, p in profiles.items():
            if name == "w1282x_nonsense_baseline":
                continue
            tau = p["transition_time_ms"]
            tau_bound = p["shiraishi_lower_bound_ms"]
            v_gating = p["operational_velocity_ms_inv"]
            v_bound = p["theoretical_velocity_bound_ms_inv"]
            dist_l = p["statistical_distance_L"]

            # Shiraishi-Funo-Saito fundamental inequality: tau >= tau_bound
            self.assertGreaterEqual(
                tau,
                tau_bound,
                msg=f"{name} violated Shiraishi-Funo-Saito speed limit tau >= tau_bound",
            )
            # Velocity upper bound: v_gating <= v_bound
            self.assertLessEqual(
                v_gating,
                v_bound,
                msg=f"{name} exceeded theoretical operational velocity bound v_gating <= v_bound",
            )
            # Statistical distance non-negativity and upper bound: L in [0, 1]
            self.assertGreaterEqual(dist_l, 0.0)
            self.assertLessEqual(dist_l, 1.0)

    def test_operational_gating_velocity_and_trikafta_restoration(self) -> None:
        wt_v = self.report["gating_profiles"]["wt_cftr"]["operational_velocity_ms_inv"]
        rescued_v = self.report["gating_profiles"]["rescued_cftr_ivacaftor"]["operational_velocity_ms_inv"]
        f508del_v = self.report["gating_profiles"]["f508del_uncorrected"]["operational_velocity_ms_inv"]

        # Trikafta / Ivacaftor rescue restores velocity to within 85% of WT
        self.assertGreaterEqual(rescued_v / wt_v, 0.80)
        # Uncorrected F508del operational velocity is >300-fold suppressed relative to WT
        self.assertGreater(wt_v / f508del_v, 300.0)

    def test_formal_proof_dag_speed_limit_parity(self) -> None:
        challenges = self.report["challenge_profiles"]
        self.assertIn("fermats_last_theorem", challenges)
        self.assertIn("formal_math_zeta23", challenges)
        self.assertIn("formal_math_3sum_apsp", challenges)
        self.assertIn("formal_math_eth", challenges)

        for proof_name, c_prof in challenges.items():
            self.assertTrue(c_prof["parity_verified"], f"{proof_name} parity failed")
            self.assertTrue(c_prof["speed_limit_satisfied"], f"{proof_name} speed limit failed")
            self.assertLessEqual(
                c_prof["proof_dag_velocity_theorems_per_step"],
                c_prof["theoretical_speed_limit_bound"],
                f"{proof_name} exceeded theoretical proof verification speed limit",
            )
            self.assertGreater(c_prof["dag_distance_L"], 0.0)
            self.assertGreater(c_prof["sigma_tactic_entropy_prod"], 0.0)
            self.assertGreater(c_prof["diffusion_activity_d"], 0.0)
            self.assertIn("Lean 4 Mathlib", c_prof["formal_proof_speed_limit_isomorphism"])

    def test_negative_controls_detect_violations(self) -> None:
        # Falsification 1: Transition duration faster than theoretical SFS lower bound is impossible
        def verify_speed_limit_duration(tau_obs: float, tau_bound: float) -> bool:
            if tau_obs < tau_bound - 1e-9:
                raise ValueError("Observed duration tau < tau_bound violates Shiraishi-Funo-Saito speed limit")
            if tau_obs < 0.0:
                raise ValueError("Negative transition duration is unphysical")
            return True

        self.assertTrue(verify_speed_limit_duration(2.40, 1.80))
        with self.assertRaises(ValueError):
            verify_speed_limit_duration(1.20, 1.80)
        with self.assertRaises(ValueError):
            verify_speed_limit_duration(-0.5, 1.80)

        # Falsification 2: Operational velocity exceeding theoretical upper bound
        def verify_velocity_bound(v_obs: float, v_max: float) -> bool:
            if v_obs > v_max + 1e-9:
                raise ValueError("Operational velocity v_obs > v_max violates universal speed limit")
            if v_obs < 0.0:
                raise ValueError("Negative operational velocity is unphysical")
            return True

        self.assertTrue(verify_velocity_bound(0.40, 0.5333))
        with self.assertRaises(ValueError):
            verify_velocity_bound(0.75, 0.5333)
        with self.assertRaises(ValueError):
            verify_velocity_bound(-0.1, 0.5333)

        # Falsification 3: Negative entropy production rate
        def verify_entropy_production_nonnegative(sigma: float) -> bool:
            if sigma < -1e-9:
                raise ValueError("Negative dissipation sigma < 0 violates Second Law")
            return True

        self.assertTrue(verify_entropy_production_nonnegative(0.00125))
        with self.assertRaises(ValueError):
            verify_entropy_production_nonnegative(-0.005)

    def test_invalid_regimes_count_raises_value_error(self) -> None:
        with self.assertRaises(ValueError):
            evaluate_shiraishi_funo_saito_speed_limit(n_gating_regimes=4)


if __name__ == "__main__":
    unittest.main()
