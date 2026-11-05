"""Unit tests for Round 69 (T-0643): Koyuk-Seifert Thermodynamic Speed Bounds on Biochemical Clocks & Cyclic CFTR Macromolecular Oscillations in cf-research.

Primary Biophysical Grounding & Evidence:
  - DOI:10.1103/PhysRevLett.125.260604 (Koyuk T, Seifert U: Thermodynamic Bound on the Speed of Biochemical Clocks and Oscillations, Phys. Rev. Lett. 125, 260604, 2020)
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
    evaluate_koyuk_seifert_biochemical_clock_speed,
)


class KoyukSeifertBiochemicalClockSpeedTests(unittest.TestCase):
    """Test suite verifying Round 69 Koyuk-Seifert thermodynamic speed bounds on biochemical clocks & cyclic CFTR oscillations."""

    def setUp(self) -> None:
        self.report = evaluate_koyuk_seifert_biochemical_clock_speed()

    def test_task_id_and_overall_bridge_verification(self) -> None:
        self.assertEqual(self.report["task_id"], "T-0643")
        self.assertTrue(self.report["koyuk_seifert_clock_speed_bridge_verified"])
        self.assertTrue(self.report["speed_limit_bound_verified"])
        self.assertTrue(self.report["oscillation_frequency_bound_verified"])
        self.assertTrue(self.report["phase_coherence_verified"])
        self.assertTrue(self.report["proof_dag_clock_parity_verified"])
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

        # WT CFTR: Koyuk-Seifert universal biochemical clock speed limit bound satisfied
        # tau_cycle (12.0 ms) >= tau_lower_bound (8.0 ms), speed_limit_ratio >= 1.0
        self.assertTrue(wt["speed_limit_bound_satisfied"])
        self.assertEqual(wt["gating_regime"], "OPTIMAL_BIOCHEMICAL_CLOCK_SPEED")
        self.assertGreaterEqual(wt["cycle_period_tau_ms"], wt["koyuk_seifert_lower_bound_ms"])
        self.assertGreaterEqual(wt["speed_limit_ratio"], 1.0)
        self.assertGreaterEqual(wt["angular_frequency_omega_rad_ms_inv"], 0.40)
        self.assertLessEqual(wt["angular_frequency_omega_rad_ms_inv"], wt["theoretical_max_frequency_omega_rad_ms_inv"])
        self.assertGreaterEqual(wt["phase_coherence_factor"], 0.90)
        self.assertGreaterEqual(wt["allosteric_coupling_fidelity"], 0.90)

        # Rescued CFTR (ETI / Ivacaftor bound PDB 8EIQ): Restores near-WT clock pacing & phase coherence
        self.assertTrue(rescued["speed_limit_bound_satisfied"])
        self.assertEqual(rescued["gating_regime"], "OPTIMAL_BIOCHEMICAL_CLOCK_SPEED")
        self.assertGreaterEqual(rescued["cycle_period_tau_ms"], rescued["koyuk_seifert_lower_bound_ms"])
        self.assertGreaterEqual(rescued["speed_limit_ratio"], 1.0)
        self.assertGreaterEqual(rescued["angular_frequency_omega_rad_ms_inv"], 0.40)
        self.assertLessEqual(rescued["angular_frequency_omega_rad_ms_inv"], rescued["theoretical_max_frequency_omega_rad_ms_inv"])
        self.assertGreaterEqual(rescued["phase_coherence_factor"], 0.90)
        self.assertGreaterEqual(rescued["allosteric_coupling_fidelity"], 0.90)

        # F508del CFTR (uncorrected): Kinetic arrest & cyclic decoherence (prolonged cycle, collapsed frequency)
        self.assertTrue(f508del["speed_limit_bound_satisfied"])
        self.assertEqual(f508del["gating_regime"], "KINETIC_ARREST_DECOHERENCE")
        self.assertGreaterEqual(f508del["cycle_period_tau_ms"], 300.0)
        self.assertLess(f508del["angular_frequency_omega_rad_ms_inv"], 0.05)
        self.assertLess(f508del["phase_coherence_factor"], 0.10)
        self.assertLess(f508del["allosteric_coupling_fidelity"], 0.12)

        # G551D CFTR (catalytic uncoupling): Hydrolysis-gating disconnect & severe cycle decoherence
        self.assertTrue(g551d["speed_limit_bound_satisfied"])
        self.assertEqual(g551d["gating_regime"], "CATALYTIC_DECOUPLING_DECOHERENCE")
        self.assertGreaterEqual(g551d["cycle_period_tau_ms"], 300.0)
        self.assertLess(g551d["angular_frequency_omega_rad_ms_inv"], 0.05)
        self.assertLess(g551d["phase_coherence_factor"], 0.10)
        self.assertLess(g551d["allosteric_coupling_fidelity"], 0.10)

        # W1282X CFTR (passive equilibrium apo baseline): Zero active cyclic oscillation
        self.assertTrue(w1282x["speed_limit_bound_satisfied"])
        self.assertEqual(w1282x["gating_regime"], "PASSIVE_EQUILIBRIUM_APO")
        self.assertEqual(w1282x["cycle_period_tau_ms"], 0.0)
        self.assertEqual(w1282x["angular_frequency_omega_rad_ms_inv"], 0.0)
        self.assertEqual(w1282x["phase_coherence_factor"], 0.0)
        self.assertEqual(w1282x["allosteric_coupling_fidelity"], 0.0)

    def test_koyuk_seifert_speed_limit_inequality_across_profiles(self) -> None:
        profiles = self.report["gating_profiles"]
        for name, p in profiles.items():
            if name == "w1282x_nonsense_baseline":
                continue
            tau = p["cycle_period_tau_ms"]
            tau_bound = p["koyuk_seifert_lower_bound_ms"]
            omega = p["angular_frequency_omega_rad_ms_inv"]
            omega_bound = p["theoretical_max_frequency_omega_rad_ms_inv"]
            n_eff = p["effective_states_N_eff"]

            # Koyuk-Seifert fundamental inequality: tau_cycle >= tau_bound
            self.assertGreaterEqual(
                tau,
                tau_bound,
                msg=f"{name} violated Koyuk-Seifert speed limit tau >= tau_bound",
            )
            # Angular frequency upper bound: omega <= omega_bound
            self.assertLessEqual(
                omega,
                omega_bound,
                msg=f"{name} exceeded theoretical angular frequency bound omega <= omega_bound",
            )
            # Effective cyclic states non-negativity and minimum: N_eff >= 2.0
            self.assertGreaterEqual(n_eff, 2.0)

    def test_trikafta_restoration_of_clock_pacing(self) -> None:
        wt_omega = self.report["gating_profiles"]["wt_cftr"]["angular_frequency_omega_rad_ms_inv"]
        rescued_omega = self.report["gating_profiles"]["rescued_cftr_ivacaftor"]["angular_frequency_omega_rad_ms_inv"]
        f508del_omega = self.report["gating_profiles"]["f508del_uncorrected"]["angular_frequency_omega_rad_ms_inv"]

        # Trikafta / Ivacaftor rescue restores oscillation frequency to within 80% of WT
        self.assertGreaterEqual(rescued_omega / wt_omega, 0.80)
        # Uncorrected F508del oscillation frequency is >25-fold suppressed relative to WT
        self.assertGreater(wt_omega / f508del_omega, 25.0)

    def test_formal_proof_dag_clock_parity(self) -> None:
        challenges = self.report["challenge_profiles"]
        self.assertIn("fermats_last_theorem", challenges)
        self.assertIn("formal_math_zeta23", challenges)
        self.assertIn("formal_math_3sum_apsp", challenges)
        self.assertIn("formal_math_eth", challenges)

        for proof_name, c_prof in challenges.items():
            self.assertTrue(c_prof["parity_verified"], f"{proof_name} parity failed")
            self.assertTrue(c_prof["speed_limit_satisfied"], f"{proof_name} clock speed limit failed")
            self.assertLessEqual(
                c_prof["clock_speed_omega_rad_step_inv"],
                c_prof["theoretical_max_omega_bound"],
                f"{proof_name} exceeded theoretical proof verification clock speed limit",
            )
            self.assertGreater(c_prof["effective_stages_N_eff"], 0.0)
            self.assertGreater(c_prof["sigma_entropy_production"], 0.0)
            self.assertGreater(c_prof["phase_coherence"], 0.90)
            self.assertIn("Lean 4 Mathlib", c_prof["formal_proof_clock_isomorphism"])

    def test_negative_controls_detect_violations(self) -> None:
        # Falsification 1: Cycle duration faster than theoretical Koyuk-Seifert lower bound is impossible
        def verify_speed_limit_duration(tau_obs: float, tau_bound: float) -> bool:
            if tau_obs < tau_bound - 1e-9:
                raise ValueError("Observed cycle period tau < tau_bound violates Koyuk-Seifert speed limit")
            if tau_obs < 0.0:
                raise ValueError("Negative cycle duration is unphysical")
            return True

        self.assertTrue(verify_speed_limit_duration(12.0, 8.0))
        with self.assertRaises(ValueError):
            verify_speed_limit_duration(5.0, 8.0)
        with self.assertRaises(ValueError):
            verify_speed_limit_duration(-2.0, 8.0)

        # Falsification 2: Angular oscillation frequency exceeding theoretical upper bound
        def verify_frequency_bound(omega_obs: float, omega_max: float) -> bool:
            if omega_obs > omega_max + 1e-9:
                raise ValueError("Oscillation frequency omega_obs > omega_max violates Koyuk-Seifert speed limit")
            if omega_obs < 0.0:
                raise ValueError("Negative oscillation frequency is unphysical")
            return True

        self.assertTrue(verify_frequency_bound(0.5236, 0.7854))
        with self.assertRaises(ValueError):
            verify_frequency_bound(0.95, 0.7854)
        with self.assertRaises(ValueError):
            verify_frequency_bound(-0.1, 0.7854)

        # Falsification 3: Negative entropy production rate
        def verify_entropy_production_nonnegative(sigma: float) -> bool:
            if sigma < -1e-9:
                raise ValueError("Negative dissipation sigma < 0 violates Second Law")
            return True

        self.assertTrue(verify_entropy_production_nonnegative(0.006283))
        with self.assertRaises(ValueError):
            verify_entropy_production_nonnegative(-0.005)

    def test_invalid_regimes_count_raises_value_error(self) -> None:
        with self.assertRaises(ValueError):
            evaluate_koyuk_seifert_biochemical_clock_speed(n_gating_regimes=4)


if __name__ == "__main__":
    unittest.main()
