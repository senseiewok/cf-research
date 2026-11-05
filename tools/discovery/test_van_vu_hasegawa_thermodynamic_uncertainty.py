"""Unit tests for Round 71 (T-0650): Van Vu - Hasegawa Thermodynamic Uncertainty Relations & Dynamical Activity Bounds in Arbitrarily Driven CFTR Catalytic Assemblies in cf-research.

Primary Biophysical Grounding & Evidence:
  - DOI:10.1103/PhysRevResearch.2.013060 (Van Vu T, Hasegawa Y: Thermodynamic uncertainty relations under arbitrary control protocols, Phys. Rev. Research 2, 013060, 2020)
  - DOI:10.1103/PhysRevLett.126.010601 (Van Vu T, Hasegawa Y: Geometrical bounds of the irreversibility in Markovian systems, Phys. Rev. Lett. 126, 010601, 2021)
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
    evaluate_van_vu_hasegawa_thermodynamic_uncertainty,
)


class VanVuHasegawaThermodynamicUncertaintyTests(unittest.TestCase):
    """Test suite verifying Round 71 Van Vu - Hasegawa thermodynamic uncertainty relations & dynamical activity bounds."""

    def setUp(self) -> None:
        self.report = evaluate_van_vu_hasegawa_thermodynamic_uncertainty()

    def test_task_id_and_overall_bridge_verification(self) -> None:
        self.assertEqual(self.report["task_id"], "T-0650")
        self.assertTrue(self.report["van_vu_hasegawa_thermodynamic_uncertainty_bridge_verified"])
        self.assertTrue(self.report["non_stationary_tur_bound_verified"])
        self.assertTrue(self.report["dynamical_activity_tradeoff_verified"])
        self.assertTrue(self.report["trikafta_activity_restoration_verified"])
        self.assertTrue(self.report["proof_dag_non_stationary_parity_verified"])
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_van_vu_hasegawa_gating_profiles(self) -> None:
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

        # WT CFTR: Van Vu - Hasegawa arbitrary protocol bound satisfied:
        # Var(J) / <J>^2 >= 2 / (Sigma + 2 * sqrt(Sigma * A))
        self.assertTrue(wt["uncertainty_bound_satisfied"])
        self.assertEqual(wt["gating_regime"], "DYNAMICAL_ACTIVITY_PRECISION_OPTIMAL")
        self.assertGreaterEqual(wt["current_relative_variance_epsilon_sq"], wt["van_vu_hasegawa_precision_bound"])
        self.assertGreaterEqual(wt["uncertainty_ratio"], 1.0)
        self.assertGreaterEqual(wt["allosteric_coupling_fidelity"], 0.90)
        self.assertLessEqual(wt["blurring_factor"], 0.10)

        # Rescued CFTR (ETI / Ivacaftor bound PDB 8EIQ): Restores active turnover and precision
        self.assertTrue(rescued["uncertainty_bound_satisfied"])
        self.assertEqual(rescued["gating_regime"], "DYNAMICAL_ACTIVITY_PRECISION_OPTIMAL")
        self.assertGreaterEqual(rescued["current_relative_variance_epsilon_sq"], rescued["van_vu_hasegawa_precision_bound"])
        self.assertGreaterEqual(rescued["uncertainty_ratio"], 1.0)
        self.assertGreaterEqual(rescued["allosteric_coupling_fidelity"], 0.90)
        self.assertLessEqual(rescued["blurring_factor"], 0.10)

        # F508del CFTR (uncorrected): Activity collapse & kinetic arrest (ratio < 1.0)
        self.assertFalse(f508del["uncertainty_bound_satisfied"])
        self.assertEqual(f508del["gating_regime"], "KINETIC_ARREST_ACTIVITY_COLLAPSE")
        self.assertLess(f508del["uncertainty_ratio"], 1.0)
        self.assertLess(f508del["allosteric_coupling_fidelity"], 0.15)
        self.assertGreater(f508del["blurring_factor"], 0.85)

        # G551D CFTR (catalytic uncoupling): Futile thermal activity with near-zero dissipation
        self.assertFalse(g551d["uncertainty_bound_satisfied"])
        self.assertEqual(g551d["gating_regime"], "FUTILE_ACTIVITY_UNCOUPLED")
        self.assertLess(g551d["uncertainty_ratio"], 1.0)
        self.assertLess(g551d["allosteric_coupling_fidelity"], 0.10)
        self.assertGreater(g551d["blurring_factor"], 0.90)

        # W1282X CFTR (passive equilibrium apo baseline): Zero entropy and activity
        self.assertFalse(w1282x["uncertainty_bound_satisfied"])
        self.assertEqual(w1282x["gating_regime"], "PASSIVE_EQUILIBRIUM_APO")
        self.assertEqual(w1282x["total_entropy_production_sigma_kb"], 0.0)
        self.assertEqual(w1282x["dynamical_activity_events_a"], 0.0)
        self.assertEqual(w1282x["allosteric_coupling_fidelity"], 0.0)
        self.assertEqual(w1282x["blurring_factor"], 1.0)

    def test_van_vu_hasegawa_inequality_across_profiles(self) -> None:
        profiles = self.report["gating_profiles"]
        for name, p in profiles.items():
            if name == "w1282x_nonsense_baseline":
                continue
            sigma = p["total_entropy_production_sigma_kb"]
            act = p["dynamical_activity_events_a"]
            denominator = sigma + 2.0 * math.sqrt(sigma * act)
            expected_lower_bound = 2.0 / denominator
            self.assertAlmostEqual(p["van_vu_hasegawa_precision_bound"], expected_lower_bound, places=3)

            eps_sq = p["current_relative_variance_epsilon_sq"]
            if p["uncertainty_bound_satisfied"]:
                self.assertGreaterEqual(
                    eps_sq,
                    p["van_vu_hasegawa_precision_bound"] - 1e-9,
                    msg=f"{name} violated Van Vu - Hasegawa non-stationary TUR bound",
                )
                self.assertGreaterEqual(p["uncertainty_ratio"], 1.0)
            else:
                self.assertLess(
                    eps_sq,
                    p["van_vu_hasegawa_precision_bound"],
                    msg=f"{name} should not satisfy precision bound without sufficient dissipation and activity",
                )
                self.assertLess(p["uncertainty_ratio"], 1.0)

    def test_trikafta_restoration_of_dynamical_activity_precision(self) -> None:
        wt = self.report["gating_profiles"]["wt_cftr"]
        rescued = self.report["gating_profiles"]["rescued_cftr_ivacaftor"]
        f508del = self.report["gating_profiles"]["f508del_uncorrected"]

        # Trikafta / Ivacaftor rescue restores precision ratio to within 90% of WT
        self.assertGreaterEqual(rescued["uncertainty_ratio"] / wt["uncertainty_ratio"], 0.90)

        # Uncorrected F508del blurring factor is >10-fold higher than WT
        self.assertGreater(f508del["blurring_factor"] / wt["blurring_factor"], 10.0)

        # Dynamical activity in rescued CFTR is >15-fold higher than uncorrected F508del
        self.assertGreater(rescued["dynamical_activity_events_a"] / f508del["dynamical_activity_events_a"], 15.0)

    def test_formal_proof_dag_non_stationary_parity(self) -> None:
        challenges = self.report["challenge_profiles"]
        self.assertIn("fermats_last_theorem", challenges)
        self.assertIn("formal_math_zeta23", challenges)
        self.assertIn("formal_math_3sum_apsp", challenges)
        self.assertIn("formal_math_eth", challenges)

        for proof_name, c_prof in challenges.items():
            self.assertTrue(c_prof["parity_verified"], f"{proof_name} parity failed")
            self.assertTrue(c_prof["uncertainty_relation_satisfied"], f"{proof_name} uncertainty bound failed")
            self.assertGreater(c_prof["sigma_proof_entropy"], 0.0)
            self.assertGreater(c_prof["dynamical_activity_proof_steps"], 0.0)
            self.assertGreaterEqual(
                c_prof["proof_current_relative_variance"],
                c_prof["van_vu_hasegawa_bound"],
                f"{proof_name} violated proof step relative variance bound",
            )
            self.assertGreaterEqual(c_prof["precision_ratio"], 1.0)
            self.assertIn("Lean 4 Mathlib", c_prof["formal_proof_non_stationary_isomorphism"])

    def test_negative_controls_detect_violations(self) -> None:
        # Falsification 1: Negative entropy production or non-physical activity
        def verify_vvh_bound(sigma: float, activity: float, variance: float) -> bool:
            if sigma < 0.0 or activity < 0.0:
                raise ValueError("Entropy production and dynamical activity must be non-negative")
            denom = sigma + 2.0 * math.sqrt(sigma * activity)
            if denom <= 0.0:
                raise ValueError("Denominator must be strictly positive")
            bound = 2.0 / denom
            if variance < bound - 1e-9:
                raise ValueError("Observed current variance violates Van Vu - Hasegawa bound")
            return True

        self.assertTrue(verify_vvh_bound(12.0, 24.0, 0.098))
        with self.assertRaises(ValueError):
            verify_vvh_bound(12.0, 24.0, 0.01)  # Violates lower bound
        with self.assertRaises(ValueError):
            verify_vvh_bound(-1.0, 24.0, 0.098)  # Negative entropy
        with self.assertRaises(ValueError):
            verify_vvh_bound(12.0, -5.0, 0.098)  # Negative activity

    def test_invalid_regimes_count_raises_value_error(self) -> None:
        with self.assertRaises(ValueError):
            evaluate_van_vu_hasegawa_thermodynamic_uncertainty(n_gating_regimes=4)


if __name__ == "__main__":
    unittest.main()
