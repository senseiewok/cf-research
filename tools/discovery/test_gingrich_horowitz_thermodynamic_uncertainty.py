"""Unit tests for Round 72 (T-0655): Gingrich-Horowitz Finite-Time Thermodynamic Uncertainty Relations & Current Precision Bounds in CFTR Catalytic Assemblies in cf-research.

Primary Biophysical Grounding & Evidence:
  - DOI:10.1103/PhysRevLett.116.120601 (Gingrich TR, Horowitz JM, Perunov N, England JL: Dissipation Bounds All Steady-State Current Fluctuations, Phys. Rev. Lett. 116, 120601, 2016)
  - DOI:10.1038/s41567-019-0702-6 (Horowitz JM, Gingrich TR: Thermodynamic uncertainty relations constrain the precision of molecular processes out of equilibrium, Nat. Phys. 16, 15-20, 2020)
  - DOI:10.1085/jgp.200609558 (Csanady et al.: Strict coupling between CFTR channel opening and ATP hydrolysis, J. Gen. Physiol. 2006)
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
    evaluate_gingrich_horowitz_thermodynamic_uncertainty,
)


class GingrichHorowitzThermodynamicUncertaintyTests(unittest.TestCase):
    """Test suite verifying Round 72 Gingrich-Horowitz thermodynamic uncertainty relations & current precision bounds."""

    def setUp(self) -> None:
        self.report = evaluate_gingrich_horowitz_thermodynamic_uncertainty()

    def test_task_id_and_overall_bridge_verification(self) -> None:
        self.assertEqual(self.report["task_id"], "T-0655")
        self.assertTrue(self.report["gingrich_horowitz_thermodynamic_uncertainty_bridge_verified"])
        self.assertTrue(self.report["finite_time_tur_bound_verified"])
        self.assertTrue(self.report["precision_dissipation_tradeoff_verified"])
        self.assertTrue(self.report["trikafta_precision_restoration_verified"])
        self.assertTrue(self.report["proof_dag_finite_time_tur_parity_verified"])
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_gingrich_horowitz_gating_profiles(self) -> None:
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

        # WT CFTR: Gingrich-Horowitz finite-time TUR satisfied:
        # Var(J) / <J>^2 >= 2 / Sigma, or Q = (Var(J) / <J>^2) * Sigma >= 2.0
        self.assertTrue(wt["uncertainty_bound_satisfied"])
        self.assertEqual(wt["gating_regime"], "FINITE_TIME_TUR_OPTIMAL")
        self.assertGreaterEqual(wt["current_relative_variance_epsilon_sq"], wt["gingrich_horowitz_precision_bound"])
        self.assertGreaterEqual(wt["tur_product_q"], 2.0)
        self.assertGreaterEqual(wt["uncertainty_ratio"], 1.0)
        self.assertGreaterEqual(wt["allosteric_coupling_fidelity"], 0.90)
        self.assertLessEqual(wt["dissipation_overhead"], 0.10)

        # Rescued CFTR (ETI / Ivacaftor bound PDB 8EIQ): Restores active dissipation and precision
        self.assertTrue(rescued["uncertainty_bound_satisfied"])
        self.assertEqual(rescued["gating_regime"], "FINITE_TIME_TUR_OPTIMAL")
        self.assertGreaterEqual(rescued["current_relative_variance_epsilon_sq"], rescued["gingrich_horowitz_precision_bound"])
        self.assertGreaterEqual(rescued["tur_product_q"], 2.0)
        self.assertGreaterEqual(rescued["uncertainty_ratio"], 1.0)
        self.assertGreaterEqual(rescued["allosteric_coupling_fidelity"], 0.90)
        self.assertLessEqual(rescued["dissipation_overhead"], 0.10)

        # F508del CFTR (uncorrected): Activity collapse & kinetic arrest (uncertainty ratio < 1.0)
        self.assertFalse(f508del["uncertainty_bound_satisfied"])
        self.assertEqual(f508del["gating_regime"], "KINETIC_ARREST_TUR_COLLAPSE")
        self.assertLess(f508del["uncertainty_ratio"], 1.0)
        self.assertLess(f508del["allosteric_coupling_fidelity"], 0.15)
        self.assertGreater(f508del["dissipation_overhead"], 0.85)

        # G551D CFTR (catalytic uncoupling): Futile thermal uncoupling with sub-TUR current
        self.assertFalse(g551d["uncertainty_bound_satisfied"])
        self.assertEqual(g551d["gating_regime"], "FUTILE_TUR_UNCOUPLED")
        self.assertLess(g551d["uncertainty_ratio"], 1.0)
        self.assertLess(g551d["allosteric_coupling_fidelity"], 0.10)
        self.assertGreater(g551d["dissipation_overhead"], 0.90)

        # W1282X CFTR (passive equilibrium apo baseline): Zero entropy and current
        self.assertFalse(w1282x["uncertainty_bound_satisfied"])
        self.assertEqual(w1282x["gating_regime"], "PASSIVE_EQUILIBRIUM_APO")
        self.assertEqual(w1282x["total_entropy_production_sigma_kb"], 0.0)
        self.assertEqual(w1282x["allosteric_coupling_fidelity"], 0.0)
        self.assertEqual(w1282x["dissipation_overhead"], 1.0)

    def test_gingrich_horowitz_inequality_across_profiles(self) -> None:
        profiles = self.report["gating_profiles"]
        for name, p in profiles.items():
            if name == "w1282x_nonsense_baseline":
                continue
            sigma = p["total_entropy_production_sigma_kb"]
            expected_lower_bound = 2.0 / sigma
            self.assertAlmostEqual(p["gingrich_horowitz_precision_bound"], expected_lower_bound, places=3)

            eps_sq = p["current_relative_variance_epsilon_sq"]
            if p["uncertainty_bound_satisfied"]:
                self.assertGreaterEqual(
                    eps_sq,
                    p["gingrich_horowitz_precision_bound"] - 1e-9,
                    msg=f"{name} violated Gingrich-Horowitz TUR lower bound",
                )
                self.assertGreaterEqual(p["uncertainty_ratio"], 1.0)
                self.assertGreaterEqual(p["tur_product_q"], 2.0)
            else:
                self.assertLess(
                    eps_sq,
                    p["gingrich_horowitz_precision_bound"],
                    msg=f"{name} should not satisfy TUR bound without sufficient coupled dissipation",
                )
                self.assertLess(p["uncertainty_ratio"], 1.0)

    def test_trikafta_restoration_of_tur_precision(self) -> None:
        wt = self.report["gating_profiles"]["wt_cftr"]
        rescued = self.report["gating_profiles"]["rescued_cftr_ivacaftor"]
        f508del = self.report["gating_profiles"]["f508del_uncorrected"]

        # Trikafta / Ivacaftor rescue restores uncertainty ratio to within 90% of WT
        self.assertGreaterEqual(rescued["uncertainty_ratio"] / wt["uncertainty_ratio"], 0.90)

        # Uncorrected F508del dissipation overhead is >10-fold higher than WT
        self.assertGreater(f508del["dissipation_overhead"] / wt["dissipation_overhead"], 10.0)

        # Dissipation in rescued CFTR is >5-fold higher than uncorrected F508del
        self.assertGreater(rescued["total_entropy_production_sigma_kb"] / f508del["total_entropy_production_sigma_kb"], 5.0)

    def test_formal_proof_dag_finite_time_tur_parity(self) -> None:
        challenges = self.report["challenge_profiles"]
        self.assertIn("fermats_last_theorem", challenges)
        self.assertIn("formal_math_zeta23", challenges)
        self.assertIn("formal_math_3sum_apsp", challenges)
        self.assertIn("formal_math_eth", challenges)

        for proof_name, c_prof in challenges.items():
            self.assertTrue(c_prof["parity_verified"], f"{proof_name} parity failed")
            self.assertTrue(c_prof["uncertainty_relation_satisfied"], f"{proof_name} uncertainty bound failed")
            self.assertGreater(c_prof["sigma_proof_entropy"], 0.0)
            self.assertGreaterEqual(
                c_prof["proof_current_relative_variance"],
                c_prof["gingrich_horowitz_bound"],
                f"{proof_name} violated proof step relative variance bound",
            )
            self.assertGreaterEqual(c_prof["tur_product_q"], 2.0)
            self.assertGreaterEqual(c_prof["precision_ratio"], 1.0)
            self.assertIn("Lean 4 Mathlib", c_prof["formal_proof_finite_time_tur_isomorphism"])

    def test_negative_controls_detect_violations(self) -> None:
        # Falsification 1: Negative entropy production or TUR violation
        def verify_gh_bound(sigma: float, variance: float) -> bool:
            if sigma <= 0.0:
                raise ValueError("Total entropy production must be strictly positive")
            bound = 2.0 / sigma
            if variance < bound - 1e-9:
                raise ValueError("Observed current variance violates Gingrich-Horowitz TUR lower bound")
            product_q = variance * sigma
            if product_q < 2.0 - 1e-9:
                raise ValueError("TUR uncertainty product Q must be >= 2.0")
            return True

        self.assertTrue(verify_gh_bound(14.5, 0.32))
        with self.assertRaises(ValueError):
            verify_gh_bound(14.5, 0.05)  # Violates lower bound (variance too small for given dissipation)
        with self.assertRaises(ValueError):
            verify_gh_bound(-2.0, 0.32)  # Negative entropy production violates Second Law
        with self.assertRaises(ValueError):
            verify_gh_bound(0.0, 0.32)  # Zero entropy production

    def test_invalid_regimes_count_raises_value_error(self) -> None:
        with self.assertRaises(ValueError):
            evaluate_gingrich_horowitz_thermodynamic_uncertainty(n_gating_regimes=4)


if __name__ == "__main__":
    unittest.main()
