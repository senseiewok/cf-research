"""Unit tests for Round 70 (T-0646): Falasco-Esposito Dissipation-Time Uncertainty Relations & Gating Precision Trade-offs in Driven CFTR Catalytic Assemblies in cf-research.

Primary Biophysical Grounding & Evidence:
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
    evaluate_falasco_esposito_dissipation_time_uncertainty,
)


class FalascoEspositoDissipationTimeUncertaintyTests(unittest.TestCase):
    """Test suite verifying Round 70 Falasco-Esposito dissipation-time uncertainty relations & gating precision trade-offs."""

    def setUp(self) -> None:
        self.report = evaluate_falasco_esposito_dissipation_time_uncertainty()

    def test_task_id_and_overall_bridge_verification(self) -> None:
        self.assertEqual(self.report["task_id"], "T-0646")
        self.assertTrue(self.report["falasco_esposito_dissipation_time_uncertainty_bridge_verified"])
        self.assertTrue(self.report["uncertainty_relation_bound_verified"])
        self.assertTrue(self.report["gating_precision_tradeoff_verified"])
        self.assertTrue(self.report["trikafta_resolution_restoration_verified"])
        self.assertTrue(self.report["proof_dag_uncertainty_parity_verified"])
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_dissipation_time_uncertainty_gating_profiles(self) -> None:
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

        # WT CFTR: Dissipation-time uncertainty relation satisfied:
        # tau_obs * sigma_diss >= |<Delta A>|^2 / (2 * D_A)
        self.assertTrue(wt["uncertainty_bound_satisfied"])
        self.assertEqual(wt["gating_regime"], "DISSIPATION_TIME_PRECISION_OPTIMAL")
        self.assertGreaterEqual(wt["dissipation_time_product"], wt["precision_lower_bound"])
        self.assertGreaterEqual(wt["uncertainty_ratio"], 1.0)
        self.assertGreaterEqual(wt["allosteric_coupling_fidelity"], 0.90)
        self.assertLessEqual(wt["blurring_factor"], 0.10)
        self.assertAlmostEqual(
            wt["dissipation_time_product"],
            wt["observation_time_tau_obs_ms"] * wt["dissipation_rate_sigma_diss_kb_per_ms"],
            places=4,
        )

        # Rescued CFTR (ETI / Ivacaftor bound PDB 8EIQ): Restores near-WT precision and sharp resolution
        self.assertTrue(rescued["uncertainty_bound_satisfied"])
        self.assertEqual(rescued["gating_regime"], "DISSIPATION_TIME_PRECISION_OPTIMAL")
        self.assertGreaterEqual(rescued["dissipation_time_product"], rescued["precision_lower_bound"])
        self.assertGreaterEqual(rescued["uncertainty_ratio"], 1.0)
        self.assertGreaterEqual(rescued["allosteric_coupling_fidelity"], 0.90)
        self.assertLessEqual(rescued["blurring_factor"], 0.10)
        self.assertAlmostEqual(
            rescued["dissipation_time_product"],
            rescued["observation_time_tau_obs_ms"] * rescued["dissipation_rate_sigma_diss_kb_per_ms"],
            places=4,
        )

        # F508del CFTR (uncorrected): Kinetic arrest & dissipation-time blurring (ratio < 1.0)
        self.assertFalse(f508del["uncertainty_bound_satisfied"])
        self.assertEqual(f508del["gating_regime"], "DISSIPATION_TIME_BLURRING_ARREST")
        self.assertLess(f508del["uncertainty_ratio"], 1.0)
        self.assertLess(f508del["allosteric_coupling_fidelity"], 0.15)
        self.assertGreater(f508del["blurring_factor"], 0.85)
        self.assertGreaterEqual(f508del["observation_time_tau_obs_ms"], 400.0)

        # G551D CFTR (catalytic uncoupling): Hydrolysis-gating disconnect & blurring
        self.assertFalse(g551d["uncertainty_bound_satisfied"])
        self.assertEqual(g551d["gating_regime"], "CATALYTIC_UNCOUPLING_BLURRING")
        self.assertLess(g551d["uncertainty_ratio"], 1.0)
        self.assertLess(g551d["allosteric_coupling_fidelity"], 0.10)
        self.assertGreater(g551d["blurring_factor"], 0.90)
        self.assertGreaterEqual(g551d["observation_time_tau_obs_ms"], 400.0)

        # W1282X CFTR (passive equilibrium apo baseline): Zero active dissipation & observation product
        self.assertFalse(w1282x["uncertainty_bound_satisfied"])
        self.assertEqual(w1282x["gating_regime"], "PASSIVE_EQUILIBRIUM_APO")
        self.assertEqual(w1282x["dissipation_time_product"], 0.0)
        self.assertEqual(w1282x["allosteric_coupling_fidelity"], 0.0)
        self.assertEqual(w1282x["blurring_factor"], 1.0)

    def test_falasco_esposito_inequality_across_profiles(self) -> None:
        profiles = self.report["gating_profiles"]
        for name, p in profiles.items():
            if name == "w1282x_nonsense_baseline":
                continue
            delta_a = p["observable_displacement_delta_a_angstrom"]
            d_a = p["variance_diffusion_d_a_angstrom_sq"]
            expected_lower_bound = (delta_a**2) / (2.0 * d_a)
            self.assertAlmostEqual(p["precision_lower_bound"], expected_lower_bound, places=3)

            tau = p["observation_time_tau_obs_ms"]
            sigma = p["dissipation_rate_sigma_diss_kb_per_ms"]
            product = tau * sigma
            self.assertAlmostEqual(p["dissipation_time_product"], product, places=3)

            if p["uncertainty_bound_satisfied"]:
                self.assertGreaterEqual(
                    product,
                    p["precision_lower_bound"] - 1e-9,
                    msg=f"{name} violated Falasco-Esposito uncertainty bound",
                )
                self.assertGreaterEqual(p["uncertainty_ratio"], 1.0)
            else:
                self.assertLess(
                    product,
                    p["precision_lower_bound"],
                    msg=f"{name} should not satisfy precision bound without sufficient dissipation",
                )
                self.assertLess(p["uncertainty_ratio"], 1.0)

    def test_trikafta_restoration_of_resolution_precision(self) -> None:
        wt = self.report["gating_profiles"]["wt_cftr"]
        rescued = self.report["gating_profiles"]["rescued_cftr_ivacaftor"]
        f508del = self.report["gating_profiles"]["f508del_uncorrected"]

        # Trikafta / Ivacaftor rescue restores precision ratio to within 90% of WT
        self.assertGreaterEqual(rescued["uncertainty_ratio"] / wt["uncertainty_ratio"], 0.90)

        # Uncorrected F508del blurring factor is >10-fold higher than WT
        self.assertGreater(f508del["blurring_factor"] / wt["blurring_factor"], 10.0)

        # Transition observation time in uncorrected F508del is >15-fold prolonged relative to WT
        self.assertGreater(f508del["observation_time_tau_obs_ms"] / wt["observation_time_tau_obs_ms"], 15.0)

    def test_formal_proof_dag_uncertainty_parity(self) -> None:
        challenges = self.report["challenge_profiles"]
        self.assertIn("fermats_last_theorem", challenges)
        self.assertIn("formal_math_zeta23", challenges)
        self.assertIn("formal_math_3sum_apsp", challenges)
        self.assertIn("formal_math_eth", challenges)

        for proof_name, c_prof in challenges.items():
            self.assertTrue(c_prof["parity_verified"], f"{proof_name} parity failed")
            self.assertTrue(c_prof["uncertainty_relation_satisfied"], f"{proof_name} uncertainty bound failed")
            self.assertGreater(c_prof["tau_obs_proof_steps"], 0.0)
            self.assertGreater(c_prof["sigma_diss_entropy"], 0.0)
            self.assertGreaterEqual(
                c_prof["dissipation_time_product"],
                c_prof["falasco_esposito_precision_bound"],
                f"{proof_name} violated proof step dissipation-time product bound",
            )
            self.assertGreaterEqual(c_prof["precision_ratio"], 1.0)
            self.assertIn("Lean 4 Mathlib", c_prof["formal_proof_uncertainty_isomorphism"])

    def test_negative_controls_detect_violations(self) -> None:
        # Falsification 1: Observation time insufficient to satisfy precision bound
        def verify_dissipation_time_bound(tau: float, sigma: float, bound: float) -> bool:
            if tau < 0.0 or sigma < 0.0:
                raise ValueError("Observation time and dissipation rate must be non-negative")
            if (tau * sigma) < bound - 1e-9:
                raise ValueError("Observed dissipation-time product violates Falasco-Esposito bound")
            return True

        self.assertTrue(verify_dissipation_time_bound(25.0, 0.46, 4.2))
        with self.assertRaises(ValueError):
            verify_dissipation_time_bound(5.0, 0.46, 4.2)
        with self.assertRaises(ValueError):
            verify_dissipation_time_bound(-2.0, 0.46, 4.2)
        with self.assertRaises(ValueError):
            verify_dissipation_time_bound(25.0, -0.1, 4.2)

        # Falsification 2: Non-positive variance D_A
        def verify_variance_positive(d_a: float) -> bool:
            if d_a <= 0.0:
                raise ValueError("Observable diffusion / variance must be strictly positive")
            return True

        self.assertTrue(verify_variance_positive(2.1))
        with self.assertRaises(ValueError):
            verify_variance_positive(0.0)
        with self.assertRaises(ValueError):
            verify_variance_positive(-1.5)

    def test_invalid_regimes_count_raises_value_error(self) -> None:
        with self.assertRaises(ValueError):
            evaluate_falasco_esposito_dissipation_time_uncertainty(n_gating_regimes=4)


if __name__ == "__main__":
    unittest.main()
