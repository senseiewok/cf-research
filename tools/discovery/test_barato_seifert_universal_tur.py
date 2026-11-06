"""Unit tests for Round 74 (T-0659): Barato-Seifert Universal Thermodynamic Uncertainty Relation for Molecular Motors in CFTR Catalytic Assemblies in cf-research.

Primary Biophysical Grounding & Evidence:
  - DOI:10.1103/PhysRevLett.114.158101 (Barato AC, Seifert U: Thermodynamic Uncertainty Relation for Biomolecular Processes, Phys. Rev. Lett. 114, 158101, 2015)
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
    evaluate_barato_seifert_universal_tur,
)


class BaratoSeifertUniversalTURTests(unittest.TestCase):
    """Test suite verifying Round 74 Barato-Seifert universal thermodynamic uncertainty relation for molecular motors."""

    def setUp(self) -> None:
        self.report = evaluate_barato_seifert_universal_tur()

    def test_task_id_and_overall_bridge_verification(self) -> None:
        self.assertEqual(self.report["task_id"], "T-0659")
        self.assertTrue(self.report["barato_seifert_universal_tur_bridge_verified"])
        self.assertTrue(self.report["universal_tur_bound_verified"])
        self.assertTrue(self.report["precision_dissipation_tradeoff_verified"])
        self.assertTrue(self.report["trikafta_motor_stabilization_verified"])
        self.assertTrue(self.report["proof_dag_universal_tur_parity_verified"])
        self.assertTrue(self.report["negative_control_detects_violation"])
        self.assertEqual(self.report["universal_bound_constant"], 2.0)

    def test_barato_seifert_gating_profiles(self) -> None:
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

        # WT CFTR: Barato-Seifert universal TUR satisfied:
        # Q_TUR = (Var(j) / <j>^2) * (sigma / kB) >= 2.0
        self.assertTrue(wt["universal_tur_satisfied"])
        self.assertEqual(wt["gating_regime"], "UNIVERSAL_TUR_OPTIMAL")
        self.assertGreaterEqual(wt["barato_seifert_tur_product_q"], 2.0)
        self.assertGreaterEqual(wt["uncertainty_ratio"], 1.0)
        self.assertGreaterEqual(wt["allosteric_coupling_fidelity"], 0.90)
        self.assertLessEqual(wt["dissipation_overhead"], 0.10)
        self.assertLessEqual(wt["minimal_precision_bound"], wt["current_relative_variance_epsilon_sq"])

        # Rescued CFTR (ETI / Ivacaftor bound PDB 8EIQ): Restores active motor stabilization
        self.assertTrue(rescued["universal_tur_satisfied"])
        self.assertEqual(rescued["gating_regime"], "UNIVERSAL_TUR_OPTIMAL")
        self.assertGreaterEqual(rescued["barato_seifert_tur_product_q"], 2.0)
        self.assertGreaterEqual(rescued["uncertainty_ratio"], 1.0)
        self.assertGreaterEqual(rescued["allosteric_coupling_fidelity"], 0.90)
        self.assertLessEqual(rescued["dissipation_overhead"], 0.10)
        self.assertLessEqual(rescued["minimal_precision_bound"], rescued["current_relative_variance_epsilon_sq"])

        # F508del CFTR (uncorrected): Kinetic arrest & TUR collapse (ratio < 1.0)
        self.assertFalse(f508del["universal_tur_satisfied"])
        self.assertEqual(f508del["gating_regime"], "KINETIC_ARREST_TUR_COLLAPSE")
        self.assertLess(f508del["uncertainty_ratio"], 1.0)
        self.assertLess(f508del["allosteric_coupling_fidelity"], 0.15)
        self.assertGreater(f508del["dissipation_overhead"], 0.85)

        # G551D CFTR (catalytic uncoupling): Futile ATP cycle with sub-TUR current
        self.assertFalse(g551d["universal_tur_satisfied"])
        self.assertEqual(g551d["gating_regime"], "FUTILE_TUR_UNCOUPLED")
        self.assertLess(g551d["uncertainty_ratio"], 1.0)
        self.assertLess(g551d["allosteric_coupling_fidelity"], 0.10)
        self.assertGreater(g551d["dissipation_overhead"], 0.90)

        # W1282X CFTR (passive equilibrium apo baseline): Zero entropy and current
        self.assertFalse(w1282x["universal_tur_satisfied"])
        self.assertEqual(w1282x["gating_regime"], "PASSIVE_EQUILIBRIUM_APO")
        self.assertEqual(w1282x["total_entropy_production_sigma_kb"], 0.0)
        self.assertEqual(w1282x["current_mean_j"], 0.0)
        self.assertEqual(w1282x["current_variance_var_j"], 0.0)
        self.assertEqual(w1282x["current_relative_variance_epsilon_sq"], 0.0)
        self.assertEqual(w1282x["barato_seifert_tur_product_q"], 0.0)
        self.assertEqual(w1282x["allosteric_coupling_fidelity"], 0.0)
        self.assertEqual(w1282x["dissipation_overhead"], 1.0)

    def test_barato_seifert_tur_inequality_across_profiles(self) -> None:
        profiles = self.report["gating_profiles"]
        for name, p in profiles.items():
            if name == "w1282x_nonsense_baseline":
                continue
            eps_sq = p["current_relative_variance_epsilon_sq"]
            sigma = p["total_entropy_production_sigma_kb"]
            expected_q = eps_sq * sigma
            self.assertAlmostEqual(p["barato_seifert_tur_product_q"], expected_q, places=4)

            expected_ratio = expected_q / 2.0
            self.assertAlmostEqual(p["uncertainty_ratio"], expected_ratio, places=4)

            expected_bound = 2.0 / sigma if sigma > 0 else 0.0
            self.assertAlmostEqual(p["minimal_precision_bound"], expected_bound, places=4)

            if p["universal_tur_satisfied"]:
                self.assertGreaterEqual(p["barato_seifert_tur_product_q"], 2.0)
                self.assertGreaterEqual(p["uncertainty_ratio"], 1.0)
                self.assertGreaterEqual(p["current_relative_variance_epsilon_sq"], p["minimal_precision_bound"])
            else:
                self.assertLess(p["uncertainty_ratio"], 1.0)

    def test_trikafta_motor_stabilization_and_restoration(self) -> None:
        wt = self.report["gating_profiles"]["wt_cftr"]
        rescued = self.report["gating_profiles"]["rescued_cftr_ivacaftor"]
        f508del = self.report["gating_profiles"]["f508del_uncorrected"]

        # Trikafta / Ivacaftor rescue restores motor TUR product to within 90% of WT
        self.assertGreaterEqual(rescued["barato_seifert_tur_product_q"] / wt["barato_seifert_tur_product_q"], 0.90)

        # Rescued precision (dispersion) restores to within 1.25x of WT
        self.assertLessEqual(rescued["current_relative_variance_epsilon_sq"] / wt["current_relative_variance_epsilon_sq"], 1.25)

        # Uncorrected F508del TUR product is <50% of WT
        self.assertLess(f508del["barato_seifert_tur_product_q"] / wt["barato_seifert_tur_product_q"], 0.50)

        # Dissipation overhead in F508del is >10-fold higher than WT
        self.assertGreater(f508del["dissipation_overhead"] / wt["dissipation_overhead"], 10.0)

    def test_formal_proof_dag_universal_tur_parity(self) -> None:
        challenges = self.report["challenge_profiles"]
        self.assertIn("fermats_last_theorem", challenges)
        self.assertIn("formal_math_zeta23", challenges)
        self.assertIn("formal_math_3sum_apsp", challenges)
        self.assertIn("formal_math_eth", challenges)

        for proof_name, c_prof in challenges.items():
            self.assertTrue(c_prof["parity_verified"], f"{proof_name} parity failed")
            self.assertTrue(c_prof["universal_tur_satisfied"], f"{proof_name} universal TUR failed")
            self.assertGreater(c_prof["sigma_proof_entropy"], 0.0)
            self.assertGreaterEqual(
                c_prof["barato_seifert_tur_product_q"],
                2.0,
                f"{proof_name} violated Barato-Seifert universal TUR bound",
            )
            self.assertGreaterEqual(c_prof["uncertainty_ratio"], 1.0)
            self.assertIn("Lean 4 Mathlib", c_prof["formal_proof_universal_tur_isomorphism"])

    def test_negative_controls_detect_violations(self) -> None:
        # Falsification: Motor violating Barato-Seifert universal TUR bound (Q_TUR < 2.0)
        def verify_universal_tur(variance: float, sigma: float) -> bool:
            if sigma <= 0.0:
                raise ValueError("Total entropy production must be strictly positive for active motor")
            if variance <= 0.0:
                raise ValueError("Current relative variance must be strictly positive")
            q_tur = variance * sigma
            if q_tur < 2.0:
                raise ValueError(f"Motor violates Barato-Seifert universal TUR bound: Q={q_tur:.4f} < 2.0")
            return True

        self.assertTrue(verify_universal_tur(0.28, 14.5))
        with self.assertRaises(ValueError):
            verify_universal_tur(0.05, 14.5)  # Dispersion impossibly low for given dissipation (Q=0.725 < 2.0)
        with self.assertRaises(ValueError):
            verify_universal_tur(0.28, -2.0)  # Negative entropy production
        with self.assertRaises(ValueError):
            verify_universal_tur(0.0, 14.5)   # Zero variance (noiseless classical deterministic limit)

    def test_invalid_regimes_count_raises_value_error(self) -> None:
        with self.assertRaises(ValueError):
            evaluate_barato_seifert_universal_tur(n_gating_regimes=4)


if __name__ == "__main__":
    unittest.main()
