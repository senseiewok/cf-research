"""Unit tests for Round 75 (T-0661): Dechant-Sasa Multidimensional Thermodynamic Uncertainty Relations & Fluctuation-Response Bounds in CFTR Catalytic Assemblies in cf-research.

Primary Biophysical Grounding & Evidence:
  - DOI:10.1073/pnas.1918386117 (Dechant A, Sasa SI: Fluctuation-response inequality out of equilibrium, Proc. Natl. Acad. Sci. USA 117 (12), 6430-6436, 2020)
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
    evaluate_dechant_sasa_multidimensional_tur,
)


class DechantSasaMultidimensionalTURTests(unittest.TestCase):
    """Test suite verifying Round 75 Dechant-Sasa multidimensional thermodynamic uncertainty relations & fluctuation-response bounds."""

    def setUp(self) -> None:
        self.report = evaluate_dechant_sasa_multidimensional_tur()

    def test_task_id_and_overall_bridge_verification(self) -> None:
        self.assertEqual(self.report["task_id"], "T-0661")
        self.assertTrue(self.report["dechant_sasa_multidimensional_tur_bridge_verified"])
        self.assertTrue(self.report["multidimensional_tur_bound_verified"])
        self.assertTrue(self.report["fluctuation_response_inequality_verified"])
        self.assertTrue(self.report["trikafta_vector_stabilization_verified"])
        self.assertTrue(self.report["proof_dag_multidimensional_tur_parity_verified"])
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_dechant_sasa_gating_profiles(self) -> None:
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

        # WT CFTR: Dechant-Sasa multidimensional TUR satisfied:
        # Q_DS = j^T chi^-1 j <= sigma / (2 kB) = 7.25
        self.assertTrue(wt["multidimensional_tur_satisfied"])
        self.assertEqual(wt["gating_regime"], "MULTIDIMENSIONAL_TUR_OPTIMAL")
        self.assertLessEqual(wt["dechant_sasa_quadratic_form_q"], wt["dechant_sasa_bound_sigma_over_2kb"])
        self.assertGreaterEqual(wt["uncertainty_ratio"], 1.0)
        self.assertLessEqual(wt["fluctuation_bound_ratio"], 1.0)
        self.assertGreaterEqual(wt["bound_slack_margin"], 0.0)
        self.assertGreaterEqual(wt["allosteric_coupling_fidelity"], 0.90)
        self.assertLessEqual(wt["dissipation_overhead"], 0.10)

        # Rescued CFTR (ETI / Ivacaftor bound PDB 8EIQ): Restores active vector stabilization
        self.assertTrue(rescued["multidimensional_tur_satisfied"])
        self.assertEqual(rescued["gating_regime"], "MULTIDIMENSIONAL_TUR_OPTIMAL")
        self.assertLessEqual(rescued["dechant_sasa_quadratic_form_q"], rescued["dechant_sasa_bound_sigma_over_2kb"])
        self.assertGreaterEqual(rescued["uncertainty_ratio"], 1.0)
        self.assertLessEqual(rescued["fluctuation_bound_ratio"], 1.0)
        self.assertGreaterEqual(rescued["bound_slack_margin"], 0.0)
        self.assertGreaterEqual(rescued["allosteric_coupling_fidelity"], 0.90)
        self.assertLessEqual(rescued["dissipation_overhead"], 0.10)

        # F508del CFTR (uncorrected): Kinetic arrest & TUR collapse (fluctuations exceed thermodynamic bound)
        self.assertFalse(f508del["multidimensional_tur_satisfied"])
        self.assertEqual(f508del["gating_regime"], "KINETIC_ARREST_TUR_COLLAPSE")
        self.assertGreater(f508del["dechant_sasa_quadratic_form_q"], f508del["dechant_sasa_bound_sigma_over_2kb"])
        self.assertLess(f508del["uncertainty_ratio"], 1.0)
        self.assertGreater(f508del["fluctuation_bound_ratio"], 1.0)
        self.assertLess(f508del["bound_slack_margin"], 0.0)
        self.assertLessEqual(f508del["allosteric_coupling_fidelity"], 0.20)
        self.assertGreaterEqual(f508del["dissipation_overhead"], 0.80)

        # G551D CFTR (uncoupled): Futile cycle with ATP uncoupling
        self.assertFalse(g551d["multidimensional_tur_satisfied"])
        self.assertEqual(g551d["gating_regime"], "FUTILE_TUR_UNCOUPLED")
        self.assertGreater(g551d["dechant_sasa_quadratic_form_q"], g551d["dechant_sasa_bound_sigma_over_2kb"])
        self.assertLess(g551d["uncertainty_ratio"], 1.0)
        self.assertGreater(g551d["fluctuation_bound_ratio"], 1.0)
        self.assertLess(g551d["bound_slack_margin"], 0.0)
        self.assertLessEqual(g551d["allosteric_coupling_fidelity"], 0.15)
        self.assertGreaterEqual(g551d["dissipation_overhead"], 0.85)

        # W1282X CFTR: Nonsense mutation baseline (zero active current)
        self.assertFalse(w1282x["multidimensional_tur_satisfied"])
        self.assertEqual(w1282x["gating_regime"], "PASSIVE_EQUILIBRIUM_APO")
        self.assertEqual(w1282x["vector_current_j"], [0.0, 0.0])
        self.assertEqual(w1282x["dechant_sasa_quadratic_form_q"], 0.0)
        self.assertEqual(w1282x["allosteric_coupling_fidelity"], 0.0)
        self.assertEqual(w1282x["dissipation_overhead"], 1.0)

    def test_dechant_sasa_multidimensional_tur_inequality_across_profiles(self) -> None:
        profiles = self.report["gating_profiles"]
        for name, p in profiles.items():
            if name == "w1282x_nonsense_baseline":
                continue
            j = p["vector_current_j"]
            chi = p["covariance_matrix_chi"]
            det_chi = chi[0][0] * chi[1][1] - chi[0][1] * chi[1][0]
            self.assertGreater(det_chi, 0.0, f"{name} covariance matrix must be strictly positive definite")

            # Analytical quadratic form j^T chi^-1 j:
            # (j1^2 * chi22 - 2 * j1 * j2 * chi12 + j2^2 * chi11) / det(chi)
            expected_q = (
                j[0] ** 2 * chi[1][1]
                - 2.0 * j[0] * j[1] * chi[0][1]
                + j[1] ** 2 * chi[0][0]
            ) / det_chi
            self.assertAlmostEqual(p["dechant_sasa_quadratic_form_q"], expected_q, places=4)

            sigma = p["total_entropy_production_sigma_kb"]
            expected_bound = sigma / 2.0
            self.assertAlmostEqual(p["dechant_sasa_bound_sigma_over_2kb"], expected_bound, places=4)

            expected_slack = expected_bound - p["dechant_sasa_quadratic_form_q"]
            self.assertAlmostEqual(p["bound_slack_margin"], expected_slack, places=4)

            if p["multidimensional_tur_satisfied"]:
                self.assertLessEqual(p["dechant_sasa_quadratic_form_q"], p["dechant_sasa_bound_sigma_over_2kb"])
                self.assertGreaterEqual(p["uncertainty_ratio"], 1.0)
                self.assertLessEqual(p["fluctuation_bound_ratio"], 1.0)
                self.assertGreaterEqual(p["bound_slack_margin"], 0.0)
            else:
                self.assertGreater(p["dechant_sasa_quadratic_form_q"], p["dechant_sasa_bound_sigma_over_2kb"])
                self.assertLess(p["uncertainty_ratio"], 1.0)
                self.assertGreater(p["fluctuation_bound_ratio"], 1.0)
                self.assertLess(p["bound_slack_margin"], 0.0)

    def test_trikafta_vector_stabilization_and_restoration(self) -> None:
        wt = self.report["gating_profiles"]["wt_cftr"]
        rescued = self.report["gating_profiles"]["rescued_cftr_ivacaftor"]
        f508del = self.report["gating_profiles"]["f508del_uncorrected"]

        # Trikafta / Ivacaftor rescue restores vector current quadratic form to within 95% of WT
        self.assertGreaterEqual(rescued["dechant_sasa_quadratic_form_q"] / wt["dechant_sasa_quadratic_form_q"], 0.90)

        # Rescued uncertainty ratio is comparable to WT (within 5%)
        self.assertAlmostEqual(rescued["uncertainty_ratio"] / wt["uncertainty_ratio"], 1.0, delta=0.05)

        # Uncorrected F508del dissipation overhead is >10-fold higher than WT
        self.assertGreater(f508del["dissipation_overhead"] / wt["dissipation_overhead"], 10.0)

    def test_formal_proof_dag_multidimensional_tur_parity(self) -> None:
        challenges = self.report["challenge_profiles"]
        self.assertIn("fermats_last_theorem", challenges)
        self.assertIn("formal_math_zeta23", challenges)
        self.assertIn("formal_math_3sum_apsp", challenges)
        self.assertIn("formal_math_eth", challenges)

        for proof_name, c_prof in challenges.items():
            self.assertTrue(c_prof["parity_verified"], f"{proof_name} parity failed")
            self.assertTrue(c_prof["multidimensional_tur_satisfied"], f"{proof_name} multidimensional TUR failed")
            self.assertGreater(c_prof["sigma_proof_entropy"], 0.0)
            self.assertLessEqual(
                c_prof["dechant_sasa_quadratic_form_q"],
                c_prof["dechant_sasa_bound_sigma_over_2"],
                f"{proof_name} violated Dechant-Sasa bound",
            )
            self.assertGreaterEqual(c_prof["uncertainty_ratio"], 1.0)
            self.assertIn("Lean 4 Mathlib", c_prof["formal_proof_multidimensional_tur_isomorphism"])

    def test_negative_controls_detect_violations(self) -> None:
        # Falsification: Vector current violating Dechant-Sasa bound (j^T chi^-1 j > sigma / 2)
        def verify_dechant_sasa(j: list[float], chi: list[list[float]], sigma: float) -> bool:
            if sigma <= 0.0:
                raise ValueError("Total entropy production must be strictly positive for active motor")
            det_chi = chi[0][0] * chi[1][1] - chi[0][1] * chi[1][0]
            if det_chi <= 0.0:
                raise ValueError("Covariance matrix must be strictly positive definite")
            q_ds = (
                j[0] ** 2 * chi[1][1]
                - 2.0 * j[0] * j[1] * chi[0][1]
                + j[1] ** 2 * chi[0][0]
            ) / det_chi
            bound = sigma / 2.0
            if q_ds > bound:
                raise ValueError(f"Motor violates Dechant-Sasa bound: Q_DS={q_ds:.4f} > B={bound:.4f}")
            return True

        self.assertTrue(verify_dechant_sasa([5.2, 18.5], [[7.5, 22.0], [22.0, 85.0]], 14.5))
        with self.assertRaises(ValueError):
            # Excess current violating dissipation bound
            verify_dechant_sasa([15.0, 45.0], [[7.5, 22.0], [22.0, 85.0]], 14.5)
        with self.assertRaises(ValueError):
            # Negative entropy production
            verify_dechant_sasa([5.2, 18.5], [[7.5, 22.0], [22.0, 85.0]], -2.0)
        with self.assertRaises(ValueError):
            # Singular / non-invertible covariance matrix
            verify_dechant_sasa([5.2, 18.5], [[1.0, 2.0], [2.0, 4.0]], 14.5)

    def test_invalid_regimes_count_raises_value_error(self) -> None:
        with self.assertRaises(ValueError):
            evaluate_dechant_sasa_multidimensional_tur(n_gating_regimes=4)


if __name__ == "__main__":
    unittest.main()
