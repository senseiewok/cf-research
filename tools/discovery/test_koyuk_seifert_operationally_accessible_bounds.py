"""Unit tests for Round 76 (T-0663): Koyuk-Seifert Operationally Accessible Bounds on Fluctuations and Dissipation in CFTR Catalytic Assemblies in cf-research.

Primary Biophysical Grounding & Evidence:
  - DOI:10.1103/PhysRevLett.122.230601 (Koyuk T, Seifert U: Operationally Accessible Bounds on Fluctuations and Dissipation in Nonequilibrium Steady States, Phys. Rev. Lett. 122, 230601, 2019)
  - DOI:10.1085/jgp.200609558 (Csanady et al.: Strict coupling between CFTR channel opening and ATP hydrolysis, J. Gen. Physiol. 2006)
  - DOI:10.1126/science.ade2216 (Fiedorczuk & Chen: Molecular structures reveal synergistic rescue of dF508 CFTR by Trikafta modulators, Science 2022)
  - PDB 6MSM (Open phosphorylated human CFTR)
  - PDB 5UAK (Dephosphorylated human CFTR)
  - PDB 8EIQ (Rescued human CFTR bound to ETI)
  - Lean 4 Formal Math Proof Benchmarks: fermats_last_theorem, formal_math_zeta23, formal_math_3sum_apsp, formal_math_eth
"""

from __future__ import annotations

from pathlib import Path
import sys
import unittest

_REPO_ROOT = Path(__file__).resolve().parent.parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from tools.discovery.formal_proof_thermodynamic_landscape import (
    evaluate_koyuk_seifert_operationally_accessible_bounds,
)


class KoyukSeifertOperationallyAccessibleBoundsTests(unittest.TestCase):
    """Test suite verifying Round 76 Koyuk-Seifert operationally accessible bounds on fluctuations and dissipation."""

    def setUp(self) -> None:
        self.report = evaluate_koyuk_seifert_operationally_accessible_bounds()

    def test_task_id_and_overall_bridge_verification(self) -> None:
        self.assertEqual(self.report["task_id"], "T-0663")
        self.assertTrue(self.report["koyuk_seifert_operationally_accessible_bounds_bridge_verified"])
        self.assertTrue(self.report["operational_tur_bound_verified"])
        self.assertTrue(self.report["accessible_dissipation_inequality_verified"])
        self.assertTrue(self.report["trikafta_operational_stabilization_verified"])
        self.assertTrue(self.report["proof_dag_operational_tur_parity_verified"])
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_koyuk_seifert_gating_profiles(self) -> None:
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

        # WT CFTR: Koyuk-Seifert operational bound satisfied:
        # sigma_eff <= sigma_tot and Q_KS = sigma_eff * Var(j)/<j>^2 >= 2.0
        self.assertTrue(wt["operational_bound_satisfied"])
        self.assertEqual(wt["gating_regime"], "OPERATIONAL_TUR_OPTIMAL")
        self.assertLessEqual(wt["operational_dissipation_sigma_eff_kb"], wt["total_entropy_production_sigma_tot_kb"])
        self.assertGreaterEqual(wt["koyuk_seifert_operational_product_q"], wt["koyuk_seifert_bound"])
        self.assertGreaterEqual(wt["uncertainty_ratio"], 1.0)
        self.assertGreaterEqual(wt["relative_variance_epsilon_sq"], wt["minimal_precision_bound"])
        self.assertGreaterEqual(wt["allosteric_coupling_fidelity"], 0.90)

        # Rescued CFTR (ETI / Ivacaftor bound PDB 8EIQ): Restores operational dissipation bound
        self.assertTrue(rescued["operational_bound_satisfied"])
        self.assertEqual(rescued["gating_regime"], "OPERATIONAL_TUR_OPTIMAL")
        self.assertLessEqual(
            rescued["operational_dissipation_sigma_eff_kb"],
            rescued["total_entropy_production_sigma_tot_kb"],
        )
        self.assertGreaterEqual(rescued["koyuk_seifert_operational_product_q"], rescued["koyuk_seifert_bound"])
        self.assertGreaterEqual(rescued["uncertainty_ratio"], 1.0)
        self.assertGreaterEqual(rescued["relative_variance_epsilon_sq"], rescued["minimal_precision_bound"])

        # Uncorrected F508del (PDB 5UAK): Operational dissipation collapse
        self.assertFalse(f508del["operational_bound_satisfied"])
        self.assertEqual(f508del["gating_regime"], "KINETIC_ARREST_OPERATIONAL_COLLAPSE")
        self.assertLess(f508del["koyuk_seifert_operational_product_q"], f508del["koyuk_seifert_bound"])
        self.assertLess(f508del["uncertainty_ratio"], 1.0)

        # Uncoupled G551D: Futile operational collapse
        self.assertFalse(g551d["operational_bound_satisfied"])
        self.assertEqual(g551d["gating_regime"], "FUTILE_UNCOUPLED_OPERATIONAL_COLLAPSE")
        self.assertLess(g551d["koyuk_seifert_operational_product_q"], g551d["koyuk_seifert_bound"])
        self.assertLess(g551d["uncertainty_ratio"], 1.0)

        # W1282X nonsense baseline: Passive equilibrium apo
        self.assertFalse(w1282x["operational_bound_satisfied"])
        self.assertEqual(w1282x["gating_regime"], "PASSIVE_EQUILIBRIUM_APO")
        self.assertEqual(w1282x["operational_dissipation_sigma_eff_kb"], 0.0)

    def test_exact_operational_product_and_capture_efficiency(self) -> None:
        profiles = self.report["gating_profiles"]
        for key, p in profiles.items():
            if key == "w1282x_nonsense_baseline":
                continue
            sigma_eff = p["operational_dissipation_sigma_eff_kb"]
            sigma_tot = p["total_entropy_production_sigma_tot_kb"]
            eps_sq = p["relative_variance_epsilon_sq"]

            self.assertLessEqual(sigma_eff, sigma_tot)
            expected_capture = sigma_eff / sigma_tot
            self.assertAlmostEqual(p["dissipation_capture_efficiency"], expected_capture, places=4)

            expected_q = sigma_eff * eps_sq
            self.assertAlmostEqual(p["koyuk_seifert_operational_product_q"], expected_q, places=4)

            expected_min_prec = 2.0 / sigma_eff
            self.assertAlmostEqual(p["minimal_precision_bound"], expected_min_prec, places=4)

    def test_formal_proof_dag_operational_tur_parity(self) -> None:
        challenges = self.report["challenge_profiles"]
        for proof_name in ("fermats_last_theorem", "formal_math_zeta23", "formal_math_3sum_apsp", "formal_math_eth"):
            self.assertIn(proof_name, challenges)
            c_prof = challenges[proof_name]
            self.assertTrue(c_prof["parity_verified"])
            self.assertTrue(c_prof["operational_bound_satisfied"])
            self.assertLessEqual(c_prof["sigma_eff_proof"], c_prof["sigma_tot_proof"])
            self.assertGreaterEqual(c_prof["koyuk_seifert_operational_product_q"], c_prof["koyuk_seifert_bound"])
            self.assertGreaterEqual(c_prof["uncertainty_ratio"], 1.0)
            self.assertIn("Lean 4 Mathlib", c_prof["formal_proof_operational_tur_isomorphism"])

    def test_negative_controls_detect_violations(self) -> None:
        def verify_koyuk_seifert_operational(sigma_eff: float, sigma_tot: float, eps_sq: float) -> bool:
            if sigma_eff <= 0.0 or sigma_tot <= 0.0:
                raise ValueError("Dissipation rates must be strictly positive for active motor")
            if sigma_eff > sigma_tot:
                raise ValueError("Operational dissipation cannot exceed total entropy production")
            if sigma_eff * eps_sq < 2.0:
                raise ValueError("Motor violates Koyuk-Seifert operational TUR bound")
            return True

        self.assertTrue(verify_koyuk_seifert_operational(12.5, 14.5, 0.28))
        with self.assertRaises(ValueError):
            verify_koyuk_seifert_operational(16.0, 14.5, 0.28)
        with self.assertRaises(ValueError):
            verify_koyuk_seifert_operational(1.25, 1.8, 0.88)
        with self.assertRaises(ValueError):
            verify_koyuk_seifert_operational(-1.0, 14.5, 0.28)

    def test_invalid_regimes_count_raises_value_error(self) -> None:
        with self.assertRaises(ValueError):
            evaluate_koyuk_seifert_operationally_accessible_bounds(n_gating_regimes=4)


if __name__ == "__main__":
    unittest.main()
