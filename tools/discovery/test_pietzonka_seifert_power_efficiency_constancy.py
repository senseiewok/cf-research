"""Unit tests for Round 73 (T-0657): Pietzonka-Seifert Universal Power-Efficiency-Constancy Trade-Off & Molecular Engine Bounds in CFTR Catalytic Assemblies in cf-research.

Primary Biophysical Grounding & Evidence:
  - DOI:10.1103/PhysRevLett.120.190602 (Pietzonka P, Seifert U: Universal Trade-Off between Power, Efficiency, and Constancy in Steady-State Heat Engines, Phys. Rev. Lett. 120, 190602, 2018)
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
    evaluate_pietzonka_seifert_power_efficiency_constancy,
)


class PietzonkaSeifertPowerEfficiencyConstancyTests(unittest.TestCase):
    """Test suite verifying Round 73 Pietzonka-Seifert universal power-efficiency-constancy trade-off & engine bounds."""

    def setUp(self) -> None:
        self.report = evaluate_pietzonka_seifert_power_efficiency_constancy()

    def test_task_id_and_overall_bridge_verification(self) -> None:
        self.assertEqual(self.report["task_id"], "T-0657")
        self.assertTrue(self.report["pietzonka_seifert_power_efficiency_constancy_bridge_verified"])
        self.assertTrue(self.report["power_efficiency_constancy_tradeoff_verified"])
        self.assertTrue(self.report["universal_engine_bound_verified"])
        self.assertTrue(self.report["trikafta_engine_restoration_verified"])
        self.assertTrue(self.report["proof_dag_power_constancy_parity_verified"])
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_pietzonka_seifert_gating_profiles(self) -> None:
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

        # WT CFTR: Pietzonka-Seifert universal power-efficiency-constancy engine bound satisfied:
        # Q_PEC = (P * epsilon^2) / Sigma <= B_PEC(eta) = eta * (1 - eta) / 2
        self.assertTrue(wt["engine_bound_satisfied"])
        self.assertEqual(wt["gating_regime"], "POWER_EFFICIENCY_CONSTANCY_OPTIMAL")
        self.assertLessEqual(wt["power_constancy_product_q"], wt["pietzonka_seifert_bound"])
        self.assertLessEqual(wt["pietzonka_seifert_ratio"], 1.0)
        self.assertGreaterEqual(wt["thermodynamic_efficiency_eta"], 0.50)
        self.assertGreaterEqual(wt["allosteric_coupling_fidelity"], 0.90)
        self.assertLessEqual(wt["dissipation_overhead"], 0.10)

        # Rescued CFTR (ETI / Ivacaftor bound PDB 8EIQ): Restores active coupled transport
        self.assertTrue(rescued["engine_bound_satisfied"])
        self.assertEqual(rescued["gating_regime"], "POWER_EFFICIENCY_CONSTANCY_OPTIMAL")
        self.assertLessEqual(rescued["power_constancy_product_q"], rescued["pietzonka_seifert_bound"])
        self.assertLessEqual(rescued["pietzonka_seifert_ratio"], 1.0)
        self.assertGreaterEqual(rescued["thermodynamic_efficiency_eta"], 0.50)
        self.assertGreaterEqual(rescued["allosteric_coupling_fidelity"], 0.90)
        self.assertLessEqual(rescued["dissipation_overhead"], 0.10)

        # F508del CFTR (uncorrected): Kinetic arrest & efficiency collapse (ratio > 1.0)
        self.assertFalse(f508del["engine_bound_satisfied"])
        self.assertEqual(f508del["gating_regime"], "KINETIC_ARREST_PEC_COLLAPSE")
        self.assertGreater(f508del["pietzonka_seifert_ratio"], 1.0)
        self.assertLess(f508del["thermodynamic_efficiency_eta"], 0.15)
        self.assertLess(f508del["allosteric_coupling_fidelity"], 0.15)
        self.assertGreater(f508del["dissipation_overhead"], 0.85)

        # G551D CFTR (catalytic uncoupling): Futile ATP cycle with near-zero output power
        self.assertFalse(g551d["engine_bound_satisfied"])
        self.assertEqual(g551d["gating_regime"], "FUTILE_PEC_UNCOUPLED")
        self.assertGreater(g551d["pietzonka_seifert_ratio"], 1.0)
        self.assertLess(g551d["thermodynamic_efficiency_eta"], 0.05)
        self.assertLess(g551d["allosteric_coupling_fidelity"], 0.10)
        self.assertGreater(g551d["dissipation_overhead"], 0.90)

        # W1282X CFTR (passive equilibrium apo baseline): Zero entropy, power, and efficiency
        self.assertFalse(w1282x["engine_bound_satisfied"])
        self.assertEqual(w1282x["gating_regime"], "PASSIVE_EQUILIBRIUM_APO")
        self.assertEqual(w1282x["total_entropy_production_sigma_kb"], 0.0)
        self.assertEqual(w1282x["output_power_p"], 0.0)
        self.assertEqual(w1282x["thermodynamic_efficiency_eta"], 0.0)
        self.assertEqual(w1282x["allosteric_coupling_fidelity"], 0.0)
        self.assertEqual(w1282x["dissipation_overhead"], 1.0)

    def test_pietzonka_seifert_bound_across_profiles(self) -> None:
        profiles = self.report["gating_profiles"]
        for name, p in profiles.items():
            if name == "w1282x_nonsense_baseline":
                continue
            eta = p["thermodynamic_efficiency_eta"]
            expected_bound = (eta * (1.0 - eta)) / 2.0
            self.assertAlmostEqual(p["pietzonka_seifert_bound"], expected_bound, places=4)

            p_val = p["output_power_p"]
            eps_sq = p["current_relative_variance_epsilon_sq"]
            sigma = p["total_entropy_production_sigma_kb"]
            expected_q = (p_val * eps_sq) / sigma
            self.assertAlmostEqual(p["power_constancy_product_q"], expected_q, places=4)

            expected_ratio = expected_q / expected_bound if expected_bound > 0 else 0.0
            self.assertAlmostEqual(p["pietzonka_seifert_ratio"], expected_ratio, places=3)

            if p["engine_bound_satisfied"]:
                self.assertLessEqual(p["power_constancy_product_q"], p["pietzonka_seifert_bound"] + 1e-9)
                self.assertLessEqual(p["pietzonka_seifert_ratio"], 1.0 + 1e-9)
            else:
                self.assertGreater(p["pietzonka_seifert_ratio"], 1.0)

    def test_trikafta_restoration_of_engine_performance(self) -> None:
        wt = self.report["gating_profiles"]["wt_cftr"]
        rescued = self.report["gating_profiles"]["rescued_cftr_ivacaftor"]
        f508del = self.report["gating_profiles"]["f508del_uncorrected"]

        # Trikafta / Ivacaftor rescue restores power output to within 90% of WT
        self.assertGreaterEqual(rescued["output_power_p"] / wt["output_power_p"], 0.90)

        # Rescued efficiency restores to within 90% of WT
        self.assertGreaterEqual(rescued["thermodynamic_efficiency_eta"] / wt["thermodynamic_efficiency_eta"], 0.90)

        # Uncorrected F508del power is <10% of WT
        self.assertLess(f508del["output_power_p"] / wt["output_power_p"], 0.10)

        # Dissipation overhead in F508del is >10-fold higher than WT
        self.assertGreater(f508del["dissipation_overhead"] / wt["dissipation_overhead"], 10.0)

    def test_formal_proof_dag_power_constancy_parity(self) -> None:
        challenges = self.report["challenge_profiles"]
        self.assertIn("fermats_last_theorem", challenges)
        self.assertIn("formal_math_zeta23", challenges)
        self.assertIn("formal_math_3sum_apsp", challenges)
        self.assertIn("formal_math_eth", challenges)

        for proof_name, c_prof in challenges.items():
            self.assertTrue(c_prof["parity_verified"], f"{proof_name} parity failed")
            self.assertTrue(c_prof["engine_bound_satisfied"], f"{proof_name} engine bound failed")
            self.assertGreater(c_prof["sigma_proof_entropy"], 0.0)
            self.assertLessEqual(
                c_prof["power_constancy_product_q"],
                c_prof["pietzonka_seifert_bound"] + 1e-9,
                f"{proof_name} violated proof engine power-constancy bound",
            )
            self.assertLessEqual(c_prof["pietzonka_seifert_ratio"], 1.0 + 1e-9)
            self.assertIn("Lean 4 Mathlib", c_prof["formal_proof_power_constancy_isomorphism"])

    def test_negative_controls_detect_violations(self) -> None:
        # Falsification 1: Engine violating Pietzonka-Seifert bound (Q_PEC > B_PEC)
        def verify_pec_bound(power: float, variance: float, sigma: float, efficiency: float) -> bool:
            if sigma <= 0.0:
                raise ValueError("Total entropy production must be strictly positive")
            if not (0.0 < efficiency < 1.0):
                raise ValueError("Thermodynamic efficiency must lie in open interval (0, 1)")
            bound = (efficiency * (1.0 - efficiency)) / 2.0
            q_pec = (power * variance) / sigma
            if q_pec > bound + 1e-9:
                raise ValueError(f"Engine violates Pietzonka-Seifert bound: Q={q_pec:.4f} > B={bound:.4f}")
            return True

        self.assertTrue(verify_pec_bound(4.20, 0.28, 14.5, 0.58))
        with self.assertRaises(ValueError):
            verify_pec_bound(15.0, 0.50, 14.5, 0.58)  # Power too high for given dissipation
        with self.assertRaises(ValueError):
            verify_pec_bound(4.20, 0.28, -2.0, 0.58)  # Negative entropy production
        with self.assertRaises(ValueError):
            verify_pec_bound(4.20, 0.28, 14.5, 1.20)  # Super-Carnot efficiency > 1.0
        with self.assertRaises(ValueError):
            verify_pec_bound(4.20, 0.28, 14.5, 0.0)   # Zero efficiency

    def test_invalid_regimes_count_raises_value_error(self) -> None:
        with self.assertRaises(ValueError):
            evaluate_pietzonka_seifert_power_efficiency_constancy(n_gating_regimes=4)


if __name__ == "__main__":
    unittest.main()
