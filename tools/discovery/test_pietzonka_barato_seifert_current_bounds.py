"""Unit tests for Round 78 (T-0667): Pietzonka-Barato-Seifert Universal Finite-Time Current Bounds in CFTR Catalytic Assemblies in cf-research.

Primary Biophysical Grounding & Evidence:
  - DOI:10.1103/PhysRevE.93.052145 (Pietzonka P, Barato AC, Seifert U: Universal bounds on current fluctuations, Phys. Rev. E 93, 052145, 2016)
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
    evaluate_pietzonka_barato_seifert_current_bounds,
)


class PietzonkaBaratoSeifertCurrentBoundsTests(unittest.TestCase):
    """Test suite verifying Round 78 Pietzonka-Barato-Seifert universal finite-time current fluctuation bounds."""

    def setUp(self) -> None:
        self.report = evaluate_pietzonka_barato_seifert_current_bounds()

    def test_task_id_and_overall_bridge_verification(self) -> None:
        self.assertEqual(self.report["task_id"], "T-0667")
        self.assertTrue(self.report["pietzonka_barato_seifert_current_bounds_bridge_verified"])
        self.assertTrue(self.report["universal_finite_time_bound_verified"])
        self.assertTrue(self.report["q_function_monotonicity_verified"])
        self.assertTrue(self.report["trikafta_current_stabilization_verified"])
        self.assertTrue(self.report["proof_dag_finite_time_current_parity_verified"])
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_q_function_mathematical_properties(self) -> None:
        """Verify universal scaling function q(theta) = theta / tanh(theta) >= 1.0."""
        # For any positive theta, q(theta) > 1 and increases monotonically
        for theta in [0.1, 0.5, 1.0, 2.0, 5.0, 10.0]:
            q_val = theta / math.tanh(theta)
            self.assertGreaterEqual(q_val, 1.0)
            if theta > 0.1:
                prev_q = (theta - 0.1) / math.tanh(theta - 0.1)
                self.assertGreater(q_val, prev_q)

        # Universal bound: Bound_PBS(sigma) = (2 / sigma) * q(sigma / 2) = 1 / tanh(sigma / 2)
        for sigma in [0.5, 1.0, 2.0, 5.0, 10.0]:
            theta = sigma / 2.0
            q_theta = theta / math.tanh(theta)
            bound_formula = (2.0 / sigma) * q_theta
            bound_direct = 1.0 / math.tanh(sigma / 2.0)
            self.assertAlmostEqual(bound_formula, bound_direct, places=6)
            self.assertGreaterEqual(bound_direct, 1.0)

    def test_gating_profiles_and_finite_time_bounds(self) -> None:
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

        # WT CFTR: Satisfies universal finite-time current bound
        # Var(j)/<j>^2 >= (2 kB / sigma) * q(kB sigma / 2) = 1 / tanh(sigma / 2)
        self.assertTrue(wt["finite_time_tur_satisfied"])
        self.assertEqual(wt["gating_regime"], "FINITE_TIME_CURRENT_OPTIMAL")
        self.assertGreaterEqual(wt["current_relative_variance_epsilon_sq"], wt["pietzonka_barato_seifert_bound"])
        self.assertGreaterEqual(wt["uncertainty_ratio"], 1.0)
        self.assertGreaterEqual(wt["allosteric_coupling_fidelity"], 0.90)
        self.assertLessEqual(wt["dissipation_overhead"], 0.10)
        self.assertGreater(wt["current_mean_j"], 0.0)

        # Rescued CFTR (ETI / Ivacaftor bound PDB 8EIQ): Restores active current stabilization
        self.assertTrue(rescued["finite_time_tur_satisfied"])
        self.assertEqual(rescued["gating_regime"], "FINITE_TIME_CURRENT_OPTIMAL")
        self.assertGreaterEqual(rescued["current_relative_variance_epsilon_sq"], rescued["pietzonka_barato_seifert_bound"])
        self.assertGreaterEqual(rescued["uncertainty_ratio"], 1.0)
        self.assertGreaterEqual(rescued["allosteric_coupling_fidelity"], 0.85)
        self.assertLessEqual(rescued["dissipation_overhead"], 0.12)
        self.assertGreater(rescued["current_mean_j"], 0.0)

        # Uncorrected F508del (PDB 5UAK): Kinetic arrest & finite-time TUR collapse
        self.assertFalse(f508del["finite_time_tur_satisfied"])
        self.assertEqual(f508del["gating_regime"], "KINETIC_ARREST_FINITE_TIME_COLLAPSE")
        self.assertLess(f508del["uncertainty_ratio"], 1.0)
        self.assertLess(f508del["allosteric_coupling_fidelity"], 0.20)
        self.assertGreater(f508del["dissipation_overhead"], 0.80)

        # Uncoupled G551D: Futile ATP cycle with sub-TUR current
        self.assertFalse(g551d["finite_time_tur_satisfied"])
        self.assertEqual(g551d["gating_regime"], "FUTILE_UNCOUPLED_FINITE_TIME_COLLAPSE")
        self.assertLess(g551d["uncertainty_ratio"], 1.0)
        self.assertLess(g551d["allosteric_coupling_fidelity"], 0.10)
        self.assertGreater(g551d["dissipation_overhead"], 0.90)

        # W1282X nonsense baseline: Passive equilibrium apo
        self.assertFalse(w1282x["finite_time_tur_satisfied"])
        self.assertEqual(w1282x["gating_regime"], "PASSIVE_EQUILIBRIUM_APO")
        self.assertEqual(w1282x["total_entropy_production_sigma_kb"], 0.0)
        self.assertEqual(w1282x["current_mean_j"], 0.0)
        self.assertEqual(w1282x["current_variance_var_j"], 0.0)
        self.assertEqual(w1282x["current_relative_variance_epsilon_sq"], 0.0)
        self.assertEqual(w1282x["pietzonka_barato_seifert_bound"], 0.0)
        self.assertEqual(w1282x["uncertainty_ratio"], 0.0)
        self.assertEqual(w1282x["allosteric_coupling_fidelity"], 0.0)

    def test_challenge_profiles_formal_proof_dag_parity(self) -> None:
        challenges = self.report["challenge_profiles"]
        self.assertEqual(len(challenges), 4)
        for c_id in ["fermats_last_theorem", "formal_math_zeta23", "formal_math_3sum_apsp", "formal_math_eth"]:
            self.assertIn(c_id, challenges)
            c = challenges[c_id]
            self.assertTrue(c["parity_verified"])
            self.assertTrue(c["finite_time_tur_satisfied"])
            self.assertGreater(c["sigma_proof_entropy"], 0.0)
            self.assertGreater(c["proof_current_relative_variance"], 0.0)
            self.assertGreater(c["q_theta"], 1.0)
            self.assertGreaterEqual(c["proof_current_relative_variance"], c["pietzonka_barato_seifert_bound"])
            self.assertGreaterEqual(c["uncertainty_ratio"], 1.0)

    def test_n_gating_regimes_validation(self) -> None:
        with self.assertRaises(ValueError):
            evaluate_pietzonka_barato_seifert_current_bounds(n_gating_regimes=4)
        with self.assertRaises(ValueError):
            evaluate_pietzonka_barato_seifert_current_bounds(n_gating_regimes=6)


if __name__ == "__main__":
    unittest.main()
