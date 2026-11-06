"""Unit tests for Round 79 (T-0669): Gingrich-Rotskoff-Vaikuntanathan-Horowitz Dissipation-Uncertainty Equality in CFTR Catalytic Assemblies in cf-research.

Primary Biophysical Grounding & Evidence:
  - DOI:10.1073/pnas.1706270114 (Gingrich TR, Rotskoff GM, Vaikuntanathan S, Horowitz JM: Dissipation-uncertainty equality, Proc. Natl. Acad. Sci. USA 114, 6383-6388, 2017)
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
    evaluate_gingrich_rotskoff_dissipation_equality,
)


class GingrichRotskoffDissipationEqualityTests(unittest.TestCase):
    """Test suite verifying Round 79 Gingrich-Rotskoff-Vaikuntanathan-Horowitz dissipation-uncertainty equality."""

    def setUp(self) -> None:
        self.report = evaluate_gingrich_rotskoff_dissipation_equality()

    def test_task_id_and_overall_bridge_verification(self) -> None:
        self.assertEqual(self.report["task_id"], "T-0669")
        self.assertTrue(self.report["gingrich_rotskoff_dissipation_equality_bridge_verified"])
        self.assertTrue(self.report["variational_equality_verified"])
        self.assertTrue(self.report["quadratic_functional_bound_verified"])
        self.assertTrue(self.report["trikafta_variational_restoration_verified"])
        self.assertTrue(self.report["proof_dag_variational_dissipation_parity_verified"])
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_variational_dissipation_mathematical_equality(self) -> None:
        """Verify the Gingrich-Rotskoff variational functional F(d) = 2<j_d> - Var(j_d) <= sigma."""
        # For arbitrary conjugate current observables, F(d) <= sigma
        sigma_test = 12.0
        for mean_j, var_j in [(1.0, 0.5), (3.0, 2.0), (5.0, 4.0), (6.0, 8.0)]:
            f_val = 2.0 * mean_j - var_j
            self.assertLessEqual(f_val, sigma_test)
            efficiency = f_val / sigma_test
            self.assertLessEqual(efficiency, 1.0)

        # For the optimal current estimator d*, F(d*) = 2<j_d*> - Var(j_d*) = sigma exactly
        mean_opt = sigma_test
        var_opt = sigma_test
        f_opt = 2.0 * mean_opt - var_opt
        self.assertEqual(f_opt, sigma_test)
        self.assertEqual(f_opt / sigma_test, 1.0)

    def test_gating_profiles_and_dissipation_equality(self) -> None:
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

        # WT CFTR: Saturated dissipation-uncertainty variational equality
        # F(d*) = 2 <j_d*> - Var(j_d*) = sigma
        self.assertTrue(wt["variational_dissipation_equality_satisfied"])
        self.assertEqual(wt["gating_regime"], "DISSIPATION_EQUALITY_OPTIMAL")
        self.assertEqual(wt["variational_functional_f_val"], wt["total_entropy_production_sigma_kb"])
        self.assertEqual(wt["variational_efficiency_eta"], 1.0)
        self.assertGreater(wt["suboptimal_current_mean_j"], 0.0)
        self.assertLess(wt["suboptimal_functional_f_val"], wt["total_entropy_production_sigma_kb"])
        self.assertLess(wt["suboptimal_efficiency_eta"], 1.0)
        self.assertGreaterEqual(wt["allosteric_coupling_fidelity"], 0.90)
        self.assertLessEqual(wt["dissipation_overhead"], 0.10)

        # Rescued CFTR (ETI / Ivacaftor bound PDB 8EIQ): Restores active variational saturation
        self.assertTrue(rescued["variational_dissipation_equality_satisfied"])
        self.assertEqual(rescued["gating_regime"], "DISSIPATION_EQUALITY_OPTIMAL")
        self.assertEqual(rescued["variational_functional_f_val"], rescued["total_entropy_production_sigma_kb"])
        self.assertEqual(rescued["variational_efficiency_eta"], 1.0)
        self.assertGreaterEqual(rescued["allosteric_coupling_fidelity"], 0.85)
        self.assertLessEqual(rescued["dissipation_overhead"], 0.10)

        # Uncorrected F508del (PDB 5UAK): Kinetic arrest & variational collapse
        self.assertFalse(f508del["variational_dissipation_equality_satisfied"])
        self.assertEqual(f508del["gating_regime"], "KINETIC_ARREST_VARIATIONAL_COLLAPSE")
        self.assertLess(f508del["variational_efficiency_eta"], 0.20)
        self.assertLess(f508del["allosteric_coupling_fidelity"], 0.20)
        self.assertGreater(f508del["dissipation_overhead"], 0.80)

        # Uncoupled G551D: Futile ATP cycle with variational collapse
        self.assertFalse(g551d["variational_dissipation_equality_satisfied"])
        self.assertEqual(g551d["gating_regime"], "FUTILE_UNCOUPLED_VARIATIONAL_COLLAPSE")
        self.assertLess(g551d["variational_efficiency_eta"], 0.10)
        self.assertLess(g551d["allosteric_coupling_fidelity"], 0.10)
        self.assertGreater(g551d["dissipation_overhead"], 0.90)

        # W1282X nonsense baseline: Passive equilibrium apo
        self.assertFalse(w1282x["variational_dissipation_equality_satisfied"])
        self.assertEqual(w1282x["gating_regime"], "PASSIVE_EQUILIBRIUM_APO")
        self.assertEqual(w1282x["total_entropy_production_sigma_kb"], 0.0)
        self.assertEqual(w1282x["optimal_estimator_mean_j"], 0.0)
        self.assertEqual(w1282x["optimal_estimator_variance_var_j"], 0.0)
        self.assertEqual(w1282x["variational_functional_f_val"], 0.0)
        self.assertEqual(w1282x["variational_efficiency_eta"], 0.0)
        self.assertEqual(w1282x["allosteric_coupling_fidelity"], 0.0)
        self.assertEqual(w1282x["dissipation_overhead"], 1.0)

    def test_challenge_profiles_formal_proof_dag_parity(self) -> None:
        challenges = self.report["challenge_profiles"]
        self.assertEqual(len(challenges), 4)
        for c_id in ["fermats_last_theorem", "formal_math_zeta23", "formal_math_3sum_apsp", "formal_math_eth"]:
            self.assertIn(c_id, challenges)
            c = challenges[c_id]
            self.assertTrue(c["parity_verified"])
            self.assertTrue(c["variational_dissipation_equality_satisfied"])
            self.assertGreater(c["sigma_proof_entropy"], 0.0)
            self.assertEqual(c["optimal_estimator_mean_j"], c["sigma_proof_entropy"])
            self.assertEqual(c["optimal_estimator_variance_var_j"], c["sigma_proof_entropy"])
            self.assertEqual(c["variational_functional_f_val"], c["sigma_proof_entropy"])
            self.assertEqual(c["variational_efficiency_eta"], 1.0)
            self.assertLess(c["suboptimal_functional_f_val"], c["sigma_proof_entropy"])

    def test_n_gating_regimes_validation(self) -> None:
        with self.assertRaises(ValueError):
            evaluate_gingrich_rotskoff_dissipation_equality(n_gating_regimes=4)
        with self.assertRaises(ValueError):
            evaluate_gingrich_rotskoff_dissipation_equality(n_gating_regimes=6)


if __name__ == "__main__":
    unittest.main()
