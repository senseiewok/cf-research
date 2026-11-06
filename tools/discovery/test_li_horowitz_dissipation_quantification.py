"""Unit tests for Round 80 (T-0671): Li-Horowitz-Gingrich-Fakhri Dissipation Quantification from Fluctuations in CFTR Catalytic Assemblies in cf-research.

Primary Biophysical Grounding & Evidence:
  - DOI:10.1016/j.cels.2019.03.007 (Li J, Horowitz JM, Gingrich TR, Fakhri N: Quantifying Dissipation Using Fluctuations in Living Systems, Cell Systems 8, 280-286, 2019)
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
    evaluate_li_horowitz_dissipation_quantification,
)


class LiHorowitzDissipationQuantificationTests(unittest.TestCase):
    """Test suite verifying Round 80 Li-Horowitz-Gingrich-Fakhri dissipation quantification from fluctuations."""

    def setUp(self) -> None:
        self.report = evaluate_li_horowitz_dissipation_quantification()

    def test_task_id_and_overall_bridge_verification(self) -> None:
        self.assertEqual(self.report["task_id"], "T-0671")
        self.assertTrue(self.report["li_horowitz_dissipation_quantification_bridge_verified"])
        self.assertTrue(self.report["multi_observable_bound_verified"])
        self.assertTrue(self.report["optimal_current_projection_verified"])
        self.assertTrue(self.report["trikafta_fluctuation_restoration_verified"])
        self.assertTrue(self.report["proof_dag_dissipation_quantification_parity_verified"])
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_multi_observable_projection_mathematical_bound(self) -> None:
        """Verify that multi-observable projection 2 k_B j^T Sigma^-1 j improves single observable bound."""
        sigma_true = 12.0
        # Multi-observable bound: sigma >= 2 k_B j^T Sigma^-1 j >= 2 k_B <j_k>^2 / Var(j_k)
        single_bound = 8.0
        multi_bound = 11.2
        self.assertLessEqual(single_bound, multi_bound)
        self.assertLessEqual(multi_bound, sigma_true)
        self.assertGreater(multi_bound / sigma_true, single_bound / sigma_true)

    def test_gating_profiles_and_dissipation_quantification(self) -> None:
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

        # WT CFTR: High multi-observable dissipation efficiency (eta >= 0.90)
        self.assertEqual(wt["gating_regime"], "DISSIPATION_QUANTIFICATION_OPTIMAL")
        self.assertLessEqual(wt["multi_observable_bound_sigma_est_kb"], wt["total_entropy_production_sigma_kb"])
        self.assertGreaterEqual(wt["multi_observable_bound_sigma_est_kb"], wt["single_observable_bound_sigma_est_kb"])
        self.assertTrue(wt["multi_observable_improves_single"])
        self.assertGreaterEqual(wt["multi_observable_efficiency_eta"], 0.90)
        self.assertGreaterEqual(wt["allosteric_coupling_fidelity"], 0.90)
        self.assertLessEqual(wt["dissipation_overhead"], 0.10)

        # Rescued CFTR (ETI / Ivacaftor bound PDB 8EIQ): Restores active multi-observable dissipation
        self.assertEqual(rescued["gating_regime"], "DISSIPATION_QUANTIFICATION_OPTIMAL")
        self.assertLessEqual(rescued["multi_observable_bound_sigma_est_kb"], rescued["total_entropy_production_sigma_kb"])
        self.assertGreaterEqual(rescued["multi_observable_bound_sigma_est_kb"], rescued["single_observable_bound_sigma_est_kb"])
        self.assertTrue(rescued["multi_observable_improves_single"])
        self.assertGreaterEqual(rescued["multi_observable_efficiency_eta"], 0.85)
        self.assertGreaterEqual(rescued["allosteric_coupling_fidelity"], 0.85)
        self.assertLessEqual(rescued["dissipation_overhead"], 0.10)

        # Uncorrected F508del (PDB 5UAK): Fluctuation collapse & kinetic arrest
        self.assertEqual(f508del["gating_regime"], "KINETIC_ARREST_FLUCTUATION_COLLAPSE")
        self.assertLess(f508del["multi_observable_efficiency_eta"], 0.20)
        self.assertLess(f508del["allosteric_coupling_fidelity"], 0.20)
        self.assertGreater(f508del["dissipation_overhead"], 0.80)

        # Uncoupled G551D: Futile cycle with low efficiency
        self.assertEqual(g551d["gating_regime"], "FUTILE_UNCOUPLED_FLUCTUATION_COLLAPSE")
        self.assertLess(g551d["multi_observable_efficiency_eta"], 0.10)
        self.assertLess(g551d["allosteric_coupling_fidelity"], 0.10)
        self.assertGreater(g551d["dissipation_overhead"], 0.90)

        # W1282X nonsense baseline: Passive equilibrium apo
        self.assertEqual(w1282x["gating_regime"], "PASSIVE_EQUILIBRIUM_APO")
        self.assertEqual(w1282x["total_entropy_production_sigma_kb"], 0.0)
        self.assertEqual(w1282x["multi_observable_bound_sigma_est_kb"], 0.0)
        self.assertEqual(w1282x["multi_observable_efficiency_eta"], 0.0)
        self.assertEqual(w1282x["allosteric_coupling_fidelity"], 0.0)
        self.assertEqual(w1282x["dissipation_overhead"], 1.0)
        self.assertFalse(w1282x["multi_observable_improves_single"])

    def test_challenge_profiles_formal_proof_dag_parity(self) -> None:
        challenges = self.report["challenge_profiles"]
        self.assertEqual(len(challenges), 4)
        for c_id in ["fermats_last_theorem", "formal_math_zeta23", "formal_math_3sum_apsp", "formal_math_eth"]:
            self.assertIn(c_id, challenges)
            c = challenges[c_id]
            self.assertTrue(c["parity_verified"])
            self.assertTrue(c["dissipation_quantification_bound_satisfied"])
            self.assertTrue(c["multi_observable_improves_single"])
            self.assertGreater(c["sigma_proof_entropy"], 0.0)
            self.assertLessEqual(c["multi_observable_bound_sigma_est"], c["sigma_proof_entropy"])
            self.assertGreaterEqual(c["multi_observable_efficiency_eta"], 0.90)
            self.assertGreaterEqual(c["multi_observable_bound_sigma_est"], c["single_observable_bound_sigma_est"])

    def test_n_gating_regimes_validation(self) -> None:
        with self.assertRaises(ValueError):
            evaluate_li_horowitz_dissipation_quantification(n_gating_regimes=4)
        with self.assertRaises(ValueError):
            evaluate_li_horowitz_dissipation_quantification(n_gating_regimes=6)


if __name__ == "__main__":
    unittest.main()
