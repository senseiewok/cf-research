"""Unit tests for Round 49 (T-0603): Sagawa-Ueda Generalized Fluctuation Theorem & Information Efficacy in Allosteric Drug Stabilization in cf-research.

Primary Biophysical Grounding & Evidence:
  - DOI:10.1038/nphys3230 (Parrondo JMR, Horowitz JM, Sagawa T: Thermodynamics of information, Nat. Phys. 11, 131–139, 2015)
  - DOI:10.1103/PhysRevLett.103.090601 (Prost J, Joanny JF, Parrondo JMR: Generalized Fluctuation-Dissipation Theorem for Steady States and Applications to Molecular Motors, PRL 103, 090601, 2009)
  - DOI:10.1088/0034-4885/75/12/126001 (Seifert U: Stochastic thermodynamics, fluctuation theorems and molecular machines, Rep. Prog. Phys. 75, 126001, 2012)
  - DOI:10.1085/jgp.200609558 (Csanády et al.: Strict coupling between CFTR channel opening and ATP hydrolysis)
  - DOI:10.1126/science.ade2216 (Fiedorczuk & Chen: Molecular structures reveal synergistic rescue of Δ508 CFTR by Trikafta modulators, Science 2022)
  - PDB 6MSM (Open phosphorylated human CFTR)
  - PDB 5UAK (Dephosphorylated human CFTR)
  - PDB 8EIQ (Rescued human CFTR bound to ETI)
  - Lean 4 Formal Math Proof Benchmarks: fermats_last_theorem, formal_math_zeta23, formal_math_3sum_apsp
"""

from __future__ import annotations

import unittest
from tools.discovery.formal_proof_thermodynamic_landscape import (
    evaluate_sagawa_ueda_fluctuation_efficiency,
)


class SagawaUedaFluctuationEfficiencyTests(unittest.TestCase):
    """Test suite verifying Round 49 Sagawa-Ueda generalized fluctuation theorem and information efficacy in allosteric drug stabilization for CFTR gating."""

    def setUp(self) -> None:
        self.report = evaluate_sagawa_ueda_fluctuation_efficiency()

    def test_task_id_and_overall_bridge_verification(self) -> None:
        self.assertEqual(self.report["task_id"], "T-0603")
        self.assertTrue(self.report["sagawa_ueda_fluctuation_efficiency_bridge_verified"])
        self.assertTrue(self.report["generalized_fluctuation_relation_verified"])
        self.assertTrue(self.report["feedback_efficacy_stratification_verified"])
        self.assertTrue(self.report["proof_dag_information_parity_verified"])
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_gating_profiles_and_feedback_efficacy_stratification(self) -> None:
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

        # WT CFTR: Generalized Sagawa-Ueda fluctuation relation <exp(-sigma - I)> = 1.0 +/- 0.05,
        # coherent feedback efficacy gamma_fb >= 0.70, I >= 1.8 kB, sigma >= 4.0 kB/s
        self.assertTrue(wt["generalized_fluctuation_relation_satisfied"])
        self.assertEqual(wt["gating_regime"], "COHERENT_FEEDBACK_COOLING")
        self.assertAlmostEqual(wt["fluctuation_integral_exp_neg_sigma_I"], 1.0, delta=0.05)
        self.assertGreaterEqual(wt["feedback_efficacy_gamma_fb"], 0.70)
        self.assertGreaterEqual(wt["mutual_information_I_kb"], 1.8)
        self.assertGreaterEqual(wt["entropy_production_sigma_kb_per_s"], 4.0)
        self.assertGreaterEqual(wt["allosteric_controller_fidelity"], 0.90)

        # Rescued CFTR (Ivacaftor / ETI): Restores generalized fluctuation balance and feedback efficacy matching WT
        self.assertTrue(rescued["generalized_fluctuation_relation_satisfied"])
        self.assertEqual(rescued["gating_regime"], "COHERENT_FEEDBACK_COOLING")
        self.assertAlmostEqual(rescued["fluctuation_integral_exp_neg_sigma_I"], 1.0, delta=0.05)
        self.assertGreaterEqual(rescued["feedback_efficacy_gamma_fb"], 0.70)
        self.assertGreaterEqual(rescued["mutual_information_I_kb"], 1.8)
        self.assertGreaterEqual(rescued["entropy_production_sigma_kb_per_s"], 3.8)
        self.assertGreaterEqual(rescued["allosteric_controller_fidelity"], 0.90)

        # F508del: Degenerate information loss (broken feedback coupling gamma_fb -> 0, futile dissipation)
        self.assertFalse(f508del["generalized_fluctuation_relation_satisfied"])
        self.assertEqual(f508del["gating_regime"], "DEGENERATE_INFORMATION_LOSS")
        self.assertLess(f508del["feedback_efficacy_gamma_fb"], 0.10)
        self.assertLess(f508del["mutual_information_I_kb"], 0.10)
        self.assertLess(f508del["allosteric_controller_fidelity"], 0.10)

        # G551D: Degenerate information loss (uncoupled futile turnover without gate control)
        self.assertFalse(g551d["generalized_fluctuation_relation_satisfied"])
        self.assertEqual(g551d["gating_regime"], "DEGENERATE_INFORMATION_LOSS")
        self.assertLess(g551d["feedback_efficacy_gamma_fb"], 0.10)
        self.assertLess(g551d["mutual_information_I_kb"], 0.10)
        self.assertLess(g551d["allosteric_controller_fidelity"], 0.10)

        # W1282X: Passive equilibrium baseline (zero information, zero dissipation, zero work extraction)
        self.assertFalse(w1282x["generalized_fluctuation_relation_satisfied"])
        self.assertEqual(w1282x["gating_regime"], "PASSIVE_EQUILIBRIUM_APO")
        self.assertEqual(w1282x["feedback_efficacy_gamma_fb"], 0.0)
        self.assertEqual(w1282x["mutual_information_I_kb"], 0.0)
        self.assertEqual(w1282x["entropy_production_sigma_kb_per_s"], 0.0)
        self.assertEqual(w1282x["extracted_work_W_kb"], 0.0)
        self.assertEqual(w1282x["allosteric_controller_fidelity"], 0.0)

    def test_formal_proof_dag_information_parity(self) -> None:
        challenges = self.report["challenge_profiles"]
        expected_benchmarks = [
            "anthropics/fermats-last-theorem",
            "anthropics/formal-math/zeta23",
            "anthropics/formal-math/3sum-apsp",
        ]
        for c_id in expected_benchmarks:
            self.assertIn(c_id, challenges)
            c_prof = challenges[c_id]
            self.assertTrue(c_prof["sagawa_ueda_parity_satisfied"])
            self.assertAlmostEqual(c_prof["fluctuation_integral_parity"], 1.00, delta=0.05)
            self.assertGreaterEqual(c_prof["feedback_pruning_efficacy"], 0.80)
            self.assertGreaterEqual(c_prof["mutual_information_bits_per_step"], 4.0)

    def test_ah10_negative_controls(self) -> None:
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_invalid_regimes_count_raises_value_error(self) -> None:
        with self.assertRaises(ValueError):
            evaluate_sagawa_ueda_fluctuation_efficiency(n_gating_regimes=3)


if __name__ == "__main__":
    unittest.main()
