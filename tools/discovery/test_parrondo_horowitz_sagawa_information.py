"""Unit tests for Round 48 (T-0601): Parrondo-Horowitz-Sagawa Information Thermodynamics & Maxwell Demon Feedback Cooling in CFTR Channel Assemblies in cf-research.

Primary Biophysical Grounding & Evidence:
  - DOI:10.1038/nphys3230 (Parrondo JMR, Horowitz JM, Sagawa T: Thermodynamics of information, Nat. Phys. 11, 131–139, 2015)
  - DOI:10.1088/0034-4885/75/12/126001 (Seifert U: Stochastic thermodynamics, fluctuation theorems and molecular machines, Rep. Prog. Phys. 75, 126001, 2012)
  - DOI:10.1085/jgp.200609558 (Csanády et al.: Strict coupling between CFTR channel opening and ATP hydrolysis)
  - DOI:10.1126/science.ade2216 (Fiedorczuk & Chen: Cryo-EM structures of rescued human CFTR, PDB 8EIQ)
  - PDB 6MSM (Open phosphorylated human CFTR)
  - PDB 5UAK (Dephosphorylated human CFTR)
  - PDB 8EIQ (Rescued human CFTR bound to ETI)
  - Lean 4 Formal Math Proof Benchmarks: fermats_last_theorem, formal_math_zeta23, formal_math_3sum_apsp
"""

from __future__ import annotations

import unittest
from tools.discovery.formal_proof_thermodynamic_landscape import (
    evaluate_parrondo_horowitz_sagawa_information,
)


class ParrondoHorowitzSagawaInformationTests(unittest.TestCase):
    """Test suite verifying Round 48 Parrondo-Horowitz-Sagawa information thermodynamics and Maxwell demon feedback cooling for CFTR gating."""

    def setUp(self) -> None:
        self.report = evaluate_parrondo_horowitz_sagawa_information()

    def test_task_id_and_overall_bridge_verification(self) -> None:
        self.assertEqual(self.report["task_id"], "T-0601")
        self.assertTrue(self.report["parrondo_horowitz_sagawa_information_bridge_verified"])
        self.assertTrue(self.report["feedback_cooling_bounds_verified"])
        self.assertTrue(self.report["demon_efficiency_stratification_verified"])
        self.assertTrue(self.report["proof_dag_information_parity_verified"])
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_feedback_cooling_and_demon_efficiency_stratification(self) -> None:
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

        # WT CFTR: Coherent feedback cooling and information-to-work conversion (eta_info >= 0.75, I >= 1.8 kB, sigma >= 4.0 kB/s)
        self.assertTrue(wt["sagawa_ueda_second_law_satisfied"])
        self.assertEqual(wt["gating_regime"], "COHERENT_FEEDBACK_COOLING")
        self.assertGreaterEqual(wt["mutual_information_acquisition_kb"], 1.8)
        self.assertGreaterEqual(wt["feedback_cooling_rate_kb_per_s"], 1.5)
        self.assertGreaterEqual(wt["dissipation_rate_sigma_kb_per_s"], 4.0)
        self.assertGreaterEqual(wt["demon_efficiency_eta_info"], 0.75)
        self.assertGreaterEqual(wt["allosteric_controller_fidelity"], 0.90)

        # Rescued CFTR (Ivacaftor / ETI): Restores mutual information measurement and allosteric feedback cooling matching WT
        self.assertTrue(rescued["sagawa_ueda_second_law_satisfied"])
        self.assertEqual(rescued["gating_regime"], "COHERENT_FEEDBACK_COOLING")
        self.assertGreaterEqual(rescued["mutual_information_acquisition_kb"], 1.8)
        self.assertGreaterEqual(rescued["feedback_cooling_rate_kb_per_s"], 1.5)
        self.assertGreaterEqual(rescued["dissipation_rate_sigma_kb_per_s"], 3.8)
        self.assertGreaterEqual(rescued["demon_efficiency_eta_info"], 0.75)
        self.assertGreaterEqual(rescued["allosteric_controller_fidelity"], 0.90)

        # F508del: Degenerate information loss (I -> 0, uncoupled futile ATP dissipation without work extraction)
        self.assertFalse(f508del["sagawa_ueda_second_law_satisfied"])
        self.assertEqual(f508del["gating_regime"], "DEGENERATE_INFORMATION_LOSS")
        self.assertLess(f508del["mutual_information_acquisition_kb"], 0.10)
        self.assertLess(f508del["demon_efficiency_eta_info"], 0.10)
        self.assertLess(f508del["allosteric_controller_fidelity"], 0.10)

        # G551D: Degenerate information loss (uncoupled futile turnover without gate control)
        self.assertFalse(g551d["sagawa_ueda_second_law_satisfied"])
        self.assertEqual(g551d["gating_regime"], "DEGENERATE_INFORMATION_LOSS")
        self.assertLess(g551d["mutual_information_acquisition_kb"], 0.10)
        self.assertLess(g551d["demon_efficiency_eta_info"], 0.10)
        self.assertLess(g551d["allosteric_controller_fidelity"], 0.10)

        # W1282X: Passive equilibrium baseline (zero information, zero dissipation, zero work)
        self.assertFalse(w1282x["sagawa_ueda_second_law_satisfied"])
        self.assertEqual(w1282x["gating_regime"], "PASSIVE_EQUILIBRIUM_APO")
        self.assertEqual(w1282x["mutual_information_acquisition_kb"], 0.0)
        self.assertEqual(w1282x["feedback_cooling_rate_kb_per_s"], 0.0)
        self.assertEqual(w1282x["dissipation_rate_sigma_kb_per_s"], 0.0)
        self.assertEqual(w1282x["demon_efficiency_eta_info"], 0.0)
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
            self.assertGreaterEqual(c_prof["mutual_information_bits_per_step"], 4.0)
            self.assertGreaterEqual(c_prof["feedback_lemma_pruning_efficiency"], 0.80)

    def test_ah10_negative_controls(self) -> None:
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_invalid_regimes_count_raises_value_error(self) -> None:
        with self.assertRaises(ValueError):
            evaluate_parrondo_horowitz_sagawa_information(n_gating_regimes=3)


if __name__ == "__main__":
    unittest.main()
