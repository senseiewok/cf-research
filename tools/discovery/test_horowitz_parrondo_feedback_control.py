"""Unit tests for Round 51 (T-0606): Horowitz-Parrondo Information Flow & Feedback Control in CFTR Catalytic Assemblies in cf-research.

Primary Biophysical Grounding & Evidence:
  - DOI:10.1038/nphys3230 (Parrondo JMR, Horowitz JM, Sagawa T: Thermodynamics of information, Nat. Phys. 11, 131–139, 2015)
  - DOI:10.1088/0034-4885/75/12/126001 (Seifert U: Stochastic thermodynamics, fluctuation theorems and molecular machines, Rep. Prog. Phys. 75, 126001, 2012)
  - DOI:10.1085/jgp.200609558 (Csanády et al.: Strict coupling between CFTR channel opening and ATP hydrolysis)
  - DOI:10.1126/science.ade2216 (Fiedorczuk & Chen: Molecular structures reveal synergistic rescue of Δ508 CFTR by Trikafta modulators, Science 2022)
  - PDB 6MSM (Open phosphorylated human CFTR)
  - PDB 5UAK (Dephosphorylated human CFTR)
  - PDB 8EIQ (Rescued human CFTR bound to ETI)
  - Lean 4 Formal Math Proof Benchmarks: fermats_last_theorem, formal_math_zeta23, formal_math_3sum_apsp
"""

from __future__ import annotations

from pathlib import Path
import sys
import unittest

_REPO_ROOT = Path(__file__).resolve().parent.parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from tools.discovery.formal_proof_thermodynamic_landscape import (
    evaluate_horowitz_parrondo_feedback_control,
)


class HorowitzParrondoFeedbackControlTests(unittest.TestCase):
    """Test suite verifying Round 51 Horowitz-Parrondo information flow and feedback control bounds for CFTR gating."""

    def setUp(self) -> None:
        self.report = evaluate_horowitz_parrondo_feedback_control()

    def test_task_id_and_overall_bridge_verification(self) -> None:
        self.assertEqual(self.report["task_id"], "T-0606")
        self.assertTrue(self.report["horowitz_parrondo_feedback_control_bridge_verified"])
        self.assertTrue(self.report["information_flow_bounds_verified"])
        self.assertTrue(self.report["feedback_fidelity_stratification_verified"])
        self.assertTrue(self.report["proof_dag_feedback_parity_verified"])
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_information_flow_bounds_and_feedback_stratification(self) -> None:
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

        # WT CFTR: Horowitz-Parrondo feedback bound satisfied (W_dot_ext <= -Delta F_dot + kB T I_dot_F)
        self.assertTrue(wt["horowitz_parrondo_bound_satisfied"])
        self.assertEqual(wt["gating_regime"], "COHERENT_FEEDBACK_CONTROL")
        self.assertGreaterEqual(wt["information_flow_rate_i_dot_f_kb_per_s"], 1.8)
        self.assertLessEqual(
            wt["extractable_work_rate_w_dot_ext_kbt_per_s"],
            -wt["free_energy_rate_delta_f_dot_kbt_per_s"] + wt["information_flow_rate_i_dot_f_kb_per_s"],
        )
        self.assertGreaterEqual(wt["allosteric_controller_fidelity"], 0.90)

        # Rescued CFTR (Ivacaftor / ETI): Restores feedback information flow matching WT
        self.assertTrue(rescued["horowitz_parrondo_bound_satisfied"])
        self.assertEqual(rescued["gating_regime"], "COHERENT_FEEDBACK_CONTROL")
        self.assertGreaterEqual(rescued["information_flow_rate_i_dot_f_kb_per_s"], 1.8)
        self.assertLessEqual(
            rescued["extractable_work_rate_w_dot_ext_kbt_per_s"],
            -rescued["free_energy_rate_delta_f_dot_kbt_per_s"] + rescued["information_flow_rate_i_dot_f_kb_per_s"],
        )
        self.assertGreaterEqual(rescued["allosteric_controller_fidelity"], 0.90)

        # F508del: Degenerate information loss (I_dot_F -> 0, uncoordinated dissipation)
        self.assertFalse(f508del["horowitz_parrondo_bound_satisfied"])
        self.assertEqual(f508del["gating_regime"], "DEGENERATE_FEEDBACK_COLLAPSE")
        self.assertLess(f508del["information_flow_rate_i_dot_f_kb_per_s"], 0.15)
        self.assertLess(f508del["extractable_work_rate_w_dot_ext_kbt_per_s"], 0.30)
        self.assertLess(f508del["allosteric_controller_fidelity"], 0.15)

        # G551D: Degenerate information loss (I_dot_F -> 0, uncoordinated catalytic turnover)
        self.assertFalse(g551d["horowitz_parrondo_bound_satisfied"])
        self.assertEqual(g551d["gating_regime"], "DEGENERATE_FEEDBACK_COLLAPSE")
        self.assertLess(g551d["information_flow_rate_i_dot_f_kb_per_s"], 0.15)
        self.assertLess(g551d["extractable_work_rate_w_dot_ext_kbt_per_s"], 0.30)
        self.assertLess(g551d["allosteric_controller_fidelity"], 0.15)

        # W1282X: Passive equilibrium baseline (W_dot_ext = Delta F_dot = I_dot_F = 0)
        self.assertFalse(w1282x["horowitz_parrondo_bound_satisfied"])
        self.assertEqual(w1282x["gating_regime"], "PASSIVE_EQUILIBRIUM_APO")
        self.assertEqual(w1282x["extractable_work_rate_w_dot_ext_kbt_per_s"], 0.0)
        self.assertEqual(w1282x["free_energy_rate_delta_f_dot_kbt_per_s"], 0.0)
        self.assertEqual(w1282x["information_flow_rate_i_dot_f_kb_per_s"], 0.0)
        self.assertEqual(w1282x["allosteric_controller_fidelity"], 0.0)

    def test_formal_proof_dag_feedback_parity(self) -> None:
        challenges = self.report["challenge_profiles"]
        expected_benchmarks = [
            "anthropics/fermats-last-theorem",
            "anthropics/formal-math/zeta23",
            "anthropics/formal-math/3sum-apsp",
        ]
        for c_id in expected_benchmarks:
            self.assertIn(c_id, challenges)
            c_prof = challenges[c_id]
            self.assertTrue(c_prof["proof_feedback_parity_satisfied"])
            self.assertGreaterEqual(c_prof["proof_controller_rate"], 3.5)
            self.assertGreaterEqual(c_prof["proof_information_flow"], 2.0)

    def test_ah10_negative_controls(self) -> None:
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_invalid_regimes_count_raises_value_error(self) -> None:
        with self.assertRaises(ValueError):
            evaluate_horowitz_parrondo_feedback_control(n_gating_regimes=3)


if __name__ == "__main__":
    unittest.main()
