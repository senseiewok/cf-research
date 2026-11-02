"""Unit tests for Round 30 (T-0574): Freidlin-Wentzell Large Deviation Theory & Quasipotential Action for Rare CFTR Opening Events.

Primary Biophysical Grounding & Evidence:
  - DOI:10.1007/978-3-642-25847-3 (Freidlin & Wentzell: Random Perturbations of Dynamical Systems)
  - DOI:10.1085/jgp.200609558 (Strict coupling between CFTR channel opening and ATP hydrolysis)
  - DOI:10.1126/science.ade2216 (Cryo-EM structures of rescued human CFTR, PDB 8EIQ)
  - PDB 6MSM (Open phosphorylated human CFTR)
  - PDB 5UAK (Dephosphorylated human CFTR)
  - Lean 4 Formal Math Proof Benchmarks: fermats_last_theorem, formal_math_zeta23, formal_math_3sum_apsp
"""

from __future__ import annotations

import unittest
from tools.discovery.formal_proof_thermodynamic_landscape import (
    evaluate_freidlin_wentzell_quasipotential,
)


class FreidlinWentzellQuasipotentialTests(unittest.TestCase):
    """Test suite verifying Round 30 Freidlin-Wentzell large deviation theory & quasipotential action in CFTR rare-event gating."""

    def setUp(self) -> None:
        self.report = evaluate_freidlin_wentzell_quasipotential()

    def test_task_id_and_overall_bridge_verification(self) -> None:
        self.assertEqual(self.report["task_id"], "T-0574")
        self.assertTrue(self.report["freidlin_wentzell_bridge_verified"])
        self.assertTrue(self.report["quasipotential_barrier_stratification_verified"])
        self.assertTrue(self.report["proof_rare_event_parity_verified"])
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_freidlin_wentzell_quasipotential_and_rare_event_stratification(self) -> None:
        profiles = self.report["gating_profiles"]
        self.assertIn("wt_cftr", profiles)
        self.assertIn("f508del_uncorrected", profiles)
        self.assertIn("g551d_uncoupled_cftr", profiles)
        self.assertIn("rescued_cftr_ivacaftor", profiles)
        self.assertIn("w1282x_nonsense_baseline", profiles)

        wt = profiles["wt_cftr"]
        f508del = profiles["f508del_uncorrected"]
        g551d = profiles["g551d_uncoupled_cftr"]
        rescued = profiles["rescued_cftr_ivacaftor"]
        w1282x = profiles["w1282x_nonsense_baseline"]

        # WT CFTR: Low quasipotential barrier (V <= 5.2 kcal/mol) under active ATP dimerization
        self.assertTrue(wt["freidlin_wentzell_criterion_satisfied"])
        self.assertEqual(wt["gating_regime"], "ACTIVE_CATALYZED_TRANSITION")
        self.assertLessEqual(wt["quasipotential_v_kcal_mol"], 5.2)
        self.assertGreaterEqual(wt["transition_rate_k_s_inv"], 1.50)
        self.assertGreaterEqual(wt["rare_event_escape_probability"], 0.85)
        self.assertGreater(wt["noise_intensity_epsilon"], 0.0)

        # F508del: High quasipotential barrier (V >= 16.0 kcal/mol), suppressing opening
        self.assertFalse(f508del["freidlin_wentzell_criterion_satisfied"])
        self.assertEqual(f508del["gating_regime"], "HIGH_QUASIPOTENTIAL_TRAP")
        self.assertGreaterEqual(f508del["quasipotential_v_kcal_mol"], 16.0)
        self.assertLess(f508del["transition_rate_k_s_inv"], 0.01)
        self.assertLess(f508del["rare_event_escape_probability"], 0.05)

        # G551D: High quasipotential barrier (V >= 16.0 kcal/mol), suppressing opening
        self.assertFalse(g551d["freidlin_wentzell_criterion_satisfied"])
        self.assertEqual(g551d["gating_regime"], "HIGH_QUASIPOTENTIAL_TRAP")
        self.assertGreaterEqual(g551d["quasipotential_v_kcal_mol"], 16.0)
        self.assertLess(g551d["transition_rate_k_s_inv"], 0.01)
        self.assertLess(g551d["rare_event_escape_probability"], 0.05)

        # Rescued CFTR (Ivacaftor / ETI): Restores minimal quasipotential barrier matching WT (V <= 5.2 kcal/mol)
        self.assertTrue(rescued["freidlin_wentzell_criterion_satisfied"])
        self.assertEqual(rescued["gating_regime"], "ACTIVE_CATALYZED_TRANSITION")
        self.assertLessEqual(rescued["quasipotential_v_kcal_mol"], 5.2)
        self.assertGreaterEqual(rescued["transition_rate_k_s_inv"], 1.50)
        self.assertGreaterEqual(rescued["rare_event_escape_probability"], 0.85)

        # W1282X: Passive equilibrium baseline
        self.assertFalse(w1282x["freidlin_wentzell_criterion_satisfied"])
        self.assertEqual(w1282x["gating_regime"], "PASSIVE_EQUILIBRIUM_APO")
        self.assertEqual(w1282x["quasipotential_v_kcal_mol"], 0.0)
        self.assertEqual(w1282x["transition_rate_k_s_inv"], 0.0)
        self.assertEqual(w1282x["rare_event_escape_probability"], 0.0)

    def test_formal_proof_dag_rare_lemma_derivation_parity(self) -> None:
        challenges = self.report["challenge_profiles"]
        expected_benchmarks = [
            "anthropics/fermats-last-theorem",
            "anthropics/formal-math/zeta23",
            "anthropics/formal-math/3sum-apsp",
        ]
        for c_id in expected_benchmarks:
            self.assertIn(c_id, challenges)
            c_prof = challenges[c_id]
            self.assertTrue(c_prof["quasipotential_proof_geodesic_satisfied"])
            self.assertTrue(c_prof["large_deviation_bound_satisfied"])
            self.assertGreaterEqual(c_prof["action_barrier_reduction_pct"], 50.0)
            self.assertLess(c_prof["quasipotential_v_dag"], 5.0)

    def test_ah10_negative_controls(self) -> None:
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_invalid_regimes_count_raises_value_error(self) -> None:
        with self.assertRaises(ValueError):
            evaluate_freidlin_wentzell_quasipotential(n_gating_regimes=2)


if __name__ == "__main__":
    unittest.main()
