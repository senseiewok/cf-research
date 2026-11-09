#!/usr/bin/env python3
"""Unit tests for Round 145: Hill-Langmuir Cooperative ATP-Binding & NBD Dimerization Invariant Bridge (T-0796)."""
from __future__ import annotations

import unittest

from tools.discovery.round145_hill_langmuir_cftr_nbd_dimerization_bridge import (
    FIVE_LAYER_LADDER,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    evaluate_falsification_negative_control,
    evaluate_round145_hill_langmuir_cftr_nbd_dimerization,
)


class Round145HillLangmuirCftrNbdDimerizationTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round145_hill_langmuir_cftr_nbd_dimerization(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Failed for {key}")
            self.assertGreaterEqual(res.open_probability, 0.0)
            self.assertLessEqual(res.open_probability, 1.0)
            self.assertGreater(res.k_close_s, 0.0)
            self.assertGreaterEqual(res.k_open_s, 0.0)
            self.assertGreaterEqual(res.epithelial_short_circuit_current_ua_cm2, 0.0)
            self.assertGreater(res.complex_step_atp_sensitivity, 0.0)

    def test_wild_type_cftr_atp_gating(self) -> None:
        wt = evaluate_round145_hill_langmuir_cftr_nbd_dimerization("wild_type_cftr_atp_gating")
        self.assertTrue(wt.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(wt.open_probability, 0.38)
        self.assertGreaterEqual(wt.fractional_saturation, 0.90)
        self.assertGreaterEqual(wt.epithelial_short_circuit_current_ua_cm2, 20.0)
        self.assertAlmostEqual(wt.fractional_wild_type_current, 1.0, delta=0.05)

    def test_cf_gating_mutant_dimer_collapse(self) -> None:
        g551d = evaluate_round145_hill_langmuir_cftr_nbd_dimerization("cf_g551d_catalytic_dimer_defect")
        f508 = evaluate_round145_hill_langmuir_cftr_nbd_dimerization("cf_delta_f508_nbd1_instability")
        for mutant in (g551d, f508):
            self.assertFalse(mutant.crosses_kahn_kalai_threshold)
            self.assertLessEqual(mutant.fractional_wild_type_current, 0.10)
        self.assertLessEqual(g551d.open_probability, 0.02)
        self.assertLessEqual(g551d.epithelial_short_circuit_current_ua_cm2, 0.50)

    def test_kahn_kalai_dimerization_boundary(self) -> None:
        kk = evaluate_round145_hill_langmuir_cftr_nbd_dimerization("kahn_kalai_dimerization_boundary")
        self.assertTrue(kk.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(kk.open_probability, 0.050)
        self.assertGreaterEqual(kk.epithelial_short_circuit_current_ua_cm2, 2.50)
        self.assertGreaterEqual(kk.fractional_wild_type_current, 0.118)

    def test_trikafta_ivacaftor_dimer_rescue(self) -> None:
        rescued = evaluate_round145_hill_langmuir_cftr_nbd_dimerization("trikafta_ivacaftor_dimer_rescue")
        self.assertTrue(rescued.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(rescued.open_probability, 0.38)
        self.assertGreaterEqual(rescued.epithelial_short_circuit_current_ua_cm2, 18.0)
        self.assertGreaterEqual(rescued.fractional_wild_type_current, 0.85)

    def test_negative_control_falsification(self) -> None:
        normal = evaluate_falsification_negative_control(poison_negative_atp=False)
        self.assertTrue(normal["passed"])
        self.assertEqual(normal["rule"], "AH-10")
        poisoned = evaluate_falsification_negative_control(poison_negative_atp=True)
        self.assertFalse(poisoned["passed"])

    def test_two_column_epistemic_hygiene_and_five_layer_ladder(self) -> None:
        self.assertIn("deductive_theorem", TWO_COLUMN_EPISTEMIC_HYGIENE)
        self.assertIn("biophysical_hypothesis", TWO_COLUMN_EPISTEMIC_HYGIENE)
        self.assertIn("falsification_condition", TWO_COLUMN_EPISTEMIC_HYGIENE["biophysical_hypothesis"])
        for layer in (
            "cool_human_fact",
            "everyday_analogy",
            "engineering_lens",
            "formal_theorem",
            "falsifiable_code_and_toy",
        ):
            self.assertTrue(len(FIVE_LAYER_LADDER.get(layer, "")) > 40)
        self.assertGreaterEqual(len(PRIMARY_CITATIONS), 5)


if __name__ == "__main__":
    unittest.main()
