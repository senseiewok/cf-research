#!/usr/bin/env python3
"""Unit tests for Round 140: Lifson-Roig Helix-Coil Transition & Statistical-Mechanical Nucleation Invariant Bridge."""
from __future__ import annotations

import unittest

from tools.discovery.round140_lifson_roig_cftr_folding_bridge import (
    FIVE_LAYER_LADDER,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    evaluate_falsification_negative_control,
    evaluate_round140_lifson_roig_cftr_folding,
)


class Round140LifsonRoigCftrFoldingTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round140_lifson_roig_cftr_folding(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Failed for {key}")
            self.assertGreater(res.complex_step_folding_sensitivity, 0.0)
            self.assertGreater(res.melting_temperature_c, 0.0)
            self.assertGreater(res.propagation_parameter_w, 0.0)
            self.assertGreater(res.nucleation_parameter_v, 0.0)
            self.assertGreaterEqual(res.maximum_eigenvalue_lambda, 1.0)
            self.assertGreaterEqual(res.helical_folded_fraction, 0.0)
            self.assertLessEqual(res.helical_folded_fraction, 1.0)

    def test_wild_type_nbd1_native_folding(self) -> None:
        wt = evaluate_round140_lifson_roig_cftr_folding("wild_type_nbd1_native")
        self.assertTrue(wt.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(wt.melting_temperature_c, 40.0)
        self.assertGreater(wt.propagation_parameter_w, 1.15)
        self.assertGreaterEqual(wt.helical_folded_fraction, 0.88)

    def test_cf_f508del_nbd1_misfolding(self) -> None:
        cf = evaluate_round140_lifson_roig_cftr_folding("cf_f508del_nbd1_misfolded")
        self.assertFalse(cf.crosses_kahn_kalai_threshold)
        self.assertLessEqual(cf.melting_temperature_c, 35.5)
        self.assertLess(cf.propagation_parameter_w, 0.90)
        self.assertLessEqual(cf.helical_folded_fraction, 0.08)

    def test_kahn_kalai_stability_threshold(self) -> None:
        kk = evaluate_round140_lifson_roig_cftr_folding("kahn_kalai_stability_threshold")
        self.assertTrue(kk.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(kk.helical_folded_fraction, 0.118)

    def test_tezacaftor_and_trikafta_stabilization(self) -> None:
        teza = evaluate_round140_lifson_roig_cftr_folding("tezacaftor_vx661_rescued")
        tri = evaluate_round140_lifson_roig_cftr_folding("trikafta_synergistic_stabilization")
        for res in (teza, tri):
            self.assertTrue(res.crosses_kahn_kalai_threshold)
            self.assertGreater(res.propagation_parameter_w, 1.05)
            self.assertGreaterEqual(res.helical_folded_fraction, 0.75)

    def test_negative_control_falsification(self) -> None:
        normal = evaluate_falsification_negative_control(poison_negative_nucleation=False)
        self.assertTrue(normal["passed"])
        self.assertEqual(normal["rule"], "AH-10")
        poisoned = evaluate_falsification_negative_control(poison_negative_nucleation=True)
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
