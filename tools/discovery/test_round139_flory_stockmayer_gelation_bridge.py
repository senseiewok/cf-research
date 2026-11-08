#!/usr/bin/env python3
"""Unit tests for Round 139: Flory-Stockmayer Gelation Percolation & Branching Function Invariant Bridge."""
from __future__ import annotations

import unittest

from tools.discovery.round139_flory_stockmayer_gelation_bridge import (
    FIVE_LAYER_LADDER,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    evaluate_falsification_negative_control,
    evaluate_round139_flory_stockmayer_gelation,
)


class Round139FloryStockmayerGelationTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round139_flory_stockmayer_gelation(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Failed for {key}")
            self.assertGreater(res.complex_step_gel_sensitivity, 0.0)
            self.assertGreaterEqual(res.branch_functionality_f, 3.0)
            self.assertGreater(res.critical_threshold_alpha_c, 0.0)
            self.assertGreaterEqual(res.sol_weight_fraction, 0.0)
            self.assertLessEqual(res.sol_weight_fraction, 1.0)
            self.assertGreaterEqual(res.gel_weight_fraction, 0.0)
            self.assertLess(res.gel_weight_fraction, 1.0)
            self.assertGreaterEqual(res.mcc_velocity_um_s, 0.0)

    def test_healthy_subgel_asl(self) -> None:
        asl = evaluate_round139_flory_stockmayer_gelation("healthy_subgel_asl")
        self.assertTrue(asl.crosses_kahn_kalai_threshold)
        self.assertFalse(asl.is_post_gel_network)
        self.assertEqual(asl.gel_weight_fraction, 0.0)
        self.assertEqual(asl.sol_weight_fraction, 1.0)
        self.assertGreaterEqual(asl.mcc_velocity_um_s, 30.0)

    def test_cf_hypercrosslinked_gel(self) -> None:
        cf = evaluate_round139_flory_stockmayer_gelation("cf_hypercrosslinked_gel")
        self.assertFalse(cf.crosses_kahn_kalai_threshold)
        self.assertTrue(cf.is_post_gel_network)
        self.assertGreaterEqual(cf.gel_weight_fraction, 0.90)
        self.assertLessEqual(cf.sol_weight_fraction, 0.10)
        self.assertLessEqual(cf.mcc_velocity_um_s, 1.0)

    def test_kahn_kalai_gelation_boundary(self) -> None:
        kk = evaluate_round139_flory_stockmayer_gelation("kahn_kalai_gelation_threshold")
        self.assertTrue(kk.crosses_kahn_kalai_threshold)
        self.assertLessEqual(kk.gel_weight_fraction, 0.15)
        self.assertGreaterEqual(kk.mcc_velocity_um_s, 5.90)

    def test_dornase_alfa_and_bicarbonate_rescue(self) -> None:
        dornase = evaluate_round139_flory_stockmayer_gelation("dornase_alfa_enzymatic_cleavage")
        bicarb = evaluate_round139_flory_stockmayer_gelation("bicarbonate_rehydration_rescue")
        for res in (dornase, bicarb):
            self.assertTrue(res.crosses_kahn_kalai_threshold)
            self.assertFalse(res.is_post_gel_network)
            self.assertEqual(res.gel_weight_fraction, 0.0)
            self.assertGreaterEqual(res.mcc_velocity_um_s, 25.0)

    def test_negative_control_falsification(self) -> None:
        normal = evaluate_falsification_negative_control(poison_negative_functionality=False)
        self.assertTrue(normal["passed"])
        self.assertEqual(normal["rule"], "AH-10")
        poisoned = evaluate_falsification_negative_control(poison_negative_functionality=True)
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
