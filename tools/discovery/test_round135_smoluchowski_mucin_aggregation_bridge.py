#!/usr/bin/env python3
"""Unit tests for Round 135: Smoluchowski Coagulation Equation & Mucin Oligomer Chain Aggregation Kernel Bridge."""
from __future__ import annotations

import unittest

from tools.discovery.round135_smoluchowski_mucin_aggregation_bridge import (
    FIVE_LAYER_LADDER,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    evaluate_falsification_negative_control,
    evaluate_round135_smoluchowski_mucin_aggregation,
)


class Round135SmoluchowskiMucinAggregationTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round135_smoluchowski_mucin_aggregation(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Failed for {key}")
            self.assertLess(res.complex_step_sensitivity, 0.0)
            self.assertGreater(res.phi_mucin, 0.0)
            self.assertGreater(res.gelation_time_hours, 0.0)
            self.assertGreaterEqual(res.clearance_velocity_um_s, 0.0)
            self.assertAlmostEqual(res.gel_fraction + res.sol_fraction, 1.0, places=4)

    def test_healthy_hydrated_asl_stability(self) -> None:
        healthy = evaluate_round135_smoluchowski_mucin_aggregation("healthy_hydrated_asl")
        self.assertTrue(healthy.crosses_kahn_kalai_threshold)
        self.assertEqual(healthy.gel_fraction, 0.0)
        self.assertEqual(healthy.sol_fraction, 1.0)
        self.assertGreater(healthy.gelation_time_hours, 100.0)
        self.assertGreaterEqual(healthy.clearance_velocity_um_s, 50.0)
        self.assertLessEqual(healthy.phi_mucin, 0.020)

    def test_cf_dehydrated_hyperconcentration_and_rapid_gelation(self) -> None:
        cf = evaluate_round135_smoluchowski_mucin_aggregation("cf_dehydrated_hyperconcentrated")
        self.assertFalse(cf.crosses_kahn_kalai_threshold)
        self.assertGreater(cf.gel_fraction, 0.70)
        self.assertLess(cf.sol_fraction, 0.30)
        self.assertLess(cf.gelation_time_hours, 1.0)
        self.assertLess(cf.clearance_velocity_um_s, 2.0)
        self.assertGreaterEqual(cf.phi_mucin, 0.080)

    def test_kahn_kalai_threshold_rescue(self) -> None:
        kk = evaluate_round135_smoluchowski_mucin_aggregation("kahn_kalai_threshold_11_8pct")
        self.assertTrue(kk.crosses_kahn_kalai_threshold)
        self.assertLessEqual(kk.gel_fraction, 0.15)
        self.assertGreaterEqual(kk.clearance_velocity_um_s, 5.0)
        self.assertLessEqual(kk.phi_mucin, 0.050)

    def test_rehydration_and_trikafta_rescue(self) -> None:
        hs = evaluate_round135_smoluchowski_mucin_aggregation("hypertonic_saline_rehydrated")
        tri = evaluate_round135_smoluchowski_mucin_aggregation("trikafta_full_mucociliary_rescue")
        for res in (hs, tri):
            self.assertTrue(res.crosses_kahn_kalai_threshold)
            self.assertEqual(res.gel_fraction, 0.0)
            self.assertEqual(res.sol_fraction, 1.0)
            self.assertGreaterEqual(res.gelation_time_hours, 30.0)
            self.assertGreaterEqual(res.clearance_velocity_um_s, 20.0)
            self.assertLessEqual(res.phi_mucin, 0.035)

    def test_negative_control_falsification(self) -> None:
        normal = evaluate_falsification_negative_control(poison_negative_kernel=False)
        self.assertTrue(normal["passed"])
        self.assertEqual(normal["rule"], "AH-10")
        poisoned = evaluate_falsification_negative_control(poison_negative_kernel=True)
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
