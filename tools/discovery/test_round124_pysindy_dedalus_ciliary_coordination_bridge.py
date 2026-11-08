#!/usr/bin/env python3
"""Unit tests for Round 124: PySINDy-Dedalus Weak-Form Sparse Continuum Dynamics & Ciliary Metachronal Wave Coordination Bridge."""
from __future__ import annotations

import unittest

from tools.discovery.round124_pysindy_dedalus_ciliary_coordination_bridge import (
    FIVE_LAYER_LADDER,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    evaluate_falsification_negative_control,
    evaluate_round124_pysindy_dedalus_ciliary_coordination,
)


class Round124PySINDyDedalusCiliaryBridgeTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round124_pysindy_dedalus_ciliary_coordination(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Failed for {key}")
            self.assertLess(res.complex_step_sensitivity, 0.0)

    def test_wild_type_healthy_crosses_jamming_safety_threshold(self) -> None:
        wt = evaluate_round124_pysindy_dedalus_ciliary_coordination("wild_type_healthy")
        self.assertTrue(wt.crosses_jamming_safety_threshold)
        self.assertGreaterEqual(wt.pcl_height_um, 7.0)
        self.assertGreaterEqual(wt.kuramoto_order_parameter, 0.85)
        self.assertGreaterEqual(wt.mcc_velocity_um_s, 25.0)

    def test_f508del_untreated_collapses_wave_and_arrests_transport(self) -> None:
        cf = evaluate_round124_pysindy_dedalus_ciliary_coordination("f508del_homozygous_untreated")
        self.assertFalse(cf.crosses_jamming_safety_threshold)
        self.assertLess(cf.pcl_height_um, 4.0)
        self.assertLess(cf.kuramoto_order_parameter, 0.35)
        self.assertLess(cf.mcc_velocity_um_s, 5.0)

    def test_kahn_kalai_threshold_rescues_coordination_and_clearance(self) -> None:
        kk = evaluate_round124_pysindy_dedalus_ciliary_coordination("kahn_kalai_threshold_11_8pct")
        self.assertAlmostEqual(kk.cftr_function_pct, 11.8, delta=0.5)
        self.assertGreaterEqual(kk.kuramoto_order_parameter, 0.65)
        self.assertGreaterEqual(kk.mcc_velocity_um_s, 10.0)

    def test_trikafta_rescued_restores_healthy_coordination(self) -> None:
        tri = evaluate_round124_pysindy_dedalus_ciliary_coordination("trikafta_rescued")
        self.assertTrue(tri.crosses_jamming_safety_threshold)
        self.assertGreaterEqual(tri.pcl_height_um, 7.0)
        self.assertGreaterEqual(tri.kuramoto_order_parameter, 0.80)
        self.assertGreaterEqual(tri.mcc_velocity_um_s, 25.0)

    def test_negative_control_falsification(self) -> None:
        normal = evaluate_falsification_negative_control(poison_break_wave=False)
        self.assertTrue(normal["passed"])
        self.assertEqual(normal["rule"], "AH-10")
        poisoned = evaluate_falsification_negative_control(poison_break_wave=True)
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
