#!/usr/bin/env python3
"""Unit tests for Round 128: Dean-Kawasaki Fluctuating Hydrodynamics & Mesoscopic Mucin Concentration Fluctuation Invariant Bridge."""
from __future__ import annotations

import unittest

from tools.discovery.round128_dean_kawasaki_mucin_fluctuations_bridge import (
    FIVE_LAYER_LADDER,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    evaluate_falsification_negative_control,
    evaluate_round128_dean_kawasaki_mucin_fluctuations,
)


class Round128DeanKawasakiMucinFluctuationsTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round128_dean_kawasaki_mucin_fluctuations(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Failed for {key}")
            self.assertGreater(res.complex_step_sensitivity, 0.0)
            self.assertGreater(res.mucus_osmotic_pressure_pa, 0.0)
            self.assertGreater(res.pcl_osmotic_pressure_pa, 0.0)
            self.assertGreater(res.pcl_brush_height_um, 0.0)
            self.assertLessEqual(res.pcl_brush_height_um, 7.0)
            self.assertGreater(res.relative_density_fluctuation_sd, 0.0)

    def test_wild_type_healthy_osmotic_balance(self) -> None:
        wt = evaluate_round128_dean_kawasaki_mucin_fluctuations("wild_type_healthy")
        self.assertTrue(wt.crosses_kahn_kalai_threshold)
        self.assertTrue(wt.ciliary_stasis_prevented)
        self.assertEqual(wt.compressive_stress_pa, 0.0)
        self.assertGreaterEqual(wt.pcl_brush_height_um, 6.8)
        self.assertLessEqual(wt.relative_density_fluctuation_sd, 0.15)

    def test_cf_homozygous_dehydrated_brush_collapse(self) -> None:
        cf = evaluate_round128_dean_kawasaki_mucin_fluctuations("cf_homozygous_dehydrated")
        self.assertFalse(cf.crosses_kahn_kalai_threshold)
        self.assertFalse(cf.ciliary_stasis_prevented)
        self.assertGreaterEqual(cf.compressive_stress_pa, 700.0)
        self.assertLess(cf.pcl_brush_height_um, 5.0)
        self.assertGreaterEqual(cf.relative_density_fluctuation_sd, 0.25)

    def test_kahn_kalai_threshold_rescue(self) -> None:
        kk = evaluate_round128_dean_kawasaki_mucin_fluctuations("kahn_kalai_threshold_11_8pct")
        self.assertTrue(kk.crosses_kahn_kalai_threshold)
        self.assertTrue(kk.ciliary_stasis_prevented)
        self.assertLessEqual(kk.compressive_stress_pa, 50.0)
        self.assertGreaterEqual(kk.pcl_brush_height_um, 5.0)

    def test_rhdnase_and_hypertonic_saline_decompression(self) -> None:
        dnase = evaluate_round128_dean_kawasaki_mucin_fluctuations("rhdnase_pulmozyme_rescued")
        saline = evaluate_round128_dean_kawasaki_mucin_fluctuations("hypertonic_saline_7pct_rehydrated")
        for res in (dnase, saline):
            self.assertTrue(res.crosses_kahn_kalai_threshold)
            self.assertTrue(res.ciliary_stasis_prevented)
            self.assertEqual(res.compressive_stress_pa, 0.0)
            self.assertGreaterEqual(res.pcl_brush_height_um, 6.0)

    def test_negative_control_falsification(self) -> None:
        normal = evaluate_falsification_negative_control(poison_negative_osmotic_pressure=False)
        self.assertTrue(normal["passed"])
        self.assertEqual(normal["rule"], "AH-10")
        poisoned = evaluate_falsification_negative_control(poison_negative_osmotic_pressure=True)
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
