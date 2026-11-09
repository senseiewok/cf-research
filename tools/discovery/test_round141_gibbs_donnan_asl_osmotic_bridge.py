#!/usr/bin/env python3
"""Unit tests for Round 141: Gibbs-Donnan Colloid Osmotic Pressure & Asymmetric Counterion Partitioning Invariant Bridge (T-0788)."""
from __future__ import annotations

import unittest

from tools.discovery.round141_gibbs_donnan_asl_osmotic_bridge import (
    FIVE_LAYER_LADDER,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    evaluate_falsification_negative_control,
    evaluate_round141_gibbs_donnan_asl_osmotic,
)


class Round141GibbsDonnanAslOsmoticTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round141_gibbs_donnan_asl_osmotic(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Failed for {key}")
            self.assertGreater(res.complex_step_osmotic_sensitivity, 0.0)
            self.assertGreaterEqual(res.donnan_partitioning_ratio_r, 1.0)
            self.assertLessEqual(res.donnan_potential_mv, 0.0)
            self.assertGreater(res.total_mucus_osmotic_pressure_pa, 0.0)
            self.assertGreater(res.pcl_height_um, 0.0)
            self.assertLessEqual(res.pcl_height_um, 7.0)
            self.assertGreaterEqual(res.mcc_velocity_um_s, 0.0)

    def test_healthy_hydrated_asl(self) -> None:
        healthy = evaluate_round141_gibbs_donnan_asl_osmotic("healthy_hydrated_asl")
        self.assertTrue(healthy.crosses_kahn_kalai_threshold)
        self.assertEqual(healthy.mucin_solids_pct, 2.0)
        self.assertLessEqual(healthy.total_mucus_osmotic_pressure_pa, 250.0)
        self.assertGreaterEqual(healthy.pcl_height_um, 6.5)
        self.assertGreaterEqual(healthy.mcc_velocity_um_s, 40.0)

    def test_cf_dehydrated_asl_collapse(self) -> None:
        cf = evaluate_round141_gibbs_donnan_asl_osmotic("cf_dehydrated_asl_collapse")
        self.assertFalse(cf.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(cf.mucin_solids_pct, 9.0)
        self.assertGreaterEqual(cf.total_mucus_osmotic_pressure_pa, 3000.0)
        self.assertLessEqual(cf.pcl_height_um, 2.80)
        self.assertLessEqual(cf.mcc_velocity_um_s, 0.5)

    def test_kahn_kalai_asl_hydration_boundary(self) -> None:
        kk = evaluate_round141_gibbs_donnan_asl_osmotic("kahn_kalai_asl_hydration_boundary")
        self.assertTrue(kk.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(kk.pcl_height_um, 3.50)
        self.assertGreaterEqual(kk.mcc_velocity_um_s, 5.90)

    def test_therapeutic_rescue_regimes(self) -> None:
        hypertonic = evaluate_round141_gibbs_donnan_asl_osmotic("hypertonic_saline_7pct_rehydration")
        trikafta = evaluate_round141_gibbs_donnan_asl_osmotic("trikafta_cftr_hydration_rescue")
        for res in (hypertonic, trikafta):
            self.assertTrue(res.crosses_kahn_kalai_threshold)
            self.assertLessEqual(res.total_mucus_osmotic_pressure_pa, 450.0)
            self.assertGreaterEqual(res.pcl_height_um, 6.0)
            self.assertGreaterEqual(res.mcc_velocity_um_s, 30.0)

    def test_negative_control_falsification(self) -> None:
        normal = evaluate_falsification_negative_control(poison_negative_fixed_charge=False)
        self.assertTrue(normal["passed"])
        self.assertEqual(normal["rule"], "AH-10")
        poisoned = evaluate_falsification_negative_control(poison_negative_fixed_charge=True)
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
