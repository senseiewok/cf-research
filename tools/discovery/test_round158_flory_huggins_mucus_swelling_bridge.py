#!/usr/bin/env python3
"""Unit tests for Round 158: Flory-Huggins Gel Swelling & Osmotic Dehydration Equilibrium Invariant Bridge (T-0819)."""
from __future__ import annotations

import math
import unittest

from tools.discovery.round158_flory_huggins_mucus_swelling_bridge import (
    FIVE_LAYER_LADDER,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    complex_step_osmotic_pressure_sensitivity,
    compute_mcc_velocity_mm_min,
    compute_mucus_osmotic_pressure,
    compute_pcl_height_um,
    evaluate_falsification_negative_control,
    evaluate_round158_flory_huggins_mucus_swelling,
)


class Round158MucusSwellingTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round158_flory_huggins_mucus_swelling(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Failed for {key}")
            self.assertGreater(res.total_mucus_osmotic_pressure_pa, 0.0)
            self.assertGreaterEqual(res.pcl_height_um, 1.50)
            self.assertLessEqual(res.pcl_height_um, 7.0)
            self.assertGreaterEqual(res.mucociliary_clearance_mm_min, 0.0)
            self.assertLessEqual(res.complex_step_relative_error, 1e-11)

    def test_healthy_euhydrated_airway_mucus(self) -> None:
        healthy = evaluate_round158_flory_huggins_mucus_swelling("healthy_euhydrated_airway_mucus")
        self.assertTrue(healthy.crosses_kahn_kalai_threshold)
        self.assertFalse(healthy.pcl_collapsed)
        self.assertLess(healthy.total_mucus_osmotic_pressure_pa, 300.0)
        self.assertGreaterEqual(healthy.pcl_height_um, 6.2)
        self.assertGreaterEqual(healthy.mucociliary_clearance_mm_min, 4.8)
        self.assertAlmostEqual(healthy.fractional_wild_type_clearance, 1.0, delta=0.10)

    def test_cf_hyperconcentrated_mucus_collapse(self) -> None:
        cf = evaluate_round158_flory_huggins_mucus_swelling("cf_hyperconcentrated_mucus_collapse")
        self.assertFalse(cf.crosses_kahn_kalai_threshold)
        self.assertTrue(cf.pcl_collapsed)
        self.assertGreater(cf.total_mucus_osmotic_pressure_pa, 3000.0)
        self.assertLess(cf.pcl_height_um, 2.5)
        self.assertLess(cf.mucociliary_clearance_mm_min, 0.20)
        self.assertLess(cf.fractional_wild_type_clearance, 0.05)

    def test_cf_extreme_mucoid_plugging(self) -> None:
        plug = evaluate_round158_flory_huggins_mucus_swelling("cf_extreme_mucoid_plugging")
        self.assertFalse(plug.crosses_kahn_kalai_threshold)
        self.assertTrue(plug.pcl_collapsed)
        self.assertGreater(plug.total_mucus_osmotic_pressure_pa, 10000.0)
        self.assertLess(plug.pcl_height_um, 1.7)
        self.assertAlmostEqual(plug.mucociliary_clearance_mm_min, 0.0, delta=0.01)

    def test_kahn_kalai_swelling_equilibrium_boundary(self) -> None:
        kk = evaluate_round158_flory_huggins_mucus_swelling("kahn_kalai_swelling_equilibrium_boundary")
        self.assertTrue(kk.crosses_kahn_kalai_threshold)
        self.assertFalse(kk.pcl_collapsed)
        self.assertGreaterEqual(kk.pcl_height_um, 4.0)
        self.assertGreaterEqual(kk.mucociliary_clearance_mm_min, 0.65)
        self.assertGreaterEqual(kk.fractional_wild_type_clearance, 0.118)

    def test_hypertonic_saline_rehydration_rescue(self) -> None:
        rescued = evaluate_round158_flory_huggins_mucus_swelling("hypertonic_saline_osmolyte_rehydration_rescue")
        self.assertTrue(rescued.crosses_kahn_kalai_threshold)
        self.assertFalse(rescued.pcl_collapsed)
        self.assertLessEqual(rescued.total_mucus_osmotic_pressure_pa, 550.0)
        self.assertGreaterEqual(rescued.pcl_height_um, 5.5)
        self.assertGreaterEqual(rescued.mucociliary_clearance_mm_min, 3.5)
        self.assertGreaterEqual(rescued.fractional_wild_type_clearance, 0.65)

    def test_monotonic_osmotic_pressure_with_solids(self) -> None:
        solids_range = [0.015, 0.025, 0.040, 0.060, 0.085, 0.120]
        pressures = [compute_mucus_osmotic_pressure(phi)[2] for phi in solids_range]
        for i in range(len(pressures) - 1):
            self.assertLess(pressures[i], pressures[i + 1])

    def test_complex_step_sensitivity(self) -> None:
        deriv_cs, deriv_ana, rel_err = complex_step_osmotic_pressure_sensitivity(0.040)
        self.assertGreater(deriv_cs, 0.0)
        self.assertLessEqual(rel_err, 1e-11)

    def test_negative_control_falsification(self) -> None:
        clean = evaluate_falsification_negative_control()
        self.assertTrue(clean["passed"])
        self.assertEqual(clean["rule"], "AH-10")

        poisoned_temp = evaluate_falsification_negative_control(poison_negative_temperature=True)
        self.assertFalse(poisoned_temp["passed"])

        poisoned_phi = evaluate_falsification_negative_control(poison_unphysical_volume_fraction=True)
        self.assertFalse(poisoned_phi["passed"])

        poisoned_mod = evaluate_falsification_negative_control(poison_negative_pcl_modulus=True)
        self.assertFalse(poisoned_mod["passed"])

    def test_two_column_epistemic_hygiene_and_five_layer_ladder(self) -> None:
        self.assertIn("deductive_theorem", TWO_COLUMN_EPISTEMIC_HYGIENE)
        self.assertIn("biophysical_hypothesis", TWO_COLUMN_EPISTEMIC_HYGIENE)
        self.assertIn("falsification_condition", TWO_COLUMN_EPISTEMIC_HYGIENE["biophysical_hypothesis"])
        self.assertIn("cool_human_fact", FIVE_LAYER_LADDER)
        self.assertIn("everyday_analogy", FIVE_LAYER_LADDER)
        self.assertIn("engineering_lens", FIVE_LAYER_LADDER)
        self.assertIn("formal_theorem", FIVE_LAYER_LADDER)
        self.assertIn("falsifiable_code_and_toy", FIVE_LAYER_LADDER)
        self.assertGreaterEqual(len(PRIMARY_CITATIONS), 3)


if __name__ == "__main__":
    unittest.main()
