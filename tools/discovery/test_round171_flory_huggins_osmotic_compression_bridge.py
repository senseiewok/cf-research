#!/usr/bin/env python3
"""Unit tests for Round 171: Flory-Huggins Polymer Osmotic Compression & Ciliary Entrapment Collapse Invariant Bridge (T-0836)."""
from __future__ import annotations

import unittest

from tools.discovery.round171_flory_huggins_osmotic_compression_bridge import (
    FIVE_LAYER_LADDER,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    complex_step_mixing_pressure_derivative,
    compute_ciliary_clearance_velocity,
    compute_ciliary_entrapment_fraction,
    compute_flory_huggins_mixing_pressure,
    compute_periciliary_brush_height,
    compute_total_mucin_osmotic_pressure,
    evaluate_falsification_negative_control,
    evaluate_round171_flory_huggins_compression,
)


class Round171FloryHugginsTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round171_flory_huggins_compression(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Failed for {key}")
            self.assertGreater(res.total_mucin_osmotic_pressure_pa, 0.0)
            self.assertGreaterEqual(res.pcl_layer_height_um, 1.80)
            self.assertLessEqual(res.pcl_layer_height_um, 7.00)
            self.assertGreaterEqual(res.ciliary_clearance_velocity_um_s, 0.0)
            self.assertLessEqual(res.ciliary_clearance_velocity_um_s, 55.0)
            self.assertGreaterEqual(res.ciliary_entrapment_fraction, 0.0)
            self.assertLessEqual(res.ciliary_entrapment_fraction, 1.0)
            self.assertLessEqual(res.complex_step_relative_error, 1e-10)

    def test_healthy_euhydrated_mucus(self) -> None:
        res = evaluate_round171_flory_huggins_compression("healthy_euhydrated_mucus")
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertLessEqual(res.total_mucin_osmotic_pressure_pa, 200.0)
        self.assertGreaterEqual(res.pcl_layer_height_um, 6.50)
        self.assertGreaterEqual(res.ciliary_clearance_velocity_um_s, 50.0)
        self.assertLessEqual(res.ciliary_entrapment_fraction, 0.10)

    def test_cf_dehydrated_osmotic_collapse(self) -> None:
        res = evaluate_round171_flory_huggins_compression("cf_dehydrated_osmotic_collapse")
        self.assertFalse(res.crosses_kahn_kalai_threshold)
        self.assertGreater(res.total_mucin_osmotic_pressure_pa, 1500.0)
        self.assertLess(res.pcl_layer_height_um, 3.00)
        self.assertLess(res.ciliary_clearance_velocity_um_s, 10.0)
        self.assertGreater(res.ciliary_entrapment_fraction, 0.75)

    def test_cf_extreme_mucoid_plugging(self) -> None:
        res = evaluate_round171_flory_huggins_compression("cf_extreme_mucoid_plugging")
        self.assertFalse(res.crosses_kahn_kalai_threshold)
        self.assertGreater(res.total_mucin_osmotic_pressure_pa, 2500.0)
        self.assertLess(res.pcl_layer_height_um, 2.00)
        self.assertLess(res.ciliary_clearance_velocity_um_s, 5.0)
        self.assertGreater(res.ciliary_entrapment_fraction, 0.90)

    def test_hypertonic_saline_7pct_rehydration(self) -> None:
        res = evaluate_round171_flory_huggins_compression("hypertonic_saline_7pct_rehydration")
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertLessEqual(res.total_mucin_osmotic_pressure_pa, 350.0)
        self.assertGreaterEqual(res.pcl_layer_height_um, 5.50)
        self.assertGreaterEqual(res.ciliary_clearance_velocity_um_s, 50.0)
        self.assertLessEqual(res.ciliary_entrapment_fraction, 0.15)

    def test_cf_intermediate_dehydration(self) -> None:
        res = evaluate_round171_flory_huggins_compression("cf_intermediate_dehydration")
        self.assertFalse(res.crosses_kahn_kalai_threshold)
        self.assertGreater(res.total_mucin_osmotic_pressure_pa, 600.0)
        self.assertLess(res.pcl_layer_height_um, 4.50)
        self.assertGreater(res.ciliary_entrapment_fraction, 0.30)

    def test_complex_step_mixing_derivative(self) -> None:
        cs_deriv, analytical, rel_err = complex_step_mixing_pressure_derivative(
            phi_muc=0.040,
        )
        self.assertLessEqual(rel_err, 1e-10)
        self.assertGreater(cs_deriv, 0.0)  # Osmotic mixing pressure increases with volume fraction

    def test_negative_control_falsification(self) -> None:
        neg1 = evaluate_falsification_negative_control(poison_negative_volume_fraction=True)
        self.assertFalse(neg1["passed"])
        self.assertIn("poisoned", neg1["error"])

        neg2 = evaluate_falsification_negative_control(poison_excessive_volume_fraction=True)
        self.assertFalse(neg2["passed"])

        neg3 = evaluate_falsification_negative_control(poison_negative_temperature=True)
        self.assertFalse(neg3["passed"])

        neg4 = evaluate_falsification_negative_control(poison_negative_pcl_pressure=True)
        self.assertFalse(neg4["passed"])

        clean = evaluate_falsification_negative_control()
        self.assertTrue(clean["passed"])

    def test_epistemic_hygiene_and_ladder_structure(self) -> None:
        self.assertIn("deductive_theorem", TWO_COLUMN_EPISTEMIC_HYGIENE)
        self.assertIn("biophysical_hypothesis", TWO_COLUMN_EPISTEMIC_HYGIENE)
        self.assertGreaterEqual(len(FIVE_LAYER_LADDER), 5)
        self.assertGreaterEqual(len(PRIMARY_CITATIONS), 4)


if __name__ == "__main__":
    unittest.main()
