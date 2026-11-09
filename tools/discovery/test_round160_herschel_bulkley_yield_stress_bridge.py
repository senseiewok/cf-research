#!/usr/bin/env python3
"""Unit tests for Round 160: Herschel-Bulkley Viscoplastic Yield-Stress Invariant Bridge (T-0825)."""
from __future__ import annotations

import math
import unittest

from tools.discovery.round160_herschel_bulkley_yield_stress_bridge import (
    FIVE_LAYER_LADDER,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    complex_step_differential_viscosity,
    compute_apparent_viscosity,
    compute_ciliary_shear_rate,
    compute_extrusion_flow_rate_mm3_s,
    compute_herschel_bulkley_shear_stress,
    compute_mcc_velocity_mm_min,
    compute_poiseuille_plug_radius_mm,
    evaluate_falsification_negative_control,
    evaluate_round160_herschel_bulkley_yield_stress,
)


class Round160HerschelBulkleyYieldStressTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round160_herschel_bulkley_yield_stress(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Failed for {key}")
            self.assertGreaterEqual(res.mcc_velocity_mm_min, 0.0)
            self.assertGreaterEqual(res.plug_radius_fraction, 0.0)
            self.assertLessEqual(res.plug_radius_fraction, 1.0)
            self.assertLessEqual(res.complex_step_relative_error, 1e-10)

    def test_healthy_euhydrated_airway_mucus(self) -> None:
        res = evaluate_round160_herschel_bulkley_yield_stress("healthy_euhydrated_airway_mucus")
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertFalse(res.ciliary_stall)
        self.assertGreaterEqual(res.mcc_velocity_mm_min, 4.50)
        self.assertGreater(res.ciliary_shear_rate_s1, 1.0)
        self.assertLess(res.yield_stress_pa, 0.20)

    def test_cf_hyperconcentrated_sputum_arrest(self) -> None:
        res = evaluate_round160_herschel_bulkley_yield_stress("cf_hyperconcentrated_sputum_arrest")
        self.assertFalse(res.crosses_kahn_kalai_threshold)
        self.assertTrue(res.ciliary_stall)
        self.assertEqual(res.mcc_velocity_mm_min, 0.0)
        self.assertEqual(res.ciliary_shear_rate_s1, 0.0)
        self.assertGreaterEqual(res.yield_stress_pa, 1.50)

    def test_cf_severe_mucoid_plugging(self) -> None:
        res = evaluate_round160_herschel_bulkley_yield_stress("cf_severe_mucoid_plugging")
        self.assertFalse(res.crosses_kahn_kalai_threshold)
        self.assertTrue(res.ciliary_stall)
        self.assertEqual(res.mcc_velocity_mm_min, 0.0)
        self.assertEqual(res.plug_radius_fraction, 1.0)
        self.assertEqual(res.extrusion_flow_rate_mm3_s, 0.0)

    def test_kahn_kalai_curative_liquefaction_boundary(self) -> None:
        res = evaluate_round160_herschel_bulkley_yield_stress("kahn_kalai_curative_liquefaction_boundary")
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertFalse(res.ciliary_stall)
        self.assertGreaterEqual(res.mcc_velocity_mm_min, 0.65)
        self.assertLessEqual(res.yield_stress_pa, 0.85)

    def test_rhdnase_enzymatic_liquefaction_rescue(self) -> None:
        res = evaluate_round160_herschel_bulkley_yield_stress("rhdnase_enzymatic_liquefaction_rescue")
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertFalse(res.ciliary_stall)
        self.assertGreaterEqual(res.mcc_velocity_mm_min, 3.80)
        self.assertGreater(res.extrusion_flow_rate_mm3_s, 0.0)

    def test_shear_stress_monotonicity(self) -> None:
        tau_1 = compute_herschel_bulkley_shear_stress(1.0, 1.0, 2.0, 0.4)
        tau_5 = compute_herschel_bulkley_shear_stress(5.0, 1.0, 2.0, 0.4)
        tau_10 = compute_herschel_bulkley_shear_stress(10.0, 1.0, 2.0, 0.4)
        self.assertGreater(tau_5, tau_1)
        self.assertGreater(tau_10, tau_5)

    def test_poiseuille_plug_radius_monotonicity(self) -> None:
        tau_w1, r_plug1 = compute_poiseuille_plug_radius_mm(1.5, 20.0, 20.0, 1.0)
        tau_w2, r_plug2 = compute_poiseuille_plug_radius_mm(1.5, 60.0, 20.0, 1.0)
        self.assertGreater(tau_w2, tau_w1)
        self.assertLessEqual(r_plug2, r_plug1)

    def test_complex_step_differential_viscosity(self) -> None:
        deriv_cs, deriv_ana, rel_err = complex_step_differential_viscosity(
            gamma_dot=2.5,
            tau_y=1.5,
            k_index=3.0,
            n_index=0.4,
        )
        self.assertGreater(deriv_cs, 0.0)
        self.assertLessEqual(rel_err, 1e-10)

    def test_negative_control_falsification(self) -> None:
        clean = evaluate_falsification_negative_control()
        self.assertTrue(clean["passed"])
        self.assertEqual(clean["rule"], "AH-10")

        p_yield = evaluate_falsification_negative_control(poison_negative_yield_stress=True)
        self.assertFalse(p_yield["passed"])

        p_const = evaluate_falsification_negative_control(poison_negative_consistency=True)
        self.assertFalse(p_const["passed"])

        p_flow = evaluate_falsification_negative_control(poison_negative_flow_index=True)
        self.assertFalse(p_flow["passed"])

        p_shear = evaluate_falsification_negative_control(poison_negative_shear_rate=True)
        self.assertFalse(p_shear["passed"])

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
