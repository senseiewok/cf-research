#!/usr/bin/env python3
"""Unit tests for Round 134: Cahn-Hilliard Phase-Field Dynamic Concentration Evolution & Mucin Spinodal Decomposition Interface Kinetics Bridge."""
from __future__ import annotations

import unittest

from tools.discovery.round134_cahn_hilliard_mucin_spinodal_interface_bridge import (
    FIVE_LAYER_LADDER,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    evaluate_falsification_negative_control,
    evaluate_round134_cahn_hilliard_mucin_spinodal,
)


class Round134CahnHilliardMucinSpinodalTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round134_cahn_hilliard_mucin_spinodal(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Failed for {key}")
            self.assertLess(res.complex_step_sensitivity, 0.0)
            self.assertGreater(res.phi_mucin, 0.0)
            self.assertGreaterEqual(res.osmotic_pressure_pa, 0.0)
            self.assertAlmostEqual(res.lifshitz_slyozov_exponent, 1.0 / 3.0, places=4)

    def test_healthy_hydrated_asl_stability(self) -> None:
        healthy = evaluate_round134_cahn_hilliard_mucin_spinodal("healthy_hydrated_asl")
        self.assertTrue(healthy.crosses_kahn_kalai_threshold)
        self.assertFalse(healthy.is_spinodal_unstable)
        self.assertFalse(healthy.cilia_compressed)
        self.assertGreater(healthy.spinodal_f_second_deriv, 0.0)
        self.assertEqual(healthy.interfacial_tension_mn_m, 0.0)
        self.assertEqual(healthy.coarsening_radius_um, 0.0)
        self.assertLessEqual(healthy.osmotic_pressure_pa, 300.0)

    def test_cf_dehydrated_spinodal_instability_and_ciliary_compression(self) -> None:
        cf = evaluate_round134_cahn_hilliard_mucin_spinodal("cf_dehydrated_spinodal")
        self.assertFalse(cf.crosses_kahn_kalai_threshold)
        self.assertTrue(cf.is_spinodal_unstable)
        self.assertTrue(cf.cilia_compressed)
        self.assertLess(cf.spinodal_f_second_deriv, 0.0)
        self.assertGreater(cf.interfacial_tension_mn_m, 1.0)
        self.assertGreaterEqual(cf.coarsening_radius_um, 1.0)
        self.assertGreaterEqual(cf.osmotic_pressure_pa, 1800.0)

    def test_kahn_kalai_threshold_rescue(self) -> None:
        kk = evaluate_round134_cahn_hilliard_mucin_spinodal("kahn_kalai_threshold_11_8pct")
        self.assertTrue(kk.crosses_kahn_kalai_threshold)
        self.assertFalse(kk.cilia_compressed)
        self.assertLessEqual(kk.osmotic_pressure_pa, 600.0)
        self.assertLessEqual(kk.phi_mucin, 0.050)

    def test_rehydration_and_trikafta_rescue(self) -> None:
        hs = evaluate_round134_cahn_hilliard_mucin_spinodal("hypertonic_saline_rehydrated")
        tri = evaluate_round134_cahn_hilliard_mucin_spinodal("trikafta_full_clearance_rescue")
        for res in (hs, tri):
            self.assertTrue(res.crosses_kahn_kalai_threshold)
            self.assertFalse(res.is_spinodal_unstable)
            self.assertFalse(res.cilia_compressed)
            self.assertLessEqual(res.osmotic_pressure_pa, 450.0)
            self.assertLessEqual(res.phi_mucin, 0.035)

    def test_negative_control_falsification(self) -> None:
        normal = evaluate_falsification_negative_control(poison_negative_gradient_energy=False)
        self.assertTrue(normal["passed"])
        self.assertEqual(normal["rule"], "AH-10")
        poisoned = evaluate_falsification_negative_control(poison_negative_gradient_energy=True)
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
