#!/usr/bin/env python3
"""Unit tests for Round 132: Flory-Huggins Polymer Solution Thermodynamics & Mucin Phase Separation/Coacervation Invariant Bridge."""
from __future__ import annotations

import unittest

from tools.discovery.round132_flory_huggins_mucin_phase_separation_bridge import (
    FIVE_LAYER_LADDER,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    evaluate_falsification_negative_control,
    evaluate_round132_flory_huggins_mucin_phase_separation,
)


class Round132FloryHugginsMucinPhaseSeparationTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round132_flory_huggins_mucin_phase_separation(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Failed for {key}")
            self.assertLess(res.complex_step_sensitivity, 0.0)
            self.assertGreater(res.critical_phi, 0.0)
            self.assertLess(res.critical_phi, 0.1)
            self.assertGreater(res.critical_chi, 0.5)
            self.assertGreater(res.spinodal_chi, 0.0)
            self.assertGreaterEqual(res.osmotic_pressure_pa, 0.0)

    def test_healthy_hydrated_asl(self) -> None:
        healthy = evaluate_round132_flory_huggins_mucin_phase_separation("healthy_hydrated_asl")
        self.assertTrue(healthy.crosses_kahn_kalai_threshold)
        self.assertFalse(healthy.is_phase_separated)
        self.assertFalse(healthy.cilia_compressed)
        self.assertLessEqual(healthy.osmotic_pressure_pa, 400.0)
        self.assertLessEqual(healthy.phi_mucin, 0.025)

    def test_cf_dehydrated_unstable_mucin(self) -> None:
        cf = evaluate_round132_flory_huggins_mucin_phase_separation("cf_dehydrated_unstable")
        self.assertFalse(cf.crosses_kahn_kalai_threshold)
        self.assertTrue(cf.is_phase_separated)
        self.assertTrue(cf.cilia_compressed)
        self.assertGreaterEqual(cf.osmotic_pressure_pa, 1500.0)
        self.assertGreaterEqual(cf.phi_mucin, 0.08)

    def test_kahn_kalai_threshold_rescue(self) -> None:
        kk = evaluate_round132_flory_huggins_mucin_phase_separation("kahn_kalai_threshold_11_8pct")
        self.assertTrue(kk.crosses_kahn_kalai_threshold)
        self.assertFalse(kk.is_phase_separated)
        self.assertFalse(kk.cilia_compressed)
        self.assertLessEqual(kk.osmotic_pressure_pa, 600.0)

    def test_hydration_and_modulator_rescue(self) -> None:
        hs = evaluate_round132_flory_huggins_mucin_phase_separation("hypertonic_saline_rescued")
        tri = evaluate_round132_flory_huggins_mucin_phase_separation("trikafta_full_hydration_rescue")
        for res in (hs, tri):
            self.assertTrue(res.crosses_kahn_kalai_threshold)
            self.assertFalse(res.is_phase_separated)
            self.assertFalse(res.cilia_compressed)
            self.assertLessEqual(res.osmotic_pressure_pa, 500.0)
            self.assertLessEqual(res.phi_mucin, 0.035)

    def test_negative_control_falsification(self) -> None:
        normal = evaluate_falsification_negative_control(poison_negative_chain_length=False)
        self.assertTrue(normal["passed"])
        self.assertEqual(normal["rule"], "AH-10")
        poisoned = evaluate_falsification_negative_control(poison_negative_chain_length=True)
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
