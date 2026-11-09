#!/usr/bin/env python3
"""Unit tests for Round 155: Acoustic Microbubble Cavitation & Rayleigh-Plesset Inertial Invariant Bridge (T-0814)."""
from __future__ import annotations

import math
import unittest

from tools.discovery.round155_rayleigh_plesset_cavitation_bridge import (
    FIVE_LAYER_LADDER,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    compute_blake_threshold_pa,
    compute_resonance_frequency_hz,
    evaluate_falsification_negative_control,
    evaluate_round155_rayleigh_plesset_cavitation,
)


class Round155CavitationTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round155_rayleigh_plesset_cavitation(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Failed for {key}")
            self.assertGreater(res.resonance_frequency_mhz, 0.0)
            self.assertGreater(res.blake_threshold_kpa, 100.0)
            self.assertGreaterEqual(res.microjet_velocity_m_s, 0.0)
            self.assertGreaterEqual(res.shear_stress_pa, 0.0)
            self.assertGreater(res.effective_mucus_diffusivity_um2_s, 0.0)
            self.assertGreaterEqual(res.basal_cell_transfection_fraction, 0.0)

    def test_unsonicated_hyperconcentrated_cf_mucus(self) -> None:
        unsonicated = evaluate_round155_rayleigh_plesset_cavitation("unsonicated_hyperconcentrated_cf_mucus")
        self.assertFalse(unsonicated.crosses_kahn_kalai_threshold)
        self.assertEqual(unsonicated.microjet_velocity_m_s, 0.0)
        self.assertEqual(unsonicated.shear_stress_pa, 0.0)
        self.assertLess(unsonicated.basal_cell_transfection_fraction, 0.01)
        self.assertLess(unsonicated.fractional_wild_type_transfection, 0.10)

    def test_diagnostic_low_intensity_ultrasound(self) -> None:
        diag = evaluate_round155_rayleigh_plesset_cavitation("diagnostic_low_intensity_ultrasound")
        self.assertFalse(diag.crosses_kahn_kalai_threshold)
        self.assertFalse(diag.is_inertial_cavitation)
        self.assertEqual(diag.microjet_velocity_m_s, 0.0)
        self.assertLess(diag.shear_stress_pa, 12.0)  # Below critical sonoporation shear
        self.assertLess(diag.basal_cell_transfection_fraction, 0.05)

    def test_kahn_kalai_inertial_cavitation_threshold(self) -> None:
        kk = evaluate_round155_rayleigh_plesset_cavitation("kahn_kalai_inertial_cavitation_threshold")
        self.assertTrue(kk.crosses_kahn_kalai_threshold)
        self.assertTrue(kk.is_inertial_cavitation)
        self.assertGreater(kk.microjet_velocity_m_s, 5.0)
        self.assertGreaterEqual(kk.basal_cell_transfection_fraction, 0.118)
        self.assertFalse(kk.exceeds_cytotoxicity_limit)

    def test_focused_ultrasound_microbubble_rescue(self) -> None:
        rescued = evaluate_round155_rayleigh_plesset_cavitation("focused_ultrasound_targeted_microbubble_rescue")
        self.assertTrue(rescued.crosses_kahn_kalai_threshold)
        self.assertTrue(rescued.is_inertial_cavitation)
        self.assertGreaterEqual(rescued.microjet_velocity_m_s, 15.0)
        self.assertGreaterEqual(rescued.shear_stress_pa, 20.0)
        self.assertGreaterEqual(rescued.basal_cell_transfection_fraction, 0.30)
        self.assertGreaterEqual(rescued.fractional_wild_type_transfection, 2.5)
        self.assertFalse(rescued.exceeds_cytotoxicity_limit)

    def test_excessive_cavitation_cytotoxicity(self) -> None:
        excess = evaluate_round155_rayleigh_plesset_cavitation("excessive_high_intensity_cavitation_cytotoxicity")
        self.assertFalse(excess.crosses_kahn_kalai_threshold)
        self.assertTrue(excess.exceeds_cytotoxicity_limit)
        self.assertGreater(excess.shear_stress_pa, 45.0)
        self.assertLess(excess.basal_cell_transfection_fraction, 0.01)

    def test_negative_control_falsification(self) -> None:
        normal = evaluate_falsification_negative_control(poison_negative_radius=False)
        self.assertTrue(normal["passed"])
        self.assertEqual(normal["rule"], "AH-10")

        poisoned_r0 = evaluate_falsification_negative_control(poison_negative_radius=True)
        self.assertFalse(poisoned_r0["passed"])

        poisoned_p0 = evaluate_falsification_negative_control(poison_negative_pressure=True)
        self.assertFalse(poisoned_p0["passed"])

    def test_blake_threshold_and_minnaert_resonance(self) -> None:
        r0 = 1.5e-6
        p0 = 101325.0
        sigma = 0.050
        p_blake = compute_blake_threshold_pa(r0, p0, sigma)
        self.assertGreater(p_blake, p0)
        self.assertAlmostEqual(p_blake, p0 + 0.77 * (sigma / r0), places=1)

        f0 = compute_resonance_frequency_hz(r0, p0, 1000.0, sigma, 1.40)
        self.assertGreater(f0, 1.5e6)  # ~ 2.25 MHz
        self.assertLess(f0, 3.5e6)

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
