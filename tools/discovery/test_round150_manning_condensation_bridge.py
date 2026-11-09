#!/usr/bin/env python3
"""Unit tests for Round 150: Manning Polyelectrolyte Counterion Condensation Invariant Bridge (T-0803)."""
from __future__ import annotations

import unittest

from tools.discovery.round150_manning_condensation_bridge import (
    FIVE_LAYER_LADDER,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    compute_manning_condensation,
    evaluate_falsification_negative_control,
    evaluate_round150_manning_condensation,
)


class Round150ManningCondensationTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round150_manning_condensation(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Failed for {key}")
            self.assertGreater(res.mucin_charge_spacing_nm, 0.0)
            self.assertGreater(res.bjerrum_length_nm, 0.0)
            self.assertGreater(res.effective_charge_fraction, 0.0)
            self.assertLessEqual(res.effective_charge_fraction, 1.0)
            self.assertGreater(res.periciliary_layer_height_um, 0.0)
            self.assertGreaterEqual(res.epithelial_short_circuit_current_ua_cm2, 0.0)

    def test_wild_type_periciliary_brush_equilibrium(self) -> None:
        wt = evaluate_round150_manning_condensation("wild_type_periciliary_brush_equilibrium")
        self.assertTrue(wt.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(wt.periciliary_layer_height_um, 6.5)
        self.assertGreaterEqual(wt.epithelial_short_circuit_current_ua_cm2, 30.0)
        self.assertAlmostEqual(wt.fractional_wild_type_current, 1.0, delta=0.05)
        self.assertEqual(wt.condensation_fraction, 0.0)
        self.assertEqual(wt.effective_charge_fraction, 1.0)
        self.assertLess(wt.manning_parameter, 1.0)

    def test_cf_mucin_hyperconcentration_and_compression(self) -> None:
        cf = evaluate_round150_manning_condensation("cf_mucin_hyperconcentration_condensation")
        self.assertFalse(cf.crosses_kahn_kalai_threshold)
        self.assertLessEqual(cf.fractional_wild_type_current, 0.05)
        self.assertLessEqual(cf.epithelial_short_circuit_current_ua_cm2, 0.5)
        self.assertLess(cf.periciliary_layer_height_um, 3.0)
        self.assertGreater(cf.manning_parameter, 1.5)
        self.assertGreater(cf.condensation_fraction, 0.40)
        self.assertLess(cf.effective_charge_fraction, 0.60)

    def test_divalent_calcium_hypercondensation(self) -> None:
        ca = evaluate_round150_manning_condensation("divalent_calcium_hypercondensation")
        self.assertFalse(ca.crosses_kahn_kalai_threshold)
        self.assertEqual(ca.counterion_valence, 2)
        self.assertGreater(ca.condensation_fraction, 0.50)
        self.assertLess(ca.effective_charge_fraction, 0.50)
        self.assertLess(ca.periciliary_layer_height_um, 3.0)

    def test_kahn_kalai_osmotic_equilibrium_boundary(self) -> None:
        kk = evaluate_round150_manning_condensation("kahn_kalai_osmotic_equilibrium_boundary")
        self.assertTrue(kk.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(kk.periciliary_layer_height_um, 5.0)
        self.assertGreaterEqual(kk.epithelial_short_circuit_current_ua_cm2, 3.75)
        self.assertGreaterEqual(kk.fractional_wild_type_current, 0.118)

    def test_trikafta_rehydration_restoration(self) -> None:
        rescued = evaluate_round150_manning_condensation("trikafta_rehydration_restoration")
        self.assertTrue(rescued.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(rescued.periciliary_layer_height_um, 6.5)
        self.assertGreaterEqual(rescued.epithelial_short_circuit_current_ua_cm2, 25.0)
        self.assertGreaterEqual(rescued.fractional_wild_type_current, 0.85)

    def test_negative_control_falsification(self) -> None:
        normal = evaluate_falsification_negative_control(poison_unphysical_charge_spacing=False)
        self.assertTrue(normal["passed"])
        self.assertEqual(normal["rule"], "AH-10")
        poisoned = evaluate_falsification_negative_control(poison_unphysical_charge_spacing=True)
        self.assertFalse(poisoned["passed"])

    def test_manning_parameter_threshold_transition(self) -> None:
        # Below threshold: spacing = 1.2 nm > l_B (~0.714 nm)
        sub_crit = compute_manning_condensation(spacing_b_nm=1.2, valence_zc=1)
        self.assertLess(sub_crit["manning_parameter"], 1.0)
        self.assertEqual(sub_crit["condensation_fraction"], 0.0)
        self.assertEqual(sub_crit["effective_charge_fraction"], 1.0)

        # Above threshold: spacing = 0.4 nm < l_B (~0.714 nm)
        super_crit = compute_manning_condensation(spacing_b_nm=0.4, valence_zc=1)
        self.assertGreater(super_crit["manning_parameter"], 1.0)
        self.assertGreater(super_crit["condensation_fraction"], 0.0)
        self.assertLess(super_crit["effective_charge_fraction"], 1.0)
        self.assertAlmostEqual(
            super_crit["condensation_fraction"] + super_crit["effective_charge_fraction"],
            1.0,
            places=6,
        )

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
