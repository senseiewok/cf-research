#!/usr/bin/env python3
"""Unit tests for Round 146: Gouy-Chapman Double-Layer Electrostatic Surface Potential Invariant Bridge (T-0797)."""
from __future__ import annotations

import unittest

from tools.discovery.round146_gouy_chapman_membrane_potential_bridge import (
    FIVE_LAYER_LADDER,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    evaluate_falsification_negative_control,
    evaluate_round146_gouy_chapman_membrane_potential,
)


class Round146GouyChapmanMembranePotentialTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round146_gouy_chapman_membrane_potential(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Failed for {key}")
            self.assertGreater(res.debye_length_nm, 0.0)
            self.assertLess(res.surface_potential_mv, 0.0)
            self.assertGreater(res.surface_chloride_concentration_mm, 0.0)
            self.assertGreater(res.effective_conductance_ps, 0.0)
            self.assertGreaterEqual(res.epithelial_short_circuit_current_ua_cm2, 0.0)
            self.assertGreater(res.complex_step_charge_sensitivity, 0.0)

    def test_wild_type_apical_membrane_bilayer(self) -> None:
        wt = evaluate_round146_gouy_chapman_membrane_potential("wild_type_apical_membrane_bilayer")
        self.assertTrue(wt.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(wt.effective_conductance_ps, 8.0)
        self.assertLessEqual(wt.surface_potential_mv, -25.0)
        self.assertGreaterEqual(wt.surface_potential_mv, -38.0)
        self.assertGreaterEqual(wt.epithelial_short_circuit_current_ua_cm2, 30.0)
        self.assertAlmostEqual(wt.fractional_wild_type_current, 1.0, delta=0.05)

    def test_cf_apical_defects_and_mucin_hyperpolarization(self) -> None:
        pip2 = evaluate_round146_gouy_chapman_membrane_potential("cf_pip2_depleted_apical_membrane")
        mucin = evaluate_round146_gouy_chapman_membrane_potential("cf_polyanionic_mucin_hyperpolarization")
        for mutant in (pip2, mucin):
            self.assertFalse(mutant.crosses_kahn_kalai_threshold)
            self.assertLessEqual(mutant.fractional_wild_type_current, 0.10)
        self.assertGreaterEqual(pip2.surface_potential_mv, -10.0)  # less negative due to depleted negative lipids
        self.assertLessEqual(mucin.surface_potential_mv, -60.0)    # hyperpolarized by anionic mucins

    def test_kahn_kalai_electrostatic_boundary(self) -> None:
        kk = evaluate_round146_gouy_chapman_membrane_potential("kahn_kalai_electrostatic_boundary")
        self.assertTrue(kk.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(kk.effective_conductance_ps, 2.0)
        self.assertGreaterEqual(kk.epithelial_short_circuit_current_ua_cm2, 3.79)
        self.assertGreaterEqual(kk.fractional_wild_type_current, 0.118)

    def test_trikafta_lipid_homeostasis_rescue(self) -> None:
        rescued = evaluate_round146_gouy_chapman_membrane_potential("trikafta_lipid_homeostasis_rescue")
        self.assertTrue(rescued.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(rescued.effective_conductance_ps, 8.0)
        self.assertLessEqual(rescued.surface_potential_mv, -25.0)
        self.assertGreaterEqual(rescued.epithelial_short_circuit_current_ua_cm2, 25.0)
        self.assertGreaterEqual(rescued.fractional_wild_type_current, 0.85)

    def test_negative_control_falsification(self) -> None:
        normal = evaluate_falsification_negative_control(poison_negative_bulk_salt=False)
        self.assertTrue(normal["passed"])
        self.assertEqual(normal["rule"], "AH-10")
        poisoned = evaluate_falsification_negative_control(poison_negative_bulk_salt=True)
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
