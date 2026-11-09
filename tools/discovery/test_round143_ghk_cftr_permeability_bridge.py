#!/usr/bin/env python3
"""Unit tests for Round 143: Goldman-Hodgkin-Katz Voltage-Dependent Current & Permeability Ratio Invariant Bridge (T-0792)."""
from __future__ import annotations

import unittest

from tools.discovery.round143_ghk_cftr_permeability_bridge import (
    FIVE_LAYER_LADDER,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    evaluate_falsification_negative_control,
    evaluate_round143_ghk_cftr_permeability,
)


class Round143GhkCftrPermeabilityTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round143_ghk_cftr_permeability(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Failed for {key}")
            self.assertGreater(res.single_channel_conductance_ps, 0.0)
            self.assertGreater(res.permeability_ratio_hco3_cl, 0.0)
            self.assertLess(res.permeability_ratio_hco3_cl, 1.0)
            self.assertLess(res.reversal_potential_mv, 0.0)
            self.assertGreater(res.epithelial_short_circuit_current_ua_cm2, 0.0)
            self.assertGreater(res.complex_step_voltage_sensitivity, 0.0)

    def test_wild_type_cftr_open_pore(self) -> None:
        wt = evaluate_round143_ghk_cftr_permeability("wild_type_cftr_open_pore")
        self.assertTrue(wt.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(wt.single_channel_conductance_ps, 8.0)
        self.assertGreaterEqual(wt.permeability_ratio_hco3_cl, 0.20)
        self.assertLessEqual(wt.reversal_potential_mv, -25.0)
        self.assertGreaterEqual(wt.epithelial_short_circuit_current_ua_cm2, 30.0)
        self.assertAlmostEqual(wt.fractional_wild_type_current, 1.0, delta=0.05)

    def test_cf_pore_mutant_conductance_collapse(self) -> None:
        r334w = evaluate_round143_ghk_cftr_permeability("cf_r334w_conductance_defect")
        r347p = evaluate_round143_ghk_cftr_permeability("cf_r347p_selectivity_collapse")
        for mutant in (r334w, r347p):
            self.assertFalse(mutant.crosses_kahn_kalai_threshold)
            self.assertLessEqual(mutant.single_channel_conductance_ps, 1.50)
            self.assertLessEqual(mutant.permeability_ratio_hco3_cl, 0.08)
            self.assertLessEqual(mutant.epithelial_short_circuit_current_ua_cm2, 2.50)
            self.assertLessEqual(mutant.fractional_wild_type_current, 0.10)

    def test_kahn_kalai_conductance_boundary(self) -> None:
        kk = evaluate_round143_ghk_cftr_permeability("kahn_kalai_conductance_boundary")
        self.assertTrue(kk.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(kk.single_channel_conductance_ps, 2.00)
        self.assertGreaterEqual(kk.epithelial_short_circuit_current_ua_cm2, 4.13)
        self.assertGreaterEqual(kk.fractional_wild_type_current, 0.118)

    def test_trikafta_ivacaftor_pore_rescue(self) -> None:
        rescued = evaluate_round143_ghk_cftr_permeability("trikafta_ivacaftor_pore_rescue")
        self.assertTrue(rescued.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(rescued.single_channel_conductance_ps, 7.8)
        self.assertGreaterEqual(rescued.permeability_ratio_hco3_cl, 0.20)
        self.assertLessEqual(rescued.reversal_potential_mv, -25.0)
        self.assertGreaterEqual(rescued.epithelial_short_circuit_current_ua_cm2, 25.0)
        self.assertGreaterEqual(rescued.fractional_wild_type_current, 0.80)

    def test_negative_control_falsification(self) -> None:
        normal = evaluate_falsification_negative_control(poison_negative_permeability=False)
        self.assertTrue(normal["passed"])
        self.assertEqual(normal["rule"], "AH-10")
        poisoned = evaluate_falsification_negative_control(poison_negative_permeability=True)
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
