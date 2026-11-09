#!/usr/bin/env python3
"""Unit tests for Round 144: Smoluchowski-Debye Diffusion-Limited Reaction Kinetics & Electrostatic Steering Invariant Bridge (T-0793)."""
from __future__ import annotations

import unittest

from tools.discovery.round144_smoluchowski_debye_cftr_capture_bridge import (
    FIVE_LAYER_LADDER,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    evaluate_falsification_negative_control,
    evaluate_round144_smoluchowski_debye_cftr_capture,
)


class Round144SmoluchowskiDebyeCftrCaptureTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round144_smoluchowski_debye_cftr_capture(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Failed for {key}")
            self.assertGreaterEqual(res.target_valence, 0.0)
            self.assertGreaterEqual(res.electrostatic_steering_factor, 1.0)
            self.assertGreater(res.neutral_diffusion_rate_m_inv_s_inv, 0.0)
            self.assertGreaterEqual(res.effective_capture_rate_m_inv_s_inv, res.neutral_diffusion_rate_m_inv_s_inv * 0.99)
            self.assertGreater(res.debye_screening_length_nm, 0.0)
            self.assertGreater(res.michaelis_constant_km_mm, 0.0)
            self.assertGreater(res.chloride_capture_turnover_rate_ions_s, 0.0)
            self.assertGreater(res.complex_step_valence_sensitivity, 0.0)

    def test_wild_type_electrostatic_steering(self) -> None:
        wt = evaluate_round144_smoluchowski_debye_cftr_capture("wild_type_electrostatic_steering")
        self.assertTrue(wt.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(wt.target_valence, 1.70)
        self.assertGreaterEqual(wt.electrostatic_steering_factor, 1.70)
        self.assertLessEqual(wt.michaelis_constant_km_mm, 30.0)
        self.assertGreaterEqual(wt.chloride_capture_turnover_rate_ions_s, 1.8e7)
        self.assertAlmostEqual(wt.fractional_wild_type_turnover, 1.0, delta=0.05)

    def test_cf_pore_mutant_steering_collapse(self) -> None:
        r334w = evaluate_round144_smoluchowski_debye_cftr_capture("cf_r334w_charge_depletion")
        r347p = evaluate_round144_smoluchowski_debye_cftr_capture("cf_r347p_vestibule_collapse")
        for mutant in (r334w, r347p):
            self.assertFalse(mutant.crosses_kahn_kalai_threshold)
            self.assertLessEqual(mutant.target_valence, 0.40)
            self.assertLessEqual(mutant.electrostatic_steering_factor, 1.20)
            self.assertGreaterEqual(mutant.michaelis_constant_km_mm, 34.0)
            self.assertLessEqual(mutant.chloride_capture_turnover_rate_ions_s, 2.5e6)
            self.assertLessEqual(mutant.fractional_wild_type_turnover, 0.13)

    def test_kahn_kalai_capture_boundary(self) -> None:
        kk = evaluate_round144_smoluchowski_debye_cftr_capture("kahn_kalai_capture_boundary")
        self.assertTrue(kk.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(kk.target_valence, 0.50)
        self.assertGreaterEqual(kk.electrostatic_steering_factor, 1.20)
        self.assertGreaterEqual(kk.chloride_capture_turnover_rate_ions_s, 2.35e6)
        self.assertGreaterEqual(kk.fractional_wild_type_turnover, 0.118)

    def test_trikafta_allosteric_rescue(self) -> None:
        rescued = evaluate_round144_smoluchowski_debye_cftr_capture("trikafta_allosteric_rescue")
        self.assertTrue(rescued.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(rescued.target_valence, 1.60)
        self.assertGreaterEqual(rescued.electrostatic_steering_factor, 1.70)
        self.assertLessEqual(rescued.michaelis_constant_km_mm, 28.0)
        self.assertGreaterEqual(rescued.chloride_capture_turnover_rate_ions_s, 1.7e7)
        self.assertGreaterEqual(rescued.fractional_wild_type_turnover, 0.85)

    def test_negative_control_falsification(self) -> None:
        normal = evaluate_falsification_negative_control(poison_negative_diffusion=False)
        self.assertTrue(normal["passed"])
        self.assertEqual(normal["rule"], "AH-10")
        poisoned = evaluate_falsification_negative_control(poison_negative_diffusion=True)
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
