#!/usr/bin/env python3
"""Unit tests for Round 125: Kramers-Moyal Finite-Sample Expansion & CFTR Patch-Clamp Stochastic Transition Bridge."""
from __future__ import annotations

import unittest

from tools.discovery.round125_kramers_moyal_cftr_patch_clamp_transition_bridge import (
    FIVE_LAYER_LADDER,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    evaluate_falsification_negative_control,
    evaluate_round125_kramers_moyal_cftr_patch_clamp,
)


class Round125KramersMoyalCFTRPatchClampTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round125_kramers_moyal_cftr_patch_clamp(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Failed for {key}")
            self.assertGreater(res.complex_step_sensitivity, 0.0)
            self.assertGreater(res.kramers_moyal_diffusion_d2, 0.0)

    def test_wild_type_healthy_open_probability_and_dwell_times(self) -> None:
        wt = evaluate_round125_kramers_moyal_cftr_patch_clamp("wild_type_healthy")
        self.assertTrue(wt.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(wt.open_probability, 0.35)
        self.assertGreaterEqual(wt.mean_open_burst_ms, 300.0)
        self.assertLessEqual(wt.mean_closed_dwell_ms, 600.0)

    def test_f508del_untreated_severe_gating_defect(self) -> None:
        cf = evaluate_round125_kramers_moyal_cftr_patch_clamp("f508del_homozygous_untreated")
        self.assertFalse(cf.crosses_kahn_kalai_threshold)
        self.assertLess(cf.open_probability, 0.05)
        self.assertLess(cf.mean_open_burst_ms, 80.0)
        self.assertGreater(cf.mean_closed_dwell_ms, 1500.0)
        self.assertGreater(cf.closing_rate_k_oc_s, 15.0)

    def test_kahn_kalai_threshold_rescues_open_probability(self) -> None:
        kk = evaluate_round125_kramers_moyal_cftr_patch_clamp("kahn_kalai_threshold_11_8pct")
        self.assertAlmostEqual(kk.cftr_function_pct, 11.8, delta=0.5)
        self.assertTrue(kk.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(kk.open_probability, 0.15)
        self.assertGreaterEqual(kk.mean_open_burst_ms, 180.0)

    def test_trikafta_rescued_burst_stabilization(self) -> None:
        tri = evaluate_round125_kramers_moyal_cftr_patch_clamp("trikafta_rescued")
        self.assertTrue(tri.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(tri.open_probability, 0.45)
        self.assertGreaterEqual(tri.mean_open_burst_ms, 700.0)
        self.assertLessEqual(tri.closing_rate_k_oc_s, 2.0)

    def test_ivacaftor_g551d_gating_rescue(self) -> None:
        iva = evaluate_round125_kramers_moyal_cftr_patch_clamp("ivacaftor_monotherapy_g551d")
        self.assertTrue(iva.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(iva.open_probability, 0.40)
        self.assertGreaterEqual(iva.mean_open_burst_ms, 600.0)

    def test_kramers_moyal_finite_sample_bias_and_positive_diffusion(self) -> None:
        res = evaluate_round125_kramers_moyal_cftr_patch_clamp("wild_type_healthy")
        self.assertGreater(res.kramers_moyal_diffusion_d2, 0.0)
        self.assertGreater(res.kramers_moyal_jump_d4, 0.0)
        self.assertGreaterEqual(res.finite_sample_bias_pct, 0.0)
        self.assertLess(res.finite_sample_bias_pct, 10.0)

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
