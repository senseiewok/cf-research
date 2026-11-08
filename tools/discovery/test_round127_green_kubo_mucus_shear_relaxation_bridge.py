#!/usr/bin/env python3
"""Unit tests for Round 127: Green-Kubo Shear Relaxation & Crosslinked Mucin Network Viscoelastic Modulus Invariant Bridge."""
from __future__ import annotations

import unittest

from tools.discovery.round127_green_kubo_mucus_shear_relaxation_bridge import (
    FIVE_LAYER_LADDER,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    evaluate_falsification_negative_control,
    evaluate_round127_green_kubo_mucus_shear_relaxation,
)


class Round127GreenKuboMucusShearRelaxationTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round127_green_kubo_mucus_shear_relaxation(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Failed for {key}")
            self.assertLess(res.complex_step_sensitivity, 0.0)
            self.assertGreater(res.storage_modulus_g_prime_pa, 0.0)
            self.assertGreater(res.loss_modulus_g_double_prime_pa, 0.0)
            self.assertGreater(res.deborah_number, 0.0)
            self.assertGreaterEqual(res.zero_shear_viscosity_pa_s, 0.001)

    def test_wild_type_healthy_rheology(self) -> None:
        wt = evaluate_round127_green_kubo_mucus_shear_relaxation("wild_type_healthy")
        self.assertTrue(wt.crosses_kahn_kalai_threshold)
        self.assertLessEqual(wt.deborah_number, 0.25)
        self.assertLessEqual(wt.plateau_modulus_gp_pa, 0.30)
        self.assertLessEqual(wt.relaxation_time_ms, 25.0)
        self.assertLessEqual(wt.zero_shear_viscosity_pa_s, 0.015)

    def test_cf_homozygous_dehydrated_mucus_jamming(self) -> None:
        cf = evaluate_round127_green_kubo_mucus_shear_relaxation("cf_homozygous_dehydrated")
        self.assertFalse(cf.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(cf.deborah_number, 5.0)
        self.assertGreaterEqual(cf.plateau_modulus_gp_pa, 7.0)
        self.assertGreaterEqual(cf.relaxation_time_ms, 500.0)
        self.assertGreaterEqual(cf.zero_shear_viscosity_pa_s, 3.0)

    def test_kahn_kalai_threshold_rescue(self) -> None:
        kk = evaluate_round127_green_kubo_mucus_shear_relaxation("kahn_kalai_threshold_11_8pct")
        self.assertTrue(kk.crosses_kahn_kalai_threshold)
        self.assertLessEqual(kk.deborah_number, 1.0)
        self.assertLessEqual(kk.plateau_modulus_gp_pa, 1.50)

    def test_rhdnase_and_hypertonic_saline_fluidization(self) -> None:
        dnase = evaluate_round127_green_kubo_mucus_shear_relaxation("rhdnase_pulmozyme_rescued")
        saline = evaluate_round127_green_kubo_mucus_shear_relaxation("hypertonic_saline_7pct_rehydrated")
        for res in (dnase, saline):
            self.assertTrue(res.crosses_kahn_kalai_threshold)
            self.assertLessEqual(res.deborah_number, 0.50)
            self.assertLessEqual(res.plateau_modulus_gp_pa, 1.0)
            self.assertLessEqual(res.relaxation_time_ms, 50.0)

    def test_negative_control_falsification(self) -> None:
        normal = evaluate_falsification_negative_control(poison_negative_plateau_modulus=False)
        self.assertTrue(normal["passed"])
        self.assertEqual(normal["rule"], "AH-10")
        poisoned = evaluate_falsification_negative_control(poison_negative_plateau_modulus=True)
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
