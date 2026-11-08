#!/usr/bin/env python3
"""Unit tests for Round 126: Zwanzig-Mori Generalized Langevin Equation & Airway Mucus Microrheology Memory Kernel Invariant Bridge."""
from __future__ import annotations

import unittest

from tools.discovery.round126_zwanzig_mori_mucus_memory_kernel_bridge import (
    FIVE_LAYER_LADDER,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    evaluate_falsification_negative_control,
    evaluate_round126_zwanzig_mori_mucus_memory_kernel,
)


class Round126ZwanzigMoriMucusMemoryKernelTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round126_zwanzig_mori_mucus_memory_kernel(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Failed for {key}")
            self.assertGreater(res.complex_step_sensitivity, 0.0)
            self.assertGreater(res.storage_modulus_g_prime_pa, 0.0)
            self.assertGreater(res.loss_modulus_g_double_prime_pa, 0.0)
            self.assertGreater(res.anomalous_diffusion_exponent_alpha, 0.0)
            self.assertLessEqual(res.anomalous_diffusion_exponent_alpha, 1.0)

    def test_wild_type_healthy_rheology(self) -> None:
        wt = evaluate_round126_zwanzig_mori_mucus_memory_kernel("wild_type_healthy")
        self.assertTrue(wt.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(wt.anomalous_diffusion_exponent_alpha, 0.85)
        self.assertLessEqual(wt.plateau_modulus_g0_pa, 0.30)
        self.assertLessEqual(wt.relaxation_time_ms, 20.0)
        self.assertLessEqual(wt.zero_shear_viscosity_pa_s, 0.01)

    def test_cf_homozygous_dehydrated_mucus_jamming(self) -> None:
        cf = evaluate_round126_zwanzig_mori_mucus_memory_kernel("cf_homozygous_dehydrated")
        self.assertFalse(cf.crosses_kahn_kalai_threshold)
        self.assertLess(cf.anomalous_diffusion_exponent_alpha, 0.25)
        self.assertGreaterEqual(cf.plateau_modulus_g0_pa, 5.0)
        self.assertGreaterEqual(cf.relaxation_time_ms, 500.0)
        self.assertGreaterEqual(cf.zero_shear_viscosity_pa_s, 2.0)

    def test_kahn_kalai_threshold_rescue(self) -> None:
        kk = evaluate_round126_zwanzig_mori_mucus_memory_kernel("kahn_kalai_threshold_11_8pct")
        self.assertTrue(kk.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(kk.anomalous_diffusion_exponent_alpha, 0.50)
        self.assertLessEqual(kk.plateau_modulus_g0_pa, 1.50)

    def test_rhdnase_and_hypertonic_saline_fluidization(self) -> None:
        dnase = evaluate_round126_zwanzig_mori_mucus_memory_kernel("rhdnase_pulmozyme_rescued")
        saline = evaluate_round126_zwanzig_mori_mucus_memory_kernel("hypertonic_saline_7pct_rehydrated")
        for res in (dnase, saline):
            self.assertTrue(res.crosses_kahn_kalai_threshold)
            self.assertGreaterEqual(res.anomalous_diffusion_exponent_alpha, 0.70)
            self.assertLessEqual(res.plateau_modulus_g0_pa, 1.0)
            self.assertLessEqual(res.relaxation_time_ms, 50.0)

    def test_negative_control_falsification(self) -> None:
        normal = evaluate_falsification_negative_control(poison_negative_storage_modulus=False)
        self.assertTrue(normal["passed"])
        self.assertEqual(normal["rule"], "AH-10")
        poisoned = evaluate_falsification_negative_control(poison_negative_storage_modulus=True)
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
