#!/usr/bin/env python3
"""Unit tests for Round 137: Rouse-Zimm Viscoelastic Hydrodynamic Relaxation & Entanglement Mesh Invariant Bridge."""
from __future__ import annotations

import unittest

from tools.discovery.round137_rouse_zimm_mucin_relaxation_bridge import (
    FIVE_LAYER_LADDER,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    evaluate_falsification_negative_control,
    evaluate_round137_rouse_zimm_mucin_relaxation,
)


class Round137RouseZimmMucinRelaxationTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round137_rouse_zimm_mucin_relaxation(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Failed for {key}")
            self.assertLess(res.complex_step_mesh_sensitivity, 0.0)
            self.assertGreater(res.phi_mucin, 0.0)
            self.assertGreater(res.mesh_correlation_length_nm, 0.0)
            self.assertGreaterEqual(res.storage_modulus_g_prime_pa, 0.0)
            self.assertGreater(res.loss_modulus_g_double_prime_pa, 0.0)
            self.assertGreater(res.nanoparticle_diffusivity_ratio, 0.0)
            self.assertLessEqual(res.nanoparticle_diffusivity_ratio, 1.0)

    def test_healthy_hydrated_asl_stability(self) -> None:
        healthy = evaluate_round137_rouse_zimm_mucin_relaxation("healthy_hydrated_asl")
        self.assertTrue(healthy.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(healthy.mesh_correlation_length_nm, 160.0)
        self.assertGreaterEqual(healthy.nanoparticle_diffusivity_ratio, 0.70)
        self.assertLessEqual(healthy.modulus_ratio_g_prime_to_g_double_prime, 0.60)
        self.assertLessEqual(healthy.phi_mucin, 0.020)

    def test_cf_dehydrated_hyperconcentration_and_mesh_entrapment(self) -> None:
        cf = evaluate_round137_rouse_zimm_mucin_relaxation("cf_dehydrated_hyperconcentrated")
        self.assertFalse(cf.crosses_kahn_kalai_threshold)
        self.assertLessEqual(cf.mesh_correlation_length_nm, 55.0)
        self.assertLessEqual(cf.nanoparticle_diffusivity_ratio, 0.050)
        self.assertGreaterEqual(cf.modulus_ratio_g_prime_to_g_double_prime, 2.0)
        self.assertGreaterEqual(cf.phi_mucin, 0.080)

    def test_kahn_kalai_threshold_rescue(self) -> None:
        kk = evaluate_round137_rouse_zimm_mucin_relaxation("kahn_kalai_delivery_threshold")
        self.assertTrue(kk.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(kk.mesh_correlation_length_nm, 85.0)
        self.assertGreaterEqual(kk.nanoparticle_diffusivity_ratio, 0.080)
        self.assertLessEqual(kk.phi_mucin, 0.050)

    def test_rehydration_and_trikafta_rescue(self) -> None:
        hs = evaluate_round137_rouse_zimm_mucin_relaxation("hypertonic_saline_rehydrated")
        tri = evaluate_round137_rouse_zimm_mucin_relaxation("trikafta_full_mucociliary_rescue")
        for res in (hs, tri):
            self.assertTrue(res.crosses_kahn_kalai_threshold)
            self.assertGreaterEqual(res.mesh_correlation_length_nm, 120.0)
            self.assertGreaterEqual(res.nanoparticle_diffusivity_ratio, 0.25)
            self.assertLessEqual(res.modulus_ratio_g_prime_to_g_double_prime, 1.20)
            self.assertLessEqual(res.phi_mucin, 0.035)

    def test_negative_control_falsification(self) -> None:
        normal = evaluate_falsification_negative_control(poison_negative_mesh=False)
        self.assertTrue(normal["passed"])
        self.assertEqual(normal["rule"], "AH-10")
        poisoned = evaluate_falsification_negative_control(poison_negative_mesh=True)
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
