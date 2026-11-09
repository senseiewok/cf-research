#!/usr/bin/env python3
"""Unit tests for Round 168: Mason-Weitz Generalized Stokes-Einstein Nanorheology & Complex Modulus Invariant Bridge (T-0833)."""
from __future__ import annotations

import unittest

from tools.discovery.round168_mason_weitz_nanorheology_bridge import (
    FIVE_LAYER_LADDER,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    complex_step_modulus_derivative,
    compute_effective_diffusivity_ratio,
    compute_mason_weitz_moduli,
    compute_mucin_mesh_pore_size,
    evaluate_falsification_negative_control,
    evaluate_round168_mason_weitz_nanorheology,
)


class Round168NanorheologyTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round168_mason_weitz_nanorheology(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Failed for {key}")
            self.assertGreaterEqual(res.subdiffusion_exponent_alpha, 0.0)
            self.assertLessEqual(res.subdiffusion_exponent_alpha, 1.0)
            self.assertGreaterEqual(res.effective_diffusivity_ratio, 0.0)
            self.assertLessEqual(res.effective_diffusivity_ratio, 1.0)
            self.assertGreaterEqual(res.storage_modulus_g_prime_pa, 0.0)
            self.assertGreaterEqual(res.loss_modulus_g_double_prime_pa, 0.0)
            self.assertLessEqual(res.complex_step_relative_error, 1e-10)

    def test_healthy_airway_mucus_gser(self) -> None:
        res = evaluate_round168_mason_weitz_nanorheology("healthy_airway_mucus_gser")
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res.subdiffusion_exponent_alpha, 0.70)
        self.assertGreaterEqual(res.effective_diffusivity_ratio, 0.70)
        self.assertGreaterEqual(res.mesh_pore_size_nm, 140.0)
        self.assertGreater(res.loss_modulus_g_double_prime_pa, res.storage_modulus_g_prime_pa)

    def test_dornase_alfa_pulmozyme_digested_cf_mucus(self) -> None:
        res = evaluate_round168_mason_weitz_nanorheology("dornase_alfa_pulmozyme_digested_cf_mucus")
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res.subdiffusion_exponent_alpha, 0.55)
        self.assertGreaterEqual(res.effective_diffusivity_ratio, 0.45)
        self.assertGreaterEqual(res.mesh_pore_size_nm, 100.0)

    def test_kahn_kalai_nanorheology_penetration_boundary(self) -> None:
        res = evaluate_round168_mason_weitz_nanorheology("kahn_kalai_nanorheology_penetration_boundary")
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res.subdiffusion_exponent_alpha, 0.50)
        self.assertGreaterEqual(res.effective_diffusivity_ratio, 0.15)

    def test_untreated_purulent_cf_mucus_steric_entrapment(self) -> None:
        res = evaluate_round168_mason_weitz_nanorheology("untreated_purulent_cf_mucus_steric_entrapment")
        self.assertFalse(res.crosses_kahn_kalai_threshold)
        self.assertLess(res.subdiffusion_exponent_alpha, 0.30)
        self.assertLess(res.effective_diffusivity_ratio, 0.05)
        self.assertLess(res.mesh_pore_size_nm, 65.0)

    def test_dilute_saline_hypertonic_flushed_asl(self) -> None:
        res = evaluate_round168_mason_weitz_nanorheology("dilute_saline_hypertonic_flushed_asl")
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res.subdiffusion_exponent_alpha, 0.80)
        self.assertGreaterEqual(res.effective_diffusivity_ratio, 0.85)

    def test_complex_step_modulus_derivative(self) -> None:
        cs_deriv, analytical, rel_err = complex_step_modulus_derivative(
            omega_rad_s=1.00,
            plateau_modulus_pa=1.00,
            alpha=0.78,
        )
        self.assertLessEqual(rel_err, 1e-10)
        self.assertGreater(cs_deriv, 0.0)

    def test_negative_control_falsification(self) -> None:
        neg1 = evaluate_falsification_negative_control(poison_negative_probe_diameter=True)
        self.assertFalse(neg1["passed"])
        self.assertIn("poisoned", neg1["error"])

        neg2 = evaluate_falsification_negative_control(poison_negative_plateau_modulus=True)
        self.assertFalse(neg2["passed"])

        neg3 = evaluate_falsification_negative_control(poison_negative_frequency=True)
        self.assertFalse(neg3["passed"])

        neg4 = evaluate_falsification_negative_control(poison_invalid_alpha=True)
        self.assertFalse(neg4["passed"])

        clean = evaluate_falsification_negative_control()
        self.assertTrue(clean["passed"])

    def test_epistemic_hygiene_and_ladder_structure(self) -> None:
        self.assertIn("deductive_theorem", TWO_COLUMN_EPISTEMIC_HYGIENE)
        self.assertIn("biophysical_hypothesis", TWO_COLUMN_EPISTEMIC_HYGIENE)
        self.assertGreaterEqual(len(FIVE_LAYER_LADDER), 5)
        self.assertGreaterEqual(len(PRIMARY_CITATIONS), 4)


if __name__ == "__main__":
    unittest.main()
