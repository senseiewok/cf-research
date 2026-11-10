#!/usr/bin/env python3
"""Unit tests for Round 173: Ogston-Brinkman Steric Sieving & Hydrodynamic Hindrance Invariant Bridge (T-0841)."""
from __future__ import annotations

import unittest

from tools.discovery.round173_ogston_brinkman_nanoparticle_sieving_bridge import (
    FIVE_LAYER_LADDER,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    complex_step_ogston_derivative,
    compute_brinkman_hydrodynamic_factor,
    compute_clearance_escape_probability,
    compute_effective_diffusivity_ratio,
    compute_mucus_translocation_time,
    compute_ogston_steric_partition,
    compute_stokes_einstein_d0,
    compute_surface_shielding_factor,
    evaluate_falsification_negative_control,
    evaluate_round173_ogston_brinkman_sieving,
)


class Round173OgstonBrinkmanTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round173_ogston_brinkman_sieving(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Failed for {key}")
            self.assertGreaterEqual(res.particle_diameter_nm, 20.0)
            self.assertLessEqual(res.particle_diameter_nm, 500.0)
            self.assertGreater(res.mucin_volume_fraction, 0.0)
            self.assertLess(res.mucin_volume_fraction, 1.0)
            self.assertGreaterEqual(res.ogston_steric_partition, 0.0)
            self.assertLessEqual(res.ogston_steric_partition, 1.0)
            self.assertGreaterEqual(res.brinkman_hydrodynamic_factor, 0.0)
            self.assertLessEqual(res.brinkman_hydrodynamic_factor, 1.0)
            self.assertGreaterEqual(res.surface_adhesive_shielding, 0.0)
            self.assertLessEqual(res.surface_adhesive_shielding, 1.0)
            self.assertGreaterEqual(res.effective_diffusivity_ratio, 0.0)
            self.assertLessEqual(res.effective_diffusivity_ratio, 1.0)
            self.assertGreater(res.unhindered_diffusivity_m2_s, 0.0)
            self.assertGreaterEqual(res.effective_diffusivity_m2_s, 0.0)
            self.assertGreater(res.mucus_translocation_time_s, 0.0)
            self.assertGreaterEqual(res.clearance_escape_probability, 0.0)
            self.assertLessEqual(res.clearance_escape_probability, 1.0)
            self.assertLessEqual(res.complex_step_relative_error, 1e-10)

    def test_dense_peg_65nm_healthy_mucus(self) -> None:
        res = evaluate_round173_ogston_brinkman_sieving("dense_peg_65nm_healthy_mucus")
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res.effective_diffusivity_ratio, 0.20)
        self.assertLess(res.mucus_translocation_time_s, 100.0)
        self.assertGreaterEqual(res.clearance_escape_probability, 0.85)

    def test_dense_peg_65nm_cf_mucus(self) -> None:
        res = evaluate_round173_ogston_brinkman_sieving("dense_peg_65nm_cf_mucus")
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res.effective_diffusivity_ratio, 0.08)
        self.assertLess(res.mucus_translocation_time_s, 300.0)
        self.assertGreaterEqual(res.clearance_escape_probability, 0.60)

    def test_standard_95nm_cf_mucus(self) -> None:
        res = evaluate_round173_ogston_brinkman_sieving("standard_95nm_cf_mucus")
        self.assertFalse(res.crosses_kahn_kalai_threshold)
        self.assertLess(res.effective_diffusivity_ratio, 0.05)
        self.assertGreater(res.mucus_translocation_time_s, 600.0)
        self.assertLess(res.clearance_escape_probability, 0.40)

    def test_unshielded_cationic_135nm_cf_mucus(self) -> None:
        res = evaluate_round173_ogston_brinkman_sieving("unshielded_cationic_135nm_cf_mucus")
        self.assertFalse(res.crosses_kahn_kalai_threshold)
        self.assertLess(res.effective_diffusivity_ratio, 0.001)
        self.assertGreater(res.mucus_translocation_time_s, 10000.0)
        self.assertLess(res.clearance_escape_probability, 0.01)

    def test_oversized_200nm_healthy_mucus(self) -> None:
        res = evaluate_round173_ogston_brinkman_sieving("oversized_200nm_healthy_mucus")
        self.assertFalse(res.crosses_kahn_kalai_threshold)
        self.assertLess(res.effective_diffusivity_ratio, 0.05)
        self.assertGreater(res.mucus_translocation_time_s, 600.0)
        self.assertLess(res.clearance_escape_probability, 0.40)

    def test_complex_step_ogston_derivative(self) -> None:
        cs_deriv, analytical, rel_err = complex_step_ogston_derivative(
            r_p_nm=32.5,
            phi_muc=0.040,
        )
        self.assertLessEqual(rel_err, 1e-10)
        self.assertLess(cs_deriv, 0.0)  # Steric transmission decreases with particle radius

    def test_negative_control_falsification(self) -> None:
        neg1 = evaluate_falsification_negative_control(poison_negative_radius=True)
        self.assertFalse(neg1["passed"])
        self.assertIn("poisoned", neg1["error"])

        neg2 = evaluate_falsification_negative_control(poison_negative_volume_fraction=True)
        self.assertFalse(neg2["passed"])

        neg3 = evaluate_falsification_negative_control(poison_excessive_volume_fraction=True)
        self.assertFalse(neg3["passed"])

        neg4 = evaluate_falsification_negative_control(poison_negative_viscosity=True)
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
