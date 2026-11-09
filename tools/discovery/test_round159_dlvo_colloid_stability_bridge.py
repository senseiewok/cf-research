#!/usr/bin/env python3
"""Unit tests for Round 159: DLVO Colloid Stability & Hamaker Van der Waals Invariant Bridge (T-0824)."""
from __future__ import annotations

import math
import unittest

from tools.discovery.round159_dlvo_colloid_stability_bridge import (
    FIVE_LAYER_LADDER,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    complex_step_colloidal_force,
    compute_debye_length_nm,
    compute_edl_potential_j,
    compute_steric_potential_j,
    compute_total_dlvo_potential_kbt,
    compute_vdw_potential_j,
    evaluate_falsification_negative_control,
    evaluate_round159_dlvo_colloid_stability,
    find_energy_barrier_kbt,
)


class Round159DLVOColloidStabilityTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round159_dlvo_colloid_stability(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Failed for {key}")
            self.assertGreater(res.debye_length_nm, 0.0)
            self.assertLessEqual(res.complex_step_relative_error, 1e-10)

    def test_dense_brush_mucus_penetrating_lnp(self) -> None:
        res = evaluate_round159_dlvo_colloid_stability("dense_brush_mucus_penetrating_lnp")
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res.energy_barrier_kbt, 18.0)
        self.assertGreaterEqual(res.mucus_diffusion_ratio, 0.55)
        self.assertGreaterEqual(res.fuchs_stability_ratio_log10, 5.0)

    def test_standard_peg_shielded_monomeric_lnp(self) -> None:
        res = evaluate_round159_dlvo_colloid_stability("standard_peg_shielded_monomeric_lnp")
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res.energy_barrier_kbt, 15.0)
        self.assertGreaterEqual(res.mucus_diffusion_ratio, 0.40)

    def test_unshielded_cationic_lnp_rapid_aggregation(self) -> None:
        res = evaluate_round159_dlvo_colloid_stability("unshielded_cationic_lnp_rapid_aggregation")
        self.assertFalse(res.crosses_kahn_kalai_threshold)
        self.assertLess(res.mucus_diffusion_ratio, 0.05)

    def test_cf_hyperconcentrated_mucus_entrapment(self) -> None:
        res = evaluate_round159_dlvo_colloid_stability("cf_hyperconcentrated_mucus_entrapment")
        self.assertFalse(res.crosses_kahn_kalai_threshold)
        self.assertLess(res.mucus_diffusion_ratio, 0.10)

    def test_kahn_kalai_colloidal_penetration_boundary(self) -> None:
        res = evaluate_round159_dlvo_colloid_stability("kahn_kalai_colloidal_penetration_boundary")
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res.energy_barrier_kbt, 12.0)
        self.assertGreaterEqual(res.mucus_diffusion_ratio, 0.25)

    def test_debye_length_ionic_strength_dependence(self) -> None:
        # Higher ionic strength -> shorter Debye screening length
        debye_low_i = compute_debye_length_nm(50.0)
        debye_phys_i = compute_debye_length_nm(150.0)
        debye_high_i = compute_debye_length_nm(300.0)
        self.assertGreater(debye_low_i, debye_phys_i)
        self.assertGreater(debye_phys_i, debye_high_i)
        # At 150 mM, Debye length should be approx 0.75 - 0.85 nm
        self.assertAlmostEqual(debye_phys_i, 0.78, delta=0.08)

    def test_vdw_potential_attractive_divergence(self) -> None:
        v_1nm = compute_vdw_potential_j(1.0, 45.0)
        v_2nm = compute_vdw_potential_j(2.0, 45.0)
        # Both must be negative and |V(1nm)| > |V(2nm)|
        self.assertLess(v_1nm, 0.0)
        self.assertLess(v_2nm, 0.0)
        self.assertLess(v_1nm, v_2nm)

    def test_steric_potential_brush_cutoff(self) -> None:
        brush_l = 3.5  # nm -> contact cutoff 2 * 3.5 = 7.0 nm
        v_inside = compute_steric_potential_j(4.0, 45.0, 0.35, brush_l)
        v_outside = compute_steric_potential_j(7.5, 45.0, 0.35, brush_l)
        self.assertGreater(v_inside, 0.0)
        self.assertEqual(v_outside, 0.0)

    def test_complex_step_force_sensitivity(self) -> None:
        force, dV_ds, rel_err = complex_step_colloidal_force(
            separation_nm=2.0,
            radius_nm=45.0,
            zeta_potential_mv=-2.5,
            peg_density_chains_nm2=0.35,
            peg_brush_length_nm=3.5,
        )
        self.assertLessEqual(rel_err, 1e-10)

    def test_negative_control_falsification(self) -> None:
        clean = evaluate_falsification_negative_control()
        self.assertTrue(clean["passed"])
        self.assertEqual(clean["rule"], "AH-10")

        poisoned_temp = evaluate_falsification_negative_control(poison_negative_temperature=True)
        self.assertFalse(poisoned_temp["passed"])

        poisoned_hamaker = evaluate_falsification_negative_control(poison_negative_hamaker=True)
        self.assertFalse(poisoned_hamaker["passed"])

        poisoned_radius = evaluate_falsification_negative_control(poison_negative_radius=True)
        self.assertFalse(poisoned_radius["passed"])

        poisoned_sep = evaluate_falsification_negative_control(poison_negative_separation=True)
        self.assertFalse(poisoned_sep["passed"])

    def test_two_column_epistemic_hygiene_and_five_layer_ladder(self) -> None:
        self.assertIn("deductive_theorem", TWO_COLUMN_EPISTEMIC_HYGIENE)
        self.assertIn("biophysical_hypothesis", TWO_COLUMN_EPISTEMIC_HYGIENE)
        self.assertIn("falsification_condition", TWO_COLUMN_EPISTEMIC_HYGIENE["biophysical_hypothesis"])
        self.assertIn("cool_human_fact", FIVE_LAYER_LADDER)
        self.assertIn("everyday_analogy", FIVE_LAYER_LADDER)
        self.assertIn("engineering_lens", FIVE_LAYER_LADDER)
        self.assertIn("formal_theorem", FIVE_LAYER_LADDER)
        self.assertIn("falsifiable_code_and_toy", FIVE_LAYER_LADDER)
        self.assertGreaterEqual(len(PRIMARY_CITATIONS), 3)


if __name__ == "__main__":
    unittest.main()
