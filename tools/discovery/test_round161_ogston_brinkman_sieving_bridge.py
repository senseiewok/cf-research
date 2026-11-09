#!/usr/bin/env python3
"""Unit tests for Round 161: Ogston-Brinkman Steric Mesh Sieving & Hydrodynamic Tortuosity Invariant Bridge (T-0826)."""
from __future__ import annotations

import math
import unittest

from tools.discovery.round161_ogston_brinkman_sieving_bridge import (
    FIVE_LAYER_LADDER,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    complex_step_diffusivity_derivative,
    compute_brinkman_hydraulic_permeability,
    compute_brinkman_hydrodynamic_factor,
    compute_effective_diffusivity_ratio,
    compute_fiber_tortuosity,
    compute_mucus_transit_probability,
    compute_ogston_partition_coefficient,
    compute_pcl_brush_entry_fraction,
    evaluate_falsification_negative_control,
    evaluate_round161_ogston_brinkman_sieving,
)


class Round161OgstonBrinkmanSievingTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round161_ogston_brinkman_sieving(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Failed for {key}")
            self.assertGreater(res.particle_radius_nm, 0.0)
            self.assertLessEqual(res.complex_step_relative_error, 1e-10)

    def test_compact_sort_lnp_relaxed_mucus(self) -> None:
        res = evaluate_round161_ogston_brinkman_sieving("compact_sort_lnp_relaxed_mucus")
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res.effective_diffusivity_ratio, 0.45)
        self.assertGreaterEqual(res.overall_transfection_index, 0.20)
        self.assertGreaterEqual(res.ogston_partition_coefficient, 0.70)

    def test_standard_mrna_lnp_healthy_airway(self) -> None:
        res = evaluate_round161_ogston_brinkman_sieving("standard_mrna_lnp_healthy_airway")
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res.effective_diffusivity_ratio, 0.30)
        self.assertGreaterEqual(res.overall_transfection_index, 0.12)

    def test_kahn_kalai_steric_filtration_boundary(self) -> None:
        res = evaluate_round161_ogston_brinkman_sieving("kahn_kalai_steric_filtration_boundary")
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res.effective_diffusivity_ratio, 0.15)
        self.assertGreaterEqual(res.overall_transfection_index, 0.118)

    def test_cf_dense_sputum_steric_mesh_cutoff(self) -> None:
        res = evaluate_round161_ogston_brinkman_sieving("cf_dense_sputum_steric_mesh_cutoff")
        self.assertFalse(res.crosses_kahn_kalai_threshold)
        self.assertLess(res.effective_diffusivity_ratio, 0.10)
        self.assertLess(res.overall_transfection_index, 0.02)

    def test_oversized_lnp_periciliary_brush_exclusion(self) -> None:
        res = evaluate_round161_ogston_brinkman_sieving("oversized_lnp_periciliary_brush_exclusion")
        self.assertFalse(res.crosses_kahn_kalai_threshold)
        self.assertLess(res.effective_diffusivity_ratio, 0.01)
        self.assertLess(res.overall_transfection_index, 0.001)

    def test_ogston_super_exponential_scaling_with_particle_radius(self) -> None:
        mesh_pore = 100.0  # nm
        f_small = compute_ogston_partition_coefficient(20.0, mesh_pore)
        f_med = compute_ogston_partition_coefficient(40.0, mesh_pore)
        f_large = compute_ogston_partition_coefficient(80.0, mesh_pore)
        self.assertGreater(f_small, f_med)
        self.assertGreater(f_med, f_large)
        self.assertGreater(f_small, 0.50)
        self.assertLess(f_large, 0.20)

    def test_brinkman_hydrodynamic_factor_decay(self) -> None:
        kp = 1000.0  # nm^2
        f_hyd_small = compute_brinkman_hydrodynamic_factor(10.0, kp)
        f_hyd_large = compute_brinkman_hydrodynamic_factor(100.0, kp)
        self.assertGreater(f_hyd_small, f_hyd_large)
        self.assertLess(f_hyd_small, 1.0)
        self.assertGreater(f_hyd_large, 0.0)

    def test_fiber_tortuosity_monotonic_increase(self) -> None:
        t_low = compute_fiber_tortuosity(0.01)
        t_high = compute_fiber_tortuosity(0.10)
        self.assertGreater(t_high, t_low)
        self.assertGreaterEqual(t_low, 1.0)

    def test_complex_step_diffusivity_sensitivity(self) -> None:
        deriv_cs, analytical, rel_err = complex_step_diffusivity_derivative(
            particle_radius_nm=30.0,
            mesh_pore_nm=120.0,
            fiber_volume_fraction=0.02,
        )
        self.assertLessEqual(rel_err, 1e-10)
        # Derivative must be negative (larger particle -> lower D_eff)
        self.assertLess(deriv_cs, 0.0)

    def test_negative_control_falsification(self) -> None:
        clean = evaluate_falsification_negative_control()
        self.assertTrue(clean["passed"])
        self.assertEqual(clean["rule"], "AH-10")

        poisoned_r = evaluate_falsification_negative_control(poison_negative_particle_radius=True)
        self.assertFalse(poisoned_r["passed"])

        poisoned_pore = evaluate_falsification_negative_control(poison_negative_mesh_pore=True)
        self.assertFalse(poisoned_pore["passed"])

        poisoned_fiber = evaluate_falsification_negative_control(poison_negative_fiber_radius=True)
        self.assertFalse(poisoned_fiber["passed"])

        poisoned_frac = evaluate_falsification_negative_control(poison_negative_fiber_fraction=True)
        self.assertFalse(poisoned_frac["passed"])

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
