#!/usr/bin/env python3
"""Unit tests for Round 152: Onsager Reciprocal Dissipation Matrix & Hydrodynamic Streaming Potential Coupling Invariant Bridge (T-0807)."""
from __future__ import annotations

import unittest

from tools.discovery.round152_onsager_reciprocal_streaming_bridge import (
    FIVE_LAYER_LADDER,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    compute_entropy_production_rate,
    compute_onsager_current_density,
    compute_onsager_volume_flux,
    evaluate_falsification_negative_control,
    evaluate_round152_onsager_reciprocal_streaming,
)


class Round152OnsagerStreamingTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round152_onsager_reciprocal_streaming(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Failed for {key}")
            self.assertGreater(res.hydraulic_permeability_l11, 0.0)
            self.assertGreater(res.electrical_conductance_l22, 0.0)
            self.assertEqual(res.onsager_reciprocal_symmetry_error, 0.0)
            self.assertGreaterEqual(res.dissipation_matrix_determinant, 0.0)
            self.assertLessEqual(abs(res.thermodynamic_coupling_efficiency), 1.0)
            self.assertGreaterEqual(res.dissipation_rate_w_m2, -1e-12)

    def test_wild_type_ciliated_epithelial_transport(self) -> None:
        wt = evaluate_round152_onsager_reciprocal_streaming("wild_type_ciliated_epithelial_transport")
        self.assertTrue(wt.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(wt.electrical_conductance_l22, 15.0)
        self.assertGreater(wt.volume_flux_um_s, 0.040)  # Positive fluid secretion
        self.assertGreaterEqual(wt.predicted_pcl_height_um, 6.5)  # Hydrated PCL propped open
        self.assertAlmostEqual(wt.fractional_wild_type_flux, 1.0, delta=0.05)
        self.assertGreaterEqual(wt.current_density_ua_cm2, 30.0)

    def test_cf_defects_and_fluid_absorption_collapse(self) -> None:
        df508 = evaluate_round152_onsager_reciprocal_streaming("cf_deltaf508_transport_collapse")
        g551d = evaluate_round152_onsager_reciprocal_streaming("cf_g551d_gating_dehydration_stasis")
        for mutant in (df508, g551d):
            self.assertFalse(mutant.crosses_kahn_kalai_threshold)
            self.assertLess(mutant.volume_flux_um_s, 0.0)  # Pathological fluid absorption
            self.assertLess(mutant.predicted_pcl_height_um, 4.5)  # Depleted PCL
            self.assertLessEqual(mutant.fractional_wild_type_flux, 0.0)
            self.assertLessEqual(mutant.current_density_ua_cm2, 5.0)

    def test_kahn_kalai_onsager_fluid_boundary(self) -> None:
        kk = evaluate_round152_onsager_reciprocal_streaming("kahn_kalai_onsager_fluid_boundary")
        self.assertTrue(kk.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(kk.electrical_conductance_l22, 4.0)
        self.assertGreaterEqual(kk.volume_flux_um_s, 0.006)
        self.assertGreaterEqual(kk.predicted_pcl_height_um, 5.0)
        self.assertGreaterEqual(kk.fractional_wild_type_flux, 0.118)

    def test_trikafta_onsager_rehydration_rescue(self) -> None:
        rescued = evaluate_round152_onsager_reciprocal_streaming("trikafta_onsager_rehydration_rescue")
        self.assertTrue(rescued.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(rescued.electrical_conductance_l22, 14.0)
        self.assertGreaterEqual(rescued.volume_flux_um_s, 0.040)
        self.assertGreaterEqual(rescued.predicted_pcl_height_um, 6.2)
        self.assertGreaterEqual(rescued.fractional_wild_type_flux, 0.85)

    def test_negative_control_falsification(self) -> None:
        normal = evaluate_falsification_negative_control(poison_broken_onsager_symmetry=False)
        self.assertTrue(normal["passed"])
        self.assertEqual(normal["rule"], "AH-10")

        broken_sym = evaluate_falsification_negative_control(poison_broken_onsager_symmetry=True)
        self.assertFalse(broken_sym["passed"])

        neg_det = evaluate_falsification_negative_control(poison_negative_determinant=True)
        self.assertFalse(neg_det["passed"])

    def test_onsager_flux_linearity_and_coupling(self) -> None:
        l11 = 2.5e-10
        l22 = 4.0
        l12 = 1.2e-6
        l21 = 1.2e-6

        # Pure pressure drive (Delta V = 0)
        jv_p = compute_onsager_volume_flux(l11, l12, delta_p_pa=100.0, delta_v_v=0.0)
        self.assertAlmostEqual(jv_p, l11 * 100.0, places=12)

        # Pure electric drive (Delta P = 0)
        jv_v = compute_onsager_volume_flux(l11, l12, delta_p_pa=0.0, delta_v_v=0.05)
        self.assertAlmostEqual(jv_v, l12 * 0.05, places=12)

        # Streaming current reciprocity
        curr_p = compute_onsager_current_density(l21, l22, delta_p_pa=100.0, delta_v_v=0.0)
        self.assertAlmostEqual(curr_p, l21 * 100.0, places=12)

        # Dissipation non-negativity
        diss = compute_entropy_production_rate(l11, l12, l21, l22, delta_p=100.0, delta_v=0.05)
        self.assertGreater(diss, 0.0)

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
