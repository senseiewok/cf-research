#!/usr/bin/env python3
"""Unit tests for Round 185: Inhaled Lentiviral Vector Mucus Penetration & Basal Stem Cell Integration Bridge (T-0862)."""
from __future__ import annotations

import unittest

from tools.discovery.round185_lentivirus_mucus_penetration_integration_bridge import (
    FIVE_LAYER_LADDER,
    PARAMETER_PROVENANCE,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    complex_step_asl_derivative,
    compute_asl_height_um,
    compute_basal_stem_cell_integration_fraction,
    compute_effective_mucus_diffusivity_um2_s,
    compute_mcc_velocity_mm_min,
    compute_mucus_penetration_fraction,
    compute_receptor_binding_fraction,
    compute_transepithelial_current_pct_wt,
    compute_unhindered_diffusivity_um2_s,
    evaluate_falsification_negative_control,
    evaluate_round185_lentivirus_mucus_penetration_integration,
)


class Round185LentivirusMucusPenetrationIntegrationTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round185_lentivirus_mucus_penetration_integration(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Invariants failed for {key}")
            self.assertGreaterEqual(res.unhindered_diffusivity_um2_s, 1.0)
            self.assertGreaterEqual(res.effective_mucus_diffusivity_um2_s, 0.0)
            self.assertGreaterEqual(res.mucus_penetration_fraction, 0.0)
            self.assertLessEqual(res.mucus_penetration_fraction, 1.0)
            self.assertGreaterEqual(res.receptor_binding_fraction, 0.0)
            self.assertLessEqual(res.receptor_binding_fraction, 1.0)
            self.assertGreaterEqual(res.basal_stem_cell_integration_fraction, 0.0)
            self.assertLessEqual(res.basal_stem_cell_integration_fraction, 1.0)
            self.assertGreaterEqual(res.transepithelial_current_pct_wt, 0.0)
            self.assertLessEqual(res.transepithelial_current_pct_wt, 100.0)
            self.assertGreaterEqual(res.airway_surface_liquid_height_um, 3.5)
            self.assertLessEqual(res.airway_surface_liquid_height_um, 8.5)
            self.assertGreaterEqual(res.mucociliary_transport_velocity_mm_min, 0.0)
            self.assertLessEqual(res.complex_step_relative_error, 1e-6)

    def test_bi_3720931_high_dose_crosses_kahn_kalai(self) -> None:
        res = evaluate_round185_lentivirus_mucus_penetration_integration("bi_3720931_therapeutic_high_dose_cf_airway")
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res.mucus_penetration_fraction, 0.15)
        self.assertGreaterEqual(res.basal_stem_cell_integration_fraction, 0.08)
        self.assertGreaterEqual(res.transepithelial_current_pct_wt, 7.0)
        self.assertGreaterEqual(res.airway_surface_liquid_height_um, 7.0)
        self.assertGreaterEqual(res.mucociliary_transport_velocity_mm_min, 1.0)

    def test_mucolytic_pretreated_crosses_kahn_kalai(self) -> None:
        res = evaluate_round185_lentivirus_mucus_penetration_integration("bi_3720931_mucolytic_pretreated_airway")
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res.mucus_penetration_fraction, 0.50)
        self.assertGreaterEqual(res.basal_stem_cell_integration_fraction, 0.15)
        self.assertGreaterEqual(res.airway_surface_liquid_height_um, 7.0)
        self.assertGreaterEqual(res.mucociliary_transport_velocity_mm_min, 2.0)

    def test_sub_curative_presets_fail_kahn_kalai(self) -> None:
        # Hyperconcentrated mucus barrier: mesh sterically blocks vector
        mucus_blocked = evaluate_round185_lentivirus_mucus_penetration_integration("hyperconcentrated_mucus_untreated_cf_barrier")
        self.assertFalse(mucus_blocked.crosses_kahn_kalai_threshold)
        self.assertLess(mucus_blocked.mucus_penetration_fraction, 0.10)
        self.assertLess(mucus_blocked.basal_stem_cell_integration_fraction, 0.05)
        self.assertEqual(mucus_blocked.mucociliary_transport_velocity_mm_min, 0.0)

        # Low dose sub-therapeutic
        low_dose = evaluate_round185_lentivirus_mucus_penetration_integration("low_dose_subtherapeutic_lentivirus")
        self.assertFalse(low_dose.crosses_kahn_kalai_threshold)
        self.assertLess(low_dose.receptor_binding_fraction, 0.15)
        self.assertLess(low_dose.basal_stem_cell_integration_fraction, 0.05)

        # Untreated CF baseline
        untreated = evaluate_round185_lentivirus_mucus_penetration_integration("untreated_cf_baseline_zero_vector")
        self.assertFalse(untreated.crosses_kahn_kalai_threshold)
        self.assertEqual(untreated.basal_stem_cell_integration_fraction, 0.0)
        self.assertAlmostEqual(untreated.airway_surface_liquid_height_um, 3.5, places=2)
        self.assertEqual(untreated.mucociliary_transport_velocity_mm_min, 0.0)

    def test_complex_step_derivative_accuracy(self) -> None:
        cs_deriv, fd_deriv, rel_err = complex_step_asl_derivative(basal_integration_fraction=0.15)
        self.assertGreater(cs_deriv, 0.0)
        self.assertGreater(fd_deriv, 0.0)
        self.assertLessEqual(rel_err, 1e-6)

    def test_poisoned_negative_controls_fail_deterministically(self) -> None:
        res_conc = evaluate_falsification_negative_control(poison_negative_conc=True)
        self.assertTrue(res_conc["passed"])
        self.assertIn("Correctly caught", res_conc["error"])

        res_mesh = evaluate_falsification_negative_control(poison_negative_mesh=True)
        self.assertTrue(res_mesh["passed"])
        self.assertIn("Correctly caught", res_mesh["error"])

        res_thick = evaluate_falsification_negative_control(poison_negative_thickness=True)
        self.assertTrue(res_thick["passed"])
        self.assertIn("Correctly caught", res_thick["error"])

        clean = evaluate_falsification_negative_control()
        self.assertTrue(clean["passed"])
        self.assertIsNone(clean["error"])

    def test_monotone_mesh_diffusivity_scaling(self) -> None:
        d0 = compute_unhindered_diffusivity_um2_s()
        d_small = compute_effective_mucus_diffusivity_um2_s(d0, mesh_pore_nm=50.0)
        d_med = compute_effective_mucus_diffusivity_um2_s(d0, mesh_pore_nm=100.0)
        d_large = compute_effective_mucus_diffusivity_um2_s(d0, mesh_pore_nm=200.0)
        self.assertLess(d_small, d_med)
        self.assertLess(d_med, d_large)
        self.assertLess(d_large, d0)

    def test_provenance_and_two_column_epistemic_hygiene_metadata(self) -> None:
        self.assertGreaterEqual(len(PRIMARY_CITATIONS), 4)
        self.assertGreaterEqual(len(PARAMETER_PROVENANCE), 5)
        for param, meta in PARAMETER_PROVENANCE.items():
            self.assertIn("value", meta)
            self.assertIn("tag", meta)
            self.assertIn("source", meta)
            self.assertIn("bounds", meta)
            self.assertIn(meta["tag"], ("empirical", "derived", "phenomenological_heuristic"))

        self.assertIn("deductive_theorem", TWO_COLUMN_EPISTEMIC_HYGIENE)
        self.assertIn("biophysical_hypothesis", TWO_COLUMN_EPISTEMIC_HYGIENE)
        self.assertIn("statement", TWO_COLUMN_EPISTEMIC_HYGIENE["biophysical_hypothesis"])
        self.assertIn("falsification_condition", TWO_COLUMN_EPISTEMIC_HYGIENE["biophysical_hypothesis"])

        self.assertIn("cool_human_fact", FIVE_LAYER_LADDER)
        self.assertIn("everyday_analogy", FIVE_LAYER_LADDER)
        self.assertIn("engineering_lens", FIVE_LAYER_LADDER)
        self.assertIn("formal_theorem", FIVE_LAYER_LADDER)
        self.assertIn("falsifiable_code_and_toy", FIVE_LAYER_LADDER)


if __name__ == "__main__":
    unittest.main()
