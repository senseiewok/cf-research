#!/usr/bin/env python3
"""Unit tests for Round 199: Basal Stem Cell Notch/Wnt Signaling & Multiciliated Lineage Differentiation Bridge (T-0888)."""
from __future__ import annotations

import unittest

from tools.discovery.round199_basal_stem_cell_differentiation_bridge import (
    FIVE_LAYER_LADDER,
    PARAMETER_PROVENANCE,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    complex_step_wnt_mcc_derivative,
    compute_lineage_simplex_fractions,
    compute_percolated_mcc_velocity,
    compute_surface_cftr_and_conductance,
    evaluate_falsification_negative_control,
    evaluate_round199_basal_stem_cell_differentiation,
)


class Round199BasalStemCellDifferentiationBridgeTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_simplex_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round199_basal_stem_cell_differentiation(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Invariants failed for {key}")
            self.assertAlmostEqual(res.f_basal + res.f_ciliated + res.f_goblet, 1.0, places=12)
            self.assertGreaterEqual(res.f_basal, 0.0)
            self.assertGreaterEqual(res.f_ciliated, 0.0)
            self.assertGreaterEqual(res.f_goblet, 0.0)
            self.assertGreaterEqual(res.i_sc_equiv_pct_wt, 0.0)
            self.assertLessEqual(res.i_sc_equiv_pct_wt, 100.0)
            self.assertLessEqual(res.complex_step_relative_error, 1e-6)

    def test_curative_presets_cross_kahn_kalai(self) -> None:
        res_cur = evaluate_round199_basal_stem_cell_differentiation("notch_modulated_wnt_edited_basal_curative")
        self.assertTrue(res_cur.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res_cur.f_ciliated, 0.40)
        self.assertLessEqual(res_cur.f_goblet, 0.30)
        self.assertGreaterEqual(res_cur.i_sc_equiv_pct_wt, 25.0)
        self.assertGreaterEqual(res_cur.v_mcc_mm_min, 1.50)

        res_high = evaluate_round199_basal_stem_cell_differentiation("high_efficiency_basal_engraftment_curative")
        self.assertTrue(res_high.crosses_kahn_kalai_threshold)
        self.assertGreater(res_high.v_mcc_mm_min, res_cur.v_mcc_mm_min)

    def test_pathology_and_negative_controls_fail_kahn_kalai(self) -> None:
        res_gob = evaluate_round199_basal_stem_cell_differentiation("chronic_notch_goblet_hyperplasia_cf")
        self.assertFalse(res_gob.crosses_kahn_kalai_threshold)
        self.assertLess(res_gob.f_ciliated, 0.10)
        self.assertGreater(res_gob.f_goblet, 0.70)
        self.assertLess(res_gob.v_mcc_mm_min, 0.05)

        res_unedit = evaluate_round199_basal_stem_cell_differentiation("unedited_cf_basal_progenitors")
        self.assertFalse(res_unedit.crosses_kahn_kalai_threshold)
        self.assertEqual(res_unedit.i_sc_equiv_pct_wt, 0.0)
        self.assertEqual(res_unedit.v_mcc_mm_min, 0.0)

        res_wnt = evaluate_round199_basal_stem_cell_differentiation("wnt_deficient_squamous_basal_arrest")
        self.assertFalse(res_wnt.crosses_kahn_kalai_threshold)
        self.assertLess(res_wnt.f_ciliated, 0.15)
        self.assertLess(res_wnt.v_mcc_mm_min, 0.05)

    def test_notch_and_wnt_monotonicity_invariants(self) -> None:
        _, fc_low_notch, fg_low_notch = compute_lineage_simplex_fractions(0.20, 1.0, 0.1)
        _, fc_high_notch, fg_high_notch = compute_lineage_simplex_fractions(1.20, 1.0, 0.1)
        self.assertLess(float(fc_high_notch), float(fc_low_notch))
        self.assertGreater(float(fg_high_notch), float(fg_low_notch))

        _, fc_low_wnt, _ = compute_lineage_simplex_fractions(0.30, 0.25, 0.1)
        _, fc_high_wnt, _ = compute_lineage_simplex_fractions(0.30, 1.25, 0.1)
        self.assertGreater(float(fc_high_wnt), float(fc_low_wnt))

        _, isc_low = compute_surface_cftr_and_conductance(0.20, 0.15, 0.65, 0.20)
        _, isc_high = compute_surface_cftr_and_conductance(0.70, 0.15, 0.65, 0.20)
        self.assertGreater(float(isc_high), float(isc_low))

        _, vm_low = compute_percolated_mcc_velocity(0.20, 0.20, 60.0)
        _, vm_high = compute_percolated_mcc_velocity(0.65, 0.20, 60.0)
        self.assertGreater(float(vm_high), float(vm_low))

    def test_complex_step_derivative_accuracy(self) -> None:
        cs_deriv, fd_deriv, rel_err = complex_step_wnt_mcc_derivative(0.25, 0.80, 0.10, 0.55)
        self.assertGreater(cs_deriv, 0.0)
        self.assertGreater(fd_deriv, 0.0)
        self.assertLessEqual(rel_err, 1e-6)

    def test_falsification_negative_controls(self) -> None:
        self.assertTrue(evaluate_falsification_negative_control(poison_negative_notch=True)["passed"])
        self.assertTrue(evaluate_falsification_negative_control(poison_negative_wnt=True)["passed"])
        self.assertTrue(evaluate_falsification_negative_control(poison_invalid_edit_fraction=True)["passed"])

    def test_provenance_and_epistemic_hygiene_metadata(self) -> None:
        self.assertGreaterEqual(len(PARAMETER_PROVENANCE), 5)
        for meta in PARAMETER_PROVENANCE.values():
            self.assertIn("value", meta)
            self.assertIn("tag", meta)
            self.assertIn("not stated", str(meta["source"]).lower())
            self.assertIn("bounds", meta)
        self.assertIn("deductive_column_a", TWO_COLUMN_EPISTEMIC_HYGIENE)
        self.assertIn("popperian_column_b", TWO_COLUMN_EPISTEMIC_HYGIENE)
        self.assertGreaterEqual(len(FIVE_LAYER_LADDER), 5)
        self.assertIn("DOI: 10.1016/j.stem.2011.04.003", PRIMARY_CITATIONS)
        self.assertIn("DOI: 10.1242/dev.177428", PRIMARY_CITATIONS)
        self.assertIn("DOI: 10.1038/s41567-019-0770-2", PRIMARY_CITATIONS)


if __name__ == "__main__":
    unittest.main()
