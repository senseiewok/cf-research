#!/usr/bin/env python3
"""Unit tests for Round 182: Airway Basal Stem Cell Clonal Persistence Bridge (T-0859)."""
from __future__ import annotations

import unittest

from tools.discovery.round182_basal_stem_cell_clonal_persistence_bridge import (
    FIVE_LAYER_LADDER,
    PARAMETER_PROVENANCE,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    complex_step_basal_derivative,
    compute_asl_and_mcc_from_ciliated_fraction,
    compute_basal_clonal_persistence,
    compute_ciliated_cftr_reconstitution,
    evaluate_falsification_negative_control,
    evaluate_round182_basal_stem_persistence,
)


class Round182BasalStemCellClonalPersistenceTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round182_basal_stem_persistence(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Invariants failed for {key}")
            self.assertGreaterEqual(res.persistent_basal_fraction, 0.0)
            self.assertLessEqual(res.persistent_basal_fraction, 1.0)
            self.assertGreaterEqual(res.ciliated_cftr_fraction, 0.0)
            self.assertLessEqual(res.ciliated_cftr_fraction, 1.0)
            self.assertGreaterEqual(res.equivalent_isc_pct_wt, 0.0)
            self.assertGreaterEqual(res.asl_height_um, 0.0)
            self.assertLessEqual(res.asl_height_um, 15.0)
            self.assertGreaterEqual(res.mcc_velocity_mm_min, 0.0)
            self.assertLessEqual(res.complex_step_relative_error, 1e-6)

    def test_basal_targeted_durable_editing_crosses_kahn_kalai(self) -> None:
        res = evaluate_round182_basal_stem_persistence("basal_targeted_durable_editing")
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res.ciliated_cftr_fraction, 0.20)
        self.assertGreaterEqual(res.asl_height_um, 7.0)
        self.assertGreaterEqual(res.mcc_velocity_mm_min, 1.0)

    def test_high_efficiency_basal_engraftment(self) -> None:
        res = evaluate_round182_basal_stem_persistence("high_efficiency_basal_engraftment")
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res.ciliated_cftr_fraction, 0.50)
        self.assertGreaterEqual(res.asl_height_um, 7.5)
        self.assertGreaterEqual(res.mcc_velocity_mm_min, 2.0)

    def test_transient_differentiated_washout_fails_at_12_months(self) -> None:
        washout = evaluate_round182_basal_stem_persistence("transient_differentiated_only_washout")
        self.assertFalse(washout.crosses_kahn_kalai_threshold)
        self.assertLess(washout.ciliated_cftr_fraction, 0.02)
        self.assertLess(washout.asl_height_um, 4.0)
        self.assertEqual(washout.mcc_velocity_mm_min, 0.0)

    def test_subthreshold_and_untreated_fail_kahn_kalai(self) -> None:
        sub = evaluate_round182_basal_stem_persistence("sub_threshold_basal_editing")
        self.assertFalse(sub.crosses_kahn_kalai_threshold)
        self.assertLess(sub.ciliated_cftr_fraction, 0.10)
        self.assertLess(sub.asl_height_um, 7.0)

        untreated = evaluate_round182_basal_stem_persistence("untreated_cf_epithelium_baseline")
        self.assertFalse(untreated.crosses_kahn_kalai_threshold)
        self.assertEqual(untreated.ciliated_cftr_fraction, 0.0)
        self.assertEqual(untreated.mcc_velocity_mm_min, 0.0)

    def test_complex_step_derivative(self) -> None:
        cs_deriv, fd_deriv, rel_err = complex_step_basal_derivative(f_basal_0=0.35)
        self.assertLessEqual(rel_err, 1e-6)
        self.assertGreater(cs_deriv, 0.8)

    def test_negative_control_falsification(self) -> None:
        neg_fb = evaluate_falsification_negative_control(poison_negative_basal_fraction=True)
        self.assertTrue(neg_fb["passed"])
        self.assertIn("Correctly caught", neg_fb["error"])

        neg_prob = evaluate_falsification_negative_control(poison_unbalanced_probabilities=True)
        self.assertTrue(neg_prob["passed"])
        self.assertIn("Correctly caught", neg_prob["error"])

        neg_hl = evaluate_falsification_negative_control(poison_negative_half_life=True)
        self.assertTrue(neg_hl["passed"])
        self.assertIn("Correctly caught", neg_hl["error"])

        neg_div = evaluate_falsification_negative_control(poison_negative_division_rate=True)
        self.assertTrue(neg_div["passed"])
        self.assertIn("Correctly caught", neg_div["error"])

        baseline = evaluate_falsification_negative_control()
        self.assertTrue(baseline["passed"])
        self.assertIsNone(baseline["error"])

    def test_provenance_and_citations(self) -> None:
        self.assertGreaterEqual(len(PRIMARY_CITATIONS), 3)
        self.assertIn("DOI: 10.1172/jci.insight.206375", PRIMARY_CITATIONS)
        self.assertIn("DOI: 10.1016/j.ymthe.2018.04.014", PRIMARY_CITATIONS)
        self.assertIn("DOI: 10.3389/fbioe.2026.1855907", PRIMARY_CITATIONS)

        self.assertGreaterEqual(len(PARAMETER_PROVENANCE), 5)
        for key, meta in PARAMETER_PROVENANCE.items():
            self.assertIn("tag", meta)
            self.assertIn("source", meta)
            self.assertIn("bounds", meta)

        self.assertIn("deductive_theorem", TWO_COLUMN_EPISTEMIC_HYGIENE)
        self.assertIn("biophysical_hypothesis", TWO_COLUMN_EPISTEMIC_HYGIENE)
        self.assertIn("cool_human_fact", FIVE_LAYER_LADDER)
        self.assertIn("falsifiable_code_and_toy", FIVE_LAYER_LADDER)


if __name__ == "__main__":
    unittest.main()
