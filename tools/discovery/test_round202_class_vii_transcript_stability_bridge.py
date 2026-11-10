#!/usr/bin/env python3
"""Unit tests for Round 202: Class VII Large Deletion vs Nonsense Transcript Stability Bridge (T-0893)."""
from __future__ import annotations

import unittest

from tools.discovery.round202_class_vii_transcript_stability_bridge import (
    FIVE_LAYER_LADDER,
    PARAMETER_PROVENANCE,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    complex_step_alpha_rt_derivative,
    compute_cftr_mrna_abundance,
    compute_surface_cftr_and_current,
    evaluate_falsification_negative_control,
    evaluate_round202_class_vii_transcript_stability,
)


class Round202ClassVIITranscriptStabilityBridgeTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round202_class_vii_transcript_stability(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Invariants failed for {key}")
            self.assertGreaterEqual(res.f_class7, 0.0)
            self.assertLessEqual(res.f_class7, 1.0)
            self.assertGreaterEqual(res.chi_ptc, 0.0)
            self.assertLessEqual(res.chi_ptc, 1.0)
            self.assertGreaterEqual(res.f_mrna_relative, 0.0)
            self.assertLessEqual(res.f_mrna_relative, 1.0)
            self.assertGreaterEqual(res.p_cftr_surface, 0.0)
            self.assertLessEqual(res.p_cftr_surface, 1.0)
            self.assertGreaterEqual(res.i_sc_equiv_pct_wt, 0.0)
            self.assertLessEqual(res.i_sc_equiv_pct_wt, 100.0)
            self.assertLessEqual(res.complex_step_relative_error, 1e-6)

    def test_curative_presets_cross_kahn_kalai(self) -> None:
        res_cur = evaluate_round202_class_vii_transcript_stability("ptc_nmd_inhibition_readthrough_curative")
        self.assertTrue(res_cur.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res_cur.f_mrna_relative, 0.40)
        self.assertGreaterEqual(res_cur.p_cftr_surface, 0.12)
        self.assertGreaterEqual(res_cur.i_sc_equiv_pct_wt, 25.0)

        res_het = evaluate_round202_class_vii_transcript_stability("compound_het_class7_ptc_rescued_curative")
        self.assertTrue(res_het.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res_het.f_mrna_relative, 0.40)
        self.assertGreaterEqual(res_het.i_sc_equiv_pct_wt, 25.0)

    def test_pathology_and_negative_controls_fail_kahn_kalai(self) -> None:
        res_null = evaluate_round202_class_vii_transcript_stability("homozygous_class7_large_deletion_null_cf")
        self.assertFalse(res_null.crosses_kahn_kalai_threshold)
        self.assertEqual(res_null.f_mrna_relative, 0.0)
        self.assertEqual(res_null.p_cftr_surface, 0.0)
        self.assertEqual(res_null.i_sc_equiv_pct_wt, 0.0)

        res_untreated = evaluate_round202_class_vii_transcript_stability("untreated_homozygous_ptc_nmd_decay_cf")
        self.assertFalse(res_untreated.crosses_kahn_kalai_threshold)
        self.assertLess(res_untreated.f_mrna_relative, 0.30)
        self.assertLess(res_untreated.i_sc_equiv_pct_wt, 10.0)

        res_no_rt = evaluate_round202_class_vii_transcript_stability("nmd_inhibited_readthrough_deficient_arrest")
        self.assertFalse(res_no_rt.crosses_kahn_kalai_threshold)
        self.assertLess(res_no_rt.i_sc_equiv_pct_wt, 25.0)

    def test_monotonicity_invariants(self) -> None:
        m_wt = compute_cftr_mrna_abundance(0.0, 0.0, 0.0)
        m_class7 = compute_cftr_mrna_abundance(1.0, 0.0, 0.0)
        self.assertEqual(float(m_wt), 1.0)
        self.assertEqual(float(m_class7), 0.0)

        m_low_nmd = compute_cftr_mrna_abundance(0.0, 1.0, 0.10)
        m_high_nmd = compute_cftr_mrna_abundance(0.0, 1.0, 0.85)
        self.assertGreater(float(m_high_nmd), float(m_low_nmd))

        _, i_low_rt = compute_surface_cftr_and_current(0.5, 1.0, 0.01)
        _, i_high_rt = compute_surface_cftr_and_current(0.5, 1.0, 0.15)
        self.assertGreater(float(i_high_rt), float(i_low_rt))

    def test_complex_step_derivative_accuracy(self) -> None:
        cs_deriv, fd_deriv, rel_err = complex_step_alpha_rt_derivative(0.10, 0.90, 0.70, 0.08)
        self.assertGreater(cs_deriv, 0.0)
        self.assertGreater(fd_deriv, 0.0)
        self.assertLessEqual(rel_err, 1e-6)

    def test_falsification_negative_controls(self) -> None:
        self.assertTrue(evaluate_falsification_negative_control(poison_negative_class7=True)["passed"])
        self.assertTrue(evaluate_falsification_negative_control(poison_exceeded_class7=True)["passed"])
        self.assertTrue(evaluate_falsification_negative_control(poison_negative_ptc=True)["passed"])
        self.assertTrue(evaluate_falsification_negative_control(poison_negative_readthrough=True)["passed"])

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
        self.assertEqual(len(PRIMARY_CITATIONS), 2)


if __name__ == "__main__":
    unittest.main()
