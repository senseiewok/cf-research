#!/usr/bin/env python3
"""Unit tests for Round 179: High-Activity Adenine Base Editor (ABE8e) Splice Rescue Bridge (T-0853)."""
from __future__ import annotations

import unittest

from tools.discovery.round179_abe8e_splice_base_editing_bridge import (
    FIVE_LAYER_LADDER,
    PARAMETER_PROVENANCE,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    complex_step_kcat_derivative,
    compute_deamination_rate,
    compute_positional_selectivity,
    compute_restored_splicing_efficiency,
    compute_target_editing_efficiency,
    evaluate_falsification_negative_control,
    evaluate_round179_abe8e_splice_base_editing,
)


class Round179Abe8eSpliceBaseEditingTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round179_abe8e_splice_base_editing(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Invariants failed for {key}")
            self.assertGreaterEqual(res.target_editing_efficiency, 0.0)
            self.assertLessEqual(res.target_editing_efficiency, 1.0)
            self.assertGreaterEqual(res.restored_splicing_efficiency, 0.0)
            self.assertLessEqual(res.restored_splicing_efficiency, 1.0)
            self.assertGreaterEqual(res.restored_cftr_synthesis_fraction, 0.0)
            self.assertLessEqual(res.restored_cftr_synthesis_fraction, 1.0)
            self.assertGreaterEqual(res.restored_isc_pct_wt, 0.0)
            self.assertGreaterEqual(res.bystander_selectivity_ratio, 1.0)
            self.assertLessEqual(res.complex_step_relative_error, 1e-6)

    def test_optimized_abe8e_lnp_mrna_sgrna(self) -> None:
        res = evaluate_round179_abe8e_splice_base_editing("optimized_abe8e_lnp_mrna_sgrna")
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res.target_editing_efficiency, 0.30)
        self.assertGreaterEqual(res.restored_splicing_efficiency, 0.40)
        self.assertGreaterEqual(res.restored_isc_pct_wt, 15.0)
        self.assertGreaterEqual(res.bystander_selectivity_ratio, 8.0)

    def test_standard_abe8e_lnp_delivery(self) -> None:
        res = evaluate_round179_abe8e_splice_base_editing("standard_abe8e_lnp_delivery")
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res.target_editing_efficiency, 0.30)
        self.assertGreaterEqual(res.restored_isc_pct_wt, 15.0)

    def test_subcurative_and_legacy_regimes(self) -> None:
        # Low dose
        low = evaluate_round179_abe8e_splice_base_editing("low_dose_abe8e_subcurative")
        self.assertFalse(low.crosses_kahn_kalai_threshold)
        self.assertLess(low.target_editing_efficiency, 0.30)

        # Legacy ABE7.10
        legacy = evaluate_round179_abe8e_splice_base_editing("legacy_abe710_unoptimized_editor")
        self.assertFalse(legacy.crosses_kahn_kalai_threshold)
        self.assertLess(legacy.restored_isc_pct_wt, 5.0)
        self.assertLess(legacy.target_editing_efficiency, 0.10)

        # Untreated baseline
        untreated = evaluate_round179_abe8e_splice_base_editing("untreated_1717_1ga_baseline")
        self.assertFalse(untreated.crosses_kahn_kalai_threshold)
        self.assertEqual(untreated.target_editing_efficiency, 0.0)
        self.assertLess(untreated.restored_isc_pct_wt, 1.0)

    def test_complex_step_kcat_derivative(self) -> None:
        cs_deriv, fd_deriv, rel_err = complex_step_kcat_derivative(
            c_abe_nm=25.0,
            k_cat_s=0.008,
            t_exp_s=300.0,
        )
        self.assertLessEqual(rel_err, 1e-6)
        self.assertGreater(cs_deriv, 0.0)

    def test_negative_control_falsification(self) -> None:
        neg_c = evaluate_falsification_negative_control(poison_negative_c_abe=True)
        self.assertTrue(neg_c["passed"])
        self.assertIn("Correctly caught", neg_c["error"])

        neg_k = evaluate_falsification_negative_control(poison_negative_k_cat=True)
        self.assertTrue(neg_k["passed"])
        self.assertIn("Correctly caught", neg_k["error"])

        neg_t = evaluate_falsification_negative_control(poison_negative_temp=True)
        self.assertTrue(neg_t["passed"])
        self.assertIn("Correctly caught", neg_t["error"])

        neg_del = evaluate_falsification_negative_control(poison_oversized_delivery=True)
        self.assertTrue(neg_del["passed"])
        self.assertIn("Correctly caught", neg_del["error"])

        baseline = evaluate_falsification_negative_control()
        self.assertTrue(baseline["passed"])
        self.assertIsNone(baseline["error"])

    def test_provenance_and_citations(self) -> None:
        self.assertGreaterEqual(len(PRIMARY_CITATIONS), 3)
        self.assertGreaterEqual(len(PARAMETER_PROVENANCE), 5)
        self.assertIn("deductive_theorem", TWO_COLUMN_EPISTEMIC_HYGIENE)
        self.assertIn("biophysical_hypothesis", TWO_COLUMN_EPISTEMIC_HYGIENE)
        self.assertEqual(len(FIVE_LAYER_LADDER), 5)


if __name__ == "__main__":
    unittest.main()
