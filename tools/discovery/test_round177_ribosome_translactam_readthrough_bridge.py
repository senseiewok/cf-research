#!/usr/bin/env python3
"""Unit tests for Round 177: Ribosome Translactam & eRF3a Degradation Invariant Bridge (T-0851)."""
from __future__ import annotations

import unittest

from tools.discovery.round177_ribosome_translactam_readthrough_bridge import (
    FIVE_LAYER_LADDER,
    PARAMETER_PROVENANCE,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    complex_step_acc_derivative,
    compute_erf3a_termination_rate,
    compute_mrna_abundance_fraction,
    compute_net_cftr_synthesis,
    compute_readthrough_efficiency,
    compute_translactam_accommodation_rate,
    evaluate_falsification_negative_control,
    evaluate_round177_ribosome_translactam_readthrough,
)


class Round177RibosomeTranslactamReadthroughTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round177_ribosome_translactam_readthrough(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Invariants failed for {key}")
            self.assertGreaterEqual(res.readthrough_efficiency, 0.0)
            self.assertLessEqual(res.readthrough_efficiency, 1.0)
            self.assertGreaterEqual(res.mrna_abundance_fraction, 0.0)
            self.assertLessEqual(res.mrna_abundance_fraction, 1.0)
            self.assertGreaterEqual(res.net_cftr_synthesis_fraction, 0.0)
            self.assertLessEqual(res.net_cftr_synthesis_fraction, 1.0)
            self.assertGreaterEqual(res.restored_isc_pct_wt, 0.0)
            self.assertLessEqual(res.erf3a_degradation_fraction, 0.85)
            self.assertLessEqual(res.complex_step_relative_error, 1e-6)

    def test_triple_combination_translactam_cc90009_smg1i(self) -> None:
        res = evaluate_round177_ribosome_translactam_readthrough(
            "triple_combination_translactam_cc90009_smg1i"
        )
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res.net_cftr_synthesis_fraction, 0.10)
        self.assertGreaterEqual(res.restored_isc_pct_wt, 10.0)
        self.assertGreaterEqual(res.readthrough_efficiency, 0.20)
        self.assertGreaterEqual(res.mrna_abundance_fraction, 0.50)
        self.assertGreaterEqual(res.erf3a_degradation_fraction, 0.65)

    def test_dual_combinations_and_monotherapy(self) -> None:
        # Without eRF3a degrader
        dual_no_deg = evaluate_round177_ribosome_translactam_readthrough(
            "dual_combination_translactam_smg1i"
        )
        self.assertFalse(dual_no_deg.crosses_kahn_kalai_threshold)
        self.assertLess(dual_no_deg.net_cftr_synthesis_fraction, 0.10)

        # Without NMD inhibitor
        dual_no_nmd = evaluate_round177_ribosome_translactam_readthrough(
            "dual_combination_translactam_cc90009"
        )
        self.assertFalse(dual_no_nmd.crosses_kahn_kalai_threshold)
        self.assertLess(dual_no_nmd.net_cftr_synthesis_fraction, 0.10)

        # Translactam monotherapy
        mono = evaluate_round177_ribosome_translactam_readthrough("monotherapy_translactam")
        self.assertFalse(mono.crosses_kahn_kalai_threshold)
        self.assertLess(mono.net_cftr_synthesis_fraction, 0.05)

    def test_untreated_w1282x_baseline(self) -> None:
        res = evaluate_round177_ribosome_translactam_readthrough("untreated_w1282x_baseline")
        self.assertFalse(res.crosses_kahn_kalai_threshold)
        self.assertLess(res.net_cftr_synthesis_fraction, 0.005)
        self.assertLess(res.readthrough_efficiency, 0.01)

    def test_complex_step_acc_derivative(self) -> None:
        cs_deriv, fd_deriv, rel_err = complex_step_acc_derivative(
            k_acc=1.48,
            k_term=3.50,
        )
        self.assertLessEqual(rel_err, 1e-6)
        self.assertGreater(cs_deriv, 0.0)

    def test_negative_control_falsification(self) -> None:
        neg_trans = evaluate_falsification_negative_control(poison_negative_translactam=True)
        self.assertTrue(neg_trans["passed"])
        self.assertIn("Correctly caught", neg_trans["error"])

        neg_deg = evaluate_falsification_negative_control(poison_negative_deg=True)
        self.assertTrue(neg_deg["passed"])
        self.assertIn("Correctly caught", neg_deg["error"])

        neg_temp = evaluate_falsification_negative_control(poison_negative_temp=True)
        self.assertTrue(neg_temp["passed"])
        self.assertIn("Correctly caught", neg_temp["error"])

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
