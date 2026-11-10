#!/usr/bin/env python3
"""Unit tests for Round 196: CRISPR-Cas12a LbCas12a Airway Epithelial Gene Correction Fidelity & NHEJ Suppressor Bridge (T-0881)."""
from __future__ import annotations

import unittest

from tools.discovery.round196_lbcas12a_hdr_gene_correction_bridge import (
    FIVE_LAYER_LADDER,
    PARAMETER_PROVENANCE,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    complex_step_template_derivative,
    compute_chloride_conductance_restoration,
    compute_cleavage_fractions,
    compute_hdr_branching_probability,
    evaluate_falsification_negative_control,
    evaluate_round196_lbcas12a_hdr_gene_correction,
)


class Round196Cas12aHDRGeneCorrectionTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round196_lbcas12a_hdr_gene_correction(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Invariants failed for {key}")
            self.assertGreaterEqual(res.rnp_nm, 0.0)
            self.assertGreaterEqual(res.s_fidelity, 0.0)
            self.assertLessEqual(res.s_fidelity, 1.0)
            self.assertGreaterEqual(res.template_hdr_nm, 0.0)
            self.assertGreaterEqual(res.inhibitor_nhej_um, 0.0)
            self.assertGreaterEqual(res.cleavage_fraction, 0.0)
            self.assertLessEqual(res.cleavage_fraction, 1.0)
            self.assertGreaterEqual(res.hdr_branching_probability, 0.0)
            self.assertLessEqual(res.hdr_branching_probability, 1.0)
            self.assertGreaterEqual(res.allele_correction_efficiency, 0.0)
            self.assertLessEqual(res.allele_correction_efficiency, 1.0)
            self.assertGreaterEqual(res.restored_chloride_current_pct_wt, 0.0)
            self.assertLessEqual(res.restored_chloride_current_pct_wt, 100.0)
            self.assertLessEqual(res.complex_step_relative_error, 1e-4)

    def test_curative_presets_cross_kahn_kalai(self) -> None:
        res1 = evaluate_round196_lbcas12a_hdr_gene_correction("enlbcas12a_plus_dna_pk_inhibitor_m3814_hdr")
        self.assertTrue(res1.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res1.allele_correction_efficiency, 0.15)
        self.assertGreaterEqual(res1.s_fidelity, 0.98)
        self.assertGreaterEqual(res1.restored_chloride_current_pct_wt, 25.0)
        self.assertLessEqual(res1.off_target_fraction, 0.02)

        res2 = evaluate_round196_lbcas12a_hdr_gene_correction("enlbcas12a_ssdna_donor_without_inhibitor")
        self.assertTrue(res2.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res2.allele_correction_efficiency, 0.15)
        self.assertGreaterEqual(res2.s_fidelity, 0.98)
        self.assertGreaterEqual(res2.restored_chloride_current_pct_wt, 25.0)

    def test_sub_curative_presets_fail_kahn_kalai(self) -> None:
        # Untreated control
        res_untreated = evaluate_round196_lbcas12a_hdr_gene_correction("untreated_nonsense_cftr_control")
        self.assertFalse(res_untreated.crosses_kahn_kalai_threshold)
        self.assertEqual(res_untreated.cleavage_fraction, 0.0)
        self.assertEqual(res_untreated.allele_correction_efficiency, 0.0)
        self.assertEqual(res_untreated.restored_chloride_current_pct_wt, 0.0)

        # Wild-type Cas12a without inhibitor (lower fidelity, higher off-target)
        res_wt = evaluate_round196_lbcas12a_hdr_gene_correction("wildtype_cas12a_without_nhej_suppression")
        self.assertFalse(res_wt.crosses_kahn_kalai_threshold)
        self.assertLess(res_wt.s_fidelity, 0.98)
        self.assertGreater(res_wt.off_target_fraction, 0.02)

        # Toxic high dose inhibitor
        res_tox = evaluate_round196_lbcas12a_hdr_gene_correction("toxic_high_dose_nhej_inhibitor_arrest")
        self.assertFalse(res_tox.crosses_kahn_kalai_threshold)
        self.assertGreater(res_tox.cytotoxicity, 0.50)
        self.assertLess(res_tox.restored_chloride_current_pct_wt, 25.0)

    def test_complex_step_derivative_accuracy(self) -> None:
        cs_deriv, fd_deriv, rel_err = complex_step_template_derivative(
            template_hdr_nm=150.0,
            inhibitor_nhej_um=1.0,
            f_cleave=0.70,
        )
        self.assertGreater(cs_deriv, 0.0, "Current derivative with respect to template must be positive")
        self.assertGreater(fd_deriv, 0.0)
        self.assertLessEqual(rel_err, 1e-4)

    def test_falsification_negative_controls(self) -> None:
        self.assertTrue(evaluate_falsification_negative_control(poison_negative_rnp=True)["passed"])
        self.assertTrue(evaluate_falsification_negative_control(poison_negative_template=True)["passed"])
        self.assertTrue(evaluate_falsification_negative_control(poison_negative_inhibitor=True)["passed"])

    def test_monotonicity_invariants(self) -> None:
        # Increasing NHEJ inhibitor concentration increases HDR branching probability
        p1 = compute_hdr_branching_probability(template_hdr_nm=100.0, inhibitor_nhej_um=0.0)
        p2 = compute_hdr_branching_probability(template_hdr_nm=100.0, inhibitor_nhej_um=1.0)
        self.assertGreater(p2, p1)

        # Increasing template concentration increases HDR branching probability
        p3 = compute_hdr_branching_probability(template_hdr_nm=50.0, inhibitor_nhej_um=0.5)
        p4 = compute_hdr_branching_probability(template_hdr_nm=200.0, inhibitor_nhej_um=0.5)
        self.assertGreater(p4, p3)

    def test_provenance_and_epistemic_hygiene_metadata(self) -> None:
        self.assertGreaterEqual(len(PARAMETER_PROVENANCE), 5)
        for p_name, p_meta in PARAMETER_PROVENANCE.items():
            self.assertIn(p_meta["tag"], ("empirical", "derived", "phenomenological_heuristic"))
            self.assertIn("source", p_meta)
            self.assertIn("bounds", p_meta)

        self.assertIn("deductive_theorem", TWO_COLUMN_EPISTEMIC_HYGIENE)
        self.assertIn("biophysical_hypothesis", TWO_COLUMN_EPISTEMIC_HYGIENE)
        self.assertIn("falsification_condition", TWO_COLUMN_EPISTEMIC_HYGIENE["biophysical_hypothesis"])

        self.assertEqual(len(FIVE_LAYER_LADDER), 5)
        for layer in ("cool_human_fact", "everyday_analogy", "engineering_lens", "formal_theorem", "falsifiable_code_and_toy"):
            self.assertIn(layer, FIVE_LAYER_LADDER)

        self.assertGreaterEqual(len(PRIMARY_CITATIONS), 2)


if __name__ == "__main__":
    unittest.main()
