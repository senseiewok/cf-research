#!/usr/bin/env python3
"""Unit tests for Round 122: Doob h-Transform Conditional Markov Survival Across 92 Class I Nonsense Alleles."""
from __future__ import annotations

import unittest

from tools.discovery.round122_doob_h_transform_92_nonsense_alleles_mrna_survival_bridge import (
    CANONICAL_92_NONSENSE_ALLELES,
    FIVE_LAYER_LADDER,
    PRIMARY_CITATIONS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    evaluate_all_92_nonsense_alleles,
    evaluate_doob_nonsense_survival,
)


class Round122DoobHTransform92NonsenseTest(unittest.TestCase):
    def test_all_92_canonical_nonsense_alleles_indexed_and_monotone(self) -> None:
        self.assertEqual(len(CANONICAL_92_NONSENSE_ALLELES), 92)
        rows = evaluate_all_92_nonsense_alleles(nmd_shield=0.0)
        self.assertEqual(len(rows), 92)
        # Early PTC Q2X (codon 2) has much lower basal Doob survival h than terminal Q1412X (codon 1412)
        q2x = rows[0]
        q1412x = [r for r in rows if r["allele"] == "Q1412X"][0]
        self.assertEqual(q2x["allele"], "Q2X")
        self.assertLess(float(q2x["doob_h"]), 0.12)
        self.assertGreater(float(q1412x["doob_h"]), 0.70)

    def test_doob_stochastic_row_conservation(self) -> None:
        for key in (
            "q1412x_terminal_ejc_escape",
            "w1282x_exon23_aso_suptrna",
            "g542x_exon12_suptrna_shield",
            "q2x_exon1_early_ptc_rescued",
            "g542x_unshielded_nmd_sink",
        ):
            res = evaluate_doob_nonsense_survival(key)
            self.assertAlmostEqual(res.stochastic_conservation_sum, 1.0, places=6)
            self.assertTrue(res.falsification_invariant_satisfied)

    def test_unshielded_g542x_fails_kahn_kalai_threshold_while_shielded_crosses(self) -> None:
        unshielded = evaluate_doob_nonsense_survival("g542x_unshielded_nmd_sink")
        shielded = evaluate_doob_nonsense_survival("g542x_exon12_suptrna_shield")
        self.assertFalse(unshielded.crosses_kahn_kalai_11_8pct)
        self.assertLess(unshielded.rescued_functional_cftr_pct, 2.0)
        self.assertTrue(shielded.crosses_kahn_kalai_11_8pct)
        self.assertGreaterEqual(shielded.rescued_functional_cftr_pct, 11.8)

    def test_terminal_exon27_evades_ejc_decay_barrier(self) -> None:
        term = evaluate_doob_nonsense_survival("q1412x_terminal_ejc_escape")
        self.assertEqual(term.exon_number, 27)
        self.assertEqual(term.downstream_ejc_count, 0.0)
        self.assertGreater(term.doob_survival_harmonic_h, 0.75)
        self.assertTrue(term.crosses_kahn_kalai_11_8pct)

    def test_codon_position_monotonicity_of_doob_survival_harmonic(self) -> None:
        early = evaluate_doob_nonsense_survival("w1282x_exon23_aso_suptrna", custom_codon=50)
        mid = evaluate_doob_nonsense_survival("w1282x_exon23_aso_suptrna", custom_codon=700)
        late = evaluate_doob_nonsense_survival("w1282x_exon23_aso_suptrna", custom_codon=1412)
        self.assertLess(early.doob_survival_harmonic_h, mid.doob_survival_harmonic_h)
        self.assertLess(mid.doob_survival_harmonic_h, late.doob_survival_harmonic_h)

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
