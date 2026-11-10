#!/usr/bin/env python3
"""Unit tests for Round 199: Basal Stem Cell Notch/Wnt Signaling & Multiciliated Cell Lineage Differentiation Invariant Bridge (T-0888)."""
from __future__ import annotations

import unittest

from tools.discovery.round199_basal_stem_cell_differentiation_bridge import (
    FIVE_LAYER_LADDER,
    PARAMETER_PROVENANCE,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    complex_step_mcc_derivative,
    compute_lineage_coverage_fractions,
    compute_mucociliary_transport_velocity,
    compute_multiciliated_commitment_probability,
    evaluate_falsification_negative_control,
    evaluate_round199_basal_stem_cell_differentiation,
)


class Round199BasalStemCellDifferentiationBridgeTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round199_basal_stem_cell_differentiation(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Invariants failed for {key}")
            self.assertGreaterEqual(res.wnt_activity, 0.0)
            self.assertLessEqual(res.wnt_activity, 1.0)
            self.assertGreaterEqual(res.notch_activity, 0.0)
            self.assertLessEqual(res.notch_activity, 1.0)
            self.assertGreaterEqual(res.p_ciliated, 0.0)
            self.assertLessEqual(res.p_ciliated, 1.0)
            self.assertGreaterEqual(res.f_ciliated, 0.0)
            self.assertLessEqual(res.f_ciliated, 1.0)
            self.assertGreaterEqual(res.f_secretory, 0.0)
            self.assertLessEqual(res.f_secretory, 1.0)
            self.assertGreaterEqual(res.f_basal_reserve, -1e-6)
            self.assertLessEqual(res.f_basal_reserve, 1.0 + 1e-6)
            self.assertAlmostEqual(res.f_ciliated + res.f_secretory + res.f_basal_reserve, 1.0, places=5)
            self.assertGreaterEqual(res.v_mcc_mm_min, 0.0)
            self.assertLessEqual(res.v_mcc_mm_min, 10.0)
            self.assertLessEqual(res.complex_step_relative_error, 1e-3)

    def test_curative_presets_cross_kahn_kalai(self) -> None:
        res_curative = evaluate_round199_basal_stem_cell_differentiation(
            "wnt_high_notch_low_curative_regeneration"
        )
        self.assertTrue(res_curative.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res_curative.p_ciliated, 0.60)
        self.assertGreaterEqual(res_curative.f_ciliated, 0.45)
        self.assertGreaterEqual(res_curative.v_mcc_mm_min, 4.0)
        self.assertLessEqual(res_curative.secretory_to_ciliated_ratio, 0.80)

        res_super = evaluate_round199_basal_stem_cell_differentiation(
            "wnt_max_notch_inhibited_super_curative"
        )
        self.assertTrue(res_super.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res_super.p_ciliated, 0.75)
        self.assertGreaterEqual(res_super.f_ciliated, 0.60)
        self.assertGreaterEqual(res_super.v_mcc_mm_min, 5.0)

    def test_sub_curative_presets_fail_kahn_kalai(self) -> None:
        res_metaplasia = evaluate_round199_basal_stem_cell_differentiation(
            "cf_inflammatory_notch_hyperactive_metaplasia"
        )
        self.assertFalse(res_metaplasia.crosses_kahn_kalai_threshold)
        self.assertLess(res_metaplasia.f_ciliated, 0.20)
        self.assertEqual(res_metaplasia.v_mcc_mm_min, 0.0)
        self.assertGreater(res_metaplasia.secretory_to_ciliated_ratio, 0.80)

        res_partial = evaluate_round199_basal_stem_cell_differentiation(
            "sub_curative_partial_differentiation"
        )
        self.assertFalse(res_partial.crosses_kahn_kalai_threshold)
        self.assertLess(res_partial.f_ciliated, 0.45)

        res_dormant = evaluate_round199_basal_stem_cell_differentiation(
            "wnt_deficient_dormant_basal_state"
        )
        self.assertFalse(res_dormant.crosses_kahn_kalai_threshold)
        self.assertGreater(res_dormant.f_basal_reserve, 0.70)
        self.assertEqual(res_dormant.v_mcc_mm_min, 0.0)

    def test_monotonicity_invariants(self) -> None:
        # Increasing Wnt increases multiciliated commitment
        p1 = float(compute_multiciliated_commitment_probability(0.3, 0.2))
        p2 = float(compute_multiciliated_commitment_probability(0.8, 0.2))
        self.assertGreater(p2, p1)

        # Increasing Notch suppresses multiciliated commitment
        p3 = float(compute_multiciliated_commitment_probability(0.8, 0.6))
        self.assertLess(p3, p2)

        # Increasing ciliated fraction increases clearance velocity
        v1 = float(compute_mucociliary_transport_velocity(0.30))
        v2 = float(compute_mucociliary_transport_velocity(0.60))
        self.assertGreater(v2, v1)

    def test_complex_step_derivative_accuracy(self) -> None:
        cs_deriv, fd_deriv, rel_err = complex_step_mcc_derivative(0.70, 0.15)
        self.assertGreater(cs_deriv, 0.0)
        self.assertLessEqual(rel_err, 1e-3)

    def test_falsification_negative_controls(self) -> None:
        neg_wnt = evaluate_falsification_negative_control(poison_negative_wnt=True)
        self.assertTrue(neg_wnt["passed"])
        self.assertIn("rejected", neg_wnt["reason"].lower())

        neg_notch = evaluate_falsification_negative_control(poison_negative_notch=True)
        self.assertTrue(neg_notch["passed"])
        self.assertIn("rejected", neg_notch["reason"].lower())

        neg_comp = evaluate_falsification_negative_control(poison_negative_competence=True)
        self.assertTrue(neg_comp["passed"])
        self.assertIn("rejected", neg_comp["reason"].lower())

    def test_provenance_and_epistemic_hygiene_metadata(self) -> None:
        self.assertGreaterEqual(len(PRIMARY_CITATIONS), 2)
        for cit in PRIMARY_CITATIONS:
            self.assertTrue(cit.startswith("DOI: 10."))

        self.assertIn("deductive_column_a", TWO_COLUMN_EPISTEMIC_HYGIENE)
        self.assertIn("popperian_column_b", TWO_COLUMN_EPISTEMIC_HYGIENE)
        self.assertIn("formal_theorem", TWO_COLUMN_EPISTEMIC_HYGIENE["deductive_column_a"])
        self.assertIn("falsifiable_hypothesis", TWO_COLUMN_EPISTEMIC_HYGIENE["popperian_column_b"])

        for layer in ("cool_human_fact", "everyday_analogy", "engineering_lens", "formal_theorem", "falsifiable_code_and_toy"):
            self.assertIn(layer, FIVE_LAYER_LADDER)
            self.assertGreater(len(FIVE_LAYER_LADDER[layer]), 20)

        self.assertGreaterEqual(len(PARAMETER_PROVENANCE), 5)
        for k, v in PARAMETER_PROVENANCE.items():
            self.assertIn("value", v)
            self.assertIn("tag", v)
            self.assertIn(v["tag"], ("empirical", "derived", "phenomenological_heuristic"))
            self.assertIn("source", v)
            self.assertIn("bounds", v)


if __name__ == "__main__":
    unittest.main()
