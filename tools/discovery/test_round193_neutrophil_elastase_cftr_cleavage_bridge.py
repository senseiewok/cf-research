#!/usr/bin/env python3
"""Unit tests for Round 193: Neutrophil Elastase CFTR Cleavage & A1AT Shielding Bridge (T-0876)."""
from __future__ import annotations

import unittest

from tools.discovery.round193_neutrophil_elastase_cftr_cleavage_bridge import (
    FIVE_LAYER_LADDER,
    PARAMETER_PROVENANCE,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    complex_step_asl_derivative,
    compute_asl_height_um,
    compute_free_elastase_concentration,
    compute_intact_cftr_fraction,
    compute_mcc_velocity_mm_min,
    compute_transepithelial_current_pct_wt,
    evaluate_falsification_negative_control,
    evaluate_round193_neutrophil_elastase_cftr_cleavage,
)


class Round193NeutrophilElastaseCleavageTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round193_neutrophil_elastase_cftr_cleavage(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Invariants failed for {key}")
            self.assertGreaterEqual(res.hne_tot_ug_ml, 0.0)
            self.assertGreaterEqual(res.a1at_tot_ug_ml, 0.0)
            self.assertGreaterEqual(res.hne_free_ug_ml, 0.0)
            self.assertGreaterEqual(res.intact_cftr_fraction, 0.0)
            self.assertLessEqual(res.intact_cftr_fraction, 1.0)
            self.assertGreaterEqual(res.transepithelial_current_pct_wt, 0.0)
            self.assertLessEqual(res.transepithelial_current_pct_wt, 100.0)
            self.assertGreaterEqual(res.airway_surface_liquid_height_um, 3.5)
            self.assertLessEqual(res.airway_surface_liquid_height_um, 8.5)
            self.assertGreaterEqual(res.mucociliary_transport_velocity_mm_min, 0.0)
            self.assertLessEqual(res.complex_step_relative_error, 1e-6)

    def test_curative_presets_cross_kahn_kalai(self) -> None:
        res1 = evaluate_round193_neutrophil_elastase_cftr_cleavage("trikafta_plus_high_dose_a1at_shielding")
        self.assertTrue(res1.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res1.intact_cftr_fraction, 0.50)
        self.assertGreaterEqual(res1.transepithelial_current_pct_wt, 12.0)
        self.assertGreaterEqual(res1.airway_surface_liquid_height_um, 7.0)
        self.assertGreaterEqual(res1.mucociliary_transport_velocity_mm_min, 1.0)
        self.assertLessEqual(res1.hne_free_ug_ml, 2.0)

        res2 = evaluate_round193_neutrophil_elastase_cftr_cleavage("trikafta_plus_stoichiometric_a1at_shielding")
        self.assertTrue(res2.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res2.intact_cftr_fraction, 0.50)
        self.assertGreaterEqual(res2.transepithelial_current_pct_wt, 12.0)
        self.assertGreaterEqual(res2.airway_surface_liquid_height_um, 7.0)
        self.assertLessEqual(res2.hne_free_ug_ml, 2.0)

    def test_sub_curative_presets_fail_kahn_kalai(self) -> None:
        # Untreated severe cleavage
        res_severe = evaluate_round193_neutrophil_elastase_cftr_cleavage("untreated_severe_elastase_cftr_cleavage")
        self.assertFalse(res_severe.crosses_kahn_kalai_threshold)
        self.assertLess(res_severe.intact_cftr_fraction, 0.10)
        self.assertLess(res_severe.transepithelial_current_pct_wt, 1.0)
        self.assertGreater(res_severe.hne_free_ug_ml, 10.0)

        # Modulator alone in inflamed niche
        res_mod_alone = evaluate_round193_neutrophil_elastase_cftr_cleavage("modulator_alone_without_antiprotease_in_inflamed_niche")
        self.assertFalse(res_mod_alone.crosses_kahn_kalai_threshold)
        self.assertLess(res_mod_alone.intact_cftr_fraction, 0.50)
        self.assertLess(res_mod_alone.transepithelial_current_pct_wt, 12.0)

        # Sub-inhibitory A1AT dose impairment
        res_sub = evaluate_round193_neutrophil_elastase_cftr_cleavage("sub_inhibitory_a1at_dose_impairment")
        self.assertFalse(res_sub.crosses_kahn_kalai_threshold)
        self.assertLess(res_sub.intact_cftr_fraction, 0.50)
        self.assertLess(res_sub.transepithelial_current_pct_wt, 12.0)
        self.assertGreater(res_sub.hne_free_ug_ml, 2.0)

    def test_complex_step_derivative_accuracy(self) -> None:
        cs_deriv, fd_deriv, rel_err = complex_step_asl_derivative(current_pct_wt=25.0)
        self.assertGreater(cs_deriv, 0.0)
        self.assertGreater(fd_deriv, 0.0)
        self.assertLessEqual(rel_err, 1e-6)

    def test_falsification_negative_controls(self) -> None:
        self.assertTrue(evaluate_falsification_negative_control(poison_negative_hne=True)["passed"])
        self.assertTrue(evaluate_falsification_negative_control(poison_negative_a1at=True)["passed"])
        self.assertTrue(evaluate_falsification_negative_control(poison_negative_po=True)["passed"])

    def test_monotonicity_invariants(self) -> None:
        # Higher A1AT reduces free elastase and increases intact CFTR
        free1 = compute_free_elastase_concentration(hne_tot_ug_ml=25.0, a1at_tot_ug_ml=10.0)
        free2 = compute_free_elastase_concentration(hne_tot_ug_ml=25.0, a1at_tot_ug_ml=20.0)
        self.assertGreater(free1, free2)

        intact1 = compute_intact_cftr_fraction(hne_free_ug_ml=free1)
        intact2 = compute_intact_cftr_fraction(hne_free_ug_ml=free2)
        self.assertLess(intact1, intact2)

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
