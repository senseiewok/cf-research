#!/usr/bin/env python3
"""Unit tests for Round 187: Engineered Nonsense Suppressor tRNA eRF1 Competition Bridge (T-0868)."""
from __future__ import annotations

import unittest

from tools.discovery.round187_suppressor_trna_readthrough_bridge import (
    FIVE_LAYER_LADDER,
    PARAMETER_PROVENANCE,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    complex_step_asl_derivative,
    compute_asl_height_um,
    compute_erf1_competition_readthrough_fraction,
    compute_mcc_velocity_mm_min,
    compute_nmd_escape_fraction,
    compute_restored_cftr_protein_fraction,
    compute_transepithelial_chloride_current_pct_wt,
    evaluate_falsification_negative_control,
    evaluate_round187_suppressor_trna_readthrough,
)


class Round187SuppressorTRNAReadthroughTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round187_suppressor_trna_readthrough(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Invariants failed for {key}")
            self.assertGreaterEqual(res.readthrough_fraction, 0.0)
            self.assertLessEqual(res.readthrough_fraction, 1.0)
            self.assertGreaterEqual(res.rescued_mrna_fraction, 0.0)
            self.assertLessEqual(res.rescued_mrna_fraction, 1.0)
            self.assertGreaterEqual(res.restored_protein_fraction, 0.0)
            self.assertLessEqual(res.restored_protein_fraction, 1.0)
            self.assertGreaterEqual(res.transepithelial_current_pct_wt, 0.0)
            self.assertLessEqual(res.transepithelial_current_pct_wt, 100.0)
            self.assertGreaterEqual(res.airway_surface_liquid_height_um, 3.5)
            self.assertLessEqual(res.airway_surface_liquid_height_um, 8.5)
            self.assertGreaterEqual(res.mucociliary_transport_velocity_mm_min, 0.0)
            self.assertLessEqual(res.complex_step_relative_error, 1e-6)

    def test_high_dose_g542x_crosses_kahn_kalai(self) -> None:
        res = evaluate_round187_suppressor_trna_readthrough("ace_trna_therapeutic_high_dose_g542x")
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res.readthrough_fraction, 0.15)
        self.assertGreaterEqual(res.restored_protein_fraction, 0.08)
        self.assertGreaterEqual(res.transepithelial_current_pct_wt, 7.0)
        self.assertGreaterEqual(res.airway_surface_liquid_height_um, 7.0)
        self.assertGreaterEqual(res.mucociliary_transport_velocity_mm_min, 1.0)

    def test_moderate_dose_w1282x_crosses_kahn_kalai(self) -> None:
        res = evaluate_round187_suppressor_trna_readthrough("ace_trna_therapeutic_moderate_dose_w1282x")
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res.readthrough_fraction, 0.50)
        self.assertGreaterEqual(res.restored_protein_fraction, 0.20)
        self.assertGreaterEqual(res.airway_surface_liquid_height_um, 7.0)
        self.assertGreaterEqual(res.mucociliary_transport_velocity_mm_min, 2.0)

    def test_sub_curative_presets_fail_kahn_kalai(self) -> None:
        # Sub-therapeutic low dose
        low_dose = evaluate_round187_suppressor_trna_readthrough("ace_trna_subtherapeutic_low_dose")
        self.assertFalse(low_dose.crosses_kahn_kalai_threshold)
        self.assertLess(low_dose.readthrough_fraction, 0.15)
        self.assertLess(low_dose.restored_protein_fraction, 0.05)
        self.assertEqual(low_dose.mucociliary_transport_velocity_mm_min, 0.0)

        # Untreated CF nonsense baseline
        untreated = evaluate_round187_suppressor_trna_readthrough("untreated_cf_nonsense_baseline")
        self.assertFalse(untreated.crosses_kahn_kalai_threshold)
        self.assertLess(untreated.readthrough_fraction, 0.02)
        self.assertLess(untreated.restored_protein_fraction, 0.01)
        self.assertAlmostEqual(untreated.airway_surface_liquid_height_um, 3.65, delta=0.2)
        self.assertEqual(untreated.mucociliary_transport_velocity_mm_min, 0.0)

        # eRF1 hyperactive termination barrier
        erf1_barrier = evaluate_round187_suppressor_trna_readthrough("erf1_hyperactive_termination_barrier")
        self.assertFalse(erf1_barrier.crosses_kahn_kalai_threshold)
        self.assertLess(erf1_barrier.readthrough_fraction, 0.15)

    def test_complex_step_derivative_accuracy(self) -> None:
        cs_deriv, fd_deriv, rel_err = complex_step_asl_derivative(protein_fraction=0.25, gamma_activity=1.00)
        self.assertGreater(cs_deriv, 0.0)
        self.assertGreater(fd_deriv, 0.0)
        self.assertLessEqual(rel_err, 1e-6)

    def test_poisoned_negative_controls_fail_deterministically(self) -> None:
        res_conc = evaluate_falsification_negative_control(poison_negative_conc=True)
        self.assertTrue(res_conc["passed"])
        self.assertIn("Correctly caught", res_conc["error"])

        res_erf1 = evaluate_falsification_negative_control(poison_negative_erf1=True)
        self.assertTrue(res_erf1["passed"])
        self.assertIn("Correctly caught", res_erf1["error"])

        res_gamma = evaluate_falsification_negative_control(poison_negative_gamma=True)
        self.assertTrue(res_gamma["passed"])
        self.assertIn("Correctly caught", res_gamma["error"])

        clean = evaluate_falsification_negative_control()
        self.assertTrue(clean["passed"])
        self.assertIsNone(clean["error"])

    def test_monotonic_trna_dose_response(self) -> None:
        eta_low = compute_erf1_competition_readthrough_fraction(0.1)
        eta_med = compute_erf1_competition_readthrough_fraction(1.0)
        eta_high = compute_erf1_competition_readthrough_fraction(5.0)
        self.assertLess(eta_low, eta_med)
        self.assertLess(eta_med, eta_high)

    def test_provenance_and_two_column_epistemic_hygiene_metadata(self) -> None:
        self.assertGreaterEqual(len(PRIMARY_CITATIONS), 4)
        self.assertGreaterEqual(len(PARAMETER_PROVENANCE), 5)
        self.assertIn("deductive_theorem", TWO_COLUMN_EPISTEMIC_HYGIENE)
        self.assertIn("biophysical_hypothesis", TWO_COLUMN_EPISTEMIC_HYGIENE)
        self.assertIn("statement", TWO_COLUMN_EPISTEMIC_HYGIENE["biophysical_hypothesis"])
        self.assertIn("falsification_condition", TWO_COLUMN_EPISTEMIC_HYGIENE["biophysical_hypothesis"])
        for layer in ["cool_human_fact", "everyday_analogy", "engineering_lens", "formal_theorem", "falsifiable_code_and_toy"]:
            self.assertIn(layer, FIVE_LAYER_LADDER)


if __name__ == "__main__":
    unittest.main()
