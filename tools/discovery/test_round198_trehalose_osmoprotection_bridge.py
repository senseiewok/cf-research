#!/usr/bin/env python3
"""Unit tests for Round 198: Inhaled Trehalose Osmoprotection & Pseudomonas Biofilm Desiccation Tolerance Invariant Bridge (T-0885)."""
from __future__ import annotations

import unittest

from tools.discovery.round198_trehalose_osmoprotection_bridge import (
    FIVE_LAYER_LADDER,
    PARAMETER_PROVENANCE,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    complex_step_viscosity_derivative,
    compute_asl_height,
    compute_biofilm_survival_fraction,
    compute_mucociliary_velocity,
    compute_sputum_complex_viscosity,
    compute_sputum_yield_stress,
    evaluate_falsification_negative_control,
    evaluate_round198_trehalose_osmoprotection,
)


class Round198TrehaloseOsmoprotectionBridgeTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round198_trehalose_osmoprotection(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Invariants failed for {key}")
            self.assertGreaterEqual(res.treh_inh_mm, 0.0)
            self.assertGreaterEqual(res.f_tobra, 0.0)
            self.assertLessEqual(res.f_tobra, 1.0)
            self.assertGreaterEqual(res.theta_hyd, 0.0)
            self.assertLessEqual(res.theta_hyd, 1.0)
            self.assertGreaterEqual(res.h_asl_um, 2.0)
            self.assertGreaterEqual(res.complex_viscosity_pa_s, 0.10)
            self.assertGreaterEqual(res.yield_stress_pa, 0.0)
            self.assertGreaterEqual(res.biofilm_survival_fraction, 0.0)
            self.assertLessEqual(res.biofilm_survival_fraction, 1.0)
            self.assertGreaterEqual(res.mucociliary_velocity_mm_min, 0.0)
            self.assertLessEqual(res.complex_step_relative_error, 1e-3)

    def test_curative_presets_cross_kahn_kalai(self) -> None:
        res_synergy = evaluate_round198_trehalose_osmoprotection("trehalose_plus_tobramycin_synergy_curative")
        self.assertTrue(res_synergy.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res_synergy.h_asl_um, 4.8)
        self.assertLessEqual(res_synergy.complex_viscosity_pa_s, 3.0)
        self.assertLessEqual(res_synergy.yield_stress_pa, 0.65)
        self.assertLessEqual(res_synergy.biofilm_survival_fraction, 0.20)
        self.assertGreaterEqual(res_synergy.mucociliary_velocity_mm_min, 1.0)

        res_mono = evaluate_round198_trehalose_osmoprotection("high_dose_trehalose_monotherapy_curative")
        self.assertTrue(res_mono.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res_mono.h_asl_um, 5.0)
        self.assertLessEqual(res_mono.complex_viscosity_pa_s, 2.0)
        self.assertLessEqual(res_mono.biofilm_survival_fraction, 0.20)
        self.assertGreaterEqual(res_mono.mucociliary_velocity_mm_min, 2.0)

    def test_sub_curative_presets_fail_kahn_kalai(self) -> None:
        # Untreated desiccated CF airway
        res_untreated = evaluate_round198_trehalose_osmoprotection("untreated_desiccated_cf_airway")
        self.assertFalse(res_untreated.crosses_kahn_kalai_threshold)
        self.assertEqual(res_untreated.h_asl_um, 2.0)
        self.assertEqual(res_untreated.complex_viscosity_pa_s, 100.0)
        self.assertEqual(res_untreated.biofilm_survival_fraction, 1.0)
        self.assertEqual(res_untreated.mucociliary_velocity_mm_min, 0.0)

        # Sub-therapeutic trehalose dose
        res_sub = evaluate_round198_trehalose_osmoprotection("sub_therapeutic_trehalose_dose")
        self.assertFalse(res_sub.crosses_kahn_kalai_threshold)
        self.assertLess(res_sub.h_asl_um, 4.8)
        self.assertGreater(res_sub.complex_viscosity_pa_s, 3.0)
        self.assertEqual(res_sub.mucociliary_velocity_mm_min, 0.0)

        # Tobramycin monotherapy without trehalose
        res_tobra_only = evaluate_round198_trehalose_osmoprotection("tobramycin_monotherapy_without_trehalose")
        self.assertFalse(res_tobra_only.crosses_kahn_kalai_threshold)
        self.assertEqual(res_tobra_only.h_asl_um, 2.0)
        self.assertGreater(res_tobra_only.complex_viscosity_pa_s, 3.0)
        self.assertEqual(res_tobra_only.mucociliary_velocity_mm_min, 0.0)

    def test_complex_step_derivative_accuracy(self) -> None:
        cs_deriv, fd_deriv, rel_err = complex_step_viscosity_derivative(treh_inh_mm=80.0)
        self.assertLess(cs_deriv, 0.0, "Viscosity derivative with respect to trehalose must be negative")
        self.assertLess(fd_deriv, 0.0)
        self.assertLessEqual(rel_err, 1e-3)

    def test_falsification_negative_controls(self) -> None:
        self.assertTrue(evaluate_falsification_negative_control(poison_negative_trehalose=True)["passed"])
        self.assertTrue(evaluate_falsification_negative_control(poison_negative_tobramycin=True)["passed"])
        self.assertTrue(evaluate_falsification_negative_control(poison_negative_viscosity=True)["passed"])

    def test_monotonicity_invariants(self) -> None:
        # Increasing trehalose increases ASL height
        h1 = float(compute_asl_height(20.0))
        h2 = float(compute_asl_height(80.0))
        self.assertGreater(h2, h1)

        # Increasing trehalose decreases sputum viscosity and yield stress
        v1 = float(compute_sputum_complex_viscosity(20.0))
        v2 = float(compute_sputum_complex_viscosity(80.0))
        self.assertGreater(v1, v2)

        tau1 = compute_sputum_yield_stress(20.0)
        tau2 = compute_sputum_yield_stress(80.0)
        self.assertGreater(tau1, tau2)

        # Increasing trehalose decreases biofilm survival fraction
        b1 = compute_biofilm_survival_fraction(20.0, f_tobra=0.5)
        b2 = compute_biofilm_survival_fraction(80.0, f_tobra=0.5)
        self.assertGreater(b1, b2)

    def test_provenance_and_epistemic_hygiene_metadata(self) -> None:
        self.assertGreaterEqual(len(PARAMETER_PROVENANCE), 5)
        for param, meta in PARAMETER_PROVENANCE.items():
            self.assertIn("value", meta)
            self.assertIn("tag", meta)
            self.assertIn(meta["tag"], ("empirical", "derived", "phenomenological_heuristic"))
            self.assertIn("source", meta)
            self.assertIn("bounds", meta)

        self.assertIn("deductive_column_a", TWO_COLUMN_EPISTEMIC_HYGIENE)
        self.assertIn("popperian_column_b", TWO_COLUMN_EPISTEMIC_HYGIENE)
        self.assertGreaterEqual(len(FIVE_LAYER_LADDER), 5)
        self.assertGreaterEqual(len(PRIMARY_CITATIONS), 2)
        self.assertIn("DOI: 10.1371/journal.pgen.1009524", PRIMARY_CITATIONS)
        self.assertIn("DOI: 10.1146/annurev.med.58.071905.105316", PRIMARY_CITATIONS)


if __name__ == "__main__":
    unittest.main()
