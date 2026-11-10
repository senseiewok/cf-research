#!/usr/bin/env python3
"""Unit tests for Round 181: Inhaled Small-Molecule ENaC Blocker (ETD001) ASL Bridge (T-0858)."""
from __future__ import annotations

import unittest

from tools.discovery.round181_etd001_enac_blocker_asl_bridge import (
    FIVE_LAYER_LADDER,
    PARAMETER_PROVENANCE,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    complex_step_enac_derivative,
    compute_asl_height_and_mcc_velocity,
    compute_enac_inhibition,
    compute_mucosal_concentration,
    compute_residual_sodium_current,
    evaluate_falsification_negative_control,
    evaluate_round181_etd001_enac_blocker,
)


class Round181Etd001EnacBlockerAslTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round181_etd001_enac_blocker(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Invariants failed for {key}")
            self.assertGreaterEqual(res.enac_inhibition_fraction, 0.0)
            self.assertLessEqual(res.enac_inhibition_fraction, 1.0)
            self.assertGreaterEqual(res.residual_na_current_ua_cm2, 0.0)
            self.assertGreaterEqual(res.asl_height_um, 0.0)
            self.assertLessEqual(res.asl_height_um, 15.0)
            self.assertGreaterEqual(res.mcc_velocity_mm_min, 0.0)
            self.assertLessEqual(res.complex_step_relative_error, 1e-6)

    def test_therapeutic_monotherapy_crosses_kahn_kalai(self) -> None:
        res = evaluate_round181_etd001_enac_blocker("etd001_therapeutic_monotherapy")
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res.asl_height_um, 7.0)
        self.assertGreaterEqual(res.mcc_velocity_mm_min, 1.0)
        self.assertGreaterEqual(res.enac_inhibition_fraction, 0.60)
        self.assertGreaterEqual(res.mucosal_concentration_um, 0.010)

    def test_peak_post_inhalation(self) -> None:
        res = evaluate_round181_etd001_enac_blocker("etd001_peak_post_inhalation")
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res.asl_height_um, 7.0)
        self.assertGreaterEqual(res.mcc_velocity_mm_min, 1.0)
        self.assertGreaterEqual(res.enac_inhibition_fraction, 0.90)

    def test_legacy_amiloride_and_untreated_fail_kahn_kalai(self) -> None:
        # Legacy amiloride rapidly washes out
        amil = evaluate_round181_etd001_enac_blocker("legacy_amiloride_rapid_washout")
        self.assertFalse(amil.crosses_kahn_kalai_threshold)
        self.assertLess(amil.asl_height_um, 7.0)
        self.assertLess(amil.mcc_velocity_mm_min, 0.5)

        # Untreated CF baseline
        untreated = evaluate_round181_etd001_enac_blocker("untreated_cf_airway_baseline")
        self.assertFalse(untreated.crosses_kahn_kalai_threshold)
        self.assertLess(untreated.asl_height_um, 5.0)
        self.assertEqual(untreated.mcc_velocity_mm_min, 0.0)

    def test_subtherapeutic_dose(self) -> None:
        low = evaluate_round181_etd001_enac_blocker("subtherapeutic_etd001_low_dose")
        self.assertFalse(low.crosses_kahn_kalai_threshold)
        self.assertLess(low.asl_height_um, 7.0)
        self.assertLess(low.mcc_velocity_mm_min, 1.0)

    def test_complex_step_derivative(self) -> None:
        cs_deriv, fd_deriv, rel_err = complex_step_enac_derivative(c_muc_um=0.010)
        self.assertLessEqual(rel_err, 1e-6)
        self.assertGreater(cs_deriv, 0.0)

    def test_negative_control_falsification(self) -> None:
        neg_c0 = evaluate_falsification_negative_control(poison_negative_c0=True)
        self.assertTrue(neg_c0["passed"])
        self.assertIn("Correctly caught", neg_c0["error"])

        neg_ic50 = evaluate_falsification_negative_control(poison_negative_ic50=True)
        self.assertTrue(neg_ic50["passed"])
        self.assertIn("Correctly caught", neg_ic50["error"])

        neg_thalf = evaluate_falsification_negative_control(poison_negative_t_half=True)
        self.assertTrue(neg_thalf["passed"])
        self.assertIn("Correctly caught", neg_thalf["error"])

        neg_finh = evaluate_falsification_negative_control(poison_unphysical_f_inh=True)
        self.assertTrue(neg_finh["passed"])
        self.assertIn("Correctly caught", neg_finh["error"])

        baseline = evaluate_falsification_negative_control()
        self.assertTrue(baseline["passed"])
        self.assertIsNone(baseline["error"])

    def test_provenance_and_citations(self) -> None:
        self.assertGreaterEqual(len(PRIMARY_CITATIONS), 3)
        self.assertIn("DOI: 10.1016/j.jcf.2024.06.002", PRIMARY_CITATIONS)
        self.assertIn("DOI: 10.1113/JP286559", PRIMARY_CITATIONS)
        self.assertIn("DOI: 10.1152/physrev.00004.2021", PRIMARY_CITATIONS)

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
