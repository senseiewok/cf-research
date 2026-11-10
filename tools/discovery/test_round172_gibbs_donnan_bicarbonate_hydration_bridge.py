#!/usr/bin/env python3
"""Unit tests for Round 172: Gibbs-Donnan Bicarbonate-Chloride Antiport & Secretory Hydration Invariant Bridge (T-0837)."""
from __future__ import annotations

import unittest

from tools.discovery.round172_gibbs_donnan_bicarbonate_hydration_bridge import (
    FIVE_LAYER_LADDER,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    complex_step_asl_ph_derivative,
    compute_amp_bactericidal_kill,
    compute_antiport_free_energy,
    compute_asl_ph,
    compute_mucin_mesh_expansion,
    evaluate_falsification_negative_control,
    evaluate_round172_gibbs_donnan_hydration,
)


class Round172GibbsDonnanTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round172_gibbs_donnan_hydration(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Failed for {key}")
            self.assertGreater(res.bicarbonate_asl_mm, 0.0)
            self.assertGreaterEqual(res.asl_ph, 5.0)
            self.assertLessEqual(res.asl_ph, 9.0)
            self.assertGreaterEqual(res.antimicrobial_kill_pct, 0.0)
            self.assertLessEqual(res.antimicrobial_kill_pct, 98.0)
            self.assertGreaterEqual(res.mucin_mesh_expansion_factor, 1.0)
            self.assertLessEqual(res.complex_step_relative_error, 1e-10)

    def test_healthy_wt_homeostasis(self) -> None:
        res = evaluate_round172_gibbs_donnan_hydration("healthy_wt_homeostasis")
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res.asl_ph, 7.40)
        self.assertGreaterEqual(res.antimicrobial_kill_pct, 90.0)
        self.assertGreaterEqual(res.mucin_mesh_expansion_factor, 9.0)

    def test_cf_untreated_acidification(self) -> None:
        res = evaluate_round172_gibbs_donnan_hydration("cf_untreated_acidification")
        self.assertFalse(res.crosses_kahn_kalai_threshold)
        self.assertLess(res.asl_ph, 6.90)
        self.assertLess(res.antimicrobial_kill_pct, 40.0)
        self.assertLess(res.mucin_mesh_expansion_factor, 3.0)

    def test_cf_extreme_acidification_atp12a(self) -> None:
        res = evaluate_round172_gibbs_donnan_hydration("cf_extreme_acidification_atp12a")
        self.assertFalse(res.crosses_kahn_kalai_threshold)
        self.assertLess(res.asl_ph, 6.70)
        self.assertLess(res.antimicrobial_kill_pct, 20.0)
        self.assertLess(res.mucin_mesh_expansion_factor, 2.0)

    def test_cftr_modulator_trikafta_rescue(self) -> None:
        res = evaluate_round172_gibbs_donnan_hydration("cftr_modulator_trikafta_rescue")
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res.asl_ph, 7.35)
        self.assertGreaterEqual(res.antimicrobial_kill_pct, 90.0)
        self.assertGreaterEqual(res.mucin_mesh_expansion_factor, 8.5)

    def test_bicarbonate_inhalation_partial(self) -> None:
        res = evaluate_round172_gibbs_donnan_hydration("bicarbonate_inhalation_partial")
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res.asl_ph, 7.20)
        self.assertGreaterEqual(res.antimicrobial_kill_pct, 80.0)
        self.assertGreaterEqual(res.mucin_mesh_expansion_factor, 6.8)

    def test_complex_step_asl_ph_derivative(self) -> None:
        cs_deriv, analytical, rel_err = complex_step_asl_ph_derivative(
            bicarbonate_asl_mm=24.0,
        )
        self.assertLessEqual(rel_err, 1e-10)
        self.assertGreater(cs_deriv, 0.0)  # ASL pH strictly increases with bicarbonate concentration

    def test_negative_control_falsification(self) -> None:
        neg1 = evaluate_falsification_negative_control(poison_negative_bicarbonate=True)
        self.assertFalse(neg1["passed"])
        self.assertIn("poisoned", neg1["error"])

        neg2 = evaluate_falsification_negative_control(poison_negative_pco2=True)
        self.assertFalse(neg2["passed"])

        neg3 = evaluate_falsification_negative_control(poison_negative_temperature=True)
        self.assertFalse(neg3["passed"])

        neg4 = evaluate_falsification_negative_control(poison_negative_chloride=True)
        self.assertFalse(neg4["passed"])

        clean = evaluate_falsification_negative_control()
        self.assertTrue(clean["passed"])

    def test_epistemic_hygiene_and_ladder_structure(self) -> None:
        self.assertIn("deductive_theorem", TWO_COLUMN_EPISTEMIC_HYGIENE)
        self.assertIn("biophysical_hypothesis", TWO_COLUMN_EPISTEMIC_HYGIENE)
        self.assertGreaterEqual(len(FIVE_LAYER_LADDER), 5)
        self.assertGreaterEqual(len(PRIMARY_CITATIONS), 4)


if __name__ == "__main__":
    unittest.main()
