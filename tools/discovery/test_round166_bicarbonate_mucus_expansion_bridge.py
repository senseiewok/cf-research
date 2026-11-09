#!/usr/bin/env python3
"""Unit tests for Round 166: Goldman-Hodgkin-Katz-Donnan Epithelial Bicarbonate/Chloride Antiport & Apical pH Regulation Invariant Bridge (T-0831)."""
from __future__ import annotations

import unittest

from tools.discovery.round166_bicarbonate_mucus_expansion_bridge import (
    FIVE_LAYER_LADDER,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    complex_step_ph_derivative,
    compute_asl_ph,
    compute_free_calcium_concentration,
    compute_ghk_anion_flux,
    compute_mucin_expansion_fold,
    evaluate_falsification_negative_control,
    evaluate_round166_bicarbonate_mucus_expansion,
)


class Round166BicarbonateMucusExpansionTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round166_bicarbonate_mucus_expansion(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Failed for {key}")
            self.assertGreaterEqual(res.asl_ph, 4.0)
            self.assertLessEqual(res.asl_ph, 9.0)
            self.assertGreaterEqual(res.mucin_expansion_fold, 1.0)
            self.assertLessEqual(res.complex_step_relative_error, 1e-10)

    def test_healthy_airway_bicarbonate_expansion(self) -> None:
        res = evaluate_round166_bicarbonate_mucus_expansion("healthy_airway_bicarbonate_expansion")
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res.asl_ph, 7.40)
        self.assertGreaterEqual(res.mucin_expansion_fold, 400.0)
        self.assertLessEqual(res.free_calcium_conc_mm, 0.80)
        self.assertGreaterEqual(res.calcium_chelation_efficiency_pct, 65.0)

    def test_trikafta_rescued_f508del_bicarbonate(self) -> None:
        res = evaluate_round166_bicarbonate_mucus_expansion("trikafta_rescued_f508del_bicarbonate")
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res.asl_ph, 7.35)
        self.assertGreaterEqual(res.mucin_expansion_fold, 350.0)

    def test_kahn_kalai_expansion_boundary(self) -> None:
        res = evaluate_round166_bicarbonate_mucus_expansion("kahn_kalai_expansion_boundary")
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res.mucin_expansion_fold, 250.0)
        self.assertGreaterEqual(res.asl_ph, 7.15)

    def test_untreated_cf_acidified_mucus_stasis(self) -> None:
        res = evaluate_round166_bicarbonate_mucus_expansion("untreated_cf_acidified_mucus_stasis")
        self.assertFalse(res.crosses_kahn_kalai_threshold)
        self.assertLess(res.mucin_expansion_fold, 160.0)
        self.assertLess(res.asl_ph, 6.90)
        self.assertGreater(res.free_calcium_conc_mm, 1.40)

    def test_severe_respiratory_acidosis_plugging(self) -> None:
        res = evaluate_round166_bicarbonate_mucus_expansion("severe_respiratory_acidosis_plugging")
        self.assertFalse(res.crosses_kahn_kalai_threshold)
        self.assertLess(res.mucin_expansion_fold, 50.0)
        self.assertLess(res.asl_ph, 6.55)

    def test_ph_monotonicity_and_buffering(self) -> None:
        ph_low = compute_asl_ph(bicarbonate_conc_mm=5.0, pco2_mmhg=40.0)
        ph_mid = compute_asl_ph(bicarbonate_conc_mm=25.0, pco2_mmhg=40.0)
        ph_high = compute_asl_ph(bicarbonate_conc_mm=50.0, pco2_mmhg=40.0)
        self.assertGreater(ph_mid, ph_low)
        self.assertGreater(ph_high, ph_mid)

    def test_calcium_chelation_monotonicity(self) -> None:
        ca_free_low_bicarb = compute_free_calcium_concentration(2.5, 5.0)
        ca_free_high_bicarb = compute_free_calcium_concentration(2.5, 30.0)
        self.assertLess(ca_free_high_bicarb, ca_free_low_bicarb)

    def test_complex_step_ph_sensitivity(self) -> None:
        cs_deriv, analytical, rel_err = complex_step_ph_derivative(
            bicarbonate_conc_mm=25.0,
            pco2_mmhg=40.0,
        )
        self.assertLessEqual(rel_err, 1e-10)
        self.assertGreater(cs_deriv, 0.0)

    def test_ghk_flux_permeability_dependence(self) -> None:
        j_wt = compute_ghk_anion_flux(permeability_cm_s=1e-6, conc_in_mm=15.0, conc_out_mm=25.0, transmembrane_voltage_mv=-40.0)
        j_mut = compute_ghk_anion_flux(permeability_cm_s=2e-7, conc_in_mm=15.0, conc_out_mm=25.0, transmembrane_voltage_mv=-40.0)
        self.assertGreater(abs(j_wt), abs(j_mut))

    def test_negative_control_falsification(self) -> None:
        neg1 = evaluate_falsification_negative_control(poison_negative_bicarbonate=True)
        self.assertFalse(neg1["passed"])
        self.assertIn("poisoned", neg1["error"])

        neg2 = evaluate_falsification_negative_control(poison_negative_pco2=True)
        self.assertFalse(neg2["passed"])

        neg3 = evaluate_falsification_negative_control(poison_negative_calcium=True)
        self.assertFalse(neg3["passed"])

        neg4 = evaluate_falsification_negative_control(poison_negative_permeability=True)
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
