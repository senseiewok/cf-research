#!/usr/bin/env python3
"""Unit tests for Round 157: Spangler-Goldman Electrodiffusion & Nernst-Planck Invariant Bridge (T-0818)."""
from __future__ import annotations

import math
import unittest

from tools.discovery.round157_goldman_electrodiffusion_bridge import (
    FIVE_LAYER_LADDER,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    complex_step_reversal_potential_sensitivity,
    compute_asl_ph,
    compute_ghk_current_density,
    compute_ghk_reversal_potential_mv,
    evaluate_falsification_negative_control,
    evaluate_round157_goldman_electrodiffusion,
)


class Round157ElectrodiffusionTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round157_goldman_electrodiffusion(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Failed for {key}")
            self.assertGreater(res.asl_ph, 6.0)
            self.assertGreater(res.rectification_ratio, 0.0)
            self.assertLessEqual(res.complex_step_relative_error, 1e-11)

    def test_wild_type_physiological_secretion(self) -> None:
        wt = evaluate_round157_goldman_electrodiffusion("wild_type_physiological_secretion")
        self.assertTrue(wt.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(wt.total_anion_current_density_ua_cm2, 28.0)
        self.assertGreaterEqual(wt.asl_ph, 7.30)
        self.assertAlmostEqual(wt.open_probability, 0.42, delta=0.01)

    def test_cf_anion_depletion_and_acidosis(self) -> None:
        df508 = evaluate_round157_goldman_electrodiffusion("cf_deltaf508_anion_depletion")
        g551d = evaluate_round157_goldman_electrodiffusion("cf_g551d_gating_blockade")

        for mutant in (df508, g551d):
            self.assertFalse(mutant.crosses_kahn_kalai_threshold)
            self.assertLess(mutant.total_anion_current_density_ua_cm2, 1.0)
            self.assertLess(mutant.asl_ph, 6.85)

    def test_kahn_kalai_rehydration_boundary(self) -> None:
        kk = evaluate_round157_goldman_electrodiffusion("kahn_kalai_rehydration_boundary")
        self.assertTrue(kk.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(kk.total_anion_current_density_ua_cm2, 3.75)
        self.assertGreaterEqual(kk.asl_ph, 7.00)

    def test_trikafta_dual_permeability_rescue(self) -> None:
        rescued = evaluate_round157_goldman_electrodiffusion("trikafta_dual_permeability_rescue")
        self.assertTrue(rescued.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(rescued.total_anion_current_density_ua_cm2, 25.0)
        self.assertGreaterEqual(rescued.asl_ph, 7.25)
        self.assertGreaterEqual(rescued.open_probability, 0.35)

    def test_symmetrical_zero_potential_invariant(self) -> None:
        # In symmetrical concentrations, V_rev must be identically 0.0 mV
        v_rev_sym = compute_ghk_reversal_potential_mv(
            c_cl_int_mm=100.0,
            c_cl_ext_mm=100.0,
            c_hco3_int_mm=25.0,
            c_hco3_ext_mm=25.0,
            p_ratio_hco3_cl=0.22,
        )
        self.assertAlmostEqual(v_rev_sym, 0.0, places=9)

    def test_eisenman_selectivity_sequence(self) -> None:
        p_cl = 1.0
        p_hco3 = 0.22
        p_scn = 3.80
        # Eisenman Sequence I: P_SCN > P_Cl > P_HCO3
        self.assertGreater(p_scn, p_cl)
        self.assertGreater(p_cl, p_hco3)

    def test_complex_step_sensitivity(self) -> None:
        deriv_cs, deriv_ana, rel_err = complex_step_reversal_potential_sensitivity(
            c_cl_int_mm=40.0,
            c_cl_ext_mm=115.0,
            c_hco3_int_mm=16.0,
            c_hco3_ext_mm=24.0,
        )
        self.assertGreater(deriv_cs, 0.0)
        self.assertLessEqual(rel_err, 1e-11)

    def test_negative_control_falsification(self) -> None:
        clean = evaluate_falsification_negative_control()
        self.assertTrue(clean["passed"])
        self.assertEqual(clean["rule"], "AH-10")

        poisoned_temp = evaluate_falsification_negative_control(poison_negative_temperature=True)
        self.assertFalse(poisoned_temp["passed"])

        poisoned_perm = evaluate_falsification_negative_control(poison_negative_permeability=True)
        self.assertFalse(poisoned_perm["passed"])

        poisoned_conc = evaluate_falsification_negative_control(poison_negative_concentration=True)
        self.assertFalse(poisoned_conc["passed"])

    def test_two_column_epistemic_hygiene_and_five_layer_ladder(self) -> None:
        self.assertIn("deductive_theorem", TWO_COLUMN_EPISTEMIC_HYGIENE)
        self.assertIn("biophysical_hypothesis", TWO_COLUMN_EPISTEMIC_HYGIENE)
        self.assertIn("falsification_condition", TWO_COLUMN_EPISTEMIC_HYGIENE["biophysical_hypothesis"])
        self.assertIn("cool_human_fact", FIVE_LAYER_LADDER)
        self.assertIn("everyday_analogy", FIVE_LAYER_LADDER)
        self.assertIn("engineering_lens", FIVE_LAYER_LADDER)
        self.assertIn("formal_theorem", FIVE_LAYER_LADDER)
        self.assertIn("falsifiable_code_and_toy", FIVE_LAYER_LADDER)
        self.assertGreaterEqual(len(PRIMARY_CITATIONS), 3)


if __name__ == "__main__":
    unittest.main()
