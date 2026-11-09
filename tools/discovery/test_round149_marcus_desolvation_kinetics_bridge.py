#!/usr/bin/env python3
"""Unit tests for Round 149: Marcus Reorganization Free Energy & Outer-Sphere Activation Invariant Bridge (T-0802)."""
from __future__ import annotations

import unittest

from tools.discovery.round149_marcus_desolvation_kinetics_bridge import (
    FIVE_LAYER_LADDER,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    evaluate_falsification_negative_control,
    evaluate_round149_marcus_desolvation_kinetics,
)


class Round149MarcusDesolvationKineticsTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round149_marcus_desolvation_kinetics(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Failed for {key}")
            self.assertGreater(res.reorganization_energy_kj_mol, 0.0)
            self.assertGreaterEqual(res.activation_free_energy_kj_mol, 0.0)
            self.assertGreater(res.chord_conductance_ps, 0.0)
            self.assertGreaterEqual(res.epithelial_short_circuit_current_ua_cm2, 0.0)

    def test_wild_type_cftr_apical_conduction(self) -> None:
        wt = evaluate_round149_marcus_desolvation_kinetics("wild_type_cftr_apical_conduction")
        self.assertTrue(wt.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(wt.chord_conductance_ps, 8.0)
        self.assertGreaterEqual(wt.epithelial_short_circuit_current_ua_cm2, 30.0)
        self.assertAlmostEqual(wt.fractional_wild_type_current, 1.0, delta=0.05)
        self.assertAlmostEqual(wt.reorganization_energy_kj_mol, 42.0, delta=0.1)
        self.assertLessEqual(wt.activation_free_energy_kj_mol, 7.0)

    def test_cf_defects_and_gating_collapse(self) -> None:
        g551d = evaluate_round149_marcus_desolvation_kinetics("cf_g551d_desolvation_collapse")
        df508 = evaluate_round149_marcus_desolvation_kinetics("cf_deltaf508_misfolded_vestibule")
        for mutant in (g551d, df508):
            self.assertFalse(mutant.crosses_kahn_kalai_threshold)
            self.assertLessEqual(mutant.fractional_wild_type_current, 0.05)
            self.assertLessEqual(mutant.epithelial_short_circuit_current_ua_cm2, 0.5)
            self.assertGreaterEqual(mutant.activation_free_energy_kj_mol, 15.0)

    def test_kahn_kalai_desolvation_boundary(self) -> None:
        kk = evaluate_round149_marcus_desolvation_kinetics("kahn_kalai_desolvation_boundary")
        self.assertTrue(kk.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(kk.chord_conductance_ps, 2.0)
        self.assertGreaterEqual(kk.epithelial_short_circuit_current_ua_cm2, 3.75)
        self.assertGreaterEqual(kk.fractional_wild_type_current, 0.118)

    def test_trikafta_allosteric_desolvation_rescue(self) -> None:
        rescued = evaluate_round149_marcus_desolvation_kinetics("trikafta_allosteric_desolvation_rescue")
        self.assertTrue(rescued.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(rescued.chord_conductance_ps, 7.0)
        self.assertGreaterEqual(rescued.epithelial_short_circuit_current_ua_cm2, 25.0)
        self.assertGreaterEqual(rescued.fractional_wild_type_current, 0.85)

    def test_negative_control_falsification(self) -> None:
        normal = evaluate_falsification_negative_control(poison_unphysical_reorganization_energy=False)
        self.assertTrue(normal["passed"])
        self.assertEqual(normal["rule"], "AH-10")
        poisoned = evaluate_falsification_negative_control(poison_unphysical_reorganization_energy=True)
        self.assertFalse(poisoned["passed"])

    def test_free_energy_sweep_and_barrierless_limit(self) -> None:
        barrierless = evaluate_round149_marcus_desolvation_kinetics(
            "wild_type_cftr_apical_conduction", custom_lambda_kj=40.0, custom_dg0_kj=-40.0
        )
        self.assertAlmostEqual(barrierless.activation_free_energy_kj_mol, 0.0, places=6)

        symmetric = evaluate_round149_marcus_desolvation_kinetics(
            "wild_type_cftr_apical_conduction", custom_lambda_kj=40.0, custom_dg0_kj=0.0
        )
        self.assertAlmostEqual(symmetric.activation_free_energy_kj_mol, 10.0, places=4)

        favorable = evaluate_round149_marcus_desolvation_kinetics(
            "wild_type_cftr_apical_conduction", custom_dg0_kj=-10.0
        )
        unfavorable = evaluate_round149_marcus_desolvation_kinetics(
            "wild_type_cftr_apical_conduction", custom_dg0_kj=+5.0
        )
        self.assertLess(favorable.activation_free_energy_kj_mol, unfavorable.activation_free_energy_kj_mol)
        self.assertGreater(favorable.chord_conductance_ps, unfavorable.chord_conductance_ps)

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
