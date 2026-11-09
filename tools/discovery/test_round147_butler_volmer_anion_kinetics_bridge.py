#!/usr/bin/env python3
"""Unit tests for Round 147: Butler-Volmer Electrochemical Transfer Coefficient & Overpotential Invariant Bridge (T-0798)."""
from __future__ import annotations

import unittest

from tools.discovery.round147_butler_volmer_anion_kinetics_bridge import (
    FIVE_LAYER_LADDER,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    evaluate_falsification_negative_control,
    evaluate_round147_butler_volmer_anion_kinetics,
)


class Round147ButlerVolmerAnionKineticsTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round147_butler_volmer_anion_kinetics(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Failed for {key}")
            self.assertGreater(res.transfer_coefficient_alpha, 0.0)
            self.assertLess(res.transfer_coefficient_alpha, 1.0)
            self.assertGreater(res.exchange_current_pa, 0.0)
            self.assertGreater(res.chord_conductance_ps, 0.0)
            self.assertGreaterEqual(res.epithelial_short_circuit_current_ua_cm2, 0.0)

    def test_wild_type_cftr_apical_conduction(self) -> None:
        wt = evaluate_round147_butler_volmer_anion_kinetics("wild_type_cftr_apical_conduction")
        self.assertTrue(wt.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(wt.chord_conductance_ps, 8.0)
        self.assertGreaterEqual(wt.epithelial_short_circuit_current_ua_cm2, 30.0)
        self.assertAlmostEqual(wt.fractional_wild_type_current, 1.0, delta=0.05)
        self.assertAlmostEqual(wt.transfer_coefficient_alpha, 0.40, delta=0.01)

    def test_cf_defects_and_gating_collapse(self) -> None:
        g551d = evaluate_round147_butler_volmer_anion_kinetics("cf_g551d_gating_rectification_collapse")
        df508 = evaluate_round147_butler_volmer_anion_kinetics("cf_deltaf508_pore_barrier_stasis")
        for mutant in (g551d, df508):
            self.assertFalse(mutant.crosses_kahn_kalai_threshold)
            self.assertLessEqual(mutant.fractional_wild_type_current, 0.10)
        self.assertLessEqual(g551d.epithelial_short_circuit_current_ua_cm2, 1.0)
        self.assertLessEqual(df508.epithelial_short_circuit_current_ua_cm2, 1.0)

    def test_kahn_kalai_electrochemical_boundary(self) -> None:
        kk = evaluate_round147_butler_volmer_anion_kinetics("kahn_kalai_electrochemical_boundary")
        self.assertTrue(kk.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(kk.chord_conductance_ps, 2.0)
        self.assertGreaterEqual(kk.epithelial_short_circuit_current_ua_cm2, 3.79)
        self.assertGreaterEqual(kk.fractional_wild_type_current, 0.118)

    def test_trikafta_rectification_restoration_rescue(self) -> None:
        rescued = evaluate_round147_butler_volmer_anion_kinetics("trikafta_rectification_restoration_rescue")
        self.assertTrue(rescued.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(rescued.chord_conductance_ps, 8.0)
        self.assertGreaterEqual(rescued.epithelial_short_circuit_current_ua_cm2, 25.0)
        self.assertGreaterEqual(rescued.fractional_wild_type_current, 0.85)

    def test_negative_control_falsification(self) -> None:
        normal = evaluate_falsification_negative_control(poison_unphysical_transfer_coeff=False)
        self.assertTrue(normal["passed"])
        self.assertEqual(normal["rule"], "AH-10")
        poisoned = evaluate_falsification_negative_control(poison_unphysical_transfer_coeff=True)
        self.assertFalse(poisoned["passed"])

    def test_overpotential_sweep_and_zero_bias(self) -> None:
        zero_bias = evaluate_round147_butler_volmer_anion_kinetics(
            "wild_type_cftr_apical_conduction", custom_eta_mv=0.0
        )
        self.assertAlmostEqual(zero_bias.single_channel_current_pa, 0.0, places=6)
        self.assertGreater(zero_bias.chord_conductance_ps, 0.0)

        small_drive = evaluate_round147_butler_volmer_anion_kinetics(
            "wild_type_cftr_apical_conduction", custom_eta_mv=-25.0
        )
        large_drive = evaluate_round147_butler_volmer_anion_kinetics(
            "wild_type_cftr_apical_conduction", custom_eta_mv=-75.0
        )
        self.assertGreater(abs(large_drive.single_channel_current_pa), abs(small_drive.single_channel_current_pa))
        self.assertGreater(
            large_drive.epithelial_short_circuit_current_ua_cm2,
            small_drive.epithelial_short_circuit_current_ua_cm2,
        )

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
