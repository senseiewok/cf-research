#!/usr/bin/env python3
"""Unit tests for Round 162: Monod-Wyman-Changeux (MWC) Allosteric Cooperativity Invariant Bridge (T-0827)."""
from __future__ import annotations

import math
import unittest

from tools.discovery.round162_mwc_allosteric_cooperativity_bridge import (
    FIVE_LAYER_LADDER,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    complex_step_atp_sensitivity,
    compute_channel_open_probability,
    compute_effective_allosteric_constant,
    compute_nbd_dimerization_free_energy_kbt,
    compute_potentiator_coupling_factor,
    evaluate_falsification_negative_control,
    evaluate_round162_mwc_allosteric_cooperativity,
)


class Round162MWCAllostericCooperativityTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round162_mwc_allosteric_cooperativity(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Failed for {key}")
            self.assertGreaterEqual(res.open_probability, 0.0)
            self.assertLessEqual(res.open_probability, 1.0)
            self.assertLessEqual(res.complex_step_relative_error, 1e-10)

    def test_wild_type_cftr_physiological_atp(self) -> None:
        res = evaluate_round162_mwc_allosteric_cooperativity("wild_type_cftr_physiological_atp")
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res.open_probability, 0.35)
        self.assertLessEqual(res.open_probability, 0.50)
        self.assertAlmostEqual(res.potentiator_coupling_factor, 1.0, places=2)

    def test_g551d_severe_gating_arrest_unpotentiated(self) -> None:
        res = evaluate_round162_mwc_allosteric_cooperativity("g551d_severe_gating_arrest_unpotentiated")
        self.assertFalse(res.crosses_kahn_kalai_threshold)
        self.assertLess(res.open_probability, 0.005)
        self.assertGreater(res.dimerization_free_energy_kbt, 5.0)

    def test_g551d_ivacaftor_synergistic_rescue(self) -> None:
        res = evaluate_round162_mwc_allosteric_cooperativity("g551d_ivacaftor_synergistic_rescue")
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res.open_probability, 0.12)
        self.assertGreaterEqual(res.potentiator_coupling_factor, 30.0)

    def test_kahn_kalai_gating_curative_boundary(self) -> None:
        res = evaluate_round162_mwc_allosteric_cooperativity("kahn_kalai_gating_curative_boundary")
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res.open_probability, 0.050)

    def test_g178r_intracellular_loop_potentiated(self) -> None:
        res = evaluate_round162_mwc_allosteric_cooperativity("g178r_intracellular_loop_potentiated")
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res.open_probability, 0.25)

    def test_potentiator_coupling_factor_monotonicity(self) -> None:
        g_0 = compute_potentiator_coupling_factor(0.0)
        g_low = compute_potentiator_coupling_factor(0.1)
        g_high = compute_potentiator_coupling_factor(1.0)
        self.assertEqual(g_0, 1.0)
        self.assertGreater(g_low, g_0)
        self.assertGreater(g_high, g_low)

    def test_atp_concentration_response_monotonicity(self) -> None:
        p_low = compute_channel_open_probability(0.05, 0.0, 24.0, 0.23)
        p_med = compute_channel_open_probability(0.50, 0.0, 24.0, 0.23)
        p_high = compute_channel_open_probability(2.00, 0.0, 24.0, 0.23)
        self.assertGreater(p_med, p_low)
        self.assertGreater(p_high, p_med)

    def test_complex_step_atp_sensitivity(self) -> None:
        deriv_cs, analytical, rel_err = complex_step_atp_sensitivity(
            atp_concentration_mm=1.0,
            potentiator_concentration_um=0.0,
            unliganded_l0=24.0,
            allosteric_ratio_c=0.23,
        )
        self.assertLessEqual(rel_err, 1e-10)
        self.assertGreater(deriv_cs, 0.0)

    def test_negative_control_falsification(self) -> None:
        clean = evaluate_falsification_negative_control()
        self.assertTrue(clean["passed"])
        self.assertEqual(clean["rule"], "AH-10")

        poisoned_atp = evaluate_falsification_negative_control(poison_negative_atp=True)
        self.assertFalse(poisoned_atp["passed"])

        poisoned_pot = evaluate_falsification_negative_control(poison_negative_potentiator=True)
        self.assertFalse(poisoned_pot["passed"])

        poisoned_l0 = evaluate_falsification_negative_control(poison_negative_l0=True)
        self.assertFalse(poisoned_l0["passed"])

        poisoned_ratio = evaluate_falsification_negative_control(poison_negative_allosteric_ratio=True)
        self.assertFalse(poisoned_ratio["passed"])

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
