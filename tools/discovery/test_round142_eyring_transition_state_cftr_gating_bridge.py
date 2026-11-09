#!/usr/bin/env python3
"""Unit tests for Round 142: Eyring Transition State Theory & Gating Energy Barrier Activation Invariant Bridge (T-0789)."""
from __future__ import annotations

import unittest

from tools.discovery.round142_eyring_transition_state_cftr_gating_bridge import (
    FIVE_LAYER_LADDER,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    evaluate_falsification_negative_control,
    evaluate_round142_eyring_transition_state_cftr_gating,
)


class Round142EyringTransitionStateCftrGatingTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round142_eyring_transition_state_cftr_gating(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Failed for {key}")
            self.assertGreater(res.complex_step_gating_sensitivity, 0.0)
            self.assertGreater(res.temperature_k, 0.0)
            self.assertGreater(res.opening_rate_k_open_s_1, 0.0)
            self.assertGreater(res.closing_rate_k_close_s_1, 0.0)
            self.assertGreater(res.closed_dwell_time_tau_closed_s, 0.0)
            self.assertGreater(res.open_dwell_time_tau_open_ms, 0.0)
            self.assertGreater(res.open_probability_p_o, 0.0)
            self.assertLess(res.open_probability_p_o, 1.0)

    def test_wild_type_physiological_gating(self) -> None:
        wt = evaluate_round142_eyring_transition_state_cftr_gating("wild_type_physiological_gating")
        self.assertTrue(wt.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(wt.open_probability_p_o, 0.25)
        self.assertGreaterEqual(wt.opening_rate_k_open_s_1, 0.70)
        self.assertGreaterEqual(wt.open_dwell_time_tau_open_ms, 300.0)
        self.assertLessEqual(wt.closed_dwell_time_tau_closed_s, 2.0)

    def test_g551d_severe_gating_defect(self) -> None:
        g551d = evaluate_round142_eyring_transition_state_cftr_gating("g551d_severe_gating_defect")
        self.assertFalse(g551d.crosses_kahn_kalai_threshold)
        self.assertLessEqual(g551d.open_probability_p_o, 0.005)
        self.assertLessEqual(g551d.opening_rate_k_open_s_1, 0.01)
        self.assertGreaterEqual(g551d.closed_dwell_time_tau_closed_s, 100.0)

    def test_kahn_kalai_gating_threshold(self) -> None:
        kk = evaluate_round142_eyring_transition_state_cftr_gating("kahn_kalai_gating_threshold")
        self.assertTrue(kk.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(kk.open_probability_p_o, 0.118)

    def test_potentiator_rescue_regimes(self) -> None:
        iva = evaluate_round142_eyring_transition_state_cftr_gating("ivacaftor_g551d_potentiator_rescue")
        tri = evaluate_round142_eyring_transition_state_cftr_gating("trikafta_f508del_dual_rescue")
        for res in (iva, tri):
            self.assertTrue(res.crosses_kahn_kalai_threshold)
            self.assertGreaterEqual(res.open_probability_p_o, 0.30)
            self.assertGreaterEqual(res.opening_rate_k_open_s_1, 0.50)

    def test_negative_control_falsification(self) -> None:
        normal = evaluate_falsification_negative_control(poison_negative_temperature=False)
        self.assertTrue(normal["passed"])
        self.assertEqual(normal["rule"], "AH-10")
        poisoned = evaluate_falsification_negative_control(poison_negative_temperature=True)
        self.assertFalse(poisoned["passed"])

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
