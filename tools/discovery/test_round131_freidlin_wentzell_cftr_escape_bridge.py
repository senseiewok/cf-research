#!/usr/bin/env python3
"""Unit tests for Round 131: Freidlin-Wentzell Large Deviation Theory & Non-Equilibrium Metastable Basin Escape Invariant Bridge for CFTR Channel Gating."""
from __future__ import annotations

import unittest

from tools.discovery.round131_freidlin_wentzell_cftr_escape_bridge import (
    FIVE_LAYER_LADDER,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    evaluate_falsification_negative_control,
    evaluate_round131_freidlin_wentzell_cftr_escape,
)


class Round131FreidlinWentzellCftrEscapeTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round131_freidlin_wentzell_cftr_escape(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Failed for {key}")
            self.assertLess(res.complex_step_sensitivity, 0.0)
            self.assertGreater(res.quasipotential_action_kb_t, 0.0)
            self.assertGreater(res.k_escape_s_inv, 0.0)
            self.assertGreater(res.k_close_s_inv, 0.0)
            self.assertGreater(res.tau_closed_ms, 0.0)
            self.assertGreater(res.tau_open_ms, 0.0)
            self.assertGreater(res.open_probability, 0.0)
            self.assertLess(res.open_probability, 1.0)

    def test_wild_type_healthy_gating(self) -> None:
        wt = evaluate_round131_freidlin_wentzell_cftr_escape("wild_type_healthy")
        self.assertTrue(wt.crosses_kahn_kalai_threshold)
        self.assertLessEqual(wt.quasipotential_action_kb_t, 4.0)
        self.assertGreaterEqual(wt.open_probability, 0.35)
        self.assertLessEqual(wt.tau_closed_ms, 600.0)
        self.assertGreaterEqual(wt.tau_open_ms, 200.0)

    def test_cf_f508del_untreated_closure_stasis(self) -> None:
        cf = evaluate_round131_freidlin_wentzell_cftr_escape("cf_f508del_untreated")
        self.assertFalse(cf.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(cf.quasipotential_action_kb_t, 8.0)
        self.assertLessEqual(cf.open_probability, 0.02)
        self.assertGreaterEqual(cf.tau_closed_ms, 8000.0)
        self.assertLessEqual(cf.tau_open_ms, 80.0)

    def test_kahn_kalai_threshold_rescue(self) -> None:
        kk = evaluate_round131_freidlin_wentzell_cftr_escape("kahn_kalai_threshold_11_8pct")
        self.assertTrue(kk.crosses_kahn_kalai_threshold)
        self.assertLessEqual(kk.quasipotential_action_kb_t, 5.5)
        self.assertGreaterEqual(kk.open_probability, 0.10)
        self.assertLessEqual(kk.tau_closed_ms, 1600.0)

    def test_potentiators_rescue_open_probability(self) -> None:
        iva = evaluate_round131_freidlin_wentzell_cftr_escape("ivacaftor_potentiated")
        tri = evaluate_round131_freidlin_wentzell_cftr_escape("trikafta_elex_tez_iva_rescued")
        for res in (iva, tri):
            self.assertTrue(res.crosses_kahn_kalai_threshold)
            self.assertLessEqual(res.quasipotential_action_kb_t, 4.0)
            self.assertGreaterEqual(res.open_probability, 0.40)
            self.assertLessEqual(res.tau_closed_ms, 600.0)
            self.assertGreaterEqual(res.tau_open_ms, 300.0)

    def test_negative_control_falsification(self) -> None:
        normal = evaluate_falsification_negative_control(poison_negative_barrier=False)
        self.assertTrue(normal["passed"])
        self.assertEqual(normal["rule"], "AH-10")
        poisoned = evaluate_falsification_negative_control(poison_negative_barrier=True)
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
