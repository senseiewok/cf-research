#!/usr/bin/env python3
"""Unit tests for Round 136: Kramers Turnover & Non-Equilibrium Gating Transition Invariant Bridge."""
from __future__ import annotations

import unittest

from tools.discovery.round136_kramers_turnover_cftr_gating_bridge import (
    FIVE_LAYER_LADDER,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    compute_kramers_turnover_gating_kinetics,
    evaluate_falsification_negative_control,
    evaluate_round136_kramers_turnover_cftr_gating,
)


class Round136KramersTurnoverCftrGatingTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round136_kramers_turnover_cftr_gating(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Failed for {key}")
            self.assertLess(res.complex_step_sensitivity, 0.0)
            self.assertGreater(res.barrier_kb_t, 0.0)
            self.assertGreater(res.damping_gamma_s1, 0.0)
            self.assertGreater(res.k_open_s1, 0.0)
            self.assertGreater(res.p_open, 0.0)
            self.assertLess(res.p_open, 1.0)
            self.assertGreater(res.gamma_eff_ps, 0.0)

    def test_healthy_cf_dimerization_at_turnover(self) -> None:
        healthy = evaluate_round136_kramers_turnover_cftr_gating("wild_type_healthy")
        self.assertTrue(healthy.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(healthy.k_open_s1, 1.2)
        self.assertGreaterEqual(healthy.p_open, 0.35)
        self.assertGreaterEqual(healthy.gamma_eff_ps, 3.0)
        self.assertLessEqual(healthy.barrier_kb_t, 4.5)

    def test_f508del_gating_collapse(self) -> None:
        cf = evaluate_round136_kramers_turnover_cftr_gating("cf_f508del_untreated")
        self.assertFalse(cf.crosses_kahn_kalai_threshold)
        self.assertLess(cf.k_open_s1, 0.05)
        self.assertLess(cf.p_open, 0.02)
        self.assertLess(cf.gamma_eff_ps, 0.20)
        self.assertGreaterEqual(cf.barrier_kb_t, 8.5)

    def test_kahn_kalai_threshold_rescue(self) -> None:
        kk = evaluate_round136_kramers_turnover_cftr_gating("kahn_kalai_threshold_11_8pct")
        self.assertTrue(kk.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(kk.gamma_eff_ps, 1.0)
        self.assertGreaterEqual(kk.p_open, 0.118)
        self.assertLessEqual(kk.barrier_kb_t, 6.5)

    def test_ivacaftor_and_trikafta_rescue(self) -> None:
        iva = evaluate_round136_kramers_turnover_cftr_gating("ivacaftor_potentiated")
        tri = evaluate_round136_kramers_turnover_cftr_gating("trikafta_full_dimerization_rescue")
        for res in (iva, tri):
            self.assertTrue(res.crosses_kahn_kalai_threshold)
            self.assertGreaterEqual(res.k_open_s1, 1.0)
            self.assertGreaterEqual(res.p_open, 0.35)
            self.assertGreaterEqual(res.gamma_eff_ps, 3.0)
            self.assertLessEqual(res.barrier_kb_t, 5.0)

    def test_kramers_turnover_asymptotics(self) -> None:
        # Non-monotonicity in friction: Rate goes down at both extremes of damping
        barrier = 4.20
        r_low = compute_kramers_turnover_gating_kinetics(barrier, gamma_s1=1e8)["r_kramers"]
        r_mid = compute_kramers_turnover_gating_kinetics(barrier, gamma_s1=1e11)["r_kramers"]
        r_high = compute_kramers_turnover_gating_kinetics(barrier, gamma_s1=1e14)["r_kramers"]
        self.assertGreater(r_mid, r_low, "Turnover peak must exceed low-damping energy diffusion rate")
        self.assertGreater(r_mid, r_high, "Turnover peak must exceed high-damping spatial diffusion rate")

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
