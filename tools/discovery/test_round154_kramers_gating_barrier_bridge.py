#!/usr/bin/env python3
"""Unit tests for Round 154: Stochastic Fluctuating Force & Kramers Barrier Crossing Invariant Bridge (T-0811)."""
from __future__ import annotations

import math
import unittest

from tools.discovery.round154_kramers_gating_barrier_bridge import (
    FIVE_LAYER_LADDER,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    compute_kramers_rate,
    evaluate_falsification_negative_control,
    evaluate_round154_kramers_gating_barrier,
)


class Round154KramersGatingTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round154_kramers_gating_barrier(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Failed for {key}")
            self.assertGreater(res.kramers_prefactor_s1, 0.0)
            self.assertGreater(res.opening_rate_s1, 0.0)
            self.assertGreater(res.closing_rate_s1, 0.0)
            self.assertGreater(res.open_probability, 0.0)
            self.assertLess(res.open_probability, 1.0)
            self.assertGreater(res.mean_burst_open_duration_ms, 0.0)
            self.assertGreater(res.mean_interburst_closed_duration_ms, 0.0)
            self.assertGreaterEqual(res.epithelial_short_circuit_current_ua_cm2, 0.0)

    def test_wild_type_phosphorylated_gating(self) -> None:
        wt = evaluate_round154_kramers_gating_barrier("wild_type_phosphorylated_gating")
        self.assertTrue(wt.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(wt.open_probability, 0.40)
        self.assertGreaterEqual(wt.epithelial_short_circuit_current_ua_cm2, 30.0)
        self.assertAlmostEqual(wt.fractional_wild_type_current, 1.0, delta=0.05)
        self.assertLess(wt.free_energy_gap_kt, 1.0)  # Moderate free-energy gap permitting physiological burst gating

    def test_cf_gating_defects(self) -> None:
        g551d = evaluate_round154_kramers_gating_barrier("cf_g551d_gating_barrier_collapse")
        df508 = evaluate_round154_kramers_gating_barrier("cf_deltaf508_misfolded_residual_gating")

        for mutant in (g551d, df508):
            self.assertFalse(mutant.crosses_kahn_kalai_threshold)
            self.assertLess(mutant.epithelial_short_circuit_current_ua_cm2, 1.0)
            self.assertLess(mutant.fractional_wild_type_current, 0.05)

        # G551D specifically has high opening barrier and low P_open
        self.assertLess(g551d.open_probability, 0.01)
        self.assertGreater(g551d.free_energy_gap_kt, 5.0)

    def test_kahn_kalai_stochastic_gating_boundary(self) -> None:
        kk = evaluate_round154_kramers_gating_barrier("kahn_kalai_stochastic_gating_boundary")
        self.assertTrue(kk.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(kk.epithelial_short_circuit_current_ua_cm2, 3.75)
        self.assertGreaterEqual(kk.fractional_wild_type_current, 0.118)

    def test_ivacaftor_potentiated_rescue(self) -> None:
        rescued = evaluate_round154_kramers_gating_barrier("ivacaftor_potentiated_g551d_rescue")
        self.assertTrue(rescued.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(rescued.open_probability, 0.35)
        self.assertGreaterEqual(rescued.epithelial_short_circuit_current_ua_cm2, 25.0)
        self.assertGreaterEqual(rescued.fractional_wild_type_current, 0.85)

    def test_negative_control_falsification(self) -> None:
        normal = evaluate_falsification_negative_control(poison_non_positive_friction=False)
        self.assertTrue(normal["passed"])
        self.assertEqual(normal["rule"], "AH-10")

        poisoned_gamma = evaluate_falsification_negative_control(poison_non_positive_friction=True)
        self.assertFalse(poisoned_gamma["passed"])

        poisoned_barrier = evaluate_falsification_negative_control(poison_negative_barrier=True)
        self.assertFalse(poisoned_barrier["passed"])

    def test_kramers_rate_scaling_and_detailed_balance(self) -> None:
        omega_min = 6.283e7
        omega_bar = 6.283e7
        gamma = 1.0e8
        temp_k = 310.15
        kb = 1.380649e-23
        kb_t = kb * temp_k

        # Lower barrier means higher rate
        r_low = compute_kramers_rate(6.0 * kb_t, omega_min, omega_bar, gamma, temp_k)
        r_high = compute_kramers_rate(12.0 * kb_t, omega_min, omega_bar, gamma, temp_k)
        self.assertGreater(r_low, r_high)
        self.assertAlmostEqual(r_low / r_high, math.exp(6.0), places=4)

        # Higher friction means lower rate
        r_visc = compute_kramers_rate(6.0 * kb_t, omega_min, omega_bar, 2.0 * gamma, temp_k)
        self.assertAlmostEqual(r_low / r_visc, 2.0, places=4)

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
            self.assertIn(layer, FIVE_LAYER_LADDER)
            self.assertGreaterEqual(len(FIVE_LAYER_LADDER[layer]), 25)
        self.assertGreaterEqual(len(PRIMARY_CITATIONS), 4)


if __name__ == "__main__":
    unittest.main()
