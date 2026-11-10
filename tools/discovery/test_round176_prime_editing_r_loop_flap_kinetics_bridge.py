#!/usr/bin/env python3
"""Unit tests for Round 176: Prime Editing pegRNA R-Loop Flap Kinetics & FEN1 Flap Cleavage Invariant Bridge (T-0847)."""
from __future__ import annotations

import unittest

from tools.discovery.round176_prime_editing_r_loop_flap_kinetics_bridge import (
    FIVE_LAYER_LADDER,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    complex_step_fen1_derivative,
    compute_flap_equilibrium_ratio,
    compute_net_editing_efficiency,
    compute_pbs_annealing_fraction,
    compute_rt_elongation_fraction,
    evaluate_falsification_negative_control,
    evaluate_round176_prime_editing_flap_kinetics,
)


class Round176PrimeEditingTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round176_prime_editing_flap_kinetics(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Invariants failed for {key}")
            self.assertGreaterEqual(res.pbs_annealing_fraction, 0.0)
            self.assertLessEqual(res.pbs_annealing_fraction, 1.0)
            self.assertGreaterEqual(res.rt_elongation_fraction, 0.0)
            self.assertLessEqual(res.rt_elongation_fraction, 1.0)
            self.assertGreaterEqual(res.flap_equilibrium_ratio, 0.0)
            self.assertLessEqual(res.flap_equilibrium_ratio, 1.0)
            self.assertGreater(res.fen1_cleavage_rate_s, 0.0)
            self.assertGreaterEqual(res.mmr_reversion_rate_s, 0.0)
            self.assertGreaterEqual(res.net_editing_efficiency, 0.0)
            self.assertLessEqual(res.net_editing_efficiency, 1.0)
            self.assertGreaterEqual(res.bystander_indel_rate, 0.0)
            self.assertLessEqual(res.bystander_indel_rate, 0.05)
            self.assertLessEqual(res.complex_step_relative_error, 1e-6)

    def test_pemax_epegRNA_cftr_w1282x_rescue(self) -> None:
        res = evaluate_round176_prime_editing_flap_kinetics("pemax_epegRNA_cftr_w1282x_rescue")
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res.net_editing_efficiency, 0.30)
        self.assertGreaterEqual(res.pbs_annealing_fraction, 0.95)
        self.assertGreaterEqual(res.flap_equilibrium_ratio, 0.70)
        self.assertLessEqual(res.bystander_indel_rate, 0.015)

    def test_standard_pe3_f508del_insertion(self) -> None:
        res = evaluate_round176_prime_editing_flap_kinetics("standard_pe3_f508del_insertion")
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res.net_editing_efficiency, 0.18)
        self.assertGreaterEqual(res.pbs_annealing_fraction, 0.90)

    def test_unengineered_pe2_active_mmr(self) -> None:
        res = evaluate_round176_prime_editing_flap_kinetics("unengineered_pe2_active_mmr")
        self.assertFalse(res.crosses_kahn_kalai_threshold)
        self.assertLess(res.net_editing_efficiency, 0.18)

    def test_suboptimal_short_pbs_7nt(self) -> None:
        res = evaluate_round176_prime_editing_flap_kinetics("suboptimal_short_pbs_7nt")
        self.assertFalse(res.crosses_kahn_kalai_threshold)
        self.assertLess(res.pbs_annealing_fraction, 0.50)
        self.assertLess(res.net_editing_efficiency, 0.10)

    def test_fen1_deficient_flaps_accumulate(self) -> None:
        res = evaluate_round176_prime_editing_flap_kinetics("fen1_deficient_flaps_accumulate")
        self.assertFalse(res.crosses_kahn_kalai_threshold)
        self.assertLess(res.net_editing_efficiency, 0.10)

    def test_complex_step_fen1_derivative(self) -> None:
        cs_deriv, fd_deriv, rel_err = complex_step_fen1_derivative(
            k_fen1=0.45,
            k_mmr=0.04,
        )
        self.assertLessEqual(rel_err, 1e-6)
        # Higher FEN1 cleavage rate increases net editing efficiency
        self.assertGreater(cs_deriv, 0.0)

    def test_negative_control_falsification(self) -> None:
        neg_fen1 = evaluate_falsification_negative_control(poison_negative_fen1=True)
        self.assertTrue(neg_fen1["passed"])
        self.assertIn("Correctly caught", neg_fen1["error"])

        neg_rt = evaluate_falsification_negative_control(poison_negative_rt=True)
        self.assertTrue(neg_rt["passed"])
        self.assertIn("Correctly caught", neg_rt["error"])

        neg_temp = evaluate_falsification_negative_control(poison_negative_temp=True)
        self.assertTrue(neg_temp["passed"])
        self.assertIn("Correctly caught", neg_temp["error"])

        baseline = evaluate_falsification_negative_control()
        self.assertTrue(baseline["passed"])
        self.assertIsNone(baseline["error"])


if __name__ == "__main__":
    unittest.main()
