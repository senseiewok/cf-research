#!/usr/bin/env python3
"""Unit tests for Round 184: CFTR Waiting-Time & Entropy Production Bridge (T-0861)."""
from __future__ import annotations

import unittest

from tools.discovery.round184_cftr_waiting_time_entropy_production_bridge import (
    FIVE_LAYER_LADDER,
    PARAMETER_PROVENANCE,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    complex_step_po_derivative,
    compute_closed_duration_ms,
    compute_coefficient_of_variation,
    compute_open_probability,
    compute_skinner_dunkel_entropy_per_cycle,
    evaluate_falsification_negative_control,
    evaluate_round184_cftr_waiting_time_entropy_production,
)


class Round184CftrWaitingTimeEntropyProductionTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round184_cftr_waiting_time_entropy_production(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Invariants failed for {key}")
            self.assertGreaterEqual(res.k_eff_shape_parameter, 1.0)
            self.assertGreaterEqual(res.coefficient_of_variation_cv, 0.0)
            self.assertLessEqual(res.coefficient_of_variation_cv, 1.0)
            self.assertGreater(res.tau_open_mean_ms, 0.0)
            self.assertGreater(res.tau_closed_mean_ms, 0.0)
            self.assertGreaterEqual(res.channel_open_probability, 0.0)
            self.assertLessEqual(res.channel_open_probability, 1.0)
            self.assertGreaterEqual(res.skinner_dunkel_entropy_per_cycle_kb, 0.0)
            self.assertGreaterEqual(res.entropy_production_rate_kb_s, 0.0)
            self.assertLessEqual(res.complex_step_relative_error, 1e-6)

    def test_wild_type_saturating_atp_crosses_kahn_kalai(self) -> None:
        res = evaluate_round184_cftr_waiting_time_entropy_production("wild_type_saturating_atp_catalytic_cycle")
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res.k_eff_shape_parameter, 1.80)
        self.assertLessEqual(res.coefficient_of_variation_cv, 0.75)
        self.assertGreaterEqual(res.skinner_dunkel_entropy_per_cycle_kb, 1.0)
        self.assertGreaterEqual(res.tau_open_mean_ms, 350.0)
        self.assertGreaterEqual(res.channel_open_probability, 0.30)

    def test_ivacaftor_potentiated_crosses_kahn_kalai(self) -> None:
        res = evaluate_round184_cftr_waiting_time_entropy_production("ivacaftor_potentiated_wt_gating")
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res.channel_open_probability, 0.70)
        self.assertGreaterEqual(res.tau_open_mean_ms, 1000.0)
        self.assertGreaterEqual(res.skinner_dunkel_entropy_per_cycle_kb, 1.0)

    def test_sub_curative_regimes_fail_kahn_kalai(self) -> None:
        # G551D mutant: decoupled gating, c_v ~ 1, P_o = 0.005
        g551d = evaluate_round184_cftr_waiting_time_entropy_production("g551d_uncorrected_catalytic_decoupled_gating")
        self.assertFalse(g551d.crosses_kahn_kalai_threshold)
        self.assertLess(g551d.channel_open_probability, 0.05)
        self.assertLess(g551d.skinner_dunkel_entropy_per_cycle_kb, 0.20)
        self.assertGreater(g551d.coefficient_of_variation_cv, 0.90)

        # ATP-depleted baseline
        low_atp = evaluate_round184_cftr_waiting_time_entropy_production("atp_depleted_sub_phys_wt_baseline")
        self.assertFalse(low_atp.crosses_kahn_kalai_threshold)
        self.assertLess(low_atp.channel_open_probability, 0.20)

        # Equilibrium AMP-PNP arrest: c_v = 1.0, Delta S = 0
        amp = evaluate_round184_cftr_waiting_time_entropy_production("equilibrium_unhydrolyzable_amp_pnp_arrest")
        self.assertFalse(amp.crosses_kahn_kalai_threshold)
        self.assertEqual(amp.skinner_dunkel_entropy_per_cycle_kb, 0.0)
        self.assertEqual(amp.coefficient_of_variation_cv, 1.0)

    def test_complex_step_derivative(self) -> None:
        cs_deriv, fd_deriv, rel_err = complex_step_po_derivative(atp_conc_mm=0.50)
        self.assertLessEqual(rel_err, 1e-6)
        self.assertGreater(cs_deriv, 0.0)

    def test_negative_control_falsification(self) -> None:
        neg_atp = evaluate_falsification_negative_control(poison_negative_atp=True)
        self.assertTrue(neg_atp["passed"])
        self.assertIn("Correctly caught", neg_atp["error"])

        neg_k = evaluate_falsification_negative_control(poison_unphysical_k_eff=True)
        self.assertTrue(neg_k["passed"])
        self.assertIn("Correctly caught", neg_k["error"])

        neg_dur = evaluate_falsification_negative_control(poison_negative_duration=True)
        self.assertTrue(neg_dur["passed"])
        self.assertIn("Correctly caught", neg_dur["error"])

        baseline = evaluate_falsification_negative_control()
        self.assertTrue(baseline["passed"])
        self.assertIsNone(baseline["error"])

    def test_provenance_and_citations(self) -> None:
        self.assertGreaterEqual(len(PRIMARY_CITATIONS), 3)
        self.assertIn("DOI: 10.1073/pnas.0911061107", PRIMARY_CITATIONS)
        self.assertIn("DOI: 10.1016/j.cell.2015.09.052", PRIMARY_CITATIONS)
        self.assertIn("DOI: 10.1103/PhysRevLett.127.198101", PRIMARY_CITATIONS)

        self.assertGreaterEqual(len(PARAMETER_PROVENANCE), 5)
        for key, meta in PARAMETER_PROVENANCE.items():
            self.assertIn("tag", meta)
            self.assertIn("source", meta)
            self.assertIn("bounds", meta)

        self.assertIn("deductive_theorem", TWO_COLUMN_EPISTEMIC_HYGIENE)
        self.assertIn("biophysical_hypothesis", TWO_COLUMN_EPISTEMIC_HYGIENE)
        self.assertIn("cool_human_fact", FIVE_LAYER_LADDER)
        self.assertIn("falsifiable_code_and_toy", FIVE_LAYER_LADDER)


if __name__ == "__main__":
    unittest.main()
