#!/usr/bin/env python3
"""Unit tests for Round 186: CFTR NBD1-NBD2 Heterodimer Dissociation Free-Energy Bridge (T-0863)."""
from __future__ import annotations

import unittest

from tools.discovery.round186_nbd_dimer_dissociation_free_energy_bridge import (
    DEFAULT_KD_SITE1_MM,
    DEFAULT_KD_SITE2_MM,
    FIVE_LAYER_LADDER,
    PARAMETER_PROVENANCE,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    complex_step_po_derivative,
    compute_channel_open_probability,
    compute_closed_duration_ms,
    compute_dimer_dissociation_barrier_kbt,
    compute_open_burst_duration_ms,
    compute_site_occupancy,
    evaluate_falsification_negative_control,
    evaluate_round186_nbd_dimer_dissociation_free_energy,
)


class Round186NbdDimerDissociationFreeEnergyTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round186_nbd_dimer_dissociation_free_energy(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Invariants failed for {key}")
            self.assertGreater(res.delta_g_dissoc_kbt, 0.0)
            self.assertGreater(res.free_energy_joules, 0.0)
            self.assertGreaterEqual(res.site1_occupancy_theta1, 0.0)
            self.assertLessEqual(res.site1_occupancy_theta1, 1.0)
            self.assertGreaterEqual(res.site2_occupancy_theta2, 0.0)
            self.assertLessEqual(res.site2_occupancy_theta2, 1.0)
            self.assertGreater(res.tau_open_mean_ms, 0.0)
            self.assertGreater(res.tau_closed_mean_ms, 0.0)
            self.assertGreaterEqual(res.channel_open_probability, 0.0)
            self.assertLessEqual(res.channel_open_probability, 1.0)
            self.assertGreaterEqual(res.catalytic_asymmetry_ratio, 5.0)
            self.assertLessEqual(res.complex_step_relative_error, 1e-6)

    def test_wild_type_saturating_atp_crosses_kahn_kalai(self) -> None:
        res = evaluate_round186_nbd_dimer_dissociation_free_energy("wild_type_saturating_atp_heterodimer")
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res.delta_g_dissoc_kbt, 14.0)
        self.assertGreaterEqual(res.tau_open_mean_ms, 350.0)
        self.assertGreaterEqual(res.channel_open_probability, 0.30)
        self.assertGreaterEqual(res.site1_occupancy_theta1, 0.95)
        self.assertGreaterEqual(res.site2_occupancy_theta2, 0.95)

    def test_trikafta_rescued_f508del_crosses_kahn_kalai(self) -> None:
        res = evaluate_round186_nbd_dimer_dissociation_free_energy("trikafta_rescued_f508del_heterodimer")
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res.delta_g_dissoc_kbt, 14.0)
        self.assertGreaterEqual(res.tau_open_mean_ms, 350.0)
        self.assertGreaterEqual(res.channel_open_probability, 0.30)

    def test_sub_curative_presets_fail_kahn_kalai(self) -> None:
        # Uncorrected F508del: collapsed barrier, premature thermal dissociation
        f508del = evaluate_round186_nbd_dimer_dissociation_free_energy("f508del_uncorrected_dimer_destabilization")
        self.assertFalse(f508del.crosses_kahn_kalai_threshold)
        self.assertLess(f508del.delta_g_dissoc_kbt, 14.0)
        self.assertLess(f508del.tau_open_mean_ms, 150.0)
        self.assertLess(f508del.channel_open_probability, 0.05)

        # ATP-depleted baseline
        low_atp = evaluate_round186_nbd_dimer_dissociation_free_energy("atp_depleted_subphys_baseline")
        self.assertFalse(low_atp.crosses_kahn_kalai_threshold)
        self.assertLess(low_atp.delta_g_dissoc_kbt, 12.0)
        self.assertLess(low_atp.channel_open_probability, 0.05)

        # G551D gating mutant
        g551d = evaluate_round186_nbd_dimer_dissociation_free_energy("g551d_catalytic_decoupled_dimer")
        self.assertFalse(g551d.crosses_kahn_kalai_threshold)
        self.assertLess(g551d.channel_open_probability, 0.01)

    def test_complex_step_po_derivative_accuracy(self) -> None:
        cs_deriv, fd_deriv, rel_err = complex_step_po_derivative(delta_g_dissoc_kbt=16.01, tau_closed_ms=683.3)
        self.assertGreater(cs_deriv, 0.0)
        self.assertGreater(fd_deriv, 0.0)
        self.assertLessEqual(rel_err, 1e-6)

    def test_poisoned_negative_controls_fail_deterministically(self) -> None:
        res_atp = evaluate_falsification_negative_control(poison_negative_atp=True)
        self.assertTrue(res_atp["passed"])
        self.assertIn("Correctly caught", res_atp["error"])

        res_temp = evaluate_falsification_negative_control(poison_negative_temperature=True)
        self.assertTrue(res_temp["passed"])
        self.assertIn("Correctly caught", res_temp["error"])

        res_barrier = evaluate_falsification_negative_control(poison_collapsed_barrier=True)
        self.assertTrue(res_barrier["passed"])
        self.assertIn("Correctly caught", res_barrier["error"])

        clean = evaluate_falsification_negative_control()
        self.assertTrue(clean["passed"])
        self.assertIsNone(clean["error"])

    def test_asymmetry_site1_vs_site2(self) -> None:
        self.assertLess(DEFAULT_KD_SITE1_MM, DEFAULT_KD_SITE2_MM)
        ratio = DEFAULT_KD_SITE2_MM / DEFAULT_KD_SITE1_MM
        self.assertGreaterEqual(ratio, 10.0)
        th1_low = compute_site_occupancy(0.010, DEFAULT_KD_SITE1_MM)
        th2_low = compute_site_occupancy(0.010, DEFAULT_KD_SITE2_MM)
        self.assertGreater(th1_low, th2_low)

    def test_provenance_and_two_column_epistemic_hygiene_metadata(self) -> None:
        self.assertGreaterEqual(len(PRIMARY_CITATIONS), 4)
        self.assertGreaterEqual(len(PARAMETER_PROVENANCE), 5)
        for param, meta in PARAMETER_PROVENANCE.items():
            self.assertIn("value", meta)
            self.assertIn("tag", meta)
            self.assertIn("source", meta)
            self.assertIn("bounds", meta)
            self.assertIn(meta["tag"], ("empirical", "derived", "phenomenological_heuristic"))

        self.assertIn("deductive_theorem", TWO_COLUMN_EPISTEMIC_HYGIENE)
        self.assertIn("biophysical_hypothesis", TWO_COLUMN_EPISTEMIC_HYGIENE)
        self.assertIn("statement", TWO_COLUMN_EPISTEMIC_HYGIENE["biophysical_hypothesis"])
        self.assertIn("falsification_condition", TWO_COLUMN_EPISTEMIC_HYGIENE["biophysical_hypothesis"])

        self.assertIn("cool_human_fact", FIVE_LAYER_LADDER)
        self.assertIn("everyday_analogy", FIVE_LAYER_LADDER)
        self.assertIn("engineering_lens", FIVE_LAYER_LADDER)
        self.assertIn("formal_theorem", FIVE_LAYER_LADDER)
        self.assertIn("falsifiable_code_and_toy", FIVE_LAYER_LADDER)


if __name__ == "__main__":
    unittest.main()
