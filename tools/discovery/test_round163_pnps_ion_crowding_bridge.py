#!/usr/bin/env python3
"""Unit tests for Round 163: Poisson-Nernst-Planck-Stern (PNPS) Finite Ion Size & Crowding Invariant Bridge (T-0828)."""
from __future__ import annotations

import unittest

from tools.discovery.round163_pnps_ion_crowding_bridge import (
    FIVE_LAYER_LADDER,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    complex_step_conductance_derivative,
    compute_anomalous_mole_fraction_effect,
    compute_maximum_packing_concentration_m3,
    compute_packing_fraction,
    compute_single_channel_current,
    compute_unitary_conductance,
    evaluate_falsification_negative_control,
    evaluate_round163_pnps_ion_crowding,
)


class Round163PNPSIonCrowdingTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round163_pnps_ion_crowding(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Failed for {key}")
            self.assertGreater(res.unitary_conductance_ps, 0.0)
            self.assertLessEqual(res.complex_step_relative_error, 1e-10)

    def test_wt_cftr_physiological_pnps(self) -> None:
        res = evaluate_round163_pnps_ion_crowding("wt_cftr_physiological_pnps")
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res.unitary_conductance_ps, 7.0)
        self.assertLess(res.unitary_conductance_ps, 9.60)
        self.assertGreater(res.packing_fraction, 0.20)
        self.assertLess(res.packing_fraction, 0.70)
        self.assertLess(res.single_channel_current_pa, 0.0)  # Inward current / negative driving force

    def test_trikafta_rescued_f508del_pnps(self) -> None:
        res = evaluate_round163_pnps_ion_crowding("trikafta_rescued_f508del_pnps")
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res.unitary_conductance_ps, 6.5)

    def test_kahn_kalai_conductance_boundary(self) -> None:
        res = evaluate_round163_pnps_ion_crowding("kahn_kalai_conductance_boundary")
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res.unitary_conductance_ps, 1.0)

    def test_untreated_r334w_pore_mutant(self) -> None:
        res = evaluate_round163_pnps_ion_crowding("untreated_r334w_pore_mutant")
        self.assertFalse(res.crosses_kahn_kalai_threshold)
        self.assertLess(res.unitary_conductance_ps, 1.0)

    def test_extreme_hyperconcentrated_crowding(self) -> None:
        res = evaluate_round163_pnps_ion_crowding("extreme_hyperconcentrated_crowding")
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res.unitary_conductance_ps, 8.5)
        self.assertGreater(res.packing_fraction, 0.70)

    def test_conductance_monotonicity_and_saturation(self) -> None:
        g_low = compute_unitary_conductance(10.0, 9.60, 38.0)
        g_mid = compute_unitary_conductance(140.0, 9.60, 38.0)
        g_high = compute_unitary_conductance(500.0, 9.60, 38.0)
        self.assertGreater(g_mid, g_low)
        self.assertGreater(g_high, g_mid)
        # Strict concavity: delta from low->mid vs mid->high
        slope_1 = (g_mid - g_low) / (140.0 - 10.0)
        slope_2 = (g_high - g_mid) / (500.0 - 140.0)
        self.assertGreater(slope_1, slope_2)

    def test_complex_step_conductance_sensitivity(self) -> None:
        cs_deriv, analytical, rel_err = complex_step_conductance_derivative(
            chloride_conc_mm=140.0,
            g_max_ps=9.60,
            k_m_mm=38.0,
        )
        self.assertLessEqual(rel_err, 1e-10)
        self.assertGreater(cs_deriv, 0.0)

    def test_anomalous_mole_fraction_effect(self) -> None:
        g_pure_cl = compute_anomalous_mole_fraction_effect(0.0)
        g_pure_comp = compute_anomalous_mole_fraction_effect(1.0)
        g_mix_half = compute_anomalous_mole_fraction_effect(0.5)
        # AMFE predicts non-linear dip at intermediate mole fraction
        linear_expected = 0.5 * (g_pure_cl + g_pure_comp)
        self.assertLess(g_mix_half, linear_expected)

    def test_high_field_current_saturation(self) -> None:
        # At very large voltage, current should saturate asymptotically rather than exploding linearly
        v_rev = -36.0
        v1 = -40.0
        v2 = -300.0
        df1 = abs(v1 - v_rev)  # 4.0 mV
        df2 = abs(v2 - v_rev)  # 264.0 mV
        linear_ratio = df2 / df1  # 66.0

        i1 = compute_single_channel_current(7.5, v1, reversal_potential_mv=v_rev)
        i2 = compute_single_channel_current(7.5, v2, reversal_potential_mv=v_rev)
        actual_ratio = abs(i2) / abs(i1)
        self.assertLess(actual_ratio, linear_ratio)

    def test_negative_control_falsification(self) -> None:
        # Poisoned parameter checks
        neg1 = evaluate_falsification_negative_control(poison_negative_chloride_conc=True)
        self.assertFalse(neg1["passed"])
        self.assertIn("poisoned", neg1["error"])

        neg2 = evaluate_falsification_negative_control(poison_negative_hydrated_diameter=True)
        self.assertFalse(neg2["passed"])

        neg3 = evaluate_falsification_negative_control(poison_negative_g_max=True)
        self.assertFalse(neg3["passed"])

        neg4 = evaluate_falsification_negative_control(poison_negative_km=True)
        self.assertFalse(neg4["passed"])

        # Unpoisoned base check
        clean = evaluate_falsification_negative_control()
        self.assertTrue(clean["passed"])

    def test_epistemic_hygiene_and_ladder_structure(self) -> None:
        self.assertIn("deductive_theorem", TWO_COLUMN_EPISTEMIC_HYGIENE)
        self.assertIn("biophysical_hypothesis", TWO_COLUMN_EPISTEMIC_HYGIENE)
        self.assertGreaterEqual(len(FIVE_LAYER_LADDER), 5)
        self.assertGreaterEqual(len(PRIMARY_CITATIONS), 4)


if __name__ == "__main__":
    unittest.main()
