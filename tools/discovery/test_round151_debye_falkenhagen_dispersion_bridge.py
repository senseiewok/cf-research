#!/usr/bin/env python3
"""Unit tests for Round 151: Debye-Falkenhagen High-Frequency Ion Atmosphere Relaxation Dispersion Invariant Bridge (T-0806)."""
from __future__ import annotations

import unittest

from tools.discovery.round151_debye_falkenhagen_dispersion_bridge import (
    FIVE_LAYER_LADDER,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    compute_debye_falkenhagen_relaxation_time,
    compute_debye_length,
    compute_falkenhagen_dispersion_conductance,
    evaluate_falsification_negative_control,
    evaluate_round151_debye_falkenhagen_dispersion,
)


class Round151DebyeFalkenhagenDispersionTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round151_debye_falkenhagen_dispersion(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Failed for {key}")
            self.assertGreater(res.debye_length_nm, 0.0)
            self.assertGreater(res.relaxation_time_ns, 0.0)
            self.assertGreater(res.dc_chord_conductance_ps, 0.0)
            self.assertGreaterEqual(res.high_freq_chord_conductance_ps, res.dc_chord_conductance_ps)
            self.assertGreaterEqual(res.epithelial_short_circuit_current_ua_cm2, 0.0)

    def test_wild_type_cftr_pore_electrokinetics(self) -> None:
        wt = evaluate_round151_debye_falkenhagen_dispersion("wild_type_cftr_pore_electrokinetics")
        self.assertTrue(wt.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(wt.dc_chord_conductance_ps, 8.0)
        self.assertGreaterEqual(wt.epithelial_short_circuit_current_ua_cm2, 30.0)
        self.assertAlmostEqual(wt.fractional_wild_type_current, 1.0, delta=0.05)
        self.assertAlmostEqual(wt.debye_length_nm, 0.895, delta=0.01)
        self.assertAlmostEqual(wt.relaxation_time_ns, 0.197, delta=0.02)
        self.assertGreater(wt.high_freq_chord_conductance_ps, wt.dc_chord_conductance_ps)

    def test_cf_defects_and_gating_collapse(self) -> None:
        g551d = evaluate_round151_debye_falkenhagen_dispersion("cf_g551d_vestibule_dispersion_collapse")
        df508 = evaluate_round151_debye_falkenhagen_dispersion("cf_deltaf508_misfolded_permeation_stasis")
        for mutant in (g551d, df508):
            self.assertFalse(mutant.crosses_kahn_kalai_threshold)
            self.assertLessEqual(mutant.fractional_wild_type_current, 0.05)
            self.assertLessEqual(mutant.epithelial_short_circuit_current_ua_cm2, 0.5)
            self.assertLessEqual(mutant.dc_chord_conductance_ps, 3.0)

    def test_kahn_kalai_electrokinetic_boundary(self) -> None:
        kk = evaluate_round151_debye_falkenhagen_dispersion("kahn_kalai_electrokinetic_boundary")
        self.assertTrue(kk.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(kk.dc_chord_conductance_ps, 4.0)
        self.assertGreaterEqual(kk.epithelial_short_circuit_current_ua_cm2, 3.75)
        self.assertGreaterEqual(kk.fractional_wild_type_current, 0.118)

    def test_trikafta_allosteric_dispersion_rescue(self) -> None:
        rescued = evaluate_round151_debye_falkenhagen_dispersion("trikafta_allosteric_dispersion_rescue")
        self.assertTrue(rescued.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(rescued.dc_chord_conductance_ps, 8.0)
        self.assertGreaterEqual(rescued.epithelial_short_circuit_current_ua_cm2, 25.0)
        self.assertGreaterEqual(rescued.fractional_wild_type_current, 0.85)

    def test_negative_control_falsification(self) -> None:
        normal = evaluate_falsification_negative_control(poison_unphysical_relaxation=False)
        self.assertTrue(normal["passed"])
        self.assertEqual(normal["rule"], "AH-10")
        poisoned = evaluate_falsification_negative_control(poison_unphysical_relaxation=True)
        self.assertFalse(poisoned["passed"])

    def test_debye_length_and_relaxation_scaling(self) -> None:
        lam_120 = compute_debye_length(0.120)
        lam_60 = compute_debye_length(0.060)
        # Lower ionic strength means larger Debye length
        self.assertGreater(lam_60, lam_120)

        tau_120 = compute_debye_falkenhagen_relaxation_time(lam_120, 2.03e-9)
        tau_60 = compute_debye_falkenhagen_relaxation_time(lam_60, 2.03e-9)
        self.assertGreater(tau_60, tau_120)

        # High frequency dispersion monotonicity
        g_dc = compute_falkenhagen_dispersion_conductance(8.40, omega_rad_s=0.0, tau_df_s=tau_120)
        g_mid = compute_falkenhagen_dispersion_conductance(8.40, omega_rad_s=1.0 / tau_120, tau_df_s=tau_120)
        g_hi = compute_falkenhagen_dispersion_conductance(8.40, omega_rad_s=1e12, tau_df_s=tau_120)
        self.assertAlmostEqual(g_dc, 8.40, places=4)
        self.assertGreater(g_mid, g_dc)
        self.assertGreater(g_hi, g_mid)
        self.assertAlmostEqual(g_hi, 8.40 * 1.10, delta=0.01)

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
