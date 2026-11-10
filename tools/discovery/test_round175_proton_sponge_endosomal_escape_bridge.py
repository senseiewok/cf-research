#!/usr/bin/env python3
"""Unit tests for Round 175: Proton-Sponge Endosomal Escape & Inverted Hexagonal H_II Phase Transition Invariant Bridge (T-0846)."""
from __future__ import annotations

import unittest

from tools.discovery.round175_proton_sponge_endosomal_escape_bridge import (
    FIVE_LAYER_LADDER,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    complex_step_escape_derivative,
    compute_critical_packing_parameter,
    compute_cytosolic_escape_fraction,
    compute_escape_rate,
    compute_inverted_hexagonal_fraction,
    compute_lipid_protonation_fraction,
    compute_osmotic_pressure_bar,
    evaluate_falsification_negative_control,
    evaluate_round175_proton_sponge_escape,
)


class Round175ProtonSpongeTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round175_proton_sponge_escape(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Invariants failed for {key}")
            self.assertGreater(res.apparent_pka, 0.0)
            self.assertLess(res.apparent_pka, 14.0)
            self.assertGreater(res.endosomal_ph, 0.0)
            self.assertLess(res.endosomal_ph, 14.0)
            self.assertGreaterEqual(res.protonation_fraction_endo, 0.0)
            self.assertLessEqual(res.protonation_fraction_endo, 1.0)
            self.assertGreaterEqual(res.protonation_fraction_extra, 0.0)
            self.assertLessEqual(res.protonation_fraction_extra, 1.0)
            self.assertGreaterEqual(res.osmotic_pressure_bar, 0.0)
            self.assertGreater(res.critical_packing_parameter, 0.0)
            self.assertGreaterEqual(res.inverted_hexagonal_fraction, 0.0)
            self.assertLessEqual(res.inverted_hexagonal_fraction, 1.0)
            self.assertGreater(res.escape_rate_s, 0.0)
            self.assertGreater(res.lysosomal_degradation_rate_s, 0.0)
            self.assertGreaterEqual(res.cytosolic_escape_fraction, 0.0)
            self.assertLessEqual(res.cytosolic_escape_fraction, 1.0)
            self.assertLessEqual(res.complex_step_relative_error, 1e-6)

    def test_optimized_ionizable_lnp_late_endosome(self) -> None:
        res = evaluate_round175_proton_sponge_escape("optimized_ionizable_lnp_late_endosome")
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res.cytosolic_escape_fraction, 0.035)
        self.assertGreaterEqual(res.critical_packing_parameter, 1.05)
        self.assertLessEqual(res.protonation_fraction_extra, 0.15)
        self.assertGreaterEqual(res.protonation_fraction_endo, 0.85)

    def test_standard_lnp_mid_endosome(self) -> None:
        res = evaluate_round175_proton_sponge_escape("standard_lnp_mid_endosome")
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res.cytosolic_escape_fraction, 0.020)
        self.assertGreaterEqual(res.critical_packing_parameter, 1.00)

    def test_suboptimal_low_pka_lnp(self) -> None:
        res = evaluate_round175_proton_sponge_escape("suboptimal_low_pka_lnp")
        self.assertFalse(res.crosses_kahn_kalai_threshold)
        self.assertLess(res.critical_packing_parameter, 1.00)

    def test_high_pka_cationic_toxic_lnp(self) -> None:
        res = evaluate_round175_proton_sponge_escape("high_pka_cationic_toxic_lnp")
        self.assertFalse(res.crosses_kahn_kalai_threshold)
        # Violates extracellular neutral stealth boundary
        self.assertGreater(res.protonation_fraction_extra, 0.20)

    def test_lbpa_deficient_endosome(self) -> None:
        res = evaluate_round175_proton_sponge_escape("lbpa_deficient_endosome")
        self.assertFalse(res.crosses_kahn_kalai_threshold)
        self.assertLess(res.cytosolic_escape_fraction, 0.020)
        self.assertLess(res.critical_packing_parameter, 1.00)

    def test_complex_step_escape_derivative(self) -> None:
        cs_deriv, fd_deriv, rel_err = complex_step_escape_derivative(
            ph_endo=5.50,
            pka=6.45,
        )
        self.assertLessEqual(rel_err, 1e-6)
        # Endosomal acidification (lower pH) increases escape fraction, so derivative wrt pH is negative
        self.assertLess(cs_deriv, 0.0)

    def test_negative_control_falsification(self) -> None:
        neg_ph = evaluate_falsification_negative_control(poison_negative_ph=True)
        self.assertTrue(neg_ph["passed"])
        self.assertIn("Correctly caught", neg_ph["error"])

        neg_pka = evaluate_falsification_negative_control(poison_negative_pka=True)
        self.assertTrue(neg_pka["passed"])
        self.assertIn("Correctly caught", neg_pka["error"])

        neg_rate = evaluate_falsification_negative_control(poison_negative_rate=True)
        self.assertTrue(neg_rate["passed"])
        self.assertIn("Correctly caught", neg_rate["error"])

        baseline = evaluate_falsification_negative_control()
        self.assertTrue(baseline["passed"])
        self.assertIsNone(baseline["error"])


if __name__ == "__main__":
    unittest.main()
