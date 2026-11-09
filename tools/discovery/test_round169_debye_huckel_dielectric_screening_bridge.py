#!/usr/bin/env python3
"""Unit tests for Round 169: Debye-Hückel-Bjerrum Ion Screening & Electrostatic Polyelectrolyte Repulsion Invariant Bridge (T-0834)."""
from __future__ import annotations

import unittest

from tools.discovery.round169_debye_huckel_dielectric_screening_bridge import (
    FIVE_LAYER_LADDER,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    complex_step_debye_derivative,
    compute_bjerrum_length,
    compute_debye_length,
    compute_donnan_osmotic_pressure,
    compute_manning_condensation,
    compute_mucus_complex_modulus,
    evaluate_falsification_negative_control,
    evaluate_round169_debye_huckel_screening,
)


class Round169DebyeHuckelTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round169_debye_huckel_screening(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Failed for {key}")
            self.assertGreater(res.debye_length_nm, 0.0)
            self.assertGreater(res.bjerrum_length_nm, 0.0)
            self.assertGreaterEqual(res.counterion_condensation_fraction, 0.0)
            self.assertLess(res.counterion_condensation_fraction, 1.0)
            self.assertGreater(res.complex_shear_modulus_pa, 0.0)
            self.assertGreater(res.softening_ratio, 0.0)
            self.assertGreater(res.clearance_velocity_recovery_ratio, 0.0)
            self.assertLessEqual(res.complex_step_relative_error, 1e-10)

    def test_hypertonic_saline_7pct_softening(self) -> None:
        res = evaluate_round169_debye_huckel_screening("hypertonic_saline_7pct_softening")
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertLessEqual(res.softening_ratio, 0.40)
        self.assertLessEqual(res.debye_length_nm, 0.35)
        self.assertGreaterEqual(res.clearance_velocity_recovery_ratio, 2.50)
        self.assertLess(res.complex_shear_modulus_pa, 2.50)

    def test_physiological_cf_mucus_baseline(self) -> None:
        res = evaluate_round169_debye_huckel_screening("physiological_cf_mucus_baseline")
        self.assertFalse(res.crosses_kahn_kalai_threshold)
        self.assertAlmostEqual(res.softening_ratio, 1.00, places=2)
        self.assertGreater(res.debye_length_nm, 0.70)
        self.assertAlmostEqual(res.clearance_velocity_recovery_ratio, 1.00, places=2)

    def test_hypertonic_saline_3pct_intermediate(self) -> None:
        res = evaluate_round169_debye_huckel_screening("hypertonic_saline_3pct_intermediate")
        self.assertFalse(res.crosses_kahn_kalai_threshold)
        self.assertGreater(res.softening_ratio, 0.40)
        self.assertGreater(res.debye_length_nm, 0.35)
        self.assertGreater(res.clearance_velocity_recovery_ratio, 1.50)

    def test_bicarbonate_hypertonic_synergy(self) -> None:
        res = evaluate_round169_debye_huckel_screening("bicarbonate_hypertonic_synergy")
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertLessEqual(res.softening_ratio, 0.25)
        self.assertGreaterEqual(res.clearance_velocity_recovery_ratio, 3.00)

    def test_hypotonic_dilution_pathology(self) -> None:
        res = evaluate_round169_debye_huckel_screening("hypotonic_dilution_pathology")
        self.assertFalse(res.crosses_kahn_kalai_threshold)
        self.assertGreater(res.softening_ratio, 1.20)
        self.assertGreater(res.debye_length_nm, 1.00)
        self.assertLess(res.clearance_velocity_recovery_ratio, 1.00)

    def test_complex_step_debye_derivative(self) -> None:
        cs_deriv, analytical, rel_err = complex_step_debye_derivative(
            ionic_strength_m=1.200,
        )
        self.assertLessEqual(rel_err, 1e-10)
        self.assertLess(cs_deriv, 0.0)  # Debye length strictly decreases with increasing ionic strength

    def test_negative_control_falsification(self) -> None:
        neg1 = evaluate_falsification_negative_control(poison_negative_ionic_strength=True)
        self.assertFalse(neg1["passed"])
        self.assertIn("poisoned", neg1["error"])

        neg2 = evaluate_falsification_negative_control(poison_negative_temperature=True)
        self.assertFalse(neg2["passed"])

        neg3 = evaluate_falsification_negative_control(poison_invalid_charge_spacing=True)
        self.assertFalse(neg3["passed"])

        neg4 = evaluate_falsification_negative_control(poison_negative_fixed_charge=True)
        self.assertFalse(neg4["passed"])

        clean = evaluate_falsification_negative_control()
        self.assertTrue(clean["passed"])

    def test_epistemic_hygiene_and_ladder_structure(self) -> None:
        self.assertIn("deductive_theorem", TWO_COLUMN_EPISTEMIC_HYGIENE)
        self.assertIn("biophysical_hypothesis", TWO_COLUMN_EPISTEMIC_HYGIENE)
        self.assertGreaterEqual(len(FIVE_LAYER_LADDER), 5)
        self.assertGreaterEqual(len(PRIMARY_CITATIONS), 4)


if __name__ == "__main__":
    unittest.main()
