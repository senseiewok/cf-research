#!/usr/bin/env python3
"""Unit tests for Round 148: Poisson-Boltzmann Nonlinear Electrostatic Salt Screening Invariant Bridge (T-0799)."""
from __future__ import annotations

import unittest

from tools.discovery.round148_poisson_boltzmann_glycocalyx_screening_bridge import (
    FIVE_LAYER_LADDER,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    evaluate_falsification_negative_control,
    evaluate_round148_poisson_boltzmann_glycocalyx_screening,
)


class Round148PoissonBoltzmannScreeningTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round148_poisson_boltzmann_glycocalyx_screening(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Failed for {key}")
            self.assertGreater(res.debye_length_nm, 0.0)
            self.assertLess(res.donnan_potential_mv, 0.0)
            self.assertGreater(res.electrostatic_transmission_factor, 0.0)
            self.assertLessEqual(res.electrostatic_transmission_factor, 1.0)
            self.assertGreater(res.effective_conductance_ps, 0.0)
            self.assertGreaterEqual(res.epithelial_short_circuit_current_ua_cm2, 0.0)

    def test_wild_type_glycocalyx_screening(self) -> None:
        wt = evaluate_round148_poisson_boltzmann_glycocalyx_screening("wild_type_glycocalyx_screening")
        self.assertTrue(wt.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(wt.effective_conductance_ps, 8.0)
        self.assertGreaterEqual(wt.epithelial_short_circuit_current_ua_cm2, 30.0)
        self.assertAlmostEqual(wt.fractional_wild_type_current, 0.97, delta=0.05)
        self.assertAlmostEqual(wt.debye_length_nm, 0.83, delta=0.02)
        self.assertGreaterEqual(wt.electrostatic_transmission_factor, 0.94)

    def test_cf_defects_and_screening_collapse(self) -> None:
        dehydrated = evaluate_round148_poisson_boltzmann_glycocalyx_screening("cf_dehydrated_mucus_hyperconcentration")
        hypotonic = evaluate_round148_poisson_boltzmann_glycocalyx_screening("cf_hypotonic_asl_screening_collapse")
        for mutant in (dehydrated, hypotonic):
            self.assertFalse(mutant.crosses_kahn_kalai_threshold)
            self.assertLessEqual(mutant.fractional_wild_type_current, 0.10)
        self.assertLessEqual(dehydrated.epithelial_short_circuit_current_ua_cm2, 1.0)
        self.assertLessEqual(hypotonic.epithelial_short_circuit_current_ua_cm2, 1.0)
        self.assertLess(dehydrated.electrostatic_transmission_factor, 0.80)
        self.assertLess(hypotonic.donnan_potential_mv, -10.0)

    def test_kahn_kalai_electrostatic_screening_boundary(self) -> None:
        kk = evaluate_round148_poisson_boltzmann_glycocalyx_screening("kahn_kalai_electrostatic_screening_boundary")
        self.assertTrue(kk.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(kk.effective_conductance_ps, 2.0)
        self.assertGreaterEqual(kk.epithelial_short_circuit_current_ua_cm2, 3.79)
        self.assertGreaterEqual(kk.fractional_wild_type_current, 0.118)

    def test_trikafta_and_hypertonic_saline_rescue(self) -> None:
        rescued = evaluate_round148_poisson_boltzmann_glycocalyx_screening("trikafta_and_hypertonic_saline_rescue")
        self.assertTrue(rescued.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(rescued.effective_conductance_ps, 8.0)
        self.assertGreaterEqual(rescued.epithelial_short_circuit_current_ua_cm2, 25.0)
        self.assertGreaterEqual(rescued.fractional_wild_type_current, 0.80)
        self.assertGreaterEqual(rescued.electrostatic_transmission_factor, 0.93)

    def test_negative_control_falsification(self) -> None:
        normal = evaluate_falsification_negative_control(poison_positive_mucin_charge=False)
        self.assertTrue(normal["passed"])
        self.assertEqual(normal["rule"], "AH-10")
        poisoned = evaluate_falsification_negative_control(poison_positive_mucin_charge=True)
        self.assertFalse(poisoned["passed"])

    def test_salt_screening_sweep(self) -> None:
        low_salt = evaluate_round148_poisson_boltzmann_glycocalyx_screening(
            "wild_type_glycocalyx_screening", custom_bulk_salt_mm=50.0
        )
        mid_salt = evaluate_round148_poisson_boltzmann_glycocalyx_screening(
            "wild_type_glycocalyx_screening", custom_bulk_salt_mm=140.0
        )
        high_salt = evaluate_round148_poisson_boltzmann_glycocalyx_screening(
            "wild_type_glycocalyx_screening", custom_bulk_salt_mm=300.0
        )
        self.assertGreater(low_salt.debye_length_nm, mid_salt.debye_length_nm)
        self.assertGreater(mid_salt.debye_length_nm, high_salt.debye_length_nm)
        self.assertLess(low_salt.electrostatic_transmission_factor, mid_salt.electrostatic_transmission_factor)
        self.assertLess(mid_salt.electrostatic_transmission_factor, high_salt.electrostatic_transmission_factor)
        self.assertLess(low_salt.donnan_potential_mv, mid_salt.donnan_potential_mv)
        self.assertLess(mid_salt.donnan_potential_mv, high_salt.donnan_potential_mv)

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
