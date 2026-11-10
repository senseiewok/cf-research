#!/usr/bin/env python3
"""Unit tests for Round 174: Eyring Transition State & NBD1 Thermal Unfolding Invariant Bridge (T-0842)."""
from __future__ import annotations

import unittest

from tools.discovery.round174_eyring_nbd1_thermostability_bridge import (
    FIVE_LAYER_LADDER,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    complex_step_eyring_derivative,
    compute_er_maturation_efficiency,
    compute_eyring_unfolding_rate,
    compute_structural_half_life,
    evaluate_falsification_negative_control,
    evaluate_round174_eyring_thermostability,
)


class Round174EyringNBD1Test(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round174_eyring_thermostability(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Failed for {key}")
            self.assertGreater(res.melting_temperature_c, 20.0)
            self.assertLess(res.melting_temperature_c, 60.0)
            self.assertGreater(res.activation_free_energy_kj_mol, 50.0)
            self.assertLess(res.activation_free_energy_kj_mol, 150.0)
            self.assertGreater(res.unfolding_rate_s_minus_1, 0.0)
            self.assertGreater(res.structural_half_life_hours, 0.0)
            self.assertGreaterEqual(res.er_maturation_efficiency, 0.0)
            self.assertLessEqual(res.er_maturation_efficiency, 1.0)
            self.assertLessEqual(res.complex_step_relative_error, 1e-10)

    def test_wild_type_cftr_nbd1(self) -> None:
        res = evaluate_round174_eyring_thermostability("wild_type_cftr_nbd1")
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res.melting_temperature_c, 42.0)
        self.assertGreaterEqual(res.activation_free_energy_kj_mol, 105.0)
        self.assertGreaterEqual(res.structural_half_life_hours, 20.0)
        self.assertGreaterEqual(res.er_maturation_efficiency, 0.90)

    def test_trikafta_synergistic_rescue(self) -> None:
        res = evaluate_round174_eyring_thermostability("trikafta_synergistic_rescue")
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res.melting_temperature_c, 41.5)
        self.assertGreaterEqual(res.activation_free_energy_kj_mol, 104.0)
        self.assertGreaterEqual(res.structural_half_life_hours, 15.0)
        self.assertGreaterEqual(res.er_maturation_efficiency, 0.85)

    def test_uncorrected_delta_f508_cftr(self) -> None:
        res = evaluate_round174_eyring_thermostability("uncorrected_delta_f508_cftr")
        self.assertFalse(res.crosses_kahn_kalai_threshold)
        self.assertLess(res.melting_temperature_c, 37.0)  # Melts below body temp!
        self.assertLess(res.activation_free_energy_kj_mol, 100.0)
        self.assertLess(res.structural_half_life_hours, 2.0)
        self.assertLess(res.er_maturation_efficiency, 0.40)

    def test_single_corrector_tezacaftor(self) -> None:
        res = evaluate_round174_eyring_thermostability("single_corrector_tezacaftor")
        self.assertFalse(res.crosses_kahn_kalai_threshold)
        self.assertLess(res.melting_temperature_c, 37.0)
        self.assertLess(res.structural_half_life_hours, 6.0)
        self.assertLess(res.er_maturation_efficiency, 0.75)

    def test_extreme_thermolabile_i507del(self) -> None:
        res = evaluate_round174_eyring_thermostability("extreme_thermolabile_i507del")
        self.assertFalse(res.crosses_kahn_kalai_threshold)
        self.assertLess(res.melting_temperature_c, 33.0)
        self.assertLess(res.structural_half_life_hours, 0.5)
        self.assertLess(res.er_maturation_efficiency, 0.15)

    def test_complex_step_eyring_derivative(self) -> None:
        cs_deriv, analytical, rel_err = complex_step_eyring_derivative(
            delta_g_act_j_mol=106500.0,
            temperature_k=310.15,
        )
        self.assertLessEqual(rel_err, 1e-10)
        self.assertLess(cs_deriv, 0.0)  # Unfolding rate decreases as activation free energy increases

    def test_negative_control_falsification(self) -> None:
        neg1 = evaluate_falsification_negative_control(poison_negative_temperature=True)
        self.assertFalse(neg1["passed"])
        self.assertIn("poisoned", neg1["error"])

        neg2 = evaluate_falsification_negative_control(poison_negative_activation_energy=True)
        self.assertFalse(neg2["passed"])

        neg3 = evaluate_falsification_negative_control(poison_negative_unfolding_rate=True)
        self.assertFalse(neg3["passed"])

        neg4 = evaluate_falsification_negative_control(poison_negative_er_transit=True)
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
