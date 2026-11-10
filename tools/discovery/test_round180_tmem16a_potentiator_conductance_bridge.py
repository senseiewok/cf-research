#!/usr/bin/env python3
"""Unit tests for Round 180: TMEM16A Anion Potentiator Conductance Bridge (T-0854)."""
from __future__ import annotations

import unittest

from tools.discovery.round180_tmem16a_potentiator_conductance_bridge import (
    FIVE_LAYER_LADDER,
    PARAMETER_PROVENANCE,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    complex_step_ca_derivative,
    compute_apparent_ec50,
    compute_asl_height_and_mcc_velocity,
    compute_epithelial_conductance_and_isc,
    compute_tmem16a_open_probability,
    evaluate_falsification_negative_control,
    evaluate_round180_tmem16a_potentiator_conductance,
)


class Round180Tmem16aPotentiatorConductanceTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round180_tmem16a_potentiator_conductance(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Invariants failed for {key}")
            self.assertGreaterEqual(res.open_probability, 0.0)
            self.assertLessEqual(res.open_probability, 1.0)
            self.assertGreaterEqual(res.conductance_ms_cm2, 0.0)
            self.assertGreaterEqual(res.equivalent_isc_pct_wt, 0.0)
            self.assertGreaterEqual(res.asl_height_um, 0.0)
            self.assertLessEqual(res.asl_height_um, 15.0)
            self.assertGreaterEqual(res.mcc_velocity_mm_min, 0.0)
            self.assertLessEqual(res.complex_step_relative_error, 1e-6)

    def test_optimal_potentiator_monotherapy(self) -> None:
        res = evaluate_round180_tmem16a_potentiator_conductance("potentiator_monotherapy_optimal")
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res.asl_height_um, 7.0)
        self.assertGreaterEqual(res.mcc_velocity_mm_min, 1.0)
        self.assertGreaterEqual(res.equivalent_isc_pct_wt, 15.0)

    def test_high_affinity_regime(self) -> None:
        res = evaluate_round180_tmem16a_potentiator_conductance("potentiator_high_affinity_regime")
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res.asl_height_um, 7.0)
        self.assertGreaterEqual(res.mcc_velocity_mm_min, 1.0)
        self.assertGreaterEqual(res.equivalent_isc_pct_wt, 30.0)

    def test_subtherapeutic_and_untreated_fail_kahn_kalai(self) -> None:
        # Subtherapeutic low dose
        low = evaluate_round180_tmem16a_potentiator_conductance("subtherapeutic_potentiator_low_dose")
        self.assertFalse(low.crosses_kahn_kalai_threshold)
        self.assertLess(low.asl_height_um, 7.0)
        self.assertLess(low.mcc_velocity_mm_min, 1.0)

        # Untreated baseline
        untreated = evaluate_round180_tmem16a_potentiator_conductance("untreated_cftr_null_baseline")
        self.assertFalse(untreated.crosses_kahn_kalai_threshold)
        self.assertLess(untreated.asl_height_um, 7.0)
        self.assertLess(untreated.mcc_velocity_mm_min, 0.1)

    def test_calcium_free_potentiated_state(self) -> None:
        ca_free = evaluate_round180_tmem16a_potentiator_conductance("calcium_free_potentiated_state")
        self.assertFalse(ca_free.crosses_kahn_kalai_threshold)
        self.assertLess(ca_free.asl_height_um, 7.0)
        self.assertGreater(ca_free.open_probability, 0.05)  # Confirms ligand-independent leak opening

    def test_complex_step_ca_derivative(self) -> None:
        cs_deriv, fd_deriv, rel_err = complex_step_ca_derivative(
            ca_i_um=0.12,
            c_pot_um=1.0,
        )
        self.assertLessEqual(rel_err, 1e-6)
        self.assertGreater(cs_deriv, 0.0)

    def test_negative_control_falsification(self) -> None:
        neg_c = evaluate_falsification_negative_control(poison_negative_c_pot=True)
        self.assertTrue(neg_c["passed"])
        self.assertIn("Correctly caught", neg_c["error"])

        neg_ca = evaluate_falsification_negative_control(poison_negative_ca_i=True)
        self.assertTrue(neg_ca["passed"])
        self.assertIn("Correctly caught", neg_ca["error"])

        neg_g = evaluate_falsification_negative_control(poison_negative_gamma=True)
        self.assertTrue(neg_g["passed"])
        self.assertIn("Correctly caught", neg_g["error"])

        neg_po = evaluate_falsification_negative_control(poison_unphysical_po=True)
        self.assertTrue(neg_po["passed"])
        self.assertIn("Correctly caught", neg_po["error"])

        baseline = evaluate_falsification_negative_control()
        self.assertTrue(baseline["passed"])
        self.assertIsNone(baseline["error"])

    def test_provenance_and_citations(self) -> None:
        self.assertGreaterEqual(len(PRIMARY_CITATIONS), 3)
        self.assertGreaterEqual(len(PARAMETER_PROVENANCE), 5)
        self.assertIn("deductive_theorem", TWO_COLUMN_EPISTEMIC_HYGIENE)
        self.assertIn("biophysical_hypothesis", TWO_COLUMN_EPISTEMIC_HYGIENE)
        self.assertEqual(len(FIVE_LAYER_LADDER), 5)


if __name__ == "__main__":
    unittest.main()
