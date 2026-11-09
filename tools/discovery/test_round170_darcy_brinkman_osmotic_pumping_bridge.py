#!/usr/bin/env python3
"""Unit tests for Round 170: Darcy-Brinkman Porous Media Hydrodynamics & Periciliary Osmotic Pumping Invariant Bridge (T-0835)."""
from __future__ import annotations

import unittest

from tools.discovery.round170_darcy_brinkman_osmotic_pumping_bridge import (
    FIVE_LAYER_LADDER,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    complex_step_clearance_derivative,
    compute_ciliary_bed_permeability,
    compute_mucociliary_clearance_velocity,
    compute_pcl_equilibrium_height,
    evaluate_falsification_negative_control,
    evaluate_round170_darcy_brinkman_osmotic_pumping,
)


class Round170DarcyBrinkmanTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round170_darcy_brinkman_osmotic_pumping(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Failed for {key}")
            self.assertGreater(res.permeability_m2, 0.0)
            self.assertGreater(res.hydraulic_conductivity_m_pa_s, 0.0)
            self.assertGreater(res.pcl_height_um, 0.0)
            self.assertGreaterEqual(res.mucociliary_clearance_velocity_um_s, 0.0)
            self.assertLessEqual(res.mucociliary_clearance_velocity_um_s, 65.0)
            self.assertLessEqual(res.complex_step_relative_error, 1e-10)

    def test_healthy_airway_pcl_homeostasis(self) -> None:
        res = evaluate_round170_darcy_brinkman_osmotic_pumping("healthy_airway_pcl_homeostasis")
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res.pcl_height_um, 6.00)
        self.assertGreaterEqual(res.mucociliary_clearance_velocity_um_s, 45.0)
        self.assertGreaterEqual(res.permeability_m2, 1.0e-14)
        self.assertLessEqual(res.permeability_m2, 1.0e-13)

    def test_cf_dehydrated_pcl_collapse(self) -> None:
        res = evaluate_round170_darcy_brinkman_osmotic_pumping("cf_dehydrated_pcl_collapse")
        self.assertFalse(res.crosses_kahn_kalai_threshold)
        self.assertLess(res.pcl_height_um, 3.50)
        self.assertLess(res.mucociliary_clearance_velocity_um_s, 2.0)

    def test_hypertonic_saline_osmotic_rehydration(self) -> None:
        res = evaluate_round170_darcy_brinkman_osmotic_pumping("hypertonic_saline_osmotic_rehydration")
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res.pcl_height_um, 6.50)
        self.assertGreaterEqual(res.mucociliary_clearance_velocity_um_s, 50.0)
        self.assertGreater(res.secretion_velocity_um_s, 0.0)

    def test_cftr_potentiator_steady_pumping(self) -> None:
        res = evaluate_round170_darcy_brinkman_osmotic_pumping("cftr_potentiator_steady_pumping")
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res.pcl_height_um, 6.00)
        self.assertGreaterEqual(res.mucociliary_clearance_velocity_um_s, 35.0)

    def test_severe_hyposecretion_enac_hyperactive(self) -> None:
        res = evaluate_round170_darcy_brinkman_osmotic_pumping("severe_hyposecretion_enac_hyperactive")
        self.assertFalse(res.crosses_kahn_kalai_threshold)
        self.assertLess(res.pcl_height_um, 2.50)
        self.assertLess(res.mucociliary_clearance_velocity_um_s, 0.5)

    def test_complex_step_clearance_derivative(self) -> None:
        cs_deriv, analytical, rel_err = complex_step_clearance_derivative(
            pcl_height_um=6.00,
        )
        self.assertLessEqual(rel_err, 1e-10)
        self.assertGreater(cs_deriv, 0.0)  # Clearance velocity strictly increases with PCL height

    def test_negative_control_falsification(self) -> None:
        neg1 = evaluate_falsification_negative_control(poison_negative_ciliary_diameter=True)
        self.assertFalse(neg1["passed"])
        self.assertIn("poisoned", neg1["error"])

        neg2 = evaluate_falsification_negative_control(poison_invalid_volume_fraction=True)
        self.assertFalse(neg2["passed"])

        neg3 = evaluate_falsification_negative_control(poison_negative_pcl_height=True)
        self.assertFalse(neg3["passed"])

        neg4 = evaluate_falsification_negative_control(poison_negative_viscosity=True)
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
