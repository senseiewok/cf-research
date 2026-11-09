#!/usr/bin/env python3
"""Unit tests for Round 165: Kuramoto-Taylor Ciliary Metachronal Synchronization & Viscous Coupling Invariant Bridge (T-0830)."""
from __future__ import annotations

import unittest

from tools.discovery.round165_ciliary_metachronal_wave_bridge import (
    FIVE_LAYER_LADDER,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    complex_step_transport_derivative,
    compute_hydrodynamic_coupling_constant,
    compute_kuramoto_order_parameter,
    compute_metachronal_wave_velocity,
    compute_mucociliary_transport_velocity,
    evaluate_falsification_negative_control,
    evaluate_round165_ciliary_metachronal_wave,
)


class Round165CiliaryMetachronalWaveTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round165_ciliary_metachronal_wave(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Failed for {key}")
            self.assertGreaterEqual(res.kuramoto_order_parameter, 0.0)
            self.assertLessEqual(res.kuramoto_order_parameter, 1.0)
            self.assertGreaterEqual(res.mucociliary_clearance_velocity_um_s, 0.0)
            self.assertLessEqual(res.complex_step_relative_error, 1e-10)

    def test_healthy_airway_metachronal_wave(self) -> None:
        res = evaluate_round165_ciliary_metachronal_wave("healthy_airway_metachronal_wave")
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res.kuramoto_order_parameter, 0.75)
        self.assertGreaterEqual(res.mucociliary_clearance_velocity_um_s, 40.0)
        self.assertGreaterEqual(res.metachronal_wave_velocity_um_s, 200.0)
        self.assertGreater(res.periciliary_shear_stress_pa, 0.05)

    def test_trikafta_rescued_f508del_cilia(self) -> None:
        res = evaluate_round165_ciliary_metachronal_wave("trikafta_rescued_f508del_cilia")
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res.kuramoto_order_parameter, 0.70)
        self.assertGreaterEqual(res.mucociliary_clearance_velocity_um_s, 30.0)

    def test_kahn_kalai_transport_boundary(self) -> None:
        res = evaluate_round165_ciliary_metachronal_wave("kahn_kalai_transport_boundary")
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res.mucociliary_clearance_velocity_um_s, 10.0)

    def test_untreated_cf_pcl_collapse_stasis(self) -> None:
        res = evaluate_round165_ciliary_metachronal_wave("untreated_cf_pcl_collapse_stasis")
        self.assertFalse(res.crosses_kahn_kalai_threshold)
        self.assertLess(res.mucociliary_clearance_velocity_um_s, 5.0)
        self.assertLess(res.kuramoto_order_parameter, 0.30)

    def test_hyperdense_synchronized_ciliary_carpet(self) -> None:
        res = evaluate_round165_ciliary_metachronal_wave("hyperdense_synchronized_ciliary_carpet")
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res.kuramoto_order_parameter, 0.90)
        self.assertGreaterEqual(res.mucociliary_clearance_velocity_um_s, 50.0)

    def test_kuramoto_order_parameter_transition(self) -> None:
        dispersion = 1.20
        # Below critical coupling: R must be strictly 0
        r_sub = compute_kuramoto_order_parameter(k_visc=1.0, frequency_dispersion_hz=dispersion)
        self.assertEqual(r_sub, 0.0)

        # Above critical coupling: R must monotonically increase towards 1
        r_mid = compute_kuramoto_order_parameter(k_visc=15.0, frequency_dispersion_hz=dispersion)
        r_high = compute_kuramoto_order_parameter(k_visc=50.0, frequency_dispersion_hz=dispersion)
        self.assertGreater(r_mid, 0.0)
        self.assertGreater(r_high, r_mid)
        self.assertLess(r_high, 1.0)

    def test_complex_step_transport_sensitivity(self) -> None:
        cs_deriv, analytical, rel_err = complex_step_transport_derivative(
            order_parameter_r=0.80,
            frequency_hz=13.50,
            pcl_height_um=7.00,
            mucus_viscosity_mpa_s=15.0,
        )
        self.assertLessEqual(rel_err, 1e-10)
        self.assertGreater(cs_deriv, 0.0)

    def test_negative_control_falsification(self) -> None:
        neg1 = evaluate_falsification_negative_control(poison_negative_frequency=True)
        self.assertFalse(neg1["passed"])
        self.assertIn("poisoned", neg1["error"])

        neg2 = evaluate_falsification_negative_control(poison_negative_pcl_height=True)
        self.assertFalse(neg2["passed"])

        neg3 = evaluate_falsification_negative_control(poison_negative_density=True)
        self.assertFalse(neg3["passed"])

        neg4 = evaluate_falsification_negative_control(poison_unphysical_order_parameter=True)
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
