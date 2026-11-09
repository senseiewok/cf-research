#!/usr/bin/env python3
"""Unit tests for Round 167: Smoluchowski-Coagulation & Shear-Induced Aerosol Droplet Coalescence Invariant Bridge (T-0832)."""
from __future__ import annotations

import unittest

from tools.discovery.round167_smoluchowski_aerosol_coalescence_bridge import (
    FIVE_LAYER_LADDER,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    complex_step_fpf_derivative,
    compute_fine_particle_fraction,
    compute_lnp_cargo_retention,
    compute_smoluchowski_shear_kernel,
    compute_tracheobronchial_deposition,
    evaluate_falsification_negative_control,
    evaluate_round167_smoluchowski_aerosol_coalescence,
)


class Round167AerosolCoalescenceTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round167_smoluchowski_aerosol_coalescence(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Failed for {key}")
            self.assertGreaterEqual(res.fine_particle_fraction_pct, 0.0)
            self.assertLessEqual(res.fine_particle_fraction_pct, 100.0)
            self.assertGreaterEqual(res.post_nebulization_retention_pct, 0.0)
            self.assertLessEqual(res.post_nebulization_retention_pct, 100.0)
            self.assertLessEqual(res.complex_step_relative_error, 1e-10)

    def test_optimal_vibrating_mesh_lnp_aerosol(self) -> None:
        res = evaluate_round167_smoluchowski_aerosol_coalescence("optimal_vibrating_mesh_lnp_aerosol")
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res.post_nebulization_retention_pct, 90.0)
        self.assertGreaterEqual(res.fine_particle_fraction_pct, 75.0)
        self.assertGreaterEqual(res.effective_delivery_efficiency, 0.20)
        self.assertGreater(res.capillary_number, 0.0)

    def test_shear_stabilized_polysarcosine_lnp(self) -> None:
        res = evaluate_round167_smoluchowski_aerosol_coalescence("shear_stabilized_polysarcosine_lnp")
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res.post_nebulization_retention_pct, 92.0)
        self.assertGreaterEqual(res.fine_particle_fraction_pct, 80.0)

    def test_kahn_kalai_aerosol_delivery_boundary(self) -> None:
        res = evaluate_round167_smoluchowski_aerosol_coalescence("kahn_kalai_aerosol_delivery_boundary")
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res.post_nebulization_retention_pct, 85.0)
        self.assertGreaterEqual(res.effective_delivery_efficiency, 0.15)

    def test_unshielded_high_shear_jet_nebulizer_rupture(self) -> None:
        res = evaluate_round167_smoluchowski_aerosol_coalescence("unshielded_high_shear_jet_nebulizer_rupture")
        self.assertFalse(res.crosses_kahn_kalai_threshold)
        self.assertLess(res.post_nebulization_retention_pct, 70.0)
        self.assertLess(res.fine_particle_fraction_pct, 50.0)
        self.assertLess(res.effective_delivery_efficiency, 0.10)

    def test_submicron_monodisperse_acoustic_aerosol(self) -> None:
        res = evaluate_round167_smoluchowski_aerosol_coalescence("submicron_monodisperse_acoustic_aerosol")
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res.fine_particle_fraction_pct, 90.0)
        self.assertGreaterEqual(res.post_nebulization_retention_pct, 93.0)

    def test_complex_step_fpf_derivative(self) -> None:
        cs_deriv, analytical, rel_err = complex_step_fpf_derivative(
            mmad_um=3.00,
            gsd=1.65,
        )
        self.assertLessEqual(rel_err, 1e-10)
        self.assertLess(cs_deriv, 0.0)  # FPF decreases with increasing MMAD

    def test_negative_control_falsification(self) -> None:
        neg1 = evaluate_falsification_negative_control(poison_negative_mmad=True)
        self.assertFalse(neg1["passed"])
        self.assertIn("poisoned", neg1["error"])

        neg2 = evaluate_falsification_negative_control(poison_invalid_gsd=True)
        self.assertFalse(neg2["passed"])

        neg3 = evaluate_falsification_negative_control(poison_negative_shear=True)
        self.assertFalse(neg3["passed"])

        neg4 = evaluate_falsification_negative_control(poison_invalid_encapsulation=True)
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
