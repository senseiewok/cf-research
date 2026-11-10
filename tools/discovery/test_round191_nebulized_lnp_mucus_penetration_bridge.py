#!/usr/bin/env python3
"""Unit tests for Round 191: Nebulized Mucus-Penetrating LNP Surface Charge & PEG Density Bridge (T-0872)."""
from __future__ import annotations

import unittest

from tools.discovery.round191_nebulized_lnp_mucus_penetration_bridge import (
    FIVE_LAYER_LADDER,
    PARAMETER_PROVENANCE,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    complex_step_asl_derivative,
    compute_adhesive_binding_factor,
    compute_airway_deposition_fraction,
    compute_asl_height_um,
    compute_effective_diffusivity_ratio,
    compute_mcc_velocity_mm_min,
    compute_mucus_penetration_fraction,
    compute_steric_obstruction_fraction,
    compute_transfection_and_chloride_current,
    evaluate_falsification_negative_control,
    evaluate_round191_nebulized_lnp,
)


class Round191LNPMucusPenetrationTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round191_nebulized_lnp(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Invariants failed for {key}")
            self.assertGreaterEqual(res.steric_fraction, 0.0)
            self.assertLessEqual(res.steric_fraction, 1.0)
            self.assertGreaterEqual(res.adhesive_binding_factor, 0.0)
            self.assertGreaterEqual(res.diffusivity_ratio, 0.0)
            self.assertLessEqual(res.diffusivity_ratio, 1.0)
            self.assertGreaterEqual(res.airway_deposition_fraction, 0.0)
            self.assertLessEqual(res.airway_deposition_fraction, 1.0)
            self.assertGreaterEqual(res.mucus_penetration_fraction, 0.0)
            self.assertLessEqual(res.mucus_penetration_fraction, 1.0)
            self.assertGreaterEqual(res.transfection_fraction, 0.0)
            self.assertLessEqual(res.transfection_fraction, 1.0)
            self.assertGreaterEqual(res.transepithelial_current_pct_wt, 0.0)
            self.assertLessEqual(res.transepithelial_current_pct_wt, 100.0)
            self.assertGreaterEqual(res.airway_surface_liquid_height_um, 3.5)
            self.assertLessEqual(res.airway_surface_liquid_height_um, 8.5)
            self.assertGreaterEqual(res.mucociliary_transport_velocity_mm_min, 0.0)
            self.assertLessEqual(res.complex_step_relative_error, 1e-6)

    def test_curative_presets_cross_kahn_kalai(self) -> None:
        res1 = evaluate_round191_nebulized_lnp("nebulized_dense_peg_neutral_lnp")
        self.assertTrue(res1.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res1.mucus_penetration_fraction, 0.30)
        self.assertGreaterEqual(res1.transfection_fraction, 0.15)
        self.assertGreaterEqual(res1.transepithelial_current_pct_wt, 12.0)
        self.assertGreaterEqual(res1.airway_surface_liquid_height_um, 7.0)
        self.assertGreaterEqual(res1.mucociliary_transport_velocity_mm_min, 1.0)

        res2 = evaluate_round191_nebulized_lnp("hyper_pegylated_high_dose_lnp")
        self.assertTrue(res2.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res2.mucus_penetration_fraction, 0.30)
        self.assertGreaterEqual(res2.transfection_fraction, 0.15)
        self.assertGreaterEqual(res2.transepithelial_current_pct_wt, 12.0)
        self.assertGreaterEqual(res2.airway_surface_liquid_height_um, 7.0)

    def test_sub_curative_presets_fail_kahn_kalai(self) -> None:
        # Cationic low-PEG entrapped
        res_cationic = evaluate_round191_nebulized_lnp("cationic_low_peg_entrapped_lnp")
        self.assertFalse(res_cationic.crosses_kahn_kalai_threshold)
        self.assertLess(res_cationic.mucus_penetration_fraction, 0.10)
        self.assertLess(res_cationic.transepithelial_current_pct_wt, 5.0)

        # Large droplet oropharyngeal impaction
        res_aero = evaluate_round191_nebulized_lnp("large_droplet_oropharyngeal_impaction")
        self.assertFalse(res_aero.crosses_kahn_kalai_threshold)
        self.assertLess(res_aero.transfection_fraction, 0.15)
        self.assertLess(res_aero.transepithelial_current_pct_wt, 12.0)

        # Severe mucus hyperconcentrated barrier
        res_mucus = evaluate_round191_nebulized_lnp("severe_mucus_hyperconcentrated_barrier")
        self.assertFalse(res_mucus.crosses_kahn_kalai_threshold)
        self.assertLess(res_mucus.diffusivity_ratio, 0.05)
        self.assertLess(res_mucus.mucus_penetration_fraction, 0.05)

    def test_complex_step_derivative_accuracy(self) -> None:
        cs_deriv, fd_deriv, rel_err = complex_step_asl_derivative(current_pct_wt=25.0)
        self.assertGreater(cs_deriv, 0.0)
        self.assertGreater(fd_deriv, 0.0)
        self.assertLessEqual(rel_err, 1e-6)

    def test_poisoned_negative_controls_fail_deterministically(self) -> None:
        res_peg = evaluate_falsification_negative_control(poison_negative_peg=True)
        self.assertTrue(res_peg["passed"])
        self.assertIn("Correctly caught", res_peg["error"])

        res_aero = evaluate_falsification_negative_control(poison_negative_aero=True)
        self.assertTrue(res_aero["passed"])
        self.assertIn("Correctly caught", res_aero["error"])

        res_pore = evaluate_falsification_negative_control(poison_negative_pore=True)
        self.assertTrue(res_pore["passed"])
        self.assertIn("Correctly caught", res_pore["error"])

        clean = evaluate_falsification_negative_control()
        self.assertTrue(clean["passed"])
        self.assertIsNone(clean["error"])

    def test_monotonic_peg_density_response(self) -> None:
        # Increasing PEG brush density monotonically reduces adhesive binding
        k_low = compute_adhesive_binding_factor(peg_density_chains_nm2=0.05, zeta_potential_mv=10.0)
        k_mid = compute_adhesive_binding_factor(peg_density_chains_nm2=0.25, zeta_potential_mv=10.0)
        k_high = compute_adhesive_binding_factor(peg_density_chains_nm2=0.50, zeta_potential_mv=10.0)
        self.assertGreater(k_low, k_mid)
        self.assertGreater(k_mid, k_high)

    def test_provenance_and_two_column_epistemic_hygiene_metadata(self) -> None:
        self.assertGreaterEqual(len(PRIMARY_CITATIONS), 4)
        self.assertGreaterEqual(len(PARAMETER_PROVENANCE), 5)
        self.assertIn("deductive_theorem", TWO_COLUMN_EPISTEMIC_HYGIENE)
        self.assertIn("biophysical_hypothesis", TWO_COLUMN_EPISTEMIC_HYGIENE)
        self.assertIn("statement", TWO_COLUMN_EPISTEMIC_HYGIENE["biophysical_hypothesis"])
        self.assertIn("falsification_condition", TWO_COLUMN_EPISTEMIC_HYGIENE["biophysical_hypothesis"])
        for layer in ["cool_human_fact", "everyday_analogy", "engineering_lens", "formal_theorem", "falsifiable_code_and_toy"]:
            self.assertIn(layer, FIVE_LAYER_LADDER)


if __name__ == "__main__":
    unittest.main()
