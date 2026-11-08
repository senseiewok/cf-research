#!/usr/bin/env python3
"""Unit tests for Round 121: Landau-Ginzburg & Privalov-Freire NBD1 Calorimetric Stabilization Bridge."""
from __future__ import annotations

import unittest

from tools.discovery.round121_landau_ginzburg_sionna_nbd1_calorimetric_stabilization_bridge import (
    FIVE_LAYER_LADDER,
    PRIMARY_CITATIONS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    complex_step_dtm_sensitivity,
    evaluate_nbd1_calorimetric_bridge,
    verify_buckingham_pi_calorimetry_dimensions,
)


class Round121LandauGinzburgSionnaNbd1Test(unittest.TestCase):
    def test_buckingham_pi_exact_rational_dimensions(self) -> None:
        dim = verify_buckingham_pi_calorimetry_dimensions()
        self.assertTrue(dim["passed"])
        self.assertTrue(dim["ligand_ratio_dimensionless"])
        self.assertEqual(dim["lhs_si_exponents"], [0, 0, 0, 0, 1, 0, 0])
        self.assertEqual(dim["rhs_si_exponents"], [0, 0, 0, 0, 1, 0, 0])

    def test_complex_step_machine_precision_derivative(self) -> None:
        cs = complex_step_dtm_sensitivity(ligand_ratio=25.0)
        self.assertLess(cs["relative_error"], 1e-12)
        self.assertGreater(cs["complex_step_deriv"], 0.0)

    def test_uncorrected_f508del_melts_near_body_temperature(self) -> None:
        uncorr = evaluate_nbd1_calorimetric_bridge("f508del_uncorrected_1ckz")
        self.assertAlmostEqual(uncorr.tm_effective_celsius, 41.2, places=2)
        self.assertEqual(uncorr.delta_tm_celsius, 0.0)
        self.assertFalse(uncorr.dual_lock_satisfied)
        self.assertLess(uncorr.band_c_maturation_pct, 15.0)

    def test_sion638_plus_glpg2737_satisfies_rabeh_mendoza_dual_lock(self) -> None:
        dual = evaluate_nbd1_calorimetric_bridge("sion638_glpg2737_dual_lock")
        eti_only = evaluate_nbd1_calorimetric_bridge("f508del_eti_interface_only_8eiq")
        nbd1_only = evaluate_nbd1_calorimetric_bridge("sion638_nbd1_monotherapy")
        self.assertTrue(dual.dual_lock_satisfied)
        self.assertGreaterEqual(dual.tm_effective_celsius, 48.0)
        self.assertGreater(dual.band_c_maturation_pct, eti_only.band_c_maturation_pct)
        self.assertGreater(dual.band_c_maturation_pct, nbd1_only.band_c_maturation_pct)

    def test_monotonic_melting_shift_with_ligand_concentration(self) -> None:
        low = evaluate_nbd1_calorimetric_bridge("sion638_glpg2737_dual_lock", custom_ligand_kd_ratio=2.0)
        mid = evaluate_nbd1_calorimetric_bridge("sion638_glpg2737_dual_lock", custom_ligand_kd_ratio=15.0)
        high = evaluate_nbd1_calorimetric_bridge("sion638_glpg2737_dual_lock", custom_ligand_kd_ratio=45.0)
        self.assertLess(low.delta_tm_celsius, mid.delta_tm_celsius)
        self.assertLess(mid.delta_tm_celsius, high.delta_tm_celsius)

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
