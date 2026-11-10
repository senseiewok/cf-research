#!/usr/bin/env python3
"""Unit tests for Round 190: CFTR Post-Translational Ubiquitination Kinetics & CHIP/Rma1 E3 Ligase Bridge (T-0871)."""
from __future__ import annotations

import unittest

from tools.discovery.round190_cftr_ubiquitination_e3_ligase_bridge import (
    FIVE_LAYER_LADDER,
    PARAMETER_PROVENANCE,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    complex_step_asl_derivative,
    compute_asl_height_um,
    compute_er_escape_and_maturation_fractions,
    compute_exposed_misfolded_fraction,
    compute_mcc_velocity_mm_min,
    compute_transepithelial_current_pct_wt,
    compute_ubiquitination_rates_and_occupancy,
    evaluate_falsification_negative_control,
    evaluate_round190_cftr_ubiquitination,
)


class Round190UbiquitinationTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round190_cftr_ubiquitination(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Invariants failed for {key}")
            self.assertGreaterEqual(res.exposed_misfolded_fraction, 0.0)
            self.assertLessEqual(res.exposed_misfolded_fraction, 1.0)
            self.assertGreaterEqual(res.ubiquitination_occupancy, 0.0)
            self.assertLessEqual(res.ubiquitination_occupancy, 1.0)
            self.assertGreaterEqual(res.er_escape_fraction, 0.0)
            self.assertLessEqual(res.er_escape_fraction, 1.0)
            self.assertGreaterEqual(res.mature_band_c_fraction, 0.0)
            self.assertLessEqual(res.mature_band_c_fraction, 1.0)
            self.assertGreaterEqual(res.transepithelial_current_pct_wt, 0.0)
            self.assertLessEqual(res.transepithelial_current_pct_wt, 100.0)
            self.assertGreaterEqual(res.airway_surface_liquid_height_um, 3.5)
            self.assertLessEqual(res.airway_surface_liquid_height_um, 8.5)
            self.assertGreaterEqual(res.mucociliary_transport_velocity_mm_min, 0.0)
            self.assertLessEqual(res.complex_step_relative_error, 1e-6)

    def test_trikafta_plus_chip_inhibition_crosses_kahn_kalai(self) -> None:
        res = evaluate_round190_cftr_ubiquitination("trikafta_plus_chip_inhibition")
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res.er_escape_fraction, 0.35)
        self.assertGreaterEqual(res.mature_band_c_fraction, 0.25)
        self.assertGreaterEqual(res.transepithelial_current_pct_wt, 12.0)
        self.assertGreaterEqual(res.airway_surface_liquid_height_um, 7.0)
        self.assertGreaterEqual(res.mucociliary_transport_velocity_mm_min, 1.0)

    def test_trikafta_standard_rescue_crosses_kahn_kalai(self) -> None:
        res = evaluate_round190_cftr_ubiquitination("trikafta_standard_rescue")
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res.er_escape_fraction, 0.35)
        self.assertGreaterEqual(res.mature_band_c_fraction, 0.25)
        self.assertGreaterEqual(res.transepithelial_current_pct_wt, 12.0)
        self.assertGreaterEqual(res.airway_surface_liquid_height_um, 7.0)
        self.assertGreaterEqual(res.mucociliary_transport_velocity_mm_min, 1.0)

    def test_sub_curative_presets_fail_kahn_kalai(self) -> None:
        # CHIP overexpression barrier
        chip_hi = evaluate_round190_cftr_ubiquitination("chip_hyperactive_overexpression")
        self.assertFalse(chip_hi.crosses_kahn_kalai_threshold)
        self.assertLess(chip_hi.er_escape_fraction, 0.35)

        # Untreated F508del baseline
        untreated = evaluate_round190_cftr_ubiquitination("untreated_f508del_rapid_erad")
        self.assertFalse(untreated.crosses_kahn_kalai_threshold)
        self.assertLess(untreated.er_escape_fraction, 0.15)
        self.assertLess(untreated.mature_band_c_fraction, 0.10)
        self.assertLess(untreated.transepithelial_current_pct_wt, 10.0)
        self.assertLess(untreated.airway_surface_liquid_height_um, 7.0)

        # Isolated Rma1 knockdown
        rma1_kd = evaluate_round190_cftr_ubiquitination("isolated_rma1_knockdown_without_corrector")
        self.assertFalse(rma1_kd.crosses_kahn_kalai_threshold)
        self.assertLess(rma1_kd.er_escape_fraction, 0.15)
        self.assertLess(rma1_kd.mature_band_c_fraction, 0.10)

    def test_complex_step_derivative_accuracy(self) -> None:
        cs_deriv, fd_deriv, rel_err = complex_step_asl_derivative(band_c_fraction=0.35, po_ratio=0.95)
        self.assertGreater(cs_deriv, 0.0)
        self.assertGreater(fd_deriv, 0.0)
        self.assertLessEqual(rel_err, 1e-6)

    def test_poisoned_negative_controls_fail_deterministically(self) -> None:
        res_chip = evaluate_falsification_negative_control(poison_negative_chip=True)
        self.assertTrue(res_chip["passed"])
        self.assertIn("Correctly caught", res_chip["error"])

        res_dub = evaluate_falsification_negative_control(poison_negative_dub=True)
        self.assertTrue(res_dub["passed"])
        self.assertIn("Correctly caught", res_dub["error"])

        res_po = evaluate_falsification_negative_control(poison_negative_po=True)
        self.assertTrue(res_po["passed"])
        self.assertIn("Correctly caught", res_po["error"])

        clean = evaluate_falsification_negative_control()
        self.assertTrue(clean["passed"])
        self.assertIsNone(clean["error"])

    def test_monotonic_corrector_response(self) -> None:
        # Increasing corrector concentration monotonically reduces exposed misfolded surface
        f_0 = compute_exposed_misfolded_fraction(0.0)
        f_low = compute_exposed_misfolded_fraction(0.5)
        f_high = compute_exposed_misfolded_fraction(5.0)
        self.assertGreater(f_0, f_low)
        self.assertGreater(f_low, f_high)

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
