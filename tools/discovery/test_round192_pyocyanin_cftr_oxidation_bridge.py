#!/usr/bin/env python3
"""Unit tests for Round 192: Pseudomonas aeruginosa Pyocyanin CFTR Inactivation Bridge (T-0873)."""
from __future__ import annotations

import unittest

from tools.discovery.round192_pyocyanin_cftr_oxidation_bridge import (
    FIVE_LAYER_LADDER,
    PARAMETER_PROVENANCE,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    complex_step_asl_derivative,
    compute_active_cftr_fraction,
    compute_asl_height_um,
    compute_mcc_velocity_mm_min,
    compute_ros_generation_flux,
    compute_steady_state_ros_concentration,
    compute_transepithelial_current_pct_wt,
    evaluate_falsification_negative_control,
    evaluate_round192_pyocyanin_cftr_oxidation,
)


class Round192PyocyaninOxidationTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round192_pyocyanin_cftr_oxidation(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Invariants failed for {key}")
            self.assertGreaterEqual(res.ros_generation_flux_um_s, 0.0)
            self.assertGreater(res.total_scavenging_rate_s, 0.0)
            self.assertGreaterEqual(res.steady_state_ros_um, 0.0)
            self.assertGreaterEqual(res.active_cftr_fraction, 0.0)
            self.assertLessEqual(res.active_cftr_fraction, 1.0)
            self.assertGreaterEqual(res.transepithelial_current_pct_wt, 0.0)
            self.assertLessEqual(res.transepithelial_current_pct_wt, 100.0)
            self.assertGreaterEqual(res.airway_surface_liquid_height_um, 3.5)
            self.assertLessEqual(res.airway_surface_liquid_height_um, 8.5)
            self.assertGreaterEqual(res.mucociliary_transport_velocity_mm_min, 0.0)
            self.assertLessEqual(res.complex_step_relative_error, 1e-6)

    def test_curative_presets_cross_kahn_kalai(self) -> None:
        res1 = evaluate_round192_pyocyanin_cftr_oxidation("trikafta_plus_high_dose_catalase_protection")
        self.assertTrue(res1.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res1.active_cftr_fraction, 0.50)
        self.assertGreaterEqual(res1.transepithelial_current_pct_wt, 12.0)
        self.assertGreaterEqual(res1.airway_surface_liquid_height_um, 7.0)
        self.assertGreaterEqual(res1.mucociliary_transport_velocity_mm_min, 1.0)
        self.assertLessEqual(res1.steady_state_ros_um, 2.0)

        res2 = evaluate_round192_pyocyanin_cftr_oxidation("trikafta_plus_moderate_antioxidant")
        self.assertTrue(res2.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res2.active_cftr_fraction, 0.50)
        self.assertGreaterEqual(res2.transepithelial_current_pct_wt, 12.0)
        self.assertGreaterEqual(res2.airway_surface_liquid_height_um, 7.0)
        self.assertLessEqual(res2.steady_state_ros_um, 2.0)

    def test_sub_curative_presets_fail_kahn_kalai(self) -> None:
        # Untreated severe shutdown
        res_severe = evaluate_round192_pyocyanin_cftr_oxidation("untreated_severe_pyocyanin_cftr_shutdown")
        self.assertFalse(res_severe.crosses_kahn_kalai_threshold)
        self.assertLess(res_severe.active_cftr_fraction, 0.15)
        self.assertLess(res_severe.transepithelial_current_pct_wt, 2.0)
        self.assertGreater(res_severe.steady_state_ros_um, 5.0)

        # Modulator alone in infected niche
        res_mod_alone = evaluate_round192_pyocyanin_cftr_oxidation("modulator_alone_without_antioxidant_in_infected_niche")
        self.assertFalse(res_mod_alone.crosses_kahn_kalai_threshold)
        self.assertLess(res_mod_alone.active_cftr_fraction, 0.50)
        self.assertLess(res_mod_alone.transepithelial_current_pct_wt, 12.0)

        # Sub-lethal chronic pyocyanin impairment
        res_sublethal = evaluate_round192_pyocyanin_cftr_oxidation("sub_lethal_chronic_pyocyanin_impairment")
        self.assertFalse(res_sublethal.crosses_kahn_kalai_threshold)
        self.assertLess(res_sublethal.active_cftr_fraction, 0.50)
        self.assertLess(res_sublethal.transepithelial_current_pct_wt, 12.0)

    def test_complex_step_derivative_accuracy(self) -> None:
        cs_deriv, fd_deriv, rel_err = complex_step_asl_derivative(current_pct_wt=20.0)
        self.assertGreater(cs_deriv, 0.0)
        self.assertGreater(fd_deriv, 0.0)
        self.assertLessEqual(rel_err, 1e-6)

    def test_poisoned_negative_controls_fail_deterministically(self) -> None:
        res_pcn = evaluate_falsification_negative_control(poison_negative_pcn=True)
        self.assertTrue(res_pcn["passed"])
        self.assertIn("Correctly caught", res_pcn["error"])

        res_cat = evaluate_falsification_negative_control(poison_negative_cat=True)
        self.assertTrue(res_cat["passed"])
        self.assertIn("Correctly caught", res_cat["error"])

        res_po = evaluate_falsification_negative_control(poison_negative_po=True)
        self.assertTrue(res_po["passed"])
        self.assertIn("Correctly caught", res_po["error"])

        clean = evaluate_falsification_negative_control()
        self.assertTrue(clean["passed"])
        self.assertIsNone(clean["error"])

    def test_monotonic_pyocyanin_dose_response(self) -> None:
        # Increasing pyocyanin monotonically increases ROS flux and decreases active CFTR fraction
        j_low = compute_ros_generation_flux(5.0)
        j_mid = compute_ros_generation_flux(20.0)
        j_high = compute_ros_generation_flux(50.0)
        self.assertLess(j_low, j_mid)
        self.assertLess(j_mid, j_high)

        _, ros_low = compute_steady_state_ros_concentration(j_low, 0.0)
        _, ros_high = compute_steady_state_ros_concentration(j_high, 0.0)
        f_act_low = compute_active_cftr_fraction(ros_low)
        f_act_high = compute_active_cftr_fraction(ros_high)
        self.assertGreater(f_act_low, f_act_high)

    def test_provenance_and_two_column_epistemic_hygiene_metadata(self) -> None:
        self.assertGreaterEqual(len(PRIMARY_CITATIONS), 4)
        self.assertIn("DOI: 10.2353/ajpath.2009.090166", PRIMARY_CITATIONS)
        self.assertIn("DOI: 10.1016/j.freeradbiomed.2008.09.011", PRIMARY_CITATIONS)
        self.assertGreaterEqual(len(PARAMETER_PROVENANCE), 5)
        for p_meta in PARAMETER_PROVENANCE.values():
            self.assertIn(p_meta["tag"], ("empirical", "derived", "phenomenological_heuristic"))
            self.assertIn("source", p_meta)
            self.assertIn("bounds", p_meta)
        self.assertIn("deductive_theorem", TWO_COLUMN_EPISTEMIC_HYGIENE)
        self.assertIn("biophysical_hypothesis", TWO_COLUMN_EPISTEMIC_HYGIENE)
        self.assertIn("statement", TWO_COLUMN_EPISTEMIC_HYGIENE["biophysical_hypothesis"])
        self.assertIn("falsification_condition", TWO_COLUMN_EPISTEMIC_HYGIENE["biophysical_hypothesis"])
        for layer in ["cool_human_fact", "everyday_analogy", "engineering_lens", "formal_theorem", "falsifiable_code_and_toy"]:
            self.assertIn(layer, FIVE_LAYER_LADDER)


if __name__ == "__main__":
    unittest.main()
