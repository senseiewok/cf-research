#!/usr/bin/env python3
"""Unit tests for Round 189: Prime Editing PegRNA Flap Annealing & Nicking Kinetics Bridge (T-0870)."""
from __future__ import annotations

import unittest

from tools.discovery.round189_prime_editing_pegrna_kinetics_bridge import (
    FIVE_LAYER_LADDER,
    PARAMETER_PROVENANCE,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    complex_step_asl_derivative,
    compute_asl_height_um,
    compute_flap_resolution_efficiency,
    compute_mcc_velocity_mm_min,
    compute_net_prime_editing_efficiency,
    compute_pbs_hybridization_fraction,
    compute_rescued_cftr_expression,
    compute_rt_extension_efficiency,
    evaluate_falsification_negative_control,
    evaluate_round189_prime_editing,
)


class Round189PrimeEditingTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round189_prime_editing(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Invariants failed for {key}")
            self.assertGreaterEqual(res.prime_editing_efficiency, 0.0)
            self.assertLessEqual(res.prime_editing_efficiency, 1.0)
            self.assertGreaterEqual(res.pbs_bound_fraction, 0.0)
            self.assertLessEqual(res.pbs_bound_fraction, 1.0)
            self.assertGreaterEqual(res.rt_elongation_efficiency, 0.0)
            self.assertLessEqual(res.rt_elongation_efficiency, 1.0)
            self.assertGreaterEqual(res.rescued_mrna_fraction, 0.0)
            self.assertLessEqual(res.rescued_mrna_fraction, 1.0)
            self.assertGreaterEqual(res.restored_protein_fraction, 0.0)
            self.assertLessEqual(res.restored_protein_fraction, 1.0)
            self.assertGreaterEqual(res.transepithelial_current_pct_wt, 0.0)
            self.assertLessEqual(res.transepithelial_current_pct_wt, 100.0)
            self.assertGreaterEqual(res.airway_surface_liquid_height_um, 3.5)
            self.assertLessEqual(res.airway_surface_liquid_height_um, 8.5)
            self.assertGreaterEqual(res.mucociliary_transport_velocity_mm_min, 0.0)
            self.assertLessEqual(res.complex_step_relative_error, 1e-6)

    def test_pe3max_optimized_w1282x_crosses_kahn_kalai(self) -> None:
        res = evaluate_round189_prime_editing("pe3max_optimized_w1282x")
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res.prime_editing_efficiency, 0.15)
        self.assertGreaterEqual(res.restored_protein_fraction, 0.12)
        self.assertGreaterEqual(res.transepithelial_current_pct_wt, 12.0)
        self.assertGreaterEqual(res.airway_surface_liquid_height_um, 7.0)
        self.assertGreaterEqual(res.mucociliary_transport_velocity_mm_min, 1.0)

    def test_pe3max_optimized_g542x_crosses_kahn_kalai(self) -> None:
        res = evaluate_round189_prime_editing("pe3max_optimized_g542x")
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res.prime_editing_efficiency, 0.15)
        self.assertGreaterEqual(res.restored_protein_fraction, 0.12)
        self.assertGreaterEqual(res.airway_surface_liquid_height_um, 7.0)
        self.assertGreaterEqual(res.mucociliary_transport_velocity_mm_min, 1.0)

    def test_sub_curative_presets_fail_kahn_kalai(self) -> None:
        # Sub-optimal PBS flap
        sub = evaluate_round189_prime_editing("pe2_suboptimal_pbs_flap")
        self.assertFalse(sub.crosses_kahn_kalai_threshold)
        self.assertLess(sub.prime_editing_efficiency, 0.15)
        self.assertLess(sub.restored_protein_fraction, 0.12)
        self.assertLess(sub.airway_surface_liquid_height_um, 7.0)

        # Untreated baseline
        untreated = evaluate_round189_prime_editing("untreated_nonsense_baseline")
        self.assertFalse(untreated.crosses_kahn_kalai_threshold)
        self.assertEqual(untreated.prime_editing_efficiency, 0.0)
        self.assertEqual(untreated.mucociliary_transport_velocity_mm_min, 0.0)

        # High MMR reversion barrier
        mmr = evaluate_round189_prime_editing("high_mmr_reversion_barrier")
        self.assertFalse(mmr.crosses_kahn_kalai_threshold)
        self.assertLess(mmr.prime_editing_efficiency, 0.15)
        self.assertLess(mmr.airway_surface_liquid_height_um, 7.0)

    def test_complex_step_derivative_accuracy(self) -> None:
        cs_deriv, fd_deriv, rel_err = complex_step_asl_derivative(protein_fraction=0.30, gamma_activity=1.00)
        self.assertGreater(cs_deriv, 0.0)
        self.assertGreater(fd_deriv, 0.0)
        self.assertLessEqual(rel_err, 1e-6)

    def test_poisoned_negative_controls_fail_deterministically(self) -> None:
        res_temp = evaluate_falsification_negative_control(poison_negative_temp=True)
        self.assertTrue(res_temp["passed"])
        self.assertIn("Correctly caught", res_temp["error"])

        res_dntp = evaluate_falsification_negative_control(poison_negative_dntp=True)
        self.assertTrue(res_dntp["passed"])
        self.assertIn("Correctly caught", res_dntp["error"])

        res_gamma = evaluate_falsification_negative_control(poison_negative_gamma=True)
        self.assertTrue(res_gamma["passed"])
        self.assertIn("Correctly caught", res_gamma["error"])

        clean = evaluate_falsification_negative_control()
        self.assertTrue(clean["passed"])
        self.assertIsNone(clean["error"])

    def test_monotonic_pbs_binding_curve(self) -> None:
        # Weaker binding (less negative Delta G) gives lower f_PBS
        f_weak = compute_pbs_hybridization_fraction(-5.0)
        f_mid = compute_pbs_hybridization_fraction(-10.0)
        f_strong = compute_pbs_hybridization_fraction(-15.0)
        self.assertLess(f_weak, f_mid)
        self.assertLess(f_mid, f_strong)

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
