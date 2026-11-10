#!/usr/bin/env python3
"""Unit tests for Round 195: Inhaled Alginate Lyase Sputum Liquefaction & Mucociliary Transport Bridge (T-0880)."""
from __future__ import annotations

import unittest

from tools.discovery.round195_alginate_lyase_sputum_liquefaction_bridge import (
    FIVE_LAYER_LADDER,
    PARAMETER_PROVENANCE,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    complex_step_viscosity_derivative,
    compute_mucociliary_velocity,
    compute_residual_alginate_fraction,
    compute_sputum_complex_viscosity,
    compute_sputum_yield_stress,
    evaluate_falsification_negative_control,
    evaluate_round195_alginate_lyase_sputum_liquefaction,
)


class Round195AlginateLyaseSputumLiquefactionTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round195_alginate_lyase_sputum_liquefaction(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Invariants failed for {key}")
            self.assertGreaterEqual(res.enz_inh_ug_ml, 0.0)
            self.assertGreaterEqual(res.phi_hydration, 0.0)
            self.assertLessEqual(res.phi_hydration, 1.0)
            self.assertGreaterEqual(res.f_dna, 0.0)
            self.assertLessEqual(res.f_dna, 1.0)
            self.assertGreaterEqual(res.f_alg, 0.0)
            self.assertLessEqual(res.f_alg, 1.0)
            self.assertGreaterEqual(res.complex_viscosity_pa_s, 0.05)
            self.assertGreaterEqual(res.yield_stress_pa, 0.0)
            self.assertGreaterEqual(res.mucociliary_velocity_mm_min, 0.0)
            self.assertLessEqual(res.complex_step_relative_error, 1e-4)

    def test_curative_presets_cross_kahn_kalai(self) -> None:
        res1 = evaluate_round195_alginate_lyase_sputum_liquefaction("alginate_lyase_triple_synergy_curative")
        self.assertTrue(res1.crosses_kahn_kalai_threshold)
        self.assertLessEqual(res1.complex_viscosity_pa_s, 3.00)
        self.assertLessEqual(res1.yield_stress_pa, 0.50)
        self.assertGreaterEqual(res1.mucociliary_velocity_mm_min, 1.00)

    def test_sub_curative_presets_fail_kahn_kalai(self) -> None:
        # Untreated mucoid CF sputum
        res_untreated = evaluate_round195_alginate_lyase_sputum_liquefaction("untreated_mucoid_pseudomonas_sputum")
        self.assertFalse(res_untreated.crosses_kahn_kalai_threshold)
        self.assertEqual(res_untreated.f_alg, 1.0)
        self.assertGreaterEqual(res_untreated.complex_viscosity_pa_s, 100.0)
        self.assertGreater(res_untreated.yield_stress_pa, 2.0)  # Exceeds ciliary stall threshold
        self.assertEqual(res_untreated.mucociliary_velocity_mm_min, 0.0)

        # Sub-therapeutic enzyme dose
        res_sub = evaluate_round195_alginate_lyase_sputum_liquefaction("sub_therapeutic_enzyme_dose")
        self.assertFalse(res_sub.crosses_kahn_kalai_threshold)
        self.assertGreater(res_sub.complex_viscosity_pa_s, 1.50)
        self.assertEqual(res_sub.mucociliary_velocity_mm_min, 0.0)

        # Inactivated enzyme degradation
        res_inact = evaluate_round195_alginate_lyase_sputum_liquefaction("inactivated_enzyme_proteolytic_degradation")
        self.assertFalse(res_inact.crosses_kahn_kalai_threshold)
        self.assertGreater(res_inact.complex_viscosity_pa_s, 1.50)
        self.assertEqual(res_inact.mucociliary_velocity_mm_min, 0.0)

    def test_complex_step_derivative_accuracy(self) -> None:
        cs_deriv, fd_deriv, rel_err = complex_step_viscosity_derivative(
            enz_inh_ug_ml=30.0,
            phi_hydration=0.50,
            f_dna=0.20,
        )
        self.assertLess(cs_deriv, 0.0, "Viscosity derivative with respect to enzyme must be negative")
        self.assertLess(fd_deriv, 0.0)
        self.assertLessEqual(rel_err, 1e-4)

    def test_falsification_negative_controls(self) -> None:
        self.assertTrue(evaluate_falsification_negative_control(poison_negative_enzyme=True)["passed"])
        self.assertTrue(evaluate_falsification_negative_control(poison_negative_viscosity=True)["passed"])
        self.assertTrue(evaluate_falsification_negative_control(poison_negative_yield_stress=True)["passed"])

    def test_monotonicity_invariants(self) -> None:
        # Increasing enzyme concentration monotonically decreases residual alginate and yield stress
        f1 = compute_residual_alginate_fraction(enz_inh_ug_ml=10.0)
        f2 = compute_residual_alginate_fraction(enz_inh_ug_ml=30.0)
        self.assertGreater(f1, f2)

        tau1 = compute_sputum_yield_stress(f_alg=f1, phi_hydration=0.2)
        tau2 = compute_sputum_yield_stress(f_alg=f2, phi_hydration=0.2)
        self.assertGreater(tau1, tau2)

        # Decreased yield stress and viscosity increases MCC velocity
        _, visc1 = compute_sputum_complex_viscosity(f1, phi_hydration=0.2, f_dna=0.5)
        _, visc2 = compute_sputum_complex_viscosity(f2, phi_hydration=0.2, f_dna=0.5)
        self.assertGreater(visc1, visc2)

    def test_provenance_and_epistemic_hygiene_metadata(self) -> None:
        self.assertGreaterEqual(len(PARAMETER_PROVENANCE), 5)
        for p_name, p_meta in PARAMETER_PROVENANCE.items():
            self.assertIn(p_meta["tag"], ("empirical", "derived", "phenomenological_heuristic"))
            self.assertIn("source", p_meta)
            self.assertIn("bounds", p_meta)

        self.assertIn("deductive_theorem", TWO_COLUMN_EPISTEMIC_HYGIENE)
        self.assertIn("biophysical_hypothesis", TWO_COLUMN_EPISTEMIC_HYGIENE)
        self.assertIn("falsification_condition", TWO_COLUMN_EPISTEMIC_HYGIENE["biophysical_hypothesis"])

        self.assertEqual(len(FIVE_LAYER_LADDER), 5)
        for layer in ("cool_human_fact", "everyday_analogy", "engineering_lens", "formal_theorem", "falsifiable_code_and_toy"):
            self.assertIn(layer, FIVE_LAYER_LADDER)

        self.assertGreaterEqual(len(PRIMARY_CITATIONS), 2)


if __name__ == "__main__":
    unittest.main()
