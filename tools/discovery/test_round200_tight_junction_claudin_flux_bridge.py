#!/usr/bin/env python3
"""Unit tests for Round 200: Bicellular Tight Junction Claudin-2/Claudin-4 Paracellular Ion Permselectivity Invariant Bridge (T-0889)."""
from __future__ import annotations

import unittest

from tools.discovery.round200_tight_junction_claudin_flux_bridge import (
    FIVE_LAYER_LADDER,
    PARAMETER_PROVENANCE,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    complex_step_teer_derivative,
    compute_asl_preservation_factor,
    compute_paracellular_conductance,
    compute_paracellular_permeabilities,
    compute_paracellular_sodium_flux,
    compute_teer,
    evaluate_falsification_negative_control,
    evaluate_round200_tight_junction_claudin_flux,
)


class Round200TightJunctionClaudinFluxBridgeTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round200_tight_junction_claudin_flux(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Invariants failed for {key}")
            self.assertGreaterEqual(res.e_cldn4, 0.0)
            self.assertGreaterEqual(res.e_cldn2, 0.0)
            self.assertGreater(res.teer_ohm_cm2, 0.0)
            self.assertGreater(res.g_para, 0.0)
            self.assertGreater(res.p_na_cm_s, 0.0)
            self.assertGreater(res.p_cl_cm_s, 0.0)
            self.assertGreater(res.permselectivity_ratio, 0.0)
            self.assertGreaterEqual(res.j_na_para_nmol_cm2_min, 0.0)
            self.assertGreaterEqual(res.phi_asl_preservation, 0.0)
            self.assertLessEqual(res.phi_asl_preservation, 1.0)
            self.assertLessEqual(res.complex_step_relative_error, 1e-3)

    def test_curative_presets_cross_kahn_kalai(self) -> None:
        res_curative = evaluate_round200_tight_junction_claudin_flux(
            "balanced_physiological_curative_junction"
        )
        self.assertTrue(res_curative.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res_curative.teer_ohm_cm2, 500.0)
        self.assertGreaterEqual(res_curative.permselectivity_ratio, 1.8)
        self.assertLessEqual(res_curative.permselectivity_ratio, 3.8)
        self.assertGreaterEqual(res_curative.j_na_para_nmol_cm2_min, 20.0)
        self.assertGreaterEqual(res_curative.phi_asl_preservation, 0.40)

        res_super = evaluate_round200_tight_junction_claudin_flux(
            "reinforced_barrier_super_curative"
        )
        self.assertTrue(res_super.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res_super.teer_ohm_cm2, 600.0)
        self.assertGreaterEqual(res_super.phi_asl_preservation, 0.45)

    def test_sub_curative_presets_fail_kahn_kalai(self) -> None:
        res_leaky = evaluate_round200_tight_junction_claudin_flux(
            "cf_inflammatory_leaky_barrier_breakdown"
        )
        self.assertFalse(res_leaky.crosses_kahn_kalai_threshold)
        self.assertLess(res_leaky.teer_ohm_cm2, 200.0)
        self.assertLess(res_leaky.phi_asl_preservation, 0.40)

        res_partial = evaluate_round200_tight_junction_claudin_flux(
            "sub_curative_attenuated_sealing"
        )
        self.assertFalse(res_partial.crosses_kahn_kalai_threshold)
        self.assertLess(res_partial.teer_ohm_cm2, 500.0)

        res_pore_def = evaluate_round200_tight_junction_claudin_flux(
            "pore_deficient_impermeable_arrest"
        )
        self.assertFalse(res_pore_def.crosses_kahn_kalai_threshold)
        self.assertLess(res_pore_def.j_na_para_nmol_cm2_min, 20.0)
        self.assertLess(res_pore_def.phi_asl_preservation, 0.40)

    def test_monotonicity_invariants(self) -> None:
        # Increasing Claudin-4 increases TEER
        g1 = float(compute_paracellular_conductance(1.0, 0.5))
        g2 = float(compute_paracellular_conductance(1.0, 2.0))
        teer1 = float(compute_teer(g1))
        teer2 = float(compute_teer(g2))
        self.assertGreater(teer2, teer1)

        # Increasing Claudin-2 increases paracellular conductance and decreases TEER
        g3 = float(compute_paracellular_conductance(3.0, 2.0))
        teer3 = float(compute_teer(g3))
        self.assertLess(teer3, teer2)

        # Increasing Claudin-2 increases paracellular Na+ permeability
        p_na_low, _, _ = compute_paracellular_permeabilities(0.5, 1.0)
        p_na_high, _, _ = compute_paracellular_permeabilities(2.5, 1.0)
        self.assertGreater(float(p_na_high), float(p_na_low))

    def test_complex_step_derivative_accuracy(self) -> None:
        cs_deriv, fd_deriv, rel_err = complex_step_teer_derivative(1.5, 1.0)
        self.assertGreater(cs_deriv, 0.0)
        self.assertLessEqual(rel_err, 1e-3)

    def test_falsification_negative_controls(self) -> None:
        neg_cldn4 = evaluate_falsification_negative_control(poison_negative_cldn4=True)
        self.assertTrue(neg_cldn4["passed"])
        self.assertIn("rejected", neg_cldn4["reason"].lower())

        neg_cldn2 = evaluate_falsification_negative_control(poison_negative_cldn2=True)
        self.assertTrue(neg_cldn2["passed"])
        self.assertIn("rejected", neg_cldn2["reason"].lower())

        neg_conc = evaluate_falsification_negative_control(poison_negative_concentration=True)
        self.assertTrue(neg_conc["passed"])
        self.assertIn("rejected", neg_conc["reason"].lower())

    def test_provenance_and_epistemic_hygiene_metadata(self) -> None:
        self.assertGreaterEqual(len(PRIMARY_CITATIONS), 2)
        for cit in PRIMARY_CITATIONS:
            self.assertTrue(cit.startswith("DOI: 10."))

        self.assertIn("deductive_column_a", TWO_COLUMN_EPISTEMIC_HYGIENE)
        self.assertIn("popperian_column_b", TWO_COLUMN_EPISTEMIC_HYGIENE)
        self.assertIn("formal_theorem", TWO_COLUMN_EPISTEMIC_HYGIENE["deductive_column_a"])
        self.assertIn("falsifiable_hypothesis", TWO_COLUMN_EPISTEMIC_HYGIENE["popperian_column_b"])

        for layer in ("cool_human_fact", "everyday_analogy", "engineering_lens", "formal_theorem", "falsifiable_code_and_toy"):
            self.assertIn(layer, FIVE_LAYER_LADDER)
            self.assertGreater(len(FIVE_LAYER_LADDER[layer]), 20)

        self.assertGreaterEqual(len(PARAMETER_PROVENANCE), 5)
        for k, v in PARAMETER_PROVENANCE.items():
            self.assertIn("value", v)
            self.assertIn("tag", v)
            self.assertIn(v["tag"], ("empirical", "derived", "phenomenological_heuristic"))
            self.assertIn("source", v)
            self.assertIn("bounds", v)


if __name__ == "__main__":
    unittest.main()
