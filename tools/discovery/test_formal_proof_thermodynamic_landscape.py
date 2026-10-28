#!/usr/bin/env python3
"""Unit tests for Formal Proof Thermodynamic Landscape & CFTR Allosteric Bridge Model (T-0475)."""
from __future__ import annotations

import sys
import unittest
from pathlib import Path

DISCOVERY_DIR = Path(__file__).resolve().parent
if str(DISCOVERY_DIR) not in sys.path:
    sys.path.insert(0, str(DISCOVERY_DIR))

from formal_proof_thermodynamic_landscape import (  # noqa: E402
    COUNTERFACTUAL_KNOCKOUT_CATALOG,
    FLT_LEAN4_STAGES,
    FORMAL_MATH_CHALLENGE_PROFILES,
    PROOF_THERMODYNAMIC_PARAMETERS,
    compute_stage_free_energy_profile,
    evaluate_falsification_negative_control,
    evaluate_two_way_bridge_isomorphism,
    simulate_counterfactual_axiom_knockout,
)


class FormalProofThermodynamicLandscapeTests(unittest.TestCase):
    def test_flt_stage_theorem_counts_sum_to_29500(self) -> None:
        total = sum(stage.theorem_count for stage in FLT_LEAN4_STAGES)
        self.assertEqual(total, 29500)
        self.assertEqual(len(FLT_LEAN4_STAGES), 5)

    def test_taylor_wiles_is_rate_limiting_transition_state(self) -> None:
        barriers = {
            s.stage_id: compute_stage_free_energy_profile(s, apply_allosteric_bridge=False)[
                "raw_activation_barrier_kcal_mol"
            ]
            for s in FLT_LEAN4_STAGES
        }
        peak_stage = max(barriers, key=barriers.get)
        self.assertEqual(peak_stage, "taylor_wiles_modularity_lift")
        self.assertGreater(barriers["taylor_wiles_modularity_lift"], barriers["frey_hellegouarch_curve"])

    def test_intact_proof_reaches_qed_and_conserves_hypothesis_flux(self) -> None:
        res = simulate_counterfactual_axiom_knockout("wild_type_intact_proof")
        self.assertTrue(res["qed_comparator_verified"])
        self.assertEqual(res["verified_theorems_count"], 29500)
        self.assertGreater(res["final_alluvial_hypothesis_flux_pct"], 60.0)
        self.assertTrue(res["parallel_cf_biophysics"]["curative_10pct_wt_threshold_met"])
        self.assertGreater(res["parallel_cf_biophysics"]["predicted_anion_efflux_pct_wt"], 60.0)

    def test_counterfactual_axiom_knockouts_fracture_pipeline(self) -> None:
        for k_id in (
            "pythagorean_exponent_p2_mutation",
            "frey_semistability_knockout",
            "taylor_wiles_p3_p5_switch_knockout",
            "unverified_sorry_axiom_injection",
        ):
            res = simulate_counterfactual_axiom_knockout(k_id)
            self.assertFalse(res["qed_comparator_verified"], f"{k_id} must not reach QED")
            self.assertLess(res["verified_theorems_count"], 29500)
            self.assertLess(res["final_alluvial_hypothesis_flux_pct"], 5.0)
            self.assertFalse(res["parallel_cf_biophysics"]["curative_10pct_wt_threshold_met"])

    def test_pythagorean_mutation_fractures_at_stage_1(self) -> None:
        res = simulate_counterfactual_axiom_knockout("pythagorean_exponent_p2_mutation")
        self.assertEqual(res["fractured_stage_id"], "frey_hellegouarch_curve")
        self.assertEqual(res["verified_theorems_count"], 0)
        self.assertIn("3^2 + 4^2 = 5^2", res["emergent_counterexample"])

    def test_two_way_bridge_and_ah11_sensitivity_stability(self) -> None:
        report = evaluate_two_way_bridge_isomorphism()
        self.assertEqual(report["answer_tier"], "T2")
        self.assertTrue(report["rate_limiting_stage_preserved_25pct"])
        self.assertTrue(report["full_stage_rank_order_preserved_25pct"])
        self.assertGreater(report["taylor_wiles_allosteric_stabilization_kcal_mol"], 4.0)
        self.assertIn("fermats_last_theorem", FORMAL_MATH_CHALLENGE_PROFILES)
        self.assertIn("formal_math_zeta23", FORMAL_MATH_CHALLENGE_PROFILES)
        self.assertIn("formal_math_3sum_apsp", FORMAL_MATH_CHALLENGE_PROFILES)

    def test_ah9_parameter_provenance_and_ah10_negative_control(self) -> None:
        self.assertGreaterEqual(len(PROOF_THERMODYNAMIC_PARAMETERS), 5)
        neg = evaluate_falsification_negative_control()
        self.assertTrue(neg["passed"])
        self.assertEqual(neg["rule"], "AH-10")

    def test_invalid_inputs_raise_value_error(self) -> None:
        with self.assertRaises(ValueError):
            simulate_counterfactual_axiom_knockout("nonexistent_knockout")
        with self.assertRaises(ValueError):
            compute_stage_free_energy_profile(FLT_LEAN4_STAGES[0], fan_in_weight=-0.5)


if __name__ == "__main__":
    unittest.main()
