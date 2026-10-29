#!/usr/bin/env python3
"""Unit tests for Formal Proof Thermodynamic Landscape & CFTR Allosteric Bridge Model (T-0479)."""
from __future__ import annotations

import sys
import unittest
from pathlib import Path

DISCOVERY_DIR = Path(__file__).resolve().parent
if str(DISCOVERY_DIR) not in sys.path:
    sys.path.insert(0, str(DISCOVERY_DIR))

from formal_proof_thermodynamic_landscape import (  # noqa: E402
    CHALLENGE_COUNTERFACTUAL_KNOCKOUTS,
    COUNTERFACTUAL_KNOCKOUT_CATALOG,
    FLT_LEAN4_STAGES,
    FORMAL_MATH_CHALLENGE_PROFILES,
    PROOF_THERMODYNAMIC_PARAMETERS,
    compute_stage_free_energy_profile,
    compute_stage_state_vector,
    evaluate_mechanistic_basis_span,
    evaluate_anthropic_proof_visual_geometries,
    evaluate_kramers_escape_and_brinkman_permeability,
    ETH_3SAT_LEAN4_STAGES,
    evaluate_scite_citation_polarity_and_kernel_parity,
    evaluate_spectral_laplacian_and_cheeger_bottleneck,
    evaluate_vector_space_isomorphism,
    evaluate_challenge_landscape,
    evaluate_continuum_percolation_phase_transition,
    evaluate_falsification_negative_control,
    evaluate_noether_conservation_and_invariant_transfer,
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
        with self.assertRaises(ValueError):
            evaluate_challenge_landscape("unknown_challenge")

    def test_all_three_anthropic_challenge_landscapes_sum_and_verify(self) -> None:
        expected_totals = {
            "fermats_last_theorem": 29500,
            "formal_math_zeta23": 1420,
            "formal_math_3sum_apsp": 980,
            "formal_math_eth": 493,
        }
        for cid, total in expected_totals.items():
            land = evaluate_challenge_landscape(cid, apply_allosteric_bridge=True)
            self.assertEqual(land["total_theorems"], total)
            self.assertEqual(len(land["stages"]), 5)
            self.assertGreater(land["final_alluvial_hypothesis_flux_pct"], 60.0)
            self.assertEqual(
                land["unbridged_peak_stage_id"],
                FORMAL_MATH_CHALLENGE_PROFILES[cid]["peak_stage_id"],
            )
        iso = evaluate_two_way_bridge_isomorphism()
        self.assertTrue(iso["all_challenges_rank_order_preserved_25pct"])




    def test_noether_kirchhoff_flux_conservation_and_invariant_transfer(self) -> None:
        gate = evaluate_noether_conservation_and_invariant_transfer()
        self.assertTrue(gate["all_invariants_passed"])
        self.assertTrue(gate["noether_kirchhoff_flux_conservation"])
        self.assertTrue(gate["monotone_allosteric_duality"])
        self.assertTrue(gate["zero_sorry_comparator_seal"])
        for cid, c_rep in gate["challenges"].items():
            self.assertTrue(c_rep["kirchhoff_flux_conservation_verified"], cid)
            self.assertLess(c_rep["flux_conservation_residual_pct"], 1e-9)
            self.assertAlmostEqual(
                c_rep["terminal_preserved_flux_pct"] + c_rep["cumulative_dissipated_flux_pct"],
                100.0,
                places=3,
            )


    def test_continuum_percolation_phase_transition(self) -> None:
        perc = evaluate_continuum_percolation_phase_transition()
        self.assertTrue(perc["phase_transition_verified"])
        self.assertAlmostEqual(perc["continuum_critical_volume_fraction_phi_c"], 0.28957, places=4)
        self.assertTrue(perc["compact_peg_lnp_in_cf_sputum"]["supercritical_percolation"])
        self.assertGreater(perc["compact_peg_lnp_in_cf_sputum"]["spanning_cluster_probability_p_inf"], 0.5)
        self.assertFalse(perc["unshielded_130nm_lnp_in_cf_sputum"]["supercritical_percolation"])
        self.assertEqual(perc["unshielded_130nm_lnp_in_cf_sputum"]["spanning_cluster_probability_p_inf"], 0.0)
        self.assertFalse(perc["cftr_allosteric_network"]["uncorrected_f508del_percolates"])
        self.assertTrue(perc["cftr_allosteric_network"]["eti_rescued_percolates"])
        with self.assertRaises(ValueError):
            evaluate_continuum_percolation_phase_transition(particle_radius_nm=0.0)


    def test_hilbert_vector_space_isomorphism_and_orthogonal_knockouts(self) -> None:
        vec_iso = evaluate_vector_space_isomorphism()
        self.assertTrue(vec_iso["pythagorean_unitary_norm_verified"])
        self.assertTrue(vec_iso["allosteric_bridge_vector_rotation_verified"])
        self.assertTrue(vec_iso["orthogonal_knockout_subspace_hierarchy_verified"])
        self.assertTrue(vec_iso["multi_axis_mechanistic_span_verified"])

        # Verify 4D non-compensatory span exposes scalar fallacy on W1282X + ETI alone
        w1282x_eti = vec_iso["mechanistic_span_matrix"]["W1282X_with_ETI_only"]
        self.assertFalse(w1282x_eti["basis_fully_spanned"])
        self.assertEqual(w1282x_eti["dominant_unspanned_basis_axis"], "e1_transcript_axiom_entry")
        self.assertGreater(w1282x_eti["scalar_fallacy_overestimate_pct"], 50.0)

        w1282x_tri = vec_iso["mechanistic_span_matrix"]["W1282X_with_trimodal"]
        self.assertTrue(w1282x_tri["basis_fully_spanned"])
        self.assertGreater(w1282x_tri["directional_cosine_alignment"], 0.99)

        with self.assertRaises(ValueError):
            evaluate_mechanistic_basis_span([0.5, 0.5], [0.5, 0.5, 0.5, 0.5])

        # Terminal QED stage in intact FLT aligns near 0 deg (cos theta > 0.99)
        flt_vec = vec_iso["challenge_vector_profiles"]["fermats_last_theorem"]
        self.assertGreater(flt_vec["terminal_qed_cos_theta"], 0.99)
        self.assertLess(flt_vec["peak_stage_bridge_rotation_deg"], -15.0)

        # Fractured stage deflects orthogonally (> 84 deg, cos theta < 0.10)
        frac_vec = compute_stage_state_vector(11.69, 3.49, is_fractured=True)
        self.assertLess(frac_vec["psi_parallel_cos_theta"], 0.10)
        self.assertGreater(frac_vec["psi_orthogonal_sin_theta"], 0.99)
        self.assertGreater(frac_vec["phase_deflection_angle_deg"], 84.0)
        self.assertAlmostEqual(frac_vec["unitary_norm_squared"], 1.0, places=9)

        with self.assertRaises(ValueError):
            compute_stage_state_vector(-1.0)



    def test_kramers_escape_and_stokes_brinkman_permeability_bridge(self) -> None:
        kb = evaluate_kramers_escape_and_brinkman_permeability()
        self.assertTrue(kb["kramers_brinkman_bridge_verified"])
        self.assertGreater(kb["darcy_permeability_collapse_fold"], 6.0)
        self.assertGreater(kb["healthy_2pct_asl"]["brinkman_screening_length_nm"], 15.0)
        self.assertLess(kb["cf_8pct_dehydrated_mucus"]["brinkman_screening_length_nm"], 8.5)

        flt_k = kb["challenge_kramers_escape_profiles"]["fermats_last_theorem"]
        self.assertGreater(flt_k["kramers_escape_speedup_fold"], 5000.0)
        self.assertLess(
            flt_k["mean_first_passage_time_ms_bridged"],
            flt_k["mean_first_passage_time_ms_unbridged"],
        )
        with self.assertRaises(ValueError):
            evaluate_kramers_escape_and_brinkman_permeability(
                healthy_solid_fraction=0.10, cf_dehydrated_solid_fraction=0.05
            )

    def test_anthropic_solved_proof_visual_geometries(self) -> None:
        geom = evaluate_anthropic_proof_visual_geometries()
        self.assertTrue(geom["all_three_geometries_verified"])
        self.assertEqual(geom["fermats_last_theorem"]["terminal_cusp_space_dim"], 0)
        self.assertEqual(geom["fermats_last_theorem"]["pythagorean_p2_lattice_witness"]["lhs"], 25)
        self.assertAlmostEqual(geom["formal_math_zeta23"]["critical_line_re_s"], 0.5)
        self.assertEqual(
            geom["formal_math_3sum_apsp"]["zero_triangle_witness_weights"]["cycle_sum"], 0
        )
        self.assertEqual(geom["formal_math_3sum_apsp"]["bit_packing_speedup_factor"], 64.0)



    def test_spectral_laplacian_and_cheeger_bottleneck_bridge(self) -> None:
        sp = evaluate_spectral_laplacian_and_cheeger_bottleneck()
        self.assertTrue(sp["spectral_cheeger_bridge_verified"])
        self.assertGreater(sp["proof_dag_bottleneck_cut"]["fiedler_connectivity_gain_fold"], 6.0)
        self.assertGreater(sp["cftr_allosteric_contact_graph"]["allosteric_fiedler_gain_fold"], 7.0)
        for cid, c_prof in sp["challenge_spectral_profiles"].items():
            self.assertTrue(c_prof["unbridged"]["cheeger_inequality_satisfied"], cid)
            self.assertTrue(c_prof["bridged"]["cheeger_inequality_satisfied"], cid)
            self.assertEqual(c_prof["knockout_disconnected_lambda_2"], 0.0)
        with self.assertRaises(ValueError):
            evaluate_spectral_laplacian_and_cheeger_bottleneck(
                unbridged_cut_conductance=0.60, bridged_cut_conductance=0.20
            )



    def test_eth_3sat_refutation_and_scite_citation_polarity_parity(self) -> None:
        self.assertEqual(sum(s.theorem_count for s in ETH_3SAT_LEAN4_STAGES), 493)
        eth_land = evaluate_challenge_landscape("formal_math_eth", apply_allosteric_bridge=True)
        self.assertEqual(eth_land["total_theorems"], 493)
        self.assertEqual(eth_land["unbridged_peak_stage_id"], "clause_vs_variable_alternation_bridge")
        self.assertGreater(eth_land["final_alluvial_hypothesis_flux_pct"], 60.0)

        geom = evaluate_anthropic_proof_visual_geometries()
        self.assertTrue(geom["all_four_geometries_verified"])
        self.assertEqual(geom["formal_math_eth"]["s3_infimum_exponent"], 0.0)
        self.assertEqual(geom["formal_math_eth"]["paper_index_theorems_count"], 493)
        self.assertAlmostEqual(geom["formal_math_zeta23"]["xiprime_simple_critical_line_proportion_quartic"], 0.86864, places=5)

        scite_rep = evaluate_scite_citation_polarity_and_kernel_parity()
        self.assertTrue(scite_rep["zero_editorial_notices_verified"])
        self.assertTrue(scite_rep["all_four_benchmarks_kernel_and_scite_verified"])
        self.assertIn("Scite ML citation classification", scite_rep["provenance_tag"])
        self.assertEqual(len(scite_rep["benchmarks"]), 4)
        for cid, ko_map in CHALLENGE_COUNTERFACTUAL_KNOCKOUTS.items():
            self.assertEqual(len(ko_map), 5)
            unbridged_peak = evaluate_challenge_landscape(cid, apply_allosteric_bridge=False)["unbridged_peak_barrier_kcal_mol"]
            self.assertIn(f"{unbridged_peak:.2f} kcal/mol", ko_map["taylor_wiles_p3_p5_switch_knockout"]["counterexample"])


if __name__ == "__main__":
    unittest.main()
