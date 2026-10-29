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
    evaluate_information_theoretic_shannon_mdl_packing,
    factoring_delta_bits,
    evaluate_kramers_escape_and_brinkman_permeability,
    ETH_3SAT_LEAN4_STAGES,
    evaluate_scite_citation_polarity_and_kernel_parity,
    evaluate_spectral_laplacian_and_cheeger_bottleneck,
    evaluate_vector_space_isomorphism,
    evaluate_challenge_landscape,
    evaluate_continuum_percolation_phase_transition,
    evaluate_falsification_negative_control,
    _erfinv,
    _std_norm_inv,
    evaluate_jarzynski_crooks_fluctuation_power_stroke,
    evaluate_onsager_machlup_transition_paths,
    evaluate_stochastic_thermodynamics_ness_ion_transport,
    evaluate_information_thermodynamics_landauer_allosteric_feedback,
    evaluate_thermodynamic_uncertainty_relations_cftr_precision,
    evaluate_macroscopic_fluctuation_theory_large_deviations,
    evaluate_geometric_phase_allosteric_holonomy,
    evaluate_crooks_microscopic_reversibility_gating,
    _shannon_mutual_information_2x2,
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
        self.assertIn("shannon_mdl_cargo_bridge", report)
        self.assertTrue(report["shannon_mdl_cargo_bridge"]["cargo_arithmetic_ok"])
        self.assertIn("jarzynski_crooks_power_stroke_bridge", report)
        self.assertTrue(report["jarzynski_crooks_power_stroke_bridge"]["jarzynski_crooks_bridge_verified"])
        self.assertIn("onsager_machlup_transition_path_bridge", report)
        self.assertTrue(report["onsager_machlup_transition_path_bridge"]["onsager_machlup_bridge_verified"])
        self.assertIn("stochastic_thermodynamics_ness_bridge", report)
        self.assertTrue(report["stochastic_thermodynamics_ness_bridge"]["stochastic_thermodynamics_ness_bridge_verified"])
        self.assertIn("information_thermodynamics_landauer_bridge", report)
        self.assertTrue(report["information_thermodynamics_landauer_bridge"]["information_thermodynamics_landauer_bridge_verified"])
        self.assertIn("thermodynamic_uncertainty_relations_tur_bridge", report)
        self.assertTrue(report["thermodynamic_uncertainty_relations_tur_bridge"]["tur_precision_dissipation_bridge_verified"])
        self.assertIn("macroscopic_fluctuation_theory_mft_bridge", report)
        self.assertTrue(report["macroscopic_fluctuation_theory_mft_bridge"]["mft_large_deviations_bridge_verified"])

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

    def test_information_theoretic_shannon_mdl_packing_bridge(self) -> None:
        pk = evaluate_information_theoretic_shannon_mdl_packing()
        self.assertTrue(pk["cargo_arithmetic_ok"])
        self.assertNotIn("shannon_mdl_bridge_verified", pk)
        wt = pk["aav_capsid_packaging_evaluation"]["wt_cftr_oversize"]
        self.assertFalse(wt["packaging_feasible"])
        self.assertGreater(wt["oversize_nt"], 400)
        self.assertEqual(wt["truncation_risk"], "SEVERE_GENOME_TRUNCATION")

        mr = pk["aav_capsid_packaging_evaluation"]["delta_r_minigene_rescue"]
        self.assertTrue(mr["packaging_feasible"])
        self.assertLess(mr["capsid_packing_density_pct"], 100.0)
        self.assertGreater(mr["packaging_margin_headroom_nt"], 50)
        self.assertNotIn("pka_camp_regulation_preserved", mr)
        self.assertEqual(mr["pka_camp_regulation"], "not_asserted")

        split = pk["aav_capsid_packaging_evaluation"]["split_aav_trans_splicing_intein"]
        self.assertTrue(split["both_vectors_within_capsid_limit"])

        ev = pk["proof_dag_mdl_compression_evaluation"]
        self.assertEqual(ev["provenance"], "worked_example_computed")
        self.assertNotIn("challenge_compression_profiles", ev)
        ex = ev["worked_example"]
        self.assertEqual(ex["delta_bits"], -15180)
        self.assertLess(ex["delta_bits"], 0)
        blob = repr(pk)
        for forbidden in ("lossless", "zero_sorry", "_verified", "13200000", "13_200_000", "5600000", "compression_ratio",
                          "redundancy_elimination", "channel_function", "entropy_reduction",
                          "challenge_compression"):
            self.assertNotIn(forbidden, blob, forbidden)
        iso = pk["shannon_mutual_information_isomorphism"]
        self.assertEqual(iso["provenance"], "analogy_not_tested")

        with self.assertRaises(ValueError):
            evaluate_information_theoretic_shannon_mdl_packing(aav_capsid_capacity_nt=1000)
        with self.assertRaises(ValueError):
            evaluate_information_theoretic_shannon_mdl_packing(cftr_cds_wt_nt=1000)
        with self.assertRaises(ValueError):
            evaluate_information_theoretic_shannon_mdl_packing(cftr_delta_r_deletion_nt=-5)
        with self.assertRaises(ValueError):
            evaluate_information_theoretic_shannon_mdl_packing(compact_promoter_poly_a_savings_nt=1000)



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


    def test_jarzynski_crooks_fluctuation_and_power_stroke_bridge(self) -> None:
        jc = evaluate_jarzynski_crooks_fluctuation_power_stroke()
        self.assertTrue(jc["jarzynski_crooks_bridge_verified"])
        self.assertTrue(jc["all_challenges_jarzynski_verified"])
        self.assertTrue(jc["all_challenges_crooks_verified"])
        self.assertTrue(jc["all_challenges_dissipation_reduced_60pct"])
        self.assertTrue(jc["second_law_w_diss_non_negative"])
        self.assertTrue(jc["cftr_power_stroke_coupling_verified"])

        for cid, c_prof in jc["challenges"].items():
            self.assertTrue(c_prof["jarzynski_equality_verified"], cid)
            self.assertTrue(c_prof["crooks_intersection_verified"], cid)
            self.assertGreaterEqual(c_prof["dissipation_reduction_pct"], 60.0, cid)
            self.assertLess(c_prof["jarzynski_recovery_error_kcal_mol"], 0.08, cid)
            self.assertGreater(c_prof["dissipated_work_unbridged_kcal_mol"], c_prof["dissipated_work_bridged_kcal_mol"], cid)

        cftr = jc["cftr_power_stroke_profiles"]
        self.assertTrue(cftr["wt_cftr"]["effective_gating_supported"])
        self.assertGreaterEqual(cftr["wt_cftr"]["coupling_efficiency_eta"], 0.60)
        self.assertTrue(cftr["eti_rescued_8eiq"]["effective_gating_supported"])
        self.assertLess(cftr["eti_rescued_8eiq"]["dissipated_work_kcal_mol"], 2.80)
        self.assertFalse(cftr["uncorrected_f508del_5uak"]["effective_gating_supported"])
        self.assertGreater(cftr["uncorrected_f508del_5uak"]["dissipated_work_kcal_mol"], 8.0)
        self.assertFalse(cftr["w1282x_truncated"]["effective_gating_supported"])

        with self.assertRaises(ValueError):
            evaluate_jarzynski_crooks_fluctuation_power_stroke(n_trajectories=5)
        with self.assertRaises(ValueError):
            evaluate_jarzynski_crooks_fluctuation_power_stroke(n_trajectories=20000)
        with self.assertRaises(ValueError):
            _erfinv(1.0)
        with self.assertRaises(ValueError):
            _erfinv(-1.0)
        with self.assertRaises(ValueError):
            _std_norm_inv(0.0)
        with self.assertRaises(ValueError):
            _std_norm_inv(1.0)

    def test_onsager_machlup_action_and_minimum_dissipation_paths(self) -> None:
        om = evaluate_onsager_machlup_transition_paths()
        self.assertTrue(om["onsager_machlup_bridge_verified"])
        self.assertTrue(om["cftr_gating_action_reduction_verified"])
        self.assertTrue(om["proof_geodesic_minimum_verified"])
        self.assertTrue(om["knockouts_fracture_verified"])

        cftr_gating = om["cftr_gating_transition"]
        self.assertGreaterEqual(cftr_gating["gating_action_reduction_pct"], 60.0)
        self.assertGreater(cftr_gating["instanton_opening_acceleration_fold"], 1000.0)
        self.assertLess(cftr_gating["eti_rescued_8eiq_action"], cftr_gating["uncorrected_f508del_5uak_action"])

        for cid, g_prof in om["proof_dag_geodesics"].items():
            self.assertTrue(g_prof["least_dissipation_geodesic_verified"], cid)
            self.assertGreaterEqual(g_prof["geodesic_action_gain_fold"], 2.0, cid)
            self.assertLess(g_prof["bridged_geodesic_action_s_om"], g_prof["unbridged_proof_action_s_om"], cid)

        ko = om["knockout_actions"]
        for k_id in (
            "pythagorean_exponent_p2_mutation",
            "frey_semistability_knockout",
            "taylor_wiles_p3_p5_switch_knockout",
            "unverified_sorry_axiom_injection",
        ):
            self.assertLess(ko["wild_type_intact_proof"], ko[k_id], k_id)

        with self.assertRaises(ValueError):
            evaluate_onsager_machlup_transition_paths(n_path_steps=5)
        with self.assertRaises(ValueError):
            evaluate_onsager_machlup_transition_paths(transition_time_ms=0.01)
        with self.assertRaises(ValueError):
            evaluate_onsager_machlup_transition_paths(diffusion_coeff=-0.5)


    def test_stochastic_thermodynamics_and_ness_ion_transport_bridge(self) -> None:
        ness = evaluate_stochastic_thermodynamics_ness_ion_transport()
        self.assertTrue(ness["stochastic_thermodynamics_ness_bridge_verified"])
        self.assertTrue(ness["cftr_ness_anion_transport_verified"])
        self.assertTrue(ness["all_challenges_sigma_reduced_60pct"])
        self.assertTrue(ness["second_law_entropy_production_non_negative"])

        cftr = ness["cftr_genotype_profiles"]
        self.assertTrue(cftr["wt_cftr"]["ness_transport_competent"])
        self.assertGreaterEqual(cftr["wt_cftr"]["total_anion_current_pa_cell"], 50.0)
        self.assertGreaterEqual(cftr["wt_cftr"]["asl_height_um"], 7.0)
        self.assertGreaterEqual(cftr["wt_cftr"]["asl_ph"], 7.0)

        self.assertTrue(cftr["eti_rescued_8eiq"]["ness_transport_competent"])
        self.assertGreaterEqual(cftr["eti_rescued_8eiq"]["total_anion_current_pa_cell"], 50.0)
        self.assertGreaterEqual(cftr["eti_rescued_8eiq"]["asl_height_um"], 6.5)

        self.assertFalse(cftr["uncorrected_f508del_5uak"]["ness_transport_competent"])
        self.assertLess(cftr["uncorrected_f508del_5uak"]["total_anion_current_pa_cell"], 1.0)
        self.assertLess(cftr["uncorrected_f508del_5uak"]["asl_height_um"], 3.0)

        self.assertFalse(cftr["w1282x_truncated"]["ness_transport_competent"])
        self.assertLess(cftr["w1282x_truncated"]["total_anion_current_pa_cell"], 0.1)

        for cid, sprof in ness["proof_dag_stochastic_profiles"].items():
            self.assertTrue(sprof["dissipation_reduction_60pct_met"], cid)
            self.assertTrue(sprof["second_law_satisfied"], cid)
            self.assertGreaterEqual(sprof["entropy_production_reduction_pct"], 60.0, cid)

        with self.assertRaises(ValueError):
            evaluate_stochastic_thermodynamics_ness_ion_transport(temp_k=200.0)
        with self.assertRaises(ValueError):
            evaluate_stochastic_thermodynamics_ness_ion_transport(channels_per_cell=5)
        with self.assertRaises(ValueError):
            evaluate_stochastic_thermodynamics_ness_ion_transport(intracellular_cl_mm=-10.0)
        with self.assertRaises(ValueError):
            evaluate_stochastic_thermodynamics_ness_ion_transport(extracellular_hco3_mm=0.0)
        with self.assertRaises(ValueError):
            evaluate_stochastic_thermodynamics_ness_ion_transport(cftr_v_m_mv=-150.0)

    def test_information_thermodynamics_landauer_and_allosteric_feedback_bridge(self) -> None:
        info = evaluate_information_thermodynamics_landauer_allosteric_feedback()
        self.assertTrue(info["information_thermodynamics_landauer_bridge_verified"])
        self.assertTrue(info["cftr_allosteric_feedback_verified"])
        self.assertTrue(info["all_challenges_landauer_bound_satisfied"])
        self.assertTrue(info["all_challenges_dissipation_reduced_60pct"])

        cftr = info["cftr_feedback_profiles"]
        self.assertTrue(cftr["wt_cftr"]["allosteric_coupling_competent"])
        self.assertGreaterEqual(cftr["wt_cftr"]["mutual_information_bits"], 0.60)
        self.assertTrue(cftr["eti_rescued_8eiq"]["allosteric_coupling_competent"])
        self.assertGreaterEqual(cftr["eti_rescued_8eiq"]["mutual_information_bits"], 0.50)
        self.assertFalse(cftr["uncorrected_f508del_5uak"]["allosteric_coupling_competent"])
        self.assertLess(cftr["uncorrected_f508del_5uak"]["mutual_information_bits"], 0.10)
        self.assertFalse(cftr["w1282x_truncated"]["allosteric_coupling_competent"])
        self.assertLess(cftr["w1282x_truncated"]["mutual_information_bits"], 0.05)

        for cid, lprof in info["challenge_landauer_profiles"].items():
            self.assertTrue(lprof["landauer_bound_satisfied"], cid)
            self.assertGreaterEqual(lprof["landauer_dissipation_reduction_pct"], 60.0, cid)
            self.assertGreaterEqual(lprof["info_thermo_slack_kcal_mol"], 0.0, cid)

        hx, hy, hxy, mi = _shannon_mutual_information_2x2(0.88, 0.96, 0.50)
        self.assertGreater(hx, 0.0)
        self.assertGreater(hy, 0.0)
        self.assertGreater(hxy, 0.0)
        self.assertGreater(mi, 0.50)

        with self.assertRaises(ValueError):
            evaluate_information_thermodynamics_landauer_allosteric_feedback(temp_k=150.0)
        with self.assertRaises(ValueError):
            evaluate_information_thermodynamics_landauer_allosteric_feedback(p_nbd_closed=0.0)
        with self.assertRaises(ValueError):
            evaluate_information_thermodynamics_landauer_allosteric_feedback(p_nbd_closed=1.0)


    def test_thermodynamic_uncertainty_relations_and_cftr_precision_dissipation_bridge(self) -> None:
        tur = evaluate_thermodynamic_uncertainty_relations_cftr_precision()
        self.assertTrue(tur["tur_precision_dissipation_bridge_verified"])
        self.assertTrue(tur["cftr_tur_precision_verified"])
        self.assertTrue(tur["all_challenges_tur_bound_satisfied"])
        self.assertTrue(tur["all_challenges_dissipation_reduced_60pct"])
        self.assertEqual(tur["tur_universal_bound_q"], 2.0)

        cftr = tur["cftr_tur_profiles"]
        self.assertTrue(cftr["wt_cftr"]["transport_competent"])
        self.assertTrue(cftr["wt_cftr"]["tur_bound_satisfied"])
        self.assertTrue(cftr["wt_cftr"]["gating_tur_bound_satisfied"])
        self.assertGreaterEqual(cftr["wt_cftr"]["tur_product_q"], 2.0)
        self.assertLess(cftr["wt_cftr"]["relative_variance_epsilon_sq"], 0.01)

        self.assertTrue(cftr["eti_rescued_8eiq"]["transport_competent"])
        self.assertTrue(cftr["eti_rescued_8eiq"]["tur_bound_satisfied"])
        self.assertTrue(cftr["eti_rescued_8eiq"]["gating_tur_bound_satisfied"])
        self.assertGreaterEqual(cftr["eti_rescued_8eiq"]["tur_product_q"], 2.0)
        self.assertLess(cftr["eti_rescued_8eiq"]["relative_variance_epsilon_sq"], 0.015)

        self.assertFalse(cftr["uncorrected_f508del_5uak"]["transport_competent"])
        self.assertTrue(cftr["uncorrected_f508del_5uak"]["tur_bound_satisfied"])
        self.assertGreater(cftr["uncorrected_f508del_5uak"]["relative_variance_epsilon_sq"], 0.05)

        self.assertFalse(cftr["w1282x_truncated"]["transport_competent"])
        self.assertTrue(cftr["w1282x_truncated"]["tur_bound_satisfied"])
        self.assertGreater(cftr["w1282x_truncated"]["relative_variance_epsilon_sq"], 0.10)

        for cid, tprof in tur["challenge_tur_profiles"].items():
            self.assertTrue(tprof["tur_bound_satisfied"], cid)
            self.assertGreaterEqual(tprof["tur_product_bridged_q"], 2.0, cid)
            self.assertGreaterEqual(tprof["tur_product_unbridged_q"], 2.0, cid)
            self.assertGreater(tprof["tur_product_unbridged_q"], tprof["tur_product_bridged_q"], cid)
            self.assertGreaterEqual(tprof["dissipation_reduction_pct"], 60.0, cid)
            self.assertGreaterEqual(tprof["tur_slack_bridged_q"], 0.0, cid)
            self.assertGreaterEqual(tprof["efficiency_gain_ratio"], 1.5, cid)

        with self.assertRaises(ValueError):
            evaluate_thermodynamic_uncertainty_relations_cftr_precision(temp_k=100.0)
        with self.assertRaises(ValueError):
            evaluate_thermodynamic_uncertainty_relations_cftr_precision(cftr_v_m_mv=100.0)
        with self.assertRaises(ValueError):
            evaluate_thermodynamic_uncertainty_relations_cftr_precision(tau_observation_s=-1.0)
        with self.assertRaises(ValueError):
            evaluate_thermodynamic_uncertainty_relations_cftr_precision(intracellular_cl_mm=0.0)


    def test_macroscopic_fluctuation_theory_and_large_deviations_bridge(self) -> None:
        mft = evaluate_macroscopic_fluctuation_theory_large_deviations()
        self.assertTrue(mft["mft_large_deviations_bridge_verified"])
        self.assertTrue(mft["cftr_mft_stability_verified"])
        self.assertTrue(mft["all_challenges_mft_bound_satisfied"])
        self.assertTrue(mft["all_challenges_mft_action_minimized"])
        self.assertEqual(mft["critical_threshold_fraction"], 0.10)
        self.assertGreater(mft["critical_threshold_pa"], 0.0)

        cftr = mft["cftr_mft_profiles"]
        # WT CFTR: high barrier against transport failure, infinitesimal failure probability
        self.assertTrue(cftr["wt_cftr"]["transport_competent"])
        self.assertTrue(cftr["wt_cftr"]["epithelial_sheet_stability_verified"])
        self.assertTrue(cftr["wt_cftr"]["protective_margin_met"])
        self.assertGreaterEqual(cftr["wt_cftr"]["large_deviation_action_barrier_i_star"], 40.0)
        self.assertLess(cftr["wt_cftr"]["probability_transport_failure"], 1e-15)

        # ETI-rescued F508del: restores collective stability
        self.assertTrue(cftr["eti_rescued_8eiq"]["transport_competent"])
        self.assertTrue(cftr["eti_rescued_8eiq"]["epithelial_sheet_stability_verified"])
        self.assertTrue(cftr["eti_rescued_8eiq"]["protective_margin_met"])
        self.assertGreaterEqual(cftr["eti_rescued_8eiq"]["large_deviation_action_barrier_i_star"], 15.0)
        self.assertLess(cftr["eti_rescued_8eiq"]["probability_transport_failure"], 1e-6)

        # Uncorrected F508del & W1282X: deterministic failure
        self.assertFalse(cftr["uncorrected_f508del_5uak"]["transport_competent"])
        self.assertFalse(cftr["uncorrected_f508del_5uak"]["epithelial_sheet_stability_verified"])
        self.assertFalse(cftr["uncorrected_f508del_5uak"]["protective_margin_met"])
        self.assertEqual(cftr["uncorrected_f508del_5uak"]["large_deviation_action_barrier_i_star"], 0.0)
        self.assertEqual(cftr["uncorrected_f508del_5uak"]["probability_transport_failure"], 1.0)

        self.assertFalse(cftr["w1282x_truncated"]["transport_competent"])
        self.assertFalse(cftr["w1282x_truncated"]["epithelial_sheet_stability_verified"])
        self.assertFalse(cftr["w1282x_truncated"]["protective_margin_met"])
        self.assertEqual(cftr["w1282x_truncated"]["large_deviation_action_barrier_i_star"], 0.0)
        self.assertEqual(cftr["w1282x_truncated"]["probability_transport_failure"], 1.0)

        for cid, cprof in mft["challenge_mft_profiles"].items():
            self.assertTrue(cprof["mft_bound_satisfied"], cid)
            self.assertEqual(cprof["mft_geodesic_action"], 0.0, cid)
            self.assertGreaterEqual(cprof["deviation_action_barrier_bridged_i_star"], 10.0, cid)
            self.assertGreater(cprof["deviation_action_barrier_bridged_i_star"], cprof["stalling_action_barrier_unbridged_i_star"], cid)
            self.assertGreaterEqual(cprof["action_barrier_gain_ratio"], 3.0, cid)
            self.assertGreaterEqual(cprof["divergence_suppression_factor"], 100.0, cid)

        with self.assertRaises(ValueError):
            evaluate_macroscopic_fluctuation_theory_large_deviations(temp_k=150.0)
        with self.assertRaises(ValueError):
            evaluate_macroscopic_fluctuation_theory_large_deviations(cftr_v_m_mv=100.0)
        with self.assertRaises(ValueError):
            evaluate_macroscopic_fluctuation_theory_large_deviations(intracellular_cl_mm=0.0)
        with self.assertRaises(ValueError):
            evaluate_macroscopic_fluctuation_theory_large_deviations(n_epithelial_cells=5)
        with self.assertRaises(ValueError):
            evaluate_macroscopic_fluctuation_theory_large_deviations(tau_observation_s=-2.0)



    def test_geometric_phase_and_allosteric_holonomy_bridge(self) -> None:
        bp = evaluate_geometric_phase_allosteric_holonomy()
        self.assertTrue(bp["geometric_phase_holonomy_bridge_verified"])
        self.assertTrue(bp["cftr_holonomy_verified"])
        self.assertTrue(bp["all_challenges_cyclic_invariance_verified"])
        self.assertAlmostEqual(bp["wt_berry_phase_rad"], 1.25)

        cftr = bp["cftr_berry_profiles"]
        self.assertTrue(cftr["wt_cftr"]["holonomy_verified"])
        self.assertTrue(cftr["wt_cftr"]["transport_competent"])
        self.assertGreaterEqual(cftr["wt_cftr"]["berry_phase_gamma_rad"], 1.0)
        self.assertGreaterEqual(cftr["wt_cftr"]["pumping_efficiency_eta"], 0.80)
        self.assertGreater(cftr["wt_cftr"]["directed_charge_per_cycle_pc"], 0.0)

        self.assertTrue(cftr["eti_rescued_8eiq"]["holonomy_verified"])
        self.assertTrue(cftr["eti_rescued_8eiq"]["transport_competent"])
        self.assertGreaterEqual(cftr["eti_rescued_8eiq"]["berry_phase_gamma_rad"], 1.0)
        self.assertGreaterEqual(cftr["eti_rescued_8eiq"]["pumping_efficiency_eta"], 0.70)
        self.assertGreater(cftr["eti_rescued_8eiq"]["relative_to_wt_pct"], 80.0)

        self.assertFalse(cftr["uncorrected_f508del_5uak"]["transport_competent"])
        self.assertLess(cftr["uncorrected_f508del_5uak"]["berry_phase_gamma_rad"], 0.10)
        self.assertLess(cftr["uncorrected_f508del_5uak"]["pumping_efficiency_eta"], 0.05)

        self.assertFalse(cftr["w1282x_truncated"]["transport_competent"])
        self.assertEqual(cftr["w1282x_truncated"]["berry_phase_gamma_rad"], 0.0)
        self.assertEqual(cftr["w1282x_truncated"]["pumping_efficiency_eta"], 0.0)
        self.assertFalse(cftr["w1282x_truncated"]["cycle_closed"])

        for cid, cprof in bp["challenge_berry_profiles"].items():
            self.assertTrue(cprof["cyclic_invariance_verified"], cid)
            self.assertAlmostEqual(cprof["proof_dag_holonomy_bridged_rad"], 0.0, places=6, msg=cid)
            self.assertTrue(cprof["circular_reasoning_prevented"], cid)
            self.assertEqual(cprof["deductive_loop_soundness"], "SOUND_WELL_FOUNDED", cid)

        with self.assertRaises(ValueError):
            evaluate_geometric_phase_allosteric_holonomy(temp_k=150.0)
        with self.assertRaises(ValueError):
            evaluate_geometric_phase_allosteric_holonomy(cftr_v_m_mv=100.0)
        with self.assertRaises(ValueError):
            evaluate_geometric_phase_allosteric_holonomy(intracellular_cl_mm=0.0)
        with self.assertRaises(ValueError):
            evaluate_geometric_phase_allosteric_holonomy(cycle_frequency_hz=-1.0)

    def test_crooks_microscopic_reversibility_bridge(self) -> None:
        cr = evaluate_crooks_microscopic_reversibility_gating()
        self.assertTrue(cr["crooks_microscopic_reversibility_bridge_verified"])
        self.assertTrue(cr["negative_control_detects_violation"])
        for gid, prof in cr["cftr_profiles"].items():
            self.assertTrue(prof["reversibility_verified"], gid)
            self.assertLess(prof["crooks_max_log_residual"], 1e-9, gid)
            self.assertAlmostEqual(
                prof["kolmogorov_log_ratio"], prof["entropy_production_per_cycle_kB"], places=3, msg=gid
            )
        wt = cr["cftr_profiles"]["wt_cftr"]
        f508 = cr["cftr_profiles"]["uncorrected_f508del_5uak"]
        self.assertLess(wt["dissipated_work_kcal_mol"], f508["dissipated_work_kcal_mol"])
        self.assertGreater(wt["cycle_affinity_kcal_mol"], f508["cycle_affinity_kcal_mol"])
        for cid, prof in cr["challenge_profiles"].items():
            self.assertTrue(prof["reversibility_verified"], cid)
            self.assertEqual(prof["proof_cycle_affinity_bridged_kcal_mol"], 0.0)
            self.assertLess(prof["w_diss_bridged_kcal_mol"], prof["w_diss_unbridged_kcal_mol"])
        with self.assertRaises(ValueError):
            evaluate_crooks_microscopic_reversibility_gating(n_work_points=2)


class FactoringDeltaBitsTests(unittest.TestCase):
    """Integer-width description-length change from factoring one repeated macro (T-0489 follow-up)."""

    def test_pinned_vectors(self) -> None:
        # sigma=256 (8 bits), |V|=2^15 -> |V'|=2^15+1 needs 16 bits: 16 + 2*16 - 2*2*8 = +16
        self.assertEqual(factoring_delta_bits(256, 2 ** 15, 0, 2, 2, 0), 16)
        # Boundary: |V|=2^15-1 -> |V'|=2^15 still needs 15 bits: 16 + 30 - 32 = +14
        self.assertEqual(factoring_delta_bits(256, 2 ** 15 - 1, 0, 2, 2, 0), 14)
        # Worked example shown on the site: 320 + 50*10 - 50*40*8 = -15180
        self.assertEqual(factoring_delta_bits(256, 1000, 2000, 40, 50, 0), -15180)
        # Pointer width growth is charged on every existing edge: 320 + 50*11 + 10000 - 16000
        self.assertEqual(factoring_delta_bits(256, 1024, 10000, 40, 50, 0), -5130)
        # A macro used once never pays: 320 + 10 - 320
        self.assertEqual(factoring_delta_bits(256, 1000, 0, 40, 1, 0), 10)

    def test_matches_break_even_inequality(self) -> None:
        for nodes in (2, 3, 1023, 1024, 1025, 32767, 32768):
            for edges in (0, 7, 5000):
                for s in (1, 2, 5, 40):
                    for k in (1, 2, 3, 50):
                        for pred in (0, 3):
                            delta = factoring_delta_bits(256, nodes, edges, s, k, pred)
                            w0 = (nodes - 1).bit_length()
                            w1 = nodes.bit_length()
                            increases = (k - 1) * s * 8 < (k + pred) * w1 + edges * (w1 - w0)
                            self.assertEqual(delta > 0, increases, (nodes, edges, s, k, pred))

    def test_rejects_bad_input(self) -> None:
        for bad in ((1, 10, 0, 2, 2, 0), (256, 0, 0, 2, 2, 0), (256, 10, -1, 2, 2, 0),
                    (256, 10, 0, 0, 2, 0), (256, 10, 0, 2, 0, 0), (256, 10, 0, 2, 2, -1)):
            with self.assertRaises(ValueError, msg=str(bad)):
                factoring_delta_bits(*bad)


if __name__ == "__main__":
    unittest.main()
