#!/usr/bin/env python3
"""Unit tests for Round 194: Bacteriophage AP-PA02 Translocation & Lytic Burst Bridge (T-0877)."""
from __future__ import annotations

import unittest

from tools.discovery.round194_phage_ap_pa02_mucus_translocation_bridge import (
    FIVE_LAYER_LADDER,
    PARAMETER_PROVENANCE,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    complex_step_asl_derivative,
    compute_active_cftr_fraction,
    compute_asl_height_um,
    compute_bacterial_lytic_reduction,
    compute_effective_diffusivity_ratio,
    compute_mcc_velocity_mm_min,
    compute_mucus_mesh_size_nm,
    compute_transepithelial_current_pct_wt,
    evaluate_falsification_negative_control,
    evaluate_round194_phage_ap_pa02_mucus_translocation,
)


class Round194PhageTranslocationTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round194_phage_ap_pa02_mucus_translocation(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Invariants failed for {key}")
            self.assertGreaterEqual(res.phage_dose_pfu_ml, 0.0)
            self.assertGreater(res.mucus_solids_pct, 0.0)
            self.assertGreater(res.mucus_mesh_size_nm, 0.0)
            self.assertGreaterEqual(res.effective_diffusivity_ratio, 0.0)
            self.assertLessEqual(res.effective_diffusivity_ratio, 1.0)
            self.assertGreaterEqual(res.effective_phage_conc_pfu_ml, 0.0)
            self.assertGreaterEqual(res.bacterial_log10_kill, 0.0)
            self.assertGreater(res.post_treatment_bacteria_cfu_ml, 0.0)
            self.assertGreaterEqual(res.active_cftr_fraction, 0.0)
            self.assertLessEqual(res.active_cftr_fraction, 1.0)
            self.assertGreaterEqual(res.transepithelial_current_pct_wt, 0.0)
            self.assertLessEqual(res.transepithelial_current_pct_wt, 100.0)
            self.assertGreaterEqual(res.airway_surface_liquid_height_um, 3.5)
            self.assertLessEqual(res.airway_surface_liquid_height_um, 8.5)
            self.assertGreaterEqual(res.mucociliary_transport_velocity_mm_min, 0.0)
            self.assertLessEqual(res.complex_step_relative_error, 1e-6)

    def test_curative_presets_cross_kahn_kalai(self) -> None:
        res1 = evaluate_round194_phage_ap_pa02_mucus_translocation("ap_pa02_high_dose_trikafta_hydrated")
        self.assertTrue(res1.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res1.bacterial_log10_kill, 2.0)
        self.assertGreaterEqual(res1.active_cftr_fraction, 0.50)
        self.assertGreaterEqual(res1.transepithelial_current_pct_wt, 12.0)
        self.assertGreaterEqual(res1.airway_surface_liquid_height_um, 7.0)
        self.assertGreaterEqual(res1.mucociliary_transport_velocity_mm_min, 1.0)
        self.assertGreaterEqual(res1.effective_diffusivity_ratio, 0.10)

        res2 = evaluate_round194_phage_ap_pa02_mucus_translocation("ap_pa02_standard_dose_trikafta")
        self.assertTrue(res2.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res2.bacterial_log10_kill, 2.0)
        self.assertGreaterEqual(res2.active_cftr_fraction, 0.50)
        self.assertGreaterEqual(res2.transepithelial_current_pct_wt, 12.0)
        self.assertGreaterEqual(res2.airway_surface_liquid_height_um, 7.0)
        self.assertGreaterEqual(res2.effective_diffusivity_ratio, 0.10)

    def test_sub_curative_presets_fail_kahn_kalai(self) -> None:
        # Untreated severe biofilm
        res_severe = evaluate_round194_phage_ap_pa02_mucus_translocation("untreated_severe_pseudomonas_biofilm")
        self.assertFalse(res_severe.crosses_kahn_kalai_threshold)
        self.assertEqual(res_severe.bacterial_log10_kill, 0.0)
        self.assertLess(res_severe.active_cftr_fraction, 0.05)
        self.assertLess(res_severe.transepithelial_current_pct_wt, 1.0)

        # Phage alone without CFTR corrector
        res_phage_alone = evaluate_round194_phage_ap_pa02_mucus_translocation("phage_alone_without_cftr_corrector")
        self.assertFalse(res_phage_alone.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res_phage_alone.bacterial_log10_kill, 2.0)
        self.assertLess(res_phage_alone.transepithelial_current_pct_wt, 12.0)

        # Dehydrated mucus entrapment
        res_entrap = evaluate_round194_phage_ap_pa02_mucus_translocation("dehydrated_mucus_phage_entrapment")
        self.assertFalse(res_entrap.crosses_kahn_kalai_threshold)
        self.assertLess(res_entrap.effective_diffusivity_ratio, 0.10)
        self.assertLess(res_entrap.bacterial_log10_kill, 2.0)
        self.assertLess(res_entrap.transepithelial_current_pct_wt, 12.0)

    def test_complex_step_derivative_accuracy(self) -> None:
        cs_deriv, fd_deriv, rel_err = complex_step_asl_derivative(current_pct_wt=20.0)
        self.assertGreater(cs_deriv, 0.0)
        self.assertGreater(fd_deriv, 0.0)
        self.assertLessEqual(rel_err, 1e-6)

    def test_falsification_negative_controls(self) -> None:
        self.assertTrue(evaluate_falsification_negative_control(poison_negative_dose=True)["passed"])
        self.assertTrue(evaluate_falsification_negative_control(poison_negative_solids=True)["passed"])
        self.assertTrue(evaluate_falsification_negative_control(poison_negative_bacteria=True)["passed"])

    def test_monotonicity_invariants(self) -> None:
        # Increasing solids decreases mesh size and effective diffusivity
        m1 = compute_mucus_mesh_size_nm(mucus_solids_pct=3.0)
        m2 = compute_mucus_mesh_size_nm(mucus_solids_pct=6.0)
        self.assertGreater(m1, m2)

        d1 = compute_effective_diffusivity_ratio(m1)
        d2 = compute_effective_diffusivity_ratio(m2)
        self.assertGreater(d1, d2)

        # Higher effective phage concentration increases log10 kill
        k1, b1 = compute_bacterial_lytic_reduction(effective_phage_conc_pfu_ml=1.0e7)
        k2, b2 = compute_bacterial_lytic_reduction(effective_phage_conc_pfu_ml=1.0e8)
        self.assertGreater(k2, k1)
        self.assertLess(b2, b1)

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

        self.assertGreaterEqual(len(PRIMARY_CITATIONS), 3)
        self.assertIn("DOI: 10.1073/pnas.1508355112", PRIMARY_CITATIONS)
        self.assertIn("DOI: 10.1136/thoraxjnl-2016-209265", PRIMARY_CITATIONS)
        self.assertIn("DOI: 10.1152/physrev.00004.2021", PRIMARY_CITATIONS)


if __name__ == "__main__":
    unittest.main()
