#!/usr/bin/env python3
"""Unit tests for Round 204: EHD Taylor-Cone Atomization, Mucin H1 Sterics, & CBF Safe Redosing Bridge (T-0897)."""
from __future__ import annotations

import unittest

from tools.discovery.round204_ehd_taylor_cbf_redosing_bridge import (
    CROSS_DOMAIN_GOVERNANCE,
    FIVE_LAYER_LADDER,
    PARAMETER_PROVENANCE,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    RESEARCH_BOUNDARY,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    complex_step_delta_t_derivative,
    compute_cbf_redosing_invariance,
    compute_ehd_atomization_and_lnp_integrity,
    compute_h1_mucus_persistent_homology_transmission,
    evaluate_falsification_negative_control,
    evaluate_round204_ehd_taylor_cbf_redosing,
)


class Round204EHDTaylorCBFRedosingBridgeTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round204_ehd_taylor_cbf_redosing(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Invariants failed for {key}")
            self.assertGreater(res.eta_intact, 0.0)
            self.assertLessEqual(res.eta_intact, 1.0)
            self.assertGreaterEqual(res.r_h_post_nm, res.r_h0_nm)
            self.assertGreaterEqual(res.f_bronch_deposition, 0.0)
            self.assertLessEqual(res.f_bronch_deposition, 1.0)
            self.assertGreaterEqual(res.t_h1_mucus_transmission, 0.0)
            self.assertLessEqual(res.t_h1_mucus_transmission, 1.0)
            self.assertGreaterEqual(res.phi_clearance_retention, 0.0)
            self.assertLessEqual(res.phi_clearance_retention, 1.0)
            self.assertGreaterEqual(res.m_cycle_maintenance, 0.0)
            self.assertLessEqual(res.m_cycle_maintenance, 1.0)
            self.assertGreaterEqual(res.j_chronic_delivery_index, 0.0)
            self.assertLessEqual(res.j_chronic_delivery_index, 1.0)
            self.assertLessEqual(res.complex_step_relative_error, 1e-6)
            self.assertEqual(res.research_boundary, RESEARCH_BOUNDARY)

    def test_curative_ehd_cbf_presets_cross_kahn_kalai(self) -> None:
        cur_72h = evaluate_round204_ehd_taylor_cbf_redosing("ehd_taylor_cbf_adaptive_72h_curative")
        self.assertTrue(cur_72h.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(cur_72h.eta_intact, 0.90)
        self.assertTrue(cur_72h.h1_steric_pass)
        self.assertLessEqual(cur_72h.r_h_post_nm, cur_72h.r_bottleneck_nm)
        self.assertGreaterEqual(cur_72h.t_h1_mucus_transmission, 0.50)
        self.assertTrue(cur_72h.cbf_forward_invariant_safe)
        self.assertGreaterEqual(cur_72h.h_cbf_margin, 0.0)
        self.assertGreaterEqual(cur_72h.delta_t_h, cur_72h.delta_t_cbf_min_h)
        self.assertGreaterEqual(cur_72h.j_chronic_delivery_index, 0.18)

        cur_48h = evaluate_round204_ehd_taylor_cbf_redosing("ehd_taylor_cbf_ultralow_tlr_48h_curative")
        self.assertTrue(cur_48h.crosses_kahn_kalai_threshold)
        self.assertTrue(cur_48h.h1_steric_pass)
        self.assertTrue(cur_48h.cbf_forward_invariant_safe)
        self.assertGreater(cur_48h.j_chronic_delivery_index, cur_72h.j_chronic_delivery_index)

    def test_shear_mucus_and_cbf_failure_presets_fail_kahn_kalai(self) -> None:
        mesh_fail = evaluate_round204_ehd_taylor_cbf_redosing("vibrating_mesh_high_shear_rupture_fail")
        self.assertFalse(mesh_fail.crosses_kahn_kalai_threshold)
        self.assertLess(mesh_fail.eta_intact, 0.10)
        self.assertFalse(mesh_fail.h1_steric_pass)
        self.assertLess(mesh_fail.j_chronic_delivery_index, 0.02)

        trapped = evaluate_round204_ehd_taylor_cbf_redosing("unshielded_280nm_lnp_h1_mucus_trapped_fail")
        self.assertFalse(trapped.crosses_kahn_kalai_threshold)
        self.assertFalse(trapped.h1_steric_pass)
        self.assertLess(trapped.t_h1_mucus_transmission, 1e-4)

        overdose = evaluate_round204_ehd_taylor_cbf_redosing("naive_fixed_12h_overdosing_cbf_violation_fail")
        self.assertFalse(overdose.crosses_kahn_kalai_threshold)
        self.assertFalse(overdose.cbf_forward_invariant_safe)
        self.assertLess(overdose.h_cbf_margin, 0.0)
        self.assertLess(overdose.phi_clearance_retention, 0.01)

    def test_monotonicity_across_ehd_h1_and_cbf_operators(self) -> None:
        eta_ehd, r_ehd, _ = compute_ehd_atomization_and_lnp_integrity(4000.0, 2.8, 1.15, 38.0)
        eta_mesh, r_mesh, _ = compute_ehd_atomization_and_lnp_integrity(180000.0, 2.8, 1.15, 38.0)
        self.assertGreater(float(eta_ehd), float(eta_mesh))
        self.assertLess(float(r_ehd), float(r_mesh))

        t_small = compute_h1_mucus_persistent_homology_transmission(38.0, 70.0, 0.92)
        t_large = compute_h1_mucus_persistent_homology_transmission(140.0, 70.0, 0.92)
        self.assertGreater(float(t_small), float(t_large))

        i_short, h_short, _, m_short = compute_cbf_redosing_invariance(0.22, 18.0)
        i_long, h_long, _, m_long = compute_cbf_redosing_invariance(0.22, 72.0)
        self.assertGreater(float(i_short), float(i_long))
        self.assertLess(float(h_short), float(h_long))
        self.assertGreater(float(m_short), float(m_long))

    def test_complex_step_derivative_accuracy(self) -> None:
        cs_deriv, fd_deriv, rel_err = complex_step_delta_t_derivative(
            72.0, 4000.0, 2.8, 1.15, 38.0, 70.0, 0.92, 0.22, h=1e-30
        )
        self.assertNotEqual(cs_deriv, 0.0)
        self.assertLessEqual(rel_err, 1e-6)

    def test_taylor_angle_and_asymptotic_overflow_safety(self) -> None:
        _, _, f_taylor = compute_ehd_atomization_and_lnp_integrity(4000.0, 2.8, 1.15, 38.0, cone_half_angle_deg=49.3)
        _, _, f_whip = compute_ehd_atomization_and_lnp_integrity(4000.0, 2.8, 1.15, 38.0, cone_half_angle_deg=62.0)
        self.assertGreater(float(f_taylor), float(f_whip))
        with self.assertRaises(ValueError):
            compute_ehd_atomization_and_lnp_integrity(4000.0, 2.8, 1.15, 38.0, cone_half_angle_deg=95.0)
        # Extreme r_h_post_nm and tiny delta_t_h must return 0.0 cleanly without OverflowError
        self.assertEqual(float(compute_h1_mucus_persistent_homology_transmission(6000.0, 70.0, 0.92)), 0.0)
        _, h_tiny, phi_tiny, m_tiny = compute_cbf_redosing_invariance(0.22, 0.05)
        self.assertLess(float(h_tiny), 0.0)
        self.assertEqual(float(phi_tiny), 0.0)
        self.assertAlmostEqual(float(m_tiny), 1.0, places=3)

    def test_falsification_negative_controls_and_unknown_regime(self) -> None:
        self.assertTrue(evaluate_falsification_negative_control(poison_negative_shear=True)["passed"])
        self.assertTrue(evaluate_falsification_negative_control(poison_invalid_gsd=True)["passed"])
        self.assertTrue(evaluate_falsification_negative_control(poison_nonpositive_interval=True)["passed"])
        self.assertTrue(evaluate_falsification_negative_control(poison_excessive_single_dose_impulse=True)["passed"])
        self.assertTrue(evaluate_falsification_negative_control(verify_h1_trapped_and_cbf_overdose_controls=True)["passed"])
        with self.assertRaises(KeyError):
            evaluate_round204_ehd_taylor_cbf_redosing("nonexistent_regime")

    def test_provenance_cross_domain_and_scite_disclosure_metadata(self) -> None:
        self.assertFalse(bool(CROSS_DOMAIN_GOVERNANCE["clinical_extrapolation_permitted"]))
        self.assertIn("inv_cf_anchor", CROSS_DOMAIN_GOVERNANCE)
        self.assertIn("inv_mech_isomorphism", CROSS_DOMAIN_GOVERNANCE)
        self.assertIn("contrasting=1", str(CROSS_DOMAIN_GOVERNANCE["inv_mech_isomorphism"]))
        self.assertGreaterEqual(len(PARAMETER_PROVENANCE), 5)
        for meta in PARAMETER_PROVENANCE.values():
            self.assertIn("value", meta)
            self.assertIn("tag", meta)
            self.assertIn("not stated", str(meta["source"]).lower())
            self.assertIn("bounds", meta)
        self.assertIn("deductive_column_a", TWO_COLUMN_EPISTEMIC_HYGIENE)
        self.assertIn("popperian_column_b", TWO_COLUMN_EPISTEMIC_HYGIENE)
        self.assertGreaterEqual(len(FIVE_LAYER_LADDER), 5)
        self.assertEqual(len(PRIMARY_CITATIONS), 4)


if __name__ == "__main__":
    unittest.main()
