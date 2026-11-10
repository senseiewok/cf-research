#!/usr/bin/env python3
"""Unit tests for Round 201: Hybrid Lipid-Polymer Nanoparticle (hLPN) Mucus Transcytosis Bridge (T-0892)."""
from __future__ import annotations

import unittest

from tools.discovery.round201_hlpn_mucus_transcytosis_bridge import (
    FIVE_LAYER_LADDER,
    PARAMETER_PROVENANCE,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    complex_step_gamma_delivery_derivative,
    compute_mucus_diffusivity_and_penetration,
    compute_net_delivery_index,
    compute_peg_shielding_factor,
    compute_uptake_and_escape,
    evaluate_falsification_negative_control,
    evaluate_round201_hlpn_mucus_transcytosis,
)


class Round201HLPNMucusTranscytosisBridgeTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round201_hlpn_mucus_transcytosis(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Invariants failed for {key}")
            self.assertGreaterEqual(res.s_peg_shielding, 0.0)
            self.assertLessEqual(res.s_peg_shielding, 1.0)
            self.assertGreaterEqual(res.d_eff_over_d0, 0.0)
            self.assertLessEqual(res.d_eff_over_d0, 1.0)
            self.assertGreaterEqual(res.p_mucus_pen, 0.0)
            self.assertLessEqual(res.p_mucus_pen, 1.0)
            self.assertGreaterEqual(res.eta_escape, 0.0)
            self.assertLessEqual(res.eta_escape, 1.0)
            self.assertGreaterEqual(res.phi_net_delivery, 0.0)
            self.assertLessEqual(res.phi_net_delivery, 1.0)
            self.assertLessEqual(res.complex_step_relative_error, 1e-6)

    def test_curative_presets_cross_kahn_kalai(self) -> None:
        res_cur = evaluate_round201_hlpn_mucus_transcytosis("hlpn_dense_peg_optimized_curative")
        self.assertTrue(res_cur.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res_cur.d_eff_over_d0, 0.15)
        self.assertGreaterEqual(res_cur.p_mucus_pen, 0.20)
        self.assertGreaterEqual(res_cur.eta_escape, 0.25)
        self.assertGreaterEqual(res_cur.phi_net_delivery, 0.040)

        res_small = evaluate_round201_hlpn_mucus_transcytosis("hlpn_ultra_small_cleavable_curative")
        self.assertTrue(res_small.crosses_kahn_kalai_threshold)
        self.assertGreater(res_small.p_mucus_pen, res_cur.p_mucus_pen)
        self.assertGreater(res_small.phi_net_delivery, res_cur.phi_net_delivery)

    def test_pathology_and_failure_presets_fail_kahn_kalai(self) -> None:
        res_unshielded = evaluate_round201_hlpn_mucus_transcytosis("unshielded_plga_hydrophobic_adhesion_cf")
        self.assertFalse(res_unshielded.crosses_kahn_kalai_threshold)
        self.assertLess(res_unshielded.d_eff_over_d0, 0.05)
        self.assertLess(res_unshielded.p_mucus_pen, 0.05)

        res_oversized = evaluate_round201_hlpn_mucus_transcytosis("oversized_mesh_excluded_aggregate_cf")
        self.assertFalse(res_oversized.crosses_kahn_kalai_threshold)
        self.assertLess(res_oversized.p_mucus_pen, 0.01)

        res_lipid_def = evaluate_round201_hlpn_mucus_transcytosis("lipid_deficient_endosomal_entrapment")
        self.assertFalse(res_lipid_def.crosses_kahn_kalai_threshold)
        self.assertLess(res_lipid_def.eta_escape, 0.10)
        self.assertLess(res_lipid_def.phi_net_delivery, 0.02)

    def test_monotonicity_and_shielding_invariants(self) -> None:
        s_low = compute_peg_shielding_factor(0.15)
        s_high = compute_peg_shielding_factor(1.50)
        self.assertGreater(float(s_high), float(s_low))

        _, d_low_peg, p_low_peg = compute_mucus_diffusivity_and_penetration(0.10, 80.0, 50.0, 35.0)
        _, d_high_peg, p_high_peg = compute_mucus_diffusivity_and_penetration(1.20, 80.0, 50.0, 35.0)
        self.assertGreater(float(d_high_peg), float(d_low_peg))
        self.assertGreater(float(p_high_peg), float(p_low_peg))

        _, d_low_conc, _ = compute_mucus_diffusivity_and_penetration(0.50, 80.0, 50.0, 20.0)
        _, d_high_conc, _ = compute_mucus_diffusivity_and_penetration(0.50, 80.0, 50.0, 70.0)
        self.assertGreater(float(d_low_conc), float(d_high_conc))

        _, _, p_small_core = compute_mucus_diffusivity_and_penetration(1.0, 60.0, 50.0, 35.0)
        _, _, p_large_core = compute_mucus_diffusivity_and_penetration(1.0, 140.0, 50.0, 35.0)
        self.assertGreater(float(p_small_core), float(p_large_core))

    def test_complex_step_derivative_accuracy(self) -> None:
        cs_deriv, fd_deriv, rel_err = complex_step_gamma_delivery_derivative(0.85, 80.0, 50.0, 35.0, 50.0, 0.35)
        self.assertGreater(cs_deriv, 0.0)
        self.assertGreater(fd_deriv, 0.0)
        self.assertLessEqual(rel_err, 1e-6)

    def test_falsification_negative_controls(self) -> None:
        self.assertTrue(evaluate_falsification_negative_control(poison_negative_gamma=True)["passed"])
        self.assertTrue(evaluate_falsification_negative_control(poison_negative_diameter=True)["passed"])
        self.assertTrue(evaluate_falsification_negative_control(poison_negative_pore=True)["passed"])
        self.assertTrue(evaluate_falsification_negative_control(poison_invalid_lipid=True)["passed"])

    def test_provenance_and_epistemic_hygiene_metadata(self) -> None:
        self.assertGreaterEqual(len(PARAMETER_PROVENANCE), 5)
        for meta in PARAMETER_PROVENANCE.values():
            self.assertIn("value", meta)
            self.assertIn("tag", meta)
            self.assertIn("not stated", str(meta["source"]).lower())
            self.assertIn("bounds", meta)
        self.assertIn("deductive_column_a", TWO_COLUMN_EPISTEMIC_HYGIENE)
        self.assertIn("popperian_column_b", TWO_COLUMN_EPISTEMIC_HYGIENE)
        self.assertGreaterEqual(len(FIVE_LAYER_LADDER), 5)
        self.assertEqual(len(PRIMARY_CITATIONS), 2)


if __name__ == "__main__":
    unittest.main()
