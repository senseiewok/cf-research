#!/usr/bin/env python3
"""Unit tests for Round 203: Chemically Modified Suppressor tRNA & Microfluidic LNP Nonsense Rescue Bridge (T-0896)."""
from __future__ import annotations

import unittest

from tools.discovery.round203_chem_mod_suptrna_lnp_bridge import (
    FIVE_LAYER_LADDER,
    PARAMETER_PROVENANCE,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    RESEARCH_BOUNDARY,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    complex_step_f_mod_derivative,
    compute_endosomal_stability_and_tlr78,
    compute_nbd_folding_and_anion_efflux,
    compute_synthetase_charging_and_decoding,
    evaluate_falsification_negative_control,
    evaluate_round203_chem_mod_suptrna_lnp,
)


class Round203ChemModSupTRNALNPBridgeTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round203_chem_mod_suptrna_lnp(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Invariants failed for {key}")
            self.assertGreaterEqual(res.t_half_trna_h, 0.99)
            self.assertGreaterEqual(res.b_trna_bioavail, 0.0)
            self.assertLessEqual(res.b_trna_bioavail, 1.0)
            self.assertGreaterEqual(res.i_tlr78_index, 0.0)
            self.assertLessEqual(res.i_tlr78_index, 1.0)
            self.assertGreaterEqual(res.eta_charge, 0.0)
            self.assertLessEqual(res.eta_charge, 1.0)
            self.assertGreaterEqual(res.alpha_suptrna_readthrough, 0.0)
            self.assertLessEqual(res.alpha_suptrna_readthrough, 1.0)
            self.assertGreaterEqual(res.p_fold_nbd, 0.0)
            self.assertLessEqual(res.p_fold_nbd, 1.0)
            self.assertGreaterEqual(res.i_sc_equiv_pct_wt, 0.0)
            self.assertLessEqual(res.i_sc_equiv_pct_wt, 100.0)
            self.assertLessEqual(res.complex_step_relative_error, 1e-6)
            self.assertEqual(res.research_boundary, RESEARCH_BOUNDARY)

    def test_curative_presets_cross_kahn_kalai(self) -> None:
        w1282x = evaluate_round203_chem_mod_suptrna_lnp("w1282x_chem_mod_suptrna_lnp_curative")
        self.assertTrue(w1282x.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(w1282x.t_half_trna_h, 4.0)
        self.assertLessEqual(w1282x.i_tlr78_index, 0.20)
        self.assertGreaterEqual(w1282x.eta_charge, 0.65)
        self.assertGreaterEqual(w1282x.i_sc_equiv_pct_wt, 15.0)

        g542x = evaluate_round203_chem_mod_suptrna_lnp("g542x_chem_mod_suptrna_lnp_curative")
        self.assertTrue(g542x.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(g542x.t_half_trna_h, 4.0)
        self.assertLessEqual(g542x.i_tlr78_index, 0.20)
        self.assertGreaterEqual(g542x.eta_charge, 0.65)
        self.assertGreaterEqual(g542x.i_sc_equiv_pct_wt, 15.0)

    def test_failure_and_negative_control_presets_fail_kahn_kalai(self) -> None:
        naked = evaluate_round203_chem_mod_suptrna_lnp("unmodified_naked_suptrna_rnase_tlr78_fail")
        self.assertFalse(naked.crosses_kahn_kalai_threshold)
        self.assertLessEqual(naked.t_half_trna_h, 1.05)
        self.assertGreater(naked.i_tlr78_index, 0.80)
        self.assertLess(naked.i_sc_equiv_pct_wt, 2.0)

        blocked = evaluate_round203_chem_mod_suptrna_lnp("anticodon_blocked_overmodified_inactivated_fail")
        self.assertFalse(blocked.crosses_kahn_kalai_threshold)
        self.assertLess(blocked.eta_charge, 0.02)
        self.assertLess(blocked.i_sc_equiv_pct_wt, 2.0)

        misfolded = evaluate_round203_chem_mod_suptrna_lnp("w1282x_nearcognate_nbd2_misfolding_unrescued_fail")
        self.assertFalse(misfolded.crosses_kahn_kalai_threshold)
        self.assertLess(misfolded.p_fold_nbd, 0.10)
        self.assertLess(misfolded.i_sc_equiv_pct_wt, 10.0)

    def test_monotonicity_and_thermodynamic_invariants(self) -> None:
        t_low, b_low, i_low = compute_endosomal_stability_and_tlr78(0.05, 0.90, 0.35)
        t_high, b_high, i_high = compute_endosomal_stability_and_tlr78(0.45, 0.90, 0.35)
        self.assertGreater(float(t_high), float(t_low))
        self.assertGreater(float(b_high), float(b_low))
        self.assertLess(float(i_high), float(i_low))

        eta_clean, alpha_clean = compute_synthetase_charging_and_decoding(0.40, 0.0, 0.95, 0.24)
        eta_blocked, alpha_blocked = compute_synthetase_charging_and_decoding(0.40, 0.65, 0.95, 0.24)
        self.assertGreater(float(eta_clean), float(eta_blocked))
        self.assertGreater(float(alpha_clean), float(alpha_blocked))

        p_cognate, _, _, i_cognate = compute_nbd_folding_and_anion_efflux(0.25, 0.30, 0.0, 0.0)
        p_nearcog, _, _, i_nearcog = compute_nbd_folding_and_anion_efflux(0.25, 0.30, 2.50, 0.0)
        self.assertGreater(float(p_cognate), float(p_nearcog))
        self.assertGreater(float(i_cognate), float(i_nearcog))

    def test_complex_step_derivative_accuracy(self) -> None:
        cs_deriv, fd_deriv, rel_err = complex_step_f_mod_derivative(
            0.35, 0.0, 0.94, 0.38, 0.95, 0.32, 0.0, 0.80, h=1e-30
        )
        self.assertGreater(cs_deriv, 0.0)
        self.assertGreater(fd_deriv, 0.0)
        self.assertLessEqual(rel_err, 1e-6)
        # Boundary f_mod = 0.0 and f_mod = 1.0 must remain accurate without central-FD clamping bias
        cs0, fd0, err0 = complex_step_f_mod_derivative(
            0.0, 0.0, 0.94, 0.38, 0.95, 0.32, 0.0, 0.80, h=1e-30
        )
        self.assertGreater(cs0, 0.0)
        self.assertGreater(fd0, 0.0)
        self.assertLessEqual(err0, 1e-6)

    def test_bounds_and_asymptotic_overflow_safety(self) -> None:
        with self.assertRaises(ValueError):
            compute_synthetase_charging_and_decoding(0.42, 0.0, 0.95, -0.1)
        with self.assertRaises(ValueError):
            compute_nbd_folding_and_anion_efflux(0.32, -0.05, 0.0, 0.0)
        # Extreme ddg_insert_kcal must asymptotically approach 0 without OverflowError
        p_fold, _, _, i_sc = compute_nbd_folding_and_anion_efflux(0.32, 0.25, 1000.0, 0.0)
        self.assertEqual(float(p_fold), 0.0)
        self.assertEqual(float(i_sc), 0.0)

    def test_falsification_negative_controls_and_unknown_regime(self) -> None:
        self.assertTrue(evaluate_falsification_negative_control(poison_negative_f_mod=True)["passed"])
        self.assertTrue(evaluate_falsification_negative_control(poison_exceeded_core_block=True)["passed"])
        self.assertTrue(evaluate_falsification_negative_control(poison_invalid_encapsulation=True)["passed"])
        self.assertTrue(evaluate_falsification_negative_control(poison_negative_ddg=True)["passed"])
        self.assertTrue(evaluate_falsification_negative_control(verify_naked_and_anticodon_controls=True)["passed"])
        with self.assertRaises(KeyError):
            evaluate_round203_chem_mod_suptrna_lnp("nonexistent_regime")

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
        self.assertEqual(len(PRIMARY_CITATIONS), 3)


if __name__ == "__main__":
    unittest.main()
