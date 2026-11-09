#!/usr/bin/env python3
"""Unit tests for Round 164: Flory-Stockmayer Branched Sol-Gel Phase Transition Invariant Bridge (T-0829)."""
from __future__ import annotations

import unittest

from tools.discovery.round164_flory_stockmayer_gelation_bridge import (
    FIVE_LAYER_LADDER,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    complex_step_gel_fraction_derivative,
    compute_critical_percolation_threshold,
    compute_equilibrium_shear_modulus,
    compute_mucociliary_velocity,
    compute_mucolytic_conversion,
    compute_sol_gel_fractions,
    compute_weight_average_functionality,
    evaluate_falsification_negative_control,
    evaluate_round164_flory_stockmayer_gelation,
    solve_extinction_probability,
)


class Round164FloryStockmayerGelationTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round164_flory_stockmayer_gelation(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Failed for {key}")
            self.assertGreater(res.weight_average_functionality, 2.0)
            self.assertAlmostEqual(res.sol_fraction + res.gel_fraction, 1.0, places=10)

    def test_healthy_mucus_sol_subcritical(self) -> None:
        res = evaluate_round164_flory_stockmayer_gelation("healthy_mucus_sol_subcritical")
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertEqual(res.gel_fraction, 0.0)
        self.assertEqual(res.sol_fraction, 1.0)
        self.assertEqual(res.equilibrium_shear_modulus_pa, 0.0)
        self.assertGreaterEqual(res.mucociliary_velocity_ums, 30.0)

    def test_cf_hyperconcentrated_mucus_gel(self) -> None:
        res = evaluate_round164_flory_stockmayer_gelation("cf_hyperconcentrated_mucus_gel")
        self.assertFalse(res.crosses_kahn_kalai_threshold)
        self.assertGreater(res.gel_fraction, 0.80)
        self.assertLess(res.sol_fraction, 0.20)
        self.assertGreater(res.equilibrium_shear_modulus_pa, 1.2)
        self.assertLess(res.mucociliary_velocity_ums, 1.0)

    def test_kahn_kalai_gelation_boundary(self) -> None:
        res = evaluate_round164_flory_stockmayer_gelation("kahn_kalai_gelation_boundary")
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res.mucociliary_velocity_ums, 5.9)

    def test_nac_mucolytic_liquefaction_rescue(self) -> None:
        res = evaluate_round164_flory_stockmayer_gelation("nac_mucolytic_liquefaction_rescue")
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertLess(res.crosslink_conversion, res.critical_threshold_pc)
        self.assertEqual(res.gel_fraction, 0.0)
        self.assertEqual(res.equilibrium_shear_modulus_pa, 0.0)
        self.assertGreaterEqual(res.mucociliary_velocity_ums, 20.0)

    def test_dtt_complete_disulfide_reduction(self) -> None:
        res = evaluate_round164_flory_stockmayer_gelation("dtt_complete_disulfide_reduction")
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertLess(res.crosslink_conversion, 0.05)
        self.assertEqual(res.gel_fraction, 0.0)
        self.assertGreaterEqual(res.mucociliary_velocity_ums, 35.0)

    def test_shear_modulus_bifurcation_at_percolation_threshold(self) -> None:
        f_w = 5.0
        p_c = compute_critical_percolation_threshold(f_w)  # 0.25
        g_sub = compute_equilibrium_shear_modulus(p_c - 0.05, f_w)
        g_at = compute_equilibrium_shear_modulus(p_c, f_w)
        g_post = compute_equilibrium_shear_modulus(p_c + 0.10, f_w)

        self.assertEqual(g_sub, 0.0)
        self.assertEqual(g_at, 0.0)
        self.assertGreater(g_post, 0.0)

    def test_sol_gel_mass_conservation_and_monotonicity(self) -> None:
        f_w = 4.5
        p_c = compute_critical_percolation_threshold(f_w)
        for p_val in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]:
            w_sol, w_gel = compute_sol_gel_fractions(p_val, f_w)
            self.assertAlmostEqual(w_sol + w_gel, 1.0, places=10)
            if p_val <= p_c:
                self.assertEqual(w_gel, 0.0)
            else:
                self.assertGreater(w_gel, 0.0)

    def test_complex_step_gel_fraction_sensitivity(self) -> None:
        cs_deriv, fd_deriv, rel_err = complex_step_gel_fraction_derivative(
            p=0.45,
            f_w=5.0,
        )
        self.assertLessEqual(rel_err, 1e-3)
        self.assertGreater(cs_deriv, 0.0)

    def test_negative_control_falsification(self) -> None:
        neg1 = evaluate_falsification_negative_control(poison_negative_functionality=True)
        self.assertFalse(neg1["passed"])
        self.assertIn("poisoned", neg1["error"])

        neg2 = evaluate_falsification_negative_control(poison_negative_conversion=True)
        self.assertFalse(neg2["passed"])

        neg3 = evaluate_falsification_negative_control(poison_excessive_conversion=True)
        self.assertFalse(neg3["passed"])

        neg4 = evaluate_falsification_negative_control(poison_negative_g_max=True)
        self.assertFalse(neg4["passed"])

        clean = evaluate_falsification_negative_control()
        self.assertTrue(clean["passed"])

    def test_epistemic_hygiene_and_ladder_structure(self) -> None:
        self.assertIn("deductive_theorem", TWO_COLUMN_EPISTEMIC_HYGIENE)
        self.assertIn("biophysical_hypothesis", TWO_COLUMN_EPISTEMIC_HYGIENE)
        self.assertGreaterEqual(len(FIVE_LAYER_LADDER), 5)
        self.assertGreaterEqual(len(PRIMARY_CITATIONS), 4)


if __name__ == "__main__":
    unittest.main()
