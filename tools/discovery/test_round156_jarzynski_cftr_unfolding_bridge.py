#!/usr/bin/env python3
"""Unit tests for Round 156: Jarzynski Free Energy Equality & Crooks Fluctuation Theorem Invariant Bridge (T-0815)."""
from __future__ import annotations

import math
import unittest

from tools.discovery.round156_jarzynski_cftr_unfolding_bridge import (
    FIVE_LAYER_LADDER,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    compute_jarzynski_free_energy,
    evaluate_falsification_negative_control,
    evaluate_round156_jarzynski_cftr_unfolding,
)


class Round156JarzynskiTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round156_jarzynski_cftr_unfolding(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Failed for {key}")
            self.assertGreaterEqual(res.dissipated_work_kt, 0.0)  # Second law: W_diss >= 0
            self.assertAlmostEqual(res.crooks_crossing_point_kt, res.jarzynski_free_energy_kt, places=4)
            self.assertGreater(res.open_probability, 0.0)
            self.assertLess(res.open_probability, 1.0)
            self.assertGreaterEqual(res.epithelial_short_circuit_current_ua_cm2, 0.0)

    def test_wild_type_cftr_dimerization(self) -> None:
        wt = evaluate_round156_jarzynski_cftr_unfolding("wild_type_cftr_nbd_dimerization")
        self.assertTrue(wt.crosses_kahn_kalai_threshold)
        self.assertAlmostEqual(wt.jarzynski_free_energy_kt, 0.30, places=2)
        self.assertAlmostEqual(wt.dimerization_free_energy_kt, 0.30, places=2)
        self.assertGreaterEqual(wt.open_probability, 0.40)
        self.assertGreaterEqual(wt.epithelial_short_circuit_current_ua_cm2, 30.0)
        self.assertAlmostEqual(wt.fractional_wild_type_current, 1.0, delta=0.05)

    def test_cf_dimerization_defects(self) -> None:
        df508 = evaluate_round156_jarzynski_cftr_unfolding("cf_deltaf508_misfolded_cavity_destabilization")
        g551d = evaluate_round156_jarzynski_cftr_unfolding("cf_g551d_catalytic_site_jamming")

        for mutant in (df508, g551d):
            self.assertFalse(mutant.crosses_kahn_kalai_threshold)
            self.assertLess(mutant.epithelial_short_circuit_current_ua_cm2, 1.0)
            self.assertLess(mutant.fractional_wild_type_current, 0.05)

        # G551D has high free energy barrier and low open probability
        self.assertGreater(g551d.dimerization_free_energy_kt, 5.0)
        self.assertLess(g551d.open_probability, 0.005)

    def test_kahn_kalai_fluctuation_stability_boundary(self) -> None:
        kk = evaluate_round156_jarzynski_cftr_unfolding("kahn_kalai_fluctuation_stability_boundary")
        self.assertTrue(kk.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(kk.epithelial_short_circuit_current_ua_cm2, 3.75)
        self.assertGreaterEqual(kk.fractional_wild_type_current, 0.118)

    def test_trikafta_allosteric_dimerization_rescue(self) -> None:
        rescued = evaluate_round156_jarzynski_cftr_unfolding("trikafta_allosteric_dimerization_rescue")
        self.assertTrue(rescued.crosses_kahn_kalai_threshold)
        self.assertLessEqual(rescued.dimerization_free_energy_kt, 0.60)
        self.assertGreaterEqual(rescued.open_probability, 0.35)
        self.assertGreaterEqual(rescued.epithelial_short_circuit_current_ua_cm2, 25.0)
        self.assertGreaterEqual(rescued.fractional_wild_type_current, 0.85)

    def test_negative_control_falsification(self) -> None:
        normal = evaluate_falsification_negative_control(poison_negative_temperature=False)
        self.assertTrue(normal["passed"])
        self.assertEqual(normal["rule"], "AH-10")

        poisoned_temp = evaluate_falsification_negative_control(poison_negative_temperature=True)
        self.assertFalse(poisoned_temp["passed"])

        poisoned_var = evaluate_falsification_negative_control(poison_negative_variance=True)
        self.assertFalse(poisoned_var["passed"])

    def test_jarzynski_cumulant_expansion_and_second_law(self) -> None:
        mean_w = 10.0
        var_w = 4.0
        delta_f, w_diss = compute_jarzynski_free_energy(mean_w, var_w)
        self.assertEqual(w_diss, 2.0)
        self.assertEqual(delta_f, 8.0)
        self.assertGreaterEqual(w_diss, 0.0)

        with self.assertRaises(ValueError):
            compute_jarzynski_free_energy(mean_w, -1.0)

    def test_two_column_epistemic_hygiene_and_five_layer_ladder(self) -> None:
        self.assertIn("deductive_theorem", TWO_COLUMN_EPISTEMIC_HYGIENE)
        self.assertIn("biophysical_hypothesis", TWO_COLUMN_EPISTEMIC_HYGIENE)
        self.assertIn("falsification_condition", TWO_COLUMN_EPISTEMIC_HYGIENE["biophysical_hypothesis"])
        self.assertIn("cool_human_fact", FIVE_LAYER_LADDER)
        self.assertIn("everyday_analogy", FIVE_LAYER_LADDER)
        self.assertIn("engineering_lens", FIVE_LAYER_LADDER)
        self.assertIn("formal_theorem", FIVE_LAYER_LADDER)
        self.assertIn("falsifiable_code_and_toy", FIVE_LAYER_LADDER)
        self.assertGreaterEqual(len(PRIMARY_CITATIONS), 3)


if __name__ == "__main__":
    unittest.main()
