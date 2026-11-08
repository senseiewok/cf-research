#!/usr/bin/env python3
"""Unit tests for Round 129: Doi-Peliti Field Theory & Reactive Mucin Network Disulfide Bond Formation/Cleavage Invariant Bridge."""
from __future__ import annotations

import unittest

from tools.discovery.round129_doi_peliti_mucin_crosslink_kinetics_bridge import (
    FIVE_LAYER_LADDER,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    evaluate_falsification_negative_control,
    evaluate_round129_doi_peliti_mucin_crosslink,
)


class Round129DoiPelitiMucinCrosslinkTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round129_doi_peliti_mucin_crosslink(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Failed for {key}")
            self.assertLess(res.complex_step_sensitivity, 0.0)
            self.assertGreaterEqual(res.plateau_modulus_gp_pa, 0.0)
            self.assertGreater(res.deborah_number, 0.0)
            self.assertGreaterEqual(res.zero_shear_viscosity_pa_s, 0.001)
            self.assertGreaterEqual(res.disulfide_conversion_p_ss, 0.0)
            self.assertLess(res.disulfide_conversion_p_ss, 1.0)

    def test_wild_type_healthy_sol(self) -> None:
        wt = evaluate_round129_doi_peliti_mucin_crosslink("wild_type_healthy")
        self.assertTrue(wt.crosses_kahn_kalai_threshold)
        self.assertFalse(wt.is_percolated_gel)
        self.assertLess(wt.disulfide_conversion_p_ss, wt.critical_percolation_p_c)
        self.assertLessEqual(wt.deborah_number, 0.35)
        self.assertLessEqual(wt.plateau_modulus_gp_pa, 0.20)
        self.assertLessEqual(wt.zero_shear_viscosity_pa_s, 0.02)

    def test_cf_homozygous_dehydrated_gel(self) -> None:
        cf = evaluate_round129_doi_peliti_mucin_crosslink("cf_homozygous_dehydrated")
        self.assertFalse(cf.crosses_kahn_kalai_threshold)
        self.assertTrue(cf.is_percolated_gel)
        self.assertGreater(cf.disulfide_conversion_p_ss, cf.critical_percolation_p_c)
        self.assertGreaterEqual(cf.deborah_number, 5.0)
        self.assertGreaterEqual(cf.plateau_modulus_gp_pa, 6.0)
        self.assertGreaterEqual(cf.zero_shear_viscosity_pa_s, 3.0)

    def test_kahn_kalai_threshold_rescue(self) -> None:
        kk = evaluate_round129_doi_peliti_mucin_crosslink("kahn_kalai_threshold_11_8pct")
        self.assertTrue(kk.crosses_kahn_kalai_threshold)
        self.assertLessEqual(kk.deborah_number, 1.0)
        self.assertLessEqual(kk.plateau_modulus_gp_pa, 1.0)

    def test_reducing_agents_dissolve_gel(self) -> None:
        nac = evaluate_round129_doi_peliti_mucin_crosslink("n_acetylcysteine_nac_cleaved")
        cysteamine = evaluate_round129_doi_peliti_mucin_crosslink("cysteamine_dtt_cleaved")
        for res in (nac, cysteamine):
            self.assertTrue(res.crosses_kahn_kalai_threshold)
            self.assertFalse(res.is_percolated_gel)
            self.assertLess(res.disulfide_conversion_p_ss, res.critical_percolation_p_c)
            self.assertLessEqual(res.deborah_number, 0.65)
            self.assertLessEqual(res.plateau_modulus_gp_pa, 0.50)

    def test_negative_control_falsification(self) -> None:
        normal = evaluate_falsification_negative_control(poison_negative_crosslinking_rate=False)
        self.assertTrue(normal["passed"])
        self.assertEqual(normal["rule"], "AH-10")
        poisoned = evaluate_falsification_negative_control(poison_negative_crosslinking_rate=True)
        self.assertFalse(poisoned["passed"])

    def test_two_column_epistemic_hygiene_and_five_layer_ladder(self) -> None:
        self.assertIn("deductive_theorem", TWO_COLUMN_EPISTEMIC_HYGIENE)
        self.assertIn("biophysical_hypothesis", TWO_COLUMN_EPISTEMIC_HYGIENE)
        self.assertIn("falsification_condition", TWO_COLUMN_EPISTEMIC_HYGIENE["biophysical_hypothesis"])
        for layer in (
            "cool_human_fact",
            "everyday_analogy",
            "engineering_lens",
            "formal_theorem",
            "falsifiable_code_and_toy",
        ):
            self.assertTrue(len(FIVE_LAYER_LADDER.get(layer, "")) > 40)
        self.assertGreaterEqual(len(PRIMARY_CITATIONS), 5)


if __name__ == "__main__":
    unittest.main()
