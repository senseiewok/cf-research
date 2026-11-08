#!/usr/bin/env python3
"""Unit tests for Round 123: PhysLean Poisson-Nernst-Planck Electro-Diffusion & Epithelial Donnan Invariant Bridge."""
from __future__ import annotations

import unittest

from tools.discovery.round123_physlean_pnp_electro_diffusion_donnan_bridge import (
    FIVE_LAYER_LADDER,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    evaluate_falsification_negative_control,
    evaluate_physlean_pnp_donnan,
)


class Round123PhysLeanPNPDonnanTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_physlean_pnp_donnan(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Failed for {key}")
            self.assertGreater(res.complex_step_sensitivity, 0.0)

    def test_wild_type_healthy_crosses_pcl_jamming_safety_threshold(self) -> None:
        wt = evaluate_physlean_pnp_donnan("wild_type_healthy")
        self.assertTrue(wt.crosses_pcl_collapse_threshold)
        self.assertGreaterEqual(wt.pcl_height_um, 7.0)
        self.assertGreater(wt.transepithelial_potential_mV, -30.0)

    def test_f508del_untreated_collapses_pcl_and_hyperpolarizes(self) -> None:
        cf = evaluate_physlean_pnp_donnan("f508del_homozygous_untreated")
        self.assertFalse(cf.crosses_pcl_collapse_threshold)
        self.assertLess(cf.pcl_height_um, 4.0)
        self.assertLess(cf.transepithelial_potential_mV, -40.0)

    def test_kahn_kalai_threshold_rescues_pcl_height(self) -> None:
        kk = evaluate_physlean_pnp_donnan("kahn_kalai_threshold_11_8pct")
        self.assertAlmostEqual(kk.cftr_function_pct, 11.8, delta=0.5)
        self.assertGreaterEqual(kk.pcl_height_um, 5.0)

    def test_trikafta_rescued_restores_physiological_potential(self) -> None:
        tri = evaluate_physlean_pnp_donnan("trikafta_rescued")
        self.assertTrue(tri.crosses_pcl_collapse_threshold)
        self.assertGreaterEqual(tri.pcl_height_um, 7.0)
        self.assertGreater(tri.transepithelial_potential_mV, -30.0)

    def test_negative_control_falsification(self) -> None:
        normal = evaluate_falsification_negative_control(poison_inverted_charge=False)
        self.assertTrue(normal["passed"])
        self.assertEqual(normal["rule"], "AH-10")
        poisoned = evaluate_falsification_negative_control(poison_inverted_charge=True)
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
