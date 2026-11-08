#!/usr/bin/env python3
"""Unit tests for Round 138: Poisson-Nernst-Planck Electrodiffusive Field & Space-Charge Asymmetry Bridge."""
from __future__ import annotations

import unittest

from tools.discovery.round138_pnp_cftr_electrodiffusion_bridge import (
    FIVE_LAYER_LADDER,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    evaluate_falsification_negative_control,
    evaluate_round138_pnp_cftr_electrodiffusion,
)


class Round138PnpCftrElectrodiffusionTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round138_pnp_cftr_electrodiffusion(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Failed for {key}")
            self.assertGreater(res.complex_step_conductance_sensitivity, 0.0)
            self.assertGreater(res.fixed_charge_density_c_m3, 0.0)
            self.assertGreater(res.potential_well_depth_mv, 0.0)
            self.assertGreater(res.chloride_conductance_ps, 0.0)
            self.assertGreaterEqual(res.anion_selectivity_ratio, 1.0)

    def test_wild_type_open_cftr_permeation(self) -> None:
        wt = evaluate_round138_pnp_cftr_electrodiffusion("wild_type_open_cftr")
        self.assertTrue(wt.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(wt.fixed_charge_density_c_m3, 2.0e8)
        self.assertGreaterEqual(wt.chloride_conductance_ps, 8.0)
        self.assertGreaterEqual(wt.anion_selectivity_ratio, 12.0)
        self.assertLessEqual(wt.reversal_potential_mv, -32.0)

    def test_cf_r334w_selectivity_collapse(self) -> None:
        cf = evaluate_round138_pnp_cftr_electrodiffusion("cf_r334w_selectivity_collapse")
        self.assertFalse(cf.crosses_kahn_kalai_threshold)
        self.assertLessEqual(cf.fixed_charge_density_c_m3, 0.50e8)
        self.assertLessEqual(cf.chloride_conductance_ps, 0.50)
        self.assertLessEqual(cf.anion_selectivity_ratio, 2.5)
        self.assertGreaterEqual(cf.reversal_potential_mv, -10.0)

    def test_kahn_kalai_threshold_rescue(self) -> None:
        kk = evaluate_round138_pnp_cftr_electrodiffusion("kahn_kalai_conductance_threshold")
        self.assertTrue(kk.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(kk.chloride_conductance_ps, 1.00)
        self.assertGreaterEqual(kk.anion_selectivity_ratio, 3.0)

    def test_ivacaftor_and_trikafta_restoration(self) -> None:
        iva = evaluate_round138_pnp_cftr_electrodiffusion("ivacaftor_potentiated_rescue")
        tri = evaluate_round138_pnp_cftr_electrodiffusion("trikafta_synergistic_restoration")
        for res in (iva, tri):
            self.assertTrue(res.crosses_kahn_kalai_threshold)
            self.assertGreaterEqual(res.chloride_conductance_ps, 7.5)
            self.assertGreaterEqual(res.anion_selectivity_ratio, 10.0)
            self.assertLessEqual(res.reversal_potential_mv, -30.0)

    def test_negative_control_falsification(self) -> None:
        normal = evaluate_falsification_negative_control(poison_negative_charge=False)
        self.assertTrue(normal["passed"])
        self.assertEqual(normal["rule"], "AH-10")
        poisoned = evaluate_falsification_negative_control(poison_negative_charge=True)
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
