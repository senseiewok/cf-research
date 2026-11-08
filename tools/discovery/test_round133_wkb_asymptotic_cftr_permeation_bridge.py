#!/usr/bin/env python3
"""Unit tests for Round 133: WKB Semiclassical Asymptotic Boundary-Layer Matching & Non-Equilibrium Barrier Penetration Bridge for CFTR Ion Permeation & Gating Pore Energetics."""
from __future__ import annotations

import unittest

from tools.discovery.round133_wkb_asymptotic_cftr_permeation_bridge import (
    FIVE_LAYER_LADDER,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    evaluate_falsification_negative_control,
    evaluate_round133_wkb_asymptotic_cftr_permeation,
)


class Round133WkbAsymptoticCftrPermeationTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round133_wkb_asymptotic_cftr_permeation(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Failed for {key}")
            self.assertLess(res.complex_step_sensitivity, 0.0)
            self.assertGreater(res.barrier_kb_t, 0.0)
            self.assertGreater(res.wkb_transmission, 0.0)
            self.assertLessEqual(res.wkb_transmission, 1.0)
            self.assertGreater(res.gamma_eff_ps, 0.0)
            self.assertGreater(res.current_cl_pa, 0.0)
            self.assertGreater(res.mucosal_flux_fraction, 0.0)

    def test_wild_type_healthy_permeation(self) -> None:
        wt = evaluate_round133_wkb_asymptotic_cftr_permeation("wild_type_healthy")
        self.assertTrue(wt.crosses_kahn_kalai_threshold)
        self.assertLessEqual(wt.barrier_kb_t, 3.50)
        self.assertGreaterEqual(wt.gamma_eff_ps, 8.0)
        self.assertGreaterEqual(wt.current_cl_pa, 0.45)
        self.assertGreaterEqual(wt.mucosal_flux_fraction, 0.90)

    def test_cf_f508del_untreated_conduction_collapse(self) -> None:
        cf = evaluate_round133_wkb_asymptotic_cftr_permeation("cf_f508del_untreated")
        self.assertFalse(cf.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(cf.barrier_kb_t, 8.0)
        self.assertLessEqual(cf.gamma_eff_ps, 0.20)
        self.assertLessEqual(cf.current_cl_pa, 0.02)
        self.assertLessEqual(cf.mucosal_flux_fraction, 0.02)

    def test_kahn_kalai_threshold_rescue(self) -> None:
        kk = evaluate_round133_wkb_asymptotic_cftr_permeation("kahn_kalai_threshold_11_8pct")
        self.assertTrue(kk.crosses_kahn_kalai_threshold)
        self.assertLessEqual(kk.barrier_kb_t, 5.0)
        self.assertGreaterEqual(kk.gamma_eff_ps, 1.0)
        self.assertGreaterEqual(kk.mucosal_flux_fraction, 0.118)

    def test_potentiators_rescue_conductance(self) -> None:
        iva = evaluate_round133_wkb_asymptotic_cftr_permeation("ivacaftor_potentiated")
        tri = evaluate_round133_wkb_asymptotic_cftr_permeation("trikafta_elex_tez_iva_rescued")
        for res in (iva, tri):
            self.assertTrue(res.crosses_kahn_kalai_threshold)
            self.assertLessEqual(res.barrier_kb_t, 3.60)
            self.assertGreaterEqual(res.gamma_eff_ps, 4.0)
            self.assertGreaterEqual(res.mucosal_flux_fraction, 0.45)

    def test_negative_control_falsification(self) -> None:
        normal = evaluate_falsification_negative_control(poison_negative_barrier=False)
        self.assertTrue(normal["passed"])
        self.assertEqual(normal["rule"], "AH-10")
        poisoned = evaluate_falsification_negative_control(poison_negative_barrier=True)
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
