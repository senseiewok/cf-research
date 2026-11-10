#!/usr/bin/env python3
"""Unit tests for Round 178: Inhaled Lung-Tropic AAV Capsid Directed Evolution Invariant Bridge (T-0852)."""
from __future__ import annotations

import unittest

from tools.discovery.round178_lung_tropic_aav_transduction_bridge import (
    FIVE_LAYER_LADDER,
    PARAMETER_PROVENANCE,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    complex_step_import_derivative,
    compute_apical_binding_fraction,
    compute_mucus_penetration_fractions,
    compute_net_transduction_efficiency,
    compute_nuclear_delivery_fraction,
    evaluate_falsification_negative_control,
    evaluate_round178_lung_tropic_aav_transduction,
)


class Round178LungTropicAAVTransductionTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round178_lung_tropic_aav_transduction(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Invariants failed for {key}")
            self.assertGreaterEqual(res.mucus_steric_fraction, 0.0)
            self.assertLessEqual(res.mucus_steric_fraction, 1.0)
            self.assertGreaterEqual(res.mucus_electrostatic_fraction, 0.0)
            self.assertLessEqual(res.mucus_electrostatic_fraction, 1.0)
            self.assertGreaterEqual(res.net_mucus_penetration, 0.0)
            self.assertLessEqual(res.net_mucus_penetration, 1.0)
            self.assertGreaterEqual(res.apical_binding_fraction, 0.0)
            self.assertLessEqual(res.apical_binding_fraction, 1.0)
            self.assertGreaterEqual(res.nuclear_delivery_fraction, 0.0)
            self.assertLessEqual(res.nuclear_delivery_fraction, 1.0)
            self.assertGreaterEqual(res.net_transduction_efficiency, 0.0)
            self.assertLessEqual(res.net_transduction_efficiency, 1.0)
            self.assertLessEqual(res.complex_step_relative_error, 1e-6)

    def test_engineered_lung_tropic_4d710_capsid(self) -> None:
        res = evaluate_round178_lung_tropic_aav_transduction(
            "engineered_lung_tropic_4d710_capsid"
        )
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res.net_transduction_efficiency, 0.25)
        self.assertGreaterEqual(res.net_mucus_penetration, 0.70)
        self.assertGreaterEqual(res.nuclear_delivery_fraction, 0.65)
        self.assertGreaterEqual(res.apical_binding_fraction, 0.75)

    def test_aav2_6_reciprocal_loop_chimera(self) -> None:
        res = evaluate_round178_lung_tropic_aav_transduction(
            "aav2_6_reciprocal_loop_chimera"
        )
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res.net_transduction_efficiency, 0.15)
        self.assertGreaterEqual(res.net_mucus_penetration, 0.65)

    def test_wild_type_aav2_cationic_entrapped(self) -> None:
        res = evaluate_round178_lung_tropic_aav_transduction(
            "wild_type_aav2_cationic_entrapped"
        )
        self.assertFalse(res.crosses_kahn_kalai_threshold)
        self.assertLess(res.net_transduction_efficiency, 0.05)
        self.assertLess(res.net_mucus_penetration, 0.10)

    def test_suboptimal_high_affinity_hspg_sink(self) -> None:
        res = evaluate_round178_lung_tropic_aav_transduction(
            "suboptimal_high_affinity_hspg_sink"
        )
        self.assertFalse(res.crosses_kahn_kalai_threshold)
        self.assertLess(res.net_transduction_efficiency, 0.10)

    def test_ubiquitination_prone_unmodified_capsid(self) -> None:
        res = evaluate_round178_lung_tropic_aav_transduction(
            "ubiquitination_prone_unmodified_capsid"
        )
        self.assertFalse(res.crosses_kahn_kalai_threshold)
        self.assertLess(res.net_transduction_efficiency, 0.15)
        self.assertLess(res.nuclear_delivery_fraction, 0.30)

    def test_complex_step_import_derivative(self) -> None:
        cs_deriv, fd_deriv, rel_err = complex_step_import_derivative(
            k_import=0.25,
            k_deg=0.08,
        )
        self.assertLessEqual(rel_err, 1e-6)
        self.assertGreater(cs_deriv, 0.0)

    def test_negative_control_falsification(self) -> None:
        neg_pore = evaluate_falsification_negative_control(poison_negative_pore=True)
        self.assertTrue(neg_pore["passed"])
        self.assertIn("Correctly caught", neg_pore["error"])

        neg_imp = evaluate_falsification_negative_control(poison_negative_import=True)
        self.assertTrue(neg_imp["passed"])
        self.assertIn("Correctly caught", neg_imp["error"])

        neg_rad = evaluate_falsification_negative_control(poison_negative_radius=True)
        self.assertTrue(neg_rad["passed"])
        self.assertIn("Correctly caught", neg_rad["error"])

        baseline = evaluate_falsification_negative_control()
        self.assertTrue(baseline["passed"])
        self.assertIsNone(baseline["error"])

    def test_provenance_and_citations(self) -> None:
        self.assertGreaterEqual(len(PRIMARY_CITATIONS), 2)
        self.assertGreaterEqual(len(PARAMETER_PROVENANCE), 5)
        self.assertIn("deductive_theorem", TWO_COLUMN_EPISTEMIC_HYGIENE)
        self.assertIn("biophysical_hypothesis", TWO_COLUMN_EPISTEMIC_HYGIENE)
        self.assertEqual(len(FIVE_LAYER_LADDER), 5)


if __name__ == "__main__":
    unittest.main()
