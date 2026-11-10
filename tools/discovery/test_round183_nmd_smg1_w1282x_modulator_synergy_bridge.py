#!/usr/bin/env python3
"""Unit tests for Round 183: NMD SMG1 Kinase Inhibition & W1282X Modulator Synergy Bridge (T-0860)."""
from __future__ import annotations

import unittest

from tools.discovery.round183_nmd_smg1_w1282x_modulator_synergy_bridge import (
    FIVE_LAYER_LADDER,
    PARAMETER_PROVENANCE,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    complex_step_smg1_derivative,
    compute_asl_height_and_mcc_velocity,
    compute_cftr_mrna_ratio,
    compute_delta_isc_ratio,
    compute_folding_efficiency,
    compute_open_probability,
    compute_smg1_inhibition,
    evaluate_falsification_negative_control,
    evaluate_round183_nmd_smg1_modulator_synergy,
)


class Round183NmdSmg1ModulatorSynergyTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round183_nmd_smg1_modulator_synergy(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Invariants failed for {key}")
            self.assertGreaterEqual(res.smg1_inhibition_fraction, 0.0)
            self.assertLessEqual(res.smg1_inhibition_fraction, 1.0)
            self.assertGreaterEqual(res.cftr_mrna_ratio, 0.0)
            self.assertLessEqual(res.cftr_mrna_ratio, 1.0)
            self.assertGreaterEqual(res.folding_maturation_efficiency, 0.0)
            self.assertLessEqual(res.folding_maturation_efficiency, 1.0)
            self.assertGreaterEqual(res.channel_open_probability, 0.0)
            self.assertLessEqual(res.channel_open_probability, 1.0)
            self.assertGreaterEqual(res.delta_isc_ratio, 0.0)
            self.assertGreaterEqual(res.asl_height_um, 0.0)
            self.assertLessEqual(res.asl_height_um, 15.0)
            self.assertGreaterEqual(res.mcc_velocity_mm_min, 0.0)
            self.assertLessEqual(res.complex_step_relative_error, 1e-6)

    def test_quad_modal_synergy_crosses_kahn_kalai(self) -> None:
        res = evaluate_round183_nmd_smg1_modulator_synergy("quad_modal_smg1i_elx_tez_iva_synergy")
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res.delta_isc_ratio, 0.150)
        self.assertGreaterEqual(res.asl_height_um, 7.0)
        self.assertGreaterEqual(res.mcc_velocity_mm_min, 1.0)
        self.assertGreaterEqual(res.cftr_mrna_ratio, 0.70)
        self.assertGreaterEqual(res.folding_maturation_efficiency, 0.25)
        self.assertGreaterEqual(res.channel_open_probability, 0.25)

    def test_sub_curative_regimes_fail_kahn_kalai(self) -> None:
        # Tri-modal without ELX
        tri = evaluate_round183_nmd_smg1_modulator_synergy("tri_modal_smg1i_tez_iva")
        self.assertFalse(tri.crosses_kahn_kalai_threshold)
        self.assertLess(tri.asl_height_um, 7.0)
        self.assertLess(tri.mcc_velocity_mm_min, 1.0)

        # ETI alone without SMG1i
        eti = evaluate_round183_nmd_smg1_modulator_synergy("eti_alone_without_smg1i")
        self.assertFalse(eti.crosses_kahn_kalai_threshold)
        self.assertLess(eti.asl_height_um, 7.0)
        self.assertLess(eti.mcc_velocity_mm_min, 1.0)

        # SMG1i monotherapy without modulators
        smg1 = evaluate_round183_nmd_smg1_modulator_synergy("smg1i_monotherapy_without_modulators")
        self.assertFalse(smg1.crosses_kahn_kalai_threshold)
        self.assertLess(smg1.delta_isc_ratio, 0.05)
        self.assertEqual(smg1.mcc_velocity_mm_min, 0.0)

        # Untreated baseline
        untreated = evaluate_round183_nmd_smg1_modulator_synergy("untreated_w1282x_baseline")
        self.assertFalse(untreated.crosses_kahn_kalai_threshold)
        self.assertLess(untreated.delta_isc_ratio, 0.01)
        self.assertEqual(untreated.mcc_velocity_mm_min, 0.0)

    def test_complex_step_derivative(self) -> None:
        cs_deriv, fd_deriv, rel_err = complex_step_smg1_derivative(c_smg1i_um=0.30)
        self.assertLessEqual(rel_err, 1e-6)
        self.assertGreater(cs_deriv, 0.0)

    def test_negative_control_falsification(self) -> None:
        neg_c = evaluate_falsification_negative_control(poison_negative_c_smg1i=True)
        self.assertTrue(neg_c["passed"])
        self.assertIn("Correctly caught", neg_c["error"])

        neg_ic50 = evaluate_falsification_negative_control(poison_negative_ic50=True)
        self.assertTrue(neg_ic50["passed"])
        self.assertIn("Correctly caught", neg_ic50["error"])

        neg_finh = evaluate_falsification_negative_control(poison_unphysical_f_inh=True)
        self.assertTrue(neg_finh["passed"])
        self.assertIn("Correctly caught", neg_finh["error"])

        neg_mod = evaluate_falsification_negative_control(poison_negative_modulator_conc=True)
        self.assertTrue(neg_mod["passed"])
        self.assertIn("Correctly caught", neg_mod["error"])

        baseline = evaluate_falsification_negative_control()
        self.assertTrue(baseline["passed"])
        self.assertIsNone(baseline["error"])

    def test_provenance_and_citations(self) -> None:
        self.assertGreaterEqual(len(PRIMARY_CITATIONS), 3)
        self.assertIn("DOI: 10.1172/jci.insight.139983", PRIMARY_CITATIONS)
        self.assertIn("DOI: 10.1016/j.jcf.2020.05.004", PRIMARY_CITATIONS)
        self.assertIn("DOI: 10.1073/pnas.1215982110", PRIMARY_CITATIONS)

        self.assertGreaterEqual(len(PARAMETER_PROVENANCE), 5)
        for key, meta in PARAMETER_PROVENANCE.items():
            self.assertIn("tag", meta)
            self.assertIn("source", meta)
            self.assertIn("bounds", meta)

        self.assertIn("deductive_theorem", TWO_COLUMN_EPISTEMIC_HYGIENE)
        self.assertIn("biophysical_hypothesis", TWO_COLUMN_EPISTEMIC_HYGIENE)
        self.assertIn("cool_human_fact", FIVE_LAYER_LADDER)
        self.assertIn("falsifiable_code_and_toy", FIVE_LAYER_LADDER)


if __name__ == "__main__":
    unittest.main()
