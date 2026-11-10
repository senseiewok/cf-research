#!/usr/bin/env python3
"""Unit tests for Round 200: Bicellular Tight Junction Claudin-2/Claudin-4 Paracellular Ion Permselectivity Bridge (T-0889)."""
from __future__ import annotations

import unittest

from tools.discovery.round200_tight_junction_claudin_flux_bridge import (
    FIVE_LAYER_LADDER,
    PARAMETER_PROVENANCE,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    bernoulli_ghk_kernel,
    complex_step_ph_transcellular_derivative,
    compute_claudin_permeabilities_and_resistance,
    compute_ghk_paracellular_fluxes_and_asl_state,
    evaluate_falsification_negative_control,
    evaluate_round200_tight_junction_claudin_flux,
)


class Round200TightJunctionClaudinFluxBridgeTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round200_tight_junction_claudin_flux(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Invariants failed for {key}")
            self.assertGreater(res.p_na_1e6_cm_s, 0.0)
            self.assertGreater(res.p_hco3_1e6_cm_s, 0.0)
            self.assertGreater(res.r_tj_ohm_cm2, 0.0)
            self.assertGreater(res.c_hco3_asl_mm, 0.0)
            self.assertLessEqual(res.complex_step_relative_error, 1e-6)

    def test_curative_presets_cross_kahn_kalai(self) -> None:
        res_cur = evaluate_round200_tight_junction_claudin_flux("claudin4_intact_barrier_cftr_rescued_curative")
        self.assertTrue(res_cur.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res_cur.r_tj_ohm_cm2, 500.0)
        self.assertGreaterEqual(res_cur.ph_asl, 7.25)
        self.assertLessEqual(abs(res_cur.j_na_para_nmol_cm2_h), 450.0)
        self.assertGreaterEqual(res_cur.h_asl_um, 5.20)

        res_enh = evaluate_round200_tight_junction_claudin_flux("claudin4_enhanced_barrier_stabilized_curative")
        self.assertTrue(res_enh.crosses_kahn_kalai_threshold)
        self.assertGreater(res_enh.r_tj_ohm_cm2, res_cur.r_tj_ohm_cm2)
        self.assertGreater(res_enh.h_asl_um, res_cur.h_asl_um)

    def test_cf_and_barrier_disruption_presets_fail_kahn_kalai(self) -> None:
        # Untreated CF: paracellular HCO3- influx > 0 tempers pH near ~6.87, but fails curative alkalinization
        res_cf = evaluate_round200_tight_junction_claudin_flux("untreated_cf_hyperpolarized_paracellular_tempering")
        self.assertFalse(res_cf.crosses_kahn_kalai_threshold)
        self.assertGreater(res_cf.j_hco3_para_nmol_cm2_h, 0.0, "Paracellular HCO3- must flow into acidic CF ASL")
        self.assertLess(res_cf.ph_asl, 7.10)
        self.assertGreater(res_cf.ph_asl, 6.50)

        # Cytokine disruption collapses R_TJ < 150 Ohm*cm2
        res_cyt = evaluate_round200_tight_junction_claudin_flux("cytokine_disrupted_leaky_tight_junction")
        self.assertFalse(res_cyt.crosses_kahn_kalai_threshold)
        self.assertLess(res_cyt.r_tj_ohm_cm2, 150.0)
        self.assertLess(res_cyt.h_asl_um, 3.0)

        # Claudin-2 overexpression creates cation-selective shunt (P_Na / P_HCO3 > 5)
        res_cldn2 = evaluate_round200_tight_junction_claudin_flux("claudin2_overexpression_cation_shunt_leak")
        self.assertFalse(res_cldn2.crosses_kahn_kalai_threshold)
        self.assertGreater(res_cldn2.cation_anion_selectivity, 5.0)
        self.assertLess(res_cldn2.r_tj_ohm_cm2, 200.0)

    def test_claudin_selectivity_and_zero_voltage_removable_singularity(self) -> None:
        pna_low, phc_low, _ = compute_claudin_permeabilities_and_resistance(0.2, 1.0, 0.0)
        pna_high, phc_high, _ = compute_claudin_permeabilities_and_resistance(2.0, 1.0, 0.0)
        self.assertGreater(float(pna_high), float(pna_low))
        self.assertGreater(float(pna_high / phc_high), float(pna_low / phc_low))

        _, _, rtj_low_c4 = compute_claudin_permeabilities_and_resistance(0.3, 0.2, 0.0)
        _, _, rtj_high_c4 = compute_claudin_permeabilities_and_resistance(0.3, 1.8, 0.0)
        self.assertGreater(float(rtj_high_c4), float(rtj_low_c4))

        # Exact removable singularity at V_te = 0.0 mV and complex V_te = 0 + 1e-30j
        self.assertAlmostEqual(float(bernoulli_ghk_kernel(0.0)), 1.0, places=14)
        chco3_0, _, jna_0, ph_0, _ = compute_ghk_paracellular_fluxes_and_asl_state(1.3, 1.0, 800.0, 50.0, 0.0)
        self.assertGreater(float(chco3_0), 24.0)
        self.assertAlmostEqual(float(jna_0), 3.6 * 1.3 * (140.0 - 105.0), places=10)
        self.assertGreater(float(ph_0), 7.3)

    def test_complex_step_derivative_accuracy(self) -> None:
        cs_deriv, fd_deriv, rel_err = complex_step_ph_transcellular_derivative(0.25, 1.35, 0.0, 56.0, -8.5)
        self.assertGreater(cs_deriv, 0.0)
        self.assertGreater(fd_deriv, 0.0)
        self.assertLessEqual(rel_err, 1e-6)

    def test_falsification_negative_controls(self) -> None:
        self.assertTrue(evaluate_falsification_negative_control(poison_negative_claudin=True)["passed"])
        self.assertTrue(evaluate_falsification_negative_control(poison_negative_cytokine=True)["passed"])
        self.assertTrue(evaluate_falsification_negative_control(poison_negative_transcellular_flux=True)["passed"])

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
        self.assertIn("DOI: 10.1113/jp280120", PRIMARY_CITATIONS)
        self.assertIn("DOI: 10.1016/j.semcdb.2015.04.009", PRIMARY_CITATIONS)
        self.assertIn("DOI: 10.1091/mbc.e02-03-0134", PRIMARY_CITATIONS)


if __name__ == "__main__":
    unittest.main()
