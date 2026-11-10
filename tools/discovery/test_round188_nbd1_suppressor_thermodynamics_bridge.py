#!/usr/bin/env python3
"""Unit tests for Round 188: NBD1 Suppressor Thermodynamics Bridge (T-0869)."""
from __future__ import annotations

import unittest

from tools.discovery.round188_nbd1_suppressor_thermodynamics_bridge import (
    FIVE_LAYER_LADDER,
    PARAMETER_PROVENANCE,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    complex_step_asl_derivative,
    compute_asl_height_um,
    compute_band_c_maturation_fraction,
    compute_gibbs_helmholtz_folding_free_energy_kj_mol,
    compute_mcc_velocity_mm_min,
    compute_nbd1_folded_fraction,
    compute_transepithelial_chloride_current_pct_wt,
    evaluate_falsification_negative_control,
    evaluate_round188_nbd1_suppressor_thermodynamics,
)


class Round188NBD1SuppressorThermodynamicsTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round188_nbd1_suppressor_thermodynamics(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Invariants failed for {key}")
            self.assertGreaterEqual(res.effective_tm_c, 35.0)
            self.assertLessEqual(res.effective_tm_c, 60.0)
            self.assertGreaterEqual(res.folded_fraction, 0.0)
            self.assertLessEqual(res.folded_fraction, 1.0)
            self.assertGreaterEqual(res.cell_surface_band_c_fraction, 0.0)
            self.assertLessEqual(res.cell_surface_band_c_fraction, 1.0)
            self.assertGreaterEqual(res.transepithelial_current_pct_wt, 0.0)
            self.assertLessEqual(res.transepithelial_current_pct_wt, 100.0)
            self.assertGreaterEqual(res.airway_surface_liquid_height_um, 3.5)
            self.assertLessEqual(res.airway_surface_liquid_height_um, 8.5)
            self.assertGreaterEqual(res.mucociliary_transport_velocity_mm_min, 0.0)
            self.assertLessEqual(res.complex_step_relative_error, 1e-6)

    def test_triple_suppressor_crosses_kahn_kalai(self) -> None:
        res = evaluate_round188_nbd1_suppressor_thermodynamics("triple_suppressor_delta_ri_r1070w_v510a")
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res.delta_tm_c, 6.0)
        self.assertLessEqual(res.delta_delta_g_supp_kj_mol, -15.0)
        self.assertGreaterEqual(res.cell_surface_band_c_fraction, 0.25)
        self.assertGreaterEqual(res.transepithelial_current_pct_wt, 7.0)
        self.assertGreaterEqual(res.airway_surface_liquid_height_um, 7.0)
        self.assertGreaterEqual(res.mucociliary_transport_velocity_mm_min, 1.0)

    def test_dual_suppressor_crosses_kahn_kalai(self) -> None:
        res = evaluate_round188_nbd1_suppressor_thermodynamics("dual_suppressor_delta_ri_r1070w")
        self.assertTrue(res.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res.delta_tm_c, 6.0)
        self.assertLessEqual(res.delta_delta_g_supp_kj_mol, -15.0)
        self.assertGreaterEqual(res.airway_surface_liquid_height_um, 7.0)
        self.assertGreaterEqual(res.mucociliary_transport_velocity_mm_min, 2.0)

    def test_sub_curative_presets_fail_kahn_kalai(self) -> None:
        # Isolated R1070W: Delta T_m (+4.5 C) < +6.0 C, Delta Delta G (-12.5) > -15.0
        r1070w = evaluate_round188_nbd1_suppressor_thermodynamics("isolated_r1070w_single_suppressor")
        self.assertFalse(r1070w.crosses_kahn_kalai_threshold)
        self.assertLess(r1070w.delta_tm_c, 6.0)

        # Isolated V510A: weak stabilizer
        v510a = evaluate_round188_nbd1_suppressor_thermodynamics("isolated_v510a_weak_suppressor")
        self.assertFalse(v510a.crosses_kahn_kalai_threshold)
        self.assertLess(v510a.delta_tm_c, 6.0)

        # Untreated F508del baseline
        untreated = evaluate_round188_nbd1_suppressor_thermodynamics("untreated_f508del_baseline")
        self.assertFalse(untreated.crosses_kahn_kalai_threshold)
        self.assertAlmostEqual(untreated.cell_surface_band_c_fraction, 0.020, places=3)
        self.assertAlmostEqual(untreated.airway_surface_liquid_height_um, 4.50, delta=0.2)
        self.assertEqual(untreated.mucociliary_transport_velocity_mm_min, 0.0)

    def test_complex_step_derivative_accuracy(self) -> None:
        cs_deriv, fd_deriv, rel_err = complex_step_asl_derivative(band_c_fraction=0.35)
        self.assertGreater(cs_deriv, 0.0)
        self.assertGreater(fd_deriv, 0.0)
        self.assertLessEqual(rel_err, 1e-6)

    def test_poisoned_negative_controls_fail_deterministically(self) -> None:
        res_tm = evaluate_falsification_negative_control(poison_negative_tm=True)
        self.assertTrue(res_tm["passed"])
        self.assertIn("Correctly caught", res_tm["error"])

        res_temp = evaluate_falsification_negative_control(poison_negative_temp=True)
        self.assertTrue(res_temp["passed"])
        self.assertIn("Correctly caught", res_temp["error"])

        res_traffick = evaluate_falsification_negative_control(poison_negative_trafficking=True)
        self.assertTrue(res_traffick["passed"])
        self.assertIn("Correctly caught", res_traffick["error"])

        clean = evaluate_falsification_negative_control()
        self.assertTrue(clean["passed"])
        self.assertIsNone(clean["error"])

    def test_monotonic_tm_stabilization(self) -> None:
        # Increasing T_m makes Delta G_fold more negative (stabilizing folded state)
        dg_low = compute_gibbs_helmholtz_folding_free_energy_kj_mol(tm_k=312.65)
        dg_high = compute_gibbs_helmholtz_folding_free_energy_kj_mol(tm_k=325.65)
        self.assertGreater(dg_low, dg_high)

    def test_provenance_and_two_column_epistemic_hygiene_metadata(self) -> None:
        self.assertGreaterEqual(len(PRIMARY_CITATIONS), 4)
        self.assertGreaterEqual(len(PARAMETER_PROVENANCE), 5)
        self.assertIn("deductive_theorem", TWO_COLUMN_EPISTEMIC_HYGIENE)
        self.assertIn("biophysical_hypothesis", TWO_COLUMN_EPISTEMIC_HYGIENE)
        self.assertIn("statement", TWO_COLUMN_EPISTEMIC_HYGIENE["biophysical_hypothesis"])
        self.assertIn("falsification_condition", TWO_COLUMN_EPISTEMIC_HYGIENE["biophysical_hypothesis"])
        for layer in ["cool_human_fact", "everyday_analogy", "engineering_lens", "formal_theorem", "falsifiable_code_and_toy"]:
            self.assertIn(layer, FIVE_LAYER_LADDER)


if __name__ == "__main__":
    unittest.main()
