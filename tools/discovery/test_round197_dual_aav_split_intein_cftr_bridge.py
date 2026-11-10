#!/usr/bin/env python3
"""Unit tests for Round 197: Dual-AAV Split-Intein Recombinant Full-Length CFTR cDNA Airway Trans-Splicing Invariant Bridge (T-0884)."""
from __future__ import annotations

import unittest

from tools.discovery.round197_dual_aav_split_intein_cftr_bridge import (
    FIVE_LAYER_LADDER,
    PARAMETER_PROVENANCE,
    PRIMARY_CITATIONS,
    REGIME_PRESETS,
    TWO_COLUMN_EPISTEMIC_HYGIENE,
    complex_step_cftr_derivative,
    compute_reconstituted_cftr_fraction,
    compute_restored_chloride_conductance,
    compute_trans_splicing_efficiency,
    compute_transduction_probability,
    evaluate_falsification_negative_control,
    evaluate_round197_dual_aav_split_intein_cftr,
)


class Round197DualAAVSplitInteinCFTRBridgeTest(unittest.TestCase):
    def test_all_presets_exist_and_satisfy_invariants(self) -> None:
        self.assertGreaterEqual(len(REGIME_PRESETS), 5)
        for key in REGIME_PRESETS:
            res = evaluate_round197_dual_aav_split_intein_cftr(key)
            self.assertTrue(res.falsification_invariant_satisfied, f"Invariants failed for {key}")
            self.assertGreaterEqual(res.moi_1, 0.0)
            self.assertGreaterEqual(res.moi_2, 0.0)
            self.assertGreaterEqual(res.p_entry_1, 0.0)
            self.assertLessEqual(res.p_entry_1, 1.0)
            self.assertGreaterEqual(res.p_entry_2, 0.0)
            self.assertLessEqual(res.p_entry_2, 1.0)
            self.assertGreaterEqual(res.f_co_trans, 0.0)
            self.assertLessEqual(res.f_co_trans, 1.0)
            self.assertGreaterEqual(res.stoichiometric_ratio, 0.0)
            self.assertLessEqual(res.stoichiometric_ratio, 1.0)
            self.assertGreaterEqual(res.eta_splice, 0.0)
            self.assertLessEqual(res.eta_splice, 1.0)
            self.assertGreaterEqual(res.f_reconstituted, 0.0)
            self.assertLessEqual(res.f_reconstituted, 1.0)
            self.assertGreaterEqual(res.toxicity, 0.0)
            self.assertLessEqual(res.toxicity, 1.0)
            self.assertGreaterEqual(res.i_sc_equiv_pct_wt, 0.0)
            self.assertLessEqual(res.i_sc_equiv_pct_wt, 100.0)
            self.assertLessEqual(res.complex_step_relative_error, 1e-3)

    def test_curative_presets_cross_kahn_kalai(self) -> None:
        res_curative = evaluate_round197_dual_aav_split_intein_cftr("dual_aav_equimolar_optimized_curative")
        self.assertTrue(res_curative.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res_curative.f_co_trans, 0.20)
        self.assertGreaterEqual(res_curative.f_reconstituted, 0.15)
        self.assertGreaterEqual(res_curative.i_sc_equiv_pct_wt, 25.0)
        self.assertLessEqual(res_curative.toxicity, 0.05)

        res_super = evaluate_round197_dual_aav_split_intein_cftr("dual_aav_high_dose_super_curative")
        self.assertTrue(res_super.crosses_kahn_kalai_threshold)
        self.assertGreaterEqual(res_super.f_co_trans, 0.50)
        self.assertGreaterEqual(res_super.f_reconstituted, 0.40)
        self.assertGreaterEqual(res_super.i_sc_equiv_pct_wt, 60.0)

    def test_sub_curative_presets_fail_kahn_kalai(self) -> None:
        # Single vector unpaired negative control
        res_unpaired = evaluate_round197_dual_aav_split_intein_cftr("single_vector_unpaired_control")
        self.assertFalse(res_unpaired.crosses_kahn_kalai_threshold)
        self.assertEqual(res_unpaired.f_co_trans, 0.0)
        self.assertEqual(res_unpaired.f_reconstituted, 0.0)
        self.assertEqual(res_unpaired.i_sc_equiv_pct_wt, 0.0)

        # Stoichiometric imbalance
        res_imbalance = evaluate_round197_dual_aav_split_intein_cftr("stoichiometric_imbalance_sub_curative")
        self.assertFalse(res_imbalance.crosses_kahn_kalai_threshold)
        self.assertLess(res_imbalance.f_reconstituted, 0.15)

        # Low-dose sub-therapeutic
        res_low = evaluate_round197_dual_aav_split_intein_cftr("low_dose_sub_therapeutic")
        self.assertFalse(res_low.crosses_kahn_kalai_threshold)
        self.assertLess(res_low.f_co_trans, 0.20)
        self.assertLess(res_low.f_reconstituted, 0.15)

    def test_complex_step_derivative_accuracy(self) -> None:
        cs_deriv, fd_deriv, rel_err = complex_step_cftr_derivative(
            moi_1=40000.0,
            moi_2=40000.0,
            k_splice=0.92,
        )
        self.assertGreater(cs_deriv, 0.0, "Current derivative with respect to MOI must be positive")
        self.assertGreater(fd_deriv, 0.0)
        self.assertLessEqual(rel_err, 1e-3)

    def test_falsification_negative_controls(self) -> None:
        self.assertTrue(evaluate_falsification_negative_control(poison_negative_moi=True)["passed"])
        self.assertTrue(evaluate_falsification_negative_control(poison_negative_splice_rate=True)["passed"])
        self.assertTrue(evaluate_falsification_negative_control(poison_negative_entry_efficiency=True)["passed"])

    def test_monotonicity_invariants(self) -> None:
        # Increasing MOI monotonically increases entry probability
        p1 = float(compute_transduction_probability(10000.0))
        p2 = float(compute_transduction_probability(50000.0))
        self.assertGreater(p2, p1)

        # Increasing splice efficiency increases restored conductance
        isc1 = float(compute_restored_chloride_conductance(0.10))
        isc2 = float(compute_restored_chloride_conductance(0.30))
        self.assertGreater(isc2, isc1)

    def test_provenance_and_epistemic_hygiene_metadata(self) -> None:
        self.assertGreaterEqual(len(PARAMETER_PROVENANCE), 5)
        for param, meta in PARAMETER_PROVENANCE.items():
            self.assertIn("value", meta)
            self.assertIn("tag", meta)
            self.assertIn(meta["tag"], ("empirical", "derived", "phenomenological_heuristic"))
            self.assertIn("source", meta)
            self.assertIn("bounds", meta)

        self.assertIn("deductive_column_a", TWO_COLUMN_EPISTEMIC_HYGIENE)
        self.assertIn("popperian_column_b", TWO_COLUMN_EPISTEMIC_HYGIENE)
        self.assertGreaterEqual(len(FIVE_LAYER_LADDER), 5)
        self.assertGreaterEqual(len(PRIMARY_CITATIONS), 2)


if __name__ == "__main__":
    unittest.main()
