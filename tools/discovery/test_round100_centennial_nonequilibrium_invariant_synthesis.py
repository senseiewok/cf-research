#!/usr/bin/env python3
"""Unit tests for Round 100 (T-0721): Centennial Unified Non-Equilibrium Invariant Synthesis."""
from __future__ import annotations

import unittest

from round100_centennial_nonequilibrium_invariant_synthesis import (
    PRIMARY_DOI_SEIFERT_RPP_2012,
    compute_centennial_synthesis_metrics,
    evaluate_centennial_nonequilibrium_invariant_synthesis,
    run_falsification,
)


class TestRound100CentennialSynthesis(unittest.TestCase):
    def test_wt_and_eti_satisfy_centennial_joint_invariant(self) -> None:
        res = evaluate_centennial_nonequilibrium_invariant_synthesis()
        self.assertEqual(res["task_id"], "T-0721")
        self.assertEqual(res["round"], 100)
        self.assertEqual(res["primary_doi"], PRIMARY_DOI_SEIFERT_RPP_2012)
        self.assertTrue(res["centennial_synthesis_verified"])
        self.assertTrue(res["joint_invariant_bound_satisfied"])
        self.assertTrue(res["eti_centennial_rescue_verified"])
        self.assertTrue(res["negative_control_detects_violation"])

        wt = res["gating_profiles"]["wt_cftr"]
        self.assertGreaterEqual(wt["joint_invariant_omega_100"], 1.0)
        self.assertTrue(wt["bound_satisfied"])

        eti = res["gating_profiles"]["rescued_cftr_eti"]
        self.assertGreaterEqual(eti["joint_invariant_omega_100"], 1.0)
        self.assertTrue(eti["bound_satisfied"])

    def test_f508del_and_g551d_fail_joint_invariant(self) -> None:
        f508 = compute_centennial_synthesis_metrics("uncorrected_f508del")
        self.assertLess(f508["joint_invariant_omega_100"], 1.0)
        self.assertFalse(f508["bound_satisfied"])

        g551d = compute_centennial_synthesis_metrics("uncoupled_g551d")
        self.assertLess(g551d["joint_invariant_omega_100"], 1.0)
        self.assertFalse(g551d["bound_satisfied"])

    def test_falsification_negative_control(self) -> None:
        neg = run_falsification()
        self.assertTrue(neg["passed"])
        self.assertEqual(neg["verdict"], "PASS")
        self.assertTrue(neg["violation_detected"])

    def test_invalid_inputs_raise_value_error(self) -> None:
        with self.assertRaises(ValueError):
            compute_centennial_synthesis_metrics("unknown_allele")
        with self.assertRaises(ValueError):
            compute_centennial_synthesis_metrics("wt_cftr", tur_margin_override=-0.5)


if __name__ == "__main__":
    unittest.main()
