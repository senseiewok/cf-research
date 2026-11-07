#!/usr/bin/env python3
"""Unit tests for Round 98: Van Vu-Hasegawa Arbitrary Drive TUR (T-0717)."""
from __future__ import annotations

import unittest

from tools.discovery.round98_van_vu_arbitrary_drive_tur import (
    compute_van_vu_arbitrary_drive_tur_metrics,
    evaluate_van_vu_arbitrary_drive_tur,
    run_falsification,
)


class Round98VanVuArbitraryDriveTURTests(unittest.TestCase):
    def test_four_cftr_genotypes_and_arbitrary_tur_bounds(self) -> None:
        wt = compute_van_vu_arbitrary_drive_tur_metrics("wt_cftr")
        self.assertTrue(wt["bound_satisfied"])
        self.assertGreaterEqual(wt["arbitrary_tur_product"], 2.0)
        self.assertEqual(wt["gating_regime"], "ARBITRARY_DRIVE_STABILIZED")

        eti = compute_van_vu_arbitrary_drive_tur_metrics("rescued_cftr_eti")
        self.assertTrue(eti["bound_satisfied"])
        self.assertGreaterEqual(eti["arbitrary_tur_product"], 2.0)
        self.assertGreaterEqual(eti["driving_fidelity_efficiency"], 0.90)

        f508 = compute_van_vu_arbitrary_drive_tur_metrics("uncorrected_f508del")
        self.assertFalse(f508["bound_satisfied"])
        self.assertLess(f508["arbitrary_tur_product"], 2.0)

        g551d = compute_van_vu_arbitrary_drive_tur_metrics("uncoupled_g551d")
        self.assertFalse(g551d["bound_satisfied"])
        self.assertLess(g551d["arbitrary_tur_product"], 2.0)

    def test_task_id_and_overall_verification(self) -> None:
        res = evaluate_van_vu_arbitrary_drive_tur()
        self.assertEqual(res["task_id"], "T-0717")
        self.assertEqual(res["round"], 98)
        self.assertTrue(res["van_vu_arbitrary_drive_tur_verified"])
        self.assertTrue(res["arbitrary_tur_bound_satisfied"])
        self.assertTrue(res["eti_arbitrary_drive_rescue_verified"])
        self.assertTrue(res["negative_control_detects_violation"])

    def test_falsification_and_invalid_inputs(self) -> None:
        neg = run_falsification()
        self.assertEqual(neg["rule"], "AH-10")
        self.assertTrue(neg["passed"])
        self.assertTrue(neg["violation_detected"])

        with self.assertRaises(ValueError):
            compute_van_vu_arbitrary_drive_tur_metrics("non_existent_genotype")

        with self.assertRaises(ValueError):
            compute_van_vu_arbitrary_drive_tur_metrics("wt_cftr", entropy_production_override=-5.0)


if __name__ == "__main__":
    unittest.main()
