#!/usr/bin/env python3
"""Unit tests for Round 99: Van Vu-Hasegawa Geometrical Irreversibility Bounds (T-0718)."""
from __future__ import annotations

import unittest

from tools.discovery.round99_van_vu_geometrical_irreversibility import (
    compute_van_vu_geometrical_irreversibility_metrics,
    evaluate_van_vu_geometrical_irreversibility,
    run_falsification,
)


class Round99VanVuGeometricalIrreversibilityTests(unittest.TestCase):
    def test_four_cftr_genotypes_and_geometrical_bounds(self) -> None:
        wt = compute_van_vu_geometrical_irreversibility_metrics("wt_cftr")
        self.assertTrue(wt["bound_satisfied"])
        self.assertGreaterEqual(wt["geometrical_irreversibility_ratio"], 1.0)
        self.assertEqual(wt["gating_regime"], "GEOMETRICAL_IRREVERSIBILITY_OPTIMAL")

        eti = compute_van_vu_geometrical_irreversibility_metrics("rescued_cftr_eti")
        self.assertTrue(eti["bound_satisfied"])
        self.assertGreaterEqual(eti["geometrical_irreversibility_ratio"], 1.0)
        self.assertGreaterEqual(eti["geometrical_coupling_efficiency"], 0.90)

        f508 = compute_van_vu_geometrical_irreversibility_metrics("uncorrected_f508del")
        self.assertFalse(f508["bound_satisfied"])
        self.assertLess(f508["geometrical_irreversibility_ratio"], 1.0)

        g551d = compute_van_vu_geometrical_irreversibility_metrics("uncoupled_g551d")
        self.assertFalse(g551d["bound_satisfied"])
        self.assertLess(g551d["geometrical_irreversibility_ratio"], 1.0)

    def test_task_id_and_overall_verification(self) -> None:
        res = evaluate_van_vu_geometrical_irreversibility()
        self.assertEqual(res["task_id"], "T-0718")
        self.assertEqual(res["round"], 99)
        self.assertTrue(res["van_vu_geometrical_irreversibility_verified"])
        self.assertTrue(res["geometrical_irreversibility_bound_satisfied"])
        self.assertTrue(res["eti_geometrical_convergence_verified"])
        self.assertTrue(res["negative_control_detects_violation"])

    def test_falsification_and_invalid_inputs(self) -> None:
        neg = run_falsification()
        self.assertEqual(neg["rule"], "AH-10")
        self.assertTrue(neg["passed"])
        self.assertTrue(neg["violation_detected"])

        with self.assertRaises(ValueError):
            compute_van_vu_geometrical_irreversibility_metrics("non_existent_genotype")

        with self.assertRaises(ValueError):
            compute_van_vu_geometrical_irreversibility_metrics("wt_cftr", entropy_production_override=-5.0)


if __name__ == "__main__":
    unittest.main()
