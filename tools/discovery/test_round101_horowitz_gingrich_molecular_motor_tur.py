"""Unit tests for Round 101 Horowitz-Gingrich Molecular Motor Precision & TUR Bridge."""

from __future__ import annotations

import unittest

from tools.discovery.round101_horowitz_gingrich_molecular_motor_tur import (
    build_round101_summary,
    evaluate_horowitz_gingrich_genotype,
    run_falsification_checks,
)


class Round101HorowitzGingrichTests(unittest.TestCase):
    def test_wt_and_eti_satisfy_motor_tur(self) -> None:
        wt = evaluate_horowitz_gingrich_genotype("wt_cftr")
        eti = evaluate_horowitz_gingrich_genotype("rescued_cftr_eti")
        self.assertTrue(wt.satisfies_horowitz_gingrich_bound)
        self.assertTrue(eti.satisfies_horowitz_gingrich_bound)
        self.assertAlmostEqual(wt.tur_uncertainty_product_kb, 2.2475, places=4)
        self.assertAlmostEqual(eti.tur_uncertainty_product_kb, 2.2008, places=4)

    def test_f508del_and_g551d_fail_motor_tur(self) -> None:
        f508 = evaluate_horowitz_gingrich_genotype("uncorrected_f508del")
        g551d = evaluate_horowitz_gingrich_genotype("uncoupled_g551d")
        self.assertFalse(f508.satisfies_horowitz_gingrich_bound)
        self.assertFalse(g551d.satisfies_horowitz_gingrich_bound)

    def test_falsification_suite_passes(self) -> None:
        res = run_falsification_checks()
        self.assertTrue(res["passed"])
        summary = build_round101_summary()
        self.assertEqual(summary["round"], 101)
        self.assertEqual(summary["primary_doi"], "10.1038/s41567-019-0702-6")


if __name__ == "__main__":
    unittest.main()
