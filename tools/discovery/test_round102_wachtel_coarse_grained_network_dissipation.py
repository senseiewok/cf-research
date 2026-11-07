"""Unit tests for Round 102 Wachtel-Rao-Esposito Coarse-Grained Network Dissipation."""

from __future__ import annotations

import unittest

from tools.discovery.round102_wachtel_coarse_grained_network_dissipation import (
    build_round102_summary,
    evaluate_wachtel_network_genotype,
    run_falsification_checks,
)


class Round102WachtelCoarseGrainedTests(unittest.TestCase):
    def test_wt_and_eti_satisfy_coarse_graining_bound(self) -> None:
        wt = evaluate_wachtel_network_genotype("wt_cftr")
        eti = evaluate_wachtel_network_genotype("rescued_cftr_eti")
        self.assertTrue(wt.satisfies_wachtel_coarse_graining_bound)
        self.assertTrue(eti.satisfies_wachtel_coarse_graining_bound)
        self.assertAlmostEqual(wt.macro_entropy_rate_kb_s, 13.3, places=3)
        self.assertAlmostEqual(eti.macro_entropy_rate_kb_s, 11.7, places=3)

    def test_f508del_and_g551d_fail_coarse_graining_bound(self) -> None:
        f508 = evaluate_wachtel_network_genotype("uncorrected_f508del")
        g551d = evaluate_wachtel_network_genotype("uncoupled_g551d")
        self.assertFalse(f508.satisfies_wachtel_coarse_graining_bound)
        self.assertFalse(g551d.satisfies_wachtel_coarse_graining_bound)

    def test_falsification_suite_passes(self) -> None:
        res = run_falsification_checks()
        self.assertTrue(res["passed"])
        summary = build_round102_summary()
        self.assertEqual(summary["round"], 102)
        self.assertEqual(summary["primary_doi"], "10.1103/PhysRevLett.120.240602")


if __name__ == "__main__":
    unittest.main()
