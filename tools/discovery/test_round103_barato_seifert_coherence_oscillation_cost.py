"""Unit tests for Round 103 Barato-Seifert Coherent Biochemical Oscillation Cost Bridge."""

from __future__ import annotations

import unittest

from tools.discovery.round103_barato_seifert_coherence_oscillation_cost import (
    build_round103_summary,
    evaluate_barato_seifert_coherence_genotype,
    run_falsification,
    run_falsification_checks,
)


class Round103BaratoSeifertCoherenceTests(unittest.TestCase):
    def test_wt_and_eti_satisfy_coherence_bound(self) -> None:
        wt = evaluate_barato_seifert_coherence_genotype("wt_cftr")
        eti = evaluate_barato_seifert_coherence_genotype("rescued_cftr_eti")
        self.assertTrue(wt.satisfies_barato_seifert_coherence_bound)
        self.assertTrue(eti.satisfies_barato_seifert_coherence_bound)
        self.assertLessEqual(wt.coherent_oscillations_n, wt.max_coherent_oscillations_bound)
        self.assertLessEqual(eti.coherent_oscillations_n, eti.max_coherent_oscillations_bound)

    def test_f508del_and_g551d_fail_coherence_bound(self) -> None:
        f508 = evaluate_barato_seifert_coherence_genotype("uncorrected_f508del")
        g551d = evaluate_barato_seifert_coherence_genotype("uncoupled_g551d")
        self.assertFalse(f508.satisfies_barato_seifert_coherence_bound)
        self.assertFalse(g551d.satisfies_barato_seifert_coherence_bound)

    def test_falsification_suite_passes(self) -> None:
        res = run_falsification_checks()
        self.assertTrue(res["passed"])
        neg = run_falsification()
        self.assertTrue(neg["passed"])
        summary = build_round103_summary()
        self.assertEqual(summary["round"], 103)
        self.assertEqual(summary["primary_doi"], "10.1103/PhysRevX.6.041053")


if __name__ == "__main__":
    unittest.main()
