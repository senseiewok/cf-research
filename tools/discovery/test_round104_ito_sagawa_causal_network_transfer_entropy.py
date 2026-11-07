"""Unit tests for Round 104 Ito-Sagawa Causal Network Transfer Entropy Bridge."""

from __future__ import annotations

import unittest

from tools.discovery.round104_ito_sagawa_causal_network_transfer_entropy import (
    build_round104_summary,
    evaluate_ito_sagawa_causal_genotype,
    run_falsification,
    run_falsification_checks,
)


class Round104ItoSagawaCausalNetworkTests(unittest.TestCase):
    def test_wt_and_eti_satisfy_causal_second_law(self) -> None:
        wt = evaluate_ito_sagawa_causal_genotype("wt_cftr")
        eti = evaluate_ito_sagawa_causal_genotype("rescued_cftr_eti")
        self.assertTrue(wt.satisfies_ito_sagawa_causal_bound)
        self.assertTrue(eti.satisfies_ito_sagawa_causal_bound)
        self.assertAlmostEqual(wt.generalized_second_law_margin_kb, 0.60, places=4)
        self.assertAlmostEqual(eti.generalized_second_law_margin_kb, 0.56, places=4)

    def test_f508del_and_g551d_fail_causal_bound(self) -> None:
        f508 = evaluate_ito_sagawa_causal_genotype("uncorrected_f508del")
        g551d = evaluate_ito_sagawa_causal_genotype("uncoupled_g551d")
        self.assertFalse(f508.satisfies_ito_sagawa_causal_bound)
        self.assertFalse(g551d.satisfies_ito_sagawa_causal_bound)

    def test_falsification_suite_passes(self) -> None:
        res = run_falsification_checks()
        self.assertTrue(res["passed"])
        neg = run_falsification()
        self.assertTrue(neg["passed"])
        summary = build_round104_summary()
        self.assertEqual(summary["round"], 104)
        self.assertEqual(summary["primary_doi"], "10.1103/PhysRevLett.111.180603")


if __name__ == "__main__":
    unittest.main()
