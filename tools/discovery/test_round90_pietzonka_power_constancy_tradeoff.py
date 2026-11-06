"""Unit tests for Round 90 (T-0697): Pietzonka-Barato-Seifert General Trade-Off Between Power, Efficiency, and Constancy.

Primary Biophysical Grounding & Evidence:
  - DOI:10.1103/PhysRevE.93.052145 (PMID:27300918 — Pietzonka P, Barato AC, Seifert U: Universal bounds on current fluctuations, Phys. Rev. E 93, 052145, 2016)
  - DOI:10.1126/science.ade2216 (Fiedorczuk & Chen: Molecular structures reveal synergistic rescue of dF508 CFTR by Trikafta modulators, Science 2022)
  - PDB 6MSM, 5UAK, 8EIQ
"""
from __future__ import annotations

from pathlib import Path
import sys
import unittest

_REPO_ROOT = Path(__file__).resolve().parent.parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from tools.discovery.round90_pietzonka_power_constancy_tradeoff import (
    compute_pietzonka_genotype_metrics,
    evaluate_pietzonka_power_constancy_tradeoff,
    run_falsification,
)


class Round90PietzonkaPowerConstancyTests(unittest.TestCase):
    """Test suite verifying Round 90 Pietzonka-Barato-Seifert power-efficiency-constancy bounds."""

    def setUp(self) -> None:
        self.report = evaluate_pietzonka_power_constancy_tradeoff()

    def test_task_id_and_overall_verification(self) -> None:
        self.assertEqual(self.report["task_id"], "T-0697")
        self.assertEqual(self.report["round"], 90)
        self.assertEqual(self.report["primary_doi"], "10.1103/PhysRevE.93.052145")
        self.assertEqual(self.report["primary_pmid"], "27300918")
        self.assertTrue(self.report["pietzonka_power_constancy_tradeoff_verified"])
        self.assertTrue(self.report["power_efficiency_constancy_bound_verified"])
        self.assertTrue(self.report["eti_power_constancy_restoration_verified"])
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_four_cftr_genotypes_and_bounds(self) -> None:
        profiles = self.report["gating_profiles"]
        for gid in ("wt_cftr", "rescued_cftr_eti", "uncorrected_f508del", "uncoupled_g551d"):
            self.assertIn(gid, profiles)

        wt = profiles["wt_cftr"]
        eti = profiles["rescued_cftr_eti"]
        f508 = profiles["uncorrected_f508del"]
        g551d = profiles["uncoupled_g551d"]

        self.assertTrue(wt["bound_satisfied"])
        self.assertAlmostEqual(wt["observed_output_power_kbt_s"], 12.8, places=4)
        self.assertAlmostEqual(wt["pietzonka_power_upper_bound_kbt_s"], 16.59375, places=4)
        self.assertLessEqual(wt["tradeoff_saturation_ratio"], 1.0)

        self.assertTrue(eti["bound_satisfied"])
        self.assertAlmostEqual(eti["observed_output_power_kbt_s"], 11.28, places=4)
        self.assertAlmostEqual(eti["pietzonka_power_upper_bound_kbt_s"], 17.045333, places=4)
        self.assertLessEqual(eti["tradeoff_saturation_ratio"], 1.0)

        self.assertFalse(f508["bound_satisfied"])
        self.assertLess(f508["observed_output_power_kbt_s"], 1.0)
        self.assertFalse(g551d["bound_satisfied"])
        self.assertLess(g551d["observed_output_power_kbt_s"], 1.0)

    def test_falsification_and_invalid_inputs(self) -> None:
        fals = run_falsification()
        self.assertEqual(fals["rule"], "AH-10")
        self.assertTrue(fals["passed"])
        self.assertEqual(fals["verdict"], "PASS")
        self.assertGreater(fals["reversible_limit_ratio"], 1.0)

        with self.assertRaises(ValueError):
            compute_pietzonka_genotype_metrics("unknown_allele")
        with self.assertRaises(ValueError):
            compute_pietzonka_genotype_metrics("wt_cftr", efficiency_override=1.05)


if __name__ == "__main__":
    unittest.main()
