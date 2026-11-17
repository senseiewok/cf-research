#!/usr/bin/env python3
"""Unit tests for tools/variant_profile/super_exon_benchmark.py (PRED-0007)."""
from __future__ import annotations

import sys
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from tools.variant_profile import super_exon_benchmark as seb  # noqa: E402


class TestSuperExonBenchmark(unittest.TestCase):
    def test_1_exon_boundary_arithmetic_and_92_allele_coverage(self) -> None:
        rep = seb.evaluate_super_exon_coverage()
        ti1 = rep["regime_coverage"]["intron1_partial_cdna_exons_2_27"]
        dist = rep["regime_coverage"]["distal_super_exon_23_27"]
        full = rep["regime_coverage"]["full_coding_exons_1_27"]

        self.assertEqual(ti1["covered_nonsense_count"], 89)
        self.assertEqual(ti1["excluded_nonsense_count"], 3)
        self.assertEqual([e["allele"] for e in ti1["excluded_examples"]], ["Q2X", "S4X", "S13X"])

        self.assertEqual(dist["codon_range"], "1240..1480")
        self.assertEqual(dist["covered_nonsense_count"], 17)
        self.assertEqual(dist["excluded_nonsense_count"], 75)

        self.assertEqual(full["covered_nonsense_count"], 92)
        self.assertEqual(full["excluded_nonsense_count"], 0)

    def test_2_negative_control_loci_and_falsification_suite(self) -> None:
        by_var = {r["variant"]: r for r in seb.verify_control_panel()}
        self.assertEqual(by_var["Q2X"]["exon"], 1)
        self.assertFalse(by_var["Q2X"]["covered_exons_2_27"])
        self.assertEqual(by_var["R1162X"]["exon"], 22)
        self.assertTrue(by_var["R1162X"]["covered_exons_2_27"])
        self.assertFalse(by_var["R1162X"]["covered_exons_23_27"])
        self.assertEqual(by_var["W1282X"]["exon"], 23)
        self.assertTrue(by_var["W1282X"]["covered_exons_23_27"])
        self.assertEqual(seb.codon_to_exon(560), 12)
        self.assertEqual(seb.codon_to_exon(1239), 22)
        self.assertEqual(seb.codon_to_exon(1240), 23)

        f_rep = seb.run_falsification_suite()
        self.assertTrue(f_rep["falsify_passed"], msg=str(f_rep))
        self.assertEqual(len(f_rep["negative_controls"]), 3)


if __name__ == "__main__":
    unittest.main()
