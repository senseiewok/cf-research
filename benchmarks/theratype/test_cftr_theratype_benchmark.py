#!/usr/bin/env python3
"""Unit tests for Open CFTR Theratype Benchmark (T-0925)."""
from __future__ import annotations

from pathlib import Path
import sys
import unittest

BENCH_DIR = Path(__file__).resolve().parent
if str(BENCH_DIR) not in sys.path:
    sys.path.insert(0, str(BENCH_DIR))

import cftr_theratype_benchmark as ctb  # noqa: E402


class TestCftrTheratypeBenchmark(unittest.TestCase):
    def test_fixed_split_passes_leakage_audit(self) -> None:
        audit = ctb.audit_split_leakage(ctb.OPEN_CFTR_THERATYPE_BENCHMARK)
        self.assertTrue(audit["leakage_audit_passed"])
        self.assertEqual(audit["codon_position_overlap"], [])
        self.assertEqual(audit["variant_family_overlap"], [])
        self.assertEqual(audit["train_count"], 6)
        self.assertEqual(audit["test_count"], 10)

    def test_evo2_and_alphagenome_cover_class_i_and_splice_beyond_alphamissense(self) -> None:
        am = ctb.evaluate_model_coverage_and_auroc(ctb.OPEN_CFTR_THERATYPE_BENCHMARK, "alphamissense")
        evo2 = ctb.evaluate_model_coverage_and_auroc(ctb.OPEN_CFTR_THERATYPE_BENCHMARK, "evo2")
        ag = ctb.evaluate_model_coverage_and_auroc(ctb.OPEN_CFTR_THERATYPE_BENCHMARK, "alphagenome_avi")
        self.assertEqual(am["test_coverage"], 0.1)
        self.assertEqual(evo2["test_coverage"], 1.0)
        self.assertEqual(ag["test_coverage"], 1.0)
        self.assertGreaterEqual(evo2["disruption_auroc"], 0.90)
        self.assertGreaterEqual(ag["disruption_auroc"], 0.90)

    def test_falsification_controls_pass(self) -> None:
        falsify = ctb.evaluate_falsification_controls()
        self.assertEqual(falsify["status"], "PASS")
        self.assertEqual(falsify["controls_passed"], 3)

    def test_cli_main_passes(self) -> None:
        self.assertEqual(ctb.main(["--check", "--falsify"]), 0)


if __name__ == "__main__":
    unittest.main()
