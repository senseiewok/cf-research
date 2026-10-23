#!/usr/bin/env python3
"""Unit tests for tools/benchmark/cf_evidence_benchmark.py."""
from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import cf_evidence_benchmark as ceb  # noqa: E402


class TestCfEvidenceBenchmark(unittest.TestCase):
    def test_build_benchmark_and_critic_accuracy(self):
        benchmark, dpo_pairs = ceb.build_benchmark_from_repo()
        summary = benchmark["summary"]
        self.assertGreaterEqual(summary["total_items"], 80)
        self.assertEqual(summary["supported_items"], summary["unsupported_items"])
        self.assertEqual(summary["dpo_pairs_count"], summary["supported_items"])

        metrics = summary["claim_support_critic_metrics"]
        self.assertEqual(metrics["accuracy"], 1.0)
        self.assertEqual(metrics["f1_on_violations"], 1.0)
        self.assertEqual(metrics["fn_missed_violations"], 0)
        self.assertEqual(metrics["fp_false_alarms"], 0)

    def test_all_six_violation_categories_detected(self):
        critic = ceb.ClaimSupportCritic()
        base = {
            "id": "test-claim-1",
            "claim": "In 2024, 92.3 percent of the total CF population in the Registry was eligible by CFTR genotype.",
            "quote": "Overall, 92.3 percent of the total CF population in the Registry was eligible by CFTR genotype for at least one CFTR modulator as of 2024.",
            "kind": "observed",
            "scope": "The CFTR Modulator Therapies chapter of the CFFPR 2024 Annual Data Report was read.",
        }
        pos_res = critic.evaluate(base["claim"], base["quote"], base["scope"], base["kind"])
        self.assertEqual(pos_res["verdict"], "PASS")

        expected_cats = {
            "R3_NUMERIC_DISTORTION",
            "R4_UNSCOPED_WIDENING",
            "R5_SILENT_EXTRAPOLATION",
            "AH3_CROSS_REGISTRY_POOLING",
            "SAFE_USE_ELIGIBILITY_DIRECTIVE",
            "AH6_FABRICATED_CITATION",
        }
        seen_cats = set()
        for i in range(6):
            neg = ceb.mutate_claim(base, i)
            seen_cats.add(neg["violation_category"])
            res = critic.evaluate(neg["claim"], neg["quote"], neg["scope"], neg["kind"])
            self.assertEqual(res["verdict"], "FAIL", f"Failed to flag {neg['violation_category']}")
            self.assertEqual(res["label"], "UNSUPPORTED")
        self.assertEqual(seen_cats, expected_cats)

    def test_cli_check_and_write(self):
        with tempfile.TemporaryDirectory() as tmp:
            out_b = Path(tmp) / "bench.json"
            out_p = Path(tmp) / "pairs.jsonl"
            rc = ceb.main(["--out-benchmark", str(out_b), "--out-pairs", str(out_p)])
            self.assertEqual(rc, 0)
            self.assertTrue(out_b.is_file())
            self.assertTrue(out_p.is_file())
            lines = out_p.read_text(encoding="utf-8").strip().splitlines()
            self.assertGreaterEqual(len(lines), 40)
            first_pair = json.loads(lines[0])
            self.assertIn("prompt", first_pair)
            self.assertIn("chosen", first_pair)
            self.assertIn("rejected", first_pair)


if __name__ == "__main__":
    unittest.main()
