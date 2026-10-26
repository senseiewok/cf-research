"""Tests for the frozen checks benchmark suite (T-0126).

Verifies cryptographic manifest integrity, task schema completeness,
clean primary source anchoring without patient data, and 100% agreement
between expected rule verdicts and benchmark runner evaluations.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys
import unittest

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from benchmark_runner import (
    BENCHMARK_DIR,
    run_benchmark,
    verify_manifest,
)

TASKS_PATH = BENCHMARK_DIR / "tasks.json"
SOURCES_PATH = BENCHMARK_DIR / "sources_fixture.json"
MANIFEST_PATH = BENCHMARK_DIR / "manifest.json"


class ChecksBenchmarkTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        with open(TASKS_PATH, "r", encoding="utf-8") as f:
            cls.tasks_data = json.load(f)
        with open(SOURCES_PATH, "r", encoding="utf-8") as f:
            cls.sources_data = json.load(f)

    def test_manifest_sha256_integrity(self) -> None:
        """All frozen benchmark files must match their SHA-256 digests in manifest.json."""
        self.assertTrue(MANIFEST_PATH.exists(), "manifest.json must exist")
        verified, details = verify_manifest(MANIFEST_PATH)
        self.assertTrue(verified, f"Manifest integrity check failed: {details}")
        for filename, result in details.items():
            self.assertEqual(result["status"], "MATCH", f"{filename} hash mismatch")

    def test_tasks_frozen_census(self) -> None:
        """The benchmark must contain between 20 and 30 frozen tasks (acceptance criteria)."""
        self.assertGreaterEqual(len(self.tasks_data), 20)
        self.assertLessEqual(len(self.tasks_data), 30)
        self.assertEqual(len(self.tasks_data), 25)

    def test_no_patient_data_and_primary_grounding(self) -> None:
        """All tasks must declare no_patient_data=True and cite verified catalog sources with source lines."""
        valid_categories = {
            "evidence_lookup",
            "quote_exactness",
            "numerical_fidelity",
            "widening_and_absence",
            "epistemic_kind_and_basis",
        }
        task_ids = set()

        for task in self.tasks_data:
            tid = task["id"]
            self.assertNotIn(tid, task_ids, f"Duplicate task ID: {tid}")
            task_ids.add(tid)

            self.assertIn(task["category"], valid_categories, f"Invalid category in {tid}")
            self.assertIs(task.get("no_patient_data"), True, f"{tid} must explicitly declare no_patient_data=True")
            self.assertIs(task.get("has_source_line"), True, f"{tid} must confirm has_source_line=True")
            self.assertTrue(task.get("source_citation"), f"{tid} must have primary source citation")
            self.assertTrue(task.get("source_line"), f"{tid} must have primary source line")

            claim = task["claim_object"]
            self.assertIn("id", claim)
            self.assertIn("claim", claim)
            self.assertIn("source", claim)
            source_id = claim["source"]
            self.assertIn(source_id, self.sources_data, f"{tid} cites source '{source_id}' not found in sources_fixture.json")

    def test_benchmark_runner_execution(self) -> None:
        """Benchmark runner must execute clean and match expected verdicts for all 25 tasks."""
        summary = run_benchmark(TASKS_PATH, SOURCES_PATH, MANIFEST_PATH)
        self.assertEqual(summary["total_tasks"], 25)
        self.assertEqual(summary["clean_tasks_count"], 12)
        self.assertEqual(summary["flawed_tasks_count"], 13)

        with_checks = summary["with_lab_checks"]
        self.assertEqual(with_checks["passed"], 12)
        self.assertEqual(with_checks["failed"], 13)
        self.assertEqual(with_checks["flaw_detection_rate"], "100.0%")
        self.assertEqual(with_checks["expected_agreement_rate"], "100.0%")

        rule_counts = with_checks["rule_violation_counts"]
        for r_key in ("R1", "R2", "R3", "R4", "R5"):
            self.assertGreaterEqual(rule_counts[r_key], 1, f"{r_key} must be triggered")

        without_checks = summary["without_lab_checks"]
        self.assertEqual(without_checks["silent_hallucination_count"], 13)
        self.assertEqual(without_checks["silent_hallucination_rate"], "100.0%")
        self.assertEqual(without_checks["flaw_detection_rate"], "0.0%")


if __name__ == "__main__":
    unittest.main()
