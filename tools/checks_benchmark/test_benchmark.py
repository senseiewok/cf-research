"""Tests for the frozen checks benchmark suite (T-0126).

Verifies cryptographic manifest integrity, task schema completeness,
clean primary source anchoring without patient data, and 100% agreement
between expected rule verdicts and benchmark runner evaluations.
"""
from __future__ import annotations

import json
from pathlib import Path
import pytest

from benchmark_runner import (
    compute_sha256,
    verify_manifest,
    run_benchmark,
    BENCHMARK_DIR,
)

TASKS_PATH = BENCHMARK_DIR / "tasks.json"
SOURCES_PATH = BENCHMARK_DIR / "sources_fixture.json"
MANIFEST_PATH = BENCHMARK_DIR / "manifest.json"


@pytest.fixture
def tasks_data():
    with open(TASKS_PATH, "r", encoding="utf-8") as f:
        return json.load(f)


@pytest.fixture
def sources_data():
    with open(SOURCES_PATH, "r", encoding="utf-8") as f:
        return json.load(f)


def test_manifest_sha256_integrity():
    """All frozen benchmark files must match their SHA-256 digests in manifest.json."""
    assert MANIFEST_PATH.exists(), "manifest.json must exist"
    verified, details = verify_manifest(MANIFEST_PATH)
    assert verified, f"Manifest integrity check failed: {details}"
    for filename, result in details.items():
        assert result["status"] == "MATCH", f"{filename} hash mismatch"


def test_tasks_frozen_census(tasks_data):
    """The benchmark must contain between 20 and 30 frozen tasks (acceptance criteria)."""
    assert 20 <= len(tasks_data) <= 30, f"Task count {len(tasks_data)} not in [20, 30]"
    assert len(tasks_data) == 25, "Expected exactly 25 tasks in frozen census"


def test_no_patient_data_and_primary_grounding(tasks_data, sources_data):
    """All tasks must declare no_patient_data=True and cite verified catalog sources with source lines."""
    valid_categories = {
        "evidence_lookup",
        "quote_exactness",
        "numerical_fidelity",
        "widening_and_absence",
        "epistemic_kind_and_basis",
    }
    task_ids = set()

    for task in tasks_data:
        tid = task["id"]
        assert tid not in task_ids, f"Duplicate task ID: {tid}"
        task_ids.add(tid)

        assert task["category"] in valid_categories, f"Invalid category in {tid}"
        assert task.get("no_patient_data") is True, f"{tid} must explicitly declare no_patient_data=True"
        assert task.get("has_source_line") is True, f"{tid} must confirm has_source_line=True"
        assert task.get("source_citation"), f"{tid} must have primary source citation"
        assert task.get("source_line"), f"{tid} must have primary source line"

        # Check claim object
        claim = task["claim_object"]
        assert "id" in claim and "claim" in claim and "source" in claim
        source_id = claim["source"]
        assert source_id in sources_data, f"{tid} cites source '{source_id}' not found in sources_fixture.json"


def test_benchmark_runner_execution():
    """Benchmark runner must execute clean and match expected verdicts for all 25 tasks."""
    summary = run_benchmark(TASKS_PATH, SOURCES_PATH, MANIFEST_PATH)
    assert summary["total_tasks"] == 25
    assert summary["clean_tasks_count"] == 12
    assert summary["flawed_tasks_count"] == 13

    # With lab checks: 100% flaw detection, 0% false positives
    with_checks = summary["with_lab_checks"]
    assert with_checks["passed"] == 12
    assert with_checks["failed"] == 13
    assert with_checks["flaw_detection_rate"] == "100.0%"
    assert with_checks["expected_agreement_rate"] == "100.0%"

    # Verify rule violation distribution covers R1 through R5
    rule_counts = with_checks["rule_violation_counts"]
    assert rule_counts["R1"] >= 1, "R1 must be triggered"
    assert rule_counts["R2"] >= 1, "R2 must be triggered"
    assert rule_counts["R3"] >= 1, "R3 must be triggered"
    assert rule_counts["R4"] >= 1, "R4 must be triggered"
    assert rule_counts["R5"] >= 1, "R5 must be triggered"

    # Without lab checks: silent hallucinations pass
    without_checks = summary["without_lab_checks"]
    assert without_checks["silent_hallucination_count"] == 13
    assert without_checks["silent_hallucination_rate"] == "100.0%"
    assert without_checks["flaw_detection_rate"] == "0.0%"
