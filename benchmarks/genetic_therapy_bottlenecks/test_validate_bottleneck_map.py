#!/usr/bin/env python3
"""Unit and planted-defect tests for benchmarks/genetic_therapy_bottlenecks/validate_bottleneck_map.py."""
from __future__ import annotations

import copy
import json
import subprocess
import sys
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from validate_bottleneck_map import (  # noqa: E402
    DEFAULT_MAP_PATH,
    DEFAULT_REPO_ROOT,
    STAGE_KEYS,
    validate_bottleneck_map_data,
)

VALIDATOR_SCRIPT = HERE / "validate_bottleneck_map.py"


class TestBottleneckMapValidation(unittest.TestCase):
    def setUp(self) -> None:
        self.base_payload = json.loads(DEFAULT_MAP_PATH.read_text(encoding="utf-8"))

    def test_shipped_map_passes_strictly(self) -> None:
        rep = validate_bottleneck_map_data(self.base_payload, repo_root=DEFAULT_REPO_ROOT)
        self.assertEqual(rep["status"], "PASS", msg=rep["errors"])
        self.assertEqual(rep["studies_validated"], 8)
        self.assertEqual(rep["errors"], [])
        self.assertEqual(len(rep["studies"]), 8)

    def test_planted_defect_1_leaderboard_score_key_rejected(self) -> None:
        bad = copy.deepcopy(self.base_payload)
        bad["studies"][0]["efficacy_score"] = 0.99
        rep = validate_bottleneck_map_data(bad, repo_root=DEFAULT_REPO_ROOT)
        self.assertEqual(rep["status"], "FAIL")
        self.assertTrue(any("efficacy_score" in e for e in rep["errors"]))

    def test_planted_defect_2_registry_hallucination_s5_rejected(self) -> None:
        bad = copy.deepcopy(self.base_payload)
        nct = next(s for s in bad["studies"] if s["source_id"] == "NCT05248230")
        nct["stages"]["S5_clinical_benefit_or_lung_function"] = {
            "status": "measured_in_system",
            "measurement_scale": "clinical_safety_registry",
            "evidence_type": "snapshot_field_pointer",
            "snapshot_nct_id": "NCT05248230",
            "snapshot_field": "overall_status",
            "expected_substring": "RECRUITING",
        }
        rep = validate_bottleneck_map_data(bad, repo_root=DEFAULT_REPO_ROOT)
        self.assertEqual(rep["status"], "FAIL")
        self.assertTrue(any("Registry hallucination guard" in e for e in rep["errors"]))

    def test_planted_defect_3_preclinical_in_vivo_overclaim_rejected(self) -> None:
        bad = copy.deepcopy(self.base_payload)
        preclin = next(s for s in bad["studies"] if s["study_id"] == "prime-editing-organoids-2024")
        preclin["stages"]["S3_in_vivo_airway_delivery"] = dict(
            preclin["stages"]["S1_cell_correction_or_expression"]
        )
        rep = validate_bottleneck_map_data(bad, repo_root=DEFAULT_REPO_ROOT)
        self.assertEqual(rep["status"], "FAIL")
        self.assertTrue(any("Preclinical overclaim guard" in e for e in rep["errors"]))

    def test_planted_defect_4_terminated_missing_reason_class_rejected(self) -> None:
        bad = copy.deepcopy(self.base_payload)
        term = next(s for s in bad["studies"] if s["source_id"] == "NCT06515002")
        del term["stages"]["S3_in_vivo_airway_delivery"]["termination_reason_class"]
        rep = validate_bottleneck_map_data(bad, repo_root=DEFAULT_REPO_ROOT)
        self.assertEqual(rep["status"], "FAIL")
        self.assertTrue(any("termination_reason_class" in e for e in rep["errors"]))

    def test_planted_defect_5_passage_sha256_or_text_drift_rejected(self) -> None:
        bad = copy.deepcopy(self.base_payload)
        cell = bad["studies"][0]["stages"]["S1_cell_correction_or_expression"]
        cell["source_file_sha256"] = "0" * 64
        cell["passage_text"] = "Altered non-verbatim paraphrase claim."
        rep = validate_bottleneck_map_data(bad, repo_root=DEFAULT_REPO_ROOT)
        self.assertEqual(rep["status"], "FAIL")
        self.assertTrue(any("source_file_sha256 drift" in e for e in rep["errors"]))
        self.assertTrue(any("not a verbatim substring" in e for e in rep["errors"]))

    def test_planted_defect_6_unverified_ledger_source_id_and_mixed_anchor_rejected(self) -> None:
        bad = copy.deepcopy(self.base_payload)
        bad["studies"][0]["source_id"] = "10.9999/unverified-synthetic-doi"
        bad["studies"][0]["stages"]["S1_cell_correction_or_expression"]["snapshot_nct_id"] = "NCT05248230"
        rep = validate_bottleneck_map_data(bad, repo_root=DEFAULT_REPO_ROOT)
        self.assertEqual(rep["status"], "FAIL")
        self.assertTrue(any("not found with verdict=PASS" in e for e in rep["errors"]))
        self.assertTrue(any("mixes snapshot field" in e for e in rep["errors"]))

    def test_planted_defect_7_path_traversal_rejected(self) -> None:
        bad = copy.deepcopy(self.base_payload)
        bad["ledger_file"] = "../outside-ledger.jsonl"
        bad["snapshot_file"] = "../outside-snapshot.json"
        rep = validate_bottleneck_map_data(bad, repo_root=DEFAULT_REPO_ROOT)
        self.assertEqual(rep["status"], "FAIL")
        self.assertTrue(any("Invalid ledger_file" in e for e in rep["errors"]))
        self.assertTrue(any("Invalid snapshot_file" in e for e in rep["errors"]))

    def test_empty_ledger_fails_closed(self) -> None:
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "empty_ledger.jsonl").write_text("", encoding="utf-8")
            (root / "empty_snapshot.json").write_text('{"studies": []}', encoding="utf-8")
            bad = copy.deepcopy(self.base_payload)
            bad["ledger_file"] = "empty_ledger.jsonl"
            bad["snapshot_file"] = "empty_snapshot.json"
            rep = validate_bottleneck_map_data(bad, repo_root=root)
            self.assertEqual(rep["status"], "FAIL")
            self.assertTrue(any("not found with verdict=PASS" in e for e in rep["errors"]))

    def test_all_five_stages_present_and_cli_json(self) -> None:
        for study in self.base_payload["studies"]:
            self.assertEqual(tuple(study["stages"].keys()), STAGE_KEYS)
        res = subprocess.run(
            [sys.executable, str(VALIDATOR_SCRIPT), "--check", "--json"],
            capture_output=True,
            text=True,
        )
        self.assertEqual(res.returncode, 0, msg=res.stderr)
        out = json.loads(res.stdout)
        self.assertEqual(out["status"], "PASS")
        self.assertEqual(out["studies_validated"], 8)


if __name__ == "__main__":
    unittest.main()
