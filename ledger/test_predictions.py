#!/usr/bin/env python3
"""Unit tests for ledger/predictions.jsonl and ledger/check_predictions.py.

Verifies:
1. The repository's ledger/predictions.jsonl is 100% valid.
2. Schema enforcement on required fields, types, and allowed outcomes.
3. Discovery Protocol invariants: unique IDs, at most one open prediction,
   outside reader review before supported claims, catalog source checking.
4. Command-line interface operation.
"""
from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from ledger.check_predictions import (
    ALLOWED_OUTCOMES,
    REQUIRED_FIELDS,
    load_catalog_source_ids,
    validate_prediction_record,
    validate_predictions,
)

REPO_ROOT = Path(__file__).resolve().parent.parent
PREDICTIONS_JSONL = REPO_ROOT / "ledger" / "predictions.jsonl"
CATALOG_YAML = REPO_ROOT / "sources" / "catalog.yaml"


class TestRepoPredictionLedger(unittest.TestCase):
    """Test that the live repository file passes validation."""

    def test_repo_predictions_file_exists(self) -> None:
        self.assertTrue(PREDICTIONS_JSONL.is_file(), "ledger/predictions.jsonl must exist")

    def test_repo_predictions_pass_validation(self) -> None:
        records, errors = validate_predictions(PREDICTIONS_JSONL, CATALOG_YAML)
        self.assertEqual(errors, [], f"Repository predictions.jsonl failed validation: {errors}")
        self.assertGreaterEqual(len(records), 1, "At least one prediction must be present")

    def test_seed_prediction_pred_0001(self) -> None:
        records, errors = validate_predictions(PREDICTIONS_JSONL, CATALOG_YAML)
        self.assertEqual(errors, [])
        pred_map = {r["id"]: r for r in records}
        self.assertIn("PRED-0001", pred_map)
        p1 = pred_map["PRED-0001"]
        self.assertEqual(p1["id"], "PRED-0001")
        self.assertEqual(p1["date"], "2026-11-14")
        self.assertIn("variability check parameter W >= 0", p1["prediction"])
        self.assertIn("Observed W >= 0", p1["refutation_observation"])
        self.assertIn("10.1098/rstb.1982.0156", p1["prior_art_dois"])
        self.assertIn("10.5281/zenodo.15476468", p1["prior_art_dois"])
        self.assertEqual(p1["measured_public_data"], "none")
        self.assertEqual(p1["outcome"], "could_not_test")
        self.assertTrue(p1["checked_by_outside_reader"])


class TestPredictionRecordValidation(unittest.TestCase):
    """Unit tests for single-record schema enforcement."""

    def setUp(self) -> None:
        self.valid_record = {
            "id": "PRED-0099",
            "date": "2026-11-15",
            "prediction": "Single test sentence predicting bounded parameter behavior.",
            "refutation_observation": "Observed value outside confidence bounds.",
            "prior_art_dois": ["10.1098/rstb.1982.0156"],
            "standard_method": "Maximum-likelihood estimation.",
            "measured_public_data": "none",
            "outcome": "could_not_test",
            "checked_by_outside_reader": True,
        }
        self.catalog_ids = {"ecfspr-adr-2024", "clinicaltrials-gov"}

    def test_valid_record(self) -> None:
        seen: set[str] = set()
        errs = validate_prediction_record(self.valid_record, 1, seen, self.catalog_ids)
        self.assertEqual(errs, [])
        self.assertIn("PRED-0099", seen)

    def test_missing_required_fields(self) -> None:
        for field in REQUIRED_FIELDS:
            with self.subTest(field=field):
                rec = dict(self.valid_record)
                del rec[field]
                seen: set[str] = set()
                errs = validate_prediction_record(rec, 1, seen, self.catalog_ids)
                self.assertTrue(
                    any(f"Missing required field '{field}'" in e for e in errs),
                    f"Expected missing field error for {field}, got {errs}",
                )

    def test_unknown_fields(self) -> None:
        rec = dict(self.valid_record)
        rec["unsupported_extra_field"] = "unexpected"
        seen: set[str] = set()
        errs = validate_prediction_record(rec, 1, seen, self.catalog_ids)
        self.assertTrue(any("Unknown fields detected" in e for e in errs))

    def test_optional_fields_allowed(self) -> None:
        rec = dict(self.valid_record)
        rec["notes"] = "Test note"
        rec["data_request_plan"] = "specs/test_plan.md"
        rec["catalog_sources_searched"] = ["clinicaltrials-gov"]
        seen: set[str] = set()
        errs = validate_prediction_record(rec, 1, seen, self.catalog_ids)
        self.assertEqual(errs, [])

    def test_invalid_id_patterns(self) -> None:
        bad_ids = ["PRED-1", "PRED-01", "PRED-001", "pred-0001", "L-0001", "PRED-12345", ""]
        for bad_id in bad_ids:
            with self.subTest(bad_id=bad_id):
                rec = dict(self.valid_record)
                rec["id"] = bad_id
                seen: set[str] = set()
                errs = validate_prediction_record(rec, 1, seen, self.catalog_ids)
                self.assertTrue(any("'id' must match pattern" in e for e in errs))

    def test_duplicate_id_detection(self) -> None:
        seen = {"PRED-0099"}
        errs = validate_prediction_record(self.valid_record, 2, seen, self.catalog_ids)
        self.assertTrue(any("Duplicate prediction id 'PRED-0099'" in e for e in errs))

    def test_invalid_dates(self) -> None:
        bad_dates = ["2026/11/15", "2026-11", "2026-02-31", "not-a-date", ""]
        for bad_date in bad_dates:
            with self.subTest(bad_date=bad_date):
                rec = dict(self.valid_record)
                rec["date"] = bad_date
                seen: set[str] = set()
                errs = validate_prediction_record(rec, 1, seen, self.catalog_ids)
                self.assertTrue(any("'date'" in e for e in errs))

    def test_invalid_outcome_enum(self) -> None:
        bad_outcomes = ["pending", "verified", "passed", "failed", "unknown"]
        for bad in bad_outcomes:
            with self.subTest(bad=bad):
                rec = dict(self.valid_record)
                rec["outcome"] = bad
                seen: set[str] = set()
                errs = validate_prediction_record(rec, 1, seen, self.catalog_ids)
                self.assertTrue(any("'outcome' must be one of" in e for e in errs))

    def test_all_allowed_outcomes_accepted(self) -> None:
        for outcome in ALLOWED_OUTCOMES:
            with self.subTest(outcome=outcome):
                rec = dict(self.valid_record)
                rec["outcome"] = outcome
                rec["checked_by_outside_reader"] = True
                seen: set[str] = set()
                errs = validate_prediction_record(rec, 1, seen, self.catalog_ids)
                self.assertEqual(errs, [])

    def test_checked_by_outside_reader_type_guard(self) -> None:
        bad_types = [1, 0, "true", "false", None, []]
        for val in bad_types:
            with self.subTest(val=val):
                rec = dict(self.valid_record)
                rec["checked_by_outside_reader"] = val
                seen: set[str] = set()
                errs = validate_prediction_record(rec, 1, seen, self.catalog_ids)
                self.assertTrue(any("'checked_by_outside_reader' must be a boolean" in e for e in errs))

    def test_supported_outcome_requires_outside_reader(self) -> None:
        rec = dict(self.valid_record)
        rec["outcome"] = "supported"
        rec["checked_by_outside_reader"] = False
        seen: set[str] = set()
        errs = validate_prediction_record(rec, 1, seen, self.catalog_ids)
        self.assertTrue(any("checked_by_outside_reader is false" in e for e in errs))

    def test_prior_art_dois_format(self) -> None:
        bad_dois = ["not-a-doi", "http://doi.org/10.1234", "10.123", ""]
        for bad in bad_dois:
            with self.subTest(bad=bad):
                rec = dict(self.valid_record)
                rec["prior_art_dois"] = [bad]
                seen: set[str] = set()
                errs = validate_prediction_record(rec, 1, seen, self.catalog_ids)
                self.assertTrue(any("is not a valid DOI format" in e for e in errs))

    def test_measured_public_data_catalog_resolution(self) -> None:
        # 'none' passes
        rec = dict(self.valid_record)
        rec["measured_public_data"] = "none"
        seen: set[str] = set()
        errs = validate_prediction_record(rec, 1, seen, self.catalog_ids)
        self.assertEqual(errs, [])

        # Valid catalog id passes
        rec["measured_public_data"] = "ecfspr-adr-2024"
        seen = set()
        errs = validate_prediction_record(rec, 1, seen, self.catalog_ids)
        self.assertEqual(errs, [])

        # Unrecognized catalog id fails
        rec["measured_public_data"] = "non-existent-source-id"
        seen = set()
        errs = validate_prediction_record(rec, 1, seen, self.catalog_ids)
        self.assertTrue(any("not recognized in sources/catalog.yaml" in e for e in errs))


class TestFileLevelInvariants(unittest.TestCase):
    """Test multi-line file-level invariants (e.g. at most one open prediction)."""

    def test_multiple_open_predictions_rejected(self) -> None:
        rec1 = {
            "id": "PRED-0010",
            "date": "2026-11-14",
            "prediction": "Prediction one.",
            "refutation_observation": "Observation one.",
            "prior_art_dois": ["10.1098/rstb.1982.0156"],
            "standard_method": "Method one.",
            "measured_public_data": "none",
            "outcome": "open",
            "checked_by_outside_reader": False,
        }
        rec2 = {
            "id": "PRED-0011",
            "date": "2026-11-15",
            "prediction": "Prediction two.",
            "refutation_observation": "Observation two.",
            "prior_art_dois": ["10.1098/rstb.1982.0156"],
            "standard_method": "Method two.",
            "measured_public_data": "none",
            "outcome": "open",
            "checked_by_outside_reader": False,
        }
        with tempfile.NamedTemporaryFile("w", encoding="utf-8", suffix=".jsonl", delete=False) as tf:
            tf.write(json.dumps(rec1) + "\n" + json.dumps(rec2) + "\n")
            temp_path = Path(tf.name)

        try:
            records, errors = validate_predictions(temp_path, CATALOG_YAML)
            self.assertEqual(len(records), 2)
            self.assertTrue(any("found 2 predictions with outcome 'open'" in e for e in errors))
        finally:
            temp_path.unlink(missing_ok=True)

    def test_single_open_with_closed_accepted(self) -> None:
        rec1 = {
            "id": "PRED-0010",
            "date": "2026-11-14",
            "prediction": "Prediction one.",
            "refutation_observation": "Observation one.",
            "prior_art_dois": ["10.1098/rstb.1982.0156"],
            "standard_method": "Method one.",
            "measured_public_data": "none",
            "outcome": "could_not_test",
            "checked_by_outside_reader": True,
        }
        rec2 = {
            "id": "PRED-0011",
            "date": "2026-11-15",
            "prediction": "Prediction two.",
            "refutation_observation": "Observation two.",
            "prior_art_dois": ["10.1098/rstb.1982.0156"],
            "standard_method": "Method two.",
            "measured_public_data": "none",
            "outcome": "open",
            "checked_by_outside_reader": False,
        }
        with tempfile.NamedTemporaryFile("w", encoding="utf-8", suffix=".jsonl", delete=False) as tf:
            tf.write(json.dumps(rec1) + "\n" + json.dumps(rec2) + "\n")
            temp_path = Path(tf.name)

        try:
            records, errors = validate_predictions(temp_path, CATALOG_YAML)
            self.assertEqual(len(records), 2)
            self.assertEqual(errors, [])
        finally:
            temp_path.unlink(missing_ok=True)


class TestCLIExecution(unittest.TestCase):
    """Test CLI invocation of check_predictions.py."""

    def test_cli_on_repo_file(self) -> None:
        cli_path = REPO_ROOT / "ledger" / "check_predictions.py"
        res = subprocess.run(
            [sys.executable, str(cli_path), "--file", str(PREDICTIONS_JSONL), "--json"],
            capture_output=True,
            text=True,
        )
        self.assertEqual(res.returncode, 0, f"CLI exited with error: {res.stderr}")
        report = json.loads(res.stdout)
        self.assertEqual(report["status"], "PASS")
        self.assertEqual(report["error_count"], 0)


if __name__ == "__main__":
    unittest.main()
