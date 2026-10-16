"""Tests for cf-registry-extract row schema and ECFSPR 2022 gold sample.

Validates that tools/registry-extract/SCHEMA.md exists, that row_schema.json enforces
the frozen schema invariants, and that the 20+ row hand-extracted gold sample from
ECFSPR 2022 page 8 passes schema validation and cryptographic source verification.
"""
import copy
import hashlib
import json
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
SCHEMA_DOC = HERE / "SCHEMA.md"
SCHEMA_JSON = HERE / "row_schema.json"
GOLD_SAMPLE = HERE / "gold_sample_ecfspr_2022_p8.json"

import sys
sys.path.insert(0, str(HERE))
from validate import validate_dataset, validate_file, validate_row


class RegistryExtractSchemaTest(unittest.TestCase):
    def setUp(self):
        self.assertTrue(SCHEMA_JSON.exists(), f"Missing {SCHEMA_JSON}")
        self.schema = json.loads(SCHEMA_JSON.read_text(encoding="utf-8"))

    def test_schema_md_exists_and_documents_required_keys(self):
        self.assertTrue(SCHEMA_DOC.exists(), f"SCHEMA.md must exist at {SCHEMA_DOC}")
        content = SCHEMA_DOC.read_text(encoding="utf-8")
        self.assertGreater(len(content), 500, "SCHEMA.md should be thorough")
        self.assertIn("source_sha256", content)
        self.assertIn("value_low", content)
        self.assertIn("value_high", content)
        self.assertIn("median_iqr", content)
        self.assertIn("ECFSPR 2022", content)

    def test_schema_json_structure(self):
        self.assertEqual(self.schema.get("type"), "object")
        required = set(self.schema.get("required", []))
        expected_required = {
            "row_id",
            "registry_id",
            "report_title",
            "report_year",
            "page_number",
            "table_id",
            "indicator_id",
            "indicator_name",
            "stratum_sex",
            "unit",
            "source_sha256",
            "raw_text_quote",
            "extraction_kind",
        }
        self.assertTrue(expected_required.issubset(required), f"Missing required keys: {expected_required - required}")
        self.assertFalse(self.schema.get("additionalProperties", True))

    def test_gold_sample_passes_schema_validation(self):
        self.assertTrue(GOLD_SAMPLE.exists(), f"Missing {GOLD_SAMPLE}")
        rows = json.loads(GOLD_SAMPLE.read_text(encoding="utf-8"))
        self.assertGreaterEqual(len(rows), 20, f"Gold sample must contain at least 20 rows, got {len(rows)}")

        errors = validate_file(GOLD_SAMPLE, SCHEMA_JSON)
        self.assertEqual(errors, [], f"Gold sample validation failed: {errors}")

    def test_gold_sample_cryptographic_source_hash(self):
        rows = json.loads(GOLD_SAMPLE.read_text(encoding="utf-8"))
        expected_hash = "22868219bcaeb5e4b2090643b0d5b1e8be91fddb1cf37ba74bf1bfaeec11506e"

        for row in rows:
            self.assertEqual(
                row.get("source_sha256"),
                expected_hash,
                f"Row {row.get('row_id')} source_sha256 mismatch",
            )

        # If PDF is present in sources/downloads/, verify hash against actual file
        pdf_path = REPO_ROOT / "sources" / "downloads" / "ECFSPR_Annual_Report_2022_v1.0.pdf"
        if pdf_path.exists():
            computed_hash = hashlib.sha256(pdf_path.read_bytes()).hexdigest()
            self.assertEqual(computed_hash, expected_hash, "Actual PDF hash must match schema declaration")

    def test_gold_sample_indicators_and_strata(self):
        rows = json.loads(GOLD_SAMPLE.read_text(encoding="utf-8"))
        row_ids = [r["row_id"] for r in rows]
        self.assertEqual(len(row_ids), len(set(row_ids)), "All row_ids in gold sample must be unique")

        # Must cover both count_and_percentage and non-parametric median_iqr rows
        count_rows = [r for r in rows if r["unit"] == "count_and_percentage"]
        median_rows = [r for r in rows if r["median_iqr"] is not None]

        self.assertGreaterEqual(len(count_rows), 10, "Must have at least 10 count/percentage rows")
        self.assertGreaterEqual(len(median_rows), 6, "Must have at least 6 median_iqr rows")

        for r in median_rows:
            self.assertIsNotNone(r["value_low"])
            self.assertIsNotNone(r["value_high"])
            self.assertEqual(r["median_iqr"]["q1"], r["value_low"])
            self.assertEqual(r["median_iqr"]["q3"], r["value_high"])

    def test_rejection_missing_required_field(self):
        sample = {
            "row_id": "TEST-01",
            "registry_id": "ecfspr",
            "report_title": "Test Report",
            "report_year": 2022,
            "page_number": 8,
            "table_id": "Test Table",
            "indicator_id": "test_ind",
            "indicator_name": "Test Indicator",
            "stratum_sex": "females",
            "unit": "count",
            # missing source_sha256
            "raw_text_quote": "123",
            "extraction_kind": "observed",
        }
        errors = validate_row(sample, self.schema)
        self.assertTrue(any("missing required field 'source_sha256'" in e for e in errors))

    def test_rejection_invalid_sha256(self):
        sample = {
            "row_id": "TEST-01",
            "registry_id": "ecfspr",
            "report_title": "Test Report",
            "report_year": 2022,
            "page_number": 8,
            "table_id": "Test Table",
            "indicator_id": "test_ind",
            "indicator_name": "Test Indicator",
            "stratum_sex": "females",
            "unit": "count",
            "source_sha256": "INVALID_HASH_NOT_64_CHARS",
            "raw_text_quote": "123",
            "extraction_kind": "observed",
        }
        errors = validate_row(sample, self.schema)
        self.assertTrue(any("does not match pattern" in e for e in errors))

    def test_rejection_percentage_out_of_bounds(self):
        sample = {
            "row_id": "TEST-01",
            "registry_id": "ecfspr",
            "report_title": "Test Report",
            "report_year": 2022,
            "page_number": 8,
            "table_id": "Test Table",
            "indicator_id": "test_ind",
            "indicator_name": "Test Indicator",
            "stratum_sex": "females",
            "percentage": 105.4,  # > 100.0
            "unit": "percentage",
            "source_sha256": "22868219bcaeb5e4b2090643b0d5b1e8be91fddb1cf37ba74bf1bfaeec11506e",
            "raw_text_quote": "105.4%",
            "extraction_kind": "observed",
        }
        errors = validate_row(sample, self.schema)
        self.assertTrue(any("greater than maximum 100.0" in e for e in errors))

    def test_rejection_unexpected_property(self):
        sample = {
            "row_id": "TEST-01",
            "registry_id": "ecfspr",
            "report_title": "Test Report",
            "report_year": 2022,
            "page_number": 8,
            "table_id": "Test Table",
            "indicator_id": "test_ind",
            "indicator_name": "Test Indicator",
            "stratum_sex": "females",
            "unit": "count",
            "source_sha256": "22868219bcaeb5e4b2090643b0d5b1e8be91fddb1cf37ba74bf1bfaeec11506e",
            "raw_text_quote": "123",
            "extraction_kind": "observed",
            "unauthorized_extra_field": "disallowed",
        }
        errors = validate_row(sample, self.schema)
        self.assertTrue(any("unexpected property 'unauthorized_extra_field'" in e for e in errors))


if __name__ == "__main__":
    unittest.main()
