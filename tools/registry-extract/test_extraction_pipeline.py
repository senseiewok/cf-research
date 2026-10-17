"""Comprehensive test suite for ECFSPR extraction pipeline (Task T-0009).

Tests:
1. Extraction pipeline across 5 report years (2020, 2021, 2022, 2023, 2024).
2. Schema validation of extracted dataset against row_schema.json.
3. 100% coverage of url and page fields across all rows.
4. Wilson 95% CI accuracy audit >= 99% on 100-cell random sample.
5. Presence and structure of disagreement_log.md and audit_100_cells.json.
6. Indicator coverage across ~15 standardized indicators.
"""
import hashlib
import json
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
SCHEMA_JSON = HERE / "row_schema.json"
DATASET_JSON = HERE / "ecfspr_extracted_dataset.json"
AUDIT_JSON = HERE / "audit_100_cells.json"
DISAGREEMENT_LOG = HERE / "disagreement_log.md"

import sys
sys.path.insert(0, str(HERE))
from extract_ecfspr import extract_all_ecfspr, ECFSPR_REPORTS
from audit_accuracy import perform_audit, wilson_95_ci
from validate import validate_dataset, validate_file


class ECFSPRExtractionPipelineTest(unittest.TestCase):
    def setUp(self):
        self.assertTrue(SCHEMA_JSON.exists(), f"Missing schema at {SCHEMA_JSON}")
        self.schema = json.loads(SCHEMA_JSON.read_text(encoding="utf-8"))

    def test_pipeline_extracts_all_five_years(self):
        rows = extract_all_ecfspr()
        self.assertGreaterEqual(len(rows), 100, f"Expected at least 100 rows, got {len(rows)}")

        years = {r["report_year"] for r in rows}
        expected_years = {2020, 2021, 2022, 2023, 2024}
        self.assertEqual(years, expected_years, f"Pipeline must extract all 5 years: {expected_years}")

        for y in expected_years:
            year_rows = [r for r in rows if r["report_year"] == y]
            self.assertGreaterEqual(len(year_rows), 10, f"Year {y} must have >= 10 rows")

    def test_schema_validates_extracted_dataset(self):
        self.assertTrue(DATASET_JSON.exists(), f"Missing dataset file {DATASET_JSON}")
        errors = validate_file(DATASET_JSON, SCHEMA_JSON)
        self.assertEqual(errors, [], f"Schema validation failed: {errors}")

    def test_one_hundred_percent_rows_carry_url_and_page(self):
        rows = json.loads(DATASET_JSON.read_text(encoding="utf-8"))
        for idx, row in enumerate(rows):
            rid = row.get("row_id", f"row_{idx}")
            self.assertIn("url", row, f"{rid} missing 'url'")
            self.assertIsInstance(row["url"], str, f"{rid} 'url' must be a string")
            self.assertTrue(row["url"].startswith("https://"), f"{rid} 'url' must be valid HTTPS URL")

            self.assertIn("page", row, f"{rid} missing 'page'")
            self.assertIsInstance(row["page"], int, f"{rid} 'page' must be an integer")
            self.assertGreaterEqual(row["page"], 1, f"{rid} 'page' must be >= 1")

            self.assertEqual(row["page"], row["page_number"], f"{rid} 'page' must equal 'page_number'")

    def test_cryptographic_source_hashes(self):
        rows = json.loads(DATASET_JSON.read_text(encoding="utf-8"))
        for row in rows:
            year = row["report_year"]
            expected_hash = ECFSPR_REPORTS[year]["sha256"]
            self.assertEqual(
                row.get("source_sha256"),
                expected_hash,
                f"Row {row['row_id']} SHA-256 hash mismatch for year {year}",
            )

    def test_indicator_diversity_covers_around_fifteen_indicators(self):
        rows = json.loads(DATASET_JSON.read_text(encoding="utf-8"))
        unique_indicators = {r["indicator_id"] for r in rows}
        self.assertGreaterEqual(len(unique_indicators), 14, f"Expected ~15 indicators, got {len(unique_indicators)}")

        # Check essential epidemiological indicators
        essential = {
            "pwcf_registered",
            "age_at_followup",
            "pwcf_adult_ge18",
            "age_at_diagnosis_months",
            "pwcf_f508del_allele",
            "pwcf_lung_transplant",
            "pwcf_liver_transplant",
            "pwcf_deceased",
            "age_at_death",
            "pwcf_seen_during_year",
            "pwcf_alive_end_of_year",
            "pwcf_mortality_cohort",
            "pwcf_uk_included",
        }
        self.assertTrue(essential.issubset(unique_indicators), f"Missing essential indicators: {essential - unique_indicators}")

    def test_hundred_cell_random_audit_accuracy_and_wilson_ci(self):
        summary = perform_audit(DATASET_JSON, seed=42, sample_size=100)
        self.assertEqual(summary["sampled_cells"], 100)
        self.assertGreaterEqual(summary["accuracy"], 0.99, f"Accuracy must be >= 99%, got {summary['accuracy']}")
        self.assertTrue(summary["meets_accuracy_threshold"])

        # Wilson 95% CI verification
        ci = summary["wilson_95_ci"]
        self.assertGreaterEqual(ci["lower"], 0.94)
        self.assertLessEqual(ci["upper"], 1.0)
        self.assertLess(ci["lower"], ci["upper"])

    def test_wilson_ci_formula_properties(self):
        lower, upper, margin = wilson_95_ci(100, 100)
        self.assertAlmostEqual(lower, 0.9630, places=3)
        self.assertAlmostEqual(upper, 1.0, places=3)

        lower_99, upper_99, margin_99 = wilson_95_ci(99, 100)
        self.assertAlmostEqual(lower_99, 0.9455, places=3)
        self.assertAlmostEqual(upper_99, 0.9982, places=3)

    def test_audit_json_and_disagreement_log_exist(self):
        self.assertTrue(AUDIT_JSON.exists(), f"Missing {AUDIT_JSON}")
        audit_data = json.loads(AUDIT_JSON.read_text(encoding="utf-8"))
        self.assertEqual(len(audit_data.get("audit_records", [])), 100)

        self.assertTrue(DISAGREEMENT_LOG.exists(), f"Missing {DISAGREEMENT_LOG}")
        log_content = DISAGREEMENT_LOG.read_text(encoding="utf-8")
        self.assertIn("Wilson 95%", log_content)
        self.assertIn("100-cell random sample", log_content)
        self.assertIn("2020", log_content)
        self.assertIn("2024", log_content)


if __name__ == "__main__":
    unittest.main()
