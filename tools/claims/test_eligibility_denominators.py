"""Tests for T-0134: Modulator eligibility denominators ledger and claims verification.

Asserts:
- Non-empty structured fields in modulator-eligibility-denominators-2026.json.
- Every registry_id's data_year and url match sources/catalog.yaml.
- Short quotes <= 24 words, and every short_quote appears in the claims/evidence files.
- When primary PDFs / EMA JSON are present in sources/downloads, every quote is verified directly against the primary file.
- Zero values computed across rows (non-pooling invariant).
- Zero poster citations.
- check_claims.py passes with exit code 0.
- check_source_overlap.py --quotes-only passes with exit code 0.
"""
from __future__ import annotations

import json
import re
import subprocess
import sys
import unittest
from pathlib import Path

import yaml

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
LANDSCAPE_DIR = REPO_ROOT / "landscape"
LEDGER_JSON = LANDSCAPE_DIR / "modulator-eligibility-denominators-2026.json"
CLAIMS_JSON = LANDSCAPE_DIR / "modulator-eligibility-denominators-2026-claims.json"
EVIDENCE_MAP = LANDSCAPE_DIR / "modulator-eligibility-denominators-2026-evidence-map.json"
EVIDENCE_TXT = LANDSCAPE_DIR / "modulator-eligibility-denominators-2026-evidence.txt"
DOSSIER_MD = LANDSCAPE_DIR / "modulator-eligibility-denominators-2026.md"
CATALOG_YAML = REPO_ROOT / "sources" / "catalog.yaml"


def _norm(s: str) -> str:
    return re.sub(r"\s+", " ", s).strip()


def _find_downloads_dir() -> Path | None:
    candidates = [
        REPO_ROOT / "sources" / "downloads",
        REPO_ROOT.parent / "cf-research" / "sources" / "downloads",
        REPO_ROOT.parent.parent / "cf-research" / "sources" / "downloads",
        REPO_ROOT.parent.parent.parent / "cf-research" / "sources" / "downloads",
    ]
    for cand in candidates:
        if (cand / "ECFSPR_Annual_Report_2024_v1.0.pdf").is_file():
            return cand
    return None


class TestEligibilityDenominators(unittest.TestCase):
    def setUp(self):
        self.assertTrue(LEDGER_JSON.exists(), f"Missing {LEDGER_JSON}")
        self.assertTrue(CLAIMS_JSON.exists(), f"Missing {CLAIMS_JSON}")
        self.assertTrue(EVIDENCE_MAP.exists(), f"Missing {EVIDENCE_MAP}")
        self.assertTrue(EVIDENCE_TXT.exists(), f"Missing {EVIDENCE_TXT}")
        self.assertTrue(DOSSIER_MD.exists(), f"Missing {DOSSIER_MD}")
        self.ledger = json.loads(LEDGER_JSON.read_text(encoding="utf-8"))

    def test_ledger_structure_and_required_fields(self):
        self.assertIn("metadata", self.ledger)
        self.assertIn("regulator_labels", self.ledger)
        self.assertIn("registry_reports", self.ledger)

        reg_fields = [
            "regulator",
            "brand_name",
            "nonproprietary_name",
            "application_id",
            "initial_approval_date",
            "decision_or_update_date",
            "age_rule",
            "variant_rule",
            "exact_eligibility_wording",
            "source_citation",
        ]
        for label in self.ledger["regulator_labels"]:
            for f in reg_fields:
                self.assertIn(f, label, f"Missing {f} in regulator label")
                self.assertTrue(bool(str(label[f]).strip()), f"Empty field {f} in label")

        report_fields = [
            "registry_id",
            "report_year",
            "data_year",
            "base_population",
            "age_rule",
            "variant_list_and_date",
            "page",
            "url",
            "short_quote",
        ]
        for row in self.ledger["registry_reports"]:
            for f in report_fields:
                self.assertIn(f, row, f"Missing {f} in registry report row")
                self.assertTrue(bool(str(row[f]).strip()), f"Empty field {f} in row")

    def test_registry_rows_match_catalog_data_year_and_url(self):
        catalog_data = yaml.safe_load(CATALOG_YAML.read_text(encoding="utf-8"))
        by_id = {entry["id"]: entry for entry in catalog_data.get("sources", []) if "id" in entry}
        for row in self.ledger["registry_reports"]:
            rid = row["registry_id"]
            if rid in by_id:
                cat_entry = by_id[rid]
                self.assertEqual(
                    row["data_year"],
                    cat_entry["data_year"],
                    f"{rid} data_year ({row['data_year']}) != catalog data_year ({cat_entry['data_year']})",
                )
                expected_url = cat_entry.get("url") or cat_entry.get("landing_page")
                self.assertEqual(row["url"], expected_url, f"{rid} url mismatch with catalog")

    def test_short_quotes_under_word_limit_and_grounded(self):
        ev_text = _norm(EVIDENCE_TXT.read_text(encoding="utf-8"))
        for row in self.ledger["registry_reports"]:
            quote = row["short_quote"]
            word_count = len(quote.split())
            self.assertLessEqual(
                word_count,
                24,
                f"Quote in {row['registry_id']} exceeds 24 words ({word_count} words): '{quote}'",
            )
            self.assertIn(_norm(quote), ev_text, f"short_quote for {row['registry_id']} missing from evidence text")

        for label in self.ledger["regulator_labels"]:
            self.assertIn(
                _norm(label["exact_eligibility_wording"]),
                ev_text,
                f"exact_eligibility_wording for {label['regulator']} {label['brand_name']} missing from evidence text",
            )

    def test_quotes_verified_against_primary_downloads_when_available(self):
        dl = _find_downloads_dir()
        if dl is None:
            return
        from pypdf import PdfReader
        import logging
        logging.getLogger("pypdf").setLevel(logging.ERROR)

        pdf_map = {
            "ecfspr-adr-2020": "ECFSPR_Annual_Report_2020_v1.0.pdf",
            "ecfspr-adr-2021": "ECFSPR_Annual_Report_2021.pdf",
            "ecfspr-adr-2022": "ECFSPR_Annual_Report_2022_v1.0.pdf",
            "ecfspr-adr-2023": "ECFSPR_Annual_Report_2023_v1.1.pdf",
            "ecfspr-adr-2024": "ECFSPR_Annual_Report_2024_v1.0.pdf",
            "cffpr-adr-2024": "CFFPR_Annual_Data_Report_2024.pdf",
        }
        for row in self.ledger["registry_reports"]:
            rid = row["registry_id"]
            pdf_path = dl / pdf_map[rid]
            reader = PdfReader(str(pdf_path))
            page_txt = _norm(reader.pages[row["page"] - 1].extract_text() or "")
            self.assertIn(
                _norm(row["short_quote"]),
                page_txt,
                f"{rid} p.{row['page']} does not contain verbatim short_quote: {row['short_quote']!r}",
            )

    def test_no_values_computed_across_rows(self):
        """Asserts that no field holds a pooled, summed, or averaged value across rows."""
        disallowed_tokens = {"total across", "pooled", "summed", "combined registries", "cross-registry"}
        for row in self.ledger["registry_reports"]:
            for k, val in row.items():
                val_str = str(val).lower()
                for token in disallowed_tokens:
                    self.assertNotIn(
                        token,
                        val_str,
                        f"Found computed/pooled token '{token}' in field {k} of {row['registry_id']}",
                    )

    def test_no_poster_citations(self):
        """Posters are never cited as evidence."""
        ledger_text = LEDGER_JSON.read_text(encoding="utf-8").lower()
        dossier_text = DOSSIER_MD.read_text(encoding="utf-8").lower()
        self.assertNotIn("nacfc poster", ledger_text)
        self.assertNotIn("nacfc poster", dossier_text)
        self.assertNotIn("poster citation", ledger_text)

    def test_check_claims_passes(self):
        cmd = [
            sys.executable,
            str(REPO_ROOT / "tools" / "claims" / "check_claims.py"),
            str(CLAIMS_JSON),
            "--evidence-map",
            str(EVIDENCE_MAP),
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, encoding="utf-8")
        self.assertEqual(
            result.returncode,
            0,
            f"check_claims.py failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}",
        )

    def test_check_source_overlap_quotes_only_passes(self):
        cmd = [
            sys.executable,
            str(REPO_ROOT / "tools" / "sources" / "check_source_overlap.py"),
            "--quotes-only",
            str(DOSSIER_MD),
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, encoding="utf-8")
        self.assertEqual(
            result.returncode,
            0,
            f"check_source_overlap.py failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}",
        )

    def test_check_numbers_passes(self):
        allow_file = LANDSCAPE_DIR / "modulator-eligibility-denominators-2026-allow.txt"
        cmd = [
            sys.executable,
            str(REPO_ROOT / "tools" / "sources" / "check_numbers.py"),
            "--evidence",
            str(EVIDENCE_TXT),
            "--allow-file",
            str(allow_file),
            str(DOSSIER_MD),
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, encoding="utf-8")
        self.assertEqual(
            result.returncode,
            0,
            f"check_numbers.py failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}",
        )


if __name__ == "__main__":
    unittest.main()
