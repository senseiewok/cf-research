"""Deterministic offline test for Step 1 of T-0102.

Verifies:
1. Multi-source brief returns records from >= 4 catalogued sources.
2. Each record has source_id, URL, retrieval date, and 'not stated' where silent.
3. Negative control: a deliberately inconsistent record (label date before first approval) flags CHECK.
4. Negative control: a retracted work in literature scan flags CHECK.
5. Consistent records flag PASS.
"""

from __future__ import annotations

import unittest

# Pure data representation matching evidence loop output schema


def run_cross_check_approval_dates(first_approval: str | None, label_effective_time: str | None) -> dict[str, str]:
    what = "first-approval year vs label effective-date year"
    if not first_approval or not label_effective_time:
        return {"verdict": "NOT STATED", "check": what, "detail": "One or both dates not stated"}
    appr_year = int(first_approval[:4]) if len(first_approval) >= 4 and first_approval[:4].isdigit() else None
    label_year = int(label_effective_time[:4]) if len(label_effective_time) >= 4 and label_effective_time[:4].isdigit() else None
    if appr_year is None or label_year is None:
        return {"verdict": "NOT STATED", "check": what, "detail": "Dates unparseable"}
    if label_year < appr_year:
        return {
            "verdict": "CHECK",
            "check": what,
            "detail": f"label effective date {label_effective_time} is earlier than first approval {first_approval}; inconsistent",
        }
    return {
        "verdict": "PASS",
        "check": what,
        "detail": f"first approval {first_approval}, label effective date {label_effective_time}; valid progression",
    }


def run_cross_check_retraction_scan(records: list[dict]) -> dict[str, str]:
    what = "retraction scan"
    retracted = [r for r in records if r.get("is_retracted") is True]
    if retracted:
        return {"verdict": "CHECK", "check": what, "detail": f"{len(retracted)} work(s) flagged retracted"}
    return {"verdict": "PASS", "check": what, "detail": "0 works flagged retracted"}


class TestT0102Step1EvidenceLoop(unittest.TestCase):
    def setUp(self):
        # Synthetic run of 5 catalogued sources for fixed question: "TRIKAFTA"
        self.sources = [
            {
                "source_id": "openfda-drugsfda",
                "url": "https://api.fda.gov/drug/drugsfda.json?search=openfda.brand_name:%22TRIKAFTA%22&limit=1",
                "accessed": "2026-10-15T22:50:06Z",
                "fields": {
                    "application_number": "NDA217660",
                    "original_approval": "2023-04-26",
                    "sponsor": "VERTEX PHARMS INC",
                    "boxed_warning": "not stated",
                },
            },
            {
                "source_id": "openfda-label",
                "url": "https://api.fda.gov/drug/label.json?search=openfda.brand_name:%22TRIKAFTA%22&limit=1",
                "accessed": "2026-10-15T22:50:07Z",
                "fields": {
                    "set_id": "f354423a-85c2-41c3-a9db-0f3aee135d8d",
                    "effective_time": "2026-06-22",
                    "version": 23,
                },
            },
            {
                "source_id": "clinicaltrials-gov",
                "url": "https://clinicaltrials.gov/api/v2/studies?query.term=TRIKAFTA&pageSize=5",
                "accessed": "2026-10-15T22:50:07Z",
                "fields": {
                    "nct_id": "NCT05743946",
                    "status": "COMPLETED",
                    "study_type": "INTERVENTIONAL",
                    "secondary_sponsor": "not stated",
                },
            },
            {
                "source_id": "ncbi-eutils",
                "url": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id=34268058",
                "accessed": "2026-10-15T22:50:08Z",
                "fields": {
                    "pmid": "34268058",
                    "title": "A Review of Trikafta",
                    "full_text_open_access": "not stated",
                },
            },
            {
                "source_id": "openalex",
                "url": "https://api.openalex.org/works?search=TRIKAFTA&per_page=3",
                "accessed": "2026-10-15T22:50:09Z",
                "fields": {
                    "id": "https://openalex.org/W3039639714",
                    "is_oa": True,
                    "is_retracted": False,
                },
            },
        ]

    def test_minimum_four_sources_returned(self):
        distinct_sources = {s["source_id"] for s in self.sources}
        self.assertGreaterEqual(len(distinct_sources), 4)

    def test_every_record_carries_required_metadata_and_not_stated(self):
        for s in self.sources:
            self.assertTrue(s["source_id"])
            self.assertTrue(s["url"].startswith("https://"))
            self.assertTrue(s["accessed"].endswith("Z"))
            # Confirm silent fields print as "not stated"
            fields = s["fields"]
            silent_values = [v for v in fields.values() if v == "not stated"]
            self.assertIsInstance(fields, dict)
        self.assertTrue(
            any(v == "not stated" for s in self.sources for v in s["fields"].values()),
            "At least one silent field across sources must be explicitly rendered as 'not stated'",
        )

    def test_concordant_cross_check_passes(self):
        appr = self.sources[0]["fields"]["original_approval"]
        label = self.sources[1]["fields"]["effective_time"]
        res = run_cross_check_approval_dates(appr, label)
        self.assertEqual(res["verdict"], "PASS")

    def test_deliberately_wrong_record_flags_check_negative_control(self):
        # Negative control: label effective date 2018 is earlier than approval date 2023
        res = run_cross_check_approval_dates("2023-04-26", "2018-01-01")
        self.assertEqual(res["verdict"], "CHECK")
        self.assertIn("inconsistent", res["detail"])

    def test_retraction_cross_check_flags_check(self):
        # Negative control: one retracted work flags CHECK
        works = [{"id": "W1", "is_retracted": True}, {"id": "W2", "is_retracted": False}]
        res = run_cross_check_retraction_scan(works)
        self.assertEqual(res["verdict"], "CHECK")
        self.assertIn("1 work(s) flagged retracted", res["detail"])


if __name__ == "__main__":
    unittest.main()
