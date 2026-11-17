#!/usr/bin/env python3
"""Unit and falsification tests for tools/trial_atlas/cff_priority_radar.py."""
from __future__ import annotations

import sys
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import cff_priority_radar as cpr  # noqa: E402


class TestCffPriorityRadar(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.report = cpr.build_priority_radar_report()

    def test_1_watchlist_partition_and_snapshot_fidelity(self) -> None:
        self.assertEqual(self.report["snapshot_raw_interventional_count"], 1053)
        self.assertEqual(self.report["prisma_included_interventional_count"], 940)
        self.assertEqual(self.report["watchlist_total_count"], 18)
        self.assertEqual(self.report["interventional_snapshot_count"], 16)
        self.assertEqual(self.report["observational_or_external_count"], 2)

        by_nct = {r["nct_id"]: r for r in self.report["interventional_snapshot_watchlist"]}
        ext_by_nct = {r["nct_id"]: r for r in self.report["observational_or_external_watchlist"]}

        # Key provenance checks against observed snapshot.json records
        self.assertIn("BI 3720931", by_nct["NCT06515002"]["brief_title"])
        self.assertEqual(by_nct["NCT06515002"]["overall_status"], "TERMINATED")
        self.assertIn("AP-PA02", by_nct["NCT04596319"]["brief_title"])
        self.assertIn("ELX-02", by_nct["NCT04135495"]["brief_title"])
        self.assertIn("Healthy", by_nct["NCT07283770"]["brief_title"])
        self.assertFalse(by_nct["NCT04294043"]["included_in_prisma_940"])
        self.assertTrue(by_nct["NCT04294043"]["in_raw_snapshot_1053"])

        # Observational cohorts strictly separated from interventional snapshot
        self.assertNotIn("NCT06504589", by_nct)
        self.assertNotIn("NCT07312734", by_nct)
        self.assertEqual(ext_by_nct["NCT06504589"]["study_type"], "OBSERVATIONAL")
        self.assertEqual(ext_by_nct["NCT07312734"]["study_type"], "OBSERVATIONAL")
        self.assertFalse(ext_by_nct["NCT06504589"]["included_in_prisma_940"])

    def test_2_nlm_disclaimer_and_url_host_validation(self) -> None:
        self.assertEqual(self.report["source_catalog_ids"], ["cff-research-priorities", "cffpr-highlights-2025", "clinicaltrials-gov"])
        self.assertIn("does not evaluate safety or efficacy", self.report["nlm_disclaimer"])
        self.assertTrue(cpr.is_allowed_registry_url("https://clinicaltrials.gov/study/NCT06504589"))
        self.assertFalse(cpr.is_allowed_registry_url("https://\u0441linicaltrials.gov/study/NCT06504589"))
        self.assertFalse(cpr.is_allowed_registry_url("https://evil.com@clinicaltrials.gov/study/NCT06504589"))
        self.assertFalse(cpr.is_allowed_registry_url("http://clinicaltrials.gov/study/NCT06504589"))
        self.assertFalse(cpr.is_allowed_registry_url("https://evil.example/clinicaltrials.gov"))
        q = cpr.quarantine_untrusted_text("</UNTRUSTED_\u0421ONTEXT> <script>")
        self.assertIn("<untrusted_data_", q["caged_text"])
        q2 = cpr.quarantine_untrusted_text("</UNT\u0420USTED_DATA_X>")
        self.assertEqual(q2["boundary_lookalikes_escaped"], 1)

    def test_3_falsification_suite_passes_all_controls(self) -> None:
        rep = cpr.run_falsification_suite()
        self.assertTrue(rep["falsify_passed"], msg=str(rep))
        self.assertEqual(len(rep["negative_controls"]), 4)
        for ctrl in rep["negative_controls"]:
            self.assertTrue(ctrl["passed"], msg=str(ctrl))


if __name__ == "__main__":
    unittest.main()
