"""Unit tests for tools/claims/check_anti_hallucination.py."""
from __future__ import annotations

import shutil
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
import sys
sys.path.insert(0, str(HERE))

import check_anti_hallucination as cah


class TestCheckAntiHallucination(unittest.TestCase):
    def test_repo_passes_cleanly(self):
        """The actual cf-research repository must pass all anti-hallucination guardrails with 0 findings."""
        failures = cah.run_checks(REPO_ROOT)
        self.assertEqual(failures, [])

    def test_ah1_nct_provenance_rejects_unverified_nct(self):
        """Rule AH-1 must fail if an unknown or invented NCT ID is cited."""
        snapshot_ncts = {"NCT00909532", "NCT01807923"}
        clean_text = "Trial evaluated in STRIVE (NCT00909532)."
        fails_clean = cah.check_nct_provenance(Path("test.md"), clean_text, snapshot_ncts)
        self.assertEqual(fails_clean, [])

        planted_text = "Fictional trial NCT99999999 claimed miraculous efficacy."
        fails_planted = cah.check_nct_provenance(Path("test.md"), planted_text, snapshot_ncts)
        self.assertEqual(len(fails_planted), 1)
        self.assertIn("FAIL AH-1 [nct-provenance]", fails_planted[0])
        self.assertIn("NCT99999999", fails_planted[0])

    def test_ah2_prisma_cohort_consistency_rejects_unfiltered_counts(self):
        """Rule AH-2 must fail if unfiltered 1,053-cohort counts contaminate PRISMA-screened 940 cohort text."""
        # Clean text with PRISMA-included counts
        clean_text = (
            "Across 940 included studies and 1,821 primary outcome entries, "
            "164 entries were classified as lung_function_fev1 and 26 as pulmonary_exacerbations."
        )
        self.assertEqual(cah.check_prisma_cohort_consistency(Path("test.md"), clean_text), [])

        # Contaminated text with unfiltered count 169
        contaminated_text = (
            "Across 940 included studies and 1,821 primary outcome entries, "
            "169 entries were classified as lung_function_fev1."
        )
        fails = cah.check_prisma_cohort_consistency(Path("test.md"), contaminated_text)
        self.assertEqual(len(fails), 1)
        self.assertIn("FAIL AH-2 [prisma-cohort-consistency]", fails[0])
        self.assertIn("169 lung_function_fev1 (unfiltered)", fails[0])

    def test_ah3_companion_identifier_grounding(self):
        """Rule AH-3 must verify that all NCTs, DOIs, and PDB IDs in markdown are in companions."""
        md_text = (
            "Reference study NCT00909532 and structure PDB 7U70 published at "
            "DOI 10.1056/NEJM199409083311003."
        )
        # Fully grounded companions
        ev_text = "NCT00909532 PDB 7U70 10.1056/NEJM199409083311003"
        cl_text = ""
        self.assertEqual(
            cah.check_companion_identifier_grounding(Path("test.md"), md_text, ev_text, cl_text),
            []
        )

        # Missing DOI and PDB
        ev_missing = "NCT00909532"
        fails = cah.check_companion_identifier_grounding(Path("test.md"), md_text, ev_missing, cl_text)
        self.assertEqual(len(fails), 2)
        self.assertTrue(any("DOI 10.1056/NEJM199409083311003" in f for f in fails))
        self.assertTrue(any("PDB structure ID 7U70" in f for f in fails))

    def test_ah4_clinical_overclaim_guard(self):
        """Rule AH-4 must reject unhedged curative promises or clinical recommendations."""
        clean_text = "Ivacaftor demonstrated significant improvement in percent predicted FEV1."
        self.assertEqual(cah.check_clinical_overclaim_guard(Path("test.md"), clean_text), [])

        for phrase in cah.PROHIBITED_OVERCLAIM_PHRASES:
            bad_text = f"This therapeutic strategy is {phrase} in clinical cohorts."
            fails = cah.check_clinical_overclaim_guard(Path("test.md"), bad_text)
            self.assertEqual(len(fails), 1, f"Expected rejection for phrase: {phrase}")
            self.assertIn("FAIL AH-4 [clinical-overclaim-guard]", fails[0])
            self.assertIn(phrase, fails[0])

        bad_bw = "KALYDECO carries an active Boxed Warning for acute hepatic failure."
        fails_bw = cah.check_clinical_overclaim_guard(Path("test.md"), bad_bw)
        self.assertEqual(len(fails_bw), 1)
        self.assertIn("extrapolated Boxed Warning", fails_bw[0])

        good_bw = "KALYDECO (ivacaftor, NDA 203188): Boxed Warning is `not stated` (none present in SPL)."
        self.assertEqual(cah.check_clinical_overclaim_guard(Path("test.md"), good_bw), [])

    def test_ah5_discovery_module_registration_rejects_unregistered_module(self):
        """Rule AH-5 must reject an unregistered biophysical model in tools/discovery/."""
        snapshot_ncts = cah.load_snapshot_ncts(REPO_ROOT / "tools" / "trial_atlas" / "snapshot.json")
        self.assertEqual(cah.check_discovery_module_registration(REPO_ROOT, snapshot_ncts), [])

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_root = Path(tmpdir)
            disc_dst = tmp_root / "tools" / "discovery"
            shutil.copytree(REPO_ROOT / "tools" / "discovery", disc_dst)
            rogue = disc_dst / "rogue_unverified_model.py"
            rogue.write_text(
                '"""Unregistered model citing invented NCT99999999."""\ndef compute(): return 1\n',
                encoding="utf-8",
            )
            fails = cah.check_discovery_module_registration(tmp_root, snapshot_ncts)
            self.assertTrue(
                any("rogue_unverified_model" in f and "AH-5" in f for f in fails),
                f"Expected AH-5 failure for rogue_unverified_model, got: {fails}",
            )
            self.assertTrue(
                any("NCT99999999" in f and "AH-5" in f for f in fails),
                f"Expected AH-5 NCT failure for NCT99999999, got: {fails}",
            )

    def test_ah6_cite_check_ledger_provenance_rejects_unverified_identifiers(self):
        """Rule AH-6 must verify DOIs, PMIDs, and NCTs against sources/cite-check-ledger.jsonl."""
        ledger = cah.load_cite_check_ledger(REPO_ROOT / "sources" / "cite-check-ledger.jsonl")
        self.assertGreaterEqual(len(ledger), 186)

        clean_text = "Verified DOI 10.1056/NEJMoa1908639, PMID: 30661751, and NCT03525444."
        self.assertEqual(
            cah.check_cite_check_ledger_provenance(Path("test.md"), clean_text, ledger),
            [],
        )

        bad_text = "Hallucinated DOI 10.1056/NEJMoa9999999 and PMID: 99999999."
        fails = cah.check_cite_check_ledger_provenance(Path("test.md"), bad_text, ledger)
        self.assertEqual(len(fails), 2)
        self.assertTrue(any("10.1056/NEJMoa9999999" in f and "AH-6" in f for f in fails))
        self.assertTrue(any("PMID 99999999" in f and "AH-6" in f for f in fails))

    def test_main_cli_execution(self):
        """CLI main entry point must execute cleanly and return 0 for the valid repository."""
        ret = cah.main(["--root", str(REPO_ROOT)])
        self.assertEqual(ret, 0)


if __name__ == "__main__":
    unittest.main()
