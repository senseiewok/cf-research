#!/usr/bin/env python3
"""Unit and falsification tests for tools/trial_atlas/rare_genotype_trial_screener.py (T-0920)."""
from __future__ import annotations

import sys
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import rare_genotype_trial_screener as rgts  # noqa: E402


class TestRareGenotypeTrialScreener(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        _, cls.studies = rgts.load_prisma_included_studies()
        cls.known_ncts = {str(s.get("nct_id")) for s in cls.studies if s.get("nct_id")}

    def test_1_screens_genotype_agnostic_and_class_i_trials_with_t1_provenance(self) -> None:
        rows = rgts.screen_rare_genotype_trials(
            self.studies,
            variant_class="genotype_agnostic_mrna_gene",
            status="ALL",
            limit=25,
        )
        self.assertGreaterEqual(len(rows), 1)
        for r in rows:
            self.assertIn(r["nct_id"], self.known_ncts)
            self.assertEqual(r["provenance_tier"], "T1")
            self.assertIn("genotype_agnostic_mrna_gene", r["matched_variant_classes"])

    def test_2_bilingual_care_team_packets_pass_safety_and_provenance(self) -> None:
        rows = rgts.screen_rare_genotype_trials(
            self.studies,
            variant_class="all",
            variant_symbol="G542X",
            status="ALL",
            limit=10,
        )
        self.assertGreaterEqual(len(rows), 1)
        pkt_en = rgts.build_care_team_packet(rows, self.known_ncts, lang="en")
        pkt_es = rgts.build_care_team_packet(rows, self.known_ncts, lang="es-MX")
        self.assertEqual(pkt_en["lang"], "en")
        self.assertEqual(pkt_es["lang"], "es-MX")
        self.assertGreaterEqual(len(pkt_en["care_team_questions"]), 3)
        self.assertGreaterEqual(len(pkt_es["care_team_questions"]), 3)
        self.assertIn("not", pkt_en["clinical_boundary"].lower())
        self.assertIn("no ", pkt_es["clinical_boundary"].lower())

    def test_3_falsification_suite_passes_all_negative_controls(self) -> None:
        rep = rgts.run_falsification_suite()
        self.assertTrue(rep["falsify_passed"], msg=str(rep))
        self.assertEqual(len(rep["negative_controls"]), 3)
        for ctrl in rep["negative_controls"]:
            self.assertTrue(ctrl["passed"], msg=str(ctrl))

    def test_4_genotype_criteria_25_trial_benchmark_and_word_boundary_guard(self) -> None:
        b_rep = rgts.evaluate_genotype_criteria_benchmark(self.studies)
        self.assertEqual(b_rep["total_trials"], 25)
        self.assertEqual(b_rep["exact_matches"], 25, msg=str(b_rep["mismatches"]))
        self.assertEqual(b_rep["exact_match_accuracy"], 1.0)
        self.assertTrue(b_rep["passed"])
        for cls_name, metrics in b_rep["class_metrics"].items():
            self.assertEqual(metrics["precision"], 1.0, msg=f"{cls_name}: {metrics}")
            self.assertEqual(metrics["recall"], 1.0, msg=f"{cls_name}: {metrics}")
        # Verify 'aso' does not match inside 'tetrasodium' or 'reasons'
        self.assertFalse(rgts._keyword_matches("denufosol tetrasodium inhalation for safety reasons", "aso"))
        self.assertTrue(rgts._keyword_matches("inhaled aso therapy for 3849+10kbc->t", "aso"))


if __name__ == "__main__":
    unittest.main()
