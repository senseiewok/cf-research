#!/usr/bin/env python3
"""Unit, benchmark, edge-case, and falsification tests for tools/trial_atlas/bilingual_care_rubric.py."""
from __future__ import annotations

import sys
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import bilingual_care_rubric as bcr  # noqa: E402
import whole_person_visit_prep as wpvp  # noqa: E402


class TestBilingualCareRubric(unittest.TestCase):
    def test_1_benchmark_meets_pred_0008_thresholds(self) -> None:
        bm = bcr.run_benchmark()
        self.assertEqual(bm["prediction_id"], "PRED-0008")
        self.assertEqual(bm["total_cases"], 16)
        self.assertGreaterEqual(bm["accuracy"], 0.90, msg=str(bm))
        self.assertGreaterEqual(bm["violation_precision"], 0.90, msg=str(bm))
        self.assertGreaterEqual(bm["violation_recall"], 0.90, msg=str(bm))

    def test_2_live_lab_guides_audit_and_honest_gaps(self) -> None:
        live = bcr.evaluate_live_lab_guides()
        self.assertEqual(live["snapshot_prisma_included"], 940)
        self.assertEqual(live["cff_research_priorities_leaks"], [])
        self.assertTrue(live["spl_boxed_warning_ground_truth"]["TRIKAFTA"]["has_boxed_warning"])
        self.assertTrue(live["spl_boxed_warning_ground_truth"]["ALYFTREK"]["has_boxed_warning"])
        self.assertFalse(live["spl_boxed_warning_ground_truth"]["KALYDECO"]["has_boxed_warning"])
        self.assertFalse(live["spl_boxed_warning_ground_truth"]["ORKAMBI"]["has_boxed_warning"])
        self.assertFalse(live["spl_boxed_warning_ground_truth"]["SYMDEKO"]["has_boxed_warning"])

        by_id = {g["guide_id"]: g for g in live["guides"]}
        self.assertTrue(by_id["whole_person_visit_prep:all:en"]["rubric_pass"])
        self.assertTrue(by_id["care_team_trials:RECRUITING:limit=8:en"]["rubric_pass"])
        self.assertTrue(by_id["care_team_trials:RECRUITING:limit=8:es-MX"]["rubric_pass"])
        # Honest readability & gap reporting:
        self.assertGreater(by_id["whole_person_visit_prep:all:en"]["readability"]["flesch_kincaid_grade"], 8.0)
        self.assertFalse(by_id["whole_person_visit_prep:all:en"]["readability"]["within_grade_8"])
        # Single-domain es-MX hepatobiliary guide omits interpreter guidance; domain coverage must not false-hit on 'fast'.
        self.assertIn(
            "es_mx_equity:missing_qualified_interpreter",
            by_id["whole_person_visit_prep:hepatobiliary_dili:es-MX"]["failed_rules"],
        )
        self.assertNotIn(
            "hepatobiliary_dili",
            by_id["whole_person_visit_prep:newborn_screening_equity:en"]["domains_covered"],
        )

    def test_3_ah4_debray_and_nct_edge_cases(self) -> None:
        # AH-4 bypasses must be flagged:
        for bad in (
            "A Boxed Warning applies to KALYDECO, ORKAMBI, and SYMDEKO.",
            "La Advertencia de Recuadro aplica a KALYDECO, ORKAMBI y SYMDEKO.",
            "KALYDECO and ORKAMBI carry Boxed Warnings for liver injury.",
            "KALYDECO has a Black Box Warning for liver injury.",
            "KALYDECO (see Debray et al. 2026) carries a Boxed Warning.",
            "KALYDECO carries a Boxed Warning (exceptions: none).",
        ):
            self.assertTrue(bcr.check_ah4_boxed_warning_isolation(bad), msg=f"Missed AH-4 bypass: {bad}")

        # AH-4 true negations must pass:
        for good in (
            "KALYDECO has no Boxed Warning.",
            "KALYDECO no tiene Advertencia de Recuadro.",
            "For ORKAMBI and SYMDEKO, a Boxed Warning is not stated in the SPL.",
        ):
            self.assertEqual(bcr.check_ah4_boxed_warning_isolation(good), [], msg=f"False AH-4 flag: {good}")

        # Debray CC BY-NC-ND 4.0 et al. & DOI checks:
        self.assertTrue(
            bcr.check_license_and_pointer_boundaries("Debray et al. (2026) adapted flowchart and reproduced table.")
        )
        self.assertTrue(
            bcr.check_license_and_pointer_boundaries("Bulk XML excerpts from 10.1097/HC9.0000000000001061.")
        )

        # Substring 'nct' in 'function' / 'conjunction' must NOT trigger NLM disclaimer:
        plain = (
            "Liver function and lung function checks are discussed in conjunction with your CF care team. "
            "If worrying, call emergency services (911 in the US). This is not a diagnosis."
        )
        ev_plain = bcr.evaluate_guide_text(plain, lang="en")
        self.assertTrue(ev_plain["rubric_pass"], msg=str(ev_plain))

        # Protective Spanish conditional must NOT trigger directive:
        base_es = wpvp.render_whole_person_markdown(wpvp.build_whole_person_packet(domain="newborn_screening_equity", lang="es-MX"))
        ev_cond = bcr.evaluate_guide_text(
            base_es + "\n\n[T1] No podemos decirle si usted es elegible; consulte a su equipo.",
            lang="es-MX",
        )
        self.assertNotIn("bilingual_safety:directive", ev_cond["failed_rules"])

    def test_4_falsification_negative_controls_all_pass(self) -> None:
        rep = bcr.run_falsification_suite()
        self.assertTrue(rep["falsify_passed"], msg=str(rep))
        self.assertEqual(len(rep["negative_controls"]), 4)
        for ctrl in rep["negative_controls"]:
            self.assertTrue(ctrl["passed"], msg=str(ctrl))

    def test_5_cli_check_and_falsify_exit_zero(self) -> None:
        self.assertEqual(bcr.main(["--check"]), 0)
        self.assertEqual(bcr.main(["--falsify"]), 0)
        self.assertEqual(bcr.main(["--benchmark"]), 0)


if __name__ == "__main__":
    unittest.main()
