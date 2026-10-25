#!/usr/bin/env python3
"""Unit tests for Counterfactual Clinical Trial Failure Forensics (T-0428)."""
from __future__ import annotations

from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parent))

from counterfactual_trial_forensics import (
    HISTORICAL_TRIAL_COHORT,
    run_counterfactual_trial_forensics,
)


class TestCounterfactualTrialForensics(unittest.TestCase):
    """Verify historical trial UNSAT decomposition and AlphaEvolve counterfactual SAT rescue."""

    def test_historical_trials_fail_invariants_and_rescued_pass(self) -> None:
        report = run_counterfactual_trial_forensics()
        self.assertEqual(report["task_id"], "T-0428")
        self.assertEqual(report["forensic_trials_count"], 4)
        self.assertTrue(report["all_historical_unsat"])
        self.assertTrue(report["all_rescued_sat"])

        for rec in report["forensic_records"]:
            self.assertFalse(rec["historical_eval"]["invariants_passed"])
            self.assertGreater(len(rec["historical_eval"]["unsat_core"]), 0)
            self.assertTrue(rec["rescued_eval"]["invariants_passed"])
            self.assertEqual(rec["rescued_eval"]["unsat_core"], [])
            self.assertGreaterEqual(rec["rescued_eval"]["delta_isc_pct_wt"], 15.0)
            self.assertGreaterEqual(rec["rescued_eval"]["asl_depth_um"], 7.0)
            self.assertGreater(rec["delta_wt_rescue_gain_pct"], 14.0)

    def test_mrt5005_mucus_pore_steric_trap_resolved(self) -> None:
        report = run_counterfactual_trial_forensics()
        mrt = next(r for r in report["forensic_records"] if r["nct_id"] == "NCT03375047")
        self.assertLess(mrt["historical_eval"]["pore_size_at_vector_nm"], 72.0)
        self.assertGreaterEqual(mrt["rescued_eval"]["pore_size_at_vector_nm"], 100.0)

    def test_prisma_cohort_numbers_locked(self) -> None:
        report = run_counterfactual_trial_forensics()
        self.assertEqual(report["prisma_context"]["screened_records"], 1053)
        self.assertEqual(report["prisma_context"]["included_studies"], 940)
        self.assertEqual(report["prisma_context"]["outcome_entries"], 1821)
        self.assertEqual(len(HISTORICAL_TRIAL_COHORT), 4)


if __name__ == "__main__":
    unittest.main()
