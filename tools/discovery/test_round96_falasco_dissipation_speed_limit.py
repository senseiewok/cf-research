"""Unit tests for Round 96 (T-0715): Falasco-Esposito Dissipation-Time Speed Limit.

Primary Biophysical Grounding & Evidence:
  - DOI:10.1103/PhysRevLett.125.120604 (Falasco G, Esposito M: Dissipation-time uncertainty relation, Phys. Rev. Lett. 125, 120604, 2020)
  - DOI:10.1085/jgp.200609558 (Csanady L et al., J. Gen. Physiol. 2006)
  - PDB 6MSM, 5UAK, 8EIQ
"""
from __future__ import annotations

from pathlib import Path
import sys
import unittest

_REPO_ROOT = Path(__file__).resolve().parent.parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from tools.discovery.round96_falasco_dissipation_speed_limit import (
    compute_falasco_dissipation_speed_limit_metrics,
    evaluate_falasco_dissipation_speed_limit,
    run_falsification,
)


class Round96FalascoDissipationSpeedLimitTests(unittest.TestCase):
    """Test suite verifying Round 96 Falasco-Esposito dissipation-time speed limit bounds."""

    def setUp(self) -> None:
        self.report = evaluate_falasco_dissipation_speed_limit()

    def test_task_id_and_overall_verification(self) -> None:
        self.assertEqual(self.report["task_id"], "T-0715")
        self.assertEqual(self.report["round"], 96)
        self.assertEqual(self.report["primary_doi"], "10.1103/PhysRevLett.125.120604")
        self.assertTrue(self.report["falasco_dissipation_speed_limit_verified"])
        self.assertTrue(self.report["speed_limit_bound_satisfied"])
        self.assertTrue(self.report["eti_speed_limit_rescue_verified"])
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_four_cftr_genotypes_and_speed_limit_bounds(self) -> None:
        profiles = self.report["gating_profiles"]
        for gid in ("wt_cftr", "rescued_cftr_eti", "uncorrected_f508del", "uncoupled_g551d"):
            self.assertIn(gid, profiles)

        wt = profiles["wt_cftr"]
        eti = profiles["rescued_cftr_eti"]
        f508 = profiles["uncorrected_f508del"]
        g551d = profiles["uncoupled_g551d"]

        self.assertTrue(wt["bound_satisfied"])
        self.assertGreaterEqual(wt["observed_transition_time_s"], wt["minimum_transition_time_bound_s"])
        self.assertTrue(eti["bound_satisfied"])
        self.assertGreaterEqual(eti["observed_transition_time_s"], eti["minimum_transition_time_bound_s"])
        self.assertFalse(f508["bound_satisfied"])
        self.assertFalse(g551d["bound_satisfied"])

    def test_falsification_and_invalid_inputs(self) -> None:
        fals = run_falsification()
        self.assertEqual(fals["rule"], "AH-10")
        self.assertTrue(fals["passed"])
        self.assertEqual(fals["verdict"], "PASS")
        self.assertTrue(fals["violation_detected"])

        with self.assertRaises(ValueError):
            compute_falasco_dissipation_speed_limit_metrics("unknown_allele")
        with self.assertRaises(ValueError):
            compute_falasco_dissipation_speed_limit_metrics("wt_cftr", entropy_rate_override=-2.0)


if __name__ == "__main__":
    unittest.main()
