"""Unit tests for Round 97 (T-0716): Koyuk-Seifert Periodic Driving Thermodynamic Uncertainty Relation (TUR).

Primary Biophysical Grounding & Evidence:
  - DOI:10.1103/PhysRevLett.125.260604 (Koyuk T, Seifert U: Thermodynamic Uncertainty Relation for Time-Dependent Driving, Phys. Rev. Lett. 125, 260604, 2020)
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

from tools.discovery.round97_koyuk_seifert_periodic_tur import (
    compute_koyuk_seifert_periodic_tur_metrics,
    evaluate_koyuk_seifert_periodic_tur,
    run_falsification,
)


class Round97KoyukSeifertPeriodicTURTests(unittest.TestCase):
    """Test suite verifying Round 97 Koyuk-Seifert periodic driving TUR bounds."""

    def setUp(self) -> None:
        self.report = evaluate_koyuk_seifert_periodic_tur()

    def test_task_id_and_overall_verification(self) -> None:
        self.assertEqual(self.report["task_id"], "T-0716")
        self.assertEqual(self.report["round"], 97)
        self.assertEqual(self.report["primary_doi"], "10.1103/PhysRevLett.125.260604")
        self.assertTrue(self.report["koyuk_seifert_periodic_tur_verified"])
        self.assertTrue(self.report["periodic_tur_bound_satisfied"])
        self.assertTrue(self.report["eti_periodic_clock_rescue_verified"])
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_four_cftr_genotypes_and_periodic_tur_bounds(self) -> None:
        profiles = self.report["gating_profiles"]
        for gid in ("wt_cftr", "rescued_cftr_eti", "uncorrected_f508del", "uncoupled_g551d"):
            self.assertIn(gid, profiles)

        wt = profiles["wt_cftr"]
        eti = profiles["rescued_cftr_eti"]
        f508 = profiles["uncorrected_f508del"]
        g551d = profiles["uncoupled_g551d"]

        self.assertTrue(wt["bound_satisfied"])
        self.assertGreaterEqual(wt["observed_relative_variance"], wt["periodic_tur_lower_bound"])
        self.assertTrue(eti["bound_satisfied"])
        self.assertGreaterEqual(eti["observed_relative_variance"], eti["periodic_tur_lower_bound"])
        self.assertFalse(f508["bound_satisfied"])
        self.assertFalse(g551d["bound_satisfied"])

    def test_falsification_and_invalid_inputs(self) -> None:
        fals = run_falsification()
        self.assertEqual(fals["rule"], "AH-10")
        self.assertTrue(fals["passed"])
        self.assertEqual(fals["verdict"], "PASS")
        self.assertTrue(fals["violation_detected"])

        with self.assertRaises(ValueError):
            compute_koyuk_seifert_periodic_tur_metrics("unknown_allele")
        with self.assertRaises(ValueError):
            compute_koyuk_seifert_periodic_tur_metrics("wt_cftr", cycle_entropy_override=-4.0)


if __name__ == "__main__":
    unittest.main()
