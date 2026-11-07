"""Unit tests for Round 95 (T-0710): Terlizzi-Baiesi Kinetic Uncertainty Relation (KUR).

Primary Biophysical Grounding & Evidence:
  - DOI:10.1088/1751-8121/ab0762 (Terlizzi ID, Baiesi M: Kinetic uncertainty relation, J. Phys. A 52, 02LT03, 2019)
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

from tools.discovery.round95_terlizzi_baiesi_kur import (
    compute_terlizzi_baiesi_kur_metrics,
    evaluate_terlizzi_baiesi_kur,
    run_falsification,
)


class Round95TerlizziBaiesiKURTests(unittest.TestCase):
    """Test suite verifying Round 95 Terlizzi-Baiesi kinetic uncertainty relation bounds."""

    def setUp(self) -> None:
        self.report = evaluate_terlizzi_baiesi_kur()

    def test_task_id_and_overall_verification(self) -> None:
        self.assertEqual(self.report["task_id"], "T-0710")
        self.assertEqual(self.report["round"], 95)
        self.assertEqual(self.report["primary_doi"], "10.1088/1751-8121/ab0762")
        self.assertTrue(self.report["terlizzi_baiesi_kur_verified"])
        self.assertTrue(self.report["kur_bound_satisfied"])
        self.assertTrue(self.report["eti_kinetic_rescue_verified"])
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_four_cftr_genotypes_and_kur_bounds(self) -> None:
        profiles = self.report["gating_profiles"]
        for gid in ("wt_cftr", "rescued_cftr_eti", "uncorrected_f508del", "uncoupled_g551d"):
            self.assertIn(gid, profiles)

        wt = profiles["wt_cftr"]
        eti = profiles["rescued_cftr_eti"]
        f508 = profiles["uncorrected_f508del"]
        g551d = profiles["uncoupled_g551d"]

        self.assertTrue(wt["bound_satisfied"])
        self.assertGreaterEqual(wt["observed_relative_variance"], wt["kur_theoretical_lower_bound"])
        self.assertTrue(eti["bound_satisfied"])
        self.assertGreaterEqual(eti["observed_relative_variance"], eti["kur_theoretical_lower_bound"])
        self.assertFalse(f508["bound_satisfied"])
        self.assertFalse(g551d["bound_satisfied"])

    def test_falsification_and_invalid_inputs(self) -> None:
        fals = run_falsification()
        self.assertEqual(fals["rule"], "AH-10")
        self.assertTrue(fals["passed"])
        self.assertEqual(fals["verdict"], "PASS")
        self.assertTrue(fals["violation_detected"])

        with self.assertRaises(ValueError):
            compute_terlizzi_baiesi_kur_metrics("unknown_allele")
        with self.assertRaises(ValueError):
            compute_terlizzi_baiesi_kur_metrics("wt_cftr", activity_override=-3.0)


if __name__ == "__main__":
    unittest.main()
