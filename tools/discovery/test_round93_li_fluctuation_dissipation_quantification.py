"""Unit tests for Round 93 (T-0704): Li-Horowitz-Gingrich-Fakhri Dissipation Quantification from Fluctuations.

Primary Biophysical Grounding & Evidence:
  - DOI:10.1016/j.cels.2019.03.007 (PMID:30904378 — Li J, Horowitz JM, Gingrich TR, Fakhri N: Quantifying Dissipation Using Fluctuations in Living Systems, Cell Syst. 8, 280-286, 2019)
  - DOI:10.1126/science.ade2216 (Fiedorczuk & Chen: Molecular structures reveal synergistic rescue of dF508 CFTR by Trikafta modulators, Science 2022)
  - PDB 6MSM, 5UAK, 8EIQ
"""
from __future__ import annotations

from pathlib import Path
import sys
import unittest

_REPO_ROOT = Path(__file__).resolve().parent.parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from tools.discovery.round93_li_fluctuation_dissipation_quantification import (
    compute_li_dissipation_metrics,
    evaluate_li_dissipation_quantification,
    run_falsification,
)


class Round93LiDissipationQuantificationTests(unittest.TestCase):
    """Test suite verifying Round 93 Li-Horowitz-Gingrich-Fakhri dissipation quantification bounds."""

    def setUp(self) -> None:
        self.report = evaluate_li_dissipation_quantification()

    def test_task_id_and_overall_verification(self) -> None:
        self.assertEqual(self.report["task_id"], "T-0704")
        self.assertEqual(self.report["round"], 93)
        self.assertEqual(self.report["primary_doi"], "10.1016/j.cels.2019.03.007")
        self.assertEqual(self.report["primary_pmid"], "30904378")
        self.assertTrue(self.report["li_dissipation_quantification_verified"])
        self.assertTrue(self.report["dissipation_lower_bound_satisfied"])
        self.assertTrue(self.report["eti_dissipation_restoration_verified"])
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_four_cftr_genotypes_and_dissipation_bounds(self) -> None:
        profiles = self.report["gating_profiles"]
        for gid in ("wt_cftr", "rescued_cftr_eti", "uncorrected_f508del", "uncoupled_g551d"):
            self.assertIn(gid, profiles)

        wt = profiles["wt_cftr"]
        eti = profiles["rescued_cftr_eti"]
        f508 = profiles["uncorrected_f508del"]
        g551d = profiles["uncoupled_g551d"]

        self.assertTrue(wt["bound_satisfied"])
        self.assertAlmostEqual(wt["dissipation_recovery_ratio"], 0.84, places=2)
        self.assertLessEqual(wt["reconstructed_lower_bound_sigma_kb_s"], wt["true_dissipation_sigma_tot_kb_s"])
        self.assertGreaterEqual(wt["dissipation_recovery_ratio"], 0.75)

        self.assertTrue(eti["bound_satisfied"])
        self.assertAlmostEqual(eti["dissipation_recovery_ratio"], 0.819697, places=5)
        self.assertLessEqual(eti["reconstructed_lower_bound_sigma_kb_s"], eti["true_dissipation_sigma_tot_kb_s"])
        self.assertGreaterEqual(eti["dissipation_recovery_ratio"], 0.75)

        self.assertFalse(f508["bound_satisfied"])
        self.assertLess(f508["dissipation_recovery_ratio"], 0.20)
        self.assertFalse(g551d["bound_satisfied"])
        self.assertLess(g551d["dissipation_recovery_ratio"], 0.20)

    def test_falsification_and_invalid_inputs(self) -> None:
        fals = run_falsification()
        self.assertEqual(fals["rule"], "AH-10")
        self.assertTrue(fals["passed"])
        self.assertEqual(fals["verdict"], "PASS")
        self.assertTrue(fals["second_law_violation_detected"])

        with self.assertRaises(ValueError):
            compute_li_dissipation_metrics("unknown_allele")
        with self.assertRaises(ValueError):
            compute_li_dissipation_metrics("wt_cftr", sigma_override=-5.0)


if __name__ == "__main__":
    unittest.main()
