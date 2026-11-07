"""Unit tests for Round 92 (T-0703): Hasegawa-Van Vu Generalized Uncertainty Relation for Biochemical Processes.

Primary Biophysical Grounding & Evidence:
  - DOI:10.1103/PhysRevLett.122.230601 (PMID:31298835 — Hasegawa Y, Van Vu T: Uncertainty Relation in Biochemical Processes, Phys. Rev. Lett. 122, 230601, 2019)
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

from tools.discovery.round92_hasegawa_biochemical_tur import (
    compute_hasegawa_genotype_metrics,
    evaluate_hasegawa_biochemical_tur,
    run_falsification,
)


class Round92HasegawaBiochemicalTURTests(unittest.TestCase):
    """Test suite verifying Round 92 Hasegawa-Van Vu biochemical TUR bounds."""

    def setUp(self) -> None:
        self.report = evaluate_hasegawa_biochemical_tur()

    def test_task_id_and_overall_verification(self) -> None:
        self.assertEqual(self.report["task_id"], "T-0703")
        self.assertEqual(self.report["round"], 92)
        self.assertEqual(self.report["primary_doi"], "10.1103/PhysRevLett.122.230601")
        self.assertEqual(self.report["primary_pmid"], "31298835")
        self.assertTrue(self.report["hasegawa_biochemical_tur_verified"])
        self.assertTrue(self.report["biochemical_bound_satisfied"])
        self.assertTrue(self.report["eti_precision_restoration_verified"])
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_four_cftr_genotypes_and_biochemical_bounds(self) -> None:
        profiles = self.report["gating_profiles"]
        for gid in ("wt_cftr", "rescued_cftr_eti", "uncorrected_f508del", "uncoupled_g551d"):
            self.assertIn(gid, profiles)

        wt = profiles["wt_cftr"]
        eti = profiles["rescued_cftr_eti"]
        f508 = profiles["uncorrected_f508del"]
        g551d = profiles["uncoupled_g551d"]

        self.assertTrue(wt["bound_satisfied"])
        self.assertAlmostEqual(wt["hasegawa_theoretical_lower_bound"], 0.064725, places=5)
        self.assertGreaterEqual(wt["observed_relative_variance"], wt["hasegawa_theoretical_lower_bound"])
        self.assertGreaterEqual(wt["precision_ratio_to_bound"], 1.0)

        self.assertTrue(eti["bound_satisfied"])
        self.assertAlmostEqual(eti["hasegawa_theoretical_lower_bound"], 0.070423, places=5)
        self.assertGreaterEqual(eti["observed_relative_variance"], eti["hasegawa_theoretical_lower_bound"])
        self.assertGreaterEqual(eti["precision_ratio_to_bound"], 1.0)

        self.assertFalse(f508["bound_satisfied"])
        self.assertLess(f508["entropy_production_sigma_kb_s"], 5.0)
        self.assertFalse(g551d["bound_satisfied"])
        self.assertLess(g551d["entropy_production_sigma_kb_s"], 5.0)

    def test_falsification_and_invalid_inputs(self) -> None:
        fals = run_falsification()
        self.assertEqual(fals["rule"], "AH-10")
        self.assertTrue(fals["passed"])
        self.assertEqual(fals["verdict"], "PASS")
        self.assertTrue(fals["violation_detected"])

        with self.assertRaises(ValueError):
            compute_hasegawa_genotype_metrics("unknown_allele")
        with self.assertRaises(ValueError):
            compute_hasegawa_genotype_metrics("wt_cftr", entropy_override=-2.0)


if __name__ == "__main__":
    unittest.main()
