"""Unit tests for Round 91 (T-0698): Dechant-Sasa Fluctuation-Response & Entropic Transport Sensitivity Bound.

Primary Biophysical Grounding & Evidence:
  - DOI:10.1073/pnas.1916060117 (PMID:32127480 — Dechant A, Sasa S: Fluctuation-response inequality out of equilibrium, PNAS 117, 6430-6436, 2020)
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

from tools.discovery.round91_dechant_sasa_fluctuation_response import (
    compute_dechant_sasa_genotype_metrics,
    evaluate_dechant_sasa_fluctuation_response,
    run_falsification,
)


class Round91DechantSasaFluctuationResponseTests(unittest.TestCase):
    """Test suite verifying Round 91 Dechant-Sasa fluctuation-response sensitivity bounds."""

    def setUp(self) -> None:
        self.report = evaluate_dechant_sasa_fluctuation_response()

    def test_task_id_and_overall_verification(self) -> None:
        self.assertEqual(self.report["task_id"], "T-0698")
        self.assertEqual(self.report["round"], 91)
        self.assertEqual(self.report["primary_doi"], "10.1073/pnas.1916060117")
        self.assertEqual(self.report["primary_pmid"], "32127480")
        self.assertTrue(self.report["dechant_sasa_fluctuation_response_verified"])
        self.assertTrue(self.report["entropic_sensitivity_bound_verified"])
        self.assertTrue(self.report["eti_sensitivity_restoration_verified"])
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_four_cftr_genotypes_and_sensitivity_bounds(self) -> None:
        profiles = self.report["gating_profiles"]
        for gid in ("wt_cftr", "rescued_cftr_eti", "uncorrected_f508del", "uncoupled_g551d"):
            self.assertIn(gid, profiles)

        wt = profiles["wt_cftr"]
        eti = profiles["rescued_cftr_eti"]
        f508 = profiles["uncorrected_f508del"]
        g551d = profiles["uncoupled_g551d"]

        self.assertTrue(wt["bound_satisfied"])
        self.assertAlmostEqual(wt["dechant_sasa_sensitivity_upper_bound"], 3.405877, places=5)
        self.assertLessEqual(wt["observed_transport_sensitivity_chi"], wt["dechant_sasa_sensitivity_upper_bound"])
        self.assertLessEqual(wt["fluctuation_response_ratio"], 1.0)

        self.assertTrue(eti["bound_satisfied"])
        self.assertAlmostEqual(eti["dechant_sasa_sensitivity_upper_bound"], 3.146427, places=5)
        self.assertLessEqual(eti["observed_transport_sensitivity_chi"], eti["dechant_sasa_sensitivity_upper_bound"])
        self.assertLessEqual(eti["fluctuation_response_ratio"], 1.0)

        self.assertFalse(f508["bound_satisfied"])
        self.assertLess(f508["observed_transport_sensitivity_chi"], 0.5)
        self.assertFalse(g551d["bound_satisfied"])
        self.assertLess(g551d["observed_transport_sensitivity_chi"], 0.5)

    def test_falsification_and_invalid_inputs(self) -> None:
        fals = run_falsification()
        self.assertEqual(fals["rule"], "AH-10")
        self.assertTrue(fals["passed"])
        self.assertEqual(fals["verdict"], "PASS")
        self.assertLess(fals["low_entropy_upper_bound"], 0.5)

        with self.assertRaises(ValueError):
            compute_dechant_sasa_genotype_metrics("unknown_allele")
        with self.assertRaises(ValueError):
            compute_dechant_sasa_genotype_metrics("wt_cftr", entropy_override=-1.0)


if __name__ == "__main__":
    unittest.main()
