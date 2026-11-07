"""Unit tests for Round 94 (T-0709): Dechant-Sasa Fluctuation-Response Inequality to Thermal Perturbations.

Primary Biophysical Grounding & Evidence:
  - DOI:10.1073/pnas.1918386117 (PMID:32220962 — Dechant A, Sasa S-i: Fluctuation-response inequality out of equilibrium, PNAS 2020)
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

from tools.discovery.round94_dechant_thermal_fri import (
    compute_dechant_thermal_fri_metrics,
    evaluate_dechant_thermal_fri,
    run_falsification,
)


class Round94DechantThermalFRITests(unittest.TestCase):
    """Test suite verifying Round 94 Dechant-Sasa thermal fluctuation-response bounds."""

    def setUp(self) -> None:
        self.report = evaluate_dechant_thermal_fri()

    def test_task_id_and_overall_verification(self) -> None:
        self.assertEqual(self.report["task_id"], "T-0709")
        self.assertEqual(self.report["round"], 94)
        self.assertEqual(self.report["primary_doi"], "10.1073/pnas.1918386117")
        self.assertEqual(self.report["primary_pmid"], "32220962")
        self.assertTrue(self.report["dechant_thermal_fri_verified"])
        self.assertTrue(self.report["thermal_fri_bound_satisfied"])
        self.assertTrue(self.report["eti_thermal_rescue_verified"])
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_four_cftr_genotypes_and_thermal_bounds(self) -> None:
        profiles = self.report["gating_profiles"]
        for gid in ("wt_cftr", "rescued_cftr_eti", "uncorrected_f508del", "uncoupled_g551d"):
            self.assertIn(gid, profiles)

        wt = profiles["wt_cftr"]
        eti = profiles["rescued_cftr_eti"]
        f508 = profiles["uncorrected_f508del"]
        g551d = profiles["uncoupled_g551d"]

        self.assertTrue(wt["bound_satisfied"])
        self.assertLessEqual(wt["response_to_fluctuation_lhs_k_inv2"], wt["thermal_fisher_bound_k_inv2"])
        self.assertTrue(eti["bound_satisfied"])
        self.assertLessEqual(eti["response_to_fluctuation_lhs_k_inv2"], eti["thermal_fisher_bound_k_inv2"])
        self.assertFalse(f508["bound_satisfied"])
        self.assertFalse(g551d["bound_satisfied"])

    def test_falsification_and_invalid_inputs(self) -> None:
        fals = run_falsification()
        self.assertEqual(fals["rule"], "AH-10")
        self.assertTrue(fals["passed"])
        self.assertEqual(fals["verdict"], "PASS")
        self.assertTrue(fals["violation_detected"])

        with self.assertRaises(ValueError):
            compute_dechant_thermal_fri_metrics("unknown_allele")
        with self.assertRaises(ValueError):
            compute_dechant_thermal_fri_metrics("wt_cftr", fisher_bound_override=-1.0)


if __name__ == "__main__":
    unittest.main()
