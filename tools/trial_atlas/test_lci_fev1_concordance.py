#!/usr/bin/env python3
"""Unit tests for Longitudinal Concordance of LCI2.5 vs FEV1 Model (T-0538)."""
from __future__ import annotations

import json
import math
import subprocess
import sys
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from tools.trial_atlas.lci_fev1_concordance import (
    GROUNDED_GENETIC_TRIALS,
    calculate_lci,
    calculate_standardized_z_score,
    self_test,
    simulate_longitudinal_concordance,
)


class TestLciCalculation(unittest.TestCase):
    """Test basic biophysical LCI calculations."""

    def test_standard_lci_ratio(self) -> None:
        lci = calculate_lci(cumulative_expired_volume_l=25.0, functional_residual_capacity_l=2.5)
        self.assertAlmostEqual(lci, 10.0, places=2)

    def test_invalid_frc_raises_value_error(self) -> None:
        with self.assertRaises(ValueError):
            calculate_lci(25.0, 0.0)
        with self.assertRaises(ValueError):
            calculate_lci(25.0, -1.5)

    def test_invalid_expired_volume_raises_value_error(self) -> None:
        with self.assertRaises(ValueError):
            calculate_lci(0.0, 2.5)


class TestZScoreCalculation(unittest.TestCase):
    """Test standardized score deviations."""

    def test_z_score_basic(self) -> None:
        z = calculate_standardized_z_score(observed=7.5, mean_baseline=6.4, sd=0.55)
        self.assertAlmostEqual(z, 2.0, places=1)

    def test_invalid_sd_raises_value_error(self) -> None:
        with self.assertRaises(ValueError):
            calculate_standardized_z_score(7.5, 6.4, 0.0)


class TestLongitudinalConcordanceSimulation(unittest.TestCase):
    """Test longitudinal trajectory and clinical sensitivity ratio."""

    def test_mild_stratum_exceeds_sensitivity_benchmark(self) -> None:
        res = simulate_longitudinal_concordance(
            disease_stratum="mild_preserved_spirometry",
            trial_id="NCT05668741",
        )
        self.assertEqual(res.disease_stratum, "mild_preserved_spirometry")
        self.assertEqual(res.trial_id, "NCT05668741")
        self.assertTrue(res.ratio_exceeds_clinical_benchmark)
        self.assertGreaterEqual(res.mean_sensitivity_ratio, 1.8)
        self.assertIsNotNone(res.earliest_lci_detection_week)
        self.assertGreater(res.lead_time_advantage_weeks, 0.0)
        self.assertEqual(len(res.trajectory), 6)

    def test_moderate_stratum_simulation(self) -> None:
        res = simulate_longitudinal_concordance(
            disease_stratum="moderate_airflow_obstruction",
            trial_id="NCT05248230",
        )
        self.assertEqual(res.disease_stratum, "moderate_airflow_obstruction")
        self.assertEqual(res.trial_id, "NCT05248230")
        self.assertTrue(res.ratio_exceeds_clinical_benchmark)
        self.assertGreaterEqual(res.mean_sensitivity_ratio, 1.8)

    def test_unknown_trial_raises_key_error(self) -> None:
        with self.assertRaises(KeyError):
            simulate_longitudinal_concordance(trial_id="NCT99999999")

    def test_unknown_stratum_raises_value_error(self) -> None:
        with self.assertRaises(ValueError):
            simulate_longitudinal_concordance(disease_stratum="unknown_cohort")

    def test_self_test_function(self) -> None:
        ret = self_test()
        self.assertEqual(ret, 0)

    def test_cli_execution_and_json(self) -> None:
        cmd = [sys.executable, str(HERE / "lci_fev1_concordance.py"), "--trial", "NCT05668741", "--json"]
        proc = subprocess.run(cmd, capture_output=True, text=True, check=True)
        data = json.loads(proc.stdout)
        self.assertEqual(data["trial_id"], "NCT05668741")
        self.assertTrue(data["ratio_exceeds_clinical_benchmark"])


if __name__ == "__main__":
    unittest.main()
