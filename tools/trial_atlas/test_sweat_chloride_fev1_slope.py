#!/usr/bin/env python3
"""Unit tests for Trial Endpoint Atlas sweat chloride vs FEV1 slope integration.

Verifies:
1. Deterministic analytical ODE solutions for sweat chloride and FEV1 trajectories.
2. Systemic modulator concordance (r >= 0.65) between sweat chloride reduction and ppFEV1 improvement.
3. Airway-restricted compartmental decoupling (|r| < 0.25) where dermal sweat chloride does not respond.
4. Trial Endpoint Atlas integration across 1,053 studies from snapshot.json.
5. Biophysical negative controls, parameter monotonicity, and NLM disclaimer compliance.
"""

from __future__ import annotations

import math
import sys
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import sweat_chloride_fev1_slope as sc_fev  # noqa: E402


class TestSweatChlorideFEV1Slope(unittest.TestCase):
    """Test suite for sweat chloride rate of change vs FEV1 slope modeling."""

    def setUp(self) -> None:
        self.sw_params = sc_fev.SweatChlorideParams(
            baseline_sweat=102.3,
            target_sweat=45.0,
            k_sweat=0.15,
        )
        self.f1_params = sc_fev.FEV1TrajectoryParams(
            baseline_fev1=61.6,
            target_fev1=76.5,
            k_fev1=0.08,
            mu_decline=1.5 / 365.25,
        )

    def test_sweat_chloride_analytical_solution(self) -> None:
        """Verifies baseline, plateau, and derivative consistency of sweat chloride ODE."""
        # Baseline at t = 0
        cl_0, rate_0 = sc_fev.compute_sweat_chloride(0.0, self.sw_params)
        self.assertAlmostEqual(cl_0, 102.3, places=4)
        # Derivative at t = 0: -k * (102.3 - 45.0) = -0.15 * 57.3 = -8.595 mmol/L/day
        self.assertAlmostEqual(rate_0, -0.15 * (102.3 - 45.0), places=4)

        # Plateau as t -> infinity
        cl_inf, rate_inf = sc_fev.compute_sweat_chloride(300.0, self.sw_params)
        self.assertAlmostEqual(cl_inf, 45.0, places=2)
        self.assertAlmostEqual(rate_inf, 0.0, places=2)

        # Numerical derivative check at day 14
        t = 14.0
        h = 1e-5
        cl_t, analytic_rate = sc_fev.compute_sweat_chloride(t, self.sw_params)
        cl_plus, _ = sc_fev.compute_sweat_chloride(t + h, self.sw_params)
        num_rate = (cl_plus - cl_t) / h
        self.assertAlmostEqual(analytic_rate, num_rate, places=4)

    def test_fev1_analytical_solution(self) -> None:
        """Verifies baseline, asymptote, and derivative consistency of FEV1 ODE."""
        # Baseline at t = 0
        f1_0, slope_0 = sc_fev.compute_fev1(0.0, self.f1_params)
        self.assertAlmostEqual(f1_0, 61.6, places=4)

        # Asymptote FEV1_inf = target - mu / k
        expected_inf = 76.5 - (1.5 / 365.25) / 0.08
        f1_inf, slope_inf = sc_fev.compute_fev1(300.0, self.f1_params)
        self.assertAlmostEqual(f1_inf, expected_inf, places=2)
        self.assertAlmostEqual(slope_inf, 0.0, places=2)

        # Numerical derivative check at day 28
        t = 28.0
        h = 1e-5
        f1_t, analytic_slope = sc_fev.compute_fev1(t, self.f1_params)
        f1_plus, _ = sc_fev.compute_fev1(t + h, self.f1_params)
        num_slope = (f1_plus - f1_t) / h
        self.assertAlmostEqual(analytic_slope, num_slope, places=4)

    def test_trajectory_simulation_points(self) -> None:
        """Tests multi-point time series trajectory generation."""
        days = [0.0, 14.0, 28.0, 56.0, 112.0, 168.0]
        points = sc_fev.simulate_patient_trajectory(days, self.sw_params, self.f1_params)
        self.assertEqual(len(points), len(days))
        for p in points:
            self.assertGreaterEqual(p.sweat_chloride, 45.0)
            self.assertLessEqual(p.sweat_chloride, 102.3)
            self.assertGreaterEqual(p.fev1, 61.6)

    def test_systemic_modulator_cohort_concordance(self) -> None:
        """Verifies that systemic modulator therapy achieves strong concordance (r >= 0.65)."""
        cohort = sc_fev.simulate_cohort_concordance("systemic_modulator", n_patients=120, seed=101)
        self.assertTrue(cohort.is_concordant)
        self.assertGreaterEqual(cohort.pearson_r, 0.65)
        self.assertGreaterEqual(cohort.r_squared, 0.42)
        self.assertLess(cohort.mean_delta_sweat, -40.0)  # Sweat chloride reduced by > 40 mmol/L
        self.assertGreater(cohort.mean_delta_fev1, 10.0)  # ppFEV1 improved by > 10%
        self.assertIn("concordance", cohort.clinical_implication.lower())

    def test_airway_restricted_compartmental_decoupling(self) -> None:
        """Verifies that airway-restricted delivery results in compartmental decoupling (|r| < 0.25)."""
        cohort = sc_fev.simulate_cohort_concordance("airway_restricted", n_patients=120, seed=202)
        self.assertTrue(cohort.is_decoupled)
        self.assertLess(abs(cohort.pearson_r), 0.25)
        # Sweat chloride change in dermal glands should be near zero (biological noise only)
        self.assertLess(abs(cohort.mean_delta_sweat), 3.0)
        # Pulmonary function still demonstrates positive therapeutic rescue
        self.assertGreater(cohort.mean_delta_fev1, 8.0)
        self.assertIn("decoupling", cohort.clinical_implication.lower())

    def test_placebo_control_negative_control(self) -> None:
        """Verifies placebo negative control: no sweat chloride improvement, chronic FEV1 decline."""
        cohort = sc_fev.simulate_cohort_concordance("placebo_control", n_patients=100, seed=303)
        self.assertLess(abs(cohort.mean_delta_sweat), 2.5)
        self.assertLess(cohort.mean_delta_fev1, 0.0)  # Negative change due to chronic decline

    def test_invalid_modality_raises(self) -> None:
        """Verifies fail-closed error on unknown modality."""
        with self.assertRaises(ValueError):
            sc_fev.simulate_cohort_concordance("unsupported_modality")

    def test_negative_time_raises(self) -> None:
        """Verifies validation on negative time inputs."""
        with self.assertRaises(ValueError):
            sc_fev.compute_sweat_chloride(-5.0, self.sw_params)
        with self.assertRaises(ValueError):
            sc_fev.compute_fev1(-1.0, self.f1_params)

    def test_trial_atlas_snapshot_integration(self) -> None:
        """Verifies integration with Trial Endpoint Atlas snapshot dataset."""
        stats = sc_fev.analyze_trial_atlas_endpoints()
        self.assertEqual(stats["total_studies_screened"], 1053)
        self.assertGreater(stats["fev1_primary_outcome_studies"], 0)
        self.assertGreater(stats["sweat_chloride_primary_outcome_studies"], 0)
        self.assertGreater(stats["co_occurring_dual_endpoint_studies"], 0)
        self.assertIn("ClinicalTrials.gov does not evaluate safety or efficacy", stats["nlm_disclaimer"])
        self.assertGreater(len(stats["co_occurring_sample_ncts"]), 0)
        for nct in stats["co_occurring_sample_ncts"]:
            self.assertTrue(nct.startswith("NCT"))

    def test_falsification_invariants(self) -> None:
        """Verifies all biophysical and mathematical invariant assertions."""
        invariants = sc_fev.verify_falsification_invariants()
        self.assertGreaterEqual(len(invariants), 5)
        for name, passed in invariants.items():
            self.assertTrue(passed, f"Invariant '{name}' failed validation")


if __name__ == "__main__":
    unittest.main()
