#!/usr/bin/env python3
"""Unit tests for tools/curvature-screen/benchmark/whole_sample/run_on_patches.py.

Verifies acceptance criteria for T-0923:
1. 0.36 / f_c default resolution calculation when filter_cutoff is supplied without dead_time.
2. Instrument-calibrated false-alarm threshold for dead time, filtering (f_c), and burst criterion (t_crit).
3. Calibration ensures instrument settings cannot fake a rejection on passive mixture / exponential gating.
4. Two-step sequential process correctly rejects mixture under calibrated threshold.
5. Error handling for negative dead time, filter cut-off, and burst criterion.
"""
import math
import os
import sys
import unittest
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
WHOLE_SAMPLE = os.path.join(HERE, "whole_sample")
sys.path.insert(0, WHOLE_SAMPLE)

import run_on_patches  # noqa: E402


class TestDeadTimeAndFilterCutoff(unittest.TestCase):
    """Test resolution derivation and parameter validations."""

    def test_filter_cutoff_sets_default_dead_time(self):
        rng = np.random.default_rng(42)
        table = {"test": {"p0": list(rng.exponential(1.0, 100))}}
        # 1000 Hz filter cut-off gives 0.36 / 1000 = 0.00036 s dead time
        res = run_on_patches.run(table, filter_cutoff=1000)
        self.assertIn("test", res["constructs"])
        self.assertAlmostEqual(res["constructs"]["test"]["resolution_s"], 0.00036, places=7)

    def test_explicit_dead_time_preserved(self):
        rng = np.random.default_rng(42)
        table = {"test": {"p0": list(rng.exponential(1.0, 100) + 0.002)}}
        res = run_on_patches.run(table, dead_time=0.002, filter_cutoff=500)
        self.assertEqual(res["constructs"]["test"]["resolution_s"], 0.002)

    def test_invalid_parameters_raise(self):
        table = {"test": {"p0": [0.1] * 60}}
        with self.assertRaises(ValueError):
            run_on_patches.run(table, dead_time=-0.001)
        with self.assertRaises(ValueError):
            run_on_patches.run(table, filter_cutoff=-100)
        with self.assertRaises(ValueError):
            run_on_patches.run(table, dead_time=0.002, burst_criterion=-0.001)


class TestInstrumentCalibratedNull(unittest.TestCase):
    """Test instrument-calibrated threshold calibration under dead time, filtering, and burst criterion."""

    def setUp(self):
        self.null = run_on_patches.Null(seed=20261113)

    def test_standard_exponential_null_reproducible(self):
        d1 = self.null.draws(100)
        d2 = self.null.draws(100)
        np.testing.assert_array_equal(d1, d2)
        self.assertGreater(d1.mean(), 0.95)

    def test_instrument_calibrated_null_generation(self):
        # Draws with dead time, filter cut-off, and burst criterion
        d_cal = self.null.draws(100, dead_time=0.002, filter_cutoff=1000, burst_criterion=0.005)
        self.assertEqual(len(d_cal), run_on_patches.REPS)
        # Critical value should be a valid float below 1.0
        crit = float(np.quantile(d_cal, run_on_patches.ALPHA))
        self.assertGreater(crit, 0.5)
        self.assertLess(crit, 0.8)

    def test_instrument_settings_do_not_fake_rejection_on_passive_channel(self):
        """Simulate passive single-channel bursts through filtering and burst criterion.

        Under a passive channel with dead time, filter cut-off, and burst criterion,
        the test must NOT fake a rejection (false alarm rate <= 2%).
        """
        rng = np.random.default_rng(101)
        td = 0.002
        tcrit = 0.005
        fc = 1000.0
        sizes = [100, 150, 80]
        # Generate passive single exponential bursts + dead time
        table = {
            "passive_wt": {f"p{i}": list(rng.exponential(0.025, n) + td) for i, n in enumerate(sizes)}
        }
        res = run_on_patches.run(table, dead_time=td, filter_cutoff=fc, burst_criterion=tcrit)
        decision = res["constructs"]["passive_wt"]["decision"]
        self.assertEqual(decision, "does_not_reject")

    def test_two_step_cascade_rejects_under_instrument_calibration(self):
        """True two-step non-passive bursts must still reject mixture under calibration."""
        rng = np.random.default_rng(202)
        td = 0.002
        sizes = [150, 200, 120]
        # Sequential two-step dwell: sum of two exponentials (CV^2 < 1)
        two_step = lambda n: rng.exponential(1 / 12, n) + rng.exponential(1 / 8, n) + td
        table = {
            "wt_twostep": {f"p{i}": list(two_step(n)) for i, n in enumerate(sizes)}
        }
        res = run_on_patches.run(table, dead_time=td, filter_cutoff=1000, burst_criterion=0.005)
        decision = res["constructs"]["wt_twostep"]["decision"]
        self.assertEqual(decision, "rejects_mixture")


if __name__ == "__main__":
    unittest.main()
