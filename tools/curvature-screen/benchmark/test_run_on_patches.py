#!/usr/bin/env python3
"""Unit and planted-defect tests for tools/curvature-screen/benchmark/whole_sample/run_on_patches.py.

Covers:
- T-0923: 0.36/f_c default resolution, instrument-calibrated threshold (dead time, f_c, t_crit),
  passive mixture non-rejection, two-step cascade rejection, invalid parameter validation.
- Bugfixes 1 & 2: quantized Null.draws unit consistency, dead-time truncation (no 0.0 clamping),
  scale-distinct cache keys, and <= 2% false-alarm rate under quantized + filtered calibration.
- Plan Step 1: CSV/JSON metadata parsing, n_channels != 1 and unstable exclusions before analysis,
  one_active_channel_verified flag and explicit note when channel counts are omitted.
- Plan Step 6: composite-null GLRT comparator (continuous & quantized), analytical gradient accuracy,
  Erlang-2 smooth limit, ll_mix >= ll_exp - 1e-6 monotonicity, cross-patch sum, resolution sensitivity.
- Plan Step 7: control power simulation (low power for small control, high power for large control)
  and Section 5 failure flags on both rejecting control and high-power comparator divergence.
"""
import json
import math
import os
import tempfile
import sys
import unittest
import numpy as np
from scipy import optimize

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
        # Uncalibrated branch preserves original 5-element RNG stream
        rng_orig = np.random.default_rng([20261113, 100, 0, 0, 0, 0])
        np.testing.assert_array_equal(d1, run_on_patches.cv2(rng_orig.exponential(1.0, (run_on_patches.REPS, 100))))

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


class TestNullQuantizationAndTruncation(unittest.TestCase):
    """Verify Bugfixes 1 & 2: quantized calibrated units, truncation (no 0 clamping), scale cache keys."""

    def test_quantized_calibrated_null_units_and_scale_distinct_cache_keys(self):
        null = run_on_patches.Null(seed=20261113)
        d_small = null.draws(
            100,
            q_over_mean=0.2,
            dead_time=0.002,
            filter_cutoff=1000.0,
            burst_criterion=0.005,
            mean_duration=0.005,
        )
        d_large = null.draws(
            100,
            q_over_mean=0.2,
            dead_time=0.002,
            filter_cutoff=1000.0,
            burst_criterion=0.005,
            mean_duration=5.0,
        )
        self.assertEqual(len(null.cache), 2)
        crit_small = float(np.quantile(d_small, run_on_patches.ALPHA))
        crit_large = float(np.quantile(d_large, run_on_patches.ALPHA))
        self.assertGreater(crit_small, 0.55)
        self.assertLess(crit_small, 0.80)
        self.assertGreater(crit_large, 0.55)
        self.assertLess(crit_large, 0.80)

    def test_dead_time_truncation_avoids_zero_clamping_inflation(self):
        # With eff_dt = 0.005 and mean_duration = 0.01 (tau_o = 0.01), ~39% of raw draws are < eff_dt.
        # Clamping short bursts to 0.0 would inflate mean(CV^2) > 1.25; true truncation keeps mean(CV^2) ~ 1.0.
        null = run_on_patches.Null(seed=20261113)
        d_trunc = null.draws(100, dead_time=0.005, filter_cutoff=72.0, burst_criterion=0.0, mean_duration=0.01)
        self.assertAlmostEqual(float(d_trunc.mean()), 1.0, delta=0.05)

    def test_quantized_instrument_calibrated_false_alarm_rate_under_2pct(self):
        rng = np.random.default_rng(303)
        null = run_on_patches.Null(seed=20261113)
        q = 0.005
        td = 0.002
        fc = 1000.0
        tc = 0.005
        n_target = 100
        trials = 250
        rej = 0
        for _ in range(trials):
            raw = rng.exponential(0.05, 250)
            raw = raw[raw >= td][:120]
            y = np.ceil(raw / q) * q
            adj = run_on_patches.adjust(y, td, q)[:n_target]
            s = run_on_patches.patch_stats(adj, q, null, dead_time=td, filter_cutoff=fc, burst_criterion=tc)
            if s["cv2"] < s["crit"]:
                rej += 1
        self.assertLessEqual(rej / trials, 0.02)


class TestPlanStep1MetadataAndExclusions(unittest.TestCase):
    """Verify Plan Step 1: one-active-channel check, unstable exclusion, CSV & JSON metadata."""

    def test_csv_exclusions_and_channel_verification_flag(self):
        rng = np.random.default_rng(404)
        with tempfile.TemporaryDirectory() as tmp:
            csv_path = os.path.join(tmp, "bursts.csv")
            with open(csv_path, "w", encoding="utf-8", newline="") as f:
                f.write("patch_id,duration,construct,n_channels,unstable,exclude_reason\n")
                for v in rng.exponential(0.1, 80) + 0.002:
                    f.write(f"wt_ok,{v:.6f},WT,1,0,\n")
                for v in rng.exponential(0.1, 80) + 0.002:
                    f.write(f"wt_two_ch,{v:.6f},WT,2,false,\n")
                for v in rng.exponential(0.1, 80) + 0.002:
                    f.write(f"wt_unstable,{v:.6f},WT,1,yes,baseline drift\n")
                for v in rng.exponential(0.5, 20) + 0.002:
                    f.write(f"ctrl_unver,{v:.6f},E1371Q,,no,\n")

            # Backwards compatibility: default read_csv returns 2-tuple
            table2, dropped2 = run_on_patches.read_csv(csv_path)
            self.assertEqual(dropped2, 0)
            self.assertIn("WT", table2)

            table, dropped, meta = run_on_patches.read_csv(csv_path, return_metadata=True)
            self.assertEqual(dropped, 0)
            res = run_on_patches.run(table, dead_time=0.002, control="E1371Q", metadata=meta, bootstrap_reps=200)

            wt = res["constructs"]["WT"]
            self.assertTrue(wt["one_active_channel_verified"])
            self.assertEqual(wt["n_patches_included"], 1)
            wt_by_id = {p["patch_id"]: p for p in wt["patches"]}
            self.assertTrue(wt_by_id["wt_ok"]["included"])
            self.assertFalse(wt_by_id["wt_two_ch"]["included"])
            self.assertTrue(wt_by_id["wt_two_ch"]["excluded_before_analysis"])
            self.assertEqual(wt_by_id["wt_two_ch"]["reason"], "n_channels != 1 (N=2)")
            self.assertFalse(wt_by_id["wt_unstable"]["included"])
            self.assertEqual(wt_by_id["wt_unstable"]["reason"], "baseline drift")

            ctrl = res["constructs"]["E1371Q"]
            self.assertFalse(ctrl["one_active_channel_verified"])
            self.assertEqual(
                ctrl["one_active_channel_note"],
                "active-channel count not supplied in input; assumed 1 per caller",
            )
            self.assertEqual(ctrl["n_patches_included"], 1)

    def test_cli_metadata_json_merge_and_active_channels_alias(self):
        rng = np.random.default_rng(405)
        with tempfile.TemporaryDirectory() as tmp:
            csv_path = os.path.join(tmp, "bursts.csv")
            meta_path = os.path.join(tmp, "meta.json")
            out_path = os.path.join(tmp, "out.json")
            with open(csv_path, "w", encoding="utf-8", newline="") as f:
                f.write("patch_id,duration,construct,active_channels\n")
                for v in rng.exponential(1 / 12, 90) + rng.exponential(1 / 8, 90) + 0.002:
                    f.write(f"p1,{v:.6f},WT,1\n")
                for v in rng.exponential(1 / 12, 90) + rng.exponential(1 / 8, 90) + 0.002:
                    f.write(f"p2,{v:.6f},WT,1\n")
            with open(meta_path, "w", encoding="utf-8") as mf:
                json.dump(
                    {
                        "p2": {
                            "n_channels": 3,
                            "unstable": True,
                            "reason": "three active channels",
                            "dead_time": 0.002,
                        }
                    },
                    mf,
                )
            rc = run_on_patches.main(
                [csv_path, "--dead-time", "0.002", "--metadata", meta_path, "--bootstrap-reps", "100", "--out", out_path]
            )
            self.assertEqual(rc, 0)
            with open(out_path, encoding="utf-8") as of:
                res = json.load(of)
            self.assertEqual(res["not_implemented"], [])
            wt = res["constructs"]["WT"]
            self.assertTrue(wt["one_active_channel_verified"])
            p_map = {p["patch_id"]: p for p in wt["patches"]}
            self.assertFalse(p_map["p2"]["included"])
            self.assertEqual(p_map["p2"]["reason"], "n_channels != 1 (N=3)")

    def test_reason_dedup_substring_safety_and_non_dict_metadata_raises(self):
        rng = np.random.default_rng(406)
        with tempfile.TemporaryDirectory() as tmp:
            csv_path = os.path.join(tmp, "bursts.csv")
            bad_meta_path = os.path.join(tmp, "bad_meta.json")
            with open(csv_path, "w", encoding="utf-8", newline="") as f:
                f.write("patch_id,duration,construct,n_channels,unstable,exclude_reason\n")
                vals = rng.exponential(0.1, 60) + 0.002
                f.write(f"p_ok,{vals[0]:.6f},WT,1,0,\n")
                for v in vals[1:]:
                    f.write(f"p_ok,{v:.6f},WT,1,0,\n")
                f.write("p_bad,0.1,WT,1,1,baseline drift\n")
                f.write("p_bad,0.2,WT,1,1,drift\n")
                f.write("p_bad,0.3,WT,1,1,drift\n")
            _, _, meta = run_on_patches.read_csv(csv_path, return_metadata=True)
            self.assertEqual(meta["WT"]["p_bad"]["reason"], "baseline drift; drift")

            with open(bad_meta_path, "w", encoding="utf-8") as mf:
                json.dump(["not", "a", "dict"], mf)
            with self.assertRaises(ValueError):
                run_on_patches.main([csv_path, "--dead-time", "0.002", "--metadata", bad_meta_path])


class TestCompositeNullGLRTComparator(unittest.TestCase):
    """Verify Plan Step 6 GLRT comparator (continuous & quantized, Erlang-2 limit, monotonicity)."""

    def test_analytical_gradients_and_erlang2_limit(self):
        rng = np.random.default_rng(501)
        x = rng.exponential(0.5, 150)
        prep_c = run_on_patches._prepare_glrt_data(x, 0.0)
        err_mix = optimize.check_grad(
            lambda p: run_on_patches._eval_mix_nll_and_grad(p, prep_c)[0],
            lambda p: run_on_patches._eval_mix_nll_and_grad(p, prep_c)[1],
            np.array([0.4, -0.5, 0.7]),
        )
        err_hypo = optimize.check_grad(
            lambda p: run_on_patches._eval_hypo_nll_continuous(p, prep_c)[0],
            lambda p: run_on_patches._eval_hypo_nll_continuous(p, prep_c)[1],
            np.array([0.5, -14.0]),
        )
        self.assertLess(err_mix, 1e-4)
        self.assertLess(err_hypo, 1e-4)

    def test_continuous_two_step_rejects_and_mixture_does_not_reject(self):
        rng = np.random.default_rng(502)
        x_two = rng.exponential(1.0 / 12.0, 200) + rng.exponential(1.0 / 8.0, 200)
        c_two = run_on_patches.comparator_patch_stats(x_two, quantum=0.0, bootstrap_reps=250)
        self.assertGreaterEqual(c_two["ll_mix"], c_two["ll_exp"] - 1e-6)
        self.assertTrue(c_two["rejects_mixture"])
        self.assertLess(c_two["p_one_sided"], 0.01)

        x_mix = np.where(rng.random(200) < 0.5, rng.exponential(1.0 / 12.0, 200), rng.exponential(1.0 / 2.0, 200))
        c_mix = run_on_patches.comparator_patch_stats(x_mix, quantum=0.0, bootstrap_reps=250)
        self.assertGreaterEqual(c_mix["ll_mix"], c_mix["ll_exp"] - 1e-6)
        self.assertFalse(c_mix["rejects_mixture"])
        self.assertLessEqual(c_mix["lambda_comp"], c_mix["lambda_exp"] + 1e-8)

    def test_quantized_two_step_rejects_and_mixture_does_not_reject(self):
        rng = np.random.default_rng(503)
        q1 = 0.025
        raw_two = rng.exponential(1.0 / 12.0, 200) + rng.exponential(1.0 / 8.0, 200)
        x_two_q = np.floor(raw_two / q1) * q1 + q1 / 2.0
        c_two = run_on_patches.comparator_patch_stats(x_two_q, quantum=q1, bootstrap_reps=250)
        self.assertGreaterEqual(c_two["ll_mix"], c_two["ll_exp"] - 1e-6)
        self.assertTrue(c_two["rejects_mixture"])

        q2 = 0.05
        raw_mix = np.where(rng.random(200) < 0.5, rng.exponential(1.0 / 12.0, 200), rng.exponential(1.0 / 2.0, 200))
        x_mix_q = np.floor(raw_mix / q2) * q2 + q2 / 2.0
        c_mix = run_on_patches.comparator_patch_stats(x_mix_q, quantum=q2, bootstrap_reps=250)
        self.assertGreaterEqual(c_mix["ll_mix"], c_mix["ll_exp"] - 1e-6)
        self.assertFalse(c_mix["rejects_mixture"])

    def test_cross_patch_comparator_and_resolution_sensitivity(self):
        rng = np.random.default_rng(504)
        td = 0.002
        table = {
            "two_step": {f"p{i}": list(rng.exponential(1 / 12, n) + rng.exponential(1 / 8, n) + td) for i, n in enumerate((100, 120))},
            "mixture": {f"m{i}": list(np.where(rng.random(n) < 0.5, rng.exponential(1 / 12, n), rng.exponential(1 / 2, n)) + td) for i, n in enumerate((100, 120))},
        }
        res = run_on_patches.run(table, dead_time=td, bootstrap_reps=200)
        c_two = res["constructs"]["two_step"]["resolution_sensitivity"]["comparator_at_1x"]["combined"]
        c_mix = res["constructs"]["mixture"]["resolution_sensitivity"]["comparator_at_1x"]["combined"]
        self.assertTrue(c_two["rejects_mixture"])
        self.assertFalse(c_mix["rejects_mixture"])


class TestControlPowerAndSection5Failure(unittest.TestCase):
    """Verify Plan Step 7 control power simulation and Section 5 failure criteria."""

    def test_control_power_low_for_small_and_high_for_large_control(self):
        rng = np.random.default_rng(601)
        td = 0.002
        table_small = {
            "wt": {"p0": list(rng.exponential(1 / 12, 120) + rng.exponential(1 / 8, 120) + td)},
            "ctrl": {f"c{i}": list(rng.exponential(0.25, n) + td) for i, n in enumerate((6, 9, 14))},
        }
        res_small = run_on_patches.run(table_small, dead_time=td, control="ctrl", bootstrap_reps=150)
        pwr_small = res_small["constructs"]["ctrl"]["control_power_simulation"]
        self.assertEqual(pwr_small["observed_patch_burst_counts"], [6, 9, 14])
        self.assertLess(pwr_small["power_at_variant_A_scale_matched"]["power_alpha_0.01"], 0.50)
        self.assertLess(pwr_small["power_at_variant_A_nominal_rates"]["power_alpha_0.01"], 0.50)
        self.assertFalse(res_small["section_5_failure_triggered"])

        table_large = {
            "wt": {"p0": list(rng.exponential(1 / 12, 120) + rng.exponential(1 / 8, 120) + td)},
            "ctrl": {f"c{i}": list(rng.exponential(0.25, n) + td) for i, n in enumerate((150, 200))},
        }
        res_large = run_on_patches.run(table_large, dead_time=td, control="ctrl", bootstrap_reps=150)
        pwr_large = res_large["constructs"]["ctrl"]["control_power_simulation"]
        self.assertEqual(pwr_large["observed_patch_burst_counts"], [150, 200])
        self.assertGreater(pwr_large["power_at_variant_A_scale_matched"]["power_alpha_0.01"], 0.80)
        self.assertGreater(pwr_large["power_at_variant_A_nominal_rates"]["power_alpha_0.01"], 0.80)

    def test_section_5_failure_triggered_when_control_rejects(self):
        rng = np.random.default_rng(602)
        td = 0.002
        # Planted defect: sequential two-step bursts placed in the control construct
        table = {
            "wt": {"p0": list(rng.exponential(1 / 12, 120) + rng.exponential(1 / 8, 120) + td)},
            "ctrl": {f"c{i}": list(rng.exponential(1 / 12, 150) + rng.exponential(1 / 8, 150) + td) for i in range(2)},
        }
        res = run_on_patches.run(table, dead_time=td, control="ctrl", bootstrap_reps=150)
        self.assertEqual(res["constructs"]["ctrl"]["decision"], "rejects_mixture")
        self.assertTrue(res["constructs"]["ctrl"]["section_5_failure_triggered"])
        self.assertTrue(res["section_5_failure_triggered"])

    def test_section_5_failure_triggered_on_test_comparator_vs_cv2_divergence(self):
        # Planted case: seed=297 with N=250 two-step (12, 4) has CV^2 just above crit (does_not_reject)
        # while the GLRT comparator rejects at alpha=0.01 and CV^2 power at fitted rates is >= 0.80.
        rng = np.random.default_rng(297)
        x = rng.exponential(1 / 12, 250) + rng.exponential(1 / 4, 250) + 0.001
        res = run_on_patches.run({"wt": {"p0": list(x)}}, dead_time=0.001, bootstrap_reps=200)
        wt = res["constructs"]["wt"]
        self.assertEqual(wt["decision"], "does_not_reject")
        self.assertTrue(wt["comparator"]["combined"]["rejects_mixture"])
        self.assertGreaterEqual(wt["test_power_at_comparator_rates"]["power_alpha_0.01"], 0.80)
        self.assertTrue(wt["section_5_failure_triggered"])
        self.assertTrue(res["section_5_failure_triggered"])


if __name__ == "__main__":
    unittest.main()
