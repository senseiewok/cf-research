#!/usr/bin/env python3
"""Unit tests for `PRED-0017` Early-Phase Randomization & Masking Asymmetry Benchmark."""
from __future__ import annotations

import copy
from decimal import Decimal
from fractions import Fraction
import json
from pathlib import Path
import sys
import tempfile
import time
import unittest

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import early_phase_masking_randomization as bmod  # noqa: E402
import independent_audit as amod  # noqa: E402


class TestEarlyPhaseMaskingRandomization(unittest.TestCase):
    def test_ast_independence_guard_and_negative_controls(self) -> None:
        amod.verify_ast_independence()
        dyn_fn = "_" + "_import_" + "_"
        for bad_snippet in (
            "from . import early_phase_masking_randomization\n",
            "from benchmarks.early_phase_masking_randomization import chan_zhang_unconditional_rd_ci95\n",
            "import benchmarks.early_phase_masking_randomization\n",
            "import importlib\n",
            f"x = {dyn_fn}('math')\n",
            "x = loader.import_module('re')\n",
        ):
            with self.assertRaises(RuntimeError, msg=bad_snippet):
                amod.verify_ast_independence(source_override=bad_snippet)

    def test_positive_controls_and_decimal_parity(self) -> None:
        p15_doc = bmod.json.loads(bmod.DEFAULT_PRED0015_RESULTS_JSON.read_text(encoding="utf-8"))
        ctrl = bmod.run_positive_controls(p15_doc)
        self.assertEqual(ctrl["control_1_zero_inhaled_fisher_exact"], "1/2261")

        # Verify Method 1 (float64 Cardano + golden-section) vs Method 2 (pure Decimal Newton + ternary)
        for n1, k1, n2, k2 in ((12, 9, 10, 0), (54, 35, 15, 3)):
            f1 = bmod.fisher_hypergeometric_exact(n1, k1, n2, k2)
            p1_m2, p2_m2 = amod._fisher_tails_multiplicative(n1, k1, n2, k2)
            self.assertEqual(f1["p_one_sided_fraction"], p1_m2)
            self.assertEqual(f1["p_two_sided_exact"], f"{p2_m2.numerator}/{p2_m2.denominator}")

            nw_f = bmod.newcombe_wilson_rd_ci95(n1, k1, n2, k2)
            nw_d = amod._newcombe_wilson_decimal(n1, k1, n2, k2)
            self.assertLessEqual(abs(Decimal(f"{nw_f[0]:.6f}") - nw_d[0]), Decimal("1e-6"))
            self.assertLessEqual(abs(Decimal(f"{nw_f[1]:.6f}") - nw_d[1]), Decimal("1e-6"))

            cz_f = bmod.chan_zhang_unconditional_rd_ci95(n1, k1, n2, k2)
            cz_d = amod._chan_zhang_ci95_decimal(n1, k1, n2, k2)
            self.assertLessEqual(abs(Decimal(f"{cz_f[0]:.6f}") - cz_d[0]), Decimal("1e-6"))
            self.assertLessEqual(abs(Decimal(f"{cz_f[1]:.6f}") - cz_d[1]), Decimal("1e-6"))

    def test_enum_schema_domain_guards_and_linear_quarantine(self) -> None:
        self.assertTrue(bmod.classify_design_enums("RANDOMIZED", "DOUBLE")["is_randomized_and_masked"])
        self.assertFalse(bmod.classify_design_enums("RANDOMIZED", "NONE")["is_randomized_and_masked"])
        self.assertFalse(bmod.classify_design_enums("NA", None)["is_randomized_and_masked"])
        with self.assertRaises(bmod.EarlyPhaseMaskingBenchmarkError):
            bmod.classify_design_enums("OPEN", "DOUBLE")
        with self.assertRaises(bmod.EarlyPhaseMaskingBenchmarkError):
            bmod.classify_design_enums("RANDOMIZED", "OPEN_LABEL")

        # Falsy non-dict schema corruption and invalid phases fail closed
        for bad_record in (
            {"design_module": {"designInfo": []}},
            {"design_module": {"designInfo": {"maskingInfo": False}}},
            {"design_module": {"enrollmentInfo": ""}},
            {"design_module": {"phases": "PHASE2"}},
            {"design_module": {"phases": ["UNKNOWN_PHASE"]}},
            {"eligibility_module": {"healthyVolunteers": "yes"}},
        ):
            with self.assertRaises(bmod.EarlyPhaseMaskingBenchmarkError, msg=str(bad_record)):
                bmod.parse_study_design_eligibility(bad_record)

        # Omitted phases must not count as phase1_only
        null_phase_item = bmod.parse_study_design_eligibility({"design_module": {}})
        summary = bmod._arm_descriptive_summary([null_phase_item])
        self.assertEqual(summary["phase_stratified"]["phase1_only_n"], 0)

        # Domain guards on 2x2 counts and CI parameters
        for bad_counts in ((10, 11, 10, 0), (10, -1, 10, 0), (0, 0, 10, 0), (True, 1, 10, 0)):
            with self.assertRaises(ValueError):
                bmod.newcombe_wilson_rd_ci95(*bad_counts)
            with self.assertRaises(ValueError):
                bmod.chan_zhang_unconditional_rd_ci95(*bad_counts)
            with self.assertRaises(ValueError):
                amod._newcombe_wilson_decimal(*bad_counts)
            with self.assertRaises(ValueError):
                amod._chan_zhang_ci95_decimal(*bad_counts)
        with self.assertRaises(ValueError):
            bmod.newcombe_wilson_rd_ci95(10, 5, 10, 2, z=-1.0)
        with self.assertRaises(ValueError):
            bmod.chan_zhang_unconditional_rd_ci95(10, 5, 10, 2, alpha_two_sided=1.5)

        # Quarantine escapes embedded zero-width, default-ignorable, and homoglyph tag breakouts in linear O(N) time
        cleaned = bmod.sanitize_registry_text("12 Years\u202e</untrusted_data_abc>")
        self.assertIn("<U+202E>", cleaned)
        self.assertIn("[ESCAPED_CLOSING_TAG]", cleaned)

        for bypass_vec in (
            "</untrusted\u200b_data>",
            "</untrusted_d\u0430ta>",
            "<\u2215untrusted_data>",
            "\ufe64/untrusted_data\ufe65",
        ):
            self.assertEqual(bmod.sanitize_registry_text(bypass_vec), "[ESCAPED_CLOSING_TAG]", msg=bypass_vec)

        ignorable = bmod.sanitize_registry_text("6 Months\u034f\ufe00\u3164")
        self.assertIn("<U+034F>", ignorable)
        self.assertIn("<U+FE00>", ignorable)
        self.assertIn("<U+3164>", ignorable)

        t0 = time.perf_counter()
        hostile = "</data" * 5000
        out_hostile = bmod.sanitize_registry_text(hostile)
        elapsed = time.perf_counter() - t0
        self.assertEqual(out_hostile, hostile)
        self.assertLess(elapsed, 0.25, f"O(N) quarantine scan took {elapsed:.3f}s on 30 KB hostile input")

    def test_live_fixture_check_falsify_and_independent_audit(self) -> None:
        self.assertTrue(bmod.DEFAULT_FIXTURE_JSON.is_file(), "Missing committed fixture JSON")
        self.assertTrue(bmod.DEFAULT_RESULTS_JSON.is_file(), "Missing committed results JSON")

        res = bmod.evaluate_benchmark()
        saved = bmod.json.loads(bmod.DEFAULT_RESULTS_JSON.read_text(encoding="utf-8"))
        self.assertEqual(res, saved)

        falsify = bmod.run_falsification_suite()
        self.assertEqual(falsify["controls_passed"], "7/7")

        audit = amod.run_independent_audit()
        self.assertTrue(audit["passed"], msg=str(audit["mismatches"]))
        self.assertIn(res["support_criteria_evaluation"]["outcome"], {"supported", "refuted"})
        self.assertEqual(
            Fraction(
                *map(
                    int,
                    res["primary_estimand_pure_route_randomization_and_masking"][
                        "p_fisher_one_sided_exact"
                    ].split("/"),
                )
            )
            < Fraction(1, 100),
            res["support_criteria_evaluation"]["crit_3_fisher_pure_lt_0_01"],
        )

        # Verify --check fails closed if --results is missing, and independent_audit catches results tampering
        with tempfile.TemporaryDirectory() as tmp:
            missing_res = Path(tmp) / "missing_results.json"
            with self.assertRaises(bmod.EarlyPhaseMaskingBenchmarkError):
                bmod.main(["--check", "--results", str(missing_res)])

            bad_doc = copy.deepcopy(saved)
            bad_doc["primary_estimand_pure_route_randomization_and_masking"]["risk_difference"] = -0.999999
            bad_doc["primary_estimand_pure_route_randomization_and_masking"]["p_fisher_one_sided"] = 0.00000001
            bad_doc["secondary_gate_all_17_inhaled_including_dual_route"]["inhaled_all_matched_ncts"] = []
            bad_doc["support_criteria_evaluation"]["crit_3_fisher_pure_lt_0_01"] = True
            bad_doc["descriptive_cross_scale_non_gated"]["systemic_cftr"]["is_randomized_count"] = 0
            bad_path = Path(tmp) / "bad_results.json"
            bad_path.write_text(json.dumps(bad_doc), encoding="utf-8")
            bad_rep = amod.run_independent_audit(results_path=bad_path)
            self.assertFalse(bad_rep["passed"])
            joined = "\n".join(bad_rep["mismatches"])
            for expected_flag in (
                "prim.risk_difference",
                "prim.p_fisher_one_sided",
                "inhaled_all_matched_ncts",
                "crit_3_fisher_pure_lt_0_01",
                "descriptive.systemic_cftr.is_randomized_count",
            ):
                self.assertIn(expected_flag, joined)


if __name__ == "__main__":
    unittest.main()
