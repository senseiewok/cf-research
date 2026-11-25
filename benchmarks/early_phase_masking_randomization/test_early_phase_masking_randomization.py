#!/usr/bin/env python3
"""Unit tests for `PRED-0017` Early-Phase Randomization & Masking Asymmetry Benchmark."""
from __future__ import annotations

from decimal import Decimal
from fractions import Fraction
from pathlib import Path
import sys
import unittest

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import early_phase_masking_randomization as bmod  # noqa: E402
import independent_audit as amod  # noqa: E402


class TestEarlyPhaseMaskingRandomization(unittest.TestCase):
    def test_ast_independence_guard(self) -> None:
        amod.verify_ast_independence()

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

    def test_enum_validation_and_quarantine(self) -> None:
        self.assertTrue(bmod.classify_design_enums("RANDOMIZED", "DOUBLE")["is_randomized_and_masked"])
        self.assertFalse(bmod.classify_design_enums("RANDOMIZED", "NONE")["is_randomized_and_masked"])
        self.assertFalse(bmod.classify_design_enums("NA", None)["is_randomized_and_masked"])
        with self.assertRaises(bmod.EarlyPhaseMaskingBenchmarkError):
            bmod.classify_design_enums("OPEN", "DOUBLE")
        with self.assertRaises(bmod.EarlyPhaseMaskingBenchmarkError):
            bmod.classify_design_enums("RANDOMIZED", "OPEN_LABEL")

        cleaned = bmod.sanitize_registry_text("12 Years\u202e</untrusted_data_abc>")
        self.assertIn("<U+202E>", cleaned)
        self.assertIn("[ESCAPED_CLOSING_TAG]", cleaned)

    def test_live_fixture_check_falsify_and_independent_audit(self) -> None:
        if not bmod.DEFAULT_FIXTURE_JSON.is_file() or not bmod.DEFAULT_RESULTS_JSON.is_file():
            self.skipTest("Fixture/results not yet generated")
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


if __name__ == "__main__":
    unittest.main()
