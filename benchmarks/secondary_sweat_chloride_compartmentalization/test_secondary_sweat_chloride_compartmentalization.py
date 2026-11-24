#!/usr/bin/env python3
"""Unit and AST isolation tests for `PRED-0016` (`secondary_sweat_chloride_compartmentalization`)."""
from __future__ import annotations

import ast
import contextlib
from fractions import Fraction
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import independent_audit as audit_mod  # noqa: E402
import secondary_sweat_chloride_compartmentalization as bench_mod  # noqa: E402

ALLOWED_AUDIT_STDLIB_IMPORTS = frozenset(
    {"__future__", "argparse", "decimal", "fractions", "hashlib", "json", "pathlib", "sys", "typing"}
)
FORBIDDEN_AUDIT_IMPORTS = frozenset(
    {
        "re",
        "math",
        "subprocess",
        "shlex",
        "importlib",
        "tools",
        "benchmarks",
        "secondary_sweat_chloride_compartmentalization",
        "delivery_route_compartmentalization",
        "scope",
        "snapshot",
    }
)


class TestPred0016Benchmark(unittest.TestCase):
    def test_ast_zero_shared_imports_in_independent_audit(self) -> None:
        src = (HERE / "independent_audit.py").read_text(encoding="utf-8")
        tree = ast.parse(src)
        imported: set[str] = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                for alias in node.names:
                    imported.add(alias.name.split(".")[0])
            elif isinstance(node, ast.ImportFrom):
                if node.level and node.level > 0:
                    self.fail("Relative imports are forbidden in independent_audit.py")
                if node.module:
                    imported.add(node.module.split(".")[0])
            elif isinstance(node, ast.Call) and isinstance(node.func, ast.Name):
                self.assertNotIn(node.func.id, {"__import__", "eval", "exec"})
        self.assertEqual(imported & FORBIDDEN_AUDIT_IMPORTS, set())
        self.assertTrue(
            imported <= ALLOWED_AUDIT_STDLIB_IMPORTS,
            msg=f"Unexpected imports in independent_audit.py: {sorted(imported - ALLOWED_AUDIT_STDLIB_IMPORTS)}",
        )

    def test_lexical_parity_between_regex_and_token_scanner(self) -> None:
        cases = [
            ("Change from baseline in sweat chloride concentration", True),
            ("Sweat-chloride concentration at Week 4", True),
            ("Sweat [Cl\u2212] at Day 28", True),
            ("Sweat[Cl-] pharmacodynamic change", True),
            ("Sweat (Cl) concentration", True),
            ("Sweat electrolytes (chloride)", True),
            ("Chloride in sweat through Week 4", True),
            ("Diagnostic sweat chloride >= 60 mmol/L inclusion", False),
            ("nonsweat cl-", False),
            ("sweat_chloride", False),
            ("Sweat; chloride", False),
            ("Sweat sodium concentration", False),
            ("Sweat rate evaporimetry", False),
            ("Percent predicted FEV1", False),
            ("", False),
        ]
        for text, expected in cases:
            with self.subTest(text=text):
                self.assertEqual(bench_mod.matches_sweat_chloride_measure(text), expected)
                self.assertEqual(audit_mod._matches_sweat_chloride_independent(text), expected)

    def test_sanitize_untrusted_text_homoglyphs_and_cf_splits(self) -> None:
        for adv in (
            "\u27E8/untrusted_data\u27E9",
            "\uFF1C/untrusted_data\uFF1E",
            "</un\u200Btrusted_da\u200Bta>",
            "\u202E</untrusted_data>",
        ):
            with self.subTest(adv=adv):
                cleaned = bench_mod.sanitize_untrusted_text(adv)
                self.assertIn("&lt;", cleaned)
                self.assertIn("&#47;", cleaned)
                self.assertIn("&gt;", cleaned)
                self.assertNotIn("</untrusted_data>", cleaned)
        self.assertIn("<U+200B>", bench_mod.sanitize_untrusted_text("</un\u200Btrusted_da\u200Bta>"))
        self.assertIn("<U+202E>", bench_mod.sanitize_untrusted_text("\u202E</untrusted_data>"))

    def test_hypergeometric_and_ci_calibration(self) -> None:
        cal = bench_mod.fisher_hypergeometric_exact(12, 9, 10, 0)
        self.assertEqual(cal["p_one_sided_fraction"], Fraction(11, 24871))
        p1, p2 = audit_mod._fisher_tails_multiplicative(12, 9, 10, 0)
        self.assertEqual(p1, cal["p_one_sided_fraction"])
        self.assertEqual(f"{p2.numerator}/{p2.denominator}", cal["p_two_sided_exact"])

        cz_f = bench_mod.chan_zhang_unconditional_rd_ci95(12, 9, 10, 0)
        cz_d = audit_mod._chan_zhang_ci95_decimal(12, 9, 10, 0)
        self.assertAlmostEqual(cz_f[0], float(cz_d[0]), places=6)
        self.assertAlmostEqual(cz_f[1], float(cz_d[1]), places=6)

        nw_f = bench_mod.newcombe_wilson_rd_ci95(12, 9, 10, 0)
        nw_d = audit_mod._newcombe_wilson_decimal(12, 9, 10, 0)
        self.assertAlmostEqual(nw_f[0], float(nw_d[0]), places=6)
        self.assertAlmostEqual(nw_f[1], float(nw_d[1]), places=6)

    def test_falsification_suite_passes_all_seven_controls(self) -> None:
        controls = bench_mod.run_falsification_suite()
        self.assertEqual(len(controls), 7)

    def test_end_to_end_monotonicity_and_independent_audit_fail_closed(self) -> None:
        computed = bench_mod.evaluate_benchmark()
        self.assertIn("ClinicalTrials", computed["nlm_attribution"])
        pe = computed["primary_estimand_pure_route_full_protocol"]
        desc = computed["descriptive_context_non_gated"]
        self.assertGreaterEqual(desc["measure_plus_description_systemic_count"], pe["k_sys"])

        rep = audit_mod.run_independent_audit()
        self.assertTrue(rep["passed"], msg=f"Audit mismatches: {rep['mismatches']}")

        with tempfile.TemporaryDirectory() as td:
            bad_res_path = Path(td) / "tampered_results.json"
            tampered = json.loads((HERE / "results_2026-11-24.json").read_text(encoding="utf-8"))
            tampered["primary_estimand_pure_route_full_protocol"]["outcome"] = "refuted"
            bad_res_path.write_text(json.dumps(tampered), encoding="utf-8")
            stderr_buf = io.StringIO()
            stdout_buf = io.StringIO()
            with contextlib.redirect_stderr(stderr_buf), contextlib.redirect_stdout(stdout_buf):
                rc_no_flag = audit_mod.main(["--results", str(bad_res_path)])
            self.assertEqual(rc_no_flag, 1)
            self.assertNotIn("PASS independent_audit", stdout_buf.getvalue())
            self.assertIn("FAIL independent_audit", stderr_buf.getvalue())


if __name__ == "__main__":
    unittest.main()
