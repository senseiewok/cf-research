"""Tests for scope.py on the SYNTHETIC snapshot (synthetic_fixtures.py), including the planted non-CF ids. No network.
Usage: python -m unittest -v (inside tools/trial_atlas)."""
import contextlib
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import scope  # noqa: E402
import snapshot as snap  # noqa: E402
import synthetic_fixtures as sf  # noqa: E402


class ScopeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._tmp = tempfile.TemporaryDirectory()
        cls.path = Path(cls._tmp.name) / "snap"
        sf.build_cf_snapshot(cls.path)
        cls.s = snap.load(cls.path)

    @classmethod
    def tearDownClass(cls):
        cls._tmp.cleanup()

    def setUp(self):
        sf.block_network(self)
        self.result = scope.apply_scope(self.s.records("condition"))
        self.excluded = {e["nct_id"]: e["rule"] for e in self.result["excluded"]}

    def test_planted_non_cf_ids_are_excluded(self):
        for nct in sf.PLANTED_NON_CF_IDS:
            self.assertIn(nct, self.excluded)
        self.assertEqual(self.excluded["NCT00000011"], "X5")       # non-CF bronchiectasis: CF named only in a negation
        self.assertEqual(self.excluded["NCT00000016"], "X4")       # asthma: CF not named

    def test_each_rule_has_its_count_and_reason(self):
        c = self.result["counts"]
        self.assertEqual(c["studies_in_route"], 16)
        self.assertEqual(c["included"], 12)
        self.assertEqual(c["excluded_by_rule"], {"X1": 1, "X2": 1, "X3": 0, "X4": 1, "X5": 1, "X6": 0})
        self.assertEqual(self.excluded["NCT00000009"], "X1")
        self.assertEqual(self.excluded["NCT00000010"], "X2")
        for e in self.result["excluded"]:
            self.assertTrue(e["reason"])

    def test_cf_only_and_cf_among_others(self):
        flags = {r["nct_id"]: r["cf_flag"] for r in self.result["included"]}
        self.assertEqual(flags["NCT00000003"], "CF among others")
        self.assertEqual(flags["NCT00000013"], "CF only")         # "Cystic Fibrosis-related Diabetes" names CF
        self.assertEqual(self.result["counts"]["cf_only"], 11)
        self.assertEqual(self.result["counts"]["cf_among_others"], 1)

    def test_planned_versus_actual_start(self):
        info = {r["nct_id"]: r for r in self.result["included"]}
        self.assertEqual(info["NCT00000004"]["start_kind"], "planned")
        self.assertEqual(info["NCT00000001"]["start_kind"], "actual")
        self.assertEqual(info["NCT00000014"]["start_kind"], "unknown")   # a date with no type is not called actual
        self.assertEqual(info["NCT00000014"]["start_year"], 2000)
        self.assertEqual(info["NCT00000012"]["start_decade"], "2000s")

    def test_condition_matching(self):
        self.assertTrue(scope.condition_names_cf("Cystic Fibrosis"))
        self.assertTrue(scope.condition_names_cf("cystic-fibrosis"))
        self.assertTrue(scope.condition_names_cf("Mucoviscidosis"))
        self.assertTrue(scope.condition_names_cf("CF"))
        self.assertIsNone(scope.condition_names_cf("Non-Cystic Fibrosis Bronchiectasis"))
        self.assertIsNone(scope.condition_names_cf("Bronchiectasis without cystic fibrosis"))
        self.assertIsNone(scope.condition_names_cf("non-CF bronchiectasis"))
        self.assertFalse(scope.condition_names_cf("Pulmonary Fibrosis"))
        self.assertFalse(scope.condition_names_cf("CF-related diabetes registry"))   # bare CF counts only as the whole condition

    def test_keywords_do_not_bring_a_study_into_scope(self):
        term = scope.apply_scope(self.s.records("term"))
        self.assertIn("NCT00000018", {e["nct_id"] for e in term["excluded"]})

    def test_manual_exclusions(self):
        r = scope.apply_scope(self.s.records("condition"), {"NCT00000001": "duplicate registration (SYNTHETIC reason)"})
        self.assertEqual(r["counts"]["excluded_by_rule"]["X6"], 1)
        self.assertIn("duplicate registration", next(e["reason"] for e in r["excluded"] if e["nct_id"] == "NCT00000001"))

    def test_route_difference(self):
        rd = scope.route_difference(self.s)
        self.assertEqual(rd["only_main"], ["NCT00000014", "NCT00000015"])
        self.assertEqual(rd["only_recall"], ["NCT00000017"])
        self.assertEqual(rd["share_unexplained"], 3 / 12)
        rd2 = scope.route_difference(self.s, explained={"NCT00000014": "x", "NCT00000015": "y"})
        self.assertEqual(rd2["unexplained"], ["NCT00000017"])

    def test_cli_prints_every_count(self):
        with contextlib.redirect_stdout(io.StringIO()) as out:
            self.assertEqual(scope.main([str(self.path)]), 0)
        text = out.getvalue()
        for rid in scope.RULES:
            self.assertIn(f"excluded {rid}:", text)
        self.assertIn("NCT00000011  X5", text)

    def test_read_id_reasons_refuses_bad_files(self):
        p = Path(self._tmp.name) / "bad.json"
        p.write_text(json.dumps({"NCT1": "too short an id"}), encoding="utf-8")
        with self.assertRaises(ValueError):
            scope.read_id_reasons(p)
        p.write_text(json.dumps({"NCT00000001": ""}), encoding="utf-8")
        with self.assertRaises(ValueError):
            scope.read_id_reasons(p)


if __name__ == "__main__":
    unittest.main()
