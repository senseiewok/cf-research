"""Tests for make_labelling_sheet.py on the SYNTHETIC snapshot: the sheet leaks no tag, rule, model result, NCT id or stratum; the
draw is reproducible for a seed; the oversamples are taken; an extension skips earlier rows; the CLI never overwrites. No network.
Usage: python -m unittest -v (inside tools/trial_atlas)."""
import contextlib
import csv
import io
import json
import re
import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import check_atlas as ca  # noqa: E402
import lexicon  # noqa: E402
import make_labelling_sheet as mls  # noqa: E402
import snapshot as snap  # noqa: E402
import synthetic_fixtures as sf  # noqa: E402

LEX = lexicon.load()


class SheetTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._tmp = tempfile.TemporaryDirectory()
        cls.root = Path(cls._tmp.name)
        sf.build_cf_snapshot(cls.root / "snap")
        cls.s = snap.load(cls.root / "snap")
        cls.tags = lexicon.tag_snapshot(cls.s, LEX)
        cls.tags_path = cls.root / "tags.json"
        cls.tags_path.write_text(json.dumps(cls.tags), encoding="utf-8")

    @classmethod
    def tearDownClass(cls):
        cls._tmp.cleanup()

    def setUp(self):
        sf.block_network(self)

    def test_the_sheet_leaks_no_tag(self):
        sheet, key, _ = mls.build(self.s, self.tags, LEX, n=20, seed=1)
        rows = list(csv.DictReader(io.StringIO(sheet)))
        self.assertEqual(list(rows[0].keys()), mls.COLUMNS)
        for r in rows:
            self.assertEqual((r["label_classes"], r["label_composite"], r["label_notes"]), ("", "", ""))
        leaks = [r["id"] for r in LEX.rules + LEX.composite] + [r["id"] for r in LEX.data["not_stated_rules"]]
        leaks += ["unclassified", "not_stated", "safety_tolerability", "by_rule", "model", "stratum", "rare class", "NCT0"]
        for word in leaks:
            self.assertNotIn(word, sheet, f"the sheet must not contain {word!r}")
        self.assertFalse(re.search(r"\d{4}s \|", sheet), "no stratum label in the sheet")
        self.assertTrue(all(r["entry_id"].startswith("NCT") for r in key["rows"]))    # the ids live only in the sealed key

    def test_the_draw_is_reproducible(self):
        a = mls.build(self.s, self.tags, LEX, n=10, seed=42)
        b = mls.build(self.s, self.tags, LEX, n=10, seed=42)
        self.assertEqual(a[0], b[0])
        self.assertEqual(a[1], b[1])
        c = mls.build(self.s, self.tags, LEX, n=10, seed=43)
        self.assertNotEqual([r["entry_id"] for r in a[1]["rows"]], [r["entry_id"] for r in c[1]["rows"]])

    def test_oversamples_are_taken_first(self):
        _, key, summary = mls.build(self.s, self.tags, LEX, n=6, seed=3)
        self.assertIn("NCT00000014:P3", [r["entry_id"] for r in key["rows"]])     # the one entry the lexicon cannot classify
        self.assertEqual(summary["unclassified"], 1)
        self.assertEqual(summary["drawn"], 6)
        self.assertEqual(len({r["entry_id"] for r in key["rows"]}), 6)

    def test_oversamples_take_at_most_half(self):
        _, _, summary = mls.build(self.s, self.tags, LEX, n=10, seed=3, rare_max=100)    # every class counts as rare here
        self.assertLessEqual(summary["unclassified"] + summary["rare_class"], 5)
        self.assertEqual(summary["stratified"], 10 - summary["unclassified"] - summary["rare_class"])

    def test_the_whole_pool_can_be_drawn(self):
        _, key, summary = mls.build(self.s, self.tags, LEX, n=50, seed=3)
        self.assertEqual(summary["drawn"], 19)            # 20 entries in scope, less the one with no primary outcome
        self.assertEqual(len(key["rows"]), 19)

    def test_an_extension_skips_earlier_rows(self):
        _, first, _ = mls.build(self.s, self.tags, LEX, n=8, seed=5)
        _, second, _ = mls.build(self.s, self.tags, LEX, n=8, seed=6, already=first)
        self.assertFalse({r["entry_id"] for r in first["rows"]} & {r["entry_id"] for r in second["rows"]})
        self.assertEqual(second["rows"][0]["sample_id"][:1], "S")
        self.assertTrue(all(int(r["sample_id"][1:]) > 8 for r in second["rows"]))
        self.assertEqual(second["extends"], first["sheet_wording_sha256"])

    def test_the_key_matches_the_sheet(self):
        sheet, key, _ = mls.build(self.s, self.tags, LEX, n=10, seed=9)
        self.assertEqual(ca.sheet_wording_sha(list(csv.DictReader(io.StringIO(sheet)))), key["sheet_wording_sha256"])

    def test_formula_like_cells_are_neutralised(self):
        self.assertEqual(mls.safe_cell("=HYPERLINK(1)"), "'=HYPERLINK(1)")
        self.assertEqual(mls.safe_cell("-5 points"), "'-5 points")
        self.assertEqual(mls.safe_cell("Change in FEV1"), "Change in FEV1")

    def test_tags_from_another_scope_are_refused(self):
        bad = {**self.tags, "entries": self.tags["entries"][:-1]}
        with self.assertRaises(ValueError):
            mls.build(self.s, bad, LEX, n=5, seed=1)

    def cli_args(self, *extra):
        sheet_dir, key_dir = Path(tempfile.mkdtemp(dir=self.root)), Path(tempfile.mkdtemp(dir=self.root))
        return sheet_dir / "s.csv", key_dir / "k.json", ["--snapshot", str(self.root / "snap"), "--tags", str(self.tags_path),
                                                         "--sheet", str(sheet_dir / "s.csv"), "--key", str(key_dir / "k.json"),
                                                         "--n", "10", "--seed", "11", *extra]

    def test_cli_writes_once_and_never_overwrites(self):
        sheet, _, args = self.cli_args()
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(mls.main(args), 0)
        first = sheet.read_bytes()
        with contextlib.redirect_stderr(io.StringIO()):
            self.assertEqual(mls.main(args), 2)
        self.assertEqual(sheet.read_bytes(), first)
        self.assertNotIn(b"\r\n", first)

    # ---- input binding, strata, encoding and the key's folder

    def test_tags_from_another_snapshot_or_lexicon_are_refused(self):
        for field in ("snapshot_sha256", "lexicon_sha256"):
            with self.subTest(field=field), self.assertRaises(ValueError):
                mls.build(self.s, {**self.tags, field: "0" * 64}, LEX, n=5, seed=1)

    def test_the_key_records_the_input_hashes(self):
        _, key, _ = mls.build(self.s, self.tags, LEX, n=5, seed=1)
        self.assertEqual(key["snapshot_sha256"], self.s.digest)
        self.assertEqual(key["lexicon_sha256"], LEX.sha256)
        self.assertEqual(key["tags_sha256"], mls.tags_digest(self.tags))

    def test_every_non_empty_stratum_gets_a_row_first(self):
        alloc = mls.allocate({"big": 100, "one": 1, "two": 2}, 5)
        self.assertEqual(sum(alloc.values()), 5)
        self.assertGreaterEqual(alloc["one"], 1)
        self.assertGreaterEqual(alloc["two"], 1)
        self.assertEqual(mls.allocate({"a": 3, "b": 1}, 10), {"a": 3, "b": 1})        # never more than a stratum holds

    def test_no_primary_outcome_entries_are_never_drawn(self):
        _, key, _ = mls.build(self.s, self.tags, LEX, n=50, seed=2)
        self.assertFalse([r for r in key["rows"] if r["entry_id"].endswith(":P0")])
        self.assertEqual(len(key["rows"]), 19)

    def test_the_sheet_is_written_as_utf8_with_a_signature(self):
        sheet, _, args = self.cli_args()
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(mls.main(args), 0)
        self.assertTrue(sheet.read_bytes().startswith(b"\xef\xbb\xbf"))

    def test_the_key_may_not_sit_beside_the_sheet(self):
        sheet, _, args = self.cli_args()
        args[args.index("--key") + 1] = str(sheet.parent / "k.json")
        with contextlib.redirect_stderr(io.StringIO()) as err:
            self.assertEqual(mls.main(args), 2)
        self.assertIn("same folder", err.getvalue())
        self.assertFalse(sheet.exists())

    def test_the_key_may_not_sit_below_the_sheets_folder(self):
        # the key's folder may not be inside the sheet's folder
        sheet, _, args = self.cli_args()
        for key in (sheet.parent / "sub" / "k.json", sheet.parent / "x" / ".." / "k3.json"):
            with self.subTest(key=str(key)):
                a = list(args)
                a[a.index("--key") + 1] = str(key)
                with contextlib.redirect_stderr(io.StringIO()) as err:
                    self.assertEqual(mls.main(a), 2)
                self.assertIn("same folder", err.getvalue())
                self.assertFalse(Path(key).exists())

    # ---- short draws and the terms file

    def test_a_short_draw_is_warned_about(self):
        _, _, args = self.cli_args()
        args[args.index("--n") + 1] = "50"
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()) as err:
            self.assertEqual(mls.main(args), 0)
        self.assertIn("WARNING: drew 19 of the 50 rows asked for", err.getvalue())
        self.assertIn("only 19 entries could be drawn", err.getvalue())

    def test_the_registry_terms_are_written_beside_the_sheet(self):
        sheet, _, args = self.cli_args()
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(mls.main(args), 0)
        terms = sheet.with_name(sheet.stem + ".TERMS.txt")
        text = terms.read_text(encoding="utf-8")
        self.assertIn("ClinicalTrials.gov", text)
        self.assertIn(ca.LICENCE_LINE, text)
        self.assertIn("for as long as the data are kept", text)
        self.assertIn(self.s.manifest["data_timestamp"], text)
        for word in ("unclassified", "fev1", "stratum", "NCT0"):
            self.assertNotIn(word, text)                                    # terms only: nothing from the rules or the key

    def test_no_aggregate_counts_without_report(self):
        _, _, args = self.cli_args()
        with contextlib.redirect_stdout(io.StringIO()) as out:
            self.assertEqual(mls.main(args), 0)
        for word in ("unclassified", "rare", "strata", "oversample"):
            self.assertNotIn(word, out.getvalue())
        _, _, args = self.cli_args("--report")
        with contextlib.redirect_stdout(io.StringIO()) as out:
            self.assertEqual(mls.main(args), 0)
        self.assertIn("unclassified", out.getvalue())

    def test_the_seed_is_required(self):
        _, _, args = self.cli_args()
        i = args.index("--seed")
        del args[i:i + 2]
        with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit) as cm:
            mls.main(args)
        self.assertEqual(cm.exception.code, 2)


if __name__ == "__main__":
    unittest.main()
