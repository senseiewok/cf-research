"""Tests for make_labelling_sheet.py on the SYNTHETIC snapshot: the sheet leaks no tag, rule, model result, NCT id or stratum; the
draw is reproducible for a seed; the oversamples are taken; an extension skips earlier rows; the CLI never overwrites. The
census-v1 design (CensusTest): coverage of the small classes and the unclassified entries, draw sizes and weights per stratum,
the partition, determinism, blindness, model tags changing the frame, the legacy design left as it was. No network.
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
from types import SimpleNamespace
from unittest import mock

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


class CensusTest(unittest.TestCase):
    """The census-v1 design. Synthetic snapshot: 19 drawable entries (the :P0 entry has no primary outcome). With small_max=1 the
    classes safety_tolerability, exacerbations and pharmacokinetics (two entries each) are large and every other class is small."""

    @classmethod
    def setUpClass(cls):
        cls._tmp = tempfile.TemporaryDirectory()
        cls.root = Path(cls._tmp.name)
        sf.build_cf_snapshot(cls.root / "snap")
        cls.s = snap.load(cls.root / "snap")
        cls.tags = lexicon.tag_snapshot(cls.s, LEX)
        cls.tags_path = cls.root / "tags.json"
        cls.tags_path.write_text(json.dumps(cls.tags), encoding="utf-8")
        cls.good = sf.good_model_tags()
        cls.good_bytes = json.dumps(cls.good).encode("utf-8")
        entries, _ = lexicon.in_scope_entries(cls.s, None)
        cls.entries = entries
        cls.drawable = sorted(e["entry_id"] for e in entries if not e.get("no_primary_outcome"))
        by_id = {t["entry_id"]: t for t in cls.tags["entries"]}
        cls.classes = {eid: [t["class"] for t in by_id[eid]["tags"]] for eid in cls.drawable}
        cls.status = {eid: by_id[eid]["status"] for eid in cls.drawable}

    @classmethod
    def tearDownClass(cls):
        cls._tmp.cleanup()

    def setUp(self):
        sf.block_network(self)

    def census(self, seed=1, **kw):
        return mls.build(self.s, self.tags, LEX, seed=seed, design="census-v1", **kw)

    def with_model_tags(self, tags, seed=1, **kw):
        raw = json.dumps(tags).encode("utf-8")
        return self.census(seed, model_tags=tags, model_tags_sha256=mls.model_tags_digest(raw), **kw)

    @staticmethod
    def by_stratum(key):
        out = {}
        for r in key["rows"]:
            out.setdefault(r["stratum"], []).append(r["entry_id"])
        return out

    def cli_args(self, *extra, seed="11"):
        sheet_dir, key_dir = Path(tempfile.mkdtemp(dir=self.root)), Path(tempfile.mkdtemp(dir=self.root))
        return sheet_dir / "s.csv", key_dir / "k.json", ["--snapshot", str(self.root / "snap"), "--tags", str(self.tags_path),
                                                         "--sheet", str(sheet_dir / "s.csv"), "--key", str(key_dir / "k.json"),
                                                         "--seed", seed, "--design", "census-v1", *extra]

    def legacy_args(self, args, *extra):
        """The census arguments with the design flag removed and the legacy ones given instead."""
        i = args.index("--design")
        return args[:i] + args[i + 2:] + list(extra)

    def run_cli(self, args):
        with contextlib.redirect_stdout(io.StringIO()) as out, contextlib.redirect_stderr(io.StringIO()) as err:
            code = mls.main(args)
        return code, out.getvalue(), err.getvalue()

    # ---- coverage, partition, sizes and weights

    def test_every_small_class_entry_and_every_unclassified_entry_is_drawn(self):
        _, key, _ = self.census(small_max=1, large_per_class=1)
        drawn = {r["entry_id"] for r in key["rows"]}
        sizes = {c: sum(c in cl for cl in self.classes.values()) for c in LEX.class_order}
        small = {c for c, k in sizes.items() if 0 < k <= 1 and c != "not_stated"}
        carrying = [eid for eid in self.drawable if small & set(self.classes[eid])]
        self.assertGreaterEqual(len(carrying), 10)                      # the check below is not vacuous
        for eid in carrying:
            self.assertIn(eid, drawn, eid)
        for eid in self.drawable:
            if self.status[eid] != "rule":
                self.assertIn(eid, drawn, eid)
        self.assertEqual(self.by_stratum(key)["census:unclassified"], ["NCT00000014:P3"])
        self.assertEqual(self.by_stratum(key)["census:not_stated"], ["NCT00000008:P1"])

    def test_strata_draw_sizes_and_weights(self):
        _, key, summary = self.census(small_max=1, large_per_class=1)
        self.assertEqual({k: (v["kind"], v["population"], v["drawn"], v["weight"]) for k, v in key["strata"].items()}, {
            "census:small": ("census", 12, 12, 1.0), "census:unclassified": ("census", 1, 1, 1.0),
            "census:not_stated": ("census", 1, 1, 1.0), "sample:safety_tolerability": ("sample", 2, 1, 2.0),
            "sample:exacerbations": ("sample", 1, 1, 1.0), "sample:pharmacokinetics": ("sample", 2, 1, 2.0)})
        self.assertEqual(summary["drawn"], 17)
        self.assertEqual(len(key["rows"]), 17)
        for name, v in key["strata"].items():
            self.assertEqual(len(self.by_stratum(key)[name]), v["drawn"], name)
            self.assertEqual(v["weight"], v["population"] / v["drawn"])
        # the draw never takes more than a stratum holds
        _, key, _ = self.census(small_max=1, large_per_class=40)
        self.assertEqual({k: v["drawn"] for k, v in key["strata"].items() if k.startswith("sample")},
                         {"sample:safety_tolerability": 2, "sample:exacerbations": 1, "sample:pharmacokinetics": 2})

    def test_with_the_default_small_max_the_synthetic_set_is_all_census(self):
        _, key, _ = self.census()
        self.assertEqual(set(key["strata"]), {"census:small", "census:unclassified", "census:not_stated"})
        self.assertEqual(sum(v["population"] for v in key["strata"].values()), 19)
        self.assertEqual({r["entry_id"] for r in key["rows"]}, set(self.drawable))

    def test_each_entry_is_in_exactly_one_stratum(self):
        for kw in ({}, {"small_max": 1, "large_per_class": 1}, {"small_max": 2, "large_per_class": 1}):
            with self.subTest(**kw):
                _, key, _ = self.census(**kw)
                ids = [r["entry_id"] for r in key["rows"]]
                self.assertEqual(len(ids), len(set(ids)))
                self.assertEqual(sum(v["population"] for v in key["strata"].values()), len(self.drawable))
                self.assertTrue(all(r["stratum"] in key["strata"] for r in key["rows"]))
                self.assertTrue(all(r["entry_id"] in self.drawable for r in key["rows"]))
                for v in key["strata"].values():
                    self.assertLessEqual(v["drawn"], v["population"])

    def test_the_key_fields(self):
        _, key, _ = self.census(small_max=1, large_per_class=2)
        self.assertEqual((key["design"], key["small_max"], key["large_per_class"], key["consensus"], key["model_tags_sha256"]),
                         ("census-v1", 1, 2, None, None))
        for field in ("seed", "snapshot_sha256", "lexicon_sha256", "tags_sha256", "sheet_wording_sha256",
                      "lexicon_version_used_for_the_draw", "strata", "rows"):
            self.assertIn(field, key)
        for field in ("n", "unclassified_share", "rare_max"):
            self.assertNotIn(field, key)
        self.assertEqual(key["tags_sha256"], mls.tags_digest(self.tags))
        self.assertEqual(sorted(key["rows"][0]), ["entry_id", "sample_id", "stratum"])
        self.assertEqual(sorted(r["sample_id"] for r in key["rows"]), [f"S{i:03d}" for i in range(1, len(key["rows"]) + 1)])

    def test_the_rarest_large_class_names_the_sample_stratum_and_ties_follow_the_lexicon(self):
        lex = SimpleNamespace(class_order=["a", "b", "c", "d", "e", "not_stated"])
        pool, tags = [], {}

        def add(classes, status="rule"):
            eid = f"E{len(pool) + 1:03d}"
            pool.append({"entry_id": eid})
            tags[eid] = {"entry_id": eid, "status": status, "tags": [{"class": c} for c in classes]}
            return eid
        for _ in range(5):
            add(["a"])
        for _ in range(50):
            add(["b"])
        both_bc = [add(["b", "c"]) for _ in range(3)]
        for _ in range(45):
            add(["c"])
        both_de = add(["d", "e"])
        for _ in range(44):
            add(["d"])
        for _ in range(44):
            add(["e"])
        unclassified, modelled, not_stated = add([], "unclassified"), add([], "unclassified"), add(["not_stated"], "not_stated")
        chosen, strata = mls.draw_census(pool, tags, lex, {modelled: [{"class": "a"}]}, seed=7, small_max=5, large_per_class=40)
        by = {}
        for c in chosen:
            by.setdefault(c["stratum"], []).append(c["entry_id"])
        # a holds 5 entries and the modelled one makes 6, so a is large; the modelled entry is then not census material
        self.assertEqual(sorted(strata), ["census:not_stated", "census:unclassified", "sample:a", "sample:b", "sample:c", "sample:d",
                                          "sample:e"])
        self.assertEqual((by["census:unclassified"], by["census:not_stated"]), ([unclassified], [not_stated]))
        self.assertEqual(strata["sample:a"], {"kind": "sample", "population": 6, "drawn": 6, "weight": 1.0})
        # b has 53 entries and c has 48: an entry of both is sampled for c
        self.assertEqual((strata["sample:b"]["population"], strata["sample:c"]["population"]), (50, 48))
        # d and e tie at 45: the entry of both is sampled for d, the earlier class in the lexicon
        self.assertEqual((strata["sample:d"]["population"], strata["sample:e"]["population"]), (45, 44))
        self.assertEqual({k: v["drawn"] for k, v in strata.items() if k in ("sample:b", "sample:c", "sample:d", "sample:e")},
                         {"sample:b": 40, "sample:c": 40, "sample:d": 40, "sample:e": 40})
        self.assertEqual((strata["sample:d"]["weight"], strata["sample:e"]["weight"]), (45 / 40, 44 / 40))
        self.assertEqual(sum(v["population"] for v in strata.values()), len(pool))
        drawn_ids = [c["entry_id"] for c in chosen]
        self.assertEqual(len(drawn_ids), len(set(drawn_ids)))
        for eid in both_bc + [both_de]:
            self.assertNotIn(eid, by.get("census:small", []))

    def test_a_small_class_censuses_its_whole_population(self):
        lex = SimpleNamespace(class_order=["a", "b", "not_stated"])
        pool = [{"entry_id": f"E{i:03d}"} for i in range(260)]
        tags = {e["entry_id"]: {"entry_id": e["entry_id"], "status": "rule", "tags": [{"class": "a" if i < 30 else "b"}]}
                for i, e in enumerate(pool)}
        chosen, strata = mls.draw_census(pool, tags, lex, {}, seed=1, small_max=30, large_per_class=40)
        self.assertEqual(strata["census:small"], {"kind": "census", "population": 30, "drawn": 30, "weight": 1.0})
        self.assertEqual(strata["sample:b"], {"kind": "sample", "population": 230, "drawn": 40, "weight": 230 / 40})
        self.assertEqual(len(chosen), 70)
        tags["E030"]["tags"] = [{"class": "a"}]               # one more entry makes the class large: "small_max or fewer" is the rule
        _, strata = mls.draw_census(pool, tags, lex, {}, seed=1, small_max=30, large_per_class=40)
        self.assertEqual(sorted(strata), ["sample:a", "sample:b"])

    def test_a_class_carried_twice_by_one_entry_counts_once_as_in_the_gate(self):
        lex = SimpleNamespace(class_order=["a", "b", "not_stated"])
        pool = [{"entry_id": f"E{i:03d}"} for i in range(40)]
        tags = {e["entry_id"]: {"entry_id": e["entry_id"], "status": "rule", "tags": [{"class": "a" if i < 3 else "b"}]}
                for i, e in enumerate(pool)}
        tags["E000"]["tags"] = [{"class": "a"}, {"class": "a"}]  # two rule tags with one class: a is carried by 3 entries, not 4
        tags["E039"]["tags"] = [{"class": "b"}, {"class": "zz_not_in_the_order"}]   # a class outside the class order: no KeyError
        chosen, strata = mls.draw_census(pool, tags, lex, {}, seed=1, small_max=3, large_per_class=40)
        self.assertEqual(strata["census:small"]["population"], 4)     # a (3 entries) and the lone unknown class are both small
        # the gate derives the same strata from the same frame
        frame = {e["entry_id"]: {"status": "rule", "classes": [t["class"] for t in tags[e["entry_id"]]["tags"]]} for e in pool}
        derived = ca.census_strata(frame, lex, 3)
        self.assertEqual({c["entry_id"]: c["stratum"] for c in chosen}, derived)
        self.assertEqual(derived["E000"], "census:small")

    # ---- model tags

    def test_model_tags_change_the_frame(self):
        _, plain, _ = self.census()
        self.assertEqual(plain["strata"]["census:unclassified"]["population"], 1)
        _, modelled, _ = self.with_model_tags(self.good)
        self.assertNotIn("census:unclassified", modelled["strata"])
        self.assertEqual(modelled["strata"]["census:small"]["population"], plain["strata"]["census:small"]["population"] + 1)
        self.assertEqual(modelled["model_tags_sha256"], mls.model_tags_digest(self.good_bytes))
        self.assertRegex(modelled["model_tags_sha256"], r"^[0-9a-f]{64}$")
        self.assertIsNone(plain["model_tags_sha256"])
        # the class a model gives counts towards that class's size: fev1 goes from one entry to two and is then large
        tags = [{"entry_id": "NCT00000014:P3", "class": "fev1", "quote": "return to school"}]
        _, shifted, _ = self.with_model_tags(tags, small_max=1, large_per_class=1)
        _, base, _ = self.census(small_max=1, large_per_class=1)
        self.assertNotIn("census:unclassified", shifted["strata"])
        self.assertEqual(shifted["strata"]["sample:fev1"], {"kind": "sample", "population": 2, "drawn": 1, "weight": 2.0})
        self.assertEqual(shifted["strata"]["census:small"]["population"], base["strata"]["census:small"]["population"] - 1)

    def test_rejected_model_tags_leave_the_frame_alone(self):
        _, plain, _ = self.census()
        bad = sf.canary_model_tags() + sf.stray_model_tags()
        _, key, summary = self.with_model_tags(bad)
        self.assertEqual(key["strata"], plain["strata"])
        self.assertEqual((summary["model_tags_accepted"], summary["model_tags_rejected"]), (0, 5))
        self.assertIsNotNone(key["model_tags_sha256"])

    def test_model_tags_for_an_entry_the_rules_decided_are_ignored(self):
        _, plain, _ = self.census()
        _, key, summary = self.with_model_tags([{"entry_id": "NCT00000001:P1", "class": "other", "quote": "FEV1"}])
        self.assertEqual(key["strata"], plain["strata"])
        self.assertEqual(summary["model_tags_rejected"], 1)

    def test_model_tags_and_their_hash_go_together(self):
        with self.assertRaises(ValueError):
            self.census(model_tags=self.good)
        with self.assertRaises(ValueError):
            self.census(model_tags_sha256="0" * 64)

    # ---- determinism and blindness

    def test_the_same_seed_gives_byte_identical_output(self):
        for kw in ({}, {"small_max": 1, "large_per_class": 1}):
            with self.subTest(**kw):
                a, b = self.census(seed=42, **kw), self.census(seed=42, **kw)
                self.assertEqual(a[0], b[0])
                self.assertEqual(json.dumps(a[1], indent=1), json.dumps(b[1], indent=1))
        a, b = self.with_model_tags(self.good, seed=42), self.with_model_tags(self.good, seed=42)
        self.assertEqual((a[0], a[1]), (b[0], b[1]))
        first, second = self.cli_args("--small-max", "1", "--large-per-class", "1"), self.cli_args("--small-max", "1", "--large-per-class", "1")
        for _, _, args in (first, second):
            self.assertEqual(self.run_cli(args)[0], 0)
        self.assertEqual(first[0].read_bytes(), second[0].read_bytes())
        self.assertEqual(first[1].read_bytes(), second[1].read_bytes())

    def test_a_different_seed_changes_only_the_sample_strata(self):
        kw = {"small_max": 1, "large_per_class": 1}
        _, first, _ = self.census(seed=1, **kw)
        base = self.by_stratum(first)
        changed = False
        for seed in range(2, 30):
            _, other, _ = self.census(seed=seed, **kw)
            now = self.by_stratum(other)
            self.assertEqual(other["strata"], first["strata"])
            for name in base:
                if name.startswith("census:"):
                    self.assertEqual(sorted(now[name]), sorted(base[name]), name)
                elif now[name] != base[name]:
                    changed = True
        self.assertTrue(changed, "no seed changed a sample stratum")
        _, shuffled, _ = self.census(seed=2, **kw)
        self.assertNotEqual([r["entry_id"] for r in shuffled["rows"]], [r["entry_id"] for r in first["rows"]])

    def test_the_sheet_is_blind(self):
        for kw in ({}, {"small_max": 1, "large_per_class": 1}):
            with self.subTest(**kw):
                sheet, key, _ = self.with_model_tags(self.good, **kw)
                rows = list(csv.DictReader(io.StringIO(sheet)))
                self.assertEqual(list(rows[0].keys()), mls.COLUMNS)
                for r in rows:
                    self.assertEqual((r["label_classes"], r["label_composite"], r["label_notes"]), ("", "", ""))
                leaks = [c for c in LEX.class_order if "_" in c] + [r["id"] for r in LEX.rules + LEX.composite]
                leaks += list(key["strata"]) + ["census", "sample:", "stratum", "unclassified", "by_rule", "model", "NCT0", "weight"]
                for word in leaks:
                    self.assertNotIn(word, sheet, f"the sheet must not contain {word!r}")
                self.assertEqual([r["sample_id"] for r in rows], [r["sample_id"] for r in key["rows"]])

    def test_the_key_matches_the_sheet_and_the_rows_are_shuffled(self):
        sheet, key, _ = self.census()
        self.assertEqual(ca.sheet_wording_sha(list(csv.DictReader(io.StringIO(sheet)))), key["sheet_wording_sha256"])
        in_row_order = [r["entry_id"] for r in key["rows"]]
        self.assertNotEqual(in_row_order, sorted(in_row_order))               # the order is the seeded shuffle, not the canonical order

    def test_tags_from_another_snapshot_or_lexicon_are_refused_in_census_mode(self):
        for field in ("snapshot_sha256", "lexicon_sha256"):
            with self.subTest(field=field), self.assertRaises(ValueError):
                mls.build(self.s, {**self.tags, field: "0" * 64}, LEX, seed=1, design="census-v1")

    # ---- the legacy design is unchanged; the cap and the options

    def test_the_legacy_key_has_no_design_field_and_the_cap_still_fires(self):
        _, key, _ = mls.build(self.s, self.tags, LEX, n=10, seed=3)
        for field in ("design", "strata", "small_max", "large_per_class", "model_tags_sha256", "consensus"):
            self.assertNotIn(field, key)
        self.assertEqual(mls.MAX_ROWS, 150)
        for n in ("151", "0"):
            _, _, args = self.cli_args()
            code, _, err = self.run_cli(self.legacy_args(args, "--n", n))
            self.assertEqual(code, 2)
            self.assertIn("--n must be 1 to 150", err)

    def test_the_cap_does_not_apply_to_census_mode(self):
        with mock.patch.object(mls, "MAX_ROWS", 5):
            _, key, args = self.cli_args()
            code, out, err = self.run_cli(args)
            self.assertEqual(code, 0, err)
            self.assertEqual(len(json.loads(key.read_text(encoding="utf-8"))["rows"]), 19)     # more than the patched cap of 5
            _, _, args = self.cli_args()
            self.assertEqual(self.run_cli(self.legacy_args(args, "--n", "10"))[0], 2)          # the legacy design still obeys it
        # and the draw function itself has no size limit
        lex = SimpleNamespace(class_order=["a", "not_stated"])
        pool = [{"entry_id": f"E{i:04d}"} for i in range(400)]
        tags = {e["entry_id"]: {"entry_id": e["entry_id"], "status": "unclassified", "tags": []} for e in pool}
        chosen, strata = mls.draw_census(pool, tags, lex, {}, seed=1, small_max=30, large_per_class=40)
        self.assertEqual((len(chosen), strata["census:unclassified"]["drawn"]), (400, 400))

    def test_census_refuses_n_and_the_legacy_only_options(self):
        for extra in (["--n", "10"], ["--unclassified-share", "0.5"], ["--rare-max", "3"], ["--already", str(self.tags_path)]):
            with self.subTest(extra=extra):
                sheet, key, args = self.cli_args(*extra)
                code, _, err = self.run_cli(args)
                self.assertEqual(code, 2)
                self.assertIn("cannot be used with --design census-v1", err)
                self.assertFalse(sheet.exists() or key.exists())
        with self.assertRaises(ValueError):
            self.census(n=10)

    def test_legacy_refuses_the_census_only_options(self):
        for extra in (["--model-tags", str(self.tags_path)], ["--small-max", "5"], ["--large-per-class", "5"]):
            with self.subTest(extra=extra):
                _, _, args = self.cli_args()
                code, _, err = self.run_cli(self.legacy_args(args, "--n", "10", *extra))
                self.assertEqual(code, 2)
                self.assertIn("cannot be used with --design legacy", err)

    def test_census_sizes_must_be_positive(self):
        for extra in (["--small-max", "0"], ["--large-per-class", "0"]):
            _, _, args = self.cli_args(*extra)
            self.assertEqual(self.run_cli(args)[0], 2)

    # ---- the CLI

    def test_cli_writes_the_census_once_and_never_overwrites(self):
        sheet, keyf, args = self.cli_args("--small-max", "1", "--large-per-class", "1")
        code, out, err = self.run_cli(args)
        self.assertEqual((code, err), (0, ""))
        self.assertIn("drew 17 rows (seed 11)", out)
        self.assertTrue(sheet.read_bytes().startswith(b"\xef\xbb\xbf"))
        self.assertTrue(sheet.with_name("s.TERMS.txt").exists())
        key = json.loads(keyf.read_text(encoding="utf-8"))
        self.assertEqual((key["design"], key["small_max"], key["large_per_class"]), ("census-v1", 1, 1))
        first = sheet.read_bytes()
        self.assertEqual(self.run_cli(args)[0], 2)
        self.assertEqual(sheet.read_bytes(), first)

    def test_cli_model_tags_file_sets_the_hash_and_changes_the_frame(self):
        mt = self.root / "model-tags.json"
        mt.write_bytes(self.good_bytes)
        _, keyf, args = self.cli_args("--model-tags", str(mt))
        self.assertEqual(self.run_cli(args)[0], 0)
        key = json.loads(keyf.read_text(encoding="utf-8"))
        self.assertEqual(key["model_tags_sha256"], mls.model_tags_digest(self.good_bytes))
        self.assertNotIn("census:unclassified", key["strata"])
        _, keyf, args = self.cli_args()
        self.assertEqual(self.run_cli(args)[0], 0)
        key = json.loads(keyf.read_text(encoding="utf-8"))
        self.assertIsNone(key["model_tags_sha256"])
        self.assertIn("census:unclassified", key["strata"])
        bad = self.root / "bad.json"                                           # a file that is not JSON is an input error
        bad.write_text("not json", encoding="utf-8")
        _, _, args = self.cli_args("--model-tags", str(bad))
        code, out, _ = self.run_cli(args)
        self.assertEqual(code, 1)
        self.assertIn("ERROR", out)

    def test_cli_consensus_report_fills_the_keys_consensus_record(self):
        mt = self.root / "kept-model-tags.json"
        mt.write_bytes(self.good_bytes)
        report = {"kind": "trial-atlas consensus report (counts only)", "rule_declared": "2026-10-11", "rule": "text",
                  "proposed": 7, "kept": 1, "dropped": 6, "dropped_disagree": 6, "labeller_family": "gemma"}
        rp = self.root / "consensus-report.json"
        rp.write_text(json.dumps(report, indent=1) + "\n", encoding="utf-8")
        _, keyf, args = self.cli_args("--model-tags", str(mt), "--consensus-report", str(rp))
        code, out, err = self.run_cli(args)
        self.assertEqual(code, 0, out + err)
        key = json.loads(keyf.read_text(encoding="utf-8"))
        self.assertEqual(key["consensus"], {"rule_declared": "2026-10-11", "proposed": 7, "kept": 1, "dropped": 6,
                                            "report_sha256": __import__("hashlib").sha256(rp.read_bytes()).hexdigest()})
        # without the option the record stays null
        _, keyf, args = self.cli_args("--model-tags", str(mt))
        self.assertEqual(self.run_cli(args)[0], 0)
        self.assertIsNone(json.loads(keyf.read_text(encoding="utf-8"))["consensus"])
        # a report lacking a field (or with a wrong type) is an input error, nothing written
        for broken in ({k: v for k, v in report.items() if k != "kept"}, {**report, "dropped": "6"}, [report]):
            with self.subTest(broken=str(broken)[:40]):
                bad = Path(tempfile.mkdtemp(dir=self.root)) / "consensus-report.json"
                bad.write_text(json.dumps(broken), encoding="utf-8")
                sheet, keyf, args = self.cli_args("--model-tags", str(mt), "--consensus-report", str(bad))
                code, out, _ = self.run_cli(args)
                self.assertEqual(code, 1, out)
                self.assertIn("ERROR", out)
                self.assertFalse(sheet.exists() or keyf.exists())
        # it needs --model-tags, and the legacy design refuses it
        _, _, args = self.cli_args("--consensus-report", str(rp))
        code, _, err = self.run_cli(args)
        self.assertEqual(code, 2)
        self.assertIn("--consensus-report goes with --model-tags", err)
        _, _, args = self.cli_args()
        code, _, err = self.run_cli(self.legacy_args(args, "--n", "10", "--consensus-report", str(rp)))
        self.assertEqual(code, 2)
        self.assertIn("cannot be used with --design legacy", err)
        with self.assertRaises(ValueError):
            self.census(consensus={"rule_declared": "2026-10-11", "proposed": 1, "kept": 1, "dropped": 0, "report_sha256": "0" * 64})

    # ---- --sample-exclude: development entries kept out of the sample strata only

    def sample_ids(self, key, name):
        return sorted(r["entry_id"] for r in key["rows"] if r["stratum"] == name)

    def test_without_exclusion_the_draw_is_unchanged(self):
        # sha256 of draw_census's output on this hand-made frame, recorded before --sample-exclude existed (2026-10-11)
        lex = SimpleNamespace(class_order=["a", "b", "c", "not_stated"])
        pool = [{"entry_id": f"E{i:03d}"} for i in range(300)]
        tags = {e["entry_id"]: {"entry_id": e["entry_id"], "status": "rule",
                                "tags": [{"class": "a" if i < 20 else ("b" if i < 150 else "c")}]} for i, e in enumerate(pool)}
        digest = lambda out: __import__("hashlib").sha256(json.dumps(list(out), sort_keys=True).encode()).hexdigest()  # noqa: E731
        before = "4763a8654d2ec569aa737faef680de135a29fbd9a781482289a765b115ad4ff2"
        self.assertEqual(digest(mls.draw_census(pool, tags, lex, {}, seed=9, small_max=30, large_per_class=40)), before)
        self.assertEqual(digest(mls.draw_census(pool, tags, lex, {}, seed=9, small_max=30, large_per_class=40,
                                                sample_exclude=frozenset())), before)
        # through build: no exclusion gives the same sheet, rows and strata, and a null sample_exclusion
        a, b = self.census(seed=5, small_max=1, large_per_class=1), self.census(seed=5, small_max=1, large_per_class=1, sample_exclude=None)
        self.assertEqual((a[0], a[1]), (b[0], b[1]))
        self.assertIsNone(a[1]["sample_exclusion"])

    def test_excluded_entries_leave_only_the_sample_strata_and_reduce_their_weight(self):
        lex = SimpleNamespace(class_order=["a", "b", "not_stated"])
        pool = [{"entry_id": f"E{i:03d}"} for i in range(260)]
        tags = {e["entry_id"]: {"entry_id": e["entry_id"], "status": "rule", "tags": [{"class": "a" if i < 30 else "b"}]}
                for i, e in enumerate(pool)}
        dev = frozenset([f"E{i:03d}" for i in range(0, 5)] + [f"E{i:03d}" for i in range(30, 60)])   # 5 census, 30 sampled
        chosen, strata = mls.draw_census(pool, tags, lex, {}, seed=1, small_max=30, large_per_class=40, sample_exclude=dev)
        self.assertEqual(strata["census:small"], {"kind": "census", "population": 30, "drawn": 30, "weight": 1.0})
        self.assertEqual(strata["sample:b"], {"kind": "sample", "population": 200, "drawn": 40, "weight": 200 / 40})
        drawn = {c["entry_id"]: c["stratum"] for c in chosen}
        self.assertTrue(all(drawn[f"E{i:03d}"] == "census:small" for i in range(5)))      # still drawn in the census stratum
        self.assertFalse(dev & {e for e, st in drawn.items() if st == "sample:b"})
        # the class sizes still count the excluded entries: with them a stays small (30); the frame is not shrunk
        self.assertNotIn("sample:a", strata)
        # a sample stratum whose whole population is excluded is left out (its weight would be undefined)
        all_b = frozenset(f"E{i:03d}" for i in range(30, 260))
        _, strata = mls.draw_census(pool, tags, lex, {}, seed=1, small_max=30, large_per_class=40, sample_exclude=all_b)
        self.assertEqual(sorted(strata), ["census:small"])

    def test_build_with_exclusion_records_it_and_the_gate_binds_it(self):
        _, base, _ = self.census(small_max=1, large_per_class=1)
        safety = sorted(e for e in self.drawable if "safety_tolerability" in self.classes[e])
        census_small = self.sample_ids(base, "census:small")
        dev = [safety[0], census_small[0], "NCT09999999:P1"]                    # one sampled, one censused, one not in the frame
        _, key, summary = self.census(small_max=1, large_per_class=1, sample_exclude=dev)
        self.assertEqual(key["sample_exclusion"], {"count": 2, "sha256": ca.sample_exclusion_digest(dev)})
        self.assertEqual(key["sample_exclusion"]["sha256"],
                         __import__("hashlib").sha256("\n".join(sorted(dev)).encode("utf-8")).hexdigest())
        self.assertEqual((summary["excluded_in_frame"], summary["excluded_not_in_frame"], summary["excluded_in_census"]), (2, 1, 1))
        self.assertEqual(key["strata"]["sample:safety_tolerability"], {"kind": "sample", "population": 1, "drawn": 1, "weight": 1.0})
        self.assertEqual(self.sample_ids(key, "sample:safety_tolerability"), [safety[1]])
        self.assertEqual(key["strata"]["census:small"], base["strata"]["census:small"])
        self.assertIn(census_small[0], self.sample_ids(key, "census:small"))
        self.assertEqual(sum(v["population"] for v in key["strata"].values()), len(self.drawable) - 1)
        # the gate re-derives the same frame with the same file, and refuses it without the file or with another one
        self.assertEqual(ca.census_binding_problems(key, self.tags["entries"], [], LEX, None, dev), [])
        self.assertIn("no --frozen-exclude file was given", " ".join(ca.census_binding_problems(key, self.tags["entries"], [], LEX, None)))
        self.assertIn("hash differs", " ".join(ca.census_binding_problems(key, self.tags["entries"], [], LEX, None, dev[:2])))
        # a key drawn without exclusion refuses an exclusion file
        self.assertIn("drawn with no sample exclusion", " ".join(ca.census_binding_problems(base, self.tags["entries"], [], LEX, None, dev)))
        # an excluded entry drawn into its sample stratum (the honest row swapped for it) is refused
        for r in key["rows"]:
            if r["entry_id"] == safety[1]:
                r["entry_id"] = safety[0]
        self.assertIn("row(s) of a sample stratum are excluded (development) entries",
                      " ".join(ca.census_binding_problems(key, self.tags["entries"], [], LEX, None, dev)))

    def test_the_same_seed_and_exclusion_give_byte_identical_output(self):
        dev = sorted(e for e in self.drawable if "pharmacokinetics" in self.classes[e])[:1]
        a = self.census(seed=42, small_max=1, large_per_class=1, sample_exclude=dev)
        b = self.census(seed=42, small_max=1, large_per_class=1, sample_exclude=list(reversed(dev + dev)))
        self.assertEqual(a[0], b[0])
        self.assertEqual(json.dumps(a[1], indent=1), json.dumps(b[1], indent=1))   # order and duplicates do not change the hash
        devf = self.root / "dev-ids.json"
        devf.write_text(json.dumps(dev), encoding="utf-8")
        first, second = (self.cli_args("--small-max", "1", "--large-per-class", "1", "--sample-exclude", str(devf)) for _ in range(2))
        for _, _, args in (first, second):
            self.assertEqual(self.run_cli(args)[0], 0)
        self.assertEqual(first[0].read_bytes(), second[0].read_bytes())
        self.assertEqual(first[1].read_bytes(), second[1].read_bytes())
        self.assertEqual(json.loads(first[1].read_text(encoding="utf-8"))["sample_exclusion"]["count"], 1)

    def test_cli_sample_exclude_inputs(self):
        devf = self.root / "dev-report.json"
        devf.write_text(json.dumps([self.drawable[0], "NCT09999999:P1", "NCT09999998:P1"]), encoding="utf-8")
        _, _, args = self.cli_args("--sample-exclude", str(devf), "--report")
        code, out, _ = self.run_cli(args)
        self.assertEqual(code, 0, out)
        self.assertIn("sample exclusion: 1 ids in the frame (1 of them in census strata, still drawn), 2 ids not in the frame (ignored)",
                      out)
        for broken in ({"ids": []}, [1, 2], "x"):
            with self.subTest(broken=broken):
                bad = Path(tempfile.mkdtemp(dir=self.root)) / "dev.json"
                bad.write_text(json.dumps(broken), encoding="utf-8")
                sheet, keyf, args = self.cli_args("--sample-exclude", str(bad))
                code, out, _ = self.run_cli(args)
                self.assertEqual(code, 1, out)
                self.assertIn("must hold a JSON list of entry ids", out)
                self.assertFalse(sheet.exists() or keyf.exists())
        _, _, args = self.cli_args()
        code, _, err = self.run_cli(self.legacy_args(args, "--n", "10", "--sample-exclude", str(devf)))
        self.assertEqual(code, 2)
        self.assertIn("cannot be used with --design legacy", err)
        with self.assertRaises(ValueError):
            mls.build(self.s, self.tags, LEX, n=10, seed=3, sample_exclude=[])

    def test_report_prints_counts_only(self):
        _, _, args = self.cli_args("--small-max", "1", "--large-per-class", "1", "--report")
        code, out, _ = self.run_cli(args)
        self.assertEqual(code, 0)
        self.assertIn("census:small: population 12, drawn 12, weight 1", out)
        self.assertIn("sample:safety_tolerability: population 2, drawn 1, weight 2", out)
        for e in self.entries:
            for field in ("measure", "description"):
                if e[field]:
                    self.assertNotIn(e[field], out)
        self.assertNotIn("NCT0", out)
        _, _, args = self.cli_args("--small-max", "1")
        code, out, _ = self.run_cli(args)
        self.assertNotIn("population", out)                                    # no counts without --report


if __name__ == "__main__":
    unittest.main()
