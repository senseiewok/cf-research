"""Tests for check_atlas.py on the SYNTHETIC snapshot: drift on a perturbed page, tags file or counts file must fail; model tags with
inexact quotes are rejected; the metrics and the Wilson interval; the publication stop rules; the negative controls, including
controls that must fail when the verifier, the scope or the lexicon is broken. No network.
Usage: python -m unittest -v (inside tools/trial_atlas)."""
import contextlib
import csv
import io
import json
import shutil
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import check_atlas as ca  # noqa: E402
import fetch_snapshot as fs  # noqa: E402
import lexicon  # noqa: E402
import make_labelling_sheet as mls  # noqa: E402
import snapshot as snap  # noqa: E402
import synthetic_fixtures as sf  # noqa: E402

LEX = lexicon.load()


def run(argv, **kw):
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        code = ca.main([str(a) for a in argv], **kw)
    return code, out.getvalue()


class Fixture(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._tmp = tempfile.TemporaryDirectory()
        cls.root = Path(cls._tmp.name)
        cls.snap_dir = cls.root / "snap"
        sf.build_cf_snapshot(cls.snap_dir)
        cls.s = snap.load(cls.snap_dir)
        cls.tags_path = cls.root / "tags.json"
        cls.tags_path.write_text(json.dumps(lexicon.tag_snapshot(cls.s, LEX)), encoding="utf-8")
        cls.model_path = cls.root / "model-tags.json"
        cls.model_path.write_text(json.dumps(sf.good_model_tags()), encoding="utf-8")
        cls.counts_path = cls.root / "counts.json"
        code, out = run(["--snapshot", cls.snap_dir, "--tags", cls.tags_path, "--model-tags", cls.model_path,
                         "--write-counts", cls.counts_path])
        assert code == 0, out
        cls.explained = cls.root / "explained.json"
        cls.explained.write_text(json.dumps({"NCT00000014": "SYNTHETIC: hand check found it", "NCT00000015": "SYNTHETIC: same",
                                             "NCT00000017": "SYNTHETIC: condition spelled differently"}), encoding="utf-8")
        cls.canary_path = cls.root / "canary.json"
        cls.canary_path.write_text(json.dumps(sf.canary_model_tags()), encoding="utf-8")
        cls.planted_path = cls.root / "planted.json"
        cls.planted_path.write_text(json.dumps(sf.PLANTED_NON_CF_IDS), encoding="utf-8")
        cls.control_dir = cls.root / "control"
        sf.build_noncf_snapshot(cls.control_dir)

    def controls(self):
        """Controls on named inputs (here the synthetic ones, named explicitly), as the publication gate requires."""
        return ["--negative-controls", "--canary", self.canary_path, "--planted", self.planted_path,
                "--control-snapshot", self.control_dir]

    @classmethod
    def tearDownClass(cls):
        cls._tmp.cleanup()

    def setUp(self):
        sf.block_network(self)
        self.work = Path(tempfile.mkdtemp(dir=self.root))

    def base(self):
        return ["--snapshot", self.snap_dir, "--tags", self.tags_path, "--counts", self.counts_path, "--model-tags", self.model_path]

    def frozen(self, label=None, n=20):
        """A filled-in frozen set. By default every row is labelled with the classes the taggers gave (a perfect labeller)."""
        tags = json.loads(self.tags_path.read_text(encoding="utf-8"))
        sheet, key, _ = mls.build(self.s, tags, LEX, n=n, seed=7)
        by_sample = {r["sample_id"]: r["entry_id"] for r in key["rows"]}
        final = {e["entry_id"]: [t["class"] for t in e["tags"]] for e in tags["entries"]}
        final["NCT00000014:P3"] = ["other"]
        rows = list(csv.DictReader(io.StringIO(sheet)))
        for r in rows:
            eid = by_sample[r["sample_id"]]
            r["label_classes"] = ";".join(label(eid, final[eid]) if label else final[eid])
        sp, kp = self.work / "sheet.csv", self.work / "key.json"
        with sp.open("w", encoding="utf-8", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=mls.COLUMNS)
            w.writeheader()
            w.writerows(rows)
        kp.write_text(json.dumps(key), encoding="utf-8")
        return ["--frozen-sheet", sp, "--frozen-key", kp]


class IntegrityTest(Fixture):
    def test_unchanged_inputs_reproduce(self):
        code, out = run(self.base() + ["--drift-only"])
        self.assertEqual(code, 0, out)
        self.assertIn("PASS integrity", out)
        self.assertIn("studies in scope: 12   outcome entries: 20", out)

    def test_a_perturbed_page_fails(self):
        copy = self.work / "snap"
        shutil.copytree(self.snap_dir, copy)
        page = copy / "condition" / "page-0001.json"
        page.write_text(page.read_text(encoding="utf-8").replace("sweat chloride", "sweat sodium"), encoding="utf-8")
        code, out = run(["--snapshot", copy, "--tags", self.tags_path, "--counts", self.counts_path, "--drift-only"])
        self.assertEqual(code, 1, out)
        self.assertIn("does not match its sha256", out)

    def test_a_perturbed_counts_file_fails(self):
        counts = json.loads(self.counts_path.read_text(encoding="utf-8"))
        counts["studies_by_class"]["fev1"] += 1
        p = self.work / "counts.json"
        p.write_text(json.dumps(counts), encoding="utf-8")
        code, out = run(["--snapshot", self.snap_dir, "--tags", self.tags_path, "--counts", p, "--model-tags", self.model_path,
                         "--drift-only"])
        self.assertEqual(code, 1, out)
        self.assertIn("drift /studies_by_class/fev1: recomputed 1, committed 2", out)

    def test_a_perturbed_tags_file_fails(self):
        tags = json.loads(self.tags_path.read_text(encoding="utf-8"))
        tags["entries"][0]["tags"][0]["class"] = "spirometry_other"
        p = self.work / "tags.json"
        p.write_text(json.dumps(tags), encoding="utf-8")
        code, out = run(["--snapshot", self.snap_dir, "--tags", p, "--counts", self.counts_path, "--model-tags", self.model_path,
                         "--drift-only"])
        self.assertEqual(code, 1, out)
        self.assertIn("does not reproduce the tags file (1 entr(ies) differ)", out)

    def test_a_tags_file_from_another_snapshot_fails(self):
        tags = json.loads(self.tags_path.read_text(encoding="utf-8"))
        tags["snapshot_sha256"] = "0" * 64
        p = self.work / "tags.json"
        p.write_text(json.dumps(tags), encoding="utf-8")
        code, out = run(["--snapshot", self.snap_dir, "--tags", p, "--counts", self.counts_path, "--model-tags", self.model_path,
                         "--drift-only"])
        self.assertEqual(code, 1, out)
        self.assertIn("different snapshot", out)

    def test_dropping_the_model_tags_is_drift(self):
        code, out = run(["--snapshot", self.snap_dir, "--tags", self.tags_path, "--counts", self.counts_path, "--drift-only"])
        self.assertEqual(code, 1, out)
        self.assertIn("/entries_by_class/other", out)

    def test_a_snapshot_without_manifest_is_refused(self):
        copy = self.work / "snap"
        shutil.copytree(self.snap_dir, copy)
        (copy / "manifest.json").unlink()
        code, out = run(["--snapshot", copy, "--tags", self.tags_path, "--counts", self.counts_path])
        self.assertEqual(code, 1, out)

    def test_counts_are_from_the_data(self):
        c = json.loads(self.counts_path.read_text(encoding="utf-8"))
        self.assertEqual(c["studies"], 12)
        self.assertEqual(c["entries"], 20)
        self.assertEqual(c["entries_by_status"], {"rule": 17, "model": 1, "not_stated": 2, "unclassified": 0})
        self.assertEqual(c["studies_by_class"]["exacerbations"], 2)       # NCT00000004 and NCT00000012
        self.assertEqual(c["studies_by_class"]["fev1"], 1)                # never the vague NCT00000007
        self.assertEqual(c["studies_by_class"]["lung_function_unspecified"], 1)
        self.assertEqual(c["flags"]["composite"], 1)
        self.assertEqual(c["studies_by_class_and_start_year"]["exacerbations"], {"2005": 1, "2027": 1})
        self.assertEqual(c["safety_subtypes"]["serious_adverse_events"], 1)
        self.assertEqual(c["scope"]["excluded_by_rule"]["X5"], 1)
        self.assertEqual(c["shares"]["other"], 0.05)

    # ---- third round (after adff5aa). Each case failed before its fix.

    def test_the_registry_terms_travel_with_the_counts(self):
        c = json.loads(self.counts_path.read_text(encoding="utf-8"))
        t = c["registry_terms"]
        self.assertEqual(t["source"], "ClinicalTrials.gov")
        self.assertEqual(t["data_processed_by_registry"], self.s.manifest["data_timestamp"])
        self.assertEqual(t["snapshot_fetched_at"], self.s.manifest["fetched_at"])
        self.assertEqual(t["terms_url"], "https://clinicaltrials.gov/about-site/terms-conditions")
        self.assertEqual(t["terms_last_updated"], "2023-01-31")
        mods = " ".join(t["modifications"])
        self.assertIn(LEX.version, mods)
        self.assertIn("X5 1", mods)                                       # each exclusion rule with its count
        self.assertIn("re-saved", mods)
        self.assertIn("clipped", mods)
        self.assertEqual(t["licence"], ca.LICENCE_LINE)
        self.assertIn("for as long as the data are kept", t["retention"])
        code, out = run(self.base() + ["--drift-only"])
        self.assertIn("source: ClinicalTrials.gov", out)
        self.assertIn(ca.LICENCE_LINE, out)

    def test_the_drift_check_covers_the_terms_block(self):
        counts = json.loads(self.counts_path.read_text(encoding="utf-8"))
        counts["registry_terms"]["licence"] = "anything goes"
        p = self.work / "counts.json"
        p.write_text(json.dumps(counts), encoding="utf-8")
        code, out = run(["--snapshot", self.snap_dir, "--tags", self.tags_path, "--counts", p, "--model-tags", self.model_path,
                         "--drift-only"])
        self.assertEqual(code, 1, out)
        self.assertIn("drift /registry_terms/licence", out)

    def test_split_counts_by_actual_and_planned_start_and_first_posted_year(self):
        c = json.loads(self.counts_path.read_text(encoding="utf-8"))
        self.assertEqual(c["studies_by_class_and_actual_start_year"]["exacerbations"], {"2005": 1})
        self.assertEqual(c["studies_by_class_and_planned_start_year"]["exacerbations"], {"2027": 1})
        self.assertEqual(c["studies_by_class_and_first_posted_year"]["exacerbations"], {"2019": 1, "2026": 1})

    # ---- round 4 (after 47adc05). Each case failed before its fix.

    def test_the_terms_block_carries_the_disclaimer_items(self):
        t = json.loads(self.counts_path.read_text(encoding="utf-8"))["registry_terms"]
        self.assertEqual(t["no_warranty"], "ClinicalTrials.gov states that the U.S. Government makes no warranties about its data and "
                                           "assumes no liability for their use.")
        self.assertIn("Study sponsors and investigators write and are responsible for their own records", t["sponsor_responsibility"])
        self.assertIn("See the registry's Disclaimer.", t["sponsor_responsibility"])
        self.assertEqual(t["disclaimer_url"], "https://clinicaltrials.gov/about-site/disclaimer")
        self.assertEqual(t["disclaimer_last_updated"], "2023-08-03")
        self.assertIn("third-party copyright", t["third_party_copyright"])
        self.assertIn("international copyright outside the United States", t["third_party_copyright"])
        self.assertIn("kept current at all times", t["keep_current"])
        text = ca.terms_text(t)
        for key in ("no_warranty", "sponsor_responsibility", "third_party_copyright", "keep_current"):
            self.assertIn(t[key], text)
        self.assertIn("https://clinicaltrials.gov/about-site/disclaimer (last updated 2023-08-03)", text)

    def test_time_zones_are_stated(self):
        c = json.loads(self.counts_path.read_text(encoding="utf-8"))
        self.assertTrue(c["registry_terms"]["snapshot_fetched_at"].endswith("Z"))
        text = ca.terms_text(c["registry_terms"])
        self.assertIn(f"snapshot fetched at: {c['registry_terms']['snapshot_fetched_at']} (UTC)", text)
        self.assertIn(f"data processed by the registry: {c['registry_terms']['data_processed_by_registry']} (as given by the registry)",
                      text)

    def test_counts_say_whether_the_data_are_synthetic(self):
        c = json.loads(self.counts_path.read_text(encoding="utf-8"))
        self.assertIs(c["synthetic"], True)
        clean = self.work / "clean"
        pages = sf.pages(sf.cf_condition_studies(), 8, "condition")
        for p in pages:
            p.pop("_synthetic")
        client = sf.FakeClient({fs.ROUTES["condition"]["query.cond"]: pages},
                               version={"apiVersion": "2.0.0-test", "dataTimestamp": "2026-10-09T09:00:00"})
        fs.run(client, clean, routes=["condition"], cfg={"page_size": 8}, contact=sf.DUMMY_CONTACT)
        s = snap.load(clean)
        self.assertIs(ca.compute_counts(s, lexicon.tag_snapshot(s, LEX), [], LEX)["synthetic"], False)

    def test_class_year_phase(self):
        c = json.loads(self.counts_path.read_text(encoding="utf-8"))["class_year_phase"]
        self.assertEqual(c["exacerbations"], {"actual": {"2005": {"PHASE3": 1}}, "planned": {"2027": {"PHASE3": 1}}})
        self.assertEqual(c["pharmacokinetics"], {"actual": {"2016": {"PHASE1": 1}}, "planned": {}})
        self.assertEqual(c["feasibility_adherence"]["actual"], {"2018": {"NA": 1}})

    def test_co_occurrence(self):
        co = json.loads(self.counts_path.read_text(encoding="utf-8"))["co_occurrence"]
        pairs = {(p["a"], p["b"]): p["studies"] for p in co["pairs"]}
        self.assertEqual(pairs[("fev1", "sweat_chloride")], 1)
        self.assertEqual(pairs[("exacerbations", "healthcare_use")], 1)
        self.assertEqual(pairs[("exacerbations", "cfqr")], 1)
        self.assertEqual(pairs[("npd", "sinus_upper_airway")], 1)
        self.assertFalse([p for p in pairs if {"other", "not_stated"} & set(p)])
        self.assertEqual(co["excluded_classes"], ["other", "not_stated"])
        self.assertEqual(co["studies_by_class"]["exacerbations"], 2)
        order = [(LEX.class_order.index(a), LEX.class_order.index(b)) for a, b in pairs]
        self.assertEqual(order, sorted(order))
        self.assertTrue(all(i < j for i, j in order))

    def test_unsorted_entries(self):
        c = json.loads(self.counts_path.read_text(encoding="utf-8"))
        self.assertEqual(c["unsorted_entries"], 0)                       # the one unclassified entry has an accepted model tag
        s = snap.load(self.snap_dir)
        self.assertEqual(ca.compute_counts(s, lexicon.tag_snapshot(s, LEX), [], LEX)["unsorted_entries"], 1)

    def test_each_new_count_is_drift_checked(self):
        for path, mutate in (("class_year_phase", lambda c: c["class_year_phase"]["exacerbations"]["actual"].update({"2005": {"PHASE3": 2}})),
                             ("co_occurrence", lambda c: c["co_occurrence"]["pairs"][0].update({"studies": 9})),
                             ("unsorted_entries", lambda c: c.update({"unsorted_entries": 5})),
                             ("synthetic", lambda c: c.update({"synthetic": False})),
                             ("registry_terms/no_warranty", lambda c: c["registry_terms"].update({"no_warranty": "x"}))):
            with self.subTest(path=path):
                counts = json.loads(self.counts_path.read_text(encoding="utf-8"))
                mutate(counts)
                p = self.work / f"counts-{path.replace('/', '-')}.json"
                p.write_text(json.dumps(counts), encoding="utf-8")
                code, out = run(["--snapshot", self.snap_dir, "--tags", self.tags_path, "--counts", p, "--model-tags", self.model_path,
                                 "--drift-only"])
                self.assertEqual(code, 1, out)
                self.assertIn(f"drift /{path}", out)

    def test_missing_input_files_are_a_clean_usage_error(self):
        nowhere = self.work / "nowhere.json"
        for flag in ("--exclude", "--explained", "--canary", "--planted", "--control-snapshot"):
            with self.subTest(flag=flag):
                extra = [flag, nowhere] + (["--negative-controls"] if flag in ("--canary", "--planted", "--control-snapshot") else [])
                with contextlib.redirect_stderr(io.StringIO()) as err:
                    code, out = run(self.base() + extra)
                self.assertEqual(code, 2, out)
                self.assertIn("ERROR", out)
                self.assertNotIn("Traceback", out + err.getvalue())
        code, out = run(["--negative-controls", "--control-snapshot", nowhere])
        self.assertEqual(code, 2, out)
        self.assertIn("ERROR", out)

    def test_confirm_real_run_refuses_a_synthetic_snapshot(self):
        copy = self.work / "snap"
        shutil.copytree(self.snap_dir, copy)
        code, out = run(["--snapshot", copy, "--confirm-real-run"])
        self.assertEqual(code, 2, out)
        self.assertIn("synthetic", out)
        self.assertFalse(json.loads((copy / "manifest.json").read_text(encoding="utf-8"))["config_verified"])

    def test_confirm_real_run_flips_config_verified(self):
        studies = sf.cf_condition_studies()
        pages = sf.pages(studies, 8, "condition")
        for p in pages:
            p.pop("_synthetic")
        version = {"apiVersion": "2.0.0-test", "dataTimestamp": "2026-10-09T09:00:00"}
        client = sf.FakeClient({fs.ROUTES["condition"]["query.cond"]: pages}, version=version)
        fs.run(client, self.work / "real", routes=["condition"], cfg={"page_size": 8}, contact=sf.DUMMY_CONTACT)
        before = snap.load(self.work / "real").digest
        code, out = run(["--snapshot", self.work / "real", "--confirm-real-run"])
        self.assertEqual(code, 0, out)
        s = snap.load(self.work / "real")
        self.assertTrue(s.manifest["config_verified"])
        self.assertTrue(s.manifest["config_verified_at"])
        self.assertNotEqual(s.digest, before)
        self.assertIn("re-run", out)
        code, out = run(["--snapshot", self.work / "real", "--confirm-real-run"])
        self.assertEqual(code, 2, out)                                    # already confirmed


class ModelTagTest(Fixture):
    def entries_and_rules(self):
        entries = ca.entry_texts(self.s)
        rules = {e["entry_id"]: e for e in json.loads(self.tags_path.read_text(encoding="utf-8"))["entries"]}
        return entries, rules

    def test_canaries_are_all_rejected_and_a_true_quote_is_accepted(self):
        entries, rules = self.entries_and_rules()
        acc, rej = ca.verify_model_tags(sf.canary_model_tags(), entries, rules, LEX)
        self.assertEqual(acc, [])
        self.assertEqual(len(rej), len(sf.canary_model_tags()))
        reasons = [r["reason"] for r in rej]
        self.assertIn("the quote matches this entry only if case is ignored", reasons)
        acc, rej = ca.verify_model_tags(sf.good_model_tags(), entries, rules, LEX)
        self.assertEqual((len(acc), rej), (1, []))
        self.assertEqual(entries["NCT00000014:P3"]["measure"][acc[0]["start"]:acc[0]["end"]], "return to school")

    def test_a_model_may_not_retag_a_rule_decision_or_use_not_stated(self):
        entries, rules = self.entries_and_rules()
        acc, rej = ca.verify_model_tags([
            {"entry_id": "NCT00000001:P1", "class": "spirometry_other", "quote": "FEV1"},
            {"entry_id": "NCT00000014:P3", "class": "not_stated", "quote": "return to school"},
            {"entry_id": "NCT00000014:P3", "class": "made_up", "quote": "return to school"},
            {"entry_id": "NCT00000014:P3", "class": "other", "quote": "Ti"},
            {"entry_id": "NCT00000014:P3", "class": "other"},
            "not an object"], entries, rules, LEX)
        self.assertEqual(acc, [])
        self.assertEqual(len(rej), 6)

    def test_verify_cli_prints_one_line_per_tag(self):
        p = self.work / "canary.json"
        p.write_text(json.dumps(sf.canary_model_tags() + sf.stray_model_tags() + sf.good_model_tags()), encoding="utf-8")
        code, out = run(["--snapshot", self.snap_dir, "--tags", self.tags_path, "--verify-model-tags", p])
        self.assertEqual(code, 0)
        self.assertIn("1 accepted, 5 rejected", out)        # the 4 canaries and the one whose entry is not in the snapshot

    def test_rejection_rate_is_printed(self):
        p = self.work / "mixed.json"
        p.write_text(json.dumps(sf.good_model_tags() + sf.canary_model_tags()[:1]), encoding="utf-8")
        code, out = run(["--snapshot", self.snap_dir, "--tags", self.tags_path, "--counts", self.counts_path, "--model-tags", p,
                         "--drift-only"])
        self.assertEqual(code, 0, out)
        self.assertIn("2 proposed, 1 accepted, 1 rejected (quote-rejection rate 50.0%)", out)


class MetricsTest(unittest.TestCase):
    def test_wilson(self):
        lo, hi = ca.wilson(9, 10)
        self.assertAlmostEqual(lo, 0.5958, places=4)
        self.assertAlmostEqual(hi, 0.9821, places=4)
        self.assertIsNone(ca.wilson(0, 0))
        self.assertEqual(ca.wilson(0, 5)[0], 0.0)

    def test_precision_recall_and_f1(self):
        gold = {f"e{i}": {"fev1"} for i in range(10)} | {"x1": {"lci"}, "x2": {"lci"}}
        pred = {f"e{i}": {"fev1"} for i in range(8)} | {"e8": set(), "e9": set(), "x1": {"fev1", "lci"}, "x2": set()}
        m = ca.metrics(gold, pred, ["fev1", "lci", "other"], 10)
        f = m["per_class"]["fev1"]
        self.assertEqual((f["tp"], f["fp"], f["fn"]), (8, 1, 2))
        self.assertAlmostEqual(f["precision"], 8 / 9)
        self.assertAlmostEqual(f["recall"], 0.8)
        self.assertFalse(f["precision_estimable"])               # 9 tagged: under 10
        self.assertTrue(f["recall_estimable"])                   # 10 labelled
        self.assertFalse(m["per_class"]["lci"]["recall_estimable"])
        self.assertEqual(m["rows_checked"], 12)
        self.assertEqual(m["macro_over_classes"], 2)
        self.assertAlmostEqual(m["micro_f1"], 2 * (9 / 10) * (9 / 12) / ((9 / 10) + (9 / 12)))


class PublicationTest(Fixture):
    def test_without_a_frozen_set_publication_is_blocked(self):
        code, out = run(self.base())
        self.assertEqual(code, 3, out)
        self.assertIn("STOP S1 not evaluated: no frozen set", out)
        self.assertIn("STOP S4 not evaluated", out)

    def test_unexplained_route_difference_blocks(self):
        code, out = run(self.base() + self.frozen() + ["--negative-controls"])
        self.assertEqual(code, 3, out)
        self.assertIn("STOP S3 the retrieval routes differ by 0.250 unexplained (3 studies)", out)

    def test_everything_in_place_passes_and_few_checked_is_said(self):
        # The synthetic inputs can pass only through the test-only function argument, which no command line can set.
        with mock.patch.object(ca, "MIN_FROZEN_ROWS", 19), mock.patch.object(ca, "MIN_CONTROL_ENTRIES", 8):
            code, out = run(self.base() + self.frozen() + self.controls() + ["--explained", self.explained],
                            allow_synthetic_for_tests=True)
        self.assertEqual(code, 0, out)
        self.assertIn("too few to estimate", out)
        self.assertIn("few checked", out)
        self.assertIn("RESULT: PASS (integrity, controls and every publication stop rule)", out)

    def test_low_precision_blocks_where_estimable(self):
        def wrong(eid, cls):
            return ["spirometry_other"] if "fev1" in cls else cls
        with mock.patch.dict(ca.THRESHOLDS, {"min_positives": 1}):
            code, out = run(self.base() + self.frozen(wrong) + ["--negative-controls", "--explained", self.explained])
        self.assertEqual(code, 3, out)
        self.assertIn("STOP S1 fev1: precision 0.000 is below 0.85", out)
        self.assertNotIn("STOP S1 spirometry_other", out)        # no study is shown under it, so it cannot block

    def test_too_much_other_blocks(self):
        with mock.patch.dict(ca.THRESHOLDS, {"other_share": 0.01}):
            code, out = run(self.base() + self.frozen() + ["--negative-controls", "--explained", self.explained])
        self.assertEqual(code, 3, out)
        self.assertIn("STOP S2 other plus unclassified are 5.0% of entries", out)

    # ---- fix round after e7ea5ed: findings 5, 6, 12, 13 and 14. Each case failed before its fix.

    def test_synthetic_controls_do_not_meet_the_gate(self):
        with mock.patch.object(ca, "MIN_FROZEN_ROWS", 19):
            code, out = run(self.base() + self.frozen() + ["--negative-controls", "--explained", self.explained])
        self.assertEqual(code, 3, out)
        self.assertIn("STOP S4 controls ran on synthetic fixtures only", out)

    def test_fewer_than_the_design_number_of_labelled_rows_blocks(self):
        self.assertEqual(ca.MIN_FROZEN_ROWS, 50)
        code, out = run(self.base() + self.frozen() + self.controls() + ["--explained", self.explained])
        self.assertEqual(code, 3, out)
        self.assertIn("STOP S1 only 19 labelled rows; the design needs at least 50", out)

    def test_an_unlabelled_row_blocks(self):
        first = []

        def skip_one(eid, cls):
            if not first:
                first.append(eid)
                return []
            return cls
        with mock.patch.object(ca, "MIN_FROZEN_ROWS", 18):
            code, out = run(self.base() + self.frozen(skip_one) + self.controls() + ["--explained", self.explained])
        self.assertEqual(code, 3, out)
        self.assertIn("STOP S1 1 row(s) of the frozen sheet are unlabelled", out)

    def test_a_shown_class_with_no_labelled_positive_blocks(self):
        def no_pk(eid, cls):
            return ["other"] if "pharmacokinetics" in cls else cls
        with mock.patch.object(ca, "MIN_FROZEN_ROWS", 19):
            code, out = run(self.base() + self.frozen(no_pk) + self.controls() + ["--explained", self.explained])
        self.assertEqual(code, 3, out)
        self.assertIn("STOP S1 shown class(es) with no labelled positive: pharmacokinetics", out)

    def test_a_bad_model_tags_file_exits_2(self):
        p = self.work / "bad.json"
        p.write_text("{not json", encoding="utf-8")
        code, out = run(["--snapshot", self.snap_dir, "--tags", self.tags_path, "--counts", self.counts_path, "--model-tags", p])
        self.assertEqual(code, 2, out)
        self.assertIn("ERROR: cannot read the proposed model tags", out)

    def test_a_bad_tags_file_exits_2(self):
        for content in ("[]", json.dumps({"entries": [{"entry_id": 1}]}), "{not json"):
            p = self.work / "tags-bad.json"
            p.write_text(content, encoding="utf-8")
            with self.subTest(content=content):
                code, out = run(["--snapshot", self.snap_dir, "--tags", p, "--counts", self.counts_path])
                self.assertEqual(code, 2, out)
                self.assertIn("ERROR:", out)

    def test_model_output_is_sanitised_when_printed(self):
        p = self.work / "evil.json"
        p.write_text(json.dumps([{"entry_id": "NCT\x1b[31mX", "class": "\x1b]0;title\x07other", "quote": "q\x1b[2Jq"}]),
                     encoding="utf-8")
        for args in (self.base()[:-2] + ["--model-tags", p, "--drift-only"],
                     ["--snapshot", self.snap_dir, "--tags", self.tags_path, "--verify-model-tags", p]):
            with self.subTest(args=args[-2:]):
                _, out = run(args)
                self.assertNotIn("\x1b", out)
                self.assertNotIn("\x07", out)

    def test_a_missing_sponsor_is_not_counted_as_a_sponsor(self):
        studies = sf.cf_condition_studies()
        for st in studies:
            if st["protocolSection"]["identificationModule"]["nctId"] == "NCT00000015":
                del st["protocolSection"]["sponsorCollaboratorsModule"]
        client = sf.FakeClient({fs.ROUTES["condition"]["query.cond"]: sf.pages(studies, 8, "condition")})
        fs.run(client, self.work / "nosponsor", routes=["condition"], cfg={"page_size": 8}, contact=sf.DUMMY_CONTACT)
        s = snap.load(self.work / "nosponsor")
        counts = ca.compute_counts(s, lexicon.tag_snapshot(s, LEX), [], LEX)
        self.assertEqual(counts["lead_sponsors_by_class"]["not_stated"], 1)          # NCT00000008's sponsor only
        self.assertEqual(counts["studies_without_lead_sponsor_by_class"]["not_stated"], 1)

    def test_an_edited_frozen_sheet_is_refused(self):
        args = self.frozen()
        sheet = Path(args[1])
        sheet.write_text(sheet.read_text(encoding="utf-8").replace("sweat chloride", "sweat sodium"), encoding="utf-8")
        code, out = run(self.base() + args)
        self.assertEqual(code, 2, out)
        self.assertIn("does not match the sealed key", out)

    # ---- second fix round (after 0401af9). Each case failed before its fix.

    def test_synthetic_inputs_never_pass_the_gate(self):
        with mock.patch.object(ca, "MIN_FROZEN_ROWS", 19), mock.patch.object(ca, "MIN_CONTROL_ENTRIES", 8):
            code, out = run(self.base() + self.frozen() + self.controls() + ["--explained", self.explained])
        self.assertEqual(code, 3, out)
        self.assertIn("STOP S0 synthetic data", out)

    def test_the_reviewers_gate_probe_construction_blocks(self):
        """gate_probe.py: 90 synthetic studies, a 50-row frozen set labelled with the rule tags, named controls. It printed PASS."""
        base = self.work
        words = [("Absolute change in ppFEV1", "Week 24"), ("Change in sweat chloride", "Week 4"),
                 ("Number of participants with adverse events", "Up to 28 days")]
        studies = [sf.study(f"NCT0000{2000 + i:04d}", f"s{i}", outcomes=[(words[i % 3][0], "", words[i % 3][1])],
                            start=f"{2001 + i % 20}-01", phases=("PHASE2", "PHASE3")[i % 2:i % 2 + 1]) for i in range(90)]
        studies += [sf.study("NCT00009001", "planted", conditions=("Asthma",), outcomes=[("Change in FEV1", "", "12 weeks")]),
                    sf.study("NCT00009002", "planted2", conditions=("COPD",), outcomes=[("Change in FEV1", "", "12 weeks")])]
        r = fs.ROUTES
        client = sf.FakeClient({r["condition"]["query.cond"]: sf.pages(studies, 50, "condition"),
                                r["term"]["query.term"]: sf.pages(studies, 50, "term")})
        fs.run(client, base / "cf", routes=["condition", "term"], cfg={"page_size": 50}, project_ua="SYNTHETIC", contact=sf.DUMMY_CONTACT)
        s = snap.load(base / "cf")
        tags = lexicon.tag_snapshot(s, LEX)
        (base / "tags.json").write_text(json.dumps(tags), encoding="utf-8")
        (base / "planted.json").write_text(json.dumps(["NCT00009001", "NCT00009002"]), encoding="utf-8")
        sheet, key, _ = mls.build(s, tags, LEX, n=50, seed=1)
        (base / "keydir").mkdir()
        (base / "keydir" / "key.json").write_text(json.dumps(key), encoding="utf-8")
        byid = {e["entry_id"]: e for e in tags["entries"]}
        samp = {row["sample_id"]: row["entry_id"] for row in key["rows"]}
        rows = list(csv.DictReader(io.StringIO(sheet)))
        for row in rows:
            row["label_classes"] = ";".join(t["class"] for t in byid[samp[row["sample_id"]]]["tags"])
        with (base / "sheet.csv").open("w", encoding="utf-8-sig", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=mls.COLUMNS, lineterminator="\n")
            w.writeheader()
            w.writerows(rows)
        args = ["--snapshot", base / "cf", "--tags", base / "tags.json"]
        self.assertEqual(run(args + ["--write-counts", base / "counts.json"])[0], 0)
        code, out = run(args + ["--counts", base / "counts.json", "--frozen-sheet", base / "sheet.csv", "--frozen-key",
                                base / "keydir" / "key.json", "--negative-controls", "--canary", self.canary_path, "--planted",
                                base / "planted.json", "--control-snapshot", self.control_dir])
        # It no longer passes: its canaries name entries that are not in its snapshot, which is now an error (exit 1). With valid
        # canaries it would still stop at S0, because every page carries the synthetic marker (test_synthetic_inputs_never_pass_the_gate).
        self.assertEqual(code, 1, out)
        self.assertIn("ERROR, not a valid canary", out)
        self.assertNotIn("RESULT: PASS", out)
        with mock.patch.object(ca, "run_negative_controls", lambda *a, **k: [("stub control", True, "stubbed")]):
            code, out = run(args + ["--counts", base / "counts.json", "--frozen-sheet", base / "sheet.csv", "--frozen-key",
                                    base / "keydir" / "key.json", "--negative-controls", "--canary", self.canary_path, "--planted",
                                    base / "planted.json", "--control-snapshot", self.control_dir])
        self.assertEqual(code, 3, out)
        self.assertIn("STOP S0 synthetic data: the snapshot carries the synthetic marker", out)
        self.assertNotIn("RESULT: PASS", out)

    def test_a_small_control_snapshot_blocks(self):
        self.assertEqual(ca.MIN_CONTROL_ENTRIES, 50)
        with mock.patch.object(ca, "MIN_FROZEN_ROWS", 19):
            code, out = run(self.base() + self.frozen() + self.controls() + ["--explained", self.explained],
                            allow_synthetic_for_tests=True)
        self.assertEqual(code, 3, out)
        self.assertIn("STOP S4 the control snapshot has 8 outcome entries; at least 50 are needed", out)

    def test_a_frozen_set_from_another_snapshot_blocks(self):
        frozen = self.frozen()                                   # drawn from snapshot A (the class fixture)
        studies = sf.cf_condition_studies()
        studies[0]["protocolSection"]["outcomesModule"]["primaryOutcomes"][1]["measure"] = "Change from baseline in sweat chloride (mmol/L)"
        client = sf.FakeClient({fs.ROUTES["condition"]["query.cond"]: sf.pages(studies, 8, "condition"),
                                fs.ROUTES["term"]["query.term"]: sf.pages(sf.cf_term_studies(), 8, "term")})
        fs.run(client, self.work / "B", routes=["condition", "term"], cfg={"page_size": 8}, contact=sf.DUMMY_CONTACT)
        sb = snap.load(self.work / "B")
        (self.work / "tagsB.json").write_text(json.dumps(lexicon.tag_snapshot(sb, LEX)), encoding="utf-8")
        args = ["--snapshot", self.work / "B", "--tags", self.work / "tagsB.json", "--model-tags", self.model_path]
        self.assertEqual(run(args + ["--write-counts", self.work / "countsB.json"])[0], 0)
        with mock.patch.object(ca, "MIN_FROZEN_ROWS", 19), mock.patch.object(ca, "MIN_CONTROL_ENTRIES", 8):
            code, out = run(args + ["--counts", self.work / "countsB.json"] + frozen + self.controls() + ["--explained", self.explained],
                            allow_synthetic_for_tests=True)
        self.assertEqual(code, 3, out)
        self.assertIn("STOP S1 the frozen set does not belong to this snapshot and lexicon", out)

    def test_a_frozen_key_from_another_lexicon_blocks(self):
        frozen = self.frozen()
        key = json.loads(Path(frozen[3]).read_text(encoding="utf-8"))
        key["lexicon_sha256"] = "0" * 64
        Path(frozen[3]).write_text(json.dumps(key), encoding="utf-8")
        with mock.patch.object(ca, "MIN_FROZEN_ROWS", 19), mock.patch.object(ca, "MIN_CONTROL_ENTRIES", 8):
            code, out = run(self.base() + frozen + self.controls() + ["--explained", self.explained], allow_synthetic_for_tests=True)
        self.assertEqual(code, 3, out)
        self.assertIn("the key was drawn with a different lexicon", out)

    def test_an_unknown_label_is_refused(self):
        code, out = run(self.base() + self.frozen(lambda e, c: ["not_a_class"]))
        self.assertEqual(code, 2, out)
        self.assertIn("unknown class", out)


class NegativeControlTest(unittest.TestCase):
    def setUp(self):
        sf.block_network(self)

    def test_controls_pass_on_the_synthetic_fixtures(self):
        code, out = run(["--negative-controls", "--controls-only"])
        self.assertEqual(code, 0, out)
        self.assertEqual(out.count("PASS control"), 3)
        self.assertIn("4 of 4 rejected for their quote", out)
        self.assertIn("2 of 2 planted ids excluded", out)
        self.assertIn("0 of 8 entries", out)

    def test_a_verifier_that_accepts_everything_fails_the_control(self):
        def accept_all(proposed, entries, rules, lex):
            return [dict(t, text=t["quote"]) for t in proposed], []
        with mock.patch.object(ca, "verify_model_tags", accept_all):
            code, out = run(["--negative-controls"])
        self.assertEqual(code, 1, out)
        self.assertIn("FAIL control: canary", out)

    def test_a_scope_that_lets_non_cf_in_fails_the_control(self):
        with mock.patch.object(ca.scope_mod, "condition_names_cf", lambda c: True):
            code, out = run(["--negative-controls"])
        self.assertEqual(code, 1, out)
        self.assertIn("FAIL control: planted", out)

    def test_a_tagger_that_finds_cf_classes_in_a_non_cf_snapshot_fails_the_control(self):
        with tempfile.TemporaryDirectory() as tmp:
            studies = sf.noncf_studies() + [sf.study("NCT00000199", "control with sweat wording", conditions=("Asthma",),
                                                     outcomes=[("Change in sweat chloride", "", "Week 4")])]
            client = sf.FakeClient({sf.NONCF_ROUTE["condition"]["query.cond"]: sf.pages(studies, 8, "condition")})
            fs.run(client, Path(tmp) / "n", routes=["condition"], cfg={"page_size": 8}, contact=sf.DUMMY_CONTACT,
                   route_params=sf.NONCF_ROUTE)
            code, out = run(["--negative-controls", "--control-snapshot", Path(tmp) / "n"])
        self.assertEqual(code, 1, out)
        self.assertIn("FAIL control: non-CF control", out)
        self.assertIn("NCT00000199:P1", out)

    def test_the_control_covers_exacerbations(self):
        with tempfile.TemporaryDirectory() as tmp:
            studies = sf.noncf_studies() + [sf.study("NCT00000198", "control with pulmonary exacerbation wording", conditions=("COPD",),
                                                     outcomes=[("Rate of pulmonary exacerbations", "", "52 weeks")])]
            client = sf.FakeClient({sf.NONCF_ROUTE["condition"]["query.cond"]: sf.pages(studies, 8, "condition")})
            fs.run(client, Path(tmp) / "n", routes=["condition"], cfg={"page_size": 8}, contact=sf.DUMMY_CONTACT,
                   route_params=sf.NONCF_ROUTE)
            code, out = run(["--negative-controls", "--control-snapshot", Path(tmp) / "n"])
        self.assertEqual(code, 1, out)
        self.assertIn("NCT00000198:P1", out)

    # ---- second fix round (after 0401af9). Each case failed before its fix.

    def test_controls_only_is_not_the_gate(self):
        code, out = run(["--negative-controls"])
        self.assertEqual(code, 4, out)
        self.assertEqual(out.strip().splitlines()[-1], "NOT THE GATE: controls only")

    def test_a_canary_whose_entry_is_not_in_the_snapshot_is_an_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            sf.build_cf_snapshot(Path(tmp) / "s")
            results = ca.run_negative_controls(LEX, snapshot=snap.load(Path(tmp) / "s"),
                                               canary=sf.canary_model_tags() + sf.stray_model_tags(), planted=sf.PLANTED_NON_CF_IDS)
        canary = next(r for r in results if r[0].startswith("canary"))
        self.assertFalse(canary[1])
        self.assertIn("NCT00000099:P1", canary[2])

    def test_a_canary_on_an_entry_the_rules_decided_is_an_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            sf.build_cf_snapshot(Path(tmp) / "s")
            results = ca.run_negative_controls(LEX, snapshot=snap.load(Path(tmp) / "s"),
                                               canary=[{"entry_id": "NCT00000001:P1", "class": "other", "quote": "SYNTHETIC false quote"}],
                                               planted=sf.PLANTED_NON_CF_IDS)
        self.assertFalse(next(r for r in results if r[0].startswith("canary"))[1])

    def test_planted_ids_found_only_by_the_recall_route_are_checked(self):
        # third round: a planted id that only the term search returned is still tested against the scope rules
        cond = [st for st in sf.cf_condition_studies()
                if st["protocolSection"]["identificationModule"]["nctId"] not in ("NCT00000011", "NCT00000016")]
        term = sf.cf_term_studies() + [st for st in sf.cf_condition_studies()
                                       if st["protocolSection"]["identificationModule"]["nctId"] == "NCT00000011"]
        with tempfile.TemporaryDirectory() as tmp:
            client = sf.FakeClient({fs.ROUTES["condition"]["query.cond"]: sf.pages(cond, 8, "condition"),
                                    fs.ROUTES["term"]["query.term"]: sf.pages(term, 8, "term")})
            fs.run(client, Path(tmp) / "s", routes=["condition", "term"], cfg={"page_size": 8}, contact=sf.DUMMY_CONTACT)
            s = snap.load(Path(tmp) / "s")
            ok_case = ca.run_negative_controls(LEX, snapshot=s, canary=sf.canary_model_tags(), planted=["NCT00000011"])
            leak_case = ca.run_negative_controls(LEX, snapshot=s, canary=sf.canary_model_tags(), planted=["NCT00000017"])
        planted = next(r for r in ok_case if r[0].startswith("planted"))
        self.assertTrue(planted[1], planted[2])
        self.assertIn("1 of 1 planted ids excluded", planted[2])
        leaked = next(r for r in leak_case if r[0].startswith("planted"))
        self.assertFalse(leaked[1])
        self.assertIn("INCLUDED", leaked[2])

    def test_a_planted_id_missing_from_the_snapshot_is_not_a_pass(self):
        lex = LEX
        with tempfile.TemporaryDirectory() as tmp:
            sf.build_cf_snapshot(Path(tmp) / "s")
            results = ca.run_negative_controls(lex, snapshot=snap.load(Path(tmp) / "s"), canary=sf.canary_model_tags(),
                                               planted=["NCT00000098"])
        planted = next(r for r in results if r[0].startswith("planted"))
        self.assertFalse(planted[1])
        self.assertIn("not in the snapshot", planted[2])


if __name__ == "__main__":
    unittest.main()
