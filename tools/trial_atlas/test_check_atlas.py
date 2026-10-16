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

    # ---- the terms block, split counts, input errors, confirming a real run

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
        self.assertEqual(c["studies_by_class_and_actual_start_year"]["exacerbations"], {})
        self.assertEqual(c["studies_by_class_and_planned_start_year"]["exacerbations"], {"2027": 1})
        self.assertEqual(c["studies_by_class_and_untyped_start_year"]["exacerbations"], {"2005": 1})   # NCT00000012: no date type
        self.assertEqual(c["studies_by_class_and_first_posted_year"]["exacerbations"], {"2019": 1, "2026": 1})

    def test_untyped_starts_are_counted_and_the_kinds_add_up(self):
        c = json.loads(self.counts_path.read_text(encoding="utf-8"))
        kinds = {k: c["scope"][k] for k in ("start_actual", "start_planned", "start_untyped", "start_no_date")}
        self.assertEqual(kinds, {"start_actual": 7, "start_planned": 1, "start_untyped": 3, "start_no_date": 1})
        self.assertEqual(sum(kinds.values()), c["studies"])
        self.assertNotIn("start_unknown", c["scope"])
        # every study with a start year is in exactly one of the three year-split counts; a study with no date is in none
        dated = sum(c["studies_by_start_year"].get(y, 0) for y in c["studies_by_start_year"] if y != "unknown")
        self.assertEqual(dated, kinds["start_actual"] + kinds["start_planned"] + kinds["start_untyped"])
        self.assertEqual(c["studies_by_start_year"].get("unknown"), kinds["start_no_date"])
        mods = " ".join(c["registry_terms"]["modifications"])
        self.assertIn("Start dates without a recorded date type are counted as their own kind; no type is assumed.", mods)

    def test_the_new_start_kind_keys_are_drift_checked(self):
        for path, mutate in (("scope/start_untyped", lambda c: c["scope"].update({"start_untyped": 4})),
                             ("scope/start_no_date", lambda c: c["scope"].update({"start_no_date": 0})),
                             ("studies_by_class_and_untyped_start_year", lambda c: c["studies_by_class_and_untyped_start_year"].update(
                                 {"exacerbations": {"2005": 2}})),
                             ("class_year_phase", lambda c: c["class_year_phase"]["exacerbations"]["untyped"].update({"2005": {"PHASE3": 2}}))):
            with self.subTest(path=path):
                counts = json.loads(self.counts_path.read_text(encoding="utf-8"))
                mutate(counts)
                p = self.work / f"counts-kind-{path.replace('/', '-')}.json"
                p.write_text(json.dumps(counts), encoding="utf-8")
                code, out = run(["--snapshot", self.snap_dir, "--tags", self.tags_path, "--counts", p, "--model-tags", self.model_path,
                                 "--drift-only"])
                self.assertEqual(code, 1, out)
                self.assertIn(f"drift /{path}", out)

    # ---- the disclaimer items, time zones, the synthetic key and the page counts

    def test_the_terms_block_carries_the_disclaimer_items(self):
        t = json.loads(self.counts_path.read_text(encoding="utf-8"))["registry_terms"]
        self.assertEqual(t["no_warranty"], "ClinicalTrials.gov states that the U.S. Government makes no warranties, expressed or "
                                           "implied, about its data and assumes no liability for any party's use of them.")
        self.assertEqual(t["sponsor_responsibility"],
                         "Study sponsors and investigators write and are responsible for their own records. The registry's Disclaimer "
                         "says the U.S. government \"does not review or approve the safety and science of all studies listed on this "
                         "website\" and that NLM staff only review study information for apparent errors, deficiencies or "
                         "inconsistencies. See the registry's Disclaimer.")
        self.assertEqual(t["disclaimer_url"], "https://clinicaltrials.gov/about-site/disclaimer")
        self.assertEqual(t["disclaimer_last_updated"], "2023-08-03")
        self.assertEqual(t["third_party_copyright"],
                         "Some registry data may be subject to third-party copyright, and the data carry an international copyright "
                         "outside the United States and its Territories or Possessions.")
        self.assertEqual(t["keep_current"],
                         "The registry says it is updated daily and that data in any publication or distribution should be kept current "
                         "at all times. This copy is dated and may be out of date; the live record is the current one.")
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

    def test_utc_is_said_only_when_the_value_says_so(self):
        base = json.loads(self.counts_path.read_text(encoding="utf-8"))["registry_terms"]
        cases = {"2026-10-10T01:02:03Z": "snapshot fetched at: 2026-10-10T01:02:03Z (UTC)",
                 "2026-10-10T01:02:03+02:00": "snapshot fetched at: 2026-10-10T01:02:03+02:00 (time zone not stated as UTC)",
                 "2026-10-10T01:02:03": "snapshot fetched at: 2026-10-10T01:02:03 (time zone not stated as UTC)",
                 None: "snapshot fetched at: not recorded\n"}
        for value, line in cases.items():
            with self.subTest(value=value):
                self.assertIn(line, ca.terms_text({**base, "snapshot_fetched_at": value}))
        text = ca.terms_text({**base, "data_processed_by_registry": None})
        self.assertIn("data processed by the registry: not recorded\n", text)
        self.assertNotIn("not recorded (as given by the registry)", text)

    def test_registry_text_is_cleaned_before_printing(self):
        base = json.loads(self.counts_path.read_text(encoding="utf-8"))["registry_terms"]
        text = ca.terms_text({**base, "data_processed_by_registry": "2026\x1b[2J\nINJECTED", "terms_last_updated": "x\x1b]0;t\x07y"})
        self.assertNotIn("\x1b", text)
        self.assertNotIn("\x07", text)
        self.assertFalse([ln for ln in text.splitlines() if ln.startswith("INJECTED")], "a newline in a value must not start a line")
        counts = json.loads(self.counts_path.read_text(encoding="utf-8"))
        counts["class_year_phase"]["exacerbations"]["untyped"]["2005"]["PHASE3\x1b[31m\nEVIL"] = 1
        p = self.work / "counts-escape.json"
        p.write_text(json.dumps(counts), encoding="utf-8")
        code, out = run(["--snapshot", self.snap_dir, "--tags", self.tags_path, "--counts", p, "--model-tags", self.model_path,
                         "--drift-only"])
        self.assertEqual(code, 1, out)
        self.assertNotIn("\x1b", out)
        self.assertFalse([ln for ln in out.splitlines() if ln.startswith("EVIL")])

    def test_drift_compares_types_not_only_values(self):
        for name, mutate in (("synthetic", lambda c: c.update({"synthetic": 1})),
                             ("unsorted_entries", lambda c: c.update({"unsorted_entries": False})),
                             ("studies", lambda c: c.update({"studies": float(c["studies"])}))):
            with self.subTest(name=name):
                counts = json.loads(self.counts_path.read_text(encoding="utf-8"))
                mutate(counts)
                p = self.work / f"counts-type-{name}.json"
                p.write_text(json.dumps(counts), encoding="utf-8")
                code, out = run(["--snapshot", self.snap_dir, "--tags", self.tags_path, "--counts", p, "--model-tags", self.model_path,
                                 "--drift-only"])
                self.assertEqual(code, 1, out)
                self.assertIn(f"drift /{name}", out)

    def test_the_readme_does_not_claim_the_synthetic_key_proves_origin(self):
        readme = (HERE / "README.md").read_text(encoding="utf-8")
        source = (HERE / "check_atlas.py").read_text(encoding="utf-8")
        for text in (readme, source):
            self.assertNotIn("without guessing", text)
        self.assertIn("not proof of origin", readme)
        self.assertIn("no synthetic marker found", readme)

    def test_co_occurrence_counts_studies_once_per_pair(self):
        studies = [
            sf.study("NCT00000301", "two entries in the same pair", start="2015-01", outcomes=[
                ("Change in FEV1", "", "Week 4"), ("Change in sweat chloride", "", "Week 4"),
                ("Change in FEV1 and sweat chloride", "", "Week 8")]),
            sf.study("NCT00000302", "three classes", start="2016-01", outcomes=[
                ("Change in FEV1", "", "Week 4"), ("Change in sweat chloride", "", "Week 4"), ("Change in HbA1c", "", "Week 4")]),
        ]
        client = sf.FakeClient({fs.ROUTES["condition"]["query.cond"]: sf.pages(studies, 8, "condition")})
        fs.run(client, self.work / "pairs", routes=["condition"], cfg={"page_size": 8}, contact=sf.DUMMY_CONTACT)
        s = snap.load(self.work / "pairs")
        pairs = {(p["a"], p["b"]): p["studies"] for p in ca.compute_counts(s, lexicon.tag_snapshot(s, LEX), [], LEX)["co_occurrence"]["pairs"]}
        self.assertEqual(pairs, {("fev1", "sweat_chloride"): 2, ("fev1", "glucose_cfrd"): 1, ("sweat_chloride", "glucose_cfrd"): 1})

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
        self.assertEqual(c["exacerbations"], {"actual": {}, "planned": {"2027": {"PHASE3": 1}}, "untyped": {"2005": {"PHASE3": 1}}})
        self.assertEqual(c["pharmacokinetics"], {"actual": {"2016": {"PHASE1": 1}}, "planned": {}, "untyped": {}})
        self.assertEqual(c["feasibility_adherence"]["actual"], {"2018": {"NA": 1}})

    def test_co_occurrence(self):
        co = json.loads(self.counts_path.read_text(encoding="utf-8"))["co_occurrence"]
        pairs = {(p["a"], p["b"]): p["studies"] for p in co["pairs"]}
        self.assertEqual(pairs[("fev1", "sweat_chloride")], 1)
        self.assertEqual(pairs[("exacerbations", "healthcare_use")], 1)
        self.assertEqual(pairs[("exacerbations", "cfqr")], 1)
        self.assertEqual(pairs[("npd", "pro_named_other")], 1)              # a named instrument (a fixture; since lexicon 0.3.1-draft SNOT-22 also gets sinus)
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
        for path, mutate in (("class_year_phase", lambda c: c["class_year_phase"]["exacerbations"]["planned"].update({"2027": {"PHASE3": 2}})),
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

    def test_model_sorted_entries_are_in_accuracy_sample(self):
        """Verifies that model-sorted entries sit inside the census-v1 frame and are drawn into the accuracy sample (T-0135)."""
        tags = json.loads(self.tags_path.read_text(encoding="utf-8"))
        good = sf.good_model_tags()
        raw = json.dumps(good).encode("utf-8")
        h = mls.model_tags_digest(raw)
        sheet, key, summary = mls.build(self.s, tags, LEX, design="census-v1", seed=42, model_tags=good, model_tags_sha256=h)
        model_eid = good[0]["entry_id"]
        # Verify model-sorted entry is in the drawn sample
        drawn_ids = {r["entry_id"] for r in key["rows"]}
        self.assertIn(model_eid, drawn_ids)
        # Verify its stratum reflects its model-assigned class
        row = next(r for r in key["rows"] if r["entry_id"] == model_eid)
        self.assertTrue(row["stratum"].startswith("census:small") or row["stratum"].startswith("sample:"))


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

    # ---- controls on named inputs, the frozen-set size, sanitised printing, input errors, sponsors

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

    # ---- the gate refuses synthetic data and unbound frozen sets

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


def tag_entry(eid, classes, status=None, no_primary=False):
    """A hand-made tags-file entry (only the fields the census code reads)."""
    status = status or ("not_stated" if classes == ["not_stated"] else ("rule" if classes else "unclassified"))
    return {"entry_id": eid, "nct_id": eid.split(":")[0], "status": status, "tags": [{"class": c} for c in classes],
            "flags": {"no_primary_outcome": no_primary}}


class CensusUnitTest(unittest.TestCase):
    """The census-v1 design on small hand-made frames and keys (proposals/2026-10-11-frozen-census-design.md)."""

    # ---- weighted metrics

    def frame(self, sampled_weight=2.0):
        gold = {"c1": {"fev1"}, "c2": {"fev1"}, "c3": {"lci"}, "s1": {"fev1"}, "s2": {"other"}}
        pred = {"c1": {"fev1"}, "c2": set(), "c3": {"fev1", "lci"}, "s1": {"fev1"}, "s2": {"fev1", "other"}}
        strata = {"c1": "census:small", "c2": "census:small", "c3": "census:small", "s1": "sample:fev1", "s2": "sample:fev1"}
        weights = {e: (sampled_weight if st.startswith("sample:") else 1.0) for e, st in strata.items()}
        return gold, pred, weights, strata

    def test_weighted_precision_and_recall_by_hand(self):
        gold, pred, w, st = self.frame()
        m = ca.metrics_weighted(gold, pred, ["fev1", "lci", "other"], w, st, {"fev1": 12, "lci": 3, "other": 40}, 10, seed=5)
        f = m["per_class"]["fev1"]
        # census: c1 TP, c2 FN, c3 FP (weight 1); sampled at weight 2: s1 TP, s2 FP
        self.assertEqual((f["tp"], f["fp"], f["fn"]), (3.0, 3.0, 1.0))
        self.assertAlmostEqual(f["precision"], 3 / 6)
        self.assertAlmostEqual(f["recall"], 3 / 4)
        self.assertEqual((f["positives"], f["predicted"]), (4.0, 6.0))
        # 12 entries in the FRAME carry fev1, but only 4 labelled rows are predicted to: under 10, so not estimable
        self.assertFalse(f["precision_estimable"])
        self.assertFalse(f["recall_estimable"])                   # 4 weighted positives from 3 labelled rows: under 10
        self.assertEqual((f["labelled_predicted_rows"], f["labelled_positive_rows"]), (4, 3))
        self.assertEqual(f["frame_tagged"], 12)
        # with a floor of 4 the frame count and the 4 predicted rows both reach it; recall's 3 positive rows do not
        f4 = ca.metrics_weighted(gold, pred, ["fev1"], w, st, {"fev1": 12}, 4, seed=5)["per_class"]["fev1"]
        self.assertEqual((f4["precision_estimable"], f4["recall_estimable"]), (True, False))
        self.assertEqual(m["resamples"], 2000)
        self.assertEqual(m["seed"], 5)
        lo, hi = f["precision_ci"]
        self.assertTrue(0.0 <= lo <= hi <= 1.0)

    def test_census_only_classes_are_exact(self):
        gold, pred, w, st = self.frame()
        m = ca.metrics_weighted(gold, pred, ["fev1", "lci"], w, st, {"fev1": 12, "lci": 3}, 10, seed=5)
        lci = m["per_class"]["lci"]
        self.assertEqual((lci["tp"], lci["fp"], lci["fn"]), (1.0, 0.0, 0.0))
        self.assertEqual((lci["precision"], lci["recall"]), (1.0, 1.0))
        self.assertEqual(lci["precision_ci"], (1.0, 1.0))         # every lci row is in a census stratum: no sampling error
        self.assertEqual(lci["recall_ci"], (1.0, 1.0))
        self.assertFalse(lci["precision_estimable"])              # 3 in the frame
        # a set drawn wholly from census strata has no bootstrap at all and every interval is the point value
        only = {e: g for e, g in gold.items() if e.startswith("c")}
        m2 = ca.metrics_weighted(only, pred, ["fev1"], w, st, {"fev1": 12}, 10, seed=5)
        f = m2["per_class"]["fev1"]
        self.assertEqual(m2["resamples"], 0)
        self.assertEqual((f["precision"], f["recall"]), (0.5, 0.5))
        self.assertEqual((f["precision_ci"], f["recall_ci"]), ((0.5, 0.5), (0.5, 0.5)))

    def test_a_sampled_row_of_weight_2_counts_double(self):
        gold, pred, w, st = self.frame(sampled_weight=1.0)
        one = ca.metrics_weighted(gold, pred, ["fev1"], w, st, {}, 10, seed=1)["per_class"]["fev1"]
        gold, pred, w, st = self.frame(sampled_weight=2.0)
        two = ca.metrics_weighted(gold, pred, ["fev1"], w, st, {}, 10, seed=1)["per_class"]["fev1"]
        self.assertEqual((one["tp"], one["fp"]), (2.0, 2.0))
        self.assertEqual((two["tp"], two["fp"]), (3.0, 3.0))      # the census rows once, s1 and s2 twice each
        self.assertEqual(two["fn"], one["fn"])

    def test_other_is_not_scored(self):
        gold, pred, w, st = self.frame()
        m = ca.metrics_weighted(gold, pred, ["fev1", "lci", "other"], w, st, {"other": 99}, 1, seed=1)
        self.assertNotIn("other", m["per_class"])
        self.assertEqual(m["macro_over_classes"], 2)

    def test_the_bootstrap_is_deterministic_for_a_seed_and_moves_with_it(self):
        gold, pred, w, st = {}, {}, {}, {}
        for i in range(12):
            e = f"s{i}"
            gold[e] = {"fev1"} if i % 3 else {"lci"}
            pred[e] = {"fev1"} if i % 4 else {"lci"}
            st[e], w[e] = ("sample:fev1", 2.5) if i < 8 else ("sample:lci", 1.5)
        gold["c0"], pred["c0"], st["c0"], w["c0"] = {"fev1"}, {"fev1"}, "census:small", 1.0

        def cis(seed):
            m = ca.metrics_weighted(gold, pred, ["fev1", "lci"], w, st, {}, 1, seed=seed, resamples=50)
            return tuple((v["precision_ci"], v["recall_ci"]) for v in m["per_class"].values())
        self.assertEqual(cis(11), cis(11))
        self.assertGreater(len({cis(s) for s in range(1, 6)}), 1)
        full = ca.metrics_weighted(gold, pred, ["fev1"], w, st, {}, 1, seed=11)
        self.assertEqual(full, ca.metrics_weighted(gold, pred, ["fev1"], w, st, {}, 1, seed=11))

    # ---- strata and the binding checks

    LEX_ENTRIES = staticmethod(lambda: [
        tag_entry("N1:f1", ["fev1"]), tag_entry("N1:f2", ["fev1"]), tag_entry("N1:f3", ["fev1"]),
        tag_entry("N2:fs", ["fev1", "sweat_chloride"]), tag_entry("N2:s1", ["sweat_chloride"]), tag_entry("N2:s2", ["sweat_chloride"]),
        tag_entry("N3:l1", ["lci"]), tag_entry("N3:fl", ["fev1", "lci"]),
        tag_entry("N4:u1", []), tag_entry("N4:u2", []),
        tag_entry("N5:n1", ["not_stated"]), tag_entry("N5:n2", ["not_stated"]), tag_entry("N5:n3", ["not_stated"]),
        tag_entry("N6:x0", ["not_stated"], no_primary=True)])
    ACCEPTED = [{"entry_id": "N4:u2", "class": "other"}]
    EXPECTED = {"N1:f1": "sample:fev1", "N1:f2": "sample:fev1", "N1:f3": "sample:fev1",
                "N2:fs": "sample:sweat_chloride", "N2:s1": "sample:sweat_chloride", "N2:s2": "sample:sweat_chloride",
                "N3:l1": "census:small", "N3:fl": "census:small", "N4:u2": "census:small",
                "N4:u1": "census:unclassified", "N5:n1": "census:not_stated", "N5:n2": "census:not_stated", "N5:n3": "census:not_stated"}

    def test_strata_follow_the_design_order(self):
        frame = ca.census_frame(self.LEX_ENTRIES(), self.ACCEPTED)
        self.assertNotIn("N6:x0", frame)                          # no primary-outcome wording: never in the frame
        self.assertEqual(ca.census_strata(frame, LEX, 2), self.EXPECTED)
        # ties between equally rare large classes go by the lexicon's class order (fev1 before lci)
        tie = ca.census_frame([tag_entry("A:1", ["lci", "fev1"]), tag_entry("A:2", ["fev1", "lci"])], [])
        self.assertEqual(ca.census_strata(tie, LEX, 0), {"A:1": "sample:fev1", "A:2": "sample:fev1"})
        # not_stated is never a small class: a lone not_stated entry is still census:not_stated
        self.assertEqual(ca.census_strata(ca.census_frame([tag_entry("B:1", ["not_stated"])], []), LEX, 30),
                         {"B:1": "census:not_stated"})

    def key(self):
        strata = {"sample:fev1": {"kind": "sample", "population": 3, "drawn": 2, "weight": 1.5},
                  "sample:sweat_chloride": {"kind": "sample", "population": 3, "drawn": 2, "weight": 1.5},
                  "census:small": {"kind": "census", "population": 3, "drawn": 3, "weight": 1},
                  "census:unclassified": {"kind": "census", "population": 1, "drawn": 1, "weight": 1},
                  "census:not_stated": {"kind": "census", "population": 3, "drawn": 3, "weight": 1}}
        # each sample stratum holds min(population, large_per_class 2) rows: N1:f3 and N2:s2 were not drawn
        rows = [{"sample_id": f"S{i:03d}", "entry_id": e, "stratum": st}
                for i, (e, st) in enumerate(sorted((e, st) for e, st in self.EXPECTED.items() if e not in ("N1:f3", "N2:s2")), 1)]
        return {"design": "census-v1", "small_max": 2, "large_per_class": 2, "seed": 4, "model_tags_sha256": None,
                "consensus": None, "strata": strata, "rows": rows}

    def problems(self, key, model_path=None):
        return ca.census_binding_problems(key, self.LEX_ENTRIES(), self.ACCEPTED, LEX, model_path)

    def test_a_key_that_describes_the_frame_binds(self):
        self.assertEqual(self.problems(self.key()), [])

    def test_populations_that_do_not_sum_are_a_problem(self):
        k = self.key()
        k["strata"]["census:unclassified"]["population"] = 2
        k["strata"]["census:unclassified"]["drawn"] = 2
        out = " | ".join(self.problems(k))
        self.assertIn("the key's strata hold 14 entries in all; this snapshot's frame has 13", out)

    def test_a_wrong_weight_is_a_problem(self):
        k = self.key()
        k["strata"]["sample:fev1"]["weight"] = 2
        self.assertEqual(self.problems(k), ["1 stratum weight(s) are not population/drawn: sample:fev1 (weight 2, population/drawn 1.5)"])
        k["strata"]["sample:fev1"]["weight"] = 1.5 + 1e-12          # within the tolerance
        self.assertEqual(self.problems(k), [])

    def test_a_row_in_another_stratum_than_the_frame_gives_is_a_problem(self):
        k = self.key()
        for r in k["rows"]:                                      # swap two rows: drawn counts and populations still add up
            if r["entry_id"] == "N1:f1":
                r["stratum"] = "sample:sweat_chloride"
            elif r["entry_id"] == "N2:s1":
                r["stratum"] = "sample:fev1"
        out = self.problems(k)
        self.assertEqual(len(out), 1, out)
        self.assertIn("2 drawn row(s) sit in a different stratum from the one this frame gives", out[0])
        self.assertIn("frame sample:fev1", out[0])

    def test_a_different_small_max_moves_the_frame(self):
        k = self.key()
        k["small_max"] = 3                                       # fev1 still large; sweat_chloride (3) becomes small
        self.assertTrue(any("different stratum" in p for p in self.problems(k)))

    def test_drawn_counts_census_and_names_are_checked(self):
        k = self.key()
        k["strata"]["census:not_stated"]["drawn"] = 2
        k["strata"]["census:not_stated"]["weight"] = 1.5
        out = " | ".join(self.problems(k))
        self.assertIn("record a drawn count that differs from their rows: census:not_stated (drawn 2, rows 3)", out)
        self.assertIn("census stratum(s) were not drawn in full: census:not_stated", out)
        k = self.key()
        k["rows"][0]["stratum"] = "sample:nowhere"
        self.assertTrue(any("name a stratum the key does not hold" in p for p in self.problems(k)))
        k = self.key()
        del k["strata"]
        self.assertEqual(self.problems(k), ["the census key has no strata"])

    def test_a_sample_stratum_must_hold_min_population_and_large_per_class(self):
        self.assertEqual(self.problems(self.key()), [])          # the honest draw: 2 of 3 in each sample stratum, large_per_class 2
        k = self.key()
        k["large_per_class"] = 3                                 # the design now says all 3 of each sample stratum
        self.assertEqual(self.problems(k), [
            "2 sample stratum(s) do not hold min(population, large_per_class 3) rows: sample:fev1 (drawn 2 of 3, the design gives 3), "
            "sample:sweat_chloride (drawn 2 of 3, the design gives 3)"])
        # a hand-edited key: large_per_class 40 but a sample stratum of 300 drawn 3 at weight 100 is refused
        k = self.key()
        k["large_per_class"] = 40
        k["strata"]["sample:fev1"] = {"kind": "sample", "population": 300, "drawn": 3, "weight": 100.0}
        out = " | ".join(self.problems(k))
        self.assertIn("sample stratum(s) do not hold min(population, large_per_class 40) rows: sample:fev1 (drawn 3 of 300, the "
                      "design gives 40)", out)

    def test_estimability_needs_labelled_rows_not_only_weight(self):
        # fev1: 50 entries in the frame, but only 4 labelled rows, each at weight 10 (weighted 40 predicted, 40 positive)
        gold = {f"s{i}": {"fev1"} for i in range(4)}
        pred = {e: {"fev1"} for e in gold}
        w = {e: 10.0 for e in gold}
        st = {e: "sample:fev1" for e in gold}
        f = ca.metrics_weighted(gold, pred, ["fev1"], w, st, {"fev1": 50}, 10, seed=2)["per_class"]["fev1"]
        self.assertEqual((f["predicted"], f["positives"]), (40.0, 40.0))
        self.assertEqual((f["labelled_predicted_rows"], f["labelled_positive_rows"]), (4, 4))
        self.assertEqual((f["precision_estimable"], f["recall_estimable"]), (False, False))
        # ten such rows are enough
        gold = {f"s{i}": {"fev1"} for i in range(10)}
        pred = {e: {"fev1"} for e in gold}
        w, st = {e: 10.0 for e in gold}, {e: "sample:fev1" for e in gold}
        f = ca.metrics_weighted(gold, pred, ["fev1"], w, st, {"fev1": 50}, 10, seed=2)["per_class"]["fev1"]
        self.assertEqual((f["precision_estimable"], f["recall_estimable"]), (True, True))

    def test_the_model_tags_hash_must_match(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / "mt.json"
            p.write_bytes(b'[{"entry_id": "N4:u2", "class": "other", "quote": "x"}]\n')
            self.assertIsNone(ca.model_tags_digest(None))
            digest = ca.model_tags_digest(p)
            self.assertEqual(digest, __import__("hashlib").sha256(p.read_bytes()).hexdigest())
            k = self.key()
            self.assertEqual(self.problems(k, p), ["the key was drawn without model tags, but --model-tags was given"])
            k["model_tags_sha256"] = digest
            self.assertEqual(self.problems(k, p), [])
            self.assertEqual(self.problems(k, None), ["the key was drawn with model tags, but no --model-tags was given"])
            k["model_tags_sha256"] = "0" * 64
            self.assertEqual(self.problems(k, p), ["the key's model_tags_sha256 does not match the model-tags file given"])


    # ---- sample exclusion (development entries kept out of the sample strata)

    def excluded_key(self, dev):
        """The honest key for this frame when the sample strata leave out `dev` (here N1:f3, so sample:fev1 holds N1:f1 and N1:f2 in
        full, and N3:l1, which stays in census:small)."""
        k = self.key()
        k["strata"]["sample:fev1"] = {"kind": "sample", "population": 2, "drawn": 2, "weight": 1.0}
        k["sample_exclusion"] = {"count": len(set(dev) & set(self.EXPECTED)), "sha256": ca.sample_exclusion_digest(dev)}
        return k

    DEV = ["N1:f3", "N3:l1", "ZZ:not-in-the-frame"]

    def test_an_honest_key_with_an_exclusion_binds(self):
        k = self.excluded_key(self.DEV)
        self.assertEqual(k["sample_exclusion"]["count"], 2)
        self.assertEqual(ca.census_binding_problems(k, self.LEX_ENTRIES(), self.ACCEPTED, LEX, None, self.DEV), [])
        # the order and duplicates of the file do not matter; its content does
        self.assertEqual(ca.census_binding_problems(k, self.LEX_ENTRIES(), self.ACCEPTED, LEX, None, self.DEV[::-1] + ["N1:f3"]), [])
        # the census stratum keeps its excluded entry: dropping N3:l1 from census:small is refused
        k2 = self.excluded_key(self.DEV)
        k2["strata"]["census:small"] = {"kind": "census", "population": 2, "drawn": 2, "weight": 1}
        k2["rows"] = [r for r in k2["rows"] if r["entry_id"] != "N3:l1"]
        out = " | ".join(ca.census_binding_problems(k2, self.LEX_ENTRIES(), self.ACCEPTED, LEX, None, self.DEV))
        self.assertIn("stratum population(s) differ from this frame: census:small (key 2, frame 3)", out)
        # the populations sum to the frame less the excluded entries that sit in sample strata (13 - 1)
        k3 = self.excluded_key(self.DEV)
        k3["strata"]["sample:fev1"] = {"kind": "sample", "population": 3, "drawn": 2, "weight": 1.5}
        out = " | ".join(ca.census_binding_problems(k3, self.LEX_ENTRIES(), self.ACCEPTED, LEX, None, self.DEV))
        self.assertIn("the key's strata hold 13 entries in all; this snapshot's frame has 13 in-scope entries with primary-outcome "
                      "wording, 1 of them excluded from the sample strata", out)

    def test_an_excluded_entry_drawn_into_a_sample_stratum_is_refused(self):
        dev = ["N1:f1"]                                         # drawn into sample:fev1 in self.key()
        k = self.key()
        k["strata"]["sample:fev1"] = {"kind": "sample", "population": 2, "drawn": 2, "weight": 1.0}
        k["sample_exclusion"] = {"count": 1, "sha256": ca.sample_exclusion_digest(dev)}
        out = ca.census_binding_problems(k, self.LEX_ENTRIES(), self.ACCEPTED, LEX, None, dev)
        self.assertEqual(out, ["1 row(s) of a sample stratum are excluded (development) entries: S001"])

    def test_a_wrong_or_missing_exclusion_file_is_refused(self):
        k = self.excluded_key(self.DEV)
        p = lambda ids: ca.census_binding_problems(k, self.LEX_ENTRIES(), self.ACCEPTED, LEX, None, ids)  # noqa: E731
        self.assertIn("the key's sampled strata exclude development entries, but no --frozen-exclude file was given", p(None))
        self.assertIn("the --frozen-exclude file's hash differs from the key's sample_exclusion sha256", p(["N1:f3"]))
        k["sample_exclusion"]["count"] = 3
        self.assertIn("the key's sample_exclusion count 3 differs from the 2 excluded id(s) in this frame", p(self.DEV))
        k["sample_exclusion"] = {"count": "2"}
        self.assertIn("the key's sample_exclusion is neither null nor a record with a whole-number count and a sha256", p(self.DEV))

    def test_an_exclusion_file_with_a_null_key_is_refused(self):
        for k in (self.key(), {**self.key(), "sample_exclusion": None}):        # an older key without the field reads as null
            self.assertEqual(ca.census_binding_problems(k, self.LEX_ENTRIES(), self.ACCEPTED, LEX, None, []),
                             ["the key was drawn with no sample exclusion, but --frozen-exclude was given"])
            self.assertEqual(ca.census_binding_problems(k, self.LEX_ENTRIES(), self.ACCEPTED, LEX, None, None), [])

    def test_the_exclusion_file_must_be_a_list_of_strings(self):
        self.assertEqual(ca.sample_exclusion_ids(["a", "b"]), ["a", "b"])
        for bad in ({"ids": []}, [1], "a", None):
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                ca.sample_exclusion_ids(bad)


class CensusGateTest(Fixture):
    """check_atlas.main with a hand-built census-v1 key on the synthetic snapshot."""

    def census(self, small_max=30, large_per_class=40, label=None, mutate=None, seed=3, exclude=None):
        tags = json.loads(self.tags_path.read_text(encoding="utf-8"))
        entries = ca.entry_texts(self.s)
        rules = {e["entry_id"]: e for e in tags["entries"]}
        accepted, _ = ca.verify_model_tags(sf.good_model_tags(), entries, rules, LEX)
        frame = ca.census_frame(tags["entries"], accepted)
        strata_of = ca.census_strata(frame, LEX, small_max)
        members: dict[str, list[str]] = {}
        for e, st in sorted(strata_of.items()):
            members.setdefault(st, []).append(e)
        rng = __import__("random").Random(seed)
        strata, chosen = {}, []
        for st, ids in sorted(members.items()):
            if exclude and st.startswith("sample:"):
                ids = [e for e in ids if e not in exclude]
                if not ids:
                    continue
            take = ids if st.startswith("census:") else rng.sample(ids, min(large_per_class, len(ids)))
            strata[st] = {"kind": "census" if st.startswith("census:") else "sample", "population": len(ids), "drawn": len(take),
                          "weight": len(ids) / len(take)}
            chosen += [(e, st) for e in take]
        rng.shuffle(chosen)
        final = {e: frame[e]["classes"] for e in frame}
        rows, key_rows = [], []
        for i, (e, st) in enumerate(chosen, 1):
            sid = f"S{i:03d}"
            rows.append({"sample_id": sid, "measure": mls.safe_cell(entries[e]["measure"]),
                         "description": mls.safe_cell(entries[e]["description"]), "time_frame": mls.safe_cell(entries[e]["time_frame"]),
                         "label_composite": "", "label_notes": "",
                         "label_classes": ";".join(label(e, final[e]) if label else final[e])})
            key_rows.append({"sample_id": sid, "entry_id": e, "stratum": st})
        key = {"kind": "trial-atlas frozen-set key", "design": "census-v1", "seed": seed, "small_max": small_max,
               "large_per_class": large_per_class, "snapshot_sha256": self.s.digest, "lexicon_sha256": LEX.sha256,
               "tags_sha256": mls.tags_digest(tags), "model_tags_sha256": ca.model_tags_digest(self.model_path), "consensus": None,
               "sheet_wording_sha256": ca.sheet_wording_sha(rows), "strata": strata, "rows": key_rows}
        if exclude is not None:
            key["sample_exclusion"] = {"count": len(set(exclude) & set(frame)), "sha256": ca.sample_exclusion_digest(exclude)}
        if mutate:
            mutate(key)
        sp, kp = self.work / "census-sheet.csv", self.work / "census-key.json"
        with sp.open("w", encoding="utf-8", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=mls.COLUMNS)
            w.writeheader()
            w.writerows(rows)
        kp.write_text(json.dumps(key), encoding="utf-8")
        return ["--frozen-sheet", sp, "--frozen-key", kp], len(rows)

    def gate(self, frozen, rows):
        with mock.patch.object(ca, "MIN_FROZEN_ROWS", rows), mock.patch.object(ca, "MIN_CONTROL_ENTRIES", 8):
            return run(self.base() + frozen + self.controls() + ["--explained", self.explained], allow_synthetic_for_tests=True)

    def test_an_all_census_set_passes_with_exact_values(self):
        frozen, rows = self.census()
        self.assertEqual(rows, 19)                               # the frame: 20 entries less the one with no primary outcome
        code, out = self.gate(frozen, rows)
        self.assertEqual(code, 0, out)
        self.assertEqual(out.count("frozen set design: census-v1 (exact over the census strata, weighted estimate over the sampled "
                                   "ones)"), 1)
        self.assertIn("frozen set strata: census:not_stated 1 of 1 (weight 1), census:small 18 of 18 (weight 1)", out)
        self.assertIn("every stratum drawn in full, so every value is exact", out)
        self.assertIn("tagged in the frame 2", out)              # e.g. pharmacokinetics
        self.assertIn("few checked (2 tagged in the frame, 2 weighted labelled positives in the frozen set; labelled rows: 2 "
                      "predicted, 2 positive)", out)
        self.assertNotRegex(out, r"\n  other ")                  # other is not scored
        self.assertIn("RESULT: PASS", out)

    def test_a_sampled_stratum_is_weighted_through_the_gate(self):
        frozen, rows = self.census(small_max=1, large_per_class=1)
        self.assertEqual(rows, 17)                               # safety and pharmacokinetics: 1 of 2 each
        code, out = self.gate(frozen, rows)
        self.assertEqual(code, 0, out)
        self.assertIn("sample:pharmacokinetics 1 of 2 (weight 2)", out)
        self.assertIn("sample:safety_tolerability 1 of 2 (weight 2)", out)
        self.assertIn("intervals: stratified bootstrap, 2000 resamples, seed 3, 95% percentile", out)
        self.assertIn("weighted positives 2", out)               # the one drawn pharmacokinetics row counts twice

    def sample_member(self, stratum, small_max=1):
        tags = json.loads(self.tags_path.read_text(encoding="utf-8"))
        accepted, _ = ca.verify_model_tags(sf.good_model_tags(), ca.entry_texts(self.s), {e["entry_id"]: e for e in tags["entries"]}, LEX)
        derived = ca.census_strata(ca.census_frame(tags["entries"], accepted), LEX, small_max)
        return sorted(e for e, st in derived.items() if st == stratum)

    def test_a_sample_exclusion_binds_through_the_gate(self):
        dev = self.sample_member("sample:safety_tolerability")[:1] + ["NCT09999999:P1"]
        devf = self.work / "dev-ids.json"
        devf.write_text(json.dumps(dev), encoding="utf-8")
        frozen, rows = self.census(small_max=1, large_per_class=1, exclude=dev)
        self.assertEqual(rows, 17)
        code, out = self.gate(frozen + ["--frozen-exclude", devf], rows)
        self.assertEqual(code, 0, out)
        self.assertIn("frozen set: sample strata exclude 1 development entries (hash checked)\n", out)
        self.assertIn("sample:safety_tolerability 1 of 1 (weight 1)", out)          # the population without the excluded entry
        # without the file, or with another file, the same key is refused
        code, out = self.gate(frozen, rows)
        self.assertEqual(code, 3, out)
        self.assertIn("no --frozen-exclude file was given", out)
        other = self.work / "other-ids.json"
        other.write_text(json.dumps(dev[:1]), encoding="utf-8")
        code, out = self.gate(frozen + ["--frozen-exclude", other], rows)
        self.assertEqual(code, 3, out)
        self.assertIn("the --frozen-exclude file's hash differs", out)
        self.assertNotIn("(hash checked)", out)
        # an excluded entry drawn into its sample stratum is refused
        frozen, rows = self.census(small_max=1, large_per_class=1, exclude=dev, mutate=lambda k: k["rows"].__setitem__(
            next(i for i, r in enumerate(k["rows"]) if r["stratum"] == "sample:safety_tolerability"),
            {**next(r for r in k["rows"] if r["stratum"] == "sample:safety_tolerability"), "entry_id": dev[0]}))
        code, out = self.gate(frozen + ["--frozen-exclude", devf], rows)
        self.assertEqual(code, 3, out)
        self.assertIn("row(s) of a sample stratum are excluded (development) entries", out)

    def test_an_exclusion_file_needs_a_census_key_that_records_one(self):
        devf = self.work / "dev-ids.json"
        devf.write_text(json.dumps(["NCT00000001:P1"]), encoding="utf-8")
        frozen, rows = self.census(small_max=1, large_per_class=1)                   # sample_exclusion absent: read as null
        code, out = self.gate(frozen + ["--frozen-exclude", devf], rows)
        self.assertEqual(code, 3, out)
        self.assertIn("the key was drawn with no sample exclusion, but --frozen-exclude was given", out)
        code, out = self.gate(self.frozen() + ["--frozen-exclude", devf], 19)            # a legacy key
        self.assertEqual(code, 3, out)
        self.assertIn("--frozen-exclude was given, but the key is not a census-v1 key", out)
        code, out = run(self.base() + ["--frozen-exclude", devf, "--drift-only"])
        self.assertEqual(code, 2, out)
        self.assertIn("--frozen-exclude goes with --frozen-sheet and --frozen-key", out)
        devf.write_text(json.dumps({"ids": []}), encoding="utf-8")
        code, out = self.gate(frozen + ["--frozen-exclude", devf], rows)
        self.assertEqual(code, 2, out)
        self.assertIn("must hold a JSON list of entry ids", out)

    def test_no_labelled_positive_still_blocks_in_census_mode(self):
        frozen, rows = self.census(label=lambda e, cls: ["other"] if "pharmacokinetics" in cls else cls)
        code, out = self.gate(frozen, rows)
        self.assertEqual(code, 3, out)
        self.assertIn("STOP S1 shown class(es) with no labelled positive: pharmacokinetics", out)

    def test_fewer_than_50_labelled_rows_still_blocks_in_census_mode(self):
        frozen, rows = self.census()
        with mock.patch.object(ca, "MIN_CONTROL_ENTRIES", 8):
            code, out = run(self.base() + frozen + self.controls() + ["--explained", self.explained], allow_synthetic_for_tests=True)
        self.assertEqual(code, 3, out)
        self.assertIn("STOP S1 only 19 labelled rows; the design needs at least 50", out)

    def test_a_key_that_does_not_describe_the_frame_blocks(self):
        def other_hash(k):
            k["model_tags_sha256"] = "0" * 64

        def bad_weight(k):
            k["strata"]["census:small"]["weight"] = 2
        for mutate, text in ((other_hash, "the key's model_tags_sha256 does not match the model-tags file given"),
                             (bad_weight, "stratum weight(s) are not population/drawn: census:small")):
            with self.subTest(text=text):
                frozen, rows = self.census(mutate=mutate)
                code, out = self.gate(frozen, rows)
                self.assertEqual(code, 3, out)
                self.assertIn("STOP S1 the frozen set does not belong to this snapshot and lexicon", out)
                self.assertIn(text, out)

    def test_a_legacy_key_is_unchanged_and_an_unknown_design_is_refused(self):
        with mock.patch.object(ca, "MIN_FROZEN_ROWS", 19), mock.patch.object(ca, "MIN_CONTROL_ENTRIES", 8):
            code, out = run(self.base() + self.frozen() + self.controls() + ["--explained", self.explained],
                            allow_synthetic_for_tests=True)
        self.assertEqual(code, 0, out)
        self.assertNotIn("frozen set design", out)
        self.assertIn("too few to estimate (", out)
        frozen = self.frozen()
        key = json.loads(Path(frozen[3]).read_text(encoding="utf-8"))
        key["design"] = "census-v9"
        Path(frozen[3]).write_text(json.dumps(key), encoding="utf-8")
        code, out = self.gate(frozen, 19)
        self.assertEqual(code, 3, out)
        self.assertIn("the key's design 'census-v9' is not one this gate knows", out)

    def test_a_census_key_with_its_design_removed_is_not_scored_as_legacy(self):
        for field in ("strata", "small_max", "large_per_class", "model_tags_sha256"):
            with self.subTest(field=field):
                def strip(k, field=field):
                    del k["design"]
                    for f in ("strata", "small_max", "large_per_class", "model_tags_sha256"):
                        if f != field:
                            k.pop(f, None)
                frozen, rows = self.census(mutate=strip)
                code, out = self.gate(frozen, rows)
                self.assertEqual(code, 3, out)
                self.assertIn(f"the key looks like a census key but names no design (it holds {field})", out)
                self.assertNotIn("micro F1", out)                # never scored

    def test_the_consensus_record_is_printed_when_the_key_holds_it(self):
        def with_record(k):
            k["consensus"] = {"rule_declared": "2026-10-11", "proposed": 9, "kept": 5, "dropped": 4, "report_sha256": "0" * 64}
        frozen, rows = self.census(mutate=with_record)
        code, out = self.gate(frozen, rows)
        self.assertEqual(code, 0, out)
        self.assertIn("consensus rule: declared 2026-10-11, proposed 9, kept 5, dropped 4\n", out)
        frozen, rows = self.census()
        code, out = self.gate(frozen, rows)
        self.assertNotIn("consensus rule:", out)                 # consensus null: no line

    def test_a_shown_class_the_metrics_do_not_hold_blocks(self):
        real = ca.metrics_weighted

        def lose_a_class(*args, **kw):                           # as if the counts named a class the lexicon's metrics lack
            m = real(*args, **kw)
            del m["per_class"]["pharmacokinetics"]
            return m
        frozen, rows = self.census()
        with mock.patch.object(ca, "metrics_weighted", side_effect=lose_a_class):
            code, out = self.gate(frozen, rows)
        self.assertEqual(code, 3, out)
        self.assertIn("STOP S1 shown class(es) the frozen-set metrics do not hold: pharmacokinetics", out)

    def test_the_plain_route_is_counted_beside_the_consensus_rule(self):
        plain = self.work / "plain.json"
        plain.write_text(json.dumps(sf.good_model_tags() + [{"entry_id": "NCT00000014:P3", "class": "pro_unnamed",
                                                            "quote": "return to school"}]), encoding="utf-8")
        code, out = run(self.base() + ["--model-tags-plain", plain, "--drift-only"])
        self.assertEqual(code, 0, out)                           # counts only: nothing changes
        self.assertIn("model route: plain 2 accepted, with the consensus rule 1 accepted, 1 dropped", out)
        code, out = run(["--snapshot", self.snap_dir, "--tags", self.tags_path, "--counts", self.counts_path,
                         "--model-tags-plain", plain])
        self.assertEqual(code, 2, out)
        self.assertIn("--model-tags-plain goes with --model-tags", out)


def no_pk(eid, cls):
    """A labeller who never gives pharmacokinetics (the synthetic frame's two pk entries are NCT00000005:P1 and P2)."""
    return ["other"] if "pharmacokinetics" in cls else cls


class WithholdTest(Fixture):
    """--withhold-failing, --record and --check-record on the synthetic census frame."""
    census = CensusGateTest.census

    def gate(self, frozen, min_rows, extra=(), withhold=True):
        with mock.patch.object(ca, "MIN_FROZEN_ROWS", min_rows), mock.patch.object(ca, "MIN_CONTROL_ENTRIES", 8):
            return run(self.base() + frozen + self.controls() + ["--explained", self.explained]
                       + (["--withhold-failing"] if withhold else []) + list(extra), allow_synthetic_for_tests=True)

    def test_a_failing_class_is_withheld_where_it_would_otherwise_block(self):
        frozen, rows = self.census(label=no_pk)
        code, out = self.gate(frozen, rows, withhold=False)
        self.assertEqual(code, 3, out)
        self.assertIn("STOP S1 shown class(es) with no labelled positive: pharmacokinetics", out)
        code, out = self.gate(frozen, rows)
        self.assertEqual(code, 0, out)
        self.assertIn("WITHHELD pharmacokinetics: no labelled positive\n", out)
        self.assertNotIn("STOP", out)
        self.assertNotIn("note: pharmacokinetics", out)                      # withheld, so not a 'few checked' class on the page
        self.assertIn("note: fev1 is shown but few checked", out)             # few checked stays shown
        # by hand: NCT00000005:P1 and P2 carry only pharmacokinetics, 2 of 20 entries; study NCT00000005 has no other class
        self.assertIn("not shown: 2 entries (10.0%) carry only withheld classes; 1 studies fall only in withheld classes\n", out)
        self.assertIn("RESULT: PUBLICATION ALLOWED WITH 1 WITHHELD CLASS(ES)\n", out)

    def test_precision_and_recall_reasons_are_computed(self):
        def label(eid, cls):
            return ["sweat_chloride", "fev1"] if eid == "NCT00000001:P2" else no_pk(eid, cls)
        frozen, rows = self.census(label=label)
        with mock.patch.dict(ca.THRESHOLDS, {"min_positives": 1}):
            code, out = self.gate(frozen, rows)
        self.assertEqual(code, 0, out)
        withheld = [ln for ln in out.splitlines() if ln.startswith("WITHHELD")]
        self.assertEqual(withheld, ["WITHHELD fev1: recall 0.500 below 0.80",
                                    "WITHHELD pharmacokinetics: no labelled positive; precision 0.000 below 0.85"])
        # NCT00000001:P1 (fev1 only) and the two pk entries; NCT00000001 keeps sweat_chloride, so one study only
        self.assertIn("not shown: 3 entries (15.0%) carry only withheld classes; 1 studies fall only in withheld classes", out)
        self.assertIn("RESULT: PUBLICATION ALLOWED WITH 2 WITHHELD CLASS(ES)", out)

    def test_every_other_stop_still_blocks(self):
        first = []

        def skip_one(eid, cls):                                  # one row left blank, besides the failing pk class
            if not first and "pharmacokinetics" not in cls:
                first.append(eid)
                return []
            return no_pk(eid, cls)
        frozen, rows = self.census(label=skip_one)
        code, out = self.gate(frozen, rows - 1)
        self.assertEqual(code, 3, out)
        self.assertIn("STOP S1 1 row(s) of the frozen sheet are unlabelled", out)
        frozen, rows = self.census(label=no_pk)
        with mock.patch.object(ca, "MIN_CONTROL_ENTRIES", 8):            # MIN_FROZEN_ROWS stays 50
            code, out = run(self.base() + frozen + self.controls() + ["--explained", self.explained, "--withhold-failing"],
                            allow_synthetic_for_tests=True)
        self.assertEqual(code, 3, out)
        self.assertIn("STOP S1 only 19 labelled rows; the design needs at least 50", out)
        frozen, rows = self.census(label=no_pk, mutate=lambda k: k["strata"]["census:small"].__setitem__("weight", 2))
        code, out = self.gate(frozen, rows)
        self.assertEqual(code, 3, out)
        self.assertIn("STOP S1 the frozen set does not belong to this snapshot and lexicon", out)
        self.assertNotIn("WITHHELD", out)                        # nothing is scored on a key that does not bind
        frozen, rows = self.census(label=no_pk)
        with mock.patch.dict(ca.THRESHOLDS, {"other_share": 0.01}):
            code, out = self.gate(frozen, rows)
        self.assertEqual(code, 3, out)
        self.assertIn("WITHHELD pharmacokinetics: no labelled positive", out)
        self.assertIn("STOP S2 other plus unclassified are 5.0% of entries", out)
        self.assertIn("RESULT: PUBLICATION BLOCKED (1 stop rule(s)); integrity passed", out)
        code, out = run(self.base() + frozen + ["--negative-controls", "--withhold-failing"], allow_synthetic_for_tests=True)
        self.assertEqual(code, 3, out)                           # S3 and S4 too
        self.assertIn("STOP S3", out)
        self.assertIn("STOP S4", out)

    def test_few_checked_classes_are_not_withheld(self):
        frozen, rows = self.census()
        code, out = self.gate(frozen, rows)
        self.assertEqual(code, 0, out)
        self.assertNotIn("WITHHELD", out)
        self.assertEqual(out.count("is shown but few checked"), 17)
        self.assertIn("not shown: 0 entries (0.0%) carry only withheld classes; 0 studies fall only in withheld classes", out)
        self.assertIn("RESULT: PASS (integrity, controls and every publication stop rule)", out)

    def test_no_class_is_withheld_because_someone_typed_it(self):
        frozen, rows = self.census(mutate=lambda k: k.update({"withheld": ["fev1"], "withhold": ["lci"]}))
        code, out = self.gate(frozen, rows)
        self.assertEqual(code, 0, out)
        self.assertNotIn("WITHHELD", out)                        # a key field naming classes is ignored
        for opt in (["--withhold-class", "fev1"], ["--withheld", "fev1"], ["--withhold-failing", "fev1"]):
            with self.subTest(opt=opt), contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit) as cm:
                self.gate(frozen, rows, extra=opt)
            self.assertEqual(cm.exception.code, 2)

    def test_usage_errors(self):
        frozen, rows = self.census(label=no_pk)
        code, out = self.gate(self.frozen(), 19)                 # a legacy key
        self.assertEqual(code, 2, out)
        self.assertIn("--withhold-failing goes with a census-v1 key only; this key's design is not named (a legacy key)", out)
        code, out = run(self.base() + ["--withhold-failing"])
        self.assertEqual(code, 2, out)
        self.assertIn("--withhold-failing needs --frozen-sheet", out)
        code, out = run(self.base() + frozen + ["--withhold-failing", "--drift-only"])
        self.assertEqual(code, 2, out)
        rec = self.work / "rec.json"
        code, out = self.gate(frozen, rows, extra=["--record", rec], withhold=False)
        self.assertEqual(code, 2, out)
        self.assertIn("--record goes with --withhold-failing", out)
        self.assertFalse(rec.exists())
        code, out = self.gate(frozen, rows, extra=["--record", self.work / "no-such-folder" / "rec.json"])
        self.assertEqual(code, 2, out)
        code, out = run(["--check-record", rec, "--snapshot", self.snap_dir])
        self.assertEqual(code, 2, out)
        self.assertIn("--check-record stands alone", out)
        code, out = run(["--check-record", self.work / "absent.json"])
        self.assertEqual(code, 2, out)

    def test_not_shown_counts_by_hand(self):
        entries = [tag_entry("A:1", ["fev1"]), tag_entry("A:2", ["lci"]),               # study A: fev1 and lci entries
                   tag_entry("B:1", ["fev1", "lci"]),                                   # both withheld: only withheld
                   tag_entry("C:1", ["fev1", "sweat_chloride"]),                        # a shown class too: not counted
                   tag_entry("D:1", []), tag_entry("D:2", ["fev1"]),                    # D: an unclassified entry has no class
                   tag_entry("E:1", []),                                                # model tag below gives it lci
                   tag_entry("F:1", [])]                                                # no class at all: never 'only withheld'
        accepted = [{"entry_id": "E:1", "class": "lci"}]
        got = ca.not_shown_counts(entries, accepted, {"fev1", "lci"})
        # entries A:1, A:2, B:1, D:2, E:1 -> 5 of 8; studies A, B, D, E (C has sweat_chloride, F has no class) -> 4
        self.assertEqual(got, {"entries": 5, "entries_share": 0.625, "studies": 4})
        self.assertEqual(ca.not_shown_counts(entries, accepted, {"fev1"}), {"entries": 2, "entries_share": 0.25, "studies": 1})
        self.assertEqual(ca.not_shown_counts(entries, accepted, set()), {"entries": 0, "entries_share": 0.0, "studies": 0})

    # ---- the run record

    def record(self, name="rec.json", label=no_pk):
        frozen, rows = self.census(label=label)
        path = self.work / name
        code, out = self.gate(frozen, rows, extra=["--record", path])
        self.assertEqual(code, 0, out)
        self.assertIn(f"run record written: {name}", out)
        return path

    def test_the_record_is_deterministic_id_free_and_never_overwritten(self):
        a, b = self.record("a.json"), self.record("b.json")
        self.assertEqual(a.read_bytes(), b.read_bytes())
        text = a.read_text(encoding="utf-8")
        self.assertTrue(text.endswith("}\n"))
        self.assertEqual(text, ca.canonical_json(json.loads(text)))
        self.assertNotRegex(text, r"NCT\d")                     # no entry or study id
        self.assertNotIn("S001", text)                           # no sample id
        for e in ca.entry_texts(self.s).values():                # no registry wording
            for f in ("measure", "description"):
                if len(e.get(f) or "") > 8:
                    self.assertNotIn(e[f], text)
        rec = json.loads(text)
        self.assertEqual(set(rec), set(ca.RECORD_FIELDS))
        self.assertEqual(rec["withheld"], ["pharmacokinetics"])
        self.assertEqual(rec["withheld_sha256"], __import__("hashlib").sha256(b"pharmacokinetics").hexdigest())
        self.assertEqual(rec["record_sha256"], ca.record_digest(rec))
        self.assertEqual(rec["verdict"], "allowed_with_withheld")
        self.assertEqual((rec["approval"], rec["reference_labels"]), (None, ca.REFERENCE_LABELS))
        self.assertEqual(rec["not_shown"], {"entries": 2, "entries_share": 0.1, "studies": 1})
        self.assertEqual(rec["key_sha256"], __import__("hashlib").sha256((self.work / "census-key.json").read_bytes()).hexdigest())
        pk = next(p for p in rec["per_class"] if p["class"] == "pharmacokinetics")
        self.assertEqual((pk["status"], pk["reason"], pk["labelled_positive_rows"]), ("withheld", "no labelled positive", 0))
        self.assertEqual({p["status"] for p in rec["per_class"] if p["class"] != "pharmacokinetics"}, {"few_checked"})
        before = a.read_bytes()
        frozen, rows = self.census(label=no_pk)
        code, out = self.gate(frozen, rows, extra=["--record", a])
        self.assertEqual(code, 2, out)
        self.assertIn("exists already; a record is never overwritten", out)
        self.assertEqual(a.read_bytes(), before)

    def check(self, path):
        return run(["--check-record", path])

    def edited(self, path, change, rehash=True):
        rec = json.loads(path.read_text(encoding="utf-8"))
        change(rec)
        if rehash:                                               # a careful forger recomputes the record hash
            rec["record_sha256"] = ca.record_digest(rec)
        out = self.work / "edited.json"
        out.write_text(ca.canonical_json(rec), encoding="utf-8")
        return out

    def test_check_record_passes_a_good_record(self):
        code, out = self.check(self.record())
        self.assertEqual(code, 0, out)
        self.assertNotIn("FAIL", out)
        self.assertIn("PASS class pharmacokinetics: recorded withheld, the numbers and thresholds give withheld (no labelled "
                      "positive)", out)

    def test_check_record_fails_each_edit(self):
        path = self.record()

        def drop_withheld(r):
            r["withheld"] = []
            r["withheld_sha256"] = ca.withheld_digest([])

        def unfail(r):                                           # pk gains labelled positives, so it no longer fails
            pk = next(p for p in r["per_class"] if p["class"] == "pharmacokinetics")
            pk["labelled_positive_rows"], pk["weighted_positives"] = 2, 2.0

        def verdict(r):
            r["verdict"] = "allowed"

        def add_withheld(r):                                     # a class typed into the list by hand
            r["withheld"] = ["fev1", "pharmacokinetics"]
            r["withheld_sha256"] = ca.withheld_digest(r["withheld"])
        cases = ((drop_withheld, True, "FAIL withheld is sorted and equals the classes whose status is withheld"),
                 (lambda r: r.update(withheld=["pharmacokinetics", "fev1"]), True, "FAIL withheld_sha256 matches"),
                 (add_withheld, True, "FAIL withheld is sorted and equals"),
                 (unfail, True, "FAIL class pharmacokinetics: recorded withheld, the numbers and thresholds give few_checked"),
                 (verdict, True, "FAIL verdict allowed agrees with"),
                 (lambda r: r.update(blocking_rules=["S2"]), True, "FAIL verdict allowed_with_withheld agrees with"),
                 (lambda r: r.update(approval="someone"), True, "FAIL reference_labels is the fixed text and approval is null"),
                 (lambda r: r.update(record_sha256="0" * 64), False, "FAIL record_sha256 matches the record"),
                 (lambda r: r["thresholds"].update(recall=0.5), True, "FAIL the record's thresholds equal this gate's"),
                 (lambda r: r.pop("verdict"), True, "FAIL format"))
        for change, rehash, text in cases:
            with self.subTest(text=text):
                code, out = self.check(self.edited(path, change, rehash))
                self.assertEqual(code, 1, out)
                self.assertIn(text, out)
        # an edit without a rehash is caught by the record hash
        code, out = self.check(self.edited(path, verdict, rehash=False))
        self.assertEqual(code, 1, out)
        self.assertIn("FAIL record_sha256 matches the record", out)


class SelfConsistencyTest(Fixture):
    """--self-consistency: the relationships a counts file must satisfy, the same rules a page generator applies."""

    def check(self, mutate=None):
        counts = json.loads(self.counts_path.read_text(encoding="utf-8"))
        if mutate:
            mutate(counts)
        p = self.work / "counts-consistency.json"
        p.write_text(json.dumps(counts), encoding="utf-8")
        return run(["--self-consistency", p])

    def test_counts_written_by_the_tools_are_consistent(self):
        code, out = self.check()
        self.assertEqual(code, 0, out)
        self.assertEqual(out.count("\nOK ") + out.startswith("OK "), len(ca.CONSISTENCY_RULES))
        self.assertNotIn("FAIL", out)

    def test_each_relationship_fails_on_a_perturbed_file(self):
        def bump_class_studies(c):
            c["studies_by_class"]["fev1"] = c["studies"] + 1

        def zero_class_entries(c):
            c["entries_by_class"] = {k: 0 for k in c["entries_by_class"]}

        def phase_total(c):
            c["class_year_phase"]["exacerbations"]["planned"]["2027"]["PHASE3"] = 2

        cases = (
            ("class studies within studies in scope", bump_class_studies),
            ("class entries cover classified entries", zero_class_entries),
            ("unsorted entries equal unclassified", lambda c: c.update({"unsorted_entries": 3})),
            ("start kinds add up to studies in scope", lambda c: c["scope"].update({"start_untyped": c["scope"]["start_untyped"] + 1})),
            ("class_year_phase totals equal year counts", phase_total),
            ("co-occurrence pairs within the smaller class", lambda c: c["co_occurrence"]["pairs"][0].update({"studies": 99})),
            ("distinct sponsors within studies", lambda c: c["lead_sponsors_by_class"].update({"exacerbations": 99})),
        )
        self.assertEqual(sorted(n for n, _ in cases), sorted(ca.CONSISTENCY_RULES))
        for name, mutate in cases:
            with self.subTest(name=name):
                code, out = self.check(mutate)
                self.assertEqual(code, 1, out)
                self.assertIn(f"FAIL {name}", out)
                self.assertEqual(out.count("FAIL "), 1, out)              # only the relationship that was broken

    def test_a_missing_key_fails_and_an_unreadable_file_is_a_usage_error(self):
        code, out = self.check(lambda c: c.pop("unsorted_entries"))
        self.assertEqual(code, 1, out)
        self.assertIn("FAIL unsorted entries equal unclassified", out)
        code, out = run(["--self-consistency", self.work / "nowhere.json"])
        self.assertEqual(code, 2, out)

    def test_write_counts_checks_before_writing(self):
        target = self.work / "counts-new.json"
        code, out = run(["--snapshot", self.snap_dir, "--tags", self.tags_path, "--model-tags", self.model_path, "--write-counts", target])
        self.assertEqual(code, 0, out)
        self.assertIn("self-consistency: OK", out)
        with mock.patch.object(ca, "consistency_results", lambda c: [("unsorted entries equal unclassified", False, "forced")]):
            target2 = self.work / "counts-refused.json"
            code, out = run(["--snapshot", self.snap_dir, "--tags", self.tags_path, "--model-tags", self.model_path,
                             "--write-counts", target2])
        self.assertEqual(code, 1, out)
        self.assertFalse(target2.exists())


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

    # ---- the gate refuses synthetic data and unbound frozen sets

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
        # a planted id that only the term search returned is still tested against the scope rules
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
