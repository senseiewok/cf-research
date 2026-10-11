"""Tests for challenge_tags.py on SYNTHETIC entries with a FAKE invoker. No model, no network.
Usage: python -m unittest -v (inside tools/trial_atlas)."""
import json
import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import challenge_tags as ct  # noqa: E402
import check_atlas as ca  # noqa: E402
import lexicon  # noqa: E402
import propose_tags as pt  # noqa: E402
import snapshot as snap  # noqa: E402
import synthetic_fixtures as sf  # noqa: E402
from test_propose_tags import build, oracle, quiet  # noqa: E402

LEX = lexicon.load()
PROPOSER = "fake-proposer:1"


class Base(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._tmp = tempfile.TemporaryDirectory()
        cls.root = Path(cls._tmp.name)
        cls.snap_dir, cls.tags_path = build(cls.root)
        cls.entries = ca.entry_texts(snap.load(cls.snap_dir))
        cls.rule_tags = {e["entry_id"]: e for e in json.loads(cls.tags_path.read_text(encoding="utf-8"))["entries"]}
        cls.prop = cls.root / "proposer"
        code, out = quiet(pt.main, ["--snapshot", str(cls.snap_dir), "--tags", str(cls.tags_path), "--out", str(cls.prop),
                                    "--model", PROPOSER], invoker=pt.ScriptedInvoker(pt.self_check_script()))
        assert code == 0, out
        cls.model_tags = json.loads((cls.prop / pt.MODEL_TAGS).read_text(encoding="utf-8"))
        assert cls.model_tags

    @classmethod
    def tearDownClass(cls):
        cls._tmp.cleanup()

    def setUp(self):
        sf.block_network(self)
        self.out = Path(tempfile.mkdtemp(dir=self.root)) / "challenge"

    def args(self, *extra, model="fake-challenger:1", proposals=None, sample="12"):
        a = ["--snapshot", str(self.snap_dir), "--tags", str(self.tags_path), "--proposals", str(proposals or self.prop),
             "--out", str(self.out), "--sample", sample, "--seed", "5", *extra]
        return a + (["--model", model] if model else [])

    def refs(self):
        out = {t["entry_id"]: [t["class"]] for t in self.model_tags}
        for e, t in self.rule_tags.items():
            if t["status"] == "rule":
                out[e] = [x["class"] for x in t["tags"]]
        return out

    def agreeing(self):
        """An oracle that answers every entry with its first reference class (its measure as the quote)."""
        return oracle({self.entries[e]["measure"]: c[0] for e, c in self.refs().items() if e in self.entries})


class RefuseTest(Base):
    def test_same_model_refused(self):
        code, out = quiet(ct.main, self.args(model=PROPOSER), invoker=self.agreeing())
        self.assertEqual(code, 2)
        self.assertIn("proposer's model", out)

    def test_same_family_name_refused(self):
        code, out = quiet(ct.main, self.args(model="fake-proposer:2-64k"), invoker=self.agreeing())
        self.assertEqual(code, 2)
        self.assertIn("same family name", out)

    def test_profile_with_the_proposers_model_refused(self):
        p = self.root / "ollama-profile.same.json"
        p.write_text(json.dumps({"model": PROPOSER}), encoding="utf-8")
        code, out = quiet(ct.main, self.args("--profile-file", str(p), model=None), invoker=self.agreeing())
        self.assertEqual(code, 2)
        self.assertIn("proposer's model", out)

    def test_proposer_tags_that_do_not_verify_refused(self):
        bad = self.root / "bad-proposer"
        bad.mkdir(exist_ok=True)
        (bad / pt.RUN).write_text((self.prop / pt.RUN).read_text(encoding="utf-8"), encoding="utf-8")
        (bad / pt.MODEL_TAGS).write_text(json.dumps(sf.canary_model_tags()), encoding="utf-8")
        code, out = quiet(ct.main, self.args(proposals=bad), invoker=self.agreeing())
        self.assertEqual(code, 2)
        self.assertIn("do not verify", out)


class BlindTest(Base):
    def test_challenger_never_receives_the_proposers_answer(self):
        """The packets sent are the same whatever the proposer said, and equal the bare packet of the entry text and class list."""
        other = self.root / "other-proposer"
        other.mkdir(exist_ok=True)
        (other / pt.RUN).write_text((self.prop / pt.RUN).read_text(encoding="utf-8"), encoding="utf-8")
        changed = [{**t, "class": "imaging" if t["class"] != "imaging" else "liver"} for t in self.model_tags]
        (other / pt.MODEL_TAGS).write_text(json.dumps(changed), encoding="utf-8")
        sent = []
        for proposals in (self.prop, other):
            fake = pt.ScriptedInvoker({})
            self.out = Path(tempfile.mkdtemp(dir=self.root)) / "c"
            code, out = quiet(ct.main, self.args(proposals=proposals), invoker=fake)
            self.assertEqual(code, 0, out)
            sent.append([(r.attempt, r.prompt) for r in fake.requests])
        self.assertEqual(sent[0], sent[1])
        classes = pt.offered_classes(LEX)
        sample = json.loads((self.out / ct.RUN).read_text(encoding="utf-8"))["sample"]
        first = [p for a, p in sent[0] if a == 1]
        self.assertEqual(first, [pt.build_packet(self.entries[e], classes) for e, _ in sample])
        for _, p in sent[0]:
            self.assertTrue(pt._boundary_ok(p))
            for word in (PROPOSER, "proposer", "rule_id", "pool", "reference"):
                self.assertNotIn(word, p)


class ResultTest(Base):
    def test_sample_mixes_pools_and_agreement_is_counted(self):
        code, out = quiet(ct.main, self.args(), invoker=self.agreeing())
        self.assertEqual(code, 0, out)
        res = json.loads((self.out / ct.OUT_JSON).read_text(encoding="utf-8"))
        n_model = len({t["entry_id"] for t in self.model_tags})
        self.assertEqual(res["sample_drawn"], {"model": n_model, "rule": 12 - n_model})
        compared = [e for e in res["entries"] if e["verified"]]
        self.assertEqual(res["unverified"], 12 - len(compared))
        self.assertTrue(all(e["agree"] for e in compared))
        o = res["agreement"]["overall"]
        self.assertEqual((o["agree"], o["compared"]), (len(compared), len(compared)))
        if len(compared) >= 10:
            lo, hi = ca.wilson(len(compared), len(compared))
            self.assertEqual(o["wilson95"], [round(lo, 4), round(hi, 4)])
        else:
            self.assertEqual(o["wilson95"], "too few to estimate")
        self.assertEqual(res["disagreements"], [])
        self.assertEqual(res["model"], "fake-challenger:1")
        self.assertEqual(res["proposer_models"], [PROPOSER])

    def test_disagreements_list_ids_and_classes_only(self):
        code, _ = quiet(ct.main, self.args(), invoker=oracle({}))     # always none: every verified answer disagrees
        self.assertEqual(code, 0)
        res = json.loads((self.out / ct.OUT_JSON).read_text(encoding="utf-8"))
        verified = [e for e in res["entries"] if e["verified"]]
        self.assertTrue(verified)
        self.assertEqual(len(res["disagreements"]), len(verified))
        self.assertEqual(res["agreement"]["overall"]["agree"], 0)
        for d in res["disagreements"]:
            self.assertEqual(set(d), {"entry_id", "pool", "reference_classes", "challenger_class"})
        text = json.dumps(res)
        for e in self.entries.values():
            if e["measure"]:
                self.assertNotIn(e["measure"], text)

    def test_rule_entries_are_quote_checked_like_the_proposers(self):
        """A true quote on a rule-tagged entry verifies; a false one is rejected."""
        rows_ok = self.run_one(lambda m: {"class_id": "fev1", "quote": m, "family": ""})
        rows_bad = self.run_one(lambda m: {"class_id": "fev1", "quote": m + " extra words here", "family": ""})
        self.assertTrue(any(r["verified"] for r in rows_ok))
        self.assertFalse(any(r["verified"] for r in rows_bad))
        self.assertTrue(all(r["reason"] in ("quote_not_substring", "quote_length") for r in rows_bad))
        self.assertTrue(all(r["role"] == "challenger" for r in rows_ok + rows_bad))

    def run_one(self, reply):
        def fake(req):
            m = req.prompt.split("measure: ", 1)[1].split("\n", 1)[0]
            return pt.InvokeResult(0, json.dumps(reply(m)))
        self.out = Path(tempfile.mkdtemp(dir=self.root)) / "c"
        code, out = quiet(ct.main, self.args(sample="6"), invoker=fake)
        self.assertEqual(code, 0, out)
        rows, _ = pt.load_rows(self.out / ct.JSONL)
        sample = dict(json.loads((self.out / ct.RUN).read_text(encoding="utf-8"))["sample"])
        return [r for r in rows if sample[r["entry_id"]] == "rule"]

    def test_rate_too_few_under_ten(self):
        self.assertEqual(ct.rate(3, 9)["wilson95"], "too few to estimate")
        self.assertEqual(ct.rate(0, 0)["rate"], None)
        lo, hi = ca.wilson(9, 10)
        self.assertEqual(ct.rate(9, 10)["wilson95"], [round(lo, 4), round(hi, 4)])

    def test_draw_tops_up_from_the_other_pool_and_is_seeded(self):
        a = ct.draw_sample(["m1"], ["r1", "r2", "r3", "r4"], 4, 1)
        self.assertEqual(sorted(p for _, p in a), ["model", "rule", "rule", "rule"])
        self.assertEqual(a, ct.draw_sample(["m1"], ["r1", "r2", "r3", "r4"], 4, 1))
        self.assertEqual(len(ct.draw_sample(["m1"], ["r1"], 10, 1)), 2)


if __name__ == "__main__":
    unittest.main()
