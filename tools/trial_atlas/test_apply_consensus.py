"""Tests for apply_consensus.py on SYNTHETIC data. No model, no network.
Usage: python -m unittest -v test_apply_consensus (inside tools/trial_atlas)."""
import contextlib
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import apply_consensus as ac  # noqa: E402
import check_atlas as ca  # noqa: E402
import lexicon  # noqa: E402
import propose_tags as pt  # noqa: E402
import snapshot as snap  # noqa: E402
from test_propose_tags import build, oracle, quiet  # noqa: E402

LEX = lexicon.load()
PROPOSER_MODEL = "fake-proposer:1"
PROPOSER_FAMILY = "fakeprop"
E = {n: f"NCT0000090{n}:P1" for n in range(1, 7)}
# the self-check studies' measures; every one is left unclassified by the rule lexicon
PROPOSALS = {  # entry number -> (measure, class the fake proposer gives)
    1: ("Time to return to school", "healthcare_use"),
    2: ("Number of nights with disturbed sleep", "pro_unnamed"),
    4: ("Change in hand dexterity score", "healthcare_use"),
    5: ("Duration of morning stiffness", "pro_unnamed"),
    6: ("Time to first school absence", "healthcare_use"),
}


def labels_file(entries, family="gemma", model="gemma4:31b-it-q4_K_M"):
    return {"labeller_family": family, "labeller_model": model, "entries": entries}


# Entry 1 agrees, 2 disagrees, 4 is invalid, 5 is missing, 6 agrees (the labeller also gave another class).
AGREEMENT = {E[1]: {"classes": ["healthcare_use"]}, E[2]: {"classes": ["fev1"]}, E[4]: {"invalid": True},
             E[6]: {"classes": ["fev1", "healthcare_use"]}, "NCT00009999:P1": {"classes": ["fev1"]}}


class Base(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._tmp = tempfile.TemporaryDirectory()
        cls.root = Path(cls._tmp.name)
        cls.snap_dir, cls.tags_path = build(cls.root)
        prop = cls.root / "proposer"
        fake = oracle({m: c for m, c in PROPOSALS.values()})
        code, out = quiet(pt.main, ["--snapshot", str(cls.snap_dir), "--tags", str(cls.tags_path), "--out", str(prop),
                                    "--model", PROPOSER_MODEL], invoker=fake)
        assert code == 0, out
        # the fixture's CF studies share a measure with entry 1, so the fake proposer also tags one of them; keep the five planted rows
        ours = {E[n] for n in PROPOSALS}
        cls.tags = [t for t in json.loads((prop / pt.MODEL_TAGS).read_text(encoding="utf-8")) if t["entry_id"] in ours]
        cls.model_tags_path = prop / "five-tags.json"
        pt.write_json(cls.model_tags_path, cls.tags)
        assert sorted(t["entry_id"] for t in cls.tags) == sorted(E[n] for n in PROPOSALS), cls.tags

    @classmethod
    def tearDownClass(cls):
        cls._tmp.cleanup()

    def setUp(self):
        self.work = Path(tempfile.mkdtemp(dir=self.root))
        self.out = self.work / "consensus"

    def labels(self, data, name="labels.json"):
        p = self.work / name
        p.write_text(json.dumps(data), encoding="utf-8")
        return p

    def run_main(self, labels_path, *, out=None, family=PROPOSER_FAMILY, tags=None):
        return quiet(ac.main, ["--model-tags", str(tags or self.model_tags_path), "--blind-labels", str(labels_path),
                               "--proposer-family", family, "--out", str(out or self.out)])

    def read(self, name, out=None):
        return json.loads(((out or self.out) / name).read_text(encoding="utf-8"))


class RuleTest(Base):
    def test_keeps_agreeing_drops_the_rest(self):
        code, text = self.run_main(self.labels(labels_file(AGREEMENT)))
        self.assertEqual(code, 0, text)
        kept = self.read(ac.MODEL_TAGS)
        self.assertEqual([t["entry_id"] for t in kept], [t["entry_id"] for t in self.tags if t["entry_id"] in (E[1], E[6])])
        for t in kept:
            self.assertIn(t, self.tags)   # kept rows are the proposer's rows, unchanged
        dropped = {d["entry_id"]: d for d in self.read(ac.DROPPED)}
        self.assertEqual({e: d["reason"] for e, d in dropped.items()},
                         {E[2]: "disagree", E[4]: "invalid_reply", E[5]: "missing"})
        self.assertEqual(dropped[E[2]]["class"], "pro_unnamed")

    def test_an_empty_class_set_drops_and_nothing_is_kept_by_default(self):
        code, _ = self.run_main(self.labels(labels_file({E[1]: {"classes": []}})))
        self.assertEqual(code, 0)
        self.assertEqual(self.read(ac.MODEL_TAGS), [])
        rep = self.read(ac.REPORT)
        self.assertEqual((rep["proposed"], rep["kept"], rep["dropped"], rep["dropped_disagree"], rep["dropped_missing"]),
                         (5, 0, 5, 1, 4))

    def test_output_loads_in_check_atlas_reader(self):
        self.run_main(self.labels(labels_file(AGREEMENT)))
        proposed = ca.read_json(self.out / ac.MODEL_TAGS, "model tags")
        s = snap.load(self.snap_dir)
        tags = json.loads(self.tags_path.read_text(encoding="utf-8"))
        rule_tags = {e["entry_id"]: e for e in tags["entries"]}
        accepted, rejected = ca.verify_model_tags(proposed, ca.entry_texts(s), rule_tags, LEX)
        self.assertEqual(len(accepted), 2)
        self.assertEqual(rejected, [])
        # the same bytes the proposer's writer produces for the same list
        ref = self.work / "ref.json"
        pt.write_json(ref, proposed)
        self.assertEqual((self.out / ac.MODEL_TAGS).read_bytes().replace(b"\r\n", b"\n"), ref.read_bytes().replace(b"\r\n", b"\n"))

    def test_report_is_counts_only_with_rule_date(self):
        self.run_main(self.labels(labels_file(AGREEMENT)))
        rep = self.read(ac.REPORT)
        self.assertEqual(ac.RULE_DECLARED, "2026-10-11")
        self.assertEqual(rep["rule_declared"], "2026-10-11")
        self.assertIn("chosen after the plain route's precision was seen", rep["rule_origin"])
        self.assertEqual((rep["proposed"], rep["kept"], rep["dropped"]), (5, 2, 3))
        self.assertEqual((rep["dropped_disagree"], rep["dropped_invalid_reply"], rep["dropped_missing"]), (1, 1, 1))
        self.assertEqual((rep["proposer_family"], rep["labeller_family"]), (PROPOSER_FAMILY, "gemma"))
        self.assertEqual(rep["proposer_models"], [PROPOSER_MODEL])
        for key, value in rep.items():   # scalars only, apart from the proposer's model ids
            if key != "proposer_models":
                self.assertNotIsInstance(value, (list, dict), key)
        # no entry id and no entry text anywhere in the report or the dropped list
        blob = (self.out / ac.REPORT).read_text(encoding="utf-8") + (self.out / ac.DROPPED).read_text(encoding="utf-8")
        for t in self.tags:
            self.assertNotIn(t["quote"], blob)
        report_text = (self.out / ac.REPORT).read_text(encoding="utf-8")
        for n in E.values():
            self.assertNotIn(n, report_text)
        for m, _ in PROPOSALS.values():
            self.assertNotIn(m, blob)

    def test_deterministic_bytes(self):
        lp = self.labels(labels_file(AGREEMENT))
        other = self.work / "second"
        self.assertEqual(self.run_main(lp)[0], 0)
        self.assertEqual(self.run_main(lp, out=other)[0], 0)
        for name in (ac.MODEL_TAGS, ac.REPORT, ac.DROPPED):
            self.assertEqual((self.out / name).read_bytes(), (other / name).read_bytes(), name)

    def test_label_order_does_not_change_outputs(self):
        items = list(AGREEMENT.items())
        a, b = self.work / "a", self.work / "b"
        self.run_main(self.labels(labels_file(dict(items)), "l1.json"), out=a)
        self.run_main(self.labels(labels_file(dict(reversed(items))), "l1.json"), out=b)
        for name in (ac.MODEL_TAGS, ac.DROPPED):
            self.assertEqual((a / name).read_bytes(), (b / name).read_bytes())


class RefusalTest(Base):
    def assert_nothing_written(self, out):
        self.assertFalse(out.exists() and any(out.iterdir()), list(out.iterdir()) if out.exists() else None)

    def test_same_family_refused_and_nothing_written(self):
        for family, model in (("fakeprop", "other-model:1"), ("FakeProp ", "other-model:1"), ("gemma", "fake-proposer:1"),
                              ("gemma", "FAKE-PROPOSER:2b"), ("fake-proposer", "gemma4:31b")):
            with self.subTest(family=family, model=model):
                code, text = self.run_main(self.labels(labels_file(AGREEMENT, family=family, model=model)))
                self.assertEqual(code, 1, text)
                self.assertIn("REFUSED", text)
                self.assert_nothing_written(self.out)

    def test_two_names_of_one_family_stem_are_refused(self):
        self.assertEqual([ac._stem(n) for n in ("qwen3.8", "qwen3-coder-next:latest", "gemma4", "claude-opus-5-5", "Qwen-Coder")],
                         ["qwen", "qwen", "gemma", "claude", "qwen"])
        for labeller_family, labeller_model in (("qwen-coder", "qwen3-coder-next:latest"),   # the review's case
                                                ("alibaba", "qwen3-coder-next:latest"),      # the model name alone gives it away
                                                ("qwen-coder", "other-model:1")):            # the family alone gives it away
            with self.subTest(labeller_family=labeller_family, labeller_model=labeller_model):
                with self.assertRaises(ac.Refused) as cm:
                    ac.refuse_same_family("qwen", ["qwen3.8:27b"], labeller_family, labeller_model)
                self.assertIn("family stem 'qwen'", str(cm.exception))
        # the stem of a proposer model id counts even when the declared proposer family names something else
        with self.assertRaises(ac.Refused):
            ac.refuse_same_family("alibaba", ["qwen3.8:27b"], "qwen-coder", "coder-next:latest")
        # gemma against qwen is still accepted
        self.assertIsNone(ac.refuse_same_family("qwen", ["qwen3.8:27b"], "gemma", "gemma4:31b-it-q4_K_M"))
        code, text = self.run_main(self.labels(labels_file(AGREEMENT, family="fake-labeller", model="fakelab:1")))
        self.assertEqual(code, 1, text)                         # through main: fake-proposer:1 and fake-labeller share "fake"
        self.assertIn("family stem 'fake'", text)
        self.assert_nothing_written(self.out)

    def test_proposer_family_is_required(self):
        with self.assertRaises(SystemExit) as cm, contextlib.redirect_stderr(io.StringIO()):
            ac.main(["--model-tags", str(self.model_tags_path), "--blind-labels", str(self.labels(labels_file(AGREEMENT))),
                     "--out", str(self.out)])
        self.assertEqual(cm.exception.code, 2)
        self.assert_nothing_written(self.out)

    def test_out_inside_a_git_repo_refused(self):
        repo = self.work / "repo"
        (repo / ".git").mkdir(parents=True)
        out = repo / "sub" / "consensus"
        code, text = self.run_main(self.labels(labels_file(AGREEMENT)), out=out)
        self.assertEqual(code, 1, text)
        self.assertIn("git working tree", text)
        self.assertFalse(out.exists())

    def test_existing_output_is_not_overwritten(self):
        self.out.mkdir()
        (self.out / ac.REPORT).write_text("keep me", encoding="utf-8")
        code, _ = self.run_main(self.labels(labels_file(AGREEMENT)))
        self.assertEqual(code, 1)
        self.assertEqual((self.out / ac.REPORT).read_text(encoding="utf-8"), "keep me")

    def test_malformed_labels_refused_and_nothing_written(self):
        bad = [
            {"labeller_model": "m", "entries": {}},                                            # no family
            {"labeller_family": "gemma", "entries": {}},                                       # no model
            {"labeller_family": "gemma", "labeller_model": "m"},                               # no entries
            labels_file({E[1]: {"classes": "fev1"}}),                                          # classes not a list
            labels_file({E[1]: {"classes": ["fev1"], "invalid": True}}),                       # both keys
            labels_file({E[1]: {"invalid": False}}),                                           # invalid not true
            labels_file({E[1]: ["fev1"]}),                                                     # not an object
            [],
        ]
        for i, data in enumerate(bad):
            with self.subTest(i=i):
                code, text = self.run_main(self.labels(data, f"bad{i}.json"))
                self.assertEqual(code, 1, text)
                self.assert_nothing_written(self.out)
        p = self.work / "notjson.json"
        p.write_text("{nope", encoding="utf-8")
        self.assertEqual(self.run_main(p)[0], 1)
        self.assert_nothing_written(self.out)

    def test_malformed_model_tags_refused(self):
        tags = self.work / "tags.json"
        for data in ({"a": 1}, [{"entry_id": "x"}], [{"entry_id": "x", "class": "fev1", "quote": "q", "model_id": 3}]):
            tags.write_text(json.dumps(data), encoding="utf-8")
            code, _ = self.run_main(self.labels(labels_file(AGREEMENT)), tags=tags)
            self.assertEqual(code, 1)
            self.assert_nothing_written(self.out)


if __name__ == "__main__":
    unittest.main()
