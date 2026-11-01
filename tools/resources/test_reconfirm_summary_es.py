"""Tests for reconfirm_summary_es.py: it re-stamps only summary_es changes and can never write confirmed_by: maintainer."""
import contextlib
import io
import shutil
import sys
import tempfile
import unittest
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
import check_resources  # noqa: E402
import reconfirm_summary_es as rc  # noqa: E402

ROW = {"name": "CF Helpline", "url": "https://www.example.org/help", "summary": "A helpline for families.",
       "who_for": "stated_countries"}


def text(rows):
    out = ["schema_version: 1\n", "resources:\n"]
    for r in rows:
        out.append(f"  - id: {r['id']}\n")
        for k in ("name", "url", "summary", "summary_es", "who_for"):
            if k in r:
                out.append(f'    {k}: "{r[k]}"\n')
        out.append(f"    confirmed_by: {r['confirmed_by']}\n")
        out.append("    confirmed_on: 2026-10-01\n")
        out.append(f"    content_hash: {r['content_hash']}\n\n")
    return "".join(out)


def make(id_, es=None, hash_with_es=False, **changes):
    r = dict(ROW, id=id_, confirmed_by="maintainer")
    r["content_hash"] = rc.h(r, False)
    if es is not None:
        r["summary_es"] = es
        if hash_with_es:
            r["content_hash"] = rc.h(r, True)
    r.update(changes)
    return r


class Reconfirm(unittest.TestCase):
    def setUp(self):
        self.dir = Path(tempfile.mkdtemp(prefix="cf-reconfirm "))
        self.addCleanup(shutil.rmtree, self.dir, ignore_errors=True)
        self.path = self.dir / "resources.yaml"

    def run_rc(self, *args):
        argv = ["--yaml", str(self.path), "--date", "2026-10-31", *args]
        with contextlib.redirect_stdout(io.StringIO()):
            rc.main(argv)
        return {r["id"]: r for r in yaml.safe_load(self.path.read_text(encoding="utf-8"))["resources"]}

    def test_council_value_is_written_with_a_correct_hash(self):
        self.path.write_text(text([make("a", es="Una línea de ayuda.")]), encoding="utf-8")
        r = self.run_rc("--confirmed-by", "council-on-behalf-of-maintainer")["a"]
        self.assertEqual(r["confirmed_by"], "council-on-behalf-of-maintainer")
        self.assertEqual(r["content_hash"], check_resources.compute_hash(r))
        self.assertEqual(str(r["confirmed_on"]), "2026-10-31")

    def test_default_leaves_confirmed_by_unchanged(self):
        self.path.write_text(text([make("a", es="Una línea de ayuda.")]), encoding="utf-8")
        self.assertEqual(self.run_rc()["a"]["confirmed_by"], "maintainer")

    def test_script_cannot_write_maintainer(self):
        self.path.write_text(text([make("a", es="Una línea de ayuda.")]), encoding="utf-8")
        before = self.path.read_bytes()
        for value in ("maintainer", "Maintainer", "model"):
            with self.subTest(value=value), contextlib.redirect_stderr(io.StringIO()):
                with self.assertRaises(SystemExit):
                    rc.main(["--yaml", str(self.path), "--date", "2026-10-31", "--confirmed-by", value])
        self.assertEqual(self.path.read_bytes(), before)
        self.assertNotIn("maintainer", rc.SCRIPT_CONFIRMED_BY)

    def test_refuses_change_outside_summary_es(self):
        r = make("a", es="Una línea de ayuda.")
        r["summary"] = "A different English summary."
        self.path.write_text(text([r]), encoding="utf-8")
        with self.assertRaises(SystemExit):
            self.run_rc("--confirmed-by", "council-on-behalf-of-maintainer")

    def test_changed_summary_es_needs_explicit_id_and_base(self):
        old = make("a", es="Texto anterior.", hash_with_es=True)
        base = self.dir / "base.yaml"
        base.write_text(text([old]), encoding="utf-8")
        new = dict(old, summary_es="Texto corregido.")
        self.path.write_text(text([new]), encoding="utf-8")
        with self.assertRaises(SystemExit):  # not named
            self.run_rc("--confirmed-by", "council-on-behalf-of-maintainer", "--base-yaml", str(base))
        with self.assertRaises(SystemExit):  # named but no base
            self.run_rc("--confirmed-by", "council-on-behalf-of-maintainer", "--ids", "a")
        r = self.run_rc("--confirmed-by", "council-on-behalf-of-maintainer", "--ids", "a", "--base-yaml", str(base))["a"]
        self.assertEqual(r["content_hash"], check_resources.compute_hash(r))
        self.assertEqual(r["confirmed_by"], "council-on-behalf-of-maintainer")

    def test_changed_summary_es_refused_when_base_fields_differ(self):
        old = make("a", es="Texto anterior.", hash_with_es=True)
        base = self.dir / "base.yaml"
        base.write_text(text([old]), encoding="utf-8")
        new = dict(old, summary_es="Texto corregido.", url="https://www.example.org/other")
        self.path.write_text(text([new]), encoding="utf-8")
        with self.assertRaises(SystemExit):
            self.run_rc("--confirmed-by", "council-on-behalf-of-maintainer", "--ids", "a", "--base-yaml", str(base))


if __name__ == "__main__":
    unittest.main()
