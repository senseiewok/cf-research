"""Tests for check_leads.py and render_leads.py. Standard library only; no network.

Run: python -m unittest discover -s tools/leads -v
"""
import contextlib
import copy
import datetime
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import check_leads  # noqa: E402
import render_leads  # noqa: E402

TODAY = datetime.date(2026, 10, 10)
CATALOG = """schema_version: 1
sources:

  - id: clinicaltrials-gov
    title: ClinicalTrials.gov API v2
    access: api

  - id: europe-pmc
    title: Europe PMC REST
    access: api

  - id: ecfspr-adr-2024
    title: ECFSPR Annual Data Report 2024
    access: fetch   # a trailing comment
"""


def good_lead(lid="L-0001", **over):
    lead = {
        "id": lid,
        "title": "A lead",
        "question": "Is this a question?",
        "why_it_matters": "It matters to a CF researcher.",
        "status": "idea",
        "sources": [{"type": "council", "ref": "a council review", "tier": "T0", "read_on": None, "note": "a short note"}],
        "data_needed": {"catalog_ids": ["clinicaltrials-gov"], "access": "api", "already_have": True, "fields_missing": []},
        "effort": {"size": "S", "estimate": True},
        "guard": {"ladder": ["computed"], "note": "Counts are computed."},
        "must_not_claim": ["Anything causal."],
        "needs_person": ["none"],
        "links": {"board": [], "proposals": [], "prs": []},
        "created": "2026-10-10",
        "updated": "2026-10-10",
        "created_by": {"kind": "model", "model_id": "test-model"},
        "visibility": "public",
    }
    lead.update(over)
    return lead


def doc(*leads):
    return {"schema_version": 1, "leads": list(leads) or [good_lead()]}


class Fixture(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)
        self.catalog_path = self.dir / "catalog.yaml"
        self.catalog_path.write_text(CATALOG, encoding="utf-8")
        self.catalog = check_leads.load_catalog_access(self.catalog_path)

    def tearDown(self):
        self.tmp.cleanup()

    def errors(self, lead):
        return check_leads.check_lead(lead, self.catalog, TODAY)

    def assertFails(self, lead, fragment):
        errs = self.errors(lead)
        self.assertTrue(any(fragment in e for e in errs), f"expected {fragment!r} in {errs}")

    def run_main(self, data, extra=()):
        path = self.dir / "leads.json"
        path.write_text(json.dumps(data), encoding="utf-8")
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            code = check_leads.main(["--leads", str(path), "--catalog", str(self.catalog_path), "--today", "2026-10-10", *extra])
        return code, buf.getvalue()


class CatalogReading(Fixture):
    def test_reads_ids_and_access(self):
        self.assertEqual(self.catalog, {"clinicaltrials-gov": "api", "europe-pmc": "api", "ecfspr-adr-2024": "fetch"})


class PassingRegister(Fixture):
    def test_good_lead_passes(self):
        self.assertEqual(self.errors(good_lead()), [])

    def test_main_exit_0(self):
        code, out = self.run_main(doc(good_lead(), good_lead("L-0002")))
        self.assertEqual(code, 0, out)
        self.assertIn("2 leads, 0 errors", out)

    def test_access_mapping_for_mixed_ids(self):
        lead = good_lead(data_needed={"catalog_ids": ["clinicaltrials-gov", "ecfspr-adr-2024"],
                                      "access": {"clinicaltrials-gov": "api", "ecfspr-adr-2024": "fetch"},
                                      "already_have": False, "fields_missing": ["x"]})
        self.assertEqual(self.errors(lead), [])

    def test_not_in_catalog(self):
        lead = good_lead(data_needed={"catalog_ids": [], "access": "not in catalog", "already_have": False,
                                      "fields_missing": [], "not_in_catalog": ["AACT"]})
        self.assertEqual(self.errors(lead), [])

    def test_done_with_pr_link(self):
        lead = good_lead(status="done", links={"board": [], "proposals": [], "prs": ["cf-research#39"]})
        self.assertEqual(self.errors(lead), [])

    def test_short_quote_and_apostrophes_pass(self):
        lead = good_lead(question="Does the sponsors' record say \"may be missing information\" and 'a short phrase'?")
        self.assertEqual(self.errors(lead), [])

    def test_shipped_register_passes(self):
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            code = check_leads.main([])
        self.assertEqual(code, 0, buf.getvalue())


class FailingRules(Fixture):
    def test_missing_field(self):
        lead = good_lead()
        del lead["guard"]
        self.assertFails(lead, "missing required field 'guard'")

    def test_unknown_field(self):
        self.assertFails(good_lead(colour="red"), "unknown field 'colour'")

    def test_bad_id(self):
        self.assertFails(good_lead("L-1"), "id must look like L-0001")

    def test_duplicate_id(self):
        code, out = self.run_main(doc(good_lead(), good_lead()))
        self.assertEqual(code, 1)
        self.assertIn("duplicate id", out)

    def test_bad_status(self):
        self.assertFails(good_lead(status="maybe"), "status 'maybe'")

    def test_bad_source_type(self):
        lead = good_lead()
        lead["sources"][0]["type"] = "blog"
        self.assertFails(lead, "type 'blog'")

    def test_bad_tier(self):
        lead = good_lead()
        lead["sources"][0]["tier"] = "T4"
        self.assertFails(lead, "tier 'T4'")

    def test_t1_needs_read_on(self):
        lead = good_lead()
        lead["sources"][0]["tier"] = "T1"
        self.assertFails(lead, "a T1 source needs read_on")

    def test_empty_sources(self):
        self.assertFails(good_lead(sources=[]), "sources must be a non-empty list")

    def test_bad_dates(self):
        self.assertFails(good_lead(created="10/10/2026"), "created must be a YYYY-MM-DD date")
        self.assertFails(good_lead(updated="2026-10-11"), "after today")
        self.assertFails(good_lead(created="2026-10-09", updated="2026-10-08"), "updated is before created")
        lead = good_lead()
        lead["sources"][0]["read_on"] = "2026-13-01"
        self.assertFails(lead, "read_on must be a YYYY-MM-DD date")

    def test_unknown_catalog_id(self):
        lead = good_lead(data_needed={"catalog_ids": ["aact"], "access": "api", "already_have": False, "fields_missing": []})
        self.assertFails(lead, "catalog id 'aact' is not in sources/catalog.yaml")

    def test_access_differs_from_catalog(self):
        lead = good_lead(data_needed={"catalog_ids": ["ecfspr-adr-2024"], "access": "api", "already_have": False,
                                      "fields_missing": []})
        self.assertFails(lead, "differs from the catalog's 'fetch'")

    def test_access_mapping_wrong(self):
        lead = good_lead(data_needed={"catalog_ids": ["clinicaltrials-gov", "ecfspr-adr-2024"],
                                      "access": {"clinicaltrials-gov": "api", "ecfspr-adr-2024": "api"},
                                      "already_have": False, "fields_missing": []})
        self.assertFails(lead, "differs from the catalog's 'fetch'")

    def test_empty_ids_need_not_in_catalog(self):
        lead = good_lead(data_needed={"catalog_ids": [], "access": "api", "already_have": False, "fields_missing": []})
        self.assertFails(lead, "must be 'not in catalog'")

    def test_done_needs_link(self):
        self.assertFails(good_lead(status="done"), "a done lead must link")
        self.assertFails(good_lead(status="in-progress"), "a in-progress lead must link")

    def test_bad_link_forms(self):
        self.assertFails(good_lead(links={"board": ["T-12"], "proposals": [], "prs": []}), "links.board")
        self.assertFails(good_lead(links={"board": [], "proposals": [], "prs": ["39"]}), "links.prs")

    def test_empty_guard_and_must_not_claim(self):
        self.assertFails(good_lead(must_not_claim=[]), "must_not_claim must be a non-empty list")
        self.assertFails(good_lead(guard={"ladder": [], "note": "x"}), "guard.ladder must be a non-empty list")
        self.assertFails(good_lead(guard={"ladder": ["vibes"], "note": "x"}), "unknown value 'vibes'")
        self.assertFails(good_lead(guard={"ladder": ["computed"], "note": ""}), "guard.note")

    def test_effort(self):
        self.assertFails(good_lead(effort={"size": "XL", "estimate": True}), "effort.size")
        self.assertFails(good_lead(effort={"size": "S", "estimate": False}), "effort.estimate must be true")

    def test_needs_person(self):
        self.assertFails(good_lead(needs_person=["statistician"]), "unknown value 'statistician'")
        self.assertFails(good_lead(needs_person=["none", "clinician-reader"]), "'none' stands alone")

    def test_created_by(self):
        self.assertFails(good_lead(created_by={"kind": "model"}), "model_id is required")
        self.assertFails(good_lead(created_by={"kind": "robot"}), "kind must be human or model")

    def test_visibility(self):
        self.assertFails(good_lead(visibility="private"), "visibility must be public")

    def test_conference_private_carries_no_detail(self):
        src = {"type": "conference-private", "ref": "a named poster", "tier": "T0", "read_on": None, "note": "detail"}
        lead = good_lead(sources=[src])
        self.assertFails(lead, "has ref 'private'")
        self.assertFails(lead, "carries no note")

    def test_email(self):
        address = "someone" + "@" + "example.org"  # built from pieces so this file holds no address
        self.assertFails(good_lead(title=f"write to {address}"), "email-like text")

    def test_windows_user_path(self):
        # built from pieces so this file itself holds no user folder path
        folder = "Us" + "ers"
        self.assertFails(good_lead(title=f"see C:\\{folder}\\someone\\notes"), "a user folder path")
        self.assertFails(good_lead(title=f"see C:/{folder}/someone/notes"), "a user folder path")
        self.assertEqual(self.errors(good_lead(title="see https://example.org/users/list")), [])

    def test_nct_id(self):
        fake_id = "NC" + "T" + "01234567"  # built from pieces so this file holds no trial id
        self.assertFails(good_lead(question=f"What about {fake_id}?"), "not on the allowlist")

    def test_long_quote(self):
        long = " ".join(["word"] * 26)
        self.assertFails(good_lead(why_it_matters=f"It says \"{long}\"."), "a quotation of 26 words")
        self.assertFails(good_lead(why_it_matters=f"It says '{long}'."), "a quotation of 26 words")
        lead = good_lead()
        lead["sources"][0]["note"] = f"\u201c{long}\u201d"
        self.assertFails(lead, "a quotation of 26 words")

    def test_top_level_and_usage(self):
        code, out = self.run_main({"schema_version": 2, "leads": [good_lead()], "extra": 1})
        self.assertEqual(code, 1)
        self.assertIn("schema_version must be 1", out)
        self.assertIn("unknown key 'extra'", out)
        code, out = self.run_main({"not_leads": []})
        self.assertEqual(code, 2)
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(io.StringIO()):
            self.assertEqual(check_leads.main(["--leads", str(self.dir / "missing.json")]), 2)
            self.assertEqual(check_leads.main(["--no-such-flag"]), 2)


class Rendering(Fixture):
    def leads(self):
        a = good_lead("L-0002", status="done", links={"board": ["T-0130"], "proposals": [], "prs": []})
        b = good_lead("L-0001")
        c = good_lead("L-0003", title="Pipe | in title")
        return [a, b, c]

    def test_deterministic_and_order_independent(self):
        first = render_leads.render(self.leads())
        again = render_leads.render(list(reversed(copy.deepcopy(self.leads()))))
        self.assertEqual(first, again)

    def test_grouped_by_status_then_id(self):
        text = render_leads.render(self.leads())
        self.assertLess(text.index("### L-0001"), text.index("### L-0003"))
        self.assertLess(text.index("## Idea"), text.index("## Done"))
        self.assertLess(text.index("### L-0003"), text.index("### L-0002"))
        self.assertIn("Pipe \\| in title", text)
        self.assertIn("S (estimate)", text)

    def test_check_mode(self):
        leads_path, out_path = self.dir / "leads.json", self.dir / "LEADS.md"
        leads_path.write_text(json.dumps(doc(*self.leads())), encoding="utf-8")
        args = ["--leads", str(leads_path), "--out", str(out_path)]
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(render_leads.main(args + ["--check"]), 1)  # missing
            self.assertEqual(render_leads.main(args), 0)
            self.assertEqual(render_leads.main(args + ["--check"]), 0)
            out_path.write_bytes(out_path.read_bytes().replace(b"\n", b"\r\n"))
            self.assertEqual(render_leads.main(args + ["--check"]), 0)  # CRLF checkout still matches
            out_path.write_text("edited by hand\n", encoding="utf-8")
            self.assertEqual(render_leads.main(args + ["--check"]), 1)
            self.assertEqual(render_leads.main(["--leads", str(self.dir / "missing.json")]), 2)

    def test_shipped_leads_md_is_current(self):
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(render_leads.main(["--check"]), 0)


if __name__ == "__main__":
    unittest.main()
