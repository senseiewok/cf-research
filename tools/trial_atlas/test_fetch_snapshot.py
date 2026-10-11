"""Tests for fetch_snapshot.py against a fake client serving SYNTHETIC pages (synthetic_fixtures.py). No network: every socket call
fails the test. The cases the packet names: a totalCount mismatch, a missing next token, a refusal by the catalog gate, a missing
contact and a budget exhaustion; plus the shapes that must stop the run. One integration case uses the evidence skill's real Client
with a fake session; it is skipped when the skill (or its `requests` dependency) is not importable.
Usage: python -m unittest -v (inside tools/trial_atlas)."""
import contextlib
import copy
import io
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import fetch_snapshot as fs  # noqa: E402
import snapshot as snap  # noqa: E402
import synthetic_fixtures as sf  # noqa: E402


class Base(unittest.TestCase):
    def setUp(self):
        sf.block_network(self)
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.out = Path(self._tmp.name) / "snap"

    def run_cf(self, client, routes=("condition", "term"), **kw):
        return fs.run(client, self.out, routes=list(routes), cfg={"page_size": 8, **kw.pop("cfg", {})}, project_ua="SYNTHETIC-agent",
                      contact=sf.DUMMY_CONTACT, **kw)

    def assertStopped(self, client, text, exc=fs.SnapshotStopped, **kw):
        with self.assertRaises(exc) as cm:
            self.run_cf(client, **kw)
        self.assertIn(text, str(cm.exception))
        self.assertFalse((self.out / "manifest.json").exists(), "a stopped run must not leave a manifest")
        if self.out.exists():
            self.assertIn("stopped:", (self.out / "INCOMPLETE.txt").read_text(encoding="utf-8"))
        return cm.exception


def client_with(cond_pages=None, term_pages=None, **kw):
    cond = cond_pages if cond_pages is not None else sf.pages(sf.cf_condition_studies(), 8, "condition")
    term = term_pages if term_pages is not None else sf.pages(sf.cf_term_studies(), 8, "term")
    return sf.FakeClient({fs.ROUTES["condition"]["query.cond"]: cond, fs.ROUTES["term"]["query.term"]: term}, **kw)


class HappyPathTest(Base):
    def test_snapshot_is_written_hashed_and_loadable(self):
        m = self.run_cf(client_with())
        self.assertEqual(m["data_timestamp"], sf.VERSION["dataTimestamp"])
        self.assertEqual(m["requests"], 5)                       # version + 2 pages + 2 pages
        self.assertEqual([(r["name"], r["total_count"], r["studies_received"], len(r["pages"])) for r in m["routes"]],
                         [("condition", 16, 16, 2), ("term", 12, 12, 2)])
        self.assertIn("query.cond=cystic+fibrosis", m["routes"][0]["query"])
        self.assertIn("pageToken=SYNTHETIC-condition-2", m["routes"][0]["pages"][1]["query"])
        self.assertFalse(m["config_verified"])
        s = snap.load(self.out)                                  # every hash checks
        self.assertEqual(s.digest, m["snapshot_sha256"])
        self.assertEqual(len(s.studies("condition")), 16)

    def test_the_contact_is_never_written(self):
        self.run_cf(client_with())
        for p in self.out.rglob("*"):
            if p.is_file():
                text = p.read_text(encoding="utf-8")
                self.assertNotIn(sf.DUMMY_CONTACT, text, p.name)
                self.assertNotIn("@", text, p.name)
        self.assertFalse(json.loads((self.out / "manifest.json").read_text(encoding="utf-8"))["contact_recorded"])

    def test_the_first_page_asks_for_the_total_and_later_pages_send_the_token(self):
        c = client_with()
        self.run_cf(c, routes=("condition",))
        studies_calls = [p for path, p in c.calls if path == "studies"]
        self.assertEqual(studies_calls[0].get("countTotal"), "true")
        self.assertNotIn("pageToken", studies_calls[0])
        self.assertEqual(studies_calls[1].get("pageToken"), "SYNTHETIC-condition-2")
        self.assertEqual(studies_calls[1].get("pageSize"), 8)

    def test_a_changed_hash_is_caught_on_load(self):
        self.run_cf(client_with())
        page = self.out / "condition" / "page-0001.json"
        page.write_text(page.read_text(encoding="utf-8").replace("SYNTHETIC modulator", "SYNTHETIC edited"), encoding="utf-8")
        with self.assertRaises(snap.SnapshotError):
            snap.load(self.out)

    def test_an_existing_folder_is_never_overwritten(self):
        self.out.mkdir()
        with self.assertRaises(fs.Refused):
            self.run_cf(client_with())


class FailClosedTest(Base):
    def test_total_count_mismatch(self):
        pages = sf.pages(sf.cf_condition_studies(), 8, "condition")
        pages[0]["totalCount"] = 17
        self.assertStopped(client_with(pages), "received 16 studies but 'totalCount' is 17", routes=("condition",))

    def test_missing_next_token(self):
        pages = sf.pages(sf.cf_condition_studies(), 8, "condition")
        del pages[0]["nextPageToken"]
        self.assertStopped(client_with(pages), "the last page named no 'nextPageToken'", routes=("condition",))

    def test_missing_total_on_the_first_page(self):
        pages = sf.pages(sf.cf_condition_studies(), 8, "condition")
        del pages[0]["totalCount"]
        self.assertStopped(client_with(pages), "no usable 'totalCount'", routes=("condition",))

    def test_missing_studies_list(self):
        pages = sf.pages(sf.cf_condition_studies(), 8, "condition")
        del pages[1]["studies"]
        self.assertStopped(client_with(pages), "has no 'studies' list", routes=("condition",))

    def test_a_study_without_an_nct_id(self):
        pages = sf.pages(sf.cf_condition_studies(), 8, "condition")
        pages[0]["studies"][2]["protocolSection"]["identificationModule"]["nctId"] = "not-an-id"
        self.assertStopped(client_with(pages), "no valid protocolSection.identificationModule.nctId", routes=("condition",))

    def test_a_total_that_changes_between_pages(self):
        pages = sf.pages(sf.cf_condition_studies(), 8, "condition")
        pages[1]["totalCount"] = 18
        self.assertStopped(client_with(pages), "changed from 16 to 18", routes=("condition",))

    def test_a_repeated_token(self):
        pages = sf.pages(sf.cf_condition_studies(), 4, "condition")
        pages[1]["nextPageToken"] = "SYNTHETIC-condition-2"
        self.assertStopped(client_with(pages), "the next token repeats", routes=("condition",))

    def test_a_study_received_twice(self):
        pages = sf.pages(sf.cf_condition_studies(), 8, "condition")
        pages[1]["studies"][0] = copy.deepcopy(pages[0]["studies"][0])
        self.assertStopped(client_with(pages), "was already received on an earlier page", routes=("condition",))

    def test_too_many_pages(self):
        self.assertStopped(client_with(sf.pages(sf.cf_condition_studies(), 2, "condition")), "more than 3 pages",
                           routes=("condition",), cfg={"max_pages": 3})

    def test_version_without_data_timestamp(self):
        self.assertStopped(client_with(version={"apiVersion": "x"}), "has no dataTimestamp")

    def test_budget_exhaustion(self):
        exc = self.assertStopped(client_with(max_requests=3), "request budget exhausted")
        self.assertIn("BudgetExhausted", str(exc))
        self.assertIn("requests made: 3", (self.out / "INCOMPLETE.txt").read_text(encoding="utf-8"))

    def test_a_refusal_by_the_catalog_gate(self):
        c = client_with(allowed=())
        self.assertStopped(c, "refused by the evidence client's gate", exc=fs.Refused)
        self.assertEqual(c.accounting.attempts, 0)
        self.assertEqual(c.calls, [])

    def test_an_error_status_is_not_retried(self):
        c = client_with()
        c.routes[fs.ROUTES["condition"]["query.cond"]] = [{"_synthetic": sf.LABEL}]
        orig = c.get

        def get(source_id, path, params=None):
            f = orig(source_id, path, params)
            return SimpleNamespace(**{**vars(f), "status": "error", "http_status": 500}) if path == "studies" else f
        c.get = get
        self.assertStopped(c, "status error (http 500); not retried", routes=("condition",))
        self.assertEqual(c.accounting.attempts, 2)


class ContactAndSkillTest(Base):
    def test_no_contact_means_no_run(self):
        with mock.patch.dict(os.environ, {"EVIDENCE_CONTACT": ""}):
            with self.assertRaises(fs.Refused) as cm:
                fs.resolve_contact()
            self.assertIn("no contact", str(cm.exception))
            with mock.patch.object(fs, "load_evidence_http", side_effect=AssertionError("must not load the client")):
                err = io.StringIO()
                with contextlib.redirect_stderr(err):
                    self.assertEqual(fs.main(["--out", str(self.out)]), 2)
                self.assertIn("no contact", err.getvalue())
        self.assertFalse(self.out.exists())

    def test_the_contact_can_come_from_the_environment(self):
        with mock.patch.dict(os.environ, {"EVIDENCE_CONTACT": sf.DUMMY_CONTACT}):
            self.assertEqual(fs.resolve_contact(), sf.DUMMY_CONTACT)

    def test_the_contact_is_not_a_command_line_option(self):
        with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit) as cm:
            fs.main(["--contact", "someone", "--out", str(self.out)])
        self.assertEqual(cm.exception.code, 2)
        self.assertFalse(self.out.exists())

    # ---- second fix round (after 0401af9): the reviewer's contact_probe.py cases. Each failed before its fix.

    FAKE = "probe.user" + "@" + "example.invalid"

    def _scenario(self, client_cls):
        ev = SimpleNamespace(PROJECT_UA="SYNTHETIC-UA")

        def fake_make_client(contact, evidence=None, *, dry_run=False):
            c = client_cls(sf.cf_client(8).routes, contact=contact)
            c.accounting.summary = lambda: f"requests: {c.accounting.attempts}"
            return c, ev
        out, err = io.StringIO(), io.StringIO()
        with mock.patch.dict(os.environ, {"EVIDENCE_CONTACT": self.FAKE}), mock.patch.object(fs, "make_client", fake_make_client), \
                contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            code = fs.main(["--out", str(self.out)])
        written = "".join(p.read_text(encoding="utf-8") for p in self.out.rglob("*") if p.is_file()) if self.out.exists() else ""
        return code, out.getvalue(), err.getvalue(), written

    def test_a_gate_error_quoting_the_user_agent_never_shows_the_contact(self):
        class BoomPerm(sf.FakeClient):
            def get(self, *a, **k):
                raise PermissionError("denied for " + self.user_agent)
        code, out, err, written = self._scenario(BoomPerm)
        self.assertEqual(code, 2)
        for text in (out, err, written):
            self.assertNotIn(self.FAKE, text)
            self.assertNotIn("@", text)
        self.assertIn("[contact]", err)
        self.assertIn("[contact]", written)

    def test_a_budget_error_quoting_the_user_agent_never_shows_the_contact(self):
        class BoomBudget(sf.FakeClient):
            def get(self, *a, **k):
                raise type("BudgetExhausted", (RuntimeError,), {})("budget " + self.user_agent)
        code, out, err, written = self._scenario(BoomBudget)
        self.assertEqual(code, 1)
        for text in (out, err, written):
            self.assertNotIn(self.FAKE, text)
        self.assertIn("[contact]", err)

    def test_a_usage_error_never_echoes_an_address(self):
        err = io.StringIO()
        with contextlib.redirect_stderr(err), self.assertRaises(SystemExit):
            fs.main(["--contact", self.FAKE])
        self.assertNotIn(self.FAKE, err.getvalue())

    # ---- third round (after adff5aa). Each case failed before its fix.

    def _main(self, args, client_factory=None, dry=False):
        ev = SimpleNamespace(PROJECT_UA="SYNTHETIC-UA")
        made = []

        def fake_make_client(contact, evidence=None, *, dry_run=False):
            c = (client_factory or (lambda: sf.cf_client(8)))()
            if dry:
                c.dry_run = True
                orig = c.get
                c.get = lambda s, p, params=None: SimpleNamespace(**{**vars(orig(s, p, params)), "dry_run": True, "status": "out_of_scope"})
            c.accounting.summary = lambda: f"requests: {c.accounting.attempts}"
            made.append(c)
            return c, ev
        out, err = io.StringIO(), io.StringIO()
        with mock.patch.dict(os.environ, {"EVIDENCE_CONTACT": sf.DUMMY_CONTACT}), mock.patch.object(fs, "make_client", fake_make_client), \
                contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            try:
                code = fs.main(args)
            except SystemExit as e:
                code = e.code
        return code, out.getvalue(), err.getvalue(), made

    def test_the_absolute_output_path_is_printed(self):
        code, out, _, _ = self._main(["--out", str(self.out), "--page-size", "8"])
        self.assertEqual(code, 0, out)
        self.assertIn(f"output folder: {self.out.resolve()}", out)

    def test_a_dry_run_prints_the_folder_and_the_page_ceiling(self):
        code, out, _, _ = self._main(["--out", str(self.out), "--max-pages", "40", "--dry-run"], dry=True)
        self.assertEqual(code, 0, out)
        self.assertIn(f"output folder (not created in a dry run): {self.out.resolve()}", out)
        self.assertIn("planned page ceiling: 40 per retrieval route x 2 routes + 1 version request = 81 of the 200-request budget", out)

    def test_a_ceiling_over_the_budget_is_refused(self):
        code, _, err, made = self._main(["--out", str(self.out), "--max-pages", "150"])
        self.assertEqual(code, 2)
        self.assertIn("200-request budget", err)
        self.assertEqual(made, [])

    def test_the_fields_can_be_overridden(self):
        code, out, _, made = self._main(["--out", str(self.out), "--page-size", "8", "--fields", "NCTId,BriefTitle"])
        self.assertEqual(code, 0, out)
        self.assertTrue(all(p.get("fields") == "NCTId,BriefTitle" for path, p in made[0].calls if path == "studies"))
        self.assertEqual(snap.load(self.out).manifest["config"]["fields"], "NCTId,BriefTitle")

    def test_a_route_query_fetches_a_control_snapshot(self):
        code, out, err, _ = self._main(["--out", str(self.out), "--page-size", "8", "--route-query",
                                        "condition=query.cond=SYNTHETIC asthma OR COPD control"], client_factory=sf.noncf_client)
        self.assertEqual(code, 0, out + err)
        s = snap.load(self.out)
        self.assertEqual(s.routes(), ["condition"])
        self.assertEqual(len(s.studies("condition")), 8)
        self.assertEqual(s.manifest["routes"][0]["params"]["query.cond"], "SYNTHETIC asthma OR COPD control")

    def test_a_bad_route_query_is_a_usage_error(self):
        for bad in ("condition", "condition=query.cond", "Bad Name=query.cond=x", "condition=pageSize=1000"):
            with self.subTest(bad=bad):
                code, _, err, made = self._main(["--out", str(self.out), "--route-query", bad])
                self.assertEqual(code, 2)
                self.assertEqual(made, [])
                self.assertFalse(self.out.exists())

    def test_scrub_replaces_the_contact_and_any_address(self):
        with mock.patch.dict(os.environ, {"EVIDENCE_CONTACT": "not-an-address-shape"}):
            self.assertEqual(fs.scrub("ua not-an-address-shape and " + self.FAKE), "ua [contact] and [contact]")

    def test_a_rejected_contact_is_never_echoed(self):
        class Ev:
            PROJECT_UA = "x"

            class Client:
                def __init__(self, contact):
                    raise ValueError(f"invalid contact {contact!r}: use a plain email address")
        with self.assertRaises(fs.Refused) as cm:
            fs.make_client("SECRET-LOOKING-VALUE", evidence=Ev)
        self.assertEqual(str(cm.exception), "the contact was refused by the evidence client: not a plain email address")

    def test_an_unexpected_client_error_leaves_incomplete(self):
        c = client_with()
        orig = c.get

        def get(source_id, path, params=None):
            if path == "studies":
                raise TypeError("SYNTHETIC unexpected failure")
            return orig(source_id, path, params)
        c.get = get
        with self.assertRaises(TypeError):
            self.run_cf(c)
        self.assertFalse((self.out / "manifest.json").exists())
        self.assertIn("TypeError", (self.out / "INCOMPLETE.txt").read_text(encoding="utf-8"))

    def test_a_missing_evidence_skill_is_a_clear_refusal(self):
        with mock.patch.object(fs, "_evidence_dirs", return_value=[Path(self._tmp.name) / "nowhere"]):
            with self.assertRaises(fs.Refused) as cm:
                fs.load_evidence_http()
        self.assertIn("evidence skill", str(cm.exception))

    def test_an_invalid_contact_is_refused_by_the_client(self):
        class Ev:
            PROJECT_UA = "x"

            class Client:
                def __init__(self, contact):
                    raise ValueError("invalid contact")
        with self.assertRaises(fs.Refused):
            fs.make_client("not an address", evidence=Ev)

    def test_a_dry_run_writes_nothing(self):
        c = client_with()
        c.dry_run = True
        orig = c.get
        c.get = lambda s, p, params=None: SimpleNamespace(**{**vars(orig(s, p, params)), "dry_run": True, "status": "out_of_scope"})
        with contextlib.redirect_stdout(io.StringIO()) as out:
            self.assertIsNone(self.run_cf(c))
        self.assertIn("planned", out.getvalue())
        self.assertFalse(self.out.exists())


def _real_evidence():
    try:
        with mock.patch.dict(os.environ, {}):
            return fs.load_evidence_http()
    except fs.Refused:
        return None


EV = _real_evidence()


@unittest.skipUnless(EV is not None, "the evidence skill (cf-skills) or its requests dependency is not importable here")
class RealClientTest(Base):
    """The evidence skill's own Client, with a fake session in place of the network and its pacing sleep switched off."""

    class Resp:
        def __init__(self, url, data, status=200):
            self.url, self.status_code, self.headers, self.raw = url, status, {}, None
            self._data = data
            self.text = json.dumps(data)

        def json(self):
            return self._data

    class Session:
        def __init__(self, routes):
            self.headers, self.routes, self.calls = {}, routes, []

        def request(self, method, url, params=None, json=None, timeout=None, allow_redirects=True, stream=False):
            self.calls.append((url, dict(params or {}), allow_redirects))
            if url.endswith("/version"):
                return RealClientTest.Resp(url, sf.VERSION)
            pages = self.routes[params.get("query.cond") or params.get("query.term")]
            tok = params.get("pageToken")
            return RealClientTest.Resp(url, pages[0 if tok is None else int(tok.rsplit("-", 1)[1]) - 1])

    def test_the_real_client_gate_pacing_and_user_agent(self):
        session = self.Session({fs.ROUTES["condition"]["query.cond"]: sf.pages(sf.cf_condition_studies(), 8, "condition"),
                                fs.ROUTES["term"]["query.term"]: sf.pages(sf.cf_term_studies(), 8, "term")})
        with mock.patch.object(EV.time, "sleep"):
            client = EV.Client(contact=sf.DUMMY_CONTACT, session=session)
            m = fs.run(client, self.out, routes=["condition", "term"], cfg={"page_size": 8}, project_ua=EV.PROJECT_UA,
                       contact=sf.DUMMY_CONTACT)
        self.assertEqual(m["requests"], 5)
        self.assertIn("mailto:" + sf.DUMMY_CONTACT, session.headers["User-Agent"])
        self.assertTrue(all(allow is False for _, _, allow in session.calls))
        self.assertTrue(all(url.startswith("https://clinicaltrials.gov/api/v2/") for url, _, _ in session.calls))
        self.assertNotIn(sf.DUMMY_CONTACT, (self.out / "manifest.json").read_text(encoding="utf-8"))

    def test_the_real_catalog_gate_refuses_a_manual_source(self):
        entries = EV.catalog.load()
        session = self.Session({})
        with mock.patch.dict(entries, {fs.SOURCE_ID: {**entries[fs.SOURCE_ID], "access": "manual"}}):
            client = EV.Client(contact=sf.DUMMY_CONTACT, session=session)
            with self.assertRaises(fs.Refused) as cm:
                fs.run(client, self.out, routes=["condition"], contact=sf.DUMMY_CONTACT)
        self.assertIn("AccessDenied", str(cm.exception))
        self.assertEqual(session.calls, [])

    def test_the_real_client_budget(self):
        session = self.Session({fs.ROUTES["condition"]["query.cond"]: sf.pages(sf.cf_condition_studies(), 2, "condition")})
        with mock.patch.object(EV.time, "sleep"), mock.patch.dict(os.environ, {"EVIDENCE_MAX_REQUESTS": "3"}):
            client = EV.Client(contact=sf.DUMMY_CONTACT, session=session)
            with self.assertRaises(fs.SnapshotStopped):
                fs.run(client, self.out, routes=["condition"], cfg={"page_size": 2}, contact=sf.DUMMY_CONTACT)
        self.assertEqual(len(session.calls), 3)


if __name__ == "__main__":
    unittest.main()
