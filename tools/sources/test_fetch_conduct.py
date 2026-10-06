"""Test for the conduct of fetch_sources.Downloader: an honest user agent, robots.txt checked at run time (AI-agent groups included), one request
per file, no retries, same-host https redirects only, a size cap, pacing. No network: a fake session answers every request.
Usage: python test_fetch_conduct.py        (standard library only; the integration test needs the evidence skill and is skipped without it)"""
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import fetch_sources as fs

UA = "senseiewok-research-evidence/1.0.0 (+https://github.com/senseiewok/cf-skills)"


def stub_robots_allows(body, path, agent):
    """Only enough of robots.txt for these tests: one `*` group and its Disallow lines."""
    disallowed = [ln.split(":", 1)[1].strip() for ln in body.splitlines() if ln.lower().startswith("disallow:")]
    return not any(d and path.startswith(d) for d in disallowed)


STUB = SimpleNamespace(robots_allows=stub_robots_allows, PROJECT_UA=UA)


class Resp:
    def __init__(self, status=200, body=b"", headers=None, text=None):
        self.status_code = status
        self.content = body
        self.text = text if text is not None else body.decode("utf-8", "replace")
        self.headers = headers or {}
        self.closed = False

    def iter_content(self, n):
        for i in range(0, len(self.content), n):
            yield self.content[i:i + n]

    def close(self):
        self.closed = True


class Session:
    def __init__(self, responses):
        self.responses = list(responses)
        self.calls = []

    def get(self, url, headers=None, timeout=None, stream=False, allow_redirects=True):
        self.calls.append(SimpleNamespace(url=url, headers=dict(headers or {}), allow_redirects=allow_redirects, stream=stream))
        item = self.responses.pop(0)
        if isinstance(item, Exception):
            raise item
        return item


class DownloaderTest(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.dest = Path(self._tmp.name) / "f.pdf"
        self.sleeps = []

    def make(self, responses, evidence=STUB, **kw):
        self.session = Session(responses)
        return fs.Downloader(session=self.session, sleep=self.sleeps.append, evidence=evidence, **kw)

    def assertRefused(self, d, url, text=""):
        with self.assertRaises(fs.Refused) as cm:
            d.download(url, self.dest)
        self.assertIn(text, str(cm.exception))
        self.assertFalse(self.dest.exists())
        self.assertFalse(self.dest.with_suffix(".pdf.part").exists())

    # --- identity and the happy path ------------------------------------------------
    def test_every_request_carries_the_project_user_agent_and_the_file_is_written(self):
        d = self.make([Resp(404), Resp(200, b"%PDF-1.7 data")])
        d.download("https://pub.example/r/a.pdf", self.dest)
        self.assertEqual(self.dest.read_bytes(), b"%PDF-1.7 data")
        self.assertEqual([c.url for c in self.session.calls], ["https://pub.example/robots.txt", "https://pub.example/r/a.pdf"])
        self.assertTrue(all(c.headers.get("User-Agent") == UA for c in self.session.calls))
        self.assertFalse(self.dest.with_suffix(".pdf.part").exists())

    def test_redirects_are_never_followed_by_the_library(self):
        d = self.make([Resp(404), Resp(200, b"x")])
        d.download("https://pub.example/a.pdf", self.dest)
        self.assertTrue(all(c.allow_redirects is False for c in self.session.calls))

    # --- robots ----------------------------------------------------------------------
    def test_a_disallowing_robots_file_refuses_and_nothing_else_is_requested(self):
        d = self.make([Resp(200, b"User-agent: *\nDisallow: /r/\n")])
        self.assertRefused(d, "https://pub.example/r/a.pdf", "robots.txt")
        self.assertEqual(len(self.session.calls), 1)

    def test_an_unreadable_robots_file_refuses(self):
        for first in (Resp(403), Resp(401), Resp(500), Resp(429), OSError("down")):
            d = self.make([first])
            self.assertRefused(d, "https://pub.example/a.pdf", "robots.txt")
            self.assertEqual(len(self.session.calls), 1)

    def test_a_missing_robots_file_allows(self):
        for status in (404, 410):
            d = self.make([Resp(status), Resp(200, b"ok")])
            d.download("https://pub.example/a.pdf", self.dest)
            self.dest.unlink()

    def test_a_huge_robots_file_refuses(self):
        d = self.make([Resp(200, b"#" * (fs.MAX_ROBOTS_BYTES + 1))])
        self.assertRefused(d, "https://pub.example/a.pdf", "robots.txt")

    def test_robots_is_fetched_once_per_host(self):
        d = self.make([Resp(404), Resp(200, b"one"), Resp(200, b"two")])
        d.download("https://pub.example/a.pdf", self.dest)
        d.download("https://pub.example/b.pdf", self.dest)
        self.assertEqual([c.url.rsplit("/", 1)[1] for c in self.session.calls], ["robots.txt", "a.pdf", "b.pdf"])

    def test_without_the_evidence_skill_nothing_is_fetched(self):
        def missing():
            raise fs.Refused("the evidence skill that holds the robots rule is not available")
        d = fs.Downloader(session=Session([]), sleep=self.sleeps.append, evidence=None, load_evidence=missing)
        self.assertRefused(d, "https://pub.example/a.pdf", "evidence skill")

    # --- the url and redirects ---------------------------------------------------------
    def test_only_https_urls(self):
        d = self.make([])
        for url in ("http://pub.example/a.pdf", "ftp://pub.example/a.pdf", "file:///etc/passwd", "https:///a.pdf", "pub.example/a.pdf"):
            self.assertRefused(d, url)
        self.assertEqual(self.session.calls, [])

    def test_a_same_host_redirect_is_followed_and_checked_against_robots(self):
        d = self.make([Resp(200, b"User-agent: *\nDisallow: /private/\n"), Resp(302, headers={"Location": "/files/a.pdf"}), Resp(200, b"ok")])
        d.download("https://pub.example/a.pdf", self.dest)
        self.assertEqual(self.session.calls[-1].url, "https://pub.example/files/a.pdf")
        d = self.make([Resp(200, b"User-agent: *\nDisallow: /private/\n"), Resp(301, headers={"Location": "https://pub.example/private/a.pdf"})])
        self.dest.unlink()
        self.assertRefused(d, "https://pub.example/a.pdf", "robots.txt")

    def test_a_redirect_to_another_host_or_to_http_is_refused(self):
        for target in ("https://other.example/a.pdf", "http://pub.example/a.pdf", "//other.example/a.pdf"):
            d = self.make([Resp(404), Resp(302, headers={"Location": target})])
            self.assertRefused(d, "https://pub.example/a.pdf", "redirect")
            self.assertEqual(len(self.session.calls), 2)

    def test_a_redirect_loop_stops_after_a_few_hops(self):
        d = self.make([Resp(404)] + [Resp(302, headers={"Location": "/a.pdf"}) for _ in range(10)])
        self.assertRefused(d, "https://pub.example/a.pdf", "redirect")
        self.assertLessEqual(len(self.session.calls), 1 + fs.MAX_REDIRECTS + 1)

    # --- no retries, no partial files, a size cap ------------------------------------------
    def test_an_error_answer_is_one_request_and_no_retry(self):
        for status in (403, 404, 429, 500, 503):
            d = self.make([Resp(404), Resp(status)])
            with self.assertRaises(fs.Refused):
                d.download("https://pub.example/a.pdf", self.dest)
            self.assertEqual(len(self.session.calls), 2, status)
            self.assertFalse(self.dest.exists())

    def test_a_declared_size_over_the_cap_is_refused_before_the_body_is_read(self):
        big = Resp(200, b"x", headers={"Content-Length": str(fs.MAX_DOWNLOAD_BYTES + 1)})
        d = self.make([Resp(404), big])
        self.assertRefused(d, "https://pub.example/a.pdf", "larger than")
        self.assertTrue(big.closed)

    def test_a_streamed_body_over_the_cap_is_refused_and_leaves_no_file(self):
        d = self.make([Resp(404), Resp(200, b"x" * 5000)], max_bytes=1000)
        self.assertRefused(d, "https://pub.example/a.pdf", "larger than")

    def test_a_transport_error_leaves_no_file(self):
        d = self.make([Resp(404), OSError("reset")])
        with self.assertRaises(OSError):
            d.download("https://pub.example/a.pdf", self.dest)
        self.assertFalse(self.dest.exists())

    # --- pacing --------------------------------------------------------------------------------
    def test_requests_to_one_host_are_paced(self):
        d = self.make([Resp(404), Resp(200, b"a"), Resp(200, b"b")])
        d.download("https://pub.example/a.pdf", self.dest)
        d.download("https://pub.example/b.pdf", self.dest)
        self.assertGreaterEqual(len(self.sleeps), 2)
        self.assertTrue(all(s >= fs.PAUSE_S for s in self.sleeps))

    def test_a_crawl_delay_longer_than_the_pause_is_honoured_and_capped(self):
        d = self.make([Resp(200, b"User-agent: *\nCrawl-delay: 7\n"), Resp(200, b"a")])
        d.download("https://pub.example/a.pdf", self.dest)
        self.assertIn(7.0, self.sleeps)
        self.dest.unlink()
        d = self.make([Resp(200, b"User-agent: *\nCrawl-delay: 9999\n"), Resp(200, b"a")])
        d.download("https://pub.example/a.pdf", self.dest)
        self.assertEqual(max(self.sleeps), fs.MAX_CRAWL_DELAY_S)

    # --- the audit: one robots.txt request per host, never a document -----------------------------------------
    def test_robots_verdict_reports_allowed_and_refused_without_touching_the_document(self):
        d = self.make([Resp(200, b"User-agent: *\nDisallow: /r/\n")])
        self.assertEqual(d.robots_verdict("https://pub.example/ok.pdf"), (True, ""))
        ok, why = d.robots_verdict("https://pub.example/r/a.pdf")
        self.assertFalse(ok)
        self.assertIn("robots.txt", why)
        self.assertEqual([c.url for c in self.session.calls], ["https://pub.example/robots.txt"])   # one request for two verdicts

    def test_robots_verdict_refuses_when_robots_cannot_be_read(self):
        d = self.make([Resp(403)])
        ok, why = d.robots_verdict("https://pub.example/a.pdf")
        self.assertFalse(ok)
        self.assertIn("robots.txt", why)

    def test_audit_lists_fetch_entries_only_and_exits_non_zero_on_a_refusal(self):
        entries = [{"id": "a", "access": "fetch", "url": "https://one.example/r/a.pdf"},
                   {"id": "b", "access": "fetch", "url": "https://two.example/b.pdf"},
                   {"id": "c", "access": "manual", "url": "https://three.example/c.pdf"},
                   {"id": "d", "access": "api", "url": "https://four.example/d"},
                   {"id": "e", "access": "fetch", "url": None}]
        d = self.make([Resp(200, b"User-agent: *\nDisallow: /r/\n"), Resp(404)])
        rows = fs.audit_robots(entries, d)
        self.assertEqual([(r[0], r[1]) for r in rows], [("a", "REFUSED"), ("b", "allowed")])
        self.assertEqual(sorted({c.url for c in self.session.calls}), ["https://one.example/robots.txt", "https://two.example/robots.txt"])

    def test_crawl_delay_parsing(self):
        self.assertEqual(fs.crawl_delay("User-agent: *\nCrawl-delay: 3\n"), 3.0)
        self.assertEqual(fs.crawl_delay("Crawl-delay: abc\n"), 0.0)
        self.assertEqual(fs.crawl_delay("# Crawl-delay: 5\n"), 0.0)
        self.assertEqual(fs.crawl_delay("Crawl-delay: -2\n"), 0.0)


if __name__ == "__main__":
    unittest.main()
