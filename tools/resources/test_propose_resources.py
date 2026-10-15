"""Offline tests for propose_resources.py: conduct, robots, redirects, size cap, data boundary, and lead_only rows.

Standard library + pytest; no real network sockets are opened during tests.
"""

from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import pytest
import yaml

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import check_resources  # noqa: E402
import propose_resources as pr  # noqa: E402


class FakeResponse:
    def __init__(self, status_code=200, content=b"", headers=None, text=None):
        self.status_code = status_code
        self.content = content
        self.text = text if text is not None else content.decode("utf-8", errors="replace")
        self.headers = headers or {}
        self.closed = False

    def iter_content(self, chunk_size=16384):
        for i in range(0, len(self.content), chunk_size):
            yield self.content[i : i + chunk_size]

    def close(self):
        self.closed = True


class FakeSession:
    def __init__(self, responses_by_url: dict[str, list[FakeResponse | Exception]] | None = None):
        self.responses_by_url = {k.lower(): list(v) for k, v in (responses_by_url or {}).items()}
        self.calls: list[SimpleNamespace] = []

    def get(self, url, headers=None, timeout=None, stream=False, allow_redirects=True):
        self.calls.append(
            SimpleNamespace(
                url=url,
                headers=dict(headers or {}),
                timeout=timeout,
                stream=stream,
                allow_redirects=allow_redirects,
            )
        )
        url_lower = url.lower()
        if url_lower in self.responses_by_url and self.responses_by_url[url_lower]:
            item = self.responses_by_url[url_lower].pop(0)
            if isinstance(item, Exception):
                raise item
            return item
        return FakeResponse(404, b"Not Found")


class ProposeResourcesConductTest(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.scratch_dir = Path(self._tmp.name) / "scratch"

        # Mock catalog
        self.catalog_path = Path(self._tmp.name) / "catalog.yaml"
        self.catalog_data = {
            "sources": [
                {
                    "id": "cftrust-resource-pages",
                    "kind": "resource_pages",
                    "access": "fetch",
                    "landing_page": "https://www.cysticfibrosis.org.uk/",
                    "publisher": "Cystic Fibrosis Trust",
                },
                {
                    "id": "cff-manual-pages",
                    "kind": "resource_pages",
                    "access": "manual",
                    "landing_page": "https://www.cff.org/support",
                    "publisher": "CFF",
                },
            ]
        }
        self.catalog_path.write_text(yaml.dump(self.catalog_data), encoding="utf-8")

    # --- 1. User Agent & Identity ----------------------------------------------------
    def test_user_agent_present_on_every_request(self):
        fake_sess = FakeSession({
            "https://www.cysticfibrosis.org.uk/robots.txt": [FakeResponse(200, b"User-agent: *\nAllow: /\n")],
            "https://www.cysticfibrosis.org.uk/support": [
                FakeResponse(200, b"<html><title>Support</title><body>Help for CF</body></html>")
            ],
        })

        def mock_worker(prompt):
            return json.dumps({
                "name": "CF Trust Support",
                "summary": "Practical information and emotional support for people living with cystic fibrosis.",
                "audience": ["adult_with_cf"],
                "category": ["daily_life"],
                "who_for": "stated_anyone",
                "format": "online",
                "countries": ["GB"],
                "lang": "en",
                "organisation": "Cystic Fibrosis Trust",
            })

        rows = pr.propose_resources(
            urls=["https://www.cysticfibrosis.org.uk/support"],
            source_id="cftrust-resource-pages",
            catalog_path=self.catalog_path,
            scratch_dir=self.scratch_dir,
            pause_s=0.0,
            session=fake_sess,
            worker_fn=mock_worker,
        )
        self.assertEqual(len(rows), 1)
        self.assertTrue(len(fake_sess.calls) >= 2)
        for call in fake_sess.calls:
            self.assertEqual(call.headers.get("User-Agent"), pr.PROJECT_UA)

    # --- 2. Robots.txt conduct -------------------------------------------------------
    def test_robots_txt_disallow_stops_fetch(self):
        fake_sess = FakeSession({
            "https://www.cysticfibrosis.org.uk/robots.txt": [
                FakeResponse(200, b"User-agent: *\nDisallow: /support\n")
            ],
            "https://www.cysticfibrosis.org.uk/support": [
                FakeResponse(200, b"<html><body>Never reached</body></html>")
            ],
        })
        with self.assertRaises(pr.Refused) as cm:
            pr.propose_resources(
                urls=["https://www.cysticfibrosis.org.uk/support"],
                source_id="cftrust-resource-pages",
                catalog_path=self.catalog_path,
                scratch_dir=self.scratch_dir,
                pause_s=0.0,
                session=fake_sess,
                worker_fn=lambda p: "",
            )
        self.assertIn("robots.txt", str(cm.exception))
        # Ensure page itself was never requested!
        called_urls = [c.url for c in fake_sess.calls]
        self.assertNotIn("https://www.cysticfibrosis.org.uk/support", called_urls)

    def test_robots_txt_unreadable_stops_run(self):
        # 401, 403, 500, network error
        for status in (401, 403, 500, 503):
            fake_sess = FakeSession({
                "https://www.cysticfibrosis.org.uk/robots.txt": [FakeResponse(status, b"Error")]
            })
            with self.assertRaises(pr.Refused) as cm:
                pr.propose_resources(
                    urls=["https://www.cysticfibrosis.org.uk/support"],
                    source_id="cftrust-resource-pages",
                    catalog_path=self.catalog_path,
                    scratch_dir=self.scratch_dir,
                    pause_s=0.0,
                    session=fake_sess,
                    worker_fn=lambda p: "",
                )
            self.assertIn("unreadable", str(cm.exception))

    def test_robots_txt_404_allows_fetch(self):
        fake_sess = FakeSession({
            "https://www.cysticfibrosis.org.uk/robots.txt": [FakeResponse(404, b"Not found")],
            "https://www.cysticfibrosis.org.uk/support": [
                FakeResponse(200, b"<html><title>Support</title><body>Help page</body></html>")
            ],
        })
        rows = pr.propose_resources(
            urls=["https://www.cysticfibrosis.org.uk/support"],
            source_id="cftrust-resource-pages",
            catalog_path=self.catalog_path,
            scratch_dir=self.scratch_dir,
            pause_s=0.0,
            session=fake_sess,
            dry_run=True,
        )
        self.assertEqual(len(rows), 1)

    def test_ai_agent_group_disallow_in_robots_stops_fetch(self):
        fake_sess = FakeSession({
            "https://www.cysticfibrosis.org.uk/robots.txt": [
                FakeResponse(200, b"User-agent: *\nAllow: /\n\nUser-agent: anthropic-ai\nDisallow: /support\n")
            ],
        })
        with self.assertRaises(pr.Refused) as cm:
            pr.propose_resources(
                urls=["https://www.cysticfibrosis.org.uk/support"],
                source_id="cftrust-resource-pages",
                catalog_path=self.catalog_path,
                scratch_dir=self.scratch_dir,
                pause_s=0.0,
                session=fake_sess,
                worker_fn=lambda p: "",
            )
        self.assertIn("robots.txt group", str(cm.exception))

    # --- 3. Redirect Handling --------------------------------------------------------
    def test_same_host_https_redirect_allowed(self):
        fake_sess = FakeSession({
            "https://www.cysticfibrosis.org.uk/robots.txt": [FakeResponse(200, b"User-agent: *\nAllow: /\n")],
            "https://www.cysticfibrosis.org.uk/support": [
                FakeResponse(301, headers={"Location": "https://www.cysticfibrosis.org.uk/support/"})
            ],
            "https://www.cysticfibrosis.org.uk/support/": [
                FakeResponse(200, b"<html><title>Support</title><body>Content here</body></html>")
            ],
        })
        rows = pr.propose_resources(
            urls=["https://www.cysticfibrosis.org.uk/support"],
            source_id="cftrust-resource-pages",
            catalog_path=self.catalog_path,
            scratch_dir=self.scratch_dir,
            pause_s=0.0,
            session=fake_sess,
            dry_run=True,
        )
        self.assertEqual(len(rows), 1)

    def test_pacing_invoked_on_redirect_hops(self):
        fake_sess = FakeSession({
            "https://www.cysticfibrosis.org.uk/robots.txt": [FakeResponse(200, b"User-agent: *\nAllow: /\n")],
            "https://www.cysticfibrosis.org.uk/support": [
                FakeResponse(301, headers={"Location": "https://www.cysticfibrosis.org.uk/support/"})
            ],
            "https://www.cysticfibrosis.org.uk/support/": [
                FakeResponse(200, b"<html><title>Support</title><body>Content here</body></html>")
            ],
        })
        sess = pr.IntakeSession(host="www.cysticfibrosis.org.uk", pause_s=1.0, session=fake_sess)
        with patch.object(sess, "_pace", wraps=sess._pace) as mock_pace:
            sess.fetch_page("https://www.cysticfibrosis.org.uk/support")
            # _pace must be called for robots.txt, first request, and redirected request (3 times total)
            self.assertEqual(mock_pace.call_count, 3)

    def test_cross_host_redirect_refused(self):
        fake_sess = FakeSession({
            "https://www.cysticfibrosis.org.uk/robots.txt": [FakeResponse(200, b"User-agent: *\nAllow: /\n")],
            "https://www.cysticfibrosis.org.uk/support": [
                FakeResponse(302, headers={"Location": "https://external.tracker.com/support"})
            ],
        })
        with self.assertRaises(pr.Refused) as cm:
            pr.propose_resources(
                urls=["https://www.cysticfibrosis.org.uk/support"],
                source_id="cftrust-resource-pages",
                catalog_path=self.catalog_path,
                scratch_dir=self.scratch_dir,
                pause_s=0.0,
                session=fake_sess,
                worker_fn=lambda p: "",
            )
        self.assertIn("cross-host redirect", str(cm.exception))

    def test_insecure_http_redirect_refused(self):
        fake_sess = FakeSession({
            "https://www.cysticfibrosis.org.uk/robots.txt": [FakeResponse(200, b"User-agent: *\nAllow: /\n")],
            "https://www.cysticfibrosis.org.uk/support": [
                FakeResponse(301, headers={"Location": "http://www.cysticfibrosis.org.uk/support"})
            ],
        })
        with self.assertRaises(pr.Refused) as cm:
            pr.propose_resources(
                urls=["https://www.cysticfibrosis.org.uk/support"],
                source_id="cftrust-resource-pages",
                catalog_path=self.catalog_path,
                scratch_dir=self.scratch_dir,
                pause_s=0.0,
                session=fake_sess,
                worker_fn=lambda p: "",
            )
        self.assertIn("insecure redirect", str(cm.exception))

    # --- 4. Size Cap -----------------------------------------------------------------
    def test_size_cap_content_length_header_refused(self):
        fake_sess = FakeSession({
            "https://www.cysticfibrosis.org.uk/robots.txt": [FakeResponse(200, b"User-agent: *\nAllow: /\n")],
            "https://www.cysticfibrosis.org.uk/support": [
                FakeResponse(200, b"", headers={"Content-Length": "99999999"})
            ],
        })
        with self.assertRaises(pr.Refused) as cm:
            pr.propose_resources(
                urls=["https://www.cysticfibrosis.org.uk/support"],
                source_id="cftrust-resource-pages",
                catalog_path=self.catalog_path,
                scratch_dir=self.scratch_dir,
                pause_s=0.0,
                max_bytes=1000,
                session=fake_sess,
                worker_fn=lambda p: "",
            )
        self.assertIn("exceeds cap", str(cm.exception))

    def test_size_cap_stream_bytes_refused(self):
        large_body = b"A" * 2000
        fake_sess = FakeSession({
            "https://www.cysticfibrosis.org.uk/robots.txt": [FakeResponse(200, b"User-agent: *\nAllow: /\n")],
            "https://www.cysticfibrosis.org.uk/support": [
                FakeResponse(200, large_body)  # no Content-Length header
            ],
        })
        with self.assertRaises(pr.Refused) as cm:
            pr.propose_resources(
                urls=["https://www.cysticfibrosis.org.uk/support"],
                source_id="cftrust-resource-pages",
                catalog_path=self.catalog_path,
                scratch_dir=self.scratch_dir,
                pause_s=0.0,
                max_bytes=1000,
                session=fake_sess,
                worker_fn=lambda p: "",
            )
        self.assertIn("exceeded size cap", str(cm.exception))

    # --- 5. Catalog Access Gate ------------------------------------------------------
    def test_catalog_access_manual_refuses_immediately(self):
        with self.assertRaises(pr.Refused) as cm:
            pr.propose_resources(
                urls=["https://www.cff.org/support"],
                source_id="cff-manual-pages",
                catalog_path=self.catalog_path,
                scratch_dir=self.scratch_dir,
            )
        self.assertIn("access 'manual'", str(cm.exception))

    def test_uncatalogued_host_refuses(self):
        with self.assertRaises(pr.Refused) as cm:
            pr.propose_resources(
                urls=["https://www.unknown-site.org/support"],
                host="www.unknown-site.org",
                catalog_path=self.catalog_path,
                scratch_dir=self.scratch_dir,
            )
        self.assertIn("no catalog entry", str(cm.exception))

    # --- 6. Data Boundary & Worker Prompting -----------------------------------------
    def test_data_boundary_structure_no_tools(self):
        prompt = pr.build_worker_prompt(
            "https://www.cysticfibrosis.org.uk/help",
            "Help Services",
            "We offer grants and confidential support for people with CF.",
        )
        self.assertIn("=== BEGIN DATA BOUNDARY ===", prompt)
        self.assertIn("=== END DATA BOUNDARY ===", prompt)
        self.assertIn("NO tools and NO network access", prompt)
        self.assertIn("Extract information ONLY from the provided text excerpt", prompt)

    # --- 7. Checker Rules & Dropping Invalid Drafts ----------------------------------
    def test_drop_draft_word_cap_exceeded(self):
        fake_sess = FakeSession({
            "https://www.cysticfibrosis.org.uk/robots.txt": [FakeResponse(200, b"User-agent: *\nAllow: /\n")],
            "https://www.cysticfibrosis.org.uk/support": [
                FakeResponse(200, b"<html><title>Support</title><body>Help</body></html>")
            ],
        })
        long_summary = " ".join(["word"] * 26)  # 26 words > 25 word cap

        def mock_worker(prompt):
            return json.dumps({
                "name": "Service Name",
                "summary": long_summary,
                "audience": ["adult_with_cf"],
                "category": ["daily_life"],
                "who_for": "stated_anyone",
                "countries": ["GB"],
            })

        rows = pr.propose_resources(
            urls=["https://www.cysticfibrosis.org.uk/support"],
            source_id="cftrust-resource-pages",
            catalog_path=self.catalog_path,
            scratch_dir=self.scratch_dir,
            pause_s=0.0,
            session=fake_sess,
            worker_fn=mock_worker,
        )
        # Invalid draft must be dropped!
        self.assertEqual(len(rows), 0)

    def test_drop_draft_with_digits(self):
        fake_sess = FakeSession({
            "https://www.cysticfibrosis.org.uk/robots.txt": [FakeResponse(200, b"User-agent: *\nAllow: /\n")],
            "https://www.cysticfibrosis.org.uk/support": [
                FakeResponse(200, b"<html><title>Support</title><body>Help</body></html>")
            ],
        })

        def mock_worker(prompt):
            return json.dumps({
                "name": "Service Name",
                "summary": "We support 500 patients every month with grant funding.",  # 500 is digit
                "audience": ["adult_with_cf"],
                "category": ["daily_life"],
                "who_for": "stated_anyone",
                "countries": ["GB"],
            })

        rows = pr.propose_resources(
            urls=["https://www.cysticfibrosis.org.uk/support"],
            source_id="cftrust-resource-pages",
            catalog_path=self.catalog_path,
            scratch_dir=self.scratch_dir,
            pause_s=0.0,
            session=fake_sess,
            worker_fn=mock_worker,
        )
        self.assertEqual(len(rows), 0)

    def test_drop_draft_with_contact_details(self):
        fake_sess = FakeSession({
            "https://www.cysticfibrosis.org.uk/robots.txt": [FakeResponse(200, b"User-agent: *\nAllow: /\n")],
            "https://www.cysticfibrosis.org.uk/support": [
                FakeResponse(200, b"<html><title>Support</title><body>Help</body></html>")
            ],
        })

        def mock_worker(prompt):
            return json.dumps({
                "name": "Service Name",
                "summary": "Contact helpline at support@example.org for direct help.",  # email address
                "audience": ["adult_with_cf"],
                "category": ["daily_life"],
                "who_for": "stated_anyone",
                "countries": ["GB"],
            })

        rows = pr.propose_resources(
            urls=["https://www.cysticfibrosis.org.uk/support"],
            source_id="cftrust-resource-pages",
            catalog_path=self.catalog_path,
            scratch_dir=self.scratch_dir,
            pause_s=0.0,
            session=fake_sess,
            worker_fn=mock_worker,
        )
        self.assertEqual(len(rows), 0)

    # --- 8. Writing lead_only Rows & confirmed_by Invariant --------------------------
    def test_tool_writes_lead_only_and_never_writes_confirmed_by(self):
        draft = {
            "name": "CF Trust Helpline",
            "summary": "Confidential support and information provided for anyone affected by cystic fibrosis.",
            "audience": ["adult_with_cf", "parent_carer_of_child"],
            "category": ["mental_health"],
            "who_for": "stated_anyone",
            "format": "online",
            "countries": ["GB"],
            "lang": "en",
            "organisation": "Cystic Fibrosis Trust",
        }
        row = pr.build_lead_only_row(draft, "https://www.cysticfibrosis.org.uk/helpline", "cftrust-resource-pages")
        self.assertEqual(row["basis"], "lead_only")
        self.assertNotIn("confirmed_by", row)
        self.assertNotIn("confirmed_on", row)
        self.assertNotIn("content_hash", row)
        self.assertNotIn("summary_attested_by", row)

        # Invariant: raising ConductViolation if worker/caller tries to inject confirmed_by
        draft_with_confirmation = dict(draft, confirmed_by="maintainer")
        with self.assertRaises(pr.ConductViolation):
            pr.build_lead_only_row(
                draft_with_confirmation,
                "https://www.cysticfibrosis.org.uk/helpline",
                "cftrust-resource-pages",
            )

    # --- 9. Scratch Storage & Cleanup ------------------------------------------------
    def test_scratch_files_saved_and_cleaned(self):
        txt_path, meta_path = pr.save_to_scratch(
            self.scratch_dir,
            "test-service",
            "Page text content",
            {"url": "https://example.org"},
        )
        self.assertTrue(txt_path.exists())
        self.assertTrue(meta_path.exists())

        # Clean scratch
        count = pr.clean_scratch(self.scratch_dir)
        self.assertEqual(count, 2)
        self.assertFalse(txt_path.exists())
        self.assertFalse(meta_path.exists())

    # --- 10. Dry Run on One Real Host ------------------------------------------------
    def test_dry_run_real_host_description(self):
        # Dry run on cftrust-resource-pages
        rows = pr.propose_resources(
            urls=["https://www.cysticfibrosis.org.uk/what-is-cystic-fibrosis/faqs"],
            source_id="cftrust-resource-pages",
            catalog_path=pr.CATALOG_PATH,
            scratch_dir=self.scratch_dir,
            dry_run=True,
        )
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["basis"], "lead_only")
        self.assertEqual(rows[0]["source_id"], "cftrust-resource-pages")
        self.assertNotIn("confirmed_by", rows[0])


if __name__ == "__main__":
    unittest.main()
