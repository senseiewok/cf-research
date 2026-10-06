"""Unit tests for fetch_sources.py. Standard library only; no network, no real catalog.

Run from this folder:  python -m unittest -v
"""

from __future__ import annotations

import contextlib
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import fetch_sources as fs

CATALOG = """\
schema_version: 1
updated: 2026-10-04
sources:
  - id: fetch-ok
    title: A permitted report
    publisher: Example Registry
    access: fetch
    filename: a.pdf
    url: https://example.invalid/a.pdf
  - id: fetch-bad
    title: Another permitted report
    publisher: Example Registry
    access: fetch
    filename: b.pdf
    url: https://example.invalid/b.pdf
  - id: manual-file
    title: A report a human must download
    publisher: Example Foundation
    access: manual
    filename: c.pdf
    url: null
    landing_page: https://example.invalid/landing
  - id: manual-nofile
    title: A site with no file
    publisher: Other Publisher
    access: manual
    landing_page: https://example.invalid/site
  - id: refused
    title: Refuses automation
    publisher: Refuser
    access: forbidden
    filename: d.pdf
    url: https://example.invalid/d.pdf
  - id: apply-first
    title: Needs an application
    publisher: Gatekeeper
    access: request
    filename: e.pdf
    url: https://example.invalid/e.pdf
"""


class FetchSourcesTest(unittest.TestCase):
    def setUp(self) -> None:
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.root = Path(self._tmp.name) / "repo"
        (self.root / "sources").mkdir(parents=True)
        self.catalog = self.root / "sources" / "catalog.yaml"
        self.catalog.write_text(CATALOG, encoding="utf-8")
        self.manifest = self.root / "sources" / "manifest.json"
        self.dl = self.root / "sources" / "downloads"
        for name, value in (("REPO_ROOT", self.root), ("CATALOG", self.catalog),
                            ("MANIFEST", self.manifest), ("DEFAULT_DEST", self.dl)):
            p = mock.patch.object(fs, name, value)
            p.start()
            self.addCleanup(p.stop)

    def run_main(self, *args: str, download=None) -> tuple[int, str, str]:
        out, err = io.StringIO(), io.StringIO()
        with mock.patch.object(sys, "argv", ["fetch_sources.py", *args]), \
                mock.patch.object(fs, "download", download or self._fake_download), \
                contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            code = fs.main()
        return code, out.getvalue(), err.getvalue()

    def _fake_download(self, url: str, dest: Path) -> None:
        self.requested.append(url)
        dest.write_bytes(b"fake pdf bytes for " + url.encode())

    @property
    def requested(self) -> list[str]:
        if not hasattr(self, "_requested"):
            self._requested: list[str] = []
        return self._requested

    # --- policy: only fetch and api entries are ever requested -------------------
    def test_fetchable_only_for_fetch_and_api_with_a_url(self) -> None:
        self.assertTrue(fs.fetchable({"access": "fetch", "url": "https://x/y"}))
        self.assertTrue(fs.fetchable({"access": "api", "url": "https://x/y"}))
        for access in ("manual", "request", "forbidden", None):
            self.assertFalse(fs.fetchable({"access": access, "url": "https://x/y"}))
        self.assertFalse(fs.fetchable({"access": "fetch", "url": None}))

    def test_default_run_requests_only_permitted_urls(self) -> None:
        code, _, _ = self.run_main()
        self.assertEqual(code, 0)
        self.assertEqual(sorted(self.requested), ["https://example.invalid/a.pdf", "https://example.invalid/b.pdf"])

    def test_forbidden_request_and_manual_entries_are_never_downloaded_even_when_named(self) -> None:
        for entry_id in ("refused", "apply-first", "manual-file", "manual-nofile"):
            self.requested.clear()
            self.run_main("--only", entry_id)
            self.assertEqual(self.requested, [], entry_id)

    # --- exit codes --------------------------------------------------------------
    def test_failed_download_exits_non_zero_and_names_the_entry(self) -> None:
        def boom(url: str, dest: Path) -> None:
            raise OSError("network down")
        code, out, err = self.run_main("--only", "fetch-ok", download=boom)
        self.assertEqual(code, 1)
        self.assertIn("FAILED (OSError)", out)
        self.assertIn("fetch-ok", err)

    def test_unknown_id_exits_2(self) -> None:
        code, _, err = self.run_main("--only", "nope")
        self.assertEqual(code, 2)
        self.assertIn("Unknown catalog id", err)

    # --- manual instructions -----------------------------------------------------
    def test_manual_listing_tolerates_entries_without_a_filename(self) -> None:
        code, out, _ = self.run_main("--manual")
        self.assertEqual(code, 0)
        self.assertIn("manual-nofile", out)
        self.assertIn("(no file; link only)", out)
        self.assertIn("manual-file", out)
        self.assertNotIn("example.invalid/a.pdf", out)  # fetchable entries are not listed
        self.assertIn("--record", out)

    # --- manifest never records a machine path -----------------------------------
    def test_dest_outside_the_repo_leaves_no_absolute_path_in_the_manifest(self) -> None:
        elsewhere = Path(self._tmp.name) / "elsewhere"
        self.run_main("--dest", str(elsewhere))
        text = self.manifest.read_text(encoding="utf-8")
        self.assertEqual(json.loads(text)["downloads_dir"], "(outside the repository)")
        self.assertNotIn(self._tmp.name, text)
        self.assertNotIn("\\\\", text)

    def test_default_dest_is_recorded_relative_with_forward_slashes(self) -> None:
        self.run_main()
        self.assertEqual(json.loads(self.manifest.read_text(encoding="utf-8"))["downloads_dir"], "sources/downloads")

    # --- record and verify -------------------------------------------------------
    def test_record_checksums_a_hand_saved_file_and_verify_then_checks_it(self) -> None:
        self.dl.mkdir(parents=True)
        (self.dl / "c.pdf").write_bytes(b"saved by a human in a browser")
        code, out, _ = self.run_main("--record")
        self.assertEqual(code, 0)
        self.assertIn("recorded 1 of 5", out)
        self.assertNotIn("still needed (save by hand): manual-file", out)  # it is on disk now
        rec = {r["id"]: r for r in json.loads(self.manifest.read_text(encoding="utf-8"))["files"]}
        self.assertEqual(rec["manual-file"]["sha256"], fs.sha256(self.dl / "c.pdf"))
        self.assertEqual(self.requested, [], "--record must never touch the network")
        self.assertEqual(self.run_main("--verify")[0], 0)
        (self.dl / "c.pdf").write_bytes(b"changed")
        code, out, _ = self.run_main("--verify")
        self.assertEqual(code, 1)
        self.assertIn("CHANGED  manual-file", out)
        (self.dl / "c.pdf").unlink()
        code, out, _ = self.run_main("--verify")
        self.assertEqual(code, 1)
        self.assertIn("MISSING  manual-file", out)

    def test_record_lists_manual_files_still_missing(self) -> None:
        code, out, _ = self.run_main("--record")
        self.assertEqual(code, 0)
        self.assertIn("still needed (save by hand): manual-file  ->  c.pdf", out)

    def test_verify_without_a_manifest_exits_1(self) -> None:
        code, out, _ = self.run_main("--verify")
        self.assertEqual(code, 1)
        self.assertIn("No manifest yet", out)

    def test_existing_files_are_not_downloaded_again_without_force(self) -> None:
        self.run_main("--only", "fetch-ok")
        self.requested.clear()
        _, out, _ = self.run_main("--only", "fetch-ok")
        self.assertEqual(self.requested, [])
        self.assertIn("have     fetch-ok", out)
        self.run_main("--only", "fetch-ok", "--force")
        self.assertEqual(len(self.requested), 1)

    # --- --verify honours --only (T-0054) ------------------------------------------
    def _two_recorded_files(self) -> None:
        self.dl.mkdir(parents=True)
        (self.dl / "a.pdf").write_bytes(b"first")
        (self.dl / "b.pdf").write_bytes(b"second")
        self.assertEqual(self.run_main("--record")[0], 0)

    def test_verify_only_checks_just_the_named_entry(self) -> None:
        self._two_recorded_files()
        code, out, _ = self.run_main("--only", "fetch-ok", "--verify")
        self.assertEqual(code, 0)
        self.assertIn("1 match, 0 changed, 0 missing", out)
        code, out, _ = self.run_main("--verify")
        self.assertEqual(code, 0)
        self.assertIn("2 match, 0 changed, 0 missing", out)

    def test_verify_only_ignores_changes_to_entries_not_named(self) -> None:
        self._two_recorded_files()
        (self.dl / "b.pdf").write_bytes(b"tampered")
        code, out, _ = self.run_main("--only", "fetch-ok", "--verify")
        self.assertEqual(code, 0)
        self.assertNotIn("fetch-bad", out)
        code, out, _ = self.run_main("--only", "fetch-bad", "--verify")
        self.assertEqual(code, 1)
        self.assertIn("CHANGED  fetch-bad", out)

    def test_verify_only_with_an_unknown_id_exits_2(self) -> None:
        self._two_recorded_files()
        code, out, err = self.run_main("--only", "no-such-id", "--verify")
        self.assertEqual(code, 2)
        self.assertIn("no-such-id", out + err)

    def test_verify_only_names_an_entry_with_no_recorded_hash_instead_of_a_silent_zero(self) -> None:
        self._two_recorded_files()  # manual-file (c.pdf) is not on disk, so no hash was recorded for it
        code, out, _ = self.run_main("--only", "manual-file", "--verify")
        self.assertEqual(code, 1)
        self.assertIn("NO HASH  manual-file", out)
        self.assertIn("0 match", out)
        # a full verify keeps skipping entries that have no hash
        self.assertEqual(self.run_main("--verify")[0], 0)


if __name__ == "__main__":
    unittest.main()
