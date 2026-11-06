"""Unit tests for tools/sources/check_repo_organization.py."""
from __future__ import annotations

import shutil
import tempfile
import unittest
from pathlib import Path

try:
    from tools.sources.check_repo_organization import REPO_ROOT, audit_repo_organization
except ModuleNotFoundError:
    from check_repo_organization import REPO_ROOT, audit_repo_organization


class CheckRepoOrganizationTests(unittest.TestCase):
    """Verify repository organization invariants on live repo and synthetic negative controls."""

    def test_live_repository_passes_all_organization_invariants(self) -> None:
        report = audit_repo_organization(REPO_ROOT)
        self.assertEqual(report["status"], "PASS", f"Expected PASS, got errors: {report['errors']}")
        self.assertEqual(report["errors"], [])
        self.assertGreaterEqual(report["catalog_total_sources"], 330)
        self.assertGreaterEqual(report["landscape_files_checked"], 55)
        self.assertGreaterEqual(report["proposals_checked"], 42)
        self.assertGreaterEqual(report["data_files_checked"], 10)
        self.assertGreaterEqual(report["tool_directories_checked"], 13)
        self.assertEqual(report["skills_checked"], 5)

    def test_negative_control_detects_readme_catalog_count_drift_and_unindexed_files(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "sources").mkdir()
            (root / "landscape").mkdir()
            (root / "proposals").mkdir()
            (root / "data").mkdir()
            (root / "tools" / "demo_tool").mkdir(parents=True)
            (root / ".claude" / "skills" / "demo-skill").mkdir(parents=True)

            (root / "sources" / "catalog.yaml").write_text(
                "updated: '2026-10-31'\nsources:\n  - id: s1\n    access: fetch\n",
                encoding="utf-8",
            )
            # Deliberately wrong count in README.md ("catalog of 999 sources")
            (root / "README.md").write_text("# Demo\ncatalog of 999 sources\n", encoding="utf-8")
            (root / "landscape" / "README.md").write_text("# Landscape\n", encoding="utf-8")
            (root / "landscape" / "unindexed.md").write_text("# Unindexed\n", encoding="utf-8")
            (root / "proposals" / "README.md").write_text("# Proposals\n", encoding="utf-8")
            (root / "proposals" / "2026-10-99-missing-status.md").write_text("# No status\n", encoding="utf-8")
            (root / "data" / "README.md").write_text("# Data\n", encoding="utf-8")

            report = audit_repo_organization(root)
            self.assertEqual(report["status"], "FAIL")
            joined = "\n".join(report["errors"])
            self.assertIn("Root README.md missing expected catalog count phrase", joined)
            self.assertIn("Unindexed files in landscape/README.md", joined)
            self.assertIn("Unindexed files in proposals/README.md", joined)
            self.assertIn("Proposals missing explicit Status header", joined)


if __name__ == "__main__":
    unittest.main()
