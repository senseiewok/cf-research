"""Deterministic grounding test suite for the 4 Track Dossiers in landscape/."""
from __future__ import annotations

import subprocess
import sys
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
LANDSCAPE_DIR = REPO_ROOT / "landscape"

FOUR_TRACK_SLUGS = (
    "rare-genotype-trial-screening-2026",
    "airway-mucus-barrier-biophysics-2026",
    "retraction-citation-integrity-2026",
    "cftr-burst-kinetics-structural-grounding-2026",
)


class TestFourTrackGroundingArticles(unittest.TestCase):
    def test_all_four_track_bundles_exist(self) -> None:
        for slug in FOUR_TRACK_SLUGS:
            for ext in (".md", "-claims.json", "-evidence.txt", "-evidence-map.json", "-allow.txt"):
                path = LANDSCAPE_DIR / f"{slug}{ext}"
                self.assertTrue(path.is_file(), f"Missing grounding artifact: {path}")

    def test_check_claims_passes_on_all_four_tracks(self) -> None:
        script = REPO_ROOT / "tools" / "claims" / "check_claims.py"
        for slug in FOUR_TRACK_SLUGS:
            claims_path = LANDSCAPE_DIR / f"{slug}-claims.json"
            map_path = LANDSCAPE_DIR / f"{slug}-evidence-map.json"
            proc = subprocess.run(
                [sys.executable, str(script), str(claims_path), "--evidence-map", str(map_path)],
                cwd=str(REPO_ROOT),
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(proc.returncode, 0, f"check_claims failed for {slug}:\n{proc.stdout}\n{proc.stderr}")
            self.assertIn("0 failures", proc.stdout)

    def test_check_claims_md_companion_mode_passes_on_all_four_tracks(self) -> None:
        script = REPO_ROOT / "tools" / "claims" / "check_claims.py"
        for slug in FOUR_TRACK_SLUGS:
            md_path = LANDSCAPE_DIR / f"{slug}.md"
            proc = subprocess.run(
                [sys.executable, str(script), str(md_path)],
                cwd=str(REPO_ROOT),
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(proc.returncode, 0, f"check_claims companion mode failed for {slug}.md:\n{proc.stdout}\n{proc.stderr}")
            self.assertIn("0 failures", proc.stdout)

    def test_check_numbers_passes_on_all_four_tracks(self) -> None:
        script = REPO_ROOT / "tools" / "sources" / "check_numbers.py"
        for slug in FOUR_TRACK_SLUGS:
            md_path = LANDSCAPE_DIR / f"{slug}.md"
            ev_path = LANDSCAPE_DIR / f"{slug}-evidence.txt"
            allow_path = LANDSCAPE_DIR / f"{slug}-allow.txt"
            proc = subprocess.run(
                [
                    sys.executable,
                    str(script),
                    "--evidence",
                    str(ev_path),
                    "--allow-file",
                    str(allow_path),
                    str(md_path),
                ],
                cwd=str(REPO_ROOT),
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(proc.returncode, 0, f"check_numbers failed for {slug}:\n{proc.stdout}\n{proc.stderr}")

    def test_check_source_overlap_quotes_passes_on_all_four_tracks(self) -> None:
        script = REPO_ROOT / "tools" / "sources" / "check_source_overlap.py"
        for slug in FOUR_TRACK_SLUGS:
            md_path = LANDSCAPE_DIR / f"{slug}.md"
            proc = subprocess.run(
                [sys.executable, str(script), "--quotes-only", str(md_path)],
                cwd=str(REPO_ROOT),
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(proc.returncode, 0, f"check_source_overlap failed for {slug}:\n{proc.stdout}\n{proc.stderr}")


if __name__ == "__main__":
    unittest.main()
