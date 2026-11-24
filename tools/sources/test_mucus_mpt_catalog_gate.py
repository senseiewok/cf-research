#!/usr/bin/env python3
"""Deterministic catalog-gate checks for Council Item #3 (mucus MPT) and Item #1 (CFTR DMS)."""
from __future__ import annotations

import re
import unittest
from pathlib import Path
from typing import Dict, List

REPO_ROOT = Path(__file__).resolve().parents[2]
CATALOG_PATH = REPO_ROOT / "sources" / "catalog.yaml"
FEASIBILITY_RECORD = REPO_ROOT / "proposals" / "2026-11-22-mucus-mpt-data-feasibility.md"
MPT_KEYWORD_RE = re.compile(r"\b(?:mucus|sputum|particle[- ]tracking|mpt)\b")


def _strip_scalar(raw_value: str) -> str:
    stripped = raw_value.strip()
    if stripped.startswith(('"', "'")):
        quote = stripped[0]
        end = stripped.find(quote, 1)
        if end != -1:
            return stripped[1:end].strip()
    val = re.split(r"(?:^|\s)#", stripped, maxsplit=1)[0].strip()
    return val


def _parse_catalog_entries_stdlib(catalog_path: Path) -> List[Dict[str, str]]:
    """Pure-stdlib parser extracting top-level fields from `sources/catalog.yaml`."""
    entries: List[Dict[str, str]] = []
    current: Dict[str, str] | None = None
    in_sources = False
    in_access_block = False
    active_block_key: str | None = None
    active_block_sep = " "

    for raw_line in catalog_path.read_text(encoding="utf-8").splitlines():
        if not raw_line.strip() or raw_line.lstrip().startswith("#"):
            continue
        indent = len(raw_line) - len(raw_line.lstrip(" "))
        stripped = raw_line.strip()

        if indent == 0:
            in_sources = stripped.startswith("sources:")
            in_access_block = False
            active_block_key = None
            continue
        if not in_sources:
            continue

        if indent == 2 and stripped.startswith("- id:"):
            if current is not None:
                entries.append(current)
            current = {"id": _strip_scalar(stripped.split(":", 1)[1])}
            in_access_block = False
            active_block_key = None
        elif current is not None:
            if indent == 4 and ":" in stripped:
                in_access_block = False
                active_block_key = None
                key, raw_rest = stripped.split(":", 1)
                key = key.strip()
                val = _strip_scalar(raw_rest)
                if key == "access" and not val:
                    in_access_block = True
                elif val in (">-", ">", "|", "|-"):
                    active_block_key = key
                    active_block_sep = "\n" if val.startswith("|") else " "
                    current[key] = ""
                else:
                    current[key] = val
            elif indent >= 6:
                if in_access_block and stripped.startswith("mode:"):
                    current["access"] = _strip_scalar(stripped.split(":", 1)[1])
                elif active_block_key is not None:
                    prev = current.get(active_block_key, "")
                    current[active_block_key] = f"{prev}{active_block_sep}{stripped}".strip() if prev else stripped

    if current is not None:
        entries.append(current)
    return entries


class TestMucusMptAndDmsCatalogGate(unittest.TestCase):
    def test_no_admitted_raw_airway_mpt_dataset_in_catalog(self) -> None:
        entries = _parse_catalog_entries_stdlib(CATALOG_PATH)
        self.assertGreaterEqual(len(entries), 300)

        admitted_datasets = [
            e for e in entries if e.get("kind") == "dataset" and e.get("access") in ("api", "fetch")
        ]
        for entry in admitted_datasets:
            blob = " ".join(entry.get(k, "") for k in ("id", "title", "notes")).lower()
            self.assertIsNone(
                MPT_KEYWORD_RE.search(blob),
                msg=f"Unexpected admitted MPT/mucus dataset {entry.get('id')!r}; update feasibility record.",
            )

    def test_dms_catalog_entries_stay_manual_and_match_verification_state(self) -> None:
        by_id = {e["id"]: e for e in _parse_catalog_entries_stdlib(CATALOG_PATH)}
        expected_labels = {
            "mavedb-dump-zenodo": "verified",
            "cftr-dms-tedman-2025": "verified",
            "cftr-dms-mckee-2023": "unverified",
        }
        for source_id, expected_label in expected_labels.items():
            self.assertIn(source_id, by_id)
            self.assertEqual(by_id[source_id].get("access"), "manual")
            self.assertEqual(by_id[source_id].get("claim_label"), expected_label)

    def test_feasibility_record_states_blocked_and_excludes_pig_gastric_gel(self) -> None:
        self.assertTrue(FEASIBILITY_RECORD.is_file())
        text = FEASIBILITY_RECORD.read_text(encoding="utf-8")
        self.assertIn("BLOCKED", text)
        self.assertIn("could not test: no data", text)
        self.assertIn("zenodo5533877_mucus_gel.npz", text)


if __name__ == "__main__":
    unittest.main()
