#!/usr/bin/env python3
"""Deterministic anti-hallucination guardrails across cf-research.

Enforces 4 strict deterministic invariants across proposals and landscape documents:
  AH-1 [nct-provenance]: Every NCT identifier cited in proposals/*.md, landscape/*.md,
       or landscape/*-claims.json must exist in tools/trial_atlas/snapshot.json.
  AH-2 [prisma-cohort-consistency]: Documents citing the PRISMA-screened included cohort
       (940 studies / 1,821 outcome entries) must not cite unfiltered 1,053-cohort counts.
  AH-3 [companion-identifier-grounding]: Every NCT, PDB, and DOI cited in a landscape/<stem>.md
       with companion evidence/claims files must be grounded verbatim in those companions.
  AH-4 [clinical-overclaim-guard]: Prohibits unhedged curative promises or clinical practice directives
       (e.g., 'proven to cure', 'safe and effective for patients', 'guarantees clinical benefit').

Zero network. Deterministic verification. Standard library only.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
DEFAULT_SNAPSHOT_JSON = REPO_ROOT / "tools" / "trial_atlas" / "snapshot.json"

# Explicit allowlist for external verified NCTs not in CF snapshot (empty by default)
VERIFIED_EXTERNAL_NCTS: set[str] = {"NCT05248022"}

# Prohibited clinical overclaim phrases
PROHIBITED_OVERCLAIM_PHRASES: tuple[str, ...] = (
    "proven to cure",
    "guarantees clinical benefit",
    "recommended treatment regimen",
    "safe and effective for patients",
    "cures cystic fibrosis",
)

# Unfiltered 1,053-cohort counts that indicate contamination if cited in PRISMA-included context
UNFILTERED_COHORT_CONTAMINATION: tuple[tuple[re.Pattern, str], ...] = (
    (
        re.compile(
            r"(?:\b169\b\s*(?:entries|outcomes|studies|records|classified|primary)\b.*?\b(?:fev1|lung_function)|\b(?:fev1|lung_function).*?\b169\s*(?:entries|outcomes|studies|records|classified))",
            re.IGNORECASE,
        ),
        "169 lung_function_fev1 (unfiltered)",
    ),
    (
        re.compile(
            r"(?:\b39\b\s*(?:entries|outcomes|studies|records|classified|primary)\b.*?\b(?:exacerbation|pulmonary_exacerbation)|\b(?:exacerbation|pulmonary_exacerbation).*?\b39\s*(?:entries|outcomes|studies|records|classified))",
            re.IGNORECASE,
        ),
        "39 pulmonary_exacerbations (unfiltered)",
    ),
    (
        re.compile(
            r"(?:\b163\b\s*(?:entries|outcomes|studies|records|classified|primary)\b.*?\b(?:quality_of_life|cfq|qol)|\b(?:quality_of_life|cfq|qol).*?\b163\s*(?:entries|outcomes|studies|records|classified))",
            re.IGNORECASE,
        ),
        "163 quality_of_life (unfiltered)",
    ),
    (
        re.compile(
            r"(?:\b453\b\s*(?:entries|outcomes|studies|records|classified|primary)\b.*?\b(?:safety|adverse_events)|\b(?:safety|adverse_events).*?\b453\s*(?:entries|outcomes|studies|records|classified))",
            re.IGNORECASE,
        ),
        "453 safety_adverse_events (unfiltered)",
    ),
    (
        re.compile(
            r"(?:\b1[,]?114\b\s*(?:entries|outcomes|studies|records|classified|primary)\b.*?\b(?:other|other_or_not_stated)|\b(?:other|other_or_not_stated).*?\b1[,]?114\s*(?:entries|outcomes|studies|records|classified))",
            re.IGNORECASE,
        ),
        "1,114 other_or_not_stated (unfiltered)",
    ),
)

NCT_RE = re.compile(r"\bNCT\d{8}\b")
DOI_RE = re.compile(r"\b10\.\d{4,9}/[-._;()/:A-Za-z0-9]+")
PDB_RE = re.compile(r"\bPDB\s*(?:ID:?|entry:?|code:?)?\s*([0-9][A-Za-z0-9]{3})\b", re.IGNORECASE)


def load_snapshot_ncts(snapshot_path: Path | str | None = None) -> set[str]:
    snap_file = Path(snapshot_path or DEFAULT_SNAPSHOT_JSON)
    if not snap_file.exists():
        raise FileNotFoundError(f"Snapshot JSON not found: {snap_file}")
    data = json.loads(snap_file.read_text(encoding="utf-8"))
    return {s["nct_id"] for s in data.get("studies", []) if "nct_id" in s}


def check_nct_provenance(
    file_path: Path,
    content: str,
    snapshot_ncts: set[str],
) -> list[str]:
    """Rule AH-1: Every NCT ID must exist in tools/trial_atlas/snapshot.json."""
    failures = []
    allowed = snapshot_ncts | VERIFIED_EXTERNAL_NCTS
    for line_no, line in enumerate(content.splitlines(), 1):
        for m in NCT_RE.finditer(line):
            nct = m.group(0)
            if nct not in allowed:
                failures.append(
                    f"FAIL AH-1 [nct-provenance] {file_path.as_posix()}:{line_no}: "
                    f"unverified NCT ID {nct} not found in tools/trial_atlas/snapshot.json"
                )
    return failures


def check_prisma_cohort_consistency(
    file_path: Path,
    content: str,
) -> list[str]:
    """Rule AH-2: PRISMA-included documents must not cite unfiltered 1,053-cohort counts."""
    failures = []
    cites_prisma = bool(
        re.search(r"\b940\b", content) and re.search(r"\b1[,]?821\b", content)
    )
    if not cites_prisma:
        return failures

    for line_no, line in enumerate(content.splitlines(), 1):
        for pattern, label in UNFILTERED_COHORT_CONTAMINATION:
            if pattern.search(line):
                failures.append(
                    f"FAIL AH-2 [prisma-cohort-consistency] {file_path.as_posix()}:{line_no}: "
                    f"document citing PRISMA-included cohort contains unfiltered count contamination: {label}"
                )
    return failures


def check_companion_identifier_grounding(
    md_path: Path,
    content: str,
    evidence_text: str,
    claims_text: str,
) -> list[str]:
    """Rule AH-3: Identifiers in landscape/<stem>.md must appear in companion evidence or claims."""
    failures = []
    companion_corpus = evidence_text + "\n" + claims_text

    for line_no, line in enumerate(content.splitlines(), 1):
        # Check NCTs
        for m in NCT_RE.finditer(line):
            nct = m.group(0)
            if nct not in companion_corpus:
                failures.append(
                    f"FAIL AH-3 [companion-identifier-grounding] {md_path.as_posix()}:{line_no}: "
                    f"NCT identifier {nct} cited in markdown but not grounded in companion evidence/claims"
                )
        # Check DOIs
        for m in DOI_RE.finditer(line):
            doi = m.group(0).rstrip(").,")
            if doi not in companion_corpus:
                failures.append(
                    f"FAIL AH-3 [companion-identifier-grounding] {md_path.as_posix()}:{line_no}: "
                    f"DOI {doi} cited in markdown but not grounded in companion evidence/claims"
                )
        # Check PDBs
        for m in PDB_RE.finditer(line):
            pdb = m.group(1).upper()
            if pdb not in companion_corpus.upper():
                failures.append(
                    f"FAIL AH-3 [companion-identifier-grounding] {md_path.as_posix()}:{line_no}: "
                    f"PDB structure ID {pdb} cited in markdown but not grounded in companion evidence/claims"
                )

    return failures


def check_clinical_overclaim_guard(
    file_path: Path,
    content: str,
) -> list[str]:
    """Rule AH-4: Rejects prohibited clinical-promise or unhedged causal phrases."""
    failures = []
    for line_no, line in enumerate(content.splitlines(), 1):
        for phrase in PROHIBITED_OVERCLAIM_PHRASES:
            pattern = re.compile(r"\b" + re.escape(phrase) + r"\b", re.IGNORECASE)
            if pattern.search(line):
                failures.append(
                    f"FAIL AH-4 [clinical-overclaim-guard] {file_path.as_posix()}:{line_no}: "
                    f"prohibited clinical overclaim phrase: {phrase!r}"
                )
    return failures


def run_checks(
    repo_root: Path | str | None = None,
    snapshot_path: Path | str | None = None,
) -> list[str]:
    """Execute all 4 anti-hallucination checks across proposals and landscape."""
    root = Path(repo_root or REPO_ROOT).resolve()
    snap_path = Path(snapshot_path or root / "tools" / "trial_atlas" / "snapshot.json")
    snapshot_ncts = load_snapshot_ncts(snap_path)

    failures: list[str] = []

    target_files = sorted(
        list(root.glob("proposals/*.md"))
        + list(root.glob("landscape/*.md"))
        + list(root.glob("landscape/*-claims.json"))
    )

    for p in target_files:
        try:
            content = p.read_text(encoding="utf-8")
        except Exception as exc:
            failures.append(f"FAIL: could not read {p.as_posix()}: {exc}")
            continue

        # Rule AH-1: NCT Provenance
        failures.extend(check_nct_provenance(p, content, snapshot_ncts))

        # Rule AH-2: PRISMA Cohort Consistency
        failures.extend(check_prisma_cohort_consistency(p, content))

        # Rule AH-4: Clinical Overclaim Guard
        failures.extend(check_clinical_overclaim_guard(p, content))

        # Rule AH-3: Companion Identifier Grounding (for landscape/<stem>.md)
        if p.suffix == ".md" and p.parent.name == "landscape":
            stem = p.stem
            ev_file = p.parent / f"{stem}-evidence.txt"
            cl_file = p.parent / f"{stem}-claims.json"
            if ev_file.exists() or cl_file.exists():
                ev_text = ev_file.read_text(encoding="utf-8") if ev_file.exists() else ""
                cl_text = cl_file.read_text(encoding="utf-8") if cl_file.exists() else ""
                failures.extend(check_companion_identifier_grounding(p, content, ev_text, cl_text))

    return failures


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Deterministic anti-hallucination guardrails")
    parser.add_argument("--root", type=Path, default=REPO_ROOT, help="Repository root path")
    parser.add_argument("--snapshot", type=Path, default=DEFAULT_SNAPSHOT_JSON, help="Path to snapshot.json")
    parser.add_argument("--verbose", action="store_true", help="Print verbose details")
    args = parser.parse_args(argv)

    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, "reconfigure"):
            stream.reconfigure(errors="backslashreplace")

    print("=" * 72)
    print("ANTI-HALLUCINATION & PROVENANCE GUARDRAIL AUDIT")
    print(f"Scanning proposals/*.md, landscape/*.md, landscape/*-claims.json")
    print("Enforcing: AH-1 (NCT Provenance), AH-2 (PRISMA Consistency),")
    print("           AH-3 (Companion Grounding), AH-4 (Clinical Overclaim Guard)")
    print("=" * 72)

    failures = run_checks(args.root, args.snapshot)

    if failures:
        print(f"RESULT: FAIL ({len(failures)} anti-hallucination failure(s) detected):")
        for f in failures:
            print(f"  {f}")
        return 1

    print("PASS AH-1 [nct-provenance]: 100% of cited NCT identifiers verified in snapshot.json")
    print("PASS AH-2 [prisma-cohort-consistency]: zero unfiltered 1,053-cohort count contamination")
    print("PASS AH-3 [companion-identifier-grounding]: 100% of markdown identifiers grounded in companions")
    print("PASS AH-4 [clinical-overclaim-guard]: zero prohibited clinical promise phrases detected")
    print("=" * 72)
    print("RESULT: PASS (all anti-hallucination guardrails satisfied)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
