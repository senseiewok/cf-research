#!/usr/bin/env python3
"""Deterministic anti-hallucination guardrails across cf-research.

Enforces 5 strict deterministic invariants across proposals, landscape documents,
and tools/discovery biophysical models:
  AH-1 [nct-provenance]: Every NCT identifier cited in proposals/*.md, landscape/*.md,
       or landscape/*-claims.json must exist in tools/trial_atlas/snapshot.json.
  AH-2 [prisma-cohort-consistency]: Documents citing the PRISMA-screened included cohort
       (940 studies / 1,821 outcome entries) must not cite unfiltered 1,053-cohort counts.
  AH-3 [companion-identifier-grounding]: Every NCT, PDB, and DOI cited in a landscape/<stem>.md
       with companion evidence/claims files must be grounded verbatim in those companions.
  AH-4 [clinical-overclaim-guard]: Prohibits unhedged curative promises or clinical practice directives
       (e.g., 'proven to cure', 'safe and effective for patients', 'guarantees clinical benefit').
  AH-5 [discovery-module-registration]: Every biophysical model in tools/discovery/*.py
       (excluding test_*.py and infrastructure scripts) must be registered in
       discovery_council.DISCOVERY_MODULE_REGISTRY with verified NCT provenance and
       zero prohibited clinical overclaims.
  AH-6 [cite-check-ledger-provenance]: Every DOI, PMID, and NCT identifier cited in
       landscape/*.md and tools/discovery/*.py must be verified with verdict=PASS in
       sources/cite-check-ledger.jsonl (offline ledger generated via cf-evidence-loop cite-check).

Zero network. Deterministic verification. Standard library only.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
DEFAULT_SNAPSHOT_JSON = REPO_ROOT / "tools" / "trial_atlas" / "snapshot.json"

# Explicit allowlist for external verified NCTs not in CF snapshot (empty by default)
VERIFIED_EXTERNAL_NCTS: set[str] = set()

# Extrapolated Boxed Warning guard (openFDA SPL Boxed Warning is `not stated` for Kalydeco, Orkambi, Symdeko)
EXTRAPOLATED_BOXED_WARNING_RE = re.compile(
    r"\b(?:kalydeco|orkambi|symdeko)\b[^.\n]*?\b(?:boxed\s+warning|black\s+box\s+warning)\b",
    re.IGNORECASE,
)
BOXED_WARNING_NEGATION_RE = re.compile(
    r"\b(?:not\s+stated|none|no\b|lacks?|without|absent|versus|vs\.?)",
    re.IGNORECASE,
)

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
PMID_RE = re.compile(r"\bPMID[:\s]+(\d{6,9})\b", re.IGNORECASE)
PDB_RE = re.compile(r"\bPDB\s*(?:ID:?|entry:?|code:?)?\s*([0-9][A-Za-z0-9]{3})\b", re.IGNORECASE)
DEFAULT_CITE_LEDGER = REPO_ROOT / "sources" / "cite-check-ledger.jsonl"


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
        for m in EXTRAPOLATED_BOXED_WARNING_RE.finditer(line):
            if not BOXED_WARNING_NEGATION_RE.search(line):
                failures.append(
                    f"FAIL AH-4 [clinical-overclaim-guard] {file_path.as_posix()}:{line_no}: "
                    f"extrapolated Boxed Warning on Kalydeco/Orkambi/Symdeko without 'not stated' boundary: {m.group(0)!r}"
                )
    return failures


def check_discovery_module_registration(
    repo_root: Path,
    snapshot_ncts: set[str],
) -> list[str]:
    """Rule AH-5: Every biophysical model in tools/discovery/*.py must be registered in discovery_council.py."""
    failures: list[str] = []
    discovery_dir = repo_root / "tools" / "discovery"
    council_file = discovery_dir / "discovery_council.py"
    if not discovery_dir.exists() or not council_file.exists():
        return failures

    spec = importlib.util.spec_from_file_location("discovery_council", council_file)
    if spec is None or spec.loader is None:
        return [f"FAIL AH-5 [discovery-module-registration]: cannot load {council_file.as_posix()}"]
    mod = importlib.util.module_from_spec(spec)
    sys.path.insert(0, str(discovery_dir))
    try:
        spec.loader.exec_module(mod)
    finally:
        if sys.path and sys.path[0] == str(discovery_dir):
            sys.path.pop(0)

    registry = getattr(mod, "DISCOVERY_MODULE_REGISTRY", {})
    infra_modules = getattr(mod, "INFRASTRUCTURE_MODULES", {"__init__", "discovery_council", "cross_domain_miner", "export_delivery_web_data"})
    audit_fn = getattr(mod, "audit_discovery_module_registry", None)
    if audit_fn is None:
        return [
            f"FAIL AH-5 [discovery-module-registration] {council_file.as_posix()}: "
            "missing audit_discovery_module_registry() or DISCOVERY_MODULE_REGISTRY"
        ]

    audit_res = audit_fn(discovery_dir)
    for mod_name in audit_res.get("unregistered_modules", []):
        failures.append(
            f"FAIL AH-5 [discovery-module-registration] {council_file.as_posix()}: "
            f"unregistered biophysical module '{mod_name}'"
        )
    for mod_name in audit_res.get("missing_test_suites", []):
        failures.append(
            f"FAIL AH-5 [discovery-module-registration] {council_file.as_posix()}: "
            f"biophysical module '{mod_name}' missing companion unit test suite"
        )

    for py_file in sorted(discovery_dir.glob("*.py")):
        if py_file.name.startswith("test_") or py_file.name.startswith("_"):
            continue
        content = py_file.read_text(encoding="utf-8")
        for nct_fail in check_nct_provenance(py_file, content, snapshot_ncts):
            failures.append(nct_fail.replace("AH-1 [nct-provenance]", "AH-5 [discovery-module-registration]"))
        for oc_fail in check_clinical_overclaim_guard(py_file, content):
            failures.append(oc_fail.replace("AH-4 [clinical-overclaim-guard]", "AH-5 [discovery-module-registration]"))

        if py_file.stem not in infra_modules and py_file.stem not in registry:
            failures.append(
                f"FAIL AH-5 [discovery-module-registration] {py_file.as_posix()}: "
                f"biophysical module '{py_file.stem}' is not registered in discovery_council.DISCOVERY_MODULE_REGISTRY"
            )

    return failures


def load_cite_check_ledger(ledger_path: Path | str | None = None) -> set[tuple[str, str]]:
    """Load verified (kind, id_lower) pairs with verdict == 'PASS' from sources/cite-check-ledger.jsonl."""
    lp = Path(ledger_path or DEFAULT_CITE_LEDGER)
    if not lp.exists():
        return set()
    verified: set[tuple[str, str]] = set()
    for line in lp.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line:
            continue
        rec = json.loads(line)
        if rec.get("verdict") == "PASS" and rec.get("kind") and rec.get("id"):
            verified.add((str(rec["kind"]).upper(), str(rec["id"]).lower()))
    return verified


def check_cite_check_ledger_provenance(
    file_path: Path,
    content: str,
    verified_ledger: set[tuple[str, str]],
) -> list[str]:
    """Rule AH-6: Every DOI, PMID, and NCT in landscape/*.md and tools/discovery/*.py must be in cite-check-ledger.jsonl."""
    failures: list[str] = []
    for line_no, line in enumerate(content.splitlines(), 1):
        for m in DOI_RE.finditer(line):
            doi = m.group(0).rstrip(").,:;")
            if ("DOI", doi.lower()) not in verified_ledger:
                failures.append(
                    f"FAIL AH-6 [cite-check-ledger-provenance] {file_path.as_posix()}:{line_no}: "
                    f"DOI {doi} not verified in sources/cite-check-ledger.jsonl"
                )
        for m in PMID_RE.finditer(line):
            pmid = m.group(1)
            if ("PMID", pmid.lower()) not in verified_ledger:
                failures.append(
                    f"FAIL AH-6 [cite-check-ledger-provenance] {file_path.as_posix()}:{line_no}: "
                    f"PMID {pmid} not verified in sources/cite-check-ledger.jsonl"
                )
        for m in NCT_RE.finditer(line):
            nct = m.group(0).upper()
            if ("NCT", nct.lower()) not in verified_ledger:
                failures.append(
                    f"FAIL AH-6 [cite-check-ledger-provenance] {file_path.as_posix()}:{line_no}: "
                    f"NCT {nct} not verified in sources/cite-check-ledger.jsonl"
                )
    return failures


def run_checks(
    repo_root: Path | str | None = None,
    snapshot_path: Path | str | None = None,
) -> list[str]:
    """Execute all 6 anti-hallucination checks across proposals, landscape, and tools/discovery."""
    root = Path(repo_root or REPO_ROOT).resolve()
    snap_path = Path(snapshot_path or root / "tools" / "trial_atlas" / "snapshot.json")
    snapshot_ncts = load_snapshot_ncts(snap_path)
    ledger_path = root / "sources" / "cite-check-ledger.jsonl"
    verified_ledger = load_cite_check_ledger(ledger_path)

    failures: list[str] = []
    if not ledger_path.exists():
        failures.append(
            f"FAIL AH-6 [cite-check-ledger-provenance]: missing {ledger_path.as_posix()}"
        )

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
            base_stem = re.sub(r"-\d{4}$", "", stem)
            ev_file = p.parent / f"{stem}-evidence.txt"
            if not ev_file.exists() and base_stem != stem:
                ev_file = p.parent / f"{base_stem}-evidence.txt"
            cl_file = p.parent / f"{stem}-claims.json"
            if not cl_file.exists() and base_stem != stem:
                cl_file = p.parent / f"{base_stem}-claims.json"
            if ev_file.exists() or cl_file.exists():
                ev_text = ev_file.read_text(encoding="utf-8") if ev_file.exists() else ""
                cl_text = cl_file.read_text(encoding="utf-8") if cl_file.exists() else ""
                failures.extend(check_companion_identifier_grounding(p, content, ev_text, cl_text))
            if verified_ledger:
                failures.extend(check_cite_check_ledger_provenance(p, content, verified_ledger))

    # Rule AH-5: Discovery Module Registration & Provenance
    failures.extend(check_discovery_module_registration(root, snapshot_ncts))

    # Rule AH-6 on tools/discovery/*.py
    if verified_ledger:
        for py_file in sorted((root / "tools" / "discovery").glob("*.py")):
            if py_file.name.startswith("test_"):
                continue
            py_content = py_file.read_text(encoding="utf-8")
            failures.extend(check_cite_check_ledger_provenance(py_file, py_content, verified_ledger))

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
    print("Scanning proposals/*.md, landscape/*.md, landscape/*-claims.json, tools/discovery/*.py")
    print("Enforcing: AH-1 (NCT Provenance), AH-2 (PRISMA Consistency),")
    print("           AH-3 (Companion Grounding), AH-4 (Clinical Overclaim Guard),")
    print("           AH-5 (Discovery Module Registration & Council Verification),")
    print("           AH-6 (Cite-Check Ledger Provenance for DOIs, PMIDs, and NCTs)")
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
    print("PASS AH-5 [discovery-module-registration]: 100% of tools/discovery models registered & verified")
    print("PASS AH-6 [cite-check-ledger-provenance]: 100% of DOIs, PMIDs, and NCTs verified in sources/cite-check-ledger.jsonl")
    print("=" * 72)
    print("RESULT: PASS (all anti-hallucination guardrails satisfied)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
