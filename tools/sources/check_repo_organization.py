#!/usr/bin/env python3
"""Deterministic repository organization, README index, and source catalog sync checker for cf-research.

Enforces six structural invariants across senseiewok/cf-research:
1. Root README.md source catalog counts (total, fetch, api, manual, forbidden, request, updated date)
   match sources/catalog.yaml exactly.
2. Every .md file and structured verification sidecar in landscape/ is indexed in landscape/README.md.
3. proposals/README.md exists, indexes every .md file in proposals/, and every proposal declares an
   explicit Status field.
4. Every file under data/ is documented in data/README.md.
5. Every tool package directory under tools/<subdir>/ has a README.md and is listed in root README.md.
6. Every repo-scoped skill under .claude/skills/<skill>/SKILL.md has valid YAML frontmatter and is
   listed in root README.md.
"""
from __future__ import annotations

import collections
import json
import re
import sys
from pathlib import Path
from typing import Any, Dict, List

import yaml

REPO_ROOT = Path(__file__).resolve().parent.parent.parent


def audit_repo_organization(repo_root: Path | None = None) -> Dict[str, Any]:
    root = Path(repo_root or REPO_ROOT)
    errors: List[str] = []

    # 1. Catalog counts vs root README.md
    catalog_path = root / "sources" / "catalog.yaml"
    readme_path = root / "README.md"
    cat_data = yaml.safe_load(catalog_path.read_text(encoding="utf-8"))
    readme_text = readme_path.read_text(encoding="utf-8")

    sources = cat_data.get("sources", [])
    updated = str(cat_data.get("updated", ""))
    mode_counts = collections.Counter(
        s.get("access", {}).get("mode") if isinstance(s.get("access"), dict) else s.get("access")
        for s in sources
    )
    total_sources = len(sources)
    expected_readme_Phrases = [
        f"catalog of {total_sources} sources",
        f"held {total_sources} sources (`updated: {updated}`)",
        f"{mode_counts.get('fetch', 0)} a tool may fetch",
        f"{mode_counts.get('api', 0)} APIs with a recorded usage policy",
        f"{mode_counts.get('manual', 0)} a person must download by hand",
        f"{mode_counts.get('forbidden', 0)} refused by their publisher",
        f"{mode_counts.get('request', 0)} that needs a written request",
    ]
    for phrase in expected_readme_Phrases:
        if phrase not in readme_text:
            errors.append(f"Root README.md missing expected catalog count phrase: {phrase!r}")

    # 2. landscape/README.md completeness
    land_dir = root / "landscape"
    land_readme = (land_dir / "README.md").read_text(encoding="utf-8")
    land_files = sorted(
        f.relative_to(land_dir).as_posix()
        for f in land_dir.rglob("*")
        if f.is_file() and f.name != "README.md"
    )
    unindexed_landscape = [
        rel for rel in land_files
        if rel not in land_readme and Path(rel).name not in land_readme
    ]
    if unindexed_landscape:
        errors.append(f"Unindexed files in landscape/README.md: {unindexed_landscape}")

    # 3. proposals/README.md completeness & status headers
    prop_dir = root / "proposals"
    prop_readme_path = prop_dir / "README.md"
    if not prop_readme_path.is_file():
        errors.append("Missing proposals/README.md index")
        prop_readme = ""
    else:
        prop_readme = prop_readme_path.read_text(encoding="utf-8")

    prop_files = sorted(f for f in prop_dir.glob("*.md") if f.name != "README.md")
    unindexed_proposals: List[str] = []
    missing_proposal_status: List[str] = []
    status_re = re.compile(
        r"(?i)(?:\*\*status\*\*\s*:|\*\*status:\*\*|^status\s*:|- \*\*status\*\*\s*:)\s*([^\n]+)",
        re.M,
    )
    for pf in prop_files:
        if pf.name not in prop_readme:
            unindexed_proposals.append(pf.name)
        txt = pf.read_text(encoding="utf-8")
        if not status_re.search(txt):
            missing_proposal_status.append(pf.name)
    if unindexed_proposals:
        errors.append(f"Unindexed files in proposals/README.md: {unindexed_proposals}")
    if missing_proposal_status:
        errors.append(f"Proposals missing explicit Status header: {missing_proposal_status}")

    # 4. data/README.md completeness
    data_dir = root / "data"
    data_readme = (data_dir / "README.md").read_text(encoding="utf-8")
    data_files = sorted(
        f.relative_to(data_dir).as_posix()
        for f in data_dir.rglob("*")
        if f.is_file() and f.name != "README.md"
    )
    unindexed_data = [
        rel for rel in data_files
        if rel not in data_readme and Path(rel).name not in data_readme
    ]
    if unindexed_data:
        errors.append(f"Unindexed files in data/README.md: {unindexed_data}")

    # 5. tools/<subdir>/ README.md and root README.md coverage
    tools_dir = root / "tools"
    tool_subdirs = sorted(
        d.name for d in tools_dir.iterdir()
        if d.is_dir() and not d.name.startswith(".") and d.name != "__pycache__"
    )
    missing_tool_readmes: List[str] = []
    unindexed_tool_dirs: List[str] = []
    for td in tool_subdirs:
        if not (tools_dir / td / "README.md").is_file():
            missing_tool_readmes.append(td)
        if f"tools/{td}/" not in readme_text:
            unindexed_tool_dirs.append(td)
    if missing_tool_readmes:
        errors.append(f"Tool directories missing README.md: {missing_tool_readmes}")
    if unindexed_tool_dirs:
        errors.append(f"Tool directories missing from root README.md Layout table: {unindexed_tool_dirs}")

    # 6. .claude/skills/<skill>/SKILL.md frontmatter and root README.md coverage
    skills_dir = root / ".claude" / "skills"
    skill_dirs = sorted(
        d.name for d in skills_dir.iterdir()
        if d.is_dir() and not d.name.startswith(".")
    )
    invalid_skills: List[str] = []
    for sd in skill_dirs:
        skill_md = skills_dir / sd / "SKILL.md"
        if not skill_md.is_file():
            invalid_skills.append(f"{sd}:missing_SKILL.md")
            continue
        stxt = skill_md.read_text(encoding="utf-8")
        if not stxt.startswith("---"):
            invalid_skills.append(f"{sd}:missing_frontmatter")
            continue
        parts = stxt.split("---", 2)
        if len(parts) < 3:
            invalid_skills.append(f"{sd}:unclosed_frontmatter")
            continue
        fm = yaml.safe_load(parts[1]) or {}
        if fm.get("name") != sd or not fm.get("description"):
            invalid_skills.append(f"{sd}:invalid_name_or_description")
        if f".claude/skills/{sd}/" not in readme_text:
            invalid_skills.append(f"{sd}:missing_in_root_README")
    if invalid_skills:
        errors.append(f"Invalid or unindexed skills: {invalid_skills}")

    return {
        "status": "PASS" if not errors else "FAIL",
        "catalog_total_sources": total_sources,
        "catalog_updated": updated,
        "catalog_access_breakdown": dict(mode_counts),
        "landscape_files_checked": len(land_files),
        "proposals_checked": len(prop_files),
        "data_files_checked": len(data_files),
        "tool_directories_checked": len(tool_subdirs),
        "skills_checked": len(skill_dirs),
        "errors": errors,
    }


def main() -> None:
    report = audit_repo_organization()
    if "--json" in sys.argv:
        print(json.dumps(report, indent=2))
    else:
        print(
            f"Repository Organization Check: {report['status']} "
            f"(sources={report['catalog_total_sources']}, landscape={report['landscape_files_checked']}, "
            f"proposals={report['proposals_checked']}, data={report['data_files_checked']}, "
            f"tools={report['tool_directories_checked']}, skills={report['skills_checked']})"
        )
        for err in report["errors"]:
            print(f"  ERROR: {err}", file=sys.stderr)
    sys.exit(0 if report["status"] == "PASS" else 1)


if __name__ == "__main__":
    main()
