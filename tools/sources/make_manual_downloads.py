#!/usr/bin/env python3
"""Generate a Markdown page listing manual-download sources from the lab catalog."""
from __future__ import annotations

import re
import sys
from pathlib import Path

import yaml


PROMPT = (
    "I have saved the manual downloads into sources/downloads/ under the file names in manual-downloads.md.\n"
    "1. Run `python tools/sources/fetch_sources.py --record`, then `python tools/sources/fetch_sources.py --verify`, and tell me which files are missing or changed. Do not fetch anything over the network.\n"
    "2. For each newly recorded file, open its title page and read it. Fill data_year, pages and text_layer in sources/catalog.yaml from what the document says, quoting the title-page text. Set claim_label to verified only for what you read.\n"
    "3. Run `python tools/sources/check_catalog_fields.py`, then sync the skill's catalog copy with `python tools/sources/sync_skill_catalog.py --source sources/catalog.yaml --target ../cf-skills/.claude/skills/cf-evidence-loop/catalog.yaml`.\n"
    "4. Show me the diff and wait for my approval before you commit anything."
)

_WHITESPACE_RUN = re.compile(r"\s+")


def _collapse(text: str) -> str:
    return _WHITESPACE_RUN.sub(" ", text).strip()


def _cell(text: str) -> str:
    """Collapse whitespace and escape pipes for a Markdown table cell."""
    collapsed = _collapse(text)
    escaped = collapsed.replace("|", "\\|")
    # remove any trailing spaces (shouldn't happen after strip, but be safe)
    return escaped.rstrip()


def _first_sentence(notes: str) -> str:
    """Collapse whitespace, cut after first '. ' occurrence, truncate to 160 chars."""
    collapsed = _WHITESPACE_RUN.sub(" ", notes).strip()
    # Find the first occurrence of ". " (full stop followed by a space)
    idx = collapsed.find(". ")
    if idx != -1:
        sentence = collapsed[: idx + 1]  # keep the full stop
    else:
        sentence = collapsed
    # Truncate to 160 characters
    if len(sentence) > 160:
        sentence = sentence[:159].rstrip() + "\u2026"
    return sentence


def _open_cell(source: dict) -> str:
    """Build the Open cell content."""
    parts = []
    landing_page = source.get("landing_page")
    url = source.get("url")
    if landing_page:
        parts.append(f"[page]({landing_page})")
    if url:
        parts.append(f"[file]({url})")
    return " ".join(parts)


def _save_as_cell(source: dict) -> str:
    """Build the Save as cell content."""
    filename = source.get("filename")
    if filename:
        return f"`{filename}`"
    return "link only"


def _notes_cell(source: dict) -> str:
    """Build the Notes cell content."""
    parts = []
    warning = source.get("warning")
    if warning:
        parts.append(f"warning: {warning}")
    notes = source.get("notes")
    if notes:
        parts.append(_first_sentence(notes))
    return " ".join(parts)


def render(catalog: dict) -> str:
    """Render the full Markdown page as a string."""
    sources = catalog.get("sources", [])

    manual_sources = [s for s in sources if s.get("access") == "manual"]
    request_sources = [s for s in sources if s.get("access") == "request"]
    forbidden_sources = [s for s in sources if s.get("access") == "forbidden"]

    N = len(manual_sources)
    M = len(request_sources)
    K = len(forbidden_sources)

    lines: list[str] = []

    # 1. Title
    lines.append("# Manual downloads")
    lines.append("")

    # 2. Generated line
    lines.append("Generated from sources/catalog.yaml by tools/sources/make_manual_downloads.py. Do not edit by hand.")
    lines.append("")

    # 3. Counts
    lines.append(f"{N} to download by hand, {M} to ask for, {K} not to download.")
    lines.append("")

    # 4. How this works
    lines.append("## How this works")
    lines.append("")
    lines.append("1. Open the link for a source in your browser.")
    lines.append("2. Save the file into `sources/downloads/` under the exact file name shown.")
    lines.append("3. When you have saved the files, paste the prompt at the bottom into your agent.")
    lines.append("")

    # 5. Manual downloads section (only if N > 0)
    if N > 0:
        lines.append("## Download in a browser")
        lines.append("")

        # Group manual sources by publisher, preserving catalog order within each group
        publisher_map: dict[str, list[dict]] = {}
        for s in manual_sources:
            pub = s.get("publisher")
            if not pub:
                key = "Unknown publisher"
            else:
                key = pub
            # Use the raw publisher name as the display key, but sort by lowercase
            if key not in publisher_map:
                publisher_map[key] = []
            publisher_map[key].append(s)

        # Sort publishers alphabetically ignoring case
        sorted_publishers = sorted(publisher_map.keys(), key=str.lower)

        for pub in sorted_publishers:
            lines.append(f"### {pub}")
            lines.append("")
            lines.append("| ID | What | Open | Save as | Notes |")
            lines.append("| --- | --- | --- | --- | --- |")
            for s in publisher_map[pub]:
                sid = _cell(str(s.get("id", "")))
                title = _cell(str(s.get("title", "")))
                open_cell = _cell(_open_cell(s))
                save_as = _cell(_save_as_cell(s))
                notes_cell = _cell(_notes_cell(s))
                lines.append(f"| {sid} | {title} | {open_cell} | {save_as} | {notes_cell} |")
            lines.append("")

    # 6. Request section (only if M > 0)
    if M > 0:
        lines.append("## Ask first (access by request)")
        lines.append("")
        for s in request_sources:
            title = _collapse(str(s.get("title", "")))
            pub = s.get("publisher")
            pub_display = pub if pub else "unknown"
            landing_page = s.get("landing_page", "")
            lines.append(f"- **{title}** ({pub_display}): [page]({landing_page}) - needs permission or an application; do not automate.")
        lines.append("")

    # 7. Forbidden section (only if K > 0)
    if K > 0:
        lines.append("## Do not download")
        lines.append("")
        for s in forbidden_sources:
            title = _collapse(str(s.get("title", "")))
            pub = s.get("publisher")
            pub_display = pub if pub else "unknown"
            landing_page = s.get("landing_page", "")
            lines.append(f"- **{title}** ({pub_display}): [page]({landing_page}) - the publisher refuses automated access.")
        lines.append("")

    # 8. Prompt section
    lines.append("## When you have saved the files")
    lines.append("")
    lines.append("```text")
    lines.extend(PROMPT.split("\n"))
    lines.append("```")

    result = "\n".join(lines) + "\n"
    # Ensure no line ends with a space (already handled by strip in _cell and collapse)
    return result


def main(argv=None) -> int:
    if argv is None:
        argv = sys.argv[1:]

    import argparse

    parser = argparse.ArgumentParser(description="Generate manual-downloads.md from the catalog.")
    # Default paths are relative to the script's location two folders up
    script_dir = Path(__file__).resolve().parent  # <repo>/tools/sources/
    repo_root = script_dir.parent.parent  # <repo>
    default_catalog = repo_root / "sources" / "catalog.yaml"
    default_out = repo_root / "sources" / "manual-downloads.md"

    parser.add_argument("--catalog", type=str, default=str(default_catalog))
    parser.add_argument("--out", type=str, default=str(default_out))
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args(argv)

    catalog_path = Path(args.catalog)
    out_path = Path(args.out)

    # Read the catalog
    with open(catalog_path, "r", encoding="utf-8") as f:
        catalog = yaml.safe_load(f)

    rendered = render(catalog)
    rendered_bytes = rendered.encode("utf-8")

    if args.check:
        if out_path.exists():
            existing = out_path.read_bytes()
            if existing == rendered_bytes:
                print("manual-downloads.md is up to date")
                return 0
        print("manual-downloads.md is out of date: run python tools/sources/make_manual_downloads.py")
        return 1

    # Write the output
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_bytes(rendered_bytes)
    print(f"wrote {out_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
