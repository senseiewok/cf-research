"""Check that every number in a draft appears in the evidence it was written from, or is allowed with a written reason.

A model that writes a note can invent a number, change a digit or do arithmetic wrongly, and the result reads well. This guard catches the
first two and forces the third into the open: a number that is not in the evidence must be listed in an allow file with the reason, for
example `564  # 127 + 437`, so a person can check the arithmetic. It checks that a number exists in the evidence, not that it is attached to
the right claim; a person or a second model still reads the claims.

Usage: python check_numbers.py --evidence PATH [--evidence PATH ...] [--allow-file FILE] [--ignore-below N] DRAFT [DRAFT ...]
Evidence is any mix of files (.txt, .md, .yaml, .json, .csv, or PDFs read with pdftotext) and folders of them. Numbers are compared without
commas and without trailing zeros (33,782 = 33782; 99.40 = 99.4). Not checked: integers below --ignore-below (default 10), dates such as
2026-10-05, identifiers such as T-0067, R117H, hg19 or NC_000007.13, numbers inside URLs and DOIs, numbers with a unit attached to a word (40s, 5T),
and anything in code fences or inline code. Exit 0: all matched. 1: some did not. 2: bad input (missing file, no evidence, an allow line with no reason)."""
import argparse
import re
import sys
from pathlib import Path

TEXT_SUFFIXES = {".txt", ".md", ".yaml", ".yml", ".json", ".csv"}
_FENCE = re.compile(r"^```.*?^```", re.S | re.M)
_INLINE_CODE = re.compile(r"`[^`\n]*`")
# a standalone number: not glued to a word, a path, a date or a version
_DRAFT_NUMBER = re.compile(r"(?<![\w.\-/])(\d[\d,]*(?:\.\d+)?)(?![\w/-]|\.\d)")
_ANY_NUMBER = re.compile(r"\d[\d,]*(?:\.\d+)?")


def normalise(token):
    s = token.replace(",", "").rstrip(",")
    if "." in s:
        s = s.rstrip("0").rstrip(".")
    return s


def evidence_numbers(text):
    return {normalise(m.group(0)) for m in _ANY_NUMBER.finditer(text)}


def unmatched_numbers(draft, evidence, allow, ignore_below=10):
    """[{"number", "line"}] for every number in `draft` that is in neither `evidence` (a string) nor `allow` (a dict of number -> reason)."""
    known = evidence_numbers(evidence) | {normalise(k) for k in allow}
    stripped = _INLINE_CODE.sub(lambda m: " " * len(m.group(0)), _FENCE.sub(lambda m: "\n" * m.group(0).count("\n"), draft))
    out = []
    for line_no, line in enumerate(stripped.split("\n"), 1):
        seen = set()
        for m in _DRAFT_NUMBER.finditer(line):
            token = m.group(1)
            n = normalise(token)
            if n in seen or n in known:
                continue
            if "." not in token and "," not in token and int(n) < ignore_below:
                continue
            seen.add(n)
            out.append({"number": n, "line": line_no})
    return out


def read_evidence(paths):
    """The text of every evidence file; PDFs go through the same reader as check_source_overlap. Raises FileNotFoundError or ValueError."""
    parts, found = [], 0
    for p in paths:
        p = Path(p)
        if not p.exists():
            raise FileNotFoundError(f"no such evidence: {p}")
        files = sorted(x for x in p.rglob("*") if x.is_file()) if p.is_dir() else [p]
        for f in files:
            if f.suffix.lower() in TEXT_SUFFIXES:
                parts.append(f.read_text(encoding="utf-8", errors="replace"))
                found += 1
            elif f.suffix.lower() == ".pdf":
                import check_source_overlap as cso  # noqa: PLC0415 (same folder; needs pdftotext)
                parts.append(cso.source_text(f))
                found += 1
    if not found:
        raise ValueError("the evidence paths hold no readable files")
    return "\n".join(parts)


def read_allow(path):
    """{number: reason} from lines like `564  # 127 + 437`. A line with no reason is an error; blank lines and # lines are skipped."""
    allow = {}
    for n, raw in enumerate(Path(path).read_text(encoding="utf-8").splitlines(), 1):
        if not raw.strip() or raw.strip().startswith("#"):
            continue
        num, _, reason = raw.partition("#")
        if not num.strip() or not reason.strip():
            raise ValueError(f"{path}:{n}: an allowed number needs a reason after '#' (for example `564  # 127 + 437`)")
        allow[normalise(num.strip())] = reason.strip()
    return allow


def main(argv=None):
    ap = argparse.ArgumentParser(description="Every number in a draft must be in the evidence or allowed with a reason.")
    ap.add_argument("drafts", nargs="+")
    ap.add_argument("--evidence", action="append", required=True, metavar="PATH", help="a file or folder of evidence; repeat the flag for more")
    ap.add_argument("--allow-file", metavar="FILE")
    ap.add_argument("--ignore-below", type=int, default=10)
    a = ap.parse_args(argv)
    for stream in (sys.stdout, sys.stderr):  # a Windows console may not show every character; never crash a report over it
        if hasattr(stream, "reconfigure"):
            stream.reconfigure(errors="backslashreplace")
    try:
        evidence = read_evidence(a.evidence)
        allow = read_allow(a.allow_file) if a.allow_file else {}
        results = []
        for d in a.drafts:
            text = Path(d).read_text(encoding="utf-8", errors="replace")
            for hit in unmatched_numbers(text, evidence, allow, a.ignore_below):
                results.append((d, hit))
    except (FileNotFoundError, ValueError) as e:
        print(f"ERROR: {e}")
        return 2
    for d, h in results:
        print(f"{d}:{h['line']}: number {h['number']} is not in the evidence and is not allowed")
    if results:
        print(f"{len(results)} number(s) not found. Check each against the source; if it is arithmetic, add it to the allow file with the sum.")
        return 1
    print(f"OK: every number in {len(a.drafts)} draft(s) is in the evidence or allowed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
