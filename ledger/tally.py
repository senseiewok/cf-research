#!/usr/bin/env python3
"""Tally ledger entries by first detector and by class.

Reads every ledger/*.md table row that starts with "| E<number> " and counts the
"Found by" and "Class" columns. The mapping from free text to a detector is a few
keywords, so the numbers are approximate; read the rows for anything that matters.

Usage:
    python tally.py              # every ledger/*.md next to this file
    python tally.py FILE ...     # specific files
    python tally.py --self-check # prove the keyword mapping classifies known strings
Exit 1 if a row cannot be classified or the self-check fails.
"""

from __future__ import annotations

import collections
import re
import sys
from pathlib import Path

# Order matters: the first matching keyword wins.
DETECTORS = [
    ("frontier review", r"frontier|blind (?:frontier )?review|public-readiness review|reviewer|same review"),
    ("local worker", r"local (?:worker|model)|thinking pass"),
    ("human", r"\bhuman\b|\buser\b"),
    ("script assertion or self-test", r"assertion|self-test|verifier case|new case that failed|mutation"),
    ("interpreter or tests", r"interpreter|py_compile|traceback|exception|running the tests|tests?\b|tool error|no such file|error text"),
    ("own re-read", r"\bme\b|\bmy\b|re-read|reading|after the fact|git show"),
]
CLASS = re.compile(r"^[A-D]\b")


def classify(text: str) -> str | None:
    t = text.lower()
    for name, pattern in DETECTORS:
        if re.search(pattern, t):
            return name
    return None


def rows(path: Path):
    for line in path.read_text(encoding="utf-8").splitlines():
        if re.match(r"^\| E\d+ ", line):
            cells = [c.strip() for c in line.strip().strip("|").split("|")]
            if len(cells) >= 6:
                yield cells[0], cells[2], cells[4]  # id, found by, class


def self_check() -> int:
    cases = {
        "interpreter, immediately": "interpreter or tests",
        "a blind frontier-model review": "frontier review",
        "me reading the output": "own re-read",
        "my own assertion that printed the tree listing": "script assertion or self-test",
        "`git show --stat` after the fact": "own re-read",
        "the same review": "frontier review",
        "a mutation proof run properly": "script assertion or self-test",
        "the script's own self-test": "script assertion or self-test",
        "the human, after publication": "human",
        "a local worker thinking pass": "local worker",
    }
    bad = [f"{k!r} -> {classify(k)!r}, wanted {v!r}" for k, v in cases.items() if classify(k) != v]
    return bad


def main() -> int:
    if "--self-check" in sys.argv:
        bad = self_check()
        for b in bad:
            print("FAIL", b)
        print("self-check:", "OK" if not bad else f"{len(bad)} failure(s)")
        return 1 if bad else 0
    paths = [Path(a) for a in sys.argv[1:]] or sorted(Path(__file__).resolve().parent.glob("*.md"))
    by_detector: collections.Counter[str] = collections.Counter()
    by_class: collections.Counter[str] = collections.Counter()
    unclassified = []
    total = 0
    for p in paths:
        for eid, found, cls in rows(p):
            total += 1
            d = classify(found)
            if d is None:
                unclassified.append(f"{p.name} {eid}: {found}")
            else:
                by_detector[d] += 1
            m = CLASS.match(cls)
            by_class[m.group(0) if m else "?"] += 1
    print(f"{total} entries in {len(paths)} file(s)")
    print("first detector:")
    for k, v in by_detector.most_common():
        print(f"  {v:3d}  {k}")
    print("class:", dict(sorted(by_class.items())))
    for u in unclassified:
        print("UNCLASSIFIED", u)
    return 1 if unclassified else 0


if __name__ == "__main__":
    sys.exit(main())
