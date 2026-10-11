#!/usr/bin/env python3
"""Ground every registry statement the atlas tools write (counts.json's registry_terms, the TERMS file, the checker's header) in a
claim with an exact quote from the registry's own pages.

registry_claims.json holds one claim per registry-attributed statement, in the format of tools/claims/check_claims.py (id, claim,
source, quote, kind, and scope where a widening word needs one). Sources: ctgov-terms (the Terms and Conditions page) and
ctgov-disclaimer (the Disclaimer page).

Part (a), coverage, needs no page text and runs everywhere, CI included:
  - the claims file parses; every claim has the checker's fields, a non-empty quote and a known source; ids are unique;
  - each statement constant (check_atlas.NO_WARRANTY, SPONSOR_RESPONSIBILITY, THIRD_PARTY_COPYRIGHT, KEEP_CURRENT, DISCLAIMER_URL,
    DISCLAIMER_LAST_UPDATED, TERMS_URL, and fetch_snapshot.TERMS["terms_last_updated"]) is covered by claims of its source:
      a sentence: the claim texts contained in it, plus the lab's own fragments listed in LAB_FRAGMENTS, must account for all of
        its words (or one claim contains the whole sentence);
      an address: a claim of its source contains it;
      a date YYYY-MM-DD: a claim of its source contains it written as "Month DD, YYYY" (as the registry's pages write it).
    So a changed word in any statement fails until its claim is changed and checked again.
Part (b), the quotes, needs the registry's page text, which is third-party text and is never committed. A person renders the Terms and
Disclaimer pages and saves their text outside the repository; then
    python check_registry_claims.py --terms PATH --disclaimer PATH
runs tools/claims/check_claims.py on the claims with those two files as evidence and fails on any rule that does not pass. When a saved
file starts with a line "URL: <address>", that address must be the one the tools emit. CI has no saved pages, so part (b) is skipped.

What this does NOT prove: that a quote entails its sentence. A person or a different model still reads each quote against its claim.

Exit 0 all checks run passed; 1 a check failed; 2 usage error (one page file without the other, or an unreadable file).
"""
from __future__ import annotations

import argparse
import contextlib
import io
import json
import re
import sys
from datetime import date
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import check_atlas as ca  # noqa: E402
import fetch_snapshot as fs  # noqa: E402

REPO_ROOT = HERE.parents[1]
CLAIMS_CHECKER = REPO_ROOT / "tools" / "claims" / "check_claims.py"
CLAIMS_PATH = HERE / "registry_claims.json"
SOURCES = ("ctgov-terms", "ctgov-disclaimer")
REQUIRED = ("id", "claim", "source", "quote", "kind")
OPTIONAL = ("computed", "scope", "basis")
KINDS = ("observed", "computed", "inferred")
# Words in the emitted sentences that are the lab's own framing, not registry statements. Each must appear in its sentence as written.
LAB_FRAGMENTS = {
    "SPONSOR_RESPONSIBILITY": ["and that", "See the registry's Disclaimer."],
    "KEEP_CURRENT": ["This copy is dated and may be out of date; the live record is the current one."],
}


def statements() -> list[tuple[str, str, str, str]]:
    """(name, kind, value, source) for every registry statement the tools emit. Read at call time, so a changed constant is seen."""
    return [
        ("NO_WARRANTY", "sentence", ca.NO_WARRANTY, "ctgov-terms"),
        ("SPONSOR_RESPONSIBILITY", "sentence", ca.SPONSOR_RESPONSIBILITY, "ctgov-disclaimer"),
        ("THIRD_PARTY_COPYRIGHT", "sentence", ca.THIRD_PARTY_COPYRIGHT, "ctgov-terms"),
        ("KEEP_CURRENT", "sentence", ca.KEEP_CURRENT, "ctgov-terms"),
        ("TERMS_URL", "address", ca.TERMS_URL, "ctgov-terms"),
        ("DISCLAIMER_URL", "address", ca.DISCLAIMER_URL, "ctgov-disclaimer"),
        ("TERMS last-updated", "date", fs.TERMS.get("terms_last_updated"), "ctgov-terms"),
        ("DISCLAIMER_LAST_UPDATED", "date", ca.DISCLAIMER_LAST_UPDATED, "ctgov-disclaimer"),
    ]


def load_claims(path) -> list[dict]:
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(data, list):
        raise ValueError("the claims file must hold a JSON list")
    return data


def _quotes(claim) -> list:
    q = claim.get("quote")
    return [q] if isinstance(q, str) else (q if isinstance(q, list) else [])


def shape_problems(claims) -> list[str]:
    problems, ids = [], []
    if not claims:
        return ["the claims file holds no claims"]
    for n, c in enumerate(claims, 1):
        if not isinstance(c, dict):
            problems.append(f"claim #{n} is not an object")
            continue
        cid = c.get("id", f"#{n}")
        ids.append(cid)
        missing = [f for f in REQUIRED if f not in c]
        unknown = sorted(set(c) - set(REQUIRED) - set(OPTIONAL))
        if missing:
            problems.append(f"claim {cid}: missing {', '.join(missing)}")
        if unknown:
            problems.append(f"claim {cid}: unknown field(s) {', '.join(unknown)}")
        quotes = _quotes(c)
        if not quotes or not all(isinstance(q, str) and q.strip() for q in quotes):
            problems.append(f"claim {cid}: no quote (a non-empty string or a list of them is needed)")
        if c.get("source") not in SOURCES:
            problems.append(f"claim {cid}: source must be one of {', '.join(SOURCES)}")
        if c.get("kind") not in KINDS:
            problems.append(f"claim {cid}: kind must be one of {', '.join(KINDS)}")
        if not isinstance(c.get("claim"), str) or not c.get("claim", "").strip():
            problems.append(f"claim {cid}: the claim text is empty")
    dupes = sorted({i for i in ids if ids.count(i) > 1})
    if dupes:
        problems.append(f"claim id(s) used twice: {', '.join(map(str, dupes))}")
    return problems


def _long_dates(iso: str) -> list[str]:
    d = date.fromisoformat(iso)
    month = d.strftime("%B")
    return [f"{month} {d.day:02d}, {d.year}", f"{month} {d.day}, {d.year}"]


def coverage_problems(claims) -> list[str]:
    """Why some emitted statement is not covered by the claims (empty when every one is)."""
    problems = shape_problems(claims)
    if problems:
        return problems
    for name, kind, value, source in statements():
        texts = sorted((c["claim"] for c in claims if c["source"] == source), key=len, reverse=True)
        if not isinstance(value, str) or not value:
            problems.append(f"{name}: no value to check")
            continue
        if kind == "address":
            # the whole address, not a prefix of a longer one: it must not run on into more address characters
            if not any(re.search(re.escape(value) + r"(?![\w/.~%-])", t) for t in texts):
                problems.append(f"{name}: no {source} claim names the address {value}")
        elif kind == "date":
            try:
                forms = _long_dates(value)
            except ValueError:
                problems.append(f"{name}: {value!r} is not a YYYY-MM-DD date")
                continue
            if not any(f in t for t in texts for f in forms):
                problems.append(f"{name}: no {source} claim says the date {value} as {forms[0]!r}")
        else:
            if any(value in t for t in texts):
                continue
            rest = value
            for t in texts:
                if t in rest:
                    rest = rest.replace(t, " ")
            for frag in LAB_FRAGMENTS.get(name, []):
                if frag not in value:
                    problems.append(f"{name}: the lab's own fragment {frag!r} is no longer in the sentence; update LAB_FRAGMENTS")
                rest = rest.replace(frag, " ")
            left = re.sub(r"[\W_]+", " ", rest).strip()
            if left:
                problems.append(f"{name}: words not covered by any {source} claim: {ca.snap.clean(left, 120)!r}")
    return problems


def saved_address_problems(path: Path, expected: str) -> list[str]:
    first = path.read_text(encoding="utf-8-sig").splitlines()[:1]
    if first and first[0].startswith("URL:"):
        named = first[0][4:].strip()
        if named != expected:
            return [f"{path.name} names another address ({ca.snap.clean(named, 80)}) than the one the tools emit ({expected})"]
    return []


def run_claims_checker(claims_path: Path, terms: Path, disclaimer: Path) -> tuple[int, str]:
    """Run tools/claims/check_claims.py (found from the repository root) and return its exit code and output."""
    if not CLAIMS_CHECKER.is_file():
        return 2, f"the claims checker is not at {CLAIMS_CHECKER.relative_to(REPO_ROOT).as_posix()}"
    sys.path.insert(0, str(CLAIMS_CHECKER.parent))
    try:
        import check_claims  # noqa: PLC0415
    finally:
        sys.path.remove(str(CLAIMS_CHECKER.parent))
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        code = check_claims.main([str(claims_path), "--evidence", f"ctgov-terms={terms}", "--evidence", f"ctgov-disclaimer={disclaimer}"])
    return code, out.getvalue()


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Ground the registry statements the atlas tools emit in exact quotes (see the module docstring).")
    ap.add_argument("--claims", type=Path, default=CLAIMS_PATH)
    ap.add_argument("--terms", type=Path, help="saved text of the registry's Terms and Conditions page (kept outside the repository)")
    ap.add_argument("--disclaimer", type=Path, help="saved text of the registry's Disclaimer page (kept outside the repository)")
    a = ap.parse_args(argv)
    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, "reconfigure"):
            stream.reconfigure(errors="backslashreplace")
    if bool(a.terms) != bool(a.disclaimer):
        print("ERROR: --terms and --disclaimer go together")
        return 2
    try:
        claims = load_claims(a.claims)
    except (OSError, ValueError) as exc:
        print(f"FAIL part (a): the claims file could not be read: {ca.snap.clean(exc, 200)}")
        return 1
    problems = coverage_problems(claims)
    for p in problems:
        print(f"FAIL part (a): {p}")
    if not problems:
        print(f"PASS part (a): {len(claims)} claims cover every registry statement the tools emit")
    if not a.terms:
        print("part (b) skipped: no saved page text given (--terms and --disclaimer); CI never has it")
        return 1 if problems else 0
    for p in (a.terms, a.disclaimer):
        if not p.is_file():
            print(f"ERROR: no such file: {p.name}")
            return 2
    failed = bool(problems)
    for p in saved_address_problems(a.terms, ca.TERMS_URL) + saved_address_problems(a.disclaimer, ca.DISCLAIMER_URL):
        print(f"FAIL part (b): {p}")
        failed = True
    code, out = run_claims_checker(a.claims, a.terms, a.disclaimer)
    fails = [ln for ln in out.splitlines() if ln.startswith("FAIL")]
    for ln in fails:
        print(ca.snap.clean(ln, 300))
    summary = next((ln for ln in out.splitlines() if re.match(r"^\d+ claims, \d+ failures$", ln)), "no summary line")
    print(f"part (b): {summary}")
    if code == 2:
        print("ERROR: the claims checker reported a usage error")
        return 2
    if code != 0 or fails:
        failed = True
    print("A pass does not prove that a quote supports its claim; a person or a different model still reads each quote against it.")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
