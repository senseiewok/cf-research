#!/usr/bin/env python3
"""The versioned rule lexicon for the trial atlas (lexicon.json) and the deterministic tagger that applies it.

A tag records the class, the rule id and the exact matched span (field, start, end, text): entry[field][start:end] == text, always.
Rules run on the measure first and on the description only when the measure gave no class (matching settings in lexicon.json).
A rule with suppressed_by is dropped when any listed class or rule id also matched in the same field group. The entry's status:
  rule          at least one class rule matched
  not_stated    no primary outcome (NS-00), a placeholder measure (NS-01), or only wording that names no measure (NS-V1, flag vague)
  unclassified  nothing matched: the remainder a model may later propose a class for, with an exact quote (check_atlas.py checks it)
Flags: composite (COMP-01 wording; the component classes are the entry's class tags), multi_class (two or more classes without
composite wording), vague (a vague rule decided the tag; vague wording is never read as FEV1), no_primary_outcome.

Usage:
    python lexicon.py try "Change in ppFEV1 from baseline"      show which rules match a piece of text
    python lexicon.py tag SNAPSHOT_DIR --out tags.json [--exclude FILE.json]
                                                             tag every outcome entry of the studies in scope (main route)
Standard library only. No network. No model.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import scope as scope_mod  # noqa: E402
import snapshot as snap  # noqa: E402

LEXICON_PATH = HERE / "lexicon.json"


def widen_spaces(pattern: str) -> str:
    """Make a pattern tolerate any whitespace where it has a literal space: ' ' becomes \\s+, or \\s inside a character class, inside a
    lookbehind (which must keep a fixed width) and before a quantifier (so ' ?' stays 'one optional whitespace'). Escapes are kept."""
    out, i, n, in_class, groups = [], 0, len(pattern), False, []
    class_start, class_has_space = 0, False
    while i < n:
        c = pattern[i]
        if c == "\\":
            out.append(pattern[i:i + 2])
            i += 2
            continue
        if in_class:
            if c == " ":
                out.append(r"\s")
                class_has_space = True
            elif c == "]":
                out.append(c)
                in_class = False
                nxt = pattern[i + 1] if i + 1 < n else ""
                if class_has_space and not any(groups) and nxt not in "+{":
                    # a class such as [- ] takes a run of whitespace: [-\s]+ ; an optional one becomes (?:[-\s]+)?
                    if nxt in "?*":
                        out[class_start:] = ["(?:", *out[class_start:], "+)"]
                    else:
                        out.append("+")
            else:
                out.append(c)
            i += 1
            continue
        if c == "[":
            class_start, class_has_space = len(out), False
            out.append(c)
            i += 1
            if i < n and pattern[i] == "^":
                out.append("^")
                i += 1
            if i < n and pattern[i] == "]":
                out.append("]")
                i += 1
            in_class = True
            continue
        if c == "(":
            groups.append(pattern.startswith("(?<=", i) or pattern.startswith("(?<!", i))
        elif c == ")" and groups:
            groups.pop()
        if c == " ":
            nxt = pattern[i + 1] if i + 1 < n else ""
            out.append(r"\s" if any(groups) or nxt in "?*+{" else r"\s+")
        else:
            out.append(c)
        i += 1
    return "".join(out)


def _compile(pattern: str, case_sensitive: bool = False):
    return re.compile(widen_spaces(pattern), 0 if case_sensitive else re.I)


class Lexicon:
    def __init__(self, data: dict):
        self.data = data
        self.version = data["lexicon_version"]
        self.match_fields = data["matching"]["match_fields"]
        self.fallback_fields = data["matching"]["fallback_fields"]
        self.fallback_excluded = set(data["matching"].get("fallback_excluded_rules", []))
        self.domain_of, self.class_order, self.labels, self.cf_specific = {}, [], {}, set()
        self.control_classes = set()
        for d in data["domains"]:
            for c in d["classes"]:
                self.domain_of[c["id"]] = d["id"]
                self.class_order.append(c["id"])
                self.labels[c["id"]] = c["label"]
                if c.get("cf_specific"):
                    self.cf_specific.add(c["id"])
                if c.get("cf_specific") or c.get("negative_control"):
                    self.control_classes.add(c["id"])      # the classes the non-CF negative control must find (near) zero of
        self.rules = []
        for r in data["rules"]:
            if r["class"] not in self.domain_of:
                raise ValueError(f"rule {r['id']}: unknown class {r['class']}")
            self.rules.append({**r, "re": _compile(r["pattern"], r.get("case_sensitive", False))})
        ids = [r["id"] for r in self.rules]
        if len(ids) != len(set(ids)):
            raise ValueError("duplicate rule ids in the lexicon")
        known = set(ids) | set(self.domain_of)
        for r in self.rules:
            bad = [s for s in r.get("suppressed_by", []) if s not in known]
            if bad:
                raise ValueError(f"rule {r['id']}: suppressed_by names unknown id(s) {bad}")
        bad = self.fallback_excluded - set(ids)
        if bad:
            raise ValueError(f"fallback_excluded_rules names unknown rule id(s) {sorted(bad)}")
        self.composite = [{**r, "re": _compile(r["pattern"])} for r in data["composite_rules"]]
        ns = {r["id"]: r for r in data["not_stated_rules"]}
        self.placeholder = _compile(ns["NS-01"]["pattern"])
        self.vague_ns = _compile(ns["NS-V1"]["pattern"])
        tf = data["timeframe"]
        self.tf = tf
        units = "|".join(sorted(map(re.escape, tf["unit_words"]), key=len, reverse=True))
        words = "|".join(sorted(map(re.escape, tf["number_words"]), key=len, reverse=True))
        self._num = rf"(?:\d+(?:\.\d+)?|\b(?:{words})\b)"
        sep = "\\s*(?:-|\u2013|\u2014|to|through|thru|and|or|&|,)\\s*"
        self._num_re = re.compile(self._num, re.I)
        # "24 weeks", "24-week", "six months", "30 minutes", "1 d"
        self._tf_after = re.compile(rf"({self._num})\s*(?:-\s*)?({units})\b", re.I)
        # "Week 24", "Weeks 4-24", "Week 4, 8, 12, 16", "Years 1 to 5", "Days 1 through 29": every number after the unit is read
        self._tf_before = re.compile(rf"\b({units})\b\.?\s*({self._num}(?:{sep}{self._num})*)", re.I)

    @property
    def sha256(self) -> str:
        """Hash of the parsed lexicon in a canonical form, so a line-ending change on checkout does not change it."""
        return hashlib.sha256(json.dumps(self.data, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")).hexdigest()

    # ------------------------------------------------------------- matching

    def rule_matches(self, rule_id: str, text: str):
        rule = next(r for r in self.rules if r["id"] == rule_id)
        return rule["re"].search(text)

    def _match_group(self, entry: dict, fields: list[str]) -> list[dict]:
        hits = []
        for f in fields:
            text = entry.get(f) or ""
            for r in self.rules:
                m = r["re"].search(text)
                if m and m.group(0).strip():
                    hits.append({"rule": r, "field": f, "start": m.start(), "end": m.end(), "text": m.group(0)})
        matched_ids = {h["rule"]["id"] for h in hits}
        kept = []
        for h in hits:
            sup = set(h["rule"].get("suppressed_by", []))
            others = {x["rule"]["class"] for x in hits if x is not h} | (matched_ids - {h["rule"]["id"]})
            if sup & others:
                continue
            if h["rule"].get("suppressed_by_any_other_class") and any(
                    x["rule"]["class"] != h["rule"]["class"] and not x["rule"].get("suppressed_by_any_other_class") for x in hits):
                continue
            kept.append(h)
        return kept

    def _fallback(self, entry: dict) -> list[dict]:
        """From the description: generic safety wording is ignored, and only the earliest remaining match is kept."""
        hits = [h for h in self._match_group(entry, self.fallback_fields) if h["rule"]["id"] not in self.fallback_excluded]
        if not hits:
            return []
        pool = [h for h in hits if not h["rule"].get("vague")] or hits
        first = min(pool, key=lambda h: (self.fallback_fields.index(h["field"]), h["start"], self.class_order.index(h["rule"]["class"])))
        return [{**first, "from_description": True}]

    @staticmethod
    def _tag(h: dict) -> dict:
        return {"class": h["rule"]["class"], "by": "rule", "rule_id": h["rule"]["id"], "field": h["field"], "start": h["start"],
                "end": h["end"], "text": h["text"], "from_description": bool(h.get("from_description"))}

    def tag_entry(self, entry: dict) -> dict:
        out = {"entry_id": entry["entry_id"], "nct_id": entry["nct_id"], "status": "unclassified", "tags": [], "safety_subtypes": [],
               "flags": {"composite": False, "multi_class": False, "vague": False,
                         "no_primary_outcome": bool(entry.get("no_primary_outcome"))}}
        out.update(self.timeframe(entry.get("time_frame") or ""))
        if entry.get("no_primary_outcome"):
            out["status"] = "not_stated"
            out["tags"] = [{"class": "not_stated", "by": "rule", "rule_id": "NS-00", "field": None, "start": None, "end": None, "text": None,
                            "from_description": False}]
            return out
        hits = self._match_group(entry, self.match_fields)
        if not any(not h["rule"].get("vague") for h in hits):
            # The measure gave no class, or only vague wording: the description is read. A measure stated there (FEV1, say) is
            # used as stated; otherwise the vague tag stands. Nothing is inferred from the vague words themselves.
            fallback = self._fallback(entry)
            if any(not h["rule"].get("vague") for h in fallback) or not hits:
                hits = fallback or hits
        if hits:
            out["status"] = "rule"
            seen = {}
            for h in hits:                                   # first hit per class, in field then rule order, is the class's tag
                seen.setdefault(h["rule"]["class"], h)
            ordered = sorted(seen.values(), key=lambda h: self.class_order.index(h["rule"]["class"]))
            out["tags"] = [self._tag(h) for h in ordered]
            out["safety_subtypes"] = sorted({h["rule"]["subtype"] for h in hits if h["rule"].get("subtype")},
                                            key=self.data["safety_subtypes"].index)
            out["flags"]["vague"] = any(seen[c]["rule"].get("vague") for c in seen)
            comp = None
            for f in self.match_fields + self.fallback_fields:
                for r in self.composite:
                    m = r["re"].search(entry.get(f) or "")
                    if m and comp is None:
                        comp = {"rule_id": r["id"], "field": f, "start": m.start(), "end": m.end(), "text": m.group(0)}
            if comp:
                out["flags"]["composite"] = True
                out["composite_span"] = comp
            elif len(seen) > 1:
                out["flags"]["multi_class"] = True
            return out
        measure = entry.get("measure") or ""
        m = self.placeholder.match(measure)
        if m:
            out["status"] = "not_stated"
            out["tags"] = [{"class": "not_stated", "by": "rule", "rule_id": "NS-01", "field": "measure", "start": m.start(),
                            "end": m.end(), "text": m.group(0), "from_description": False}]
            return out
        for f in self.match_fields + self.fallback_fields:
            m = self.vague_ns.search(entry.get(f) or "")
            if m:
                out["status"] = "not_stated"
                out["flags"]["vague"] = True
                out["tags"] = [{"class": "not_stated", "by": "rule", "rule_id": "NS-V1", "field": f, "start": m.start(), "end": m.end(),
                                "text": m.group(0), "from_description": f != self.match_fields[0]}]
                return out
        return out

    # ------------------------------------------------------------- time frame

    def timeframe(self, text: str) -> dict:
        if not text.strip():
            return {"timeframe_bucket": self.tf["empty_bucket"], "timeframe_days": None}
        days = []
        for m in self._tf_after.finditer(text):
            days.append(self._days(m.group(1), m.group(2)))
        for m in self._tf_before.finditer(text):
            days.extend(self._days(n, m.group(1)) for n in self._num_re.findall(m.group(2)))
        if not days:
            return {"timeframe_bucket": self.tf["unparseable_bucket"], "timeframe_days": None}
        longest = max(days)
        for b in self.tf["buckets"]:
            if b["max_days"] is None or longest <= b["max_days"]:
                return {"timeframe_bucket": b["id"], "timeframe_days": round(longest, 3)}
        raise AssertionError("the last bucket must have max_days null")

    def _days(self, number: str, unit: str) -> float:
        words, units = self.tf["number_words"], self.tf["unit_words"]
        n = float(words[number.lower()]) if number.lower() in words else float(number)
        return n * float(self.tf["units_in_days"][units[unit.lower()]])


def load(path: Path | str = LEXICON_PATH) -> Lexicon:
    return Lexicon(json.loads(Path(path).read_text(encoding="utf-8")))


def in_scope_entries(snapshot: snap.Snapshot, manual: dict | None = None) -> tuple[list[dict], dict]:
    """The outcome entries of the main route's in-scope studies, and the scope result."""
    result = scope_mod.apply_scope(snapshot.records(scope_mod.MAIN_ROUTE), manual)
    entries = [e for r in result["included"] for e in snap.outcome_entries(r)]
    return entries, result


def tag_snapshot(snapshot: snap.Snapshot, lex: Lexicon, manual: dict | None = None) -> dict:
    entries, _ = in_scope_entries(snapshot, manual)
    return {"kind": "trial-atlas rule tags", "lexicon_version": lex.version, "lexicon_sha256": lex.sha256,
            "snapshot_sha256": snapshot.digest, "entries": [lex.tag_entry(e) for e in entries]}


def tag_all_entries(snapshot: snap.Snapshot, lex: Lexicon, route: str = scope_mod.MAIN_ROUTE) -> list[dict]:
    """Every entry of a route, scope ignored: for the negative control over a non-CF snapshot."""
    return [lex.tag_entry(e) for r in snapshot.records(route) for e in snap.outcome_entries(r)]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Apply the atlas rule lexicon.")
    sub = ap.add_subparsers(dest="cmd", required=True)
    t = sub.add_parser("try", help="show the rules that match a piece of text (as a measure)")
    t.add_argument("text")
    t.add_argument("--time-frame", default="")
    g = sub.add_parser("tag", help="tag the in-scope entries of a snapshot")
    g.add_argument("snapshot", type=Path)
    g.add_argument("--out", type=Path, required=True)
    g.add_argument("--exclude", type=Path)
    a = ap.parse_args(argv)
    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, "reconfigure"):
            stream.reconfigure(errors="backslashreplace")
    lex = load()
    if a.cmd == "try":
        res = lex.tag_entry({"entry_id": "TRY:P1", "nct_id": "TRY", "measure": a.text, "description": "", "time_frame": a.time_frame})
        print(json.dumps(res, indent=2, ensure_ascii=False))
        return 0
    try:
        s = snap.load(a.snapshot)
        manual = scope_mod.read_id_reasons(a.exclude) if a.exclude else {}
    except (snap.SnapshotError, OSError, ValueError) as exc:
        print(f"ERROR: {exc}")
        return 1
    tags = tag_snapshot(s, lex, manual)
    a.out.write_text(json.dumps(tags, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
    n = len(tags["entries"])
    by = {k: sum(1 for e in tags["entries"] if e["status"] == k) for k in ("rule", "not_stated", "unclassified")}
    print(f"{n} entries tagged with lexicon {lex.version}: by rule {by['rule']}, not stated {by['not_stated']}, "
          f"unclassified {by['unclassified']}  ->  {a.out.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
