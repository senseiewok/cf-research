#!/usr/bin/env python3
"""The versioned rule lexicon for the trial atlas (lexicon.json) and the deterministic tagger that applies it.

A tag records the class, the rule id and the exact matched span (field, start, end, text): entry[field][start:end] == text, always.
Rules run on the measure first and on the description only when the measure gave no class (matching settings in lexicon.json),
with one declared exception (matching.description_upgrades): a measure that gave only an unnamed patient-reported class may have the
instrument named in its description, and only the named-instrument rules listed there are read on the description for it.
A rule with suppressed_by is dropped when any listed class or rule id also matched in the same field group. The entry's status:
  rule          at least one class rule matched
  not_stated    no primary outcome (NS-00), a placeholder measure (NS-01), or only wording that names no measure (NS-V1, flag vague)
  unclassified  nothing matched: the remainder a model may later propose a class for, with an exact quote (check_atlas.py checks it)
Flags: composite (COMP-01 wording; the component classes are the entry's class tags), multi_class (two or more classes without
composite wording), vague (a vague rule decided the tag; vague wording is never read as FEV1), no_primary_outcome.

Usage:
    python lexicon.py try "Change in ppFEV1 from baseline"      show which rules match a piece of text
    python lexicon.py tag SNAPSHOT_DIR --out tags.json [--exclude FILE.json]
                                                             tag every outcome entry of the studies in scope (main retrieval
                                                             route); the output file must not exist yet (exit 2 otherwise)
    python lexicon.py classes                                every class id, its domain and a one-line gloss for the labeller
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

# The 8 canonical outcome classes per PRISMA 2020 CF trial census consensus
CANONICAL_8_CLASSES: tuple[str, ...] = (
    "lung_function_fev1",
    "sweat_chloride",
    "pulmonary_exacerbations",
    "quality_of_life",
    "microbiology_cfu",
    "nutritional_bmi",
    "safety_adverse_events",
    "other_or_not_stated",
)

CANONICAL_8_MAP: dict[str, str] = {
    "fev1": "lung_function_fev1",
    "sweat_chloride": "sweat_chloride",
    "exacerbations": "pulmonary_exacerbations",
    "cfqr": "quality_of_life",
    "pro_named_other": "quality_of_life",
    "pro_unnamed": "quality_of_life",
    "microbiology_sputum": "microbiology_cfu",
    "nutrition_growth": "nutritional_bmi",
    "safety_tolerability": "safety_adverse_events",
}


def to_canonical_8(cls: str) -> str:
    """Map a granular trial atlas class to one of the 8 canonical consensus outcome classes."""
    return CANONICAL_8_MAP.get(cls, "other_or_not_stated")


def map_tags_to_canonical_8(tags: list[dict]) -> list[str]:
    """Return deduplicated, ordered list of 8 canonical classes for a list of tag dictionaries."""
    if not tags:
        return ["other_or_not_stated"]
    canonical: list[str] = []
    for t in tags:
        c = to_canonical_8(t.get("class", ""))
        if c not in canonical:
            canonical.append(c)
    return canonical or ["other_or_not_stated"]


def canonical_8_spans(entry: dict, tags: list[dict]) -> list[dict]:
    """Generate verbatim character span tuples [start, end] and verify 100% exact substring match."""
    spans = []
    for t in tags:
        field = t.get("field")
        start = t.get("start")
        end = t.get("end")
        text = t.get("text")
        c_cls = to_canonical_8(t.get("class", ""))
        span_tuple = [start, end] if start is not None and end is not None else None

        # Verify verbatim substring invariant when field, span, and text are present
        if field and start is not None and end is not None and text is not None:
            raw_field = entry.get(field) or ""
            extracted = raw_field[start:end]
            if extracted != text:
                raise ValueError(
                    f"Verbatim span mismatch for rule {t.get('rule_id')}: "
                    f"raw[{start}:{end}]='{extracted}' != '{text}'"
                )
        spans.append({
            "canonical_class": c_cls,
            "granular_class": t.get("class"),
            "rule_id": t.get("rule_id"),
            "field": field,
            "start": start,
            "end": end,
            "span": span_tuple,
            "text": text,
            "from_description": bool(t.get("from_description")),
        })
    return spans


_WS_RUN = re.compile(r"\s+")
# One character to one character, so spans still map back exactly: subscript digits to ASCII digits, and the hyphen-like dashes
# (hyphen, non-breaking hyphen, figure dash, en dash, em dash, minus sign) to '-'.
_FOLD = str.maketrans({**{chr(0x2080 + i): str(i) for i in range(10)},
                       **{c: "-" for c in "\u2010\u2011\u2012\u2013\u2014\u2212"}})
MAX_TIMEFRAME_CHARS = 2000    # a longer time frame is reported "unparseable" instead of being parsed


def normalise_whitespace(text: str) -> tuple[str, list[tuple[int, int]]]:
    """(text with every run of whitespace, newlines included, collapsed to one space, and subscript digits and dashes folded one to one
    (_FOLD); for each character of it, the (start, end) it covers in the original). Rules match the normalised text, and spans are
    mapped back, so a reported span is always an exact piece of the ORIGINAL text, double spaces, line breaks and subscripts included."""
    norm, spans, pos = [], [], 0
    for m in _WS_RUN.finditer(text):
        for i in range(pos, m.start()):
            norm.append(text[i])
            spans.append((i, i + 1))
        norm.append(" ")
        spans.append((m.start(), m.end()))
        pos = m.end()
    for i in range(pos, len(text)):
        norm.append(text[i])
        spans.append((i, i + 1))
    return "".join(norm).translate(_FOLD), spans


def _mapped(m, spans, text):
    """A match on the collapsed text as (start, end, exact original text), or None for an empty match."""
    if m is None or m.end() <= m.start():
        return None
    start, end = spans[m.start()][0], spans[m.end() - 1][1]
    return start, end, text[start:end]


_VERBOSE_FLAG = re.compile(r"\(\?[a-zA-Z]*x[a-zA-Z]*[):]")


def _has_quantifier(inner: str) -> bool:
    i, in_class = 0, False
    while i < len(inner):
        c = inner[i]
        if c == "\\":
            i += 2
            continue
        if in_class:
            in_class = c != "]"
        elif c == "[":
            in_class = True
            if inner[i + 1:i + 2] == "^":
                i += 1
            if inner[i + 1:i + 2] == "]":
                i += 1
        elif c in "*+{":
            return True
        i += 1
    return False


def risky_pattern(pattern: str) -> str | None:
    """Why a pattern is refused at load, or None. Refused: an inline verbose flag such as (?x), which changes what a space means; an
    escaped space; and a quantified group whose content is itself quantified, such as (a[bc]*)* or (x+)+, which can take time that
    grows exponentially with the text (catastrophic backtracking)."""
    if _VERBOSE_FLAG.search(pattern):
        return "an inline verbose flag such as (?x)"
    stack, i, n, in_class = [], 0, len(pattern), False
    while i < n:
        c = pattern[i]
        if c == "\\":
            if pattern[i + 1:i + 2] == " ":
                return "an escaped space"
            i += 2
            continue
        if in_class:
            in_class = c != "]"
        elif c == "[":
            in_class = True
            if pattern[i + 1:i + 2] == "^":
                i += 1
            if pattern[i + 1:i + 2] == "]":
                i += 1
        elif c == "(":
            stack.append(i)
        elif c == ")" and stack:
            start = stack.pop()
            if pattern[i + 1:i + 2] in ("*", "+", "{") and _has_quantifier(pattern[start + 1:i]):
                return "a nested quantifier (a quantified group whose content is quantified)"
        i += 1
    return None


def _compile(pattern: str, case_sensitive: bool = False, rule_id: str = "?"):
    why = risky_pattern(pattern)
    if why:
        raise ValueError(f"rule {rule_id}: the pattern is refused: {why}")
    return re.compile(pattern, 0 if case_sensitive else re.I)


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
            rule = {**r, "re": _compile(r["pattern"], r.get("case_sensitive", False), r["id"])}
            if r.get("suppressed_by_pattern"):
                rule["sup_re"] = _compile(r["suppressed_by_pattern"], False, r["id"])
            self.rules.append(rule)
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
        self.fallback_restrictions = []
        for f in data["matching"].get("fallback_restrictions", []):
            bad = set(f["classes"]) - set(self.domain_of)
            if bad:
                raise ValueError(f"fallback restriction {f['id']}: unknown class(es) {sorted(bad)}")
            self.fallback_restrictions.append({**f, "re": _compile(f["measure_pattern"], False, f["id"])})
        # description_upgrades: when the measure gave ONLY the listed classes, the listed rules are also read on the description; a
        # match there adds its class (from_description) and removes the classes in "replaces". Nothing else reads the description then.
        self.description_upgrades = []
        for u in data["matching"].get("description_upgrades", []):
            bad = (set(u["measure_classes"]) | set(u["replaces"])) - set(self.domain_of)
            if bad:
                raise ValueError(f"description upgrade {u['id']}: unknown class(es) {sorted(bad)}")
            bad = set(u["description_rules"]) - set(ids)
            if bad or not u["description_rules"]:
                raise ValueError(f"description upgrade {u['id']}: unknown or missing rule id(s) {sorted(bad)}")
            if not set(u["replaces"]) <= set(u["measure_classes"]):
                raise ValueError(f"description upgrade {u['id']}: replaces must be among measure_classes")
            self.description_upgrades.append(u)
        # fallback_companions: groups of classes that one description may give together (two distinct marker families named side by
        # side). When the earliest description match is in a group, the earliest match of each other class of that group is kept too.
        self.fallback_companions = []
        for g in data["matching"].get("fallback_companions", []):
            bad = set(g["classes"]) - set(self.domain_of)
            if bad or len(set(g["classes"])) < 2:
                raise ValueError(f"fallback companion group {g['id']}: needs two or more known classes (unknown: {sorted(bad)})")
            self.fallback_companions.append(g)
        self.composite = [{**r, "re": _compile(r["pattern"], False, r["id"])} for r in data["composite_rules"]]
        ns = {r["id"]: r for r in data["not_stated_rules"]}
        self.placeholder = _compile(ns["NS-01"]["pattern"], False, "NS-01")
        self.vague_ns = _compile(ns["NS-V1"]["pattern"], False, "NS-V1")
        tf = data["timeframe"]
        self.tf = tf
        units = "|".join(sorted(map(re.escape, tf["unit_words"]), key=len, reverse=True))
        words = "|".join(sorted(map(re.escape, tf["number_words"]), key=len, reverse=True))
        self._num = rf"(?:\d+(?:\.\d+)?|\b(?:{words})\b)"
        sep = "\\s*(?:-|\u2013|\u2014|to|through|thru|and|or|&|,)\\s*"
        self._num_re = re.compile(self._num, re.I)
        # "24 weeks", "24-week", "six months", "30 minutes", "1 d"
        # The lookbehind starts a number only where no digit or dot precedes it, so a long digit run is not retried at every position.
        self._tf_after = re.compile(rf"(?<![\d.])({self._num})\s*(?:-\s*)?({units})\b", re.I)
        # "Week 24", "Weeks 4-24", "Week 4, 8, 12, 16", "Years 1 to 5", "Days 1 through 29": every number after the unit is read
        self._tf_before = re.compile(rf"\b({units})\b\.?\s*({self._num}(?:{sep}{self._num})*)", re.I)

    @property
    def sha256(self) -> str:
        """Hash of the parsed lexicon in a canonical form, so a line-ending change on checkout does not change it."""
        return hashlib.sha256(json.dumps(self.data, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")).hexdigest()

    # ------------------------------------------------------------- matching

    @staticmethod
    def find(rx, text: str):
        """(start, end, exact original text) of the first match of rx in text with its whitespace collapsed, or None."""
        norm, spans = normalise_whitespace(text)
        return _mapped(rx.search(norm), spans, text)

    def rule_matches(self, rule_id: str, text: str):
        rule = next(r for r in self.rules if r["id"] == rule_id)
        return self.find(rule["re"], text)

    def _match_group(self, entry: dict, fields: list[str]) -> list[dict]:
        hits = []
        for f in fields:
            text = entry.get(f) or ""
            norm, spans = normalise_whitespace(text)
            for r in self.rules:
                found = _mapped(r["re"].search(norm), spans, text)
                if found and r.get("sup_re") and r["sup_re"].search(norm):
                    continue                                  # suppressed_by_pattern: the same field names a measure this rule must not claim
                if found and found[2].strip():
                    hits.append({"rule": r, "field": f, "start": found[0], "end": found[1], "text": found[2]})
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
        """From the description: generic safety wording is ignored, a fallback restriction whose measure pattern matches limits the
        classes the description may give, and only the earliest remaining match is kept, together with the earliest match of each
        other class in its fallback companion group (matching.fallback_companions), if any."""
        hits = [h for h in self._match_group(entry, self.fallback_fields) if h["rule"]["id"] not in self.fallback_excluded]
        measure, _ = normalise_whitespace(entry.get(self.match_fields[0]) or "")
        for f in self.fallback_restrictions:
            if f["re"].search(measure):                       # the measure names an outcome; the description may only confirm its class
                hits = [h for h in hits if h["rule"]["class"] in f["classes"]]
        if not hits:
            return []
        pool = [h for h in hits if not h["rule"].get("vague")] or hits
        order = lambda h: (self.fallback_fields.index(h["field"]), h["start"], self.class_order.index(h["rule"]["class"]))  # noqa: E731
        first = min(pool, key=order)
        kept = [first]
        for g in self.fallback_companions:                    # matching.fallback_companions: named side by side, both classes stay
            if first["rule"]["class"] in g["classes"]:
                for c in g["classes"]:
                    same = [h for h in pool if h["rule"]["class"] == c and c != first["rule"]["class"]]
                    if same and all(k["rule"]["class"] != c for k in kept):
                        kept.append(min(same, key=order))
        return [{**h, "from_description": True} for h in kept]

    def _description_upgrade(self, entry: dict, hits: list[dict]) -> list[dict]:
        """matching.description_upgrades: if every class the measure gave is in an upgrade's measure_classes, its description_rules
        are read on the description (with the usual suppression there); their matches are added as from_description tags and the
        classes in 'replaces' are dropped. It applies only when the measure gave at least one class in 'replaces'. Any other class
        in the measure, or no match in the description, leaves the hits as they are."""
        measure_classes = {h["rule"]["class"] for h in hits}
        for u in self.description_upgrades:
            if not measure_classes or not measure_classes <= set(u["measure_classes"]):
                continue
            if not measure_classes & set(u["replaces"]):          # the upgrade only ever replaces a class the measure actually gave
                continue
            allowed = set(u["description_rules"])
            found = [h for h in self._match_group(entry, self.fallback_fields) if h["rule"]["id"] in allowed]
            if found:
                kept = [h for h in hits if h["rule"]["class"] not in set(u["replaces"])]
                return kept + [{**h, "from_description": True} for h in found]
        return hits

    @staticmethod
    def _tag(h: dict) -> dict:
        return {"class": h["rule"]["class"], "by": "rule", "rule_id": h["rule"]["id"], "field": h["field"], "start": h["start"],
                "end": h["end"], "text": h["text"], "from_description": bool(h.get("from_description"))}

    def _finalize_entry(self, out: dict, entry: dict) -> dict:
        out["canonical_8_classes"] = map_tags_to_canonical_8(out.get("tags", []))
        out["canonical_8_spans"] = canonical_8_spans(entry, out.get("tags", []))
        return out

    def tag_entry(self, entry: dict) -> dict:
        out = {"entry_id": entry["entry_id"], "nct_id": entry["nct_id"], "status": "unclassified", "tags": [], "safety_subtypes": [],
               "flags": {"composite": False, "multi_class": False, "vague": False,
                         "no_primary_outcome": bool(entry.get("no_primary_outcome"))}}
        out.update(self.timeframe(entry.get("time_frame") or ""))
        if entry.get("no_primary_outcome"):
            out["status"] = "not_stated"
            out["tags"] = [{"class": "not_stated", "by": "rule", "rule_id": "NS-00", "field": None, "start": None, "end": None, "text": None,
                            "from_description": False}]
            return self._finalize_entry(out, entry)
        hits = self._match_group(entry, self.match_fields)
        if not any(not h["rule"].get("vague") for h in hits):
            # The measure gave no class, or only vague wording: the description is read. A measure stated there (FEV1, say) is
            # used as stated; otherwise the vague tag stands. Nothing is inferred from the vague words themselves.
            fallback = self._fallback(entry)
            if any(not h["rule"].get("vague") for h in fallback) or not hits:
                hits = fallback or hits
        else:
            hits = self._description_upgrade(entry, hits)
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
                    found = self.find(r["re"], entry.get(f) or "")
                    if found and comp is None:
                        comp = {"rule_id": r["id"], "field": f, "start": found[0], "end": found[1], "text": found[2]}
            if comp:
                out["flags"]["composite"] = True
                out["composite_span"] = comp
            elif len(seen) > 1:
                out["flags"]["multi_class"] = True
            return self._finalize_entry(out, entry)
        measure = entry.get("measure") or ""
        m = self.placeholder.match(measure)
        if m:
            out["status"] = "not_stated"
            out["tags"] = [{"class": "not_stated", "by": "rule", "rule_id": "NS-01", "field": "measure", "start": m.start(),
                            "end": m.end(), "text": m.group(0), "from_description": False}]
            return self._finalize_entry(out, entry)
        for f in self.match_fields + self.fallback_fields:
            found = self.find(self.vague_ns, entry.get(f) or "")
            if found:
                out["status"] = "not_stated"
                out["flags"]["vague"] = True
                out["tags"] = [{"class": "not_stated", "by": "rule", "rule_id": "NS-V1", "field": f, "start": found[0], "end": found[1],
                                "text": found[2], "from_description": f != self.match_fields[0]}]
                return self._finalize_entry(out, entry)
        return self._finalize_entry(out, entry)

    # ------------------------------------------------------------- time frame

    def timeframe(self, text: str) -> dict:
        if not text.strip():
            return {"timeframe_bucket": self.tf["empty_bucket"], "timeframe_days": None}
        if len(text) > MAX_TIMEFRAME_CHARS:
            return {"timeframe_bucket": self.tf["unparseable_bucket"], "timeframe_days": None}
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
    sub.add_parser("classes", help="list every class id with its domain and a one-line gloss (for the person labelling)")
    sub.add_parser("canonical-classes", help="list the 8 canonical outcome classes and their mapped granular classes")
    a = ap.parse_args(argv)
    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, "reconfigure"):
            stream.reconfigure(errors="backslashreplace")
    lex = load()
    if a.cmd == "try":
        res = lex.tag_entry({"entry_id": "TRY:P1", "nct_id": "TRY", "measure": a.text, "description": "", "time_frame": a.time_frame})
        print(json.dumps(res, indent=2, ensure_ascii=False))
        return 0
    if a.cmd == "classes":
        for d in lex.data["domains"]:
            for c in d["classes"]:
                print(f"{c['id']:28} {d['id']:24} {c.get('gloss', '')}")
        return 0
    if a.cmd == "canonical-classes":
        for canon in CANONICAL_8_CLASSES:
            granular = [k for k, v in CANONICAL_8_MAP.items() if v == canon]
            mapped_str = ", ".join(granular) if granular else "all other physiological/unclassified/not-stated measures"
            print(f"{canon:28} -> {mapped_str}")
        return 0
    if a.out.exists():
        print(f"ERROR: {a.out.name} exists; a tags file is never overwritten (choose a new name)", file=sys.stderr)
        return 2
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
