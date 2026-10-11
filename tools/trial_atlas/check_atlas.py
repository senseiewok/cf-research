#!/usr/bin/env python3
"""Check the trial atlas before anything is published: recompute every count, verify model quotes, measure the taggers on the frozen
set, run the negative controls and apply the publication stop rules from the design (decision 5).

What it does, in order:
  1. Integrity. The snapshot's files match its manifest; the tags file was made from this snapshot and from the current lexicon, and
     re-running the lexicon reproduces it exactly; the counts recomputed from snapshot + tags + accepted model tags equal the committed
     counts file. Any difference is drift and fails the run (exit 1), with each differing count printed.
  2. Model tags (optional file of proposed tags {entry_id, class, quote}). A proposed tag is rejected unless its entry is in scope and
     was left unclassified by the rules, its class is a lexicon class other than not_stated, and its quote is an exact substring
     (case and whitespace included) of that entry's measure, description or time frame. Only accepted tags are counted.
  3. Metrics on the frozen set (a labelling sheet a person filled in, and its sealed key): per-class precision and recall with the
     number of checked rows and a Wilson 95% interval ("too few to estimate" under 10 positives), micro and macro F1, the shares of
     entries tagged by rule, by model, as other, as not stated and left unclassified, and the model quote-rejection rate.
  4. Negative controls (--negative-controls): canary model tags with false quotes must all be rejected; planted non-CF ids must all be
     excluded by scope; the tagger over a non-CF snapshot must give near zero CF-specific classes and near zero of the classes the
     lexicon marks negative_control (exacerbations). By default they run on the synthetic fixtures; --canary, --planted and
     --control-snapshot point them at other inputs. Synthetic controls never satisfy the publication gate. A canary counts only when
     its entry is in scope, unclassified, and the tag is rejected for its quote; any other rejection makes the canary an error.
  5. Publication stop rules (skipped with --drift-only). Publication is blocked (exit 3) when:
       S0 the snapshot or the control snapshot carries the synthetic marker (`_synthetic` on a page or in version.json, or an
          apiVersion naming SYNTHETIC); only main(allow_synthetic_for_tests=True), which no option sets, skips this rule;
       S1 the frozen set is missing or not bound to this snapshot and lexicon (key hashes, and each row's wording against this
          snapshot's text), has fewer than MIN_FROZEN_ROWS (50) labelled rows, has any unlabelled row, leaves a shown class with
          no labelled positive, or a shown class has precision below 0.85 or recall below 0.80 where estimable (10 or more);
       S2 "other" plus unclassified exceed 15% of entries (stricter than the design's "other" alone);
       S3 the two retrieval routes differ by more than 10% unexplained, or there is no recall route;
       S4 the controls were not run, ran on synthetic fixtures instead of --canary, --planted and --control-snapshot, or the control
          snapshot holds fewer than MIN_CONTROL_ENTRIES (50) outcome entries.
     A failed control or a rerun that does not reproduce counts and hashes is an integrity failure (exit 1). The thresholds are the
     reviewer's judgement, not a standard; they may be changed before the frozen set is labelled and not after.
Text that came from a model, a file or the registry is printed with control characters removed.

No model is called. No network. Every number printed comes from the files given.

Usage:
  python check_atlas.py --snapshot DIR --tags tags.json --counts counts.json [--model-tags FILE] [--exclude FILE] [--explained FILE]
                        [--frozen-sheet SHEET.csv --frozen-key KEY.json] [--negative-controls] [--drift-only]
  python check_atlas.py --snapshot DIR --tags tags.json --write-counts counts.json [...]   write the counts file (no check)
  python check_atlas.py --verify-model-tags FILE --snapshot DIR --tags tags.json          only the quote check, one line per tag
  python check_atlas.py --negative-controls [--controls-only]                            only the controls, on the synthetic fixtures
  python check_atlas.py --snapshot DIR --confirm-real-run      a person confirms a REAL snapshot's fetch configuration (config_verified)
counts.json carries a registry_terms block (source, processing and fetch dates, terms pointer, modifications, licence line, retention
note); the check prints it first and the drift check covers it.
Exit 0 all checks passed (and publication allowed unless --drift-only); 1 integrity failure or a failed control; 2 usage error;
3 integrity passed but publication is blocked; 4 the controls-only run passed but is not the gate (0 with --controls-only).
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import lexicon as lex_mod  # noqa: E402
import scope as scope_mod  # noqa: E402
import snapshot as snap  # noqa: E402

THRESHOLDS = {"precision": 0.85, "recall": 0.80, "other_share": 0.15, "route_unexplained": 0.10, "min_positives": 10,
              "near_zero_cf_specific": 0.01}
Z95 = 1.959963984540054
MIN_QUOTE_CHARS = 3
MODEL_FORBIDDEN_CLASSES = {"not_stated"}
MIN_FROZEN_ROWS = 50          # the design's number of person-labelled outcome entries (decision 5)
MIN_CONTROL_ENTRIES = 50      # the fewest outcome entries a non-CF control snapshot may hold for the gate (a reviewer's number)
TAG_STATUSES = {"rule", "not_stated", "unclassified"}
CO_OCCURRENCE_EXCLUDED = ("other", "not_stated")

# The registry's terms travel with every count (counts.json, the printed header, the TERMS file beside a labelling sheet).
REGISTRY_SOURCE = "ClinicalTrials.gov"
# The Terms and Conditions page as given in the third review round; it was not opened during this offline build.
TERMS_URL = "https://clinicaltrials.gov/about-site/terms-conditions"
LICENCE_LINE = ("The lab's licence covers its own tags, code and counts only; registry text and fields remain ClinicalTrials.gov data "
                "under its terms.")
RETENTION_NOTE = "The registry's terms apply for as long as the data are kept, in any copy, file or page made from them."
# Round 4: wording quoted or closely paraphrased from the registry's pages as read on 2026-10-10 by the controlling agent (not by
# these tools, which open no URL). The Disclaimer page said "Last updated on August 03, 2023".
NO_WARRANTY = ("ClinicalTrials.gov states that the U.S. Government makes no warranties about its data and assumes no liability for "
               "their use.")
SPONSOR_RESPONSIBILITY = ("Study sponsors and investigators write and are responsible for their own records; the U.S. Government does "
                          "not review or approve the safety and science of all studies listed. See the registry's Disclaimer.")
DISCLAIMER_URL = "https://clinicaltrials.gov/about-site/disclaimer"
DISCLAIMER_LAST_UPDATED = "2023-08-03"
THIRD_PARTY_COPYRIGHT = ("Some registry data may be subject to third-party copyright, and the data carry an international copyright "
                         "outside the United States.")
KEEP_CURRENT = ("The registry asks that data in any publication or distribution be kept current at all times; this copy is dated and "
                "the live record is current.")


def registry_terms(snapshot, lex, scope_counts: dict) -> dict:
    """The block that must accompany any count or file made from the snapshot: source, the registry's processing date, the fetch
    date, the terms pointer, every modification the lab made, the licence line and the retention note."""
    excluded = ", ".join(f"{rid} {n}" for rid, n in scope_counts.get("excluded_by_rule", {}).items())
    terms = snapshot.manifest.get("terms") if isinstance(snapshot.manifest.get("terms"), dict) else {}
    return {
        "source": REGISTRY_SOURCE,
        "data_processed_by_registry": snapshot.manifest.get("data_timestamp"),
        "snapshot_fetched_at": snapshot.manifest.get("fetched_at"),
        "terms_url": TERMS_URL,
        "terms_last_updated": terms.get("terms_last_updated"),
        "modifications": [
            f"Primary-outcome wording was classified into classes by the versioned rule lexicon {lex.version} (sha256 {lex.sha256[:12]}); "
            "a model-proposed class is counted only when its quote is an exact piece of the registry text.",
            f"Studies were left out by stated scope rules, with these counts: {excluded} (X1 to X6 are explained in "
            "tools/trial_atlas/scope.py; X6 is a person's exclusion list).",
            "No registry wording was clipped by these tools; any clipping on a published page must be added to this list.",
            "Registry pages were re-saved as UTF-8 JSON after parsing; values and key order are kept, whitespace is not.",
        ],
        "licence": LICENCE_LINE,
        "retention": RETENTION_NOTE,
        "no_warranty": NO_WARRANTY,
        "sponsor_responsibility": SPONSOR_RESPONSIBILITY,
        "disclaimer_url": DISCLAIMER_URL,
        "disclaimer_last_updated": DISCLAIMER_LAST_UPDATED,
        "third_party_copyright": THIRD_PARTY_COPYRIGHT,
        "keep_current": KEEP_CURRENT,
    }


def terms_text(block: dict) -> str:
    """The terms block as plain lines (for the printed header and the TERMS file beside a CSV). The fetch time is UTC (the manifest
    writes it with a trailing Z); the registry's processing time is printed exactly as the registry gave it."""
    lines = [f"source: {block['source']}",
             f"data processed by the registry: {block['data_processed_by_registry']} (as given by the registry)",
             f"snapshot fetched at: {block['snapshot_fetched_at']} (UTC)",
             f"terms: {block['terms_url']} (last updated {block['terms_last_updated'] or 'not recorded in the manifest'})",
             f"disclaimer: {block['disclaimer_url']} (last updated {block['disclaimer_last_updated']})",
             "modifications made by the lab:"]
    lines += [f"  - {m}" for m in block["modifications"]]
    lines += [f"no warranty: {block['no_warranty']}", f"sponsor responsibility: {block['sponsor_responsibility']}",
              f"copyright: {block['third_party_copyright']}", f"keep current: {block['keep_current']}",
              f"licence: {block['licence']}", f"retention: {block['retention']}"]
    return "\n".join(lines) + "\n"


class UsageError(Exception):
    pass


def validate_tags(tags, lex) -> dict:
    """The tags file's shape, checked before anything reads it, so a wrong file is a usage error and not a traceback."""
    def bad(why):
        raise UsageError(f"the tags file is not a tags file written by lexicon.py ({why})")
    if not isinstance(tags, dict) or not isinstance(tags.get("entries"), list):
        bad("it must be an object with an 'entries' list")
    for i, e in enumerate(tags["entries"]):
        if not isinstance(e, dict):
            bad(f"entry {i} is not an object")
        if not isinstance(e.get("entry_id"), str) or not isinstance(e.get("nct_id"), str):
            bad(f"entry {i} has no string entry_id and nct_id")
        if e.get("status") not in TAG_STATUSES or not isinstance(e.get("tags"), list) or not isinstance(e.get("flags"), dict):
            bad(f"entry {i} has no valid status, tags list or flags")
        if not all(isinstance(t, dict) and t.get("class") in lex.domain_of for t in e["tags"]):
            bad(f"entry {i} has a tag with an unknown class")
        if not isinstance(e.get("timeframe_bucket"), str) or not isinstance(e.get("safety_subtypes"), list):
            bad(f"entry {i} has no time-frame bucket or safety subtypes")
    return tags


# ------------------------------------------------------------------ small helpers

def wilson(k: int, n: int, z: float = Z95):
    """Wilson score interval for k successes in n trials; None when n is 0."""
    if n <= 0:
        return None
    p = k / n
    denom = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
    return max(0.0, centre - half), min(1.0, centre + half)


def read_json(path, what):
    try:
        return json.loads(Path(path).read_text(encoding="utf-8-sig"))
    except (OSError, ValueError) as exc:
        raise UsageError(f"cannot read the {what} {path}: {exc}") from None


def entry_texts(snapshot: snap.Snapshot, manual=None) -> dict[str, dict]:
    entries, _ = lex_mod.in_scope_entries(snapshot, manual)
    return {e["entry_id"]: e for e in entries}


# ------------------------------------------------------------------ model tags: the exact-quote verifier

def verify_model_tags(proposed, entries: dict[str, dict], rule_tags: dict[str, dict], lex) -> tuple[list[dict], list[dict]]:
    """(accepted, rejected). Each rejected item is {"tag": the proposed tag, "reason": why, "code": shape|entry|class|decided|short|
    quote}; "quote" means the tag reached the quote check and its quote was not found. An accepted tag records the field and span
    where the quote was found. A quote is accepted only as an exact substring: no case folding, no whitespace normalising."""
    accepted, rejected = [], []
    if not isinstance(proposed, list):
        return [], [{"tag": proposed, "reason": "the proposed-tags file must hold a JSON list", "code": "shape"}]
    for t in proposed:
        def no(why, code, t=t):
            rejected.append({"tag": t, "reason": why, "code": code})
        if not isinstance(t, dict) or not all(isinstance(t.get(k), str) for k in ("entry_id", "class", "quote")):
            no("a proposed tag needs string fields entry_id, class and quote", "shape")
            continue
        eid, cls, quote = t["entry_id"], t["class"], t["quote"]
        if eid not in entries:
            no(f"{eid} is not an in-scope outcome entry of this snapshot", "entry")
            continue
        if cls not in lex.domain_of or cls in MODEL_FORBIDDEN_CLASSES:
            no(f"'{cls}' is not a class a model may propose", "class")
            continue
        if rule_tags.get(eid, {}).get("status") != "unclassified":
            no(f"{eid} was already decided by the rules ({rule_tags.get(eid, {}).get('status', 'no tag')}); a model proposes only for "
               "the remainder", "decided")
            continue
        if len(quote.strip()) < MIN_QUOTE_CHARS:
            no(f"the quote is shorter than {MIN_QUOTE_CHARS} characters", "short")
            continue
        e = entries[eid]
        where = next(((f, e[f].find(quote)) for f in ("measure", "description", "time_frame") if quote in (e.get(f) or "")), None)
        if where is None:
            low = any(quote.lower() in (e.get(f) or "").lower() for f in ("measure", "description", "time_frame"))
            no("the quote matches this entry only if case is ignored" if low else
               "the quote is not an exact substring of this entry's measure, description or time frame", "quote")
            continue
        f, start = where
        accepted.append({"entry_id": eid, "class": cls, "by": "model", "rule_id": None, "field": f, "start": start,
                         "end": start + len(quote), "text": quote, "model_id": t.get("model_id")})
    return accepted, rejected


# ------------------------------------------------------------------ counts

def final_classes(tag_entry: dict, model_by_entry: dict[str, list[dict]]) -> list[str]:
    classes = [t["class"] for t in tag_entry["tags"]]
    if tag_entry["status"] == "unclassified":
        classes = list(dict.fromkeys(t["class"] for t in model_by_entry.get(tag_entry["entry_id"], [])))
    return classes


def phase_label(phases) -> str:
    return "/".join(phases) if phases else "none"


def compute_counts(snapshot: snap.Snapshot, tags: dict, accepted_model: list[dict], lex, manual=None) -> dict:
    result = scope_mod.apply_scope(snapshot.records(scope_mod.MAIN_ROUTE), manual)
    studies = {r["nct_id"]: r for r in result["included"]}
    model_by_entry: dict[str, list[dict]] = {}
    for t in accepted_model:
        model_by_entry.setdefault(t["entry_id"], []).append(t)
    order = lex.class_order
    entries_by_class = {c: 0 for c in order}
    study_classes: dict[str, set] = {n: set() for n in studies}
    status = {"rule": 0, "model": 0, "not_stated": 0, "unclassified": 0}
    flags = {"composite": 0, "multi_class": 0, "vague": 0, "no_primary_outcome": 0}
    buckets: dict[str, int] = {b["id"]: 0 for b in lex.tf["buckets"]}
    buckets.update({lex.tf["unparseable_bucket"]: 0, lex.tf["empty_bucket"]: 0})
    subtypes = {s: 0 for s in lex.data["safety_subtypes"]}
    other_entries = 0
    for e in tags["entries"]:
        classes = final_classes(e, model_by_entry)
        st = e["status"]
        if st == "unclassified" and classes:
            st = "model"
        status[st] += 1
        for c in classes:
            entries_by_class[c] += 1
            study_classes[e["nct_id"]].add(c)
        if "other" in classes:
            other_entries += 1
        for k in flags:
            flags[k] += 1 if e["flags"].get(k) else 0
        buckets[e["timeframe_bucket"]] += 1
        for s in e.get("safety_subtypes", []):
            subtypes[s] += 1
    by_year, by_phase, cf_only, sponsors, no_sponsor = {}, {}, {}, {}, {}
    by_kind, by_posted, cyp = {"actual": {}, "planned": {}}, {}, {}
    studies_by_year: dict[str, int] = {}
    for nct, r in studies.items():
        y = str(r["start_year"]) if r["start_year"] else "unknown"
        fp = str(r["first_posted_year"]) if r.get("first_posted_year") else "unknown"
        studies_by_year[y] = studies_by_year.get(y, 0) + 1
        for c in study_classes[nct]:
            by_year.setdefault(c, {})[y] = by_year.setdefault(c, {}).get(y, 0) + 1
            if r["start_kind"] in by_kind:
                by_kind[r["start_kind"]].setdefault(c, {})[y] = by_kind[r["start_kind"]].setdefault(c, {}).get(y, 0) + 1
            by_posted.setdefault(c, {})[fp] = by_posted.setdefault(c, {}).get(fp, 0) + 1
            p = phase_label(r["phases"])
            by_phase.setdefault(c, {})[p] = by_phase.setdefault(c, {}).get(p, 0) + 1
            if r["start_kind"] in ("actual", "planned"):
                cell = cyp.setdefault(c, {"actual": {}, "planned": {}})[r["start_kind"]].setdefault(y, {})
                cell[p] = cell.get(p, 0) + 1
            if r["cf_flag"] == "CF only":
                cf_only[c] = cf_only.get(c, 0) + 1
            name = (r.get("sponsor") or "").strip() if isinstance(r.get("sponsor"), str) else ""
            if name:
                sponsors.setdefault(c, set()).add(name)
            else:                                            # a missing sponsor is counted apart, never as one more sponsor
                no_sponsor[c] = no_sponsor.get(c, 0) + 1
    n = len(tags["entries"])

    def share(k):
        return round(k / n, 4) if n else None

    # Studies (not entries) that name both classes of a pair; "other" and "not_stated" are left out of pairs. Pairs are listed in the
    # lexicon's class order, only when at least one study names both; the per-class totals are the denominators for a share.
    paired = [c for c in order if c not in CO_OCCURRENCE_EXCLUDED]
    pair_counts: dict[tuple[str, str], int] = {}
    for classes in study_classes.values():
        present = [c for c in paired if c in classes]
        for i, a in enumerate(present):
            for b in present[i + 1:]:
                pair_counts[(a, b)] = pair_counts.get((a, b), 0) + 1
    co_occurrence = {
        "unit": "studies",
        "excluded_classes": list(CO_OCCURRENCE_EXCLUDED),
        "pairs": [{"a": a, "b": b, "studies": pair_counts[(a, b)]}
                  for a, b in sorted(pair_counts, key=lambda p: (order.index(p[0]), order.index(p[1])))],
        "studies_by_class": {c: sum(1 for s in study_classes.values() if c in s) for c in paired},
    }

    def nested_sorted(d):
        return {k: nested_sorted(v) if isinstance(v, dict) else v for k, v in sorted(d.items())}

    rd = scope_mod.route_difference(snapshot, manual)      # raw differences; explanations are applied only by stop rule S3
    return {
        "kind": "trial-atlas counts",
        "units": "studies are counted once per class; entries are primary-outcome entries; a study can name several classes",
        # true when the snapshot carries the synthetic marker; a consumer can refuse synthetic data without guessing
        "synthetic": bool(synthetic_marks(snapshot)),
        # entries left with no class by a rule or by an accepted model tag (shown as "left unsorted"); not the class "other",
        # which a model or a person assigns to a measure that fits none of the classes
        "unsorted_entries": status["unclassified"],
        # studies per class x start year x phase label, with ACTUAL and planned (ESTIMATED) start types kept apart
        "class_year_phase": {c: {"actual": nested_sorted(cyp.get(c, {}).get("actual", {})),
                                 "planned": nested_sorted(cyp.get(c, {}).get("planned", {}))} for c in order},
        "co_occurrence": co_occurrence,
        "snapshot_sha256": snapshot.digest,
        "data_timestamp": snapshot.manifest.get("data_timestamp"),
        "lexicon_version": lex.version,
        "lexicon_sha256": lex.sha256,
        "scope": result["counts"],
        "routes": None if rd is None else {k: rd[k] for k in ("main_in_scope", "recall_in_scope", "only_main", "only_recall")},
        "studies": len(studies),
        "entries": n,
        "entries_by_status": status,
        "entries_by_class": entries_by_class,
        "studies_by_class": {c: sum(1 for s in study_classes.values() if c in s) for c in order},
        "studies_by_class_cf_only": {c: cf_only.get(c, 0) for c in order},
        "lead_sponsors_by_class": {c: len(sponsors.get(c, ())) for c in order},
        "studies_without_lead_sponsor_by_class": {c: no_sponsor.get(c, 0) for c in order},
        "studies_by_start_year": dict(sorted(studies_by_year.items())),
        "studies_by_class_and_start_year": {c: dict(sorted(by_year.get(c, {}).items())) for c in order},
        # Split by the registry's start type: ACTUAL and ESTIMATED (planned). A start with no date or no type is in neither.
        "studies_by_class_and_actual_start_year": {c: dict(sorted(by_kind["actual"].get(c, {}).items())) for c in order},
        "studies_by_class_and_planned_start_year": {c: dict(sorted(by_kind["planned"].get(c, {}).items())) for c in order},
        # The first-posted year, for banding registration eras (design, decision 2).
        "studies_by_class_and_first_posted_year": {c: dict(sorted(by_posted.get(c, {}).items())) for c in order},
        "studies_by_class_and_phase": {c: dict(sorted(by_phase.get(c, {}).items())) for c in order},
        "registry_terms": registry_terms(snapshot, lex, result["counts"]),
        "flags": flags,
        "timeframe_buckets": buckets,
        "safety_subtypes": subtypes,
        "shares": {"by_rule": share(status["rule"]), "by_model": share(status["model"]), "other": share(other_entries),
                   "not_stated": share(status["not_stated"]), "unclassified": share(status["unclassified"])},
        "other_entries": other_entries,
    }


def diff_paths(a, b, path="") -> list[str]:
    """Every path where two JSON values differ."""
    if isinstance(a, dict) and isinstance(b, dict):
        out = []
        for k in sorted(set(a) | set(b), key=str):
            if k not in a or k not in b:
                out.append(f"{path}/{k}: {'missing in committed' if k not in b else 'not recomputed'}")
            else:
                out.extend(diff_paths(a[k], b[k], f"{path}/{k}"))
        return out
    if a != b:
        return [f"{path}: recomputed {json.dumps(a)[:80]}, committed {json.dumps(b)[:80]}"]
    return []


# ------------------------------------------------------------------ frozen set and metrics

def sheet_wording_sha(rows: list[dict]) -> str:
    """Hash of the wording columns only, in sample order, so a filled-in sheet can be matched to the blank one it came from."""
    h = hashlib.sha256()
    for r in sorted(rows, key=lambda r: r["sample_id"]):
        h.update(json.dumps([r["sample_id"], r["measure"], r["description"], r["time_frame"]], ensure_ascii=False).encode("utf-8"))
        h.update(b"\n")
    return h.hexdigest()


def read_sheet(path) -> list[dict]:
    try:
        with Path(path).open(encoding="utf-8-sig", newline="") as fh:
            rows = list(csv.DictReader(fh))
    except OSError as exc:
        raise UsageError(f"cannot read the frozen sheet {path}: {exc}") from None
    need = {"sample_id", "measure", "description", "time_frame", "label_classes"}
    if not rows or not need <= set(rows[0]):
        raise UsageError(f"the frozen sheet {path} needs the columns {', '.join(sorted(need))}")
    return rows


def frozen_binding_problems(sheet_path, key_path, snapshot, lex, entries: dict[str, dict]) -> list[str]:
    """Why the frozen set does not belong to this snapshot and lexicon (empty when it does): the key must name this snapshot's hash
    and this lexicon's hash, every key row's entry must be in scope here, and every sheet row's wording must equal this snapshot's text
    for that entry (as the sheet tool writes it). No allowance flag exists."""
    import make_labelling_sheet as mls  # noqa: PLC0415  (imported here: it imports this module)
    rows = read_sheet(sheet_path)
    key = read_json(key_path, "frozen-set key")
    problems = []
    if key.get("snapshot_sha256") != snapshot.digest:
        problems.append("the key was drawn from a different snapshot")
    if key.get("lexicon_sha256") != lex.sha256:
        problems.append("the key was drawn with a different lexicon")
    by_sample = {r.get("sample_id"): r.get("entry_id") for r in key.get("rows", []) if isinstance(r, dict)}
    differ = []
    for r in rows:
        e = entries.get(by_sample.get(r["sample_id"]))
        if e is None:
            differ.append(r["sample_id"])
            continue
        if any(r[f] != mls.safe_cell(e[f]) for f in ("measure", "description", "time_frame")):
            differ.append(r["sample_id"])
    if differ:
        problems.append(f"{len(differ)} sheet row(s) do not match this snapshot's wording for their entry: {', '.join(snap.clean(d, 10) for d in differ[:5])}")
    return problems


def load_frozen(sheet_path, key_path, lex) -> tuple[dict[str, set], dict]:
    """{entry_id: set of gold classes} for every labelled row, and a summary. Fails (UsageError) when the sheet's wording does not
    match the key's hash, a sample id is unknown, or a label is not a lexicon class."""
    rows = read_sheet(sheet_path)
    key = read_json(key_path, "frozen-set key")
    by_sample = {r["sample_id"]: r["entry_id"] for r in key.get("rows", [])}
    if sheet_wording_sha(rows) != key.get("sheet_wording_sha256"):
        raise UsageError("the sheet's wording does not match the sealed key (rows edited, added or removed)")
    gold, unlabelled = {}, 0
    for r in rows:
        if r["sample_id"] not in by_sample:
            raise UsageError(f"sample {r['sample_id']} is not in the key")
        labels = [x.strip() for x in (r.get("label_classes") or "").replace(",", ";").split(";") if x.strip()]
        if not labels:
            unlabelled += 1
            continue
        bad = [x for x in labels if x not in lex.domain_of]
        if bad:
            raise UsageError(f"sample {r['sample_id']}: unknown class(es) {bad}; use the lexicon's class ids")
        gold[by_sample[r["sample_id"]]] = set(labels)
    return gold, {"rows": len(rows), "labelled": len(gold), "unlabelled": unlabelled}


def metrics(gold: dict[str, set], predicted: dict[str, set], classes: list[str], min_pos: int) -> dict:
    per, tp_all, fp_all, fn_all, f1s = {}, 0, 0, 0, []
    for c in classes:
        tp = sum(1 for e, g in gold.items() if c in g and c in predicted.get(e, set()))
        fp = sum(1 for e, g in gold.items() if c not in g and c in predicted.get(e, set()))
        fn = sum(1 for e, g in gold.items() if c in g and c not in predicted.get(e, set()))
        tp_all, fp_all, fn_all = tp_all + tp, fp_all + fp, fn_all + fn
        prec = tp / (tp + fp) if tp + fp else None
        rec = tp / (tp + fn) if tp + fn else None
        f1 = (2 * prec * rec / (prec + rec)) if prec and rec else (0.0 if (tp + fp or tp + fn) else None)
        if f1 is not None:
            f1s.append(f1)
        per[c] = {"tp": tp, "fp": fp, "fn": fn, "positives": tp + fn, "predicted": tp + fp, "precision": prec, "recall": rec,
                  "precision_ci": wilson(tp, tp + fp), "recall_ci": wilson(tp, tp + fn), "f1": f1,
                  "precision_estimable": tp + fp >= min_pos, "recall_estimable": tp + fn >= min_pos}
    mp = tp_all / (tp_all + fp_all) if tp_all + fp_all else 0.0
    mr = tp_all / (tp_all + fn_all) if tp_all + fn_all else 0.0
    return {"per_class": per, "rows_checked": len(gold), "micro_f1": (2 * mp * mr / (mp + mr)) if mp + mr else 0.0,
            "macro_f1": (sum(f1s) / len(f1s)) if f1s else None, "macro_over_classes": len(f1s)}


# ------------------------------------------------------------------ negative controls

def run_negative_controls(lex, *, snapshot=None, canary=None, planted=None, control_snapshot=None, manual=None,
                          near_zero=THRESHOLDS["near_zero_cf_specific"]) -> list[tuple[str, bool, str]]:
    """[(name, passed, detail)]. With no inputs given, everything runs on the synthetic fixtures in a temporary folder."""
    import synthetic_fixtures as sf  # noqa: PLC0415
    results = []
    with tempfile.TemporaryDirectory() as tmp:
        if snapshot is None:
            sf.build_cf_snapshot(Path(tmp) / "cf")
            snapshot = snap.load(Path(tmp) / "cf")
            canary = sf.canary_model_tags() if canary is None else canary
            planted = sf.PLANTED_NON_CF_IDS if planted is None else planted
        if control_snapshot is None:
            sf.build_noncf_snapshot(Path(tmp) / "noncf")
            control_snapshot = snap.load(Path(tmp) / "noncf")
        entries = entry_texts(snapshot, manual)
        rule_tags = {e["entry_id"]: e for e in lex_mod.tag_snapshot(snapshot, lex, manual)["entries"]}
        if canary:
            # A canary tests the quote check: it passes only when its entry exists in scope, is unclassified, and the tag is rejected
            # for its quote. A canary rejected for any other reason (entry missing, entry already decided, bad class) is an error.
            acc, rej = verify_model_tags(canary, entries, rule_tags, lex)
            by_quote = [r for r in rej if r.get("code") == "quote"]
            wrong = [r for r in rej if r.get("code") != "quote"]
            detail = f"{len(by_quote)} of {len(canary)} rejected for their quote"
            if acc:
                detail += f"; ACCEPTED: {[snap.clean(a['entry_id'], 40) + ' ' + snap.clean(repr(a['text']), 60) for a in acc]}"
            if wrong:
                detail += ("; ERROR, not a valid canary (rejected before the quote check): "
                           + "; ".join(f"{snap.clean((r['tag'] or {}).get('entry_id', '?') if isinstance(r['tag'], dict) else '?', 40)} "
                                       f"({r.get('code')})" for r in wrong))
            results.append(("canary model tags with false quotes are all rejected", not acc and not wrong, detail))
        else:
            results.append(("canary model tags with false quotes are all rejected", False, "no canary tags were given"))
        if planted:
            # A planted id is looked for in the main retrieval route and then in the recall route: an id only the term search
            # returned is still tested against the scope rules.
            records = {}
            for route in [scope_mod.MAIN_ROUTE] + [r for r in snapshot.routes() if r != scope_mod.MAIN_ROUTE]:
                for r in snapshot.records(route):
                    records.setdefault(r["nct_id"], r)
            missing = [i for i in planted if i not in records]
            included = {r["nct_id"] for r in scope_mod.apply_scope(list(records.values()), manual)["included"]}
            leaked = [i for i in planted if i in included]
            ok = not missing and not leaked
            detail = f"{len(planted) - len(missing) - len(leaked)} of {len(planted)} planted ids excluded"
            if missing:
                detail += f"; not in the snapshot, so not tested: {missing}"
            if leaked:
                detail += f"; INCLUDED: {leaked}"
            results.append(("planted non-CF ids are all excluded by scope", ok, detail))
        else:
            results.append(("planted non-CF ids are all excluded by scope", False, "no planted ids were given"))
        control = lex_mod.tag_all_entries(control_snapshot, lex)
        watched = lex.control_classes                       # the CF-specific classes plus classes marked negative_control
        hits = [e for e in control if any(t["class"] in watched for t in e["tags"])]
        share = len(hits) / len(control) if control else None
        ok = bool(control) and share <= near_zero
        results.append((f"non-CF control: CF-specific and watched classes ({', '.join(sorted(watched - lex.cf_specific))}) "
                        f"at most {near_zero:.0%} of entries", ok,
                        f"{len(hits)} of {len(control)} entries" + (f" ({share:.1%})" if control else "; the control has no entries")
                        + (f"; tagged: {[snap.clean(e['entry_id'], 40) for e in hits]}" if hits else "")))
    return results


# ------------------------------------------------------------------ the run

def fmt(x, digits=3):
    return "n/a" if x is None else f"{x:.{digits}f}"


def synthetic_marks(snapshot) -> list[str]:
    """Where a snapshot says it is synthetic: a page or version.json carrying `_synthetic`, or an apiVersion naming SYNTHETIC."""
    marks = []
    if isinstance(snapshot.version, dict) and ("_synthetic" in snapshot.version or "SYNTHETIC" in str(snapshot.version.get("apiVersion", "")).upper()):
        marks.append("version.json")
    for route, plist in snapshot.pages.items():
        if any(isinstance(p, dict) and "_synthetic" in p for p in plist):
            marks.append(f"{route} pages")
    return marks


def load_input_snapshot(path, what: str):
    """A snapshot given as an input other than the one under check: a missing or broken one is a usage error (exit 2)."""
    try:
        return snap.load(path)
    except snap.SnapshotError as exc:
        raise UsageError(f"cannot read the {what} {path}: {exc}") from None


def confirm_real_run(path) -> int:
    """Set config_verified in a real snapshot's manifest after a person has checked its fetch configuration against the registry
    (page size, token keys, query parameters, field names). Refuses a synthetic or an already confirmed snapshot. The manifest hash
    is recomputed, so every tags file, sheet and counts file made before must be made again."""
    if not path:
        print("ERROR: --confirm-real-run needs --snapshot")
        return 2
    try:
        s = snap.load(path)
    except snap.SnapshotError as exc:
        print(f"ERROR: {snap.clean(exc, 300)}")
        return 2
    marks = synthetic_marks(s)
    if marks:
        print(f"ERROR: this snapshot carries the synthetic marker ({', '.join(marks)}); only a real snapshot can be confirmed")
        return 2
    if s.manifest.get("config_verified") is True:
        print("ERROR: this snapshot's configuration is already confirmed")
        return 2
    m = dict(s.manifest)
    m["config_verified"] = True
    m["config_verified_at"] = datetime.now(timezone.utc).strftime("%Y-%m-%d")
    m["snapshot_sha256"] = snap.manifest_digest(m)
    (Path(path) / snap.MANIFEST).write_text(json.dumps(m, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    snap.load(path)
    print(f"config_verified set; the snapshot sha256 is now {m['snapshot_sha256']}")
    print("re-run lexicon.py tag, make_labelling_sheet.py and check_atlas.py --write-counts: anything made before refers to the old hash")
    return 0


def main(argv=None, *, allow_synthetic_for_tests: bool = False) -> int:
    """allow_synthetic_for_tests exists only for the unit tests: no command-line option sets it, so synthetic data can never pass
    the gate from a command line."""
    ap = argparse.ArgumentParser(description="Check the trial atlas: drift, model quotes, metrics, negative controls, publication gate.")
    ap.add_argument("--snapshot", type=Path)
    ap.add_argument("--tags", type=Path)
    ap.add_argument("--counts", type=Path, help="the committed counts file to check against")
    ap.add_argument("--write-counts", type=Path, help="write the recomputed counts here and stop (no check)")
    ap.add_argument("--model-tags", type=Path, help="proposed model tags: a JSON list of {entry_id, class, quote}")
    ap.add_argument("--verify-model-tags", type=Path, help="only verify these proposed tags and print one line each")
    ap.add_argument("--exclude", type=Path)
    ap.add_argument("--explained", type=Path)
    ap.add_argument("--frozen-sheet", type=Path)
    ap.add_argument("--frozen-key", type=Path)
    ap.add_argument("--negative-controls", action="store_true")
    ap.add_argument("--canary", type=Path, help="canary model tags for the controls (default: the synthetic ones)")
    ap.add_argument("--planted", type=Path, help="a JSON list of planted non-CF NCT ids (default: the synthetic ones)")
    ap.add_argument("--control-snapshot", type=Path, help="a non-CF snapshot for the tagger control (default: the synthetic one)")
    ap.add_argument("--drift-only", action="store_true", help="check integrity only; do not apply the publication stop rules")
    ap.add_argument("--controls-only", action="store_true",
                    help="with --negative-controls and no --snapshot: exit 0 when the controls pass (otherwise that run exits 4: it is not the gate)")
    ap.add_argument("--confirm-real-run", action="store_true",
                    help="with --snapshot only: a person confirms that this real snapshot's fetch configuration was checked; sets "
                         "config_verified in its manifest (the snapshot hash changes, so tag, draw and count again afterwards)")
    a = ap.parse_args(argv)
    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, "reconfigure"):
            stream.reconfigure(errors="backslashreplace")
    lex = lex_mod.load()
    failures, blocked, summary = [], [], None
    try:
        if a.confirm_real_run:
            return confirm_real_run(a.snapshot)
        manual = scope_mod.read_id_reasons(a.exclude) if a.exclude else {}
        explained = scope_mod.read_id_reasons(a.explained) if a.explained else {}
        if bool(a.frozen_sheet) != bool(a.frozen_key):
            raise UsageError("--frozen-sheet and --frozen-key go together")
        if not a.snapshot:
            if not a.negative_controls or any((a.tags, a.counts, a.write_counts, a.model_tags, a.verify_model_tags, a.frozen_sheet)):
                raise UsageError("give --snapshot and --tags (or only --negative-controls)")
            results = run_negative_controls(lex, control_snapshot=load_input_snapshot(a.control_snapshot, "control snapshot")
                                            if a.control_snapshot else None)
            for name, ok, detail in results:
                print(f"{'PASS' if ok else 'FAIL'} control: {name}: {detail}")
            print("NOT THE GATE: controls only")
            if not all(ok for _, ok, _ in results):
                return 1
            return 0 if a.controls_only else 4
        if not a.tags:
            raise UsageError("--tags is required with --snapshot")
        s = snap.load(a.snapshot)
        tags = validate_tags(read_json(a.tags, "tags file"), lex)
        mt = a.model_tags or a.verify_model_tags
        proposed = read_json(mt, "proposed model tags") if mt else []
    except (UsageError, ValueError, OSError) as exc:
        print(f"ERROR: {snap.clean(exc, 300)}")
        return 2
    except snap.SnapshotError as exc:
        print(f"FAIL integrity: {snap.clean(exc, 300)}")
        return 1

    entries = entry_texts(s, manual)
    rule_tags = {e["entry_id"]: e for e in tags.get("entries", [])}
    accepted, rejected = verify_model_tags(proposed, entries, rule_tags, lex)
    if a.verify_model_tags:
        for t in accepted:
            print(f"ACCEPT {snap.clean(t['entry_id'], 40)} {snap.clean(t['class'], 40)}: found in {t['field']} at {t['start']}")
        for r in rejected:
            t = r["tag"] if isinstance(r["tag"], dict) else {}
            print(f"REJECT {snap.clean(t.get('entry_id', '?'), 40)} {snap.clean(t.get('class', '?'), 40)}: {snap.clean(r['reason'], 200)}")
        print(f"{len(accepted)} accepted, {len(rejected)} rejected")
        return 0

    # 1. integrity
    if tags.get("snapshot_sha256") != s.digest:
        failures.append("the tags file was made from a different snapshot")
    if tags.get("lexicon_sha256") != lex.sha256:
        failures.append(f"the tags file was made with a different lexicon ({tags.get('lexicon_version')}); re-run lexicon.py tag")
    retag = lex_mod.tag_snapshot(s, lex, manual)
    if retag["entries"] != tags.get("entries"):
        n_diff = sum(1 for x, y in zip(retag["entries"], tags.get("entries", [])) if x != y) + abs(len(retag["entries"]) - len(tags.get("entries", [])))
        failures.append(f"re-running the lexicon does not reproduce the tags file ({n_diff} entr(ies) differ)")
    counts = compute_counts(s, retag, accepted, lex, manual)
    if a.write_counts:
        if failures:
            for f in failures:
                print(f"FAIL integrity: {f}")
            return 1
        a.write_counts.write_text(json.dumps(counts, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
        print(f"wrote {a.write_counts.name}: {counts['studies']} studies, {counts['entries']} entries")
        return 0
    if not a.counts:
        failures.append("no committed counts file was given (--counts); nothing to compare the recomputed counts with")
    else:
        committed = read_json(a.counts, "counts file")
        drift = diff_paths(counts, committed)
        if drift:
            failures.append(f"{len(drift)} count(s) differ from the committed counts file")
            for d in drift[:40]:
                print(f"  drift {d}")
    print(f"snapshot {s.digest[:12]}  dataTimestamp {s.manifest.get('data_timestamp')}  lexicon {lex.version}")
    print(terms_text(counts["registry_terms"]), end="")
    print(f"studies in scope: {counts['studies']}   outcome entries: {counts['entries']}")
    for f in failures:
        print(f"FAIL integrity: {f}")
    if not failures:
        print("PASS integrity: snapshot hashes, tags and every committed count reproduce")

    # 2. model tags
    if proposed:
        rate = len(rejected) / len(proposed)
        print(f"model tags: {len(proposed)} proposed, {len(accepted)} accepted, {len(rejected)} rejected (quote-rejection rate {rate:.1%})")
        for r in rejected[:20]:
            t = r["tag"] if isinstance(r["tag"], dict) else {}
            print(f"  rejected {snap.clean(t.get('entry_id', '?'), 40)} {snap.clean(t.get('class', '?'), 40)}: {snap.clean(r['reason'], 200)}")
    else:
        print("model tags: none proposed (quote-rejection rate not measured)")
    print("challenger disagreement rate: not measured (no challenger review in this step)")
    sh = counts["shares"]
    print(f"shares of entries: by rule {fmt(sh['by_rule'])}, by model {fmt(sh['by_model'])}, other {fmt(sh['other'])}, "
          f"not stated {fmt(sh['not_stated'])}, unclassified {fmt(sh['unclassified'])}; vague wording {counts['flags']['vague']} entr(ies)")

    # 3. metrics
    m, frozen_problems = None, []
    if a.frozen_sheet:
        try:
            gold, summary = load_frozen(a.frozen_sheet, a.frozen_key, lex)
            frozen_problems = frozen_binding_problems(a.frozen_sheet, a.frozen_key, s, lex, entries)
        except UsageError as exc:
            print(f"ERROR: {snap.clean(exc, 300)}")
            return 2
        for p in frozen_problems:
            print(f"frozen set: {p}")
    if a.frozen_sheet and not frozen_problems:
        model_by = {}
        for t in accepted:
            model_by.setdefault(t["entry_id"], []).append(t)
        predicted = {e["entry_id"]: set(final_classes(e, model_by)) for e in retag["entries"]}
        missing = [e for e in gold if e not in predicted]
        if missing:
            failures.append(f"{len(missing)} frozen-set entr(ies) are not in this snapshot's scope: {missing[:5]}")
        m = metrics({e: g for e, g in gold.items() if e in predicted}, predicted, lex.class_order, THRESHOLDS["min_positives"])
        print(f"frozen set: {summary['rows']} rows, {summary['labelled']} labelled, {summary['unlabelled']} unlabelled")
        for c, v in m["per_class"].items():
            if v["positives"] == 0 and v["predicted"] == 0:
                continue
            parts = []
            for name, k, n, est, ci in (("precision", v["tp"], v["predicted"], v["precision_estimable"], v["precision_ci"]),
                                        ("recall", v["tp"], v["positives"], v["recall_estimable"], v["recall_ci"])):
                if not est:
                    parts.append(f"{name} too few to estimate ({k}/{n} checked)")
                else:
                    parts.append(f"{name} {fmt(k / n)} ({k}/{n}, 95% CI {fmt(ci[0])}-{fmt(ci[1])})")
            print(f"  {c:28} " + "   ".join(parts))
        print(f"  micro F1 {fmt(m['micro_f1'])}   macro F1 {fmt(m['macro_f1'])} over {m['macro_over_classes']} class(es)   "
              f"rows checked {m['rows_checked']}")

    # 4. negative controls
    controls, control = None, None
    if a.negative_controls:
        try:
            canary = read_json(a.canary, "canary tags") if a.canary else None
            planted = read_json(a.planted, "planted ids") if a.planted else None
            if canary is not None and not (isinstance(canary, list) and canary and all(isinstance(t, dict) for t in canary)):
                raise UsageError("the canary file must be a non-empty JSON list of {entry_id, class, quote} objects")
            if planted is not None and not (isinstance(planted, list) and planted and all(snap.valid_nct(i) for i in planted)):
                raise UsageError("the planted-ids file must be a non-empty JSON list of NCT ids (NCT followed by 8 digits)")
            control = snap.load(a.control_snapshot) if a.control_snapshot else None
        except (UsageError, snap.SnapshotError) as exc:
            print(f"ERROR: {snap.clean(exc, 300)}")
            return 2
        on_real = canary is not None or planted is not None
        controls = run_negative_controls(lex, snapshot=s if on_real else None, canary=canary, planted=planted,
                                         control_snapshot=control, manual=manual)
        for name, ok, detail in controls:
            print(f"{'PASS' if ok else 'FAIL'} control: {name}: {detail}")
            if not ok:
                failures.append(f"negative control failed: {name}")

    if failures:
        print(f"RESULT: FAIL ({len(failures)} integrity or control failure(s))")
        return 1
    if a.drift_only:
        print("RESULT: PASS (integrity only; publication stop rules not applied)")
        return 0

    # 5. publication stop rules
    if not allow_synthetic_for_tests:
        for label, snp in (("the snapshot", s), ("the control snapshot", control)):
            marks = synthetic_marks(snp) if snp is not None else []
            if marks:
                blocked.append(f"S0 synthetic data: {label} carries the synthetic marker ({', '.join(marks)}); synthetic data can "
                               "never be published")
    shown = [c for c, n in counts["studies_by_class"].items() if n > 0]
    if frozen_problems:
        blocked.append("S1 the frozen set does not belong to this snapshot and lexicon: " + "; ".join(frozen_problems))
    elif m is None:
        blocked.append("S1 not evaluated: no frozen set was given (publication needs the person-labelled set)")
    else:
        if summary["labelled"] < MIN_FROZEN_ROWS:
            blocked.append(f"S1 only {summary['labelled']} labelled rows; the design needs at least {MIN_FROZEN_ROWS}")
        if summary["unlabelled"]:
            blocked.append(f"S1 {summary['unlabelled']} row(s) of the frozen sheet are unlabelled; every drawn row must be labelled")
        unlabelled_classes = [c for c in shown if m["per_class"][c]["positives"] == 0]
        if unlabelled_classes:
            blocked.append(f"S1 shown class(es) with no labelled positive: {', '.join(unlabelled_classes)} (nothing checks their recall)")
        for c in shown:
            v = m["per_class"][c]
            if not (v["precision_estimable"] and v["recall_estimable"]):
                print(f"  note: {c} is shown but few checked ({v['predicted']} tagged, {v['positives']} labelled in the frozen set); "
                      "the page must say 'few checked'")
            if v["precision_estimable"] and v["precision"] is not None and v["precision"] < THRESHOLDS["precision"]:
                blocked.append(f"S1 {c}: precision {fmt(v['precision'])} is below {THRESHOLDS['precision']}")
            if v["recall_estimable"] and v["recall"] is not None and v["recall"] < THRESHOLDS["recall"]:
                blocked.append(f"S1 {c}: recall {fmt(v['recall'])} is below {THRESHOLDS['recall']}")
    other_share = ((counts["other_entries"] + counts["entries_by_status"]["unclassified"]) / counts["entries"]) if counts["entries"] else 1.0
    if other_share > THRESHOLDS["other_share"]:
        blocked.append(f"S2 other plus unclassified are {other_share:.1%} of entries, over {THRESHOLDS['other_share']:.0%}")
    rd = scope_mod.route_difference(s, manual, explained)
    if rd is None:
        blocked.append("S3 not evaluated: the snapshot has no recall route")
    elif rd["share_unexplained"] is None or rd["share_unexplained"] > THRESHOLDS["route_unexplained"]:
        blocked.append(f"S3 the retrieval routes differ by {fmt(rd['share_unexplained'])} unexplained "
                       f"({len(rd['unexplained'])} studies), over {THRESHOLDS['route_unexplained']:.0%}")
    if controls is None:
        blocked.append("S4 not evaluated: the negative controls were not run (--negative-controls)")
    elif not (a.canary and a.planted and a.control_snapshot):
        blocked.append("S4 controls ran on synthetic fixtures only; the gate needs --canary, --planted and --control-snapshot "
                       "(canary tags and planted ids for this snapshot, and a fetched non-CF control snapshot)")
    else:
        n_control = len(lex_mod.tag_all_entries(control, lex))
        if n_control < MIN_CONTROL_ENTRIES:
            blocked.append(f"S4 the control snapshot has {n_control} outcome entries; at least {MIN_CONTROL_ENTRIES} are needed")
    for b in blocked:
        print(f"STOP {b}")
    if blocked:
        print(f"RESULT: PUBLICATION BLOCKED ({len(blocked)} stop rule(s)); integrity passed")
        return 3
    print("RESULT: PASS (integrity, controls and every publication stop rule)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
