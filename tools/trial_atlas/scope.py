#!/usr/bin/env python3
"""Which studies of a snapshot the atlas counts, which it leaves out and why (design, decision 3).

Inclusion: an interventional study (designModule.studyType INTERVENTIONAL) whose condition list names cystic fibrosis, any status and
any date. Each exclusion rule has an id, a reason and a printed count:

  X1  observational study                      X2  expanded access
  X3  study type missing or another type       X4  the condition list does not name CF
  X5  the condition list names CF only in a negation (for example "Non-Cystic Fibrosis Bronchiectasis")
  X6  listed by hand in an exclusions file (--exclude), with the reason written there

Flags on each included study: "CF only" (every condition names CF) or "CF among others"; start "actual", "planned" (the registry's
ESTIMATED start) or "unknown" (no date, or a date with no type); start year and decade; first-posted year.

Only the condition list decides; keywords and titles do not (a keyword is how the recall route finds studies, not proof of scope).
The second retrieval route (a term search) is compared with the main route after the same rules, and the difference is printed.

Usage: python scope.py SNAPSHOT_DIR [--exclude FILE.json] [--explained FILE.json] [--json]
  --exclude    a JSON object {"NCT...": "reason"}: studies a person excluded by hand (X6)
  --explained  a JSON object {"NCT...": "reason"}: route differences a person explained after a hand check
Exit 0 always after a successful read (it reports; check_atlas.py applies the stop rules); 1 the snapshot could not be read.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import snapshot as snap  # noqa: E402

MAIN_ROUTE = "condition"
RECALL_ROUTE = "term"

# A condition names CF when one of these matches and the match is not inside a negation.
CF_TERMS = re.compile(r"\bcystic[\s-]+fibrosis\b|\bmucoviscidosis\b|^\s*CF\s*$", re.I)
CF_NEGATION = re.compile(r"\bnon[\s-]*(?:cystic[\s-]+fibrosis|CF)\b|\bwithout\s+(?:cystic[\s-]+fibrosis|CF)\b|"
                         r"\bother\s+than\s+(?:cystic[\s-]+fibrosis|CF)\b", re.I)

RULES = {
    "X1": "observational study (out of scope: the atlas counts interventional studies)",
    "X2": "expanded access (out of scope: not a trial with a registered primary outcome to compare)",
    "X3": "study type missing or not interventional",
    "X4": "the condition list does not name cystic fibrosis",
    "X5": "the condition list names cystic fibrosis only in a negation (for example non-CF bronchiectasis)",
    "X6": "excluded by hand (reason in the exclusions file)",
}


def condition_names_cf(condition: str) -> bool | None:
    """True: names CF. None: names CF only inside a negation. False: does not name CF."""
    if CF_TERMS.search(CF_NEGATION.sub(" ", condition)):
        return True
    return None if CF_NEGATION.search(condition) else False


def start_info(record: dict) -> dict:
    date = record.get("start_date") or ""
    m = re.match(r"^(\d{4})", date)
    year = int(m.group(1)) if m else None
    stype = (record.get("start_type") or "").upper()
    kind = "unknown" if year is None or not stype else ("actual" if stype == "ACTUAL" else "planned" if stype == "ESTIMATED" else "unknown")
    fp = re.match(r"^(\d{4})", record.get("first_posted") or "")
    return {"start_year": year, "start_decade": f"{year // 10 * 10}s" if year else "unknown", "start_kind": kind,
            "first_posted_year": int(fp.group(1)) if fp else None}


def classify(record: dict, manual: dict[str, str] | None = None) -> tuple[str | None, str]:
    """(rule id, reason) when excluded, (None, flag) when included; flag is 'CF only' or 'CF among others'."""
    manual = manual or {}
    if record["nct_id"] in manual:
        return "X6", f"{RULES['X6']}: {manual[record['nct_id']]}"
    st = (record.get("study_type") or "").upper()
    if st == "OBSERVATIONAL":
        return "X1", RULES["X1"]
    if st == "EXPANDED_ACCESS":
        return "X2", RULES["X2"]
    if st != "INTERVENTIONAL":
        return "X3", f"{RULES['X3']} ({st or 'missing'})"
    verdicts = [condition_names_cf(c) for c in record.get("conditions", [])]
    if True not in verdicts:
        if None in verdicts:
            return "X5", RULES["X5"]
        return "X4", RULES["X4"]
    return None, "CF only" if all(v is True for v in verdicts) else "CF among others"


def apply_scope(records: list[dict], manual: dict[str, str] | None = None) -> dict:
    included, excluded = [], []
    for r in records:
        rule, why = classify(r, manual)
        if rule:
            excluded.append({"nct_id": r["nct_id"], "rule": rule, "reason": why})
        else:
            included.append({**r, "cf_flag": why, **start_info(r)})
    counts = {rid: sum(1 for e in excluded if e["rule"] == rid) for rid in RULES}
    return {
        "included": included,
        "excluded": excluded,
        "counts": {
            "studies_in_route": len(records),
            "included": len(included),
            "excluded_by_rule": counts,
            "cf_only": sum(1 for r in included if r["cf_flag"] == "CF only"),
            "cf_among_others": sum(1 for r in included if r["cf_flag"] == "CF among others"),
            "start_actual": sum(1 for r in included if r["start_kind"] == "actual"),
            "start_planned": sum(1 for r in included if r["start_kind"] == "planned"),
            "start_unknown": sum(1 for r in included if r["start_kind"] == "unknown"),
        },
    }


def route_difference(snapshot: snap.Snapshot, manual=None, explained=None) -> dict | None:
    """The in-scope studies of the two routes compared after the same rules. None when the snapshot has no recall route.
    share_unexplained = (differing ids not in the explained file) / (in-scope studies of the main route)."""
    if RECALL_ROUTE not in snapshot.routes():
        return None
    explained = explained or {}
    main_ids = {r["nct_id"] for r in apply_scope(snapshot.records(MAIN_ROUTE), manual)["included"]}
    recall_ids = {r["nct_id"] for r in apply_scope(snapshot.records(RECALL_ROUTE), manual)["included"]}
    only_main, only_recall = sorted(main_ids - recall_ids), sorted(recall_ids - main_ids)
    unexplained = [i for i in only_main + only_recall if i not in explained]
    return {"main_in_scope": len(main_ids), "recall_in_scope": len(recall_ids), "only_main": only_main, "only_recall": only_recall,
            "unexplained": unexplained, "share_unexplained": (len(unexplained) / len(main_ids)) if main_ids else None}


def read_id_reasons(path) -> dict[str, str]:
    data = json.loads(Path(path).read_text(encoding="utf-8-sig"))
    if not isinstance(data, dict) or not all(isinstance(k, str) and snap.NCT_ID.match(k) and isinstance(v, str) and v.strip()
                                             for k, v in data.items()):
        raise ValueError(f"{path} must be a JSON object of NCT id -> a non-empty reason")
    return data


def report(snapshot: snap.Snapshot, manual=None, explained=None) -> dict:
    result = apply_scope(snapshot.records(MAIN_ROUTE), manual)
    return {"counts": result["counts"], "excluded": result["excluded"], "routes": route_difference(snapshot, manual, explained)}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Apply the atlas scope rules to a snapshot and print every count with its reason.")
    ap.add_argument("snapshot", type=Path)
    ap.add_argument("--exclude", type=Path)
    ap.add_argument("--explained", type=Path)
    ap.add_argument("--json", action="store_true")
    a = ap.parse_args(argv)
    try:
        s = snap.load(a.snapshot)
        manual = read_id_reasons(a.exclude) if a.exclude else {}
        explained = read_id_reasons(a.explained) if a.explained else {}
    except (snap.SnapshotError, OSError, ValueError) as exc:
        print(f"ERROR: {exc}")
        return 1
    rep = report(s, manual, explained)
    if a.json:
        print(json.dumps(rep, indent=2))
        return 0
    c = rep["counts"]
    print(f"route '{MAIN_ROUTE}': {c['studies_in_route']} studies in the snapshot, {c['included']} included")
    for rid, n in c["excluded_by_rule"].items():
        print(f"  excluded {rid}: {n}  {RULES[rid]}")
    print(f"  CF only: {c['cf_only']}   CF among others: {c['cf_among_others']}")
    print(f"  start actual: {c['start_actual']}   planned: {c['start_planned']}   unknown: {c['start_unknown']}")
    for e in rep["excluded"]:
        print(f"  {e['nct_id']}  {e['rule']}  {e['reason']}")
    rd = rep["routes"]
    if rd is None:
        print(f"no '{RECALL_ROUTE}' route in this snapshot: the recall check was not run")
    else:
        share = "n/a" if rd["share_unexplained"] is None else f"{rd['share_unexplained']:.1%}"
        print(f"routes: main {rd['main_in_scope']} in scope, recall {rd['recall_in_scope']} in scope; only in main "
              f"{len(rd['only_main'])}, only in recall {len(rd['only_recall'])}; unexplained {len(rd['unexplained'])} ({share} of main)")
        for i in rd["unexplained"]:
            print(f"  unexplained difference: {i}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
