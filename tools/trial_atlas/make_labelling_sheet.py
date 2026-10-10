#!/usr/bin/env python3
"""Draw the frozen set for the trial atlas and write a blind labelling sheet plus a sealed key (design, decision 5).

The draw is stratified by start decade and phase, with two oversamples taken first: entries the lexicon left unclassified (a share of
the sample, default 20%) and up to two entries of each rare class (fewer entries than --rare-max, default 2% of entries but at least
2), rarest first. The two oversamples together take at most half the sample. The rest is allocated to the strata in proportion to
their size (largest remainder). These proportions are a first draft for the maintainer to confirm before the draw. A fixed seed makes the draw reproducible: the same snapshot, tags and
seed give the same sheet and key, byte for byte.

The sheet shows only the registry wording (measure, description, time frame) and empty columns for the person: label_classes (class
ids from lexicon.json, separated by ';'), label_composite (yes/no) and label_notes. It holds NO rule result, NO model result, no NCT
id and no stratum, and its rows are shuffled, so the person labels blind. A cell that starts with = + - @ or a tab gets a leading
apostrophe so a spreadsheet does not run it as a formula (the only change to the wording; it is written into the key's notes).

The key (keep it closed until labelling is finished) maps each sample id to its entry id and stratum, and holds the seed, the snapshot
hash and the hash of the sheet's wording, so check_atlas.py can tell a filled-in sheet from an edited one. Rows drawn here stay out of
rule development. --already KEY.json skips the rows of an earlier key, to extend the set (up to 150 in all).

Usage: python make_labelling_sheet.py --snapshot DIR --tags tags.json --sheet SHEET.csv --key KEY.json [--n 50] [--seed 20261010]
       [--unclassified-share 0.2] [--rare-max N] [--already KEY.json] [--exclude FILE.json]
Exit 0 written; 1 the inputs could not be read or do not match; 2 usage error.
"""
from __future__ import annotations

import argparse
import csv
import io
import json
import random
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import check_atlas  # noqa: E402
import lexicon as lex_mod  # noqa: E402
import scope as scope_mod  # noqa: E402
import snapshot as snap  # noqa: E402

COLUMNS = ["sample_id", "measure", "description", "time_frame", "label_classes", "label_composite", "label_notes"]
MAX_ROWS = 150
DEFAULT_SEED = 20261010
FORMULA_START = ("=", "+", "-", "@", "\t", "\r")


def safe_cell(text: str) -> str:
    return "'" + text if text.startswith(FORMULA_START) else text


def phase_label(phases) -> str:
    return "/".join(phases) if phases else "none"


def draw(entries: list[dict], tags_by_id: dict[str, dict], study: dict[str, dict], lex, *, n: int, seed: int,
         unclassified_share: float, rare_max: int, skip: set[str]) -> tuple[list[dict], dict]:
    """[{entry_id, stratum, reason}] in drawn order, and the allocation summary (counts only)."""
    rng = random.Random(seed)
    pool = sorted((e for e in entries if e["entry_id"] not in skip), key=lambda e: e["entry_id"])
    chosen, taken = [], set()

    def stratum(e):
        r = study[e["nct_id"]]
        return f"{r['start_decade']} | {phase_label(r['phases'])}"

    def take(cands, k, reason):
        cands = [e for e in cands if e["entry_id"] not in taken]
        for e in rng.sample(cands, min(k, len(cands), n - len(chosen))):
            chosen.append({"entry_id": e["entry_id"], "stratum": stratum(e), "reason": reason})
            taken.add(e["entry_id"])

    take([e for e in pool if tags_by_id[e["entry_id"]]["status"] == "unclassified"], round(n * unclassified_share), "unclassified")
    class_count = {c: 0 for c in lex.class_order}
    for e in pool:
        for t in tags_by_id[e["entry_id"]]["tags"]:
            class_count[t["class"]] += 1
    # Rarest first; the two oversamples together never take more than half the sample, so the strata keep the other half.
    rare = sorted((c for c in lex.class_order if 0 < class_count[c] < rare_max), key=lambda c: (class_count[c], lex.class_order.index(c)))
    for c in rare:
        budget = n // 2 - len(chosen)
        if budget <= 0:
            break
        take([e for e in pool if any(t["class"] == c for t in tags_by_id[e["entry_id"]]["tags"])], min(2, budget), "rare class")
    left = n - len(chosen)
    rest = [e for e in pool if e["entry_id"] not in taken]
    strata: dict[str, list[dict]] = {}
    for e in rest:
        strata.setdefault(stratum(e), []).append(e)
    if left > 0 and rest:
        total = len(rest)
        exact = {s: left * len(v) / total for s, v in strata.items()}
        alloc = {s: int(x) for s, x in exact.items()}
        for s in sorted(exact, key=lambda s: (-(exact[s] - alloc[s]), s))[: left - sum(alloc.values())]:
            alloc[s] += 1
        for s in sorted(strata):
            take(strata[s], alloc[s], "stratified")
    summary = {"drawn": len(chosen), "unclassified": sum(1 for c in chosen if c["reason"] == "unclassified"),
               "rare_class": sum(1 for c in chosen if c["reason"] == "rare class"), "rare_classes": len(rare),
               "stratified": sum(1 for c in chosen if c["reason"] == "stratified"), "strata": len(strata), "pool": len(pool)}
    return chosen, summary


def build(snapshot: snap.Snapshot, tags: dict, lex, *, n: int, seed: int, unclassified_share: float = 0.2, rare_max: int | None = None,
          already: dict | None = None, manual=None) -> tuple[str, dict, dict]:
    """(sheet CSV text, key, summary)."""
    entries, result = lex_mod.in_scope_entries(snapshot, manual)
    tags_by_id = {e["entry_id"]: e for e in tags["entries"]}
    if set(tags_by_id) != {e["entry_id"] for e in entries}:
        raise ValueError("the tags file does not cover exactly this snapshot's in-scope entries; re-run lexicon.py tag")
    study = {r["nct_id"]: r for r in result["included"]}
    skip = {r["entry_id"] for r in (already or {}).get("rows", [])}
    start_no = len((already or {}).get("rows", [])) + 1
    if rare_max is None:
        rare_max = max(2, round(0.02 * len(entries)))
    chosen, summary = draw(entries, tags_by_id, study, lex, n=n, seed=seed, unclassified_share=unclassified_share,
                           rare_max=rare_max, skip=skip)
    rng = random.Random(seed + 1)
    rng.shuffle(chosen)                                   # the sheet's order says nothing about why a row was drawn
    by_id = {e["entry_id"]: e for e in entries}
    rows = []
    for i, c in enumerate(chosen, start_no):
        e = by_id[c["entry_id"]]
        rows.append({"sample_id": f"S{i:03d}", "measure": safe_cell(e["measure"]), "description": safe_cell(e["description"]),
                     "time_frame": safe_cell(e["time_frame"]), "label_classes": "", "label_composite": "", "label_notes": ""})
    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=COLUMNS, lineterminator="\n")
    w.writeheader()
    w.writerows(rows)
    key = {
        "kind": "trial-atlas frozen-set key",
        "sealed": "Keep closed until labelling is finished. These rows stay out of rule development.",
        "seed": seed, "n": n, "unclassified_share": unclassified_share, "rare_max": rare_max,
        "snapshot_sha256": snapshot.digest, "lexicon_version_used_for_the_draw": lex.version,
        "extends": (already or {}).get("sheet_wording_sha256"),
        "sheet_wording_sha256": check_atlas.sheet_wording_sha(rows),
        "notes": "A cell starting with = + - @ or a tab was given a leading apostrophe; nothing else in the wording was changed.",
        "rows": [{"sample_id": r["sample_id"], "entry_id": c["entry_id"], "stratum": c["stratum"]} for r, c in zip(rows, chosen)],
    }
    return buf.getvalue(), key, summary


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Draw the frozen set and write a blind labelling sheet and a sealed key.")
    ap.add_argument("--snapshot", type=Path, required=True)
    ap.add_argument("--tags", type=Path, required=True)
    ap.add_argument("--sheet", type=Path, required=True)
    ap.add_argument("--key", type=Path, required=True)
    ap.add_argument("--n", type=int, default=50)
    ap.add_argument("--seed", type=int, default=DEFAULT_SEED)
    ap.add_argument("--unclassified-share", type=float, default=0.2)
    ap.add_argument("--rare-max", type=int)
    ap.add_argument("--already", type=Path)
    ap.add_argument("--exclude", type=Path)
    a = ap.parse_args(argv)
    if not 1 <= a.n <= MAX_ROWS or not 0 <= a.unclassified_share <= 1:
        print(f"error: --n must be 1 to {MAX_ROWS} and --unclassified-share 0 to 1", file=sys.stderr)
        return 2
    for p in (a.sheet, a.key):
        if p.exists():
            print(f"error: {p.name} exists; a sheet or key is never overwritten", file=sys.stderr)
            return 2
    try:
        s = snap.load(a.snapshot)
        tags = json.loads(a.tags.read_text(encoding="utf-8"))
        already = json.loads(a.already.read_text(encoding="utf-8")) if a.already else None
        manual = scope_mod.read_id_reasons(a.exclude) if a.exclude else {}
        if already and len(already.get("rows", [])) + a.n > MAX_ROWS:
            raise ValueError(f"the earlier key holds {len(already['rows'])} rows; {a.n} more would pass {MAX_ROWS}")
        sheet, key, summary = build(s, tags, lex_mod.load(), n=a.n, seed=a.seed, unclassified_share=a.unclassified_share,
                                    rare_max=a.rare_max, already=already, manual=manual)
    except (snap.SnapshotError, OSError, ValueError) as exc:
        print(f"ERROR: {exc}")
        return 1
    a.sheet.write_text(sheet, encoding="utf-8", newline="")
    a.key.write_text(json.dumps(key, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"drew {summary['drawn']} of {summary['pool']} entries (seed {a.seed}): {summary['unclassified']} unclassified oversample, "
          f"{summary['rare_class']} from {summary['rare_classes']} rare class(es), {summary['stratified']} across {summary['strata']} strata")
    print(f"sheet: {a.sheet.name} (blind: wording only)   key: {a.key.name} (sealed)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
