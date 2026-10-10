#!/usr/bin/env python3
"""Draw the frozen set for the trial atlas and write a blind labelling sheet plus a sealed key (design, decision 5).

Run by the person who prepares the set, NOT by the person who labels it: with --report it prints how many rows each oversample took,
which is a rule result.

The draw is stratified by start decade and phase, with two oversamples taken first: entries the lexicon left unclassified (a share of
the sample, default 20%) and up to two entries of each rare class (fewer entries than --rare-max, default 2% of entries but at least
2), rarest first. The two oversamples together take at most half the sample. The rest goes first one row to every non-empty stratum
(the largest strata first when there are more strata than rows), then to the strata in proportion to what is left in them (largest
remainder). Entries of studies that entered no primary outcome are never drawn: there is no wording to label. These proportions are a
first draft for the maintainer to confirm before the draw. The seed must be given (--seed): the same snapshot, tags and seed give the
same sheet and key, byte for byte.

The sheet shows only the registry wording (measure, description, time frame) and empty columns for the person: label_classes (class
ids from lexicon.json, separated by ';'), label_composite (yes/no) and label_notes. It holds NO rule result, NO model result, no NCT
id and no stratum, and its rows are shuffled, so the person labels blind. It is written as UTF-8 with a signature so a spreadsheet
opens it correctly; save it back as "CSV UTF-8". A cell that starts with = + - @ or a tab gets a leading apostrophe so a spreadsheet
does not run it as a formula (the only change to the wording; it is written into the key's notes).

The key (keep it closed until labelling is finished; it may not be written into the sheet's folder) maps each sample id to its entry
id and stratum, and holds the seed, the hashes of the snapshot, the lexicon and the tags, and the hash of the sheet's wording, so
check_atlas.py can tell a filled-in sheet from an edited one. Tags made from another snapshot or another lexicon are refused. Rows
drawn here stay out of rule development. --already KEY.json skips the rows of an earlier key, to extend the set (up to 150 in all).

Usage: python make_labelling_sheet.py --snapshot DIR --tags tags.json --sheet SHEET.csv --key KEY.json --seed N [--n 50]
       [--unclassified-share 0.2] [--rare-max N] [--already KEY.json] [--exclude FILE.json] [--report]
Exit 0 written; 1 the inputs could not be read or do not match; 2 usage error.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import os
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
FORMULA_START = ("=", "+", "-", "@", "\t", "\r")


def safe_cell(text: str) -> str:
    return "'" + text if text.startswith(FORMULA_START) else text


def phase_label(phases) -> str:
    return "/".join(phases) if phases else "none"


def key_beside_sheet(key: Path, sheet: Path) -> bool:
    """True when the key's folder is the sheet's folder or inside it (paths resolved; letter case ignored where the system ignores it)."""
    key_dir = Path(os.path.normcase(str(Path(key).resolve().parent)))
    sheet_dir = Path(os.path.normcase(str(Path(sheet).resolve().parent)))
    return key_dir == sheet_dir or key_dir.is_relative_to(sheet_dir)


def tags_digest(tags: dict) -> str:
    """sha256 of the tags in a canonical form (sorted keys, no spaces), so a line-ending change does not change it."""
    return hashlib.sha256(json.dumps(tags, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")).hexdigest()


def allocate(sizes: dict[str, int], left: int) -> dict[str, int]:
    """Rows per stratum: one for every non-empty stratum first (largest strata first if rows run short), then the rest in proportion
    to what each stratum still holds (largest remainder). Never more rows than a stratum holds."""
    alloc = {s: 0 for s in sizes}
    for s in sorted(sizes, key=lambda s: (-sizes[s], s)):
        if sum(alloc.values()) >= left:
            break
        if sizes[s] > 0:
            alloc[s] = 1
    cap = {s: sizes[s] - alloc[s] for s in sizes}
    remaining = min(left - sum(alloc.values()), sum(cap.values()))
    if remaining > 0:
        total = sum(cap.values())
        exact = {s: remaining * cap[s] / total for s in sizes}
        add = {s: int(x) for s, x in exact.items()}
        for s in sorted(exact, key=lambda s: (-(exact[s] - add[s]), s))[: remaining - sum(add.values())]:
            add[s] += 1
        for s in sizes:
            alloc[s] += add[s]
    return alloc


def draw(entries: list[dict], tags_by_id: dict[str, dict], study: dict[str, dict], lex, *, n: int, seed: int,
         unclassified_share: float, rare_max: int, skip: set[str]) -> tuple[list[dict], dict]:
    """[{entry_id, stratum, reason}] in drawn order, and the allocation summary (counts only)."""
    rng = random.Random(seed)
    pool = sorted((e for e in entries if e["entry_id"] not in skip and not e.get("no_primary_outcome")), key=lambda e: e["entry_id"])
    chosen, taken = [], set()

    def stratum(e):
        r = study[e["nct_id"]]
        return f"{r['start_decade']} | {phase_label(r['phases'])}"

    def take(cands, k, reason):
        cands = [e for e in cands if e["entry_id"] not in taken]
        for e in rng.sample(cands, max(0, min(k, len(cands), n - len(chosen)))):
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
        alloc = allocate({s: len(v) for s, v in strata.items()}, left)
        for s in sorted(strata):
            take(strata[s], alloc[s], "stratified")
    summary = {"drawn": len(chosen), "unclassified": sum(1 for c in chosen if c["reason"] == "unclassified"),
               "rare_class": sum(1 for c in chosen if c["reason"] == "rare class"), "rare_classes": len(rare),
               "stratified": sum(1 for c in chosen if c["reason"] == "stratified"), "strata": len(strata), "pool": len(pool)}
    return chosen, summary


def build(snapshot: snap.Snapshot, tags: dict, lex, *, n: int, seed: int, unclassified_share: float = 0.2, rare_max: int | None = None,
          already: dict | None = None, manual=None) -> tuple[str, dict, dict]:
    """(sheet CSV text, key, summary). Raises ValueError when the tags were not made from this snapshot with this lexicon."""
    if not isinstance(tags, dict) or not isinstance(tags.get("entries"), list):
        raise ValueError("the tags file is not a tags file (an object with an 'entries' list)")
    if tags.get("snapshot_sha256") != snapshot.digest:
        raise ValueError("the tags were made from a different snapshot; re-run lexicon.py tag on this one")
    if tags.get("lexicon_sha256") != lex.sha256:
        raise ValueError("the tags were made with a different lexicon; re-run lexicon.py tag")
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
        "snapshot_sha256": snapshot.digest, "lexicon_sha256": lex.sha256, "tags_sha256": tags_digest(tags),
        "lexicon_version_used_for_the_draw": lex.version,
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
    ap.add_argument("--seed", type=int, required=True, help="required: choose it when the draw is made, and write it down")
    ap.add_argument("--n", type=int, default=50)
    ap.add_argument("--unclassified-share", type=float, default=0.2)
    ap.add_argument("--rare-max", type=int)
    ap.add_argument("--already", type=Path)
    ap.add_argument("--exclude", type=Path)
    ap.add_argument("--report", action="store_true", help="print how many rows each oversample took (a rule result; not for the labeller)")
    a = ap.parse_args(argv)
    if not 1 <= a.n <= MAX_ROWS or not 0 <= a.unclassified_share <= 1:
        print(f"error: --n must be 1 to {MAX_ROWS} and --unclassified-share 0 to 1", file=sys.stderr)
        return 2
    if key_beside_sheet(a.key, a.sheet):
        print("error: the sealed key may not be written into the same folder as the sheet, or into a folder inside it; give the "
              "labeller only the sheet's folder", file=sys.stderr)
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
        print(f"ERROR: {snap.clean(exc, 300)}")
        return 1
    a.sheet.write_text(sheet, encoding="utf-8-sig", newline="")
    a.key.write_text(json.dumps(key, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"drew {summary['drawn']} rows (seed {a.seed})")
    if a.report:
        print(f"of {summary['pool']} eligible entries: {summary['unclassified']} unclassified oversample, {summary['rare_class']} from "
              f"{summary['rare_classes']} rare class(es), {summary['stratified']} across {summary['strata']} strata")
    print(f"sheet: {a.sheet.name} (blind: wording only)   key: {a.key.name} (sealed; keep it away from the labeller)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
