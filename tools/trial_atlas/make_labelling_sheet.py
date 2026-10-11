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

--design census-v1 (proposals/2026-10-11-frozen-census-design.md) draws a census instead: every entry that carries a class with
--small-max (default 30) or fewer entries, every unclassified entry and every not-stated entry, plus a seeded random sample of
--large-per-class (default 40) entries for each larger class (an entry is sampled for its rarest large class). The size follows from the
design, so --n, --unclassified-share, --rare-max and --already are refused with it, and MAX_ROWS does not apply. The classes an entry
carries are its rule tags plus, for an entry the rules left unclassified, its accepted model tag (--model-tags FILE, the proposed-tags
file check_atlas.py reads; each tag is checked by check_atlas.verify_model_tags). The key then records the design, both sizes, the
model-tags hash and, per stratum, its population, the number drawn and the weight (population / drawn). The sheet is as blind as before.
With --consensus-report FILE (the consensus-report.json apply_consensus.py wrote for the --model-tags file) the key's `consensus`
holds the rule's date, the proposed, kept and dropped counts and the report's sha256; otherwise it is null.
--sample-exclude FILE (census-v1 only; a JSON list of entry ids, e.g. the development entries the lexicon was tuned on): those entries
still count for class sizes, stratum assignment and the census strata (an excluded entry in a census stratum is still drawn), but
are never drawn into a sample:<class> stratum and are not part of its population, so its weight is (population without them) /
drawn. The key's `sample_exclusion` is null, or {count: ids of the file in the frame, sha256: check_atlas.sample_exclusion_digest};
check_atlas.py then needs the same file (--frozen-exclude). Ids not in the frame are ignored and counted in --report.

Usage: python make_labelling_sheet.py --snapshot DIR --tags tags.json --sheet SHEET.csv --key KEY.json --seed N [--n 50]
       [--unclassified-share 0.2] [--rare-max N] [--already KEY.json] [--exclude FILE.json] [--report]
       python make_labelling_sheet.py ... --seed N --design census-v1 [--model-tags FILE.json [--consensus-report FILE.json]] [--small-max 30] [--large-per-class 40]
       [--sample-exclude IDS.json]
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
MAX_ROWS = 150                                           # the legacy design only; census-v1's size follows from its strata
DESIGNS = ("legacy", "census-v1")
SMALL_MAX_DEFAULT, LARGE_PER_CLASS_DEFAULT = 30, 40
NOT_STATED = "not_stated"                                # treated through its own census stratum, never as a small or large class
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


def model_tags_digest(raw: bytes) -> str:
    """sha256 of the model-tags file's bytes. check_atlas.py has no hash helper for a model-tags file (it only reads and verifies the
    tags), so the plain file hash is used; any change to the file, a line-ending change included, gives a new hash."""
    return hashlib.sha256(raw).hexdigest()


CONSENSUS_FIELDS = ("rule_declared", "proposed", "kept", "dropped")


def consensus_record(raw: bytes) -> dict:
    """The key's `consensus` record from the bytes of the consensus-report.json that apply_consensus.py writes: only rule_declared,
    proposed, kept and dropped, copied, plus report_sha256 (sha256 of the file's bytes). ValueError when a field is missing or of
    the wrong type."""
    try:
        report = json.loads(raw.decode("utf-8-sig"))
    except (UnicodeDecodeError, ValueError) as exc:
        raise ValueError(f"the consensus report is not JSON: {exc}") from None
    if not isinstance(report, dict):
        raise ValueError("the consensus report must hold a JSON object")
    bad = [f for f in CONSENSUS_FIELDS if not (isinstance(report.get(f), str) and report[f].strip() if f == "rule_declared" else
                                               isinstance(report.get(f), int) and not isinstance(report.get(f), bool) and report[f] >= 0)]
    if bad:
        raise ValueError(f"the consensus report lacks a valid {', '.join(bad)} (is it the consensus-report.json apply_consensus.py writes?)")
    return {**{f: report[f] for f in CONSENSUS_FIELDS}, "report_sha256": hashlib.sha256(raw).hexdigest()}


def draw_census(pool: list[dict], tags_by_id: dict[str, dict], lex, model_by_entry: dict[str, list[dict]], *, seed: int,
                small_max: int, large_per_class: int, sample_exclude=frozenset()) -> tuple[list[dict], dict]:
    """([{entry_id, stratum}] in canonical order, {stratum name: {kind, population, drawn, weight}}) for the census-v1 design.
    Strata are assigned in the order census:small, census:unclassified, census:not_stated, sample:<class>. An empty stratum is left
    out (its weight would be undefined). The class `not_stated` is never small or large: its entries are the census:not_stated stratum.
    sample_exclude (entry ids, e.g. the development entries the rules were tuned on) still counts for the class sizes and the stratum
    assignment, and an excluded entry in a census stratum is drawn as before; it is only taken out of the sample strata, from both
    their draw and their population. With no exclusion the draw is the same as before the option existed."""
    final = {e["entry_id"]: [c for c in check_atlas.final_classes(tags_by_id[e["entry_id"]], model_by_entry) if c != NOT_STATED]
             for e in pool}
    # Counted as check_atlas.frame_class_counts counts them: once per entry (a class an entry carries twice counts once), and a class
    # outside the lexicon's class order counted with .get, so the sheet and the gate agree on every class's size.
    count = {c: 0 for c in lex.class_order}
    for classes in final.values():
        for c in set(classes):
            count[c] = count.get(c, 0) + 1
    order = {c: i for i, c in enumerate(lex.class_order)}
    small = {c for c, k in count.items() if 0 < k <= small_max}
    members: dict[str, list[dict]] = {}
    for e in pool:                                       # the pool is sorted by entry id, so every stratum is too
        classes = final[e["entry_id"]]
        if any(c in small for c in classes):
            name = "census:small"
        elif tags_by_id[e["entry_id"]]["status"] == "not_stated":
            name = "census:not_stated"
        elif not classes:
            name = "census:unclassified"
        else:
            name = "sample:" + min(classes, key=lambda c: (count[c], order.get(c, len(order))))
        members.setdefault(name, []).append(e)
    rng = random.Random(seed)
    order = ["census:small", "census:unclassified", "census:not_stated"] + [f"sample:{c}" for c in lex.class_order]
    order += sorted(n for n in members if n not in order)   # a class outside the class order still keeps its entries
    chosen, strata = [], {}
    for name in order:
        pop = members.get(name, [])
        if not pop:
            continue
        kind = name.split(":")[0]
        if kind == "sample" and sample_exclude:
            pop = [e for e in pop if e["entry_id"] not in sample_exclude]
            if not pop:
                continue
        if kind == "sample" and len(pop) > large_per_class:
            picked = sorted(rng.sample(pop, large_per_class), key=lambda e: e["entry_id"])
        else:
            picked = pop
        chosen += [{"entry_id": e["entry_id"], "stratum": name} for e in picked]
        strata[name] = {"kind": kind, "population": len(pop), "drawn": len(picked), "weight": len(pop) / len(picked)}
    return chosen, strata


def build(snapshot: snap.Snapshot, tags: dict, lex, *, n: int | None = None, seed: int, unclassified_share: float = 0.2,
          rare_max: int | None = None, already: dict | None = None, manual=None, design: str = "legacy",
          small_max: int = SMALL_MAX_DEFAULT, large_per_class: int = LARGE_PER_CLASS_DEFAULT, model_tags: list | None = None,
          model_tags_sha256: str | None = None, consensus: dict | None = None,
          sample_exclude: list | None = None) -> tuple[str, dict, dict]:
    """(sheet CSV text, key, summary). Raises ValueError when the tags were not made from this snapshot with this lexicon.
    design "legacy" needs n; design "census-v1" takes no n and no earlier key, and may take the model tags (the proposed-tags list as
    check_atlas reads it) together with the sha256 of their file, and then the consensus record (consensus_record()) of the rule that
    filtered them. sample_exclude (census-v1 only) is a list of entry ids the sample strata may not draw from (see draw_census); the
    key's sample_exclusion records how many of them sit in the frame and the hash of the list (check_atlas.sample_exclusion_digest)."""
    if design not in DESIGNS:
        raise ValueError(f"unknown design {design!r}")
    census = design == "census-v1"
    if consensus is not None and (not census or model_tags is None):
        raise ValueError("a consensus record goes with census-v1 and its model tags only")
    if census and (n is not None or already):
        raise ValueError("census-v1 sizes itself: it takes neither n nor an earlier key")
    if not census and n is None:
        raise ValueError("the legacy design needs n")
    if not census and model_tags is not None:
        raise ValueError("model tags are used by the census-v1 design only")
    if sample_exclude is not None:
        if not census:
            raise ValueError("a sample exclusion is used by the census-v1 design only")
        check_atlas.sample_exclusion_ids(sample_exclude)
    if census and (small_max < 1 or large_per_class < 1):
        raise ValueError("small_max and large_per_class must each be at least 1")
    if (model_tags is None) != (model_tags_sha256 is None):
        raise ValueError("give the model tags and their file hash together")
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
    by_id = {e["entry_id"]: e for e in entries}
    if census:
        pool = sorted((e for e in entries if not e.get("no_primary_outcome")), key=lambda e: e["entry_id"])
        accepted, rejected = check_atlas.verify_model_tags(model_tags, by_id, tags_by_id, lex) if model_tags is not None else ([], [])
        model_by_entry: dict[str, list[dict]] = {}
        for t in accepted:
            model_by_entry.setdefault(t["entry_id"], []).append(t)
        pool_ids = {e["entry_id"] for e in pool}
        excluded = set(sample_exclude or ()) & pool_ids
        chosen, strata = draw_census(pool, tags_by_id, lex, model_by_entry, seed=seed, small_max=small_max,
                                     large_per_class=large_per_class, sample_exclude=frozenset(excluded))
        if any(c["entry_id"] in excluded and c["stratum"].startswith("sample:") for c in chosen):
            raise ValueError("internal error: an excluded entry was drawn into a sample stratum")
        # census strata are drawn in full, so the excluded entries not among their rows are the ones the sample strata left out
        in_census = sum(1 for c in chosen if c["entry_id"] in excluded and c["stratum"].startswith("census:"))
        if sum(v["population"] for v in strata.values()) != len(pool) - (len(excluded) - in_census):
            raise ValueError("internal error: the strata do not partition the drawable entries")
        if not all(c["stratum"] in strata for c in chosen) or len(chosen) != sum(v["drawn"] for v in strata.values()):
            raise ValueError("internal error: a drawn row is outside the strata")
        summary = {"drawn": len(chosen), "pool": len(pool), "strata": strata, "model_tags_accepted": len(accepted),
                   "model_tags_rejected": len(rejected)}
        sample_exclusion = None
        if sample_exclude is not None:
            sample_exclusion = {"count": len(excluded), "sha256": check_atlas.sample_exclusion_digest(sample_exclude)}
            summary.update({"excluded_in_frame": len(excluded), "excluded_not_in_frame": len(set(sample_exclude) - pool_ids),
                            "excluded_in_census": in_census})
    else:
        if rare_max is None:
            rare_max = max(2, round(0.02 * len(entries)))
        chosen, summary = draw(entries, tags_by_id, study, lex, n=n, seed=seed, unclassified_share=unclassified_share,
                               rare_max=rare_max, skip=skip)
    rng = random.Random(seed + 1)
    rng.shuffle(chosen)                                   # the sheet's order says nothing about why a row was drawn
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
        "seed": seed,
    }
    if census:
        key["design"] = "census-v1"
        key.update({"small_max": small_max, "large_per_class": large_per_class})
    else:
        key.update({"n": n, "unclassified_share": unclassified_share, "rare_max": rare_max})
    key.update({"snapshot_sha256": snapshot.digest, "lexicon_sha256": lex.sha256, "tags_sha256": tags_digest(tags)})
    if census:
        key["model_tags_sha256"] = model_tags_sha256
    key.update({
        "lexicon_version_used_for_the_draw": lex.version,
        "extends": (already or {}).get("sheet_wording_sha256"),
        "sheet_wording_sha256": check_atlas.sheet_wording_sha(rows),
        "notes": "A cell starting with = + - @ or a tab was given a leading apostrophe; nothing else in the wording was changed.",
    })
    if census:
        key["consensus"] = dict(consensus) if consensus is not None else None
        key["sample_exclusion"] = sample_exclusion
        key["strata"] = summary["strata"]
    key["rows"] = [{"sample_id": r["sample_id"], "entry_id": c["entry_id"], "stratum": c["stratum"]} for r, c in zip(rows, chosen)]
    return buf.getvalue(), key, summary


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Draw the frozen set and write a blind labelling sheet and a sealed key.")
    ap.add_argument("--snapshot", type=Path, required=True)
    ap.add_argument("--tags", type=Path, required=True)
    ap.add_argument("--sheet", type=Path, required=True)
    ap.add_argument("--key", type=Path, required=True)
    ap.add_argument("--seed", type=int, required=True, help="required: choose it when the draw is made, and write it down")
    ap.add_argument("--design", choices=DESIGNS, default="legacy", help="legacy: stratified sample of --n rows; census-v1: a census of "
                    "the small classes, the unclassified and not-stated entries, plus a sample per large class (size follows from the data)")
    ap.add_argument("--n", type=int, help=f"legacy design only (default 50, at most {MAX_ROWS})")
    ap.add_argument("--unclassified-share", type=float, help="legacy design only (default 0.2)")
    ap.add_argument("--rare-max", type=int, help="legacy design only")
    ap.add_argument("--model-tags", type=Path, help="census-v1 only: the proposed model-tags file check_atlas.py reads with --model-tags")
    ap.add_argument("--consensus-report", type=Path, help="census-v1 with --model-tags only: the consensus-report.json apply_consensus.py "
                    "wrote for those tags; the key records its rule date and counts")
    ap.add_argument("--small-max", type=int, help=f"census-v1 only: a class with this many entries or fewer is censused (default {SMALL_MAX_DEFAULT})")
    ap.add_argument("--large-per-class", type=int, help=f"census-v1 only: rows sampled per larger class (default {LARGE_PER_CLASS_DEFAULT})")
    ap.add_argument("--sample-exclude", type=Path, help="census-v1 only: a JSON list of entry ids (the development entries) the "
                    "sample strata may not draw from; they still count for class sizes and census strata")
    ap.add_argument("--already", type=Path)
    ap.add_argument("--exclude", type=Path)
    ap.add_argument("--report", action="store_true", help="print how many rows each oversample took (a rule result; not for the labeller)")
    a = ap.parse_args(argv)
    census = a.design == "census-v1"
    wrong = ([f for f, v in (("--n", a.n), ("--unclassified-share", a.unclassified_share), ("--rare-max", a.rare_max),
                             ("--already", a.already)) if v is not None] if census else
             [f for f, v in (("--model-tags", a.model_tags), ("--consensus-report", a.consensus_report), ("--small-max", a.small_max),
                             ("--large-per-class", a.large_per_class), ("--sample-exclude", a.sample_exclude)) if v is not None])
    if wrong:
        print(f"error: {', '.join(wrong)} cannot be used with --design {a.design}", file=sys.stderr)
        return 2
    if a.consensus_report and not a.model_tags:
        print("error: --consensus-report goes with --model-tags (the tags the consensus rule kept)", file=sys.stderr)
        return 2
    if census:
        a.small_max = SMALL_MAX_DEFAULT if a.small_max is None else a.small_max
        a.large_per_class = LARGE_PER_CLASS_DEFAULT if a.large_per_class is None else a.large_per_class
        if a.small_max < 1 or a.large_per_class < 1:
            print("error: --small-max and --large-per-class must each be at least 1", file=sys.stderr)
            return 2
    else:
        a.n = 50 if a.n is None else a.n
        a.unclassified_share = 0.2 if a.unclassified_share is None else a.unclassified_share
        if not 1 <= a.n <= MAX_ROWS or not 0 <= a.unclassified_share <= 1:
            print(f"error: --n must be 1 to {MAX_ROWS} and --unclassified-share 0 to 1", file=sys.stderr)
            return 2
    if key_beside_sheet(a.key, a.sheet):
        print("error: the sealed key may not be written into the same folder as the sheet, or into a folder inside it; give the "
              "labeller only the sheet's folder", file=sys.stderr)
        return 2
    terms_path = a.sheet.with_name(a.sheet.stem + ".TERMS.txt")
    for p in (a.sheet, a.key, terms_path):
        if p.exists():
            print(f"error: {p.name} exists; a sheet, key or terms file is never overwritten", file=sys.stderr)
            return 2
    try:
        s = snap.load(a.snapshot)
        tags = json.loads(a.tags.read_text(encoding="utf-8"))
        already = json.loads(a.already.read_text(encoding="utf-8")) if a.already else None
        manual = scope_mod.read_id_reasons(a.exclude) if a.exclude else {}
        if already and len(already.get("rows", [])) + a.n > MAX_ROWS:
            raise ValueError(f"the earlier key holds {len(already['rows'])} rows; {a.n} more would pass {MAX_ROWS}")
        lex = lex_mod.load()
        if census:
            raw = a.model_tags.read_bytes() if a.model_tags else None
            proposed = json.loads(raw.decode("utf-8-sig")) if raw is not None else None
            consensus = consensus_record(a.consensus_report.read_bytes()) if a.consensus_report else None
            sample_exclude = (check_atlas.sample_exclusion_ids(json.loads(a.sample_exclude.read_text(encoding="utf-8-sig")))
                              if a.sample_exclude else None)
            sheet, key, summary = build(s, tags, lex, seed=a.seed, design="census-v1", small_max=a.small_max,
                                        large_per_class=a.large_per_class, model_tags=proposed, manual=manual,
                                        model_tags_sha256=model_tags_digest(raw) if raw is not None else None, consensus=consensus,
                                        sample_exclude=sample_exclude)
        else:
            sheet, key, summary = build(s, tags, lex, n=a.n, seed=a.seed, unclassified_share=a.unclassified_share,
                                        rare_max=a.rare_max, already=already, manual=manual)
    except (snap.SnapshotError, OSError, ValueError) as exc:
        print(f"ERROR: {snap.clean(exc, 300)}")
        return 1
    a.sheet.write_text(sheet, encoding="utf-8-sig", newline="")
    a.key.write_text(json.dumps(key, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
    # The sheet holds registry text, so the registry's terms go beside it (a CSV has no room for a header comment).
    block = check_atlas.registry_terms(s, lex, scope_mod.apply_scope(s.records(scope_mod.MAIN_ROUTE), manual)["counts"])
    terms_path.write_text(f"Terms for {a.sheet.name}: it holds registry text from {check_atlas.REGISTRY_SOURCE}.\n"
                          + check_atlas.terms_text(block), encoding="utf-8")
    print(f"drew {summary['drawn']} rows (seed {a.seed})")
    if not census and summary["drawn"] < a.n:
        print(f"WARNING: drew {summary['drawn']} of the {a.n} rows asked for: only {summary['pool']} entries could be drawn "
              "(in-scope entries with a primary outcome, less the rows of any earlier key)", file=sys.stderr)
    if a.report and census:                               # counts only: never an entry's text
        print(f"of {summary['pool']} eligible entries in {len(summary['strata'])} strata (model tags: {summary['model_tags_accepted']} "
              f"accepted, {summary['model_tags_rejected']} rejected)")
        if "excluded_in_frame" in summary:
            print(f"sample exclusion: {summary['excluded_in_frame']} ids in the frame ({summary['excluded_in_census']} of them in census "
                  f"strata, still drawn), {summary['excluded_not_in_frame']} ids not in the frame (ignored)")
        for name, v in summary["strata"].items():
            print(f"  {name}: population {v['population']}, drawn {v['drawn']}, weight {v['weight']:.4g}")
    elif a.report:
        print(f"of {summary['pool']} eligible entries: {summary['unclassified']} unclassified oversample, {summary['rare_class']} from "
              f"{summary['rare_classes']} rare class(es), {summary['stratified']} across {summary['strata']} strata")
    print(f"sheet: {a.sheet.name} (blind: wording only)   terms: {terms_path.name}   key: {a.key.name} (sealed; keep it away from "
          "the labeller)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
