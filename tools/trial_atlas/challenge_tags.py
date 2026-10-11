#!/usr/bin/env python3
"""A BLIND second opinion on the atlas tags from a model of a different family.

A seeded sample is drawn from two pools: entries the proposer's model route accepted (model-tags.json in the proposer's output
folder) and entries the rule lexicon tagged. The challenger sees exactly what the proposer saw: the entry's measure, description and
time frame inside the data boundary, and the class list. It never sees the proposer's class, the rule's class or quote, or which
pool the entry came from (the packet is built by propose_tags.build_packet, which takes only the entry text and the class list). It
answers like the proposer (one class id or none, plus an exact quote), and its quote is verified the same way (check_atlas's
verifier; a rule-tagged entry is not "unclassified", so for the challenger's quote check only it is treated as such). An entry whose
challenger answer never verifies is left out of the agreement rate and counted as unverified.

Agreement: the challenger's verified class is one of the reference classes (the proposer's class, or the rule tag's classes; a rule
entry can carry several). A verified none disagrees. The rate is given overall and per pool, with a Wilson 95% interval from
check_atlas.wilson, or "too few to estimate" under 10 compared entries.

Outputs in --out (outside any git working tree):
  challenge.jsonl   append-only, one row per attempt (as proposals.jsonl, with the challenger's model)
  challenge.json    per entry: pool, reference classes, challenger class, verified, agree; the rates; a disagreement list (entry
                    ids and classes only, for a person or a council to review). No registry text.
  challenge-run.json  what the folder is bound to, and the drawn sample

Usage:
  python challenge_tags.py --snapshot SNAP --tags TAGS --proposals PROPOSER_OUT --out DIR (--model NAME | --profile-file P.json
                           [--thinking-profile-file T.json]) --lab-repo LAB [--invoker PATH] --sample N [--seed N] [--resume]
                           [--think-after-fast] [--temperature T] [--timeout SEC] [--exclude FILE]
The challenger's model must differ from every model the proposer used, also in the part before ':' (qwen3.8:27b and
qwen3.8:27b-64k are one family here). Exit 0 finished; 2 usage error; 130 interrupted (run again with --resume).

Standard library only. No network of its own; the invoker talks to a local Ollama only. Tests use a fake invoker.
"""
from __future__ import annotations

import argparse
import json
import random
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import check_atlas as ca  # noqa: E402
import lexicon as lex_mod  # noqa: E402
import propose_tags as pt  # noqa: E402
import snapshot as snap  # noqa: E402

JSONL, OUT_JSON, RUN = "challenge.jsonl", "challenge.json", "challenge-run.json"
MIN_FOR_INTERVAL = 10


def base_name(model: str) -> str:
    return (model or "").split(":", 1)[0].strip().lower()


def refuse_same_model(challenger_models: list[str], proposer_models: list[str]) -> None:
    for c in challenger_models:
        for p in proposer_models:
            if not p:
                continue
            if c.strip().lower() == p.strip().lower():
                raise pt.UsageError(f"the challenger model {snap.clean(c, 60)} is the proposer's model; a challenger must be a "
                                    "different model family")
            if base_name(c) == base_name(p):
                raise pt.UsageError(f"the challenger model {snap.clean(c, 60)} has the same family name as the proposer's "
                                    f"{snap.clean(p, 60)}; choose a model of a different family")


def draw_sample(model_pool: list[str], rule_pool: list[str], n: int, seed: int) -> list[tuple[str, str]]:
    """[(entry id, pool)]: half (rounded up) from the model pool and the rest from the rule pool; a short pool is topped up from the
    other. The order is shuffled so the pools are mixed."""
    rng = random.Random(seed)
    m, r = sorted(model_pool), sorted(rule_pool)
    want_m = min(len(m), (n + 1) // 2)
    want_r = min(len(r), n - want_m)
    want_m = min(len(m), n - want_r)
    picked = [(e, "model") for e in rng.sample(m, want_m)] + [(e, "rule") for e in rng.sample(r, want_r)]
    rng.shuffle(picked)
    return picked


def rate(agree: int, n: int) -> dict:
    if n == 0:
        return {"agree": 0, "compared": 0, "rate": None, "wilson95": "too few to estimate"}
    ci = ca.wilson(agree, n) if n >= MIN_FOR_INTERVAL else None
    return {"agree": agree, "compared": n, "rate": agree / n,
            "wilson95": [round(ci[0], 4), round(ci[1], 4)] if ci else "too few to estimate"}


def result(rows: list[dict], sample: list[list[str]], refs: dict[str, list[str]], binding: dict) -> dict:
    fv = pt.first_verified(rows)
    per, counts = [], {"all": [0, 0], "model": [0, 0], "rule": [0, 0]}
    for eid, pool in sample:
        r = fv.get(eid)
        item = {"entry_id": eid, "pool": pool, "reference_classes": refs[eid], "challenger_class": r["class_id"] if r else None,
                "verified": bool(r), "agree": None}
        if r:
            item["agree"] = r["class_id"] in refs[eid]
            for k in ("all", pool):
                counts[k][0] += 1 if item["agree"] else 0
                counts[k][1] += 1
        per.append(item)
    return {"kind": "trial-atlas challenger review (private)", **{k: binding[k] for k in
            ("snapshot_sha256", "lexicon_version", "proposer_models", "model", "profile", "think_model", "think_profile",
             "template", "seed", "sample_requested")},
            "sample_drawn": {"model": sum(1 for _, p in sample if p == "model"), "rule": sum(1 for _, p in sample if p == "rule")},
            "unverified": sum(1 for i in per if not i["verified"]),
            "agreement": {"overall": rate(*counts["all"]), "model": rate(*counts["model"]), "rule": rate(*counts["rule"])},
            "note": "Agreement between two models is not correctness; the frozen set, labelled by a person, measures that.",
            "entries": per,
            "disagreements": [{k: i[k] for k in ("entry_id", "pool", "reference_classes", "challenger_class")}
                              for i in per if i["agree"] is False]}


def main(argv=None, *, invoker: pt.Invoker | None = None) -> int:
    ap = argparse.ArgumentParser(description="Blind challenger review of the trial atlas tags (a different model family).")
    ap.add_argument("--snapshot", type=Path)
    ap.add_argument("--tags", type=Path)
    ap.add_argument("--proposals", type=Path, help="the proposer's output folder (run.json and model-tags.json)")
    ap.add_argument("--out", type=Path)
    pt.add_model_options(ap)
    ap.add_argument("--lab-repo", type=Path)
    ap.add_argument("--invoker", type=Path)
    ap.add_argument("--exclude", type=Path)
    ap.add_argument("--sample", type=int)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--resume", action="store_true")
    ap.add_argument("--think-after-fast", action="store_true")
    ap.add_argument("--temperature", type=float, default=0.0)
    ap.add_argument("--timeout", type=int, default=pt.DEFAULT_TIMEOUT)
    a = ap.parse_args(argv)
    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, "reconfigure"):
            stream.reconfigure(errors="backslashreplace")
    lex = lex_mod.load()
    try:
        if not (a.snapshot and a.tags and a.proposals and a.out and a.sample):
            raise pt.UsageError("--snapshot, --tags, --proposals, --out and --sample are required, with --model or --profile-file")
        if a.sample < 1:
            raise pt.UsageError("--sample must be 1 or more")
        model, profile, think_model, think_profile = pt.model_choice(a)
        prop_run = ca.read_json(a.proposals / pt.RUN, "proposer's run file")
        proposer_models = [prop_run.get("model"), prop_run.get("think_model")]
        refuse_same_model([model, think_model], proposer_models)
        s, entries, rule_tags, _ = pt.load_inputs(a.snapshot, a.tags, a.exclude, lex)
        if prop_run.get("snapshot_sha256") != s.digest or prop_run.get("lexicon_sha256") != lex.sha256:
            raise pt.UsageError("the proposer's run was made from a different snapshot or lexicon")
        proposed = ca.read_json(a.proposals / pt.MODEL_TAGS, "proposer's model tags")
        accepted, rejected = ca.verify_model_tags(proposed, entries, rule_tags, lex)
        if rejected:
            raise pt.UsageError(f"{len(rejected)} of the proposer's model tags do not verify against this snapshot and tags file")
        refs: dict[str, list[str]] = {}
        for t in accepted:
            refs.setdefault(t["entry_id"], []).append(t["class"])
        rule_pool = [e for e, t in rule_tags.items() if t["status"] == "rule" and e in entries]
        for e in rule_pool:
            refs[e] = [t["class"] for t in rule_tags[e]["tags"]]
        sample = [list(x) for x in draw_sample(list(dict.fromkeys(t["entry_id"] for t in accepted)), rule_pool, a.sample, a.seed)]
        if invoker is None:
            invoker = pt.PwshInvoker(pt.resolve_invoker(a.invoker, a.lab_repo))
        binding = {"kind": "trial-atlas challenger", "snapshot_sha256": s.digest, "lexicon_version": lex.version,
                   "lexicon_sha256": lex.sha256, "tags_sha256": pt.sha256_bytes(a.tags),
                   "model_tags_sha256": pt.sha256_bytes(a.proposals / pt.MODEL_TAGS),
                   "proposer_models": [m for m in dict.fromkeys(proposer_models) if m], "model": model,
                   "profile": Path(profile).name if profile else None, "think_model": think_model,
                   "think_profile": Path(think_profile).name if think_profile else None, "template": pt.TEMPLATE_VERSION,
                   "seed": a.seed, "temperature": a.temperature, "sample_requested": a.sample, "sample": sample}
        pt.prepare_out(a.out, a.resume, binding, run_name=RUN, jsonl_name=JSONL)
    except (pt.UsageError, ca.UsageError, OSError, ValueError) as exc:
        print(f"ERROR: {snap.clean(exc, 300)}")
        return 2
    except snap.SnapshotError as exc:
        print(f"ERROR: the snapshot was refused: {snap.clean(exc, 300)}")
        return 2
    # The challenger's quote check: the same verifier, with each sampled entry treated as unclassified (a rule-tagged entry is
    # otherwise refused as "already decided"). This view is used for nothing else.
    open_view = {eid: {"status": "unclassified"} for eid, _ in sample}
    cfg = pt.Settings(model=model, profile=profile, think_model=think_model, think_profile=think_profile, temperature=a.temperature,
                      seed=a.seed, timeout_sec=a.timeout, think_after_fast=a.think_after_fast, extra={"role": "challenger"})
    print(f"challenger {snap.clean(model, 60)}: {len(sample)} entr(ies) drawn "
          f"({sum(1 for _, p in sample if p == 'model')} model, {sum(1 for _, p in sample if p == 'rule')} rule)")
    code = 0
    try:
        pt.run_entries([e for e, _ in sample], entries, open_view, lex, invoker, a.out / JSONL, cfg)
    except KeyboardInterrupt:
        print("interrupted: challenge.jsonl holds every finished attempt; run again with --resume")
        code = 130
    rows, _ = pt.load_rows(a.out / JSONL)
    res = result(rows, sample, refs, binding)
    pt.write_json(a.out / OUT_JSON, res)
    o = res["agreement"]["overall"]
    shown = "not measured" if o["rate"] is None else f"{o['rate']:.3f} ({o['agree']} of {o['compared']}; Wilson 95% {o['wilson95']})"
    print(f"agreement {shown}; unverified {res['unverified']}; disagreements {len(res['disagreements'])}; see {OUT_JSON}")
    return code


if __name__ == "__main__":
    sys.exit(main())
