#!/usr/bin/env python3
"""Apply the frozen consensus rule to the model route's proposals.

The rule (declared 2026-10-11, see proposals/2026-10-11-frozen-census-design.md): a model proposal for an unclassified entry is kept
only when a blind model of a DIFFERENT family, labelling the entry without seeing the proposal, also gave that class. A proposal
whose entry has an invalid blind reply, or is absent from the blind labels, is dropped: nothing is kept by default.

Inputs
  --model-tags FILE      the plain model-tags.json that propose_tags.py writes: a JSON list of {entry_id, class, quote, model_id}.
                         That file records each proposal's model_id and no family, so --proposer-family is always required.
  --proposer-family NAME the family of the proposing model (for example "qwen"). Required.
  --blind-labels FILE    JSON, made by the person who runs this tool from a blind labelling of the same entries by a model of
                         another family (JSON only, no CSV: the family and model are given once at the top level):
        {"labeller_family": "gemma",                         non-empty string, once
         "labeller_model":  "gemma4:31b-it-q4_K_M",          non-empty string, once
         "entries": {"NCT00000001:P1": {"classes": ["fev1", "pro_unnamed"]},     the set of class ids the labeller gave
                     "NCT00000002:P1": {"classes": []},                           the labeller answered none: a valid reply
                     "NCT00000003:P1": {"invalid": true}}}                        the reply was unusable
        Each entry value holds exactly one key: "classes" (a list of strings) or "invalid" (true). Any other shape is refused.
  --out DIR              a folder outside every git working tree (refused otherwise, as propose_tags.py does). The three
                         files below must not already exist there.

Outputs in --out
  model-tags.json        the kept proposals, in the input order, in exactly the format propose_tags.py writes (indent 1, UTF-8,
                         one trailing newline), so check_atlas.py --model-tags reads it unchanged. It carries each kept proposal's
                         quote (that format requires it), so it is private like proposals.jsonl.
  consensus-report.json  counts only: proposed, kept, dropped, dropped by reason, the two families and models, the rule text,
                         RULE_DECLARED and the sentence that the rule was chosen after the plain route's precision was seen, and the
                         hashes of the two inputs. No entry id, no entry text.
  dropped.json           the dropped proposals: entry id, class, reason (disagree | invalid_reply | missing). No entry text.

Refusals (exit 1, nothing written): the labeller's family or model is the proposer's (compared after lower-casing, and by the part
before ':' of each model id, as challenge_tags.py does); an unreadable or malformed input; --out inside a git working tree or
already holding an output file. Exit 2 is a usage error (a missing argument).

Deterministic: no clock, no randomness; the same inputs give byte-identical outputs. Standard library only. No network, no model.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import check_atlas as ca  # noqa: E402
import propose_tags as pt  # noqa: E402

RULE_DECLARED = "2026-10-11"
RULE_TEXT = ("A model proposal for an unclassified entry is kept only when a blind model of a different family, labelling the "
             "entry without seeing the proposal, also gave that class; an entry with an invalid blind reply, or missing from the "
             "blind labels, is dropped.")
RULE_ORIGIN = ("The rule was chosen after the plain route's precision was seen, so a figure measured on entries labelled before "
               "the choice is a development figure; only entries labelled after it, scored under it, confirm the rule.")
MODEL_TAGS, REPORT, DROPPED = "model-tags.json", "consensus-report.json", "dropped.json"
DISAGREE, INVALID, MISSING = "disagree", "invalid_reply", "missing"


class Refused(Exception):
    pass


def _family(text) -> str:
    return " ".join(str(text or "").split()).casefold()


def _base(model) -> str:
    return _family(model).split(":", 1)[0].strip()


def load_proposals(path: Path) -> list[dict]:
    try:
        data = ca.read_json(path, "model-tags file")
    except ca.UsageError as exc:
        raise Refused(str(exc)) from None
    if not isinstance(data, list):
        raise Refused("the model-tags file must hold a JSON list")
    for t in data:
        if not isinstance(t, dict) or not all(isinstance(t.get(k), str) and t[k] for k in ("entry_id", "class", "quote")):
            raise Refused("every model tag needs non-empty string fields entry_id, class and quote")
        if "model_id" in t and t["model_id"] is not None and not isinstance(t["model_id"], str):
            raise Refused("a model tag's model_id must be a string")
    return data


def load_labels(path: Path) -> dict:
    """{"family": str, "model": str, "entries": {entry_id: set of class ids, or None for an invalid reply}}"""
    try:
        data = ca.read_json(path, "blind-labels file")
    except ca.UsageError as exc:
        raise Refused(str(exc)) from None
    if not isinstance(data, dict):
        raise Refused("the blind-labels file must hold a JSON object")
    for k in ("labeller_family", "labeller_model"):
        if not isinstance(data.get(k), str) or not data[k].strip():
            raise Refused(f"the blind-labels file needs a non-empty string {k}")
    raw = data.get("entries")
    if not isinstance(raw, dict):
        raise Refused("the blind-labels file needs an entries object keyed by entry id")
    entries: dict[str, set | None] = {}
    for eid, v in raw.items():
        if not isinstance(v, dict) or len(v) != 1:
            raise Refused('each blind-labels entry must be {"classes": [...]} or {"invalid": true}')
        if "invalid" in v:
            if v["invalid"] is not True:
                raise Refused('"invalid" must be true')
            entries[eid] = None
        elif isinstance(v.get("classes"), list) and all(isinstance(c, str) and c for c in v["classes"]):
            entries[eid] = set(v["classes"])
        else:
            raise Refused('each blind-labels entry must be {"classes": [...]} or {"invalid": true}')
    return {"family": data["labeller_family"].strip(), "model": data["labeller_model"].strip(), "entries": entries}


def refuse_same_family(proposer_family: str, proposer_models: list[str], labeller_family: str, labeller_model: str) -> None:
    pf, lf = _family(proposer_family), _family(labeller_family)
    if not pf:
        raise Refused("--proposer-family must not be empty")
    if pf == lf:
        raise Refused(f"the blind labeller's family ({labeller_family!r}) is the proposer's; the rule needs a different family")
    for m in proposer_models:
        if _family(m) == _family(labeller_model) or _base(m) == _base(labeller_model) or _base(m) == lf:
            raise Refused("the blind labeller's model matches a model that made the proposals; the rule needs a different family")
    if _base(labeller_model) == pf:
        raise Refused("the blind labeller's model name carries the proposer's family; the rule needs a different family")


def apply_rule(proposals: list[dict], labels: dict[str, set | None]) -> tuple[list[dict], list[dict]]:
    kept, dropped = [], []
    for t in proposals:
        eid, cls = t["entry_id"], t["class"]
        if eid not in labels:
            reason = MISSING
        elif labels[eid] is None:
            reason = INVALID
        elif cls not in labels[eid]:
            reason = DISAGREE
        else:
            kept.append(t)
            continue
        dropped.append({"entry_id": eid, "class": cls, "reason": reason})
    return kept, dropped


def _sha(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def make_report(proposals, kept, dropped, proposer_family, labeller, tags_path, labels_path) -> dict:
    by = {r: sum(1 for d in dropped if d["reason"] == r) for r in (DISAGREE, INVALID, MISSING)}
    return {"kind": "trial-atlas consensus report (counts only)",
            "rule_declared": RULE_DECLARED,
            "rule": RULE_TEXT,
            "rule_origin": RULE_ORIGIN,
            "proposer_family": proposer_family.strip(),
            "proposer_models": sorted({t["model_id"] for t in proposals if t.get("model_id")}),
            "labeller_family": labeller["family"],
            "labeller_model": labeller["model"],
            "proposed": len(proposals),
            "kept": len(kept),
            "dropped": len(dropped),
            "dropped_disagree": by[DISAGREE],
            "dropped_invalid_reply": by[INVALID],
            "dropped_missing": by[MISSING],
            "model_tags_sha256": _sha(tags_path),
            "blind_labels_sha256": _sha(labels_path)}


def _dump(data) -> bytes:
    """The bytes propose_tags.write_json writes (indent 1, UTF-8, trailing newline), with \\n on every platform."""
    return (json.dumps(data, indent=1, ensure_ascii=False) + "\n").encode("utf-8")


def run(model_tags: Path, blind_labels: Path, out: Path, proposer_family: str) -> dict:
    """Check everything first; write only when every check has passed."""
    proposals = load_proposals(model_tags)
    labeller = load_labels(blind_labels)
    refuse_same_family(proposer_family, [t["model_id"] for t in proposals if t.get("model_id")], labeller["family"], labeller["model"])
    if pt.inside_git_tree(out):
        raise Refused("--out lies inside a git working tree; the kept tags carry registry quotes and stay outside every repository")
    if any((out / n).exists() for n in (MODEL_TAGS, REPORT, DROPPED)):
        raise Refused(f"{out.name} already holds an output file; choose a new folder")
    kept, dropped = apply_rule(proposals, labeller["entries"])
    dropped.sort(key=lambda d: (d["entry_id"], d["class"]))
    report = make_report(proposals, kept, dropped, proposer_family, labeller, model_tags, blind_labels)
    out.mkdir(parents=True, exist_ok=True)
    for name, data in ((MODEL_TAGS, kept), (REPORT, report), (DROPPED, dropped)):
        tmp = out / (name + ".tmp")
        tmp.write_bytes(_dump(data))
        os.replace(tmp, out / name)
    return report


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Apply the frozen consensus rule (blind labeller of another family) to model tags.")
    ap.add_argument("--model-tags", type=Path, required=True)
    ap.add_argument("--blind-labels", type=Path, required=True)
    ap.add_argument("--proposer-family", required=True,
                    help="the proposing model's family; model-tags.json records only model_id, so this is always required")
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(argv)
    try:
        r = run(a.model_tags, a.blind_labels, a.out, a.proposer_family)
    except Refused as exc:
        print(f"REFUSED: {str(exc)[:300]}")
        return 1
    print(f"proposed {r['proposed']}, kept {r['kept']}, dropped {r['dropped']} (disagree {r['dropped_disagree']}, invalid reply "
          f"{r['dropped_invalid_reply']}, missing {r['dropped_missing']}); rule declared {RULE_DECLARED}; see {REPORT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
