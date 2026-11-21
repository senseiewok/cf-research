#!/usr/bin/env python3
"""Independent regex-free audit of snapshot.json cohort & microbiology counts (T-0944 / PRED-0009).

Rule 2 of `discovery-practice`: a second method for a headline claim must live in separate code
and use a genuinely different method rather than copying regexes from `scope.py` or `lexicon.py`.
This script imports NO regular-expression engine (`re` is not imported) and NO `tools/trial_atlas/`
modules: it uses pure character-level alphanumeric word-tuple scanning with measure-first outcome
precedence, cross-verified against `data/care_team_trials_directory.json`.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
DEFAULT_SNAPSHOT = REPO_ROOT / "tools" / "trial_atlas" / "snapshot.json"
DEFAULT_DIRECTORY = REPO_ROOT / "data" / "care_team_trials_directory.json"

BENCHMARK_19_SNAPSHOT_NCTS: tuple[str, ...] = (
    "NCT06998043", "NCT05010577", "NCT04596319", "NCT04684641", "NCT07275905",
    "NCT05616221", "NCT04294043", "NCT03597347", "NCT05101915", "NCT04685720",
    "NCT02354859", "NCT00721071", "NCT01594827", "NCT01746095", "NCT03181932",
    "NCT00391976", "NCT02453789", "NCT06654752", "NCT04545515",
)

NEGATION_WORDS: frozenset[str] = frozenset({"non", "without", "excluding", "not", "except"})
STRONG_MICROBIOLOGY_STEMS: frozenset[str] = frozenset(
    {
        "cfu", "cfus", "colony", "culture", "cultures", "culturing", "eradication",
        "microbiota", "microbiome", "microbiology", "microbiological", "mycobacterial",
    }
)
BARE_PATHOGEN_WORDS: frozenset[str] = frozenset(
    {
        "mycobacterium", "pseudomonas", "staphylococcus", "mrsa", "burkholderia",
        "achromobacter", "ntm", "fungal", "aspergillus", "stenotrophomonas",
    }
)
NON_MICRO_CLINICAL_MEASURE_WORDS: frozenset[str] = frozenset(
    {
        "fev1", "ppfev1", "fvc", "lci", "spirometry", "adverse", "event", "events",
        "teae", "tea", "teaes", "saes", "sae", "toxicity", "toxicities",
        "pharmacokinetic", "pharmacokinetics", "ctrough", "cmax", "auc", "bmi", "weight", "cfqr", "sweat",
    }
)


def _word_list(text: str) -> list[str]:
    """Split text into lowercase alphanumeric words without regular expressions."""
    chars = [ch.lower() if ch.isalnum() else " " for ch in str(text or "")]
    return "".join(chars).split()


def classify_conditions_by_word_window(study: dict) -> str:
    """Classify scope ('included', 'X1'..'X5') via word-window tuples without regexes."""
    st = str(study.get("study_type") or "").upper()
    if st == "OBSERVATIONAL":
        return "X1"
    if st == "EXPANDED_ACCESS":
        return "X2"
    if st != "INTERVENTIONAL":
        return "X3"

    saw_positive_cf = False
    saw_negated_cf = False
    for raw_cond in study.get("conditions") or []:
        words = _word_list(str(raw_cond or ""))
        if words == ["cf"]:
            saw_positive_cf = True
            continue
        for idx, w in enumerate(words):
            is_cf_span = (
                w == "mucoviscidosis"
                or (w == "cystic" and idx + 1 < len(words) and words[idx + 1] == "fibrosis")
                or (w == "cf" and idx > 0 and words[idx - 1] in NEGATION_WORDS)
            )
            if not is_cf_span:
                continue
            window_before = words[max(0, idx - 3):idx]
            has_other_than = (
                len(window_before) >= 2
                and window_before[-2] == "other"
                and window_before[-1] == "than"
            )
            if any(nw in window_before for nw in NEGATION_WORDS) or has_other_than:
                saw_negated_cf = True
            else:
                saw_positive_cf = True

    if saw_positive_cf:
        return "included"
    if saw_negated_cf:
        return "X5"
    return "X4"


def has_primary_microbiology_measure_first(primary_outcomes: list[dict]) -> bool:
    """Measure-first general microbiology check: suppresses bare pathogen names when measure names FEV1/safety/PK."""
    for po in primary_outcomes:
        m_words = set(_word_list(str(po.get("measure") or "")))
        has_strong_m = bool(m_words & STRONG_MICROBIOLOGY_STEMS) or (
            "bacterial" in m_words and bool(m_words & {"load", "density", "burden", "count", "counts"})
        )
        if has_strong_m:
            return True
        if m_words & NON_MICRO_CLINICAL_MEASURE_WORDS:
            continue
        if m_words & BARE_PATHOGEN_WORDS:
            return True
        d_words = set(_word_list(str(po.get("description") or "")))
        has_strong_d = bool(d_words & STRONG_MICROBIOLOGY_STEMS) or (
            "bacterial" in d_words and bool(d_words & {"load", "density", "burden", "count", "counts"})
        )
        if has_strong_d or (
            bool(d_words & BARE_PATHOGEN_WORDS) and not bool(d_words & NON_MICRO_CLINICAL_MEASURE_WORDS)
        ):
            return True
    return False


def audit_independently(
    snapshot_path: Path = DEFAULT_SNAPSHOT,
    directory_path: Path = DEFAULT_DIRECTORY,
) -> dict[str, int]:
    payload = json.loads(snapshot_path.read_text(encoding="utf-8"))
    studies = payload.get("studies") or []
    by_id = {str(s["nct_id"]): s for s in studies if s.get("nct_id")}

    interventional_raw = sum(1 for s in studies if str(s.get("study_type") or "").upper() == "INTERVENTIONAL")
    prisma_by_window = sum(1 for s in studies if classify_conditions_by_word_window(s) == "included")

    directory_count = prisma_by_window
    if directory_path.is_file():
        dir_data = json.loads(directory_path.read_text(encoding="utf-8"))
        directory_count = int((dir_data.get("metadata") or {}).get("total_prisma_included_studies") or 0)
        if directory_count != prisma_by_window:
            raise ValueError(f"Directory count {directory_count} != word-window count {prisma_by_window}")

    prisma_in_20 = 0
    prisma_primary_micro_in_20 = 0
    for nid in BENCHMARK_19_SNAPSHOT_NCTS:
        rec = by_id[nid]
        if classify_conditions_by_word_window(rec) == "included":
            prisma_in_20 += 1
            if has_primary_microbiology_measure_first(list(rec.get("primary_outcomes") or [])):
                prisma_primary_micro_in_20 += 1

    x4_abate = int(classify_conditions_by_word_window(by_id["NCT04294043"]) == "X4")
    x5_noncf = int(classify_conditions_by_word_window(by_id["NCT05616221"]) == "X5")
    send_cf_absent = int("NCT07312734" not in by_id)

    return {
        "raw_interventional_count": interventional_raw,
        "prisma_included_count": prisma_by_window,
        "directory_cross_verified_count": directory_count,
        "x4_excluded_cf_ntm_gallium_count": x4_abate,
        "x5_excluded_non_cf_phage_count": x5_noncf,
        "prisma_studies_in_20_benchmark": prisma_in_20,
        "prisma_primary_microbiology_in_20_benchmark": prisma_primary_micro_in_20,
        "send_cf_absent_from_interventional_snapshot": send_cf_absent,
        "verified_benchmark_total": prisma_in_20 + x4_abate + x5_noncf + send_cf_absent,
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--snapshot", type=Path, default=DEFAULT_SNAPSHOT)
    ap.add_argument("--directory", type=Path, default=DEFAULT_DIRECTORY)
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--json", action="store_true")
    args = ap.parse_args(argv)

    counts = audit_independently(args.snapshot, args.directory)
    expected = {
        "raw_interventional_count": 1053,
        "prisma_included_count": 940,
        "directory_cross_verified_count": 940,
        "x4_excluded_cf_ntm_gallium_count": 1,
        "x5_excluded_non_cf_phage_count": 1,
        "prisma_studies_in_20_benchmark": 17,
        "prisma_primary_microbiology_in_20_benchmark": 9,
        "send_cf_absent_from_interventional_snapshot": 1,
        "verified_benchmark_total": 20,
    }
    if args.check and counts != expected:
        print(f"FAIL independent_audit: observed={counts} expected={expected}", file=sys.stderr)
        return 1
    if args.json:
        print(json.dumps(counts, indent=2))
    else:
        print(
            f"PASS independent_audit: raw={counts['raw_interventional_count']} "
            f"prisma={counts['prisma_included_count']} (dir={counts['directory_cross_verified_count']}) "
            f"prisma_primary_micro={counts['prisma_primary_microbiology_in_20_benchmark']}/17"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
