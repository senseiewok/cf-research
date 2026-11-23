#!/usr/bin/env python3
"""Independent evaluation harness for CF ML readiness extraction (Council Item #4).

Implements design-council invariants for evaluating `tools/readiness/extractor.py`:
1. Keeps `cf_ml_benchmark.json` untouched as the frozen regression baseline.
2. Evaluates factual extraction (`cohort_size`, `split_type`, `external_validation`,
   `headline_metric`, verbatim sentence grounding) separately from heuristic
   `readiness_tier` (`T0`-`T3`).
3. Enforces explicit per-field matching policies (exact numeric, closed enum, and
   normalized metric string / accepted alias list) without free-text substring guessing.
4. Reports positive-assertion precision with exact 95% Wilson score confidence
   intervals, statistical power status (`N >= 73` for zero-error >= 0.95 lower bound),
   recall on stated ground truth, answer coverage, abstention rate, verbatim grounding
   rate, and structured error taxonomy.
5. Provides a version-pinned function-to-function comparator against the pre-PR #388
   `extract_external_validation` logic (`d8acab74d94c6d5bfdf1b4c2cbbf6769da4bb830`).
6. Provides a `readiness-eval-v1` schema validator, dev-set study-group/DOI/PMID leakage
   checker, and dual-annotator agreement / nominal Cohen's kappa comparison tool.

Pure Python standard library only; zero network calls.
"""
from __future__ import annotations

import argparse
import hashlib
import inspect
import json
import math
import re
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from tools.readiness.extractor import (  # noqa: E402
    ReadinessExtractor,
    determine_readiness_tier,
    extract_external_validation,
    sanitize_input_text,
    tokenize_sentences,
)

READINESS_DIR = Path(__file__).resolve().parent
DEFAULT_STRESS_PATH = READINESS_DIR / "stress_eval_synthetic.json"
DEFAULT_DEV_BENCHMARK_PATH = READINESS_DIR / "cf_ml_benchmark.json"

SCHEMA_VERSION = "readiness-eval-v1"
FACTUAL_FIELDS: Tuple[str, ...] = (
    "cohort_size",
    "split_type",
    "external_validation",
    "headline_metric",
)
ALL_EVAL_FIELDS: Tuple[str, ...] = FACTUAL_FIELDS + ("readiness_tier",)

SENTENCE_FIELDS: Dict[str, str] = {
    "cohort_size": "cohort_sentence",
    "split_type": "split_sentence",
    "external_validation": "external_validation_sentence",
    "headline_metric": "metric_sentence",
}

MATCHING_POLICIES: Dict[str, str] = {
    "cohort_size": "exact_int_or_not_stated",
    "split_type": "exact_closed_enum",
    "external_validation": "exact_closed_enum",
    "headline_metric": "normalized_text_or_accepted_aliases",
    "readiness_tier": "exact_closed_enum",
}

CLOSED_ENUM_FIELDS: Dict[str, Tuple[str, ...]] = {
    "split_type": (
        "patient_level",
        "sample_or_record_level",
        "unsupervised_or_zero_shot",
        "not_stated",
    ),
    "external_validation": (
        "external_multicenter_or_independent",
        "absent_single_center_or_internal_only",
        "not_stated",
    ),
    "readiness_tier": (
        "T0_PRECLINICAL_BENCHMARK",
        "T1_RETROSPECTIVE_OR_UNVALIDATED",
        "T2_INTERNALLY_VALIDATED",
        "T3_EXTERNALLY_VALIDATED",
    ),
}

TIER_ORDER: Dict[str, int] = {
    "T0_PRECLINICAL_BENCHMARK": 0,
    "T1_RETROSPECTIVE_OR_UNVALIDATED": 1,
    "T2_INTERNALLY_VALIDATED": 2,
    "T3_EXTERNALLY_VALIDATED": 3,
}

# Exact pre-PR #388 commit on senseiewok/cf-research (`16c9ff67...~1`).
BASELINE_GIT_COMMIT = "d8acab74d94c6d5bfdf1b4c2cbbf6769da4bb830"
BASELINE_FUNCTION_NAME = "extract_external_validation"
WILSON_Z_95 = 1.959963984540054
PRECISION_TARGET = 0.95
MIN_KAPPA_N = 30
ABSTAINED_TOKENS = frozenset({"", "not stated", "not_stated"})
UNICODE_DASH_RE = re.compile(r"[\u2012\u2013\u2014\u2015\u2212]")


def _wilson_bounds_raw(
    successes: int, total: int, z: float = WILSON_Z_95
) -> Optional[Tuple[float, float]]:
    if total < 0 or successes < 0 or successes > total:
        raise ValueError(f"Invalid counts for Wilson interval: successes={successes}, total={total}")
    if z <= 0:
        raise ValueError(f"z must be positive, got {z}")
    if total == 0:
        return None
    p = successes / total
    z2 = z * z
    denom = 1.0 + z2 / total
    center = (p + z2 / (2.0 * total)) / denom
    radicand = (p * (1.0 - p) + z2 / (4.0 * total)) / total
    margin = (z * math.sqrt(max(0.0, radicand))) / denom
    return (max(0.0, center - margin), min(1.0, center + margin))


def wilson_score_interval(
    successes: int, total: int, z: float = WILSON_Z_95
) -> Optional[Tuple[float, float]]:
    """Compute exact Wilson score confidence interval (rounded to 4 decimals)."""
    raw = _wilson_bounds_raw(successes, total, z=z)
    if raw is None:
        return None
    return (round(raw[0], 4), round(raw[1], 4))


def min_n_for_wilson_lower_bound(
    target: float = PRECISION_TARGET, z: float = WILSON_Z_95
) -> int:
    """Return smallest N where 0 errors (k=N) achieves Wilson lower bound >= target."""
    if not (0.0 < target < 1.0):
        raise ValueError(f"target must be in (0, 1), got {target}")
    if z <= 0:
        raise ValueError(f"z must be positive, got {z}")
    approx = max(1, math.ceil((z * z * target) / (1.0 - target)))
    # Adjust locally around the analytical closed form to eliminate float boundary error.
    n = max(1, approx - 2)
    while True:
        bounds = _wilson_bounds_raw(n, n, z=z)
        if bounds is not None and bounds[0] >= target:
            return n
        n += 1


MIN_N_PRECISION_95 = min_n_for_wilson_lower_bound(PRECISION_TARGET, WILSON_Z_95)


def is_abstained(field: str, value: Any) -> bool:
    """Return True if value represents explicit abstention / silence for the given field."""
    if value is None:
        return True
    if field == "cohort_size" and isinstance(value, int) and not isinstance(value, bool):
        return False
    return str(value).strip().lower() in ABSTAINED_TOKENS


def normalize_metric_text(value: Any) -> str:
    """Normalize headline metric strings deterministically without fuzzy substring matching."""
    if value is None:
        return ""
    text = str(value).strip().lower()
    if text in ABSTAINED_TOKENS:
        return "not stated"
    text = UNICODE_DASH_RE.sub("-", text)
    text = re.sub(r"\s*±\s*", "±", text)
    text = re.sub(r"\s+", " ", text)
    return text.strip()


def values_match(
    field: str,
    got: Any,
    expected: Any,
    accepted_aliases: Optional[List[str]] = None,
) -> bool:
    """Check value equivalence under the field's declared matching policy."""
    if field == "cohort_size":
        if is_abstained(field, got) and is_abstained(field, expected):
            return True
        if isinstance(got, bool) or isinstance(expected, bool):
            return False
        if isinstance(got, int) and isinstance(expected, int):
            return got == expected
        return False

    if field in CLOSED_ENUM_FIELDS:
        return str(got) == str(expected)

    if field == "headline_metric":
        got_abs = is_abstained(field, got)
        exp_abs = is_abstained(field, expected)
        if got_abs and exp_abs:
            return True
        if got_abs or exp_abs:
            return False
        norm_got = normalize_metric_text(got)
        candidates = [normalize_metric_text(expected)]
        if accepted_aliases:
            candidates.extend(normalize_metric_text(a) for a in accepted_aliases)
        return norm_got in candidates

    raise ValueError(f"Unsupported field for matching: {field!r}")


def legacy_pre_pr388_extract_external_validation(sentences: List[str]) -> Tuple[str, str]:
    """Frozen pre-PR #388 `extract_external_validation` (commit d8acab74d94c6d5bfdf1b4c2cbbf6769da4bb830).

    Kept verbatim in structure for regression/delta comparison against the hardened PR #388 logic.
    """
    for sent in sentences:
        low = sent.lower()
        if (
            "externally validated" in low
            or "independent data" in low
            or "independent validation" in low
            or "external validation cohort" in low
        ) and "needed" not in low and "requires" not in low:
            return "external_multicenter_or_independent", sent

        if "single-center" in low or "single center" in low:
            return "absent_single_center_or_internal_only", sent
        if "external validation" in low and (
            "needed" in low
            or "lacks" in low
            or "required" in low
            or "future" in low
            or "before clinical use" in low
        ):
            return "absent_single_center_or_internal_only", sent

    return "not_stated", ""


BASELINE_FUNCTION_SHA256 = hashlib.sha256(
    inspect.getsource(legacy_pre_pr388_extract_external_validation).encode("utf-8")
).hexdigest()


def _verify_ext_val_contract(result: Any) -> None:
    if (
        not isinstance(result, tuple)
        or len(result) != 2
        or result[0] not in CLOSED_ENUM_FIELDS["external_validation"]
        or not isinstance(result[1], str)
    ):
        raise ValueError(f"Invalid external_validation extractor contract return: {result!r}")


def _unpack_dataset(raw: Any) -> Tuple[Dict[str, Any], List[Dict[str, Any]]]:
    if isinstance(raw, dict) and isinstance(raw.get("items"), list):
        meta = {
            "schema_version": raw.get("schema_version", "unspecified"),
            "dataset_kind": raw.get("dataset_kind", "unspecified"),
            "is_unseen_published_papers": bool(raw.get("is_unseen_published_papers", False)),
            "annotator_id": raw.get("annotator_id"),
            "epistemic_notice": raw.get("epistemic_notice"),
        }
        return meta, raw["items"]
    if isinstance(raw, list):
        meta = {
            "schema_version": "legacy-list",
            "dataset_kind": "legacy_list_fixture",
            "is_unseen_published_papers": False,
            "annotator_id": None,
            "epistemic_notice": None,
        }
        return meta, raw
    raise ValueError("Dataset must be a 'readiness-eval-v1' dict with 'items' or a list of items.")


def _normalize_doi(raw: Any) -> str:
    text = str(raw or "").strip().lower()
    for prefix in ("https://doi.org/", "http://doi.org/", "https://dx.doi.org/", "http://dx.doi.org/", "doi:"):
        if text.startswith(prefix):
            text = text[len(prefix):].strip()
    return text


def _normalize_pmid(raw: Any) -> str:
    text = str(raw or "").strip()
    return re.sub(r"^pmid\s*:\s*", "", text, flags=re.IGNORECASE).strip()


def check_dev_leakage(
    items: List[Dict[str, Any]],
    dev_benchmark_path: Optional[Path] = DEFAULT_DEV_BENCHMARK_PATH,
) -> Dict[str, Any]:
    """Verify no paper_id, DOI, PMID, or study_group_id overlaps with the frozen dev benchmark."""
    if dev_benchmark_path is None or not Path(dev_benchmark_path).is_file():
        return {
            "checked": False,
            "dev_benchmark_path": str(dev_benchmark_path) if dev_benchmark_path else None,
            "overlap_detected": False,
            "overlapping_paper_ids": [],
            "overlapping_dois": [],
            "overlapping_pmids": [],
            "overlapping_study_groups": [],
        }

    dev_raw = json.loads(Path(dev_benchmark_path).read_text(encoding="utf-8"))
    _, dev_items = _unpack_dataset(dev_raw)

    dev_ids = {str(it.get("paper_id") or it.get("id") or "").strip() for it in dev_items} - {""}
    dev_dois = {_normalize_doi(it.get("doi")) for it in dev_items} - {""}
    dev_pmids = {_normalize_pmid(it.get("pmid")) for it in dev_items} - {""}
    dev_groups = {
        str(it.get("study_group_id") or it.get("paper_id") or it.get("id") or "").strip()
        for it in dev_items
    } - {""}

    overlap_ids: List[str] = []
    overlap_dois: List[str] = []
    overlap_pmids: List[str] = []
    overlap_groups: List[str] = []

    for it in items:
        pid = str(it.get("paper_id") or it.get("id") or "").strip()
        doi = _normalize_doi(it.get("doi"))
        pmid = _normalize_pmid(it.get("pmid"))
        gid = str(it.get("study_group_id") or "").strip()
        if pid and pid in dev_ids:
            overlap_ids.append(pid)
        if doi and doi in dev_dois:
            overlap_dois.append(doi)
        if pmid and pmid in dev_pmids:
            overlap_pmids.append(pmid)
        if gid and gid in dev_groups:
            overlap_groups.append(gid)

    overlap_detected = bool(overlap_ids or overlap_dois or overlap_pmids or overlap_groups)
    return {
        "checked": True,
        "dev_benchmark_path": str(dev_benchmark_path),
        "dev_item_count": len(dev_items),
        "overlap_detected": overlap_detected,
        "overlapping_paper_ids": sorted(set(overlap_ids)),
        "overlapping_dois": sorted(set(overlap_dois)),
        "overlapping_pmids": sorted(set(overlap_pmids)),
        "overlapping_study_groups": sorted(set(overlap_groups)),
    }


def _build_field_tracker() -> Dict[str, Any]:
    return {
        "positive_assertions": 0,
        "correct_positive_assertions": 0,
        "stated_ground_truth": 0,
        "recalled_correct": 0,
        "abstention_count": 0,
        "true_abstentions": 0,
        "sentences_emitted": 0,
        "sentences_grounded": 0,
        "error_taxonomy": {
            "false_positive_assertion": 0,
            "false_negative_miss": 0,
            "value_mismatch": 0,
            "ungrounded_quote": 0,
            "grounding_span_mismatch": 0,
        },
    }


def _finalize_field_stats(
    field: str, tracker: Dict[str, Any], total_items: int
) -> Dict[str, Any]:
    pos = tracker["positive_assertions"]
    corr_pos = tracker["correct_positive_assertions"]
    stated = tracker["stated_ground_truth"]
    rec_corr = tracker["recalled_correct"]
    abst = tracker["abstention_count"]
    true_abst = tracker["true_abstentions"]
    sent_emitted = tracker["sentences_emitted"]
    sent_grounded = tracker["sentences_grounded"]

    precision = round(corr_pos / pos, 4) if pos > 0 else None
    ci = wilson_score_interval(corr_pos, pos) if pos > 0 else None
    if pos == 0:
        power_status = "no_positive_assertions"
    elif pos < MIN_N_PRECISION_95:
        power_status = f"underpowered_n_lt_{MIN_N_PRECISION_95}"
    elif ci is not None and ci[0] >= PRECISION_TARGET:
        power_status = "meets_lower_bound_target"
    else:
        power_status = "below_lower_bound_target"

    return {
        "matching_policy": MATCHING_POLICIES[field],
        "positive_assertions": pos,
        "correct_positive_assertions": corr_pos,
        "positive_assertion_precision": precision,
        "wilson_95_ci": list(ci) if ci is not None else None,
        "precision_target": PRECISION_TARGET,
        "min_n_for_95_lower_bound": MIN_N_PRECISION_95,
        "precision_target_status": power_status,
        "stated_ground_truth": stated,
        "recalled_correct": rec_corr,
        "recall_on_stated": round(rec_corr / stated, 4) if stated > 0 else None,
        "answer_coverage": round(pos / total_items, 4) if total_items > 0 else None,
        "abstention_count": abst,
        "abstention_rate": round(abst / total_items, 4) if total_items > 0 else None,
        "true_abstentions": true_abst,
        "mutual_abstention_rate": round(true_abst / total_items, 4) if total_items > 0 else None,
        "sentences_emitted": sent_emitted,
        "sentences_grounded": sent_grounded,
        "verbatim_grounding_rate": (
            round(sent_grounded / sent_emitted, 4) if sent_emitted > 0 else None
        ),
        "error_taxonomy": dict(tracker["error_taxonomy"]),
    }


def evaluate_dataset(
    dataset: Any,
    *,
    compare_baseline: bool = False,
    dev_benchmark_path: Optional[Path] = DEFAULT_DEV_BENCHMARK_PATH,
    extractor: Optional[ReadinessExtractor] = None,
) -> Dict[str, Any]:
    """Evaluate factual extraction separately from heuristic readiness_tier scoring."""
    meta, items = _unpack_dataset(dataset)
    if not items:
        raise ValueError("Cannot evaluate an empty dataset.")

    seen_ids = set()
    for idx, item in enumerate(items):
        pid = str(item.get("paper_id") or item.get("id") or "").strip()
        if not pid:
            raise ValueError(f"Item at index {idx} missing paper_id/id.")
        if pid in seen_ids:
            raise ValueError(f"Duplicate paper_id detected: {pid!r}")
        seen_ids.add(pid)
        if "text" not in item or not isinstance(item.get("expected"), dict):
            raise ValueError(f"Item {pid!r} must include 'text' and an 'expected' dict.")
        exp_obj = item["expected"]
        for req in ALL_EVAL_FIELDS:
            if req not in exp_obj:
                raise ValueError(f"Item {pid!r} missing expected.{req!r}.")
        for enum_field, allowed in CLOSED_ENUM_FIELDS.items():
            if exp_obj[enum_field] not in allowed:
                raise ValueError(
                    f"Item {pid!r} expected.{enum_field}={exp_obj[enum_field]!r} not in {allowed}."
                )

    ext = extractor or ReadinessExtractor()
    trackers = {f: _build_field_tracker() for f in FACTUAL_FIELDS}

    tier_total = 0
    tier_exact = 0
    tier_overclaim = 0
    tier_t3_overclaim = 0
    tier_underclaim = 0

    baseline_trackers = _build_field_tracker()
    baseline_t3_overclaim = 0
    current_t3_overclaim = 0
    comparator_item_deltas: List[Dict[str, Any]] = []
    per_item: List[Dict[str, Any]] = []

    for item in items:
        pid = str(item.get("paper_id") or item.get("id")).strip()
        text = str(item["text"])
        sanitized = sanitize_input_text(text)
        rec = ext.extract(text, paper_id=pid)
        rec_dict = rec.to_dict()
        exp = item["expected"]
        accepted_metrics = item.get("accepted_headline_metrics") or []
        grounding_substrings = item.get("expected_grounding_substrings") or {}

        item_factual_match = True
        item_errors: Dict[str, str] = {}

        for field in FACTUAL_FIELDS:
            tr = trackers[field]
            got_val = rec_dict[field]
            exp_val = exp[field]
            aliases = accepted_metrics if field == "headline_metric" else None
            got_abs = is_abstained(field, got_val)
            exp_abs = is_abstained(field, exp_val)
            matched = values_match(field, got_val, exp_val, accepted_aliases=aliases)

            if got_abs:
                tr["abstention_count"] += 1
                if exp_abs:
                    tr["true_abstentions"] += 1
            else:
                tr["positive_assertions"] += 1
                if matched:
                    tr["correct_positive_assertions"] += 1

            if not exp_abs:
                tr["stated_ground_truth"] += 1
                if not got_abs and matched:
                    tr["recalled_correct"] += 1

            if got_abs and not exp_abs:
                tr["error_taxonomy"]["false_negative_miss"] += 1
                item_factual_match = False
                item_errors[field] = "false_negative_miss"
            elif not got_abs and exp_abs:
                tr["error_taxonomy"]["false_positive_assertion"] += 1
                item_factual_match = False
                item_errors[field] = "false_positive_assertion"
            elif not got_abs and not exp_abs and not matched:
                tr["error_taxonomy"]["value_mismatch"] += 1
                item_factual_match = False
                item_errors[field] = "value_mismatch"

            sent_key = SENTENCE_FIELDS[field]
            got_sent = str(rec_dict.get(sent_key) or "")
            if got_sent:
                tr["sentences_emitted"] += 1
                if got_sent in sanitized:
                    tr["sentences_grounded"] += 1
                else:
                    tr["error_taxonomy"]["ungrounded_quote"] += 1
                    item_factual_match = False
                    item_errors[f"{field}_grounding"] = "ungrounded_quote"
                expected_sub = grounding_substrings.get(field)
                if expected_sub and expected_sub.lower() not in got_sent.lower():
                    tr["error_taxonomy"]["grounding_span_mismatch"] += 1
                    item_factual_match = False
                    item_errors[f"{field}_span"] = "grounding_span_mismatch"
            elif not got_abs:
                # Asserted field without any source sentence is treated as ungrounded span failure.
                tr["sentences_emitted"] += 1
                tr["error_taxonomy"]["ungrounded_quote"] += 1
                item_factual_match = False
                item_errors[f"{field}_grounding"] = "ungrounded_quote"

        # Evaluate heuristic readiness_tier in its own isolated section.
        tier_total += 1
        got_tier = rec.readiness_tier
        exp_tier = exp["readiness_tier"]
        if got_tier == exp_tier:
            tier_exact += 1
        elif TIER_ORDER[got_tier] > TIER_ORDER[exp_tier]:
            tier_overclaim += 1
            if got_tier == "T3_EXTERNALLY_VALIDATED" and exp_tier != "T3_EXTERNALLY_VALIDATED":
                tier_t3_overclaim += 1
        else:
            tier_underclaim += 1

        if compare_baseline:
            sentences = tokenize_sentences(text)
            base_res = legacy_pre_pr388_extract_external_validation(sentences)
            _verify_ext_val_contract(base_res)
            cur_res = extract_external_validation(sentences)
            _verify_ext_val_contract(cur_res)

            base_label, base_sent = base_res
            cur_label, _ = cur_res
            exp_ext = exp["external_validation"]

            b_got_abs = is_abstained("external_validation", base_label)
            b_exp_abs = is_abstained("external_validation", exp_ext)
            b_match = values_match("external_validation", base_label, exp_ext)

            if b_got_abs:
                baseline_trackers["abstention_count"] += 1
                if b_exp_abs:
                    baseline_trackers["true_abstentions"] += 1
            else:
                baseline_trackers["positive_assertions"] += 1
                if b_match:
                    baseline_trackers["correct_positive_assertions"] += 1

            if not b_exp_abs:
                baseline_trackers["stated_ground_truth"] += 1
                if not b_got_abs and b_match:
                    baseline_trackers["recalled_correct"] += 1

            if b_got_abs and not b_exp_abs:
                baseline_trackers["error_taxonomy"]["false_negative_miss"] += 1
            elif not b_got_abs and b_exp_abs:
                baseline_trackers["error_taxonomy"]["false_positive_assertion"] += 1
            elif not b_got_abs and not b_exp_abs and not b_match:
                baseline_trackers["error_taxonomy"]["value_mismatch"] += 1

            if base_sent:
                baseline_trackers["sentences_emitted"] += 1
                if base_sent in sanitized:
                    baseline_trackers["sentences_grounded"] += 1
                else:
                    baseline_trackers["error_taxonomy"]["ungrounded_quote"] += 1

            base_tier = determine_readiness_tier(rec.cohort_size, rec.split_type, base_label)
            if base_tier == "T3_EXTERNALLY_VALIDATED" and exp_tier != "T3_EXTERNALLY_VALIDATED":
                baseline_t3_overclaim += 1
            if got_tier == "T3_EXTERNALLY_VALIDATED" and exp_tier != "T3_EXTERNALLY_VALIDATED":
                current_t3_overclaim += 1

            if base_label != cur_label or base_tier != got_tier:
                comparator_item_deltas.append(
                    {
                        "paper_id": pid,
                        "expected_external_validation": exp_ext,
                        "baseline_external_validation": base_label,
                        "current_external_validation": cur_label,
                        "expected_tier": exp_tier,
                        "baseline_tier": base_tier,
                        "current_tier": got_tier,
                    }
                )

        per_item.append(
            {
                "paper_id": pid,
                "study_group_id": item.get("study_group_id"),
                "category": item.get("category"),
                "factual_all_match": item_factual_match,
                "factual_errors": item_errors,
                "predicted": rec_dict,
                "expected": exp,
            }
        )

    n_items = len(items)
    factual_fields_report = {
        field: _finalize_field_stats(field, trackers[field], n_items)
        for field in FACTUAL_FIELDS
    }

    all_grounded_emitted = sum(tr["sentences_emitted"] for tr in trackers.values())
    all_grounded_ok = sum(tr["sentences_grounded"] for tr in trackers.values())

    report: Dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "dataset_kind": meta["dataset_kind"],
        "is_unseen_published_papers": meta["is_unseen_published_papers"],
        "epistemic_notice": meta["epistemic_notice"],
        "items_evaluated": n_items,
        "dev_leakage_check": check_dev_leakage(items, dev_benchmark_path),
        "factual_extraction": {
            "evaluated_independently_from_readiness_tier": True,
            "fields": factual_fields_report,
            "overall_sentences_emitted": all_grounded_emitted,
            "overall_sentences_grounded": all_grounded_ok,
            "overall_verbatim_grounding_rate": (
                round(all_grounded_ok / all_grounded_emitted, 4)
                if all_grounded_emitted > 0
                else None
            ),
        },
        "heuristic_readiness_tier": {
            "heuristic_derived_from_factual_fields": True,
            "total": tier_total,
            "exact_matches": tier_exact,
            "exact_match_rate": round(tier_exact / tier_total, 4),
            "overclaim_count": tier_overclaim,
            "overclaim_rate": round(tier_overclaim / tier_total, 4),
            "t3_overclaim_count": tier_t3_overclaim,
            "underclaim_count": tier_underclaim,
            "underclaim_rate": round(tier_underclaim / tier_total, 4),
        },
        "items": per_item,
    }

    if compare_baseline:
        base_ext_stats = _finalize_field_stats("external_validation", baseline_trackers, n_items)
        cur_ext_stats = factual_fields_report["external_validation"]
        base_prec = base_ext_stats["positive_assertion_precision"]
        cur_prec = cur_ext_stats["positive_assertion_precision"]
        prec_delta = (
            round(cur_prec - base_prec, 4)
            if (base_prec is not None and cur_prec is not None)
            else None
        )
        report["baseline_comparison"] = {
            "comparison_scope": "function_to_function_external_validation_on_identical_sentences",
            "baseline_git_commit": BASELINE_GIT_COMMIT,
            "baseline_function_name": BASELINE_FUNCTION_NAME,
            "baseline_function_sha256": BASELINE_FUNCTION_SHA256,
            "baseline_external_validation": base_ext_stats,
            "current_external_validation": cur_ext_stats,
            "positive_assertion_precision_delta": prec_delta,
            "baseline_t3_overclaims": baseline_t3_overclaim,
            "current_t3_overclaims": current_t3_overclaim,
            "t3_overclaim_reduction": baseline_t3_overclaim - current_t3_overclaim,
            "items_changed_count": len(comparator_item_deltas),
            "item_deltas": comparator_item_deltas,
        }

    return report


def validate_annotation_payload(payload: Any) -> List[str]:
    """Validate a dual-annotator `readiness-eval-v1` file and return structural error messages."""
    errors: List[str] = []
    if not isinstance(payload, dict):
        return ["Top-level JSON payload must be an object."]
    if payload.get("schema_version") != SCHEMA_VERSION:
        errors.append(
            f"schema_version must be {SCHEMA_VERSION!r}, got {payload.get('schema_version')!r}"
        )
    if not isinstance(payload.get("dataset_kind"), str) or not payload["dataset_kind"].strip():
        errors.append("Missing non-empty 'dataset_kind'.")
    if not isinstance(payload.get("annotator_id"), str) or not payload["annotator_id"].strip():
        errors.append("Missing non-empty 'annotator_id'.")

    items = payload.get("items")
    if not isinstance(items, list) or not items:
        errors.append("'items' must be a non-empty list.")
        return errors

    seen_pids = set()
    for idx, item in enumerate(items):
        if not isinstance(item, dict):
            errors.append(f"Item {idx} must be an object.")
            continue
        pid = item.get("paper_id")
        gid = item.get("study_group_id")
        if not isinstance(pid, str) or not pid.strip():
            errors.append(f"Item {idx} missing non-empty 'paper_id'.")
        elif pid in seen_pids:
            errors.append(f"Duplicate 'paper_id': {pid!r}")
        else:
            seen_pids.add(pid)
        if not isinstance(gid, str) or not gid.strip():
            errors.append(f"Item {idx} ({pid!r}) missing non-empty 'study_group_id'.")

        exp = item.get("expected")
        if not isinstance(exp, dict):
            errors.append(f"Item {idx} ({pid!r}) missing 'expected' dict.")
            continue

        for req_field in ALL_EVAL_FIELDS:
            if req_field not in exp:
                errors.append(f"Item {idx} ({pid!r}) missing expected.{req_field}.")
        if "cohort_size" in exp:
            cs = exp["cohort_size"]
            valid_cs = (
                (isinstance(cs, int) and not isinstance(cs, bool) and cs >= 0)
                or (isinstance(cs, str) and cs.strip().lower() == "not stated")
            )
            if not valid_cs:
                errors.append(
                    f"Item {idx} ({pid!r}) expected.cohort_size must be non-negative int or 'not stated', got {cs!r}."
                )
        for enum_field, allowed in CLOSED_ENUM_FIELDS.items():
            if enum_field in exp and exp[enum_field] not in allowed:
                errors.append(
                    f"Item {idx} ({pid!r}) expected.{enum_field}={exp[enum_field]!r} not in {allowed}."
                )
    return errors


def compute_nominal_cohens_kappa(
    labels_a: List[str], labels_b: List[str]
) -> Dict[str, Any]:
    """Compute nominal Cohen's kappa with explicit N < 30 power and zero-variance guards."""
    if len(labels_a) != len(labels_b):
        raise ValueError("Label sequences must have identical length.")
    n = len(labels_a)
    if n == 0:
        raise ValueError("Cannot compute Cohen's kappa on 0 items.")

    agree = sum(1 for a, b in zip(labels_a, labels_b) if a == b)
    p_o = agree / n

    counts_a: Dict[str, int] = {}
    counts_b: Dict[str, int] = {}
    for a, b in zip(labels_a, labels_b):
        counts_a[a] = counts_a.get(a, 0) + 1
        counts_b[b] = counts_b.get(b, 0) + 1
    categories = set(counts_a.keys()) | set(counts_b.keys())
    p_e = sum((counts_a.get(cat, 0) / n) * (counts_b.get(cat, 0) / n) for cat in categories)

    if abs(1.0 - p_e) <= 1e-12:
        return {
            "n": n,
            "agreement_scope": "all_rows",
            "observed_agreement": round(p_o, 4),
            "expected_chance_agreement": round(p_e, 4),
            "kappa": None,
            "kappa_status": "undefined_zero_variance",
            "min_n_recommended": MIN_KAPPA_N,
        }

    kappa = (p_o - p_e) / (1.0 - p_e)
    status = f"underpowered_n_lt_{MIN_KAPPA_N}" if n < MIN_KAPPA_N else "computed"
    return {
        "n": n,
        "agreement_scope": "all_rows",
        "observed_agreement": round(p_o, 4),
        "expected_chance_agreement": round(p_e, 4),
        "kappa": round(kappa, 4),
        "kappa_status": status,
        "min_n_recommended": MIN_KAPPA_N,
    }


def compare_annotators(
    payload_a: Dict[str, Any],
    payload_b: Dict[str, Any],
    *,
    dev_benchmark_path: Optional[Path] = DEFAULT_DEV_BENCHMARK_PATH,
) -> Dict[str, Any]:
    """Compare two independent `readiness-eval-v1` annotation files for agreement & leakage."""
    errs_a = validate_annotation_payload(payload_a)
    if errs_a:
        raise ValueError(f"Annotator A payload invalid: {'; '.join(errs_a)}")
    errs_b = validate_annotation_payload(payload_b)
    if errs_b:
        raise ValueError(f"Annotator B payload invalid: {'; '.join(errs_b)}")

    items_a = {it["paper_id"]: it for it in payload_a["items"]}
    items_b = {it["paper_id"]: it for it in payload_b["items"]}
    keys_a = set(items_a.keys())
    keys_b = set(items_b.keys())
    if keys_a != keys_b:
        missing_in_b = sorted(keys_a - keys_b)
        missing_in_a = sorted(keys_b - keys_a)
        raise ValueError(
            f"Annotator item sets differ: missing_in_b={missing_in_b}, missing_in_a={missing_in_a}"
        )

    ordered_pids = [it["paper_id"] for it in payload_a["items"]]
    study_group_mismatches: List[Dict[str, str]] = []
    for pid in ordered_pids:
        ga = items_a[pid]["study_group_id"]
        gb = items_b[pid]["study_group_id"]
        if ga != gb:
            study_group_mismatches.append(
                {"paper_id": pid, "annotator_a_group": ga, "annotator_b_group": gb}
            )

    leakage_a = check_dev_leakage(payload_a["items"], dev_benchmark_path)
    leakage_b = check_dev_leakage(payload_b["items"], dev_benchmark_path)

    field_agreement: Dict[str, Any] = {}
    disagreements: List[Dict[str, Any]] = []

    for field in ALL_EVAL_FIELDS:
        all_matches = 0
        both_asserted_total = 0
        both_asserted_matches = 0
        labels_a: List[str] = []
        labels_b: List[str] = []

        for pid in ordered_pids:
            va = items_a[pid]["expected"][field]
            vb = items_b[pid]["expected"][field]
            matched = values_match(field, va, vb)
            if matched:
                all_matches += 1
            else:
                disagreements.append(
                    {
                        "paper_id": pid,
                        "study_group_id": items_a[pid]["study_group_id"],
                        "field": field,
                        "annotator_a_value": va,
                        "annotator_b_value": vb,
                    }
                )

            if not is_abstained(field, va) and not is_abstained(field, vb):
                both_asserted_total += 1
                if matched:
                    both_asserted_matches += 1

            if field in CLOSED_ENUM_FIELDS:
                labels_a.append(str(va))
                labels_b.append(str(vb))

        n = len(ordered_pids)
        if field in CLOSED_ENUM_FIELDS:
            kappa_block = compute_nominal_cohens_kappa(labels_a, labels_b)
        else:
            kappa_block = {
                "n": n,
                "agreement_scope": "all_rows",
                "kappa": None,
                "kappa_status": "not_applicable_non_nominal",
                "min_n_recommended": MIN_KAPPA_N,
            }

        field_agreement[field] = {
            "matching_policy": MATCHING_POLICIES[field],
            "total_items": n,
            "raw_agreement_all_rows": round(all_matches / n, 4),
            "both_asserted_count": both_asserted_total,
            "raw_agreement_both_asserted": (
                round(both_asserted_matches / both_asserted_total, 4)
                if both_asserted_total > 0
                else None
            ),
            "cohens_kappa": kappa_block,
        }

    return {
        "schema_version": SCHEMA_VERSION,
        "annotator_a_id": payload_a["annotator_id"],
        "annotator_b_id": payload_b["annotator_id"],
        "items_compared": len(ordered_pids),
        "min_n_for_95_precision_lower_bound": MIN_N_PRECISION_95,
        "min_n_for_kappa_report": MIN_KAPPA_N,
        "study_group_id_mismatches": study_group_mismatches,
        "dev_leakage_annotator_a": leakage_a,
        "dev_leakage_annotator_b": leakage_b,
        "field_agreement": field_agreement,
        "disagreements": disagreements,
    }


def _print_text_report(report: Dict[str, Any]) -> None:
    print(
        f"Readiness Evaluation ({report['dataset_kind']}, items={report['items_evaluated']}, "
        f"unseen_published_papers={report['is_unseen_published_papers']})"
    )
    if report.get("epistemic_notice"):
        print(f"Notice: {report['epistemic_notice']}")
    print("Factual Extraction:")
    for field, stats in report["factual_extraction"]["fields"].items():
        prec_s = (
            f"{stats['positive_assertion_precision'] * 100:.1f}%"
            if stats["positive_assertion_precision"] is not None
            else "n/a"
        )
        ci = stats["wilson_95_ci"]
        ci_s = f"[{ci[0]:.4f}, {ci[1]:.4f}]" if ci else "n/a"
        rec_s = (
            f"{stats['recall_on_stated'] * 100:.1f}%"
            if stats["recall_on_stated"] is not None
            else "n/a"
        )
        print(
            f"  - {field}: prec={prec_s} CI95={ci_s} ({stats['precision_target_status']}) "
            f"recall={rec_s} cov={stats['answer_coverage'] * 100:.1f}% "
            f"abst={stats['abstention_rate'] * 100:.1f}%"
        )
    tier = report["heuristic_readiness_tier"]
    print(
        f"Heuristic Readiness Tier (kept separate): exact={tier['exact_matches']}/{tier['total']} "
        f"({tier['exact_match_rate'] * 100:.1f}%), t3_overclaims={tier['t3_overclaim_count']}"
    )
    if "baseline_comparison" in report:
        bc = report["baseline_comparison"]
        print(
            f"Pre-PR#388 Baseline ({bc['baseline_git_commit'][:10]}) Comparison: "
            f"t3_overclaims {bc['baseline_t3_overclaims']} -> {bc['current_t3_overclaims']} "
            f"(changed={bc['items_changed_count']})"
        )


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--stress",
        action="store_true",
        help="Evaluate on the labeled synthetic stress suite (stress_eval_synthetic.json).",
    )
    group.add_argument(
        "--file",
        type=Path,
        help="Evaluate a readiness-eval-v1 JSON dataset.",
    )
    group.add_argument(
        "--compare-annotators",
        nargs=2,
        type=Path,
        metavar=("ANNOTATOR_A", "ANNOTATOR_B"),
        help="Compare two independent readiness-eval-v1 annotation JSON files.",
    )
    parser.add_argument(
        "--compare-baseline",
        action="store_true",
        help="Include version-pinned pre-PR #388 extract_external_validation comparison.",
    )
    parser.add_argument(
        "--dev-benchmark",
        type=Path,
        default=DEFAULT_DEV_BENCHMARK_PATH,
        help="Frozen dev benchmark JSON to check for cohort/paper leakage.",
    )
    parser.add_argument("--json", action="store_true", help="Emit JSON output.")

    args = parser.parse_args(argv)

    if args.compare_annotators:
        path_a, path_b = args.compare_annotators
        payload_a = json.loads(path_a.read_text(encoding="utf-8"))
        payload_b = json.loads(path_b.read_text(encoding="utf-8"))
        comp = compare_annotators(
            payload_a, payload_b, dev_benchmark_path=args.dev_benchmark
        )
        if args.json:
            print(json.dumps(comp, indent=2))
        else:
            print(
                f"Annotator Comparison ({comp['annotator_a_id']} vs {comp['annotator_b_id']}, "
                f"items={comp['items_compared']}, disagreements={len(comp['disagreements'])})"
            )
            for field, stats in comp["field_agreement"].items():
                kp = stats["cohens_kappa"]
                print(
                    f"  - {field}: agree_all={stats['raw_agreement_all_rows']:.4f} "
                    f"kappa={kp.get('kappa')} ({kp['kappa_status']})"
                )
        return 0

    dataset_path = DEFAULT_STRESS_PATH if args.stress else args.file
    dataset = json.loads(dataset_path.read_text(encoding="utf-8"))
    report = evaluate_dataset(
        dataset,
        compare_baseline=args.compare_baseline,
        dev_benchmark_path=args.dev_benchmark,
    )
    if args.json:
        print(json.dumps(report, indent=2))
    else:
        _print_text_report(report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
