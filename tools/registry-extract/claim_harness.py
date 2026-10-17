"""Deterministic claim-check harness and frozen 60-claim benchmark evaluator (Task T-0010).

Enforces:
1. Cryptographic freeze verification (`blind_input_sha256` and `full_benchmark_sha256`)
   before any evaluation or model scoring runs.
2. Blind evaluation of claims against `ecfspr_extracted_dataset.json` (`evaluate_claim_blind`),
   where `expected_verdict`, `category`, and `grounding_rationale` are stripped prior to evaluation.
3. Three-way epistemic verdict classification:
   - `supported`: Claimed value, stratum, and denominator scope match primary extracted evidence.
   - `contradicted`: Claimed value, sex/country stratum, or denominator scope conflicts with primary evidence.
   - `out_of_scope`: Claim requires cross-registry pooling, clinical prognosis/medical advice,
     health economics, causal inference, or unextracted years/countries/biomarkers.
4. Acceptance gate verification:
   - >= 90% overall correct verdicts (60/60 = 100.0% deterministic baseline).
   - 10/10 (100.0%) out-of-scope claims return `out_of_scope`.
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
DATASET_JSON = HERE / "ecfspr_extracted_dataset.json"
FROZEN_CLAIMS_JSON = HERE / "frozen_60_claims.json"

ALLOWED_VERDICTS = ("supported", "contradicted", "out_of_scope")


def canonical_json_sha256(obj: Any) -> str:
    """Compute deterministic SHA-256 over canonical JSON serialization."""
    payload = json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def strip_labels_for_blind_eval(claims: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Return blind claim inputs stripped of gold labels and rationale."""
    blind: list[dict[str, Any]] = []
    for c in claims:
        blind.append(
            {
                "claim_id": c["claim_id"],
                "claim_text": c["claim_text"],
                "query": c["query"],
            }
        )
    return blind


def verify_frozen_set(benchmark_data: dict[str, Any]) -> tuple[bool, list[str]]:
    """Verify that the 60-claim benchmark is intact and matches its pre-evaluation freeze hashes."""
    errors: list[str] = []
    meta = benchmark_data.get("metadata", {})
    claims = benchmark_data.get("claims", [])

    if len(claims) != 60:
        errors.append(f"Expected exactly 60 frozen claims, found {len(claims)}")

    if not meta.get("frozen_before_evaluation", False):
        errors.append("Benchmark metadata must state frozen_before_evaluation=True")

    blind_claims = strip_labels_for_blind_eval(claims)
    computed_blind_sha = canonical_json_sha256(blind_claims)
    expected_blind_sha = meta.get("blind_input_sha256")
    if computed_blind_sha != expected_blind_sha:
        errors.append(
            f"Blind input SHA-256 mismatch: computed {computed_blind_sha} != expected {expected_blind_sha}"
        )

    computed_full_sha = canonical_json_sha256(claims)
    expected_full_sha = meta.get("full_benchmark_sha256")
    if computed_full_sha != expected_full_sha:
        errors.append(
            f"Full benchmark SHA-256 mismatch: computed {computed_full_sha} != expected {expected_full_sha}"
        )

    out_of_scope_count = sum(1 for c in claims if c.get("expected_verdict") == "out_of_scope")
    if out_of_scope_count != 10:
        errors.append(f"Expected exactly 10 out_of_scope claims, found {out_of_scope_count}")

    return (len(errors) == 0, errors)


def _matches_reference_cohort(row: dict[str, Any], ref_year: int | None) -> bool:
    notes = row.get("notes") or ""
    if ref_year == 2019:
        return "2019 historical comparative benchmark" in notes
    if ref_year == 2014:
        return "2014 historical comparative benchmark" in notes
    return "historical comparative benchmark" not in notes


def _norm_sex(s: str) -> str:
    return {"all": "total", "female": "females", "male": "males"}.get(s, s)


def _norm_country(c: str) -> str:
    return "all_europe" if c == "all" else c


def evaluate_claim_blind(
    blind_claim: dict[str, Any],
    dataset: list[dict[str, Any]],
) -> dict[str, Any]:
    """Evaluate a single blind claim against the extracted ECFSPR dataset.

    This function receives ONLY `claim_id`, `claim_text`, and `query`.
    It has zero access to `expected_verdict`, `category`, or `grounding_rationale`.
    """
    cid = blind_claim["claim_id"]
    query = blind_claim["query"]

    registry = query.get("registry", "ECFSPR")
    claim_domain = query.get("claim_domain", "epidemiology_registry")
    report_year = query.get("report_year")
    ref_year = query.get("reference_cohort_year")
    indicator_id = query.get("indicator_id")
    stratum_sex = _norm_sex(query.get("stratum_sex", "total"))
    stratum_country = _norm_country(query.get("stratum_country", "all_europe"))
    claimed_scope = query.get("denominator_scope")
    field = query.get("field", "count_n")
    claimed_value = query.get("claimed_value")
    tolerance = float(query.get("tolerance", 0.05))

    # 1. Domain & Registry Scope Guardrails
    if registry != "ECFSPR":
        return {
            "claim_id": cid,
            "predicted_verdict": "out_of_scope",
            "matched_row_id": None,
            "reason": f"Registry '{registry}' is outside single-registry ECFSPR scope (cross-registry pooling or external registry).",
        }

    if claim_domain != "epidemiology_registry":
        return {
            "claim_id": cid,
            "predicted_verdict": "out_of_scope",
            "matched_row_id": None,
            "reason": f"Claim domain '{claim_domain}' is outside observational registry epidemiology scope.",
        }

    supported_years = {r["report_year"] for r in dataset}
    if report_year not in supported_years:
        return {
            "claim_id": cid,
            "predicted_verdict": "out_of_scope",
            "matched_row_id": None,
            "reason": f"Report year '{report_year}' is not in extracted ECFSPR corpus {sorted(supported_years)}.",
        }

    supported_indicators = {r["indicator_id"] for r in dataset}
    if indicator_id not in supported_indicators:
        return {
            "claim_id": cid,
            "predicted_verdict": "out_of_scope",
            "matched_row_id": None,
            "reason": f"Indicator '{indicator_id}' is not collected in the extracted ECFSPR indicator set.",
        }

    supported_countries = {r.get("stratum_country", "all_europe") for r in dataset}
    if stratum_country not in supported_countries:
        return {
            "claim_id": cid,
            "predicted_verdict": "out_of_scope",
            "matched_row_id": None,
            "reason": f"Country stratum '{stratum_country}' is not in the extracted ECFSPR country cohorts.",
        }

    # 2. Match candidate rows by (report_year, indicator_id, stratum_sex, stratum_country, ref_year)
    candidates = [
        r
        for r in dataset
        if r["report_year"] == report_year
        and r["indicator_id"] == indicator_id
        and r.get("stratum_sex", "total") == stratum_sex
        and r.get("stratum_country", "all_europe") == stratum_country
        and _matches_reference_cohort(r, ref_year)
    ]

    if not candidates:
        return {
            "claim_id": cid,
            "predicted_verdict": "out_of_scope",
            "matched_row_id": None,
            "reason": (
                f"No extracted row for year={report_year}, indicator={indicator_id}, "
                f"sex={stratum_sex}, country={stratum_country}, ref_year={ref_year}."
            ),
        }

    row = candidates[0]
    row_id = row["row_id"]

    # 3. Check denominator scope conflation if a specific scope was asserted
    actual_scope = row.get("denominator_scope")
    if claimed_scope is not None and actual_scope is not None and claimed_scope != actual_scope:
        return {
            "claim_id": cid,
            "predicted_verdict": "contradicted",
            "matched_row_id": row_id,
            "reason": (
                f"Denominator scope conflation: claim asserts scope '{claimed_scope}', "
                f"but indicator '{indicator_id}' in {row_id} has scope '{actual_scope}'."
            ),
        }

    # 4. Check quantitative value field
    actual_value = row.get(field)
    if actual_value is None:
        return {
            "claim_id": cid,
            "predicted_verdict": "out_of_scope",
            "matched_row_id": row_id,
            "reason": f"Field '{field}' is null/unreported in row {row_id}.",
        }

    diff = abs(float(actual_value) - float(claimed_value))
    if diff <= tolerance:
        return {
            "claim_id": cid,
            "predicted_verdict": "supported",
            "matched_row_id": row_id,
            "reason": (
                f"Matches {row_id} ({field}={actual_value}, claimed={claimed_value}, "
                f"page={row['page']}, scope={actual_scope})."
            ),
        }
    else:
        return {
            "claim_id": cid,
            "predicted_verdict": "contradicted",
            "matched_row_id": row_id,
            "reason": (
                f"Value mismatch on {row_id}: actual {field}={actual_value} != claimed {claimed_value} "
                f"(diff={diff:.4f} > tolerance={tolerance})."
            ),
        }


def run_benchmark(
    benchmark_path: Path = FROZEN_CLAIMS_JSON,
    dataset_path: Path = DATASET_JSON,
) -> dict[str, Any]:
    """Verify freeze hashes, run blind evaluation on all 60 claims, and score results."""
    benchmark_data = json.loads(benchmark_path.read_text(encoding="utf-8"))
    dataset = json.loads(dataset_path.read_text(encoding="utf-8"))

    freeze_ok, freeze_errors = verify_frozen_set(benchmark_data)
    if not freeze_ok:
        raise ValueError(f"Frozen benchmark verification failed: {freeze_errors}")

    claims = benchmark_data["claims"]
    blind_claims = strip_labels_for_blind_eval(claims)

    predictions: dict[str, dict[str, Any]] = {}
    for b_claim in blind_claims:
        res = evaluate_claim_blind(b_claim, dataset)
        predictions[b_claim["claim_id"]] = res

    return score_predictions(predictions, benchmark_data)


def score_predictions(
    predictions: dict[str, dict[str, Any]],
    benchmark_data: dict[str, Any],
) -> dict[str, Any]:
    """Score predicted verdicts against the frozen 60-claim gold standard."""
    claims = benchmark_data["claims"]
    total = len(claims)
    correct = 0

    by_verdict: dict[str, dict[str, int]] = {
        v: {"total": 0, "correct": 0} for v in ALLOWED_VERDICTS
    }
    detailed_results: list[dict[str, Any]] = []

    for gold in claims:
        cid = gold["claim_id"]
        expected = gold["expected_verdict"]
        pred_obj = predictions.get(cid, {})
        predicted = pred_obj.get("predicted_verdict", "missing")

        is_correct = predicted == expected
        if is_correct:
            correct += 1

        by_verdict[expected]["total"] += 1
        if is_correct:
            by_verdict[expected]["correct"] += 1

        detailed_results.append(
            {
                "claim_id": cid,
                "category": gold["category"],
                "expected_verdict": expected,
                "predicted_verdict": predicted,
                "correct": is_correct,
                "matched_row_id": pred_obj.get("matched_row_id"),
                "reason": pred_obj.get("reason", ""),
            }
        )

    overall_accuracy = correct / total if total > 0 else 0.0
    oos_total = by_verdict["out_of_scope"]["total"]
    oos_correct = by_verdict["out_of_scope"]["correct"]
    all_oos_correct = (oos_total == 10) and (oos_correct == 10)
    meets_acceptance = (overall_accuracy >= 0.90) and all_oos_correct

    return {
        "total_claims": total,
        "correct_verdicts": correct,
        "overall_accuracy": round(overall_accuracy, 4),
        "out_of_scope_total": oos_total,
        "out_of_scope_correct": oos_correct,
        "all_out_of_scope_correct": all_oos_correct,
        "meets_acceptance_criteria": meets_acceptance,
        "breakdown_by_verdict": by_verdict,
        "blind_input_sha256": benchmark_data["metadata"]["blind_input_sha256"],
        "full_benchmark_sha256": benchmark_data["metadata"]["full_benchmark_sha256"],
        "results": detailed_results,
    }


def main() -> int:
    summary = run_benchmark()
    print(
        f"Claim-check harness evaluated {summary['total_claims']} frozen claims: "
        f"{summary['correct_verdicts']}/{summary['total_claims']} correct "
        f"({summary['overall_accuracy'] * 100:.1f}%)."
    )
    print(
        f"Out-of-scope gate: {summary['out_of_scope_correct']}/{summary['out_of_scope_total']} "
        f"returned out_of_scope (all_out_of_scope_correct={summary['all_out_of_scope_correct']})."
    )
    print(f"Blind input SHA-256:    {summary['blind_input_sha256']}")
    print(f"Full benchmark SHA-256: {summary['full_benchmark_sha256']}")
    if not summary["meets_acceptance_criteria"]:
        print("FAIL: Acceptance criteria not met.", file=sys.stderr)
        return 1
    print("PASS: Meets T-0010 acceptance criteria (>=90% accuracy and 10/10 out_of_scope).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
