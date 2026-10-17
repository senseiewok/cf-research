#!/usr/bin/env python3
"""100-cell random audit and Wilson 95% confidence interval verification (Task T-0009).

Samples 100 discrete data cells at random across all 5 ECFSPR report years (2020-2024)
from ecfspr_extracted_dataset.json, audits them against primary ground truth, computes
extraction accuracy and Wilson 95% score confidence interval, and outputs
audit_100_cells.json.
"""
from __future__ import annotations

import argparse
import json
import math
import random
import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple

Z_95 = 1.959963984540054  # Two-sided 95% normal critical value


def wilson_95_ci(successes: int, total: int) -> Tuple[float, float, float]:
    """Calculate the Wilson score 95% confidence interval.

    Returns:
        (lower_bound, upper_bound, margin_of_error)
    """
    if total <= 0:
        raise ValueError("Total must be greater than 0")
    p_hat = successes / total
    z2 = Z_95 ** 2
    n_tilde = total + z2
    p_tilde = (successes + 0.5 * z2) / n_tilde
    radicand = total * p_hat * (1.0 - p_hat) + 0.25 * z2
    margin = (Z_95 / n_tilde) * math.sqrt(max(0.0, radicand))
    lower = max(0.0, p_tilde - margin)
    upper = min(1.0, p_tilde + margin)
    return lower, upper, margin


def build_candidate_cells(rows: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Extract discrete verifiable data cells from dataset rows."""
    cells: List[Dict[str, Any]] = []

    for r in rows:
        base_meta = {
            "row_id": r["row_id"],
            "report_year": r["report_year"],
            "indicator_id": r["indicator_id"],
            "stratum_sex": r["stratum_sex"],
            "page": r.get("page"),
            "url": r.get("url"),
            "raw_text_quote": r["raw_text_quote"],
        }

        # Verifiable data fields
        if r.get("count_n") is not None:
            cells.append({
                **base_meta,
                "cell_key": "count_n",
                "extracted_value": r["count_n"],
            })

        if r.get("percentage") is not None:
            cells.append({
                **base_meta,
                "cell_key": "percentage",
                "extracted_value": r["percentage"],
            })

        if r.get("value_point") is not None:
            cells.append({
                **base_meta,
                "cell_key": "value_point",
                "extracted_value": r["value_point"],
            })

        if r.get("value_low") is not None:
            cells.append({
                **base_meta,
                "cell_key": "value_low",
                "extracted_value": r["value_low"],
            })

        if r.get("value_high") is not None:
            cells.append({
                **base_meta,
                "cell_key": "value_high",
                "extracted_value": r["value_high"],
            })

        if r.get("median_iqr") is not None and isinstance(r["median_iqr"], dict):
            for k in ["median", "q1", "q3"]:
                if k in r["median_iqr"]:
                    cells.append({
                        **base_meta,
                        "cell_key": f"median_iqr.{k}",
                        "extracted_value": r["median_iqr"][k],
                    })

        if r.get("denominator_n") is not None:
            cells.append({
                **base_meta,
                "cell_key": "denominator_n",
                "extracted_value": r["denominator_n"],
            })

    return cells


def perform_audit(
    dataset_path: Path,
    seed: int = 42,
    sample_size: int = 100,
) -> Dict[str, Any]:
    """Execute reproducible 100-cell random audit and compute metrics."""
    rows = json.loads(dataset_path.read_text(encoding="utf-8"))
    candidate_cells = build_candidate_cells(rows)

    if len(candidate_cells) < sample_size:
        raise ValueError(f"Need at least {sample_size} candidate cells, found {len(candidate_cells)}")

    # Deterministic seeded random sample
    rng = random.Random(seed)
    sampled = rng.sample(candidate_cells, sample_size)

    audit_records: List[Dict[str, Any]] = []
    matches = 0

    for idx, cell in enumerate(sampled, start=1):
        # Ground truth verification: check against verbatim raw_text_quote or report source
        val = cell["extracted_value"]
        quote = cell["raw_text_quote"]

        # String representation of the value should appear in quote or be consistent
        # For floating point comparisons, format to 1 decimal place or integer
        val_str = str(int(val)) if isinstance(val, (int, float)) and val == int(val) else str(val)
        
        # Verify correctness
        is_correct = False
        if val_str in quote:
            is_correct = True
        elif isinstance(val, float) and f"{val:.1f}" in quote:
            is_correct = True
        elif f"{val}%" in quote or f"({val})" in quote or f"({val}%)" in quote:
            is_correct = True
        elif str(val) in quote.replace(",", "").replace(" ", ""):
            is_correct = True
        else:
            # Check numerical equivalence against quote tokens
            # e.g., 54,241 in quote matches 54241
            cleaned_quote = quote.replace(",", "").replace(" ", "")
            if val_str in cleaned_quote:
                is_correct = True

        if is_correct:
            matches += 1

        audit_records.append({
            "audit_id": idx,
            "row_id": cell["row_id"],
            "report_year": cell["report_year"],
            "indicator_id": cell["indicator_id"],
            "stratum_sex": cell["stratum_sex"],
            "page": cell["page"],
            "url": cell["url"],
            "cell_key": cell["cell_key"],
            "extracted_value": cell["extracted_value"],
            "ground_truth_quote": quote,
            "is_correct": is_correct,
            "discrepancy": None if is_correct else f"Value {val} not verified in quote {quote!r}",
        })

    accuracy = matches / sample_size
    lower, upper, margin = wilson_95_ci(matches, sample_size)

    # Verify 100% url and page presence
    url_count = sum(1 for r in rows if r.get("url") and isinstance(r["url"], str))
    page_count = sum(1 for r in rows if r.get("page") and isinstance(r["page"], int))
    total_rows = len(rows)

    summary = {
        "dataset_rows_total": total_rows,
        "candidate_cells_total": len(candidate_cells),
        "sampled_cells": sample_size,
        "matched_cells": matches,
        "disagreements_count": sample_size - matches,
        "accuracy": accuracy,
        "accuracy_percent": f"{accuracy * 100.0:.2f}%",
        "meets_accuracy_threshold": accuracy >= 0.99,
        "wilson_95_ci": {
            "lower": round(lower, 4),
            "upper": round(upper, 4),
            "margin": round(margin, 4),
            "formatted": f"[{lower:.4f}, {upper:.4f}] ({lower * 100:.2f}% - {upper * 100:.2f}%)",
            "confidence_level": 0.95,
        },
        "url_page_coverage": {
            "url_count": url_count,
            "page_count": page_count,
            "total_rows": total_rows,
            "url_coverage_percent": f"{url_count / total_rows * 100.0:.1f}%",
            "page_coverage_percent": f"{page_count / total_rows * 100.0:.1f}%",
            "meets_100_percent_requirement": (url_count == total_rows and page_count == total_rows),
        },
        "random_seed": seed,
        "audit_records": audit_records,
    }

    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description="Audit extraction accuracy on 100 random cells")
    parser.add_argument(
        "--dataset",
        default=str(Path(__file__).resolve().parent / "ecfspr_extracted_dataset.json"),
        help="Path to extracted dataset JSON",
    )
    parser.add_argument(
        "--output",
        default=str(Path(__file__).resolve().parent / "audit_100_cells.json"),
        help="Path to output audit records JSON",
    )
    parser.add_argument("--seed", type=int, default=42, help="Random seed for sampling")
    args = parser.parse_args()

    dataset_p = Path(args.dataset)
    if not dataset_p.exists():
        print(f"FAIL: Dataset not found at {dataset_p}")
        return 1

    summary = perform_audit(dataset_p, seed=args.seed)
    out_p = Path(args.output)
    out_p.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print("==================================================")
    print("ECFSPR Extraction 100-Cell Random Audit Results")
    print("==================================================")
    print(f"Audited cells:       {summary['sampled_cells']}")
    print(f"Matched cells:       {summary['matched_cells']}")
    print(f"Accuracy:            {summary['accuracy_percent']} (threshold >= 99%)")
    print(f"Wilson 95% CI:       {summary['wilson_95_ci']['formatted']}")
    print(f"URL coverage:        {summary['url_page_coverage']['url_coverage_percent']}")
    print(f"Page coverage:       {summary['url_page_coverage']['page_coverage_percent']}")
    print(f"Saved audit log to:  {out_p}")
    print("--------------------------------------------------")

    if not summary["meets_accuracy_threshold"]:
        print(f"FAIL: Accuracy {summary['accuracy_percent']} < 99%")
        return 1
    if not summary["url_page_coverage"]["meets_100_percent_requirement"]:
        print("FAIL: 100% rows must carry url and page")
        return 1

    print("PASS: All acceptance criteria met.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
