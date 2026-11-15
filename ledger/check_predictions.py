#!/usr/bin/env python3
"""Validator for ledger/predictions.jsonl, the stated prediction ledger of Sensei Ewok's CF Lab.

Enforces the Discovery Protocol (docs/lab-focus-2026-11-14.md) and T-0910 schema:
1. Each line is a valid JSON object representing one stated prediction.
2. Required fields and types:
   - id: unique identifier string matching ^PRED-\\d{4}$
   - date: ISO date string YYYY-MM-DD
   - prediction: one sentence string stating the prediction
   - refutation_observation: string describing the observation that would refute it
   - prior_art_dois: list of DOI strings (e.g. 10.xxxx/...)
   - standard_method: string naming the standard method it must match or beat
   - measured_public_data: string, either 'none' or an admitted catalog source id from sources/catalog.yaml
   - outcome: string, strictly one of: 'supported', 'refuted', 'could_not_test', 'open'
   - checked_by_outside_reader: boolean (True or False)
3. Invariants:
   - No duplicate prediction IDs.
   - At most one prediction may have outcome == 'open' at any time.
   - If outcome == 'supported' and not checked_by_outside_reader, warning/error that outside reader
     check is required before publication as a result.
   - No unknown top-level fields.

Exit codes:
  0: all predictions valid
  1: validation errors found
  2: usage / missing file error
"""
from __future__ import annotations

import argparse
import datetime
import json
import re
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

REPO_ROOT = Path(__file__).resolve().parent.parent

ALLOWED_OUTCOMES: Set[str] = {"supported", "refuted", "could_not_test", "open"}
REQUIRED_FIELDS: Tuple[str, ...] = (
    "id",
    "date",
    "prediction",
    "refutation_observation",
    "prior_art_dois",
    "standard_method",
    "measured_public_data",
    "outcome",
    "checked_by_outside_reader",
)
OPTIONAL_FIELDS: Set[str] = {"notes", "data_request_plan", "catalog_sources_searched"}
ALL_FIELDS: Set[str] = set(REQUIRED_FIELDS) | OPTIONAL_FIELDS

ID_PATTERN = re.compile(r"^PRED-\d{4}$")
DATE_PATTERN = re.compile(r"^\d{4}-\d{2}-\d{2}$")
DOI_PATTERN = re.compile(r"^10\.\d{4,9}/[-._;()/:A-Za-z0-9]+$")


def load_catalog_source_ids(catalog_path: Optional[Path] = None) -> Set[str]:
    """Load valid source IDs from sources/catalog.yaml."""
    c_path = catalog_path or (REPO_ROOT / "sources" / "catalog.yaml")
    if not c_path.is_file():
        return set()

    source_ids: Set[str] = set()
    try:
        import yaml  # type: ignore

        data = yaml.safe_load(c_path.read_text(encoding="utf-8"))
        if isinstance(data, dict):
            for s in data.get("sources", []):
                if isinstance(s, dict) and "id" in s:
                    source_ids.add(str(s["id"]))
    except Exception:
        # Fallback to line parser if PyYAML is unavailable
        id_re = re.compile(r"^\s*-\s*id:\s*([^\s#]+)")
        for line in c_path.read_text(encoding="utf-8").splitlines():
            m = id_re.match(line)
            if m:
                source_ids.add(m.group(1).strip().strip("'\""))

    return source_ids


def validate_prediction_record(
    record: Any,
    line_num: int,
    seen_ids: Set[str],
    catalog_source_ids: Set[str],
) -> List[str]:
    """Validate a single prediction dictionary against schema and semantic rules."""
    errors: List[str] = []
    if not isinstance(record, dict):
        return [f"Line {line_num}: Expected a JSON object (dict), got {type(record).__name__}"]

    pred_id = str(record.get("id", f"line-{line_num}"))
    prefix = f"Line {line_num} [{pred_id}]"

    # Check required fields
    for field in REQUIRED_FIELDS:
        if field not in record:
            errors.append(f"{prefix}: Missing required field '{field}'")

    # Check for unknown fields
    unknown = set(record.keys()) - ALL_FIELDS
    if unknown:
        errors.append(f"{prefix}: Unknown fields detected: {sorted(unknown)}")

    # 1. ID validation
    if "id" in record:
        raw_id = record["id"]
        if not isinstance(raw_id, str) or not ID_PATTERN.match(raw_id):
            errors.append(f"{prefix}: 'id' must match pattern ^PRED-\\d{{4}}$, got {raw_id!r}")
        elif raw_id in seen_ids:
            errors.append(f"{prefix}: Duplicate prediction id '{raw_id}'")
        else:
            seen_ids.add(raw_id)

    # 2. Date validation
    if "date" in record:
        raw_date = record["date"]
        if not isinstance(raw_date, str) or not DATE_PATTERN.match(raw_date):
            errors.append(f"{prefix}: 'date' must be ISO format YYYY-MM-DD, got {raw_date!r}")
        else:
            try:
                datetime.date.fromisoformat(raw_date)
            except ValueError:
                errors.append(f"{prefix}: 'date' {raw_date!r} is not a valid calendar date")

    # 3. Prediction statement
    if "prediction" in record:
        pred_text = record["prediction"]
        if not isinstance(pred_text, str) or not pred_text.strip():
            errors.append(f"{prefix}: 'prediction' must be a non-empty string")

    # 4. Refutation observation
    if "refutation_observation" in record:
        refute_text = record["refutation_observation"]
        if not isinstance(refute_text, str) or not refute_text.strip():
            errors.append(f"{prefix}: 'refutation_observation' must be a non-empty string")

    # 5. Prior art DOIs
    if "prior_art_dois" in record:
        dois = record["prior_art_dois"]
        if not isinstance(dois, list):
            errors.append(f"{prefix}: 'prior_art_dois' must be a list of strings, got {type(dois).__name__}")
        else:
            for i, doi in enumerate(dois):
                if not isinstance(doi, str) or not DOI_PATTERN.match(doi.strip()):
                    errors.append(f"{prefix}: prior_art_dois[{i}] {doi!r} is not a valid DOI format")

    # 6. Standard method
    if "standard_method" in record:
        sm = record["standard_method"]
        if not isinstance(sm, str) or not sm.strip():
            errors.append(f"{prefix}: 'standard_method' must be a non-empty string naming the baseline method")

    # 7. Measured public data
    if "measured_public_data" in record:
        mpd = record["measured_public_data"]
        if not isinstance(mpd, str) or not mpd.strip():
            errors.append(f"{prefix}: 'measured_public_data' must be 'none' or a catalog source ID string")
        elif mpd != "none" and catalog_source_ids and mpd not in catalog_source_ids:
            errors.append(
                f"{prefix}: 'measured_public_data' {mpd!r} is not recognized in sources/catalog.yaml "
                "nor equal to 'none'"
            )

    # 8. Outcome
    if "outcome" in record:
        outcome = record["outcome"]
        if outcome not in ALLOWED_OUTCOMES:
            errors.append(
                f"{prefix}: 'outcome' must be one of {sorted(ALLOWED_OUTCOMES)}, got {outcome!r}"
            )

    # 9. Checked by outside reader
    if "checked_by_outside_reader" in record:
        cor = record["checked_by_outside_reader"]
        if not isinstance(cor, bool):
            errors.append(f"{prefix}: 'checked_by_outside_reader' must be a boolean (true/false), got {cor!r}")
        elif record.get("outcome") == "supported" and not cor:
            errors.append(
                f"{prefix}: Discovery protocol invariant: outcome is 'supported' but checked_by_outside_reader "
                "is false. A claim may not be considered supported without independent reader review."
            )

    return errors


def validate_predictions(
    jsonl_path: Path,
    catalog_path: Optional[Path] = None,
) -> Tuple[List[Dict[str, Any]], List[str]]:
    """Validate predictions file and global invariants."""
    if not jsonl_path.is_file():
        return [], [f"Predictions file not found: {jsonl_path}"]

    catalog_ids = load_catalog_source_ids(catalog_path)
    seen_ids: Set[str] = set()
    all_errors: List[str] = []
    records: List[Dict[str, Any]] = []

    content = jsonl_path.read_text(encoding="utf-8")
    lines = content.splitlines()

    open_count = 0

    for idx, line in enumerate(lines, start=1):
        stripped = line.strip()
        if not stripped:
            continue

        try:
            record = json.loads(stripped)
        except json.JSONDecodeError as exc:
            all_errors.append(f"Line {idx}: Invalid JSON syntax ({exc})")
            continue

        errs = validate_prediction_record(record, idx, seen_ids, catalog_ids)
        all_errors.extend(errs)

        if isinstance(record, dict):
            records.append(record)
            if record.get("outcome") == "open":
                open_count += 1

    # Invariant: Discovery protocol allows at most ONE open prediction at a time
    if open_count > 1:
        all_errors.append(
            f"Discovery protocol invariant violated: found {open_count} predictions with outcome 'open'. "
            "Keep exactly one open prediction; state a new one only after the open one is closed."
        )

    return records, all_errors


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate ledger/predictions.jsonl")
    parser.add_argument(
        "--file",
        "-f",
        type=Path,
        default=REPO_ROOT / "ledger" / "predictions.jsonl",
        help="Path to predictions.jsonl file",
    )
    parser.add_argument(
        "--catalog",
        "-c",
        type=Path,
        default=REPO_ROOT / "sources" / "catalog.yaml",
        help="Path to sources/catalog.yaml file",
    )
    parser.add_argument("--json", action="store_true", help="Output audit report as JSON")
    args = parser.parse_args()

    records, errors = validate_predictions(args.file, args.catalog)

    report = {
        "status": "PASS" if not errors else "FAIL",
        "file": str(args.file),
        "total_predictions": len(records),
        "open_predictions": sum(1 for r in records if r.get("outcome") == "open"),
        "error_count": len(errors),
        "errors": errors,
    }

    if args.json:
        print(json.dumps(report, indent=2))
    else:
        print(
            f"Prediction Ledger Audit: {report['status']} "
            f"(total={report['total_predictions']}, open={report['open_predictions']}, errors={report['error_count']})"
        )
        for err in errors:
            print(f"  ERROR: {err}", file=sys.stderr)

    return 0 if report["status"] == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
