#!/usr/bin/env python3
"""Deterministic JSON Schema validator for cf-registry-extract row datasets.

Implements lightweight, standard-library-only JSON Schema validation against
tools/registry-extract/row_schema.json without external dependencies.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Any, Dict, List

DEFAULT_SCHEMA_PATH = Path(__file__).resolve().parent / "row_schema.json"


def _check_type(val: Any, expected_type: Any) -> bool:
    if isinstance(expected_type, list):
        return any(_check_type(val, t) for t in expected_type)
    if expected_type == "string":
        return isinstance(val, str)
    if expected_type == "integer":
        return isinstance(val, int) and not isinstance(val, bool)
    if expected_type == "number":
        return (isinstance(val, (int, float)) and not isinstance(val, bool))
    if expected_type == "object":
        return isinstance(val, dict)
    if expected_type == "null":
        return val is None
    if expected_type == "boolean":
        return isinstance(val, bool)
    if expected_type == "array":
        return isinstance(val, list)
    return True


def validate_row(row: Dict[str, Any], schema: Dict[str, Any], path_prefix: str = "row") -> List[str]:
    """Validate a single registry extract row dictionary against the schema."""
    errors: List[str] = []

    if not isinstance(row, dict):
        return [f"{path_prefix}: expected object, got {type(row).__name__}"]

    # 1. Required fields
    required = schema.get("required", [])
    for field in required:
        if field not in row:
            errors.append(f"{path_prefix}: missing required field '{field}'")

    # 2. Properties
    properties = schema.get("properties", {})
    additional_props = schema.get("additionalProperties", True)

    if not additional_props:
        for key in row:
            if key not in properties:
                errors.append(f"{path_prefix}: unexpected property '{key}' not permitted by schema")

    for key, spec in properties.items():
        if key not in row:
            continue
        val = row[key]

        # Type check
        expected_type = spec.get("type")
        if expected_type is not None and not _check_type(val, expected_type):
            errors.append(f"{path_prefix}.{key}: expected type {expected_type}, got {type(val).__name__} ({val!r})")
            continue

        if val is None:
            continue

        # Enum check
        if "enum" in spec and val not in spec["enum"]:
            errors.append(f"{path_prefix}.{key}: value {val!r} not in allowed enum {spec['enum']}")

        # Pattern check
        if "pattern" in spec and isinstance(val, str):
            if not re.search(spec["pattern"], val):
                errors.append(f"{path_prefix}.{key}: value {val!r} does not match pattern {spec['pattern']}")

        # String minLength
        if "minLength" in spec and isinstance(val, str):
            if len(val) < spec["minLength"]:
                errors.append(f"{path_prefix}.{key}: string length {len(val)} is less than minLength {spec['minLength']}")

        # Number/integer minimum / maximum
        if isinstance(val, (int, float)):
            if "minimum" in spec and val < spec["minimum"]:
                errors.append(f"{path_prefix}.{key}: value {val} is less than minimum {spec['minimum']}")
            if "maximum" in spec and val > spec["maximum"]:
                errors.append(f"{path_prefix}.{key}: value {val} is greater than maximum {spec['maximum']}")

        # Nested object properties (e.g. median_iqr)
        if isinstance(val, dict) and "properties" in spec:
            nested_errors = validate_row(val, spec, path_prefix=f"{path_prefix}.{key}")
            errors.extend(nested_errors)

    return errors


def validate_dataset(rows: List[Dict[str, Any]], schema: Dict[str, Any]) -> List[str]:
    """Validate a list of row dictionaries."""
    errors: List[str] = []
    if not isinstance(rows, list):
        return [f"Dataset root must be a list of rows, got {type(rows).__name__}"]

    row_ids = set()
    for idx, row in enumerate(rows):
        prefix = f"row[{idx}]"
        if isinstance(row, dict) and "row_id" in row:
            rid = row["row_id"]
            if rid in row_ids:
                errors.append(f"{prefix}: duplicate row_id {rid!r}")
            row_ids.add(rid)
            prefix = f"row[{idx}:{rid}]"

        row_errors = validate_row(row, schema, path_prefix=prefix)
        errors.extend(row_errors)

    return errors


def validate_file(file_path: Path, schema_path: Path = DEFAULT_SCHEMA_PATH) -> List[str]:
    """Load JSON file and validate against schema."""
    if not file_path.exists():
        return [f"File not found: {file_path}"]
    if not schema_path.exists():
        return [f"Schema not found: {schema_path}"]

    try:
        schema = json.loads(schema_path.read_text(encoding="utf-8"))
    except Exception as e:
        return [f"Failed to parse schema JSON from {schema_path}: {e}"]

    try:
        data = json.loads(file_path.read_text(encoding="utf-8"))
    except Exception as e:
        return [f"Failed to parse dataset JSON from {file_path}: {e}"]

    if isinstance(data, list):
        return validate_dataset(data, schema)
    elif isinstance(data, dict):
        return validate_row(data, schema)
    else:
        return [f"Unsupported root JSON type: {type(data).__name__}"]


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate registry extract dataset against JSON Schema.")
    parser.add_argument("dataset", nargs="?", default=None, help="Path to dataset JSON file to validate.")
    parser.add_argument("--schema", default=str(DEFAULT_SCHEMA_PATH), help="Path to row_schema.json.")
    args = parser.parse_args()

    schema_p = Path(args.schema)
    if not args.dataset:
        default_gold = Path(__file__).resolve().parent / "gold_sample_ecfspr_2022_p8.json"
        dataset_p = default_gold
    else:
        dataset_p = Path(args.dataset)

    print(f"Validating {dataset_p} against {schema_p}...")
    errors = validate_file(dataset_p, schema_p)
    if errors:
        print(f"FAIL: {len(errors)} validation errors encountered:")
        for err in errors:
            print(f"  - {err}")
        return 1

    print(f"PASS: {dataset_p} successfully validated against schema (0 errors).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
