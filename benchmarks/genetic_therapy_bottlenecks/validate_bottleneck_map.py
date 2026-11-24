#!/usr/bin/env python3
"""Deterministic validator for the CF genetic-therapy bottleneck evidence map (Council Item #5).

Enforces strict epistemic and provenance invariants over `evidence_map.json`:
1. Rejects any ranking, leaderboard, or composite efficacy scoring keys anywhere in the JSON tree.
2. Enforces closed status and termination-reason enums so `TERMINATED` never quietly becomes failed efficacy
   and clinical registry listings never masquerade as measured lung-function benefit.
3. Enforces typed evidence anchors for every populated cell (`verbatim_passage` vs `snapshot_field_pointer`),
   verifying repo-relative file LF-normalized SHA-256 digests and exact substrings on disk.
4. Verifies every study's `source_id` (DOI/NCT) exists with `verdict == "PASS"` in `sources/cite-check-ledger.jsonl`.
5. Forbids clinical/delivery `measured_in_system` cells on preclinical studies and forbids `S5` clinical benefit
   `measured_in_system` cells on unposted ClinicalTrials.gov registry snapshots.

Pure Python standard library only; zero network calls.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys
from typing import Any, Dict, List, Optional, Set, Tuple

HERE = Path(__file__).resolve().parent
DEFAULT_MAP_PATH = HERE / "evidence_map.json"
DEFAULT_REPO_ROOT = HERE.parents[1]

SCHEMA_VERSION = "gt-bottleneck-map-v1"

STAGE_KEYS: Tuple[str, ...] = (
    "S1_cell_correction_or_expression",
    "S2_restored_channel_function",
    "S3_in_vivo_airway_delivery",
    "S4_durability_or_redosing_persistence",
    "S5_clinical_benefit_or_lung_function",
)

ALLOWED_STAGE_STATUSES: Tuple[str, ...] = (
    "measured_in_system",
    "registered_endpoint_not_reported",
    "study_terminated_unreported",
    "untested_in_study",
    "not_stated",
)

POPULATED_EVIDENCE_STATUSES: Tuple[str, ...] = (
    "measured_in_system",
    "registered_endpoint_not_reported",
    "study_terminated_unreported",
)

ALLOWED_TERMINATION_REASONS: Tuple[str, ...] = (
    "safety_signal",
    "funding_or_sponsor",
    "enrollment",
    "manufacturing",
    "unspecified_in_snapshot",
)

ALLOWED_STUDY_DESIGNS: Tuple[str, ...] = (
    "preclinical_in_vitro_or_organoid",
    "clinical_trial_registry_snapshot",
    "clinical_randomized_phase2b",
)

ALLOWED_MEASUREMENT_SCALES: Tuple[str, ...] = (
    "molecular_and_cellular",
    "clinical_safety_registry",
    "clinical_spirometry_rct",
)

ALLOWED_SNAPSHOT_FIELDS: Tuple[str, ...] = (
    "overall_status",
    "phases",
    "interventions",
    "primary_outcomes",
)

FORBIDDEN_LEADERBOARD_KEYS: frozenset[str] = frozenset(
    {
        "score",
        "rank",
        "ranking",
        "efficacy_score",
        "leaderboard",
        "overall_rating",
        "best_therapy",
    }
)


def _lf_sha256(path: Path) -> str:
    """Compute SHA-256 over LF-normalized file bytes so Windows CRLF and Linux LF checkouts agree."""
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def _resolve_repo_path(repo_root: Path, rel_path: str) -> Path:
    """Resolve a repo-relative path and reject absolute paths or parent escapes."""
    candidate = Path(rel_path)
    if candidate.is_absolute():
        raise ValueError(f"Path must be repo-relative, got {rel_path!r}")
    root_resolved = repo_root.resolve()
    resolved = (root_resolved / candidate).resolve()
    resolved.relative_to(root_resolved)
    return resolved


def _scan_forbidden_keys(obj: Any, path_prefix: str = "$") -> List[str]:
    violations: List[str] = []
    if isinstance(obj, dict):
        for k, v in obj.items():
            kl = str(k).strip().lower()
            cur = f"{path_prefix}.{k}"
            if kl in FORBIDDEN_LEADERBOARD_KEYS:
                violations.append(f"Forbidden ranking/leaderboard key at {cur!r}")
            violations.extend(_scan_forbidden_keys(v, cur))
    elif isinstance(obj, list):
        for idx, item in enumerate(obj):
            violations.extend(_scan_forbidden_keys(item, f"{path_prefix}[{idx}]"))
    return violations


def _load_verified_ledger_ids(ledger_path: Path) -> Set[str]:
    if not ledger_path.is_file():
        raise FileNotFoundError(f"Missing cite-check ledger: {ledger_path}")
    verified: Set[str] = set()
    for line in ledger_path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        rec = json.loads(line)
        if str(rec.get("verdict") or "").upper() == "PASS" and rec.get("id"):
            verified.add(str(rec["id"]).strip())
    return verified


def _load_snapshot_by_nct(snapshot_path: Path) -> Dict[str, Dict[str, Any]]:
    if not snapshot_path.is_file():
        raise FileNotFoundError(f"Missing ClinicalTrials.gov snapshot: {snapshot_path}")
    data = json.loads(snapshot_path.read_text(encoding="utf-8"))
    studies = data.get("studies") or []
    return {str(s["nct_id"]).strip(): s for s in studies if s.get("nct_id")}


def validate_bottleneck_map_data(
    payload: Any,
    *,
    repo_root: Path = DEFAULT_REPO_ROOT,
) -> Dict[str, Any]:
    """Validate bottleneck evidence map structure, provenance hashes, snapshot pointers, and anti-overclaim rules."""
    errors: List[str] = []
    stage_status_counts: Dict[str, Dict[str, int]] = {
        s: {st: 0 for st in ALLOWED_STAGE_STATUSES} for s in STAGE_KEYS
    }

    if not isinstance(payload, dict):
        return {
            "status": "FAIL",
            "errors": ["Top-level payload must be a JSON object."],
            "studies_validated": 0,
            "stage_status_counts": stage_status_counts,
        }

    errors.extend(_scan_forbidden_keys(payload))

    if payload.get("schema_version") != SCHEMA_VERSION:
        errors.append(f"schema_version must be {SCHEMA_VERSION!r}, got {payload.get('schema_version')!r}")
    if payload.get("non_comparable_across_scales") is not True:
        errors.append("non_comparable_across_scales must be true.")
    if payload.get("no_ranking_or_leaderboard") is not True:
        errors.append("no_ranking_or_leaderboard must be true.")
    if not isinstance(payload.get("research_boundary"), str) or not payload["research_boundary"].strip():
        errors.append("Missing non-empty 'research_boundary' statement.")

    ledger_rel = str(payload.get("ledger_file") or "sources/cite-check-ledger.jsonl")
    snapshot_rel = str(payload.get("snapshot_file") or "tools/trial_atlas/snapshot.json")
    ledger_loaded = False
    try:
        verified_ids = _load_verified_ledger_ids(_resolve_repo_path(repo_root, ledger_rel))
        ledger_loaded = True
    except Exception as exc:
        errors.append(f"Invalid ledger_file {ledger_rel!r}: {exc}")
        verified_ids = set()

    try:
        snapshot_by_nct = _load_snapshot_by_nct(_resolve_repo_path(repo_root, snapshot_rel))
    except Exception as exc:
        errors.append(f"Invalid snapshot_file {snapshot_rel!r}: {exc}")
        snapshot_by_nct = {}

    studies = payload.get("studies")
    if not isinstance(studies, list) or not studies:
        errors.append("'studies' must be a non-empty list.")
        return {
            "status": "FAIL",
            "errors": errors,
            "studies_validated": 0,
            "stage_status_counts": stage_status_counts,
        }

    seen_study_ids: Set[str] = set()
    file_cache: Dict[Path, Tuple[str, str]] = {}
    per_study_summary: List[Dict[str, Any]] = []

    for idx, study in enumerate(studies):
        if not isinstance(study, dict):
            errors.append(f"Study index {idx} must be an object.")
            continue

        sid = str(study.get("study_id") or "").strip()
        if not sid:
            errors.append(f"Study index {idx} missing non-empty 'study_id'.")
            sid = f"index-{idx}"
        elif sid in seen_study_ids:
            errors.append(f"Duplicate study_id: {sid!r}")
        seen_study_ids.add(sid)

        source_id = str(study.get("source_id") or "").strip()
        if not source_id:
            errors.append(f"Study {sid!r} missing 'source_id'.")
        elif ledger_loaded and source_id not in verified_ids:
            errors.append(f"Study {sid!r} source_id {source_id!r} not found with verdict=PASS in {ledger_rel}")

        design = study.get("study_design")
        if design not in ALLOWED_STUDY_DESIGNS:
            errors.append(f"Study {sid!r} invalid study_design={design!r}")
        scale = study.get("measurement_scale")
        if scale not in ALLOWED_MEASUREMENT_SCALES:
            errors.append(f"Study {sid!r} invalid measurement_scale={scale!r}")

        for req_header in ("experimental_system", "controls", "follow_up"):
            if not isinstance(study.get(req_header), str) or not study[req_header].strip():
                errors.append(f"Study {sid!r} missing header field {req_header!r}.")

        stages = study.get("stages")
        if not isinstance(stages, dict):
            errors.append(f"Study {sid!r} missing 'stages' dict.")
            continue

        missing_stages = [k for k in STAGE_KEYS if k not in stages]
        extra_stages = [k for k in stages if k not in STAGE_KEYS]
        if missing_stages:
            errors.append(f"Study {sid!r} missing stage keys: {missing_stages}")
        if extra_stages:
            errors.append(f"Study {sid!r} unknown stage keys: {extra_stages}")

        cell_statuses: Dict[str, str] = {}
        for skey in STAGE_KEYS:
            cell = stages.get(skey)
            if not isinstance(cell, dict):
                errors.append(f"Study {sid!r} stage {skey!r} must be a dict.")
                continue

            status = cell.get("status")
            if status not in ALLOWED_STAGE_STATUSES:
                errors.append(f"Study {sid!r} stage {skey!r} invalid status={status!r}")
                continue
            cell_statuses[skey] = status
            stage_status_counts[skey][status] += 1

            # Hard guardrail 1: Registry snapshot cannot claim S5 clinical benefit measured.
            if (
                design == "clinical_trial_registry_snapshot"
                and skey == "S5_clinical_benefit_or_lung_function"
                and status == "measured_in_system"
            ):
                errors.append(
                    f"Registry hallucination guard: {sid!r} has study_design={design!r} "
                    "and cannot mark S5 clinical benefit as 'measured_in_system'."
                )

            # Hard guardrail 2: Preclinical study cannot claim in vivo clinical/delivery measured stages.
            if (
                design == "preclinical_in_vitro_or_organoid"
                and skey in (
                    "S3_in_vivo_airway_delivery",
                    "S4_durability_or_redosing_persistence",
                    "S5_clinical_benefit_or_lung_function",
                )
                and status == "measured_in_system"
            ):
                errors.append(
                    f"Preclinical overclaim guard: {sid!r} has study_design={design!r} "
                    f"and cannot mark {skey!r} as 'measured_in_system'."
                )

            # Hard guardrail 3: Terminated studies require explicit termination_reason_class.
            if status == "study_terminated_unreported":
                trc = cell.get("termination_reason_class")
                if trc not in ALLOWED_TERMINATION_REASONS:
                    errors.append(
                        f"Study {sid!r} stage {skey!r} has status='study_terminated_unreported' "
                        f"without valid termination_reason_class in {ALLOWED_TERMINATION_REASONS}; got {trc!r}"
                    )

            # Hard guardrail 4: Typed evidence anchors for populated statuses.
            if status in POPULATED_EVIDENCE_STATUSES:
                etype = cell.get("evidence_type")
                if etype == "verbatim_passage":
                    for bad_key in ("snapshot_nct_id", "snapshot_field", "expected_substring"):
                        if bad_key in cell:
                            errors.append(
                                f"Study {sid!r} stage {skey!r} mixes snapshot field {bad_key!r} "
                                "with evidence_type='verbatim_passage'."
                            )
                    rel_src = str(cell.get("source_file") or "").strip()
                    expected_sha = str(cell.get("source_file_sha256") or "").strip()
                    locator = str(cell.get("locator") or "").strip()
                    passage = str(cell.get("passage_text") or "")
                    if not rel_src or not expected_sha or not locator or not passage.strip():
                        errors.append(
                            f"Study {sid!r} stage {skey!r} verbatim_passage missing source_file/source_file_sha256/locator/passage_text."
                        )
                        continue
                    try:
                        src_path = _resolve_repo_path(repo_root, rel_src)
                    except ValueError as exc:
                        errors.append(f"Study {sid!r} stage {skey!r} escapes repo root ({rel_src!r}): {exc}")
                        continue
                    if not src_path.is_file():
                        errors.append(f"Study {sid!r} stage {skey!r} missing source_file: {rel_src!r}")
                        continue
                    if src_path not in file_cache:
                        file_cache[src_path] = (_lf_sha256(src_path), src_path.read_text(encoding="utf-8"))
                    actual_sha, actual_text = file_cache[src_path]
                    if actual_sha != expected_sha:
                        errors.append(
                            f"Study {sid!r} stage {skey!r} source_file_sha256 drift on {rel_src!r}: "
                            f"expected {expected_sha}, got {actual_sha}"
                        )
                    if passage not in actual_text:
                        errors.append(
                            f"Study {sid!r} stage {skey!r} passage_text is not a verbatim substring of {rel_src!r}"
                        )
                elif etype == "snapshot_field_pointer":
                    for bad_key in ("source_file", "source_file_sha256", "passage_text"):
                        if bad_key in cell:
                            errors.append(
                                f"Study {sid!r} stage {skey!r} mixes verbatim field {bad_key!r} "
                                "with evidence_type='snapshot_field_pointer'."
                            )
                    nct_id = str(cell.get("snapshot_nct_id") or "").strip()
                    sfield = str(cell.get("snapshot_field") or "").strip()
                    expected_sub = str(cell.get("expected_substring") or "")
                    if not nct_id or sfield not in ALLOWED_SNAPSHOT_FIELDS or not expected_sub.strip():
                        errors.append(
                            f"Study {sid!r} stage {skey!r} invalid snapshot_field_pointer fields."
                        )
                        continue
                    if nct_id != source_id:
                        errors.append(
                            f"Study {sid!r} stage {skey!r} snapshot_nct_id={nct_id!r} != study source_id={source_id!r}"
                        )
                    rec = snapshot_by_nct.get(nct_id)
                    if rec is None:
                        errors.append(f"Study {sid!r} stage {skey!r} NCT {nct_id!r} absent from {snapshot_rel}")
                        continue
                    field_blob = json.dumps(rec.get(sfield), ensure_ascii=False)
                    if expected_sub not in field_blob:
                        errors.append(
                            f"Study {sid!r} stage {skey!r} expected_substring {expected_sub!r} not found in "
                            f"{nct_id}.{sfield}"
                        )
                else:
                    errors.append(
                        f"Study {sid!r} stage {skey!r} with status={status!r} requires valid evidence_type "
                        "('verbatim_passage' or 'snapshot_field_pointer')."
                    )

        per_study_summary.append(
            {
                "study_id": sid,
                "source_id": source_id,
                "study_design": design,
                "measurement_scale": scale,
                "stage_statuses": cell_statuses,
            }
        )

    return {
        "status": "PASS" if not errors else "FAIL",
        "schema_version": payload.get("schema_version"),
        "studies_validated": len(studies),
        "errors": errors,
        "stage_status_counts": stage_status_counts,
        "studies": per_study_summary,
    }


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--map", type=Path, default=DEFAULT_MAP_PATH, help="Path to evidence_map.json")
    parser.add_argument("--repo-root", type=Path, default=DEFAULT_REPO_ROOT, help="Repository root path")
    parser.add_argument("--check", action="store_true", help="Exit non-zero if validation fails")
    parser.add_argument("--json", action="store_true", help="Print JSON output (allow_nan=False)")
    args = parser.parse_args(argv)

    payload = json.loads(args.map.read_text(encoding="utf-8"))
    report = validate_bottleneck_map_data(payload, repo_root=args.repo_root)
    if args.json:
        print(json.dumps(report, indent=2, allow_nan=False))
    else:
        print(
            f"Genetic-Therapy Bottleneck Map Check: {report['status']} "
            f"(studies={report['studies_validated']}, errors={len(report['errors'])})"
        )
        for s in report.get("studies", []):
            summary = ", ".join(f"{k.split('_', 1)[0]}={v}" for k, v in s["stage_statuses"].items())
            print(f"  - {s['study_id']} [{s['source_id']}]: {summary}")
        for err in report["errors"]:
            print(f"  ERROR: {err}", file=sys.stderr)

    if args.check and report["status"] != "PASS":
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
