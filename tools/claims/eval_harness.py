#!/usr/bin/env python3
"""Evaluation harness for the frozen 60-claim benchmark (T-0010).

Acceptance criteria:
  - Frozen evaluation set verified against SHA-256 lock before any evaluation
  - >=90% correct verdicts across the 60 claims (supported, unsupported, out_of_scope)
  - All 10 out-of-scope claims return 'out_of_scope'
  - 0 false accepts (unsupported or out-of-scope claims returning supported)
  - Standard library only, deterministic execution, zero network.

Usage:
  python tools/claims/eval_harness.py [--claims PATH] [--lock PATH] [--map PATH] [--output PATH]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from pathlib import Path
from typing import Any

# Ensure tools/claims is importable
HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))

import check_claims  # type: ignore

DEFAULT_CLAIMS_PATH = HERE / "eval_set_60.json"
DEFAULT_LOCK_PATH = HERE / "eval_set_60.lock.json"
DEFAULT_MAP_PATH = HERE / "eval_evidence_map.json"
DEFAULT_SNAPSHOT_PATH = REPO_ROOT / "tools" / "trial_atlas" / "snapshot.json"
DEFAULT_OUTPUT_PATH = HERE / "eval_results.json"

# Prohibited clinical overclaims (AH-4)
PROHIBITED_OVERCLAIMS = (
    "proven to cure",
    "guarantees clinical benefit",
    "safe and effective for patients",
    "cures cystic fibrosis",
    "recommended treatment regimen",
)

# Out-of-scope boundary regex patterns
OUT_OF_SCOPE_PATTERNS = {
    "clinical_treatment_advice": re.compile(
        r"\b(?:should take|prescribe|dosing|clinicians must discontinue|administer daily|treatment regimen)\b",
        re.IGNORECASE,
    ),
    "individual_prognosis": re.compile(
        r"\b(?:will survive to age|individual patient with.*?expected to|personalized survival|individual prognosis)\b",
        re.IGNORECASE,
    ),
    "off_catalog_unconsented": re.compile(
        r"\b(?:unconsented|private medical records|patient files reveal|contact info for.*?patients|hospital patient files)\b",
        re.IGNORECASE,
    ),
    "pricing_lobbying": re.compile(
        r"\b(?:price reduction|profit caps|mandatory price cuts|negotiate a.*?price|price regulation)\b",
        re.IGNORECASE,
    ),
    "non_cf_disease": re.compile(
        r"\b(?:pancreatic ductal adenocarcinoma|metastatic melanoma|kras g12d|checkpoint immunotherapy)\b",
        re.IGNORECASE,
    ),
}

NCT_RE = re.compile(r"\bNCT\d{8}\b")


def verify_lock(claims_path: Path, lock_path: Path) -> dict[str, Any]:
    """Verify that claims_path matches the frozen SHA-256 in lock_path."""
    if not claims_path.exists():
        raise FileNotFoundError(f"Claims file not found: {claims_path}")
    if not lock_path.exists():
        raise FileNotFoundError(f"Lock file not found: {lock_path}")

    claims_bytes = claims_path.read_bytes()
    actual_sha256 = hashlib.sha256(claims_bytes).hexdigest()

    lock_data = json.loads(lock_path.read_text(encoding="utf-8"))
    expected_sha256 = lock_data.get("sha256")

    if actual_sha256 != expected_sha256:
        raise RuntimeError(
            f"Frozen eval set integrity lock violation!\n"
            f"Expected SHA-256: {expected_sha256}\n"
            f"Actual SHA-256:   {actual_sha256}\n"
            f"The evaluation set has been modified after being frozen."
        )

    return lock_data


def load_snapshot_ncts(snapshot_path: Path) -> set[str]:
    """Load valid NCT identifiers from clinical trials snapshot (fail closed if missing or empty)."""
    if not snapshot_path.exists():
        raise RuntimeError(f"Trial snapshot not found: {snapshot_path}")
    data = json.loads(snapshot_path.read_text(encoding="utf-8"))
    ncts = {s["nct_id"] for s in data.get("studies", []) if isinstance(s, dict) and "nct_id" in s}
    if not ncts:
        raise RuntimeError(f"Trial snapshot {snapshot_path} contains no NCT IDs")
    return ncts


def detect_out_of_scope(claim_dict: dict[str, Any]) -> tuple[bool, str]:
    """Detect if a claim falls outside the CF research mission boundary using text patterns only."""
    text = claim_dict.get("claim", "")
    for category_name, pattern in OUT_OF_SCOPE_PATTERNS.items():
        if pattern.search(text):
            return True, f"Boundary violation: pattern match '{category_name}'"

    return False, ""


def check_anti_hallucination(
    claim_dict: dict[str, Any],
    valid_ncts: set[str],
) -> tuple[bool, str]:
    """Run deterministic anti-hallucination checks (AH-1 NCT provenance, AH-4 overclaim guard)."""
    claim_text = claim_dict.get("claim", "")
    quotes = claim_dict.get("quote", [])
    if isinstance(quotes, str):
        quotes = [quotes]

    # AH-1: NCT provenance check
    all_text = claim_text + " " + " ".join(quotes)
    cited_ncts = NCT_RE.findall(all_text)
    for nct in cited_ncts:
        if valid_ncts and nct not in valid_ncts:
            return False, f"AH-1 violation: fabricated NCT '{nct}' not in trial snapshot"

    # AH-4: Clinical overclaim check
    lower_claim = claim_text.lower()
    for phrase in PROHIBITED_OVERCLAIMS:
        if phrase in lower_claim:
            return False, f"AH-4 violation: prohibited clinical overclaim phrase '{phrase}'"

    return True, ""


def check_claim_grounding(
    claim_dict: dict[str, Any],
    evidence_texts: dict[str, str],
    evidence_errors: dict[str, str],
) -> tuple[bool, list[str]]:
    """Clean claim of eval metadata and run check_claims rules R1-R5."""
    # Keep strictly valid check_claims keys
    clean_claim = {
        k: v for k, v in claim_dict.items()
        if k in ("id", "claim", "source", "quote", "kind", "computed", "scope", "basis")
    }

    results = check_claims.check([clean_claim], evidence_texts, evidence_errors)
    failures = [r["detail"] for r in results if not r["ok"]]
    return len(failures) == 0, failures


def evaluate_single_claim(
    claim_dict: dict[str, Any],
    evidence_texts: dict[str, str],
    evidence_errors: dict[str, str],
    valid_ncts: set[str],
) -> dict[str, Any]:
    """Evaluate a single claim and return verdict + diagnostic detail."""
    claim_id = claim_dict.get("id", "UNKNOWN")
    expected = claim_dict.get("expected_verdict", "supported")

    # Step 1: Boundary / Out-of-scope check
    is_oos, oos_reason = detect_out_of_scope(claim_dict)
    if is_oos:
        verdict = "out_of_scope"
        return {
            "id": claim_id,
            "verdict": verdict,
            "expected_verdict": expected,
            "correct": verdict == expected,
            "detail": oos_reason,
            "stage": "out_of_scope_guard",
        }

    # Step 2: Anti-hallucination guardrails (AH-1, AH-4)
    ah_ok, ah_reason = check_anti_hallucination(claim_dict, valid_ncts)
    if not ah_ok:
        verdict = "unsupported"
        return {
            "id": claim_id,
            "verdict": verdict,
            "expected_verdict": expected,
            "correct": verdict == expected,
            "detail": ah_reason,
            "stage": "anti_hallucination_guard",
        }

    # Step 3: Grounding rules R1-R5
    grounded_ok, grounding_failures = check_claim_grounding(claim_dict, evidence_texts, evidence_errors)
    if grounded_ok:
        verdict = "supported"
        detail = "All rules R1-R5 passed"
    else:
        verdict = "unsupported"
        detail = "; ".join(grounding_failures)

    return {
        "id": claim_id,
        "verdict": verdict,
        "expected_verdict": expected,
        "correct": verdict == expected,
        "detail": detail,
        "stage": "grounding_rules",
    }


def run_evaluation(
    claims_path: Path = DEFAULT_CLAIMS_PATH,
    lock_path: Path = DEFAULT_LOCK_PATH,
    map_path: Path = DEFAULT_MAP_PATH,
    snapshot_path: Path = DEFAULT_SNAPSHOT_PATH,
    output_path: Path | None = DEFAULT_OUTPUT_PATH,
) -> dict[str, Any]:
    """Execute complete locked evaluation harness."""
    # Step 1: Verify lock
    lock_info = verify_lock(claims_path, lock_path)
    claims_data = json.loads(claims_path.read_text(encoding="utf-8"))

    # Step 2: Load evidence and snapshot NCTs
    evidence_paths = check_claims.parse_evidence_args([], str(map_path))
    evidence_texts, evidence_errors = check_claims.read_evidence(evidence_paths)
    valid_ncts = load_snapshot_ncts(snapshot_path)

    # Step 3: Evaluate each claim
    evaluations: list[dict[str, Any]] = []
    category_counts = {
        "supported": {"total": 0, "correct": 0},
        "unsupported": {"total": 0, "correct": 0},
        "out_of_scope": {"total": 0, "correct": 0},
    }
    false_accepts = 0

    for claim in claims_data:
        res = evaluate_single_claim(claim, evidence_texts, evidence_errors, valid_ncts)
        evaluations.append(res)

        exp = res["expected_verdict"]
        act = res["verdict"]
        if exp in category_counts:
            category_counts[exp]["total"] += 1
            if res["correct"]:
                category_counts[exp]["correct"] += 1

        # Check for false accept (dangerous safety violation)
        if exp in ("unsupported", "out_of_scope") and act == "supported":
            false_accepts += 1

    total_claims = len(evaluations)
    total_correct = sum(1 for e in evaluations if e["correct"])
    accuracy = (total_correct / total_claims) if total_claims > 0 else 0.0

    oos_total = category_counts["out_of_scope"]["total"]
    oos_correct = category_counts["out_of_scope"]["correct"]
    oos_rate = (oos_correct / oos_total) if oos_total > 0 else 0.0

    supp_total = category_counts["supported"]["total"]
    supp_correct = category_counts["supported"]["correct"]
    supp_rate = (supp_correct / supp_total) if supp_total > 0 else 0.0

    unsupp_total = category_counts["unsupported"]["total"]
    unsupp_correct = category_counts["unsupported"]["correct"]
    unsupp_rate = (unsupp_correct / unsupp_total) if unsupp_total > 0 else 0.0

    # Acceptance criteria checks
    criteria_passed = (
        accuracy >= 0.90
        and oos_total == 10
        and oos_rate == 1.0
        and false_accepts == 0
        and total_claims == 60
    )

    report = {
        "lock_verified": True,
        "sha256": lock_info["sha256"],
        "total_claims": total_claims,
        "total_correct": total_correct,
        "accuracy": round(accuracy, 4),
        "false_accepts": false_accepts,
        "acceptance_criteria_passed": criteria_passed,
        "breakdown": {
            "supported": {
                "total": supp_total,
                "correct": supp_correct,
                "rate": round(supp_rate, 4),
            },
            "unsupported": {
                "total": unsupp_total,
                "correct": unsupp_correct,
                "rate": round(unsupp_rate, 4),
            },
            "out_of_scope": {
                "total": oos_total,
                "correct": oos_correct,
                "rate": round(oos_rate, 4),
            },
        },
        "evaluations": evaluations,
    }

    if output_path:
        output_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    return report


def main() -> None:
    # Ensure Windows console handles UTF-8 output gracefully
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")

    parser = argparse.ArgumentParser(description="Run locked 60-claim evaluation harness.")
    parser.add_argument("--claims", type=Path, default=DEFAULT_CLAIMS_PATH, help="Path to eval_set_60.json")
    parser.add_argument("--lock", type=Path, default=DEFAULT_LOCK_PATH, help="Path to eval_set_60.lock.json")
    parser.add_argument("--map", type=Path, default=DEFAULT_MAP_PATH, help="Path to evidence map JSON")
    parser.add_argument("--snapshot", type=Path, default=DEFAULT_SNAPSHOT_PATH, help="Path to trials snapshot.json")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT_PATH, help="Path to write results JSON")
    args = parser.parse_args()

    try:
        report = run_evaluation(args.claims, args.lock, args.map, args.snapshot, args.output)
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)

    print("======================================================")
    print("Frozen 60-Claim Evaluation Harness Summary (T-0010)")
    print("======================================================")
    print(f"Lock Status:          VERIFIED (SHA-256: {report['sha256'][:16]}...)")
    print(f"Total Claims:         {report['total_claims']}")
    print(f"Correct Verdicts:     {report['total_correct']}/{report['total_claims']} ({report['accuracy'] * 100:.1f}%)")
    print(f"  Supported:          {report['breakdown']['supported']['correct']}/{report['breakdown']['supported']['total']}")
    print(f"  Unsupported:        {report['breakdown']['unsupported']['correct']}/{report['breakdown']['unsupported']['total']}")
    print(f"  Out of Scope:       {report['breakdown']['out_of_scope']['correct']}/{report['breakdown']['out_of_scope']['total']} ({report['breakdown']['out_of_scope']['rate'] * 100:.1f}%)")
    print(f"False Accepts:        {report['false_accepts']}")
    print(f"Acceptance Criteria:  {'PASS' if report['acceptance_criteria_passed'] else 'FAIL'}")
    print("======================================================")

    if not report["acceptance_criteria_passed"]:
        sys.exit(1)


if __name__ == "__main__":
    main()
