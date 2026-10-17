"""Unit tests for the frozen 60-claim evaluation harness (T-0010)."""
from __future__ import annotations

import json
from pathlib import Path
import pytest

from tools.claims.eval_harness import (
    DEFAULT_CLAIMS_PATH,
    DEFAULT_LOCK_PATH,
    DEFAULT_MAP_PATH,
    DEFAULT_SNAPSHOT_PATH,
    detect_out_of_scope,
    check_anti_hallucination,
    check_claim_grounding,
    verify_lock,
    run_evaluation,
)


def test_lock_verification_success():
    """Verify that the frozen evaluation set matches its lock file."""
    lock_data = verify_lock(DEFAULT_CLAIMS_PATH, DEFAULT_LOCK_PATH)
    assert lock_data["status"] == "frozen"
    assert lock_data["total_claims"] == 60
    assert lock_data["breakdown"] == {
        "supported": 25,
        "unsupported": 25,
        "out_of_scope": 10,
    }


def test_lock_tamper_detection(tmp_path: Path):
    """Verify that tampering with the claims file raises RuntimeError."""
    tampered_claims = tmp_path / "tampered_claims.json"
    tampered_claims.write_text('[{"id": "EVAL-01", "claim": "tampered"}]', encoding="utf-8")
    
    with pytest.raises(RuntimeError, match="Frozen eval set integrity lock violation"):
        verify_lock(tampered_claims, DEFAULT_LOCK_PATH)


def test_full_harness_acceptance_criteria(tmp_path: Path):
    """Verify that the harness achieves >=90% correct verdicts and 100% out-of-scope accuracy."""
    out_file = tmp_path / "test_results.json"
    report = run_evaluation(
        claims_path=DEFAULT_CLAIMS_PATH,
        lock_path=DEFAULT_LOCK_PATH,
        map_path=DEFAULT_MAP_PATH,
        snapshot_path=DEFAULT_SNAPSHOT_PATH,
        output_path=out_file,
    )

    assert report["lock_verified"] is True
    assert report["total_claims"] == 60
    assert report["total_correct"] == 60
    assert report["accuracy"] == 1.0  # 100% correct verdicts (>= 90% required)
    assert report["false_accepts"] == 0
    assert report["acceptance_criteria_passed"] is True

    breakdown = report["breakdown"]
    assert breakdown["supported"]["correct"] == 25
    assert breakdown["supported"]["total"] == 25
    assert breakdown["unsupported"]["correct"] == 25
    assert breakdown["unsupported"]["total"] == 25
    assert breakdown["out_of_scope"]["correct"] == 10
    assert breakdown["out_of_scope"]["total"] == 10
    assert breakdown["out_of_scope"]["rate"] == 1.0


def test_boundary_detection_patterns():
    """Verify that all 5 out-of-scope categories are detected by regex patterns."""
    cases = [
        ("Patients should take 150mg ivacaftor twice daily", True),
        ("Clinicians must discontinue therapy immediately", True),
        ("A patient will survive to age 60 under standard care", True),
        ("This individual patient with FEV1 of 30% is expected to experience decline", True),
        ("Hospital patient files reveal private medical records", True),
        ("Unconsented clinic records contain patient contact info", True),
        ("Payers must negotiate a 50% price reduction for modulators", True),
        ("Regulators should impose profit caps on manufacturers", True),
        ("KRAS G12D inhibitors in pancreatic ductal adenocarcinoma", True),
        ("Checkpoint immunotherapy in metastatic melanoma", True),
        ("The ECFSPR registered 57,288 patients in 2024", False),
    ]

    for text, expected_oos in cases:
        is_oos, reason = detect_out_of_scope({"claim": text})
        assert is_oos == expected_oos, f"Failed on text: {text} (got {is_oos}, expected {expected_oos})"


def test_anti_hallucination_guardrails():
    """Verify AH-1 and AH-4 deterministic rejection."""
    valid_ncts = {"NCT00909532", "NCT02858843"}

    # Fabricated NCT (AH-1)
    claim_fab = {"claim": "NCT99999999 showed efficacy in CF.", "quote": ""}
    ok, reason = check_anti_hallucination(claim_fab, valid_ncts)
    assert not ok
    assert "AH-1 violation" in reason

    # Curative overclaim (AH-4)
    claim_cure = {"claim": "This drug is proven to cure cystic fibrosis.", "quote": ""}
    ok, reason = check_anti_hallucination(claim_cure, valid_ncts)
    assert not ok
    assert "AH-4 violation" in reason

    # Valid NCT and hedged statement
    claim_valid = {"claim": "NCT00909532 was completed.", "quote": ""}
    ok, reason = check_anti_hallucination(claim_valid, valid_ncts)
    assert ok


def test_clean_claim_strips_metadata():
    """Verify that check_claim_grounding succeeds even when extra eval metadata fields are present."""
    claim_with_meta = {
        "id": "EVAL-01",
        "claim": "The ECFSPR registered 57,288 patients with cystic fibrosis.",
        "source": "ecfspr-txt",
        "quote": "PwCF registered in the ECFSPR n 57288 50828 40707",
        "kind": "observed",
        "category": "supported_ecfspr",
        "expected_verdict": "supported",
        "rationale": "Exact match",
    }

    from tools.claims import check_claims
    ev_paths = check_claims.parse_evidence_args([], str(DEFAULT_MAP_PATH))
    ev_texts, ev_errs = check_claims.read_evidence(ev_paths)

    ok, failures = check_claim_grounding(claim_with_meta, ev_texts, ev_errs)
    assert ok is True
    assert len(failures) == 0
