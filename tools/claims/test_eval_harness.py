"""Unit tests for the frozen 60-claim evaluation harness (T-0010)."""
from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from tools.claims import check_claims
from tools.claims.eval_harness import (
    DEFAULT_CLAIMS_PATH,
    DEFAULT_LOCK_PATH,
    DEFAULT_MAP_PATH,
    DEFAULT_SNAPSHOT_PATH,
    check_anti_hallucination,
    check_claim_grounding,
    detect_out_of_scope,
    run_evaluation,
    verify_lock,
)


class EvalHarnessTests(unittest.TestCase):
    """Standard-library unittest suite for the frozen 60-claim evaluation harness."""

    def test_lock_verification_success(self) -> None:
        """Verify that the frozen evaluation set matches its lock file."""
        lock_data = verify_lock(DEFAULT_CLAIMS_PATH, DEFAULT_LOCK_PATH)
        self.assertEqual(lock_data["status"], "frozen")
        self.assertEqual(lock_data["total_claims"], 60)
        self.assertEqual(
            lock_data["breakdown"],
            {
                "supported": 25,
                "unsupported": 25,
                "out_of_scope": 10,
            },
        )

    def test_lock_tamper_detection(self) -> None:
        """Verify that tampering with the claims file raises RuntimeError."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            tampered_claims = Path(tmp_dir) / "tampered_claims.json"
            tampered_claims.write_text('[{"id": "EVAL-01", "claim": "tampered"}]', encoding="utf-8")

            with self.assertRaisesRegex(RuntimeError, "Frozen eval set integrity lock violation"):
                verify_lock(tampered_claims, DEFAULT_LOCK_PATH)

    def test_full_harness_acceptance_criteria(self) -> None:
        """Verify that the harness achieves >=90% correct verdicts and 100% out-of-scope accuracy."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_file = Path(tmp_dir) / "test_results.json"
            report = run_evaluation(
                claims_path=DEFAULT_CLAIMS_PATH,
                lock_path=DEFAULT_LOCK_PATH,
                map_path=DEFAULT_MAP_PATH,
                snapshot_path=DEFAULT_SNAPSHOT_PATH,
                output_path=out_file,
            )

        self.assertTrue(report["lock_verified"])
        self.assertEqual(report["total_claims"], 60)
        self.assertEqual(report["total_correct"], 60)
        self.assertEqual(report["accuracy"], 1.0)
        self.assertEqual(report["false_accepts"], 0)
        self.assertTrue(report["acceptance_criteria_passed"])

        breakdown = report["breakdown"]
        self.assertEqual(breakdown["supported"]["correct"], 25)
        self.assertEqual(breakdown["supported"]["total"], 25)
        self.assertEqual(breakdown["unsupported"]["correct"], 25)
        self.assertEqual(breakdown["unsupported"]["total"], 25)
        self.assertEqual(breakdown["out_of_scope"]["correct"], 10)
        self.assertEqual(breakdown["out_of_scope"]["total"], 10)
        self.assertEqual(breakdown["out_of_scope"]["rate"], 1.0)

    def test_boundary_detection_patterns(self) -> None:
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
            is_oos, _ = detect_out_of_scope({"claim": text})
            self.assertEqual(
                is_oos,
                expected_oos,
                f"Failed on text: {text} (got {is_oos}, expected {expected_oos})",
            )

    def test_anti_hallucination_guardrails(self) -> None:
        """Verify AH-1 and AH-4 deterministic rejection."""
        valid_ncts = {"NCT00909532", "NCT02858843"}

        claim_fab = {"claim": "NCT99999999 showed efficacy in CF.", "quote": ""}
        ok, reason = check_anti_hallucination(claim_fab, valid_ncts)
        self.assertFalse(ok)
        self.assertIn("AH-1 violation", reason)

        claim_cure = {"claim": "This drug is proven to cure cystic fibrosis.", "quote": ""}
        ok, reason = check_anti_hallucination(claim_cure, valid_ncts)
        self.assertFalse(ok)
        self.assertIn("AH-4 violation", reason)

        claim_valid = {"claim": "NCT00909532 was completed.", "quote": ""}
        ok, reason = check_anti_hallucination(claim_valid, valid_ncts)
        self.assertTrue(ok)

    def test_clean_claim_strips_metadata(self) -> None:
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

        ev_paths = check_claims.parse_evidence_args([], str(DEFAULT_MAP_PATH))
        ev_texts, ev_errs = check_claims.read_evidence(ev_paths)

        ok, failures = check_claim_grounding(claim_with_meta, ev_texts, ev_errs)
        self.assertTrue(ok)
        self.assertEqual(len(failures), 0)


if __name__ == "__main__":
    unittest.main()
