"""Unit tests for the ECFSPR claim-check harness and frozen 60-claim benchmark (Task T-0010)."""
from __future__ import annotations

import copy
import json
import sys
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from claim_harness import (
    DATASET_JSON,
    FROZEN_CLAIMS_JSON,
    evaluate_claim_blind,
    run_benchmark,
    strip_labels_for_blind_eval,
    verify_frozen_set,
)


class ClaimHarnessTest(unittest.TestCase):
    def setUp(self):
        self.assertTrue(FROZEN_CLAIMS_JSON.exists(), f"Missing {FROZEN_CLAIMS_JSON}")
        self.assertTrue(DATASET_JSON.exists(), f"Missing {DATASET_JSON}")
        self.benchmark = json.loads(FROZEN_CLAIMS_JSON.read_text(encoding="utf-8"))
        self.dataset = json.loads(DATASET_JSON.read_text(encoding="utf-8"))

    def test_frozen_set_cryptographic_integrity(self):
        ok, errors = verify_frozen_set(self.benchmark)
        self.assertTrue(ok, f"Frozen set verification failed: {errors}")
        self.assertEqual(len(self.benchmark["claims"]), 60)
        self.assertEqual(
            self.benchmark["metadata"]["blind_input_sha256"],
            "d306e5575cd6151c1ad431aa7ce83b30e118c3a4e745dd97b48923742dcd4a21",
        )
        self.assertEqual(
            self.benchmark["metadata"]["full_benchmark_sha256"],
            "12df3288f5c735a57554e89f2f85eb5d21bfb8df9eed5031d892888382389d5e",
        )

    def test_tampering_with_frozen_set_is_detected(self):
        tampered = copy.deepcopy(self.benchmark)
        tampered["claims"][0]["query"]["claimed_value"] = 99999
        ok, errors = verify_frozen_set(tampered)
        self.assertFalse(ok)
        self.assertTrue(any("SHA-256 mismatch" in e for e in errors))

    def test_blind_strip_removes_gold_labels(self):
        blind = strip_labels_for_blind_eval(self.benchmark["claims"])
        self.assertEqual(len(blind), 60)
        for item in blind:
            self.assertNotIn("expected_verdict", item)
            self.assertNotIn("category", item)
            self.assertNotIn("grounding_rationale", item)
            self.assertEqual(set(item.keys()), {"claim_id", "claim_text", "query"})

    def test_benchmark_exceeds_ninety_percent_and_ten_of_ten_out_of_scope(self):
        summary = run_benchmark(FROZEN_CLAIMS_JSON, DATASET_JSON)
        self.assertEqual(summary["total_claims"], 60)
        self.assertGreaterEqual(summary["overall_accuracy"], 0.90)
        self.assertEqual(summary["overall_accuracy"], 1.0)
        self.assertEqual(summary["out_of_scope_total"], 10)
        self.assertEqual(summary["out_of_scope_correct"], 10)
        self.assertTrue(summary["all_out_of_scope_correct"])
        self.assertTrue(summary["meets_acceptance_criteria"])

    def test_denominator_conflation_claims_return_contradicted(self):
        denom_claims = [
            c for c in self.benchmark["claims"] if c["category"] == "contradicted_denominator_conflation"
        ]
        self.assertEqual(len(denom_claims), 5)
        for c in denom_claims:
            blind_input = strip_labels_for_blind_eval([c])[0]
            res = evaluate_claim_blind(blind_input, self.dataset)
            self.assertEqual(res["predicted_verdict"], "contradicted", f"Failed on {c['claim_id']}")


if __name__ == "__main__":
    unittest.main()
