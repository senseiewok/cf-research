#!/usr/bin/env python3
"""Unit tests for tools/benchmark/train_claim_support_critic.py."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from tools.benchmark.train_claim_support_critic import (
    FEATURE_NAMES,
    extract_feature_vector,
    main,
    predict_proba,
    train_and_export_critic,
)


class TestTrainClaimSupportCritic(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.artifact = train_and_export_critic()
        cls.weights = [cls.artifact["feature_weights"][k] for k in FEATURE_NAMES]
        cls.bias = cls.artifact["bias"]

    def test_cross_validation_and_dpo_preference_accuracy(self) -> None:
        cv = self.artifact["cross_validation_5fold"]
        dpo = self.artifact["dpo_preference_ranking"]
        self.assertGreaterEqual(cv["oof_f1_on_violations"], 0.98)
        self.assertGreaterEqual(cv["oof_accuracy"], 0.98)
        self.assertEqual(dpo["preference_accuracy"], 1.0)
        self.assertGreater(dpo["mean_probability_margin"], 0.70)

    def test_scores_grounded_vs_adversarial_claims(self) -> None:
        quote = (
            "Overall, 92.3 percent of the total CF population in the Registry was eligible "
            "by CFTR genotype for at least one CFTR modulator as of 2024."
        )
        scope = "CFFPR 2024 Annual Data Report page 8."
        x_good = extract_feature_vector(
            "In 2024, 92.3 percent of the CFFPR population was eligible by genotype for at least one modulator.",
            quote,
            scope,
            "observed",
        )
        x_bad_num = extract_feature_vector(
            "In 2024, 98.7 percent of the CFFPR population was eligible by genotype for at least one modulator.",
            quote,
            scope,
            "observed",
        )
        x_bad_safe = extract_feature_vector(
            "In 2024, 92.3 percent was eligible; therefore you are eligible and should take this recommended dose for you.",
            quote,
            scope,
            "observed",
        )

        p_good = predict_proba(x_good, self.weights, self.bias)
        p_bad_num = predict_proba(x_bad_num, self.weights, self.bias)
        p_bad_safe = predict_proba(x_bad_safe, self.weights, self.bias)

        self.assertLess(p_good, 0.25)
        self.assertGreater(p_bad_num, 0.80)
        self.assertGreater(p_bad_safe, 0.85)

    def test_cli_write_and_check(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out_file = Path(tmp) / "weights.json"
            rc = main(["--out", str(out_file)])
            self.assertEqual(rc, 0)
            self.assertTrue(out_file.exists())
            loaded = json.loads(out_file.read_text(encoding="utf-8"))
            self.assertEqual(loaded["model_id"], "cf-claim-support-critic-l2-v1")


if __name__ == "__main__":
    unittest.main()
