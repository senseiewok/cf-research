#!/usr/bin/env python3
"""Unit and planted-defect tests for tools/readiness/evaluate_readiness.py."""
from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from tools.readiness.evaluate_readiness import (
    BASELINE_GIT_COMMIT,
    MIN_KAPPA_N,
    MIN_N_PRECISION_95,
    SCHEMA_VERSION,
    compare_annotators,
    compute_nominal_cohens_kappa,
    evaluate_dataset,
    legacy_pre_pr388_extract_external_validation,
    min_n_for_wilson_lower_bound,
    normalize_metric_text,
    validate_annotation_payload,
    values_match,
    wilson_score_interval,
)
from tools.readiness.extractor import ReadinessRecord

READINESS_DIR = Path(__file__).resolve().parent
STRESS_FILE = READINESS_DIR / "stress_eval_synthetic.json"
BENCHMARK_FILE = READINESS_DIR / "cf_ml_benchmark.json"


class _StubExtractor:
    """Deterministic stub extractor for planted-defect harness tests."""

    def __init__(self, record: ReadinessRecord) -> None:
        self._record = record

    def extract(self, text: str, paper_id: str = "") -> ReadinessRecord:
        return self._record


def _make_item(
    paper_id: str,
    study_group_id: str,
    *,
    cohort_size: int | str = 120,
    split_type: str = "patient_level",
    external_validation: str = "not_stated",
    headline_metric: str = "AUC-ROC of 0.80",
    readiness_tier: str = "T2_INTERNALLY_VALIDATED",
    doi: str = "",
    pmid: str = "",
) -> dict:
    item = {
        "paper_id": paper_id,
        "study_group_id": study_group_id,
        "text": "Placeholder text.",
        "expected": {
            "cohort_size": cohort_size,
            "split_type": split_type,
            "external_validation": external_validation,
            "headline_metric": headline_metric,
            "readiness_tier": readiness_tier,
        },
    }
    if doi:
        item["doi"] = doi
    if pmid:
        item["pmid"] = pmid
    return item


class TestWilsonAndPower(unittest.TestCase):
    def test_wilson_known_answers(self) -> None:
        self.assertIsNone(wilson_score_interval(0, 0))
        ci_6 = wilson_score_interval(6, 6)
        self.assertIsNotNone(ci_6)
        self.assertAlmostEqual(ci_6[0], 0.6097, places=3)
        self.assertEqual(ci_6[1], 1.0)

        ci_72 = wilson_score_interval(72, 72)
        self.assertIsNotNone(ci_72)
        self.assertLess(ci_72[0], 0.95)

        ci_73 = wilson_score_interval(73, 73)
        self.assertIsNotNone(ci_73)
        self.assertGreaterEqual(ci_73[0], 0.95)

        ci_95 = wilson_score_interval(95, 100)
        self.assertIsNotNone(ci_95)
        self.assertAlmostEqual(ci_95[0], 0.8882, places=3)
        self.assertAlmostEqual(ci_95[1], 0.9785, places=3)

    def test_min_n_for_95_lower_bound(self) -> None:
        self.assertEqual(min_n_for_wilson_lower_bound(0.95), 73)
        self.assertEqual(MIN_N_PRECISION_95, 73)

    def test_invalid_wilson_inputs_raise(self) -> None:
        with self.assertRaises(ValueError):
            wilson_score_interval(7, 6)
        with self.assertRaises(ValueError):
            wilson_score_interval(-1, 5)


class TestMatchingPoliciesAndDefects(unittest.TestCase):
    def test_cohort_size_rejects_bool_and_str(self) -> None:
        self.assertTrue(values_match("cohort_size", 200, 200))
        self.assertFalse(values_match("cohort_size", True, 1))
        self.assertFalse(values_match("cohort_size", "200", 200))
        self.assertTrue(values_match("cohort_size", "not stated", "not_stated"))

    def test_headline_metric_normalization_and_no_partial_substring_cheat(self) -> None:
        self.assertEqual(normalize_metric_text("-1.40 ± 1.53 points"), "-1.40±1.53 points")
        self.assertTrue(values_match("headline_metric", "AUCs of 0.75–0.80", "AUCs of 0.75-0.80"))
        self.assertTrue(
            values_match(
                "headline_metric",
                "C-index of 0.81",
                "concordance index of 0.81",
                accepted_aliases=["C-index of 0.81"],
            )
        )
        # Planted defect: a truncated substring must NOT be rewarded as a match.
        self.assertFalse(values_match("headline_metric", "AUC", "AUC-ROC of 0.89"))
        self.assertFalse(values_match("headline_metric", "AUC-ROC of 0.89", "AUC"))

    def test_planted_ungrounded_quote_and_span_mismatch_detected(self) -> None:
        fabricated = ReadinessRecord(
            paper_id="planted-1",
            cohort_size=100,
            cohort_sentence="Fabricated sentence absent from source abstract.",
            split_type="patient_level",
            split_sentence="Patient-level splitting was applied.",
            external_validation="not_stated",
            external_validation_sentence="",
            headline_metric="AUC-ROC of 0.82",
            metric_sentence="The model achieved an AUC-ROC of 0.82.",
            readiness_tier="T2_INTERNALLY_VALIDATED",
            anti_overclaim_flags=[],
        )
        dataset = {
            "schema_version": SCHEMA_VERSION,
            "dataset_kind": "synthetic_defect_probe",
            "is_unseen_published_papers": False,
            "items": [
                {
                    "paper_id": "planted-1",
                    "study_group_id": "grp-planted-1",
                    "text": (
                        "We evaluated 100 patients with CF. "
                        "Patient-level splitting was applied. "
                        "The model achieved an AUC-ROC of 0.82."
                    ),
                    "expected": {
                        "cohort_size": 100,
                        "split_type": "patient_level",
                        "external_validation": "not_stated",
                        "headline_metric": "AUC-ROC of 0.82",
                        "readiness_tier": "T2_INTERNALLY_VALIDATED",
                    },
                    "expected_grounding_substrings": {
                        "split_type": "wrong-span-token"
                    },
                }
            ],
        }
        report = evaluate_dataset(dataset, extractor=_StubExtractor(fabricated))
        cs = report["factual_extraction"]["fields"]["cohort_size"]
        st = report["factual_extraction"]["fields"]["split_type"]
        self.assertEqual(cs["error_taxonomy"]["ungrounded_quote"], 1)
        self.assertEqual(cs["verbatim_grounding_rate"], 0.0)
        self.assertEqual(st["error_taxonomy"]["grounding_span_mismatch"], 1)


class TestSyntheticStressEvaluation(unittest.TestCase):
    def test_stress_suite_separates_factual_from_tier_and_beats_baseline(self) -> None:
        dataset = json.loads(STRESS_FILE.read_text(encoding="utf-8"))
        report = evaluate_dataset(dataset, compare_baseline=True)

        self.assertEqual(report["dataset_kind"], "synthetic_stress_fixture")
        self.assertFalse(report["is_unseen_published_papers"])
        self.assertEqual(report["items_evaluated"], 18)
        self.assertFalse(report["dev_leakage_check"]["overlap_detected"])
        self.assertTrue(report["factual_extraction"]["evaluated_independently_from_readiness_tier"])
        self.assertTrue(report["heuristic_readiness_tier"]["heuristic_derived_from_factual_fields"])

        fields = report["factual_extraction"]["fields"]
        for fname in ("cohort_size", "split_type", "external_validation", "headline_metric"):
            self.assertEqual(fields[fname]["precision_target_status"], "underpowered_n_lt_73")
            self.assertEqual(fields[fname]["positive_assertion_precision"], 1.0)
            self.assertEqual(fields[fname]["verbatim_grounding_rate"], 1.0)

        # Novel metric phrasing (`synth-recall-01`) must surface honestly as recall < 1.0.
        hm = fields["headline_metric"]
        self.assertLess(hm["recall_on_stated"], 1.0)
        self.assertEqual(hm["error_taxonomy"]["false_negative_miss"], 1)

        # Zero T3 overclaims under current hardened extractor vs multiple under pre-PR #388 baseline.
        tier = report["heuristic_readiness_tier"]
        self.assertEqual(tier["t3_overclaim_count"], 0)
        self.assertEqual(
            report["factual_extraction"]["overall_verbatim_grounding_rate"], 1.0
        )

        bc = report["baseline_comparison"]
        self.assertEqual(bc["baseline_git_commit"], BASELINE_GIT_COMMIT)
        self.assertGreater(bc["baseline_t3_overclaims"], 0)
        self.assertEqual(bc["current_t3_overclaims"], 0)
        self.assertGreater(bc["t3_overclaim_reduction"], 0)
        self.assertGreater(bc["positive_assertion_precision_delta"], 0.0)
        self.assertEqual(legacy_pre_pr388_extract_external_validation([])[0], "not_stated")


class TestDualAnnotatorAndLeakage(unittest.TestCase):
    def test_schema_validator_rejects_bad_payloads(self) -> None:
        errs = validate_annotation_payload({"schema_version": "wrong", "items": []})
        self.assertTrue(any("schema_version" in e for e in errs))
        self.assertTrue(any("annotator_id" in e for e in errs))

        bad_enum = {
            "schema_version": SCHEMA_VERSION,
            "dataset_kind": "unseen_review",
            "annotator_id": "rev-a",
            "items": [
                _make_item("u1", "grp-1", cohort_size=True, split_type="bad_split")
            ],
        }
        errs2 = validate_annotation_payload(bad_enum)
        self.assertTrue(any("cohort_size" in e for e in errs2))
        self.assertTrue(any("split_type" in e for e in errs2))

    def test_cohens_kappa_known_answers_and_guards(self) -> None:
        self.assertEqual(compute_nominal_cohens_kappa(["a", "b"], ["a", "b"])["kappa"], 1.0)
        zero_var = compute_nominal_cohens_kappa(["a", "a", "a"], ["a", "a", "a"])
        self.assertIsNone(zero_var["kappa"])
        self.assertEqual(zero_var["kappa_status"], "undefined_zero_variance")

        # N = 40: 15 (a,a), 15 (b,b), 5 (a,b), 5 (b,a) => p_o=0.75, p_e=0.50, kappa=0.50 exact.
        la = ["a"] * 15 + ["b"] * 15 + ["a"] * 5 + ["b"] * 5
        lb = ["a"] * 15 + ["b"] * 15 + ["b"] * 5 + ["a"] * 5
        k_40 = compute_nominal_cohens_kappa(la, lb)
        self.assertEqual(k_40["n"], 40)
        self.assertGreaterEqual(k_40["n"], MIN_KAPPA_N)
        self.assertEqual(k_40["kappa_status"], "computed")
        self.assertAlmostEqual(k_40["kappa"], 0.5, places=4)

        k_small = compute_nominal_cohens_kappa(["a", "b", "a", "b"], ["a", "b", "b", "b"])
        self.assertEqual(k_small["kappa_status"], f"underpowered_n_lt_{MIN_KAPPA_N}")

    def test_compare_annotators_symmetry_and_dev_leakage_detection(self) -> None:
        # Includes a leaked DOI matching Jiang et al. 2026 in cf_ml_benchmark.json via https://doi.org/ prefix.
        items_a = [
            _make_item("unseen-1", "grp-u1", doi="https://doi.org/10.1002/MCO2.70894", headline_metric="AUC-ROC of 0.85"),
            _make_item("unseen-2", "grp-u2", external_validation="absent_single_center_or_internal_only", readiness_tier="T1_RETROSPECTIVE_OR_UNVALIDATED"),
        ]
        items_b = [
            _make_item("unseen-1", "grp-u1", doi="https://doi.org/10.1002/mco2.70894", headline_metric="AUCs of 0.70"),
            _make_item("unseen-2", "grp-u2", external_validation="not_stated", readiness_tier="T2_INTERNALLY_VALIDATED"),
        ]
        pay_a = {
            "schema_version": SCHEMA_VERSION,
            "dataset_kind": "unseen_review",
            "annotator_id": "reviewer-a",
            "items": items_a,
        }
        pay_b = {
            "schema_version": SCHEMA_VERSION,
            "dataset_kind": "unseen_review",
            "annotator_id": "reviewer-b",
            "items": items_b,
        }
        res_ab = compare_annotators(pay_a, pay_b, dev_benchmark_path=BENCHMARK_FILE)
        res_ba = compare_annotators(pay_b, pay_a, dev_benchmark_path=BENCHMARK_FILE)

        self.assertTrue(res_ab["dev_leakage_annotator_a"]["overlap_detected"])
        self.assertIn("10.1002/mco2.70894", res_ab["dev_leakage_annotator_a"]["overlapping_dois"])
        self.assertEqual(len(res_ab["disagreements"]), 3)
        self.assertEqual(
            res_ab["field_agreement"]["cohort_size"]["cohens_kappa"]["kappa_status"],
            "not_applicable_non_nominal",
        )
        self.assertEqual(
            res_ab["field_agreement"]["headline_metric"]["cohens_kappa"]["kappa_status"],
            "not_applicable_non_nominal",
        )
        self.assertEqual(
            res_ab["field_agreement"]["external_validation"]["cohens_kappa"]["kappa_status"],
            f"underpowered_n_lt_{MIN_KAPPA_N}",
        )
        for f in ("cohort_size", "split_type", "external_validation", "headline_metric", "readiness_tier"):
            self.assertEqual(
                res_ab["field_agreement"][f]["raw_agreement_all_rows"],
                res_ba["field_agreement"][f]["raw_agreement_all_rows"],
            )
            self.assertEqual(
                res_ab["field_agreement"][f]["raw_agreement_both_asserted"],
                res_ba["field_agreement"][f]["raw_agreement_both_asserted"],
            )


class TestEvaluatorCLI(unittest.TestCase):
    def test_cli_stress_and_compare_annotators(self) -> None:
        cli_script = READINESS_DIR / "evaluate_readiness.py"
        res = subprocess.run(
            [sys.executable, str(cli_script), "--stress", "--compare-baseline", "--json"],
            capture_output=True,
            text=True,
        )
        self.assertEqual(res.returncode, 0, msg=res.stderr)
        payload = json.loads(res.stdout)
        self.assertEqual(payload["dataset_kind"], "synthetic_stress_fixture")
        self.assertFalse(payload["is_unseen_published_papers"])
        self.assertEqual(payload["heuristic_readiness_tier"]["t3_overclaim_count"], 0)
        self.assertGreater(payload["baseline_comparison"]["t3_overclaim_reduction"], 0)

        with tempfile.TemporaryDirectory() as tmp:
            p_a = Path(tmp) / "a.json"
            p_b = Path(tmp) / "b.json"
            base_pay = {
                "schema_version": SCHEMA_VERSION,
                "dataset_kind": "unseen_review",
                "annotator_id": "rev-1",
                "items": [_make_item("u-cli-1", "grp-cli-1")],
            }
            p_a.write_text(json.dumps(base_pay), encoding="utf-8")
            base_pay["annotator_id"] = "rev-2"
            p_b.write_text(json.dumps(base_pay), encoding="utf-8")
            res2 = subprocess.run(
                [sys.executable, str(cli_script), "--compare-annotators", str(p_a), str(p_b), "--json"],
                capture_output=True,
                text=True,
            )
            self.assertEqual(res2.returncode, 0, msg=res2.stderr)
            comp = json.loads(res2.stdout)
            self.assertEqual(comp["items_compared"], 1)
            self.assertEqual(comp["disagreements"], [])


if __name__ == "__main__":
    unittest.main()
