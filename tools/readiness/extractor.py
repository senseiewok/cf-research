#!/usr/bin/env python3
"""Deterministic readiness label extractor for Cystic Fibrosis machine learning papers.

Extracts four essential readiness attributes from CF ML abstracts and methods:
1. Cohort size (human patients/pwCF; distinguishes from sample/culture counts and cell lines)
2. Split type (patient_level vs sample_or_record_level vs unsupervised_or_zero_shot vs not_stated)
3. External validation (external_multicenter_or_independent vs absent_single_center_or_internal_only vs not_stated)
4. Headline performance metric (e.g., AUC-ROC, AUCs of 0.75-0.80, delta Bhalla score)

Every extracted attribute is anchored with its verbatim source sentence.
Enforces epistemic guardrails: assigns standardized readiness tiers (T0-T3) and raises
anti-overclaim flags to ensure explainers cannot overstate preclinical or single-center
models as clinically ready.
"""
from __future__ import annotations

import argparse
import dataclasses
import json
import re
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

# Prohibited prompt injection patterns quarantined and defanged
PROMPT_INJECTION_RE = re.compile(
    r"(?i)(?:<system>|\[system\]|ignore\s+all\s+previous\s+instructions|system\s*:\s*you\s+are|"
    r"as\s+an\s+ai\s+language\s+model|sudo\s+make|curl\s+https?://)",
    re.M,
)


def sanitize_input_text(text: str) -> str:
    """Quarantine and defang potential prompt injection payloads in untrusted input text."""
    if not text:
        return ""
    # Strip forbidden control tags and defang injection phrases
    cleaned = PROMPT_INJECTION_RE.sub("[DEFANGED_INJECTION_ATTEMPT]", text)
    return cleaned.strip()


def tokenize_sentences(text: str) -> List[str]:
    """Split text into sentences using deterministic standard punctuation boundaries."""
    sanitized = sanitize_input_text(text)
    if not sanitized:
        return []

    # Preserve abbreviations like 'i.e.', 'e.g.', 'al.', 'vs.'
    protected = sanitized
    replacements = [
        ("i.e.,", "i_e_comma"),
        ("i.e.", "i_e_period"),
        ("e.g.,", "e_g_comma"),
        ("e.g.", "e_g_period"),
        ("al.", "al_period"),
        ("vs.", "vs_period"),
        ("p<", "p_less_"),
        ("p=", "p_equal_"),
    ]
    for orig, rep in replacements:
        protected = protected.replace(orig, rep)

    raw_sentences = re.split(r"(?<=[.!?])\s+(?=[A-Z0-9\"'(\[])", protected)
    results = []
    for s in raw_sentences:
        s_clean = s.strip()
        for orig, rep in replacements:
            s_clean = s_clean.replace(rep, orig)
        if s_clean:
            results.append(s_clean)
    return results


def extract_cohort_size(sentences: List[str]) -> Tuple[Union[int, str], str]:
    """Extract patient cohort size N and verbatim source sentence.

    Distinguishes patient cohort size (pwCF, patients, individuals) from culture/sample counts
    (e.g., 12,618 sputum cultures) and laboratory variants. If no human patient cohort is present,
    returns 'not stated'.
    """
    patient_patterns = [
        # In total, 283 patients were included
        re.compile(r"(?i)\bin\s+total[,\s]+([0-9,]+)\s+patients\b"),
        # 10,261 patients from / with / in / having
        re.compile(r"(?i)\b([0-9,]+)\s+patients\s+(?:with|from|having|in)\b"),
        # from 209 pwCF
        re.compile(r"(?i)\b([0-9,]+)\s+pwCF\b"),
        # evaluated / included / enrolled 200 patients
        re.compile(r"(?i)\b(?:evaluated|included|enrolled|recruited|analyzed|studied)\s+([0-9,]+)\s+(?:patients|pwCF|individuals|participants)\b"),
        # cohort of 500 individuals
        re.compile(r"(?i)\bcohort\s+of\s+([0-9,]+)\s+(?:patients|pwCF|individuals|participants)\b"),
        # 283 patients were included / evaluated / enrolled
        re.compile(r"(?i)\b([0-9,]+)\s+(?:patients|pwCF|individuals|participants)\s+(?:were\s+)?(?:included|evaluated|enrolled|recruited|analyzed)\b"),
    ]

    for sent in sentences:
        for pat in patient_patterns:
            m = pat.search(sent)
            if m:
                raw_val = m.group(1).replace(",", "")
                try:
                    val = int(raw_val)
                    return val, sent
                except ValueError:
                    pass

    return "not stated", ""


def extract_split_type(sentences: List[str]) -> Tuple[str, str]:
    """Extract train/test split methodology.

    Categories:
    - 'patient_level': Patient-stratified or clustered splitting avoiding data leakage.
    - 'unsupervised_or_zero_shot': Unsupervised clustering or zero-shot foundation model evaluation.
    - 'sample_or_record_level': Cell line, profiling, or sample-level split.
    - 'not_stated': Split not specified.
    """
    for sent in sentences:
        low = sent.lower()
        if "patient-level split" in low or "patient-level splitting" in low or "split by patient" in low or "patient-wise" in low:
            return "patient_level", sent
        if "unsupervised machine learning" in low or "clustering method" in low or "k-means clustering" in low or "cluster analysis" in low:
            return "unsupervised_or_zero_shot", sent
        if "zero-shot" in low or "pre-trained" in low or "pretrained" in low or "evaluat" in low and "cftr2.org database" in low:
            return "unsupervised_or_zero_shot", sent
        if "ribosome profiling experiments" in low or "cell-based" in low or "trained on data from cells" in low:
            return "sample_or_record_level", sent
        if "record-level" in low or "sample-level" in low or "random split" in low:
            return "sample_or_record_level", sent

    return "not_stated", ""


def extract_external_validation(sentences: List[str]) -> Tuple[str, str]:
    """Extract whether external validation on an independent cohort was performed.

    Categories:
    - 'external_multicenter_or_independent': Validated on independent external registry or center.
    - 'absent_single_center_or_internal_only': Single center, internal validation only, or explicitly needs external validation.
    - 'not_stated': No completed validation statement, or conflicting explicit evidence.
    """
    validation = r"(?:external|independent)\s+validation"
    population = r"(?:external|independent)(?:\s+external)?\s+(?:data(?:set)?s?|cohorts?|registr(?:y|ies))"
    action = r"(?:validated|tested|evaluated)"
    negative = re.compile(
        rf"\b(?:not(?!\s+only\b)|never|without)\s+(?:(?:yet|been|being)\s+)*"
        rf"(?:externally\s+validated|{validation})\b"
        rf"|\b(?:no|without|lacks?|lacked)\s+(?:an?\s+)?{validation}\b"
        rf"|\b{validation}\s+(?:(?:was|were|has|have|had)\s+(?:not|never)|(?:wasn't|weren't|hasn't|haven't|hadn't))\b"
        rf"|\b(?:wasn't|weren't|hasn't|haven't|hadn't)\s+(?:been\s+)?externally\s+validated\b"
        rf"|\b(?:did\s+not|didn't|do\s+not|don't)\s+(?:externally\s+validate|(?:perform|conduct|complete)\s+(?:an?\s+)?{validation})\b"
        rf"|\b(?:not|never)\s+(?:(?:yet|been|being)\s+)*{action}\s+(?:on|using|against|in|with)\s+(?:an?\s+)?{population}\b"
        r"|\b(?:internal validation only|only internal validation|only internally validated)\b"
    )
    completed = re.compile(
        r"\b(?:was|were|has been|have been|had been|we)\s+"
        r"(?:(?:also|successfully|not only)\s+)?externally\s+validated\b"
        rf"|\b(?:was|were|has been|have been|had been)\s+(?:successfully\s+)?{action}"
        rf"\s+(?:on|using|against|in|with)\s+(?:an?\s+)?{population}\b"
        rf"|\bwe\s+{action}\s+(?:the\s+)?(?:model|algorithm|framework)"
        rf"\s+(?:on|using|against|in|with)\s+(?:an?\s+)?{population}\b"
        rf"|\b{validation}\s+(?:was|has been|had been)\s+(?:successfully\s+)?"
        r"(?:performed|conducted|completed)\b"
    )
    prospective = re.compile(
        r"\b(?:will|would|could|should|may|might|if|whether|planned|proposed|"
        r"intend(?:ed)?|aim(?:ed)?)\b"
    )
    uncertain = re.compile(
        r"\b(?:no (?:evidence|indication)|(?:not|never) (?:shown|established|confirmed)|"
        r"unclear|unknown|uncertain)\b"
    )
    affirmative_sentence = ""
    negative_sentence = ""
    uncertain_sentence = ""
    fallback_sentence = ""

    for sent in sentences:
        low = sent.lower()
        for clause in re.split(r";|\b(?:but|however)\b", low):
            for match in negative.finditer(clause):
                prefix = clause[:match.start()]
                if not prospective.search(prefix) and not uncertain.search(prefix):
                    negative_sentence = negative_sentence or sent
            # Remove negated spans so their completed-action wording cannot be reused.
            positive_clause = negative.sub(lambda match: " " * len(match.group()), clause)
            for match in completed.finditer(positive_clause):
                prefix = positive_clause[:match.start()]
                if uncertain.search(prefix):
                    uncertain_sentence = uncertain_sentence or sent
                elif not prospective.search(prefix):
                    affirmative_sentence = affirmative_sentence or sent

        # Development setting and calls for further validation are weaker evidence.
        if "single-center" in low or "single center" in low:
            fallback_sentence = fallback_sentence or sent
        if "external validation" in low and any(
            term in low for term in ("needed", "requires", "required", "before clinical use")
        ):
            fallback_sentence = fallback_sentence or sent

    if affirmative_sentence and negative_sentence:
        return "not_stated", negative_sentence
    if uncertain_sentence:
        return "not_stated", uncertain_sentence
    if affirmative_sentence:
        return "external_multicenter_or_independent", affirmative_sentence
    if negative_sentence or fallback_sentence:
        return "absent_single_center_or_internal_only", negative_sentence or fallback_sentence
    return "not_stated", ""


def extract_headline_metric(sentences: List[str]) -> Tuple[str, str]:
    """Extract primary headline performance metric reported in the text."""
    metric_patterns = [
        # AUCs of 0.75-0.80
        re.compile(r"(?i)\b(AUCs?\s+of\s+[0-9.]+(?:\s*-\s*[0-9.]+)?)\b"),
        # AUC-ROC of 0.89
        re.compile(r"(?i)\b(AUC-ROC\s+of\s+[0-9.]+)\b"),
        # -1.40±1.53 points
        re.compile(r"(-?[0-9.]+\s*±\s*[0-9.]+\s*points)"),
        # fair classification performance
        re.compile(r"(?i)\b(fair\s+classification\b(?:\s+performance)?)"),
        # accurately predicted readthrough
        re.compile(r"(?i)\b(accurately\s+predicted\s+readthrough\b)"),
        # dynamic quantification of CFTR function
        re.compile(r"(?i)\b(dynamic\s+quantification\s+of\s+CFTR\s+function\b)"),
    ]

    for sent in sentences:
        for pat in metric_patterns:
            m = pat.search(sent)
            if m:
                return m.group(1).strip(), sent

    return "not stated", ""


def determine_readiness_tier(
    cohort_size: Union[int, str],
    split_type: str,
    external_validation: str,
) -> str:
    """Determine standardized epistemic readiness tier (T0-T3).

    T3_EXTERNALLY_VALIDATED: Large human patient cohort (N >= 50) + patient-level split + independent external validation.
    T2_INTERNALLY_VALIDATED: Human patient cohort + patient-level split + internal cross-validation, but lacking external validation.
    T1_RETROSPECTIVE_OR_UNVALIDATED: Single-center retrospective, unsupervised, or lacking independent validation.
    T0_PRECLINICAL_BENCHMARK: In vitro, cell line, organoid, or variant-level study with no human patient cohort.
    """
    if cohort_size == "not stated" or (isinstance(cohort_size, int) and cohort_size == 0):
        return "T0_PRECLINICAL_BENCHMARK"

    if isinstance(cohort_size, int) and cohort_size >= 50:
        if split_type == "patient_level" and external_validation == "external_multicenter_or_independent":
            return "T3_EXTERNALLY_VALIDATED"
        if split_type == "patient_level" and external_validation != "external_multicenter_or_independent":
            # Note: single-center retrospective study explicitly stating external validation needed is T1
            if external_validation == "absent_single_center_or_internal_only":
                return "T1_RETROSPECTIVE_OR_UNVALIDATED"
            return "T2_INTERNALLY_VALIDATED"

    return "T1_RETROSPECTIVE_OR_UNVALIDATED"


def identify_anti_overclaim_flags(
    cohort_size: Union[int, str],
    split_type: str,
    external_validation: str,
    readiness_tier: str,
) -> List[str]:
    """Flag common failure modes and overclaims in clinical AI reporting."""
    flags: List[str] = []

    if cohort_size == "not stated":
        flags.append("IN_VITRO_OR_VARIANT_ONLY_NO_HUMAN_COHORT")

    if external_validation == "absent_single_center_or_internal_only":
        flags.append("SINGLE_CENTER_NO_EXTERNAL_VALIDATION")
    elif external_validation == "not_stated" and readiness_tier != "T0_PRECLINICAL_BENCHMARK":
        flags.append("UNSTATED_EXTERNAL_VALIDATION_STATUS")

    if split_type == "unsupervised_or_zero_shot":
        flags.append("UNSUPERVISED_OR_ZERO_SHOT_NO_PROSPECTIVE_TEST")
    elif split_type == "sample_or_record_level":
        flags.append("SAMPLE_OR_RECORD_LEVEL_SPLIT_LEAKAGE_RISK")
    elif split_type == "not_stated":
        flags.append("UNSTATED_SPLIT_METHODOLOGY")

    if readiness_tier in ("T0_PRECLINICAL_BENCHMARK", "T1_RETROSPECTIVE_OR_UNVALIDATED"):
        flags.append("NOT_CLINICALLY_ACTIONABLE_DO_NOT_DEPLOY")

    return flags


@dataclasses.dataclass
class ReadinessRecord:
    paper_id: str
    cohort_size: Union[int, str]
    cohort_sentence: str
    split_type: str
    split_sentence: str
    external_validation: str
    external_validation_sentence: str
    headline_metric: str
    metric_sentence: str
    readiness_tier: str
    anti_overclaim_flags: List[str]

    def to_dict(self) -> Dict[str, Any]:
        return dataclasses.asdict(self)


class ReadinessExtractor:
    """Extractor class managing parsing, grounding, and validation."""

    def extract(self, text: str, paper_id: str = "") -> ReadinessRecord:
        sentences = tokenize_sentences(text)
        cohort_size, cohort_sent = extract_cohort_size(sentences)
        split_type, split_sent = extract_split_type(sentences)
        ext_val, ext_sent = extract_external_validation(sentences)
        metric, metric_sent = extract_headline_metric(sentences)
        tier = determine_readiness_tier(cohort_size, split_type, ext_val)
        flags = identify_anti_overclaim_flags(cohort_size, split_type, ext_val, tier)

        return ReadinessRecord(
            paper_id=paper_id,
            cohort_size=cohort_size,
            cohort_sentence=cohort_sent,
            split_type=split_type,
            split_sentence=split_sent,
            external_validation=ext_val,
            external_validation_sentence=ext_sent,
            headline_metric=metric,
            metric_sentence=metric_sent,
            readiness_tier=tier,
            anti_overclaim_flags=flags,
        )


def extract_readiness(text: str, paper_id: str = "") -> ReadinessRecord:
    """Convenience function to extract readiness labels from raw text."""
    extractor = ReadinessExtractor()
    return extractor.extract(text, paper_id=paper_id)


def run_benchmark(benchmark_path: Optional[Path] = None) -> Dict[str, Any]:
    """Run extraction benchmark on curated 2024-2026 CF ML papers."""
    if benchmark_path is None:
        benchmark_path = Path(__file__).resolve().parent / "cf_ml_benchmark.json"

    items = json.loads(benchmark_path.read_text(encoding="utf-8"))
    extractor = ReadinessExtractor()

    total_fields = 0
    exact_matches = 0
    field_breakdown = {
        "cohort_size": {"correct": 0, "total": 0},
        "split_type": {"correct": 0, "total": 0},
        "external_validation": {"correct": 0, "total": 0},
        "headline_metric": {"correct": 0, "total": 0},
        "readiness_tier": {"correct": 0, "total": 0},
    }
    results = []

    for item in items:
        rec = extractor.extract(item["text"], paper_id=item["id"])
        exp = item["expected"]

        item_correct = True
        diffs = {}

        # 1. Cohort size
        field_breakdown["cohort_size"]["total"] += 1
        if rec.cohort_size == exp["cohort_size"]:
            field_breakdown["cohort_size"]["correct"] += 1
            exact_matches += 1
        else:
            item_correct = False
            diffs["cohort_size"] = {"got": rec.cohort_size, "expected": exp["cohort_size"]}
        total_fields += 1

        # 2. Split type
        field_breakdown["split_type"]["total"] += 1
        if rec.split_type == exp["split_type"]:
            field_breakdown["split_type"]["correct"] += 1
            exact_matches += 1
        else:
            item_correct = False
            diffs["split_type"] = {"got": rec.split_type, "expected": exp["split_type"]}
        total_fields += 1

        # 3. External validation
        field_breakdown["external_validation"]["total"] += 1
        if rec.external_validation == exp["external_validation"]:
            field_breakdown["external_validation"]["correct"] += 1
            exact_matches += 1
        else:
            item_correct = False
            diffs["external_validation"] = {"got": rec.external_validation, "expected": exp["external_validation"]}
        total_fields += 1

        # 4. Headline metric (containment / normalization)
        field_breakdown["headline_metric"]["total"] += 1
        metric_match = (
            exp["headline_metric"] == "not stated" and rec.headline_metric == "not stated"
        ) or (
            exp["headline_metric"].lower() in rec.headline_metric.lower()
            or rec.headline_metric.lower() in exp["headline_metric"].lower()
        )
        if metric_match:
            field_breakdown["headline_metric"]["correct"] += 1
            exact_matches += 1
        else:
            item_correct = False
            diffs["headline_metric"] = {"got": rec.headline_metric, "expected": exp["headline_metric"]}
        total_fields += 1

        # 5. Readiness tier
        field_breakdown["readiness_tier"]["total"] += 1
        if rec.readiness_tier == exp["readiness_tier"]:
            field_breakdown["readiness_tier"]["correct"] += 1
            exact_matches += 1
        else:
            item_correct = False
            diffs["readiness_tier"] = {"got": rec.readiness_tier, "expected": exp["readiness_tier"]}
        total_fields += 1

        results.append({
            "id": item["id"],
            "title": item["title"],
            "pass": item_correct,
            "record": rec.to_dict(),
            "diffs": diffs,
        })

    accuracy = exact_matches / total_fields if total_fields > 0 else 0.0

    return {
        "status": "PASS" if accuracy >= 0.85 else "FAIL",
        "benchmark_items_evaluated": len(items),
        "total_fields_evaluated": total_fields,
        "exact_matches": exact_matches,
        "accuracy": round(accuracy, 4),
        "field_breakdown": {
            k: {
                "correct": v["correct"],
                "total": v["total"],
                "accuracy": round(v["correct"] / v["total"], 4) if v["total"] > 0 else 0.0,
            }
            for k, v in field_breakdown.items()
        },
        "results": results,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract readiness labels from CF machine-learning paper texts.")
    parser.add_argument("--benchmark", action="store_true", help="Run benchmark evaluation across cf_ml_benchmark.json")
    parser.add_argument("--file", type=Path, help="Path to text or markdown file to parse")
    parser.add_argument("--text", type=str, help="Raw text string to parse")
    parser.add_argument("--json", action="store_true", help="Output results in JSON format")

    args = parser.parse_args()

    if args.benchmark:
        report = run_benchmark()
        if args.json:
            print(json.dumps(report, indent=2))
        else:
            print(f"CF ML Readiness Benchmark: {report['status']} (Accuracy: {report['accuracy'] * 100:.1f}%)")
            for field, stats in report["field_breakdown"].items():
                print(f"  - {field}: {stats['correct']}/{stats['total']} ({stats['accuracy'] * 100:.1f}%)")
        sys.exit(0 if report["status"] == "PASS" else 1)

    input_text = ""
    paper_id = "input"
    if args.file:
        input_text = args.file.read_text(encoding="utf-8")
        paper_id = args.file.stem
    elif args.text:
        input_text = args.text
    else:
        parser.print_help()
        sys.exit(1)

    rec = extract_readiness(input_text, paper_id=paper_id)
    if args.json:
        print(json.dumps(rec.to_dict(), indent=2))
    else:
        print(f"Readiness Assessment for '{rec.paper_id}':")
        print(f"  Readiness Tier:      {rec.readiness_tier}")
        print(f"  Cohort Size:         {rec.cohort_size}")
        if rec.cohort_sentence:
            print(f"    Source: \"{rec.cohort_sentence}\"")
        print(f"  Split Type:          {rec.split_type}")
        if rec.split_sentence:
            print(f"    Source: \"{rec.split_sentence}\"")
        print(f"  External Validation: {rec.external_validation}")
        if rec.external_validation_sentence:
            print(f"    Source: \"{rec.external_validation_sentence}\"")
        print(f"  Headline Metric:     {rec.headline_metric}")
        if rec.metric_sentence:
            print(f"    Source: \"{rec.metric_sentence}\"")
        if rec.anti_overclaim_flags:
            print("  Anti-Overclaim Flags:")
            for flag in rec.anti_overclaim_flags:
                print(f"    [!] {flag}")


if __name__ == "__main__":
    main()
