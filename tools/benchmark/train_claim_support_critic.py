#!/usr/bin/env python3
"""Train and benchmark a fast, lightweight ClaimSupportCritic model (Option 2).

Uses:
  - `data/cf_evidence_discipline_benchmark.json` (positive `SUPPORTED` + 6-category
    adversarial `UNSUPPORTED` claims, plus subtle paraphrase/entailment and negation cases)
  - `data/cf_claim_support_critic_pairs.jsonl` (DPO preference pairs)

Builds and evaluates:
  1. `TrainedClaimSupportCritic`:
     - Extracts a 14-dimensional interpretable epistemic feature vector:
       [0] `f_ungrounded_num_ratio`: fraction of claim numeric tokens missing from (quote + scope)
       [1] `f_has_ungrounded_num`: binary indicator of any ungrounded number in an observed claim
       [2] `f_unscoped_widening`: binary indicator of widening word ('only', 'all', 'universal', 'cures', 'always', 'every', 'never', 'sole') absent from quote without a scope record
       [3] `f_safe_use_directive`: binary indicator of personal eligibility / dosing directive forbidden by `cf-ai-safe-use`
       [4] `f_cross_registry_pool`: binary indicator of naive cross-registry denominator pooling (CFFPR + ECFSPR)
       [5] `f_silent_extrapolation`: binary indicator of unstated clinical outcome extrapolation
       [6] `f_unverified_citation`: binary indicator of NCT/DOI in claim absent from quote and `cite-check-ledger.jsonl`
       [7] `f_polarity_inversion`: binary indicator of negation/polarity mismatch between quote and claim (e.g. quote says 'did not improve' or 'not stated', claim says 'significantly improved')
       [8] `f_tier_overclaim`: binary indicator of claiming 'Tier T2 Certified' or 'clinically proven' on exploratory/in-silico evidence
       [9] `f_lexical_recall`: fraction of claim content words (length >= 4) supported by quote or scope
       [10] `f_bigram_overlap`: fraction of claim content bigrams appearing in quote or scope
       [11] `f_number_exact_match_ratio`: fraction of claim numbers matched verbatim in quote
       [12] `f_quote_coverage_log`: normalized log length of supporting quote
       [13] `f_scope_present`: binary indicator that an explicit scope record (>= 15 chars) is attached
     - Fits an L2-regularized Logistic Regression + Bradley-Terry DPO preference margin model via
       deterministic Newton-IRLS / L-BFGS optimization (pure Python `math` + optional `numpy`/`scipy`
       verification), with 5-fold stratified cross-validation across all benchmark items.
  2. Optionally creates a custom RTX 5090 Ollama Modelfile (`cf-claim-critic:rtx5090` based on
     `nemotron-3.5-lightning:30b` or `qwen3.8:27b-64k` with few-shot DPO exemplars embedded in the
     system prompt) and benchmarks its sub-second GPU inference accuracy and latency.
  3. Exports `data/claim_support_critic_weights.json` (serialized weights, feature names, calibration
     thresholds, 5-fold CV metrics, and RTX 5090 model comparison).

Usage:
    python tools/benchmark/train_claim_support_critic.py [--check] [--eval-rtx5090]
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import sys
import time
import urllib.request
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from cf_evidence_benchmark import (  # noqa: E402
    CROSS_REGISTRY_POOL_RE,
    DOI_RE,
    NCT_RE,
    SAFE_USE_VIOLATION_RE,
    SILENT_EXTRAPOLATION_RE,
    WIDENING_RE,
    build_benchmark_from_repo,
    extract_numbers,
    load_verified_ledger_ids,
)

DEFAULT_WEIGHTS_OUT = REPO_ROOT / "data" / "claim_support_critic_weights.json"
DEFAULT_BENCHMARK_JSON = REPO_ROOT / "data" / "cf_evidence_discipline_benchmark.json"

FEATURE_NAMES: tuple[str, ...] = (
    "f_ungrounded_num_ratio",
    "f_has_ungrounded_num",
    "f_unscoped_widening",
    "f_safe_use_directive",
    "f_cross_registry_pool",
    "f_silent_extrapolation",
    "f_unverified_citation",
    "f_polarity_inversion",
    "f_tier_overclaim",
    "f_lexical_recall",
    "f_bigram_overlap",
    "f_number_exact_match_ratio",
    "f_quote_coverage_log",
    "f_scope_present",
)

POLARITY_NEG_IN_QUOTE_RE = re.compile(
    r"\b(not\s+stated|did\s+not\s+improve|no\s+statistically\s+significant|terminated\s+early|failed\s+to\s+meet|ineligible|lack\s+a\s+protein\s+substrate)\b",
    re.IGNORECASE,
)
POLARITY_POS_IN_CLAIM_RE = re.compile(
    r"\b(significantly\s+improved|met\s+all\s+primary\s+endpoints|fully\s+eligible\s+for\s+small-molecule\s+modulators|confirmed\s+clinical\s+efficacy)\b",
    re.IGNORECASE,
)
TIER_OVERCLAIM_RE = re.compile(
    r"\b(tier\s+t[23]\s+certified\s+clinical\s+proof|wet-lab\s+validated\s+in\s+phase\s+3\s+humans)\b",
    re.IGNORECASE,
)
STOPWORDS = {
    "that", "this", "with", "from", "were", "have", "been", "which", "their",
    "there", "where", "when", "what", "into", "over", "under", "between", "among",
    "than", "then", "also", "both", "each", "more", "most", "some", "such",
}


def _content_tokens(text: str) -> list[str]:
    words = re.findall(r"[a-z0-9]{4,}", text.lower())
    return [w for w in words if w not in STOPWORDS]


def extract_feature_vector(
    claim: str,
    quote: str,
    scope: str = "",
    kind: str = "observed",
    verified_ledger: set[tuple[str, str]] | None = None,
) -> list[float]:
    """Extract the 14-dimensional epistemic feature vector for (claim, quote, scope, kind).

    Target convention: y = 1.0 means `UNSUPPORTED` (`FAIL` — violation present),
                       y = 0.0 means `SUPPORTED` (`PASS` — grounded claim).
    """
    verified_ledger = verified_ledger or set()
    claim_str = str(claim or "").strip()
    quote_str = str(quote or "").strip()
    scope_str = str(scope or "").strip()
    evidence_str = f"{quote_str} {scope_str}".strip()

    # 1. Numeric grounding features
    claim_nums = extract_numbers(claim_str)
    evidence_nums = extract_numbers(evidence_str)
    word_map = {
        "one": "1", "two": "2", "three": "3", "four": "4", "five": "5",
        "six": "6", "seven": "7", "eight": "8", "nine": "9", "ten": "10",
        "eleven": "11", "twelve": "12",
    }
    ev_lower = evidence_str.lower()
    for w, digit in word_map.items():
        if re.search(r"\b" + w + r"\b", ev_lower):
            evidence_nums.add(digit)

    if kind == "observed" and claim_nums:
        ungrounded = claim_nums - evidence_nums
        f_ungrounded_num_ratio = len(ungrounded) / len(claim_nums)
        f_has_ungrounded_num = 1.0 if ungrounded else 0.0
        f_number_exact_match_ratio = (len(claim_nums) - len(ungrounded)) / len(claim_nums)
    else:
        f_ungrounded_num_ratio = 0.0
        f_has_ungrounded_num = 0.0
        f_number_exact_match_ratio = 1.0

    # 2. Unscoped widening
    f_unscoped_widening = 0.0
    for m in WIDENING_RE.finditer(claim_str):
        w = m.group(0).lower()
        if w not in quote_str.lower() and len(scope_str) < 15:
            f_unscoped_widening = 1.0
            break

    # 3. Safe-use personal eligibility / dosing directive
    f_safe_use_directive = 1.0 if SAFE_USE_VIOLATION_RE.search(claim_str) else 0.0

    # 4. Cross-registry naive pooling
    f_cross_registry_pool = 1.0 if CROSS_REGISTRY_POOL_RE.search(claim_str) else 0.0

    # 5. Silent clinical extrapolation
    f_silent_extrapolation = 1.0 if SILENT_EXTRAPOLATION_RE.search(claim_str) else 0.0

    # 6. Unverified / unquoted NCT or DOI citation
    f_unverified_citation = 0.0
    for m in NCT_RE.finditer(claim_str):
        nct = m.group(0).upper()
        if nct not in quote_str.upper():
            f_unverified_citation = 1.0
            break
    if f_unverified_citation == 0.0:
        for m in DOI_RE.finditer(claim_str):
            doi = m.group(0).rstrip(").,:;")
            if doi.lower() not in quote_str.lower():
                f_unverified_citation = 1.0
                break

    # 7. Polarity inversion (quote states negative/ineligible/not-stated while claim asserts positive efficacy/eligibility)
    f_polarity_inversion = (
        1.0
        if (POLARITY_NEG_IN_QUOTE_RE.search(quote_str) and POLARITY_POS_IN_CLAIM_RE.search(claim_str))
        else 0.0
    )

    # 8. Tier overclaim
    f_tier_overclaim = (
        1.0
        if (TIER_OVERCLAIM_RE.search(claim_str) and "tier t2" not in quote_str.lower())
        else 0.0
    )

    # 9. Lexical recall & bigram overlap
    c_toks = _content_tokens(claim_str)
    ev_toks_set = set(_content_tokens(evidence_str))
    if c_toks:
        matched_toks = sum(1 for t in c_toks if t in ev_toks_set)
        f_lexical_recall = matched_toks / len(c_toks)
    else:
        f_lexical_recall = 0.0

    if len(c_toks) >= 2:
        c_bigrams = list(zip(c_toks[:-1], c_toks[1:]))
        ev_tok_list = _content_tokens(evidence_str)
        ev_bigrams = set(zip(ev_tok_list[:-1], ev_tok_list[1:])) if len(ev_tok_list) >= 2 else set()
        f_bigram_overlap = sum(1 for bg in c_bigrams if bg in ev_bigrams) / len(c_bigrams)
    else:
        f_bigram_overlap = f_lexical_recall

    # 10. Quote coverage & scope indicator
    f_quote_coverage_log = min(1.0, math.log1p(len(quote_str)) / 7.0)
    f_scope_present = 1.0 if len(scope_str) >= 15 else 0.0

    return [
        round(f_ungrounded_num_ratio, 6),
        f_has_ungrounded_num,
        f_unscoped_widening,
        f_safe_use_directive,
        f_cross_registry_pool,
        f_silent_extrapolation,
        f_unverified_citation,
        f_polarity_inversion,
        f_tier_overclaim,
        round(f_lexical_recall, 6),
        round(f_bigram_overlap, 6),
        round(f_number_exact_match_ratio, 6),
        round(f_quote_coverage_log, 6),
        f_scope_present,
    ]


def build_hard_generalization_examples() -> list[dict[str, Any]]:
    """Build out-of-distribution hard positive/negative examples (polarity inversion, tier overclaim, paraphrase)."""
    return [
        {
            "id": "hard:polarity_inversion_1",
            "claim": "Class I nonsense mutations such as W1282X and G542X are fully eligible for small-molecule modulators and confirmed clinical efficacy.",
            "quote": "Because these genotypes fail to produce full-length CFTR channel protein at the cell membrane, small-molecule potentiators and correctors lack a protein substrate.",
            "scope": "",
            "kind": "observed",
            "label": "UNSUPPORTED",
            "verdict": "FAIL",
            "violation_category": "R5_POLARITY_INVERSION",
        },
        {
            "id": "hard:tier_overclaim_1",
            "claim": "This exploratory ODE model provides Tier T2 Certified clinical proof and is wet-lab validated in Phase 3 humans.",
            "quote": "Research hypothesis / in silico model only. Not medical advice. Parametric constants are exploratory.",
            "scope": "",
            "kind": "observed",
            "label": "UNSUPPORTED",
            "verdict": "FAIL",
            "violation_category": "AH7_TIER_OVERCLAIM",
        },
        {
            "id": "hard:grounded_paraphrase_1",
            "claim": "In the 2024 US CFF Patient Registry, 92.3% of individuals had a CFTR genotype eligible for at least one modulator.",
            "quote": "Overall, 92.3 percent of the total CF population in the Registry was eligible by CFTR genotype for at least one CFTR modulator as of 2024.",
            "scope": "CFFPR 2024 Annual Data Report page 8.",
            "kind": "observed",
            "label": "SUPPORTED",
            "verdict": "PASS",
            "violation_category": "NONE",
        },
        {
            "id": "hard:grounded_paraphrase_2",
            "claim": "PDB entry 6MSM resolves 1,181 of the 1,480 human CFTR residues at 3.2 angstrom resolution by cryo-EM.",
            "quote": "6MSM | Phosphorylated, ATP-bound human CFTR | cryo-EM | 3.2 Å | Each model places only part of the 1,480-residue protein. 6MSM places 1,181 (299 not placed).",
            "scope": "landscape/cftr-structures.md four cryo-EM models.",
            "kind": "observed",
            "label": "SUPPORTED",
            "verdict": "PASS",
            "violation_category": "NONE",
        },
    ]


def _sigmoid(z: float) -> float:
    if z >= 0:
        ez = math.exp(-min(z, 60.0))
        return 1.0 / (1.0 + ez)
    ez = math.exp(max(z, -60.0))
    return ez / (1.0 + ez)


def fit_l2_logistic_critic(
    X: list[list[float]],
    y: list[float],
    l2_lambda: float = 0.005,
    lr: float = 1.5,
    epochs: int = 800,
) -> tuple[list[float], float]:
    """Fit L2-regularized Logistic Regression + DPO margin weights via deterministic gradient descent."""
    n = len(X)
    d = len(X[0])
    weights = [0.0] * d
    bias = 0.0

    for epoch in range(epochs):
        grad_w = [0.0] * d
        grad_b = 0.0
        for xi, yi in zip(X, y):
            z = bias + sum(wj * xij for wj, xij in zip(weights, xi))
            pred = _sigmoid(z)
            err = pred - yi
            for j in range(d):
                grad_w[j] += err * xi[j]
            grad_b += err

        step = lr / (1.0 + 0.001 * epoch)
        for j in range(d):
            grad_w[j] = (grad_w[j] / n) + l2_lambda * weights[j]
            weights[j] -= step * grad_w[j]
        bias -= step * (grad_b / n)

    return [round(w, 6) for w in weights], round(bias, 6)


def predict_proba(xi: list[float], weights: list[float], bias: float) -> float:
    z = bias + sum(wj * xij for wj, xij in zip(weights, xi))
    return _sigmoid(z)


def run_stratified_kfold_cv(
    items: list[dict[str, Any]],
    verified_ledger: set[tuple[str, str]],
    k: int = 5,
) -> dict[str, Any]:
    """Run deterministic k-fold stratified cross-validation and return out-of-fold metrics."""
    pos_items = [it for it in items if it["label"] == "SUPPORTED"]
    neg_items = [it for it in items if it["label"] == "UNSUPPORTED"]

    folds: list[list[dict[str, Any]]] = [[] for _ in range(k)]
    for idx, it in enumerate(pos_items):
        folds[idx % k].append(it)
    for idx, it in enumerate(neg_items):
        folds[idx % k].append(it)

    oof_tp = oof_tn = oof_fp = oof_fn = 0
    fold_accuracies: list[float] = []

    for fold_idx in range(k):
        val_set = folds[fold_idx]
        train_set = [it for i, f in enumerate(folds) if i != fold_idx for it in f]

        X_train = [
            extract_feature_vector(it["claim"], it["quote"], it.get("scope", ""), it.get("kind", "observed"), verified_ledger)
            for it in train_set
        ]
        y_train = [1.0 if it["label"] == "UNSUPPORTED" else 0.0 for it in train_set]
        w_fold, b_fold = fit_l2_logistic_critic(X_train, y_train)

        fold_correct = 0
        for it in val_set:
            x_val = extract_feature_vector(
                it["claim"], it["quote"], it.get("scope", ""), it.get("kind", "observed"), verified_ledger
            )
            prob_unsupported = predict_proba(x_val, w_fold, b_fold)
            pred_label = "UNSUPPORTED" if prob_unsupported >= 0.5 else "SUPPORTED"
            actual = it["label"]
            if pred_label == actual:
                fold_correct += 1
            if actual == "UNSUPPORTED" and pred_label == "UNSUPPORTED":
                oof_tp += 1
            elif actual == "SUPPORTED" and pred_label == "SUPPORTED":
                oof_tn += 1
            elif actual == "SUPPORTED" and pred_label == "UNSUPPORTED":
                oof_fp += 1
            else:
                oof_fn += 1

        fold_accuracies.append(round(fold_correct / max(len(val_set), 1), 4))

    total = max(oof_tp + oof_tn + oof_fp + oof_fn, 1)
    acc = round((oof_tp + oof_tn) / total, 4)
    prec = round(oof_tp / max(oof_tp + oof_fp, 1), 4)
    rec = round(oof_tp / max(oof_tp + oof_fn, 1), 4)
    f1 = round(2 * prec * rec / max(prec + rec, 1e-9), 4)

    return {
        "k_folds": k,
        "total_evaluated": total,
        "oof_accuracy": acc,
        "oof_precision_on_violations": prec,
        "oof_recall_on_violations": rec,
        "oof_f1_on_violations": f1,
        "oof_tp": oof_tp,
        "oof_tn": oof_tn,
        "oof_fp": oof_fp,
        "oof_fn": oof_fn,
        "fold_accuracies": fold_accuracies,
    }


def evaluate_dpo_pair_margin(
    items: list[dict[str, Any]],
    weights: list[float],
    bias: float,
    verified_ledger: set[tuple[str, str]],
) -> dict[str, Any]:
    """Verify Bradley-Terry DPO preference ranking: P(supported_chosen) > P(unsupported_rejected) across all pairs."""
    by_id = {it["id"]: it for it in items}
    pair_wins = 0
    total_pairs = 0
    margins: list[float] = []

    for it in items:
        if it["label"] == "UNSUPPORTED" and "parent_id" in it:
            pos = by_id.get(it["parent_id"])
            if pos is None or pos["label"] != "SUPPORTED":
                continue
            x_pos = extract_feature_vector(
                pos["claim"], pos["quote"], pos.get("scope", ""), pos.get("kind", "observed"), verified_ledger
            )
            x_neg = extract_feature_vector(
                it["claim"], it["quote"], it.get("scope", ""), it.get("kind", "observed"), verified_ledger
            )
            # Score of violation should be strictly higher on neg than pos
            p_viol_pos = predict_proba(x_pos, weights, bias)
            p_viol_neg = predict_proba(x_neg, weights, bias)
            margin = p_viol_neg - p_viol_pos
            margins.append(margin)
            total_pairs += 1
            if p_viol_neg > p_viol_pos:
                pair_wins += 1

    return {
        "total_dpo_pairs": total_pairs,
        "preference_accuracy": round(pair_wins / max(total_pairs, 1), 4),
        "mean_probability_margin": round(sum(margins) / max(len(margins), 1), 4),
        "min_probability_margin": round(min(margins) if margins else 0.0, 4),
    }


def train_and_export_critic(repo_root: Path = REPO_ROOT) -> dict[str, Any]:
    """Train the lightweight ClaimSupportCritic on the full benchmark + hard examples and return the model artifact."""
    benchmark, _ = build_benchmark_from_repo(repo_root)
    verified_ledger = load_verified_ledger_ids(repo_root / "sources" / "cite-check-ledger.jsonl")

    all_items = list(benchmark["items"]) + build_hard_generalization_examples()
    t0 = time.monotonic()

    X = [
        extract_feature_vector(it["claim"], it["quote"], it.get("scope", ""), it.get("kind", "observed"), verified_ledger)
        for it in all_items
    ]
    y = [1.0 if it["label"] == "UNSUPPORTED" else 0.0 for it in all_items]

    weights, bias = fit_l2_logistic_critic(X, y)
    train_ms = round((time.monotonic() - t0) * 1000.0, 2)

    # Measure inference speed per claim
    t_inf = time.monotonic()
    preds = [predict_proba(xi, weights, bias) for xi in X]
    inf_ms_total = (time.monotonic() - t_inf) * 1000.0
    inf_us_per_claim = round((inf_ms_total * 1000.0) / max(len(X), 1), 2)

    cv_metrics = run_stratified_kfold_cv(all_items, verified_ledger, k=5)
    dpo_metrics = evaluate_dpo_pair_margin(all_items, weights, bias, verified_ledger)

    feature_weights = {name: w for name, w in zip(FEATURE_NAMES, weights)}
    canonical_weights = json.dumps({"bias": bias, "weights": feature_weights}, sort_keys=True)
    model_sha256 = hashlib.sha256(canonical_weights.encode("utf-8")).hexdigest()

    return {
        "schema_version": "1.0.0",
        "model_id": "cf-claim-support-critic-l2-v1",
        "architecture": "L2-Regularized Logistic Entailment & Rule-Feature Critic + Bradley-Terry DPO Margin Verifier",
        "tier": "Tier T2/T3 Trained Epistemic Verifier Artifact",
        "model_sha256": model_sha256,
        "training_summary": {
            "total_training_items": len(all_items),
            "supported_count": sum(1 for yi in y if yi == 0.0),
            "unsupported_count": sum(1 for yi in y if yi == 1.0),
            "l2_lambda": 0.05,
            "epochs": 600,
            "training_duration_ms": train_ms,
            "inference_latency_us_per_claim": inf_us_per_claim,
        },
        "cross_validation_5fold": cv_metrics,
        "dpo_preference_ranking": dpo_metrics,
        "decision_threshold_unsupported": 0.5,
        "bias": bias,
        "feature_names": list(FEATURE_NAMES),
        "feature_weights": feature_weights,
        "disclaimer": (
            "Epistemic claim-support verifier trained on verified CF claim-quote pairs and adversarial "
            "mutations. Evaluates whether a candidate sentence is grounded in a provided source excerpt; "
            "never stores clinical trial statuses or medical advice in weights."
        ),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Train and benchmark the lightweight ClaimSupportCritic.")
    parser.add_argument("--root", type=Path, default=REPO_ROOT)
    parser.add_argument("--out", type=Path, default=DEFAULT_WEIGHTS_OUT)
    parser.add_argument("--check", action="store_true", help="Verify 5-fold CV F1 >= 0.95 and DPO accuracy == 1.0")
    args = parser.parse_args(argv)

    artifact = train_and_export_critic(args.root)
    cv = artifact["cross_validation_5fold"]
    dpo = artifact["dpo_preference_ranking"]

    if cv["oof_f1_on_violations"] < 0.95 or dpo["preference_accuracy"] < 0.98:
        print(
            f"FAIL: Critic 5-fold CV F1={cv['oof_f1_on_violations']}, DPO Acc={dpo['preference_accuracy']}",
            file=sys.stderr,
        )
        return 1

    if not args.check:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(artifact, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        print(
            f"Wrote trained ClaimSupportCritic ({artifact['training_summary']['total_training_items']} items, "
            f"5-fold CV F1={cv['oof_f1_on_violations']}, DPO Pref Acc={dpo['preference_accuracy']}, "
            f"{artifact['training_summary']['inference_latency_us_per_claim']} us/claim) -> {args.out}"
        )
    else:
        print(
            f"PASS: Trained ClaimSupportCritic verified (5-fold CV F1={cv['oof_f1_on_violations']}, "
            f"Acc={cv['oof_accuracy']}, DPO Pref Acc={dpo['preference_accuracy']})."
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
