#!/usr/bin/env python3
"""Run real tasks against the local worker (qwen3.8:27b-64k) with and without lab checks.

Captures actual local worker generations, tests them against the lab's deterministic
checkers, records failures and passes, and demonstrates the practical necessity of R1-R5.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys
import os
import urllib.error
import urllib.request
from typing import Any, Dict, List

BENCHMARK_DIR = Path(__file__).resolve().parent
TOOLS_DIR = BENCHMARK_DIR.parent
CLAIMS_TOOL_DIR = TOOLS_DIR / "claims"

if str(CLAIMS_TOOL_DIR) not in sys.path:
    sys.path.insert(0, str(CLAIMS_TOOL_DIR))

import check_claims


DEFAULT_HOST = os.environ.get("OLLAMA_HOST", "http://127.0.0.1:11434").rstrip("/")
OLLAMA_URL = f"{DEFAULT_HOST}/api/generate"
DEFAULT_MODEL = os.environ.get("LOCAL_WORKER_MODEL", "qwen3.8:27b-64k")


def query_local_worker(prompt: str, model: str = DEFAULT_MODEL, timeout: int = 60) -> str:
    """Send prompt to local Ollama worker via loopback HTTP with graceful fallback."""
    payload = json.dumps({
        "model": model,
        "prompt": prompt,
        "stream": False,
        "options": {
            "temperature": 0.2,
            "num_predict": 1024,
        }
    }).encode("utf-8")
    
    req = urllib.request.Request(
        OLLAMA_URL,
        data=payload,
        headers={"Content-Type": "application/json"}
    )
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            data = json.loads(resp.read().decode("utf-8"))
            return data.get("response", "").strip()
    except (urllib.error.URLError, TimeoutError) as e:
        sys.stderr.write(f"Warning: Could not contact local model at {OLLAMA_URL}: {e}\n")
        return ""


def run_worker_evaluation(output_path: Path | None = None) -> Dict[str, Any]:
    sources_path = BENCHMARK_DIR / "sources_fixture.json"
    with open(sources_path, "r", encoding="utf-8") as f:
        sources = json.load(f)

    ctx = {
        "evidence": {eid: check_claims.normalise_text(t) for eid, t in sources.items()},
        "evidence_errors": {},
        "id_counts": {},
        "max_quote_words": 25,
        "widening": check_claims.DEFAULT_WIDENING,
        "absence": check_claims.DEFAULT_ABSENCE,
    }

    test_scenarios = [
        {
            "id": "WORKER-001",
            "title": "Unconstrained summary vs verbatim extraction of patient trial count",
            "source_id": "ctgov-nct05033080",
            "prompt": (
                "You are writing a Cystic Fibrosis research note. Based on this text:\n"
                "'NCT05033080: A Phase 3, Randomized, Double-Blind, Controlled Study Evaluating Efficacy and Safety. A total of 510 patients were randomized across 115 clinical sites.'\n\n"
                "Write a 1-sentence observed claim summarizing the patient enrollment."
            ),
            "expected_evidence_numbers": ["510", "115"]
        },
        {
            "id": "WORKER-002",
            "title": "Modulator uptake percentage extraction from registry report",
            "source_id": "ecfspr-2024-modulators",
            "prompt": (
                "Based on this registry statement:\n"
                "'In 2024, approximately 82.4 percent of eligible individuals aged 6 years and older were receiving CFTR modulator therapy in reporting centres.'\n\n"
                "Provide a one-sentence claim stating CFTR modulator therapy uptake."
            ),
            "expected_evidence_numbers": ["2024", "82.4", "6"]
        },
        {
            "id": "WORKER-003",
            "title": "Strict structured claim JSON generation with exact quote instruction",
            "source_id": "trikafta-spl-warnings",
            "prompt": (
                "From this FDA label passage:\n"
                "'Elevations in liver transaminases (ALT and AST) have been reported in patients treated with TRIKAFTA.'\n\n"
                "Return a JSON object with keys 'claim', 'quote', 'kind'.\n"
                "The 'quote' MUST be an exact copy of the text without modifying any words or punctuation.\n"
                "'kind' must be 'observed'. Return JSON only."
            ),
            "expected_evidence_numbers": []
        },
        {
            "id": "WORKER-004",
            "title": "Extrapolation and number hallucination without full evidence",
            "source_id": "ctgov-nct05033080",
            "prompt": (
                "From this text only: 'A total of 510 patients were randomized across 115 clinical sites.'\n"
                "State how many patients completed the full 24-week study period."
            ),
            "expected_evidence_numbers": []
        },
        {
            "id": "WORKER-005",
            "title": "Unhedged inference on biological restoration",
            "source_id": "ctgov-sweat-chloride",
            "prompt": (
                "From this statement: 'Sweat chloride concentration decreased by a mean of 42.8 mmol/L at week 24 in the active intervention group.'\n"
                "State what this implies for channel opening. State it directly as a fact without hedging words."
            ),
            "expected_evidence_numbers": []
        }
    ]

    results = []

    for scenario in test_scenarios:
        sid = scenario["id"]
        source_id = scenario["source_id"]
        prompt = scenario["prompt"]
        source_text = sources[source_id]

        worker_reply = query_local_worker(prompt)

        # Parse claim object from worker reply or construct one for evaluation
        claim_obj = {
            "id": f"claim-{sid.lower()}",
            "source": source_id,
            "kind": "observed",
        }

        # Check if worker returned JSON
        parsed_json = None
        if "{" in worker_reply and "}" in worker_reply:
            try:
                start = worker_reply.find("{")
                end = worker_reply.rfind("}") + 1
                parsed_json = json.loads(worker_reply[start:end])
            except Exception:
                parsed_json = None

        if parsed_json and isinstance(parsed_json, dict) and "quote" in parsed_json:
            claim_obj["claim"] = parsed_json.get("claim", worker_reply)
            claim_obj["quote"] = parsed_json.get("quote", "")
            claim_obj["kind"] = parsed_json.get("kind", "observed")
        else:
            claim_obj["claim"] = worker_reply.strip('"\n ')
            # Test whether worker's generation itself constitutes an exact verbatim quote from evidence:
            claim_obj["quote"] = worker_reply.strip('"\n ')

        ctx["id_counts"][claim_obj["id"]] = 1

        # Evaluate WITH lab checks
        r1_ok, r1_msg = check_claims.rule_r1(claim_obj, ctx)
        rule_evals = {"R1": {"ok": r1_ok, "msg": r1_msg}}
        failed_rules = []
        if not r1_ok:
            failed_rules.append("R1")
        else:
            for r_name in ("R2", "R3", "R4", "R5"):
                ok, msg = check_claims.RULES[r_name](claim_obj, ctx)
                rule_evals[r_name] = {"ok": ok, "msg": msg}
                if not ok:
                    failed_rules.append(r_name)

        status_with_checks = "FAIL" if failed_rules else "PASS"

        # Evaluate WITHOUT lab checks (baseline unverified pipeline accepts any non-empty claim)
        status_without_checks = "PASS" if bool(claim_obj.get("claim")) else "FAIL"

        result_entry = {
            "scenario_id": sid,
            "title": scenario["title"],
            "model": DEFAULT_MODEL,
            "source_passage": source_text,
            "prompt": prompt,
            "worker_raw_reply": worker_reply,
            "evaluated_claim_object": claim_obj,
            "without_checks": {
                "verdict": status_without_checks,
                "notes": "Unconstrained baseline model accepts generation without checking exact substring quote or number matches"
            },
            "with_checks": {
                "verdict": status_with_checks,
                "failed_rules": failed_rules,
                "rule_evaluations": rule_evals,
                "notes": "Lab checkers enforce exact lexical containment, number grounding, and scope requirements"
            }
        }
        results.append(result_entry)

    summary = {
        "evaluation_title": "Local Worker (qwen3.8:27b-64k) Verification Analysis",
        "model": DEFAULT_MODEL,
        "scenarios_tested": len(results),
        "results": results,
    }

    if output_path:
        with open(output_path, "w", encoding="utf-8") as f:
            json.dump(summary, f, indent=2)

    return summary


if __name__ == "__main__":
    out = BENCHMARK_DIR / "local_worker_results.json"
    print(f"Running evaluation against local worker {DEFAULT_MODEL}...")
    res = run_worker_evaluation(out)
    print(f"Done. Evaluated {len(res['results'])} scenarios. Results saved to {out}")
