#!/usr/bin/env python3
"""Two-Stage AGILE + Mucus Barrier LNP Candidate Ranker for sup-tRNA & mRNA Delivery.

100% Pure Python Standard Library Runtime.
Stage 1: Physical CF Sputum Mucus Barrier Gate via `mucus_barrier_screener.py` (T-0941 / PRED-0004):
  Rejects any LNP formulation whose transit time exceeds the 900 s mucociliary clearance window
  in 140 nm untreated CF sputum.
Stage 2: Joint AGILE-style ionizable lipid airway transfection score + Round 203 m1A sup-tRNA
  bioavailability ranking (`DOI: 10.1126/science.aeb0054`, `DOI: 10.1021/acsnano.5c15613`).

Research simulation only; not medical advice.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
from pathlib import Path
import sys
from typing import Any, Dict, List, Sequence

TOOLS_DIR = Path(__file__).resolve().parent
if str(TOOLS_DIR) not in sys.path:
    sys.path.insert(0, str(TOOLS_DIR))

import mucus_barrier_screener as mbs  # noqa: E402

RESEARCH_BOUNDARY = "Research simulation only; not medical advice."


@dataclass(frozen=True)
class LnpCandidateSpec:
    candidate_id: str
    diameter_nm: float
    zeta_mv: float
    peg_density_chains_nm2: float
    mucus_thickness_um: float
    agile_airway_potency_score: float  # Normalized [0, 1] AGILE-predicted airway transfection potency
    cargo_modality: str                # "m1a_suptrna" | "cftr_mrna" | "control"


DEFAULT_LNP_CANDIDATES: Tuple[LnpCandidateSpec, ...] = (
    LnpCandidateSpec("HLPN-PEG-80nm-m1A-suptRNA", 80.0, -2.0, 0.28, 15.0, 0.88, "m1a_suptrna"),
    LnpCandidateSpec("Hybrid-LNP-100nm-mRNA", 100.0, -3.0, 0.28, 15.0, 0.82, "cftr_mrna"),
    LnpCandidateSpec("Cationic-Lipoplex-110nm-Trap", 110.0, 28.0, 0.0, 20.0, 0.95, "control"),
    LnpCandidateSpec("Oversized-LNP-300nm-Steric", 300.0, -4.0, 0.25, 20.0, 0.79, "control"),
)


def rank_lnp_candidates(
    candidates: Sequence[LnpCandidateSpec] = DEFAULT_LNP_CANDIDATES,
    pore_diameter_nm: float = 140.0,
) -> Dict[str, Any]:
    """Apply Stage 1 SEOB mucus barrier gate and Stage 2 AGILE potency ranking."""
    evaluated: List[Dict[str, Any]] = []
    for cand in candidates:
        if not (0.0 <= cand.agile_airway_potency_score <= 1.0):
            raise ValueError(f"agile_airway_potency_score must be in [0, 1], got {cand.agile_airway_potency_score}")
        screen = mbs.evaluate_mucus_barrier(
            mbs.VectorBiophysicalInput(
                diameter_nm=cand.diameter_nm,
                zeta_mv=cand.zeta_mv,
                peg_density_chains_nm2=cand.peg_density_chains_nm2,
                pore_diameter_nm=pore_diameter_nm,
                mucus_thickness_um=cand.mucus_thickness_um,
            )
        )
        metrics = screen["computed_metrics"]
        verdict = metrics["crosses_before_mcc_window"]["verdict"]
        d_eff_over_dw = metrics["d_eff_over_d_water"]["value"]
        transit_time_s = metrics["mean_transit_time_s"]["value"]
        mucus_pass = verdict == "FEASIBLE_MUCUS_PENETRATING_CANDIDATE"
        composite_score = round(d_eff_over_dw * cand.agile_airway_potency_score * 100.0, 3) if mucus_pass else 0.0
        evaluated.append(
            {
                "candidate_id": cand.candidate_id,
                "cargo_modality": cand.cargo_modality,
                "mucus_gate_passed": mucus_pass,
                "d_eff_over_dw": d_eff_over_dw,
                "transit_time_s": transit_time_s,
                "agile_airway_potency_score": cand.agile_airway_potency_score,
                "composite_rank_score": composite_score,
            }
        )

    evaluated.sort(key=lambda row: (-int(row["mucus_gate_passed"]), -row["composite_rank_score"]))
    feasible = [r for r in evaluated if r["mucus_gate_passed"]]
    rejected = [r for r in evaluated if not r["mucus_gate_passed"]]
    return {
        "total_candidates": len(evaluated),
        "feasible_count": len(feasible),
        "rejected_count": len(rejected),
        "top_candidate": feasible[0]["candidate_id"] if feasible else None,
        "ranked": evaluated,
    }


def evaluate_falsification_controls() -> Dict[str, Any]:
    """Verify high-potency cationic & oversized traps are rejected by the Stage 1 mucus gate."""
    controls: List[Dict[str, Any]] = []
    rep = rank_lnp_candidates()
    by_id = {r["candidate_id"]: r for r in rep["ranked"]}

    # 1. Cationic lipoplex has highest raw AGILE potency (0.95) but MUST be rejected by Stage 1 mucus gate
    controls.append(
        {"id": "cationic_high_potency_trap_rejected", "falsified": not by_id["Cationic-Lipoplex-110nm-Trap"]["mucus_gate_passed"]}
    )
    # 2. Oversized 300 nm LNP MUST be rejected by steric sieving
    controls.append(
        {"id": "oversized_300nm_steric_rejected", "falsified": not by_id["Oversized-LNP-300nm-Steric"]["mucus_gate_passed"]}
    )
    # 3. Out-of-range AGILE score must raise ValueError
    try:
        rank_lnp_candidates([LnpCandidateSpec("BAD", 80.0, -2.0, 0.9, 15.0, 1.5, "control")])
        controls.append({"id": "invalid_agile_score_rejected", "falsified": False})
    except ValueError:
        controls.append({"id": "invalid_agile_score_rejected", "falsified": True})

    all_passed = all(c["falsified"] for c in controls)
    return {
        "status": "PASS" if all_passed else "FAIL",
        "controls_passed": sum(1 for c in controls if c["falsified"]),
        "total_controls": len(controls),
        "controls": controls,
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Two-Stage AGILE + Mucus Barrier LNP Candidate Ranker")
    parser.add_argument("--falsify", action="store_true", help="Run negative-control falsification suite")
    parser.add_argument("--json", action="store_true", help="Emit JSON output")
    args = parser.parse_args(list(argv) if argv is not None else None)

    ranked = rank_lnp_candidates()
    falsify = evaluate_falsification_controls() if args.falsify else {"status": "SKIPPED", "controls_passed": 0, "total_controls": 0}
    ok = (
        ranked["feasible_count"] == 2
        and ranked["rejected_count"] == 2
        and ranked["top_candidate"] == "HLPN-PEG-80nm-m1A-suptRNA"
        and (not args.falsify or falsify["status"] == "PASS")
    )
    payload = {"status": "PASS" if ok else "FAIL", "boundary": RESEARCH_BOUNDARY, "ranking": ranked, "falsification": falsify}
    if args.json:
        print(json.dumps(payload, indent=2))
    else:
        print(
            f"AGILE + Mucus LNP Ranker: {payload['status']} "
            f"(feasible={ranked['feasible_count']}/{ranked['total_candidates']}, "
            f"top={ranked['top_candidate']}, falsify={falsify['controls_passed']}/{falsify['total_controls']})"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
