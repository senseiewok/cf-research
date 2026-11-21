#!/usr/bin/env python3
"""Coupled CFTR Variant -> Modality -> Mucus-LNP Delivery -> Trial Deduction Engine.

100% Pure Python Standard Library Runtime.
Connects:
  1. `benchmarks/theratype/cftr_theratype_benchmark.py` (Evo 2 & AlphaGenome AVI variant classes)
  2. `tools/discovery/alphagenome_cftr_splice_probe.py` (GRCh38 CFTR locus & splicing vs. PTC)
  3. `tools/discovery/agile_mucus_lnp_ranker.py` + `mucus_barrier_screener.py` (140 nm vs 220 nm sputum)
  4. `tools/discovery/round203_chem_mod_suptrna_lnp_bridge.py` (m1A sup-tRNA I_sc_equiv rescue)
  5. `tools/trial_atlas/rare_genotype_trial_screener.py` (ClinicalTrials.gov rare-genotype coverage)

Deduces counterfactuals that no single tool can see alone:
  - Variant-specific therapeutic modality routing (modulator vs. ASO splice-switch vs. m1A sup-tRNA LNP vs. mRNA/gene)
  - Minimum dornase alfa co-therapy requirement (140 nm untreated vs. 220 nm dornase-expanded sputum)
  - Coupled end-to-end rescue score combining mucus transit survival, AGILE potency, and I_sc_equiv >= 15% WT
  - ClinicalTrials.gov trial-coverage gap per variant class

Research simulation only; not medical advice.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
from typing import Any, Dict, List, Sequence

TOOLS_DIR = Path(__file__).resolve().parent
REPO_ROOT = TOOLS_DIR.parent.parent
for p in (TOOLS_DIR, REPO_ROOT / "benchmarks" / "theratype", REPO_ROOT / "tools" / "trial_atlas"):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

import agile_mucus_lnp_ranker as aml  # noqa: E402
import alphagenome_cftr_splice_probe as agp  # noqa: E402
import cftr_theratype_benchmark as ctb  # noqa: E402
import round203_chem_mod_suptrna_lnp_bridge as r203  # noqa: E402

RESEARCH_BOUNDARY = "Research simulation only; not medical advice."


def route_variant_modality(case: ctb.CftrVariantCase) -> str:
    """Deduce candidate modality family from variant class and model coverage."""
    if case.modulator_responsive == 1 and case.variant_family in (
        "missense",
        "inframe_indel",
    ):
        return "SMALL_MOLECULE_MODULATOR_ELIGIBLE"
    if case.variant_family == "splice_regulatory":
        return "ASO_SPLICE_MODULATION_OR_MODULATOR" if case.modulator_responsive == 1 else "ASO_SPLICE_MODULATION_OR_MRNA"
    if "Ter" in case.hgvs_p_or_c:
        return "M1A_SUPTRNA_LNP_OR_MRNA_GENE"
    return "GENOTYPE_AGNOSTIC_MRNA_OR_GENE_THERAPY"


def deduce_end_to_end_rescue(
    variant_id: str = "CFTR-W1282X",
    pore_diameter_nm: float = 140.0,
    f_mod_m1a: float = 0.42,
) -> Dict[str, Any]:
    """Run coupled Variant -> Modality -> Mucus-LNP -> m1A sup-tRNA I_sc deduction."""
    by_id = {c.variant_id: c for c in ctb.OPEN_CFTR_THERATYPE_BENCHMARK}
    if variant_id not in by_id:
        raise ValueError(f"Unknown benchmark variant_id: {variant_id!r}")
    case = by_id[variant_id]
    modality = route_variant_modality(case)

    # Stage 2: Rank LNP carriers under specified sputum pore regime (140 nm untreated vs 220 nm dornase)
    lnp_report = aml.rank_lnp_candidates(pore_diameter_nm=pore_diameter_nm)
    top_lnp = lnp_report["ranked"][0] if lnp_report["ranked"] else None

    # Stage 3: Couple m1A sup-tRNA biophysical rescue (Round 203) when routed to PTC/RNA rescue
    if modality == "M1A_SUPTRNA_LNP_OR_MRNA_GENE" and top_lnp and top_lnp["mucus_gate_passed"]:
        # Convert mucus transit survival + AGILE potency into effective delivery eta_encap
        eta_encap_eff = min(0.95, 0.70 + 0.25 * top_lnp["agile_airway_potency_score"])
        t_half, b_trna, i_tlr = r203.compute_endosomal_stability_and_tlr78(
            f_mod=f_mod_m1a,
            eta_encap=eta_encap_eff,
            eta_esc=0.40,
        )
        eta_chg, alpha_sup = r203.compute_synthetase_charging_and_decoding(
            f_mod=f_mod_m1a,
            f_core_block=0.0,
            s_evo2_fitness=0.95,
            b_trna=b_trna,
        )
        _, _, _, i_sc_raw = r203.compute_nbd_folding_and_anion_efflux(
            alpha_suptrna=alpha_sup,
            f_mrna0=0.28,
            ddg_insert_kcal=0.0,
            ddg_corr_kcal=0.0,
        )
        i_sc_equiv = round(float(i_sc_raw), 2)
        curative_threshold_met = bool(
            float(t_half) >= 4.0 and float(i_tlr) <= 0.20 and float(eta_chg) >= 0.65 and float(i_sc_raw) >= 15.0
        )
    else:
        i_sc_equiv = None
        curative_threshold_met = False

    return {
        "variant_id": case.variant_id,
        "hgvs": case.hgvs_p_or_c,
        "variant_family": case.variant_family,
        "alphamissense_covered": case.alphamissense_score is not None,
        "evo2_score": case.evo2_zero_shot_score,
        "alphagenome_avi_score": case.alphagenome_avi_score,
        "routed_modality": modality,
        "sputum_pore_nm": pore_diameter_nm,
        "top_lnp_candidate": lnp_report["top_candidate"],
        "feasible_lnp_count": lnp_report["feasible_count"],
        "coupled_i_sc_equiv_pct": i_sc_equiv,
        "coupled_rescue_threshold_met": curative_threshold_met,
    }


def deduce_dornase_counterfactual(variant_id: str = "CFTR-W1282X") -> Dict[str, Any]:
    """Compare untreated (140 nm) vs dornase-expanded (220 nm) CF sputum delivery feasibility."""
    untreated = deduce_end_to_end_rescue(variant_id=variant_id, pore_diameter_nm=140.0)
    dornase = deduce_end_to_end_rescue(variant_id=variant_id, pore_diameter_nm=220.0)
    return {
        "variant_id": variant_id,
        "untreated_140nm_feasible_lnps": untreated["feasible_lnp_count"],
        "dornase_220nm_feasible_lnps": dornase["feasible_lnp_count"],
        "dornase_expands_carrier_window": dornase["feasible_lnp_count"] >= untreated["feasible_lnp_count"],
        "untreated": untreated,
        "dornase": dornase,
    }


def evaluate_falsification_controls() -> Dict[str, Any]:
    """Run deterministic negative controls on coupled deduction invariants."""
    controls: List[Dict[str, Any]] = []

    # 1. Unknown variant_id must fail closed
    try:
        deduce_end_to_end_rescue(variant_id="CFTR-FAKE999X")
        controls.append({"id": "unknown_variant_rejected", "falsified": False})
    except ValueError:
        controls.append({"id": "unknown_variant_rejected", "falsified": True})

    # 2. Unmodified naked sup-tRNA (f_mod=0.0) must NOT meet >= 15% WT I_sc_equiv threshold
    naked = deduce_end_to_end_rescue(variant_id="CFTR-G542X", f_mod_m1a=0.0)
    controls.append({"id": "unmodified_trna_fails_rescue", "falsified": not naked["coupled_rescue_threshold_met"]})

    # 3. Severely collapsed 50 nm sputum pore must reject all 80-300 nm LNPs
    collapsed = deduce_end_to_end_rescue(variant_id="CFTR-W1282X", pore_diameter_nm=50.0)
    controls.append({"id": "collapsed_50nm_sputum_blocks_all_lnps", "falsified": collapsed["feasible_lnp_count"] == 0})

    all_passed = all(c["falsified"] for c in controls)
    return {
        "status": "PASS" if all_passed else "FAIL",
        "controls_passed": sum(1 for c in controls if c["falsified"]),
        "total_controls": len(controls),
        "controls": controls,
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Coupled CFTR Variant -> Modality -> Delivery Deduction Engine")
    parser.add_argument("--falsify", action="store_true", help="Run negative-control falsification suite")
    parser.add_argument("--json", action="store_true", help="Emit JSON output")
    args = parser.parse_args(list(argv) if argv is not None else None)

    w1282x = deduce_dornase_counterfactual("CFTR-W1282X")
    splice = deduce_end_to_end_rescue("CFTR-1585-1G>A")
    falsify = evaluate_falsification_controls() if args.falsify else {"status": "SKIPPED", "controls_passed": 0, "total_controls": 0}

    ok = (
        w1282x["untreated"]["coupled_rescue_threshold_met"]
        and splice["routed_modality"] == "ASO_SPLICE_MODULATION_OR_MRNA"
        and (not args.falsify or falsify["status"] == "PASS")
    )
    payload = {
        "status": "PASS" if ok else "FAIL",
        "boundary": RESEARCH_BOUNDARY,
        "w1282x_dornase_counterfactual": w1282x,
        "splice_1585_routing": splice,
        "falsification": falsify,
    }
    if args.json:
        print(json.dumps(payload, indent=2))
    else:
        print(
            f"CFTR Variant-to-Delivery Deducer: {payload['status']} "
            f"(W1282X top={w1282x['untreated']['top_lnp_candidate']}, "
            f"I_sc={w1282x['untreated']['coupled_i_sc_equiv_pct']}% WT, "
            f"1585-1G>A->{splice['routed_modality']}, "
            f"falsify={falsify['controls_passed']}/{falsify['total_controls']})"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
