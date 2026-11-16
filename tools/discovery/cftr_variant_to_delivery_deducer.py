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
  - Cargo-matched mucus-gate feasibility and dornase alfa co-therapy expansion (140 nm vs. 220 nm sputum)
  - Variant-specific NMD/NBD2 kinetics (exon 23 W1282X vs. internal exon 12 G542X) coupled to LNP potency
  - ClinicalTrials.gov recruiting trial-coverage count per variant class

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
import rare_genotype_trial_screener as rgts  # noqa: E402
import round203_chem_mod_suptrna_lnp_bridge as r203  # noqa: E402

RESEARCH_BOUNDARY = "Research simulation only; not medical advice."


def _is_nonsense_ptc(case: ctb.CftrVariantCase) -> bool:
    """Detect premature termination codon variants across Ter/* /X HGVS forms."""
    if case.variant_family != "class_i_ptc_frameshift":
        return False
    return (
        "Ter" in case.hgvs_p_or_c
        or "*" in case.hgvs_p_or_c
        or case.variant_id.endswith("X")
    )


def route_variant_modality(case: ctb.CftrVariantCase) -> str:
    """Deduce candidate modality family from variant class and model coverage."""
    if case.modulator_responsive == 1 and case.variant_family in (
        "inframe_folding",
        "gating_missense",
        "gating_missense_heldout",
    ):
        return "SMALL_MOLECULE_MODULATOR_ELIGIBLE"
    if case.variant_family == "splice_regulatory":
        return "ASO_SPLICE_MODULATION_OR_MODULATOR" if case.modulator_responsive == 1 else "ASO_SPLICE_MODULATION_OR_MRNA"
    if _is_nonsense_ptc(case):
        return "M1A_SUPTRNA_LNP_OR_MRNA_GENE"
    return "GENOTYPE_AGNOSTIC_MRNA_OR_GENE_THERAPY"


def _expected_cargo_modality(routed_modality: str) -> str | None:
    """Map routed therapeutic modality to agile_mucus_lnp_ranker cargo_modality."""
    if routed_modality == "M1A_SUPTRNA_LNP_OR_MRNA_GENE":
        return "m1a_suptrna"
    if routed_modality in ("GENOTYPE_AGNOSTIC_MRNA_OR_GENE_THERAPY", "ASO_SPLICE_MODULATION_OR_MRNA"):
        return "cftr_mrna"
    return None


def _lookup_grch38_locus_context(case: ctb.CftrVariantCase) -> Dict[str, Any] | None:
    """Match benchmark variant against GRCh38 AlphaGenome splice/PTC panel."""
    short_id = case.variant_id.replace("CFTR-", "").lower()
    for item in agp.CFTR_SPLICE_AND_PTC_PANEL:
        item_norm = item.label.replace(" ", "").lower()
        if short_id == item_norm or item.hgvs_c.lower() == case.hgvs_p_or_c.lower():
            agp.validate_cftr_grch38_variant(item.grch38_variant)
            return {
                "grch38_variant": item.grch38_variant,
                "expected_modality": item.expected_modality,
                "context": item.aso_or_readthrough_context,
            }
    return None


_CACHED_PRISMA_STUDIES: List[Dict[str, Any]] | None = None


def _get_prisma_studies() -> List[Dict[str, Any]]:
    global _CACHED_PRISMA_STUDIES
    if _CACHED_PRISMA_STUDIES is None:
        _, _CACHED_PRISMA_STUDIES = rgts.load_prisma_included_studies(rgts.DEFAULT_SNAPSHOT_JSON)
    return _CACHED_PRISMA_STUDIES


def _lookup_trial_coverage(case: ctb.CftrVariantCase, routed_modality: str) -> Dict[str, Any]:
    """Screen PRISMA-included ClinicalTrials.gov snapshot for matching variant-class trials."""
    if not rgts.DEFAULT_SNAPSHOT_JSON.exists():
        return {"recruiting_trial_count": 0, "sample_nct_ids": [], "trial_class_queried": "snapshot_unavailable"}

    if _is_nonsense_ptc(case):
        trial_class = "class_i_stop_nonsense"
    elif case.variant_family == "splice_regulatory":
        trial_class = "splice_modulation"
    elif case.variant_family == "class_i_ptc_frameshift" or routed_modality == "GENOTYPE_AGNOSTIC_MRNA_OR_GENE_THERAPY":
        trial_class = "genotype_agnostic_mrna_gene"
    else:
        trial_class = "all"

    studies = _get_prisma_studies()
    rows = rgts.screen_rare_genotype_trials(
        studies,
        variant_class=trial_class,
        status="RECRUITING",
        limit=50,
    )
    return {
        "trial_class_queried": trial_class,
        "recruiting_trial_count": len(rows),
        "sample_nct_ids": [r["nct_id"] for r in rows[:3]],
    }


def _variant_ptc_kinetics(case: ctb.CftrVariantCase) -> Dict[str, float]:
    """Select variant-specific NMD baseline, NBD2 corrector delta, and Evo2 fitness."""
    # Exon 23 / 3'-terminal PTC (e.g. W1282X) partially escapes NMD but destabilizes NBD2
    if case.codon_position >= 1200:
        return {
            "f_mrna0": 0.32,
            "ddg_insert_kcal": 0.0,
            "ddg_corr_kcal": 0.80,
            "s_evo2_fitness": 0.95,
        }
    # Internal PTCs (e.g. exon 12 G542X, R553X) undergo stronger NMD without NBD2 cavity penalty
    return {
        "f_mrna0": 0.22,
        "ddg_insert_kcal": 0.0,
        "ddg_corr_kcal": 0.0,
        "s_evo2_fitness": min(0.96, max(0.85, case.evo2_zero_shot_score)),
    }


def deduce_end_to_end_rescue(
    variant_id: str = "CFTR-W1282X",
    pore_diameter_nm: float = 140.0,
    f_mod_m1a: float = 0.42,
) -> Dict[str, Any]:
    """Run coupled Variant -> Modality -> Mucus-LNP -> m1A sup-tRNA I_sc & Trial deduction."""
    if pore_diameter_nm <= 0.0:
        raise ValueError(f"pore_diameter_nm must be positive, got {pore_diameter_nm}")
    if not (0.0 <= f_mod_m1a <= 1.0):
        raise ValueError(f"f_mod_m1a must be in [0.0, 1.0], got {f_mod_m1a}")

    by_id = {c.variant_id: c for c in ctb.OPEN_CFTR_THERATYPE_BENCHMARK}
    if variant_id not in by_id:
        raise ValueError(f"Unknown benchmark variant_id: {variant_id!r}")
    case = by_id[variant_id]
    modality = route_variant_modality(case)
    target_cargo = _expected_cargo_modality(modality)

    # Stage 2: Rank LNP carriers under specified sputum pore regime & filter by mucus gate + cargo
    lnp_report = aml.rank_lnp_candidates(pore_diameter_nm=pore_diameter_nm)
    feasible_all = [c for c in lnp_report["ranked"] if c["mucus_gate_passed"]]
    feasible_cargo = (
        [c for c in feasible_all if c["cargo_modality"] == target_cargo]
        if target_cargo is not None
        else feasible_all
    )
    top_feasible_lnp = feasible_cargo[0] if feasible_cargo else None

    # Stage 3: Couple variant-specific m1A sup-tRNA kinetics (Round 203) for PTC rescue
    kinetics = _variant_ptc_kinetics(case) if modality == "M1A_SUPTRNA_LNP_OR_MRNA_GENE" else None
    if modality == "M1A_SUPTRNA_LNP_OR_MRNA_GENE" and top_feasible_lnp is not None and kinetics is not None:
        # Keep formulation encapsulation anchored (0.94); scale endosomal escape by AGILE airway potency
        eta_encap_formulation = 0.94
        eta_esc_eff = min(0.45, 0.20 + 0.25 * top_feasible_lnp["agile_airway_potency_score"])
        t_half, b_trna, i_tlr = r203.compute_endosomal_stability_and_tlr78(
            f_mod=f_mod_m1a,
            eta_encap=eta_encap_formulation,
            eta_esc=eta_esc_eff,
        )
        eta_chg, alpha_sup = r203.compute_synthetase_charging_and_decoding(
            f_mod=f_mod_m1a,
            f_core_block=0.0,
            s_evo2_fitness=kinetics["s_evo2_fitness"],
            b_trna=b_trna,
        )
        _, _, _, i_sc_raw = r203.compute_nbd_folding_and_anion_efflux(
            alpha_suptrna=alpha_sup,
            f_mrna0=kinetics["f_mrna0"],
            ddg_insert_kcal=kinetics["ddg_insert_kcal"],
            ddg_corr_kcal=kinetics["ddg_corr_kcal"],
        )
        i_sc_equiv = round(float(i_sc_raw), 2)
        curative_threshold_met = bool(
            float(t_half) >= 4.0 and float(i_tlr) <= 0.20 and float(eta_chg) >= 0.65 and float(i_sc_raw) >= 15.0
        )
    else:
        i_sc_equiv = None
        curative_threshold_met = False

    locus_ctx = _lookup_grch38_locus_context(case)
    trial_cov = _lookup_trial_coverage(case, modality)

    return {
        "variant_id": case.variant_id,
        "hgvs": case.hgvs_p_or_c,
        "codon_position": case.codon_position,
        "variant_family": case.variant_family,
        "alphamissense_covered": case.alphamissense_score is not None,
        "evo2_score": case.evo2_zero_shot_score,
        "alphagenome_avi_score": case.alphagenome_avi_score,
        "grch38_locus_context": locus_ctx,
        "routed_modality": modality,
        "target_cargo_modality": target_cargo,
        "sputum_pore_nm": pore_diameter_nm,
        "top_lnp_candidate": top_feasible_lnp["candidate_id"] if top_feasible_lnp else None,
        "feasible_lnp_count": len(feasible_all),
        "cargo_matched_feasible_count": len(feasible_cargo),
        "f_mrna0_used": kinetics["f_mrna0"] if kinetics else None,
        "ddg_corr_kcal_used": kinetics["ddg_corr_kcal"] if kinetics else None,
        "coupled_i_sc_equiv_pct": i_sc_equiv,
        "coupled_rescue_threshold_met": curative_threshold_met,
        "trial_coverage": trial_cov,
    }


def deduce_dornase_counterfactual(variant_id: str = "CFTR-W1282X") -> Dict[str, Any]:
    """Compare untreated (140 nm) vs dornase-expanded (220 nm) CF sputum delivery feasibility."""
    untreated = deduce_end_to_end_rescue(variant_id=variant_id, pore_diameter_nm=140.0)
    dornase = deduce_end_to_end_rescue(variant_id=variant_id, pore_diameter_nm=220.0)
    delta = dornase["feasible_lnp_count"] - untreated["feasible_lnp_count"]
    return {
        "variant_id": variant_id,
        "untreated_140nm_feasible_lnps": untreated["feasible_lnp_count"],
        "dornase_220nm_feasible_lnps": dornase["feasible_lnp_count"],
        "carrier_window_expansion_delta": delta,
        "dornase_expands_carrier_window": delta > 0,
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
    controls.append(
        {
            "id": "collapsed_50nm_sputum_blocks_all_lnps",
            "falsified": collapsed["feasible_lnp_count"] == 0 and collapsed["top_lnp_candidate"] is None,
        }
    )

    # 4. W1282X (exon 23) vs G542X (exon 12) must exhibit distinct NMD/folding kinetics & I_sc
    w1282x = deduce_end_to_end_rescue(variant_id="CFTR-W1282X")
    g542x = deduce_end_to_end_rescue(variant_id="CFTR-G542X")
    controls.append(
        {
            "id": "w1282x_vs_g542x_kinetics_differ",
            "falsified": (
                w1282x["f_mrna0_used"] != g542x["f_mrna0_used"]
                and w1282x["coupled_i_sc_equiv_pct"] != g542x["coupled_i_sc_equiv_pct"]
            ),
        }
    )

    all_passed = all(c["falsified"] for c in controls)
    return {
        "status": "PASS" if all_passed else "FAIL",
        "controls_passed": sum(1 for c in controls if c["falsified"]),
        "total_controls": len(controls),
        "controls": controls,
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Coupled CFTR Variant -> Modality -> Delivery Deduction Engine")
    parser.add_argument("--variant", default="CFTR-W1282X", help="Benchmark variant_id (default: CFTR-W1282X)")
    parser.add_argument("--pore-diameter-nm", type=float, default=140.0, help="Sputum pore diameter in nm")
    parser.add_argument("--f-mod", type=float, default=0.42, help="m1A modification fraction in [0, 1]")
    parser.add_argument("--falsify", action="store_true", help="Run negative-control falsification suite")
    parser.add_argument("--json", action="store_true", help="Emit JSON output")
    args = parser.parse_args(list(argv) if argv is not None else None)

    primary = deduce_end_to_end_rescue(
        variant_id=args.variant,
        pore_diameter_nm=args.pore_diameter_nm,
        f_mod_m1a=args.f_mod,
    )
    counterfactual = deduce_dornase_counterfactual(args.variant)
    splice = deduce_end_to_end_rescue("CFTR-1585-1G>A")
    falsify = (
        evaluate_falsification_controls()
        if args.falsify
        else {"status": "SKIPPED", "controls_passed": 0, "total_controls": 0}
    )

    ok = not args.falsify or falsify["status"] == "PASS"
    payload = {
        "status": "PASS" if ok else "FAIL",
        "boundary": RESEARCH_BOUNDARY,
        "primary_deduction": primary,
        "dornase_counterfactual": counterfactual,
        "splice_1585_routing": splice,
        "falsification": falsify,
    }
    if args.json:
        print(json.dumps(payload, indent=2))
    else:
        print(
            f"CFTR Variant-to-Delivery Deducer: {payload['status']} "
            f"({args.variant} top={primary['top_lnp_candidate']}, "
            f"I_sc={primary['coupled_i_sc_equiv_pct']}% WT, "
            f"trials={primary['trial_coverage']['recruiting_trial_count']}, "
            f"1585-1G>A->{splice['routed_modality']}, "
            f"falsify={falsify['controls_passed']}/{falsify['total_controls']})"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
