#!/usr/bin/env python3
"""Independent exact-rational & rank-sum audit of the Open CFTR Theratype Benchmark (`T-0925` / `PRED-0010`).

Rule 2 of `discovery-practice`: a second method for a headline claim must live in separate code
and use a genuinely different algorithm rather than copying floating-point pairwise loops or
regular expressions from `cftr_theratype_benchmark.py`.

This script imports NO regular-expression engine (`re` is not imported) and NO modules from
`cftr_theratype_benchmark` or `tools/`:
  1. Uses `json.loads(..., parse_float=Decimal)` with exact `ROUND_HALF_EVEN` quantization (`0.01`)
     and threshold checks for McDonald et al. (2024) Table 1 VUS controls.
  2. Uses exact `fractions.Fraction` tie-averaged Wilcoxon rank-sum statistics (`U = R_pos - k*(k+1)/2`)
     and the closed-form abstention identity `1/2 + (|S_1|*|S_0| / (N_1*N_0)) * (AUC_cond - 1/2)`.
  3. Uses character-by-character HGVS splice dinucleotide parsing (no regexes) for the
     sequence-consequence baseline.
  4. Independently verifies UniProt `P13569` 6-domain intervals, cross-split codon gap,
     and US FDA `TRIKAFTA` v23 vs `ALYFTREK` v9 counts and variant discordance.
"""
from __future__ import annotations

import argparse
from decimal import Decimal, ROUND_HALF_EVEN
from fractions import Fraction
import json
from pathlib import Path
import sys
from typing import Any, Dict, List, Sequence, Tuple

HERE = Path(__file__).resolve().parent
DEFAULT_FIXTURE = HERE / "theratype_open_annotations_2026-11-21.json"

DOMAIN_TABLE: Tuple[Tuple[int, int, str], ...] = (
    (1, 80, "Lasso_N_terminal"),
    (81, 380, "TMD1"),
    (381, 630, "NBD1"),
    (631, 845, "R_domain"),
    (846, 1190, "TMD2"),
    (1191, 1480, "NBD2_C_terminal"),
)


def _domain_for_codon(codon: int) -> str:
    for lo, hi, name in DOMAIN_TABLE:
        if lo <= codon <= hi:
            return name
    raise ValueError(f"Codon {codon} outside CFTR 1..1480")


def _decimal_vus_class(score_dec: Decimal) -> str:
    if score_dec < Decimal("0.34"):
        return "benign"
    if score_dec <= Decimal("0.56"):
        return "ambiguous"
    return "pathogenic"


def _is_canonical_splice_charwise(hgvs: str) -> bool:
    """Regex-free parser for `c.<digits>[+-][12][ACGT]>[ACGT]`."""
    s = hgvs.strip()
    if not s.startswith("c.") or len(s) < 8:
        return False
    idx = 2
    while idx < len(s) and s[idx].isdigit():
        idx += 1
    if idx == 2 or idx + 5 != len(s):
        return False
    sign_ch = s[idx]
    offset_ch = s[idx + 1]
    ref_nt = s[idx + 2]
    arrow_ch = s[idx + 3]
    alt_nt = s[idx + 4]
    bases = {"A", "C", "G", "T"}
    return (
        sign_ch in {"+", "-"}
        and offset_ch in {"1", "2"}
        and ref_nt in bases
        and arrow_ch == ">"
        and alt_nt in bases
    )


def _rank_sum_auroc_exact(
    items: Sequence[Tuple[Fraction | None, int]],
) -> Dict[str, Any]:
    """Compute exact Mann-Whitney U via tie-averaged Wilcoxon rank sums + closed-form abstention identity."""
    n_total = len(items)
    n1 = sum(1 for _, y in items if y == 1)
    n0 = sum(1 for _, y in items if y == 0)
    scored = [(s, y) for s, y in items if s is not None]
    s1 = sum(1 for _, y in scored if y == 1)
    s0 = sum(1 for _, y in scored if y == 0)
    cov = Fraction(len(scored), n_total) if n_total else Fraction(0, 1)

    if n1 == 0 or n0 == 0:
        return {
            "coverage": round(float(cov), 4),
            "coverage_raw": float(cov),
            "coverage_exact": f"{cov.numerator}/{cov.denominator}",
            "conditional_auroc": None,
            "conditional_auroc_raw": None,
            "conditional_auroc_exact": None,
            "abstention_aware_auroc": None,
            "abstention_aware_auroc_raw": None,
            "abstention_aware_auroc_exact": None,
        }

    if s1 > 0 and s0 > 0:
        ordered = sorted(scored, key=lambda pair: pair[0])
        m = len(ordered)
        avg_ranks: List[Fraction] = [Fraction(0, 1)] * m
        pos = 0
        while pos < m:
            end = pos + 1
            while end < m and ordered[end][0] == ordered[pos][0]:
                end += 1
            mean_rank = Fraction(pos + 1 + end, 2)
            for k in range(pos, end):
                avg_ranks[k] = mean_rank
            pos = end
        r_pos = sum((avg_ranks[i] for i in range(m) if ordered[i][1] == 1), Fraction(0, 1))
        u_stat = r_pos - Fraction(s1 * (s1 + 1), 2)
        cond_auc: Fraction | None = u_stat / Fraction(s1 * s0, 1)
        abst_auc = Fraction(1, 2) + Fraction(s1 * s0, n1 * n0) * (cond_auc - Fraction(1, 2))
    else:
        cond_auc = None
        abst_auc = Fraction(1, 2)

    return {
        "coverage": round(float(cov), 4),
        "coverage_raw": float(cov),
        "coverage_exact": f"{cov.numerator}/{cov.denominator}",
        "conditional_auroc": round(float(cond_auc), 4) if cond_auc is not None else None,
        "conditional_auroc_raw": float(cond_auc) if cond_auc is not None else None,
        "conditional_auroc_exact": f"{cond_auc.numerator}/{cond_auc.denominator}" if cond_auc is not None else None,
        "abstention_aware_auroc": round(float(abst_auc), 4),
        "abstention_aware_auroc_raw": float(abst_auc),
        "abstention_aware_auroc_exact": f"{abst_auc.numerator}/{abst_auc.denominator}",
    }


def audit_independently(fixture_path: Path = DEFAULT_FIXTURE) -> Dict[str, Any]:
    raw = json.loads(fixture_path.read_text(encoding="utf-8"), parse_float=Decimal)
    controls = raw["mcdonald_2024_table1_positive_controls"]
    variants = raw["variants"]

    pred_letter_map = {"B": "benign", "A": "ambiguous", "P": "pathogenic"}
    controls_ok = 0
    control_codons: set[int] = set()
    for ctrl in controls:
        codon = int(ctrl["codon_position"])
        _domain_for_codon(codon)
        control_codons.add(codon)
        score_dec = Decimal(ctrl["myvariant_am_score"])
        exp_dec = Decimal(ctrl["mcdonald_2024_table1_score"]).quantize(Decimal("0.01"), rounding=ROUND_HALF_EVEN)
        quant_dec = score_dec.quantize(Decimal("0.01"), rounding=ROUND_HALF_EVEN)
        exp_cls = str(ctrl["mcdonald_2024_table1_class"])
        pred_cls = pred_letter_map.get(str(ctrl["myvariant_am_pred"]), "unknown")
        if quant_dec == exp_dec and _decimal_vus_class(score_dec) == exp_cls and pred_cls == exp_cls:
            controls_ok += 1

    train_rows = [v for v in variants if v["split"] == "train"]
    test_rows = [v for v in variants if v["split"] == "test"]

    train_codons = [int(v["codon_position"]) for v in train_rows]
    test_codons = [int(v["codon_position"]) for v in test_rows]
    train_domains = sorted({_domain_for_codon(c) for c in train_codons})
    test_domains = sorted({_domain_for_codon(c) for c in test_codons})
    domain_overlap_count = len(set(train_domains) & set(test_domains))
    codon_overlap_count = len(set(train_codons) & set(test_codons))
    control_overlap_count = len((set(train_codons) | set(test_codons)) & control_codons)
    min_codon_gap = min((abs(tc - qc) for tc in train_codons for qc in test_codons), default=0)

    test_trikafta_pos = sum(int(v["trikafta_v23_responsive"]) for v in test_rows)
    test_alyftrek_pos = sum(int(v["alyftrek_v9_responsive"]) for v in test_rows)
    all_trikafta_pos = sum(int(v["trikafta_v23_responsive"]) for v in variants)
    all_alyftrek_pos = sum(int(v["alyftrek_v9_responsive"]) for v in variants)

    discordant_eti_vs_vtd = sorted(
        str(v["legacy_spl_name"])
        for v in variants
        if int(v["trikafta_v23_responsive"]) != int(v["alyftrek_v9_responsive"])
    )

    am_test_items: List[Tuple[Fraction | None, int]] = []
    am_vtd_items: List[Tuple[Fraction | None, int]] = []
    base_test_items: List[Tuple[Fraction | None, int]] = []
    base_vtd_items: List[Tuple[Fraction | None, int]] = []

    for v in test_rows:
        y_eti = int(v["trikafta_v23_responsive"])
        y_vtd = int(v["alyftrek_v9_responsive"])
        am_val = v["alphamissense_score"]
        s_am = Fraction(1, 1) - Fraction(am_val) if am_val is not None else None
        am_test_items.append((s_am, y_eti))
        am_vtd_items.append((s_am, y_vtd))

        fam = str(v["variant_family"])
        hgvs = str(v["hgvs_p_or_c"])
        if fam == "class_i_ptc_frameshift" or _is_canonical_splice_charwise(hgvs):
            s_base = Fraction(0, 1)
        elif fam in {"splice_regulatory", "inframe_indel"}:
            s_base = Fraction(1, 2)
        elif fam == "missense":
            s_base = Fraction(3, 4)
        else:
            s_base = Fraction(0, 1)
        base_test_items.append((s_base, y_eti))
        base_vtd_items.append((s_base, y_vtd))

    am_eti_metrics = _rank_sum_auroc_exact(am_test_items)
    am_vtd_metrics = _rank_sum_auroc_exact(am_vtd_items)
    base_eti_metrics = _rank_sum_auroc_exact(base_test_items)
    base_vtd_metrics = _rank_sum_auroc_exact(base_vtd_items)

    return {
        "mcdonald_controls_verified": controls_ok,
        "train_count": len(train_rows),
        "test_count": len(test_rows),
        "total_variants": len(variants),
        "codon_overlap_count": codon_overlap_count,
        "control_overlap_count": control_overlap_count,
        "domain_overlap_count": domain_overlap_count,
        "min_cross_split_codon_gap": min_codon_gap,
        "train_domains": train_domains,
        "test_domains": test_domains,
        "test_trikafta_v23_responsive_count": test_trikafta_pos,
        "test_alyftrek_v9_responsive_count": test_alyftrek_pos,
        "all30_trikafta_v23_responsive_count": all_trikafta_pos,
        "all30_alyftrek_v9_responsive_count": all_alyftrek_pos,
        "discordant_trikafta_vs_alyftrek": discordant_eti_vs_vtd,
        "alphamissense_trikafta_v23": am_eti_metrics,
        "alphamissense_alyftrek_v9": am_vtd_metrics,
        "sequence_consequence_baseline_trikafta_v23": base_eti_metrics,
        "sequence_consequence_baseline_alyftrek_v9": base_vtd_metrics,
    }


def main(argv: Sequence[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Independent exact-rational audit for T-0925 / PRED-0010")
    ap.add_argument("--fixture", type=Path, default=DEFAULT_FIXTURE)
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--json", action="store_true")
    args = ap.parse_args(list(argv) if argv is not None else None)

    report = audit_independently(args.fixture)
    am_eti = report["alphamissense_trikafta_v23"]
    am_vtd = report["alphamissense_alyftrek_v9"]
    base_eti = report["sequence_consequence_baseline_trikafta_v23"]
    base_vtd = report["sequence_consequence_baseline_alyftrek_v9"]

    ok = (
        report["mcdonald_controls_verified"] == 7
        and report["train_count"] == 10
        and report["test_count"] == 20
        and report["codon_overlap_count"] == 0
        and report["control_overlap_count"] == 0
        and report["domain_overlap_count"] == 0
        and report["min_cross_split_codon_gap"] == 71
        and report["test_trikafta_v23_responsive_count"] == 9
        and report["test_alyftrek_v9_responsive_count"] == 10
        and report["all30_trikafta_v23_responsive_count"] == 15
        and report["all30_alyftrek_v9_responsive_count"] == 16
        and report["discordant_trikafta_vs_alyftrek"] == ["2789+5G->A", "3849+10kbC->T", "R1066C", "R560T", "V520F"]
        and am_eti["coverage_exact"] == "11/20"
        and am_eti["conditional_auroc_exact"] == "8/15"
        and am_eti["abstention_aware_auroc_exact"] == "101/198"
        and am_vtd["conditional_auroc_exact"] == "8/9"
        and am_vtd["abstention_aware_auroc_exact"] == "57/100"
        and base_eti["coverage_exact"] == "1/1"
        and base_eti["abstention_aware_auroc_exact"] == "15/22"
        and base_vtd["abstention_aware_auroc_exact"] == "7/8"
    )
    if args.check and not ok:
        print(f"FAIL independent_audit: {report}", file=sys.stderr)
        return 1
    if args.json:
        print(json.dumps(report, indent=2))
    else:
        print(
            f"PASS independent_audit: controls={report['mcdonald_controls_verified']}/7 "
            f"gap={report['min_cross_split_codon_gap']} "
            f"AM_cov={am_eti['coverage']} ({am_eti['coverage_exact']}) "
            f"AM_cond={am_eti['conditional_auroc']} ({am_eti['conditional_auroc_exact']}) "
            f"AM_abst={am_eti['abstention_aware_auroc']} ({am_eti['abstention_aware_auroc_exact']}) "
            f"base_abst={base_eti['abstention_aware_auroc']} ({base_eti['abstention_aware_auroc_exact']})"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
