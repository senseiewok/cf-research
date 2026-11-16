#!/usr/bin/env python3
"""AlphaGenome CFTR Splicing & Regulatory Variant Probe (GRCh38 / GENCODE v46).

100% Pure Python Standard Library offline verifier + AlphaGenome Atlas AVI JSON parser.
Enforces strict molecular & functional framing (research use only; no clinical diagnosis).
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
import re
import sys
from typing import Any, Dict, List, Sequence, Tuple

RESEARCH_BOUNDARY = (
    "Research use only; AlphaGenome AVI reflects predicted molecular/functional impact "
    "and must not be used for clinical diagnosis or medical advice."
)

# GRCh38 CFTR locus bounds (ENSG00000001626, chr7:117,465,784-117,715,971)
CFTR_CHR = "chr7"
CFTR_START_GRCH38 = 117_465_784
CFTR_END_GRCH38 = 117_715_971

_VARIANT_RE = re.compile(r"^(chr[0-9XYM]+):(\d+):([ACGT]+)>([ACGT]+)$")
_DIAGNOSTIC_OVERCLAIM_RE = re.compile(
    r"\b(diagnos(e|is)|prescribe|treat(ment)?\s+plan|clinically\s+actionable\s+cure)\b", re.I
)


@dataclass(frozen=True)
class CftrRegulatoryVariant:
    label: str
    hgvs_c: str
    grch38_variant: str
    expected_modality: str  # "splicing" | "coding_ptc"
    aso_or_readthrough_context: str


CFTR_SPLICE_AND_PTC_PANEL: Tuple[CftrRegulatoryVariant, ...] = (
    CftrRegulatoryVariant(
        "3849+10kb C>T",
        "c.3718-2477C>T",
        "chr7:117642553:C>T",
        "splicing",
        "Cryptic donor activation in intron 22; SPL84 ASO target (NCT06429176)",
    ),
    CftrRegulatoryVariant(
        "2789+5G>A",
        "c.2657+5G>A",
        "chr7:117611653:G>A",
        "splicing",
        "Exon 16 splice donor weakening",
    ),
    CftrRegulatoryVariant(
        "3272-26A>G",
        "c.3140-26A>G",
        "chr7:117614674:A>G",
        "splicing",
        "Cryptic splice acceptor activation upstream of exon 20",
    ),
    CftrRegulatoryVariant(
        "c.1585-1G>A",
        "c.1585-1G>A",
        "chr7:117587747:G>A",
        "splicing",
        "Canonical splice acceptor loss at exon 12",
    ),
    CftrRegulatoryVariant(
        "G542X",
        "c.1624G>T",
        "chr7:117587806:G>T",
        "coding_ptc",
        "Class I premature termination codon in NBD1 (ELX-02 / sup-tRNA target)",
    ),
)


def validate_cftr_grch38_variant(variant_str: str) -> Tuple[str, int, str, str]:
    """Validate chr:pos:ref>alt syntax and verify coordinate lies inside GRCh38 CFTR locus."""
    m = _VARIANT_RE.match(variant_str.strip())
    if not m:
        raise ValueError(f"Invalid chr:pos:ref>alt format: {variant_str!r}")
    chrom, pos_s, ref, alt = m.group(1), m.group(2), m.group(3), m.group(4)
    pos = int(pos_s)
    if ref == alt:
        raise ValueError(f"Reference and alternate alleles are identical: {variant_str!r}")
    if chrom != CFTR_CHR or not (CFTR_START_GRCH38 <= pos <= CFTR_END_GRCH38):
        raise ValueError(f"Coordinate {chrom}:{pos} lies outside GRCh38 CFTR locus")
    return chrom, pos, ref, alt


def parse_alphagenome_avi_json(payload: Dict[str, Any]) -> Dict[str, Any]:
    """Parse an AlphaGenome Atlas AVI query JSON record and enforce molecular-only framing."""
    variant = str(payload.get("variant") or "")
    chrom, pos, _, _ = validate_cftr_grch38_variant(variant)
    summary = str(payload.get("molecular_summary") or "")
    if _DIAGNOSTIC_OVERCLAIM_RE.search(summary):
        raise ValueError(f"Clinical diagnostic overclaim detected in AVI summary: {summary!r}")
    avi_phred = float(payload.get("avi_phred", 0.0))
    if avi_phred < 0.0:
        raise ValueError(f"avi_phred must be non-negative, got {avi_phred}")
    return {
        "variant": variant,
        "chrom": chrom,
        "pos": pos,
        "avi_phred": round(avi_phred, 2),
        "top_modality": str(payload.get("top_modality") or "not_stated"),
        "molecular_summary": summary or "Predicted molecular/functional impact only.",
    }


def evaluate_panel() -> Dict[str, Any]:
    """Validate all 5 panel variants against GRCh38 CFTR coordinates."""
    validated: List[Dict[str, Any]] = []
    for item in CFTR_SPLICE_AND_PTC_PANEL:
        chrom, pos, ref, alt = validate_cftr_grch38_variant(item.grch38_variant)
        validated.append(
            {
                "label": item.label,
                "hgvs_c": item.hgvs_c,
                "chrom": chrom,
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "expected_modality": item.expected_modality,
            }
        )
    return {"panel_size": len(validated), "variants": validated}


def evaluate_falsification_controls() -> Dict[str, Any]:
    """Run deterministic negative controls for coordinate bounds, syntax, and clinical framing."""
    controls: List[Dict[str, Any]] = []

    try:
        validate_cftr_grch38_variant("7-117642553-C-T")
        controls.append({"id": "malformed_variant_rejected", "falsified": False})
    except ValueError:
        controls.append({"id": "malformed_variant_rejected", "falsified": True})

    try:
        validate_cftr_grch38_variant("chr1:100000:A>G")
        controls.append({"id": "non_cftr_locus_rejected", "falsified": False})
    except ValueError:
        controls.append({"id": "non_cftr_locus_rejected", "falsified": True})

    try:
        parse_alphagenome_avi_json(
            {
                "variant": "chr7:117642553:C>T",
                "avi_phred": 24.5,
                "top_modality": "splicing",
                "molecular_summary": "Diagnose CF and prescribe treatment plan.",
            }
        )
        controls.append({"id": "clinical_overclaim_rejected", "falsified": False})
    except ValueError:
        controls.append({"id": "clinical_overclaim_rejected", "falsified": True})

    all_passed = all(c["falsified"] for c in controls)
    return {
        "status": "PASS" if all_passed else "FAIL",
        "controls_passed": sum(1 for c in controls if c["falsified"]),
        "total_controls": len(controls),
        "controls": controls,
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="AlphaGenome CFTR Splicing & Regulatory Variant Probe")
    parser.add_argument("--falsify", action="store_true", help="Run negative-control falsification suite")
    parser.add_argument("--json", action="store_true", help="Emit JSON output")
    args = parser.parse_args(list(argv) if argv is not None else None)

    panel = evaluate_panel()
    falsify = evaluate_falsification_controls() if args.falsify else {"status": "SKIPPED", "controls_passed": 0, "total_controls": 0}
    ok = panel["panel_size"] == 5 and (not args.falsify or falsify["status"] == "PASS")
    payload = {"status": "PASS" if ok else "FAIL", "boundary": RESEARCH_BOUNDARY, "panel": panel, "falsification": falsify}
    if args.json:
        print(json.dumps(payload, indent=2))
    else:
        print(
            f"AlphaGenome CFTR Splice Probe: {payload['status']} "
            f"(panel={panel['panel_size']} GRCh38 CFTR variants, "
            f"falsify={falsify['controls_passed']}/{falsify['total_controls']})"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
