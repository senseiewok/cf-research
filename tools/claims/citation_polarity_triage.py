#!/usr/bin/env python3
"""Section-Weighted Citation Polarity & Retraction Surveillance Model for CF Discovery Literature (scite-api).

Evaluates section-weighted citation tallies (/tallies/cited-by-sections, /tallies/aggregate,
/papers/sources/{target_doi}) across cited CF discovery DOIs:
1. Section-Weighted Dispute Index: Distinguishes empirical/methodological contrasting citations
   occurring in Results or Methods sections (EMPIRICAL_OR_METHODOLOGICAL_ATTENTION) from
   rhetorical framing in Introduction or Discussion (DISCURSIVE_DIFFERENCE).
2. Retracted Citing-Source Contamination Filter: Identifies whether supporting or contrasting
   tallies for a target CF study include retracted citing publications
   (retracted_citing_sources_counted) and computes an adjusted net-support ratio.
3. Zero-Quota Batch Ledger Audit: Verifies that every cited DOI in landscape/ and proposals/
   carries explicit [observed: Scite ML citation classification, accessed YYYY-MM-DD] provenance
   when polarity counts are stated (AH-8).

Primary Literature Grounding:
- Nicholson JM et al. (2021) Quantitative Science Studies, DOI: 10.1162/qss_a_00146
- Middleton PG et al. (2019) N Engl J Med, DOI: 10.1056/NEJMoa1908639
- Keating D et al. (2023) N Engl J Med, DOI: 10.1056/NEJMoa2209892

Clinical Registry Grounding:
- NCT04452240 (Phase 3 ETI long-term safety and efficacy)
- NCT05033080 (Phase 3 Vanzacaftor/Tezacaftor/Deutivacaftor vs ETI)
"""
from __future__ import annotations

import re
from dataclasses import asdict, dataclass
from typing import Any, Dict, Iterable, List


EMPIRICAL_SECTIONS = frozenset({"results", "methods", "materials and methods"})
DISCURSIVE_SECTIONS = frozenset({"introduction", "discussion", "other"})

SCITE_UNLABELLED_POLARITY_RE = re.compile(
    r"\bScite\b[^\n]*?\b(?:supporting|contrasting)\s+(?:citations?|statements?|studies|papers)\b",
    re.IGNORECASE,
)
SCITE_PROVENANCE_LABEL_RE = re.compile(
    r"\[observed:\s*Scite\s+ML\s+citation\s+classification|Scite\s+ML\s+classifier|machine-learning\s+citation",
    re.IGNORECASE,
)


@dataclass(frozen=True)
class SectionCitationTally:
    """Section-level citation polarity counts for a target DOI."""
    section: str
    supporting: int = 0
    contrasting: int = 0
    mentioning: int = 0


@dataclass(frozen=True)
class CitingSourceRecord:
    """Metadata for an individual citing publication referencing a target CF study."""
    citing_doi: str
    polarity: str  # "supporting", "contrasting", or "mentioning"
    section: str = "other"
    is_retracted: bool = False
    editorial_notices: int = 0


@dataclass(frozen=True)
class CitationPolarityEvaluation:
    """Quantitative citation-polarity and methodological-dispute triage outcome."""
    target_doi: str
    total_supporting: int
    total_contrasting: int
    total_mentioning: int
    empirical_contrasting: int
    discursive_contrasting: int
    section_weighted_dispute_index: float
    dispute_classification: str
    retracted_citing_sources_counted: int
    adjusted_supporting: int
    adjusted_contrasting: int
    adjusted_net_support_ratio: float
    provenance_label: str

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


def compute_section_weighted_dispute_index(
    section_tallies: Iterable[SectionCitationTally],
    empirical_weight: float = 3.0,
    discursive_weight: float = 1.0,
) -> Dict[str, Any]:
    """Compute section-weighted methodological dispute index from section tallies.

    Contrasting citations in Results or Methods sections carry higher empirical weight than
    contrasting citations in Introduction or Discussion sections.
    """
    if empirical_weight <= 0 or discursive_weight <= 0:
        raise ValueError("Section weights must be positive")

    total_sup = 0
    total_con = 0
    total_men = 0
    emp_con = 0
    disc_con = 0

    for item in section_tallies:
        if item.supporting < 0 or item.contrasting < 0 or item.mentioning < 0:
            raise ValueError("Citation tallies cannot be negative")
        sec_norm = item.section.strip().lower()
        total_sup += item.supporting
        total_con += item.contrasting
        total_men += item.mentioning
        if sec_norm in EMPIRICAL_SECTIONS:
            emp_con += item.contrasting
        else:
            disc_con += item.contrasting

    weighted_con = empirical_weight * emp_con + discursive_weight * disc_con
    denom = weighted_con + float(total_sup) + 1.0
    dispute_index = round(weighted_con / denom, 4)

    if emp_con >= 2 or (emp_con >= 1 and dispute_index >= 0.20):
        classification = "EMPIRICAL_OR_METHODOLOGICAL_ATTENTION"
    elif total_con > 0:
        classification = "DISCURSIVE_DIFFERENCE"
    else:
        classification = "CONSENSUS_OR_UNCONTESTED"

    return {
        "total_supporting": total_sup,
        "total_contrasting": total_con,
        "total_mentioning": total_men,
        "empirical_contrasting": emp_con,
        "discursive_contrasting": disc_con,
        "section_weighted_dispute_index": dispute_index,
        "dispute_classification": classification,
    }


def filter_retracted_citing_sources(
    total_supporting: int,
    total_contrasting: int,
    citing_sources: Iterable[CitingSourceRecord] = (),
) -> Dict[str, Any]:
    """Filter retracted citing publications from supporting/contrasting tallies."""
    if total_supporting < 0 or total_contrasting < 0:
        raise ValueError("Tallies cannot be negative")

    retracted_sup = 0
    retracted_con = 0
    retracted_total = 0

    for src in citing_sources:
        if src.is_retracted or src.editorial_notices > 0:
            retracted_total += 1
            pol = src.polarity.strip().lower()
            if pol == "supporting":
                retracted_sup += 1
            elif pol == "contrasting":
                retracted_con += 1

    adj_sup = max(0, total_supporting - retracted_sup)
    adj_con = max(0, total_contrasting - retracted_con)
    polarized_total = adj_sup + adj_con
    net_support_ratio = round((adj_sup - adj_con) / polarized_total, 4) if polarized_total > 0 else 0.0

    return {
        "retracted_citing_sources_counted": retracted_total,
        "retracted_supporting_removed": retracted_sup,
        "retracted_contrasting_removed": retracted_con,
        "adjusted_supporting": adj_sup,
        "adjusted_contrasting": adj_con,
        "adjusted_net_support_ratio": net_support_ratio,
    }


def evaluate_citation_polarity_triage(
    target_doi: str,
    section_tallies: Iterable[SectionCitationTally],
    citing_sources: Iterable[CitingSourceRecord] = (),
    accessed_date: str = "2026-10-24",
) -> CitationPolarityEvaluation:
    """Evaluate section-weighted dispute classification and retracted citing-source contamination."""
    sec_res = compute_section_weighted_dispute_index(section_tallies)
    ret_res = filter_retracted_citing_sources(
        sec_res["total_supporting"],
        sec_res["total_contrasting"],
        citing_sources,
    )
    return CitationPolarityEvaluation(
        target_doi=target_doi,
        total_supporting=sec_res["total_supporting"],
        total_contrasting=sec_res["total_contrasting"],
        total_mentioning=sec_res["total_mentioning"],
        empirical_contrasting=sec_res["empirical_contrasting"],
        discursive_contrasting=sec_res["discursive_contrasting"],
        section_weighted_dispute_index=sec_res["section_weighted_dispute_index"],
        dispute_classification=sec_res["dispute_classification"],
        retracted_citing_sources_counted=ret_res["retracted_citing_sources_counted"],
        adjusted_supporting=ret_res["adjusted_supporting"],
        adjusted_contrasting=ret_res["adjusted_contrasting"],
        adjusted_net_support_ratio=ret_res["adjusted_net_support_ratio"],
        provenance_label=f"[observed: Scite ML citation classification, accessed {accessed_date}]",
    )


def audit_scite_polarity_provenance_text(text: str) -> List[str]:
    """Verify that any statement citing Scite supporting/contrasting counts carries AH-8 ML provenance."""
    violations: List[str] = []
    lines = text.splitlines()
    for idx, line in enumerate(lines):
        window = " ".join(lines[max(0, idx - 1) : min(len(lines), idx + 2)])
        if SCITE_UNLABELLED_POLARITY_RE.search(line) and not SCITE_PROVENANCE_LABEL_RE.search(window):
            violations.append(
                f"Line {idx + 1}: Scite supporting/contrasting polarity cited without "
                "[observed: Scite ML citation classification] provenance label"
            )
    return violations
