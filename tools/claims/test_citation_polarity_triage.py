"""Unit tests for Section-Weighted Citation Polarity & Retraction Surveillance Model (T-0421)."""
from __future__ import annotations

import pytest
from citation_polarity_triage import (
    CitingSourceRecord,
    SectionCitationTally,
    audit_scite_polarity_provenance_text,
    compute_section_weighted_dispute_index,
    evaluate_citation_polarity_triage,
    filter_retracted_citing_sources,
)


def test_section_weighted_dispute_index_empirical_vs_discursive():
    """Contrasting citations in Results/Methods must trigger EMPIRICAL_OR_METHODOLOGICAL_ATTENTION."""
    empirical_tallies = [
        SectionCitationTally(section="results", supporting=12, contrasting=3, mentioning=45),
        SectionCitationTally(section="discussion", supporting=8, contrasting=1, mentioning=30),
    ]
    res_emp = compute_section_weighted_dispute_index(empirical_tallies)
    assert res_emp["empirical_contrasting"] == 3
    assert res_emp["discursive_contrasting"] == 1
    assert res_emp["dispute_classification"] == "EMPIRICAL_OR_METHODOLOGICAL_ATTENTION"
    assert res_emp["section_weighted_dispute_index"] > 0.25

    discursive_tallies = [
        SectionCitationTally(section="introduction", supporting=25, contrasting=1, mentioning=80),
        SectionCitationTally(section="discussion", supporting=15, contrasting=1, mentioning=40),
        SectionCitationTally(section="results", supporting=10, contrasting=0, mentioning=20),
    ]
    res_disc = compute_section_weighted_dispute_index(discursive_tallies)
    assert res_disc["empirical_contrasting"] == 0
    assert res_disc["discursive_contrasting"] == 2
    assert res_disc["dispute_classification"] == "DISCURSIVE_DIFFERENCE"

    consensus_tallies = [
        SectionCitationTally(section="results", supporting=18, contrasting=0, mentioning=50),
    ]
    res_cons = compute_section_weighted_dispute_index(consensus_tallies)
    assert res_cons["dispute_classification"] == "CONSENSUS_OR_UNCONTESTED"


def test_filter_retracted_citing_sources():
    """Verify retracted citing publications are subtracted from supporting and contrasting tallies."""
    sources = [
        CitingSourceRecord(citing_doi="10.1000/retracted1", polarity="contrasting", section="results", is_retracted=True),
        CitingSourceRecord(citing_doi="10.1000/retracted2", polarity="supporting", section="discussion", editorial_notices=1),
        CitingSourceRecord(citing_doi="10.1000/clean1", polarity="supporting", section="results", is_retracted=False),
    ]
    filtered = filter_retracted_citing_sources(total_supporting=10, total_contrasting=2, citing_sources=sources)
    assert filtered["retracted_citing_sources_counted"] == 2
    assert filtered["adjusted_supporting"] == 9
    assert filtered["adjusted_contrasting"] == 1
    assert filtered["adjusted_net_support_ratio"] == 0.8


def test_evaluate_citation_polarity_triage_end_to_end():
    """Verify end-to-end citation polarity evaluation with AH-8 provenance label."""
    tallies = [
        SectionCitationTally(section="methods", supporting=4, contrasting=2, mentioning=15),
        SectionCitationTally(section="discussion", supporting=16, contrasting=1, mentioning=40),
    ]
    sources = [
        CitingSourceRecord(citing_doi="10.1002/path.1597", polarity="contrasting", section="methods", is_retracted=True),
    ]
    eval_res = evaluate_citation_polarity_triage(
        "10.1056/NEJMoa1908639",
        section_tallies=tallies,
        citing_sources=sources,
        accessed_date="2026-10-24",
    )
    assert eval_res.target_doi == "10.1056/NEJMoa1908639"
    assert eval_res.dispute_classification == "EMPIRICAL_OR_METHODOLOGICAL_ATTENTION"
    assert eval_res.retracted_citing_sources_counted == 1
    assert eval_res.adjusted_supporting == 20
    assert eval_res.adjusted_contrasting == 2
    assert "[observed: Scite ML citation classification, accessed 2026-10-24]" in eval_res.provenance_label


def test_audit_scite_polarity_provenance_text():
    """Verify AH-8 provenance check catches unlabelled Scite polarity claims."""
    bad_text = "Scite reports 42 supporting citations and 2 contrasting statements for this trial."
    assert len(audit_scite_polarity_provenance_text(bad_text)) == 1

    good_text = (
        "Scite reports 42 supporting citations and 2 contrasting statements for this trial "
        "[observed: Scite ML citation classification, accessed 2026-10-24]."
    )
    assert audit_scite_polarity_provenance_text(good_text) == []

    with pytest.raises(ValueError):
        compute_section_weighted_dispute_index([SectionCitationTally(section="results", supporting=-1)])
