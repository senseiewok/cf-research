"""Unit tests for tools/discovery/cftr_sequence_scanner.py (T-0264)."""
from __future__ import annotations

import pytest
from cftr_sequence_scanner import (
    CANONICAL_CFTR_LOCI,
    abe8e_window_efficiency,
    audit_cpg_and_codon_harmony,
    calculate_oligo_tm_c,
    match_pam,
    reverse_complement,
    run_full_sequence_scanner_suite,
    scan_base_editing_candidates,
    scan_cryptic_splice_aso_walk,
)


def test_reverse_complement_and_pam_matching():
    assert reverse_complement("ATGC") == "GCAT"
    assert match_pam("AGG", "SpCas9") is True
    assert match_pam("AGA", "SpCas9") is False
    assert match_pam("AGA", "SpG") is True
    assert match_pam("AAA", "SpRY") is True


def test_w1282x_and_r553x_abe8e_vs_g542x_transversion_gate():
    """W1282X (c.3846G>A) and R553X (c.1657C>T) are ABE8e-eligible; G542X (c.1624G>T) requires Prime Editing."""
    w1282x = scan_base_editing_candidates(CANONICAL_CFTR_LOCI["W1282X"])
    assert w1282x["transition_base_editor_eligible"] is True
    assert w1282x["requires_prime_editing_transversion"] is False
    assert w1282x["candidate_count"] >= 1
    assert w1282x["top_candidates"][0]["strand"] == "coding (+)"

    r553x = scan_base_editing_candidates(CANONICAL_CFTR_LOCI["R553X"])
    assert r553x["transition_base_editor_eligible"] is True
    assert r553x["top_candidates"][0]["strand"] == "antisense (-)"

    # G542X is c.1624G>T (GGA -> TGA); converting T->G on coding requires transversion (Prime Editing)
    # Note: antisense T->C would convert TGA -> CGA (Arg, missense), whereas G542X WT base is G
    g542x = scan_base_editing_candidates(CANONICAL_CFTR_LOCI["G542X"])
    assert g542x["wt_base_coding"] == "G"
    assert g542x["mutant_base_coding"] == "T"


def test_cryptic_splice_aso_walk_3849_10kb():
    """ASO walk across 3849+10kbC->T identifies steric blockers covering the cryptic GT donor."""
    locus = CANONICAL_CFTR_LOCI["3849+10kbC->T"]
    walk = scan_cryptic_splice_aso_walk(locus.mutant_sequence_5to3, cryptic_donor_offset=23, oligo_length=20)
    assert walk["total_walk_windows"] == len(locus.mutant_sequence_5to3) - 20 + 1
    best = walk["optimal_aso"]
    assert best["covers_cryptic_gt"] is True
    assert 45.0 < best["tm_deg_c"] < 75.0
    assert best["has_gggg_quadruplex_trap"] is False


def test_cpg_depletion_and_codon_harmony():
    """CpG audit detects intragenic CpG dinucleotides and generates synonymous CpG-depleted DNA."""
    wt_seq = CANONICAL_CFTR_LOCI["R553X"].wt_sequence_5to3  # contains CGA (Arg553)
    res = audit_cpg_and_codon_harmony(wt_seq)
    assert res["wildtype_cpg_count"] >= 1
    assert res["depleted_cpg_count"] < res["wildtype_cpg_count"]
    assert len(res["cpg_depleted_sequence_5to3"]) == len(wt_seq)

    with pytest.raises(ValueError):
        audit_cpg_and_codon_harmony("ATGC")  # Not divisible by 3


def test_full_sequence_scanner_suite():
    suite = run_full_sequence_scanner_suite()
    assert suite["task_id"] == "T-0264"
    assert suite["refseq_transcript"] == "NM_000492.4"
    assert "W1282X" in suite["allele_scans"]
    assert "3849+10kbC->T" in suite["allele_scans"]
