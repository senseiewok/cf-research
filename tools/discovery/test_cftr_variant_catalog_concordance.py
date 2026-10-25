"""Unit tests for CFTR Missense & Splice Variant Cross-Catalog Concordance Model (T-0420)."""
from __future__ import annotations

import pytest
from cftr_variant_catalog_concordance import (
    CANONICAL_ENSEMBL_TRANSCRIPT,
    CANONICAL_REFSEQ_TRANSCRIPT,
    CANONICAL_UNIPROT_ACCESSION,
    PRIORITY_CFTR_VARIANT_CATALOG,
    evaluate_variant_catalog_concordance,
    map_codon_to_cftr_domain,
    summarize_catalog_concordance_panel,
    verify_spdi_hgvs_coordinate_shift,
)


def test_spdi_hgvs_coordinate_shift_f508del_and_snvs():
    """Verify 0-based SPDI vs 1-based GRCh38 coordinate normalization."""
    f508 = verify_spdi_hgvs_coordinate_shift("NC_000007.14:117559592:CTT:", 117559593)
    assert f508["is_concordant"] is True
    assert f508["coordinate_shift"] == 1
    assert f508["deleted_length"] == 3
    assert f508["inserted_length"] == 0

    g551d = verify_spdi_hgvs_coordinate_shift("NC_000007.14:117587805:G:A", 117587806)
    assert g551d["is_concordant"] is True
    assert g551d["deleted_length"] == 1
    assert g551d["inserted_length"] == 1

    with pytest.raises(ValueError):
        verify_spdi_hgvs_coordinate_shift("invalid_spdi", 100)


def test_map_codon_to_cftr_domain_boundaries():
    """Verify UniProtKB P13569 and PDB 6MSM structural domain mapping across CFTR."""
    assert map_codon_to_cftr_domain(117) == "TMD1"
    assert map_codon_to_cftr_domain(508) == "NBD1"
    assert map_codon_to_cftr_domain(542) == "NBD1"
    assert map_codon_to_cftr_domain(551) == "NBD1"
    assert map_codon_to_cftr_domain(750) == "R"
    assert map_codon_to_cftr_domain(1000) == "TMD2"
    assert map_codon_to_cftr_domain(1282) == "NBD2"
    assert map_codon_to_cftr_domain(1303) == "NBD2"
    assert map_codon_to_cftr_domain(None) == "INTRON_22_CRYPTIC_SPLICE"

    with pytest.raises(ValueError):
        map_codon_to_cftr_domain(0)
    with pytest.raises(ValueError):
        map_codon_to_cftr_domain(1481)


def test_evaluate_variant_catalog_concordance_theratypes():
    """Verify modulator eligibility vs refractory nucleic acid pathways across priority alleles."""
    f508 = evaluate_variant_catalog_concordance("F508del")
    assert f508["concordance_verdict"] == "CONCORDANT"
    assert f508["structural_domain"] == "NBD1"
    assert f508["fda_eti_responsive"] is True
    assert f508["fda_ivacaftor_monotherapy_responsive"] is False
    assert f508["regulatory_pathway"] == "TRIPLE_MODULATOR_ELIGIBLE"
    assert f508["transcript_refseq"] == CANONICAL_REFSEQ_TRANSCRIPT
    assert f508["transcript_ensembl"] == CANONICAL_ENSEMBL_TRANSCRIPT
    assert f508["uniprot_accession"] == CANONICAL_UNIPROT_ACCESSION

    g551d = evaluate_variant_catalog_concordance("G551D")
    assert g551d["regulatory_pathway"] == "MODULATOR_MONOTHERAPY_AND_TRIPLE_ELIGIBLE"

    w1282x = evaluate_variant_catalog_concordance("W1282X")
    assert w1282x["fda_eti_responsive"] is False
    assert w1282x["requires_genetic_replacement_or_aso"] is True
    assert w1282x["regulatory_pathway"] == "MODULATOR_REFRACTORY_NUCLEIC_ACID_PATHWAY"
    assert "NCT05248230" in w1282x["active_clinical_trials"]
    assert "NCT05668741" in w1282x["active_clinical_trials"]

    with pytest.raises(ValueError):
        evaluate_variant_catalog_concordance("UNKNOWN_ALLELE")


def test_summarize_catalog_concordance_panel():
    """Verify full 7-variant panel concordance summary."""
    summary = summarize_catalog_concordance_panel()
    assert summary["total_variants_evaluated"] == len(PRIORITY_CFTR_VARIANT_CATALOG) == 7
    assert summary["all_concordant"] is True
    assert summary["concordant_variants"] == 7
    assert summary["modulator_eligible_variants"] == 4
    assert summary["modulator_refractory_variants"] == 3
