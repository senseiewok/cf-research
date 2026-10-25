"""Unit tests for CFTR Missense & Splice Variant Cross-Catalog Concordance Model (T-0420)."""
from __future__ import annotations

import sys
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent.parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(REPO_ROOT))

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


class TestCftrVariantCatalogConcordance(unittest.TestCase):
    def test_spdi_hgvs_coordinate_shift_f508del_and_snvs(self):
        """Verify 0-based SPDI vs 1-based GRCh38 coordinate normalization."""
        f508 = verify_spdi_hgvs_coordinate_shift("NC_000007.14:117559592:CTT:", 117559593)
        self.assertTrue(f508["is_concordant"])
        self.assertEqual(f508["coordinate_shift"], 1)
        self.assertEqual(f508["deleted_length"], 3)
        self.assertEqual(f508["inserted_length"], 0)

        g551d = verify_spdi_hgvs_coordinate_shift("NC_000007.14:117587805:G:A", 117587806)
        self.assertTrue(g551d["is_concordant"])
        self.assertEqual(g551d["deleted_length"], 1)
        self.assertEqual(g551d["inserted_length"], 1)

        with self.assertRaises(ValueError):
            verify_spdi_hgvs_coordinate_shift("invalid_spdi", 100)

    def test_map_codon_to_cftr_domain_boundaries(self):
        """Verify UniProtKB P13569 and PDB 6MSM structural domain mapping across CFTR."""
        self.assertEqual(map_codon_to_cftr_domain(117), "TMD1")
        self.assertEqual(map_codon_to_cftr_domain(508), "NBD1")
        self.assertEqual(map_codon_to_cftr_domain(542), "NBD1")
        self.assertEqual(map_codon_to_cftr_domain(551), "NBD1")
        self.assertEqual(map_codon_to_cftr_domain(750), "R")
        self.assertEqual(map_codon_to_cftr_domain(1000), "TMD2")
        self.assertEqual(map_codon_to_cftr_domain(1282), "NBD2")
        self.assertEqual(map_codon_to_cftr_domain(1303), "NBD2")
        self.assertEqual(map_codon_to_cftr_domain(None), "INTRON_22_CRYPTIC_SPLICE")

        with self.assertRaises(ValueError):
            map_codon_to_cftr_domain(0)
        with self.assertRaises(ValueError):
            map_codon_to_cftr_domain(1481)

    def test_evaluate_variant_catalog_concordance_theratypes(self):
        """Verify modulator eligibility vs refractory nucleic acid pathways across priority alleles."""
        f508 = evaluate_variant_catalog_concordance("F508del")
        self.assertEqual(f508["concordance_verdict"], "CONCORDANT")
        self.assertEqual(f508["structural_domain"], "NBD1")
        self.assertTrue(f508["fda_eti_responsive"])
        self.assertFalse(f508["fda_ivacaftor_monotherapy_responsive"])
        self.assertEqual(f508["regulatory_pathway"], "TRIPLE_MODULATOR_ELIGIBLE")
        self.assertEqual(f508["transcript_refseq"], CANONICAL_REFSEQ_TRANSCRIPT)
        self.assertEqual(f508["transcript_ensembl"], CANONICAL_ENSEMBL_TRANSCRIPT)
        self.assertEqual(f508["uniprot_accession"], CANONICAL_UNIPROT_ACCESSION)

        g551d = evaluate_variant_catalog_concordance("G551D")
        self.assertEqual(g551d["regulatory_pathway"], "MODULATOR_MONOTHERAPY_AND_TRIPLE_ELIGIBLE")

        w1282x = evaluate_variant_catalog_concordance("W1282X")
        self.assertFalse(w1282x["fda_eti_responsive"])
        self.assertTrue(w1282x["requires_genetic_replacement_or_aso"])
        self.assertEqual(w1282x["regulatory_pathway"], "MODULATOR_REFRACTORY_NUCLEIC_ACID_PATHWAY")
        self.assertIn("NCT05248230", w1282x["active_clinical_trials"])
        self.assertIn("NCT05668741", w1282x["active_clinical_trials"])

        with self.assertRaises(ValueError):
            evaluate_variant_catalog_concordance("UNKNOWN_ALLELE")

    def test_summarize_catalog_concordance_panel(self):
        """Verify full 7-variant panel concordance summary."""
        summary = summarize_catalog_concordance_panel()
        self.assertEqual(summary["total_variants_evaluated"], len(PRIORITY_CFTR_VARIANT_CATALOG))
        self.assertEqual(summary["total_variants_evaluated"], 7)
        self.assertTrue(summary["all_concordant"])
        self.assertEqual(summary["concordant_variants"], 7)
        self.assertEqual(summary["modulator_eligible_variants"], 4)
        self.assertEqual(summary["modulator_refractory_variants"], 3)


if __name__ == "__main__":
    unittest.main()
