#!/usr/bin/env python3
"""Unit tests for Rare CFTR Frameshift, Large Deletion & Splicing Theratyping Matrix (T-0201)."""
from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(REPO_ROOT))

from rare_allele_matrix import (
    DEFAULT_SNAPSHOT_JSON,
    MODALITIES,
    NLM_ATTRIBUTION,
    RARE_ALLELE_CATALOG,
    RESEARCH_BOUNDARY,
    build_theratyping_matrix,
    evaluate_modality_feasibility,
    load_snapshot_nct_ids,
    self_test,
)


class TestRareAlleleMatrix(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.snapshot_ncts = load_snapshot_nct_ids(DEFAULT_SNAPSHOT_JSON)

    def test_catalog_completeness(self):
        """Verify catalog contains all 12 key rare Class I/V/VII alleles."""
        expected_alleles = {
            "G542X",
            "W1282X",
            "R553X",
            "3905insT",
            "2184insA",
            "CFTRdele2,3",
            "1717-1G->A",
            "621+1G->T",
            "3849+10kb C->T",
            "2789+5G->A",
            "1811+1.6kbA->G",
            "N1303K",
        }
        self.assertTrue(expected_alleles.issubset(set(RARE_ALLELE_CATALOG.keys())))
        for name in expected_alleles:
            prof = RARE_ALLELE_CATALOG[name]
            self.assertTrue(prof.cdna_name)
            self.assertTrue(prof.protein_name)
            self.assertTrue(prof.variant_class)
            self.assertTrue(prof.mutation_type)

    def test_frameshift_feasibility_invariants(self):
        """Verify frameshifts (3905insT, 2184insA) score 0 on modulators/readthrough/ASO and high on mRNA/gene addition."""
        for fs in ("3905insT", "2184insA"):
            prof = RARE_ALLELE_CATALOG[fs]
            self.assertEqual(prof.mutation_type, "frameshift")

            res_mod = evaluate_modality_feasibility(prof, "CFTR_MODULATOR", self.snapshot_ncts)
            self.assertEqual(res_mod.feasibility_score, 0.0)
            self.assertEqual(res_mod.feasibility_verdict, "Incompatible")

            res_rt = evaluate_modality_feasibility(prof, "READTHROUGH_NMD", self.snapshot_ncts)
            self.assertEqual(res_rt.feasibility_score, 0.0)
            self.assertEqual(res_rt.feasibility_verdict, "Incompatible")

            res_aso = evaluate_modality_feasibility(prof, "SPLICE_ASO", self.snapshot_ncts)
            self.assertEqual(res_aso.feasibility_score, 0.0)
            self.assertEqual(res_aso.feasibility_verdict, "Incompatible")

            res_mrna = evaluate_modality_feasibility(prof, "MRNA_LNP", self.snapshot_ncts)
            self.assertGreaterEqual(res_mrna.feasibility_score, 0.85)

            res_gene = evaluate_modality_feasibility(prof, "VIRAL_GENE", self.snapshot_ncts)
            self.assertGreaterEqual(res_gene.feasibility_score, 0.85)

    def test_large_deletion_cftrdele23(self):
        """Verify CFTRdele2,3 scores 0 on modulators, readthrough, ASO, and low on point editing."""
        prof = RARE_ALLELE_CATALOG["CFTRdele2,3"]
        self.assertEqual(prof.variant_class, "Class VII")
        self.assertEqual(prof.mutation_type, "large_deletion")

        res_mod = evaluate_modality_feasibility(prof, "CFTR_MODULATOR", self.snapshot_ncts)
        self.assertEqual(res_mod.feasibility_score, 0.0)

        res_rt = evaluate_modality_feasibility(prof, "READTHROUGH_NMD", self.snapshot_ncts)
        self.assertEqual(res_rt.feasibility_score, 0.0)

        res_aso = evaluate_modality_feasibility(prof, "SPLICE_ASO", self.snapshot_ncts)
        self.assertEqual(res_aso.feasibility_score, 0.0)

        res_edit = evaluate_modality_feasibility(prof, "GENE_EDITING", self.snapshot_ncts)
        self.assertLessEqual(res_edit.feasibility_score, 0.20)

        res_mrna = evaluate_modality_feasibility(prof, "MRNA_LNP", self.snapshot_ncts)
        self.assertGreaterEqual(res_mrna.feasibility_score, 0.85)

        res_gene = evaluate_modality_feasibility(prof, "VIRAL_GENE", self.snapshot_ncts)
        self.assertGreaterEqual(res_gene.feasibility_score, 0.85)

    def test_cryptic_splice_3849_10kb(self):
        """Verify 3849+10kb C->T scores highest on splice-switching ASO and links to NCT06429176."""
        prof = RARE_ALLELE_CATALOG["3849+10kb C->T"]
        self.assertEqual(prof.variant_class, "Class V")
        self.assertTrue(prof.protein_substrate_present)

        res_aso = evaluate_modality_feasibility(prof, "SPLICE_ASO", self.snapshot_ncts)
        self.assertGreaterEqual(res_aso.feasibility_score, 0.90)
        self.assertIn("NCT06429176", res_aso.relevant_nct_ids)
        self.assertEqual(res_aso.feasibility_verdict, "Clinically Validated")

        res_mod = evaluate_modality_feasibility(prof, "CFTR_MODULATOR", self.snapshot_ncts)
        self.assertGreaterEqual(res_mod.feasibility_score, 0.65)

    def test_canonical_splice_mutations(self):
        """Verify canonical splice variants (1717-1G->A, 621+1G->T) score 0 on splice ASO."""
        for cs in ("1717-1G->A", "621+1G->T"):
            prof = RARE_ALLELE_CATALOG[cs]
            self.assertEqual(prof.mutation_type, "canonical_splice")

            res_aso = evaluate_modality_feasibility(prof, "SPLICE_ASO", self.snapshot_ncts)
            self.assertEqual(res_aso.feasibility_score, 0.0)
            self.assertEqual(res_aso.feasibility_verdict, "Incompatible")

            res_mod = evaluate_modality_feasibility(prof, "CFTR_MODULATOR", self.snapshot_ncts)
            self.assertEqual(res_mod.feasibility_score, 0.0)

    def test_nonsense_mutations_readthrough(self):
        """Verify nonsense alleles (G542X, W1282X, R553X) score high on readthrough and mRNA."""
        for ns in ("G542X", "W1282X", "R553X"):
            prof = RARE_ALLELE_CATALOG[ns]
            self.assertEqual(prof.mutation_type, "nonsense")

            res_rt = evaluate_modality_feasibility(prof, "READTHROUGH_NMD", self.snapshot_ncts)
            self.assertGreaterEqual(res_rt.feasibility_score, 0.70)

            res_mrna = evaluate_modality_feasibility(prof, "MRNA_LNP", self.snapshot_ncts)
            self.assertGreaterEqual(res_mrna.feasibility_score, 0.85)

    def test_missense_misfolding_n1303k(self):
        """Verify N1303K scores moderately on modulators and 0 on readthrough."""
        prof = RARE_ALLELE_CATALOG["N1303K"]
        self.assertEqual(prof.variant_class, "Class II")

        res_mod = evaluate_modality_feasibility(prof, "CFTR_MODULATOR", self.snapshot_ncts)
        self.assertGreaterEqual(res_mod.feasibility_score, 0.60)

        res_rt = evaluate_modality_feasibility(prof, "READTHROUGH_NMD", self.snapshot_ncts)
        self.assertEqual(res_rt.feasibility_score, 0.0)

        res_aso = evaluate_modality_feasibility(prof, "SPLICE_ASO", self.snapshot_ncts)
        self.assertEqual(res_aso.feasibility_score, 0.0)

    def test_nct_provenance_guard(self):
        """Rule AH-1: Ensure 100% of cited NCT IDs exist in snapshot.json."""
        if not self.snapshot_ncts:
            self.skipTest("snapshot.json not found")
        for prof in RARE_ALLELE_CATALOG.values():
            for mod in MODALITIES:
                res = evaluate_modality_feasibility(prof, mod, self.snapshot_ncts)
                for nct in res.relevant_nct_ids:
                    self.assertIn(
                        nct,
                        self.snapshot_ncts,
                        f"NCT {nct} cited for {prof.allele_name} not in snapshot.json",
                    )

    def test_registry_non_pooling_imperative(self):
        """Verify registry notes keep national/source counts unpooled."""
        for prof in RARE_ALLELE_CATALOG.values():
            ecfspr_note = prof.registry_notes.get("ecfspr_2024", "")
            self.assertNotIn("global pooled percentage", ecfspr_note.lower())
            self.assertTrue(len(prof.registry_notes) >= 2)

    def test_nlm_attribution_and_research_boundary(self):
        """Verify NLM attribution and research boundary disclaimers."""
        self.assertIn("National Library of Medicine", NLM_ATTRIBUTION)
        self.assertIn("Research, not medical advice", RESEARCH_BOUNDARY)

    def test_no_clinical_overclaim_phrases(self):
        """Rule AH-4: Verify no prohibited clinical overclaim phrases in descriptions."""
        prohibited = (
            "proven to cure",
            "guarantees clinical benefit",
            "recommended treatment regimen",
            "safe and effective for patients",
            "cures cystic fibrosis",
        )
        source_code = Path(__file__).resolve().parent.joinpath("rare_allele_matrix.py").read_text(encoding="utf-8")
        for p in prohibited:
            self.assertNotIn(p.lower(), source_code.lower())

    def test_internal_self_test(self):
        """Verify self_test() returns True."""
        self.assertTrue(self_test())

    def test_cli_execution_and_json(self):
        """Verify CLI execution with --self-test and --json."""
        cmd_self = [sys.executable, "-m", "tools.discovery.rare_allele_matrix", "--self-test"]
        proc_self = subprocess.run(cmd_self, capture_output=True, text=True, cwd=str(REPO_ROOT))
        self.assertEqual(proc_self.returncode, 0)
        self.assertIn("PASS", proc_self.stdout)

        cmd_json = [sys.executable, "-m", "tools.discovery.rare_allele_matrix", "--allele", "G542X", "--json"]
        proc_json = subprocess.run(cmd_json, capture_output=True, text=True, cwd=str(REPO_ROOT))
        self.assertEqual(proc_json.returncode, 0)
        parsed = json.loads(proc_json.stdout)
        self.assertIn("matrix", parsed)
        self.assertIn("G542X", parsed["matrix"])


if __name__ == "__main__":
    unittest.main()
