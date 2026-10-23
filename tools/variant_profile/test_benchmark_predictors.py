#!/usr/bin/env python3
"""Unit tests for tools/variant_profile/benchmark_predictors.py."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from tools.variant_profile.benchmark_predictors import (
    NOT_STATED,
    build_variant_domain_benchmark,
    classify_molecular_mechanism,
    evaluate_predictor_concordance,
    extract_aa_position,
    main,
    map_residue_to_domain,
    map_residue_to_exon,
)


class TestBenchmarkPredictors(unittest.TestCase):
    def test_domain_and_exon_mapping_boundaries(self) -> None:
        self.assertEqual(map_residue_to_domain(31)["domain_id"], "Lasso_N_terminal")
        self.assertEqual(map_residue_to_domain(117)["domain_id"], "TMD1")
        self.assertEqual(map_residue_to_domain(508)["domain_id"], "NBD1")
        self.assertEqual(map_residue_to_domain(700)["domain_id"], "R_domain")
        self.assertEqual(map_residue_to_domain(1152)["domain_id"], "TMD2")
        self.assertEqual(map_residue_to_domain(1303)["domain_id"], "NBD2_C_terminal")
        self.assertEqual(map_residue_to_domain(None)["domain_id"], NOT_STATED)
        self.assertEqual(map_residue_to_domain(9999)["domain_id"], NOT_STATED)

        self.assertEqual(map_residue_to_exon(31)["exon_label"], "2/27")
        self.assertEqual(map_residue_to_exon(508)["exon_label"], "11/27")
        self.assertEqual(map_residue_to_exon(1282)["exon_label"], "23/27")
        self.assertEqual(map_residue_to_exon(1303)["exon_label"], "24/27")
        self.assertTrue(map_residue_to_exon(55)["near_splice_junction_le_3aa"])
        self.assertEqual(map_residue_to_exon(None)["exon_label"], NOT_STATED)

    def test_extract_aa_position_and_mechanism(self) -> None:
        self.assertEqual(extract_aa_position("R31L"), 31)
        self.assertEqual(extract_aa_position("F508del"), 508)
        self.assertEqual(extract_aa_position("W1282*"), 1282)
        self.assertEqual(extract_aa_position("R104fs"), 104)

        self.assertEqual(
            classify_molecular_mechanism("F508del", "Deletion", "p.Phe508del")["consequence_class"],
            "inframe_deletion",
        )
        self.assertEqual(
            classify_molecular_mechanism("W1282*", "single nucleotide variant", "p.Trp1282Ter")["consequence_class"],
            "nonsense_premature_termination",
        )
        self.assertEqual(
            classify_molecular_mechanism("R117H", "Haplotype", "c.[350G>A;1210-12[5]]")["consequence_class"],
            "complex_haplotype",
        )

    def test_predictor_concordance_detects_r31l_split(self) -> None:
        res = evaluate_predictor_concordance("Uncertain significance", 22.9, "benign", "tolerated")
        self.assertEqual(res["status"], "evaluated_from_fixture")
        self.assertTrue(res["cadd_ge_20_deleterious"])
        self.assertFalse(res["polyphen_damaging"])
        self.assertFalse(res["sift_deleterious"])
        self.assertTrue(res["predictors_disagree_internally"])

    def test_build_benchmark_completeness_and_not_stated_discipline(self) -> None:
        bench = build_variant_domain_benchmark()
        sm = bench["summary_metrics"]
        self.assertEqual(sm["total_clinvar_records_evaluated"], 17)
        self.assertGreaterEqual(sm["domains_represented"], 4)
        self.assertIn("N1303K", bench["multi_nucleotide_degeneracy_audit"])
        self.assertIn("F508del", bench["multi_nucleotide_degeneracy_audit"])

        by_uid = {r["clinvar_variation_id"]["value"]: r for r in bench["variant_records"]}
        r31l = by_uid["54087"]
        self.assertEqual(r31l["in_silico_annotations"]["cadd_phred"]["value"], 22.9)
        self.assertEqual(r31l["structural_domain"]["domain_id"], "Lasso_N_terminal")
        self.assertEqual(r31l["refseq_exon"]["exon_label"], "2/27")

        d1152h = by_uid["35867"]
        self.assertEqual(d1152h["clinvar_classification"]["description"]["value"], "drug response")
        self.assertEqual(d1152h["in_silico_annotations"]["cadd_phred"]["value"], NOT_STATED)
        self.assertEqual(d1152h["refseq_exon"]["exon_label"], "21/27")

    def test_cli_check_and_write(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out_path = Path(tmp) / "domain_bench.json"
            rc = main(["--out", str(out_path)])
            self.assertEqual(rc, 0)
            self.assertTrue(out_path.exists())
            loaded = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertEqual(loaded["summary_metrics"]["total_clinvar_records_evaluated"], 17)


if __name__ == "__main__":
    unittest.main()
