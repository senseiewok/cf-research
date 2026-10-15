"""Unit tests for the variant layer (myvariant and variant_layer modules).

Tests run offline against saved fixtures and mocks. No network calls.
"""

from __future__ import annotations

import json
from pathlib import Path
import unittest
from unittest.mock import patch

from tools.variant_profile import myvariant, variant_layer, variant_profile

FIXTURES_DIR = Path(__file__).parent / "fixtures"


class TestMyVariantParsing(unittest.TestCase):
    def setUp(self):
        fixture_path = FIXTURES_DIR / "myvariant_CFTR_R31L.json"
        self.raw_data = json.loads(fixture_path.read_text(encoding="utf-8"))
        self.hit = self.raw_data["hits"][0]

    def test_parse_tied_fields(self):
        parsed = myvariant.parse_myvariant_record(self.hit)

        # 1. Genomic ID
        val, src = parsed["variant_id_genomic"]
        self.assertEqual(val, "chr7:g.117144345G>T")
        self.assertEqual(src, "[myvariant:_id]")

        # 2. CADD Phred
        val, src = parsed["cadd_phred"]
        self.assertEqual(val, 22.9)
        self.assertEqual(src, "[myvariant:cadd.phred]")

        # 3. CADD rawscore
        val, src = parsed["cadd_rawscore"]
        self.assertEqual(val, 3.36173)
        self.assertEqual(src, "[myvariant:cadd.rawscore]")

        # 4. PolyPhen-2
        val_cat, src_cat = parsed["polyphen_category"]
        val_score, _ = parsed["polyphen_score"]
        self.assertEqual(val_cat, "benign")
        self.assertEqual(val_score, 0.131)
        self.assertEqual(src_cat, "[myvariant:cadd.polyphen.cat]")

        # 5. SIFT
        val_cat, src_cat = parsed["sift_category"]
        val_score, _ = parsed["sift_score"]
        self.assertEqual(val_cat, "tolerated")
        self.assertEqual(val_score, 0.13)
        self.assertEqual(src_cat, "[myvariant:cadd.sift.cat]")

        # 6. SnpEff
        val_eff, src_eff = parsed["snpeff_effect"]
        val_impact, src_impact = parsed["snpeff_putative_impact"]
        self.assertEqual(val_eff, "missense_variant")
        self.assertEqual(val_impact, "MODERATE")
        self.assertEqual(src_eff, "[myvariant:snpeff.ann.effect]")
        self.assertEqual(src_impact, "[myvariant:snpeff.ann.putative_impact]")

    def test_missing_fields_render_not_stated(self):
        empty_hit = {}
        parsed = myvariant.parse_myvariant_record(empty_hit)
        for key, (val, src) in parsed.items():
            self.assertEqual(val, variant_layer.NOT_STATED, f"Key {key} should be not stated")
            self.assertTrue(src.startswith("[myvariant:"), f"Key {key} should have source tag")


class TestCrossCheckAndAlphaMissense(unittest.TestCase):
    def test_cross_check_agree_passes(self):
        with patch.object(variant_profile, "run") as mock_cv, patch.object(myvariant, "query_myvariant") as mock_mv:
            mock_cv.return_value = ({"match_status": "single"}, [{"clinvar_variation_id": 54087}])
            mock_mv.return_value = {"total": 1, "hits": [{"clinvar": {"variant_id": 54087}}]}

            res = variant_layer.run_variant_layer("CFTR", "R31L")
            cc = res["cross_checks"][0]
            self.assertEqual(cc["verdict"], "PASS")
            self.assertIn("54087", cc["detail"])

    def test_cross_check_disagree_flags_check(self):
        # Negative control: deliberately mismatched IDs must flag CHECK
        with patch.object(variant_profile, "run") as mock_cv, patch.object(myvariant, "query_myvariant") as mock_mv:
            mock_cv.return_value = ({"match_status": "single"}, [{"clinvar_variation_id": 54087}])
            mock_mv.return_value = {"total": 1, "hits": [{"clinvar": {"variant_id": 99999}}]}

            res = variant_layer.run_variant_layer("CFTR", "R31L")
            cc = res["cross_checks"][0]
            self.assertEqual(cc["verdict"], "CHECK")
            self.assertIn("differs", cc["detail"])

    def test_alphamissense_catalog_boundary_enforced(self):
        with patch.object(variant_profile, "run") as mock_cv, patch.object(myvariant, "query_myvariant") as mock_mv:
            mock_cv.return_value = ({"match_status": "single"}, [{}])
            mock_mv.return_value = {"total": 1, "hits": [{}]}

            res = variant_layer.run_variant_layer("CFTR", "R31L")
            am = res["alphamissense"]
            self.assertEqual(am["catalog_id"], "alphamissense-predictions")
            self.assertEqual(am["access"], "manual")
            self.assertFalse(am["admitted_for_api"])
            self.assertIn("manual", am["detail"])


class TestInputSafety(unittest.TestCase):
    def test_patient_like_input_refused(self):
        with self.assertRaises(variant_profile.Refused):
            variant_profile.validate("CFTR", "JohnDoe")
        with self.assertRaises(variant_profile.Refused):
            variant_profile.validate("CFTR", "1990-01-01")


if __name__ == "__main__":
    unittest.main()
