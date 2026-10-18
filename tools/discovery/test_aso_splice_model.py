#!/usr/bin/env python3
"""Deterministic verification test suite for aso_splice_model.py (T-0183).

Dual compatibility:
- pytest tools/discovery/test_aso_splice_model.py
- python -m unittest discover -s tools/discovery -p "test_*.py"
"""
from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))

import aso_splice_model as aso


class TestASOHybridizationThermodynamics(unittest.TestCase):
    """Verifies nearest-neighbor hybridization thermodynamics and chemistry modifications."""

    def test_moe_vs_rna_stabilization(self):
        seq = "AAAGGGUAAGUAGCAU"
        res_rna = aso.calculate_hybridization_thermodynamics(seq, chemistry="RNA")
        res_moe = aso.calculate_hybridization_thermodynamics(seq, chemistry="2MOE")

        # 2MOE should have significantly more negative Delta G and higher Tm
        self.assertLess(res_moe.delta_g_37_kcal, res_rna.delta_g_37_kcal)
        self.assertGreater(res_moe.melting_temperature_c, res_rna.melting_temperature_c)
        self.assertEqual(res_moe.chemistry, "2MOE")
        self.assertEqual(res_rna.chemistry, "RNA")

    def test_pmo_vs_rna_stabilization(self):
        seq = "AAAGGGUAAGUAGCAU"
        res_rna = aso.calculate_hybridization_thermodynamics(seq, chemistry="RNA")
        res_pmo = aso.calculate_hybridization_thermodynamics(seq, chemistry="PMO")

        # PMO uncharged backbone confers enhanced affinity
        self.assertLess(res_pmo.delta_g_37_kcal, res_rna.delta_g_37_kcal)
        self.assertGreater(res_pmo.melting_temperature_c, res_rna.melting_temperature_c)
        self.assertEqual(res_pmo.chemistry, "PMO")

    def test_pmo_salt_independence(self):
        seq = "AAAGGGUAAGUAGCAU"
        # Low salt vs physiological salt
        res_moe_low = aso.calculate_hybridization_thermodynamics(seq, chemistry="2MOE", na_m=0.020)
        res_moe_high = aso.calculate_hybridization_thermodynamics(seq, chemistry="2MOE", na_m=0.150)
        delta_s_moe_shift = abs(res_moe_high.delta_s_cal_k - res_moe_low.delta_s_cal_k)

        res_pmo_low = aso.calculate_hybridization_thermodynamics(seq, chemistry="PMO", na_m=0.020)
        res_pmo_high = aso.calculate_hybridization_thermodynamics(seq, chemistry="PMO", na_m=0.150)
        delta_s_pmo_shift = abs(res_pmo_high.delta_s_cal_k - res_pmo_low.delta_s_cal_k)

        # PMO salt shift should be much smaller than charged MOE
        self.assertLess(delta_s_pmo_shift, delta_s_moe_shift)

    def test_gc_content_thermodynamic_effect(self):
        seq_at = "AUAUAUAUAUAUAUAU"
        seq_gc = "GCGCGCGCGCGCGCGC"
        res_at = aso.calculate_hybridization_thermodynamics(seq_at, chemistry="RNA")
        res_gc = aso.calculate_hybridization_thermodynamics(seq_gc, chemistry="RNA")

        self.assertGreater(res_gc.gc_content_pct, res_at.gc_content_pct)
        self.assertGreater(res_gc.melting_temperature_c, res_at.melting_temperature_c)
        self.assertLess(res_gc.delta_g_37_kcal, res_at.delta_g_37_kcal)

    def test_invalid_chemistry_or_sequence(self):
        with self.assertRaises(ValueError):
            aso.calculate_hybridization_thermodynamics("AUGC", chemistry="INVALID")
        with self.assertRaises(ValueError):
            aso.calculate_hybridization_thermodynamics("XYZ123", chemistry="2MOE")


class TestASODoseResponseAndSplicing(unittest.TestCase):
    """Verifies competitive spliceosome blocking partition function and Hill dose-response."""

    def setUp(self):
        self.summary_3849 = aso.evaluate_aso_splice_switching_design("3849+10kb C->T", chemistry="2MOE")

    def test_dose_response_monotonicity(self):
        points = self.summary_3849.dose_response
        for i in range(len(points) - 1):
            p_curr = points[i]
            p_next = points[i + 1]
            # Concentration increases
            self.assertLessEqual(p_curr.concentration_nm, p_next.concentration_nm)
            # Bound fraction increases monotonically
            self.assertLessEqual(p_curr.fraction_target_bound, p_next.fraction_target_bound)
            # Canonical transcript increases monotonically
            self.assertLessEqual(p_curr.fraction_canonical_transcript, p_next.fraction_canonical_transcript)
            # Chloride current increases monotonically
            self.assertLessEqual(p_curr.restored_chloride_current_pct_wt, p_next.restored_chloride_current_pct_wt)

    def test_pseudoexon_exclusion_monotonicity(self):
        points = self.summary_3849.dose_response
        for i in range(len(points) - 1):
            p_curr = points[i]
            p_next = points[i + 1]
            # Pseudoexon inclusion fraction decreases monotonically
            self.assertGreaterEqual(p_curr.fraction_pseudoexon_included, p_next.fraction_pseudoexon_included)

    def test_zero_dose_matches_baseline(self):
        zero_pt = self.summary_3849.dose_response[0]
        self.assertEqual(zero_pt.concentration_nm, 0.0)
        self.assertEqual(zero_pt.fraction_target_bound, 0.0)
        self.assertEqual(zero_pt.fraction_pseudoexon_included, 0.85)
        self.assertEqual(zero_pt.fraction_canonical_transcript, 0.15)
        self.assertEqual(zero_pt.restored_chloride_current_pct_wt, 6.5)
        self.assertFalse(zero_pt.meets_therapeutic_threshold)

    def test_saturating_dose_reaches_plateau(self):
        sat_pt = self.summary_3849.dose_response[-1]
        self.assertEqual(sat_pt.concentration_nm, 1000.0)
        self.assertGreaterEqual(sat_pt.fraction_target_bound, 0.98)
        self.assertLessEqual(sat_pt.fraction_pseudoexon_included, 0.05)
        self.assertGreaterEqual(sat_pt.restored_chloride_current_pct_wt, 45.0)
        self.assertTrue(sat_pt.meets_therapeutic_threshold)

    def test_therapeutic_threshold_crossing(self):
        # Meets therapeutic threshold (>= 10% WT current) at achievable nM dose
        self.assertIsNotNone(self.summary_3849.minimal_curative_dose_nm)
        self.assertLessEqual(self.summary_3849.minimal_curative_dose_nm, 25.0)
        self.assertEqual(self.summary_3849.verdict, "POTENT_CURATIVE_RESCUE")


class TestSpliceVariantCatalogAndGrounding(unittest.TestCase):
    """Verifies variant catalog completeness and strict anti-hallucination identifier invariants."""

    def test_catalog_entries_present(self):
        catalog = aso.SPLICE_VARIANT_CATALOG
        self.assertIn("3849+10kb C->T", catalog)
        self.assertIn("2789+5G->A", catalog)
        self.assertIn("1811+1.6kbA->G", catalog)

    def test_variant_3849_metadata(self):
        target = aso.SPLICE_VARIANT_CATALOG["3849+10kb C->T"]
        self.assertEqual(target.cdna_nomenclature, "c.3718-2477C>T")
        self.assertEqual(target.clinical_trial_id, "NCT06429176")
        self.assertIn("RCSB PDB 6MSM", target.evidence_citations)

    def test_all_variants_evaluable(self):
        for vname in aso.SPLICE_VARIANT_CATALOG:
            res = aso.evaluate_aso_splice_switching_design(vname)
            self.assertGreater(res.predicted_max_cftr_current_pct_wt, 10.0)
            self.assertIn(res.verdict, ("POTENT_CURATIVE_RESCUE", "VIABLE_THERAPEUTIC_RESCUE"))

    def test_anti_hallucination_provenance_guardrails(self):
        # AH-1: All cited NCTs must exist in snapshot.json
        snapshot_file = REPO_ROOT / "tools" / "trial_atlas" / "snapshot.json"
        self.assertTrue(snapshot_file.exists())
        with open(snapshot_file, encoding="utf-8") as f:
            snap_data = json.load(f)
        valid_ncts = {s["nct_id"] for s in snap_data.get("studies", []) if "nct_id" in s}

        for vname, target in aso.SPLICE_VARIANT_CATALOG.items():
            if target.clinical_trial_id:
                self.assertIn(target.clinical_trial_id, valid_ncts)
            for cit in target.evidence_citations:
                if cit.startswith("NCT"):
                    self.assertIn(cit, valid_ncts)

        # Prohibited clinical overclaim check
        source_code = (HERE / "aso_splice_model.py").read_text(encoding="utf-8").lower()
        prohibited = ["proven to cure", "guarantees clinical benefit", "recommended treatment regimen"]
        for p in prohibited:
            self.assertNotIn(p, source_code)


class TestCLIExecution(unittest.TestCase):
    """Verifies command-line execution, self-test, and JSON export."""

    def test_cli_self_test(self):
        cmd = [sys.executable, str(HERE / "aso_splice_model.py"), "--self-test"]
        proc = subprocess.run(cmd, capture_output=True, text=True, check=True)
        self.assertEqual(proc.returncode, 0)
        self.assertIn("PASS: aso_splice_model self-test passed.", proc.stdout)

    def test_cli_json_export(self):
        export_path = HERE / "scratch_test_aso_export.json"
        try:
            cmd = [
                sys.executable,
                str(HERE / "aso_splice_model.py"),
                "--variant",
                "all",
                "--export",
                str(export_path),
                "--json",
            ]
            proc = subprocess.run(cmd, capture_output=True, text=True, check=True)
            self.assertEqual(proc.returncode, 0)
            self.assertTrue(export_path.exists())

            data = json.loads(export_path.read_text(encoding="utf-8"))
            self.assertIn("3849+10kb C->T", data)
            self.assertIn("2789+5G->A", data)
            self.assertIn("1811+1.6kbA->G", data)
            self.assertEqual(data["3849+10kb C->T"]["verdict"], "POTENT_CURATIVE_RESCUE")
        finally:
            if export_path.exists():
                export_path.unlink()


if __name__ == "__main__":
    unittest.main()
