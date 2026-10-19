#!/usr/bin/env python3
"""Unit tests for CFTR Delta-R Minigene Packaging & Patch-Clamp Simulator (T-0207)."""
from __future__ import annotations

import json
import math
import subprocess
import sys
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(REPO_ROOT))

from cftr_minigene_model import (
    AAV_HARD_PACKAGING_LIMIT_NT,
    DEFAULT_SNAPSHOT_JSON,
    NLM_ATTRIBUTION,
    R_DOMAIN_END_RESIDUE,
    R_DOMAIN_START_RESIDUE,
    RESEARCH_BOUNDARY,
    UNITARY_CONDUCTANCE_PS,
    evaluate_minigene_packaging,
    self_test,
    simulate_cftr_minigene,
    simulate_patch_clamp_electrophysiology,
)


class TestCFTRMinigeneModel(unittest.TestCase):
    def test_packaging_capacity_deltar_vs_wt(self):
        """DeltaR(708-759) with compact elements fits single AAV, while WT CFTR exceeds capacity."""
        pkg_dr = evaluate_minigene_packaging("CFTR_DeltaR_708_759", promoter_key="tg83", polya_key="synPolyA_compact")
        self.assertTrue(pkg_dr.fits_single_aav)
        self.assertLessEqual(pkg_dr.total_genome_nt, AAV_HARD_PACKAGING_LIMIT_NT)
        self.assertGreaterEqual(pkg_dr.packaging_headroom_nt, 0)

        pkg_wt = evaluate_minigene_packaging("CFTR_WT", promoter_key="tg83", polya_key="synPolyA_compact")
        self.assertFalse(pkg_wt.fits_single_aav)
        self.assertGreater(pkg_wt.total_genome_nt, AAV_HARD_PACKAGING_LIMIT_NT)
        self.assertLess(pkg_wt.packaging_headroom_nt, 0)

    def test_promoter_size_impact_on_packaging(self):
        """DeltaR(708-759) fits with compact promoters (tg83) but fails with large promoters (CBA_full)."""
        pkg_tg83 = evaluate_minigene_packaging("CFTR_DeltaR_708_759", promoter_key="tg83", polya_key="synPolyA_compact")
        self.assertTrue(pkg_tg83.fits_single_aav)

        pkg_cba = evaluate_minigene_packaging("CFTR_DeltaR_708_759", promoter_key="CBA_full", polya_key="synPolyA_compact")
        self.assertFalse(pkg_cba.fits_single_aav)
        self.assertGreater(pkg_cba.total_genome_nt, AAV_HARD_PACKAGING_LIMIT_NT)

    def test_unitary_conductance_invariance(self):
        """Unitary single-channel conductance gamma = 8.5 pS is invariant between WT and DeltaR designs."""
        ep_wt = simulate_patch_clamp_electrophysiology("CFTR_WT")
        ep_dr = simulate_patch_clamp_electrophysiology("CFTR_DeltaR_708_759")
        ep_ext = simulate_patch_clamp_electrophysiology("CFTR_DeltaR_708_835")

        self.assertAlmostEqual(ep_wt.unitary_conductance_ps, 8.5)
        self.assertAlmostEqual(ep_dr.unitary_conductance_ps, 8.5)
        self.assertAlmostEqual(ep_ext.unitary_conductance_ps, 8.5)

    def test_atp_gating_requirement(self):
        """Strict ATP dependence: zero ATP yields Po = 0.0 regardless of PKA phosphorylation."""
        ep_no_atp_stim = simulate_patch_clamp_electrophysiology("CFTR_DeltaR_708_759", pka_stimulated=True, atp_concentration_mm=0.0)
        ep_no_atp_basal = simulate_patch_clamp_electrophysiology("CFTR_DeltaR_708_759", pka_stimulated=False, atp_concentration_mm=0.0)

        self.assertEqual(ep_no_atp_stim.open_probability_po, 0.0)
        self.assertEqual(ep_no_atp_basal.open_probability_po, 0.0)
        self.assertEqual(ep_no_atp_stim.macroscopic_current_pa, 0.0)

    def test_atp_concentration_dependence(self):
        """Po increases monotonically with ATP concentration."""
        ep_low = simulate_patch_clamp_electrophysiology("CFTR_DeltaR_708_759", atp_concentration_mm=0.05)
        ep_mid = simulate_patch_clamp_electrophysiology("CFTR_DeltaR_708_759", atp_concentration_mm=0.5)
        ep_high = simulate_patch_clamp_electrophysiology("CFTR_DeltaR_708_759", atp_concentration_mm=2.0)

        self.assertLess(ep_low.open_probability_po, ep_mid.open_probability_po)
        self.assertLess(ep_mid.open_probability_po, ep_high.open_probability_po)

    def test_basal_vs_stimulated_po(self):
        """DeltaR(708-759) exhibits elevated basal Po (~0.075) vs WT (~0.008), with full stimulation ~0.42."""
        ep_wt_basal = simulate_patch_clamp_electrophysiology("CFTR_WT", pka_stimulated=False, atp_concentration_mm=1.0)
        ep_dr_basal = simulate_patch_clamp_electrophysiology("CFTR_DeltaR_708_759", pka_stimulated=False, atp_concentration_mm=1.0)
        ep_dr_stim = simulate_patch_clamp_electrophysiology("CFTR_DeltaR_708_759", pka_stimulated=True, atp_concentration_mm=1.0)

        self.assertLess(ep_wt_basal.open_probability_po, 0.02)
        self.assertGreater(ep_dr_basal.open_probability_po, 0.05)
        self.assertGreater(ep_dr_stim.open_probability_po, ep_dr_basal.open_probability_po)
        self.assertGreaterEqual(ep_dr_stim.open_probability_po, 0.38)

    def test_delta264_truncation_er_failure(self):
        """CFTR_Delta264 N-terminal truncation fails ER export (f_Band_C < 0.05) and has zero current."""
        ep_d264 = simulate_patch_clamp_electrophysiology("CFTR_Delta264")
        self.assertLess(ep_d264.fraction_band_c_mature, 0.05)
        self.assertEqual(ep_d264.open_probability_po, 0.0)
        self.assertEqual(ep_d264.macroscopic_current_pa, 0.0)

    def test_driving_force_and_current_calculation(self):
        """Verify driving force (Vm - E_Cl) and single-channel current i = gamma * (Vm - E_Cl)."""
        ep = simulate_patch_clamp_electrophysiology(
            "CFTR_DeltaR_708_759",
            v_membrane_mv=-60.0,
            e_cl_mv=-20.0,
        )
        self.assertEqual(ep.driving_force_mv, -40.0)
        # i = 8.5 pS * -40 mV = 8.5e-12 S * -0.040 V = -0.34 pA
        self.assertAlmostEqual(ep.single_channel_current_pa, -0.34, places=2)

    def test_r_domain_structural_boundaries(self):
        """Deleted residues 708-759 lie strictly within cryo-EM disordered R-domain (654-831)."""
        del_start, del_end = 708, 759
        self.assertGreaterEqual(del_start, R_DOMAIN_START_RESIDUE)
        self.assertLessEqual(del_end, R_DOMAIN_END_RESIDUE)

    def test_nct_provenance_guard(self):
        """Rule AH-1: Ensure NCT05248230 is present in tools/trial_atlas/snapshot.json."""
        if not DEFAULT_SNAPSHOT_JSON.is_file():
            self.skipTest("snapshot.json not found")
        data = json.loads(DEFAULT_SNAPSHOT_JSON.read_text(encoding="utf-8"))
        snapshot_ncts = {s["nct_id"] for s in data.get("studies", []) if "nct_id" in s}
        self.assertIn("NCT05248230", snapshot_ncts)

    def test_no_clinical_overclaims(self):
        """Rule AH-4: Verify zero prohibited clinical overclaim phrases."""
        prohibited = (
            "proven to cure",
            "guarantees clinical benefit",
            "recommended treatment regimen",
            "safe and effective for patients",
            "cures cystic fibrosis",
        )
        source = Path(__file__).resolve().parent.joinpath("cftr_minigene_model.py").read_text(encoding="utf-8")
        for p in prohibited:
            self.assertNotIn(p.lower(), source.lower())

    def test_disclaimers_present(self):
        """Verify NLM attribution and research boundary texts are present."""
        self.assertIn("National Library of Medicine", NLM_ATTRIBUTION)
        self.assertIn("Research, not medical advice", RESEARCH_BOUNDARY)

    def test_internal_self_test(self):
        """Verify self_test() returns True."""
        self.assertTrue(self_test())

    def test_cli_execution_and_json(self):
        """Verify CLI execution with --self-test and --json."""
        cmd_self = [sys.executable, "-m", "tools.discovery.cftr_minigene_model", "--self-test"]
        proc_self = subprocess.run(cmd_self, capture_output=True, text=True, cwd=str(REPO_ROOT))
        self.assertEqual(proc_self.returncode, 0)
        self.assertIn("PASS", proc_self.stdout)

        cmd_json = [sys.executable, "-m", "tools.discovery.cftr_minigene_model", "--construct", "CFTR_DeltaR_708_759", "--json"]
        proc_json = subprocess.run(cmd_json, capture_output=True, text=True, cwd=str(REPO_ROOT))
        self.assertEqual(proc_json.returncode, 0)
        parsed = json.loads(proc_json.stdout)
        self.assertIn("packaging", parsed)
        self.assertIn("electrophysiology", parsed)
        self.assertEqual(parsed["packaging"]["construct_key"], "CFTR_DeltaR_708_759")


if __name__ == "__main__":
    unittest.main()
