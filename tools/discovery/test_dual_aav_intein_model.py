#!/usr/bin/env python3
"""Unit tests for Dual-AAV Split-Intein Trans-Splicing CFTR Biophysical Reconstitution Model (T-0206)."""
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

from dual_aav_intein_model import (
    AAV_PACKAGING_HARD_LIMIT_NT,
    DEFAULT_SNAPSHOT_JSON,
    INTEIN_SYSTEMS,
    NLM_ATTRIBUTION,
    RESEARCH_BOUNDARY,
    calculate_poisson_co_transduction,
    check_single_aav_packaging_capacity,
    design_split_cftr_vectors,
    self_test,
    simulate_dual_aav_cftr_reconstitution,
    simulate_second_order_trans_splicing_kinetics,
)


class TestDualAAVInteinModel(unittest.TestCase):
    def test_single_aav_packaging_failure(self):
        """Single-AAV full-length CFTR cDNA + expression elements must exceed 4,700 nt."""
        budget = check_single_aav_packaging_capacity()
        self.assertFalse(budget.fits_packaging_limit)
        self.assertGreater(budget.total_genome_nt, AAV_PACKAGING_HARD_LIMIT_NT)
        self.assertGreaterEqual(budget.total_genome_nt, 5200)
        self.assertGreater(budget.payload_excess_nt, 500)

    def test_split_aav_packaging_within_limits(self):
        """Both split vectors (N and C halves) must fit strictly within 4,700 nt."""
        design = design_split_cftr_vectors(split_residue=744, intein_key="Cfa", junction_c_plus1_aa="S")
        self.assertTrue(design.n_vector.fits_packaging_limit)
        self.assertTrue(design.c_vector.fits_packaging_limit)
        self.assertLessEqual(design.n_vector.total_genome_nt, AAV_PACKAGING_HARD_LIMIT_NT)
        self.assertLessEqual(design.c_vector.total_genome_nt, AAV_PACKAGING_HARD_LIMIT_NT)
        self.assertEqual(design.n_vector.payload_excess_nt, 0)
        self.assertEqual(design.c_vector.payload_excess_nt, 0)

    def test_r_domain_disordered_split_targeting(self):
        """Verify split sites inside disordered R-domain (654-831) vs folded domains."""
        # Residue 744, 783, 800 are inside R-domain
        for res in (700, 744, 783, 820):
            d = design_split_cftr_vectors(split_residue=res, intein_key="Cfa", junction_c_plus1_aa="S")
            self.assertTrue(d.split_in_disordered_r_domain, f"Residue {res} should be in R-domain")
            self.assertGreater(d.predicted_splicing_efficiency, 0.30)

        # Splitting in NBD1 (e.g. 508) or TMD2 (e.g. 900) should be flagged and penalized
        d_nbd1 = design_split_cftr_vectors(split_residue=508, intein_key="Cfa", junction_c_plus1_aa="S")
        self.assertFalse(d_nbd1.split_in_disordered_r_domain)
        self.assertLess(d_nbd1.predicted_splicing_efficiency, 0.15)

    def test_catalytic_plus1_junction_rules(self):
        """Verify C-extein +1 residue constraints (Cfa accepts C/S, Npu accepts C only)."""
        # Cfa accepts C and S
        d_cfa_s = design_split_cftr_vectors(split_residue=744, intein_key="Cfa", junction_c_plus1_aa="S")
        d_cfa_c = design_split_cftr_vectors(split_residue=744, intein_key="Cfa", junction_c_plus1_aa="C")
        self.assertTrue(d_cfa_s.is_valid_junction)
        self.assertTrue(d_cfa_c.is_valid_junction)

        # Cfa rejects Proline or Glycine
        d_cfa_p = design_split_cftr_vectors(split_residue=744, intein_key="Cfa", junction_c_plus1_aa="P")
        self.assertFalse(d_cfa_p.is_valid_junction)
        self.assertEqual(d_cfa_p.predicted_splicing_efficiency, 0.0)

        # Npu DnaE strictly requires C
        d_npu_c = design_split_cftr_vectors(split_residue=744, intein_key="Npu_DnaE", junction_c_plus1_aa="C")
        self.assertTrue(d_npu_c.is_valid_junction)
        d_npu_s = design_split_cftr_vectors(split_residue=744, intein_key="Npu_DnaE", junction_c_plus1_aa="S")
        self.assertFalse(d_npu_s.is_valid_junction)
        self.assertEqual(d_npu_s.predicted_splicing_efficiency, 0.0)

    def test_kinetics_cfa_vs_npu_dnae(self):
        """Verify Cfa has higher k_splice and shorter half-life than Npu DnaE."""
        k_cfa = INTEIN_SYSTEMS["Cfa"]["k_splice_s_inv"]
        k_npu = INTEIN_SYSTEMS["Npu_DnaE"]["k_splice_s_inv"]
        t_cfa = INTEIN_SYSTEMS["Cfa"]["t_half_seconds"]
        t_npu = INTEIN_SYSTEMS["Npu_DnaE"]["t_half_seconds"]

        self.assertGreater(k_cfa, k_npu)
        self.assertLess(t_cfa, t_npu)
        self.assertAlmostEqual(k_cfa, 0.015, places=4)
        self.assertAlmostEqual(k_npu, 0.005, places=4)

    def test_second_order_kinetics_progress(self):
        """Verify time-dependent trans-splicing yield increases monotonically."""
        k0 = simulate_second_order_trans_splicing_kinetics(time_seconds=0.0, intein_key="Cfa")
        k10 = simulate_second_order_trans_splicing_kinetics(time_seconds=10.0, intein_key="Cfa")
        k60 = simulate_second_order_trans_splicing_kinetics(time_seconds=60.0, intein_key="Cfa")
        k300 = simulate_second_order_trans_splicing_kinetics(time_seconds=300.0, intein_key="Cfa")

        self.assertEqual(k0["fraction_spliced"], 0.0)
        self.assertLess(k0["fraction_spliced"], k10["fraction_spliced"])
        self.assertLess(k10["fraction_spliced"], k60["fraction_spliced"])
        self.assertLess(k60["fraction_spliced"], k300["fraction_spliced"])
        self.assertLess(k300["fraction_spliced"], 1.0)

    def test_poisson_co_transduction_scaling(self):
        """Verify independent Poisson entry: P_dual = (1 - e^-MOI)^2."""
        for moi in (0.5, 1.0, 2.0, 3.0, 5.0, 10.0):
            res = calculate_poisson_co_transduction(moi_n=moi, moi_c=moi)
            expected = (1.0 - math.exp(-moi)) ** 2
            self.assertAlmostEqual(res.p_dual_co_transduction, expected, places=3)
            self.assertAlmostEqual(res.percent_dual_positive_cells, expected * 100.0, places=1)

    def test_moi_zero_and_negative_handling(self):
        """Verify MOI 0 gives zero co-transduction, negative MOI raises ValueError."""
        res_zero = calculate_poisson_co_transduction(moi_n=0.0, moi_c=0.0)
        self.assertEqual(res_zero.p_dual_co_transduction, 0.0)
        with self.assertRaises(ValueError):
            calculate_poisson_co_transduction(moi_n=-1.0, moi_c=2.0)

    def test_restored_chloride_transport_matches_benchmark(self):
        """At MOI 3 with Cfa, restored transepithelial Cl- transport reaches ~42% WT (DOI 10.3389/fgeed.2023.1271813)."""
        sim = simulate_dual_aav_cftr_reconstitution(moi=3.0, intein_key="Cfa", junction_c_plus1_aa="S")
        self.assertAlmostEqual(sim.restored_chloride_transport_pct_wt, 42.0, delta=2.0)
        self.assertGreater(sim.reconstituted_mature_cftr_yield_pct_wt, 30.0)

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
        source = Path(__file__).resolve().parent.joinpath("dual_aav_intein_model.py").read_text(encoding="utf-8")
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
        cmd_self = [sys.executable, "-m", "tools.discovery.dual_aav_intein_model", "--self-test"]
        proc_self = subprocess.run(cmd_self, capture_output=True, text=True, cwd=str(REPO_ROOT))
        self.assertEqual(proc_self.returncode, 0)
        self.assertIn("PASS", proc_self.stdout)

        cmd_json = [sys.executable, "-m", "tools.discovery.dual_aav_intein_model", "--moi", "3.0", "--json"]
        proc_json = subprocess.run(cmd_json, capture_output=True, text=True, cwd=str(REPO_ROOT))
        self.assertEqual(proc_json.returncode, 0)
        parsed = json.loads(proc_json.stdout)
        self.assertIn("intein_design", parsed)
        self.assertIn("co_transduction", parsed)
        self.assertIn("reconstituted_mature_cftr_yield_pct_wt", parsed)


if __name__ == "__main__":
    unittest.main()
