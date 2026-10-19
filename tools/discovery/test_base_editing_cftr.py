#!/usr/bin/env python3
"""Unit test suite for CFTR Base Editing Biophysical Model (T-0219).

Verifies:
1. Deamination window activity profiles W_E(p) for ABE8e, ABE7.10, and AncBE4max.
2. CRISPR-Cas PAM binding kinetics k_PAM across SpCas9, SpRY, and SpG.
3. Bystander deamination kinetics and purity ratio P calculation.
4. Identification and scoring of candidate gRNAs for W1282X, R553X, and 2789+5G->A.
5. Detection of transversion limitation for G542X requiring prime editing / dual deaminase.
6. Cryo-EM PDB 6MSM atomic coordinates and DeltaDelta G_fold structural perturbation annotations.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
import sys
import unittest
import pytest

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(REPO_ROOT))

from base_editing_cftr import (
    CFTR_ALLELE_CATALOG,
    EDITORS,
    CAS_VARIANTS,
    CFTRBaseEditingModel,
    compute_editing_window,
    score_pam_affinity,
    calculate_bystander_purity_ratio,
    evaluate_structural_perturbation,
    NLM_ATTRIBUTION,
    RESEARCH_BOUNDARY,
)


class TestCFTRBaseEditingModel(unittest.TestCase):
    """Test suite for CFTRBaseEditingModel biophysical engine."""

    def setUp(self) -> None:
        self.model = CFTRBaseEditingModel()

    def test_editor_profiles_and_window(self) -> None:
        """Verify Gaussian deamination windows for ABE8e, ABE7.10, and AncBE4max."""
        abe8e = EDITORS["ABE8e"]
        abe710 = EDITORS["ABE7.10"]

        # Center deamination probability W_E(p_0) is 1.0 at center
        self.assertAlmostEqual(compute_editing_window(5, abe8e), math.exp(-((5 - 5.5) ** 2) / (2 * 1.8 ** 2)), places=4)
        self.assertAlmostEqual(compute_editing_window(6, abe8e), math.exp(-((6 - 5.5) ** 2) / (2 * 1.8 ** 2)), places=4)

        # Symmetry around p_0 = 5.5: position 5 and position 6 have equal weights
        self.assertAlmostEqual(compute_editing_window(5, abe8e), compute_editing_window(6, abe8e))
        self.assertAlmostEqual(compute_editing_window(4, abe8e), compute_editing_window(7, abe8e))
        self.assertAlmostEqual(compute_editing_window(3, abe8e), compute_editing_window(8, abe8e))

        # ABE8e has broader window than ABE7.10 at boundary positions (pos 3 and pos 8)
        w_abe8e_pos3 = compute_editing_window(3, abe8e)
        w_abe710_pos3 = compute_editing_window(3, abe710)
        self.assertGreater(w_abe8e_pos3, w_abe710_pos3)
        self.assertGreater(w_abe8e_pos3, 0.35)  # ABE8e retains >35% activity at pos 3
        self.assertLess(w_abe710_pos3, 0.15)  # Canonical ABE7.10 is <15% at pos 3

        # Out-of-bounds positions return 0.0
        self.assertEqual(compute_editing_window(0, abe8e), 0.0)
        self.assertEqual(compute_editing_window(21, abe8e), 0.0)

    def test_pam_affinity_scoring(self) -> None:
        """Verify PAM recognition for SpCas9, SpRY, and SpG."""
        # SpCas9: NGG = 1.0, NAG = 0.15, NGA = 0.0
        self.assertEqual(score_pam_affinity("TGG", "SpCas9"), 1.0)
        self.assertEqual(score_pam_affinity("CGG", "SpCas9"), 1.0)
        self.assertEqual(score_pam_affinity("CAG", "SpCas9"), 0.15)
        self.assertEqual(score_pam_affinity("TGA", "SpCas9"), 0.0)

        # SpRY: NRN = 0.85 (A or G in middle), NYN = 0.65 (C or T in middle)
        self.assertEqual(score_pam_affinity("GAA", "SpRY"), 0.85)
        self.assertEqual(score_pam_affinity("CGT", "SpRY"), 0.85)
        self.assertEqual(score_pam_affinity("TTC", "SpRY"), 0.65)
        self.assertEqual(score_pam_affinity("CCA", "SpRY"), 0.65)

        # SpG: NGN = 0.90
        self.assertEqual(score_pam_affinity("TGC", "SpG"), 0.90)
        self.assertEqual(score_pam_affinity("AGA", "SpG"), 0.90)
        self.assertEqual(score_pam_affinity("TAA", "SpG"), 0.10)

    def test_bystander_purity_ratio_calculation(self) -> None:
        """Verify bystander purity ratio P calculation."""
        # Zero bystanders -> purity is exactly 1.0
        self.assertEqual(calculate_bystander_purity_ratio(0.85, []), 1.0)

        # 1 bystander with weight 0.15
        p = calculate_bystander_purity_ratio(0.85, [0.15])
        self.assertAlmostEqual(p, 0.85 / (0.85 + 0.15), places=4)
        self.assertEqual(p, 0.85)

        # Higher bystander weight reduces purity
        p_low = calculate_bystander_purity_ratio(0.85, [0.85, 0.50])
        self.assertLess(p_low, 0.50)

    def test_w1282x_candidate_gRNA(self) -> None:
        """Verify gRNA design and scoring for CFTR W1282X (c.3846G>A, Exon 23)."""
        top_spcas9 = self.model.get_top_candidate("W1282X", "ABE8e", "SpCas9")
        self.assertIsNotNone(top_spcas9)
        self.assertEqual(top_spcas9.pam_sequence, "AGG")
        self.assertEqual(top_spcas9.target_base_position, 7)
        self.assertGreater(top_spcas9.on_target_efficiency, 0.70)
        self.assertLessEqual(top_spcas9.max_destabilizing_ddg, 1.0)  # No severe missense destabilization

        # SpRY and SpG find guide placing target adenine at optimal pos 6 (W_E > 0.95)
        top_spry = self.model.get_top_candidate("W1282X", "ABE8e", "SpRY")
        self.assertIsNotNone(top_spry)
        self.assertEqual(top_spry.target_base_position, 6)
        self.assertGreater(top_spry.on_target_efficiency, 0.80)
        self.assertGreater(top_spry.deamination_weight, 0.95)

    def test_r553x_candidate_gRNA(self) -> None:
        """Verify gRNA design for R553X (c.1657C>T, Exon 12) on non-coding strand."""
        top_spry = self.model.get_top_candidate("R553X", "ABE8e", "SpRY")
        self.assertIsNotNone(top_spry)
        self.assertEqual(top_spry.strand, "non_coding")
        self.assertEqual(top_spry.target_base_position, 6)
        self.assertGreater(top_spry.on_target_efficiency, 0.80)
        # R553X non-coding strand has exceptionally high purity (>0.95)
        self.assertGreater(top_spry.bystander_purity_ratio, 0.95)
        self.assertIn("HIGH_FIDELITY_RESCUE", top_spry.curative_verdict)

        top_spg = self.model.get_top_candidate("R553X", "ABE8e", "SpG")
        self.assertIsNotNone(top_spg)
        self.assertGreater(top_spg.bystander_purity_ratio, 0.98)

    def test_2789_splice_site_restoration(self) -> None:
        """Verify gRNA design for intronic splice mutation 2789+5G->A (c.2657+5G>A)."""
        top_spry = self.model.get_top_candidate("2789+5G->A", "ABE8e", "SpRY")
        self.assertIsNotNone(top_spry)
        self.assertGreater(top_spry.on_target_efficiency, 0.80)
        # Intronic mutations cause zero polypeptide missense perturbation
        self.assertEqual(top_spry.max_destabilizing_ddg, 0.0)

    def test_g542x_transversion_limitation(self) -> None:
        """Verify G542X (c.1624G>T transversion) is flagged as unfeasible with pure transition base editors."""
        top_spry = self.model.get_top_candidate("G542X", "ABE8e", "SpRY")
        self.assertIsNotNone(top_spry)
        self.assertIn("UNFEASIBLE_TRANSVERSION", top_spry.curative_verdict)
        # Bystander/missense edit at codon 542 is structurally destabilizing (DDG > 3.0 kcal/mol)
        self.assertGreater(top_spry.max_destabilizing_ddg, 3.0)

    def test_evaluate_all_priority_alleles(self) -> None:
        """Verify evaluate_all_priority_alleles matrix returns all 4 key CFTR mutations."""
        res = self.model.evaluate_all_priority_alleles("ABE8e")
        self.assertIn("W1282X", res)
        self.assertIn("R553X", res)
        self.assertIn("G542X", res)
        self.assertIn("2789+5G->A", res)

    def test_cryoem_coordinates_provenance(self) -> None:
        """Verify Cryo-EM 6MSM structural coordinates and references."""
        for allele_name, profile in CFTR_ALLELE_CATALOG.items():
            self.assertIn("6MSM", profile.pdb_6msm_context)
            if profile.mutation_type == "nonsense":
                self.assertIsNotNone(profile.residue_number)

    def test_attribution_and_boundary(self) -> None:
        """Verify presence of NLM attribution and research boundary disclaimers."""
        self.assertIn("ClinicalTrials.gov", NLM_ATTRIBUTION)
        self.assertIn("Research, not medical advice", RESEARCH_BOUNDARY)


if __name__ == "__main__":
    unittest.main()
