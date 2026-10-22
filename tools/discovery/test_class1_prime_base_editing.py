#!/usr/bin/env python3
"""Unit tests for Class I Nonsense Allele Base Editing & Prime Editing Model (T-0255).

Verifies:
1. Trial provenance: 100% of cited NCT identifiers (NCT05668741, NCT06237335) verified in
   tools/trial_atlas/snapshot.json.
2. Canonical Class I allele catalog invariants (G542X, W1282X, R553X, R1162X) with Cryo-EM coordinates.
3. ABE8e Gaussian deamination window (positions 4-8) and PAM affinity profiles (SpCas9, SpCas9-NG, SpRY).
4. Direct wild-type restoration for W1282X (TGG Trp) and antisense deamination for R553X and R1162X (CGA Arg).
5. G542X transversion roadblock and requirement for Prime Editing (PE3).
6. Prime editing pegRNA design: PBS Tm thermodynamics, RTT flap hybridization, and 100% sequence fidelity.
7. Basal stem cell (KRT5+/TP63+) delivery dynamics, mucus barrier penetration, and 10-year durability.
8. Anti-hallucination compliance, NLM disclaimer, and research boundary enforcement.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path
import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
DISCOVERY_DIR = Path(__file__).resolve().parent

if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(DISCOVERY_DIR) not in sys.path:
    sys.path.insert(0, str(DISCOVERY_DIR))

from class1_prime_base_editing import (
    CLASS1_ALLELE_CATALOG,
    CURATIVE_THRESHOLD_PCT_WT,
    NLM_REGISTRY_DISCLAIMER,
    RESEARCH_BOUNDARY,
    compute_abe8e_activity_weight,
    compute_dna_tm_and_delta_g,
    design_prime_editing_allele,
    evaluate_class1_nonsense_landscape,
    get_pam_affinity,
    model_abe8e_allele,
    model_basal_stem_cell_delivery,
    reverse_complement,
)


class TestClass1PrimeBaseEditing:
    """Comprehensive test suite for Class I CFTR nonsense allele editing and stem cell delivery."""

    def test_01_trial_provenance_and_snapshot_grounding(self) -> None:
        """Rule AH-1: Every cited clinical trial must exist in tools/trial_atlas/snapshot.json."""
        snapshot_path = REPO_ROOT / "tools" / "trial_atlas" / "snapshot.json"
        assert snapshot_path.exists(), f"Snapshot file not found: {snapshot_path}"

        with open(snapshot_path, "r", encoding="utf-8") as f:
            snap = json.load(f)

        study_ncts = {s["nct_id"] for s in snap.get("studies", []) if "nct_id" in s}

        expected_ncts = ["NCT05668741", "NCT06237335"]
        for nct in expected_ncts:
            assert nct in study_ncts, f"NCT {nct} not grounded in snapshot.json"

    def test_02_class1_allele_catalog_invariants(self) -> None:
        """Verify the 4 canonical Class I nonsense alleles with biophysical and structural parameters."""
        assert len(CLASS1_ALLELE_CATALOG) == 4
        assert set(CLASS1_ALLELE_CATALOG.keys()) == {"G542X", "W1282X", "R553X", "R1162X"}

        for name, allele in CLASS1_ALLELE_CATALOG.items():
            assert allele.mutant_codon == "TGA", f"{name} must have UGA/TGA stop codon"
            assert allele.residue_number in (542, 1282, 553, 1162)
            assert len(allele.coding_context) >= 50, f"{name} context sequence too short"
            assert len(allele.cryoem_6msm_coords) == 3, f"{name} missing 3D Cryo-EM coordinates"

        # Check allele-specific details
        g542x = CLASS1_ALLELE_CATALOG["G542X"]
        assert g542x.cdna == "c.1624G>T"
        assert g542x.wt_codon == "GGA"
        assert g542x.transversion_roadblock is True

        w1282x = CLASS1_ALLELE_CATALOG["W1282X"]
        assert w1282x.cdna == "c.3846G>A"
        assert w1282x.wt_codon == "TGG"
        assert w1282x.transversion_roadblock is False

        r553x = CLASS1_ALLELE_CATALOG["R553X"]
        assert r553x.cdna == "c.1657C>T"
        assert r553x.wt_codon == "CGA"
        assert r553x.target_strand_abe == "non_coding"

        r1162x = CLASS1_ALLELE_CATALOG["R1162X"]
        assert r1162x.cdna == "c.3484C>T"
        assert r1162x.wt_codon == "CGA"
        assert r1162x.target_strand_abe == "non_coding"

    def test_03_abe8e_gaussian_window_and_pam_affinity(self) -> None:
        """Verify Gaussian activity profile and PAM variant affinity scaling."""
        # Positions 5 and 6 should be near peak weight (~0.9-1.0)
        w5 = compute_abe8e_activity_weight(5)
        w6 = compute_abe8e_activity_weight(6)
        assert w5 > 0.90
        assert w6 > 0.90
        assert w5 == w6  # Symmetric around 5.5

        # Position 1 and 20 should have negligible deamination
        w1 = compute_abe8e_activity_weight(1)
        w20 = compute_abe8e_activity_weight(20)
        assert w1 < 0.05
        assert w20 < 0.001

        # Positions outside 1-20 should be 0.0
        assert compute_abe8e_activity_weight(0) == 0.0
        assert compute_abe8e_activity_weight(21) == 0.0

        # PAM affinities
        assert get_pam_affinity("SpCas9", "AGG") == 1.0
        assert get_pam_affinity("SpCas9", "AGC") < 0.20
        assert get_pam_affinity("SpCas9-NG", "AGA") >= 0.85
        assert get_pam_affinity("SpRY", "AGA") >= 0.75

    def test_04_w1282x_exact_wildtype_restoration(self) -> None:
        """Verify that ABE8e directly restores wild-type Trp1282 (TGG) for W1282X."""
        w1282x = CLASS1_ALLELE_CATALOG["W1282X"]
        result = model_abe8e_allele(w1282x, "SpCas9")

        assert result.is_exact_wildtype is True
        assert result.repaired_amino_acid == "W"
        assert result.target_base_position == 7
        assert result.on_target_efficiency > 0.50
        assert result.bystander_count > 0
        assert result.bystander_purity_ratio > 0.40
        assert result.curative_verdict in ("HIGH_FIDELITY", "MODERATE_FIDELITY")

    def test_05_g542x_transversion_roadblock(self) -> None:
        """Verify that G542X is blocked from clean ABE reversion due to G>T transversion."""
        g542x = CLASS1_ALLELE_CATALOG["G542X"]
        result = model_abe8e_allele(g542x, "SpCas9")

        assert result.is_exact_wildtype is False
        assert result.curative_verdict == "TRANSVERSION_BLOCKED"
        assert result.overall_feasibility_score <= 30.0

    def test_06_r553x_and_r1162x_antisense_abe8e(self) -> None:
        """Verify antisense strand deamination restores wild-type Arg (CGA) for R553X and R1162X."""
        r553x = CLASS1_ALLELE_CATALOG["R553X"]
        res_r553x = model_abe8e_allele(r553x, "SpCas9")
        assert res_r553x.is_exact_wildtype is True
        assert res_r553x.repaired_amino_acid == "R"
        assert res_r553x.strand == "non_coding"

        r1162x = CLASS1_ALLELE_CATALOG["R1162X"]
        res_r1162x = model_abe8e_allele(r1162x, "SpRY")
        assert res_r1162x.is_exact_wildtype is True
        assert res_r1162x.repaired_amino_acid == "R"
        assert res_r1162x.strand == "non_coding"

    def test_07_prime_editing_pegRNA_design_and_100pct_fidelity(self) -> None:
        """Verify pegRNA thermodynamic design, zero bystander edits, and 100% sequence fidelity."""
        for name, allele in CLASS1_ALLELE_CATALOG.items():
            pe = design_prime_editing_allele(allele, epegRNA_motif="mpknot", use_pe3b=True)

            assert pe.bystander_amino_acid_alterations == 0, f"{name} must have 0 bystander mutations"
            assert pe.sequence_fidelity_pct == 100.0
            assert pe.curative_verdict == "OPTIMAL_FIDELITY"
            assert pe.overall_feasibility_score >= 60.0
            assert pe.indel_hazard_pct < 1.0  # PE3b protects from high indels

            # Thermodynamic PBS parameters
            assert 30.0 <= pe.pbs_tm_c <= 42.0, f"{name} PBS Tm {pe.pbs_tm_c} C out of expected range"
            assert pe.pbs_delta_g_kcal_mol < 0.0, f"{name} PBS Delta G must be negative (spontaneous binding)"

    def test_08_basal_stem_cell_delivery_and_10year_durability(self) -> None:
        """Verify in vivo mucus penetration, basal cell transduction, and 10-year persistence floor."""
        sort_lnp = model_basal_stem_cell_delivery("SORT-LNP")
        ev = model_basal_stem_cell_delivery("Targeted-EV")
        aav = model_basal_stem_cell_delivery("Dual-AAV")

        # SORT-LNP PEG brush ensures high mucus diffusion
        assert sort_lnp.peg_surface_density_mol_pct >= 5.0
        assert sort_lnp.mucus_diffusion_ratio_d_d0 > 0.50
        assert sort_lnp.net_basal_transduction_efficiency > 0.10

        # Basal stem cell persistence sustains >= 10% WT CFTR across 10 years
        assert sort_lnp.sustains_curative_threshold_10yr is True
        assert sort_lnp.functional_cftr_pct_wt_year10 >= CURATIVE_THRESHOLD_PCT_WT
        assert ev.sustains_curative_threshold_10yr is True
        assert ev.functional_cftr_pct_wt_year10 >= CURATIVE_THRESHOLD_PCT_WT

    def test_09_anti_hallucination_and_research_boundaries(self) -> None:
        """Rule AH-4: Verify absence of prohibited clinical overclaims and presence of disclaimers."""
        assert "ClinicalTrials.gov does not evaluate safety or efficacy" in NLM_REGISTRY_DISCLAIMER
        assert "Research, not medical advice" in RESEARCH_BOUNDARY

        # Test full module text
        mod_path = DISCOVERY_DIR / "class1_prime_base_editing.py"
        assert mod_path.exists()
        content = mod_path.read_text(encoding="utf-8").lower()

        prohibited = [
            "proven to cure",
            "guarantees clinical benefit",
            "recommended treatment regimen",
            "safe and effective for patients",
            "cures cystic fibrosis",
        ]
        for phrase in prohibited:
            assert phrase not in content, f"Prohibited phrase '{phrase}' found in module"

    def test_10_full_landscape_evaluation_and_decision_matrix(self) -> None:
        """Verify end-to-end evaluation pipeline and consensus decision matrix."""
        landscape = evaluate_class1_nonsense_landscape()

        assert "metadata" in landscape
        assert "alleles" in landscape
        assert "decision_matrix" in landscape
        assert len(landscape["decision_matrix"]) == 4

        # G542X must recommend Prime Editing
        g542_row = next(r for r in landscape["decision_matrix"] if r["allele"] == "G542X")
        assert "PRIME_EDITING_MANDATORY" in g542_row["consensus_curative_recommendation"]

        # W1282X should note ABE or PE
        w1282_row = next(r for r in landscape["decision_matrix"] if r["allele"] == "W1282X")
        assert ("ABE8e" in w1282_row["consensus_curative_recommendation"] or
                "PRIME_EDITING" in w1282_row["consensus_curative_recommendation"])
