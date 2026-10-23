#!/usr/bin/env python3
"""Deterministic verification suite for prime_editing_cftr.py (T-0184).

Verifies:
1. pegRNA Primer Binding Site (PBS) melting temperature bounds (30-37C).
2. Reverse Transcriptase Template (RTT) length constraints and processivity limits (10-35 nt).
3. PEmax, epegRNA, and delivery vector packaging limits (Single AAV vs Dual AAV vs LNP).
4. CFTR F508del and W1282X therapeutic prime editing feasibility.
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import prime_editing_cftr as pe


class TestPBSMeltingTemperatureBounds:
    """Verifies PBS thermodynamic calculation and 30-37C optimal window bounds."""

    def test_pbs_optimal_window(self):
        # 13 nt PBS with optimal Tm
        seq = "TATCATTGGTGTT"
        res = pe.score_pbs(seq)
        assert res.within_optimal_bounds is True
        assert 30.0 <= res.melting_temperature_c <= 37.0
        assert res.score >= 0.85
        assert "Optimal PBS duplex stability" in res.details

    def test_pbs_suboptimal_low_tm(self):
        # Very AT-rich short oligo resulting in Tm < 30 C
        seq = "ATATATATAA"
        res = pe.score_pbs(seq)
        assert res.within_optimal_bounds is False
        assert res.melting_temperature_c < 30.0
        assert res.score < 0.75
        assert "Suboptimal low PBS Tm" in res.details

    def test_pbs_suboptimal_high_tm(self):
        # GC-rich oligo resulting in Tm > 37 C
        seq = "CCGCGCCGCGGC"
        res = pe.score_pbs(seq)
        assert res.within_optimal_bounds is False
        assert res.melting_temperature_c > 37.0
        assert res.score < 0.70
        assert "High PBS Tm" in res.details

    def test_pbs_length_bounds(self):
        # Out of bounds length (< 10 nt)
        seq_short = "ATGC"
        res = pe.score_pbs(seq_short)
        assert res.within_optimal_bounds is False

    def test_thermodynamic_properties(self):
        # Monotonicity with salt: higher Na+ stabilizes duplex -> higher Tm
        seq = "ACACCAAAGATG"
        tm_low_salt = pe.calculate_melting_temperature_nn(seq, na_m=0.020)
        tm_high_salt = pe.calculate_melting_temperature_nn(seq, na_m=0.100)
        assert tm_high_salt > tm_low_salt


class TestRTTLengthConstraints:
    """Verifies RTT length limits, processivity constraints, and distance penalties."""

    def test_rtt_optimal_range(self):
        seq = "ACACCAAAGATGATA"  # 15 nt
        res = pe.score_rtt(seq, nick_to_edit_distance_nt=4)
        assert res.within_constraints is True
        assert pe.RTT_LENGTH_OPTIMAL_MIN_NT <= res.length_nt <= pe.RTT_LENGTH_OPTIMAL_MAX_NT
        assert res.score >= 0.85
        assert "Optimal RTT length" in res.details

    def test_rtt_minimum_homology_constraint(self):
        seq = "ACG"  # 3 nt (< 10 nt)
        res = pe.score_rtt(seq, nick_to_edit_distance_nt=2)
        assert res.within_constraints is False
        assert res.score <= 0.20
        assert "Inviable RTT length" in res.details

    def test_rtt_maximum_processivity_limit(self):
        # RTT exceeding 35 nt processivity threshold
        seq = "A" * 42  # 42 nt
        res = pe.score_rtt(seq, nick_to_edit_distance_nt=10)
        assert res.within_constraints is False
        assert res.score <= 0.30
        assert "Exceeded processivity limit" in res.details

    def test_rtt_acceptable_extended_range(self):
        # 30 nt is acceptable but penalized for potential RT stalling
        seq = "ACGT" * 7 + "AA"  # 30 nt
        res = pe.score_rtt(seq, nick_to_edit_distance_nt=8)
        assert res.within_constraints is True
        assert res.score < 0.85


class TestPackagingAndEpegRNALimits:
    """Verifies delivery packaging capacity (AAV vs Dual AAV vs LNP) and epegRNA motifs."""

    def test_single_aav_packaging_exceeded(self):
        pack = pe.evaluate_packaging_limits(system_cds_bp=pe.PEMAX_CDS_BP)
        assert pack.single_aav_feasible is False
        assert pack.pemax_size_bp > pe.AAV_PACKAGING_LIMIT_BP
        assert "Single AAV packaging exceeded" in pack.recommendation

    def test_dual_aav_split_intein_feasibility(self):
        pack = pe.evaluate_packaging_limits(system_cds_bp=pe.PEMAX_CDS_BP)
        assert pack.dual_aav_split_intein_feasible is True

    def test_lnp_mrna_feasibility(self):
        pack = pe.evaluate_packaging_limits(system_cds_bp=pe.PEMAX_CDS_BP)
        assert pack.lnp_mrna_feasible is True
        total_cargo = pe.PEMAX_CDS_BP + pe.EPEGRNA_AVERAGE_NT + pe.NICKING_SGRNA_NT
        assert total_cargo <= pe.LNP_MRNA_PACKAGING_LIMIT_NT

    def test_epegrna_structural_motif_enhancement(self):
        # Engineered motifs vs linear
        mpknot = pe.score_epegRNA_motif("mpknot")
        tevpre = pe.score_epegRNA_motif("tevpre")
        linear = pe.score_epegRNA_motif("unmodified")

        assert mpknot.is_engineered is True
        assert tevpre.is_engineered is True
        assert linear.is_engineered is False
        assert mpknot.stability_retention_score > linear.stability_retention_score
        assert mpknot.stability_retention_score == 1.0


class TestCFTRPrimeEditingDesigns:
    """Validates complete biophysical designs for F508del and W1282X."""

    def test_f508del_cftr_design(self):
        design = pe.build_cftr_f508del_prime_editing_design()
        assert design.target_id == "T-0184-F508DEL"
        assert "F508del" in design.cftr_mutation
        assert "6MSM" in design.cryoem_pdb_context
        # PBS within 30-37C
        assert design.pbs.within_optimal_bounds is True
        assert 30.0 <= design.pbs.melting_temperature_c <= 37.0
        # RTT length constraints
        assert design.rtt.within_constraints is True
        assert design.rtt.length_nt <= 35
        # Packaging
        assert design.packaging.single_aav_feasible is False
        assert design.packaging.lnp_mrna_feasible is True
        # Score & Verdict
        assert design.composite_feasibility_score >= 85.0
        assert design.verdict == "OPTIMAL_FEASIBILITY"
        assert "DOI 10.1016/j.cell.2021.09.018" in design.evidence_citations

    def test_w1282x_cftr_design(self):
        design = pe.build_cftr_w1282x_prime_editing_design()
        assert design.target_id == "T-0184-W1282X"
        assert "W1282X" in design.cftr_mutation
        assert "6MSM" in design.cryoem_pdb_context
        # PBS within 30-37C
        assert design.pbs.within_optimal_bounds is True
        assert 30.0 <= design.pbs.melting_temperature_c <= 37.0
        # RTT length constraints
        assert design.rtt.within_constraints is True
        assert design.rtt.length_nt <= 35
        # Packaging
        assert design.packaging.single_aav_feasible is False
        assert design.packaging.lnp_mrna_feasible is True
        # Score & Verdict
        assert design.composite_feasibility_score >= 85.0
        assert design.verdict == "OPTIMAL_FEASIBILITY"
        assert "DOI 10.1016/j.cell.2021.09.018" in design.evidence_citations


class TestNickingSgRNAConstraints:
    """Verifies PE3 nicking sgRNA spacing and indel hazard prevention."""

    def test_indel_hazard_under_threshold(self):
        # Spacing of 20 bp without PE3b creates indel hazard
        res = pe.score_nicking_guide("ACGTACGTACGTACGTACGT", relative_distance_bp=20, is_pe3b=False)
        assert res.indel_hazard_flag is True
        assert res.score <= 0.30
        assert "High indel hazard" in res.details

    def test_pe3b_exemption(self):
        # Spacing of 20 bp with PE3b is exempt from indel hazard
        res = pe.score_nicking_guide("ACGTACGTACGTACGTACGT", relative_distance_bp=20, is_pe3b=True)
        assert res.indel_hazard_flag is False
        assert res.is_pe3b is True
        assert res.score == 1.0


class TestCLIExecution:
    """Verifies command-line execution and JSON export."""

    def test_cli_json_export(self, tmp_path):
        export_file = tmp_path / "pe_export.json"
        script_path = HERE / "prime_editing_cftr.py"
        cmd = [sys.executable, str(script_path), "--target", "all", "--export", str(export_file), "--json"]
        proc = subprocess.run(cmd, capture_output=True, text=True, check=True)
        assert proc.returncode == 0
        assert export_file.exists()

        with open(export_file, encoding="utf-8") as f:
            data = json.load(f)
        assert "F508del" in data
        assert "W1282X" in data
        assert data["F508del"]["composite_feasibility_score"] >= 85.0
        assert data["W1282X"]["composite_feasibility_score"] >= 85.0
