#!/usr/bin/env python3
"""Unit tests for Suppressor tRNA PTC Readthrough & Transcriptome Orthogonality Model (T-0394).

Verifies ribosomal A-site competition, context hierarchy, transcriptome orthogonality index (OI >= 10^2),
NMD mRNA decay kinetics, cognate folding fidelity, and electrophysiological channel restoration.
"""
import math
import subprocess
import sys
from pathlib import Path
import pytest

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(REPO_ROOT))

from suppressor_trna_orthogonality import (
    CLINICAL_CURATIVE_THRESHOLD_PCT,
    IVACAFTOR_POTENTIATED_PO,
    ORTHOGONALITY_SAFETY_THRESHOLD,
    STOP_CODON_CONTEXTS,
    SUPPRESSOR_TRNA_PROFILES,
    WT_CAMP_OPEN_PROBABILITY,
    WT_CFTR_CONDUCTANCE_PS,
    SuppressorTRNAOrthogonalityModel,
    run_comprehensive_matrix,
)


class TestSuppressorTRNAOrthogonalityModel:
    @pytest.fixture
    def model(self):
        return SuppressorTRNAOrthogonalityModel()

    def test_context_hierarchy_ptc_readthrough(self, model):
        """UGA-C (permissive) > UGA-A (W1282X) > UGA-G (G542X) at identical tRNA concentrations."""
        rt_c, _ = model.compute_ptc_readthrough("UGA_C_PERMISSIVE", "tRNA_Trp_UCA", 30.0)
        rt_a, _ = model.compute_ptc_readthrough("W1282X_UGA_A", "tRNA_Trp_UCA", 30.0)
        rt_g, _ = model.compute_ptc_readthrough("G542X_UGA_G", "tRNA_Trp_UCA", 30.0)

        assert rt_c > rt_a > rt_g
        assert rt_c > 0.08  # >8% readthrough
        assert rt_g > 0.02  # >2% readthrough

    def test_orthogonality_index_safety_threshold(self, model):
        """Transcriptome Orthogonality Index (OI) must be >= 10^2 (100.0) safeguarding natural stop codons."""
        for conc in [10.0, 30.0, 60.0]:
            sim_w1282x = model.simulate("W1282X_UGA_A", "tRNA_Trp_UCA", conc, 0.0)
            assert sim_w1282x.orthogonality_index >= ORTHOGONALITY_SAFETY_THRESHOLD
            assert sim_w1282x.orthogonality_pass is True

            # Off-target readthrough across NTCs must remain low (< 0.5%)
            assert sim_w1282x.transcriptome_off_target_pct < 0.50

    def test_trna_concentration_monotonicity_and_saturation(self, model):
        """Readthrough increases monotonically with tRNA concentration and saturates within physical limit (<= 40%)."""
        rt_5, _ = model.compute_ptc_readthrough("W1282X_UGA_A", "tRNA_Trp_UCA", 5.0)
        rt_15, _ = model.compute_ptc_readthrough("W1282X_UGA_A", "tRNA_Trp_UCA", 15.0)
        rt_30, _ = model.compute_ptc_readthrough("W1282X_UGA_A", "tRNA_Trp_UCA", 30.0)
        rt_120, _ = model.compute_ptc_readthrough("W1282X_UGA_A", "tRNA_Trp_UCA", 120.0)

        assert rt_5 < rt_15 < rt_30 < rt_120
        assert rt_120 <= 0.40  # Physical ceiling

    def test_cognate_vs_non_cognate_folding_superiority(self, model):
        """Cognate tRNA-Trp yields higher mature Band C protein than non-cognate tRNA-Arg insertion."""
        sim_cog = model.simulate("W1282X_UGA_A", "tRNA_Trp_UCA", 30.0, 50.0)
        sim_non_cog = model.simulate("W1282X_UGA_A", "tRNA_Arg_UCA", 30.0, 50.0)

        # Cognate folding eta is 0.92 vs 0.32
        assert sim_cog.mature_band_c_pct > 2.0 * sim_non_cog.mature_band_c_pct

    def test_nmd_inhibition_synergy(self, model):
        """NMD inhibition rescues transcript stability (tau_mRNA) and synergistically elevates mature CFTR."""
        # Basal (uninhibited NMD)
        sim_basal = model.simulate("G542X_UGA_G", "tRNA_Gly_UCA", 30.0, 0.0)
        # 80% NMD inhibition
        sim_nmd_inh = model.simulate("G542X_UGA_G", "tRNA_Gly_UCA", 30.0, 80.0)

        assert sim_nmd_inh.tau_mrna_hours > 2.5 * sim_basal.tau_mrna_hours
        assert sim_nmd_inh.mrna_abundance_pct > 3.0 * sim_basal.mrna_abundance_pct
        assert sim_nmd_inh.mature_band_c_pct > 3.0 * sim_basal.mature_band_c_pct

    def test_ivacaftor_potentiation_coupling(self, model):
        """Ivacaftor potentiation scales open probability from 0.38 to 0.72 (1.89x)."""
        sim = model.simulate("W1282X_UGA_A", "tRNA_Trp_UCA", 30.0, 50.0)
        ratio = sim.delta_isc_potentiated_pct / sim.delta_isc_basal_pct
        assert pytest.approx(ratio, rel=1e-3) == (IVACAFTOR_POTENTIATED_PO / WT_CAMP_OPEN_PROBABILITY)

    def test_clinical_curative_threshold_attainment(self, model):
        """Combined suppressor tRNA + NMD inhibition + ivacaftor achieves curative threshold (>= 10% WT Delta Isc)."""
        sim_combined = model.simulate("W1282X_UGA_A", "tRNA_Trp_UCA", 45.0, 60.0)
        assert sim_combined.delta_isc_potentiated_pct >= CLINICAL_CURATIVE_THRESHOLD_PCT
        assert sim_combined.meets_clinical_curative_threshold is True

    def test_comprehensive_matrix_execution(self, model):
        """Comprehensive matrix runs cleanly and produces all expected simulation combinations."""
        res = run_comprehensive_matrix(model)
        assert res["status"] == "PASS"
        assert res["total_simulations"] == 108  # 3 contexts * 3 tRNAs * 6 concs * 2 NMD states
        assert len(res["results"]) == 108

    def test_cli_verification_flag(self):
        """CLI --verify flag passes and exits 0."""
        script_path = HERE / "suppressor_trna_orthogonality.py"
        proc = subprocess.run(
            [sys.executable, str(script_path), "--verify"],
            capture_output=True,
            text=True,
            check=True,
        )
        assert "PASS: All deterministic biophysical invariants verified successfully." in proc.stdout
