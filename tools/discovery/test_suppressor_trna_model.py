#!/usr/bin/env python3
"""Unit tests for Suppressor tRNA & SMG1 NMD Inhibition Biophysical Model (T-0216).

Verifies ribosomal A-site competition, tetranucleotide context hierarchy,
NMD mRNA decay suppression, supra-additive synergy, and electrophysiological restoration.
"""
import sys
from pathlib import Path
import pytest

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(REPO_ROOT))

from suppressor_trna_model import (
    CLINICAL_CURATIVE_THRESHOLD_PCT,
    IVACAFTOR_POTENTIATED_PO,
    STOP_CODON_CONTEXTS,
    SUPPRESSOR_TRNA_PROFILES,
    WT_CAMP_OPEN_PROBABILITY,
    WT_CFTR_CONDUCTANCE_PS,
    SuppressorTRNAModel,
    run_comprehensive_matrix,
)


class TestSuppressorTRNAModel:
    @pytest.fixture
    def model(self):
        return SuppressorTRNAModel()

    def test_context_hierarchy_readthrough(self, model):
        """UGA-C (permissive) > UGA-A (W1282X) > UGA-G (G542X) at identical tRNA concentrations."""
        rt_c, _ = model.compute_readthrough("UGA_C_PERMISSIVE", "tRNA_Trp_UCA", 30.0)
        rt_a, _ = model.compute_readthrough("W1282X_UGA_A", "tRNA_Trp_UCA", 30.0)
        rt_g, _ = model.compute_readthrough("G542X_UGA_G", "tRNA_Gly_UCA", 30.0)

        assert rt_c > rt_a > rt_g
        assert rt_c > 20.0
        assert rt_g > 10.0

    def test_smg1_mrna_tau_extension(self, model):
        """SMG1 inhibition extends G542X mRNA tau from ~2.1h to ~8.4h."""
        # Uninhibited G542X
        tau_0, mrna_0, inh_0 = model.compute_mrna_dynamics("G542X_UGA_G", 0.0)
        assert pytest.approx(tau_0, 0.1) == 2.1
        assert pytest.approx(mrna_0, 0.5) == 14.5
        assert inh_0 == 0.0

        # Saturating SMG1 inhibition (250 nM)
        tau_sat, mrna_sat, inh_sat = model.compute_mrna_dynamics("G542X_UGA_G", 250.0)
        assert tau_sat > 7.5
        assert mrna_sat > 55.0
        assert inh_sat > 75.0

    def test_trna_concentration_saturation(self, model):
        """Readthrough efficiency increases monotonically and saturates with tRNA concentration."""
        rt_low, _ = model.compute_readthrough("W1282X_UGA_A", "tRNA_Trp_UCA", 5.0)
        rt_med, _ = model.compute_readthrough("W1282X_UGA_A", "tRNA_Trp_UCA", 25.0)
        rt_high, _ = model.compute_readthrough("W1282X_UGA_A", "tRNA_Trp_UCA", 100.0)

        assert rt_low < rt_med < rt_high
        assert rt_high <= 40.0  # Physical ceiling

    def test_cognate_vs_non_cognate_insertion(self, model):
        """Cognate tRNA-Trp yields WT 8.5 pS conductance and high folding; non-cognate tRNA-Arg is lower."""
        sim_cognate = model.simulate("W1282X_UGA_A", "tRNA_Trp_UCA", 30.0, 100.0, True)
        sim_non_cognate = model.simulate("W1282X_UGA_A", "tRNA_Arg_UCA", 30.0, 100.0, True)

        assert sim_cognate.single_channel_conductance_pS == WT_CFTR_CONDUCTANCE_PS
        assert sim_non_cognate.single_channel_conductance_pS < WT_CFTR_CONDUCTANCE_PS
        assert sim_cognate.band_c_pct_wt > sim_non_cognate.band_c_pct_wt
        assert sim_cognate.delta_isc_pct_wt > sim_non_cognate.delta_isc_pct_wt

    def test_supra_additive_synergy(self, model):
        """Bimodal combination yields supra-additive synergy quotient (> 1.5x) over monotherapies."""
        synergy_g542 = model.evaluate_synergy("G542X_UGA_G", "tRNA_Gly_UCA", 30.0, 100.0, True)
        assert synergy_g542["is_supra_additive"] is True
        assert synergy_g542["synergy_quotient"] >= 2.0
        assert synergy_g542["curative_threshold_exceeded"] is True
        assert synergy_g542["bimodal_synergy_delta_isc"] >= CLINICAL_CURATIVE_THRESHOLD_PCT

    def test_curative_threshold_w1282x(self, model):
        """W1282X with bimodal therapy exceeds 10% curative threshold."""
        synergy_w1282 = model.evaluate_synergy("W1282X_UGA_A", "tRNA_Trp_UCA", 30.0, 100.0, True)
        assert synergy_w1282["curative_threshold_exceeded"] is True
        assert synergy_w1282["bimodal_synergy_delta_isc"] >= 20.0

    def test_ivacaftor_potentiator_amplification(self, model):
        """Ivacaftor potentiator amplifies Delta Isc by ~1.89x through Po elevation."""
        sim_unpot = model.simulate("W1282X_UGA_A", "tRNA_Trp_UCA", 30.0, 100.0, False)
        sim_pot = model.simulate("W1282X_UGA_A", "tRNA_Trp_UCA", 30.0, 100.0, True)

        ratio = sim_pot.delta_isc_pct_wt / sim_unpot.delta_isc_pct_wt
        expected_ratio = IVACAFTOR_POTENTIATED_PO / WT_CAMP_OPEN_PROBABILITY
        assert pytest.approx(ratio, 0.05) == expected_ratio

    def test_invalid_context_and_profile_raises(self, model):
        """Invalid context or tRNA profile raises ValueError."""
        with pytest.raises(ValueError):
            model.compute_mrna_dynamics("INVALID_CONTEXT", 50.0)
        with pytest.raises(ValueError):
            model.compute_readthrough("W1282X_UGA_A", "INVALID_TRNA", 30.0)

    def test_comprehensive_matrix_runs(self):
        """Comprehensive matrix runs cleanly across all alleles."""
        matrix = run_comprehensive_matrix()
        assert len(matrix) >= 4
        for k, item in matrix.items():
            assert "synergy_quotient" in item
            assert item["bimodal_synergy_delta_isc"] > 0.0
