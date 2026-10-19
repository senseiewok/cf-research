#!/usr/bin/env python3
"""Unit tests for CFTR Cryo-EM PTC Readthrough Allosteric Rescue Model (T-0217).

Verifies atomistic C-alpha coordinate extraction, functional landmark distances,
internal cavity volume calculations, thermal destabilization (Delta Tm),
allosteric modulator rescue, and electrophysiological conductance.
"""
from __future__ import annotations

import sys
from pathlib import Path
import pytest

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(REPO_ROOT))

from ptc_cryoem_allostery import (
    CLINICAL_CURATIVE_THRESHOLD_PCT,
    IVACAFTOR_POTENTIATED_PO,
    NLM_ATTRIBUTION,
    RESEARCH_BOUNDARY,
    WT_CAMP_OPEN_PROBABILITY,
    WT_CFTR_CONDUCTANCE_PS,
    PTCCryoEMAllosteryModel,
)


class TestPTCCryoEMAllostery:
    @pytest.fixture
    def model(self):
        return PTCCryoEMAllosteryModel()

    def test_pdb_calpha_extraction(self, model):
        """Extracts valid C-alpha coordinates for Trp1282 and Gly542 in 6MSM and 8EIQ."""
        c1282_6msm = model.get_calpha_coordinates("6MSM", 1282)
        c542_6msm = model.get_calpha_coordinates("6MSM", 542)
        assert c1282_6msm is not None
        assert c542_6msm is not None
        assert len(c1282_6msm) == 3
        assert len(c542_6msm) == 3
        assert pytest.approx(c1282_6msm[0], 0.5) == -11.7
        assert pytest.approx(c542_6msm[0], 0.5) == -10.6

        c1282_8eiq = model.get_calpha_coordinates("8EIQ", 1282)
        c542_8eiq = model.get_calpha_coordinates("8EIQ", 542)
        assert c1282_8eiq is not None
        assert c542_8eiq is not None

    def test_landmark_distances(self, model):
        """Computes Euclidean distances to functional landmarks in Cryo-EM coordinates."""
        dist_1282 = model.compute_landmark_distances("6MSM", 1282)
        assert "Walker_B_1371" in dist_1282
        assert "ABC_sig_1347" in dist_1282
        assert "ATP_site_2" in dist_1282
        assert "TMD1_pocket" in dist_1282
        assert "TMD2_pocket" in dist_1282

        # Walker B catalytic Glu1371 is proximal (~21.6 A) in NBD2
        assert 18.0 <= dist_1282["Walker_B_1371"] <= 25.0
        # TMD1 Tezacaftor pocket is > 50 A away (long-range allosteric coupling)
        assert dist_1282["TMD1_pocket"] > 45.0

        dist_542 = model.compute_landmark_distances("6MSM", 542)
        # Gly542 is in NBD1, close to TMD2 pocket (~14.9 A) and ICL4 coupling (~18.0 A)
        assert dist_542["TMD2_pocket"] < 20.0
        assert dist_542["ICL4_coupling"] < 25.0

    def test_w1282x_vacuum_cavity_quantification(self, model):
        """Near-cognate Cys creates ~119-120 A^3 vacuum cavity and thermal destabilization."""
        # Cognate Trp
        cog = model.compute_substitution_energetics("W1282X", "W")
        assert cog.readthrough_type == "cognate"
        assert cog.cavity_nature == "isosteric_native"
        assert cog.volume_delta_A3 == 0.0
        assert cog.delta_tm_celsius == 0.0
        assert cog.delta_delta_g_destabilization_kcal_mol == 0.0
        assert cog.single_channel_conductance_pS == WT_CFTR_CONDUCTANCE_PS
        assert cog.erad_targeting_relative_rate == 1.0

        # Near-cognate Cys
        nc_cys = model.compute_substitution_energetics("W1282X", "C")
        assert nc_cys.readthrough_type == "near_cognate"
        assert nc_cys.cavity_nature == "vacuum_cavity"
        # Volume delta: 108.5 - 227.8 = -119.3 A^3 (~120 A^3 vacuum cavity)
        assert -120.0 <= nc_cys.volume_delta_A3 <= -118.0
        assert nc_cys.delta_tm_celsius < -4.0
        assert nc_cys.delta_delta_g_destabilization_kcal_mol > 2.5
        assert nc_cys.single_channel_conductance_pS < WT_CFTR_CONDUCTANCE_PS
        assert nc_cys.folding_efficiency_band_c_pct < cog.folding_efficiency_band_c_pct
        assert nc_cys.erad_targeting_relative_rate > 3.0

    def test_w1282x_charge_burial_destabilization(self, model):
        """Near-cognate Arg substitution causes severe charge burial destabilization."""
        nc_arg = model.compute_substitution_energetics("W1282X", "R")
        assert nc_arg.readthrough_type == "near_cognate"
        assert nc_arg.cavity_nature == "charge_burial"
        assert nc_arg.delta_tm_celsius <= -6.0
        assert nc_arg.folding_efficiency_band_c_pct < 15.0
        assert nc_arg.erad_targeting_relative_rate > 5.0

    def test_g542x_steric_clash_quantification(self, model):
        """Near-cognate substitutions at Gly542 produce steric clashes in tight NBD1 core."""
        cog_gly = model.compute_substitution_energetics("G542X", "G")
        assert cog_gly.readthrough_type == "cognate"
        assert cog_gly.volume_delta_A3 == 0.0
        assert cog_gly.delta_tm_celsius == 0.0

        # Cys has +48.4 A^3 steric excess over Gly
        nc_cys = model.compute_substitution_energetics("G542X", "C")
        assert nc_cys.readthrough_type == "near_cognate"
        assert nc_cys.cavity_nature == "steric_clash"
        assert nc_cys.volume_delta_A3 > 40.0
        assert nc_cys.delta_tm_celsius < -3.5

        # Trp has massive +167.7 A^3 steric clash
        nc_trp = model.compute_substitution_energetics("G542X", "W")
        assert nc_trp.cavity_nature == "extreme_steric_clash"
        assert nc_trp.volume_delta_A3 > 150.0
        assert nc_trp.delta_tm_celsius < -7.5

    def test_allosteric_rescue_trikafta_cooperativity(self, model):
        """Trikafta confers cooperative thermal stabilization (+5.8 C) and elevates Delta Isc."""
        # Cognate Trp1282 + Trikafta + Ivacaftor
        rescue_cog = model.simulate_allosteric_rescue("W1282X", "W", corrector="trikafta", potentiated=True)
        assert rescue_cog.net_delta_tm_celsius > 5.0
        assert rescue_cog.net_band_c_pct_wt >= 85.0
        assert rescue_cog.exceeds_curative_threshold is True
        assert rescue_cog.delta_isc_pct_wt > 50.0  # Robust cure-level conductance

        # Near-cognate Cys1282 unpotentiated monotherapy
        rescue_nc_unpot = model.simulate_allosteric_rescue("W1282X", "C", corrector="none", potentiated=False)
        assert rescue_nc_unpot.net_delta_tm_celsius < -4.0
        assert rescue_nc_unpot.exceeds_curative_threshold is False
        assert rescue_nc_unpot.delta_isc_pct_wt < 5.0  # Sub-therapeutic

    def test_potentiation_gating_amplification(self, model):
        """Ivacaftor potentiation elevates effective Po and amplifies Delta Isc."""
        unpot = model.simulate_allosteric_rescue("W1282X", "W", corrector="none", potentiated=False)
        pot = model.simulate_allosteric_rescue("W1282X", "W", corrector="none", potentiated=True)

        assert unpot.effective_po == WT_CAMP_OPEN_PROBABILITY
        assert pot.effective_po == IVACAFTOR_POTENTIATED_PO
        ratio = pot.delta_isc_pct_wt / unpot.delta_isc_pct_wt
        assert pytest.approx(ratio, 0.05) == IVACAFTOR_POTENTIATED_PO / WT_CAMP_OPEN_PROBABILITY

    def test_full_report_structure(self, model):
        """generate_full_allostery_report emits all required epistemic sections and citations."""
        rep = model.generate_full_allostery_report()
        assert rep["curative_track"] == "Track A (Biophysical Science & Cure Discovery)"
        assert rep["answer_tier"] == "T2"
        assert "6MSM" in rep["reference_structures"]
        assert "8EIQ" in rep["reference_structures"]
        assert len(rep["primary_citations"]) >= 3
        assert "W1282X" in rep["substitutions"]
        assert "G542X" in rep["substitutions"]
        assert len(rep["allosteric_rescue_scenarios"]) > 10
        assert NLM_ATTRIBUTION in rep["nlm_attribution"]
        assert RESEARCH_BOUNDARY in rep["research_boundary"]

    def test_invalid_parameters_raise(self, model):
        """Invalid variant, amino acid, or corrector raises appropriate ValueError."""
        with pytest.raises(ValueError):
            model.compute_substitution_energetics("UNKNOWN_VARIANT", "W")

        with pytest.raises(ValueError):
            model.compute_substitution_energetics("W1282X", "Z")

        with pytest.raises(ValueError):
            model.simulate_allosteric_rescue("W1282X", "W", corrector="unknown_drug")

        with pytest.raises(KeyError):
            model.get_calpha_coordinates("NONEXISTENT_PDB", 1282)
