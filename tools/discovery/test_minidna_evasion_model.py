#!/usr/bin/env python3
"""Unit tests for Nonviral Minimal DNA Vector (minidnacftr) CpG Depletion & TLR9 Evasion Model (T-0233).

Verifies CpG density calculation, TLR9 activation flux, cytokine reactogenicity,
nuclear persistence dynamics, and comparison against wild-type bacterial plasmid baselines.
"""
from __future__ import annotations

import sys
from pathlib import Path
import pytest

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(REPO_ROOT))

from minidna_evasion_model import (
    MiniDNAEvasionModel,
    compute_cpg_density,
    compute_tlr9_flux,
    compute_cytokine_reactogenicity,
    VECTOR_CATALOG,
)


class TestMiniDNAEvasionModel:
    @pytest.fixture
    def model(self) -> MiniDNAEvasionModel:
        return MiniDNAEvasionModel()

    def test_cpg_density_calculation(self) -> None:
        """Verify CpG density calculation rho_CpG = N_CpG / L_vector."""
        assert compute_cpg_density(0, 1000) == 0.0
        assert pytest.approx(compute_cpg_density(3, 1200), 1e-6) == 0.0025
        assert pytest.approx(compute_cpg_density(310, 5000), 1e-6) == 0.062

        with pytest.raises(ValueError, match="positive"):
            compute_cpg_density(10, 0)
        with pytest.raises(ValueError, match="non-negative"):
            compute_cpg_density(-1, 1000)

    def test_tlr9_flux_and_saturation(self) -> None:
        """TLR9 activation flux scales with CpG count and saturates at high endosomal dose."""
        flux_zero = compute_tlr9_flux(0, 50.0)
        flux_low = compute_tlr9_flux(3, 50.0)
        flux_high = compute_tlr9_flux(310, 50.0)

        assert flux_zero == 0.0
        assert flux_low > 0.0
        assert flux_high > flux_low * 50.0

        # Michaelis-Menten saturation check
        flux_50nm = compute_tlr9_flux(10, 50.0)
        flux_500nm = compute_tlr9_flux(10, 500.0)
        flux_5000nm = compute_tlr9_flux(10, 5000.0)
        assert flux_50nm < flux_500nm < flux_5000nm
        # Saturation asymptote should approach k_tlr9 * N_CpG = 0.42 * 10 = 4.2
        assert flux_5000nm < 4.2

    def test_cytokine_reactogenicity_index_bounds(self) -> None:
        """R_cytokine must be in [0, 1] and strictly monotone with flux."""
        r_zero = compute_cytokine_reactogenicity(0.0)
        r_low = compute_cytokine_reactogenicity(0.74)
        r_high = compute_cytokine_reactogenicity(76.5)

        assert r_zero == 0.0
        assert 0.0 < r_low < 0.20  # Clinically tolerable tier
        assert r_high > 0.99       # Severe cytokine reactivity
        assert r_low < r_high

    def test_plasmid_vs_minidna_reactogenicity_and_persistence(self, model: MiniDNAEvasionModel) -> None:
        """minidnacftr shows dramatic reduction in reactogenicity and enhanced day-90 persistence vs pDNA."""
        res_plasmid = model.evaluate_vector("bacterial_plasmid", endosomal_dose_nm=50.0)
        res_minidna = model.evaluate_vector("minidnacftr_4x", endosomal_dose_nm=50.0)

        assert res_plasmid["safety_tier"] == "SEVERE_TLR9_CYTOKINE_REACTIVITY"
        assert res_minidna["safety_tier"] == "CLINICALLY_TOLERABLE_LOW_REACTOGENICITY"

        # Reactogenicity reduction > 5-fold
        fold_reduction = res_plasmid["r_cytokine_index"] / res_minidna["r_cytokine_index"]
        assert fold_reduction > 5.0

        # Day 90 nuclear persistence: plasmid is nearly completely cleared by silencing/dilution,
        # while minidnacftr retains substantial nuclear copy
        assert res_plasmid["nuclear_vector_day90_nm"] < 0.01
        assert res_minidna["nuclear_vector_day90_nm"] > 3.0
        assert res_minidna["retention_day90_vs_peak_pct"] > 25.0

        # Curative benchmark: minidnacftr meets >= 10% WT CFTR current at day 30, plasmid fails
        assert res_minidna["meets_curative_10pct_wt_day30"] is True
        assert res_plasmid["meets_curative_10pct_wt_day30"] is False

    def test_zero_cpg_doggybone_dna(self, model: MiniDNAEvasionModel) -> None:
        """Synthetic dbDNA with zero CpGs exhibits zero TLR9 flux and highest 90-day retention."""
        res_dbdna = model.evaluate_vector("dbdna_zero_cpg", endosomal_dose_nm=50.0)
        assert res_dbdna["j_tlr9_flux"] == 0.0
        assert res_dbdna["r_cytokine_index"] == 0.0
        assert res_dbdna["safety_tier"] == "CLINICALLY_TOLERABLE_LOW_REACTOGENICITY"
        assert res_dbdna["retention_day90_vs_peak_pct"] > 40.0
        assert res_dbdna["restored_cftr_current_day30_pct_wt"] > 30.0

    def test_compare_all_vectors_catalog(self, model: MiniDNAEvasionModel) -> None:
        """compare_all_vectors evaluates all 5 cataloged vector formats."""
        comparison = model.compare_all_vectors(endosomal_dose_nm=50.0)
        assert len(comparison["vectors"]) == len(VECTOR_CATALOG)
        assert comparison["minidnacftr_vs_plasmid_reactogenicity_fold_reduction"] > 5.0
        assert comparison["minidnacftr_vs_plasmid_day90_nuclear_ratio"] > 1000.0

    def test_evidence_grounding_and_boundary_tags(self, model: MiniDNAEvasionModel) -> None:
        """Each vector evaluation retains NLM attribution, surveillance NCT, and research boundary."""
        res = model.evaluate_vector("minidnacftr_4x")
        assert "NLM" in res["nlm_attribution"]
        assert res["clinical_surveillance_nct"] == "NCT05668741"
        assert "research, not medical advice" in res["research_boundary"]
        assert res["citation_doi"] == "10.1093/nar/gkag082"
