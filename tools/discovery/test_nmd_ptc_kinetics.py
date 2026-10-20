#!/usr/bin/env python3
"""Unit tests for Nonsense-Mediated mRNA Decay UPF1-SMG1 Branch Kinetics & PTC Escape Model (T-0241).

Verifies:
  - Ribosome dwell times and +4 flanking nucleotide context sensitivity (W1282X vs G542X)
  - SMG1 kinase inhibition kinetics (Ki = 25.0 nM)
  - UPF1 phosphorylation cycle and EJC boundary condition
  - Dual-branch NMD degradation (SMG6 endonucleolytic vs SMG5/7 exonucleolytic)
  - CFTR transcript half-life and steady-state abundance
  - Supra-additive synergy of SMG1i + Readthrough + Corrector exceeding the 10.0% WT curative Delta I_sc threshold
  - Benchmarks against ELX-02 Phase 2 trial surveillance (NCT04135495)
"""

from __future__ import annotations

import math
import sys
from pathlib import Path
import pytest

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(REPO_ROOT))

from nmd_ptc_kinetics import (
    NMDPtcKineticsModel,
    CFTR_PTC_CATALOG,
    PTCVariantSpec,
    compute_dwell_time,
    compute_readthrough_probability,
    compute_smg1_activity,
    compute_upf1_phosphorylation_level,
    compute_nmd_decay_flux,
    compute_cftr_mrna_kinetics,
    compute_functional_restoration_isc,
    format_markdown_report,
    DEFAULT_CURATIVE_ISC_THRESHOLD_PCT,
    DEFAULT_SMG1_KI_NM,
    DEFAULT_WT_HALF_LIFE_HR,
    DEFAULT_PTC_UNINHIBITED_HALF_LIFE_HR,
    NLM_ATTRIBUTION,
    RESEARCH_BOUNDARY,
)


class TestNMDPtcKineticsModel:
    @pytest.fixture
    def model(self) -> NMDPtcKineticsModel:
        return NMDPtcKineticsModel()

    def test_dwell_time_and_flanking_plus4_context(self) -> None:
        """+4 G (W1282X) exhibits lower termination rate and longer dwell time than +4 A (G542X)."""
        w1282x = CFTR_PTC_CATALOG["W1282X"]
        g542x = CFTR_PTC_CATALOG["G542X"]
        y122x = CFTR_PTC_CATALOG["Y122X"]

        tau_w1282x = compute_dwell_time(k_term=w1282x.k_term_per_s, k_read=0.0)
        tau_g542x = compute_dwell_time(k_term=g542x.k_term_per_s, k_read=0.0)
        tau_y122x = compute_dwell_time(k_term=y122x.k_term_per_s, k_read=0.0)

        # W1282X (UGA-G) dwell time ~ 0.95 s > G542X (UGA-A) dwell time ~ 0.38 s > Y122X (UAA-U) ~ 0.28 s
        assert tau_w1282x > tau_g542x
        assert tau_g542x > tau_y122x
        assert 0.80 <= tau_w1282x <= 1.20
        assert 0.30 <= tau_g542x <= 0.50

        # Adding readthrough tRNA activity decreases ribosome stall dwell time
        tau_rescued = compute_dwell_time(k_term=w1282x.k_term_per_s, k_read=0.30)
        assert tau_rescued < tau_w1282x

        # Invalid inputs raise ValueError
        with pytest.raises(ValueError, match="positive"):
            compute_dwell_time(k_term=-1.0, k_read=0.0)
        with pytest.raises(ValueError, match="positive"):
            compute_dwell_time(k_term=1.0, k_read=0.0, erf_conc=-0.5)

    def test_readthrough_probability_bounds_and_monotony(self) -> None:
        """Readthrough probability is bounded in [0, 1] and increases with tRNA recruitment."""
        # Zero readthrough tRNA produces 0.0 probability
        assert compute_readthrough_probability(k_term=1.5, k_read=0.0) == 0.0

        # Equal competition gives 0.50
        assert pytest.approx(compute_readthrough_probability(k_term=2.0, k_read=2.0), 1e-5) == 0.50

        # High readthrough gives asymptotically 1.0
        assert compute_readthrough_probability(k_term=1.0, k_read=999.0) > 0.99

        # Increases monotonically
        p1 = compute_readthrough_probability(k_term=2.0, k_read=0.1)
        p2 = compute_readthrough_probability(k_term=2.0, k_read=0.5)
        p3 = compute_readthrough_probability(k_term=2.0, k_read=1.0)
        assert p1 < p2 < p3

    def test_smg1_kinase_inhibition_kinetics(self) -> None:
        """SMG1 kinase activity follows competitive inhibition with Ki = 25.0 nM."""
        # No inhibitor -> 100% activity
        assert compute_smg1_activity(0.0, ki_nm=DEFAULT_SMG1_KI_NM) == 1.0

        # Exactly at Ki -> 50% activity
        assert pytest.approx(compute_smg1_activity(25.0, ki_nm=25.0), 1e-5) == 0.50

        # At 100 nM (4x Ki) -> 20% activity
        assert pytest.approx(compute_smg1_activity(100.0, ki_nm=25.0), 1e-5) == 0.20

        # At 250 nM (10x Ki) -> ~9.09% activity
        act_250 = compute_smg1_activity(250.0, ki_nm=25.0)
        assert pytest.approx(act_250, 1e-3) == 0.091

        # Monotonically decreasing
        assert compute_smg1_activity(10.0) > compute_smg1_activity(50.0) > compute_smg1_activity(200.0)

        with pytest.raises(ValueError, match="negative"):
            compute_smg1_activity(-5.0)

    def test_upf1_phosphorylation_and_ejc_boundary(self) -> None:
        """UPF1 phosphorylation requires premature termination >50 nt upstream of an EJC."""
        # Non-NMD substrate produces 0.0 UPF1 phosphorylation
        upf1_non_nmd = compute_upf1_phosphorylation_level(
            smg1_activity=1.0,
            dwell_time_s=0.50,
            readthrough_prob=0.01,
            is_nmd_substrate=False,
        )
        assert upf1_non_nmd == 0.0

        # Active NMD substrate produces positive phosphorylation
        upf1_active = compute_upf1_phosphorylation_level(
            smg1_activity=1.0,
            dwell_time_s=0.38,
            readthrough_prob=0.002,
            is_nmd_substrate=True,
        )
        assert 0.90 <= upf1_active <= 1.10

        # SMG1 inhibition directly reduces UPF1 phosphorylation
        upf1_inhibited = compute_upf1_phosphorylation_level(
            smg1_activity=0.20,
            dwell_time_s=0.38,
            readthrough_prob=0.002,
            is_nmd_substrate=True,
        )
        assert pytest.approx(upf1_inhibited / upf1_active, 1e-3) == 0.20

    def test_dual_branch_nmd_degradation_flux(self) -> None:
        """Dual-branch NMD degradation is partitioned into SMG6 (75%) and SMG5/7 (25%)."""
        decay_flux = compute_nmd_decay_flux(upf1_p_level=1.0)

        # Sum of branches equals total flux
        assert pytest.approx(decay_flux["j_smg6_per_hr"] + decay_flux["j_smg5_7_per_hr"], 1e-5) == decay_flux["j_total_nmd_per_hr"]

        # Branching ratios: SMG6 ~ 75%, SMG5/7 ~ 25%
        assert pytest.approx(decay_flux["smg6_branch_fraction"], 1e-3) == 0.75
        assert pytest.approx(decay_flux["smg5_7_branch_fraction"], 1e-3) == 0.25

        # Uninhibited rate matches empirical CFTR NMD decay (~0.288 hr^-1)
        assert 0.25 <= decay_flux["j_total_nmd_per_hr"] <= 0.35

        # Zero UPF1-P yields zero NMD degradation
        zero_flux = compute_nmd_decay_flux(upf1_p_level=0.0)
        assert zero_flux["j_total_nmd_per_hr"] == 0.0

    def test_mrna_kinetics_and_half_life(self) -> None:
        """Uninhibited PTC mRNA is severely depleted (<15% WT) and restored by NMD inhibition."""
        decay_uninhibited = compute_nmd_decay_flux(upf1_p_level=1.0)
        mrna_uninhibited = compute_cftr_mrna_kinetics(j_nmd_per_hr=decay_uninhibited["j_total_nmd_per_hr"])

        # Half-life of uninhibited G542X mRNA matches nominal range (1.5 - 3.0 hr)
        assert 1.8 <= mrna_uninhibited["t_half_hr"] <= 2.4
        # Depleted steady-state mRNA matches nominal range (< 15% WT)
        assert 10.0 <= mrna_uninhibited["mrna_pct_of_wt"] <= 15.0

        # With 250 nM SMG1i, NMD decay drops ~11-fold, restoring t1/2 > 8 hr and mRNA > 50% WT
        decay_inhibited = compute_nmd_decay_flux(upf1_p_level=0.091)
        mrna_inhibited = compute_cftr_mrna_kinetics(j_nmd_per_hr=decay_inhibited["j_total_nmd_per_hr"])

        assert mrna_inhibited["t_half_hr"] >= 8.0
        assert mrna_inhibited["mrna_pct_of_wt"] >= 50.0

        # Completely uninhibited wild-type (zero NMD flux) matches physiological stability
        mrna_wt = compute_cftr_mrna_kinetics(j_nmd_per_hr=0.0)
        assert pytest.approx(mrna_wt["t_half_hr"], 1e-3) == DEFAULT_WT_HALF_LIFE_HR
        assert pytest.approx(mrna_wt["mrna_pct_of_wt"], 1e-3) == 100.0

    def test_cftr_functional_restoration_and_curative_threshold(self, model: NMDPtcKineticsModel) -> None:
        """Tri-modal therapy synergistically exceeds the 10.0% WT curative Delta I_sc threshold."""
        # 1. Baseline untreated W1282X: severely sub-curative (< 0.1% WT)
        base = model.simulate_condition("W1282X", smg1i_conc_nm=0.0, readthrough_induced=False)
        assert base.delta_isc_pct_wt < 0.10
        assert not base.curative_threshold_met

        # 2. Readthrough monotherapy (ELX-02 benchmark on W1282X without SMG1i):
        # Despite 18.5% readthrough efficiency, depleted mRNA substrate caps Delta I_sc at ~2.0% WT (< 10%)
        rt_only = model.simulate_condition("W1282X", smg1i_conc_nm=0.0, readthrough_induced=True)
        assert 1.2 <= rt_only.delta_isc_pct_wt <= 2.8
        assert not rt_only.curative_threshold_met

        # 3. SMG1i monotherapy (250 nM SMG1i on W1282X without readthrough):
        # Despite mRNA stabilized to ~56%, basal readthrough (0.35%) caps Delta I_sc at ~0.17% WT (< 10%)
        smg1i_only = model.simulate_condition("W1282X", smg1i_conc_nm=250.0, readthrough_induced=False)
        assert smg1i_only.delta_isc_pct_wt < 0.50
        assert not smg1i_only.curative_threshold_met

        # 4. Dual Combination Therapy (250 nM SMG1i + Readthrough on W1282X):
        # Substantial rescue approaching the threshold (> 8.0% WT)
        dual = model.simulate_condition("W1282X", smg1i_conc_nm=250.0, readthrough_induced=True)
        assert dual.delta_isc_pct_wt >= 8.0

        # 5. Tri-Modal Combination (250 nM SMG1i + Readthrough + Corrector on W1282X):
        # Boosts trafficking to 95%, delivering Delta I_sc >= 10.0% WT curative threshold!
        tri = model.simulate_condition(
            "W1282X",
            smg1i_conc_nm=250.0,
            readthrough_induced=True,
            corrector_co_administered=True,
        )
        assert tri.delta_isc_pct_wt >= DEFAULT_CURATIVE_ISC_THRESHOLD_PCT
        assert tri.curative_threshold_met
        assert tri.delta_isc_pct_wt > dual.delta_isc_pct_wt

    def test_g542x_strict_context_response(self, model: NMDPtcKineticsModel) -> None:
        """G542X (+4 A) has lower readthrough permissiveness than W1282X (+4 G)."""
        g542x_dual = model.simulate_condition("G542X", smg1i_conc_nm=250.0, readthrough_induced=True)
        w1282x_dual = model.simulate_condition("W1282X", smg1i_conc_nm=250.0, readthrough_induced=True)

        # W1282X (+4 G) produces greater Delta I_sc than G542X (+4 A) due to higher readthrough permissiveness
        assert w1282x_dual.delta_isc_pct_wt > g542x_dual.delta_isc_pct_wt
        assert g542x_dual.effective_readthrough_prob == 0.095
        assert w1282x_dual.effective_readthrough_prob == 0.185

        # With tri-modal therapy and higher corrector dosing, G542X approaches the therapeutic window
        g542x_tri = model.simulate_condition(
            "G542X",
            smg1i_conc_nm=250.0,
            readthrough_induced=True,
            corrector_co_administered=True,
        )
        assert g542x_tri.delta_isc_pct_wt > 3.0

    def test_escape_matrix_and_reporting(self, model: NMDPtcKineticsModel) -> None:
        """Full escape matrix across all 5 variants generates grounded outputs."""
        matrix = model.generate_escape_matrix()
        assert len(matrix) >= 30

        report = format_markdown_report(matrix)
        assert "NCT04135495" in report
        assert "10.1038/s41580-019-0126-2" in report
        assert "W1282X" in report
        assert "G542X" in report
        assert NLM_ATTRIBUTION in report
        assert RESEARCH_BOUNDARY in report

    def test_time_course_ode_simulation(self, model: NMDPtcKineticsModel) -> None:
        """Time course simulation shows smooth asymptotic mRNA accumulation over 48 hr."""
        time_points = model.simulate_time_course(
            variant_name="W1282X",
            smg1i_conc_nm=250.0,
            readthrough_induced=True,
            t_max_hr=24.0,
            dt_hr=1.0,
        )

        assert len(time_points) == 25  # 0 to 24 hr inclusive
        assert time_points[0]["time_hr"] == 0.0
        assert time_points[-1]["time_hr"] == 24.0

        # Monotonically increasing mRNA level toward steady state
        mrna_levels = [p["mrna_pct_wt"] for p in time_points]
        for i in range(len(mrna_levels) - 1):
            assert mrna_levels[i + 1] >= mrna_levels[i]

        # Initial mRNA is low (~10%), final mRNA is high (> 45%)
        assert mrna_levels[0] < 20.0
        assert mrna_levels[-1] > 45.0
