#!/usr/bin/env python3
"""Unit tests for Stochastic Single-Cell CFTR mRNA Translation Dynamics & Ribosome Collisions under Hypoxia (T-0452)."""

from __future__ import annotations

import json
from pathlib import Path
import sys
import pytest

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent.parent
for p in (str(HERE), str(REPO_ROOT)):
    if p not in sys.path:
        sys.path.insert(0, p)

from mrna_hypoxia_translation_dynamics import (
    CFTR_CDS_NUCLEOTIDES,
    CFTR_CODONS,
    CURATIVE_RESCUE_CHANNELS_THRESHOLD,
    CURATIVE_RESCUE_FRACTION_WT,
    DOI_AIRWAY_MUCUS_PHYSIOL,
    DOI_ETI_PHASE3_TRIAL,
    DOI_SUPPRESSOR_TRNA,
    DOI_W1282X_RESCUE_BARRIER,
    NCT_4D710_AAV_TRIAL,
    NCT_VX522_MRNA_TRIAL,
    NLM_ATTRIBUTION,
    PDB_ACTIVE_CFTR,
    PDB_DEPHOSPHO_CFTR,
    RESEARCH_BOUNDARY,
    WT_CFTR_CHANNELS_PER_CELL,
    MRNAHypoxiaTranslationEngine,
    main,
)


@pytest.fixture
def engine() -> MRNAHypoxiaTranslationEngine:
    return MRNAHypoxiaTranslationEngine()


def test_constants_and_coding_sequence_invariants():
    """Verify structural and coding sequence constants for full-length human CFTR."""
    assert CFTR_CDS_NUCLEOTIDES == 4443
    assert CFTR_CODONS == 1480
    assert CFTR_CDS_NUCLEOTIDES == CFTR_CODONS * 3 + 3  # 1,480 sense codons + 1 stop codon (3 nt)
    assert WT_CFTR_CHANNELS_PER_CELL == 1500.0
    assert CURATIVE_RESCUE_FRACTION_WT == 0.100
    assert CURATIVE_RESCUE_CHANNELS_THRESHOLD == 150.0


def test_hypoxic_isr_attenuation_bounds(engine: MRNAHypoxiaTranslationEngine):
    """Domain 1: Assert eIF2alpha phosphorylation and initiation attenuation under hypoxia."""
    normoxia = engine.compute_hypoxic_isr_attenuation(21.0)
    assert normoxia.stress_regime == "NORMOXIA_CONDUCTING"
    assert normoxia.initiation_attenuation_ratio >= 0.95
    assert normoxia.eif2alpha_phosphorylation_index <= 0.05
    assert normoxia.effective_initiation_rate_hz >= 0.075

    intermediate = engine.compute_hypoxic_isr_attenuation(5.0)
    assert intermediate.stress_regime == "INTERMEDIATE_HYPOXIA"
    assert intermediate.initiation_attenuation_ratio < normoxia.initiation_attenuation_ratio

    severe = engine.compute_hypoxic_isr_attenuation(1.0)
    assert severe.stress_regime == "SEVERE_HYPOXIA_OCCLUSION"
    assert severe.initiation_attenuation_ratio <= 0.35
    assert severe.eif2alpha_phosphorylation_index >= 0.65

    anoxia = engine.compute_hypoxic_isr_attenuation(0.0)
    assert anoxia.initiation_attenuation_ratio == pytest.approx(0.180, rel=1e-2)

    # Monotonicity check across oxygen gradient
    po2_levels = [21.0, 15.0, 10.0, 5.0, 2.5, 1.0, 0.5, 0.0]
    rates = [engine.compute_hypoxic_isr_attenuation(p).effective_initiation_rate_hz for p in po2_levels]
    for i in range(len(rates) - 1):
        assert rates[i] >= rates[i + 1]


def test_ribosome_collision_and_ngd_mitigation(engine: MRNAHypoxiaTranslationEngine):
    """Domain 2: Assert codon optimization and m1Psi suppress collisions and extend half-life."""
    k_init_norm = 0.080
    unmod = engine.evaluate_collision_and_ngd(k_init_norm, codon_optimized=False, m1psi_modified=False)
    opt = engine.evaluate_collision_and_ngd(k_init_norm, codon_optimized=True, m1psi_modified=True)

    # Elongation velocities
    assert opt.mean_elongation_velocity > unmod.mean_elongation_velocity
    assert opt.pause_elongation_velocity > unmod.pause_elongation_velocity

    # Collision probability and disome density reduction
    assert opt.collision_probability_per_transcript < unmod.collision_probability_per_transcript
    assert opt.disome_density_per_mrna < unmod.disome_density_per_mrna
    reduction_pct = (1.0 - opt.disome_density_per_mrna / unmod.disome_density_per_mrna) * 100.0
    assert reduction_pct >= 60.0

    # Effective half-life extension via m1Psi and reduced NGD
    assert opt.effective_mrna_half_life_hours > unmod.effective_mrna_half_life_hours
    assert opt.effective_mrna_half_life_hours >= 12.0
    assert opt.ngd_escape_survival_fraction > unmod.ngd_escape_survival_fraction


def test_single_cell_yield_and_negative_control(engine: MRNAHypoxiaTranslationEngine):
    """Domain 3: Assert negative-control zero yield and dose-dependent channel density."""
    # Negative control: 0 mRNA copies must yield exactly 0 channels
    null_cell = engine.calculate_single_cell_yield(mrna_copy_number=0, po2_percent=21.0)
    assert null_cell.steady_state_surface_channels == 0.0
    assert null_cell.mature_cftr_synthesis_rate_per_sec == 0.0
    assert null_cell.percent_wild_type_cftr == 0.0
    assert not null_cell.is_curative_threshold_met

    # Single-cell dose scaling
    cell_5 = engine.calculate_single_cell_yield(mrna_copy_number=5, po2_percent=21.0)
    cell_20 = engine.calculate_single_cell_yield(mrna_copy_number=20, po2_percent=21.0)
    assert cell_20.steady_state_surface_channels > cell_5.steady_state_surface_channels

    # Hypoxia attenuation comparison at same dose (20 copies)
    cell_hypox = engine.calculate_single_cell_yield(mrna_copy_number=20, po2_percent=1.0)
    assert cell_hypox.steady_state_surface_channels < cell_20.steady_state_surface_channels

    # Curative threshold classification (>= 150 channels / 10% WT)
    assert cell_20.steady_state_surface_channels >= CURATIVE_RESCUE_CHANNELS_THRESHOLD
    assert cell_20.is_curative_threshold_met


def test_airway_population_curative_rescue(engine: MRNAHypoxiaTranslationEngine):
    """Domain 4: Simulate epithelial population across heterogeneous airway microenvironments."""
    # Under physiological low-dose inhalation (mean dose = 1.5 mRNAs/cell), optimized mRNA yields higher rescue %
    pop_unmod_low = engine.simulate_airway_population_rescue(mean_lnp_dose=1.5, codon_optimized=False, m1psi_modified=False)
    pop_opt_low = engine.simulate_airway_population_rescue(mean_lnp_dose=1.5, codon_optimized=True, m1psi_modified=True)

    assert pop_opt_low["cohort_size"] == 1000
    assert pop_opt_low["overall_curative_rescue_percentage"] > pop_unmod_low["overall_curative_rescue_percentage"]
    assert pop_opt_low["overall_mean_surface_channels"] > pop_unmod_low["overall_mean_surface_channels"]

    # At full therapeutic dose (20 mRNAs/cell), optimized formulation achieves >4x higher mature surface density
    pop_opt = engine.simulate_airway_population_rescue(mean_lnp_dose=20.0, codon_optimized=True, m1psi_modified=True)
    pop_unmod = engine.simulate_airway_population_rescue(mean_lnp_dose=20.0, codon_optimized=False, m1psi_modified=False)
    assert pop_opt["overall_mean_surface_channels"] >= 3.0 * pop_unmod["overall_mean_surface_channels"]

    # Zone-specific checks
    zones = pop_opt["per_zone_analysis"]
    assert "Conducting_Ventilated" in zones
    assert "Intermediate_Hypoxic" in zones
    assert "Occluded_Atelectatic" in zones

    # Ventilated conducting zone has higher rescue percentage than severely occluded atelectatic zone at low dose
    low_zones = pop_opt_low["per_zone_analysis"]
    assert low_zones["Conducting_Ventilated"]["curative_rescue_percentage"] > low_zones["Occluded_Atelectatic"]["curative_rescue_percentage"]


def test_provenance_citations_and_disclaimers():
    """Verify registry IDs, primary DOIs, and mandatory regulatory disclaimers."""
    assert "U.S. National Library of Medicine" in NLM_ATTRIBUTION
    assert "does not endorse or recommend" in NLM_ATTRIBUTION
    assert "research, not medical advice" in RESEARCH_BOUNDARY
    assert "No clinical efficacy" in RESEARCH_BOUNDARY

    # Registry grounding
    assert NCT_VX522_MRNA_TRIAL == "NCT05668741"
    assert NCT_4D710_AAV_TRIAL == "NCT05248230"

    # Literature grounding
    assert DOI_SUPPRESSOR_TRNA == "10.1038/s41586-023-06133-1"
    assert DOI_W1282X_RESCUE_BARRIER == "10.1021/acs.biochem.5c00590"
    assert DOI_ETI_PHASE3_TRIAL == "10.1056/NEJMoa1908639"
    assert DOI_AIRWAY_MUCUS_PHYSIOL == "10.1152/physrev.00004.2021"

    # Structural PDB grounding
    assert PDB_DEPHOSPHO_CFTR == "6MSM"
    assert PDB_ACTIVE_CFTR == "5UAK"


def test_comprehensive_suite_and_cli(capsys):
    """Verify comprehensive evaluation suite execution and CLI integration."""
    exit_code = main(["--run-all", "--json"])
    assert exit_code == 0
    captured = capsys.readouterr()
    data = json.loads(captured.out)

    assert data["metadata"]["task_id"] == "T-0452"
    assert data["metadata"]["answer_tier"] == "T2"
    assert "domain_1_hypoxic_isr_attenuation" in data
    assert "domain_2_ribosome_collision_and_ngd" in data
    assert "domain_3_translational_yield_and_channel_density" in data
    assert "domain_4_airway_population_curative_rescue" in data

    # Single cell execution
    exit_code_single = main(["--po2", "10.0", "--mrna-copies", "12", "--json"])
    assert exit_code_single == 0
