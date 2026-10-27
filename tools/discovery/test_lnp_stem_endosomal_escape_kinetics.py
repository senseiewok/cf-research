#!/usr/bin/env python3
"""Unit tests for Biophysical Model of Inhaled CFTR mRNA-LNP Cytosolic Delivery, Endosomal Escape & Basal Stem Cell Kinetics (T-0455)."""

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

from lnp_stem_endosomal_escape_kinetics import (
    CFTR_CDS_NUCLEOTIDES,
    CFTR_CODONS,
    CURATIVE_RESCUE_CHANNELS_THRESHOLD,
    CURATIVE_RESCUE_FRACTION_WT,
    DEFAULT_IONIZABLE_LIPID_PKA,
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
    LNPStemEndosomalEscapeEngine,
    main,
)


@pytest.fixture
def engine() -> LNPStemEndosomalEscapeEngine:
    return LNPStemEndosomalEscapeEngine()


def test_constants_and_coding_sequence_invariants():
    """Verify structural and coding sequence constants for human CFTR."""
    assert CFTR_CDS_NUCLEOTIDES == 4443
    assert CFTR_CODONS == 1480
    assert CFTR_CDS_NUCLEOTIDES == CFTR_CODONS * 3 + 3
    assert WT_CFTR_CHANNELS_PER_CELL == 1500.0
    assert CURATIVE_RESCUE_FRACTION_WT == 0.100
    assert CURATIVE_RESCUE_CHANNELS_THRESHOLD == 150.0
    assert "NCT05668741" in NCT_VX522_MRNA_TRIAL
    assert "NCT05248230" in NCT_4D710_AAV_TRIAL
    assert PDB_DEPHOSPHO_CFTR == "6MSM"
    assert PDB_ACTIVE_CFTR == "5UAK"


def test_endosomal_acidification_and_protonation(engine: LNPStemEndosomalEscapeEngine):
    """Domain 1: Assert pH acidification kinetics and Henderson-Hasselbalch lipid protonation."""
    # Extracellular neutrality
    ph_0 = engine.calculate_endosomal_ph(0.0)
    assert ph_0 == pytest.approx(7.40, abs=1e-2)
    alpha_0 = engine.calculate_protonation_fraction(ph_0)
    assert alpha_0 < 0.15  # Neutral at extracellular pH 7.4

    # Acidification at 1 hr and 3 hr
    ph_1 = engine.calculate_endosomal_ph(1.0)
    ph_3 = engine.calculate_endosomal_ph(3.0)
    assert ph_0 > ph_1 > ph_3
    assert ph_3 <= 5.20

    # Protonation at pKa
    alpha_pka = engine.calculate_protonation_fraction(engine.pka)
    assert alpha_pka == pytest.approx(0.500, abs=1e-3)

    # Lysosomal acidic core
    alpha_lyso = engine.calculate_protonation_fraction(4.80)
    assert alpha_lyso > 0.95

    # Stage evaluation
    early_state = engine.evaluate_acidification_state(0.1)
    assert early_state.compartment_stage == "EARLY_ENDOSOME"

    late_state = engine.evaluate_acidification_state(0.8)
    assert late_state.compartment_stage == "LATE_ENDOSOME_ESCAPE_WINDOW"
    assert late_state.h_ii_phase_propensity > early_state.h_ii_phase_propensity
    assert late_state.membrane_disruption_rate_per_hr > early_state.membrane_disruption_rate_per_hr


def test_helper_lipid_phase_propensity_comparison():
    """Domain 1: Verify DOPE (inverted cone) promotes H_II phase compared to DSPC (lamellar stabilizer)."""
    engine_dope = LNPStemEndosomalEscapeEngine(dope_molar_fraction=0.40, dspc_molar_fraction=0.00)
    engine_dspc = LNPStemEndosomalEscapeEngine(dope_molar_fraction=0.00, dspc_molar_fraction=0.40)

    # In late endosome at pH 5.5
    h_ii_dope = engine_dope.compute_h_ii_phase_propensity(5.5)
    h_ii_dspc = engine_dspc.compute_h_ii_phase_propensity(5.5)

    assert h_ii_dope > h_ii_dspc
    assert (h_ii_dope / h_ii_dspc) > 1.10


def test_cytosolic_cargo_delivery_yield(engine: LNPStemEndosomalEscapeEngine):
    """Domain 2: Assert compartmental ODE trafficking and intact mRNA cytosolic yield."""
    initial = 25.0
    yield_res = engine.simulate_intracellular_trafficking(initial_mRNAs=initial, total_time_hr=24.0)

    assert yield_res.initial_internalized_copies == initial
    assert yield_res.cytosolic_intact_copies > 0.5
    assert 0.03 <= yield_res.escape_efficiency_fraction <= 0.25
    assert 0.5 <= yield_res.time_to_peak_cytosolic_hr <= 8.0
    assert yield_res.lysosomal_degraded > 0.0
    assert yield_res.recycled_exocytosed > 0.0


def test_subpopulation_endocytic_kinetics(engine: LNPStemEndosomalEscapeEngine):
    """Domain 3: Assert cell-type differences across luminal, active basal, and quiescent basal stem cells."""
    lum = engine.evaluate_subpopulation_kinetics("ciliated_luminal")
    active_b = engine.evaluate_subpopulation_kinetics("active_basal_stem")
    quiesc_b = engine.evaluate_subpopulation_kinetics("quiescent_basal_stem")

    # Volume differences
    assert lum.relative_cell_volume_um3 > active_b.relative_cell_volume_um3 > quiesc_b.relative_cell_volume_um3

    # Endocytic and translation differences
    assert lum.paracellular_access_fraction == 1.000
    assert active_b.paracellular_access_fraction < 1.000
    assert active_b.translation_initiation_rate_hz > quiesc_b.translation_initiation_rate_hz
    assert active_b.endocytic_uptake_rate_per_hr > quiesc_b.endocytic_uptake_rate_per_hr

    # Curative threshold evaluation
    assert lum.achieves_curative_threshold is True
    assert active_b.achieves_curative_threshold is True


def test_paracellular_tight_junction_barrier_modulation(engine: LNPStemEndosomalEscapeEngine):
    """Domain 3: Assert TEER reduction (transient tight junction opening) increases basal stem access."""
    intact_teer = engine.evaluate_subpopulation_kinetics("active_basal_stem", paracellular_ter_ohm_cm2=1500.0)
    permeabilized_teer = engine.evaluate_subpopulation_kinetics("active_basal_stem", paracellular_ter_ohm_cm2=450.0)

    assert permeabilized_teer.paracellular_access_fraction > intact_teer.paracellular_access_fraction
    assert permeabilized_teer.peak_channels_per_cell > intact_teer.peak_channels_per_cell


def test_epithelial_turnover_and_curative_longevity(engine: LNPStemEndosomalEscapeEngine):
    """Domain 4: Assert basal stem cell transduction extends curative longevity across renewal cycles."""
    trajectory = engine.simulate_epithelial_turnover_trajectory(days_total=90.0)

    assert len(trajectory.days) == 91
    assert trajectory.curative_threshold_channels == 150.0

    # Luminal-only channels decay rapidly
    assert trajectory.luminal_therapeutic_duration_days <= 14.0

    # Basal-anchored channels persist across epithelial renewal cycles (30-90 days)
    assert trajectory.basal_anchored_therapeutic_duration_days >= 30.0

    # Stem cell longevity benefit ratio
    assert trajectory.curative_longevity_ratio >= 2.0


def test_full_simulation_and_cli(engine: LNPStemEndosomalEscapeEngine, monkeypatch):
    """Assert full simulation runner and CLI interface."""
    sim = engine.run_full_simulation()

    assert sim["model"] == "lnp_stem_endosomal_escape_kinetics"
    assert sim["task"] == "T-0455"
    assert sim["answer_tier"] == "T2"
    assert len(sim["primary_evidence"]["dois"]) == 4
    assert len(sim["primary_evidence"]["clinical_trials"]) == 2
    assert len(sim["primary_evidence"]["structural_pdbs"]) == 2

    # CLI test with --json
    monkeypatch.setattr(sys, "argv", ["lnp_stem_endosomal_escape_kinetics.py", "--json"])
    assert main() == 0

    # CLI test default output
    monkeypatch.setattr(sys, "argv", ["lnp_stem_endosomal_escape_kinetics.py"])
    assert main() == 0


def test_invalid_parameters_raise(engine: LNPStemEndosomalEscapeEngine):
    """Assert parameter bounds checking and exception handling."""
    with pytest.raises(ValueError):
        LNPStemEndosomalEscapeEngine(ionizable_pka=4.5)  # Below 5.0

    with pytest.raises(ValueError):
        engine.calculate_endosomal_ph(-1.0)

    with pytest.raises(ValueError):
        engine.simulate_intracellular_trafficking(initial_mRNAs=0.0)

    with pytest.raises(ValueError):
        engine.evaluate_subpopulation_kinetics("invalid_cell_type")
