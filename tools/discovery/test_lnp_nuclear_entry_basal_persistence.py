#!/usr/bin/env python3
"""Unit tests for LNP Nuclear Entry and Basal Stem Cell Genomic Persistence Kinetics Model (T-0478)."""

from __future__ import annotations

import json
import sys
from pathlib import Path
import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
DISCOVERY_DIR = Path(__file__).resolve().parent
if str(DISCOVERY_DIR) not in sys.path:
    sys.path.insert(0, str(DISCOVERY_DIR))

try:
    from tools.discovery.lnp_nuclear_entry_basal_persistence import (
        CURATIVE_THRESHOLD_CFTR_PCT,
        DEFAULT_RNP_MOLECULAR_WEIGHT_KDA,
        NLM_ATTRIBUTION,
        PRIMARY_DOI_BUTTON,
        PRIMARY_DOI_COONEY,
        PRIMARY_DOI_MIDDLETON,
        CLINICAL_TRIAL_VX522,
        CLINICAL_TRIAL_VANZACAFTOR,
        PDB_STRUCTURE_CFTR,
        RESEARCH_BOUNDARY,
        EpithelialNicheParameters,
        LnpNuclearEntryBasalPersistenceModel,
        NpcTranslocationParameters,
        TurnoverKineticsParameters,
    )
except ImportError:
    from lnp_nuclear_entry_basal_persistence import (
        CURATIVE_THRESHOLD_CFTR_PCT,
        DEFAULT_RNP_MOLECULAR_WEIGHT_KDA,
        NLM_ATTRIBUTION,
        PRIMARY_DOI_BUTTON,
        PRIMARY_DOI_COONEY,
        PRIMARY_DOI_MIDDLETON,
        CLINICAL_TRIAL_VX522,
        CLINICAL_TRIAL_VANZACAFTOR,
        PDB_STRUCTURE_CFTR,
        RESEARCH_BOUNDARY,
        EpithelialNicheParameters,
        LnpNuclearEntryBasalPersistenceModel,
        NpcTranslocationParameters,
        TurnoverKineticsParameters,
    )


@pytest.fixture
def default_model() -> LnpNuclearEntryBasalPersistenceModel:
    return LnpNuclearEntryBasalPersistenceModel()


def test_stokes_hydrodynamic_radius_calculation(default_model: LnpNuclearEntryBasalPersistenceModel) -> None:
    """Verifies Stokes-Einstein hydrodynamic radius calculation and physical scaling."""
    # Standard Cas9 / base editor RNP (~160 kDa)
    r_h = default_model.compute_stokes_hydrodynamic_radius_nm(160.0)
    assert 4.0 <= r_h <= 4.8, f"Expected Stokes radius ~4.4 nm for 160 kDa RNP, got {r_h}"

    # Small protein (~27 kDa GFP)
    r_gfp = default_model.compute_stokes_hydrodynamic_radius_nm(27.0)
    assert 2.0 <= r_gfp <= 2.6, f"Expected Stokes radius ~2.3 nm for GFP, got {r_gfp}"

    # Error handling for invalid mass
    with pytest.raises(ValueError):
        default_model.compute_stokes_hydrodynamic_radius_nm(-10.0)


def test_npc_passive_vs_active_translocation(default_model: LnpNuclearEntryBasalPersistenceModel) -> None:
    """Tests passive exclusion of large RNPs and importin-facilitated active translocation."""
    # Active import with NLS
    active_res = default_model.solve_npc_translocation(
        NpcTranslocationParameters(
            molecular_weight_kda=DEFAULT_RNP_MOLECULAR_WEIGHT_KDA,
            has_nls=True,
            nls_kd_nm=25.0,
        )
    )
    assert not active_res.is_passive_permissive
    assert active_res.importin_bound_fraction > 0.85
    assert active_res.active_flux_molecules_per_s > 1000.0
    assert active_res.nuclear_cytoplasmic_steady_state_ratio > 10.0
    assert "Importin-alpha/beta" in active_res.translocation_mechanism

    # Cargo without NLS
    passive_res = default_model.solve_npc_translocation(
        NpcTranslocationParameters(
            molecular_weight_kda=DEFAULT_RNP_MOLECULAR_WEIGHT_KDA,
            has_nls=False,
        )
    )
    assert not passive_res.is_passive_permissive
    assert passive_res.importin_bound_fraction == 0.0
    assert passive_res.active_flux_molecules_per_s == 0.0
    assert passive_res.nuclear_cytoplasmic_steady_state_ratio < 1.0
    assert "Passive diffusion" in passive_res.translocation_mechanism


def test_paracellular_basal_niche_delivery(default_model: LnpNuclearEntryBasalPersistenceModel) -> None:
    """Evaluates tight junction barrier and permeabilized paracellular delivery into basal stem cells."""
    # Intact epithelium (TEER 600 Ohm*cm^2)
    intact_res = default_model.solve_paracellular_basal_niche_delivery(
        EpithelialNicheParameters(teer_ohm_cm2=600.0, is_permeabilized=False)
    )
    assert intact_res.paracellular_penetration_fraction < 0.03
    assert intact_res.apical_columnar_uptake_fraction > 0.97
    assert "Intact Epithelial Barrier" in intact_res.permeabilization_state

    # Permeabilized epithelium (e.g. transient tight junction disassembly)
    perm_res = default_model.solve_paracellular_basal_niche_delivery(
        EpithelialNicheParameters(teer_ohm_cm2=80.0, is_permeabilized=True)
    )
    assert 0.15 <= perm_res.paracellular_penetration_fraction <= 0.38
    assert perm_res.basal_stem_uptake_fraction >= 0.15
    assert perm_res.basal_stem_dose_per_cell > intact_res.basal_stem_dose_per_cell * 5.0
    assert "Transiently Permeabilized" in perm_res.permeabilization_state


def test_lineage_turnover_genomic_persistence(default_model: LnpNuclearEntryBasalPersistenceModel) -> None:
    """Verifies multi-year persistence of genomic edits in basal stem cells vs rapid decay of episomal mRNA."""
    trajectory = default_model.simulate_turnover_and_persistence(duration_days=3650.0)  # 10 years
    assert len(trajectory) > 50

    p0 = trajectory[0]
    p_30d = min(trajectory, key=lambda p: abs(p.day - 30.0))
    p_90d = min(trajectory, key=lambda p: abs(p.day - 90.0))
    p_1yr = min(trajectory, key=lambda p: abs(p.day - 365.0))
    p_5yr = min(trajectory, key=lambda p: abs(p.day - 1825.0))
    p_10yr = min(trajectory, key=lambda p: abs(p.day - 3650.0))

    # Day 0: both modalities show strong initial expression
    assert p0.episomal_cftr_pct > 50.0
    assert p0.genomic_total_cftr_pct > 25.0

    # Day 30: episomal mRNA decays significantly
    assert p_30d.episomal_cftr_pct < 5.0

    # Day 90: episomal mRNA vanishes (< 0.1%), dropping below curative threshold
    assert p_90d.episomal_cftr_pct < 0.1
    assert not p_90d.episomal_above_curative

    # Multi-year genomic persistence: basal stem cell edits sustain the lineage
    # Asymptotic floor is stably maintained >= 10% WT curative threshold
    assert p_1yr.genomic_total_cftr_pct >= CURATIVE_THRESHOLD_CFTR_PCT
    assert p_5yr.genomic_total_cftr_pct >= CURATIVE_THRESHOLD_CFTR_PCT
    assert p_10yr.genomic_total_cftr_pct >= CURATIVE_THRESHOLD_CFTR_PCT
    assert p_10yr.genomic_above_curative


def test_durability_comparison_structure(default_model: LnpNuclearEntryBasalPersistenceModel) -> None:
    """Verifies modality durability comparison output and clinical attributes."""
    trajectory = default_model.simulate_turnover_and_persistence(duration_days=3650.0)
    summaries = default_model.generate_durability_comparison(trajectory)
    assert len(summaries) == 2

    genomic_mod, episomal_mod = summaries[0], summaries[1]

    # Genomic editing checks
    assert "Base / Prime Editing" in genomic_mod.modality_name
    assert not genomic_mod.requires_indefinite_redosing
    assert genomic_mod.curative_duration_days >= 3650.0
    assert genomic_mod.asymptotic_homeostatic_floor_pct >= CURATIVE_THRESHOLD_CFTR_PCT
    assert PRIMARY_DOI_COONEY in genomic_mod.clinical_reference

    # Episomal mRNA checks
    assert "Episomal mRNA" in episomal_mod.modality_name
    assert episomal_mod.requires_indefinite_redosing
    assert episomal_mod.curative_duration_days < 30.0
    assert episomal_mod.asymptotic_homeostatic_floor_pct == 0.0
    assert CLINICAL_TRIAL_VX522 in episomal_mod.clinical_reference


def test_full_simulation_report_and_disclaimers(default_model: LnpNuclearEntryBasalPersistenceModel) -> None:
    """Verifies full simulation report structure, primary groundings, and mandatory disclaimers."""
    report = default_model.run_full_simulation_report(duration_days=3650.0)

    # Metadata & primary grounding
    assert report["model_metadata"]["task_id"] == "T-0478"
    assert report["model_metadata"]["track"] == "Track A"
    assert report["model_metadata"]["discipline"] == "airway_delivery"
    assert report["model_metadata"]["answer_tier"] == "T2"

    grounding = report["primary_grounding"]
    dois = [lit["doi"] for lit in grounding["literature"]]
    assert PRIMARY_DOI_BUTTON in dois
    assert PRIMARY_DOI_COONEY in dois
    assert PRIMARY_DOI_MIDDLETON in dois

    ncts = [ct["nct"] for ct in grounding["clinical_trials"]]
    assert CLINICAL_TRIAL_VX522 in ncts
    assert CLINICAL_TRIAL_VANZACAFTOR in ncts
    assert grounding["pdb_structural_coordinate"] == PDB_STRUCTURE_CFTR

    # Disclaimers
    assert report["disclaimers"]["nlm_attribution"] == NLM_ATTRIBUTION
    assert report["disclaimers"]["research_boundary"] == RESEARCH_BOUNDARY


def test_cli_execution_json(monkeypatch, capsys) -> None:
    """Tests CLI entrypoint with JSON output."""
    from tools.discovery.lnp_nuclear_entry_basal_persistence import main
    monkeypatch.setattr("sys.argv", ["lnp_nuclear_entry_basal_persistence.py", "--json", "--duration-days", "180"])
    exit_code = main()
    assert exit_code == 0
    captured = capsys.readouterr()
    data = json.loads(captured.out)
    assert data["model_metadata"]["task_id"] == "T-0478"
