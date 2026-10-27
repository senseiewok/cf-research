#!/usr/bin/env python3
"""Unit tests for Cross-Scale Biophysical Modeling of ASL Acidification and Antimicrobial Potency (T-0469).

Validates:
1. Grounding identifiers (DOIs, NCT trials, PDB structures, NLM/Research disclaimers).
2. CFTR GHK electrodiffusion and Henderson-Hasselbalch ASL pH regulation.
3. pH-dependent protonation states and net charge for hBD-1 and LL-37.
4. Bicarbonate-dependent alpha-helical folding of LL-37 and ionic screening attenuation.
5. Pseudomonas aeruginosa bactericidal kill kinetics and AH-10 negative control.
6. Small-molecule modulator (ETI / Trikafta, VTD) and mRNA LNP (VX-522) rescue.
7. Topical aerosolized bicarbonate buffer clearance kinetics.
8. Discovery dossier schema and CLI execution.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys
import pytest

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import asl_ph_antimicrobial_potency as ap


@pytest.fixture
def engine() -> ap.AslPhAntimicrobialPotencyEngine:
    return ap.AslPhAntimicrobialPotencyEngine()


def test_grounding_identifiers_and_disclaimers():
    """Verify primary DOIs, clinical registry NCTs, PDB codes, and mandatory attributions."""
    assert ap.CLINICAL_TRIAL_VX522 == "NCT05668741"
    assert ap.CLINICAL_TRIAL_VANZACAFTOR == "NCT05033080"
    assert ap.PDB_DEPHOSPHO_CFTR == "6MSM"
    assert ap.PDB_ACTIVE_CFTR == "5UAK"

    assert ap.DOI_AIRWAY_MUCUS == "10.1152/physrev.00004.2021"
    assert ap.DOI_ASL_REGULATION == "10.1164/rccm.202102-0461oc"
    assert ap.DOI_ETI_TRIAL == "10.1056/nejmoa1908639"

    assert "ClinicalTrials.gov does not evaluate safety or efficacy" in ap.NLM_ATTRIBUTION
    assert "research, not medical advice" in ap.RESEARCH_BOUNDARY
    assert ap.CURATIVE_LOG10_CFU_REDUCTION_BENCHMARK == 3.0


def test_ghk_and_henderson_hasselbalch_ph_regulation(engine: ap.AslPhAntimicrobialPotencyEngine):
    """Verify GHK bicarbonate electrodiffusion and steady-state ASL pH in healthy vs CF."""
    wt_state = engine.calculate_steady_state_asl_ph(cftr_activity_fraction=1.0)
    cf_state = engine.calculate_steady_state_asl_ph(cftr_activity_fraction=0.01)

    # Healthy ASL pH: 7.20 to 7.55; CF ASL pH: 6.50 to 6.80
    assert 7.20 <= wt_state.asl_ph <= 7.55
    assert 6.50 <= cf_state.asl_ph <= 6.80

    # Bicarbonate concentration
    assert wt_state.asl_bicarbonate_mm >= 20.0
    assert cf_state.asl_bicarbonate_mm <= 6.0

    # GHK flux difference
    assert wt_state.hco3_efflux_pmol_cm2_s > cf_state.hco3_efflux_pmol_cm2_s * 10.0

    # ASL height and ionic strength
    assert wt_state.asl_height_microns > cf_state.asl_height_microns * 1.8
    assert cf_state.asl_ionic_strength_mm > wt_state.asl_ionic_strength_mm


def test_peptide_protonation_and_charge(engine: ap.AslPhAntimicrobialPotencyEngine):
    """Verify protonation states and net charge for hBD-1 and LL-37 across physiological pH."""
    hbd1_74 = engine.calculate_hbd1_protonation(asl_ph=7.40)
    hbd1_66 = engine.calculate_hbd1_protonation(asl_ph=6.60)

    assert 4.0 <= hbd1_74.nominal_charge_z <= 6.0
    assert hbd1_74.amphipathic_helix_fraction == 0.0  # beta-sheet defensin
    assert hbd1_74.effective_bactericidal_rate_k > hbd1_66.effective_bactericidal_rate_k * 2.0

    ll37_74 = engine.calculate_ll37_protonation(asl_ph=7.40, bicarbonate_mm=24.5)
    ll37_66 = engine.calculate_ll37_protonation(asl_ph=6.60, bicarbonate_mm=4.2)

    assert 5.0 <= ll37_74.nominal_charge_z <= 7.0
    # LL-37 helical folding is bicarbonate-dependent
    assert ll37_74.amphipathic_helix_fraction > 0.65
    assert ll37_66.amphipathic_helix_fraction < 0.35


def test_ionic_strength_and_screening(engine: ap.AslPhAntimicrobialPotencyEngine):
    """Verify ionic screening and MBC_99.9 shift between healthy (110 mM) and CF (155 mM) salt."""
    hbd1_healthy = engine.calculate_hbd1_protonation(asl_ph=7.40, ionic_strength_mm=110.0)
    hbd1_cf = engine.calculate_hbd1_protonation(asl_ph=6.60, ionic_strength_mm=155.0)

    assert hbd1_healthy.ionic_screening_factor > hbd1_cf.ionic_screening_factor
    assert hbd1_healthy.mbc999_ug_ml < 12.0
    assert hbd1_cf.mbc999_ug_ml > 25.0

    ll37_healthy = engine.calculate_ll37_protonation(asl_ph=7.40, bicarbonate_mm=24.5, ionic_strength_mm=110.0)
    ll37_cf = engine.calculate_ll37_protonation(asl_ph=6.60, bicarbonate_mm=4.2, ionic_strength_mm=155.0)

    assert ll37_healthy.mbc999_ug_ml < 10.0
    assert ll37_cf.mbc999_ug_ml > 20.0


def test_bactericidal_clearance_and_ah10_negative_control(engine: ap.AslPhAntimicrobialPotencyEngine):
    """Rule AH-10: Untreated CF must fail bactericidal threshold; healthy control must pass."""
    cf_hbd1 = engine.evaluate_bactericidal_potency("hBD-1", cftr_activity_fraction=0.01)
    cf_ll37 = engine.evaluate_bactericidal_potency("LL-37", cftr_activity_fraction=0.01)

    # AH-10 Negative Control: In untreated CF, endogenous peptides fail to achieve 3-log kill
    assert cf_hbd1.bactericidal_threshold_met is False
    assert cf_hbd1.log10_cfu_reduction_60min < 1.5
    assert cf_ll37.bactericidal_threshold_met is False
    assert cf_ll37.log10_cfu_reduction_60min < 1.5

    # Healthy positive control
    wt_hbd1 = engine.evaluate_bactericidal_potency("hBD-1", cftr_activity_fraction=1.0)
    wt_ll37 = engine.evaluate_bactericidal_potency("LL-37", cftr_activity_fraction=1.0)

    assert wt_hbd1.bactericidal_threshold_met is True
    assert wt_hbd1.log10_cfu_reduction_60min >= 3.0
    assert wt_ll37.bactericidal_threshold_met is True
    assert wt_ll37.log10_cfu_reduction_60min >= 3.0


def test_modulator_regimens_rescue(engine: ap.AslPhAntimicrobialPotencyEngine):
    """Verify small-molecule modulators (Trikafta/ETI, VTD) and mRNA LNP rescue."""
    res = engine.simulate_modulator_regimens()

    assert "UNTREATED_CF" in res
    assert "TRIKAFTA_ETI" in res
    assert "VANZACAFTOR_VTD" in res
    assert "INHALED_MRNA_VX522" in res
    assert "HEALTHY_CONTROL" in res

    # Untreated CF fails
    assert res["UNTREATED_CF"]["curative_clearance_supported"] is False

    # Trikafta (ETI) recovers pH >= 7.20 and achieves threshold
    eti = res["TRIKAFTA_ETI"]
    assert eti["asl_ph"] >= 7.20
    assert eti["hbd1_threshold_met"] is True
    assert eti["ll37_threshold_met"] is True
    assert eti["curative_clearance_supported"] is True

    # Next-generation VTD achieves higher pH than Trikafta
    vtd = res["VANZACAFTOR_VTD"]
    assert vtd["asl_ph"] > eti["asl_ph"]
    assert vtd["curative_clearance_supported"] is True

    # mRNA LNP achieves therapeutic recovery
    mrna = res["INHALED_MRNA_VX522"]
    assert mrna["asl_ph"] >= 7.15


def test_topical_aerosol_bicarbonate_kinetics(engine: ap.AslPhAntimicrobialPotencyEngine):
    """Verify acute alkaline peak and exponential re-acidification kinetics of topical NaHCO3."""
    traj = engine.simulate_topical_bicarbonate_aerosol(
        peak_bicarb_mm=48.0,
        clearance_half_life_min=42.0,
        duration_min=180.0,
        time_step_min=15.0,
    )

    assert len(traj) == 13  # 0 to 180 min by 15 min steps

    # Initial time point t = 0 min
    t0 = traj[0]
    assert t0.time_minutes == 0.0
    assert t0.asl_ph >= 7.60
    assert t0.bactericidal_threshold_met is True

    # After 90 min, clearance causes pH drop
    t90 = next(pt for pt in traj if pt.time_minutes == 90.0)
    assert t90.asl_ph < t0.asl_ph

    # By 180 min, re-acidification returns near CF baseline
    t180 = traj[-1]
    assert t180.asl_ph <= 6.85
    assert t180.bactericidal_threshold_met is False


def test_discovery_dossier_and_cli():
    """Verify discovery dossier structure and CLI exit code."""
    dossier = ap.build_discovery_dossier()
    assert dossier["task_id"] == "T-0469"
    assert dossier["status"] == "APPROVED"
    assert dossier["answer_tier"] == "T2"

    assert len(dossier["primary_evidence_citations"]) >= 5
    assert "steady_state_asl_ph" in dossier["biophysical_findings"]
    assert "antimicrobial_peptide_inactivation" in dossier["biophysical_findings"]

    # Verify CLI execution
    exit_code = ap.main(["--json"])
    assert exit_code == 0
