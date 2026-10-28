#!/usr/bin/env python3
"""Unit tests for TMEM16A potentiation and SLC26A9 coupling conductance model (T-0476).

Tests allosteric calcium/voltage gating of TMEM16A, SLC26A9 constitutive halide conductance
and PDZ scaffolding kinetics, transepithelial Kirchhoff circuit equilibrium,
ASL fluid dynamics and mucociliary clearance threshold restoration, Henderson-Hasselbalch
pH buffering, scenario benchmarking, and AH-10 falsification negative controls.
"""

from __future__ import annotations

import math
import pytest

from tmem16a_slc26a9_conductance_coupling import (
    ASL_CLEARANCE_THRESHOLD_UM,
    ASL_HEIGHT_COLLAPSED_CF_UM,
    ASL_TARGET_HEALTHY_UM,
    CLINICAL_TRIAL_VANZACAFTOR,
    CLINICAL_TRIAL_VX522,
    DEFAULT_CFTR_DENSITY_WT,
    DEFAULT_CFTR_PO_PKA_ACTIVE,
    DEFAULT_ENAC_DENSITY_CF,
    DEFAULT_ENAC_PO_CF,
    DEFAULT_SLC_DENSITY_CFTR_NULL,
    DEFAULT_SLC_DENSITY_RESCUED,
    DEFAULT_SLC_PO_CONSTITUTIVE,
    DEFAULT_TMEM_DENSITY_TOTAL,
    DEFAULT_TMEM_EC50_BASAL_UM,
    DEFAULT_TMEM_HILL_NH,
    NLM_ATTRIBUTION,
    PDB_STRUCTURE_CFTR,
    PRIMARY_DOI_BUTTON,
    PRIMARY_DOI_MIDDLETON,
    RESEARCH_BOUNDARY,
    benchmark_bypass_regimens,
    calculate_slc26a9_apical_fraction,
    calculate_tmem16a_ec50,
    calculate_tmem16a_open_probability,
    evaluate_falsification_negative_control,
    simulate_asl_hydration,
    simulate_scenario,
    solve_epithelial_circuit,
)


def test_grounding_identifiers_and_disclaimers():
    """Verify primary DOIs, clinical trial IDs, PDB structure, and mandatory disclaimers."""
    assert PRIMARY_DOI_BUTTON == "10.1152/physrev.00004.2021"
    assert PRIMARY_DOI_MIDDLETON == "10.1056/nejmoa1908639"
    assert CLINICAL_TRIAL_VANZACAFTOR == "NCT05033080"
    assert CLINICAL_TRIAL_VX522 == "NCT05668741"
    assert PDB_STRUCTURE_CFTR == "6MSM"
    assert "ClinicalTrials.gov does not evaluate safety or efficacy" in NLM_ATTRIBUTION
    assert "research, not medical advice" in RESEARCH_BOUNDARY


def test_tmem16a_ec50_voltage_and_potentiator_shifts():
    """Verify TMEM16A apparent EC50 modulation by membrane voltage and small-molecule potentiator."""
    ec50_zero_mv = calculate_tmem16a_ec50(membrane_potential_mv=0.0, potentiator_conc_nm=0.0)
    assert math.isclose(ec50_zero_mv, DEFAULT_TMEM_EC50_BASAL_UM, rel_tol=1e-3)

    # Depolarization (positive Vm) left-shifts EC50 to lower concentrations
    ec50_depol = calculate_tmem16a_ec50(membrane_potential_mv=40.0, potentiator_conc_nm=0.0)
    assert ec50_depol < ec50_zero_mv

    # Hyperpolarization (negative Vm) right-shifts EC50 to higher concentrations
    ec50_hyper = calculate_tmem16a_ec50(membrane_potential_mv=-40.0, potentiator_conc_nm=0.0)
    assert ec50_hyper > ec50_zero_mv

    # Small-molecule potentiator (100 nM) substantially left-shifts apparent EC50
    ec50_potentiated = calculate_tmem16a_ec50(membrane_potential_mv=-40.0, potentiator_conc_nm=100.0)
    assert ec50_potentiated < ec50_hyper
    assert ec50_potentiated <= 0.55 * ec50_hyper


def test_tmem16a_hill_open_probability():
    """Verify TMEM16A allosteric Hill activation kinetics across Ca2+ concentrations."""
    # Zero calcium must yield zero open probability
    po_zero = calculate_tmem16a_open_probability(cytosolic_ca_um=0.0, membrane_potential_mv=-40.0)
    assert po_zero == 0.0

    # Basal resting Ca2+ (0.1 uM = 100 nM) gives low basal Po (< 0.10) without potentiator
    po_basal = calculate_tmem16a_open_probability(cytosolic_ca_um=0.10, membrane_potential_mv=-40.0, potentiator_conc_nm=0.0)
    assert 0.01 < po_basal < 0.10

    # Potentiator (100 nM) elevates Po at basal Ca2+ by > 2.5-fold
    po_potentiated = calculate_tmem16a_open_probability(cytosolic_ca_um=0.10, membrane_potential_mv=-40.0, potentiator_conc_nm=100.0)
    assert po_potentiated > po_basal * 2.5
    assert po_potentiated >= 0.12

    # High micromolar calcium (5.0 uM) saturates open probability (> 0.90)
    po_saturating = calculate_tmem16a_open_probability(cytosolic_ca_um=5.0, membrane_potential_mv=-40.0)
    assert po_saturating > 0.90


def test_slc26a9_constitutive_and_pdz_trafficking():
    """Verify SLC26A9 constitutive gating and PDZ scaffold-dependent apical retention."""
    # Constitutive open probability is fixed at ~0.85
    assert DEFAULT_SLC_PO_CONSTITUTIVE == 0.85

    # In CFTR-null baseline without PDZ scaffold, apical retention is basal (~12%)
    f_null = calculate_slc26a9_apical_fraction(pdz_scaffold_conc_um=0.0, has_cftr_scaffold=False)
    assert math.isclose(f_null, 0.12, rel_tol=1e-3)

    # Scaffolding with PDZ chaperone (2.5 uM) rescues apical retention towards maximal
    f_rescued = calculate_slc26a9_apical_fraction(pdz_scaffold_conc_um=2.5, has_cftr_scaffold=False)
    assert f_rescued > 0.50
    assert f_rescued > f_null * 4.0

    # CFTR-anchored scaffold automatically provides maximal stabilization
    f_cftr = calculate_slc26a9_apical_fraction(has_cftr_scaffold=True)
    assert f_cftr >= 0.85


def test_epithelial_circuit_kirchhoff_solution():
    """Verify analytical Kirchhoff circuit equilibrium and current conservation."""
    physio, tmem_po, slc_eff = solve_epithelial_circuit(
        cytosolic_ca_um=0.10,
        potentiator_conc_nm=50.0,
        slc_channel_density=1000.0,
        pdz_scaffold_conc_um=1.5,
    )

    # Membrane potential must fall within physiological range (-65 mV to -15 mV)
    assert -65.0 <= physio.membrane_potential_mv <= -15.0

    # Conductances must be strictly positive
    assert physio.g_cl_total_ps > 0.0
    assert physio.g_hco3_total_ps > 0.0
    assert physio.g_na_total_ps > 0.0
    assert physio.g_apical_total_ps > 0.0

    # ENaC sodium current is inward (negative pA)
    assert physio.i_na_apical_pa < 0.0


def test_asl_fluid_dynamics_collapsing_in_cf():
    """Verify that untreated CF airway epithelia collapse ASL height below 3.5 um."""
    cf_res = simulate_scenario(
        "Untreated CF",
        cftr_density=0.0,
        enac_density=DEFAULT_ENAC_DENSITY_CF,
        enac_p_open=DEFAULT_ENAC_PO_CF,
        slc_channel_density=DEFAULT_SLC_DENSITY_CFTR_NULL,
        pdz_scaffold_conc_um=0.0,
        potentiator_conc_nm=0.0,
        initial_asl_height_um=ASL_HEIGHT_COLLAPSED_CF_UM,
    )

    # Collapsed ASL height < 3.5 um
    assert cf_res.asl_hydration.asl_height_um < 3.5
    # Ciliary clearance incompetent (mucus tethering)
    assert not cf_res.asl_hydration.ciliary_clearance_competent
    # Net solute flux is negative (absorption dominant)
    assert cf_res.electrophysiology.j_osm_solute_pmol_s_cm2 < 0.0
    # Acidic ASL pH < 6.85
    assert cf_res.asl_hydration.asl_ph < 6.85


def test_dual_channel_restoration_above_threshold():
    """Verify that dual-channel coupling (TMEM16A + SLC26A9) restores ASL height >= 7.0 um."""
    dual_res = simulate_scenario(
        "Dual-Channel Bypass",
        cftr_density=0.0,
        enac_density=DEFAULT_ENAC_DENSITY_CF,
        slc_channel_density=DEFAULT_SLC_DENSITY_RESCUED,
        pdz_scaffold_conc_um=2.5,
        potentiator_conc_nm=100.0,
    )

    # ASL height restored above clearance threshold (>= 7.0 um)
    assert dual_res.asl_hydration.asl_height_um >= ASL_CLEARANCE_THRESHOLD_UM
    assert dual_res.asl_hydration.ciliary_clearance_competent
    # Ciliary beat frequency restored (> 8 Hz)
    assert dual_res.asl_hydration.ciliary_beat_frequency_hz >= 8.0
    # Net solute flux is outward secretion (positive)
    assert dual_res.electrophysiology.j_osm_solute_pmol_s_cm2 > 0.0
    # Alkaline buffering restored (pH >= 7.15)
    assert dual_res.asl_hydration.asl_ph >= 7.15


def test_henderson_hasselbalch_asl_ph_buffering():
    """Verify ASL bicarbonate concentration and Henderson-Hasselbalch pH calculations."""
    cf_res = simulate_scenario("CF Baseline", cftr_density=0.0, potentiator_conc_nm=0.0)
    dual_res = simulate_scenario("Dual Channel", cftr_density=0.0, potentiator_conc_nm=100.0, pdz_scaffold_conc_um=2.5, slc_channel_density=2500.0)

    # CF baseline ASL HCO3- is depleted (~6.8 mM, pH ~6.84)
    assert cf_res.asl_hydration.asl_hco3_concentration_mm <= 10.0
    assert cf_res.asl_hydration.asl_ph < 6.90

    # Dual channel bypass secretes HCO3-, raising luminal concentration to > 15 mM and pH >= 7.15
    assert dual_res.asl_hydration.asl_hco3_concentration_mm >= 15.0
    assert dual_res.asl_hydration.asl_ph >= 7.15


def test_multi_scenario_benchmark_comparisons():
    """Verify comparative benchmark across all 6 scenarios."""
    bench = benchmark_bypass_regimens()
    assert len(bench) == 6

    # Order of ASL heights: Untreated CF < TMEM16A alone < Dual Bypass <= Curative Triple Synergy
    h_cf = bench["untreated_cf_baseline"]["asl_height_um"]
    h_tmem = bench["tmem16a_potentiated_alone"]["asl_height_um"]
    h_dual = bench["dual_channel_coupling"]["asl_height_um"]
    h_triple = bench["curative_triple_synergy"]["asl_height_um"]
    h_wt = bench["healthy_wild_type"]["asl_height_um"]

    assert h_cf < h_tmem < h_dual <= h_triple
    assert h_cf < 3.5
    assert h_dual >= 7.0
    assert h_wt >= 7.5


def test_falsification_negative_controls_pass():
    """Verify that AH-10 dynamic falsification suite returns status=PASS."""
    res = evaluate_falsification_negative_control()
    assert res["status"] == "PASS"
    assert res["all_negative_controls_satisfied"] is True
    assert res["cf_collapsed_asl_um"] < 3.5
    assert res["zero_ca_tmem16a_po"] == 0.0
    assert res["dual_bypass_asl_height_um"] >= ASL_CLEARANCE_THRESHOLD_UM
    assert res["ciliary_clearance_restored"] is True
