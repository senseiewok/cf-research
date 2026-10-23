#!/usr/bin/env python3
"""Comprehensive unit test suite for CFTR PDI Redox Remodeling & S-Nitrosylation Model (T-0281)."""
from __future__ import annotations

import json
import pytest
from pathlib import Path

from pdi_redox_remodeling import (
    CFTR_CYSTEINE_CATALOG,
    CFTR_VARIANT_CATALOG,
    MODULATOR_PROFILES,
    NLM_ATTRIBUTION,
    RESEARCH_BOUNDARY,
    calculate_gsno_occupancy,
    calculate_pdi_remodeling_rate,
    estimate_sweat_chloride_equivalent,
    simulate_pdi_redox_remodeling,
    simulate_dose_response_curve,
    simulate_combination_synergy,
    run_full_evaluation,
    PDIRedoxState,
)


def test_catalogs_integrity():
    """Verify primary literature catalogs and disclaimers."""
    assert len(CFTR_CYSTEINE_CATALOG) >= 5
    assert "Cys524" in CFTR_CYSTEINE_CATALOG
    assert "Cys1355" in CFTR_CYSTEINE_CATALOG
    assert CFTR_CYSTEINE_CATALOG["Cys524"].is_cryptic is True
    assert CFTR_CYSTEINE_CATALOG["Cys1355"].is_cryptic is True

    assert "F508del" in CFTR_VARIANT_CATALOG
    assert "N1303K" in CFTR_VARIANT_CATALOG
    assert "G85E" in CFTR_VARIANT_CATALOG
    assert "WT" in CFTR_VARIANT_CATALOG

    # Verify primary DOIs and NCTs
    for var_name, profile in CFTR_VARIANT_CATALOG.items():
        assert profile.primary_doi.startswith("10.")
        assert profile.clinical_nct.startswith("NCT")

    # Mandatory disclaimers
    assert "ClinicalTrials.gov does not evaluate safety or efficacy" in NLM_ATTRIBUTION
    assert "Preclinical in-silico" in RESEARCH_BOUNDARY


def test_gsno_occupancy_hill_equation():
    """Verify Hill-Langmuir binding equation properties."""
    # Zero dose -> 0 occupancy
    assert calculate_gsno_occupancy(0.0) == 0.0
    assert calculate_gsno_occupancy(-5.0) == 0.0

    # Monotonicity
    doses = [0.1, 0.5, 1.0, 2.5, 5.0, 10.0, 50.0, 100.0]
    occupancies = [calculate_gsno_occupancy(d) for d in doses]
    for i in range(len(occupancies) - 1):
        assert occupancies[i] < occupancies[i + 1]

    # At Kd = 2.5 uM, occupancy should be exactly 50%
    occ_kd = calculate_gsno_occupancy(2.5, kd_um=2.5)
    assert pytest.approx(occ_kd, abs=0.01) == 0.50

    # Saturation at high doses
    assert pytest.approx(calculate_gsno_occupancy(1000.0), abs=0.01) == 1.0


def test_pdi_remodeling_kinetics():
    """Verify Michaelis-Menten kinetics for PDI disulfide remodeling."""
    pdi_state = PDIRedoxState(
        pdi_concentration_um=20.0,
        active_oxidoreductase_fraction=0.90,
        sno_pdi_fraction=0.30,
        kcat_isomerization_s_inv=1.20,
        km_substrate_um=15.0,
    )
    # Zero substrate -> 0 flux
    assert calculate_pdi_remodeling_rate(pdi_state, substrate_concentration_um=0.0) == 0.0

    # Non-zero flux
    flux = calculate_pdi_remodeling_rate(pdi_state, substrate_concentration_um=5.0)
    assert flux > 0.0
    # Michaelis-Menten: v = Vmax * [S] / (Km + [S]) = 1.2 * 18 * 5 / 20 = 5.4 uM/s
    assert pytest.approx(flux, abs=0.1) == 5.4


def test_f508del_rescue_simulation():
    """Verify GSNO rescue of F508del CFTR per Zaman et al. 2001."""
    # Baseline (no GSNO)
    res_base = simulate_pdi_redox_remodeling("F508del", gsno_dose_um=0.0)
    assert res_base.theta_sno_occupancy == 0.0
    assert pytest.approx(res_base.surface_rescue_fraction, abs=0.005) == 0.012  # ~1.2%
    assert res_base.delta_sweat_chloride_mmol_l == 0.0
    assert res_base.fold_maturation_enhancement == 1.0

    # Therapeutic dose 10 uM GSNO
    res_treated = simulate_pdi_redox_remodeling("F508del", gsno_dose_um=10.0)
    assert res_treated.theta_sno_occupancy > 0.70
    # Zaman et al. 2001 reported ~3-5 fold Band C increase
    assert res_treated.fold_maturation_enhancement >= 2.5
    assert res_treated.surface_rescue_fraction > res_base.surface_rescue_fraction
    assert res_treated.effective_k_erad_s_inv < res_base.effective_k_erad_s_inv
    assert res_treated.effective_k_fold_s_inv > res_base.effective_k_fold_s_inv
    # Significant sweat chloride normalization
    assert res_treated.delta_sweat_chloride_mmol_l <= -5.0


def test_multi_variant_profiles():
    """Verify comparative response across WT, F508del, N1303K, G85E."""
    for variant in ["WT", "F508del", "N1303K", "G85E"]:
        res = simulate_pdi_redox_remodeling(variant, gsno_dose_um=10.0)
        assert res.variant == variant
        assert 0.0 <= res.surface_rescue_fraction <= 1.0
        assert res.effective_k_fold_s_inv > 0.0
        assert res.effective_k_erad_s_inv > 0.0

    # G85E is known to be more refractory than F508del to chaperone remodeling
    res_g85e = simulate_pdi_redox_remodeling("G85E", gsno_dose_um=10.0)
    res_f508 = simulate_pdi_redox_remodeling("F508del", gsno_dose_um=10.0)
    assert res_g85e.surface_rescue_fraction < res_f508.surface_rescue_fraction


def test_dose_response_curve():
    """Verify dose-response curve generation and monotonic surface rescue."""
    curve = simulate_dose_response_curve("F508del", doses_um=[0.0, 1.0, 5.0, 25.0])
    assert len(curve) == 4
    for i in range(len(curve) - 1):
        assert curve[i].surface_rescue_fraction <= curve[i + 1].surface_rescue_fraction
        assert curve[i].delta_sweat_chloride_mmol_l >= curve[i + 1].delta_sweat_chloride_mmol_l


def test_combination_synergy():
    """Verify Bliss independence and synergy with CFTR correctors."""
    # Test combination of GSNO with Tezacaftor (VX-661)
    syn = simulate_combination_synergy("F508del", gsno_dose_um=10.0, modulator_key="tezacaftor_vx661")
    assert syn.is_synergistic is True
    assert syn.combination_surface_fraction > syn.single_agent_gsno_surface_fraction
    assert syn.combination_surface_fraction > syn.single_agent_modulator_surface_fraction
    assert syn.delta_sweat_chloride_mmol_l <= -15.0  # Clinically meaningful sweat chloride drop


def test_sweat_chloride_transfer_function():
    """Verify transfer function boundaries."""
    # Baseline -> 0 change
    assert estimate_sweat_chloride_equivalent(0.012, baseline_surface_fraction=0.012) == 0.0
    # Lower than baseline -> 0 change (no worsening beyond baseline modeled)
    assert estimate_sweat_chloride_equivalent(0.005, baseline_surface_fraction=0.012) == 0.0
    # Full WT restoration -> ~ -58 mmol/L
    sw_full = estimate_sweat_chloride_equivalent(0.28, baseline_surface_fraction=0.012, wt_surface_fraction=0.28)
    assert pytest.approx(sw_full, abs=1.0) == -58.0


def test_run_full_evaluation():
    """Verify full system evaluation generates valid report data."""
    output = run_full_evaluation()
    assert output["model"] == "pdi_redox_remodeling"
    assert output["task_id"] == "T-0281"
    assert "F508del" in output["variants"]
    assert len(output["synergy_combinations"]) >= 3
    # Check JSON serializability
    serialized = json.dumps(output)
    assert len(serialized) > 500


def test_invalid_inputs_raise_errors():
    """Verify proper error handling for invalid variants and modulators."""
    with pytest.raises(ValueError, match="Unknown CFTR variant"):
        simulate_pdi_redox_remodeling("INVALID_MUTATION", gsno_dose_um=10.0)

    with pytest.raises(ValueError, match="Unknown modulator"):
        simulate_combination_synergy("F508del", gsno_dose_um=10.0, modulator_key="INVALID_DRUG")
