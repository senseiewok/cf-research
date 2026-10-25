#!/usr/bin/env python3
"""Unit tests for tools/discovery/enac_proteolysis_asl_regulation.py (T-0417)."""
import math
import pytest

from enac_proteolysis_asl_regulation import (
    CRITICAL_CILIARY_HEIGHT_UM,
    PO_ALPHA_FURIN_CLEAVED,
    PO_GAMMA_ELASTASE_ACTIVE,
    PO_UNCLEAVED_SILENT,
    ENaCASLParameters,
    calculate_enac_open_probability,
    calculate_gamma_cleavage_rate,
    compare_enac_therapeutic_regimens,
    simulate_enac_asl_regulation,
)


def test_michaelis_menten_elastase_cleavage_saturation():
    """Verify gamma-ENaC cleavage rate increases monotonically and saturates at high elastase."""
    k_zero = calculate_gamma_cleavage_rate(neutrophil_elastase_nm=0.0, prostasin_activity_fraction=0.5)
    k_mid = calculate_gamma_cleavage_rate(neutrophil_elastase_nm=80.0, prostasin_activity_fraction=0.5, elastase_km_nm=80.0, kcat_elastase_min_inv=0.20)
    k_high = calculate_gamma_cleavage_rate(neutrophil_elastase_nm=8000.0, prostasin_activity_fraction=0.5, elastase_km_nm=80.0, kcat_elastase_min_inv=0.20)

    assert k_zero > 0
    assert k_mid > k_zero
    # At [NE] == Km, elastase contribution equals 0.5 * kcat = 0.10
    assert math.isclose(k_mid - k_zero, 0.10, rel_tol=1e-4)
    assert k_high > k_mid
    assert k_high < k_zero + 0.20


def test_open_probability_bounded_between_silent_and_constitutive_states():
    """Verify ensemble Po stays strictly within [PO_UNCLEAVED_SILENT, PO_GAMMA_ELASTASE_ACTIVE]."""
    f_low, po_low = calculate_enac_open_probability(gamma_cleavage_rate_min_inv=0.0)
    f_high, po_high = calculate_enac_open_probability(gamma_cleavage_rate_min_inv=2.0)

    assert f_low == 0.0
    assert PO_UNCLEAVED_SILENT < po_low < PO_ALPHA_FURIN_CLEAVED
    assert f_high > 0.80
    assert po_low < po_high < PO_GAMMA_ELASTASE_ACTIVE


def test_neutrophil_elastase_depletes_cf_asl_below_critical_7um_threshold():
    """Verify untreated CF airway under neutrophil elastase load collapses below 7.0 um."""
    comp = compare_enac_therapeutic_regimens(neutrophil_elastase_nm=200.0)
    healthy = comp["healthy"]
    cf = comp["untreated_cf"]

    assert healthy.final_asl_height_um >= CRITICAL_CILIARY_HEIGHT_UM
    assert healthy.ciliary_clearance_preserved is True
    assert cf.effective_open_probability > healthy.effective_open_probability
    assert cf.final_asl_height_um < CRITICAL_CILIARY_HEIGHT_UM
    assert cf.ciliary_clearance_preserved is False
    assert cf.mucociliary_velocity_fraction < healthy.mucociliary_velocity_fraction


def test_enac_inhibition_restores_asl_depth_in_non_responsive_cf():
    """Verify ASO/siRNA ENaC knockdown restores ASL height and mucociliary clearance."""
    comp = compare_enac_therapeutic_regimens(neutrophil_elastase_nm=200.0)
    cf = comp["untreated_cf"]
    rescued = comp["enac_aso_rescue"]

    assert rescued.functional_enac_conductance_fraction < cf.functional_enac_conductance_fraction
    assert rescued.sodium_absorption_flux_um_min < cf.sodium_absorption_flux_um_min
    assert rescued.steady_state_asl_height_um > cf.steady_state_asl_height_um
    assert rescued.steady_state_asl_height_um >= CRITICAL_CILIARY_HEIGHT_UM


def test_cftr_chloride_secretion_synergy_with_asl_height():
    """Verify higher CFTR chloride secretory drive increases steady-state ASL height."""
    low_cftr = simulate_enac_asl_regulation(ENaCASLParameters(cftr_chloride_flux_fraction=0.02))
    high_cftr = simulate_enac_asl_regulation(ENaCASLParameters(cftr_chloride_flux_fraction=0.85))

    assert high_cftr.chloride_secretion_flux_um_min > low_cftr.chloride_secretion_flux_um_min
    assert high_cftr.steady_state_asl_height_um > low_cftr.steady_state_asl_height_um


def test_invalid_inputs_raise_value_errors():
    """Verify physical bounds and parameter validation."""
    with pytest.raises(ValueError, match="Neutrophil elastase concentration cannot be negative"):
        calculate_gamma_cleavage_rate(neutrophil_elastase_nm=-5.0)

    with pytest.raises(ValueError, match="Prostasin activity fraction must be in"):
        calculate_gamma_cleavage_rate(neutrophil_elastase_nm=10.0, prostasin_activity_fraction=1.5)

    with pytest.raises(ValueError, match="Endocytic retrieval rate must be positive"):
        calculate_enac_open_probability(gamma_cleavage_rate_min_inv=0.1, k_endocytic_retrieval_min_inv=0.0)

    with pytest.raises(ValueError, match="ENaC inhibitor efficacy must be in"):
        simulate_enac_asl_regulation(ENaCASLParameters(enac_inhibitor_efficacy=1.05))
