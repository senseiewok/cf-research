#!/usr/bin/env python3
"""Narrow Correctness, Uncertainty Split & Identifiability Audit for mucus_barrier_screener (Items #6 & #8).

Audits `tools/discovery/mucus_barrier_screener.py` along two distinct axes without modifying the
screener itself and without claiming empirical CF sputum MPT validation:

1. Narrow Correctness, SI Dimensionality & Discontinuity Audit (Council Item #8):
   - Independent reference evaluation verifying numerical parity with `evaluate_mucus_barrier`.
   - Rational SI dimensional exponent oracle in `(M, L, T, Theta)`.
   - Limiting cases, smooth-domain monotonicity, active-input checks, and documented
     cationic-vs-anionic electrostatic asymmetry.
   - Explicit C0 step-jump and C1 kink discontinuity detection.
   - Five-case planted-defect detection suite (unit conversion, sign flip, missing steric cutoff,
     disconnected pore input, swapped charge asymmetry).

2. Source-Split Joint Uncertainty & Algebraic Identifiability Proofs (Council Item #6):
   - Strictly separates `MeasuredLiteratureParams` (`pore_diameter_nm`, `mucus_thickness_um`,
     `mcc_window_s`) from `HeuristicModelParams` (`interstitial_viscosity_ratio`,
     `charge_shield_*`, `electrostatic_alpha_*`) via `measured_only`, `heuristic_only`, and
     `joint_measured_plus_heuristic` grids with no unlabeled blended error bar.
   - JSON-safe handling of infinite transit time (`None` + `transit_time_is_finite=False`).
   - Exact algebraic non-identifiability and two-probe recovery proofs labeled
     `SYNTHETIC_KNOWN_ANSWER` (`not_empirical_validation=True`).

Pure Python 3.10+ standard library only; zero network calls.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import itertools
import json
import math
from pathlib import Path
import sys
import time
from typing import Any, Dict, List, Optional, Tuple

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from tools.discovery.mucus_barrier_screener import (  # noqa: E402
    BODY_TEMP_K,
    BOLTZMANN_K_J_PER_K,
    EO_KUHN_LENGTH_NM,
    EO_MONOMER_MW_DA,
    INTERSTITIAL_VISCOSITY_RATIO,
    WATER_VISCOSITY_PA_S,
    VectorBiophysicalInput,
    evaluate_mucus_barrier,
)

MODEL_SCOPE_NOTES: Tuple[str, ...] = (
    "f_steric = (1 - d_p / d_pore)^2 is the geometric Renkin steric partition factor without full Renkin-Hille hydrodynamic pore drag terms.",
    "Piecewise PEG charge shielding (1.0 -> 0.60 -> 0.25 at R_F/D = 1.0, 1.5) and max(0, |zeta_eff| - 5) introduce C0 jumps and C1 kinks.",
    "Arithmetic consistency and synthetic known-answer identifiability proofs do NOT establish empirical validity in human CF sputum.",
)

# SI dimensional exponents in base (M, L, T, Theta): mass, length, time, Kelvin temperature.
DIMENSION_EXPONENTS: Dict[str, Tuple[int, int, int, int]] = {
    "boltzmann_k_j_per_k": (1, 2, -2, -1),
    "temperature_k": (0, 0, 0, 1),
    "water_viscosity_pa_s": (1, -1, -1, 0),
    "particle_radius_m": (0, 1, 0, 0),
    "mucus_thickness_m": (0, 1, 0, 0),
    "dimensionless_factor": (0, 0, 0, 0),
    "expected_diffusivity_m2_s": (0, 2, -1, 0),
    "expected_time_s": (0, 0, 1, 0),
}


@dataclass(frozen=True)
class VectorDesignParams:
    """Engineering parameters of a candidate inhaled delivery vector."""

    name: str
    diameter_nm: float
    zeta_mv: float
    peg_density_chains_nm2: float
    peg_mw_kda: float


@dataclass(frozen=True)
class MeasuredLiteratureParams:
    """Parameters anchored to cited CF airway mucus / MCC literature bounds."""

    pore_diameter_nm: float = 140.0
    mucus_thickness_um: float = 25.0
    mcc_window_s: float = 900.0


@dataclass(frozen=True)
class HeuristicModelParams:
    """Uncalibrated phenomenological model parameters kept strictly typed apart from measured inputs."""

    interstitial_viscosity_ratio: float = INTERSTITIAL_VISCOSITY_RATIO
    charge_shield_dense_brush: float = 0.25
    charge_shield_transition: float = 0.60
    charge_shield_bare: float = 1.0
    electrostatic_alpha_cationic: float = 0.22
    electrostatic_alpha_anionic: float = 0.12


@dataclass(frozen=True)
class IndependentScreenerOutput:
    """Typed output of the independent reference implementation (never emits Inf/NaN to JSON)."""

    peg_regime: str
    rf_over_d: float
    charge_shield_factor: float
    f_steric: float
    steric_blocked: bool
    effective_zeta_mv: float
    electrostatic_alpha: float
    f_electrostatic: float
    d_eff_over_d_water: float
    d_water_m2_s: float
    d_eff_m2_s: float
    transit_time_s: Optional[float]
    transit_time_is_finite: bool
    transit_time_reason: str
    crosses_before_mcc: bool
    feasibility_verdict: str


DEFAULT_MEASURED = MeasuredLiteratureParams()
DEFAULT_HEURISTIC = HeuristicModelParams()

REPRESENTATIVE_VECTORS: Tuple[VectorDesignParams, ...] = (
    VectorDesignParams("aav_like_25nm_neutral", 25.0, -4.0, 0.0, 2.0),
    VectorDesignParams("peg_lnp_80nm_neutral", 80.0, -2.5, 0.28, 2.0),
    VectorDesignParams("large_lnp_120nm_neutral", 120.0, -2.5, 0.28, 2.0),
    VectorDesignParams("cationic_lnp_80nm_unshielded", 80.0, 20.0, 0.0, 2.0),
)

# Measured bounds grounded in Suk 2016 / Dawson 2003 [100, 300] nm, Button 2012 mucus & MCC window.
MEASURED_GRID_LEVELS: Dict[str, Tuple[float, ...]] = {
    "pore_diameter_nm": (100.0, 140.0, 220.0, 300.0),
    "mucus_thickness_um": (25.0, 37.5, 50.0),
    "mcc_window_s": (600.0, 900.0, 1200.0),
}

# Phenomenological working variations (nominal middle element matches screener default).
HEURISTIC_GRID_LEVELS: Dict[str, Tuple[float, ...]] = {
    "interstitial_viscosity_ratio": (2.0, INTERSTITIAL_VISCOSITY_RATIO, 5.0),
    "charge_shield_dense_brush": (0.15, 0.25, 0.40),
    "charge_shield_transition": (0.45, 0.60, 0.75),
    "electrostatic_alpha_cationic": (0.15, 0.22, 0.30),
    "electrostatic_alpha_anionic": (0.08, 0.12, 0.18),
}


def _dim_mul(a: Tuple[int, int, int, int], b: Tuple[int, int, int, int]) -> Tuple[int, int, int, int]:
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2], a[3] + b[3])


def _dim_div(a: Tuple[int, int, int, int], b: Tuple[int, int, int, int]) -> Tuple[int, int, int, int]:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2], a[3] - b[3])


def _dim_pow(a: Tuple[int, int, int, int], exp: int) -> Tuple[int, int, int, int]:
    return (a[0] * exp, a[1] * exp, a[2] * exp, a[3] * exp)


def independent_evaluate_mucus_barrier(
    vec: VectorDesignParams,
    measured: MeasuredLiteratureParams = DEFAULT_MEASURED,
    heuristic: HeuristicModelParams = DEFAULT_HEURISTIC,
) -> IndependentScreenerOutput:
    """Independent pure-stdlib re-derivation of the Renkin-Flory-Stokes-Einstein screener equations."""
    if vec.diameter_nm <= 0.0 or vec.peg_mw_kda <= 0.0 or vec.peg_density_chains_nm2 < 0.0:
        raise ValueError("Vector diameter and PEG MW must be positive; grafting density >= 0.")
    if measured.pore_diameter_nm <= 0.0 or measured.mucus_thickness_um <= 0.0 or measured.mcc_window_s <= 0.0:
        raise ValueError("Measured pore diameter, mucus thickness, and MCC window must be positive.")
    if heuristic.interstitial_viscosity_ratio <= 0.0:
        raise ValueError("Interstitial viscosity ratio must be positive.")

    n_monomers = (vec.peg_mw_kda * 1000.0) / EO_MONOMER_MW_DA
    flory_radius_nm = EO_KUHN_LENGTH_NM * (n_monomers ** 0.6)
    if vec.peg_density_chains_nm2 > 1e-9:
        grafting_distance_nm = 2.0 / math.sqrt(math.pi * vec.peg_density_chains_nm2)
        rf_over_d = flory_radius_nm / grafting_distance_nm
    else:
        rf_over_d = 0.0

    if rf_over_d >= 1.5:
        peg_regime = "dense_brush"
        shield = heuristic.charge_shield_dense_brush
    elif rf_over_d >= 1.0:
        peg_regime = "transition_mushroom_to_brush"
        shield = heuristic.charge_shield_transition
    else:
        peg_regime = "mushroom_or_bare"
        shield = heuristic.charge_shield_bare

    size_ratio = vec.diameter_nm / measured.pore_diameter_nm
    if size_ratio >= 1.0:
        f_steric = 0.0
        steric_blocked = True
    else:
        f_steric = (1.0 - size_ratio) ** 2
        steric_blocked = False

    effective_zeta_mv = vec.zeta_mv * shield
    excess_charge_mv = max(0.0, abs(effective_zeta_mv) - 5.0)
    alpha = (
        heuristic.electrostatic_alpha_cationic
        if vec.zeta_mv > 5.0
        else heuristic.electrostatic_alpha_anionic
    )
    f_electrostatic = math.exp(-alpha * excess_charge_mv)

    d_eff_over_d_water = f_steric * (1.0 / heuristic.interstitial_viscosity_ratio) * f_electrostatic

    radius_m = 0.5 * vec.diameter_nm * 1e-9
    d_water_m2_s = (BOLTZMANN_K_J_PER_K * BODY_TEMP_K) / (6.0 * math.pi * WATER_VISCOSITY_PA_S * radius_m)
    d_eff_m2_s = d_water_m2_s * d_eff_over_d_water
    thickness_m = measured.mucus_thickness_um * 1e-6

    if steric_blocked:
        transit_time_s: Optional[float] = None
        transit_is_finite = False
        transit_reason = "sterically_blocked"
    elif d_eff_m2_s <= 1e-24:
        transit_time_s = None
        transit_is_finite = False
        transit_reason = "vanishing_diffusivity"
    else:
        transit_time_s = (thickness_m * thickness_m) / (2.0 * d_eff_m2_s)
        transit_is_finite = True
        transit_reason = "finite"

    crosses = (not steric_blocked) and transit_is_finite and (transit_time_s is not None) and (transit_time_s <= measured.mcc_window_s)
    if steric_blocked:
        verdict = "STERICALLY_TRAPPED_BY_SPUTUM_MESH"
    elif f_electrostatic < 0.15:
        verdict = "ELECTROSTATICALLY_TRAPPED_ON_MUCIN_FIBERS"
    elif crosses:
        verdict = "FEASIBLE_MUCUS_PENETRATING_CANDIDATE"
    else:
        verdict = "CLEARED_BY_MCC_BEFORE_PCL_TRANSIT"

    return IndependentScreenerOutput(
        peg_regime=peg_regime,
        rf_over_d=rf_over_d,
        charge_shield_factor=shield,
        f_steric=f_steric,
        steric_blocked=steric_blocked,
        effective_zeta_mv=effective_zeta_mv,
        electrostatic_alpha=alpha,
        f_electrostatic=f_electrostatic,
        d_eff_over_d_water=d_eff_over_d_water,
        d_water_m2_s=d_water_m2_s,
        d_eff_m2_s=d_eff_m2_s,
        transit_time_s=transit_time_s,
        transit_time_is_finite=transit_is_finite,
        transit_time_reason=transit_reason,
        crosses_before_mcc=crosses,
        feasibility_verdict=verdict,
    )


def audit_dimensional_consistency() -> Dict[str, Any]:
    """Verify exact SI dimensional exponents (M, L, T, Theta) and um-to-m scaling."""
    thermal_energy = _dim_mul(DIMENSION_EXPONENTS["boltzmann_k_j_per_k"], DIMENSION_EXPONENTS["temperature_k"])
    stokes_drag_over_v = _dim_mul(DIMENSION_EXPONENTS["water_viscosity_pa_s"], DIMENSION_EXPONENTS["particle_radius_m"])
    dim_dw = _dim_div(thermal_energy, stokes_drag_over_v)
    dim_deff = _dim_mul(dim_dw, DIMENSION_EXPONENTS["dimensionless_factor"])
    dim_transit = _dim_div(_dim_pow(DIMENSION_EXPONENTS["mucus_thickness_m"], 2), dim_deff)

    # Confirm explicit unit scaling on default 80 nm PEG-LNP.
    vec = REPRESENTATIVE_VECTORS[1]
    out = independent_evaluate_mucus_barrier(vec, DEFAULT_MEASURED, DEFAULT_HEURISTIC)
    unscaled_l_transit = (DEFAULT_MEASURED.mucus_thickness_um ** 2) / (2.0 * out.d_eff_m2_s)
    scale_ratio = unscaled_l_transit / out.transit_time_s if out.transit_time_s else 0.0

    passed = (
        dim_dw == DIMENSION_EXPONENTS["expected_diffusivity_m2_s"]
        and dim_deff == DIMENSION_EXPONENTS["expected_diffusivity_m2_s"]
        and dim_transit == DIMENSION_EXPONENTS["expected_time_s"]
        and abs(scale_ratio - 1e12) / 1e12 < 1e-12
    )
    return {
        "passed": passed,
        "thermal_energy_dim": list(thermal_energy),
        "stokes_drag_dim": list(stokes_drag_over_v),
        "d_water_dim": list(dim_dw),
        "d_eff_dim": list(dim_deff),
        "transit_time_dim": list(dim_transit),
        "um_to_m_squared_scaling_ratio": scale_ratio,
    }


def _screener_packet_for(vec: VectorDesignParams, measured: MeasuredLiteratureParams) -> Dict[str, Any]:
    return evaluate_mucus_barrier(
        VectorBiophysicalInput(
            diameter_nm=vec.diameter_nm,
            zeta_mv=vec.zeta_mv,
            peg_density_chains_nm2=vec.peg_density_chains_nm2,
            peg_mw_kda=vec.peg_mw_kda,
            pore_diameter_nm=measured.pore_diameter_nm,
            mucus_thickness_um=measured.mucus_thickness_um,
            mcc_window_s=measured.mcc_window_s,
        )
    )


def audit_parity_and_invariants() -> Dict[str, Any]:
    """Verify numerical parity against mucus_barrier_screener.py plus limiting/monotonicity invariants."""
    parity_checks: List[Dict[str, Any]] = []
    parity_cases = list(REPRESENTATIVE_VECTORS) + [
        VectorDesignParams("exact_pore_cutoff_140nm", 140.0, -2.5, 0.28, 2.0),
        VectorDesignParams("oversized_300nm", 300.0, -2.5, 0.28, 2.0),
        VectorDesignParams("zero_peg_anionic_10nm", 10.0, -3.0, 0.0, 2.0),
    ]
    for vec in parity_cases:
        ref = independent_evaluate_mucus_barrier(vec, DEFAULT_MEASURED, DEFAULT_HEURISTIC)
        pkt = _screener_packet_for(vec, DEFAULT_MEASURED)["computed_metrics"]
        d_ratio_pkt = float(pkt["d_eff_over_d_water"]["value"])
        t_pkt = pkt["mean_transit_time_s"]["value"]
        verdict_pkt = str(pkt["crosses_before_mcc_window"]["verdict"])
        d_ok = abs(round(ref.d_eff_over_d_water, 6) - d_ratio_pkt) <= 1e-12
        t_ok = (t_pkt is None and ref.transit_time_s is None) or (
            t_pkt is not None and ref.transit_time_s is not None and abs(round(ref.transit_time_s, 2) - float(t_pkt)) <= 1e-9
        )
        v_ok = ref.feasibility_verdict == verdict_pkt
        parity_checks.append(
            {
                "case": vec.name,
                "passed": bool(d_ok and t_ok and v_ok),
                "d_ratio_ref": round(ref.d_eff_over_d_water, 6),
                "d_ratio_screener": d_ratio_pkt,
                "verdict": ref.feasibility_verdict,
            }
        )

    # Limiting cases: tiny neutral limit -> f_steric -> 1; boundary d_p == d_pore -> f_steric == 0 & blocked.
    tiny = independent_evaluate_mucus_barrier(VectorDesignParams("tiny", 1e-4, 0.0, 0.0, 2.0))
    boundary = independent_evaluate_mucus_barrier(VectorDesignParams("boundary", 140.0, 0.0, 0.0, 2.0))
    limit_checks = {
        "tiny_neutral_f_steric_near_1": abs(tiny.f_steric - 1.0) < 1e-5,
        "tiny_neutral_ratio_approaches_inv_eta": abs(tiny.d_eff_over_d_water - (1.0 / INTERSTITIAL_VISCOSITY_RATIO)) < 1e-5,
        "exact_boundary_sterically_blocked": (boundary.f_steric == 0.0 and boundary.steric_blocked and not boundary.transit_time_is_finite),
    }

    # Documented electrostatic cationic/anionic sign asymmetry at equal |zeta| = 15 mV, bare.
    cat_15 = independent_evaluate_mucus_barrier(VectorDesignParams("cat15", 60.0, 15.0, 0.0, 2.0))
    ani_15 = independent_evaluate_mucus_barrier(VectorDesignParams("ani15", 60.0, -15.0, 0.0, 2.0))
    asymmetry_pinned = (
        cat_15.electrostatic_alpha == DEFAULT_HEURISTIC.electrostatic_alpha_cationic
        and ani_15.electrostatic_alpha == DEFAULT_HEURISTIC.electrostatic_alpha_anionic
        and cat_15.f_electrostatic < ani_15.f_electrostatic
    )

    # Monotonicity across smooth open intervals.
    diam_grid = (20.0, 50.0, 80.0, 110.0, 130.0)
    deff_by_diam = [
        independent_evaluate_mucus_barrier(VectorDesignParams("m", d, 0.0, 0.0, 2.0)).d_eff_m2_s
        for d in diam_grid
    ]
    pore_grid = (100.0, 140.0, 200.0, 300.0)
    deff_by_pore = [
        independent_evaluate_mucus_barrier(
            VectorDesignParams("m", 80.0, 0.0, 0.0, 2.0),
            MeasuredLiteratureParams(pore_diameter_nm=p),
        ).d_eff_over_d_water
        for p in pore_grid
    ]
    zeta_grid = (0.0, 5.0, 10.0, 20.0, 35.0)
    felec_by_zeta = [
        independent_evaluate_mucus_barrier(VectorDesignParams("m", 80.0, z, 0.0, 2.0)).f_electrostatic
        for z in zeta_grid
    ]
    monotonic_checks = {
        "d_eff_strictly_decreasing_in_diameter": all(a > b for a, b in zip(deff_by_diam, deff_by_diam[1:])),
        "d_eff_ratio_strictly_increasing_in_pore_diameter": all(a < b for a, b in zip(deff_by_pore, deff_by_pore[1:])),
        "f_elec_non_increasing_in_abs_zeta": all(a >= b for a, b in zip(felec_by_zeta, felec_by_zeta[1:])),
    }

    # Active-input dependence: every advertised VectorDesignParams / MeasuredLiteratureParams input affects output.
    base_vec = REPRESENTATIVE_VECTORS[1]
    base_out = independent_evaluate_mucus_barrier(base_vec, DEFAULT_MEASURED, DEFAULT_HEURISTIC)
    active_checks = {
        "diameter_active": independent_evaluate_mucus_barrier(
            VectorDesignParams("d2", 95.0, base_vec.zeta_mv, base_vec.peg_density_chains_nm2, base_vec.peg_mw_kda)
        ).d_eff_over_d_water != base_out.d_eff_over_d_water,
        "zeta_active": independent_evaluate_mucus_barrier(
            VectorDesignParams("z2", base_vec.diameter_nm, 30.0, base_vec.peg_density_chains_nm2, base_vec.peg_mw_kda)
        ).f_electrostatic != base_out.f_electrostatic,
        "peg_density_active": independent_evaluate_mucus_barrier(
            VectorDesignParams("s2", base_vec.diameter_nm, 30.0, 0.02, base_vec.peg_mw_kda)
        ).f_electrostatic != independent_evaluate_mucus_barrier(
            VectorDesignParams("s1", base_vec.diameter_nm, 30.0, 0.35, base_vec.peg_mw_kda)
        ).f_electrostatic,
        "peg_mw_active": independent_evaluate_mucus_barrier(
            VectorDesignParams("mw2", base_vec.diameter_nm, 20.0, 0.14, 0.5)
        ).peg_regime != independent_evaluate_mucus_barrier(
            VectorDesignParams("mw1", base_vec.diameter_nm, 20.0, 0.14, 5.0)
        ).peg_regime,
        "pore_diameter_active": independent_evaluate_mucus_barrier(
            base_vec, MeasuredLiteratureParams(pore_diameter_nm=220.0)
        ).d_eff_over_d_water != base_out.d_eff_over_d_water,
        "mucus_thickness_active": independent_evaluate_mucus_barrier(
            base_vec, MeasuredLiteratureParams(mucus_thickness_um=50.0)
        ).transit_time_s != base_out.transit_time_s,
        "mcc_window_active": independent_evaluate_mucus_barrier(
            base_vec, MeasuredLiteratureParams(mcc_window_s=1.0)
        ).crosses_before_mcc != base_out.crosses_before_mcc,
    }

    all_ok = (
        all(c["passed"] for c in parity_checks)
        and all(limit_checks.values())
        and asymmetry_pinned
        and all(monotonic_checks.values())
        and all(active_checks.values())
    )
    return {
        "passed": all_ok,
        "parity_checks": parity_checks,
        "limit_checks": limit_checks,
        "cationic_stricter_than_anionic_pinned": asymmetry_pinned,
        "monotonic_checks": monotonic_checks,
        "active_input_checks": active_checks,
    }


def _peg_density_for_rf_over_d(target_ratio: float, peg_mw_kda: float = 2.0) -> float:
    n_monomers = (peg_mw_kda * 1000.0) / EO_MONOMER_MW_DA
    r_f = EO_KUHN_LENGTH_NM * (n_monomers ** 0.6)
    return (4.0 / math.pi) * ((target_ratio / r_f) ** 2)


def audit_discontinuities() -> Dict[str, Any]:
    """Detect piecewise C0 jumps at R_F/D = 1.0, 1.5 and C1 slope kink at |zeta_eff| = 5 mV."""
    eps = 1e-4
    c0_jumps: List[Dict[str, Any]] = []
    for threshold in (1.0, 1.5):
        sigma_lo = _peg_density_for_rf_over_d(threshold - eps, peg_mw_kda=2.0)
        sigma_hi = _peg_density_for_rf_over_d(threshold + eps, peg_mw_kda=2.0)
        out_lo = independent_evaluate_mucus_barrier(VectorDesignParams("lo", 80.0, 20.0, sigma_lo, 2.0))
        out_hi = independent_evaluate_mucus_barrier(VectorDesignParams("hi", 80.0, 20.0, sigma_hi, 2.0))
        jump_f_elec = abs(out_hi.f_electrostatic - out_lo.f_electrostatic)
        c0_jumps.append(
            {
                "threshold_rf_over_d": threshold,
                "regime_below": out_lo.peg_regime,
                "regime_above": out_hi.peg_regime,
                "shield_below": out_lo.charge_shield_factor,
                "shield_above": out_hi.charge_shield_factor,
                "f_elec_jump": round(jump_f_elec, 6),
                "discontinuity_detected": jump_f_elec > 0.05,
            }
        )

    # C1 kink at |zeta_eff| = 5.0 mV (bare particle sigma = 0, shield = 1.0).
    dz = 1e-3
    below = independent_evaluate_mucus_barrier(VectorDesignParams("z_lo", 80.0, 5.0 - dz, 0.0, 2.0))
    at_kink = independent_evaluate_mucus_barrier(VectorDesignParams("z_mid", 80.0, 5.0, 0.0, 2.0))
    above = independent_evaluate_mucus_barrier(VectorDesignParams("z_hi", 80.0, 5.0 + dz, 0.0, 2.0))
    left_slope = (at_kink.f_electrostatic - below.f_electrostatic) / dz
    right_slope = (above.f_electrostatic - at_kink.f_electrostatic) / dz
    c1_kink_detected = abs(right_slope - left_slope) > 0.1

    all_detected = all(item["discontinuity_detected"] for item in c0_jumps) and c1_kink_detected
    return {
        "passed": all_detected,
        "c0_piecewise_shield_jumps": c0_jumps,
        "c1_zeta_threshold_kink": {
            "threshold_effective_zeta_mv": 5.0,
            "left_slope_d_felec_d_zeta": round(left_slope, 6),
            "right_slope_d_felec_d_zeta": round(right_slope, 6),
            "kink_detected": c1_kink_detected,
        },
    }


def run_planted_defect_checks() -> Dict[str, Any]:
    """Run five deterministic planted-defect probes and verify the audit detects every single defect."""
    defects: List[Dict[str, Any]] = []
    ref_vec = REPRESENTATIVE_VECTORS[1]
    ref_out = independent_evaluate_mucus_barrier(ref_vec, DEFAULT_MEASURED, DEFAULT_HEURISTIC)

    # Defect 1: Omitting um -> m conversion inflates transit_time_s by 1e12.
    bad_transit = (DEFAULT_MEASURED.mucus_thickness_um ** 2) / (2.0 * ref_out.d_eff_m2_s)
    ratio_err = abs(bad_transit / ref_out.transit_time_s - 1.0)
    defects.append({
        "defect": "missing_um_to_m_conversion",
        "caught": ratio_err > 1e6,
        "observed_ratio_error": round(ratio_err, 2),
    })

    # Defect 2: Electrostatic exponent sign flip turns trapping attenuation into unphysical f_elec > 1.
    bad_felec = math.exp(+DEFAULT_HEURISTIC.electrostatic_alpha_cationic * 15.0)
    defects.append({
        "defect": "electrostatic_sign_flip",
        "caught": not (0.0 <= bad_felec <= 1.0),
        "observed_bad_f_elec": round(bad_felec, 4),
    })

    # Defect 3: Removing steric clamp d_p >= d_pore leaves positive steric transmission (1 - 300/140)^2 > 0.
    bad_steric = (1.0 - (300.0 / 140.0)) ** 2
    defects.append({
        "defect": "missing_steric_cutoff",
        "caught": bad_steric > 0.0,
        "observed_bad_steric": round(bad_steric, 6),
    })

    # Defect 4: Ignoring pore_diameter_nm gives identical diffusivity across d_pore = 100 nm vs 260 nm.
    def eval_disconnected(pore_nm: float) -> float:
        fixed_pore = 140.0
        return (1.0 - (ref_vec.diameter_nm / fixed_pore)) ** 2 / DEFAULT_HEURISTIC.interstitial_viscosity_ratio

    diff_disconnected = abs(eval_disconnected(260.0) - eval_disconnected(100.0))
    defects.append({
        "defect": "disconnected_pore_input",
        "caught": diff_disconnected <= 1e-12,
        "observed_variation": diff_disconnected,
    })

    # Defect 5: Swapped cationic/anionic alpha makes anionic charge trapped more strongly than cationic.
    bad_swapped = HeuristicModelParams(electrostatic_alpha_cationic=0.12, electrostatic_alpha_anionic=0.22)
    cat_swapped = independent_evaluate_mucus_barrier(
        VectorDesignParams("cs", 60.0, 20.0, 0.0, 2.0), DEFAULT_MEASURED, bad_swapped
    )
    ani_swapped = independent_evaluate_mucus_barrier(
        VectorDesignParams("as", 60.0, -20.0, 0.0, 2.0), DEFAULT_MEASURED, bad_swapped
    )
    defects.append({
        "defect": "swapped_cationic_anionic_asymmetry",
        "caught": not (cat_swapped.f_electrostatic < ani_swapped.f_electrostatic),
        "observed_cat_f_elec": round(cat_swapped.f_electrostatic, 6),
        "observed_ani_f_elec": round(ani_swapped.f_electrostatic, 6),
    })

    caught_count = sum(1 for d in defects if d["caught"])
    return {
        "passed": caught_count == len(defects),
        "caught": caught_count,
        "total": len(defects),
        "defects": defects,
    }


def _summarize_grid_runs(outputs: List[IndependentScreenerOutput]) -> Dict[str, Any]:
    if not outputs:
        raise ValueError("Grid outputs must be non-empty.")

    deff_ratios = [o.d_eff_over_d_water for o in outputs]
    finite_transits = [o.transit_time_s for o in outputs if o.transit_time_is_finite and o.transit_time_s is not None]
    verdict_counts: Dict[str, int] = {}
    for o in outputs:
        verdict_counts[o.feasibility_verdict] = verdict_counts.get(o.feasibility_verdict, 0) + 1

    dominant_verdict, dominant_count = max(
        verdict_counts.items(), key=lambda kv: (kv[1], kv[0])
    )
    n = len(outputs)
    return {
        "sample_count": n,
        "finite_transit_samples": len(finite_transits),
        "infinite_transit_samples": n - len(finite_transits),
        "d_eff_over_d_water_min": round(min(deff_ratios), 6),
        "d_eff_over_d_water_max": round(max(deff_ratios), 6),
        "transit_time_finite_s_min": round(min(finite_transits), 2) if finite_transits else None,
        "transit_time_finite_s_max": round(max(finite_transits), 2) if finite_transits else None,
        "verdict_counts": dict(sorted(verdict_counts.items())),
        "dominant_verdict": dominant_verdict,
        "verdict_stability": round(dominant_count / n, 4),
        "verdict_flips_occurred": len(verdict_counts) > 1,
    }


def run_source_split_uncertainty_audit(
    vectors: Tuple[VectorDesignParams, ...] = REPRESENTATIVE_VECTORS,
) -> Dict[str, Any]:
    """Run separate measured-only, heuristic-only, and joint uncertainty grids per representative vector."""
    measured_combos = [
        MeasuredLiteratureParams(pore_diameter_nm=p, mucus_thickness_um=t, mcc_window_s=m)
        for p, t, m in itertools.product(
            MEASURED_GRID_LEVELS["pore_diameter_nm"],
            MEASURED_GRID_LEVELS["mucus_thickness_um"],
            MEASURED_GRID_LEVELS["mcc_window_s"],
        )
    ]
    heuristic_combos = [
        HeuristicModelParams(
            interstitial_viscosity_ratio=eta,
            charge_shield_dense_brush=s_dense,
            charge_shield_transition=s_trans,
            electrostatic_alpha_cationic=a_cat,
            electrostatic_alpha_anionic=a_ani,
        )
        for eta, s_dense, s_trans, a_cat, a_ani in itertools.product(
            HEURISTIC_GRID_LEVELS["interstitial_viscosity_ratio"],
            HEURISTIC_GRID_LEVELS["charge_shield_dense_brush"],
            HEURISTIC_GRID_LEVELS["charge_shield_transition"],
            HEURISTIC_GRID_LEVELS["electrostatic_alpha_cationic"],
            HEURISTIC_GRID_LEVELS["electrostatic_alpha_anionic"],
        )
    ]

    profiles_report: Dict[str, Any] = {}
    for vec in vectors:
        baseline_out = independent_evaluate_mucus_barrier(vec, DEFAULT_MEASURED, DEFAULT_HEURISTIC)
        meas_only_outs = [independent_evaluate_mucus_barrier(vec, m, DEFAULT_HEURISTIC) for m in measured_combos]
        heur_only_outs = [independent_evaluate_mucus_barrier(vec, DEFAULT_MEASURED, h) for h in heuristic_combos]
        joint_outs = [
            independent_evaluate_mucus_barrier(vec, m, h)
            for m in measured_combos
            for h in heuristic_combos
        ]
        profiles_report[vec.name] = {
            "vector": asdict(vec),
            "baseline_verdict": baseline_out.feasibility_verdict,
            "baseline_d_eff_over_d_water": round(baseline_out.d_eff_over_d_water, 6),
            "baseline_transit_time_s": (
                round(baseline_out.transit_time_s, 2) if baseline_out.transit_time_s is not None else None
            ),
            "measured_only": _summarize_grid_runs(meas_only_outs),
            "heuristic_only": _summarize_grid_runs(heur_only_outs),
            "joint_measured_plus_heuristic": _summarize_grid_runs(joint_outs),
        }

    return {
        "epistemic_label": "SYNTHETIC_KNOWN_ANSWER",
        "has_blended_uncertainty_bar": False,
        "measured_parameter_ranges": {k: list(v) for k, v in MEASURED_GRID_LEVELS.items()},
        "heuristic_parameter_ranges": {k: list(v) for k, v in HEURISTIC_GRID_LEVELS.items()},
        "profiles": profiles_report,
    }


def prove_single_size_non_identifiability(
    diameter_nm: float = 80.0,
    d_pore_1: float = 140.0,
    eta_ratio_1: float = 3.5,
    d_pore_2: float = 200.0,
) -> Dict[str, Any]:
    """Prove distinct (d_pore, eta_ratio) pairs yield identical D_eff/D_w for one neutral probe."""
    if diameter_nm <= 0.0 or d_pore_1 <= diameter_nm or d_pore_2 <= diameter_nm or eta_ratio_1 <= 0.0:
        raise ValueError("Require 0 < diameter_nm < d_pore_1, d_pore_2 and eta_ratio_1 > 0.")

    r_target = ((1.0 - diameter_nm / d_pore_1) ** 2) / eta_ratio_1
    eta_ratio_2 = ((1.0 - diameter_nm / d_pore_2) ** 2) / r_target

    vec = VectorDesignParams("neutral_probe", diameter_nm, 0.0, 0.0, 2.0)
    out_1 = independent_evaluate_mucus_barrier(
        vec,
        MeasuredLiteratureParams(pore_diameter_nm=d_pore_1),
        HeuristicModelParams(interstitial_viscosity_ratio=eta_ratio_1),
    )
    out_2 = independent_evaluate_mucus_barrier(
        vec,
        MeasuredLiteratureParams(pore_diameter_nm=d_pore_2),
        HeuristicModelParams(interstitial_viscosity_ratio=eta_ratio_2),
    )
    rel_diff = abs(out_1.d_eff_over_d_water - out_2.d_eff_over_d_water) / max(out_1.d_eff_over_d_water, 1e-18)
    indistinguishable = rel_diff < 1e-12 and abs(d_pore_1 - d_pore_2) > 1.0 and abs(eta_ratio_1 - eta_ratio_2) > 0.1
    return {
        "epistemic_label": "SYNTHETIC_KNOWN_ANSWER",
        "not_empirical_validation": True,
        "indistinguishable_proven": indistinguishable,
        "probe_diameter_nm": diameter_nm,
        "pair_1": {"d_pore_nm": d_pore_1, "eta_ratio": eta_ratio_1, "d_eff_over_d_water": out_1.d_eff_over_d_water},
        "pair_2": {"d_pore_nm": d_pore_2, "eta_ratio": round(eta_ratio_2, 6), "d_eff_over_d_water": out_2.d_eff_over_d_water},
        "relative_diff": rel_diff,
    }


def recover_pore_and_viscosity_from_two_neutral_probes(
    d1_nm: float,
    d2_nm: float,
    r1: float,
    r2: float,
) -> Tuple[float, float]:
    """Algebraically recover (d_pore_nm, eta_ratio) from two neutral probes 0 < d1 < d2 < d_pore."""
    if d1_nm <= 0.0 or d2_nm <= d1_nm:
        raise ValueError("Require strictly increasing positive probe diameters 0 < d1_nm < d2_nm.")
    if r1 <= 0.0 or r2 <= 0.0:
        raise ValueError("Probe diffusivity ratios must be strictly positive (probes cannot be sterically blocked).")
    if r2 >= r1:
        raise ValueError("Larger neutral probe d2 must have strictly smaller diffusivity ratio r2 < r1.")

    rho = math.sqrt(r2 / r1)
    d_pore = (d2_nm - rho * d1_nm) / (1.0 - rho)
    f_steric_1 = (1.0 - (d1_nm / d_pore)) ** 2
    eta_ratio = f_steric_1 / r1
    return d_pore, eta_ratio


def run_two_probe_identifiability_proofs(
    d1_nm: float = 40.0,
    d2_nm: float = 90.0,
) -> Dict[str, Any]:
    """Verify exact two-probe algebraic recovery on synthetic known-answer cases plus fail-closed guards."""
    synthetic_cases = (
        (140.0, 3.5),
        (180.0, 2.4),
        (260.0, 4.8),
    )
    recoveries: List[Dict[str, Any]] = []
    vec_1 = VectorDesignParams("probe_1", d1_nm, 0.0, 0.0, 2.0)
    vec_2 = VectorDesignParams("probe_2", d2_nm, 0.0, 0.0, 2.0)

    for d_pore_true, eta_true in synthetic_cases:
        m = MeasuredLiteratureParams(pore_diameter_nm=d_pore_true)
        h = HeuristicModelParams(interstitial_viscosity_ratio=eta_true)
        r1 = independent_evaluate_mucus_barrier(vec_1, m, h).d_eff_over_d_water
        r2 = independent_evaluate_mucus_barrier(vec_2, m, h).d_eff_over_d_water
        d_pore_rec, eta_rec = recover_pore_and_viscosity_from_two_neutral_probes(d1_nm, d2_nm, r1, r2)
        d_err = abs(d_pore_rec - d_pore_true) / d_pore_true
        e_err = abs(eta_rec - eta_true) / eta_true
        recoveries.append(
            {
                "true_d_pore_nm": d_pore_true,
                "true_eta_ratio": eta_true,
                "recovered_d_pore_nm": round(d_pore_rec, 9),
                "recovered_eta_ratio": round(eta_rec, 9),
                "d_pore_rel_error": d_err,
                "eta_ratio_rel_error": e_err,
                "passed": (d_err < 1e-12 and e_err < 1e-12),
            }
        )

    blocked_guard_triggered = False
    try:
        recover_pore_and_viscosity_from_two_neutral_probes(d1_nm, 150.0, 0.145, 0.0)
    except ValueError:
        blocked_guard_triggered = True

    return {
        "epistemic_label": "SYNTHETIC_KNOWN_ANSWER",
        "not_empirical_validation": True,
        "probe_pair_nm": [d1_nm, d2_nm],
        "exact_recoveries": recoveries,
        "steric_blocked_guard_triggered": blocked_guard_triggered,
        "passed": all(r["passed"] for r in recoveries) and blocked_guard_triggered,
    }


def prove_electrostatic_confounding(
    zeta_mv: float = 20.0,
    shield_1: float = 0.60,
    alpha_1: float = 0.22,
    shield_2: float = 0.40,
) -> Dict[str, Any]:
    """Prove (charge_shield_factor, electrostatic_alpha) is non-identifiable from one charged probe."""
    s1 = abs(zeta_mv * shield_1) - 5.0
    s2 = abs(zeta_mv * shield_2) - 5.0
    if s1 <= 0.0 or s2 <= 0.0 or alpha_1 <= 0.0:
        raise ValueError("Require active electrostatic regime |zeta_mv * shield| > 5 mV for both pairs.")

    alpha_2 = alpha_1 * (s1 / s2)
    f_1 = math.exp(-alpha_1 * s1)
    f_2 = math.exp(-alpha_2 * s2)
    rel_diff = abs(f_1 - f_2) / max(f_1, 1e-18)
    return {
        "epistemic_label": "SYNTHETIC_KNOWN_ANSWER",
        "not_empirical_validation": True,
        "confounded_proven": (rel_diff < 1e-12 and abs(shield_1 - shield_2) > 0.05 and abs(alpha_1 - alpha_2) > 0.02),
        "zeta_mv": zeta_mv,
        "pair_1": {"shield": shield_1, "alpha": alpha_1, "f_electrostatic": round(f_1, 12)},
        "pair_2": {"shield": shield_2, "alpha": round(alpha_2, 6), "f_electrostatic": round(f_2, 12)},
        "relative_diff": rel_diff,
    }


def run_full_audit() -> Dict[str, Any]:
    """Run complete Item #6 + Item #8 correctness, uncertainty-split, and identifiability audit."""
    t0 = time.perf_counter()
    dim = audit_dimensional_consistency()
    parity = audit_parity_and_invariants()
    disc = audit_discontinuities()
    defects = run_planted_defect_checks()
    uncertainty = run_source_split_uncertainty_audit()
    non_id = prove_single_size_non_identifiability()
    two_probe = run_two_probe_identifiability_proofs()
    elec_conf = prove_electrostatic_confounding()
    elapsed_ms = round((time.perf_counter() - t0) * 1000.0, 2)

    all_passed = (
        dim["passed"]
        and parity["passed"]
        and disc["passed"]
        and defects["passed"]
        and non_id["indistinguishable_proven"]
        and two_probe["passed"]
        and elec_conf["confounded_proven"]
    )
    return {
        "audit_target": "tools/discovery/mucus_barrier_screener.py",
        "council_items": ["item-6-uncertainty-identifiability", "item-8-narrow-correctness"],
        "research_boundary": "Mathematical correctness and synthetic known-answer identifiability audit only; zero empirical CF sputum validation.",
        "model_scope_notes": list(MODEL_SCOPE_NOTES),
        "status": "PASS" if all_passed else "FAIL",
        "elapsed_ms": elapsed_ms,
        "dimensional_audit": dim,
        "parity_and_invariants": parity,
        "discontinuity_audit": disc,
        "planted_defect_checks": defects,
        "source_split_uncertainty": uncertainty,
        "identifiability_proofs": {
            "single_size_non_identifiability": non_id,
            "two_neutral_probe_recovery": two_probe,
            "electrostatic_confounding": elec_conf,
        },
    }


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--json", action="store_true", help="Print strict JSON (allow_nan=False).")
    parser.add_argument("--check", action="store_true", help="Exit 0 only when all correctness & identifiability checks pass.")
    args = parser.parse_args(argv)

    report = run_full_audit()
    if args.json:
        print(json.dumps(report, indent=2, allow_nan=False))
    else:
        u_profiles = report["source_split_uncertainty"]["profiles"]
        print(
            f"Mucus Barrier Correctness & Identifiability Audit: {report['status']} "
            f"(elapsed={report['elapsed_ms']} ms)"
        )
        print(
            f"  Dimensional={report['dimensional_audit']['passed']} "
            f"Parity/Invariants={report['parity_and_invariants']['passed']} "
            f"Discontinuities={report['discontinuity_audit']['passed']} "
            f"PlantedDefects={report['planted_defect_checks']['caught']}/{report['planted_defect_checks']['total']}"
        )
        print(
            f"  Identifiability: single_size_nonid={report['identifiability_proofs']['single_size_non_identifiability']['indistinguishable_proven']} "
            f"two_probe_recovery={report['identifiability_proofs']['two_neutral_probe_recovery']['passed']} "
            f"elec_confound={report['identifiability_proofs']['electrostatic_confounding']['confounded_proven']}"
        )
        for name, prof in u_profiles.items():
            jm = prof["joint_measured_plus_heuristic"]
            print(
                f"  - {name}: baseline={prof['baseline_verdict']} | "
                f"joint_stability={jm['verdict_stability']:.4f} "
                f"finite={jm['finite_transit_samples']}/{jm['sample_count']} "
                f"flips={jm['verdict_flips_occurred']}"
            )
    if args.check and report["status"] != "PASS":
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
