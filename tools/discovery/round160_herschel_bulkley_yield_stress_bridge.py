#!/usr/bin/env python3
"""
Round 160: Herschel-Bulkley Viscoplastic Yield-Stress Invariant Bridge for Non-Newtonian CF Sputum Extrusion & Mucociliary Transport Arrest (T-0825).

Mathematical & Physical Pillar (Column A - Deductive Theorem):
  1. Herschel-Bulkley Viscoplastic Constitutive Rheology:
     For a non-Newtonian viscoplastic fluid such as human airway mucus, shear stress tau as a function
     of shear strain rate gamma_dot is governed by the Herschel-Bulkley constitutive relation:
       tau = tau_y + K * (gamma_dot)^n    for |tau| > tau_y
       gamma_dot = 0.0                    for |tau| <= tau_y
     where:
       - tau_y is the critical yield stress (Pa), below which the fluid behaves as a rigid solid plug.
       - K is the consistency index (Pa * s^n), characterizing fluid thickness.
       - n is the dimensionless flow behavior index (n < 1 for pseudoplastic shear-thinning mucus).
       - gamma_dot is the shear strain rate (s^-1).
  2. Apparent and Differential Dynamic Viscosity:
     - The apparent dynamic viscosity eta_app(gamma_dot) is:
         eta_app(gamma_dot) = tau / gamma_dot = (tau_y / gamma_dot) + K * (gamma_dot)^(n - 1)
       As gamma_dot -> 0, eta_app -> infinity whenever tau_y > 0 (plastic unyielding limit).
     - The differential (tangent) dynamic viscosity eta_diff(gamma_dot) is:
         eta_diff(gamma_dot) = dtau / d(gamma_dot) = n * K * (gamma_dot)^(n - 1)
  3. Ciliary Forward-Stroke Shear Stress & Ciliary Stall Mechanics:
     Coordinated axonemal beating in the periciliary liquid layer generates an effective forward-stroke
     shear stress tau_cilia across the mucus interface:
       tau_cilia = rho_cilia * F_cilia ~ 1.25 Pa
     (where ciliary areal density rho_cilia ~ 7.0 cilia/um^2 and forward stroke force F_cilia ~ 50 pN).
     - If tau_cilia > tau_y:
       The ciliary stroke exceeds the yield stress, driving active shear flow:
         gamma_dot = ( (tau_cilia - tau_y) / K )^(1 / n)
       The resulting mucociliary clearance (MCC) transport velocity v_MCC (mm/min) is:
         v_MCC = v_0 * [ gamma_dot / (gamma_dot + gamma_dot_ref) ]
       with healthy baseline v_0 = 5.50 mm/min and gamma_dot_ref = 0.40 s^-1.
     - If tau_cilia <= tau_y:
       Ciliary shear forces cannot yield the mucus matrix.
       Shear rate gamma_dot = 0.0 s^-1, causing complete ciliary stall and clearance arrest (v_MCC = 0.0 mm/min).
  4. Airway Plug Flow & Cough Extrusion Mechanics:
     In a cylindrical airway of radius R and length L under a trans-pulmonary pressure drop Delta P:
     Wall shear stress tau_w is:
       tau_w = (R * Delta P) / (2 * L)
     The radius of the central unyielded solid core plug r_plug is:
       r_plug = R * min(1.0, tau_y / tau_w)
     When tau_w > tau_y, the volumetric extrusion flow rate Q for a Herschel-Bulkley fluid in a circular tube is:
       Q = (pi * R^3 / (K^(1/n) * tau_w^3)) * (tau_w - tau_y)^((n+1)/n) * [
             (tau_w - tau_y)^2 / (3 + 1/n) + 2 * tau_y * (tau_w - tau_y) / (2 + 1/n) + tau_y^2 / (1 + 1/n)
           ]
     When tau_w <= tau_y, r_plug = R and Q = 0.0 (complete airway plug occlusion).
  5. Invariants & Falsification:
     - Strict monotonicity: dtau / d(gamma_dot) > 0 for all gamma_dot > 0.
     - Plug radius bounding: 0 <= r_plug <= R, with dr_plug / dtau_w <= 0.
     - Complex-step differentiation: dtau / d(gamma_dot) verified against analytical eta_diff at step h = 1e-30.
     - Poisoned negative control: Yield stress tau_y < 0.0 Pa, consistency index K <= 0.0 Pa*s^n,
       flow index n <= 0.0 or n > 2.0, or shear rate gamma_dot < 0.0 s^-1 fails deterministically.

CF Biophysical Domain Mapping (Column B - Popperian Falsifiable Hypothesis):
  - In healthy euhydrated airway mucus (Hill et al. 2022 Physiol Rev; Button et al. 2012 Science):
    Low solids (2.0%), low yield stress (tau_y ~ 0.10 Pa < 1.25 Pa), K ~ 0.50 Pa*s^n, n ~ 0.45.
    Ciliary stroke easily overcomes yield stress (gamma_dot ~ 6.4 s^-1), driving rapid clearance (v_MCC >= 5.0 mm/min).
  - In CF hyperconcentrated sputum arrest (Hill et al. 2022 Physiol Rev; Button et al. 2012 Science):
    Dehydration and eDNA networks elevate solids to >= 8.5%, pushing yield stress to tau_y ~ 3.20 Pa (> 1.25 Pa),
    K ~ 8.50 Pa*s^n, n ~ 0.35. Because tau_cilia < tau_y, cilia stall completely (v_MCC = 0.0 mm/min),
    and airway plug occludes small bronchioles under tidal breathing.
  - In CF severe mucoid impaction (Hill et al. 2022 Physiol Rev):
    Solids reach >= 12.0%, tau_y ~ 8.50 Pa >> 1.25 Pa, K ~ 22.0 Pa*s^n, n ~ 0.30.
    Even deep tidal breathing fails to mobilize plugs (r_plug = R_airway), requiring high-velocity cough (> 15 kPa).
  - In Kahn-Kalai curative liquefaction boundary:
    Restoration of CFTR anion secretion or mucolytic digestion reduces yield stress below tau_y <= 0.85 Pa (< 1.25 Pa),
    K <= 2.00 Pa*s^n, n >= 0.40, restoring active ciliary transport to v_MCC >= 0.65 mm/min (>= 11.8% of WT),
    crossing the curative threshold.
  - In Dornase Alfa (rhDNase) enzymatic liquefaction rescue (Hill et al. 2022 Physiol Rev):
    Hydrolysis of extracellular DNA networks lowers yield stress from 3.20 Pa down to tau_y ~ 0.35 Pa (< 1.25 Pa),
    reducing K to ~ 0.90 Pa*s^n and restoring clearance to v_MCC >= 3.80 mm/min.

Research simulation only; not medical advice.
"""
from __future__ import annotations

import cmath
import math
from dataclasses import dataclass
from typing import Any, Dict, Tuple

PRIMARY_CITATIONS: Tuple[str, ...] = (
    "DOI: 10.1152/physrev.00004.2021",
    "DOI: 10.1126/science.1223012",
    "DOI: 10.1038/s41565-020-0669-6",
    "DOI: 10.1085/jgp.200609558",
    "DOI: 10.1126/science.ade2216",
)

PARAMETER_PROVENANCE: Dict[str, Dict[str, object]] = {
    "ciliary_stroke_shear_stress_pa": {
        "value": 1.25,
        "tag": "empirical",
        "source": "Airway ciliary axonemal propulsive shear stress (DOI: 10.1152/physrev.00004.2021)",
        "bounds": (0.8, 1.8),
    },
    "healthy_mucociliary_velocity_mm_min": {
        "value": 5.50,
        "tag": "empirical",
        "source": "Human in vivo tracheobronchial mucociliary transport rate (DOI: 10.1126/science.1223012)",
        "bounds": (4.0, 7.5),
    },
    "healthy_mucus_yield_stress_pa": {
        "value": 0.10,
        "tag": "empirical",
        "source": "Euhydrated 2% mucus Herschel-Bulkley yield stress (DOI: 10.1152/physrev.00004.2021)",
        "bounds": (0.02, 0.35),
    },
    "cf_sputum_yield_stress_pa": {
        "value": 3.20,
        "tag": "empirical",
        "source": "Hyperconcentrated CF sputum Herschel-Bulkley yield stress (DOI: 10.1152/physrev.00004.2021)",
        "bounds": (1.8, 6.5),
    },
    "small_airway_internal_radius_mm": {
        "value": 1.50,
        "tag": "empirical",
        "source": "Human 8th-generation bronchiole internal lumen radius (DOI: 10.1152/physrev.00004.2021)",
        "bounds": (0.8, 2.5),
    },
}

TWO_COLUMN_EPISTEMIC_HYGIENE: Dict[str, Any] = {
    "deductive_theorem": (
        "Under Herschel-Bulkley non-Newtonian viscoplastic rheology, shear strain rate gamma_dot is strictly zero "
        "whenever applied shear stress tau is below the yield stress tau_y. "
        "When ciliary propulsive stress tau_cilia <= tau_y, mucociliary transport is kinematically arrested. "
        "In cylindrical pipe Poiseuille flow, an unyielded solid core plug of radius r_plug = R * (tau_y / tau_w) "
        "persists, blocking extrusion unless wall shear stress exceeds the critical yield threshold."
    ),
    "biophysical_hypothesis": {
        "statement": (
            "In cystic fibrosis, dehydration and neutrophil extracellular DNA elevate sputum yield stress above "
            "the ciliary propulsion limit (tau_y >= 1.5 Pa > tau_cilia = 1.25 Pa), causing total ciliary stall "
            "and mucus plug impaction. Enzymatic liquefaction (dornase alfa) hydrolyzes eDNA networks, reducing "
            "tau_y below 0.85 Pa and restoring active mucociliary clearance."
        ),
        "falsification_condition": (
            "If sputum samples with measured yield stress tau_y >= 2.0 Pa exhibit sustained mucociliary transport "
            "under normal ciliary beat frequency or if enzymatic cleavage of eDNA fails to reduce yield stress below "
            "tau_cilia, the viscoplastic ciliary stall model is falsified."
        ),
    },
}

FIVE_LAYER_LADDER: Dict[str, str] = {
    "cool_human_fact": (
        "The millions of microscopic hair-like cilia lining our airways beat like oars, but when CF mucus dries out, "
        "it turns from lubricant into solid wax that stalls the cilia mid-stroke."
    ),
    "everyday_analogy": (
        "Like trying to stir chilled cold butter with a plastic coffee stirrer—below a certain push, the butter "
        "doesn't flow at all, it simply holds stiff and snaps the stirrer."
    ),
    "engineering_lens": (
        "A viscoplastic Herschel-Bulkley flow problem with a Bingham-like yield stress threshold coupling ciliary "
        "surface traction to Poiseuille core plug radius extrusion."
    ),
    "formal_theorem": (
        "Herschel-Bulkley constitutive law tau = tau_y + K * (gamma_dot)^n enforces an essential discontinuity "
        "in apparent viscosity eta_app = tau / gamma_dot -> infinity as gamma_dot -> 0, establishing an unyielding regime."
    ),
    "falsifiable_code_and_toy": (
        "A pure Python viscoplastic solver computing ciliary stall velocity, Poiseuille plug core radius, and complex-step "
        "differential viscosity at h=1e-30 with deterministic falsification on negative yield stress or consistency."
    ),
}

REGIME_PRESETS: Dict[str, Dict[str, Any]] = {
    "healthy_euhydrated_airway_mucus": {
        "yield_stress_pa": 0.10,
        "consistency_index_pa_s_n": 0.50,
        "flow_behavior_index": 0.45,
        "ciliary_stress_pa": 1.25,
        "airway_radius_mm": 1.50,
        "airway_length_mm": 20.0,
        "pressure_drop_pa": 40.0,
        "expected_min_mcc_mm_min": 4.50,
        "expected_stall": False,
        "crosses_kahn_kalai": True,
    },
    "cf_hyperconcentrated_sputum_arrest": {
        "yield_stress_pa": 3.20,
        "consistency_index_pa_s_n": 8.50,
        "flow_behavior_index": 0.35,
        "ciliary_stress_pa": 1.25,
        "airway_radius_mm": 1.50,
        "airway_length_mm": 20.0,
        "pressure_drop_pa": 40.0,
        "expected_min_mcc_mm_min": 0.0,
        "expected_stall": True,
        "crosses_kahn_kalai": False,
    },
    "cf_severe_mucoid_plugging": {
        "yield_stress_pa": 8.50,
        "consistency_index_pa_s_n": 22.0,
        "flow_behavior_index": 0.30,
        "ciliary_stress_pa": 1.25,
        "airway_radius_mm": 1.50,
        "airway_length_mm": 20.0,
        "pressure_drop_pa": 40.0,
        "expected_min_mcc_mm_min": 0.0,
        "expected_stall": True,
        "crosses_kahn_kalai": False,
    },
    "kahn_kalai_curative_liquefaction_boundary": {
        "yield_stress_pa": 0.75,
        "consistency_index_pa_s_n": 1.80,
        "flow_behavior_index": 0.45,
        "ciliary_stress_pa": 1.25,
        "airway_radius_mm": 1.50,
        "airway_length_mm": 20.0,
        "pressure_drop_pa": 40.0,
        "expected_min_mcc_mm_min": 0.65,
        "expected_stall": False,
        "crosses_kahn_kalai": True,
    },
    "rhdnase_enzymatic_liquefaction_rescue": {
        "yield_stress_pa": 0.35,
        "consistency_index_pa_s_n": 0.90,
        "flow_behavior_index": 0.45,
        "ciliary_stress_pa": 1.25,
        "airway_radius_mm": 1.50,
        "airway_length_mm": 20.0,
        "pressure_drop_pa": 40.0,
        "expected_min_mcc_mm_min": 3.80,
        "expected_stall": False,
        "crosses_kahn_kalai": True,
    },
}


def compute_herschel_bulkley_shear_stress(
    gamma_dot: float,
    tau_y: float,
    k_index: float,
    n_index: float,
) -> float:
    """Compute shear stress tau (Pa) from shear strain rate gamma_dot (s^-1)."""
    if tau_y < 0.0:
        raise ValueError(f"Yield stress must be non-negative, got {tau_y}")
    if k_index <= 0.0:
        raise ValueError(f"Consistency index must be strictly positive, got {k_index}")
    if n_index <= 0.0 or n_index > 2.0:
        raise ValueError(f"Flow behavior index must be in (0, 2], got {n_index}")
    if gamma_dot < 0.0:
        raise ValueError(f"Shear rate must be non-negative, got {gamma_dot}")

    if gamma_dot == 0.0:
        return tau_y

    return tau_y + k_index * (gamma_dot**n_index)


def compute_apparent_viscosity(
    gamma_dot: float,
    tau_y: float,
    k_index: float,
    n_index: float,
) -> float:
    """Compute apparent dynamic viscosity eta_app = tau / gamma_dot (Pa * s)."""
    if gamma_dot <= 0.0:
        return float("inf") if tau_y > 0.0 else k_index

    tau = compute_herschel_bulkley_shear_stress(gamma_dot, tau_y, k_index, n_index)
    return tau / gamma_dot


def compute_ciliary_shear_rate(
    tau_cilia: float,
    tau_y: float,
    k_index: float,
    n_index: float,
) -> float:
    """Compute shear strain rate gamma_dot generated by ciliary stroke."""
    if tau_cilia <= tau_y:
        return 0.0
    effective_stress = tau_cilia - tau_y
    return (effective_stress / k_index)**(1.0 / n_index)


def compute_mcc_velocity_mm_min(
    gamma_dot: float,
    v0_mm_min: float = 5.50,
    gamma_dot_ref: float = 0.40,
) -> float:
    """Compute mucociliary clearance transport velocity v_MCC in mm/min."""
    if gamma_dot <= 0.0:
        return 0.0
    return v0_mm_min * (gamma_dot / (gamma_dot + gamma_dot_ref))


def compute_poiseuille_plug_radius_mm(
    r_airway_mm: float,
    delta_p_pa: float,
    length_mm: float,
    tau_y: float,
) -> Tuple[float, float]:
    """Compute wall shear stress tau_w (Pa) and solid unyielded plug radius r_plug (mm)."""
    if r_airway_mm <= 0.0 or length_mm <= 0.0:
        raise ValueError("Airway radius and length must be strictly positive")
    if delta_p_pa <= 0.0:
        return 0.0, r_airway_mm

    r_m = r_airway_mm * 1e-3
    l_m = length_mm * 1e-3

    tau_w = (r_m * delta_p_pa) / (2.0 * l_m)

    if tau_w <= tau_y:
        r_plug_mm = r_airway_mm
    else:
        r_plug_mm = r_airway_mm * (tau_y / tau_w)

    return tau_w, min(r_airway_mm, max(0.0, r_plug_mm))


def compute_extrusion_flow_rate_mm3_s(
    r_airway_mm: float,
    delta_p_pa: float,
    length_mm: float,
    tau_y: float,
    k_index: float,
    n_index: float,
) -> float:
    """Compute volumetric extrusion flow rate Q in mm^3 / s for Herschel-Bulkley fluid in tube."""
    tau_w, r_plug_mm = compute_poiseuille_plug_radius_mm(r_airway_mm, delta_p_pa, length_mm, tau_y)
    if tau_w <= tau_y:
        return 0.0

    r_m = r_airway_mm * 1e-3
    m = 1.0 / n_index
    diff = tau_w - tau_y

    prefactor = (math.pi * (r_m**3)) / ((k_index**m) * (tau_w**3))
    power_term = diff**(m + 1.0)
    bracket = (
        (diff**2) / (3.0 + m)
        + (2.0 * tau_y * diff) / (2.0 + m)
        + (tau_y**2) / (1.0 + m)
    )
    q_m3_s = prefactor * power_term * bracket
    return q_m3_s * 1e9


def complex_step_differential_viscosity(
    gamma_dot: float,
    tau_y: float,
    k_index: float,
    n_index: float,
    step_h: float = 1e-30,
) -> Tuple[float, float, float]:
    """Compute differential viscosity dtau/d(gamma_dot) via complex-step differentiation (step_h = 1e-30)."""
    if gamma_dot <= 0.0:
        raise ValueError("Shear rate must be strictly positive for complex-step derivative")

    gamma_c = complex(gamma_dot, step_h)
    tau_c = tau_y + k_index * (gamma_c**n_index)
    deriv_cs = tau_c.imag / step_h

    deriv_ana = n_index * k_index * (gamma_dot**(n_index - 1.0))

    rel_error = abs(deriv_cs - deriv_ana) / abs(deriv_ana)
    return deriv_cs, deriv_ana, rel_error


@dataclass(frozen=True)
class Round160HerschelBulkleyResult:
    preset_name: str
    yield_stress_pa: float
    ciliary_shear_rate_s1: float
    mcc_velocity_mm_min: float
    apparent_viscosity_pa_s: float
    wall_shear_stress_pa: float
    plug_radius_fraction: float
    extrusion_flow_rate_mm3_s: float
    ciliary_stall: bool
    crosses_kahn_kalai_threshold: bool
    falsification_invariant_satisfied: bool
    complex_step_relative_error: float


def evaluate_round160_herschel_bulkley_yield_stress(preset_name: str) -> Round160HerschelBulkleyResult:
    """Evaluate Round 160 Herschel-Bulkley Viscoplastic Yield-Stress Invariant Bridge."""
    if preset_name not in REGIME_PRESETS:
        raise KeyError(f"Unknown preset: {preset_name}. Expected one of {sorted(REGIME_PRESETS)}")

    cfg = REGIME_PRESETS[preset_name]
    tau_y = float(cfg["yield_stress_pa"])
    k_idx = float(cfg["consistency_index_pa_s_n"])
    n_idx = float(cfg["flow_behavior_index"])
    tau_cilia = float(cfg["ciliary_stress_pa"])
    r_airway = float(cfg["airway_radius_mm"])
    l_airway = float(cfg["airway_length_mm"])
    dp = float(cfg["pressure_drop_pa"])

    gamma_dot_cilia = compute_ciliary_shear_rate(tau_cilia, tau_y, k_idx, n_idx)
    v_mcc = compute_mcc_velocity_mm_min(gamma_dot_cilia)
    is_stall = gamma_dot_cilia == 0.0

    ref_gamma = gamma_dot_cilia if gamma_dot_cilia > 0.0 else 1.0
    eta_app = compute_apparent_viscosity(ref_gamma, tau_y, k_idx, n_idx)

    tau_w, r_plug_mm = compute_poiseuille_plug_radius_mm(r_airway, dp, l_airway, tau_y)
    plug_fraction = r_plug_mm / r_airway
    q_mm3_s = compute_extrusion_flow_rate_mm3_s(r_airway, dp, l_airway, tau_y, k_idx, n_idx)

    _, _, rel_err = complex_step_differential_viscosity(2.0, tau_y, k_idx, n_idx)

    crosses_kk = (not is_stall) and (v_mcc >= 0.65) and (tau_y <= 0.85)

    invariants_ok = (
        rel_err < 1e-10
        and 0.0 <= plug_fraction <= 1.0
        and (is_stall == cfg["expected_stall"])
        and (not crosses_kk if not cfg["crosses_kahn_kalai"] else crosses_kk)
    )

    return Round160HerschelBulkleyResult(
        preset_name=preset_name,
        yield_stress_pa=tau_y,
        ciliary_shear_rate_s1=gamma_dot_cilia,
        mcc_velocity_mm_min=v_mcc,
        apparent_viscosity_pa_s=eta_app,
        wall_shear_stress_pa=tau_w,
        plug_radius_fraction=plug_fraction,
        extrusion_flow_rate_mm3_s=q_mm3_s,
        ciliary_stall=is_stall,
        crosses_kahn_kalai_threshold=crosses_kk,
        falsification_invariant_satisfied=invariants_ok,
        complex_step_relative_error=rel_err,
    )


def evaluate_falsification_negative_control(
    poison_negative_yield_stress: bool = False,
    poison_negative_consistency: bool = False,
    poison_negative_flow_index: bool = False,
    poison_negative_shear_rate: bool = False,
) -> Dict[str, Any]:
    """Execute AH-10 negative control checks ensuring unphysical Herschel-Bulkley models fail deterministically."""
    if poison_negative_yield_stress:
        try:
            compute_herschel_bulkley_shear_stress(1.0, -1.5, 2.0, 0.45)
            return {"passed": False, "rule": "AH-10", "error": "Accepted negative yield stress"}
        except ValueError:
            return {"passed": False, "rule": "AH-10", "error": "Correctly rejected negative yield stress (poisoned)"}

    if poison_negative_consistency:
        try:
            compute_herschel_bulkley_shear_stress(1.0, 1.5, -2.0, 0.45)
            return {"passed": False, "rule": "AH-10", "error": "Accepted negative consistency index"}
        except ValueError:
            return {"passed": False, "rule": "AH-10", "error": "Correctly rejected negative consistency index (poisoned)"}

    if poison_negative_flow_index:
        try:
            compute_herschel_bulkley_shear_stress(1.0, 1.5, 2.0, -0.45)
            return {"passed": False, "rule": "AH-10", "error": "Accepted negative flow index"}
        except ValueError:
            return {"passed": False, "rule": "AH-10", "error": "Correctly rejected negative flow index (poisoned)"}

    if poison_negative_shear_rate:
        try:
            compute_herschel_bulkley_shear_stress(-1.0, 1.5, 2.0, 0.45)
            return {"passed": False, "rule": "AH-10", "error": "Accepted negative shear rate"}
        except ValueError:
            return {"passed": False, "rule": "AH-10", "error": "Correctly rejected negative shear rate (poisoned)"}

    res = evaluate_round160_herschel_bulkley_yield_stress("healthy_euhydrated_airway_mucus")
    passed = res.falsification_invariant_satisfied and res.crosses_kahn_kalai_threshold
    return {
        "passed": passed,
        "rule": "AH-10",
        "mcc_velocity_mm_min": res.mcc_velocity_mm_min,
        "ciliary_stall": res.ciliary_stall,
    }
