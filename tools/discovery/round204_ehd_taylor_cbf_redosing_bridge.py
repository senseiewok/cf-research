#!/usr/bin/env python3
"""
Round 204: Electrohydrodynamic (EHD) Taylor-Cone Aerosol Atomization, Mucin H1 Persistent-Homology Sterics, & Control Barrier Function (CBF) Safe Redosing Bridge (T-0897).

Cross-Domain Isomorphism & CF Anchor (INV-CF-ANCHOR + INV-MECH-ISOMORPHISM):
  Couples three Informatico formal operators to chronic nonviral inhaled site-specific m1A LNP-sup-tRNA therapy for
  cystic fibrosis nonsense genotypes (Chen, Cui, Li et al., Science 2026, DOI: 10.1126/science.aeb0054; CFFPR 2024:
  8.0% of registered individuals genotype-ineligible for CFTR modulators):
  1. Informatico Folio P.09 — Electrohydrodynamic (EHD) Taylor-Cone Monodisperse Atomization:
     Under Maxwell electrostatic-capillary stress equilibrium at Taylor semi-vertical angle theta_T = 49.3 deg
     (Taylor 1964, DOI: 10.1098/rspa.1964.0151), electrospray jet breakup generates monodisperse respirable droplets
     (d_aero ~ 2.6-2.8 um, GSD ~ 1.12-1.15) at low shear gamma_dot ~ 3.2e3-4.0e3 s^-1, avoiding high-shear bilayer
     rupture imposed by vibrating-mesh orifices (gamma_dot ~ 1.8e5 s^-1) and off-equilibrium cone-angle jet whipping:
       eta_intact(gamma_dot) = 1.0 / (1.0 + (gamma_dot / gamma_crit)^2)
       r_h_post(gamma_dot) = r_h0_nm * (1.0 + beta_agg * (1.0 - eta_intact(gamma_dot)))
       gsd_eff(gsd, theta_deg) = gsd * (1.0 + 0.015 * (theta_deg - theta_T)^2)
       f_bronch(d_aero_um, gsd_eff) = f_bronch_max * exp(-((d_aero_um - d_opt_um)^2) / (2.0 * sigma_d_um^2))
                                      / (1.0 + lambda_gsd * (gsd_eff - 1.0)^2)
     where gamma_crit = 4.5e4 s^-1, beta_agg = 0.85, f_bronch_max = 0.72, d_opt_um = 2.8 um, sigma_d_um = 1.35 um,
     and lambda_gsd = 0.90.

  2. Informatico Folio M.10 — CF Mucin Mesh Persistent-Homology H1 Bottleneck Sterics:
     Modeling entangled hyperconcentrated CF MUC5B/MUC5AC + eDNA fibers via alpha-complex persistent homology H1
     1-cycle filtration (Edelsbrunner, Letscher, Zomorodian 2002, DOI: 10.1007/s00454-002-2885-2; Scite records 1 contrasting citation [observed: Scite ML citation classification, accessed 2026-11-11: supporting=8, contrasting=1, mentioning=2230 across 2416 citingPublications])
     defines a topological bottleneck pore radius r_bottleneck_nm. Percolation requires r_h_post <= r_bottleneck_nm
     paired with muco-inert PEG/zwitterionic shielding s_shield in [0, 1]:
       T_h1 = (s_shield / (1.0 + exp((r_h_post - r_bottleneck_nm) / w_pore_nm))) * exp(-((r_h_post / r_steric_nm)^2))
     where w_pore_nm = 8.0 nm and r_steric_nm = 95.0 nm.

  3. Informatico Folio E.12 — Control Barrier Function (CBF) Forward-Invariant Safe Redosing:
     Chronic inhaled LNP-sup-tRNA redosing at interval delta_t_h > 0 drives periodic innate TLR7/8 + anti-PEG
     immune state peaks I_tlr_peak(delta_t_h) = delta_i_dose / (1.0 - exp(-k_resolve * delta_t_h)) with decay rate
     k_resolve = 0.032 h^-1. Zeroing Control Barrier Function safety (Ames et al. 2017, DOI: 10.1109/TAC.2016.2638961)
     enforces forward invariance of C = {x : h(x) >= 0} via:
       h_cbf(delta_t_h) = I_max - I_tlr_peak(delta_t_h) >= 0
     where I_max = 0.35. Violating h_cbf < 0 triggers accelerated macrophage/anti-PEG clearance:
       phi_clearance = 1.0 / (1.0 + exp(kappa_cbf * (I_tlr_peak(delta_t_h) - I_max)))
     while intracellular functional sup-tRNA + apical CFTR cycle maintenance (k_trna_eff = 0.018 h^-1) obeys:
       M_cycle(delta_t_h) = (1.0 - exp(-k_trna_eff * delta_t_h)) / (k_trna_eff * delta_t_h).
     Net chronic sustained apical delivery index:
       J_chronic = eta_intact * f_bronch * T_h1 * phi_clearance * M_cycle.

  4. Kahn-Kalai Curative Threshold & Deterministic Invariants:
     - Post-atomization LNP integrity eta_intact >= 0.90 and steric condition r_h_post <= r_bottleneck_nm.
     - Topological mucus H1 transmission T_h1 >= 0.50.
     - Forward-invariant CBF margin h_cbf >= 0.0 (equivalently I_tlr_peak <= I_max = 0.35).
     - Chronic sustained therapeutic delivery index J_chronic >= 0.18.
     - Complex-step sensitivity d(J_chronic)/d(delta_t_h) at h = 1e-30 matches central finite differences (rel_err <= 1e-6).
     - Poisoned negative controls (AH-10): unshielded >250 nm (d=280 nm, r_h0=140 nm) LNPs are topologically trapped in
       CF sputum H1 pores (T_h1 < 1e-4), naive 12 h fixed-interval overdosing violates h_cbf < 0, and non-physical
       negative inputs raise ValueError deterministically.

CF Biophysical Domain Mapping (Column B - Popperian Falsifiable Hypothesis):
  - EHD Taylor-Cone + CBF Adaptive 72 h Curative Regimen (gamma_dot=4000 s^-1, r_h0=38 nm, delta_t=72 h):
    eta_intact = 0.9922 >= 0.90, r_h_post = 38.25 nm <= 70.0 nm, T_h1 = 0.7678 >= 0.50, h_cbf = +0.1056 >= 0,
    J_chronic = 0.2621 >= 0.18 (Pass).
  - EHD Taylor-Cone + Ultra-Low TLR 48 h Curative Regimen (gamma_dot=3200 s^-1, r_h0=34 nm, delta_t=48 h):
    eta_intact = 0.9950 >= 0.90, r_h_post = 34.15 nm <= 68.0 nm, T_h1 = 0.8229 >= 0.50, h_cbf = +0.1716 >= 0,
    J_chronic = 0.3686 >= 0.18 (Pass).
  - Vibrating-Mesh High-Shear Rupture Failure Regimen (gamma_dot=1.8e5 s^-1, d_aero=4.5 um, gsd=1.85):
    eta_intact = 0.0588 < 0.90, r_h_post = 75.60 nm > 70.0 nm, J_chronic = 0.0009 < 0.18 (Fail).
  - Unshielded 280 nm LNP H1 Mucus Entrapment Control (r_h0=140 nm > 125 nm, s_shield=0.12):
    r_h_post = 140.93 nm > 70.0 nm, T_h1 = 1.87e-6 < 0.50, J_chronic = 6.40e-7 < 0.18 (Fail).
  - Naive Fixed 12 h Overdosing CBF Violation Control (delta_i_dose=0.24, delta_t=12 h < delta_t_min=36.17 h):
    I_tlr_peak = 0.7527 > 0.35, h_cbf = -0.4027 < 0, phi_clearance = 0.0007, J_chronic = 0.0003 < 0.18 (Fail).

Research simulation only; not medical advice.
"""
from __future__ import annotations

import cmath
import math
from dataclasses import dataclass
from typing import Any, Dict, Tuple

RESEARCH_BOUNDARY: str = "Research simulation only; not medical advice."

PRIMARY_CITATIONS: Tuple[str, ...] = (
    "DOI: 10.1126/science.aeb0054",
    "DOI: 10.1098/rspa.1964.0151",
    "DOI: 10.1007/s00454-002-2885-2",
    "DOI: 10.1109/TAC.2016.2638961",
)

CROSS_DOMAIN_GOVERNANCE: Dict[str, object] = {
    "inv_cf_anchor": (
        "Anchored to chronic nonviral airway delivery of chemically modified suppressor tRNA LNPs for CFTR nonsense "
        "variants W1282X/G542X (Chen et al., Science 2026, DOI: 10.1126/science.aeb0054; CFFPR 2024: 8.0% modulator-ineligible)."
    ),
    "inv_mech_isomorphism": (
        "Couples Informatico Folio P.09 (Taylor 1964 EHD cone atomization, DOI: 10.1098/rspa.1964.0151), Folio M.10 "
        "(Edelsbrunner et al. 2002 H1 persistent-homology bottleneck filtration, DOI: 10.1007/s00454-002-2885-2; Scite records 1 contrasting citation [observed: Scite ML citation classification, accessed 2026-11-11: supporting=8, contrasting=1, mentioning=2230 across 2416 citingPublications]), "
        "and Folio E.12 (Ames et al. 2017 Control Barrier Function forward invariance, DOI: 10.1109/TAC.2016.2638961)."
    ),
    "clinical_extrapolation_permitted": False,
}

# Biophysical & Control Constants
DEFAULT_TAYLOR_CONE_HALF_ANGLE_DEG: float = 49.3
DEFAULT_GAMMA_CRIT_S_INV: float = 45000.0
DEFAULT_BETA_SHEAR_AGG: float = 0.85
DEFAULT_F_BRONCH_MAX: float = 0.72
DEFAULT_D_OPT_AERO_UM: float = 2.8
DEFAULT_SIGMA_AERO_UM: float = 1.35
DEFAULT_LAMBDA_GSD_DISPERSION: float = 0.90
DEFAULT_W_PORE_TRANSITION_NM: float = 8.0
DEFAULT_R_STERIC_CHARACTERISTIC_NM: float = 95.0
DEFAULT_K_IMMUNE_RESOLVE_PER_H: float = 0.032
DEFAULT_I_TLR_MAX_CEILING: float = 0.35
DEFAULT_KAPPA_CBF_CLEARANCE: float = 18.0
DEFAULT_K_TRNA_PERSIST_PER_H: float = 0.018

PARAMETER_PROVENANCE: Dict[str, Dict[str, object]] = {
    "lnp_bilayer_critical_shear_rate_s_inv": {
        "value": DEFAULT_GAMMA_CRIT_S_INV,
        "tag": "phenomenological_heuristic",
        "source": (
            "Critical nebulization shear rate gamma_crit governing LNP bilayer rupture vs Taylor-cone EHD preservation "
            "(exact 4.5e4 s^-1 scalar not stated in Taylor 1964, DOI: 10.1098/rspa.1964.0151, or Chen et al., Science 2026, "
            "DOI: 10.1126/science.aeb0054)"
        ),
        "bounds": (1.5e4, 1.2e5),
    },
    "shear_induced_lnp_fusion_inflation_beta": {
        "value": DEFAULT_BETA_SHEAR_AGG,
        "tag": "phenomenological_heuristic",
        "source": (
            "Fractional hydrodynamic radius inflation factor beta_agg upon high-shear orifice nebulization "
            "(exact 0.85 dimensionless scalar not stated in Taylor 1964, DOI: 10.1098/rspa.1964.0151)"
        ),
        "bounds": (0.20, 1.80),
    },
    "mucin_h1_bottleneck_transition_width_nm": {
        "value": DEFAULT_W_PORE_TRANSITION_NM,
        "tag": "phenomenological_heuristic",
        "source": (
            "Sigmoid transition width w_pore_nm around the persistent-homology H1 bottleneck pore radius in CF sputum "
            "(exact 8.0 nm scalar not stated in Edelsbrunner et al. 2002, DOI: 10.1007/s00454-002-2885-2; Scite records 1 contrasting citation [observed: Scite ML citation classification, accessed 2026-11-11: supporting=8, contrasting=1, mentioning=2230 across 2416 citingPublications])"
        ),
        "bounds": (3.0, 20.0),
    },
    "tlr78_innate_resolution_rate_per_h": {
        "value": DEFAULT_K_IMMUNE_RESOLVE_PER_H,
        "tag": "phenomenological_heuristic",
        "source": (
            "First-order airway innate TLR7/8 and anti-PEG priming resolution rate k_resolve between inhaled doses "
            "(exact 0.032 h^-1 scalar not stated in Chen et al., Science 2026, DOI: 10.1126/science.aeb0054)"
        ),
        "bounds": (0.012, 0.080),
    },
    "cbf_innate_safety_ceiling_imax": {
        "value": DEFAULT_I_TLR_MAX_CEILING,
        "tag": "phenomenological_heuristic",
        "source": (
            "Control Barrier Function safe-set ceiling I_max for cumulative TLR7/8 / anti-PEG activation "
            "(exact 0.35 dimensionless ceiling not stated in Ames et al. 2017, DOI: 10.1109/TAC.2016.2638961)"
        ),
        "bounds": (0.15, 0.60),
    },
    "chronic_suptrna_rescue_decay_rate_per_h": {
        "value": DEFAULT_K_TRNA_PERSIST_PER_H,
        "tag": "phenomenological_heuristic",
        "source": (
            "Effective intracellular persistence decay rate k_trna_eff of rescued CFTR + modified sup-tRNA activity "
            "(exact 0.018 h^-1 scalar not stated in Chen et al., Science 2026, DOI: 10.1126/science.aeb0054)"
        ),
        "bounds": (0.008, 0.045),
    },
}

TWO_COLUMN_EPISTEMIC_HYGIENE: Dict[str, Dict[str, str]] = {
    "deductive_column_a": {
        "formal_theorem": (
            "For any single-dose innate impulse 0 < delta_i_dose < I_max, the periodic impulse-decay map admits a "
            "unique forward-invariant Control Barrier Function boundary delta_t_min = -ln(1 - delta_i_dose / I_max) / k_resolve "
            "such that h_cbf(delta_t_h) = I_max - I_tlr_peak(delta_t_h) >= 0 iff delta_t_h >= delta_t_min, while "
            "J_chronic(delta_t_h) is holomorphic on Re(delta_t_h) > 0."
        ),
        "invariants": (
            "0 < eta_intact <= 1; r_h_post >= r_h0_nm; 0 <= f_bronch, T_h1, phi_clearance, M_cycle, J_chronic <= 1; "
            "d(eta_intact)/d(gamma_dot) < 0; d(T_h1)/d(r_h0_nm) < 0; d(h_cbf)/d(delta_t_h) > 0; "
            "complex-step relative error <= 1e-6."
        ),
    },
    "popperian_column_b": {
        "falsifiable_hypothesis": (
            "Low-shear EHD Taylor-cone atomization of sub-80 nm muco-inert LNP-sup-tRNAs coupled to a CBF-constrained "
            "48-72 h redosing schedule satisfies r_h_post <= r_bottleneck_nm and h_cbf >= 0 in CF sputum models "
            "(J_chronic >= 0.18), whereas vibrating-mesh high-shear nebulization, unshielded >250 nm LNPs, or naive "
            "12 h fixed-interval overdosing deterministically fail."
        ),
        "empirical_falsification_threshold": (
            "Falsified if piezoelectric vibrating-mesh nebulization at gamma_dot >= 1.5e5 s^-1 preserves >= 90% LNP "
            "encapsulation without radius inflation, if unshielded >250 nm LNPs percolate through CF sputum H1 pores "
            "(Edelsbrunner et al. 2002, DOI: 10.1007/s00454-002-2885-2; Scite records 1 contrasting citation [observed: Scite ML citation classification, accessed 2026-11-11: supporting=8, contrasting=1, mentioning=2230 across 2416 citingPublications]) "
            "with T_h1 >= 0.50, or if 12 h redosing avoids cumulative TLR7/8 / anti-PEG clearance escalation."
        ),
    },
}

FIVE_LAYER_LADDER: Dict[str, str] = {
    "cool_human_fact": (
        "Because suppressor tRNAs rescue CFTR mRNA transiently without permanently editing DNA, turning a one-time "
        "laboratory rescue for the 8.0% modulator-ineligible CF cohort into lifelong breathing support requires gentle "
        "inhaled redosing that neither shreds the lipid carrier nor provokes airway immunity."
    ),
    "everyday_analogy": (
        "Instead of forcing delicate soap-bubble couriers through a violent kitchen sieve every few hours until the "
        "building's alarm locks the doors, a quiet electric cone guides uniform mist droplets through the widest mesh "
        "tunnels on a paced schedule that lets the alarm reset cleanly."
    ),
    "engineering_lens": (
        "Electrohydrodynamic Taylor-cone atomization keeps shear below gamma_crit to preserve sub-80 nm LNP radius "
        "inside the CF mucin H1 bottleneck pore, while a zeroing Control Barrier Function paces dose intervals so "
        "cumulative TLR7/8 and anti-PEG priming never cross I_max."
    ),
    "formal_theorem": (
        "The composite delivery functional J_chronic = eta_intact * f_bronch * T_h1 * phi_clearance * M_cycle attains "
        "its constrained interior maximum inside the CBF forward-invariant safe set delta_t_h >= delta_t_min."
    ),
    "falsifiable_code_and_toy": (
        "evaluate_round204_ehd_taylor_cbf_redosing(regime_key) evaluates EHD vs mesh shear rupture, H1 mucus sterics, "
        "and CBF redosing invariance across 5 curative and negative-control regimens."
    ),
}

REGIME_PRESETS: Dict[str, Dict[str, float]] = {
    "ehd_taylor_cbf_adaptive_72h_curative": {
        "gamma_dot_s_inv": 4000.0,
        "d_aero_um": 2.8,
        "gsd": 1.15,
        "r_h0_nm": 38.0,
        "r_bottleneck_nm": 70.0,
        "s_shield": 0.92,
        "delta_i_dose": 0.22,
        "delta_t_h": 72.0,
    },
    "ehd_taylor_cbf_ultralow_tlr_48h_curative": {
        "gamma_dot_s_inv": 3200.0,
        "d_aero_um": 2.6,
        "gsd": 1.12,
        "r_h0_nm": 34.0,
        "r_bottleneck_nm": 68.0,
        "s_shield": 0.95,
        "delta_i_dose": 0.14,
        "delta_t_h": 48.0,
    },
    "vibrating_mesh_high_shear_rupture_fail": {
        "gamma_dot_s_inv": 180000.0,
        "d_aero_um": 4.5,
        "gsd": 1.85,
        "r_h0_nm": 42.0,
        "r_bottleneck_nm": 70.0,
        "s_shield": 0.88,
        "delta_i_dose": 0.22,
        "delta_t_h": 72.0,
    },
    "unshielded_280nm_lnp_h1_mucus_trapped_fail": {
        "gamma_dot_s_inv": 4000.0,
        "d_aero_um": 2.8,
        "gsd": 1.15,
        "r_h0_nm": 140.0,  # 280 nm diameter (> 250 nm unshielded LNP negative control)
        "r_bottleneck_nm": 70.0,
        "s_shield": 0.12,
        "delta_i_dose": 0.22,
        "delta_t_h": 72.0,
    },
    "naive_fixed_12h_overdosing_cbf_violation_fail": {
        "gamma_dot_s_inv": 4000.0,
        "d_aero_um": 2.8,
        "gsd": 1.15,
        "r_h0_nm": 38.0,
        "r_bottleneck_nm": 70.0,
        "s_shield": 0.92,
        "delta_i_dose": 0.24,
        "delta_t_h": 12.0,
    },
}


@dataclass(frozen=True)
class Round204EHDTaylorCBFResult:
    regime_key: str
    gamma_dot_s_inv: float
    d_aero_um: float
    gsd: float
    r_h0_nm: float
    r_h_post_nm: float
    r_bottleneck_nm: float
    s_shield: float
    delta_i_dose: float
    delta_t_h: float
    delta_t_cbf_min_h: float
    eta_intact: float
    f_bronch_deposition: float
    t_h1_mucus_transmission: float
    h1_steric_pass: bool
    i_tlr_peak: float
    h_cbf_margin: float
    cbf_forward_invariant_safe: bool
    phi_clearance_retention: float
    m_cycle_maintenance: float
    j_chronic_delivery_index: float
    crosses_kahn_kalai_threshold: bool
    falsification_invariant_satisfied: bool
    complex_step_relative_error: float
    research_boundary: str = RESEARCH_BOUNDARY


def _exp_num(x: float | complex) -> float | complex:
    return cmath.exp(x) if isinstance(x, complex) else math.exp(x)


def _logistic_neg(z: float | complex) -> float | complex:
    """Overflow-safe evaluation of 1 / (1 + exp(z)) on Re(z) in (-inf, +inf)."""
    real_part = z.real if isinstance(z, complex) else float(z)
    if real_part >= 0.0:
        e_neg = _exp_num(-z)
        return e_neg / (1.0 + e_neg)
    return 1.0 / (1.0 + _exp_num(z))


def _one_minus_exp_neg(z: float | complex) -> float | complex:
    """Cancellation-safe evaluation of 1 - exp(-z) for small positive z while preserving complex-step."""
    if isinstance(z, complex):
        return 1.0 - cmath.exp(-z)
    return -math.expm1(-float(z))


def compute_ehd_atomization_and_lnp_integrity(
    gamma_dot_s_inv: float | complex,
    d_aero_um: float | complex,
    gsd: float | complex,
    r_h0_nm: float | complex,
    cone_half_angle_deg: float | complex = DEFAULT_TAYLOR_CONE_HALF_ANGLE_DEG,
) -> Tuple[float | complex, float | complex, float | complex]:
    """Computes intact LNP retention eta_intact, post-atomization radius r_h_post_nm, and bronchial deposition f_bronch."""
    if isinstance(gamma_dot_s_inv, (int, float)) and gamma_dot_s_inv < 0.0:
        raise ValueError(f"gamma_dot_s_inv must be non-negative, got {gamma_dot_s_inv}")
    if isinstance(d_aero_um, (int, float)) and d_aero_um <= 0.0:
        raise ValueError(f"d_aero_um must be positive, got {d_aero_um}")
    if isinstance(gsd, (int, float)) and gsd < 1.0:
        raise ValueError(f"gsd must be >= 1.0, got {gsd}")
    if isinstance(r_h0_nm, (int, float)) and r_h0_nm <= 0.0:
        raise ValueError(f"r_h0_nm must be positive, got {r_h0_nm}")
    if isinstance(cone_half_angle_deg, (int, float)) and (cone_half_angle_deg <= 0.0 or cone_half_angle_deg >= 90.0):
        raise ValueError(f"cone_half_angle_deg must be in (0, 90), got {cone_half_angle_deg}")

    shear_ratio_sq = (gamma_dot_s_inv / DEFAULT_GAMMA_CRIT_S_INV) ** 2
    eta_intact = 1.0 / (1.0 + shear_ratio_sq)
    r_h_post = r_h0_nm * (1.0 + DEFAULT_BETA_SHEAR_AGG * (1.0 - eta_intact))
    aero_gauss = _exp_num(-((d_aero_um - DEFAULT_D_OPT_AERO_UM) ** 2) / (2.0 * (DEFAULT_SIGMA_AERO_UM ** 2)))
    gsd_eff = gsd * (1.0 + 0.015 * ((cone_half_angle_deg - DEFAULT_TAYLOR_CONE_HALF_ANGLE_DEG) ** 2))
    gsd_penalty = 1.0 + DEFAULT_LAMBDA_GSD_DISPERSION * ((gsd_eff - 1.0) ** 2)
    f_bronch = DEFAULT_F_BRONCH_MAX * aero_gauss / gsd_penalty
    return eta_intact, r_h_post, f_bronch


def compute_h1_mucus_persistent_homology_transmission(
    r_h_post_nm: float | complex,
    r_bottleneck_nm: float | complex,
    s_shield: float | complex,
) -> float | complex:
    """Computes CF mucin H1 bottleneck steric-topological transmission T_h1."""
    if isinstance(r_h_post_nm, (int, float)) and r_h_post_nm <= 0.0:
        raise ValueError(f"r_h_post_nm must be positive, got {r_h_post_nm}")
    if isinstance(r_bottleneck_nm, (int, float)) and r_bottleneck_nm <= 0.0:
        raise ValueError(f"r_bottleneck_nm must be positive, got {r_bottleneck_nm}")
    if isinstance(s_shield, (int, float)) and (s_shield < 0.0 or s_shield > 1.0):
        raise ValueError(f"s_shield must be in [0, 1], got {s_shield}")

    bottleneck_gate = _logistic_neg((r_h_post_nm - r_bottleneck_nm) / DEFAULT_W_PORE_TRANSITION_NM)
    steric_decay = _exp_num(-((r_h_post_nm / DEFAULT_R_STERIC_CHARACTERISTIC_NM) ** 2))
    return s_shield * bottleneck_gate * steric_decay


def compute_cbf_redosing_invariance(
    delta_i_dose: float | complex,
    delta_t_h: float | complex,
) -> Tuple[float | complex, float | complex, float | complex, float | complex]:
    """Computes periodic peak TLR7/8 state, CBF safety margin h_cbf, clearance retention, and cycle maintenance."""
    if isinstance(delta_i_dose, (int, float)) and (delta_i_dose <= 0.0 or delta_i_dose >= DEFAULT_I_TLR_MAX_CEILING):
        raise ValueError(f"delta_i_dose must be in (0, {DEFAULT_I_TLR_MAX_CEILING}), got {delta_i_dose}")
    if isinstance(delta_t_h, (int, float)) and delta_t_h <= 0.0:
        raise ValueError(f"delta_t_h must be positive, got {delta_t_h}")

    denom_decay = _one_minus_exp_neg(DEFAULT_K_IMMUNE_RESOLVE_PER_H * delta_t_h)
    i_tlr_peak = delta_i_dose / denom_decay
    h_cbf = DEFAULT_I_TLR_MAX_CEILING - i_tlr_peak
    phi_clearance = _logistic_neg(DEFAULT_KAPPA_CBF_CLEARANCE * (i_tlr_peak - DEFAULT_I_TLR_MAX_CEILING))
    trna_arg = DEFAULT_K_TRNA_PERSIST_PER_H * delta_t_h
    m_cycle = _one_minus_exp_neg(trna_arg) / trna_arg
    return i_tlr_peak, h_cbf, phi_clearance, m_cycle


def _eval_j_chronic_from_delta_t(
    delta_t_h: float | complex,
    gamma_dot_s_inv: float,
    d_aero_um: float,
    gsd: float,
    r_h0_nm: float,
    r_bottleneck_nm: float,
    s_shield: float,
    delta_i_dose: float,
) -> float | complex:
    eta_intact, r_h_post, f_bronch = compute_ehd_atomization_and_lnp_integrity(
        gamma_dot_s_inv, d_aero_um, gsd, r_h0_nm
    )
    t_h1 = compute_h1_mucus_persistent_homology_transmission(r_h_post, r_bottleneck_nm, s_shield)
    _, _, phi_clear, m_cycle = compute_cbf_redosing_invariance(delta_i_dose, delta_t_h)
    return eta_intact * f_bronch * t_h1 * phi_clear * m_cycle


def complex_step_delta_t_derivative(
    delta_t_h: float,
    gamma_dot_s_inv: float,
    d_aero_um: float,
    gsd: float,
    r_h0_nm: float,
    r_bottleneck_nm: float,
    s_shield: float,
    delta_i_dose: float,
    h: float = 1e-30,
) -> Tuple[float, float, float]:
    """Verifies d(J_chronic)/d(delta_t_h) via complex-step perturbation f(x + 1j*1e-30).imag / 1e-30."""
    j_c = _eval_j_chronic_from_delta_t(
        complex(delta_t_h, h),
        gamma_dot_s_inv,
        d_aero_um,
        gsd,
        r_h0_nm,
        r_bottleneck_nm,
        s_shield,
        delta_i_dose,
    )
    cs_deriv = complex(j_c).imag / h

    delta = 1e-5
    j_plus = float(
        _eval_j_chronic_from_delta_t(
            delta_t_h + delta,
            gamma_dot_s_inv,
            d_aero_um,
            gsd,
            r_h0_nm,
            r_bottleneck_nm,
            s_shield,
            delta_i_dose,
        )
    )
    j_minus = float(
        _eval_j_chronic_from_delta_t(
            delta_t_h - delta,
            gamma_dot_s_inv,
            d_aero_um,
            gsd,
            r_h0_nm,
            r_bottleneck_nm,
            s_shield,
            delta_i_dose,
        )
    )
    fd_deriv = (j_plus - j_minus) / (2.0 * delta)
    denom = max(abs(cs_deriv), abs(fd_deriv), 1e-12)
    rel_err = abs(cs_deriv - fd_deriv) / denom
    return cs_deriv, fd_deriv, rel_err


def evaluate_round204_ehd_taylor_cbf_redosing(regime_key: str) -> Round204EHDTaylorCBFResult:
    """Evaluates EHD Taylor-cone atomization, mucin H1 sterics, and CBF redosing safety for a named regimen."""
    if regime_key not in REGIME_PRESETS:
        raise KeyError(f"Unknown regime '{regime_key}'. Valid regimes: {sorted(REGIME_PRESETS)}")

    p = REGIME_PRESETS[regime_key]
    gamma_dot = p["gamma_dot_s_inv"]
    d_aero = p["d_aero_um"]
    gsd = p["gsd"]
    r_h0 = p["r_h0_nm"]
    r_bot = p["r_bottleneck_nm"]
    s_sh = p["s_shield"]
    d_i = p["delta_i_dose"]
    d_t = p["delta_t_h"]

    eta_v, r_post_v, f_br_v = compute_ehd_atomization_and_lnp_integrity(gamma_dot, d_aero, gsd, r_h0)
    t_h1_v = compute_h1_mucus_persistent_homology_transmission(r_post_v, r_bot, s_sh)
    i_pk_v, h_cbf_v, phi_cl_v, m_cyc_v = compute_cbf_redosing_invariance(d_i, d_t)

    eta_int = float(eta_v)
    r_post = float(r_post_v)
    f_br = float(f_br_v)
    t_h1 = float(t_h1_v)
    i_pk = float(i_pk_v)
    h_cbf = float(h_cbf_v)
    phi_cl = float(phi_cl_v)
    m_cyc = float(m_cyc_v)

    j_chronic = eta_int * f_br * t_h1 * phi_cl * m_cyc
    dt_cbf_min = -math.log(1.0 - (d_i / DEFAULT_I_TLR_MAX_CEILING)) / DEFAULT_K_IMMUNE_RESOLVE_PER_H

    _, _, rel_err = complex_step_delta_t_derivative(
        d_t, gamma_dot, d_aero, gsd, r_h0, r_bot, s_sh, d_i
    )

    h1_pass = r_post <= r_bot
    cbf_safe = h_cbf >= 0.0

    inv_satisfied = (
        0.0 < eta_int <= 1.0
        and r_post >= r_h0
        and 0.0 <= f_br <= 1.0
        and 0.0 <= t_h1 <= 1.0
        and 0.0 <= phi_cl <= 1.0
        and 0.0 <= m_cyc <= 1.0
        and 0.0 <= j_chronic <= 1.0
        and rel_err <= 1e-6
    )

    crosses_threshold = (
        eta_int >= 0.90
        and h1_pass
        and t_h1 >= 0.50
        and cbf_safe
        and j_chronic >= 0.18
    )

    return Round204EHDTaylorCBFResult(
        regime_key=regime_key,
        gamma_dot_s_inv=gamma_dot,
        d_aero_um=d_aero,
        gsd=gsd,
        r_h0_nm=r_h0,
        r_h_post_nm=r_post,
        r_bottleneck_nm=r_bot,
        s_shield=s_sh,
        delta_i_dose=d_i,
        delta_t_h=d_t,
        delta_t_cbf_min_h=dt_cbf_min,
        eta_intact=eta_int,
        f_bronch_deposition=f_br,
        t_h1_mucus_transmission=t_h1,
        h1_steric_pass=h1_pass,
        i_tlr_peak=i_pk,
        h_cbf_margin=h_cbf,
        cbf_forward_invariant_safe=cbf_safe,
        phi_clearance_retention=phi_cl,
        m_cycle_maintenance=m_cyc,
        j_chronic_delivery_index=j_chronic,
        crosses_kahn_kalai_threshold=crosses_threshold,
        falsification_invariant_satisfied=inv_satisfied,
        complex_step_relative_error=rel_err,
    )


def evaluate_falsification_negative_control(
    poison_negative_shear: bool = False,
    poison_invalid_gsd: bool = False,
    poison_nonpositive_interval: bool = False,
    poison_excessive_single_dose_impulse: bool = False,
    verify_h1_trapped_and_cbf_overdose_controls: bool = False,
) -> Dict[str, Any]:
    """Tests poisoned negative controls (ValueError + AH-10 >250 nm H1 trapping and CBF h(x)<0 violation)."""
    reasons: list[str] = []
    if poison_negative_shear:
        try:
            compute_ehd_atomization_and_lnp_integrity(-500.0, 2.8, 1.15, 38.0)
            return {"passed": False, "reason": "Negative gamma_dot_s_inv accepted unexpectedly"}
        except ValueError:
            reasons.append("Deterministically rejected negative gamma_dot_s_inv")

    if poison_invalid_gsd:
        try:
            compute_ehd_atomization_and_lnp_integrity(4000.0, 2.8, 0.85, 38.0)
            return {"passed": False, "reason": "GSD < 1.0 accepted unexpectedly"}
        except ValueError:
            reasons.append("Deterministically rejected GSD < 1.0")

    if poison_nonpositive_interval:
        try:
            compute_cbf_redosing_invariance(0.20, 0.0)
            return {"passed": False, "reason": "Non-positive delta_t_h accepted unexpectedly"}
        except ValueError:
            reasons.append("Deterministically rejected non-positive delta_t_h")

    if poison_excessive_single_dose_impulse:
        try:
            compute_cbf_redosing_invariance(0.45, 72.0)
            return {"passed": False, "reason": "Single-dose impulse >= I_max accepted unexpectedly"}
        except ValueError:
            reasons.append("Deterministically rejected delta_i_dose >= I_max")

    if verify_h1_trapped_and_cbf_overdose_controls:
        trapped = evaluate_round204_ehd_taylor_cbf_redosing("unshielded_280nm_lnp_h1_mucus_trapped_fail")
        overdosed = evaluate_round204_ehd_taylor_cbf_redosing("naive_fixed_12h_overdosing_cbf_violation_fail")
        if trapped.crosses_kahn_kalai_threshold or trapped.h1_steric_pass or trapped.t_h1_mucus_transmission >= 1e-4:
            return {"passed": False, "reason": "Unshielded >250 nm LNP failed H1 pore trapping negative control"}
        if overdosed.crosses_kahn_kalai_threshold or overdosed.cbf_forward_invariant_safe or overdosed.h_cbf_margin >= 0.0:
            return {"passed": False, "reason": "Naive 12 h overdosing failed CBF h(x) < 0 violation negative control"}
        reasons.append("Verified >250 nm H1 mucus entrapment and 12 h overdosing CBF h(x) < 0 violation")

    return {"passed": True, "reason": "; ".join(reasons) if reasons else "Nominal baseline"}
