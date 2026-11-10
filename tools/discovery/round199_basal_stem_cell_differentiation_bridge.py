#!/usr/bin/env python3
"""
Round 199: Basal Stem Cell Notch/Wnt Signaling & Multiciliated Cell Lineage Differentiation Invariant Bridge (T-0888).

Mathematical & Physical Pillar (Column A - Deductive Theorem):
  1. Three-State Epithelial Lineage Partition Simplex (Notch/Wnt Signaling Balance):
     Adult human conducting airway KRT5+/TP63+ basal stem cells self-renew and differentiate along competing
     mucociliary trajectories governed by Notch intracellular domain activity (S_notch >= 0) and canonical
     Wnt/beta-catenin + FOXJ1 ciliogenic drive (S_wnt >= 0) under mucosal inflammatory bias (I_inf >= 0).
     Sustained Notch signaling promotes goblet/club mucin-lineage commitment while suppressing multiciliated cell
     differentiation; conversely, Notch attenuation paired with Wnt/FOXJ1 activation drives deuterosomal
     multiciliated specification (Rock et al., Cell Stem Cell 2011, DOI: 10.1016/j.stem.2011.04.003;
     Ruiz Garcia et al., Development 2019, DOI: 10.1242/dev.177428).
     On the stationary epithelial lineage simplex with partition weights:
       W_basal = w_b0 = 0.28
       W_cil   = (alpha_cil * S_wnt) / (1.0 + (S_notch / K_notch)**2 + gamma_inf * I_inf)
       W_gob   = w_g0 + alpha_gob * (S_notch**2 / (K_gob**2 + S_notch**2)) * (1.0 + beta_inf * I_inf)
     where alpha_cil = 1.60, K_notch = 0.45, gamma_inf = 0.80, w_g0 = 0.12, alpha_gob = 1.35, K_gob = 0.50,
     and beta_inf = 0.90, the partition sum Z = W_basal + W_cil + W_gob > 0 yields exact simplex fractions:
       f_basal = W_basal / Z,   f_ciliated = W_cil / Z,   f_goblet = W_gob / Z,
     satisfying f_basal + f_ciliated + f_goblet = 1.0 identically.

  2. Clonally Inherited Surface CFTR Conductance & Hydrodynamic Ciliary Percolation:
     Given an edited/corrected basal stem cell fraction eta_edit in [0, 1], progeny populate differentiated
     luminal surface cells while goblet-cell hyperplasia dilutes effective apical ion-transporting surface area:
       f_cftr_surface = eta_edit * (1.0 - f_basal) * (f_ciliated / (f_ciliated + K_arch * f_goblet))
       I_sc_equiv_pct_wt = 100.0 * (f_cftr_surface / (f_cftr_surface + K_05_cftr))
     where K_arch = 0.65 and K_05_cftr = 0.15.
     Long-range mucociliary transport velocity v_mcc (mm/min) requires multiciliated coverage exceeding the
     hydrodynamic metachronal percolation threshold phi_crit = 0.32 (Loiseau et al., Nat Phys 2020,
     DOI: 10.1038/s41567-019-0770-2):
       Phi_perc = f_ciliated**4 / (phi_crit**4 + f_ciliated**4)
       v_mcc = v_max * Phi_perc * (I_sc_equiv_pct_wt / (I_sc_equiv_pct_wt + K_05_isc)) / (1.0 + lambda_gob * f_goblet**2)
     where v_max = 4.2 mm/min, K_05_isc = 18.0 (% WT), and lambda_gob = 2.2.

  3. Kahn-Kalai Curative Threshold for Basal Lineage Regeneration & Ciliary Percolation:
     - Multiciliated surface fraction f_ciliated >= 0.40 (>= 40.0% coverage, exceeding phi_crit = 0.32).
     - Goblet/mucin-lineage fraction f_goblet <= 0.30 (<= 30.0%, preventing goblet hyperplasia).
     - Restored transepithelial chloride current I_sc_equiv_pct_wt >= 25.0% WT.
     - Coordinated mucociliary transport velocity v_mcc >= 1.50 mm/min.

  4. Invariants & Deterministic Bounds:
     - Simplex conservation: 0.0 <= f_basal, f_ciliated, f_goblet <= 1.0 and |f_basal + f_ciliated + f_goblet - 1.0| <= 1e-12.
     - Directional signaling derivatives: d(f_ciliated)/d(S_notch) < 0; d(f_goblet)/d(S_notch) > 0; d(f_ciliated)/d(S_wnt) > 0.
     - Complex-step differentiation: d(v_mcc)/d(S_wnt) verified at step h = 1e-30 against central finite differences
       with relative error <= 1e-6.
     - Poisoned negative controls: negative Notch/Wnt signaling or out-of-bounds eta_edit fail deterministically with ValueError.

CF Biophysical Domain Mapping (Column B - Popperian Falsifiable Hypothesis):
  - In Notch-Modulated Wnt-Driven Edited Basal Curative Regimen (S_notch=0.25, S_wnt=1.15, I_inf=0.10, eta_edit=0.55):
    f_basal = 0.1387, f_ciliated = 0.6562 >= 0.40, f_goblet = 0.2052 <= 0.30, f_cftr_surface = 0.3937,
    I_sc_equiv = 72.41% WT >= 25.0%, Phi_perc = 0.9465, v_mcc = 2.91 mm/min >= 1.50 mm/min (Pass).
  - In High-Efficiency Basal Engraftment Curative Regimen (S_notch=0.20, S_wnt=1.30, I_inf=0.05, eta_edit=0.78):
    f_basal = 0.1231, f_ciliated = 0.7387, f_goblet = 0.1383, f_cftr_surface = 0.6098,
    I_sc_equiv = 80.26% WT, Phi_perc = 0.9660, v_mcc = 3.18 mm/min (Pass).
  - In Chronic Notch Hyperactivation Goblet Hyperplasia CF Regimen (S_notch=1.25, S_wnt=0.55, I_inf=1.40, eta_edit=0.45):
    f_basal = 0.0898, f_ciliated = 0.0287 < 0.40, f_goblet = 0.8816 > 0.30, I_sc_equiv = 11.52% WT,
    Phi_perc = 0.0001, v_mcc = 0.00 mm/min, failing due to Notch-driven goblet hyperplasia and ciliary loss.
  - In Unedited CF Basal Progenitor Negative Control (S_notch=0.30, S_wnt=1.05, I_inf=0.20, eta_edit=0.00):
    f_ciliated = 0.5603, f_goblet = 0.2899, f_cftr_surface = 0.0000, I_sc_equiv = 0.00% WT, v_mcc = 0.00 mm/min (Fail).
  - In Wnt-Deficient Ciliogenesis Arrest Regimen (S_notch=0.35, S_wnt=0.12, I_inf=0.60, eta_edit=0.60):
    f_basal = 0.2381, f_ciliated = 0.0783 < 0.40, f_goblet = 0.6835 > 0.30, I_sc_equiv = 31.35% WT,
    Phi_perc = 0.0036, v_mcc = 0.00 mm/min, failing hydrodynamic percolation despite 60% basal editing.

Research simulation only; not medical advice.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Tuple

PRIMARY_CITATIONS: Tuple[str, ...] = (
    "DOI: 10.1016/j.stem.2011.04.003",
    "DOI: 10.1242/dev.177428",
    "DOI: 10.1038/s41567-019-0770-2",
)

# Lineage Simplex & Hydrodynamic Percolation Constants
DEFAULT_W_BASAL_BASE: float = 0.28
DEFAULT_ALPHA_CIL: float = 1.60
DEFAULT_K_NOTCH: float = 0.45
DEFAULT_GAMMA_INF: float = 0.80
DEFAULT_W_GOB_BASE: float = 0.12
DEFAULT_ALPHA_GOB: float = 1.35
DEFAULT_K_GOB: float = 0.50
DEFAULT_BETA_INF: float = 0.90
DEFAULT_K_ARCH: float = 0.65
DEFAULT_K_05_CFTR: float = 0.15
DEFAULT_PHI_CRIT_PERC: float = 0.32
DEFAULT_V_MAX_MM_MIN: float = 4.2
DEFAULT_K_05_ISC: float = 18.0
DEFAULT_LAMBDA_GOB: float = 2.2

PARAMETER_PROVENANCE: Dict[str, Dict[str, object]] = {
    "notch_ciliated_repression_constant": {
        "value": DEFAULT_K_NOTCH,
        "tag": "phenomenological_heuristic",
        "source": (
            "Analytical half-repression Notch signaling constant K_notch modeling Notch lateral inhibition of "
            "multiciliated cell fate in adult airway basal progenitors (exact 0.45 dimensionless scalar not stated "
            "in Rock et al., Cell Stem Cell 2011, DOI: 10.1016/j.stem.2011.04.003)"
        ),
        "bounds": (0.20, 0.90),
    },
    "wnt_foxj1_ciliogenic_gain": {
        "value": DEFAULT_ALPHA_CIL,
        "tag": "phenomenological_heuristic",
        "source": (
            "Analytical ciliogenic partition gain alpha_cil for deuterosomal multiciliated differentiation trajectories "
            "(exact 1.60 dimensionless weight not stated in Ruiz Garcia et al., Development 2019, DOI: 10.1242/dev.177428)"
        ),
        "bounds": (0.80, 2.80),
    },
    "notch_goblet_commitment_gain": {
        "value": DEFAULT_ALPHA_GOB,
        "tag": "phenomenological_heuristic",
        "source": (
            "Analytical goblet/club lineage partition weight alpha_gob under sustained Notch activation "
            "(exact 1.35 dimensionless scalar not stated in Rock et al., Cell Stem Cell 2011, DOI: 10.1016/j.stem.2011.04.003)"
        ),
        "bounds": (0.60, 2.50),
    },
    "ciliary_percolation_threshold": {
        "value": DEFAULT_PHI_CRIT_PERC,
        "tag": "phenomenological_heuristic",
        "source": (
            "Analytical Hill percolation midpoint phi_crit parameterizing ciliated cell coverage required for long-range "
            "mucociliary transport (exact 0.32 Hill midpoint scalar not stated in Loiseau et al., Nat Phys 2020, "
            "DOI: 10.1038/s41567-019-0770-2)"
        ),
        "bounds": (0.20, 0.50),
    },
    "cftr_surface_conductance_half_saturation": {
        "value": DEFAULT_K_05_CFTR,
        "tag": "phenomenological_heuristic",
        "source": (
            "Analytical half-saturation constant K_05_cftr relating differentiated corrected luminal CFTR fraction to "
            "apical short-circuit current % WT (exact 0.15 scalar not stated in Ruiz Garcia et al., Development 2019, "
            "DOI: 10.1242/dev.177428)"
        ),
        "bounds": (0.08, 0.25),
    },
}

TWO_COLUMN_EPISTEMIC_HYGIENE: Dict[str, Dict[str, str]] = {
    "deductive_column_a": {
        "formal_theorem": (
            "On a positive three-state Boltzmann-Hill partition simplex Z = W_basal + W_cil(S_wnt, S_notch) + "
            "W_gob(S_notch), the stationary lineage fractions f_i = W_i / Z satisfy exact simplex conservation "
            "sum_i f_i = 1, strict interior positivity for finite positive weights, d(f_ciliated)/d(S_notch) < 0, "
            "and d(f_ciliated)/d(S_wnt) > 0, coupling clonally inherited CFTR expression to fourth-order hydrodynamic "
            "percolation of mucociliary clearance."
        ),
        "invariants": (
            "0 <= f_basal, f_ciliated, f_goblet <= 1; |f_basal + f_ciliated + f_goblet - 1| <= 1e-12; "
            "0 <= I_sc_equiv_pct_wt <= 100%; d(f_ciliated)/d(S_notch) < 0; d(f_goblet)/d(S_notch) > 0; "
            "d(f_ciliated)/d(S_wnt) > 0."
        ),
    },
    "popperian_column_b": {
        "falsifiable_hypothesis": (
            "Gene correction of >= 50% of airway KRT5+ basal stem cells restores long-range mucociliary clearance "
            "(v_mcc >= 1.5 mm/min) only when Notch/Wnt signaling maintains multiciliated surface coverage >= 40% "
            "and restricts goblet-cell hyperplasia to <= 30%; sustained Notch hyperactivation abolishes mucociliary "
            "transport despite high progenitor editing."
        ),
        "empirical_falsification_threshold": (
            "Falsified if Notch-hyperactivated goblet-hyperplastic epithelia (f_ciliated < 10%, f_goblet > 75%) "
            "exhibit coordinated long-range mucociliary clearance >= 1.5 mm/min without ciliated coverage recovery, "
            "or if Notch attenuation fails to increase multiciliated lineage commitment in human airway basal cultures."
        ),
    },
}

FIVE_LAYER_LADDER: Dict[str, str] = {
    "cool_human_fact": (
        "Editing the CFTR gene inside airway basal stem cells is only half the victory: those stem cells still use "
        "two ancient developmental dials—Notch and Wnt—to decide whether their daughters become microscopic broom "
        "cells with 200 beating cilia or mucus-making goblet cells."
    ),
    "everyday_analogy": (
        "Think of repairing the engines in a city's bus depot: even if half the fleet has brand-new clean engines, "
        "traffic still stalls if the dispatch switch keeps sending out gravel trucks instead of passenger buses "
        "or leaves large stretches of highway disconnected."
    ),
    "engineering_lens": (
        "Airway epithelial regeneration couples a three-state lineage partition system to a 2D hydrodynamic "
        "percolation network: Notch lateral inhibition shifts probability mass from FOXJ1+ multiciliated nodes "
        "to goblet nodes, dropping ciliated coverage below the percolation threshold phi_crit."
    ),
    "formal_theorem": (
        "For partition weights W_basal > 0, W_cil(S_wnt, S_notch), and W_gob(S_notch), the normalized simplex "
        "f = W / (W_basal + W_cil + W_gob) preserves sum_i f_i = 1 everywhere in C^n and yields analytic "
        "complex-step sensitivities coupling progenitor editing fraction eta_edit to metachronal transport v_mcc."
    ),
    "falsifiable_code_and_toy": (
        "evaluate_round199_basal_stem_cell_differentiation(regime_key) evaluates lineage simplex fractions, "
        "clonally inherited surface CFTR current, and hydrodynamic ciliary percolation across 5 airway regimens."
    ),
}

REGIME_PRESETS: Dict[str, Dict[str, float]] = {
    "notch_modulated_wnt_edited_basal_curative": {
        "s_notch": 0.25,
        "s_wnt": 1.15,
        "i_inflammation": 0.10,
        "eta_basal_edit": 0.55,
    },
    "high_efficiency_basal_engraftment_curative": {
        "s_notch": 0.20,
        "s_wnt": 1.30,
        "i_inflammation": 0.05,
        "eta_basal_edit": 0.78,
    },
    "chronic_notch_goblet_hyperplasia_cf": {
        "s_notch": 1.25,
        "s_wnt": 0.55,
        "i_inflammation": 1.40,
        "eta_basal_edit": 0.45,
    },
    "unedited_cf_basal_progenitors": {
        "s_notch": 0.30,
        "s_wnt": 1.05,
        "i_inflammation": 0.20,
        "eta_basal_edit": 0.00,
    },
    "wnt_deficient_squamous_basal_arrest": {
        "s_notch": 0.35,
        "s_wnt": 0.12,
        "i_inflammation": 0.60,
        "eta_basal_edit": 0.60,
    },
}


@dataclass(frozen=True)
class Round199BasalDifferentiationResult:
    regime_key: str
    s_notch: float
    s_wnt: float
    i_inflammation: float
    eta_basal_edit: float
    f_basal: float
    f_ciliated: float
    f_goblet: float
    f_cftr_surface: float
    i_sc_equiv_pct_wt: float
    percolation_factor: float
    v_mcc_mm_min: float
    crosses_kahn_kalai_threshold: bool
    falsification_invariant_satisfied: bool
    complex_step_relative_error: float
    research_boundary: str = "Research simulation only; not medical advice."


def compute_lineage_simplex_fractions(
    s_notch: float | complex,
    s_wnt: float | complex,
    i_inflammation: float | complex = 0.0,
) -> Tuple[float | complex, float | complex, float | complex]:
    """Computes normalized stationary fractions (f_basal, f_ciliated, f_goblet) on the airway lineage simplex."""
    for name, val in (("s_notch", s_notch), ("s_wnt", s_wnt), ("i_inflammation", i_inflammation)):
        if isinstance(val, (int, float)) and val < 0.0:
            raise ValueError(f"{name} must be non-negative, got {val}")

    w_basal = DEFAULT_W_BASAL_BASE
    w_cil = (DEFAULT_ALPHA_CIL * s_wnt) / (
        1.0 + (s_notch / DEFAULT_K_NOTCH) ** 2 + DEFAULT_GAMMA_INF * i_inflammation
    )
    w_gob = DEFAULT_W_GOB_BASE + DEFAULT_ALPHA_GOB * (
        (s_notch ** 2) / (DEFAULT_K_GOB ** 2 + s_notch ** 2)
    ) * (1.0 + DEFAULT_BETA_INF * i_inflammation)

    z_sum = w_basal + w_cil + w_gob
    return w_basal / z_sum, w_cil / z_sum, w_gob / z_sum


def compute_surface_cftr_and_conductance(
    eta_basal_edit: float | complex,
    f_basal: float | complex,
    f_ciliated: float | complex,
    f_goblet: float | complex,
) -> Tuple[float | complex, float | complex]:
    """Computes differentiated surface CFTR fraction and restored equivalent chloride short-circuit current (% WT)."""
    if isinstance(eta_basal_edit, (int, float)) and (eta_basal_edit < 0.0 or eta_basal_edit > 1.0):
        raise ValueError(f"eta_basal_edit must be in [0, 1], got {eta_basal_edit}")
    for name, val in (("f_basal", f_basal), ("f_ciliated", f_ciliated), ("f_goblet", f_goblet)):
        if isinstance(val, (int, float)) and val < 0.0:
            raise ValueError(f"{name} must be non-negative, got {val}")

    arch_denom = f_ciliated + DEFAULT_K_ARCH * f_goblet + 1e-15
    f_cftr_surface = eta_basal_edit * (1.0 - f_basal) * (f_ciliated / arch_denom)
    i_sc_pct_wt = 100.0 * (f_cftr_surface / (f_cftr_surface + DEFAULT_K_05_CFTR))
    return f_cftr_surface, i_sc_pct_wt


def compute_percolated_mcc_velocity(
    f_ciliated: float | complex,
    f_goblet: float | complex,
    i_sc_equiv_pct_wt: float | complex,
) -> Tuple[float | complex, float | complex]:
    """Computes fourth-order hydrodynamic ciliary percolation factor and mucociliary transport velocity (mm/min)."""
    for name, val in (("f_ciliated", f_ciliated), ("f_goblet", f_goblet), ("i_sc_equiv_pct_wt", i_sc_equiv_pct_wt)):
        if isinstance(val, (int, float)) and val < 0.0:
            raise ValueError(f"{name} must be non-negative, got {val}")

    phi_perc = (f_ciliated ** 4) / (DEFAULT_PHI_CRIT_PERC ** 4 + f_ciliated ** 4)
    hydration_factor = i_sc_equiv_pct_wt / (i_sc_equiv_pct_wt + DEFAULT_K_05_ISC)
    goblet_drag = 1.0 + DEFAULT_LAMBDA_GOB * (f_goblet ** 2)
    v_mcc = DEFAULT_V_MAX_MM_MIN * phi_perc * hydration_factor / goblet_drag
    return phi_perc, v_mcc


def complex_step_wnt_mcc_derivative(
    s_notch: float,
    s_wnt: float,
    i_inflammation: float,
    eta_basal_edit: float,
    step_h: float = 1e-30,
) -> Tuple[float, float, float]:
    """Verifies d(v_mcc)/d(S_wnt) via complex-step differentiation at h = 1e-30 against central finite differences."""
    def _eval_vmcc(wnt_val: float | complex) -> float | complex:
        fb, fc, fg = compute_lineage_simplex_fractions(s_notch, wnt_val, i_inflammation)
        _, isc = compute_surface_cftr_and_conductance(eta_basal_edit, fb, fc, fg)
        _, vm = compute_percolated_mcc_velocity(fc, fg, isc)
        return vm

    cs_val = _eval_vmcc(complex(s_wnt, step_h))
    cs_deriv = cs_val.imag / step_h

    fd_h = 1e-6
    fd_deriv = float(_eval_vmcc(s_wnt + fd_h) - _eval_vmcc(s_wnt - fd_h)) / (2.0 * fd_h)
    rel_err = abs(cs_deriv - fd_deriv) / max(abs(cs_deriv), 1e-12)
    return cs_deriv, fd_deriv, rel_err


def evaluate_round199_basal_stem_cell_differentiation(regime_key: str) -> Round199BasalDifferentiationResult:
    """Evaluates the Round 199 basal stem cell Notch/Wnt differentiation and ciliary percolation bridge."""
    if regime_key not in REGIME_PRESETS:
        raise KeyError(f"Unknown regime '{regime_key}'. Valid regimes: {list(REGIME_PRESETS.keys())}")

    cfg = REGIME_PRESETS[regime_key]
    s_notch = cfg["s_notch"]
    s_wnt = cfg["s_wnt"]
    i_inf = cfg["i_inflammation"]
    eta_edit = cfg["eta_basal_edit"]

    fb_r, fc_r, fg_r = compute_lineage_simplex_fractions(s_notch, s_wnt, i_inf)
    f_basal = float(fb_r)
    f_ciliated = float(fc_r)
    f_goblet = float(fg_r)

    fcf_r, isc_r = compute_surface_cftr_and_conductance(eta_edit, f_basal, f_ciliated, f_goblet)
    f_cftr_surface = float(fcf_r)
    i_sc_pct_wt = float(isc_r)

    phi_r, vmcc_r = compute_percolated_mcc_velocity(f_ciliated, f_goblet, i_sc_pct_wt)
    phi_perc = float(phi_r)
    v_mcc = float(vmcc_r)

    eval_edit = max(eta_edit, 0.25)
    eval_wnt = max(s_wnt, 0.20)
    _, _, rel_err = complex_step_wnt_mcc_derivative(s_notch, eval_wnt, i_inf, eval_edit)

    simplex_sum = f_basal + f_ciliated + f_goblet
    inv_satisfied = (
        0.0 <= f_basal <= 1.0
        and 0.0 <= f_ciliated <= 1.0
        and 0.0 <= f_goblet <= 1.0
        and abs(simplex_sum - 1.0) <= 1e-12
        and 0.0 <= f_cftr_surface <= 1.0
        and 0.0 <= i_sc_pct_wt <= 100.0
        and 0.0 <= phi_perc <= 1.0
        and 0.0 <= v_mcc <= DEFAULT_V_MAX_MM_MIN
        and rel_err <= 1e-6
    )

    crosses_threshold = (
        f_ciliated >= 0.40
        and f_goblet <= 0.30
        and i_sc_pct_wt >= 25.0
        and v_mcc >= 1.50
    )

    return Round199BasalDifferentiationResult(
        regime_key=regime_key,
        s_notch=s_notch,
        s_wnt=s_wnt,
        i_inflammation=i_inf,
        eta_basal_edit=eta_edit,
        f_basal=f_basal,
        f_ciliated=f_ciliated,
        f_goblet=f_goblet,
        f_cftr_surface=f_cftr_surface,
        i_sc_equiv_pct_wt=i_sc_pct_wt,
        percolation_factor=phi_perc,
        v_mcc_mm_min=v_mcc,
        crosses_kahn_kalai_threshold=crosses_threshold,
        falsification_invariant_satisfied=inv_satisfied,
        complex_step_relative_error=rel_err,
    )


def evaluate_falsification_negative_control(
    poison_negative_notch: bool = False,
    poison_negative_wnt: bool = False,
    poison_invalid_edit_fraction: bool = False,
) -> Dict[str, Any]:
    """Tests poisoned negative controls that must deterministically raise ValueError."""
    reasons: list[str] = []
    if poison_negative_notch:
        try:
            compute_lineage_simplex_fractions(-0.2, 1.0, 0.0)
            return {"passed": False, "reason": "Negative Notch signaling accepted unexpectedly"}
        except ValueError:
            reasons.append("Deterministically rejected negative Notch signaling")

    if poison_negative_wnt:
        try:
            compute_lineage_simplex_fractions(0.3, -0.5, 0.0)
            return {"passed": False, "reason": "Negative Wnt signaling accepted unexpectedly"}
        except ValueError:
            reasons.append("Deterministically rejected negative Wnt signaling")

    if poison_invalid_edit_fraction:
        try:
            compute_surface_cftr_and_conductance(1.25, 0.2, 0.6, 0.2)
            return {"passed": False, "reason": "Edit fraction > 1.0 accepted unexpectedly"}
        except ValueError:
            reasons.append("Deterministically rejected out-of-bounds edit fraction")

    return {"passed": True, "reason": "; ".join(reasons) if reasons else "Nominal baseline"}
