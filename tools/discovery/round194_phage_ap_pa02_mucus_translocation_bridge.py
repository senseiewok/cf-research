#!/usr/bin/env python3
"""
Round 194: Mucus-Penetrating Inhaled Bacteriophage AP-PA02 Translocation & Lytic Burst Invariant Bridge (T-0877).

Mathematical & Physical Pillar (Column A - Deductive Theorem):
  1. Aerosolized Bacteriophage Sputum Translocation Dynamics:
     AP-PA02 is an investigational cocktail of virulent obligately lytic bacteriophages targeting Pseudomonas aeruginosa,
     evaluated in cystic fibrosis chronic respiratory infections (SWARM-P.a. trial, NCT04596319).
     Aerosolized phage virions (representative capsid hydrodynamic diameter d_phage = 70.0 nm; mucosal phage subdiffusion grounded in
     Barr et al., Proc Natl Acad Sci USA 2015, DOI: 10.1073/pnas.1508355112) diffuse through dense,
     negatively charged CF mucus whose steric mesh pore size xi_mesh scales inversely with mucin solids weight fraction w_solid:
       xi_mesh = xi_0 * ( w_solid / w_0 )^-0.75
     where xi_0 = 140.0 nm and w_0 = 3.0% solids (Hill et al., Physiol Rev 2022, DOI: 10.1152/physrev.00004.2021;
     chronic respiratory P. aeruginosa phage lysis grounded in Waters et al., Thorax 2017, DOI: 10.1136/thoraxjnl-2016-209265).
     Steric hydrogel obstruction follows the Ogston-Amsden macromolecular pore permeability scaling:
       ratio = ( d_phage + d_fiber ) / ( xi_mesh + d_fiber )
       D_eff / D_0 = exp( - ( pi / 4 ) * ratio^2 * c_mesh )
     where d_fiber = 10.0 nm and c_mesh = 3.50 represents the multi-scale steric and viscoelastic hindrance factor.
     Effective phage virion concentration reaching bacterial biofilm microcolonies:
       [P]_eff = [P]_dose * ( D_eff / D_0 )

  2. Receptor Adsorption & Lytic Amplification Kinetics:
     Phages adsorb to outer membrane LPS core and Type IV pili receptors on P. aeruginosa. Across the therapeutic exposure interval,
     lytic infection cycles (burst size b ~ 100 virions/cell) reduce chronic bacterial density [B]_0 = 1.0e8 CFU/mL according to
     saturable multi-round lytic kinetics:
       Delta_log10_B = min( 4.0, max( 0.0, k_lytic * log10( 1.0 + [P]_eff / P_half ) ) )
       [B]_post = [B]_0 * 10^( -Delta_log10_B )
     where k_lytic = 1.85 and P_half = 1.0e7 PFU/mL.

  3. Bacterial Virulence Inactivation & CFTR Channel Functional Preservation:
     P. aeruginosa exotoxins (pyocyanin, elastase, alkaline protease) inactivate CFTR gating and cleave surface channels (PDB 8EIQ).
     Functional CFTR channel preservation follows cooperative Hill attenuation driven by post-treatment bacterial burden:
       f_active = 1.0 / [ 1.0 + ( [B]_post / B_crit )^n_H ]
     where B_crit = 1.5e6 CFU/mL and n_H = 1.50.
     When bacterial burden is suppressed ([B]_post <= 5.0e5 CFU/mL), f_active >= 0.85.

  4. Macroscopic Transepithelial Chloride Current & Airway Surface Liquid Hydration:
     In the presence of mature cell-surface CFTR (fraction f_Band_C restored by correctors) and potentiator ivacaftor:
       I_sc_equiv_pct_wt = f_active * f_Band_C * 100.0 * ( P_o / P_o_WT )
     where P_o / P_o_WT = 0.95.
     Airway surface liquid (ASL) height and mucociliary clearance (MCC) velocity:
       h_ASL = h_ASL_CF + ( h_ASL_WT - h_ASL_CF ) * [ I_sc_equiv_pct_wt / ( I_sc_equiv_pct_wt + K_05_ASL ) ]
       v_MCC = v_MCC_WT * max(0.0, ( h_ASL - h_cilia_min ) / ( h_ASL_WT - h_cilia_min ))^1.5

  5. Kahn-Kalai Curative Threshold for Bacteriophage AP-PA02 Therapy:
     - Bacterial burden reduction: Delta_log10_B >= 2.0 log10 CFU/mL.
     - Active uninhibited CFTR fraction: f_active >= 0.50 (50.0%).
     - Restored equivalent chloride current: I_sc_equiv_pct_wt >= 12.0% WT.
     - Steady-state ASL height: h_ASL >= 7.0 um.
     - Sustained mucociliary transport velocity: v_MCC >= 1.0 mm/min.
     - Mucus penetration efficiency: D_eff / D_0 >= 0.10.

  6. Invariants & Deterministic Bounds:
     - Strict bounds: 0.0 <= D_eff/D_0 <= 1.0; 0.0 <= f_active <= 1.0; 0.0 <= I_sc_equiv_pct_wt <= 100.0%.
     - Monotonicity: d(D_eff/D_0)/d(w_solid) < 0; d(Delta_log10_B)/d([P]_dose) >= 0; d(f_active)/d([P]_eff) >= 0.
     - Complex-step differentiation: d(h_ASL)/d(I_sc) verified at step h = 1e-30 against central finite differences.
     - Poisoned negative controls: negative phage dose, negative bacterial density, or negative solids fraction fail deterministically.

CF Biophysical Domain Mapping (Column B - Popperian Falsifiable Hypothesis):
  - In AP-PA02 High-Dose Trikafta Hydrated Regimen (NCT04596319; DOI: 10.1073/pnas.1508355112; DOI: 10.1136/thoraxjnl-2016-209265):
    [P]_dose = 1.0e9 PFU/mL, w_solid = 3.0%, f_Band_C = 0.35 (Trikafta):
    D_eff/D_0 = 0.458, [P]_eff = 4.58e8 PFU/mL, Delta_log10_B = 3.09 log10, [B]_post = 8.1e4 CFU/mL,
    f_active = 0.988, I_sc_equiv = 32.85% WT, h_ASL = 7.52 um >= 7.0 um, v_MCC = 2.44 mm/min >= 1.0 mm/min,
    robustly crossing Kahn-Kalai.
  - In AP-PA02 Standard-Dose Trikafta Regimen:
    [P]_dose = 5.0e8 PFU/mL, w_solid = 3.5%:
    D_eff/D_0 = 0.379, [P]_eff = 1.89e8 PFU/mL, Delta_log10_B = 2.41 log10, [B]_post = 3.9e5 CFU/mL,
    f_active = 0.883, I_sc_equiv = 29.36% WT, h_ASL = 7.43 um, v_MCC = 2.31 mm/min, crossing threshold.
  - In Untreated Severe Pseudomonas Biofilm Regimen:
    [P]_dose = 0, w_solid = 7.5%, f_Band_C = 0.06 (untreated F508del):
    Delta_log10_B = 0, [B]_post = 1.0e8 CFU/mL, f_active = 0.0018, I_sc_equiv = 0.01% WT,
    h_ASL = 3.51 um < 7.0 um, v_MCC = 0.0 mm/min, failing threshold.
  - In Phage Alone Without CFTR Corrector:
    [P]_dose = 1.0e9 PFU/mL, w_solid = 3.0%, but untreated F508del (f_Band_C = 0.06):
    Bacterial burden cleared, but low mature channel density yields I_sc = 5.63% WT < 12.0%,
    failing threshold; proves antimicrobial clearance alone is insufficient without modulator rescue.
  - In Dehydrated Mucus Phage Entrapment Regimen:
    [P]_dose = 2.5e8 PFU/mL, dehydrated w_solid = 8.5%:
    D_eff/D_0 = 0.041 < 0.10 (phage virions entrapped in hyper-concentrated mucus mesh),
    Delta_log10_B = 0.56 < 2.0 log10, [B]_post = 2.8e7 CFU/mL, f_active = 0.013,
    I_sc_equiv = 0.43% WT, failing threshold.

Research simulation only; not medical advice.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Dict, Tuple

PRIMARY_CITATIONS: Tuple[str, ...] = (
    "DOI: 10.1073/pnas.1508355112",
    "DOI: 10.1136/thoraxjnl-2016-209265",
    "DOI: 10.1152/physrev.00004.2021",
    "NCT04596319",
    "PDB: 8EIQ",
)

# Physical & Kinetic Constants
DEFAULT_D_PHAGE_NM: float = 70.0                # AP-PA02 hydrodynamic capsid diameter (nm)
DEFAULT_D_FIBER_NM: float = 10.0                # Mucin polymer fiber diameter (nm)
DEFAULT_XI_0_NM: float = 140.0                  # Baseline mucus mesh pore size at 3.0% solids (nm)
DEFAULT_W_0_PCT: float = 3.0                    # Reference mucus solids percentage (%)
DEFAULT_C_MESH: float = 3.50                    # Steric hindrance scaling exponent
DEFAULT_K_LYTIC: float = 1.85                   # Multi-round lytic burst log10 killing potency coefficient
DEFAULT_P_HALF_PFU_ML: float = 1.0e7            # Half-maximal effective phage concentration (PFU/mL)
DEFAULT_B_0_CFU_ML: float = 1.0e8               # Baseline chronic airway P. aeruginosa density (CFU/mL)
DEFAULT_B_CRIT_CFU_ML: float = 1.5e6            # Critical bacterial burden causing 50% CFTR functional inhibition (CFU/mL)
DEFAULT_HILL_N: float = 1.50                    # Hill coefficient for microbial virulence CFTR inactivation
DEFAULT_PO_RATIO: float = 0.95                  # Channel open probability ratio with ivacaftor (P_o / P_o_WT)
DEFAULT_H_ASL_CF_UM: float = 3.5                # Dehydrated CF airway surface liquid height (um)
DEFAULT_H_ASL_WT_UM: float = 8.5                # Healthy wild-type ASL height (um)
DEFAULT_K_05_ASL_PCT: float = 8.0               # Half-maximal ASL hydration chloride conductance (% WT)
DEFAULT_V_MCC_WT_MM_MIN: float = 4.0            # Healthy wild-type mucociliary transport velocity (mm/min)
DEFAULT_H_CILIA_MIN_UM: float = 5.0             # Minimum ASL height required for active ciliary stroke (um)

PARAMETER_PROVENANCE: Dict[str, Dict[str, object]] = {
    "ap_pa02_phage_diameter": {
        "value": 70.0,
        "tag": "derived",
        "source": "Representative caudoviricete capsid hydrodynamic diameter parameterization for mucosal phage diffusion (exact AP-PA02 virion diameter not stated in NCT04596319 registry record; mucosal phage subdiffusion grounded in Barr et al., Proc Natl Acad Sci USA 2015, DOI: 10.1073/pnas.1508355112)",
        "bounds": (50.0, 100.0),
    },
    "mucus_pore_scaling": {
        "value": 140.0,
        "tag": "derived",
        "source": "Reference CF sputum mesh spacing (~140 nm at 3% solids) with polymer gel correlation scaling xi ~ w^-0.75 (Hill et al., Physiol Rev 2022, DOI: 10.1152/physrev.00004.2021; Suk et al., Biomaterials 2009, DOI: 10.1016/j.biomaterials.2008.12.076)",
        "bounds": (80.0, 220.0),
    },
    "phage_lytic_potency": {
        "value": 1.85,
        "tag": "phenomenological_heuristic",
        "source": "Analytical multi-round lytic log10 reduction coefficient parameterized to chronic P. aeruginosa respiratory phage killing (Waters et al., Thorax 2017, DOI: 10.1136/thoraxjnl-2016-209265; clinical trial context NCT04596319, efficacy constant not stated in registry)",
        "bounds": (1.00, 2.50),
    },
    "bacterial_virulence_inactivation_threshold": {
        "value": 1.5e6,
        "tag": "phenomenological_heuristic",
        "source": "Analytical critical P. aeruginosa CFU/mL burden parameterization for exotoxin-driven CFTR inhibition (exact CFU/mL half-maximal constant not stated in Waters et al., Thorax 2017, DOI: 10.1136/thoraxjnl-2016-209265 or Barr et al., Proc Natl Acad Sci USA 2015, DOI: 10.1073/pnas.1508355112)",
        "bounds": (5.0e5, 5.0e6),
    },
    "cftr_gating_hill_cooperativity": {
        "value": 1.50,
        "tag": "phenomenological_heuristic",
        "source": "Analytical Hill cooperativity parameterization for bacterial exotoxin-mediated CFTR gating attenuation (structural reference PDB: 8EIQ)",
        "bounds": (1.0, 2.5),
    },
}

TWO_COLUMN_EPISTEMIC_HYGIENE: Dict[str, Any] = {
    "deductive_theorem": (
        "Inhaled bacteriophage cocktail AP-PA02 (NCT04596319) penetrates CF airway sputum with effective diffusivity ratio "
        "D_eff/D_0 determined by steric pore exclusion xi_mesh = xi_0 * (w_solid / w_0)^-0.75 and Ogston scaling. "
        "Effective phage concentration [P]_eff = [P]_dose * (D_eff/D_0) drives mass-action receptor binding and multi-round lytic burst cycles, "
        "reducing Pseudomonas aeruginosa density Delta_log10_B = k_lytic * log10(1 + [P]_eff / P_half). Suppressing bacterial burden "
        "below critical virulence threshold B_crit halts microbial CFTR inactivation (PDB 8EIQ), rescuing transepithelial chloride current, "
        "airway surface liquid height, and mucociliary clearance."
    ),
    "biophysical_hypothesis": {
        "statement": (
            "In Pseudomonas aeruginosa colonized CF airways ([B]_0 = 1.0e8 CFU/mL), inhaled AP-PA02 cocktail at "
            "[P]_dose >= 5.0e8 PFU/mL with mucus solids w_solid <= 3.5% achieves >= 2.0 log10 bacterial burden reduction, "
            "preserving active CFTR channel fraction f_active >= 0.75, transepithelial chloride current I_sc >= 12.0% WT, "
            "ASL height h_ASL >= 7.0 um, and mucociliary transport v_MCC >= 1.0 mm/min, crossing the Kahn-Kalai curative threshold."
        ),
        "falsification_condition": (
            "If nebulized AP-PA02 at [P]_dose >= 5.0e8 PFU/mL in sputum with w_solid <= 3.5% fails to achieve at least 1.8 log10 "
            "bacterial reduction or fails to maintain steady-state ASL height h_ASL >= 6.8 um in the presence of CFTR correctors, "
            "the AP-PA02 mucus translocation and lytic burst invariant bridge is falsified."
        ),
    },
}

FIVE_LAYER_LADDER: Dict[str, str] = {
    "cool_human_fact": (
        "Tiny natural virus hunters (phages) swim through airway mucus like microscopic torpedoes, latch onto Pseudomonas bacteria, "
        "and multiply inside them until the bacteria burst, saving the lung's chloride pumps from being poisoned."
    ),
    "everyday_analogy": (
        "Like deploying targeted robotic weed-cutters that weave through tall dense reeds to destroy invasive weeds, allowing the "
        "drainage pipes beneath to stay clear and carry fresh water."
    ),
    "engineering_lens": (
        "Steric pore-diffusion filtration kinetics in a viscoelastic polymer mesh coupled to nonlinear autocatalytic lytic "
        "predator-prey population dynamics."
    ),
    "formal_theorem": (
        "The effective bacterial log10 clearance functional Delta_log10_B([P]_dose, w_solid) is strictly monotonically increasing "
        "with phage dose and strictly monotonically decreasing with mucus solids content."
    ),
    "falsifiable_code_and_toy": (
        "A pure Python biophysical bridge with complex-step derivative parity at h=1e-30 and deterministic falsification on "
        "negative phage doses, unphysical mucus solids fractions, or negative bacterial densities."
    ),
}

REGIME_PRESETS: Dict[str, Dict[str, Any]] = {
    "ap_pa02_high_dose_trikafta_hydrated": {
        "phage_dose_pfu_ml": 1.0e9,
        "mucus_solids_pct": 3.0,
        "band_c_fraction": 0.35,
        "expected_crosses_kahn_kalai": True,
        "description": "1.0e9 PFU/mL AP-PA02 + 3.0% solids + Trikafta: 3.09 log10 kill, 98.8% active CFTR, crosses curative threshold",
    },
    "ap_pa02_standard_dose_trikafta": {
        "phage_dose_pfu_ml": 5.0e8,
        "mucus_solids_pct": 3.5,
        "band_c_fraction": 0.35,
        "expected_crosses_kahn_kalai": True,
        "description": "5.0e8 PFU/mL AP-PA02 + 3.5% solids + Trikafta: 2.41 log10 kill, 88.3% active CFTR, crosses curative threshold",
    },
    "untreated_severe_pseudomonas_biofilm": {
        "phage_dose_pfu_ml": 0.0,
        "mucus_solids_pct": 7.5,
        "band_c_fraction": 0.06,
        "expected_crosses_kahn_kalai": False,
        "description": "0 phage, 7.5% solids, untreated F508del: 0 kill, 0.2% active CFTR, 0.01% current, fails curative threshold",
    },
    "phage_alone_without_cftr_corrector": {
        "phage_dose_pfu_ml": 1.0e9,
        "mucus_solids_pct": 3.0,
        "band_c_fraction": 0.06,
        "expected_crosses_kahn_kalai": False,
        "description": "1.0e9 PFU/mL AP-PA02 clears bacteria, but untreated F508del yields I_sc=5.6% < 12%, fails curative threshold",
    },
    "dehydrated_mucus_phage_entrapment": {
        "phage_dose_pfu_ml": 2.5e8,
        "mucus_solids_pct": 8.5,
        "band_c_fraction": 0.35,
        "expected_crosses_kahn_kalai": False,
        "description": "2.5e8 PFU/mL AP-PA02 in 8.5% dehydrated mucus: virions entrapped (D/D0=0.041), kill only 0.56 log10, fails threshold",
    },
}


@dataclass(frozen=True)
class Round194PhageTranslocationResult:
    """Evaluation result for Round 194 AP-PA02 Phage Translocation & Lytic Burst Bridge."""

    preset_name: str
    phage_dose_pfu_ml: float
    mucus_solids_pct: float
    mucus_mesh_size_nm: float
    effective_diffusivity_ratio: float
    effective_phage_conc_pfu_ml: float
    bacterial_log10_kill: float
    post_treatment_bacteria_cfu_ml: float
    active_cftr_fraction: float
    transepithelial_current_pct_wt: float
    airway_surface_liquid_height_um: float
    mucociliary_transport_velocity_mm_min: float
    crosses_kahn_kalai_threshold: bool
    falsification_invariant_satisfied: bool
    complex_step_relative_error: float


def compute_mucus_mesh_size_nm(
    mucus_solids_pct: float,
    xi_0: float = DEFAULT_XI_0_NM,
    w_0: float = DEFAULT_W_0_PCT,
) -> float:
    """Compute steric pore mesh size of airway mucus as a function of solids percentage."""
    if mucus_solids_pct <= 0.0 or xi_0 <= 0.0 or w_0 <= 0.0:
        raise ValueError("Solids percentage and reference parameters must be strictly positive.")

    ratio = mucus_solids_pct / w_0
    return xi_0 * (ratio ** -0.75)


def compute_effective_diffusivity_ratio(
    mucus_mesh_size_nm: float,
    d_phage: float = DEFAULT_D_PHAGE_NM,
    d_fiber: float = DEFAULT_D_FIBER_NM,
    c_mesh: float = DEFAULT_C_MESH,
) -> float:
    """Compute effective diffusivity ratio D_eff / D_0 using the Ogston-Amsden steric obstruction model."""
    if mucus_mesh_size_nm <= 0.0 or d_phage <= 0.0 or d_fiber < 0.0 or c_mesh <= 0.0:
        raise ValueError("Mesh size, diameters, and mesh factor must be positive.")

    steric_ratio = (d_phage + d_fiber) / (mucus_mesh_size_nm + d_fiber)
    exponent = -(math.pi / 4.0) * (steric_ratio ** 2) * c_mesh
    return min(1.0, max(0.0, math.exp(exponent)))


def compute_bacterial_lytic_reduction(
    effective_phage_conc_pfu_ml: float,
    b_0: float = DEFAULT_B_0_CFU_ML,
    k_lytic: float = DEFAULT_K_LYTIC,
    p_half: float = DEFAULT_P_HALF_PFU_ML,
) -> Tuple[float, float]:
    """Compute log10 bacterial killing and surviving bacterial burden (CFU/mL)."""
    if effective_phage_conc_pfu_ml < 0.0 or b_0 <= 0.0 or k_lytic <= 0.0 or p_half <= 0.0:
        raise ValueError("Phage concentration and parameters must be non-negative with positive baseline and rates.")

    if effective_phage_conc_pfu_ml == 0.0:
        return 0.0, b_0

    dose_ratio = effective_phage_conc_pfu_ml / p_half
    log_reduction = min(4.0, max(0.0, k_lytic * math.log10(1.0 + dose_ratio)))
    post_burden = max(1.0, b_0 * (10.0 ** -log_reduction))
    return log_reduction, post_burden


def compute_active_cftr_fraction(
    post_treatment_bacteria_cfu_ml: float,
    b_crit: float = DEFAULT_B_CRIT_CFU_ML,
    hill_n: float = DEFAULT_HILL_N,
) -> float:
    """Compute preserved functional uninhibited CFTR channel fraction."""
    if post_treatment_bacteria_cfu_ml < 0.0 or b_crit <= 0.0 or hill_n <= 0.0:
        raise ValueError("Bacterial density must be non-negative with positive critical density and Hill exponent.")

    if post_treatment_bacteria_cfu_ml == 0.0:
        return 1.0

    ratio = post_treatment_bacteria_cfu_ml / b_crit
    power_term = ratio ** hill_n
    return min(1.0, max(0.0, 1.0 / (1.0 + power_term)))


def compute_transepithelial_current_pct_wt(
    active_cftr_fraction: float,
    band_c_fraction: float,
    po_ratio: float = DEFAULT_PO_RATIO,
) -> float:
    """Compute macroscopic transepithelial equivalent chloride current (% WT)."""
    if active_cftr_fraction < 0.0 or active_cftr_fraction > 1.0:
        raise ValueError("Active CFTR fraction must be in [0.0, 1.0].")
    if band_c_fraction < 0.0 or band_c_fraction > 1.0:
        raise ValueError("Mature Band C fraction must be in [0.0, 1.0].")
    if po_ratio <= 0.0 or po_ratio > 1.50:
        raise ValueError("Po ratio must be positive and <= 1.50.")

    return min(100.0, max(0.0, active_cftr_fraction * band_c_fraction * 100.0 * po_ratio))


def compute_asl_height_um(
    current_pct_wt: Any,
    h_asl_cf_um: float = DEFAULT_H_ASL_CF_UM,
    h_asl_wt_um: float = DEFAULT_H_ASL_WT_UM,
    k_05_asl_pct: float = DEFAULT_K_05_ASL_PCT,
) -> Any:
    """Compute steady-state airway surface liquid height (um), supporting complex numbers."""
    if isinstance(current_pct_wt, (int, float)) and current_pct_wt < 0.0:
        raise ValueError("Chloride current cannot be negative.")
    delta_h = h_asl_wt_um - h_asl_cf_um
    hydration_ratio = current_pct_wt / (current_pct_wt + k_05_asl_pct)
    return h_asl_cf_um + delta_h * hydration_ratio


def compute_mcc_velocity_mm_min(
    asl_height_um: float,
    h_cilia_min_um: float = DEFAULT_H_CILIA_MIN_UM,
    v_mcc_wt_mm_min: float = DEFAULT_V_MCC_WT_MM_MIN,
    h_asl_wt_um: float = DEFAULT_H_ASL_WT_UM,
) -> float:
    """Compute mucociliary clearance velocity (mm/min) above ciliary stroke floor."""
    if asl_height_um <= h_cilia_min_um:
        return 0.0
    effective_stroke = (asl_height_um - h_cilia_min_um) / (h_asl_wt_um - h_cilia_min_um)
    return v_mcc_wt_mm_min * (min(1.5, max(0.0, effective_stroke)) ** 1.5)


def complex_step_asl_derivative(
    current_pct_wt: float,
    h: float = 1e-30,
) -> Tuple[float, float, float]:
    """Verify sensitivity d(h_ASL)/d(I_sc) using complex-step differentiation against central finite differences."""
    curr_complex = complex(current_pct_wt, h)
    h_asl_complex = compute_asl_height_um(curr_complex)
    cs_derivative = h_asl_complex.imag / h

    delta = 1e-6
    curr_plus = min(100.0, current_pct_wt + delta)
    h_plus = compute_asl_height_um(curr_plus)

    curr_minus = max(0.0, current_pct_wt - delta)
    h_minus = compute_asl_height_um(curr_minus)

    fd_derivative = (h_plus - h_minus) / (2.0 * delta)
    denom = abs(fd_derivative) if abs(fd_derivative) > 1e-12 else 1.0
    rel_error = abs(cs_derivative - fd_derivative) / denom
    return cs_derivative, fd_derivative, rel_error


def evaluate_round194_phage_ap_pa02_mucus_translocation(preset_key: str) -> Round194PhageTranslocationResult:
    """Evaluate a named biophysical regime preset for Round 194."""
    if preset_key not in REGIME_PRESETS:
        raise KeyError(f"Unknown preset key: {preset_key}")

    p = REGIME_PRESETS[preset_key]
    dose = float(p["phage_dose_pfu_ml"])
    w_solids = float(p["mucus_solids_pct"])
    band_c = float(p["band_c_fraction"])

    mesh_nm = compute_mucus_mesh_size_nm(w_solids)
    d_ratio = compute_effective_diffusivity_ratio(mesh_nm)
    p_eff = dose * d_ratio
    log_kill, b_post = compute_bacterial_lytic_reduction(p_eff)
    f_act = compute_active_cftr_fraction(b_post)
    curr_pct = compute_transepithelial_current_pct_wt(f_act, band_c)
    h_asl = compute_asl_height_um(curr_pct)
    v_mcc = compute_mcc_velocity_mm_min(h_asl)

    _, _, rel_err = complex_step_asl_derivative(max(0.1, curr_pct))

    # Kahn-Kalai Curative Threshold for AP-PA02 Bacteriophage Therapy:
    # Delta_log10_B >= 2.0 log10, f_active >= 0.50, I_sc >= 12.0% WT, h_ASL >= 7.0 um, v_MCC >= 1.0 mm/min, D_eff/D_0 >= 0.10
    crosses_threshold = (
        (log_kill >= 2.0)
        and (f_act >= 0.50)
        and (curr_pct >= 12.0)
        and (h_asl >= 7.0)
        and (v_mcc >= 1.0)
        and (d_ratio >= 0.10)
    )
    expected_cross = bool(p["expected_crosses_kahn_kalai"])
    invariant_satisfied = (crosses_threshold == expected_cross) and (rel_err <= 1e-6)

    return Round194PhageTranslocationResult(
        preset_name=preset_key,
        phage_dose_pfu_ml=dose,
        mucus_solids_pct=w_solids,
        mucus_mesh_size_nm=mesh_nm,
        effective_diffusivity_ratio=d_ratio,
        effective_phage_conc_pfu_ml=p_eff,
        bacterial_log10_kill=log_kill,
        post_treatment_bacteria_cfu_ml=b_post,
        active_cftr_fraction=f_act,
        transepithelial_current_pct_wt=curr_pct,
        airway_surface_liquid_height_um=h_asl,
        mucociliary_transport_velocity_mm_min=v_mcc,
        crosses_kahn_kalai_threshold=crosses_threshold,
        falsification_invariant_satisfied=invariant_satisfied,
        complex_step_relative_error=rel_err,
    )


def evaluate_falsification_negative_control(
    poison_negative_dose: bool = False,
    poison_negative_solids: bool = False,
    poison_negative_bacteria: bool = False,
) -> Dict[str, Any]:
    """Execute poisoned negative controls to verify deterministic falsification."""
    if poison_negative_dose:
        try:
            compute_bacterial_lytic_reduction(effective_phage_conc_pfu_ml=-1.0e8)
            return {"passed": False, "error": "Failed to catch negative phage concentration"}
        except ValueError as e:
            return {"passed": True, "error": f"Correctly caught: {e}"}

    if poison_negative_solids:
        try:
            compute_mucus_mesh_size_nm(mucus_solids_pct=-2.0)
            return {"passed": False, "error": "Failed to catch negative mucus solids"}
        except ValueError as e:
            return {"passed": True, "error": f"Correctly caught: {e}"}

    if poison_negative_bacteria:
        try:
            compute_active_cftr_fraction(post_treatment_bacteria_cfu_ml=-5.0e6)
            return {"passed": False, "error": "Failed to catch negative bacterial density"}
        except ValueError as e:
            return {"passed": True, "error": f"Correctly caught: {e}"}

    return {"passed": True, "error": "No negative controls triggered"}


if __name__ == "__main__":
    print("=== Round 194: Bacteriophage AP-PA02 Translocation & Lytic Burst Bridge ===")
    for key in REGIME_PRESETS:
        res = evaluate_round194_phage_ap_pa02_mucus_translocation(key)
        status = "PASS" if res.crosses_kahn_kalai_threshold else "FAIL (Sub-curative)"
        print(f"[{status}] {res.preset_name}: D/D0={res.effective_diffusivity_ratio:.3f}, Kill={res.bacterial_log10_kill:.2f} log10, Post-B={res.post_treatment_bacteria_cfu_ml:.1e} CFU/mL, Intact={res.active_cftr_fraction*100:.1f}%, I_sc={res.transepithelial_current_pct_wt:.2f}% WT, ASL={res.airway_surface_liquid_height_um:.2f} um, MCC={res.mucociliary_transport_velocity_mm_min:.2f} mm/min")
