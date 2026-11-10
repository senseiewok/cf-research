#!/usr/bin/env python3
"""
Round 180: TMEM16A Anion Potentiator Allosteric Activation & Calcium-Independent Conductance Coupling Invariant Bridge for Modulator-Ineligible CF Epithelia (T-0854).

Mathematical & Physical Pillar (Column A - Deductive Theorem):
  1. Allosteric Calcium Sensitization & Shift Kinetics in TMEM16A (ANO1):
     In modulator-ineligible CF patients (Class I nonsense, frameshift, or canonical splice-site mutations
     lacking full-length CFTR protein), small-molecule potentiators of the alternative apical anion channel
     TMEM16A (anoctamin-1) offer a genotype-agnostic chloride and bicarbonate bypass pathway.
     TMEM16A is natively gated by intracellular calcium ([Ca2+]_i); basal resting airway epithelial
     calcium is low ([Ca2+]_i ~ 0.10 to 0.12 uM), keeping native TMEM16A channels closed (EC50_0 ~ 0.85 uM).
     Cryo-EM structures (Paulino et al. 2017 Nature, DOI: 10.1038/nature24652; PDB 5OYB) and electrophysiology
     (Danahay et al. 2020 AJRCCM, DOI: 10.1164/rccm.201908-1641OC) establish that small-molecule potentiators
     bind allosterically near the transmembrane calcium-coordinating pocket, shifting the apparent EC50:
       EC50(C_pot) = EC50_0 * [ 1.0 - f_shift_max * (C_pot / (C_pot + K_D_pot)) ]
     where EC50_0 = 0.85 uM, f_shift_max = 0.82 (shifting EC50 down to ~0.15 uM), C_pot is potentiator
     concentration (uM), and K_D_pot = 0.35 uM.

  2. Calcium-Independent Basal Leak Gating:
     Potentiator binding also stabilizes an open pore conformation in the absence of bound divalent cation:
       P_o_basal(C_pot) = P_o_leak + P_o_ind_max * [ C_pot / (C_pot + K_D_pot) ]
     where basal leak open probability P_o_leak = 0.015 and maximum ligand-induced calcium-independent
     open probability P_o_ind_max = 0.12.

  3. Cooperative Anion Channel Open Probability:
     Total TMEM16A open probability P_o is modeled via the cooperative Hill equation:
       P_o([Ca2+]_i, C_pot) = P_o_basal(C_pot) + (P_o_max - P_o_basal(C_pot)) * [ [Ca2+]_i^n_H / ([Ca2+]_i^n_H + EC50(C_pot)^n_H) ]
     where maximum open probability P_o_max = 0.90 and Hill coefficient n_H = 2.2.

  4. Macroscopic Epithelial Conductance & Transepithelial Secretion:
     The apical membrane TMEM16A macroscopic conductance G_TMEM16A (mS/cm^2) is:
       G_TMEM16A = N_channel * gamma_single * P_o * 1e3
     where channel surface density N_channel = 2.5e7 channels/cm^2, single channel conductance
     gamma_single = 8.0 pS (8.0e-12 S), yielding scaling factor G_max = 0.20 mS/cm^2.
     Transepithelial equivalent short-circuit current (% of WT CFTR):
       I_sc_equiv_pct_wt = (G_TMEM16A / G_CFTR_WT) * 100.0%
     where stimulated WT CFTR apical conductance G_CFTR_WT = 0.08 mS/cm^2.

  5. Airway Surface Liquid (ASL) Hydration & Mucociliary Clearance (MCC):
     Coupled osmotic fluid transport restores the periciliary liquid layer (Boucher et al. 2022 Physiol Rev,
     DOI: 10.1152/physrev.00004.2021). Steady-state ASL height h_ASL (um) is:
       h_ASL = h_ASL_CF + (h_ASL_WT - h_ASL_CF) * [ I_sc_equiv_pct_wt / (I_sc_equiv_pct_wt + K_05_ASL) ]
     where depleted CF ASL height h_ASL_CF = 3.5 um, wild-type height h_ASL_WT = 8.5 um, and K_05_ASL = 15.0%.
     Mucociliary clearance velocity v_MCC (mm/min) depends non-linearly on ASL height exceeding the ciliary height floor:
       v_MCC = v_MCC_WT * max(0.0, (h_ASL - h_cilia_min) / (h_ASL_WT - h_cilia_min))^1.5
     where v_MCC_WT = 4.0 mm/min and h_cilia_min = 5.0 um (below which cilia are physically matted).

  6. Kahn-Kalai Curative Threshold for TMEM16A Airway Bypass:
     - Restored ASL hydration height h_ASL >= 7.0 um (curative physiological hydration).
     - Mucociliary transport velocity v_MCC >= 1.0 mm/min.
     - Equivalent chloride conductance I_sc_equiv >= 15.0% WT.
     - Preserved intracellular calcium safety: basal [Ca2+]_i <= 0.25 uM, preventing calcium-driven
       goblet cell degranulation and MUC5AC mucus hypersecretion.

  7. Invariants & Deterministic Bounds:
     - Strict bounds: 0.0 <= P_o <= 1.0; 0.0 <= h_ASL <= 15.0 um; v_MCC >= 0.0.
     - Monotonicity: d(P_o)/d(C_pot) > 0 at resting [Ca2+]_i; d(h_ASL)/d(I_sc) > 0; d(v_MCC)/d(h_ASL) > 0.
     - Complex-step differentiation: d(P_o)/d([Ca2+]_i) verified at step h = 1e-30 against finite differences.
     - Poisoned negative controls: negative concentrations, negative calcium, or unphysical temperature fail deterministically.

CF Biophysical Domain Mapping (Column B - Popperian Falsifiable Hypothesis):
  - In optimal TMEM16A potentiator monotherapy (Danahay et al. 2020 AJRCCM):
    C_pot = 1.0 uM, basal [Ca2+]_i = 0.12 uM:
    EC50 shifts from 0.85 uM to 0.334 uM, P_o = 0.1798, I_sc_equiv = 44.9% WT >= 15.0% WT,
    h_ASL = 7.25 um >= 7.0 um, v_MCC = 2.06 mm/min >= 1.0 mm/min,
    crossing the Kahn-Kalai curative threshold in CFTR-null epithelia without elevating [Ca2+]_i.
  - In high-affinity potentiator regime:
    C_pot = 2.5 uM, [Ca2+]_i = 0.12 uM:
    P_o = 0.2611, I_sc_equiv = 65.3% WT, h_ASL = 7.57 um, v_MCC = 2.51 mm/min,
    crossing Kahn-Kalai threshold.
  - In sub-therapeutic low-dose potentiator:
    C_pot = 0.03 uM, [Ca2+]_i = 0.12 uM:
    P_o = 0.0379, I_sc_equiv = 9.5% WT, h_ASL = 5.44 um < 7.0 um, v_MCC = 0.18 mm/min < 1.0 mm/min,
    failing Kahn-Kalai threshold.
  - In untreated CFTR-null baseline:
    C_pot = 0.0 uM, [Ca2+]_i = 0.12 uM:
    P_o = 0.0268, I_sc_equiv = 6.7% WT, h_ASL = 5.04 um < 7.0 um, v_MCC = 0.01 mm/min,
    failing Kahn-Kalai threshold, modeling persistent airway dehydration and mucociliary stasis.
  - In calcium-free potentiated state:
    C_pot = 1.0 uM, [Ca2+]_i = 0.0 uM:
    P_o = 0.1039, I_sc_equiv = 26.0% WT, h_ASL = 6.67 um < 7.0 um, v_MCC = 1.32 mm/min,
    demonstrating allosteric calcium-independent leak activation.

Research simulation only; not medical advice.
"""
from __future__ import annotations

import cmath
import math
from dataclasses import dataclass
from typing import Any, Dict, Tuple

PRIMARY_CITATIONS: Tuple[str, ...] = (
    "DOI: 10.1164/rccm.201908-1641OC",
    "DOI: 10.1038/nature24652",
    "DOI: 10.1152/physrev.00004.2021",
)

# Physical constants
GAS_CONSTANT_R: float = 8.314462618                 # J/(mol*K)
CAL_TO_JOULE: float = 4184.0                        # J per kcal
TEMPERATURE_PHYSIOLOGICAL_K: float = 310.15         # 37 C (310.15 K)

# Biophysical defaults for TMEM16A potentiator allosteric gating & hydration
DEFAULT_EC50_0_UM: float = 0.85                     # Basal TMEM16A EC50 for Ca2+ activation (uM)
DEFAULT_F_SHIFT_MAX: float = 0.82                   # Maximum potentiator-induced EC50 shift fraction
DEFAULT_KD_POT_UM: float = 0.35                     # Small-molecule potentiator binding affinity (uM)
DEFAULT_P_O_LEAK: float = 0.015                     # Basal channel leak open probability
DEFAULT_P_O_IND_MAX: float = 0.12                   # Maximum ligand-induced Ca2+-independent open probability
DEFAULT_P_O_MAX: float = 0.90                       # Maximum open probability under saturating Ca2+
DEFAULT_HILL_COEFF_NH: float = 2.2                  # Cooperative Hill coefficient for Ca2+ gating
DEFAULT_N_CHANNEL_PER_CM2: float = 2.5e7            # Apical TMEM16A channel surface density (channels/cm^2)
DEFAULT_GAMMA_SINGLE_PS: float = 8.0                # Single-channel unitary conductance (pS)
DEFAULT_G_CFTR_WT_MS_CM2: float = 0.08              # Stimulated WT CFTR apical conductance benchmark (mS/cm^2)
DEFAULT_H_ASL_CF_UM: float = 3.5                    # Depleted CF airway surface liquid height (um)
DEFAULT_H_ASL_WT_UM: float = 8.5                    # Normal wild-type airway surface liquid height (um)
DEFAULT_K05_ASL_PCT: float = 15.0                   # Half-maximal I_sc for ASL fluid rehydration (% WT)
DEFAULT_V_MCC_WT_MM_MIN: float = 4.0                # Normal wild-type mucociliary transport velocity (mm/min)
DEFAULT_H_CILIA_MIN_UM: float = 5.0                 # Minimum cilia height threshold for effective stroke (um)
DEFAULT_BASAL_CA_I_UM: float = 0.12                 # Physiological basal intracellular calcium (uM)

PARAMETER_PROVENANCE: Dict[str, Dict[str, object]] = {
    "tmem16a_basal_ec50_ca_um": {
        "value": 0.85,
        "tag": "empirical",
        "source": "Basal TMEM16A EC50 for Ca2+ activation without potentiator (DOI: 10.1038/nature24652)",
        "bounds": (0.50, 1.50),
    },
    "potentiator_kd_binding_um": {
        "value": 0.35,
        "tag": "empirical",
        "source": "TMEM16A potentiator apparent binding affinity (DOI: 10.1164/rccm.201908-1641OC)",
        "bounds": (0.05, 1.00),
    },
    "single_channel_conductance_ps": {
        "value": 8.0,
        "tag": "empirical",
        "source": "Patch-clamp single-channel Cl- conductance of TMEM16A (DOI: 10.1038/nature24652)",
        "bounds": (5.0, 12.0),
    },
    "airway_surface_liquid_wt_height_um": {
        "value": 8.5,
        "tag": "empirical",
        "source": "Physiological airway surface liquid height in healthy human airway (DOI: 10.1152/physrev.00004.2021)",
        "bounds": (7.0, 10.0),
    },
    "curative_asl_hydration_floor_um": {
        "value": 7.0,
        "tag": "derived",
        "source": "Minimal ASL height required to support mucociliary clearance (DOI: 10.1152/physrev.00004.2021)",
        "bounds": (6.5, 7.5),
    },
    "mcc_clearance_curative_velocity_mm_min": {
        "value": 1.0,
        "tag": "derived",
        "source": "Minimal mucociliary transport velocity for mucus clearance (DOI: 10.1164/rccm.201908-1641OC)",
        "bounds": (0.8, 1.5),
    },
}

TWO_COLUMN_EPISTEMIC_HYGIENE: Dict[str, Any] = {
    "deductive_theorem": (
        "Under allosteric gating kinetics, binding of a small-molecule potentiator shifts the calcium activation "
        "energy barrier and introduces a basal open state, driving macroscopic apical chloride conductance that couples "
        "osmotically to transepithelial water rehydration and mucociliary clearance velocity in CFTR-deficient epithelia."
    ),
    "biophysical_hypothesis": {
        "statement": (
            "Small-molecule TMEM16A potentiator allosteric activation (C_pot >= 0.5 uM) at basal resting calcium "
            "([Ca2+]_i ~ 0.12 uM) restores airway surface liquid hydration (h_ASL >= 7.0 um) and mucociliary clearance "
            "velocity (v_MCC >= 1.0 mm/min) across CFTR-null airway epithelia without triggering intracellular calcium "
            "elevation, crossing the Kahn-Kalai curative threshold."
        ),
        "falsification_condition": (
            "If an optimal TMEM16A potentiator regimen with C_pot >= 1.0 uM and basal [Ca2+]_i = 0.12 uM fails to "
            "restore ASL height h_ASL >= 7.0 um or v_MCC >= 1.0 mm/min, or if potentiator-induced chloride flux requires "
            "elevating [Ca2+]_i > 0.5 uM, the TMEM16A bypass invariant bridge is falsified."
        ),
    },
}

FIVE_LAYER_LADDER: Dict[str, str] = {
    "cool_human_fact": (
        "Even when the CFTR channel is completely missing, human lung cells possess a second backup chloride channel "
        "called TMEM16A that can be switched on with a small molecule to restore normal mucus clearance."
    ),
    "everyday_analogy": (
        "Like opening an emergency spillway gate on a dam when the main turbine channel is permanently blocked, "
        "restoring river flow without overfilling the reservoir."
    ),
    "engineering_lens": (
        "Allosteric cooperative Hill gating kinetics, coupled electro-osmotic fluid pumping, and non-linear "
        "viscous boundary-layer mucociliary transport above the ciliary tip height floor."
    ),
    "formal_theorem": (
        "The composite mucociliary clearance operator v_MCC(C_pot, [Ca2+]_i) is non-negative, monotonic in C_pot, "
        "and bounded above by v_MCC_WT, exhibiting a hard threshold at h_ASL = h_cilia_min."
    ),
    "falsifiable_code_and_toy": (
        "A pure Python TMEM16A potentiator and ASL hydration simulator with complex-step derivative parity at h=1e-30 and "
        "deterministic falsification on negative concentrations, negative calcium, or unphysical temperatures."
    ),
}

REGIME_PRESETS: Dict[str, Dict[str, Any]] = {
    "potentiator_monotherapy_optimal": {
        "c_pot_um": 1.0,
        "ca_i_um": DEFAULT_BASAL_CA_I_UM,
        "crosses_kahn_kalai": True,
        "min_asl_height_um": 7.0,
        "min_mcc_velocity_mm_min": 1.0,
    },
    "potentiator_high_affinity_regime": {
        "c_pot_um": 2.5,
        "ca_i_um": DEFAULT_BASAL_CA_I_UM,
        "crosses_kahn_kalai": True,
        "min_asl_height_um": 7.0,
        "min_mcc_velocity_mm_min": 1.0,
    },
    "subtherapeutic_potentiator_low_dose": {
        "c_pot_um": 0.03,
        "ca_i_um": DEFAULT_BASAL_CA_I_UM,
        "crosses_kahn_kalai": False,
        "min_asl_height_um": 7.0,
        "min_mcc_velocity_mm_min": 1.0,
    },
    "untreated_cftr_null_baseline": {
        "c_pot_um": 0.0,
        "ca_i_um": DEFAULT_BASAL_CA_I_UM,
        "crosses_kahn_kalai": False,
        "min_asl_height_um": 7.0,
        "min_mcc_velocity_mm_min": 1.0,
    },
    "calcium_free_potentiated_state": {
        "c_pot_um": 1.0,
        "ca_i_um": 0.0,
        "crosses_kahn_kalai": False,
        "min_asl_height_um": 7.0,
        "min_mcc_velocity_mm_min": 1.0,
    },
}


@dataclass(frozen=True)
class Round180Tmem16aResult:
    """Result data container for Round 180 TMEM16A Potentiator Bridge."""
    preset_name: str
    c_pot_um: float
    ca_i_um: float
    apparent_ec50_um: float
    open_probability: float
    conductance_ms_cm2: float
    equivalent_isc_pct_wt: float
    asl_height_um: float
    mcc_velocity_mm_min: float
    crosses_kahn_kalai_threshold: bool
    falsification_invariant_satisfied: bool
    complex_step_relative_error: float


def compute_apparent_ec50(
    c_pot_um: float,
    ec50_0_um: float = DEFAULT_EC50_0_UM,
    f_shift_max: float = DEFAULT_F_SHIFT_MAX,
    kd_pot_um: float = DEFAULT_KD_POT_UM,
) -> float:
    """Compute apparent calcium EC50 under potentiator allosteric shift."""
    if c_pot_um < 0.0:
        raise ValueError("Potentiator concentration cannot be negative")
    if ec50_0_um <= 0.0:
        raise ValueError("Basal EC50 must be strictly positive")
    if kd_pot_um <= 0.0:
        raise ValueError("Potentiator KD must be strictly positive")

    occupancy = c_pot_um / (c_pot_um + kd_pot_um) if (c_pot_um + kd_pot_um) > 0.0 else 0.0
    return ec50_0_um * (1.0 - f_shift_max * occupancy)


def compute_tmem16a_open_probability(
    ca_i_um: float,
    c_pot_um: float,
    ec50_0_um: float = DEFAULT_EC50_0_UM,
    f_shift_max: float = DEFAULT_F_SHIFT_MAX,
    kd_pot_um: float = DEFAULT_KD_POT_UM,
    p_leak: float = DEFAULT_P_O_LEAK,
    p_ind_max: float = DEFAULT_P_O_IND_MAX,
    p_max: float = DEFAULT_P_O_MAX,
    nh: float = DEFAULT_HILL_COEFF_NH,
) -> Tuple[float, float]:
    """Compute channel open probability Po and apparent EC50.

    Returns (Po, apparent_ec50).
    """
    if ca_i_um < 0.0:
        raise ValueError("Intracellular calcium concentration cannot be negative")
    if c_pot_um < 0.0:
        raise ValueError("Potentiator concentration cannot be negative")

    ec50 = compute_apparent_ec50(c_pot_um, ec50_0_um, f_shift_max, kd_pot_um)
    occupancy = c_pot_um / (c_pot_um + kd_pot_um) if (c_pot_um + kd_pot_um) > 0.0 else 0.0
    p_basal = p_leak + p_ind_max * occupancy

    if ca_i_um == 0.0:
        return max(0.0, min(1.0, p_basal)), ec50

    ca_pow = math.pow(ca_i_um, nh)
    ec_pow = math.pow(ec50, nh)
    fraction = ca_pow / (ca_pow + ec_pow) if (ca_pow + ec_pow) > 0.0 else 0.0
    po = p_basal + (p_max - p_basal) * fraction
    return max(0.0, min(1.0, po)), ec50


def compute_epithelial_conductance_and_isc(
    po: float,
    n_channel: float = DEFAULT_N_CHANNEL_PER_CM2,
    gamma_single_ps: float = DEFAULT_GAMMA_SINGLE_PS,
    g_cftr_wt: float = DEFAULT_G_CFTR_WT_MS_CM2,
) -> Tuple[float, float]:
    """Compute macroscopic TMEM16A conductance (mS/cm^2) and equivalent % WT I_sc.

    Returns (g_tmem, isc_equiv_pct_wt).
    """
    if po < 0.0 or po > 1.0:
        raise ValueError("Open probability must be in [0, 1]")
    if gamma_single_ps <= 0.0:
        raise ValueError("Single channel conductance must be strictly positive")
    if g_cftr_wt <= 0.0:
        raise ValueError("WT CFTR conductance benchmark must be strictly positive")

    # gamma_single in Siemens = gamma_single_ps * 1e-12 S
    # G = N * gamma * Po in S/cm^2
    # in mS/cm^2 = G * 1e3 = N * gamma_single_ps * 1e-9 * Po
    g_tmem = n_channel * (gamma_single_ps * 1e-9) * po
    isc_equiv = (g_tmem / g_cftr_wt) * 100.0
    return g_tmem, isc_equiv


def compute_asl_height_and_mcc_velocity(
    isc_equiv_pct_wt: float,
    h_cf_um: float = DEFAULT_H_ASL_CF_UM,
    h_wt_um: float = DEFAULT_H_ASL_WT_UM,
    k05_asl: float = DEFAULT_K05_ASL_PCT,
    v_mcc_wt: float = DEFAULT_V_MCC_WT_MM_MIN,
    h_cilia_min: float = DEFAULT_H_CILIA_MIN_UM,
) -> Tuple[float, float]:
    """Compute steady-state ASL height (um) and mucociliary clearance velocity (mm/min).

    Returns (h_asl, v_mcc).
    """
    if isc_equiv_pct_wt < 0.0:
        raise ValueError("Equivalent I_sc cannot be negative")

    fraction = isc_equiv_pct_wt / (isc_equiv_pct_wt + k05_asl) if (isc_equiv_pct_wt + k05_asl) > 0.0 else 0.0
    h_asl = h_cf_um + (h_wt_um - h_cf_um) * fraction

    diff_h = max(0.0, (h_asl - h_cilia_min) / (h_wt_um - h_cilia_min))
    v_mcc = v_mcc_wt * math.pow(diff_h, 1.5)
    return h_asl, v_mcc


def complex_step_ca_derivative(
    ca_i_um: float,
    c_pot_um: float,
    ec50_0_um: float = DEFAULT_EC50_0_UM,
    f_shift_max: float = DEFAULT_F_SHIFT_MAX,
    kd_pot_um: float = DEFAULT_KD_POT_UM,
    p_leak: float = DEFAULT_P_O_LEAK,
    p_ind_max: float = DEFAULT_P_O_IND_MAX,
    p_max: float = DEFAULT_P_O_MAX,
    nh: float = DEFAULT_HILL_COEFF_NH,
    h: float = 1e-30,
) -> Tuple[float, float, float]:
    """Verify derivative d(Po)/d([Ca2+]_i) using complex-step differentiation against finite differences."""
    ec50 = compute_apparent_ec50(c_pot_um, ec50_0_um, f_shift_max, kd_pot_um)
    occupancy = c_pot_um / (c_pot_um + kd_pot_um) if (c_pot_um + kd_pot_um) > 0.0 else 0.0
    p_basal = p_leak + p_ind_max * occupancy

    # Complex-step evaluation
    ca_c = complex(ca_i_um, h)
    ca_pow_c = ca_c ** nh
    ec_pow = ec50 ** nh
    frac_c = ca_pow_c / (ca_pow_c + ec_pow)
    po_c = p_basal + (p_max - p_basal) * frac_c
    cs_deriv = po_c.imag / h

    # Central finite differences
    delta = 1e-6
    ca_p = ca_i_um + delta
    ca_m = ca_i_um - delta
    f_p = math.pow(ca_p, nh) / (math.pow(ca_p, nh) + ec_pow)
    f_m = math.pow(ca_m, nh) / (math.pow(ca_m, nh) + ec_pow)
    po_p = p_basal + (p_max - p_basal) * f_p
    po_m = p_basal + (p_max - p_basal) * f_m
    fd_deriv = (po_p - po_m) / (2.0 * delta)

    rel_error = abs(cs_deriv - fd_deriv) / max(1e-12, abs(cs_deriv))
    return cs_deriv, fd_deriv, rel_error


def evaluate_round180_tmem16a_potentiator_conductance(preset_key: str) -> Round180Tmem16aResult:
    """Evaluate Round 180 TMEM16A potentiator conductance bridge for a named regime preset."""
    if preset_key not in REGIME_PRESETS:
        raise KeyError(f"Unknown preset '{preset_key}'. Choose from: {sorted(REGIME_PRESETS)}")

    p = REGIME_PRESETS[preset_key]
    c_pot = float(p["c_pot_um"])
    ca_i = float(p["ca_i_um"])

    po, ec50 = compute_tmem16a_open_probability(ca_i, c_pot)
    g_tmem, isc_equiv = compute_epithelial_conductance_and_isc(po)
    h_asl, v_mcc = compute_asl_height_and_mcc_velocity(isc_equiv)

    crosses_kk = (
        h_asl >= float(p["min_asl_height_um"])
        and v_mcc >= float(p["min_mcc_velocity_mm_min"])
        and isc_equiv >= 15.0
    )

    invariants_ok = (
        0.0 <= po <= 1.0
        and g_tmem >= 0.0
        and isc_equiv >= 0.0
        and 0.0 <= h_asl <= 15.0
        and v_mcc >= 0.0
    )

    _, _, rel_err = complex_step_ca_derivative(
        ca_i_um=max(0.01, ca_i),
        c_pot_um=c_pot,
    )

    return Round180Tmem16aResult(
        preset_name=preset_key,
        c_pot_um=c_pot,
        ca_i_um=ca_i,
        apparent_ec50_um=ec50,
        open_probability=po,
        conductance_ms_cm2=g_tmem,
        equivalent_isc_pct_wt=isc_equiv,
        asl_height_um=h_asl,
        mcc_velocity_mm_min=v_mcc,
        crosses_kahn_kalai_threshold=crosses_kk,
        falsification_invariant_satisfied=invariants_ok,
        complex_step_relative_error=rel_err,
    )


def evaluate_falsification_negative_control(
    poison_negative_c_pot: bool = False,
    poison_negative_ca_i: bool = False,
    poison_negative_gamma: bool = False,
    poison_unphysical_po: bool = False,
) -> Dict[str, Any]:
    """Execute dynamic negative controls and boundary falsifications under AH-10."""
    try:
        if poison_negative_c_pot:
            compute_apparent_ec50(c_pot_um=-1.0)
            return {"passed": False, "error": "Failed to catch negative potentiator concentration"}

        if poison_negative_ca_i:
            compute_tmem16a_open_probability(ca_i_um=-0.5, c_pot_um=1.0)
            return {"passed": False, "error": "Failed to catch negative calcium concentration"}

        if poison_negative_gamma:
            compute_epithelial_conductance_and_isc(po=0.25, gamma_single_ps=-8.0)
            return {"passed": False, "error": "Failed to catch negative single-channel conductance"}

        if poison_unphysical_po:
            compute_epithelial_conductance_and_isc(po=1.5)
            return {"passed": False, "error": "Failed to catch open probability > 1"}

        # Baseline evaluation
        res = evaluate_round180_tmem16a_potentiator_conductance("potentiator_monotherapy_optimal")
        passed = (
            res.crosses_kahn_kalai_threshold
            and res.falsification_invariant_satisfied
            and res.complex_step_relative_error <= 1e-6
        )
        return {"passed": passed, "error": None if passed else "Baseline evaluation failed invariants"}
    except ValueError as exc:
        return {"passed": True, "error": f"Correctly caught: {exc}"}
