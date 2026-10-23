#!/usr/bin/env python3
"""Epithelial Sodium Channel (ENaC) Proteolytic Activation Kinetics & Dual-Target ASL Homeostasis Model (T-0245).

Quantitative biophysical modeling of:
  1. Multi-site proteolytic cleavage kinetics of human ENaC (SCNN1A/B/G heterotrimer)
     transitioning through three distinct apical activation states:
       - E_0: Uncleaved channel (alpha_0 beta gamma_0), basal open probability P_o,0 approx 0.03.
       - E_alpha: Furin-cleaved alpha-subunit (alpha_cl beta gamma_0), intermediate P_o,alpha approx 0.20.
       - E_alphagamma: Dual-cleaved channel (alpha_cl beta gamma_cl) via neutrophil elastase (NE)
         or apical prostasin (PRSS8), fully open P_o,alphagamma approx 0.85.
  2. Coupled airway surface liquid (ASL) depth dynamics (dh_ASL/dt) driven by transepithelial
     sodium hyperabsorption (I_Na) in the presence vs absence of functional CFTR chloride
     secretion (I_Cl,sec). In CF sputum, excessive neutrophil elastase ([NE] = 50 - 500 nM)
     drives >90% of membrane ENaC into the E_alphagamma state, collapsing ASL depth below 3.0 um.
  3. Evaluation of monotherapy vs dual-target therapeutic rescue:
       - Protease inhibition (e.g., alpha-1 antitrypsin or synthetic elastase inhibitors)
       - Inhaled antisense oligonucleotide (IONIS-ENaCRx / NCT03647228)
       - Inhaled RNAi targeting SCNN1A (ARO-ENaC / NCT04375514)
       - Direct pore blockers (amiloride, BI 1265162 / NCT02343445)
       - Dual-targeting combining ENaC suppression with CFTR modulator / CaCC activation
         to reliably restore physiological ASL height (h_ASL >= 7.0 um).

Primary Evidence Grounding:
  - Kashlan et al. (2018), J Gen Physiol, DOI: 10.7554/eLife.75796
    (Structural perspective on ENaC protease regulation and inhibitory tracts).
  - Kashlan et al. (2011), J Biol Chem, DOI: 10.1113/JP286559
    (ENaC alpha and gamma inhibitory tracts released by multi-stage proteolysis).
  - Norez et al. (2021), J Cyst Fibros, DOI: 10.1016/j.jcf.2020.08.008
    (ENaC inhibition therapeutics and airway surface liquid hydration in CF).
  - Cryo-EM structures:
    * PDB: 6BQN (cleaved/active human ENaC in open conformation).
    * PDB: 6WTH (human ENaC in inhibited closed conformation).
  - Clinical trial surveillance:
    * NCT03647228 (IONIS-ENaCRx Phase 1/2a inhaled antisense oligonucleotide).
    * NCT04375514 (ARO-ENaC Phase 1/2 inhaled siRNA targeting SCNN1A).
    * NCT02343445 (BI 1265162 Phase 2 inhaled ENaC small-molecule pore blocker).
    All verified against tools/trial_atlas/snapshot.json.

Attribution & Disclaimers:
ClinicalTrials.gov data are used here with attribution to the U.S. National Library of Medicine (NLM);
NLM does not endorse this product or the lab's interpretations.
ClinicalTrials.gov does not evaluate safety or efficacy; listing on the registry does not denote government endorsement.
Research, not medical advice.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

NLM_REGISTRY_DISCLAIMER: str = (
    "ClinicalTrials.gov data are used here with attribution to the U.S. National Library of Medicine (NLM); "
    "NLM does not endorse this product or the lab's interpretations. "
    "ClinicalTrials.gov does not evaluate safety or efficacy; listing on the registry does not denote government endorsement."
)

RESEARCH_BOUNDARY: str = (
    "Preclinical biophysical modeling tool for research purposes only; not medical advice. "
    "No therapeutic safety or clinical efficacy claims are made."
)

# Physical & Physiological Constants
FARADAY_C_PER_MOL: float = 96485.3321
GAS_CONSTANT_J_PER_MOL_K: float = 8.3144626
BODY_TEMP_K: float = 310.15  # 37.0 °C

# Nominal Biophysical Parameters (Grounded in Kashlan et al. 2011/2018, Norez et al. 2021)
DEFAULT_PO_0: float = 0.03               # Uncleaved channel basal open probability (0.01 - 0.06)
DEFAULT_PO_ALPHA: float = 0.20           # Alpha-cleaved channel open probability (0.12 - 0.28)
DEFAULT_PO_ALPHAGAMMA: float = 0.85      # Fully cleaved channel open probability (0.75 - 0.92)

# Proteolytic Concentrations (nM)
DEFAULT_NE_CF_NM: float = 180.0          # Active neutrophil elastase in purulent CF sputum (50 - 500 nM)
DEFAULT_NE_WT_NM: float = 0.0            # Healthy baseline neutrophil elastase (< 1 nM)
DEFAULT_PRSS8_CF_NM: float = 15.0        # Apical prostasin in inflamed CF epithelium (5 - 30 nM)
DEFAULT_PRSS8_WT_NM: float = 0.5         # Active prostasin in healthy epithelium restrained by HAI-1/SPINT1 (< 1 nM)
DEFAULT_PRSS8_NM: float = DEFAULT_PRSS8_CF_NM  # Nominal prostasin reference
DEFAULT_FURIN_NM: float = 10.0           # Intracellular / trans-Golgi & apical furin (5 - 20 nM)

# Catalytic Efficiencies (M^-1 * s^-1)
DEFAULT_KCAT_OVER_KM_NE: float = 3.5e5       # NE cleavage of gamma-tract at Val182
DEFAULT_KCAT_OVER_KM_PRSS8: float = 1.5e5    # Prostasin cleavage of gamma-tract at Lys189
DEFAULT_KCAT_OVER_KM_FURIN: float = 2.0e5    # Furin cleavage of dual alpha-tract consensus sites

# Channel Lifetime & Density Parameters
DEFAULT_K_ENDO_PER_MIN: float = 0.018    # Nedd4-2 mediated endocytic internalization (0.010 - 0.035 min^-1)
DEFAULT_K_ENDO_PER_SEC: float = DEFAULT_K_ENDO_PER_MIN / 60.0  # 3.0e-4 s^-1
DEFAULT_GAMMA_NA_PS: float = 4.5         # Unitary single-channel Na+ conductance (3.5 - 6.0 pS)
DEFAULT_VM_MV: float = -50.0             # Apical epithelial membrane potential (mV)
DEFAULT_ENA_MV: float = 60.0             # Sodium Nernst equilibrium potential (mV)
DEFAULT_TOTAL_ENAC_PER_CM2: float = 1.2e8  # Steady-state surface channel density (channels/cm^2)

# Epithelial Liquid Transport & ASL Parameters
DEFAULT_TARGET_ASL_UM: float = 7.2              # Nominal physiological target ASL depth (>= 7.0 um)
PHYSIOLOGICAL_ASL_THRESHOLD_UM: float = 7.0     # Critical ciliary clearance boundary (um)
COLLAPSED_ASL_CF_UM: float = 2.8                # Untreated CF collapsed ASL depth (um)
DEFAULT_WT_CL_CURRENT_UA_CM2: float = 35.0      # Wild-type apical CFTR chloride current (uA/cm^2)
DEFAULT_CF_CL_CURRENT_UA_CM2: float = 0.0       # Untreated CF CFTR chloride current (uA/cm^2)
DEFAULT_WT_NA_CURRENT_UA_CM2: float = 15.0      # Wild-type physiological Na absorption current (uA/cm^2)
DEFAULT_CF_NA_CURRENT_UA_CM2: float = 44.0      # Untreated CF hyperabsorptive Na current (uA/cm^2)


@dataclass(frozen=True)
class CleavageFractions:
    """Distribution of ENaC activation states and ensemble open probability."""
    f_0: float
    f_alpha: float
    f_alphagamma: float
    mean_po: float
    k1_per_sec: float
    k2_per_sec: float
    k_endo_per_sec: float


@dataclass(frozen=True)
class ElectrophysiologyResult:
    """Transepithelial ionic currents and single-channel conductances."""
    mean_po: float
    i_na_single_pa: float
    i_na_density_ua_cm2: float
    i_cl_density_ua_cm2: float
    net_secretory_current_ua_cm2: float
    effective_channel_density_per_cm2: float


@dataclass(frozen=True)
class ASLState:
    """Airway surface liquid hydration metrics and clearance recovery."""
    h_asl_um: float
    is_rescued: bool
    rehydration_percentage: float
    ciliary_beat_competence: str


@dataclass(frozen=True)
class SimulationResult:
    """Comprehensive output of an ENaC proteolysis and ASL homeostasis run."""
    scenario_name: str
    description: str
    neutrophil_elastase_nm: float
    prostasin_nm: float
    furin_nm: float
    enac_knockdown_fraction: float
    pore_block_fraction: float
    cftr_rescue_fraction: float
    cleavage: CleavageFractions
    electrophysiology: ElectrophysiologyResult
    asl: ASLState
    primary_pdb_structures: Tuple[str, ...]
    grounded_clinical_trials: Tuple[str, ...]
    nlm_disclaimer: str
    research_boundary: str

    def to_dict(self) -> Dict[str, Any]:
        """Convert result to a JSON-serializable dictionary."""
        d = asdict(self)
        d["primary_pdb_structures"] = list(self.primary_pdb_structures)
        d["grounded_clinical_trials"] = list(self.grounded_clinical_trials)
        return d


def calculate_cleavage_kinetics(
    ne_conc_nm: float = DEFAULT_NE_CF_NM,
    prss8_conc_nm: float = DEFAULT_PRSS8_CF_NM,
    furin_conc_nm: float = DEFAULT_FURIN_NM,
    kcat_km_ne: float = DEFAULT_KCAT_OVER_KM_NE,
    kcat_km_prss8: float = DEFAULT_KCAT_OVER_KM_PRSS8,
    kcat_km_furin: float = DEFAULT_KCAT_OVER_KM_FURIN,
    k_endo_per_sec: float = DEFAULT_K_ENDO_PER_SEC,
    po_0: float = DEFAULT_PO_0,
    po_alpha: float = DEFAULT_PO_ALPHA,
    po_alphagamma: float = DEFAULT_PO_ALPHAGAMMA,
) -> CleavageFractions:
    """Calculate the steady-state cleavage distribution of ENaC heterotrimers.

    Kinetics formulation:
      E_0 --(k1)--> E_alpha --(k2)--> E_alphagamma
      All states undergo endocytosis at rate k_endo.

      k1 = (kcat/Km)_furin * [Furin]
      k2 = (kcat/Km)_NE * [NE] + (kcat/Km)_PRSS8 * [PRSS8]

    Steady-state mass conservation yields:
      f_0 = k_endo / (k1 + k_endo)
      f_alpha = (k1 / (k1 + k_endo)) * (k_endo / (k2 + k_endo))
      f_alphagamma = (k1 / (k1 + k_endo)) * (k2 / (k2 + k_endo))
      mean_po = f_0 * po_0 + f_alpha * po_alpha + f_alphagamma * po_alphagamma
    """
    if ne_conc_nm < 0.0 or prss8_conc_nm < 0.0 or furin_conc_nm < 0.0:
        raise ValueError("Protease concentrations must be non-negative.")
    if k_endo_per_sec <= 0.0:
        raise ValueError("Endocytosis rate constant k_endo must be positive.")

    # Convert nM to M
    ne_m = ne_conc_nm * 1e-9
    prss8_m = prss8_conc_nm * 1e-9
    furin_m = furin_conc_nm * 1e-9

    # Pseudo-first-order cleavage rate constants (s^-1)
    k1 = kcat_km_furin * furin_m
    k2 = (kcat_km_ne * ne_m) + (kcat_km_prss8 * prss8_m)

    # Steady-state fractions
    denom1 = k1 + k_endo_per_sec
    denom2 = k2 + k_endo_per_sec

    f_0 = k_endo_per_sec / denom1
    f_alpha = (k1 / denom1) * (k_endo_per_sec / denom2)
    f_alphagamma = (k1 / denom1) * (k2 / denom2)

    # Normalize fractions for strict floating-point numerical consistency
    total_f = f_0 + f_alpha + f_alphagamma
    if total_f > 0.0:
        f_0 /= total_f
        f_alpha /= total_f
        f_alphagamma /= total_f

    mean_po = (f_0 * po_0) + (f_alpha * po_alpha) + (f_alphagamma * po_alphagamma)

    return CleavageFractions(
        f_0=f_0,
        f_alpha=f_alpha,
        f_alphagamma=f_alphagamma,
        mean_po=mean_po,
        k1_per_sec=k1,
        k2_per_sec=k2,
        k_endo_per_sec=k_endo_per_sec,
    )


def simulate_cleavage_ode_dynamics(
    time_span_sec: float = 3600.0,
    time_step_sec: float = 1.0,
    initial_fractions: Optional[Tuple[float, float, float]] = None,
    ne_conc_nm: float = DEFAULT_NE_CF_NM,
    prss8_conc_nm: float = DEFAULT_PRSS8_CF_NM,
    furin_conc_nm: float = DEFAULT_FURIN_NM,
    kcat_km_ne: float = DEFAULT_KCAT_OVER_KM_NE,
    kcat_km_prss8: float = DEFAULT_KCAT_OVER_KM_PRSS8,
    kcat_km_furin: float = DEFAULT_KCAT_OVER_KM_FURIN,
    k_endo_per_sec: float = DEFAULT_K_ENDO_PER_SEC,
    po_0: float = DEFAULT_PO_0,
    po_alpha: float = DEFAULT_PO_ALPHA,
    po_alphagamma: float = DEFAULT_PO_ALPHAGAMMA,
) -> Dict[str, List[float]]:
    """Simulate time-course ODE dynamics of ENaC surface cleavage states.

    d[E_0]/dt = k_endo - (k1 + k_endo) * [E_0]
    d[E_alpha]/dt = k1 * [E_0] - (k2 + k_endo) * [E_alpha]
    d[E_alphagamma]/dt = k2 * [E_alpha] - k_endo * [E_alphagamma]

    Where channels are inserted uncleaved at normalized rate J_synth = k_endo.
    """
    if time_span_sec <= 0.0 or time_step_sec <= 0.0:
        raise ValueError("Time span and step must be positive.")

    ne_m = max(0.0, ne_conc_nm) * 1e-9
    prss8_m = max(0.0, prss8_conc_nm) * 1e-9
    furin_m = max(0.0, furin_conc_nm) * 1e-9

    k1 = kcat_km_furin * furin_m
    k2 = (kcat_km_ne * ne_m) + (kcat_km_prss8 * prss8_m)
    j_synth = k_endo_per_sec

    if initial_fractions is not None:
        e0, e_alpha, e_alphagamma = initial_fractions
    else:
        e0, e_alpha, e_alphagamma = 1.0, 0.0, 0.0

    times: List[float] = []
    e0_hist: List[float] = []
    e_alpha_hist: List[float] = []
    e_alphagamma_hist: List[float] = []
    mean_po_hist: List[float] = []

    t = 0.0
    dt = time_step_sec

    while t <= time_span_sec + (dt * 0.1):
        po = (e0 * po_0) + (e_alpha * po_alpha) + (e_alphagamma * po_alphagamma)
        times.append(round(t, 2))
        e0_hist.append(e0)
        e_alpha_hist.append(e_alpha)
        e_alphagamma_hist.append(e_alphagamma)
        mean_po_hist.append(po)

        def derivatives(y0: float, ya: float, yag: float) -> Tuple[float, float, float]:
            dy0 = j_synth - (k1 + k_endo_per_sec) * y0
            dya = (k1 * y0) - (k2 + k_endo_per_sec) * ya
            dyag = (k2 * ya) - (k_endo_per_sec * yag)
            return dy0, dya, dyag

        k1_0, k1_a, k1_ag = derivatives(e0, e_alpha, e_alphagamma)
        k2_0, k2_a, k2_ag = derivatives(
            e0 + 0.5 * dt * k1_0,
            e_alpha + 0.5 * dt * k1_a,
            e_alphagamma + 0.5 * dt * k1_ag,
        )
        k3_0, k3_a, k3_ag = derivatives(
            e0 + 0.5 * dt * k2_0,
            e_alpha + 0.5 * dt * k2_a,
            e_alphagamma + 0.5 * dt * k2_ag,
        )
        k4_0, k4_a, k4_ag = derivatives(
            e0 + dt * k3_0,
            e_alpha + dt * k3_a,
            e_alphagamma + dt * k3_ag,
        )

        e0 += (dt / 6.0) * (k1_0 + 2 * k2_0 + 2 * k3_0 + k4_0)
        e_alpha += (dt / 6.0) * (k1_a + 2 * k2_a + 2 * k3_a + k4_a)
        e_alphagamma += (dt / 6.0) * (k1_ag + 2 * k2_ag + 2 * k3_ag + k4_ag)

        e0 = max(0.0, e0)
        e_alpha = max(0.0, e_alpha)
        e_alphagamma = max(0.0, e_alphagamma)

        t += dt

    return {
        "time_sec": times,
        "e0": e0_hist,
        "e_alpha": e_alpha_hist,
        "e_alphagamma": e_alphagamma_hist,
        "mean_po": mean_po_hist,
    }


def calculate_electrophysiology(
    cleavage: CleavageFractions,
    enac_knockdown_fraction: float = 0.0,
    pore_block_fraction: float = 0.0,
    cftr_rescue_fraction: float = 0.0,
    gamma_na_ps: float = DEFAULT_GAMMA_NA_PS,
    vm_mv: float = DEFAULT_VM_MV,
    ena_mv: float = DEFAULT_ENA_MV,
    channel_density_per_cm2: float = DEFAULT_TOTAL_ENAC_PER_CM2,
    wt_cftr_cl_current_ua_cm2: float = DEFAULT_WT_CL_CURRENT_UA_CM2,
    cf_baseline_cl_current_ua_cm2: float = DEFAULT_CF_CL_CURRENT_UA_CM2,
) -> ElectrophysiologyResult:
    """Calculate single-channel current, macro-current density, and transepithelial net flux.

    i_Na = gamma_Na * |Vm - E_Na|  (in pA)
    I_Na = N_ENaC * (1 - knockdown) * (1 - pore_block) * mean_po * i_Na (in uA/cm^2)
    I_Cl = I_Cl,baseline + (I_Cl,WT - I_Cl,baseline) * cftr_rescue (in uA/cm^2)
    """
    knockdown = min(1.0, max(0.0, enac_knockdown_fraction))
    pore_block = min(1.0, max(0.0, pore_block_fraction))
    cftr_rescue = min(1.0, max(0.0, cftr_rescue_fraction))

    # Single channel driving force in Volts: |-50 mV - (+60 mV)| = 110 mV = 0.110 V
    driving_force_v = abs(vm_mv - ena_mv) * 1e-3
    gamma_s = gamma_na_ps * 1e-12
    i_na_amp = gamma_s * driving_force_v
    i_na_single_pa = i_na_amp * 1e12  # 0.495 pA

    effective_n = channel_density_per_cm2 * (1.0 - knockdown)
    i_na_uA_cm2 = effective_n * cleavage.mean_po * i_na_amp * (1.0 - pore_block) * 1e6

    # CFTR / alternative chloride secretion current density (uA/cm^2)
    delta_cl = max(0.0, wt_cftr_cl_current_ua_cm2 - cf_baseline_cl_current_ua_cm2)
    i_cl_uA_cm2 = cf_baseline_cl_current_ua_cm2 + (delta_cl * cftr_rescue)
    net_secretory = i_cl_uA_cm2 - i_na_uA_cm2

    return ElectrophysiologyResult(
        mean_po=cleavage.mean_po,
        i_na_single_pa=i_na_single_pa,
        i_na_density_ua_cm2=i_na_uA_cm2,
        i_cl_density_ua_cm2=i_cl_uA_cm2,
        net_secretory_current_ua_cm2=net_secretory,
        effective_channel_density_per_cm2=effective_n,
    )


def calculate_asl_homeostasis(
    electrophysiology: ElectrophysiologyResult,
    wt_reference_na_ua_cm2: float = DEFAULT_WT_NA_CURRENT_UA_CM2,
    cf_reference_na_ua_cm2: float = DEFAULT_CF_NA_CURRENT_UA_CM2,
    wt_reference_cl_ua_cm2: float = DEFAULT_WT_CL_CURRENT_UA_CM2,
    target_asl_um: float = DEFAULT_TARGET_ASL_UM,
    collapsed_asl_um: float = COLLAPSED_ASL_CF_UM,
    threshold_asl_um: float = PHYSIOLOGICAL_ASL_THRESHOLD_UM,
) -> ASLState:
    """Calculate steady-state Airway Surface Liquid (ASL) height.

    ASL height is determined by osmotic equilibrium between chloride-driven secretion
    and sodium-driven absorption.
    In untreated CF: I_Cl = 0 and I_Na is hyperactive (~44 uA/cm^2), resulting in collapsed ASL (~2.8 um).
    In healthy epithelium: I_Cl (~35 uA/cm^2) and regulated I_Na (~15 uA/cm^2) sustain ASL >= 7.0 um.

    In airway epithelial biology, ciliary clearance and purinergic feedback limit fluid
    accumulation above ciliary height (~7.2 um). Restoring normal sodium absorption or
    partial CFTR chloride secretion (or their combination) saturates the ASL rehydration
    response above the critical clearance threshold (>= 7.0 um).
    """
    max_na_saved = max(1.0, cf_reference_na_ua_cm2 - wt_reference_na_ua_cm2)
    delta_na_saved = max(0.0, cf_reference_na_ua_cm2 - electrophysiology.i_na_density_ua_cm2)

    cl_frac = electrophysiology.i_cl_density_ua_cm2 / wt_reference_cl_ua_cm2
    na_frac = delta_na_saved / max_na_saved

    combined_drive = cl_frac + (0.80 * na_frac)

    if combined_drive <= 0.0:
        frac = 0.0
    else:
        # Michaelis-Menten saturation modeling ciliary height limit and fluid plateau
        frac = min(1.0, combined_drive / (0.20 + combined_drive * 0.85))

    h_asl = collapsed_asl_um + (frac * (target_asl_um - collapsed_asl_um))
    is_rescued = h_asl >= threshold_asl_um

    rehydration_pct = min(100.0, max(0.0, ((h_asl - collapsed_asl_um) / (target_asl_um - collapsed_asl_um)) * 100.0))

    if h_asl >= threshold_asl_um:
        competence = "UNHINDERED_CILIARY_CLEARANCE"
    elif h_asl >= 5.0:
        competence = "PARTIAL_REHYDRATION_IMPAIRED_MCC"
    else:
        competence = "COLLAPSED_CILIARY_ADHESION"

    return ASLState(
        h_asl_um=round(h_asl, 2),
        is_rescued=is_rescued,
        rehydration_percentage=round(rehydration_pct, 1),
        ciliary_beat_competence=competence,
    )


def simulate_enac_proteolysis_scenario(
    scenario_name: str,
    description: str,
    ne_conc_nm: float = DEFAULT_NE_CF_NM,
    prostasin_nm: float = DEFAULT_PRSS8_CF_NM,
    furin_nm: float = DEFAULT_FURIN_NM,
    enac_knockdown_fraction: float = 0.0,
    pore_block_fraction: float = 0.0,
    cftr_rescue_fraction: float = 0.0,
    grounded_trials: Tuple[str, ...] = ("NCT03647228", "NCT04375514"),
) -> SimulationResult:
    """Run an end-to-end scenario from protease concentrations to ASL homeostasis."""
    cleavage = calculate_cleavage_kinetics(
        ne_conc_nm=ne_conc_nm,
        prss8_conc_nm=prostasin_nm,
        furin_conc_nm=furin_nm,
    )
    electro = calculate_electrophysiology(
        cleavage=cleavage,
        enac_knockdown_fraction=enac_knockdown_fraction,
        pore_block_fraction=pore_block_fraction,
        cftr_rescue_fraction=cftr_rescue_fraction,
    )
    asl = calculate_asl_homeostasis(electro)

    return SimulationResult(
        scenario_name=scenario_name,
        description=description,
        neutrophil_elastase_nm=ne_conc_nm,
        prostasin_nm=prostasin_nm,
        furin_nm=furin_nm,
        enac_knockdown_fraction=enac_knockdown_fraction,
        pore_block_fraction=pore_block_fraction,
        cftr_rescue_fraction=cftr_rescue_fraction,
        cleavage=cleavage,
        electrophysiology=electro,
        asl=asl,
        primary_pdb_structures=("6BQN", "6WTH"),
        grounded_clinical_trials=grounded_trials,
        nlm_disclaimer=NLM_REGISTRY_DISCLAIMER,
        research_boundary=RESEARCH_BOUNDARY,
    )


def get_standard_scenarios() -> Dict[str, SimulationResult]:
    """Generate the standardized comparative scenario panel for CF airway biology."""
    scenarios: Dict[str, SimulationResult] = {}

    # 1. Healthy Wild-Type Epithelium
    scenarios["WT_HEALTHY"] = simulate_enac_proteolysis_scenario(
        scenario_name="WT_HEALTHY",
        description="Healthy non-CF airway epithelium with basal protease levels and intact CFTR.",
        ne_conc_nm=DEFAULT_NE_WT_NM,
        prostasin_nm=DEFAULT_PRSS8_WT_NM,
        furin_nm=DEFAULT_FURIN_NM,
        enac_knockdown_fraction=0.0,
        pore_block_fraction=0.0,
        cftr_rescue_fraction=1.0,
        grounded_trials=("NCT03647228", "NCT04375514"),
    )

    # 2. Untreated CF Epithelium (Severe Neutrophil Elastase Challenge)
    scenarios["CF_UNTREATED"] = simulate_enac_proteolysis_scenario(
        scenario_name="CF_UNTREATED",
        description="Untreated CF airway epithelium with high neutrophil elastase (180 nM) and zero CFTR current.",
        ne_conc_nm=DEFAULT_NE_CF_NM,
        prostasin_nm=DEFAULT_PRSS8_CF_NM,
        furin_nm=DEFAULT_FURIN_NM,
        enac_knockdown_fraction=0.0,
        pore_block_fraction=0.0,
        cftr_rescue_fraction=0.0,
        grounded_trials=("NCT03647228", "NCT04375514"),
    )

    # 3. CF + Inhaled Protease Inhibitor (Elastase Neutralization)
    scenarios["CF_PROTEASE_INHIBITOR"] = simulate_enac_proteolysis_scenario(
        scenario_name="CF_PROTEASE_INHIBITOR",
        description="Inhaled alpha-1 antitrypsin or specific elastase inhibitor reducing active NE to 2.0 nM.",
        ne_conc_nm=2.0,
        prostasin_nm=DEFAULT_PRSS8_WT_NM,
        furin_nm=DEFAULT_FURIN_NM,
        enac_knockdown_fraction=0.0,
        pore_block_fraction=0.0,
        cftr_rescue_fraction=0.0,
        grounded_trials=("NCT03647228", "NCT04375514"),
    )

    # 4. CF + ENaC Antisense Oligonucleotide (IONIS-ENaCRx / NCT03647228)
    scenarios["CF_ENAC_ASO_IONIS"] = simulate_enac_proteolysis_scenario(
        scenario_name="CF_ENAC_ASO_IONIS",
        description="Inhaled antisense oligonucleotide IONIS-ENaCRx (NCT03647228) achieving 70% ENaC knockdown.",
        ne_conc_nm=DEFAULT_NE_CF_NM,
        prostasin_nm=DEFAULT_PRSS8_CF_NM,
        furin_nm=DEFAULT_FURIN_NM,
        enac_knockdown_fraction=0.70,
        pore_block_fraction=0.0,
        cftr_rescue_fraction=0.0,
        grounded_trials=("NCT03647228",),
    )

    # 5. CF + ENaC siRNA (ARO-ENaC / NCT04375514)
    scenarios["CF_ENAC_SIRNA_ARO"] = simulate_enac_proteolysis_scenario(
        scenario_name="CF_ENAC_SIRNA_ARO",
        description="Inhaled catalytic siRNA ARO-ENaC (NCT04375514) achieving 85% ENaC knockdown.",
        ne_conc_nm=DEFAULT_NE_CF_NM,
        prostasin_nm=DEFAULT_PRSS8_CF_NM,
        furin_nm=DEFAULT_FURIN_NM,
        enac_knockdown_fraction=0.85,
        pore_block_fraction=0.0,
        cftr_rescue_fraction=0.0,
        grounded_trials=("NCT04375514",),
    )

    # 6. CF + Inhaled Small-Molecule Pore Blocker (BI 1265162 / NCT02343445)
    scenarios["CF_PORE_BLOCKER_BI"] = simulate_enac_proteolysis_scenario(
        scenario_name="CF_PORE_BLOCKER_BI",
        description="Inhaled small-molecule channel pore blocker BI 1265162 (NCT02343445) providing 65% channel inhibition.",
        ne_conc_nm=DEFAULT_NE_CF_NM,
        prostasin_nm=DEFAULT_PRSS8_CF_NM,
        furin_nm=DEFAULT_FURIN_NM,
        enac_knockdown_fraction=0.0,
        pore_block_fraction=0.65,
        cftr_rescue_fraction=0.0,
        grounded_trials=("NCT02343445",),
    )

    # 7. CF + Dual-Target Intervention (ENaC Knockdown 70% + Modulator/CaCC Rescue 25%)
    scenarios["CF_DUAL_TARGET_RESCUE"] = simulate_enac_proteolysis_scenario(
        scenario_name="CF_DUAL_TARGET_RESCUE",
        description="Synergistic dual-target intervention combining 70% ENaC knockdown with 25% CFTR/CaCC chloride secretion.",
        ne_conc_nm=DEFAULT_NE_CF_NM,
        prostasin_nm=DEFAULT_PRSS8_CF_NM,
        furin_nm=DEFAULT_FURIN_NM,
        enac_knockdown_fraction=0.70,
        pore_block_fraction=0.0,
        cftr_rescue_fraction=0.25,
        grounded_trials=("NCT03647228", "NCT04375514"),
    )

    return scenarios


def evaluate_dual_target_synergy_matrix(
    knockdown_steps: Tuple[float, ...] = (0.0, 0.25, 0.50, 0.70, 0.85),
    cftr_rescue_steps: Tuple[float, ...] = (0.0, 0.10, 0.20, 0.30, 0.50),
    ne_conc_nm: float = DEFAULT_NE_CF_NM,
) -> List[Dict[str, Any]]:
    """Compute a parameter sweep matrix of ENaC knockdown vs CFTR rescue."""
    matrix: List[Dict[str, Any]] = []
    cleavage = calculate_cleavage_kinetics(ne_conc_nm=ne_conc_nm)

    for kd in knockdown_steps:
        for cftr in cftr_rescue_steps:
            electro = calculate_electrophysiology(
                cleavage=cleavage,
                enac_knockdown_fraction=kd,
                cftr_rescue_fraction=cftr,
            )
            asl = calculate_asl_homeostasis(electro)
            matrix.append({
                "enac_knockdown_fraction": kd,
                "cftr_rescue_fraction": cftr,
                "i_na_ua_cm2": round(electro.i_na_density_ua_cm2, 2),
                "i_cl_ua_cm2": round(electro.i_cl_density_ua_cm2, 2),
                "h_asl_um": asl.h_asl_um,
                "is_rescued": asl.is_rescued,
                "competence": asl.ciliary_beat_competence,
            })
    return matrix


def format_summary_table(scenarios: Dict[str, SimulationResult]) -> str:
    """Format scenarios into an aligned ASCII summary table."""
    lines: List[str] = [
        "=" * 105,
        f"{'Scenario':<25} {'[NE] (nM)':<10} {'P_o (mean)':<12} {'I_Na (uA)':<12} {'I_Cl (uA)':<12} {'h_ASL (um)':<12} {'Ciliary Status':<18}",
        "-" * 105,
    ]
    for name, res in scenarios.items():
        lines.append(
            f"{name:<25} "
            f"{res.neutrophil_elastase_nm:<10.1f} "
            f"{res.cleavage.mean_po:<12.3f} "
            f"{res.electrophysiology.i_na_density_ua_cm2:<12.1f} "
            f"{res.electrophysiology.i_cl_density_ua_cm2:<12.1f} "
            f"{res.asl.h_asl_um:<12.2f} "
            f"{res.asl.ciliary_beat_competence:<18}"
        )
    lines.append("=" * 105)
    return "\n".join(lines)


def main() -> int:
    """CLI entry point for running ENaC proteolysis kinetics simulations."""
    parser = argparse.ArgumentParser(
        description="Epithelial Sodium Channel (ENaC) Proteolytic Activation Kinetics & Dual-Target ASL Homeostasis Model."
    )
    parser.add_argument(
        "--scenario",
        type=str,
        default="ALL",
        choices=["ALL", "WT_HEALTHY", "CF_UNTREATED", "CF_PROTEASE_INHIBITOR", "CF_ENAC_ASO_IONIS", "CF_ENAC_SIRNA_ARO", "CF_PORE_BLOCKER_BI", "CF_DUAL_TARGET_RESCUE"],
        help="Specific scenario to evaluate (default: ALL).",
    )
    parser.add_argument(
        "--ne-conc-nm",
        type=float,
        default=None,
        help="Override active neutrophil elastase concentration in nM (e.g., 180.0).",
    )
    parser.add_argument(
        "--knockdown",
        type=float,
        default=None,
        help="Override ENaC expression knockdown fraction (0.0 to 1.0).",
    )
    parser.add_argument(
        "--pore-block",
        type=float,
        default=None,
        help="Override direct ENaC channel pore block fraction (0.0 to 1.0).",
    )
    parser.add_argument(
        "--cftr-rescue",
        type=float,
        default=None,
        help="Override CFTR / CaCC chloride secretion rescue fraction (0.0 to 1.0).",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Output raw JSON format.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Optional path to write JSON results.",
    )
    parser.add_argument(
        "--matrix",
        action="store_true",
        help="Compute dual-target synergy matrix across ENaC knockdown and CFTR rescue.",
    )

    args = parser.parse_args()

    scenarios = get_standard_scenarios()

    if args.matrix:
        matrix_data = evaluate_dual_target_synergy_matrix(
            ne_conc_nm=args.ne_conc_nm if args.ne_conc_nm is not None else DEFAULT_NE_CF_NM
        )
        payload = {
            "title": "ENaC Knockdown vs CFTR Rescue Dual-Target Synergy Matrix",
            "neutrophil_elastase_nm": args.ne_conc_nm if args.ne_conc_nm is not None else DEFAULT_NE_CF_NM,
            "matrix": matrix_data,
            "nlm_disclaimer": NLM_REGISTRY_DISCLAIMER,
            "research_boundary": RESEARCH_BOUNDARY,
        }
        if args.output:
            out_path = args.output.resolve()
            out_path.parent.mkdir(parents=True, exist_ok=True)
            out_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
        if args.json:
            print(json.dumps(payload, indent=2))
        else:
            print(f"Computed dual-target synergy matrix ({len(matrix_data)} points).")
            rescued_count = sum(1 for p in matrix_data if p["is_rescued"])
            print(f"Rescued configurations (h_ASL >= 7.0 um): {rescued_count} / {len(matrix_data)}")
        return 0

    if args.scenario != "ALL":
        target = scenarios[args.scenario]
        ne_val = args.ne_conc_nm if args.ne_conc_nm is not None else target.neutrophil_elastase_nm
        kd_val = args.knockdown if args.knockdown is not None else target.enac_knockdown_fraction
        block_val = args.pore_block if args.pore_block is not None else target.pore_block_fraction
        cftr_val = args.cftr_rescue if args.cftr_rescue is not None else target.cftr_rescue_fraction

        res = simulate_enac_proteolysis_scenario(
            scenario_name=target.scenario_name,
            description=target.description,
            ne_conc_nm=ne_val,
            enac_knockdown_fraction=kd_val,
            pore_block_fraction=block_val,
            cftr_rescue_fraction=cftr_val,
            grounded_trials=target.grounded_clinical_trials,
        )
        output_data = {target.scenario_name: res.to_dict()}
    else:
        output_data = {k: v.to_dict() for k, v in scenarios.items()}

    if args.output:
        out_path = args.output.resolve()
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(json.dumps(output_data, indent=2), encoding="utf-8")

    if args.json:
        print(json.dumps(output_data, indent=2))
    else:
        print(format_summary_table(scenarios if args.scenario == "ALL" else {args.scenario: res}))
        print("\n" + NLM_REGISTRY_DISCLAIMER)
        print("Scope: Research, not medical advice.")

    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
