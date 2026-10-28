#!/usr/bin/env python3
"""Thermodynamic and Conductance Modeling of Alternative Anion Secretion (TMEM16A Potentiation & SLC26A9 Coupling) in CFTR-Deficient Airway Epithelia (T-0476).

Quantitative biophysical, electrophysiological, and airway surface liquid (ASL) fluid dynamics model:
1. TMEM16A / ANO1 Calcium-Activated Anion Conductance:
   - Formulates allosteric Hill intracellular calcium activation (Hill coefficient nH ~ 2.0, basal EC50 ~ 0.4 uM).
   - Voltage-dependent gating charge modulation (zc ~ 0.35) and allosteric left-shift by small-molecule potentiators
     (reducing apparent EC50 from 0.4 uM to <= 0.15 uM at resting cytosolic [Ca2+]i ~ 100 nM).
   - Dimeric architecture (unitary conductance gamma_TMEM_Cl ~ 8.3 pS, P_HCO3 / P_Cl ~ 0.12).
2. SLC26A9 Constitutive Halide Transport & Macromolecular PDZ Scaffolding:
   - Uncoupled constitutive chloride and bicarbonate transport independent of elevated [Ca2+]i (Po ~ 0.85).
   - Unitary conductance gamma_SLC_Cl ~ 3.0 pS, P_HCO3 / P_Cl ~ 0.18.
   - Apical membrane trafficking regulation: in CFTR-deficient epithelia, loss of macromolecular PDZ scaffolding
     (NHERF1 / EBP50 association via C-terminal -TRL motif) and kinase dysregulation downregulates apical density
     (CF baseline: ~300 channels/cell; rescued/stabilized: ~2500 channels/cell).
3. Coupled ENaC Sodium Absorption:
   - Hyperactivated sodium absorption in CFTR-null epithelia (Po ~ 0.85, N_ENaC ~ 3000 channels/cell, gamma_Na ~ 4.5 pS).
   - Models therapeutic catalytic knockdown (e.g., inhaled siRNA or channel block) relieving fluid absorption.
4. Transepithelial Kirchhoff Circuit Equilibrium & Apical Membrane Potential (Vm):
   - Closed-form analytical equilibrium solution balancing all ionic currents:
     I_Cl,TMEM + I_HCO3,TMEM + I_Cl,SLC + I_HCO3,SLC + I_Na,ENaC + I_K,apical + I_basal = 0.
5. ASL Fluid Dynamics & Mucociliary Clearance Threshold Restoration:
   - Dynamic osmotic water flux driven by net osmolyte secretion and ENaC relief.
   - Dynamic ASL height balance: restores collapsed CF baseline (< 3.5 um, baseline ~3.2 um) above the critical
     ciliary beat clearance threshold (>= 7.0 um, target 7.5 - 8.5 um, Button et al. 2021).
   - Henderson-Hasselbalch ASL pH buffering: restores acidic CF ASL (pH < 6.85) to physiological alkaline buffering
     (pH >= 7.15) via combined TMEM16A and SLC26A9 bicarbonate secretion.

Primary Evidence Grounding:
- Primary Literature:
  - Button B et al. (2021) Physiol Rev, DOI: 10.1152/physrev.00004.2021
    (Physiology and pathophysiology of human airway mucus)
  - Middleton PG et al. (2019) N Engl J Med, DOI: 10.1056/nejmoa1908639
    (Elexacaftor-Tezacaftor-Ivacaftor for Cystic Fibrosis with a Single Phe508del Allele)
- Clinical Registry Grounding:
  - NCT05033080: Phase 3 Study of Vanzacaftor-Tezacaftor-Deutivacaftor in Cystic Fibrosis
  - NCT05668741: Phase 1/2 Study of Inhaled mRNA LNP VX-522 in Participants With Cystic Fibrosis
- Structural Grounding:
  - RCSB PDB 6MSM (Human CFTR dephosphorylated ATP-free conformation, 3.9 A Cryo-EM)

Epistemic Status:
- Answer Tier T2 (Open-access peer-reviewed electrophysiology, ClinicalTrials.gov API v2).
- Preclinical in-silico computational model; research, not medical advice.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
import math
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

SCRIPT_PATH = Path(__file__).resolve()

# Attribution & Mandatory Disclaimers
NLM_ATTRIBUTION: str = (
    "This product uses publicly available data from the U.S. National Library of Medicine (NLM), "
    "National Institutes of Health, Department of Health and Human Services; NLM is not responsible "
    "for the product and does not endorse or recommend this or any other product. "
    "ClinicalTrials.gov does not evaluate safety or efficacy; listing on the registry does not denote government endorsement."
)

RESEARCH_BOUNDARY: str = (
    "Preclinical in-silico coupled TMEM16A-SLC26A9 alternative anion conductance and ASL hydration model; "
    "research, not medical advice. No clinical efficacy or therapeutic safety claim is made."
)

# Reference Grounding Identifiers (Verified in sources/cite-check-ledger.jsonl & snapshot.json)
PRIMARY_DOI_BUTTON: str = "10.1152/physrev.00004.2021"
PRIMARY_DOI_MIDDLETON: str = "10.1056/nejmoa1908639"
CLINICAL_TRIAL_VANZACAFTOR: str = "NCT05033080"
CLINICAL_TRIAL_VX522: str = "NCT05668741"
PDB_STRUCTURE_CFTR: str = "6MSM"

# Physical and Thermodynamic Constants (37 C / 310.15 K)
FARADAY_CONSTANT_C_PER_MOL: float = 96485.3321
UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K: float = 8.3144626
BODY_TEMPERATURE_K: float = 310.15
THERMAL_VOLTAGE_MV: float = (UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * BODY_TEMPERATURE_K / FARADAY_CONSTANT_C_PER_MOL) * 1000.0  # ~26.726 mV

# Physiological Nernst Equilibrium Potentials (mV)
DEFAULT_E_CL_MV: float = -32.0          # Apical chloride equilibrium potential
DEFAULT_E_HCO3_MV: float = -18.0        # Apical bicarbonate equilibrium potential
DEFAULT_E_NA_MV: float = 60.0           # Apical sodium equilibrium potential (ENaC)
DEFAULT_E_K_MV: float = -85.0           # Apical potassium leak equilibrium potential
DEFAULT_E_BASAL_MV: float = -70.0       # Basolateral resting/shunt potential

# Conductances and Cell Geometry
DEFAULT_APICAL_AREA_UM2: float = 100.0  # Human airway columnar ciliated cell apical surface area (100 um^2)
DEFAULT_G_BASAL_PS: float = 15000.0     # Basolateral equivalent conductance (pS/cell)
DEFAULT_G_K_APICAL_PS: float = 1500.0   # Apical K+ leak conductance (pS/cell)

# TMEM16A / ANO1 Biophysical Parameters
DEFAULT_TMEM_GAMMA_CL_PS: float = 8.3       # Unitary Cl- single-channel conductance (pS)
DEFAULT_TMEM_PERM_RATIO_HCO3_CL: float = 0.12  # Permeability / conductance ratio P_HCO3 / P_Cl
DEFAULT_TMEM_DENSITY_TOTAL: float = 2500.0  # Functional apical channels per airway cell
DEFAULT_TMEM_HILL_NH: float = 2.0           # Hill cooperativity coefficient
DEFAULT_TMEM_EC50_BASAL_UM: float = 0.40    # Half-maximal Ca2+ activation at 0 mV without potentiator (uM)
DEFAULT_TMEM_GATING_ZC: float = 0.35        # Effective gating charge for voltage-dependent Ca2+ sensitivity
DEFAULT_TMEM_POTENTIATOR_KD_NM: float = 50.0  # Potentiator binding affinity (nM)
DEFAULT_TMEM_POTENTIATOR_ALPHA: float = 0.75  # Maximal fractional reduction in apparent EC50 (shifts 0.4 uM -> 0.1 uM)

# SLC26A9 Biophysical Parameters
DEFAULT_SLC_GAMMA_CL_PS: float = 3.0        # Unitary Cl- single-channel conductance (pS)
DEFAULT_SLC_PERM_RATIO_HCO3_CL: float = 0.18  # Permeability / conductance ratio P_HCO3 / P_Cl
DEFAULT_SLC_PO_CONSTITUTIVE: float = 0.85   # Ca2+-independent constitutive open probability
DEFAULT_SLC_DENSITY_CFTR_NULL: float = 300.0  # Basal apical channels in untreated CFTR-null epithelia
DEFAULT_SLC_DENSITY_RESCUED: float = 2500.0   # Rescued apical channels under PDZ stabilization / trafficking
DEFAULT_SLC_PDZ_KD_UM: float = 1.2          # Binding affinity for apical PDZ scaffolding (uM)
DEFAULT_SLC_PDZ_MAX_STABILIZATION: float = 0.88 # Maximal apical retention fraction under full PDZ coupling

# ENaC Biophysical Parameters
DEFAULT_ENAC_GAMMA_NA_PS: float = 4.5       # Unitary Na+ conductance (pS)
DEFAULT_ENAC_DENSITY_CF: float = 3000.0     # Apical ENaC channel density in CF (uninhibited by CFTR)
DEFAULT_ENAC_PO_CF: float = 0.85            # Hyperactivated open probability in CFTR-deficient epithelia
DEFAULT_ENAC_DENSITY_HEALTHY: float = 1200.0  # Apical ENaC density under tonic CFTR inhibition
DEFAULT_ENAC_PO_HEALTHY: float = 0.35       # Normal ENaC open probability in healthy epithelia

# Wild-Type CFTR Reference Parameters
DEFAULT_CFTR_GAMMA_CL_PS: float = 8.0       # CFTR unitary Cl- conductance (pS)
DEFAULT_CFTR_PERM_RATIO_HCO3_CL: float = 0.25  # P_HCO3 / P_Cl for wild-type CFTR
DEFAULT_CFTR_DENSITY_WT: float = 1800.0     # Apical CFTR channels/cell in healthy airway ciliated cells
DEFAULT_CFTR_PO_PKA_ACTIVE: float = 0.45    # Phosphorylated open probability

# ASL Fluid & Acid-Base Parameters (Button et al. 2021)
ASL_HEIGHT_COLLAPSED_CF_UM: float = 3.2           # Collapsed CF baseline ASL height (< 3.5 um)
ASL_CLEARANCE_THRESHOLD_UM: float = 7.0           # Ciliary beat clearance threshold (>= 7.0 um)
ASL_TARGET_HEALTHY_UM: float = 8.0                # Healthy wild-type steady-state ASL height (um)
DEFAULT_K_HYDR_PER_HR: float = 0.35               # ASL fluid volume equilibration rate constant (1/h)

# Henderson-Hasselbalch Buffer Constants
CARBONIC_ACID_PKA: float = 6.10
CO2_SOLUBILITY_MM_PER_MMHG: float = 0.0307        # Henry's law constant for CO2 in water at 37 C
DEFAULT_PCO2_MMHG: float = 40.0                   # Physiological alveolar / tissue PCO2 (mmHg)
DEFAULT_BASAL_CF_HCO3_MM: float = 6.8             # Deficient ASL HCO3- in CFTR-null (pH ~6.84)
DEFAULT_NORMAL_AIRWAY_HCO3_MM: float = 23.5       # Restored ASL HCO3- in healthy / bypassed epithelia (pH ~7.38)


@dataclass(frozen=True)
class EpithelialElectrophysiology:
    """Cellular electrophysiological state and ionic conductances."""
    membrane_potential_mv: float
    g_cl_total_ps: float
    g_hco3_total_ps: float
    g_na_total_ps: float
    g_apical_total_ps: float
    i_cl_apical_pa: float
    i_hco3_apical_pa: float
    i_na_apical_pa: float
    i_net_apical_pa: float
    j_osm_solute_pmol_s_cm2: float


@dataclass(frozen=True)
class ASLHydrationState:
    """Airway surface liquid (ASL) height, fluid dynamics, and acid-base buffering."""
    asl_height_um: float
    ciliary_clearance_competent: bool
    net_water_flux_um_h: float
    asl_hco3_concentration_mm: float
    asl_ph: float
    ciliary_beat_frequency_hz: float


@dataclass(frozen=True)
class SimulationResult:
    """Comprehensive multi-scale biophysical output for an airway epithelial state."""
    scenario_name: str
    electrophysiology: EpithelialElectrophysiology
    asl_hydration: ASLHydrationState
    tmem16a_p_open: float
    slc26a9_effective_channels: float
    grounding: Dict[str, str]


def calculate_tmem16a_ec50(
    membrane_potential_mv: float,
    potentiator_conc_nm: float = 0.0,
    ec50_basal_um: float = DEFAULT_TMEM_EC50_BASAL_UM,
    gating_zc: float = DEFAULT_TMEM_GATING_ZC,
    potentiator_kd_nm: float = DEFAULT_TMEM_POTENTIATOR_KD_NM,
    potentiator_alpha: float = DEFAULT_TMEM_POTENTIATOR_ALPHA,
) -> float:
    """Calculate apparent EC50 for TMEM16A Ca2+ activation accounting for voltage and potentiator.

    Voltage shift: depolarizing membrane potentials (positive Vm) left-shift EC50:
      EC50(Vm) = EC50_basal * exp(-zc * F * Vm / (2 * R * T))
    Potentiator shift: allosteric potentiator further shifts apparent EC50 to lower concentrations:
      factor = 1.0 - alpha * ([pot] / (Kd + [pot]))
    """
    v_norm = (gating_zc * membrane_potential_mv) / (2.0 * THERMAL_VOLTAGE_MV)
    voltage_factor = math.exp(-v_norm)

    if potentiator_conc_nm > 0.0:
        pot_fraction = potentiator_conc_nm / (potentiator_kd_nm + potentiator_conc_nm)
        pot_factor = max(0.10, 1.0 - (potentiator_alpha * pot_fraction))
    else:
        pot_factor = 1.0

    return max(0.01, ec50_basal_um * voltage_factor * pot_factor)


def calculate_tmem16a_open_probability(
    cytosolic_ca_um: float,
    membrane_potential_mv: float,
    potentiator_conc_nm: float = 0.0,
    hill_nh: float = DEFAULT_TMEM_HILL_NH,
) -> float:
    """Compute steady-state open probability P_open of TMEM16A via allosteric Hill kinetics."""
    if cytosolic_ca_um <= 0.0:
        return 0.0

    ec50 = calculate_tmem16a_ec50(membrane_potential_mv, potentiator_conc_nm)
    ratio = cytosolic_ca_um / ec50
    ratio_pow = math.pow(ratio, hill_nh)
    return ratio_pow / (1.0 + ratio_pow)


def calculate_slc26a9_apical_fraction(
    pdz_scaffold_conc_um: float = 0.0,
    pdz_kd_um: float = DEFAULT_SLC_PDZ_KD_UM,
    max_stabilization: float = DEFAULT_SLC_PDZ_MAX_STABILIZATION,
    has_cftr_scaffold: bool = False,
) -> float:
    """Compute fraction of SLC26A9 retained at the apical membrane governed by PDZ scaffolding.

    In healthy airways, CFTR forms a macromolecular complex with NHERF1, anchoring SLC26A9.
    In CFTR-null cells, SLC26A9 is rapidly internalized unless exogenous PDZ stabilization
    or trafficking rescue is active.
    """
    basal_retention = 0.12  # Baseline retention fraction in unrescued CFTR-null epithelia
    if has_cftr_scaffold:
        return max_stabilization

    if pdz_scaffold_conc_um <= 0.0:
        return basal_retention

    pdz_occupancy = pdz_scaffold_conc_um / (pdz_kd_um + pdz_scaffold_conc_um)
    return basal_retention + (max_stabilization - basal_retention) * pdz_occupancy


def solve_epithelial_circuit(
    cytosolic_ca_um: float = 0.10,  # Resting cytosolic Ca2+ (100 nM = 0.10 uM)
    potentiator_conc_nm: float = 0.0,
    tmem_channel_density: float = DEFAULT_TMEM_DENSITY_TOTAL,
    slc_channel_density: float = DEFAULT_SLC_DENSITY_CFTR_NULL,
    pdz_scaffold_conc_um: float = 0.0,
    enac_density: float = DEFAULT_ENAC_DENSITY_CF,
    enac_p_open: float = DEFAULT_ENAC_PO_CF,
    cftr_density: float = 0.0,
    cftr_p_open: float = 0.0,
    has_cftr_scaffold: bool = False,
    e_cl_mv: float = DEFAULT_E_CL_MV,
    e_hco3_mv: float = DEFAULT_E_HCO3_MV,
    e_na_mv: float = DEFAULT_E_NA_MV,
    e_k_mv: float = DEFAULT_E_K_MV,
    e_basal_mv: float = DEFAULT_E_BASAL_MV,
    g_basal_ps: float = DEFAULT_G_BASAL_PS,
    g_k_apical_ps: float = DEFAULT_G_K_APICAL_PS,
) -> Tuple[EpithelialElectrophysiology, float, float]:
    """Solve Kirchhoff circuit equilibrium for the apical membrane potential Vm and conductances.

    Performs fixed-point iteration to account for the voltage-dependence of TMEM16A gating.
    Returns:
      (EpithelialElectrophysiology, tmem_p_open, slc_effective_channels)
    """
    # 1. SLC26A9 constitutive conductance
    slc_f_apical = calculate_slc26a9_apical_fraction(pdz_scaffold_conc_um, has_cftr_scaffold=has_cftr_scaffold)
    slc_eff_channels = slc_channel_density * slc_f_apical
    g_cl_slc = slc_eff_channels * DEFAULT_SLC_PO_CONSTITUTIVE * DEFAULT_SLC_GAMMA_CL_PS
    g_hco3_slc = g_cl_slc * DEFAULT_SLC_PERM_RATIO_HCO3_CL

    # 2. CFTR conductance (if present)
    g_cl_cftr = cftr_density * cftr_p_open * DEFAULT_CFTR_GAMMA_CL_PS
    g_hco3_cftr = g_cl_cftr * DEFAULT_CFTR_PERM_RATIO_HCO3_CL

    # 3. ENaC sodium conductance
    g_na_enac = enac_density * enac_p_open * DEFAULT_ENAC_GAMMA_NA_PS

    # Fixed-point iteration for Vm (converges within 3-4 steps)
    vm = -35.0  # Initial guess
    tmem_po = 0.0

    for _ in range(12):
        tmem_po = calculate_tmem16a_open_probability(cytosolic_ca_um, vm, potentiator_conc_nm)
        g_cl_tmem = tmem_channel_density * tmem_po * DEFAULT_TMEM_GAMMA_CL_PS
        g_hco3_tmem = g_cl_tmem * DEFAULT_TMEM_PERM_RATIO_HCO3_CL

        g_cl_tot = g_cl_slc + g_cl_cftr + g_cl_tmem
        g_hco3_tot = g_hco3_slc + g_hco3_cftr + g_hco3_tmem
        g_na_tot = g_na_enac

        # Total conductance and driving electromotive forces
        sum_g_e = (
            (g_cl_tot * e_cl_mv)
            + (g_hco3_tot * e_hco3_mv)
            + (g_na_tot * e_na_mv)
            + (g_k_apical_ps * e_k_mv)
            + (g_basal_ps * e_basal_mv)
        )
        sum_g = g_cl_tot + g_hco3_tot + g_na_tot + g_k_apical_ps + g_basal_ps

        vm_new = sum_g_e / max(1.0, sum_g)
        if abs(vm_new - vm) < 0.005:
            vm = vm_new
            break
        vm = 0.5 * (vm + vm_new)

    # Compute steady-state ionic currents (pA per cell)
    i_cl = g_cl_tot * (vm - e_cl_mv) * 1e-3
    i_hco3 = g_hco3_tot * (vm - e_hco3_mv) * 1e-3
    i_na = g_na_tot * (vm - e_na_mv) * 1e-3  # Negative for inward Na+ absorption
    i_net = i_cl + i_hco3 + i_na

    # Net solute secretion rate:
    # Anion conductance provides outward osmotic driving force, while ENaC sodium conductance
    # drives hyperabsorption. In CFTR-null, ENaC dominates, causing net negative flux (fluid absorption).
    total_anion_g = g_cl_tot + g_hco3_tot
    cf_baseline_g_anion = (
        (DEFAULT_SLC_DENSITY_CFTR_NULL * 0.12 * DEFAULT_SLC_PO_CONSTITUTIVE * DEFAULT_SLC_GAMMA_CL_PS * 1.18)
        + (DEFAULT_TMEM_DENSITY_TOTAL * 0.045 * DEFAULT_TMEM_GAMMA_CL_PS * 1.12)
    )  # ~1200 pS
    cf_baseline_g_na = DEFAULT_ENAC_DENSITY_CF * DEFAULT_ENAC_PO_CF * DEFAULT_ENAC_GAMMA_NA_PS  # 11475 pS

    anion_drive = (total_anion_g - cf_baseline_g_anion) / 300.0
    na_load = (g_na_tot / cf_baseline_g_na) * 20.0
    j_osm_pmol_s_cm2 = anion_drive - na_load

    g_apical_tot = g_cl_tot + g_hco3_tot + g_na_tot + g_k_apical_ps

    physio = EpithelialElectrophysiology(
        membrane_potential_mv=round(vm, 2),
        g_cl_total_ps=round(g_cl_tot, 1),
        g_hco3_total_ps=round(g_hco3_tot, 1),
        g_na_total_ps=round(g_na_tot, 1),
        g_apical_total_ps=round(g_apical_tot, 1),
        i_cl_apical_pa=round(i_cl, 3),
        i_hco3_apical_pa=round(i_hco3, 3),
        i_na_apical_pa=round(i_na, 3),
        i_net_apical_pa=round(i_net, 3),
        j_osm_solute_pmol_s_cm2=round(j_osm_pmol_s_cm2, 2),
    )
    return physio, tmem_po, slc_eff_channels


def simulate_asl_hydration(
    physio: EpithelialElectrophysiology,
    initial_asl_height_um: float = ASL_HEIGHT_COLLAPSED_CF_UM,
    simulation_hours: float = 24.0,
    k_hydr: float = DEFAULT_K_HYDR_PER_HR,
) -> ASLHydrationState:
    """Simulate dynamic 24-hour Airway Surface Liquid (ASL) height and pH buffering.

    1. ASL Target Height:
       - ENaC hyperabsorption suppression raises baseline from 3.2 um up to ~5.5 um.
       - Apical anion secretion (TMEM16A + SLC26A9) provides osmotic expansion up to 9.5 um.
    2. Dynamic differential relaxation:
       dh_ASL/dt = k_hydr * (h_target - h_ASL)
    3. Bicarbonate buffering:
       Calculates luminal [HCO3-] and Henderson-Hasselbalch pH.
    """
    cf_baseline_genac = DEFAULT_ENAC_DENSITY_CF * DEFAULT_ENAC_GAMMA_NA_PS * DEFAULT_ENAC_PO_CF
    enac_suppression = max(0.0, min(1.0, 1.0 - (physio.g_na_total_ps / cf_baseline_genac)))

    # Height baseline under ENaC modulation
    h_baseline_um = ASL_HEIGHT_COLLAPSED_CF_UM + (5.5 - ASL_HEIGHT_COLLAPSED_CF_UM) * enac_suppression

    # Apical anion osmotic expansion
    total_anion_g = physio.g_cl_total_ps + physio.g_hco3_total_ps
    cf_baseline_anion_g = 1200.0  # pS
    delta_anion_g = max(0.0, total_anion_g - cf_baseline_anion_g)
    k_anion_half = 4200.0  # pS for half-maximal osmotic hydration
    anion_gain = (9.5 - h_baseline_um) * (delta_anion_g / (k_anion_half + delta_anion_g))
    h_target = h_baseline_um + anion_gain

    # Numerical integration over simulation period
    h_asl = initial_asl_height_um
    dt = 0.25  # 15 minutes
    steps = int(simulation_hours / dt)
    for _ in range(steps):
        dh = k_hydr * (h_target - h_asl) * dt
        h_asl = max(1.5, min(12.0, h_asl + dh))

    # Ciliary beat competence requires h_ASL >= 7.0 um (Button et al. 2021)
    ciliary_competent = (h_asl >= ASL_CLEARANCE_THRESHOLD_UM)

    # Ciliary beat frequency (CBF): optimal 10-12 Hz in hydrated ASL, drops to < 3 Hz in collapsed CF mucus
    if h_asl >= 7.0:
        cbf_hz = min(12.0, 8.5 + (h_asl - 7.0) * 1.5)
    elif h_asl >= 4.5:
        cbf_hz = 5.0 + (h_asl - 4.5) * 1.4
    else:
        cbf_hz = max(1.0, 2.5 * (h_asl / 4.5))

    # Bicarbonate secretion & Henderson-Hasselbalch pH buffering
    # Basal CF: deficient HCO3- (6.8 mM) -> pH 6.84
    # Rescued: HCO3- rises toward 23.5 mM -> pH 7.38
    cf_baseline_hco3_g = 150.0  # pS
    delta_hco3_g = max(0.0, physio.g_hco3_total_ps - cf_baseline_hco3_g)
    k_hco3_half = 750.0  # pS/cell half-maximal bicarbonate saturation
    asl_hco3_mm = DEFAULT_BASAL_CF_HCO3_MM + (DEFAULT_NORMAL_AIRWAY_HCO3_MM - DEFAULT_BASAL_CF_HCO3_MM) * (
        delta_hco3_g / (k_hco3_half + delta_hco3_g)
    )

    # Henderson-Hasselbalch: pH = 6.10 + log10([HCO3-] / (s * PCO2))
    dissolved_co2_mm = CO2_SOLUBILITY_MM_PER_MMHG * DEFAULT_PCO2_MMHG  # ~1.228 mM
    asl_ph = CARBONIC_ACID_PKA + math.log10(asl_hco3_mm / dissolved_co2_mm)

    net_water_flux_rate = k_hydr * (h_target - h_asl)

    return ASLHydrationState(
        asl_height_um=round(h_asl, 2),
        ciliary_clearance_competent=ciliary_competent,
        net_water_flux_um_h=round(net_water_flux_rate, 3),
        asl_hco3_concentration_mm=round(asl_hco3_mm, 2),
        asl_ph=round(asl_ph, 2),
        ciliary_beat_frequency_hz=round(cbf_hz, 1),
    )


def simulate_scenario(
    scenario_name: str,
    cytosolic_ca_um: float = 0.10,
    potentiator_conc_nm: float = 0.0,
    pdz_scaffold_conc_um: float = 0.0,
    slc_channel_density: float = DEFAULT_SLC_DENSITY_CFTR_NULL,
    tmem_channel_density: float = DEFAULT_TMEM_DENSITY_TOTAL,
    enac_density: float = DEFAULT_ENAC_DENSITY_CF,
    enac_p_open: float = DEFAULT_ENAC_PO_CF,
    cftr_density: float = 0.0,
    cftr_p_open: float = 0.0,
    has_cftr_scaffold: bool = False,
    initial_asl_height_um: float = ASL_HEIGHT_COLLAPSED_CF_UM,
) -> SimulationResult:
    """Run full biophysical simulation for a named scenario."""
    physio, tmem_po, slc_eff = solve_epithelial_circuit(
        cytosolic_ca_um=cytosolic_ca_um,
        potentiator_conc_nm=potentiator_conc_nm,
        tmem_channel_density=tmem_channel_density,
        slc_channel_density=slc_channel_density,
        pdz_scaffold_conc_um=pdz_scaffold_conc_um,
        enac_density=enac_density,
        enac_p_open=enac_p_open,
        cftr_density=cftr_density,
        cftr_p_open=cftr_p_open,
        has_cftr_scaffold=has_cftr_scaffold,
    )

    asl = simulate_asl_hydration(physio, initial_asl_height_um=initial_asl_height_um)

    grounding = {
        "button_physiol_rev_2021": PRIMARY_DOI_BUTTON,
        "middleton_nejm_2019": PRIMARY_DOI_MIDDLETON,
        "vanzacaftor_phase3_trial": CLINICAL_TRIAL_VANZACAFTOR,
        "vx522_mrna_lnp_trial": CLINICAL_TRIAL_VX522,
        "pdb_cftr_structure": PDB_STRUCTURE_CFTR,
    }

    return SimulationResult(
        scenario_name=scenario_name,
        electrophysiology=physio,
        asl_hydration=asl,
        tmem16a_p_open=round(tmem_po, 4),
        slc26a9_effective_channels=round(slc_eff, 1),
        grounding=grounding,
    )


def benchmark_bypass_regimens() -> Dict[str, Any]:
    """Benchmark comparative biophysical efficacy across 6 standard airway scenarios:
    1. Healthy Wild-Type: Full CFTR conductance, tonic ENaC inhibition.
    2. Untreated CF Baseline: CFTR null, hyperactive ENaC, downregulated SLC26A9, basal Ca2+.
    3. TMEM16A Potentiation Alone: 100 nM small-molecule potentiator, resting Ca2+ (0.1 uM).
    4. SLC26A9 Scaffolding Rescue Alone: PDZ stabilization (2.5 uM) restoring surface SLC26A9.
    5. Dual-Channel Bypass Coupling: Combined TMEM16A potentiator (100 nM) + SLC26A9 PDZ rescue.
    6. Curative Triple Synergy: Dual-channel bypass + 60% catalytic ENaC inhibition.
    """
    scenarios = {
        "healthy_wild_type": simulate_scenario(
            "Healthy Wild-Type",
            cftr_density=DEFAULT_CFTR_DENSITY_WT,
            cftr_p_open=DEFAULT_CFTR_PO_PKA_ACTIVE,
            enac_density=DEFAULT_ENAC_DENSITY_HEALTHY,
            enac_p_open=DEFAULT_ENAC_PO_HEALTHY,
            slc_channel_density=DEFAULT_SLC_DENSITY_RESCUED,
            has_cftr_scaffold=True,
            initial_asl_height_um=ASL_TARGET_HEALTHY_UM,
        ),
        "untreated_cf_baseline": simulate_scenario(
            "Untreated CF Baseline",
            cftr_density=0.0,
            cftr_p_open=0.0,
            enac_density=DEFAULT_ENAC_DENSITY_CF,
            enac_p_open=DEFAULT_ENAC_PO_CF,
            slc_channel_density=DEFAULT_SLC_DENSITY_CFTR_NULL,
            pdz_scaffold_conc_um=0.0,
            potentiator_conc_nm=0.0,
            initial_asl_height_um=ASL_HEIGHT_COLLAPSED_CF_UM,
        ),
        "tmem16a_potentiated_alone": simulate_scenario(
            "TMEM16A Potentiated Alone",
            cftr_density=0.0,
            enac_density=DEFAULT_ENAC_DENSITY_CF,
            slc_channel_density=DEFAULT_SLC_DENSITY_CFTR_NULL,
            pdz_scaffold_conc_um=0.0,
            potentiator_conc_nm=100.0,
        ),
        "slc26a9_rescued_alone": simulate_scenario(
            "SLC26A9 Rescued Alone",
            cftr_density=0.0,
            enac_density=DEFAULT_ENAC_DENSITY_CF,
            slc_channel_density=DEFAULT_SLC_DENSITY_RESCUED,
            pdz_scaffold_conc_um=2.5,
            potentiator_conc_nm=0.0,
        ),
        "dual_channel_coupling": simulate_scenario(
            "Dual-Channel Bypass (TMEM16A + SLC26A9)",
            cftr_density=0.0,
            enac_density=DEFAULT_ENAC_DENSITY_CF,
            slc_channel_density=DEFAULT_SLC_DENSITY_RESCUED,
            pdz_scaffold_conc_um=2.5,
            potentiator_conc_nm=100.0,
        ),
        "curative_triple_synergy": simulate_scenario(
            "Curative Triple Synergy (Dual Bypass + ENaC Knockdown)",
            cftr_density=0.0,
            enac_density=DEFAULT_ENAC_DENSITY_CF * 0.40,  # 60% ENaC knockdown
            enac_p_open=DEFAULT_ENAC_PO_CF * 0.70,
            slc_channel_density=DEFAULT_SLC_DENSITY_RESCUED,
            pdz_scaffold_conc_um=2.5,
            potentiator_conc_nm=100.0,
        ),
    }

    results = {}
    for key, sc in scenarios.items():
        results[key] = {
            "name": sc.scenario_name,
            "membrane_potential_mv": sc.electrophysiology.membrane_potential_mv,
            "g_cl_total_ps": sc.electrophysiology.g_cl_total_ps,
            "g_hco3_total_ps": sc.electrophysiology.g_hco3_total_ps,
            "j_osm_pmol_s_cm2": sc.electrophysiology.j_osm_solute_pmol_s_cm2,
            "asl_height_um": sc.asl_hydration.asl_height_um,
            "ciliary_clearance_competent": sc.asl_hydration.ciliary_clearance_competent,
            "asl_ph": sc.asl_hydration.asl_ph,
            "cbf_hz": sc.asl_hydration.ciliary_beat_frequency_hz,
            "tmem16a_p_open": sc.tmem16a_p_open,
            "slc26a9_eff_channels": sc.slc26a9_effective_channels,
        }
    return results


def evaluate_falsification_negative_control() -> Dict[str, Any]:
    """Execute AH-10 dynamic falsification and negative control assertions:
    1. Untreated CF Negative Control: In absent CFTR and unrescued bypass channels,
       net solute flux must be negative (net absorption), and steady-state ASL height
       must collapse below ciliary clearance threshold (< 3.5 um).
    2. Zero-Calcium TMEM16A Negation: At [Ca2+]i = 0 uM without potentiator, TMEM16A
       open probability must equal 0.0, producing zero CaCC conductance.
    3. Dual-Channel Synergistic Rescue: TMEM16A potentiation (100 nM) combined with
       SLC26A9 PDZ stabilization must restore ASL height above the ciliary clearance
       threshold (>= 7.0 um) and physiological pH >= 7.15.
    """
    bench = benchmark_bypass_regimens()
    cf = bench["untreated_cf_baseline"]
    dual = bench["dual_channel_coupling"]

    # Assertion 1: CF collapsed baseline
    cf_collapsed = (cf["asl_height_um"] < 3.5) and (not cf["ciliary_clearance_competent"]) and (cf["j_osm_pmol_s_cm2"] < 0.0)

    # Assertion 2: Zero calcium negation
    zero_ca_po = calculate_tmem16a_open_probability(cytosolic_ca_um=0.0, membrane_potential_mv=-40.0)
    zero_ca_passed = (zero_ca_po == 0.0)

    # Assertion 3: Dual channel synergy
    dual_passed = (dual["asl_height_um"] >= ASL_CLEARANCE_THRESHOLD_UM) and dual["ciliary_clearance_competent"] and (dual["asl_ph"] >= 7.15)

    all_passed = cf_collapsed and zero_ca_passed and dual_passed
    return {
        "status": "PASS" if all_passed else "FAIL",
        "cf_collapsed_asl_um": cf["asl_height_um"],
        "cf_net_solute_flux_pmol_s_cm2": cf["j_osm_pmol_s_cm2"],
        "zero_ca_tmem16a_po": zero_ca_po,
        "dual_bypass_asl_height_um": dual["asl_height_um"],
        "dual_bypass_asl_ph": dual["asl_ph"],
        "ciliary_clearance_restored": dual["ciliary_clearance_competent"],
        "all_negative_controls_satisfied": all_passed,
    }


def main() -> None:
    """Command-line interface for TMEM16A potentiation and SLC26A9 coupling model."""
    parser = argparse.ArgumentParser(
        description="Thermodynamic and Conductance Modeling of Alternative Anion Secretion (T-0476)."
    )
    parser.add_argument("--ca-conc", type=float, default=0.10, help="Cytosolic [Ca2+] in uM (default: 0.10 uM / 100 nM)")
    parser.add_argument("--pot-conc", type=float, default=0.0, help="TMEM16A potentiator concentration in nM (default: 0.0)")
    parser.add_argument("--pdz-density", type=float, default=0.0, help="Apical PDZ scaffolding concentration in uM (default: 0.0)")
    parser.add_argument("--enac-knockdown", type=float, default=0.0, help="Fractional catalytic ENaC knockdown 0.0 to 0.8 (default: 0.0)")
    parser.add_argument("--benchmark", action="store_true", help="Execute comparative multi-scenario bypass benchmark")
    parser.add_argument("--falsification", action="store_true", help="Run AH-10 dynamic falsification negative controls")
    parser.add_argument("--json", action="store_true", help="Emit raw JSON output")

    args = parser.parse_args()

    if args.falsification:
        fals = evaluate_falsification_negative_control()
        print(json.dumps(fals, indent=2))
        return

    if args.benchmark:
        bench_data = benchmark_bypass_regimens()
        if args.json:
            print(json.dumps(bench_data, indent=2))
        else:
            print("=" * 85)
            print("ALTERNATIVE ANION SECRETION BYPASS BENCHMARK (TMEM16A & SLC26A9)")
            print("=" * 85)
            for k, v in bench_data.items():
                print(f"Scenario: {v['name']}")
                print(f"  - Membrane Potential:   {v['membrane_potential_mv']:>6.2f} mV")
                print(f"  - Total Cl- Conductance: {v['g_cl_total_ps']:>6.1f} pS/cell | HCO3-: {v['g_hco3_total_ps']:>5.1f} pS/cell")
                print(f"  - Net Solute Flux:      {v['j_osm_pmol_s_cm2']:>6.2f} pmol/(s*cm^2)")
                print(f"  - ASL Steady Height:    {v['asl_height_um']:>6.2f} um (Target >= 7.0 um)")
                print(f"  - Ciliary Competence:   {'RESTORED' if v['ciliary_clearance_competent'] else 'COLLAPSED'}")
                print(f"  - ASL pH & CBF:         pH {v['asl_ph']:.2f} | {v['cbf_hz']} Hz\n")
        return

    # Single simulation run
    sim_res = simulate_scenario(
        scenario_name="Custom User Parameter Sweep",
        cytosolic_ca_um=args.ca_conc,
        potentiator_conc_nm=args.pot_conc,
        pdz_scaffold_conc_um=args.pdz_density,
        slc_channel_density=DEFAULT_SLC_DENSITY_RESCUED if args.pdz_density > 0.5 else DEFAULT_SLC_DENSITY_CFTR_NULL,
        enac_density=DEFAULT_ENAC_DENSITY_CF * (1.0 - max(0.0, min(0.9, args.enac_knockdown))),
    )

    if args.json:
        print(json.dumps(asdict(sim_res), indent=2))
    else:
        print("=" * 80)
        print("TMEM16A & SLC26A9 CONDUCTANCE COUPLING SIMULATION (T-0476)")
        print("=" * 80)
        print(f"Parameters: [Ca2+]i = {args.ca_conc:.2f} uM | Potentiator = {args.pot_conc:.1f} nM | PDZ = {args.pdz_density:.2f} uM")
        print(f"  - Apical Vm:            {sim_res.electrophysiology.membrane_potential_mv} mV")
        print(f"  - G_Cl / G_HCO3:        {sim_res.electrophysiology.g_cl_total_ps} pS / {sim_res.electrophysiology.g_hco3_total_ps} pS")
        print(f"  - Net Solute Flux:      {sim_res.electrophysiology.j_osm_solute_pmol_s_cm2} pmol/(s*cm^2)")
        print(f"  - ASL Height:           {sim_res.asl_hydration.asl_height_um} um ({'COMPETENT' if sim_res.asl_hydration.ciliary_clearance_competent else 'SUB-THRESHOLD'})")
        print(f"  - ASL pH / CBF:         pH {sim_res.asl_hydration.asl_ph} | {sim_res.asl_hydration.ciliary_beat_frequency_hz} Hz")
        print("=" * 80)


if __name__ == "__main__":
    main()
