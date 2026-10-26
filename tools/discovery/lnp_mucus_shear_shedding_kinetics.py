#!/usr/bin/env python3
"""Biophysical Simulation of Aerosolized LNP Mucus Fluidic Shear & Charge-Neutralized PEG-Shedding Kinetics (T-0439).

Simulates vibrating mesh nebulization fluidic shear stresses, membrane lysis tension thresholds,
Ogston-Renkin steric mesh filtration, DLVO electrostatic double-layer interactions,
and acyl-chain desorption kinetics for deep airway basal stem cell penetration.

Primary Evidence Grounding:
  - Adv Mater (2012), DOI: 10.1002/adma.201201800
  - Adv Drug Deliv Rev (2009), DOI: 10.1016/j.addr.2008.11.002
  - Adv Drug Deliv Rev (2015), DOI: 10.1016/j.addr.2015.09.012
  - Science (2012), DOI: 10.1126/science.1223012
  - Physiol Rev (2022), DOI: 10.1152/physrev.00004.2021
  - Sci Transl Med (2021), DOI: 10.1126/scitranslmed.adw8886

Clinical Trials:
  - NCT05668741 (VX-522 Phase 1/2 inhaled mRNA LNP in CF)
  - NCT05248230 (4D-710 Phase 1/2 aerosolized AAV in CF)

Structures:
  - RCSB PDB 6MSM, 5UAK
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Dict, List, Tuple

# Mandatory ClinicalTrials.gov Attribution
NLM_ATTRIBUTION = (
    "ClinicalTrials.gov does not evaluate safety or efficacy; "
    "listing on the registry does not denote government endorsement."
)

CLINICAL_TRIAL_VX522 = "NCT05668741"
CLINICAL_TRIAL_4D710 = "NCT05248230"
PDB_DEPHOSPHO_CFTR = "6MSM"
PDB_ACTIVE_CFTR = "5UAK"

PRIMARY_DOI_AEROSOL_MRNA = "10.1002/adma.201201800"
PRIMARY_DOI_MUCUS_PENETRATION = "10.1016/j.addr.2008.11.002"
PRIMARY_DOI_PEG_SHEDDING = "10.1016/j.addr.2015.09.012"
PRIMARY_DOI_PERICILIARY_BRUSH = "10.1126/science.1223012"
PRIMARY_DOI_AIRWAY_MUCUS_PHYS = "10.1152/physrev.00004.2021"
PRIMARY_DOI_GENE_EDITING_LUNG = "10.1126/scitranslmed.adw8886"

KB_J_PER_K = 1.380649e-23  # Boltzmann constant (J/K)
T_KELVIN = 310.15          # Body temperature (37 C in K)
KB_T = KB_J_PER_K * T_KELVIN


@dataclass(frozen=True)
class NebulizerShearParameters:
    """Parameters for vibrating mesh nebulizer shear and LNP membrane stability."""
    fluid_viscosity_pa_s: float = 0.001       # Aqueous suspension viscosity (Pa*s)
    critical_lysis_tension_mn_m: float = 20.0 # Bilayer critical rupture tension (mN/m)
    activation_area_nm2: float = 0.40         # Area expansion per lipid at transition (nm^2)
    k0_destab_s1: float = 1.0e-4              # Baseline thermal destabilization rate (1/s)


def calculate_nebulizer_shear_stress(
    shear_rate_s1: float,
    viscosity_pa_s: float = 0.001,
) -> float:
    """Compute local fluidic shear stress tau = eta * gamma_dot (Pa)."""
    if shear_rate_s1 < 0.0:
        raise ValueError("Shear rate must be non-negative.")
    return viscosity_pa_s * shear_rate_s1


def calculate_membrane_tension(
    radius_lnp_nm: float,
    shear_stress_pa: float,
    extensional_amplification: float = 1000.0,
) -> float:
    """Compute bilayer tension sigma = 0.5 * r * tau_ext (mN/m).
    
    Includes extensional flow acceleration / converging aperture stress concentration
    (amplification ~ 1000x in vibrating mesh micro-pores).
    """
    if radius_lnp_nm <= 0.0:
        raise ValueError("LNP radius must be positive.")
    # r in meters = radius_lnp_nm * 1e-9
    # tau_ext in Pa = shear_stress_pa * extensional_amplification
    # tension in N/m = 0.5 * (radius_lnp_nm * 1e-9) * (shear_stress_pa * extensional_amplification)
    # tension in mN/m = tension_N_m * 1000
    tension_mn_m = 0.5 * (radius_lnp_nm * 1e-9) * (shear_stress_pa * extensional_amplification) * 1000.0
    return tension_mn_m


def evaluate_lnp_nebulization_survival(
    radius_lnp_nm: float,
    shear_rate_s1: float,
    params: NebulizerShearParameters | None = None,
) -> Dict[str, Any]:
    """Evaluate LNP structural integrity under vibrating mesh nebulizer shear."""
    p = params or NebulizerShearParameters()
    tau_pa = calculate_nebulizer_shear_stress(shear_rate_s1, p.fluid_viscosity_pa_s)
    sigma_mn_m = calculate_membrane_tension(radius_lnp_nm, tau_pa)

    # Tension ratio relative to critical lysis threshold
    tension_ratio = sigma_mn_m / p.critical_lysis_tension_mn_m

    # Mechanical rupture when membrane tension approaches or exceeds critical lysis threshold
    if sigma_mn_m >= p.critical_lysis_tension_mn_m:
        excess = (sigma_mn_m - p.critical_lysis_tension_mn_m) / p.critical_lysis_tension_mn_m
        survival_fraction = math.exp(-5.0 * excess)
        k_destab = 1.0e5 * (1.0 + excess)
    else:
        # Sub-critical tension: thermal Eyring fluctuation rate
        delta_sigma_n_m = (sigma_mn_m - p.critical_lysis_tension_mn_m) * 1e-3
        area_m2 = p.activation_area_nm2 * 1e-18
        work_j = delta_sigma_n_m * area_m2
        exponent = max(-30.0, min(30.0, work_j / KB_T))
        k_destab = p.k0_destab_s1 * math.exp(exponent)
        t_passage_s = 1.0e-5
        rupture_prob = 1.0 - math.exp(-k_destab * t_passage_s)
        survival_fraction = max(0.0, min(1.0, 1.0 - rupture_prob))

    return {
        "shear_rate_s1": shear_rate_s1,
        "shear_stress_pa": tau_pa,
        "membrane_tension_mn_m": sigma_mn_m,
        "tension_ratio": tension_ratio,
        "k_destabilization_s1": k_destab,
        "survival_fraction": survival_fraction,
        "verdict": "INTACT" if survival_fraction >= 0.90 else "LYSED",
    }


def calculate_ogston_renkin_diffusivity(
    diameter_lnp_nm: float,
    mesh_size_nm: float = 120.0,
    fiber_radius_nm: float = 3.5,
) -> float:
    """Compute effective steric diffusivity ratio D_steric / D_0 using Ogston-Renkin model."""
    r_lnp = diameter_lnp_nm / 2.0
    if r_lnp <= 0.0 or mesh_size_nm <= 0.0:
        return 0.0

    # Ogston steric partitioning: exp(-( (r_lnp + r_fiber) / mesh_size )^2 )
    steric_exponent = ((r_lnp + fiber_radius_nm) / mesh_size_nm) ** 2
    k_ogston = math.exp(-steric_exponent)

    # Hydrodynamic Renkin drag factor lambda = r_lnp / mesh_size
    lam = min(1.0, r_lnp / mesh_size_nm)
    f_renkin = (1.0 - lam) ** 2 * (
        1.0 - 2.104 * lam + 2.09 * (lam ** 3) - 0.95 * (lam ** 5)
    )
    f_renkin = max(0.0, f_renkin)

    return max(0.0, min(1.0, k_ogston * f_renkin))


def calculate_dlvo_electrostatic_factor(
    zeta_potential_mv: float,
    mucus_potential_mv: float = -15.0,
    debye_length_nm: float = 0.80,
) -> float:
    """Compute electrostatic hindrance/attraction transmission factor E_elec in [0.0, 1.0]."""
    # Trapping energy in units of k_B * T
    # Near-neutral PEGylation (-2 to +2 mV) avoids adhesive or repulsive trapping
    # Strongly cationic (> +5 mV) binds negatively charged sialylated mucin (mucus_potential ~ -15 mV)
    delta_zeta = abs(zeta_potential_mv)
    # Phenomenological electrostatic trapping potential
    u_trap_kbt = 0.15 * (zeta_potential_mv - 0.0) ** 2
    if zeta_potential_mv > 0.0:
        # Cationic adhesive cross-linking with negative mucin fibers
        u_trap_kbt *= 1.8

    e_factor = 1.0 / (1.0 + u_trap_kbt)
    return max(0.0, min(1.0, e_factor))


@dataclass(frozen=True)
class PEGSheddingKinetics:
    """PEG-lipid shedding parameters by lipid anchor composition."""
    anchor_type: str = "DMG_PEG2000"     # DMG (C14) vs DSPE (C18)
    k_shed_per_hour: float = 0.30        # Shedding rate (1/hour)
    initial_peg_mol_pct: float = 2.0     # Initial surface PEGylation (mol%)
    critical_endocytic_peg_pct: float = 0.50 # Threshold for cellular uptake

    def peg_density_at_time(self, t_hours: float) -> float:
        """Compute remaining surface PEG mol% at time t."""
        if t_hours <= 0.0:
            return self.initial_peg_mol_pct
        return self.initial_peg_mol_pct * math.exp(-self.k_shed_per_hour * t_hours)

    def is_endocytosis_competent(self, t_hours: float) -> bool:
        """Check if PEG has shed sufficiently for epithelial/stem cell endocytosis."""
        return self.peg_density_at_time(t_hours) <= self.critical_endocytic_peg_pct


def simulate_mucus_transit_and_basal_delivery(
    diameter_lnp_nm: float = 65.0,
    zeta_potential_mv: float = -1.2,
    anchor_type: str = "DMG_PEG2000",
    k_shed_per_hour: float = 0.30,
    mucus_thickness_um: float = 25.0,
    pcl_thickness_um: float = 7.0,
    mesh_size_nm: float = 120.0,
    duration_hours: float = 4.0,
    dt_hours: float = 0.05,
) -> Dict[str, Any]:
    """Simulate coupled transient transport through mucus and PCL to the basal stem cell niche."""
    # Steric diffusivity
    d_steric_ratio = calculate_ogston_renkin_diffusivity(diameter_lnp_nm, mesh_size_nm)
    e_elec = calculate_dlvo_electrostatic_factor(zeta_potential_mv)
    d_eff_ratio = d_steric_ratio * e_elec

    # Water diffusivity D_0 (Stokes-Einstein) for d=65nm in m^2/s
    # D0 = kBT / (6*pi*eta*r) = (4.28e-21) / (6 * pi * 0.001 * 32.5e-9) = 7.0e-12 m^2/s = 2.52e4 um^2/h
    d0_um2_h = (KB_T / (6.0 * math.pi * 0.001 * (diameter_lnp_nm * 0.5 * 1e-9))) * 1e12 * 3600.0
    d_mucus_um2_h = d0_um2_h * d_eff_ratio

    # Mucociliary clearance constant: k_mcc ~ 0.5 / hour (half-life ~ 1.4 hours)
    k_mcc_per_hour = 0.50

    shedder = PEGSheddingKinetics(
        anchor_type=anchor_type,
        k_shed_per_hour=k_shed_per_hour,
        initial_peg_mol_pct=2.0,
        critical_endocytic_peg_pct=0.50,
    )

    total_depth_um = mucus_thickness_um + pcl_thickness_um
    # Effective transit timescale tau_diff = L^2 / (2 * D)
    tau_transit_h = (total_depth_um ** 2) / (2.0 * max(d_mucus_um2_h, 1.0))

    # Fractional breakthrough reaching basal layer before mucociliary clearance
    # Breakthrough ~ exp(-sqrt(k_mcc * L^2 / D))
    penetration_factor = math.exp(-math.sqrt(k_mcc_per_hour * (total_depth_um ** 2) / max(d_mucus_um2_h, 1.0)))

    # Evaluate endocytic competence upon basal arrival
    endocytosis_competent = shedder.is_endocytosis_competent(min(duration_hours, tau_transit_h))

    return {
        "diameter_lnp_nm": diameter_lnp_nm,
        "zeta_potential_mv": zeta_potential_mv,
        "steric_diffusivity_ratio": d_steric_ratio,
        "electrostatic_factor": e_elec,
        "effective_diffusivity_ratio": d_eff_ratio,
        "transit_half_time_hours": tau_transit_h,
        "delivered_fraction_at_basal": penetration_factor,
        "endocytosis_competent_at_arrival": endocytosis_competent,
        "curative_delivery_feasible": penetration_factor >= 0.10 and d_eff_ratio >= 0.05,
    }
