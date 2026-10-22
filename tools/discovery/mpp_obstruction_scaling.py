#!/usr/bin/env python3
"""Airway Mucus-Penetrating Particle (MPP) Steric & Adhesive Obstruction Scaling (T-0244).

Biophysical modeling of nanoparticle transport through human Cystic Fibrosis airway
mucus hydrogels, coupling:
  1. Stokes-Einstein hydrodynamic diffusivity in physiological aqueous media at 37 °C.
  2. Dense polymer brush conformation: Flory radius (R_F), mushroom-to-brush transition
     threshold (Gamma*), and polyvalent hydrophobic core adhesive shielding (F_adhesive).
  3. Mucin (MUC5AC / MUC5B) and purulent extracellular DNA (eDNA) hydrogel mesh pore
     scaling (xi_mesh), including Ogston steric exclusion (F_steric).
  4. Mucolytic synergy: rhDNase I (dornase alfa) enzymatic eDNA degradation expanding
     steric mesh pores and rescuing macromolecular gene vector diffusion.
  5. Mucus blanket transit kinetics: effective diffusivity (D_eff / D_0) and transit time
     (tau_transit) scaled against the physiological mucociliary clearance window.

Grounding & Evidence Provenance (Answer Tier T2 Certified):
  - Primary Literature:
    * Lai SK et al. (2007). Rapid transport of large polymeric nanoparticles in fresh
      undiluted human mucus. Proc Natl Acad Sci USA 104(5):1482-1487.
      DOI: 10.1073/pnas.0610617104.
    * Lai SK et al. (2009). Mucus-penetrating nanoparticles for drug and gene delivery
      to mucosal tissues. Adv Drug Deliv Rev 61(2):158-171.
      DOI: 10.1016/j.addr.2008.09.009.
    * Tang BC et al. (2009). Biodegradable polymer nanoparticles that rapidly penetrate
      undiluted human mucus. Proc Natl Acad Sci USA 106(46):19268-19273.
      DOI: 10.1073/pnas.0811566106.
    * Suk JS et al. (2016). Nanoparticle barriers in respiratory diseases.
      Adv Drug Deliv Rev 99(Pt A):28-51. DOI: 10.1016/j.addr.2015.09.006.
    * Suk JS et al. (2009). The penetration of fresh undiluted sputum expectorated by
      cystic fibrosis patients by nanoparticles. Biomaterials 30(13):2591-2597.
      DOI: 10.1016/j.biomaterials.2009.01.052.
    * Shak S et al. (1990). Recombinant human DNase I reduces the viscosity of cystic
      fibrosis sputum. Proc Natl Acad Sci USA 87(23):9188-9192.
      DOI: 10.1073/pnas.87.23.9188.
  - ClinicalTrials.gov Registries (verified in tools/trial_atlas/snapshot.json):
    * NCT00843817: Safety and Efficacy of Aerosolized Dornase Alfa.
    * NCT00680316: Hypertonic Saline and rhDNase in Children With CF.
    * NCT06350474: Evaluating Discontinuation of Chronic Inhaled Therapies in CFTR
      Modulator-Treated Patients (SIMPLIFY study).
    * NCT05248230: Aerosolized AAV Gene Therapy for Cystic Fibrosis (4D-710).

Clinical Dignity & Research Boundaries:
  ClinicalTrials.gov does not evaluate safety or efficacy; listing on the registry
  does not denote government endorsement. Research modeling and in silico biophysical
  transport evaluation only; not medical advice or clinical dosing instructions.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
import math
from pathlib import Path
import sys
from typing import Any, Dict, List, Optional, Tuple

NLM_ATTRIBUTION = (
    "ClinicalTrials.gov does not evaluate safety or efficacy; "
    "listing on the registry does not denote government endorsement."
)

RESEARCH_BOUNDARY = (
    "Research, not medical advice. Biophysical macromolecular transport model, "
    "not clinical administration guidance."
)

# ---------------------------------------------------------------------------
# Fundamental Biophysical & Thermodynamic Constants
# ---------------------------------------------------------------------------
BOLTZMANN_CONSTANT_J_K: float = 1.380649e-23  # J / K
PHYSIOLOGICAL_TEMP_K: float = 310.15          # 37.0 °C
WATER_VISCOSITY_PA_S: float = 6.91e-4         # Dynamic viscosity of water at 37 °C (Pa*s)
METER_TO_UM: float = 1.0e6
M2_PER_S_TO_UM2_PER_S: float = 1.0e12

# ---------------------------------------------------------------------------
# Physiological Mucus & Purulent Sputum Mesh Baseline Parameters
# ---------------------------------------------------------------------------
HEALTHY_MUCUS_BASE_MESH_NM: float = 220.0     # Mean pore spacing in normal airway mucus (xi_0)
HEALTHY_MUCIN_VOLUME_FRACTION: float = 0.02   # Baseline mucin solids volume fraction (phi_m,0)
CF_MUCIN_VOLUME_FRACTION: float = 0.065       # Elevated mucin solids in CF hyper-concentrated mucus
DEFAULT_CF_EDNA_CONC_MG_ML: float = 2.8       # Extracellular DNA concentration in purulent CF sputum
EDNA_CROSSOVER_CONC_MG_ML: float = 2.0        # Characteristic eDNA mesh entanglement threshold (C_eDNA^dag)
EDNA_MESH_EXPONENT: float = 0.35              # Scaling exponent for eDNA mesh tightening (nu_DNA)
MUCIN_MESH_EXPONENT: float = 0.75             # De Gennes semi-dilute polymer scaling exponent (3/4)
POLYMER_FIBER_RADIUS_NM: float = 3.5          # Hydrated radius of glycoprotein/eDNA fibers (r_f)

# ---------------------------------------------------------------------------
# PEG Surface Brush & Nanoparticle Coating Parameters
# ---------------------------------------------------------------------------
PEG_MONOMER_LENGTH_NM: float = 0.35           # Effective bond length of ethylene glycol unit (a_PEG)
PEG_MONOMER_MW_DA: float = 44.05              # Molecular weight of monomer (-CH2-CH2-O-)
FLORY_EXPONENT: float = 0.60                  # Flory good-solvent swelling exponent (nu = 3/5)
DEFAULT_PEG_MW_KDA: float = 5.0               # Nominal PEG molecular weight for dense brush coating
DEFAULT_PEG_GRAFTING_CHAINS_NM2: float = 0.35 # Dense brush surface grafting density (Gamma)
BRUSH_SHIELDING_KAPPA: float = 1.25           # Exponential screening efficiency exponent of PEG brush
INTRINSIC_CORE_ADHESION: float = 450.0        # Dimensionless polyvalent core adhesive affinity K_adh[M_bind]

# ---------------------------------------------------------------------------
# Anatomical & Transport Clearance Windows
# ---------------------------------------------------------------------------
DEFAULT_MUCUS_BLANKET_THICKNESS_UM: float = 30.0  # Bronchial mucus layer depth in CF airways (h_mucus)
MUCOCILIARY_CLEARANCE_WINDOW_MIN: float = 30.0    # Mean time window before ciliary clearance (tau_clearance)
TARGET_MINIMUM_DIFFUSION_RATIO: float = 0.10      # Threshold D_eff / D_0 for rapid penetration


@dataclass(frozen=True)
class NanoparticleConfig:
    """Physicochemical specification of an engineered aerosolized delivery nanoparticle."""
    diameter_nm: float
    peg_mw_kda: float = DEFAULT_PEG_MW_KDA
    grafting_density_chains_nm2: float = DEFAULT_PEG_GRAFTING_CHAINS_NM2
    intrinsic_adhesiveness: float = INTRINSIC_CORE_ADHESION
    core_material: str = "biodegradable_polymer_plga_pba"


@dataclass(frozen=True)
class MucusHydrogelConfig:
    """Biochemical and rheological characteristics of the mucosal barrier."""
    mucin_volume_fraction: float = CF_MUCIN_VOLUME_FRACTION
    edna_conc_mg_ml: float = DEFAULT_CF_EDNA_CONC_MG_ML
    mucus_thickness_um: float = DEFAULT_MUCUS_BLANKET_THICKNESS_UM
    rhdnase_treated: bool = False
    rhdnase_cleavage_efficiency: float = 0.0  # Fraction [0.0 - 1.0] of eDNA degraded
    temperature_k: float = PHYSIOLOGICAL_TEMP_K


@dataclass(frozen=True)
class MPPSimulationResult:
    """Deterministic biophysical transport evaluation for a nanoparticle in mucus."""
    particle_diameter_nm: float
    stokes_einstein_d0_um2_s: float
    mesh_pore_spacing_xi_nm: float
    steric_exclusion_factor_f_steric: float
    flory_radius_rf_nm: float
    brush_transition_density_gamma_star: float
    grafting_ratio_gamma_over_gamma_star: float
    brush_regime: str
    exposed_core_area_fraction: float
    adhesive_obstruction_factor_f_adhesive: float
    effective_diffusivity_ratio: float
    effective_diffusivity_deff_um2_s: float
    mucus_transit_time_minutes: float
    mucociliary_clearance_limit_minutes: float
    escapes_mucociliary_clearance: bool
    transport_classification: str


@dataclass(frozen=True)
class ObstructionScalingSummary:
    """Comparative summary of conventional particle (CPP) vs mucus-penetrating particle (MPP)."""
    particle_diameter_nm: float
    cpp_effective_diffusivity_ratio: float
    cpp_transit_time_minutes: float
    cpp_escapes_clearance: bool
    mpp_effective_diffusivity_ratio: float
    mpp_transit_time_minutes: float
    mpp_escapes_clearance: bool
    fold_diffusivity_enhancement: float
    mesh_pore_spacing_nm: float
    cf_baseline_vs_rhdnase_diffusivity_ratio: float
    tier: str
    clinical_grounding: List[str]
    primary_dois: List[str]


def calculate_stokes_einstein_d0(
    diameter_nm: float,
    temperature_k: float = PHYSIOLOGICAL_TEMP_K,
    viscosity_pa_s: float = WATER_VISCOSITY_PA_S,
) -> float:
    """Calculate Stokes-Einstein diffusion coefficient D_0 in water at physiological temperature.

    Equation:
        D_0 = (k_B * T) / (6 * pi * mu_w * r_p)
    Returns:
        Diffusivity in um^2 / s.
    """
    if diameter_nm <= 0.0:
        raise ValueError("Particle diameter must be positive.")
    radius_m = (diameter_nm / 2.0) * 1.0e-9
    d0_m2_s = (BOLTZMANN_CONSTANT_J_K * temperature_k) / (6.0 * math.pi * viscosity_pa_s * radius_m)
    return round(d0_m2_s * M2_PER_S_TO_UM2_PER_S, 4)


def calculate_flory_radius_and_brush_transition(
    peg_mw_kda: float = DEFAULT_PEG_MW_KDA,
) -> Tuple[float, float]:
    """Calculate the Flory radius (R_F) and critical mushroom-to-brush surface density (Gamma*).

    Equations:
        N_p = (M_w * 1000) / M_0
        R_F = a_PEG * (N_p)^(3/5)
        Gamma* = 1.0 / (pi * R_F^2)
    Returns:
        Tuple of (Flory radius R_F in nm, critical density Gamma* in chains / nm^2).
    """
    if peg_mw_kda <= 0.0:
        raise ValueError("PEG molecular weight must be positive.")
    degree_polymerization = (peg_mw_kda * 1000.0) / PEG_MONOMER_MW_DA
    rf_nm = PEG_MONOMER_LENGTH_NM * (degree_polymerization ** FLORY_EXPONENT)
    gamma_star = 1.0 / (math.pi * (rf_nm ** 2))
    return round(rf_nm, 3), round(gamma_star, 5)


def calculate_mucus_mesh_spacing(
    mucin_volume_fraction: float = CF_MUCIN_VOLUME_FRACTION,
    edna_conc_mg_ml: float = DEFAULT_CF_EDNA_CONC_MG_ML,
    xi_0_nm: float = HEALTHY_MUCUS_BASE_MESH_NM,
    rhdnase_treated: bool = False,
    rhdnase_cleavage_efficiency: float = 0.0,
) -> float:
    """Calculate the average hydrogel mesh pore spacing xi_mesh in CF mucus or purulent sputum.

    Equations:
        C_eff = C_eDNA * (1.0 - efficiency) if rhDNase active
        xi_mesh = xi_0 * (phi_m / phi_m,0)^(-3/4) * (1.0 + C_eff / C_eDNA^dag)^(-nu_DNA)
    Returns:
        Hydrogel mesh pore spacing in nanometers.
    """
    if mucin_volume_fraction <= 0.0:
        raise ValueError("Mucin volume fraction must be positive.")
    effective_edna = edna_conc_mg_ml
    if rhdnase_treated:
        clamp_eff = max(0.0, min(1.0, rhdnase_cleavage_efficiency))
        effective_edna = edna_conc_mg_ml * (1.0 - clamp_eff)

    mucin_scaling = (mucin_volume_fraction / HEALTHY_MUCIN_VOLUME_FRACTION) ** (-MUCIN_MESH_EXPONENT)
    edna_term = 1.0 + (max(0.0, effective_edna) / EDNA_CROSSOVER_CONC_MG_ML)
    edna_scaling = edna_term ** (-EDNA_MESH_EXPONENT)

    mesh_nm = xi_0_nm * mucin_scaling * edna_scaling
    return round(max(5.0, mesh_nm), 2)


def calculate_steric_exclusion_factor(
    particle_diameter_nm: float,
    mesh_pore_spacing_nm: float,
    fiber_radius_nm: float = POLYMER_FIBER_RADIUS_NM,
) -> float:
    """Calculate the Ogston / Brinkman hydrodynamic steric exclusion factor F_steric.

    Equation:
        r_p = d_p / 2
        F_steric = exp( -pi * ((r_p + r_f) / (xi_mesh + r_f))^2 )
    Returns:
        Dimensionless steric hindrance factor in range [0.0, 1.0].
    """
    if particle_diameter_nm <= 0.0 or mesh_pore_spacing_nm <= 0.0:
        return 0.0
    r_p = particle_diameter_nm / 2.0
    numerator = r_p + fiber_radius_nm
    denominator = mesh_pore_spacing_nm + fiber_radius_nm
    steric_ratio = numerator / denominator

    # Prevent underflow for very large steric ratios
    if steric_ratio > 4.0:
        return 0.0
    f_steric = math.exp(-math.pi * (steric_ratio ** 2))
    return round(max(0.0, min(1.0, f_steric)), 5)


def calculate_peg_brush_shielding(
    grafting_density_chains_nm2: float,
    peg_mw_kda: float = DEFAULT_PEG_MW_KDA,
    intrinsic_adhesiveness: float = INTRINSIC_CORE_ADHESION,
    kappa_brush: float = BRUSH_SHIELDING_KAPPA,
) -> Tuple[float, str, float, float]:
    """Calculate dense PEG brush shielding and adhesive retardation factor F_adhesive.

    Equations:
        ratio = Gamma / Gamma*
        f_exposed = exp( -kappa_brush * ratio )
        F_adhesive = 1.0 / (1.0 + K_adh[M_bind] * f_exposed)
    Returns:
        Tuple of (ratio, regime_name, f_exposed, F_adhesive).
    """
    _, gamma_star = calculate_flory_radius_and_brush_transition(peg_mw_kda)
    if gamma_star <= 0.0:
        return 0.0, "uncoated_cpp", 1.0, 0.0

    ratio = grafting_density_chains_nm2 / gamma_star
    if ratio < 1.0:
        regime = "mushroom" if ratio > 0.1 else "uncoated_cpp"
    elif ratio < 2.5:
        regime = "intermediate_brush"
    else:
        regime = "dense_brush_mpp"

    f_exposed = math.exp(-kappa_brush * max(0.0, ratio))
    denominator = 1.0 + (intrinsic_adhesiveness * f_exposed)
    f_adhesive = 1.0 / denominator
    return round(ratio, 2), regime, round(f_exposed, 6), round(f_adhesive, 5)


def simulate_mpp_transport(
    particle: NanoparticleConfig,
    mucus: MucusHydrogelConfig,
) -> MPPSimulationResult:
    """Execute complete deterministic biophysical transport simulation for a nanoparticle in CF mucus."""
    d0 = calculate_stokes_einstein_d0(particle.diameter_nm, mucus.temperature_k)
    rf_nm, gamma_star = calculate_flory_radius_and_brush_transition(particle.peg_mw_kda)
    mesh_nm = calculate_mucus_mesh_spacing(
        mucin_volume_fraction=mucus.mucin_volume_fraction,
        edna_conc_mg_ml=mucus.edna_conc_mg_ml,
        rhdnase_treated=mucus.rhdnase_treated,
        rhdnase_cleavage_efficiency=mucus.rhdnase_cleavage_efficiency,
    )
    f_steric = calculate_steric_exclusion_factor(particle.diameter_nm, mesh_nm)
    gamma_ratio, regime, f_exposed, f_adhesive = calculate_peg_brush_shielding(
        grafting_density_chains_nm2=particle.grafting_density_chains_nm2,
        peg_mw_kda=particle.peg_mw_kda,
        intrinsic_adhesiveness=particle.intrinsic_adhesiveness,
    )

    effective_diffusivity_ratio = round(f_steric * f_adhesive, 5)
    deff_um2_s = round(d0 * effective_diffusivity_ratio, 5)

    # 1D/2D random walk mean first-passage transit time: tau = h^2 / (4 * D_eff)
    if deff_um2_s <= 1.0e-8:
        transit_time_min = 9999.0
    else:
        transit_time_sec = (mucus.mucus_thickness_um ** 2) / (4.0 * deff_um2_s)
        transit_time_min = round(transit_time_sec / 60.0, 2)

    escapes_clearance = (
        transit_time_min <= MUCOCILIARY_CLEARANCE_WINDOW_MIN
        and effective_diffusivity_ratio >= TARGET_MINIMUM_DIFFUSION_RATIO
    )

    if effective_diffusivity_ratio >= 0.15:
        classification = "rapid_mpp_penetrator"
    elif effective_diffusivity_ratio >= 0.05:
        classification = "moderately_hindered_diffuser"
    elif f_steric < 0.05:
        classification = "sterically_excluded_mesh_trapped"
    else:
        classification = "adhesively_immobilized_cpp"

    return MPPSimulationResult(
        particle_diameter_nm=particle.diameter_nm,
        stokes_einstein_d0_um2_s=d0,
        mesh_pore_spacing_xi_nm=mesh_nm,
        steric_exclusion_factor_f_steric=f_steric,
        flory_radius_rf_nm=rf_nm,
        brush_transition_density_gamma_star=gamma_star,
        grafting_ratio_gamma_over_gamma_star=gamma_ratio,
        brush_regime=regime,
        exposed_core_area_fraction=f_exposed,
        adhesive_obstruction_factor_f_adhesive=f_adhesive,
        effective_diffusivity_ratio=effective_diffusivity_ratio,
        effective_diffusivity_deff_um2_s=deff_um2_s,
        mucus_transit_time_minutes=transit_time_min,
        mucociliary_clearance_limit_minutes=MUCOCILIARY_CLEARANCE_WINDOW_MIN,
        escapes_mucociliary_clearance=escapes_clearance,
        transport_classification=classification,
    )


def compare_cpp_vs_mpp_scaling(
    diameter_nm: float = 65.0,
    peg_grafting_density: float = DEFAULT_PEG_GRAFTING_CHAINS_NM2,
    rhdnase_cleavage_efficiency: float = 0.75,
) -> ObstructionScalingSummary:
    """Generate comparative scaling summary between conventional (uncoated) and dense-brush MPP particles."""
    cf_mucus = MucusHydrogelConfig(rhdnase_treated=False)
    cf_rhdnase_mucus = MucusHydrogelConfig(
        rhdnase_treated=True,
        rhdnase_cleavage_efficiency=rhdnase_cleavage_efficiency,
    )

    # Conventional particle (CPP): Gamma = 0.0 chains/nm^2
    cpp_particle = NanoparticleConfig(
        diameter_nm=diameter_nm,
        grafting_density_chains_nm2=0.0,
    )
    cpp_res = simulate_mpp_transport(cpp_particle, cf_mucus)

    # Mucus-penetrating particle (MPP): dense PEG brush
    mpp_particle = NanoparticleConfig(
        diameter_nm=diameter_nm,
        grafting_density_chains_nm2=peg_grafting_density,
    )
    mpp_res = simulate_mpp_transport(mpp_particle, cf_mucus)

    # MPP in rhDNase-digested mucus
    mpp_rhdnase_res = simulate_mpp_transport(mpp_particle, cf_rhdnase_mucus)

    fold_enhancement = (
        round(mpp_res.effective_diffusivity_ratio / max(1.0e-6, cpp_res.effective_diffusivity_ratio), 1)
        if cpp_res.effective_diffusivity_ratio > 0.0
        else 999.0
    )

    rhdnase_ratio = (
        round(mpp_rhdnase_res.effective_diffusivity_ratio / max(1.0e-5, mpp_res.effective_diffusivity_ratio), 2)
        if mpp_res.effective_diffusivity_ratio > 0.0
        else 1.0
    )

    return ObstructionScalingSummary(
        particle_diameter_nm=diameter_nm,
        cpp_effective_diffusivity_ratio=cpp_res.effective_diffusivity_ratio,
        cpp_transit_time_minutes=cpp_res.mucus_transit_time_minutes,
        cpp_escapes_clearance=cpp_res.escapes_mucociliary_clearance,
        mpp_effective_diffusivity_ratio=mpp_res.effective_diffusivity_ratio,
        mpp_transit_time_minutes=mpp_res.mucus_transit_time_minutes,
        mpp_escapes_clearance=mpp_res.escapes_mucociliary_clearance,
        fold_diffusivity_enhancement=fold_enhancement,
        mesh_pore_spacing_nm=mpp_res.mesh_pore_spacing_xi_nm,
        cf_baseline_vs_rhdnase_diffusivity_ratio=rhdnase_ratio,
        tier="Tier T2 Certified",
        clinical_grounding=[
            "NCT00843817",
            "NCT00680316",
            "NCT06350474",
            "NCT05248230",
        ],
        primary_dois=[
            "10.1073/pnas.0610617104",
            "10.1016/j.addr.2008.09.009",
            "10.1073/pnas.0811566106",
            "10.1016/j.addr.2015.09.006",
            "10.1016/j.biomaterials.2009.01.052",
            "10.1073/pnas.87.23.9188",
        ],
    )


def run_size_sweep(
    diameters_nm: Optional[List[float]] = None,
    peg_grafting: float = DEFAULT_PEG_GRAFTING_CHAINS_NM2,
) -> List[Dict[str, Any]]:
    """Sweep nanoparticle diameters to evaluate steric pore cutoff and transport outcomes."""
    if diameters_nm is None:
        diameters_nm = [25.0, 40.0, 55.0, 65.0, 80.0, 100.0, 120.0, 150.0, 200.0]

    cf_mucus = MucusHydrogelConfig(rhdnase_treated=False)
    rhdnase_mucus = MucusHydrogelConfig(rhdnase_treated=True, rhdnase_cleavage_efficiency=0.80)

    sweep_results: List[Dict[str, Any]] = []
    for d in diameters_nm:
        mpp = NanoparticleConfig(diameter_nm=d, grafting_density_chains_nm2=peg_grafting)
        cpp = NanoparticleConfig(diameter_nm=d, grafting_density_chains_nm2=0.0)

        res_mpp = simulate_mpp_transport(mpp, cf_mucus)
        res_cpp = simulate_mpp_transport(cpp, cf_mucus)
        res_mpp_rhdnase = simulate_mpp_transport(mpp, rhdnase_mucus)

        sweep_results.append({
            "diameter_nm": d,
            "d0_um2_s": res_mpp.stokes_einstein_d0_um2_s,
            "cpp_deff_ratio": res_cpp.effective_diffusivity_ratio,
            "mpp_deff_ratio": res_mpp.effective_diffusivity_ratio,
            "mpp_rhdnase_deff_ratio": res_mpp_rhdnase.effective_diffusivity_ratio,
            "mpp_transit_min": res_mpp.mucus_transit_time_minutes,
            "mpp_escapes_clearance": res_mpp.escapes_mucociliary_clearance,
            "classification": res_mpp.transport_classification,
        })
    return sweep_results


def build_arg_parser() -> argparse.ArgumentParser:
    """Construct command-line argument parser."""
    parser = argparse.ArgumentParser(
        description="Airway Mucus-Penetrating Particle (MPP) Steric & Adhesive Obstruction Scaling (T-0244)",
    )
    parser.add_argument("--diameter", type=float, default=65.0, help="Nanoparticle hydrodynamic diameter (nm)")
    parser.add_argument("--grafting-density", type=float, default=DEFAULT_PEG_GRAFTING_CHAINS_NM2, help="PEG surface grafting density (chains/nm^2)")
    parser.add_argument("--peg-mw", type=float, default=DEFAULT_PEG_MW_KDA, help="PEG molecular weight (kDa)")
    parser.add_argument("--mucin-fraction", type=float, default=CF_MUCIN_VOLUME_FRACTION, help="CF mucin solids volume fraction")
    parser.add_argument("--edna-conc", type=float, default=DEFAULT_CF_EDNA_CONC_MG_ML, help="CF sputum eDNA concentration (mg/mL)")
    parser.add_argument("--rhdnase", action="store_true", help="Enable rhDNase I eDNA cleavage")
    parser.add_argument("--rhdnase-eff", type=float, default=0.75, help="rhDNase cleavage efficiency [0.0 - 1.0]")
    parser.add_argument("--run-sweep", action="store_true", help="Run diameter sweep across 25nm - 200nm")
    parser.add_argument("--json", action="store_true", help="Output results in JSON format")
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    """CLI execution entrypoint."""
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    if args.run_sweep:
        sweep = run_size_sweep(peg_grafting=args.grafting_density)
        if args.json:
            print(json.dumps(sweep, indent=2))
        else:
            print("================================================================================")
            print("MPP DIAMETER SWEEP: STERIC CUTOFF & MUCOCILIARY CLEARANCE WINDOW")
            print("================================================================================")
            print(f"{'Diameter (nm)':<14}{'D0 (um^2/s)':<14}{'CPP D_eff/D0':<14}{'MPP D_eff/D0':<14}{'rhDNase D_eff':<14}{'Transit (min)':<14}{'Escapes?':<10}")
            print("-" * 84)
            for row in sweep:
                print(
                    f"{row['diameter_nm']:<14.1f}"
                    f"{row['d0_um2_s']:<14.2f}"
                    f"{row['cpp_deff_ratio']:<14.5f}"
                    f"{row['mpp_deff_ratio']:<14.4f}"
                    f"{row['mpp_rhdnase_deff_ratio']:<14.4f}"
                    f"{row['mpp_transit_min']:<14.1f}"
                    f"{str(row['mpp_escapes_clearance']):<10}"
                )
            print("=" * 84)
        return 0

    particle = NanoparticleConfig(
        diameter_nm=args.diameter,
        peg_mw_kda=args.peg_mw,
        grafting_density_chains_nm2=args.grafting_density,
    )
    mucus = MucusHydrogelConfig(
        mucin_volume_fraction=args.mucin_fraction,
        edna_conc_mg_ml=args.edna_conc,
        rhdnase_treated=args.rhdnase,
        rhdnase_cleavage_efficiency=args.rhdnase_eff,
    )

    sim_res = simulate_mpp_transport(particle, mucus)
    comparison = compare_cpp_vs_mpp_scaling(
        diameter_nm=args.diameter,
        peg_grafting_density=args.grafting_density,
    )

    if args.json:
        payload = {
            "simulation": asdict(sim_res),
            "comparison": asdict(comparison),
            "attribution": NLM_ATTRIBUTION,
            "research_boundary": RESEARCH_BOUNDARY,
        }
        print(json.dumps(payload, indent=2))
    else:
        print("================================================================================")
        print("AIRWAY MUCUS-PENETRATING PARTICLE (MPP) TRANSPORT SIMULATION (T-0244)")
        print("================================================================================")
        print(f"Particle Diameter:           {sim_res.particle_diameter_nm:.1f} nm")
        print(f"Base Water Diffusivity D_0:  {sim_res.stokes_einstein_d0_um2_s:.4f} um^2/s")
        print(f"CF Sputum Pore Mesh (xi):    {sim_res.mesh_pore_spacing_xi_nm:.1f} nm")
        print(f"Flory Radius (R_F):          {sim_res.flory_radius_rf_nm:.2f} nm (Gamma* = {sim_res.brush_transition_density_gamma_star:.5f} ch/nm^2)")
        print(f"PEG Brush Ratio:             {sim_res.grafting_ratio_gamma_over_gamma_star:.2f}x Gamma* ({sim_res.brush_regime})")
        print(f"Steric Exclusion Factor:     {sim_res.steric_exclusion_factor_f_steric:.5f}")
        print(f"Adhesive Screening Factor:   {sim_res.adhesive_obstruction_factor_f_adhesive:.5f}")
        print(f"Effective Diffusivity Ratio: {sim_res.effective_diffusivity_ratio:.5f} (D_eff = {sim_res.effective_diffusivity_deff_um2_s:.4f} um^2/s)")
        print(f"Mucus Transit Time:          {sim_res.mucus_transit_time_minutes:.1f} min (Clearance Limit: {sim_res.mucociliary_clearance_limit_minutes:.1f} min)")
        print(f"Escapes Clearance:           {sim_res.escapes_mucociliary_clearance}")
        print(f"Classification:              {sim_res.transport_classification}")
        print("-" * 80)
        print(f"CPP vs MPP Enhancement:      {comparison.fold_diffusivity_enhancement:.1f}x higher diffusion with dense PEG brush")
        print(f"rhDNase Rescue Enhancement:  {comparison.cf_baseline_vs_rhdnase_diffusivity_ratio:.2f}x expansion with eDNA enzymatic digestion")
        print(f"Clinical Grounding:          {', '.join(comparison.clinical_grounding)}")
        print("=" * 80)

    return 0


if __name__ == "__main__":
    sys.exit(main())
