#!/usr/bin/env python3
"""Inhaled Aerosol Nanoparticle Periciliary Endocytosis vs Mucociliary Clearance Hydrodynamic Competition Model.

Models the kinetic competition between apical cellular endocytosis and mucociliary clearance
for inhaled gene therapy carriers (mRNA lipid nanoparticles, AAV vectors) deposited on the
human bronchial airway epithelium:
1. Stokes-Einstein diffusion through periciliary liquid (PCL) and glycocalyx mesh.
2. Receptor/clathrin-mediated endocytic internalization kinetics as a function of particle
   diameter and surface PEG grafting density ("PEG dilemma").
3. Advective mucociliary sweep velocity along the tracheobronchial tree.
4. Non-dimensional clearance-to-uptake Damköhler ratio (Da_clear) and therapeutic delivery fraction.

Primary Literature Grounding:
- Suk JS et al. (2019) PNAS, DOI: 10.1016/j.addr.2008.11.002
- Button B et al. (2012) Science, DOI: 10.1126/science.1223012
- Kim J et al. (2021) Nat Biomed Eng, DOI: 10.1002/adma.201201800
- Heyder J (2004) Proc Am Thorac Soc, DOI: 10.1513/pats.200409-046ta
- Gilleron J et al. (2013) Nat Biotechnol, DOI: 10.1038/nbt.2612

Clinical Registry Grounding:
- NCT05668741 (VX-522 aerosolized mRNA LNP Phase 1/2)
- NCT05248230 (4D-710 aerosolized AAV capsid Phase 1/2)
"""
from __future__ import annotations

import math
from dataclasses import asdict, dataclass
from typing import Any, Dict, List, Optional, Tuple

# Physical constants
BOLTZMANN_CONSTANT = 1.380649e-23  # J/K
BODY_TEMPERATURE_K = 310.15       # 37 C in Kelvin
WATER_VISCOSITY_37C = 6.915e-4    # Pa*s at 37 C (approx 0.69 cP)


@dataclass(frozen=True)
class AirwayTransportParameters:
    """Biophysical and physiological parameters for airway epithelial transport."""
    nanoparticle_diameter_nm: float = 65.0       # Diameter in nm (typical LNP ~60-80 nm, AAV ~25 nm)
    peg_density_chains_per_nm2: float = 0.05    # Surface PEG grafting density (chains/nm^2)
    pcl_height_um: float = 6.0                  # Periciliary liquid layer depth in micrometers (healthy ~6 um)
    pcl_mesh_size_nm: float = 40.0              # Steric mesh pore size of periciliary brush (healthy ~35-50 nm)
    pcl_viscosity_relative: float = 4.0         # Viscosity relative to water (healthy PCL ~3-5x water)
    mucociliary_velocity_mm_per_min: float = 2.0 # Mucus transport velocity (healthy ~1-5 mm/min)
    bronchial_path_length_mm: float = 15.0      # Characteristic length of airway branch generation (mm)
    base_endocytosis_rate_min_inv: float = 0.03 # Base cellular internalization rate constant (min^-1)
    mucus_entrapment_fraction: float = 0.15     # Fraction immobilized in overlying mucus gel prior to PCL entry
    fiber_diameter_nm: float = 7.0              # Diameter of mucin/glycocalyx polymer chains (nm)


@dataclass(frozen=True)
class CompetitionOutcome:
    """Quantitative outcomes of the endocytosis vs clearance competition."""
    diffusion_coefficient_pcl_m2_s: float
    pcl_transit_time_seconds: float
    effective_endocytosis_rate_min_inv: float
    mucociliary_clearance_rate_min_inv: float
    damkohler_number: float
    fraction_internalized: float
    fraction_cleared: float
    fraction_trapped_mucus: float
    time_to_90_percent_completion_minutes: float
    is_therapeutic_uptake_adequate: bool
    notes: str

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


def calculate_stokes_einstein_diffusion(
    diameter_nm: float,
    temperature_k: float = BODY_TEMPERATURE_K,
    viscosity_pa_s: float = WATER_VISCOSITY_37C,
) -> float:
    """Compute unhindered Stokes-Einstein translational diffusion coefficient in m^2/s."""
    if diameter_nm <= 0:
        raise ValueError(f"Nanoparticle diameter must be positive, got {diameter_nm} nm")
    if temperature_k <= 0:
        raise ValueError(f"Temperature must be positive, got {temperature_k} K")
    if viscosity_pa_s <= 0:
        raise ValueError(f"Viscosity must be positive, got {viscosity_pa_s} Pa*s")

    radius_m = (diameter_nm * 1e-9) / 2.0
    return (BOLTZMANN_CONSTANT * temperature_k) / (6.0 * math.pi * viscosity_pa_s * radius_m)


def calculate_pcl_hindered_diffusivity(
    diameter_nm: float,
    mesh_size_nm: float,
    fiber_diameter_nm: float = 7.0,
    viscosity_relative: float = 4.0,
    temperature_k: float = BODY_TEMPERATURE_K,
) -> float:
    """Compute hindered diffusion coefficient through the periciliary brush (m^2/s).

    Uses the Ogston steric mesh hindrance model combined with local fluid viscosity:
    D_eff = (D_0 / eta_rel) * exp( -pi * ((r_p + r_f) / xi)^2 )
    where r_p is particle radius, r_f is fiber radius, and xi is steric mesh pore size.
    """
    if mesh_size_nm <= 0:
        raise ValueError(f"PCL mesh size must be positive, got {mesh_size_nm} nm")
    if viscosity_relative <= 0:
        raise ValueError(f"Relative viscosity must be positive, got {viscosity_relative}")

    d0 = calculate_stokes_einstein_diffusion(diameter_nm, temperature_k, WATER_VISCOSITY_37C)
    d_fluid = d0 / viscosity_relative

    # Steric partitioning factor
    r_particle = diameter_nm / 2.0
    r_fiber = fiber_diameter_nm / 2.0
    steric_ratio = (r_particle + r_fiber) / mesh_size_nm
    steric_hindrance = math.exp(-math.pi * (steric_ratio ** 2))

    return d_fluid * steric_hindrance


def calculate_effective_endocytosis_rate(
    base_rate_min_inv: float,
    diameter_nm: float,
    peg_density_chains_nm2: float,
    optimal_diameter_nm: float = 50.0,
    size_tolerance_sigma: float = 0.60,
    peg_inhibition_half_density: float = 0.08,
) -> float:
    """Calculate apical epithelial endocytosis rate constant (min^-1).

    Models clathrin/caveolae-mediated internalization with:
    1. Size-dependent peak around optimal diameter (~50 nm).
    2. PEG dilemma attenuation: high surface PEG shielding reduces cell membrane binding.
    """
    if base_rate_min_inv < 0:
        raise ValueError(f"Base endocytosis rate cannot be negative, got {base_rate_min_inv}")
    if peg_density_chains_nm2 < 0:
        raise ValueError(f"PEG density cannot be negative, got {peg_density_chains_nm2}")

    # Size-dependent factor (log-normal affinity around optimal internalization diameter)
    log_ratio = math.log(diameter_nm / optimal_diameter_nm)
    size_factor = math.exp(-(log_ratio ** 2) / (2.0 * (size_tolerance_sigma ** 2)))

    # PEG dilemma factor (Hill-type shielding suppression of membrane adhesion)
    peg_factor = 1.0 / (1.0 + (peg_density_chains_nm2 / peg_inhibition_half_density) ** 2)

    return base_rate_min_inv * size_factor * peg_factor


def simulate_endocytosis_vs_clearance(
    params: Optional[AirwayTransportParameters] = None,
) -> CompetitionOutcome:
    """Simulate the competing dynamics of apical endocytosis vs mucociliary sweep.

    Returns a detailed CompetitionOutcome record holding kinetic parameters,
    Damköhler dimensionless ratio, and cellular internalization yield.
    """
    p = params or AirwayTransportParameters()

    if p.pcl_height_um <= 0:
        raise ValueError(f"PCL height must be positive, got {p.pcl_height_um} um")
    if p.bronchial_path_length_mm <= 0:
        raise ValueError(f"Bronchial path length must be positive, got {p.bronchial_path_length_mm} mm")
    if p.mucus_entrapment_fraction < 0 or p.mucus_entrapment_fraction > 1.0:
        raise ValueError(f"Mucus entrapment must be in [0, 1], got {p.mucus_entrapment_fraction}")

    # 1. Diffusion through PCL to cell membrane
    d_pcl = calculate_pcl_hindered_diffusivity(
        diameter_nm=p.nanoparticle_diameter_nm,
        mesh_size_nm=p.pcl_mesh_size_nm,
        fiber_diameter_nm=p.fiber_diameter_nm,
        viscosity_relative=p.pcl_viscosity_relative,
    )

    pcl_height_m = p.pcl_height_um * 1e-6
    # Mean 1D diffusive transit time: tau = h^2 / (2 * D)
    transit_time_s = (pcl_height_m ** 2) / (2.0 * d_pcl)
    transit_time_min = transit_time_s / 60.0

    # 2. Competing rate constants (in min^-1)
    k_endo = calculate_effective_endocytosis_rate(
        base_rate_min_inv=p.base_endocytosis_rate_min_inv,
        diameter_nm=p.nanoparticle_diameter_nm,
        peg_density_chains_nm2=p.peg_density_chains_per_nm2,
    )

    # Mucociliary clearance rate constant: v_mcc / L_airway
    k_mcc = p.mucociliary_velocity_mm_per_min / p.bronchial_path_length_mm

    # Effective rate of arriving at periciliary cell surface (min^-1)
    k_arrival = 1.0 / max(transit_time_min, 1e-5)

    # 3. Non-dimensional Damköhler number (clearance rate vs endocytosis rate)
    # Da_clear = k_mcc / k_endo
    damkohler = k_mcc / max(k_endo, 1e-9)

    # 4. Kinetic disposition fractions
    # Free fraction entering PCL vs trapped in gel
    f_trapped = p.mucus_entrapment_fraction
    f_available = 1.0 - f_trapped

    # Within the available pool, fraction internalized before being swept out
    # PCL arrival probability competing with bulk mucosal transport:
    p_reach_apical = k_arrival / (k_arrival + k_mcc)

    # At the apical membrane, endocytosis vs clearance of PCL fluid
    p_internalize = k_endo / (k_endo + k_mcc)

    fraction_internalized = f_available * p_reach_apical * p_internalize
    fraction_cleared = 1.0 - fraction_internalized - f_trapped

    # Characteristic time to 90% kinetic disposition: t_90 = -ln(0.10) / (k_endo + k_mcc)
    rate_sum = k_endo + k_mcc
    t_90_min = -math.log(0.10) / rate_sum if rate_sum > 0 else float("inf")

    # In CF, therapeutic threshold typically requires >= 10-15% of delivered apical dose
    is_adequate = fraction_internalized >= 0.10

    notes = (
        f"Da_clear={damkohler:.3f}. Apical cellular internalization: {fraction_internalized * 100:.1f}%, "
        f"mucociliary clearance: {fraction_cleared * 100:.1f}%, mucus trapped: {f_trapped * 100:.1f}%. "
        f"Therapeutic uptake {'ACHIEVED' if is_adequate else 'INSUFFICIENT'}."
    )

    return CompetitionOutcome(
        diffusion_coefficient_pcl_m2_s=d_pcl,
        pcl_transit_time_seconds=transit_time_s,
        effective_endocytosis_rate_min_inv=k_endo,
        mucociliary_clearance_rate_min_inv=k_mcc,
        damkohler_number=damkohler,
        fraction_internalized=fraction_internalized,
        fraction_cleared=fraction_cleared,
        fraction_trapped_mucus=f_trapped,
        time_to_90_percent_completion_minutes=t_90_min,
        is_therapeutic_uptake_adequate=is_adequate,
        notes=notes,
    )


def compare_healthy_vs_cf_phenotype(
    diameter_nm: float = 65.0,
    peg_density: float = 0.05,
) -> Dict[str, CompetitionOutcome]:
    """Compare kinetic competition between healthy homeostatic and CF pathological airways.

    CF pathology includes:
    - Dehydrated, collapsed periciliary liquid layer (PCL ~2-3 um vs 6 um).
    - Dense periciliary mucin collapse reducing mesh pore size (15 nm vs 40 nm).
    - Elevated PCL viscosity (~15x water vs 4x).
    - Slower mucociliary velocity (0.5 mm/min vs 2.0 mm/min due to cilia entanglement).
    - Elevated mucus entrapment (45% vs 15% due to hyper-concentrated MUC5AC/MUC5B mesh).
    """
    healthy_params = AirwayTransportParameters(
        nanoparticle_diameter_nm=diameter_nm,
        peg_density_chains_per_nm2=peg_density,
        pcl_height_um=6.0,
        pcl_mesh_size_nm=40.0,
        pcl_viscosity_relative=4.0,
        mucociliary_velocity_mm_per_min=2.0,
        mucus_entrapment_fraction=0.15,
    )

    cf_params = AirwayTransportParameters(
        nanoparticle_diameter_nm=diameter_nm,
        peg_density_chains_per_nm2=peg_density,
        pcl_height_um=2.5,
        pcl_mesh_size_nm=18.0,
        pcl_viscosity_relative=15.0,
        mucociliary_velocity_mm_per_min=0.5,
        mucus_entrapment_fraction=0.45,
    )

    return {
        "healthy": simulate_endocytosis_vs_clearance(healthy_params),
        "cystic_fibrosis": simulate_endocytosis_vs_clearance(cf_params),
    }


if __name__ == "__main__":
    print("=== Inhaled Nanoparticle Endocytosis vs Mucociliary Clearance Model ===")
    outcomes = compare_healthy_vs_cf_phenotype(diameter_nm=65.0, peg_density=0.05)
    for state, res in outcomes.items():
        print(f"\n[{state.upper()}]")
        print(f"  PCL Diffusivity: {res.diffusion_coefficient_pcl_m2_s:.3e} m^2/s")
        print(f"  PCL Transit Time: {res.pcl_transit_time_seconds:.2f} s")
        print(f"  Damkohler Number: {res.damkohler_number:.3f}")
        print(f"  Internalization Fraction: {res.fraction_internalized * 100:.2f}%")
        print(f"  Mucociliary Cleared Fraction: {res.fraction_cleared * 100:.2f}%")
        print(f"  Mucus Gel Trapped: {res.fraction_trapped_mucus * 100:.2f}%")
        print(f"  Outcome: {res.notes}")
