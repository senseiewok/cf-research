#!/usr/bin/env python3
"""Aerosol Droplet Aerodynamic Size Distribution & Mucociliary Impaction Mechanics for Inhaled Gene Therapies (T-0408).

Quantitative biophysical model of aerosol fluid dynamics, droplet aerodynamic size distributions,
and tracheobronchial multi-compartment airway deposition for inhaled cystic fibrosis gene therapies:
  1. Multi-Compartment Airway Deposition Fractions:
     Implements Weibel symmetric lung morphology (generations 0 to 16) modeling inertial impaction
     in upper conducting airways (G0-G6), gravitational sedimentation in small conducting
     bronchioles (G7-G12), and Brownian diffusion in terminal bronchioles (G13-G16).
  2. Log-Normal Droplet Aerodynamic Sizing:
     Analyzes Mass Median Aerodynamic Diameter (MMAD in [1.0, 5.0] um) and Geometric Standard
     Deviation (GSD in [1.5, 2.5]) to maximize target bronchial deposition while minimizing
     alveolar dissipation and oropharyngeal impaction loss.
  3. Inspiratory Flow & Airway Narrowing Impact:
     Simulates deposition shifts under CF pathological airway narrowing (luminal plugging,
     bronchiectasis, mucosal thickening) across variable inspiratory flow rates (Q_insp in [15, 60] L/min).
  4. Targeted Delivery Fraction Metric:
     Derives the Bronchial Delivery Index (BDI = Fraction_bronchial / Fraction_total_deposited)
     benchmarking delivery platforms (vibrating mesh nebulizers vs jet nebulizers vs DPI).

Primary Evidence Grounding:
  - Discovery Brief: cf-research/landscape/mission-grounding-2026-10-24.md (Curative Discovery Council Track A Wave 18).
  - Primary Literature:
    - Kim J et al. (2021) Nat Biomed Eng, DOI: 10.1038/s41551-021-00724-5
    - Front Bioeng Biotechnol (2026), DOI: 10.3389/fbioe.2026.1855907
    - Heyder J (2004) Proc Am Thorac Soc, DOI: 10.1513/pats.200409-046ta
  - Clinical Registry Grounding:
    - NCT05668741 (Phase 1/2 aerosolized mRNA LNP VX-522 in CF)
    - NCT05248230 (Phase 1/2 aerosolized AAV vector 4D-710 in CF)

Epistemic Status:
  - Answer Tier T2 (Open-access primary literature, ClinicalTrials.gov API v2).
  - Preclinical computational model; research, not medical advice.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

SCRIPT_PATH = Path(__file__).resolve()

NLM_ATTRIBUTION = (
    "This product uses publicly available data from the U.S. National Library of Medicine (NLM), "
    "National Institutes of Health, Department of Health and Human Services; NLM is not responsible "
    "for the product and does not endorse or recommend this or any other product. "
    "ClinicalTrials.gov does not evaluate safety or efficacy; listing on the registry does not denote government endorsement."
)

RESEARCH_BOUNDARY = (
    "Preclinical in-silico aerosol droplet aerodynamic size distribution & mucociliary impaction model; "
    "research, not medical advice. No clinical efficacy or therapeutic safety claim is made."
)

# Reference Grounding Identifiers
PRIMARY_DOI_KIM = "10.1038/s41551-021-00724-5"
PRIMARY_DOI_BIOENG = "10.3389/fbioe.2026.1855907"
PRIMARY_DOI_HEYDER = "10.1513/pats.200409-046ta"
CLINICAL_TRIAL_VX522 = "NCT05668741"
CLINICAL_TRIAL_4D710 = "NCT05248230"

# Physical and Physiological Constants
BOLTZMANN_CONSTANT = 1.380649e-23              # J / K
DEFAULT_TEMPERATURE_K = 310.15                 # 37 C (Physiological body temperature)
AIR_DYNAMIC_VISCOSITY_PA_S = 1.87e-5           # Dynamic viscosity of air at 37 C (Pa * s)
AIR_DENSITY_KG_M3 = 1.14                       # Density of humidified air at 37 C (kg / m^3)
UNIT_DROPLET_DENSITY_KG_M3 = 1000.0            # Standard unit droplet density (1 g / cm^3 = 1000 kg / m^3)
GRAVITATIONAL_ACCEL_M_S2 = 9.80665             # Gravitational acceleration (m / s^2)
AIR_MEAN_FREE_PATH_M = 6.65e-8                 # Mean free path of air molecules at 37 C (m)

# Standard Weibel Morphological Model A parameters for Generations 0 to 16
# Format: (generation, number_of_airways, diameter_mm, length_mm, branching_angle_deg)
WEIBEL_CONDUCTING_AIRWAYS: Tuple[Tuple[int, int, float, float, float], ...] = (
    (0, 1, 18.0, 120.0, 0.0),       # Trachea (G0)
    (1, 2, 12.2, 47.6, 35.0),       # Main bronchi (G1)
    (2, 4, 8.3, 19.0, 35.0),        # Lobar bronchi (G2)
    (3, 8, 5.6, 7.6, 35.0),         # Segmental bronchi (G3)
    (4, 16, 4.5, 12.7, 35.0),       # Subsegmental bronchi (G4)
    (5, 32, 3.5, 10.7, 35.0),       # Small bronchi (G5)
    (6, 64, 2.8, 9.0, 35.0),        # Small bronchi (G6)
    (7, 128, 2.3, 7.6, 35.0),       # Conducting bronchioles (G7)
    (8, 256, 1.86, 6.4, 35.0),      # Conducting bronchioles (G8)
    (9, 512, 1.54, 5.4, 35.0),      # Conducting bronchioles (G9)
    (10, 1024, 1.30, 4.6, 35.0),    # Conducting bronchioles (G10)
    (11, 2048, 1.09, 3.9, 35.0),    # Conducting bronchioles (G11)
    (12, 4096, 0.95, 3.3, 35.0),    # Conducting bronchioles (G12)
    (13, 8192, 0.82, 2.7, 35.0),    # Terminal bronchioles (G13)
    (14, 16384, 0.74, 2.3, 35.0),   # Terminal bronchioles (G14)
    (15, 32768, 0.66, 2.0, 35.0),   # Terminal bronchioles (G15)
    (16, 65536, 0.60, 1.7, 35.0),   # Terminal bronchioles (G16)
)


@dataclass(frozen=True)
class AirwayGeneration:
    """Biophysical characteristics of an airway generation."""

    generation: int
    num_airways: int
    diameter_m: float
    length_m: float
    branching_angle_rad: float
    cross_sectional_area_m2: float
    volume_m3: float


@dataclass(frozen=True)
class DropletDistribution:
    """Log-normal droplet aerodynamic size distribution parameters."""

    name: str
    mmad_um: float
    gsd: float
    description: str


@dataclass(frozen=True)
class GenerationDeposition:
    """Fractional deposition and mechanics for a specific airway generation."""

    generation: int
    airway_velocity_m_s: float
    residence_time_s: float
    stokes_number: float
    p_impaction: float
    p_sedimentation: float
    p_diffusion: float
    p_combined: float
    fraction_deposited: float
    fraction_entering: float


@dataclass(frozen=True)
class DepositionProfile:
    """Comprehensive regional lung deposition profile."""

    distribution: DropletDistribution
    inspiratory_flow_l_min: float
    cf_narrowing_factor: float
    extrathoracic_fraction: float
    upper_conducting_fraction: float      # G0-G6
    small_conducting_fraction: float      # G7-G12
    terminal_bronchiolar_fraction: float  # G13-G16
    tracheobronchial_total_fraction: float # G0-G16
    alveolar_fraction: float              # G17-G23 (peripheral loss)
    exhaled_fraction: float               # Non-deposited
    bronchial_delivery_index: float       # BDI = tracheobronchial / total deposited
    generation_breakdown: List[GenerationDeposition]


def cunningham_slip_correction(diameter_m: float) -> float:
    """Calculates the Cunningham slip correction factor Cc for a submicron/micron droplet.

    Cc = 1 + (2 * lambda / d) * [1.257 + 0.400 * exp(-1.10 * d / (2 * lambda))]
    """
    if diameter_m <= 0.0:
        return 1.0
    knudsen_ratio = (2.0 * AIR_MEAN_FREE_PATH_M) / diameter_m
    exponent = -0.55 * diameter_m / AIR_MEAN_FREE_PATH_M
    return 1.0 + knudsen_ratio * (1.257 + 0.400 * math.exp(exponent))


def droplet_relaxation_time(diameter_m: float, density_kg_m3: float = UNIT_DROPLET_DENSITY_KG_M3) -> float:
    """Calculates the aerodynamic particle relaxation time tau.

    tau = (rho_p * d^2 * Cc) / (18 * mu_air)
    """
    cc = cunningham_slip_correction(diameter_m)
    return (density_kg_m3 * (diameter_m ** 2) * cc) / (18.0 * AIR_DYNAMIC_VISCOSITY_PA_S)


def droplet_settling_velocity(diameter_m: float, density_kg_m3: float = UNIT_DROPLET_DENSITY_KG_M3) -> float:
    """Calculates gravitational settling velocity v_sed in still air.

    v_sed = tau * g
    """
    tau = droplet_relaxation_time(diameter_m, density_kg_m3)
    return tau * GRAVITATIONAL_ACCEL_M_S2


def droplet_brownian_diffusivity(diameter_m: float, temp_k: float = DEFAULT_TEMPERATURE_K) -> float:
    """Calculates Brownian diffusion coefficient D for a spherical droplet in air.

    D = (k_B * T * Cc) / (3 * pi * mu_air * d)
    """
    cc = cunningham_slip_correction(diameter_m)
    return (BOLTZMANN_CONSTANT * temp_k * cc) / (3.0 * math.pi * AIR_DYNAMIC_VISCOSITY_PA_S * diameter_m)


def build_weibel_tree(cf_narrowing_factor: float = 1.0) -> List[AirwayGeneration]:
    """Instantiates the Weibel Model A conducting airway morphology with optional CF narrowing."""
    narrowing = max(0.2, min(1.0, cf_narrowing_factor))
    tree: List[AirwayGeneration] = []

    for gen, n_airways, diam_mm, len_mm, angle_deg in WEIBEL_CONDUCTING_AIRWAYS:
        effective_diam_m = (diam_mm * 1e-3) * narrowing
        length_m = len_mm * 1e-3
        angle_rad = math.radians(angle_deg)
        single_area_m2 = math.pi * (effective_diam_m / 2.0) ** 2
        total_area_m2 = single_area_m2 * n_airways
        volume_m3 = total_area_m2 * length_m

        tree.append(
            AirwayGeneration(
                generation=gen,
                num_airways=n_airways,
                diameter_m=effective_diam_m,
                length_m=length_m,
                branching_angle_rad=angle_rad,
                cross_sectional_area_m2=total_area_m2,
                volume_m3=volume_m3,
            )
        )
    return tree


def extrathoracic_deposition_efficiency(
    diameter_ae_um: float,
    flow_l_min: float,
) -> float:
    """Calculates oropharyngeal / extrathoracic impaction efficiency based on Heyder empirical data.

    P_extra = 1 - 1 / (1 + 0.00022 * (d_ae^2 * Q)^1.5)
    """
    impaction_param = (diameter_ae_um ** 2) * flow_l_min
    if impaction_param <= 0.0:
        return 0.0
    denom = 1.0 + 0.00022 * (impaction_param ** 1.5)
    p_dep = 1.0 - (1.0 / denom)
    return max(0.0, min(0.999, p_dep))


def single_size_generation_mechanics(
    generation: AirwayGeneration,
    diameter_m: float,
    flow_m3_s: float,
) -> Tuple[float, float, float, float, float, float]:
    """Computes impaction, sedimentation, and diffusion probabilities for a single droplet size in a generation.

    Returns:
      (velocity, residence_time, stokes_num, p_imp, p_sed, p_diff)
    """
    v_air = flow_m3_s / generation.cross_sectional_area_m2
    t_res = generation.length_m / max(v_air, 1e-6)

    # 1. Inertial Impaction (Stokes number based)
    tau = droplet_relaxation_time(diameter_m)
    stokes_num = (tau * v_air) / generation.diameter_m
    sin_theta = math.sin(generation.branching_angle_rad)
    if sin_theta > 0.0:
        # Cai & Yu / Landahl empirical impaction parameter
        p_imp = 1.0 - math.exp(-2.0 * stokes_num * sin_theta)
    else:
        p_imp = 0.0
    p_imp = max(0.0, min(1.0, p_imp))

    # 2. Gravitational Sedimentation
    v_sed = droplet_settling_velocity(diameter_m)
    sed_param = (4.0 * v_sed * t_res) / (math.pi * generation.diameter_m)
    p_sed = 1.0 - math.exp(-sed_param)
    p_sed = max(0.0, min(1.0, p_sed))

    # 3. Brownian Diffusion (Gormley & Kennedy cylindrical tube deposition)
    diffusivity = droplet_brownian_diffusivity(diameter_m)
    q_single_tube = flow_m3_s / generation.num_airways
    delta = (math.pi * diffusivity * generation.length_m) / (4.0 * max(q_single_tube, 1e-12))
    if delta < 0.02:
        p_diff = 2.56 * (delta ** (2.0 / 3.0)) - 1.2 * delta - 0.177 * (delta ** (4.0 / 3.0))
    else:
        p_diff = 1.0 - 0.819 * math.exp(-3.657 * delta) - 0.097 * math.exp(-22.3 * delta)
    p_diff = max(0.0, min(1.0, p_diff))

    return v_air, t_res, stokes_num, p_imp, p_sed, p_diff


def log_normal_size_bins(
    mmad_um: float,
    gsd: float,
    num_bins: int = 40,
    d_min_um: float = 0.2,
    d_max_um: float = 12.0,
) -> List[Tuple[float, float]]:
    """Generates discretized droplet diameter bins and corresponding mass fractions for a log-normal distribution.

    Returns list of (diameter_um, mass_fraction).
    """
    ln_gsd = math.log(max(1.01, gsd))
    ln_mmad = math.log(max(0.1, mmad_um))

    ln_min = math.log(d_min_um)
    ln_max = math.log(d_max_um)
    step = (ln_max - ln_min) / num_bins

    bins: List[Tuple[float, float]] = []
    total_weight = 0.0

    for i in range(num_bins):
        ln_center = ln_min + (i + 0.5) * step
        d_center = math.exp(ln_center)
        # Log-normal probability density
        diff = ln_center - ln_mmad
        prob_density = (1.0 / (math.sqrt(2.0 * math.pi) * ln_gsd)) * math.exp(-(diff ** 2) / (2.0 * (ln_gsd ** 2)))
        weight = prob_density * step
        bins.append((d_center, weight))
        total_weight += weight

    # Normalize weights so sum equals 1.0
    normalized_bins = [(d, w / total_weight) for d, w in bins]
    return normalized_bins


def simulate_deposition(
    distribution: DropletDistribution,
    inspiratory_flow_l_min: float = 30.0,
    cf_narrowing_factor: float = 1.0,
) -> DepositionProfile:
    """Simulates full respiratory airway deposition across conducting generations G0-G16."""
    flow_l_min = max(5.0, min(90.0, inspiratory_flow_l_min))
    flow_m3_s = (flow_l_min * 1e-3) / 60.0
    tree = build_weibel_tree(cf_narrowing_factor)
    size_bins = log_normal_size_bins(distribution.mmad_um, distribution.gsd)

    # Accumulators across droplet size bins
    total_extrathoracic = 0.0
    gen_deposited_mass = [0.0] * len(tree)
    gen_entering_mass = [0.0] * len(tree)
    gen_avg_vel = [0.0] * len(tree)
    gen_avg_tres = [0.0] * len(tree)
    gen_avg_stk = [0.0] * len(tree)
    gen_avg_pimp = [0.0] * len(tree)
    gen_avg_psed = [0.0] * len(tree)
    gen_avg_pdiff = [0.0] * len(tree)
    gen_avg_pcomb = [0.0] * len(tree)
    total_alveolar_loss = 0.0

    for d_um, mass_frac in size_bins:
        d_m = d_um * 1e-6

        # Extrathoracic deposition
        p_extra = extrathoracic_deposition_efficiency(d_um, flow_l_min)
        dep_extra = mass_frac * p_extra
        total_extrathoracic += dep_extra

        remaining_mass = mass_frac - dep_extra

        # Traverse conducting airways G0 to G16
        for idx, gen in enumerate(tree):
            gen_entering_mass[idx] += remaining_mass
            v_air, t_res, stk, p_imp, p_sed, p_diff = single_size_generation_mechanics(
                gen, d_m, flow_m3_s
            )
            # Independent probabilities combined
            p_comb = 1.0 - (1.0 - p_imp) * (1.0 - p_sed) * (1.0 - p_diff)
            p_comb = max(0.0, min(1.0, p_comb))

            dep_gen = remaining_mass * p_comb
            gen_deposited_mass[idx] += dep_gen
            remaining_mass -= dep_gen

            # Weighted averages
            gen_avg_vel[idx] += v_air * mass_frac
            gen_avg_tres[idx] += t_res * mass_frac
            gen_avg_stk[idx] += stk * mass_frac
            gen_avg_pimp[idx] += p_imp * mass_frac
            gen_avg_psed[idx] += p_sed * mass_frac
            gen_avg_pdiff[idx] += p_diff * mass_frac
            gen_avg_pcomb[idx] += p_comb * mass_frac

        # Penetration into acinar/alveolar zone (generations 17-23)
        # Deep lung deposition parameterized via Heyder retention kinetics
        alveolar_retention_prob = 1.0 - math.exp(-0.45 * (d_um ** 0.8)) if d_um >= 0.5 else 0.2
        dep_alveolar = remaining_mass * alveolar_retention_prob
        total_alveolar_loss += dep_alveolar

    upper_conducting = sum(gen_deposited_mass[0:7])   # G0-G6
    small_conducting = sum(gen_deposited_mass[7:13])  # G7-G12
    terminal_conducting = sum(gen_deposited_mass[13:17]) # G13-G16
    tracheobronchial_total = upper_conducting + small_conducting + terminal_conducting

    total_deposited = total_extrathoracic + tracheobronchial_total + total_alveolar_loss
    exhaled = max(0.0, 1.0 - total_deposited)

    # Bronchial Delivery Index (BDI) = target conducting airways / total deposited
    bdi = (tracheobronchial_total / total_deposited) if total_deposited > 1e-6 else 0.0

    breakdown: List[GenerationDeposition] = []
    for idx, gen in enumerate(tree):
        breakdown.append(
            GenerationDeposition(
                generation=gen.generation,
                airway_velocity_m_s=round(gen_avg_vel[idx], 4),
                residence_time_s=round(gen_avg_tres[idx], 5),
                stokes_number=round(gen_avg_stk[idx], 6),
                p_impaction=round(gen_avg_pimp[idx], 4),
                p_sedimentation=round(gen_avg_psed[idx], 4),
                p_diffusion=round(gen_avg_pdiff[idx], 4),
                p_combined=round(gen_avg_pcomb[idx], 4),
                fraction_deposited=round(gen_deposited_mass[idx], 5),
                fraction_entering=round(gen_entering_mass[idx], 5),
            )
        )

    return DepositionProfile(
        distribution=distribution,
        inspiratory_flow_l_min=flow_l_min,
        cf_narrowing_factor=cf_narrowing_factor,
        extrathoracic_fraction=round(total_extrathoracic, 4),
        upper_conducting_fraction=round(upper_conducting, 4),
        small_conducting_fraction=round(small_conducting, 4),
        terminal_bronchiolar_fraction=round(terminal_conducting, 4),
        tracheobronchial_total_fraction=round(tracheobronchial_total, 4),
        alveolar_fraction=round(total_alveolar_loss, 4),
        exhaled_fraction=round(exhaled, 4),
        bronchial_delivery_index=round(bdi, 4),
        generation_breakdown=breakdown,
    )


# Standard Preset Delivery Platforms
PRESET_PLATFORMS: Dict[str, DropletDistribution] = {
    "vx522_mesh_lnp": DropletDistribution(
        name="VX-522 Inhaled mRNA LNP Mesh Aerosol",
        mmad_um=3.2,
        gsd=1.65,
        description="Vibrating mesh aerosolized lipid nanoparticle formulation for CFTR mRNA (NCT05668741).",
    ),
    "4d710_mesh_aav": DropletDistribution(
        name="4D-710 Inhaled AAV Vector Mesh Aerosol",
        mmad_um=2.8,
        gsd=1.75,
        description="Aerosolized targeted AAV capsid delivery for CFTR gene addition (NCT05248230).",
    ),
    "standard_jet_nebulizer": DropletDistribution(
        name="Standard Clinical Jet Nebulizer",
        mmad_um=4.6,
        gsd=2.2,
        description="High-velocity polydisperse pneumatic jet aerosolization with high oropharyngeal loss.",
    ),
    "optimized_monodisperse_mesh": DropletDistribution(
        name="Monodisperse Target Bronchial Aerosol",
        mmad_um=2.5,
        gsd=1.35,
        description="Narrow-distribution vibrating mesh tailored for small bronchiole conductance rescue.",
    ),
}


def evaluate_cf_obstruction_sweep(
    distribution: DropletDistribution,
    flow_l_min: float = 30.0,
    narrowing_steps: Optional[List[float]] = None,
) -> List[Dict[str, Any]]:
    """Simulates the impact of progressive CF airway narrowing (1.0 healthy down to 0.4 severe obstruction)."""
    if narrowing_steps is None:
        narrowing_steps = [1.0, 0.85, 0.70, 0.55, 0.40]

    sweep_results: List[Dict[str, Any]] = []
    for factor in narrowing_steps:
        profile = simulate_deposition(
            distribution=distribution,
            inspiratory_flow_l_min=flow_l_min,
            cf_narrowing_factor=factor,
        )
        sweep_results.append({
            "cf_narrowing_factor": factor,
            "luminal_occlusion_pct": round((1.0 - (factor ** 2)) * 100.0, 1),
            "extrathoracic_pct": round(profile.extrathoracic_fraction * 100.0, 2),
            "upper_conducting_pct": round(profile.upper_conducting_fraction * 100.0, 2),
            "small_conducting_pct": round(profile.small_conducting_fraction * 100.0, 2),
            "terminal_bronchiolar_pct": round(profile.terminal_bronchiolar_fraction * 100.0, 2),
            "tracheobronchial_total_pct": round(profile.tracheobronchial_total_fraction * 100.0, 2),
            "alveolar_pct": round(profile.alveolar_fraction * 100.0, 2),
            "bdi": profile.bronchial_delivery_index,
        })
    return sweep_results


def benchmark_delivery_platforms(
    flow_l_min: float = 30.0,
    cf_narrowing: float = 0.70,
) -> Dict[str, Any]:
    """Benchmarks all preset delivery platforms under standardized CF pathological conditions."""
    results: Dict[str, Any] = {}
    for key, platform in PRESET_PLATFORMS.items():
        profile = simulate_deposition(
            distribution=platform,
            inspiratory_flow_l_min=flow_l_min,
            cf_narrowing_factor=cf_narrowing,
        )
        results[key] = {
            "name": platform.name,
            "mmad_um": platform.mmad_um,
            "gsd": platform.gsd,
            "extrathoracic_pct": round(profile.extrathoracic_fraction * 100.0, 2),
            "tracheobronchial_pct": round(profile.tracheobronchial_total_fraction * 100.0, 2),
            "small_bronchioles_pct": round(profile.small_conducting_fraction * 100.0, 2),
            "alveolar_pct": round(profile.alveolar_fraction * 100.0, 2),
            "exhaled_pct": round(profile.exhaled_fraction * 100.0, 2),
            "bdi": profile.bronchial_delivery_index,
        }
    return results


def main() -> None:
    """Command-line interface for the aerosol droplet deposition mechanics model."""
    parser = argparse.ArgumentParser(
        description="Aerosol Droplet Aerodynamic Size Distribution & Mucociliary Impaction Mechanics (T-0408)."
    )
    parser.add_argument("--platform", choices=list(PRESET_PLATFORMS.keys()), default="vx522_mesh_lnp")
    parser.add_argument("--mmad", type=float, default=None, help="Custom MMAD in micrometers (1.0 to 6.0)")
    parser.add_argument("--gsd", type=float, default=None, help="Custom GSD (1.2 to 2.5)")
    parser.add_argument("--flow", type=float, default=30.0, help="Inspiratory flow rate in L/min (15 to 60)")
    parser.add_argument("--cf-narrowing", type=float, default=1.0, help="Airway narrowing factor (0.4 to 1.0)")
    parser.add_argument("--benchmark", action="store_true", help="Run comparative platform benchmark")
    parser.add_argument("--sweep", action="store_true", help="Run CF obstruction parameter sweep")
    parser.add_argument("--json", action="store_true", help="Emit raw JSON output")

    args = parser.parse_args()

    platform = PRESET_PLATFORMS[args.platform]
    mmad = args.mmad if args.mmad is not None else platform.mmad_um
    gsd = args.gsd if args.gsd is not None else platform.gsd
    dist = DropletDistribution(
        name=f"Custom ({platform.name})" if args.mmad or args.gsd else platform.name,
        mmad_um=mmad,
        gsd=gsd,
        description=platform.description,
    )

    if args.benchmark:
        bench_data = benchmark_delivery_platforms(flow_l_min=args.flow, cf_narrowing=args.cf_narrowing)
        if args.json:
            print(json.dumps(bench_data, indent=2))
        else:
            print("=" * 78)
            print("DELIVERY PLATFORM REGIONAL DEPOSITION BENCHMARK (CF Obstruction Narrowing = {:.2f})".format(args.cf_narrowing))
            print("=" * 78)
            for k, v in bench_data.items():
                print(f"Platform: {v['name']} (MMAD={v['mmad_um']} um, GSD={v['gsd']})")
                print(f"  - Extrathoracic Loss:   {v['extrathoracic_pct']:>5.2f}%")
                print(f"  - Tracheobronchial:     {v['tracheobronchial_pct']:>5.2f}% (Small Bronchioles: {v['small_bronchioles_pct']}%)")
                print(f"  - Alveolar Dissipation: {v['alveolar_pct']:>5.2f}%")
                print(f"  - Exhaled Fraction:     {v['exhaled_pct']:>5.2f}%")
                print(f"  - Bronchial Delivery Index (BDI): {v['bdi']:.3f}\n")
        return

    if args.sweep:
        sweep_data = evaluate_cf_obstruction_sweep(dist, flow_l_min=args.flow)
        if args.json:
            print(json.dumps(sweep_data, indent=2))
        else:
            print("=" * 78)
            print(f"CF PROGRESSIVE AIRWAY OBSTRUCTION SWEEP — {dist.name}")
            print("=" * 78)
            for row in sweep_data:
                print(f"Narrowing Factor: {row['cf_narrowing_factor']:.2f} (Occlusion: {row['luminal_occlusion_pct']:>4.1f}%)")
                print(f"  Upper Conducting: {row['upper_conducting_pct']:>5.2f}% | Small Bronchioles: {row['small_conducting_pct']:>5.2f}% | BDI: {row['bdi']:.3f}")
        return

    profile = simulate_deposition(dist, inspiratory_flow_l_min=args.flow, cf_narrowing_factor=args.cf_narrowing)

    if args.json:
        print(json.dumps(asdict(profile), indent=2))
    else:
        print("=" * 78)
        print(f"AEROSOL DROPLET DEPOSITION PROFILE: {dist.name}")
        print(f"MMAD: {dist.mmad_um:.2f} um | GSD: {dist.gsd:.2f} | Flow: {args.flow:.1f} L/min | CF Narrowing: {args.cf_narrowing:.2f}")
        print("=" * 78)
        print(f"Extrathoracic Impaction Loss:        {profile.extrathoracic_fraction * 100:>6.2f}%")
        print(f"Upper Conducting Airways (G0-G6):    {profile.upper_conducting_fraction * 100:>6.2f}%")
        print(f"Small Conducting Bronchioles (G7-G12): {profile.small_conducting_fraction * 100:>6.2f}%")
        print(f"Terminal Bronchioles (G13-G16):      {profile.terminal_bronchiolar_fraction * 100:>6.2f}%")
        print(f"Tracheobronchial Total (Target):     {profile.tracheobronchial_total_fraction * 100:>6.2f}%")
        print(f"Alveolar Gas-Exchange Loss (G17-23): {profile.alveolar_fraction * 100:>6.2f}%")
        print(f"Exhaled Non-Deposited Fraction:      {profile.exhaled_fraction * 100:>6.2f}%")
        print(f"Bronchial Delivery Index (BDI):      {profile.bronchial_delivery_index:>6.3f}")
        print("=" * 78)
        print(NLM_ATTRIBUTION)


if __name__ == "__main__":
    main()
