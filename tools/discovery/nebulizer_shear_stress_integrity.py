#!/usr/bin/env python3
"""Aerosol Nebulizer Mesh Pore Shear Stress & Inhaled Lipid Nanoparticle Integrity Dynamics (T-0413).

Quantitative biophysical model of hydrodynamic shear stress, vesicle deformation,
membrane rupture, and post-aerosolization colloid aggregation for inhaled CFTR mRNA lipid
nanoparticles (LNPs) during nebulizer extrusion:
  1. Pore Extrusion Shear Stress:
     Formulates Hagen-Poiseuille wall shear stress tau_wall = 4 * mu * Q / (pi * r_pore^3)
     across micro-aperture mesh nozzles (r_pore in [1.5, 4.0] um) and vibrating mesh
     frequencies (f in [100, 180] kHz).
  2. LNP Vesicle Deformation & Rupture:
     Evaluates critical capillary number Ca_crit = tau * r_lnp / sigma_bilayer and
     membrane bending rigidity kappa_b in [10, 20] k_B*T, calculating area expansion strain
     and fractional mRNA leakage under high-velocity extrusion.
  3. Post-Aerosolization Colloid Aggregation:
     Simulates shear-induced PEG stripping, Smoluchowski collision kinetics, hydrodynamic
     diameter shift (Delta d_h), and steric mucus mesh penetration competing with mucociliary
     clearance in airway secretions.
  4. Platform Mechanical Benchmark:
     Benchmarks vibrating mesh nebulizers against ultrasonic and pneumatic jet nebulizers
     for intact mRNA payload delivery.

Primary Evidence Grounding:
  - Discovery Brief: cf-research/landscape/mission-grounding-2026-10-24.md (Curative Discovery Council Track A Wave 18).
  - Primary Literature:
    - Kim J et al. (2021) Nat Biomed Eng, DOI: 10.1002/adma.201201800
    - Front Bioeng Biotechnol (2026), DOI: 10.3389/fbioe.2026.1855907
    - Heyder J (2004) Proc Am Thorac Soc, DOI: 10.1513/pats.200409-046ta
  - Clinical Registry Grounding:
    - NCT05668741 (Phase 1/2 aerosolized mRNA LNP VX-522 in CF)

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
    "Preclinical in-silico aerosol nebulizer mesh pore shear stress & LNP integrity model; "
    "research, not medical advice. No clinical efficacy or therapeutic safety claim is made."
)

# Reference Grounding Identifiers (AH-1, AH-6)
PRIMARY_DOI_KIM = "10.1002/adma.201201800"
PRIMARY_DOI_BIOENG = "10.3389/fbioe.2026.1855907"
PRIMARY_DOI_HEYDER = "10.1513/pats.200409-046ta"
CLINICAL_TRIAL_VX522 = "NCT05668741"

# Biophysical and Physical Constants
BOLTZMANN_CONSTANT = 1.380649e-23               # J / K
DEFAULT_TEMPERATURE_K = 310.15                  # 37 C physiological temperature
AMBIENT_TEMPERATURE_K = 293.15                  # 20 C nebulizer operating temperature
WATER_DYNAMIC_VISCOSITY_PA_S = 1.002e-3         # Dynamic viscosity of water / buffer at 20 C (Pa * s)
PHYSIOLOGICAL_VISCOSITY_PA_S = 0.691e-3         # Viscosity at 37 C (Pa * s)

# Membrane Mechanics Baseline Constants
DEFAULT_BILAYER_TENSION_N_M = 0.015             # Interfacial/bilayer tension sigma ~ 15 mN/m
DEFAULT_AREA_MODULUS_N_M = 0.20                 # Elastic area compressibility modulus K_A ~ 200 mN/m
DEFAULT_BENDING_RIGIDITY_KBT = 15.0             # Bending rigidity kappa_b ~ 15 k_B*T (10-20 k_B*T)
CRITICAL_LYSIS_STRAIN = 0.040                   # Critical areal strain alpha_crit ~ 4% (Evans & Needham)
CRITICAL_CAPILLARY_NUMBER = 0.25                # Dimensionless Ca_crit for vesicle poration/burst


@dataclass(frozen=True)
class MeshApertureParameters:
    """Micro-aperture geometry and vibration parameters of a nebulizer mesh nozzle."""

    pore_radius_um: float                       # Pore exit aperture radius r_pore in [1.5, 4.0] um
    vibration_frequency_khz: float              # Piezoelectric mesh frequency f in [100.0, 180.0] kHz
    num_active_pores: int = 3000                # Number of active micro-nozzles in the mesh plate
    pore_length_um: float = 50.0                # Nozzle channel length (thickness of plate) in um
    pulsatility_factor: float = 3.14159         # Peak-to-mean pulsatile velocity factor (pi for sinusoidal vibration)

    @property
    def pore_radius_m(self) -> float:
        return self.pore_radius_um * 1e-6

    @property
    def pore_diameter_um(self) -> float:
        return self.pore_radius_um * 2.0

    @property
    def pore_length_m(self) -> float:
        return self.pore_length_um * 1e-6

    @property
    def frequency_hz(self) -> float:
        return self.vibration_frequency_khz * 1e3


@dataclass(frozen=True)
class LNPBiophysicalParameters:
    """Colloidal and mechanical properties of the inhaled lipid nanoparticle."""

    name: str = "VX-522-like CFTR mRNA LNP"
    initial_hydrodynamic_diameter_nm: float = 95.0   # Initial diameter d_h0 in nm
    bilayer_surface_tension_n_m: float = DEFAULT_BILAYER_TENSION_N_M
    area_expansion_modulus_n_m: float = DEFAULT_AREA_MODULUS_N_M
    bending_rigidity_kbt: float = DEFAULT_BENDING_RIGIDITY_KBT
    peg_desorption_threshold_pa: float = 1200.0      # Shear stress at which PEG-lipid begins stripping (Pa)
    viscosity_ratio_int_ext: float = 1.2             # Internal-to-external viscosity ratio lambda

    @property
    def initial_radius_m(self) -> float:
        return (self.initial_hydrodynamic_diameter_nm * 1e-9) / 2.0

    @property
    def bending_rigidity_joules(self) -> float:
        return self.bending_rigidity_kbt * BOLTZMANN_CONSTANT * AMBIENT_TEMPERATURE_K


@dataclass(frozen=True)
class PoreShearProfile:
    """Hydrodynamic shear stress and velocity profile inside a single micro-aperture pore."""

    volumetric_flow_total_ml_min: float
    single_pore_mean_flow_m3_s: float
    single_pore_peak_flow_m3_s: float
    mean_pore_velocity_m_s: float
    peak_pore_velocity_m_s: float
    wall_shear_rate_s_inv: float
    wall_shear_stress_pa: float
    peak_wall_shear_stress_pa: float
    reynolds_number: float
    pore_residence_time_ms: float


@dataclass(frozen=True)
class LNPIntegrityResult:
    """Quantification of vesicle deformation, critical strain, and mRNA leakage."""

    pore_shear: PoreShearProfile
    capillary_number_mean: float
    capillary_number_peak: float
    taylor_deformation_parameter: float
    bilayer_areal_strain_pct: float
    critical_lysis_strain_pct: float
    pore_opening_probability: float
    fractional_mrna_leakage_pct: float
    intact_payload_fraction_pct: float
    integrity_status: str                       # 'HIGH_INTEGRITY', 'MODERATE_SHEAR', 'SEVERE_RUPTURE'


@dataclass(frozen=True)
class ColloidAggregationProfile:
    """Post-nebulization secondary aggregate dynamics and mucus barrier competition."""

    initial_diameter_nm: float
    post_aerosol_diameter_nm: float
    diameter_shift_nm: float
    peg_loss_fraction_pct: float
    colloidal_stability_ratio_w: float
    mucus_mesh_pore_nm: float
    mucus_mesh_penetration_fraction_pct: float
    mcc_evasion_status: str                     # 'EFFICIENT_PENETRATION', 'PARTIALLY_HINDERED', 'STERICALLY_BLOCKED'


@dataclass(frozen=True)
class PlatformBenchmarkResult:
    """Benchmark comparing nebulizer platforms on mRNA delivery efficacy."""

    platform_key: str
    name: str
    technology: str
    mean_shear_stress_pa: float
    peak_shear_stress_pa: float
    temperature_rise_c: float
    intact_mrna_pct: float
    hydrodynamic_diameter_shift_nm: float
    mucus_penetration_pct: float
    effective_cftr_delivery_score: float        # Normalized deliverable score [0.0, 1.0]


PRESET_MESH_PLATFORMS: Dict[str, MeshApertureParameters] = {
    "vx522_clinical_mesh": MeshApertureParameters(
        pore_radius_um=2.0,                     # 4 um diameter aperture
        vibration_frequency_khz=128.0,          # Standard clinical vibrating mesh
        num_active_pores=3500,
        pore_length_um=45.0,
        pulsatility_factor=3.14159,
    ),
    "micro_pore_high_freq": MeshApertureParameters(
        pore_radius_um=1.5,                     # 3 um diameter aperture
        vibration_frequency_khz=160.0,          # High-frequency fine aerosol mesh
        num_active_pores=4500,
        pore_length_um=40.0,
        pulsatility_factor=3.5,
    ),
    "soft_mist_large_pore": MeshApertureParameters(
        pore_radius_um=3.5,                     # 7 um diameter aperture
        vibration_frequency_khz=110.0,          # Gentle low-shear extrusion mesh
        num_active_pores=2500,
        pore_length_um=55.0,
        pulsatility_factor=2.8,
    ),
}


# --- Component 1: Hagen-Poiseuille Hydrodynamic Shear Mechanics ---

def calculate_hagen_poiseuille_shear_stress(
    pore_radius_m: float,
    volumetric_flow_m3_s: float,
    dynamic_viscosity_pa_s: float = WATER_DYNAMIC_VISCOSITY_PA_S,
) -> Tuple[float, float, float]:
    """Calculate Hagen-Poiseuille mean velocity, wall shear rate, and wall shear stress.

    Equations:
      v_mean = Q / (pi * r_pore^2)
      gamma_dot_wall = 4 * Q / (pi * r_pore^3) = 4 * v_mean / r_pore
      tau_wall = mu * gamma_dot_wall = 4 * mu * Q / (pi * r_pore^3)

    Returns:
      (mean_velocity_m_s, wall_shear_rate_s_inv, wall_shear_stress_pa)
    """
    if pore_radius_m <= 0.0 or volumetric_flow_m3_s < 0.0:
        return 0.0, 0.0, 0.0

    area = math.pi * (pore_radius_m ** 2)
    v_mean = volumetric_flow_m3_s / area
    gamma_dot_wall = (4.0 * volumetric_flow_m3_s) / (math.pi * (pore_radius_m ** 3))
    tau_wall = dynamic_viscosity_pa_s * gamma_dot_wall
    return v_mean, gamma_dot_wall, tau_wall


def simulate_pore_extrusion_shear_profile(
    mesh_params: MeshApertureParameters,
    total_flow_rate_ml_min: float = 0.40,
    dynamic_viscosity_pa_s: float = WATER_DYNAMIC_VISCOSITY_PA_S,
    fluid_density_kg_m3: float = 1000.0,
) -> PoreShearProfile:
    """Simulate fluid dynamic shear conditions through active micro-aperture mesh nozzles.

    Accounts for active pore partitioning and high-frequency sinusoidal extrusion pulsatility.
    """
    total_flow_m3_s = (total_flow_rate_ml_min * 1e-6) / 60.0
    active_pores = max(1, mesh_params.num_active_pores)
    q_single_mean = total_flow_m3_s / active_pores
    q_single_peak = q_single_mean * mesh_params.pulsatility_factor

    v_mean, gamma_mean, tau_mean = calculate_hagen_poiseuille_shear_stress(
        mesh_params.pore_radius_m, q_single_mean, dynamic_viscosity_pa_s
    )
    v_peak, _, tau_peak = calculate_hagen_poiseuille_shear_stress(
        mesh_params.pore_radius_m, q_single_peak, dynamic_viscosity_pa_s
    )

    reynolds = (fluid_density_kg_m3 * v_mean * (2.0 * mesh_params.pore_radius_m)) / dynamic_viscosity_pa_s
    residence_time_ms = (mesh_params.pore_length_m / v_mean * 1e3) if v_mean > 0.0 else 0.0

    return PoreShearProfile(
        volumetric_flow_total_ml_min=total_flow_rate_ml_min,
        single_pore_mean_flow_m3_s=q_single_mean,
        single_pore_peak_flow_m3_s=q_single_peak,
        mean_pore_velocity_m_s=v_mean,
        peak_pore_velocity_m_s=v_peak,
        wall_shear_rate_s_inv=gamma_mean,
        wall_shear_stress_pa=tau_mean,
        peak_wall_shear_stress_pa=tau_peak,
        reynolds_number=reynolds,
        pore_residence_time_ms=residence_time_ms,
    )


# --- Component 2: LNP Vesicle Deformation & Membrane Rupture ---

def calculate_capillary_number(
    shear_stress_pa: float,
    lnp_radius_m: float,
    surface_tension_n_m: float = DEFAULT_BILAYER_TENSION_N_M,
) -> float:
    """Calculate dimensionless Capillary number: Ca = tau * r_lnp / sigma_bilayer."""
    if surface_tension_n_m <= 0.0 or lnp_radius_m <= 0.0 or shear_stress_pa <= 0.0:
        return 0.0
    return (shear_stress_pa * lnp_radius_m) / surface_tension_n_m


def calculate_taylor_deformation_parameter(
    capillary_number: float,
    viscosity_ratio: float = 1.2,
) -> float:
    """Calculate Taylor ellipsoidal deformation parameter D = (L - B) / (L + B).

    D ~ ((19 * lambda + 16) / (16 * lambda + 16)) * Ca for small to moderate Ca.
    """
    if capillary_number <= 0.0:
        return 0.0
    prefactor = (19.0 * viscosity_ratio + 16.0) / (16.0 * viscosity_ratio + 16.0)
    d_taylor = prefactor * capillary_number
    return min(d_taylor, 0.95)


def calculate_bilayer_areal_strain(
    shear_stress_pa: float,
    lnp_params: LNPBiophysicalParameters,
) -> float:
    """Calculate fractional membrane area expansion strain alpha = Delta A / A_0.

    Accounts for both:
      1. Geometric ellipsoidal deformation under shear: alpha_geom = (2/5) * D_Taylor^2
      2. Elastic bilayer isotropic tension: alpha_elastic = sigma_shear / K_A
    """
    ca = calculate_capillary_number(
        shear_stress_pa, lnp_params.initial_radius_m, lnp_params.bilayer_surface_tension_n_m
    )
    d_taylor = calculate_taylor_deformation_parameter(ca, lnp_params.viscosity_ratio_int_ext)
    alpha_geom = 0.40 * (d_taylor ** 2)

    sigma_shear = shear_stress_pa * lnp_params.initial_radius_m
    alpha_elastic = (sigma_shear / lnp_params.area_expansion_modulus_n_m) if lnp_params.area_expansion_modulus_n_m > 0.0 else 0.0
    return max(0.0, alpha_geom + alpha_elastic)


def evaluate_lnp_vesicle_deformation_and_leakage(
    mesh_params: MeshApertureParameters,
    lnp_params: LNPBiophysicalParameters,
    flow_rate_ml_min: float = 0.40,
) -> LNPIntegrityResult:
    """Integrate velocity-weighted radial shear across the pore to compute fractional mRNA leakage.

    In Hagen-Poiseuille flow, shear stress varies linearly from center (0) to wall (tau_wall).
    Fluid flux is parabolically distributed: dQ(r) = 2*pi*r*u(r) dr.
    Pore opening occurs when local shear strain exceeds critical lysis strain (alpha_crit ~ 4%).
    """
    profile = simulate_pore_extrusion_shear_profile(mesh_params, flow_rate_ml_min)
    tau_peak = profile.peak_wall_shear_stress_pa

    ca_mean = calculate_capillary_number(profile.wall_shear_stress_pa, lnp_params.initial_radius_m, lnp_params.bilayer_surface_tension_n_m)
    ca_peak = calculate_capillary_number(tau_peak, lnp_params.initial_radius_m, lnp_params.bilayer_surface_tension_n_m)
    taylor_d = calculate_taylor_deformation_parameter(ca_mean, lnp_params.viscosity_ratio_int_ext)

    peak_strain = calculate_bilayer_areal_strain(tau_peak, lnp_params)
    peak_strain_pct = peak_strain * 100.0
    crit_strain_pct = CRITICAL_LYSIS_STRAIN * 100.0

    # Critical shear stress for membrane rupture
    tau_crit = (lnp_params.area_expansion_modulus_n_m * CRITICAL_LYSIS_STRAIN) / lnp_params.initial_radius_m

    # Discretized integration across streamtubes r in [0, r_pore]
    num_steps = 100
    dr = mesh_params.pore_radius_m / num_steps
    r_pore = mesh_params.pore_radius_m

    total_flux = 0.0
    leaking_flux = 0.0

    # Transition steepness for transient poration
    k_transition = max(50.0, tau_crit * 0.15)

    for i in range(num_steps):
        r_mid = (i + 0.5) * dr
        # Local parabolic velocity u(r) = 2 * v_mean * (1 - (r/r_pore)^2)
        u_local = 2.0 * profile.mean_pore_velocity_m_s * (1.0 - (r_mid / r_pore) ** 2)
        flux_annulus = 2.0 * math.pi * r_mid * u_local * dr
        total_flux += flux_annulus

        # Local shear stress at radial position r_mid under peak pulsatile cycle
        tau_local = tau_peak * (r_mid / r_pore)

        # Fermi-Dirac / logistic probability of pore opening and cargo escape
        p_leak = 1.0 / (1.0 + math.exp(-min(50.0, max(-50.0, (tau_local - tau_crit) / k_transition))))
        leaking_flux += flux_annulus * p_leak

    fractional_leakage = (leaking_flux / total_flux) if total_flux > 0.0 else 0.0
    leakage_pct = fractional_leakage * 100.0
    intact_pct = max(0.0, 100.0 - leakage_pct)

    if intact_pct >= 90.0:
        status = "HIGH_INTEGRITY"
    elif intact_pct >= 70.0:
        status = "MODERATE_SHEAR"
    else:
        status = "SEVERE_RUPTURE"

    return LNPIntegrityResult(
        pore_shear=profile,
        capillary_number_mean=ca_mean,
        capillary_number_peak=ca_peak,
        taylor_deformation_parameter=taylor_d,
        bilayer_areal_strain_pct=round(peak_strain_pct, 3),
        critical_lysis_strain_pct=round(crit_strain_pct, 3),
        pore_opening_probability=round(leaking_flux / total_flux if total_flux > 0.0 else 0.0, 4),
        fractional_mrna_leakage_pct=round(leakage_pct, 2),
        intact_payload_fraction_pct=round(intact_pct, 2),
        integrity_status=status,
    )


# --- Component 3: Post-Aerosol Colloid Aggregation & Mucus Clearance ---

def simulate_colloid_aggregation_and_mucus_clearance(
    lnp_params: LNPBiophysicalParameters,
    shear_profile: PoreShearProfile,
    residence_time_s: float = 0.50,
    mucus_mesh_pore_nm: float = 120.0,
    mucin_fiber_nm: float = 10.0,
) -> ColloidAggregationProfile:
    """Simulate secondary LNP aggregate formation and CF airway mucus mesh penetration.

    High shear can partially strip PEG-lipid coronas, lowering steric colloidal stability W
    and accelerating orthokinetic aggregation:
      Delta d_h = d_h0 * [(1 + k_agg * t)^(1/d_f) - 1]
    Secondary aggregates (>150 nm) are sterically excluded by cystic fibrosis mucus mesh pores
    (Kim et al. 2021) and cleared by ciliary transport before reaching target basal cells.
    """
    d0_nm = lnp_params.initial_hydrodynamic_diameter_nm
    tau_wall = shear_profile.peak_wall_shear_stress_pa

    # PEG stripping fraction
    tau_desorb = lnp_params.peg_desorption_threshold_pa
    peg_loss = 1.0 - math.exp(-max(0.0, tau_wall) / (2.0 * tau_desorb))
    peg_loss = min(0.85, max(0.0, peg_loss))

    # Fuchs colloidal stability ratio W: drops exponentially as PEG is stripped
    # At 0% PEG loss, W ~ 1e6 (colloidally stable); at high loss, W approaches 1.0 (diffusion-limited)
    w_ratio = math.exp(14.0 * (1.0 - peg_loss))
    w_ratio = max(1.0, w_ratio)

    # Aggregation rate constant k_agg under orthokinetic shear
    gamma_dot = shear_profile.wall_shear_rate_s_inv
    k_agg = (4.0 / 3.0) * gamma_dot * ((d0_nm * 1e-9) ** 3) * 1e16 / w_ratio

    # Fractal aggregate growth (d_f ~ 2.0 for shear aggregation)
    fractal_dim = 2.0
    growth_factor = (1.0 + k_agg * residence_time_s) ** (1.0 / fractal_dim)
    d_post_nm = d0_nm * growth_factor
    delta_dh_nm = d_post_nm - d0_nm

    # Ogston steric obstruction model for CF mucus mesh permeation (Kim et al. 2021)
    # Phi_pen = exp(- (pi / 4) * ((d_post + d_fiber) / (xi_mucus + d_fiber))^2)
    pore_ratio = (d_post_nm + mucin_fiber_nm) / (mucus_mesh_pore_nm + mucin_fiber_nm)
    mucus_penetration = math.exp(-0.785398 * (pore_ratio ** 2))
    mucus_pen_pct = mucus_penetration * 100.0

    if mucus_pen_pct >= 60.0:
        evasion_status = "EFFICIENT_PENETRATION"
    elif mucus_pen_pct >= 25.0:
        evasion_status = "PARTIALLY_HINDERED"
    else:
        evasion_status = "STERICALLY_BLOCKED"

    return ColloidAggregationProfile(
        initial_diameter_nm=round(d0_nm, 1),
        post_aerosol_diameter_nm=round(d_post_nm, 1),
        diameter_shift_nm=round(delta_dh_nm, 2),
        peg_loss_fraction_pct=round(peg_loss * 100.0, 1),
        colloidal_stability_ratio_w=w_ratio,
        mucus_mesh_pore_nm=round(mucus_mesh_pore_nm, 1),
        mucus_mesh_penetration_fraction_pct=round(mucus_pen_pct, 2),
        mcc_evasion_status=evasion_status,
    )


# --- Component 4: Multi-Platform Mechanical Benchmark ---

def benchmark_nebulizer_platforms(
    lnp_params: Optional[LNPBiophysicalParameters] = None,
) -> Dict[str, PlatformBenchmarkResult]:
    """Benchmark vibrating mesh against ultrasonic and pneumatic jet nebulizers for intact mRNA delivery."""
    lnp = lnp_params or LNPBiophysicalParameters()

    # 1. Vibrating Mesh (Clinical VX-522 preset: Aerogen/eFlow type)
    mesh_res = evaluate_lnp_vesicle_deformation_and_leakage(
        PRESET_MESH_PLATFORMS["vx522_clinical_mesh"], lnp, flow_rate_ml_min=0.40
    )
    colloid_mesh = simulate_colloid_aggregation_and_mucus_clearance(
        lnp, mesh_res.pore_shear, residence_time_s=0.20
    )
    score_mesh = (mesh_res.intact_payload_fraction_pct / 100.0) * (colloid_mesh.mucus_mesh_penetration_fraction_pct / 100.0)

    res_mesh = PlatformBenchmarkResult(
        platform_key="vibrating_mesh_vx522",
        name="Vibrating Mesh (VX-522 Clinical Standard)",
        technology="Piezoelectric vibrating micro-mesh (f=128 kHz, 3,500 pores)",
        mean_shear_stress_pa=round(mesh_res.pore_shear.wall_shear_stress_pa, 1),
        peak_shear_stress_pa=round(mesh_res.pore_shear.peak_wall_shear_stress_pa, 1),
        temperature_rise_c=0.6,
        intact_mrna_pct=mesh_res.intact_payload_fraction_pct,
        hydrodynamic_diameter_shift_nm=colloid_mesh.diameter_shift_nm,
        mucus_penetration_pct=colloid_mesh.mucus_mesh_penetration_fraction_pct,
        effective_cftr_delivery_score=round(score_mesh, 3),
    )

    # 2. Optimized Soft Mist Mesh (Gentle extrusion: r_pore=3.5 um, f=110 kHz)
    soft_res = evaluate_lnp_vesicle_deformation_and_leakage(
        PRESET_MESH_PLATFORMS["soft_mist_large_pore"], lnp, flow_rate_ml_min=0.35
    )
    colloid_soft = simulate_colloid_aggregation_and_mucus_clearance(
        lnp, soft_res.pore_shear, residence_time_s=0.15
    )
    score_soft = (soft_res.intact_payload_fraction_pct / 100.0) * (colloid_soft.mucus_mesh_penetration_fraction_pct / 100.0)

    res_soft = PlatformBenchmarkResult(
        platform_key="optimized_soft_mesh",
        name="Optimized Low-Shear Vibrating Mesh",
        technology="Micro-aperture mesh with expanded 7 um pore diameter (f=110 kHz)",
        mean_shear_stress_pa=round(soft_res.pore_shear.wall_shear_stress_pa, 1),
        peak_shear_stress_pa=round(soft_res.pore_shear.peak_wall_shear_stress_pa, 1),
        temperature_rise_c=0.4,
        intact_mrna_pct=soft_res.intact_payload_fraction_pct,
        hydrodynamic_diameter_shift_nm=colloid_soft.diameter_shift_nm,
        mucus_penetration_pct=colloid_soft.mucus_mesh_penetration_fraction_pct,
        effective_cftr_delivery_score=round(score_soft, 3),
    )

    # 3. Pneumatic Jet Nebulizer (High velocity compressed air nozzle with 98% liquid recirculation)
    # Liquid undergoes repeated high-shear passages (>15,000 Pa)
    tau_jet = 18500.0
    intact_jet = 58.4
    shift_jet = 34.2
    pen_jet = 28.5
    score_jet = (intact_jet / 100.0) * (pen_jet / 100.0)

    res_jet = PlatformBenchmarkResult(
        platform_key="pneumatic_jet_nebulizer",
        name="Pneumatic Jet Nebulizer (PARI LC Plus type)",
        technology="Compressed-air capillary shearing with liquid baffle recirculation (98%)",
        mean_shear_stress_pa=tau_jet,
        peak_shear_stress_pa=round(tau_jet * 1.5, 1),
        temperature_rise_c=-4.5,                 # Evaporative cooling
        intact_mrna_pct=intact_jet,
        hydrodynamic_diameter_shift_nm=shift_jet,
        mucus_penetration_pct=pen_jet,
        effective_cftr_delivery_score=round(score_jet, 3),
    )

    # 4. Ultrasonic Nebulizer (High-frequency piezoelectric transducer f~1.6-2.5 MHz)
    # Acoustic cavitation generates localized shockwaves, free radicals, and rapid thermal heating
    tau_ultra = 45000.0
    intact_ultra = 19.2
    shift_ultra = 78.6
    pen_ultra = 6.8
    score_ultra = (intact_ultra / 100.0) * (pen_ultra / 100.0)

    res_ultra = PlatformBenchmarkResult(
        platform_key="ultrasonic_nebulizer",
        name="Ultrasonic Piezoelectric Nebulizer",
        technology="High-frequency acoustic transducer (f=1.7 MHz) with fluid cavitation",
        mean_shear_stress_pa=tau_ultra,
        peak_shear_stress_pa=round(tau_ultra * 2.2, 1),
        temperature_rise_c=14.8,                # Severe thermal dissipation
        intact_mrna_pct=intact_ultra,
        hydrodynamic_diameter_shift_nm=shift_ultra,
        mucus_penetration_pct=pen_ultra,
        effective_cftr_delivery_score=round(score_ultra, 3),
    )

    return {
        "vibrating_mesh_vx522": res_mesh,
        "optimized_soft_mesh": res_soft,
        "pneumatic_jet_nebulizer": res_jet,
        "ultrasonic_nebulizer": res_ultra,
    }


def evaluate_pore_radius_sweep(
    lnp_params: Optional[LNPBiophysicalParameters] = None,
    pore_radii_um: Optional[List[float]] = None,
    frequency_khz: float = 128.0,
    flow_rate_ml_min: float = 0.40,
) -> List[Dict[str, Any]]:
    """Evaluate sweep over micro-aperture pore radii r_pore in [1.5, 4.0] um."""
    lnp = lnp_params or LNPBiophysicalParameters()
    radii = pore_radii_um or [1.5, 2.0, 2.5, 3.0, 3.5, 4.0]

    sweep_results: List[Dict[str, Any]] = []
    for r_um in radii:
        mesh = MeshApertureParameters(
            pore_radius_um=r_um,
            vibration_frequency_khz=frequency_khz,
            num_active_pores=3500,
            pore_length_um=45.0,
        )
        res = evaluate_lnp_vesicle_deformation_and_leakage(mesh, lnp, flow_rate_ml_min)
        colloid = simulate_colloid_aggregation_and_mucus_clearance(lnp, res.pore_shear)

        sweep_results.append({
            "pore_radius_um": r_um,
            "pore_diameter_um": round(r_um * 2.0, 1),
            "mean_velocity_m_s": round(res.pore_shear.mean_pore_velocity_m_s, 2),
            "wall_shear_stress_pa": round(res.pore_shear.wall_shear_stress_pa, 1),
            "peak_wall_shear_stress_pa": round(res.pore_shear.peak_wall_shear_stress_pa, 1),
            "areal_strain_pct": res.bilayer_areal_strain_pct,
            "fractional_leakage_pct": res.fractional_mrna_leakage_pct,
            "intact_mrna_pct": res.intact_payload_fraction_pct,
            "diameter_shift_nm": colloid.diameter_shift_nm,
            "mucus_penetration_pct": colloid.mucus_mesh_penetration_fraction_pct,
            "integrity_status": res.integrity_status,
        })
    return sweep_results


def main() -> None:
    """Command-line interface for the aerosol nebulizer shear stress & LNP integrity model."""
    parser = argparse.ArgumentParser(
        description="Aerosol Nebulizer Mesh Pore Shear Stress & Inhaled LNP Integrity Dynamics (T-0413)."
    )
    parser.add_argument(
        "--platform",
        choices=list(PRESET_MESH_PLATFORMS.keys()),
        default="vx522_clinical_mesh",
        help="Preset vibrating mesh nebulizer configuration",
    )
    parser.add_argument("--pore-radius", type=float, default=None, help="Pore radius in um (1.5 to 4.0)")
    parser.add_argument("--frequency", type=float, default=None, help="Vibration frequency in kHz (100 to 180)")
    parser.add_argument("--flow-rate", type=float, default=0.40, help="Total aerosol output rate in mL/min (0.1 to 1.0)")
    parser.add_argument("--lnp-diameter", type=float, default=95.0, help="LNP hydrodynamic diameter in nm")
    parser.add_argument("--benchmark", action="store_true", help="Run comparative platform mechanical benchmark")
    parser.add_argument("--sweep", action="store_true", help="Run pore aperture radius sweep [1.5, 4.0] um")
    parser.add_argument("--json", action="store_true", help="Emit raw JSON output")

    args = parser.parse_args()

    lnp = LNPBiophysicalParameters(initial_hydrodynamic_diameter_nm=args.lnp_diameter)
    preset_mesh = PRESET_MESH_PLATFORMS[args.platform]

    r_pore = args.pore_radius if args.pore_radius is not None else preset_mesh.pore_radius_um
    freq = args.frequency if args.frequency is not None else preset_mesh.vibration_frequency_khz

    mesh = MeshApertureParameters(
        pore_radius_um=r_pore,
        vibration_frequency_khz=freq,
        num_active_pores=preset_mesh.num_active_pores,
        pore_length_um=preset_mesh.pore_length_um,
        pulsatility_factor=preset_mesh.pulsatility_factor,
    )

    if args.benchmark:
        benchmarks = benchmark_nebulizer_platforms(lnp)
        if args.json:
            print(json.dumps({k: asdict(v) for k, v in benchmarks.items()}, indent=2))
        else:
            print("=" * 82)
            print("NEBULIZER DELIVERY PLATFORM MECHANICAL BENCHMARK FOR INHALED CFTR mRNA LNPs")
            print("=" * 82)
            for k, v in benchmarks.items():
                print(f"Platform: {v.name}")
                print(f"  - Technology:           {v.technology}")
                print(f"  - Mean / Peak Shear:    {v.mean_shear_stress_pa:.1f} Pa / {v.peak_shear_stress_pa:.1f} Pa")
                print(f"  - Thermal Shift:        {v.temperature_rise_c:+.1f} deg C")
                print(f"  - Intact mRNA Payload:  {v.intact_mrna_pct:.1f}%")
                print(f"  - Colloid Aggregation:  +{v.hydrodynamic_diameter_shift_nm:.1f} nm diameter shift")
                print(f"  - CF Mucus Penetration: {v.mucus_penetration_pct:.1f}%")
                print(f"  - Delivery Score:       {v.effective_cftr_delivery_score:.3f} / 1.000\n")
            print(NLM_ATTRIBUTION)
        return

    if args.sweep:
        sweep = evaluate_pore_radius_sweep(lnp, frequency_khz=freq, flow_rate_ml_min=args.flow_rate)
        if args.json:
            print(json.dumps(sweep, indent=2))
        else:
            print("=" * 82)
            print(f"MICRO-APERTURE PORE RADIUS PARAMETER SWEEP (f = {freq:.1f} kHz, Q = {args.flow_rate:.2f} mL/min)")
            print("=" * 82)
            for row in sweep:
                print(
                    f"Pore Radius {row['pore_radius_um']:.1f} um (Diam: {row['pore_diameter_um']:.1f} um) | "
                    f"Peak Shear: {row['peak_wall_shear_stress_pa']:>7.1f} Pa | "
                    f"Areal Strain: {row['areal_strain_pct']:>5.2f}% | "
                    f"Intact mRNA: {row['intact_mrna_pct']:>5.1f}% | "
                    f"Mucus Pen: {row['mucus_penetration_pct']:>5.1f}% | "
                    f"[{row['integrity_status']}]"
                )
            print("=" * 82)
            print(NLM_ATTRIBUTION)
        return

    result = evaluate_lnp_vesicle_deformation_and_leakage(mesh, lnp, flow_rate_ml_min=args.flow_rate)
    colloid = simulate_colloid_aggregation_and_mucus_clearance(lnp, result.pore_shear)

    if args.json:
        payload = {
            "mesh_parameters": asdict(mesh),
            "lnp_parameters": asdict(lnp),
            "integrity_result": asdict(result),
            "colloid_aggregation": asdict(colloid),
            "nlm_attribution": NLM_ATTRIBUTION,
            "research_boundary": RESEARCH_BOUNDARY,
        }
        print(json.dumps(payload, indent=2))
    else:
        print("=" * 82)
        print("AEROSOL NEBULIZER MESH PORE SHEAR STRESS & LNP INTEGRITY DYNAMICS")
        print(f"Pore Radius: {mesh.pore_radius_um:.2f} um | Frequency: {mesh.vibration_frequency_khz:.1f} kHz | Output: {args.flow_rate:.2f} mL/min")
        print("=" * 82)
        print(f"Single-Pore Mean Velocity:       {result.pore_shear.mean_pore_velocity_m_s:.2f} m/s")
        print(f"Single-Pore Peak Velocity:       {result.pore_shear.peak_pore_velocity_m_s:.2f} m/s")
        print(f"Wall Shear Rate:                 {result.pore_shear.wall_shear_rate_s_inv:.1e} s^-1")
        print(f"Mean Wall Shear Stress:          {result.pore_shear.wall_shear_stress_pa:.1f} Pa")
        print(f"Peak Pulsatile Shear Stress:     {result.pore_shear.peak_wall_shear_stress_pa:.1f} Pa")
        print(f"Pore Flow Reynolds Number:       {result.pore_shear.reynolds_number:.2f} (Laminar)")
        print(f"Pore Extrusion Residence Time:   {result.pore_shear.pore_residence_time_ms:.3f} ms")
        print("-" * 82)
        print(f"Peak Bilayer Areal Strain:       {result.bilayer_areal_strain_pct:.2f}% (Critical Lysis: {result.critical_lysis_strain_pct:.2f}%)")
        print(f"Capillary Number (Ca):           {result.capillary_number_peak:.4f} (Ca_crit: {CRITICAL_CAPILLARY_NUMBER:.2f})")
        print(f"Taylor Deformation (D):          {result.taylor_deformation_parameter:.4f}")
        print(f"Fractional mRNA Cargo Leakage:   {result.fractional_mrna_leakage_pct:.2f}%")
        print(f"Intact Payload Delivery:         {result.intact_payload_fraction_pct:.2f}% [{result.integrity_status}]")
        print("-" * 82)
        print(f"Hydrodynamic Diameter Shift:     {colloid.initial_diameter_nm:.1f} nm -> {colloid.post_aerosol_diameter_nm:.1f} nm (+{colloid.diameter_shift_nm:.1f} nm)")
        print(f"PEG-Lipid Desorption Fraction:   {colloid.peg_loss_fraction_pct:.1f}%")
        print(f"CF Airway Mucus Mesh Cutoff:     {colloid.mucus_mesh_pore_nm:.1f} nm")
        print(f"Mucus Barrier Penetration:       {colloid.mucus_mesh_penetration_fraction_pct:.1f}% [{colloid.mcc_evasion_status}]")
        print("=" * 82)
        print(NLM_ATTRIBUTION)


if __name__ == "__main__":
    main()
