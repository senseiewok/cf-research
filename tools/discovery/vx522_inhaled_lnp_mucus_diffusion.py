#!/usr/bin/env python3
"""Inhaled mRNA Lipid Nanoparticle Mucus Mesh Diffusion & Surface Charge Neutralization Model (T-0400).

Quantitative biophysical model of aerosolized CFTR mRNA lipid nanoparticle (LNP) transport,
mucin mesh steric sieving, surface charge electrostatic neutralization, PEG-lipid shedding,
and pH-dependent endosomal escape for inhaled cystic fibrosis gene therapies:
  1. Mucus Mesh Steric Obstruction & Ogston/Renkin Sieving:
     Formulates combined Ogston steric sieving and Renkin hydrodynamic hindrance scaling
     relating LNP diameter (d_lnp approx 80 nm) to heterogeneous airway mucus pore sizes
     (d_pore in [100, 200] nm):
       D_0 = (k_B * T) / (6 * pi * eta_0 * r_lnp)
       D_steric / D_0 = exp(-pi * (r_lnp + r_f)^2 / xi_mesh^2) * F(lambda)
       F(lambda) = (1 - lambda)^2 * (1 - 2.104*lambda + 2.09*lambda^3 - 0.95*lambda^5)
  2. Surface Charge Neutralization & DLVO Electrostatic Trap Energy:
     Simulates double-layer electrostatic interaction energy between LNPs and negatively
     charged sialylated/sulfated mucin fibers (psi_mucus approx -15 mV). Near-neutral or
     zwitterionic surface charge (zeta in [-2.0, +2.0] mV) maintains |U_elec| < 1.0 k_B*T,
     preventing adhesive entrapment:
       U_elec(h) = 2 * pi * epsilon_r * epsilon_0 * r_lnp * psi_mucus * zeta_lnp * exp(-kappa * h)
  3. PEG-Lipid Shedding Kinetics:
     Models shedding of hydrophilic PEG-lipids (k_shed in [0.05, 0.20] h^-1) to permit
     apolipoprotein adsorption and endocytic uptake following mucus transit:
       [PEG_bound](t) = [PEG_0] * exp(-k_shed * t)
  4. Ionizable Lipid pKa Tuning & Endosomal Escape:
     Formulates pH-dependent protonation alpha_prot(pH) switching from < 15% at extracellular
     pH 7.2 to > 85% at late endosomal pH 5.5 for ionizable lipids with pKa in [6.2, 6.5],
     driving inverted hexagonal H_II phase transition and endosomal escape.

Primary Evidence Grounding:
  - Discovery Brief: cf-research/landscape/mission-grounding-2026-10-23-wave17.md
  - Primary Literature:
    - Suk JS et al. (2019) PNAS, DOI: 10.1073/pnas.1901667116
    - Kim J et al. (2021) Nat Biomed Eng, DOI: 10.1038/s41551-021-00724-5
  - Clinical Registry Grounding:
    - NCT06237335: Phase 1/2 study of aerosolized CFTR mRNA LNP RCT2100 in healthy participants and CF

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
    "Preclinical in-silico inhaled mRNA LNP mucus diffusion & endosomal escape biophysical model; "
    "research, not medical advice. No clinical efficacy or therapeutic safety claim is made."
)

# Reference Grounding Identifiers
PRIMARY_DOI_SUK = "10.1073/pnas.1901667116"       # Suk et al. PNAS 2019
PRIMARY_DOI_KIM = "10.1038/s41551-021-00724-5"     # Kim et al. Nat Biomed Eng 2021
CLINICAL_TRIAL_REGISTRY = "NCT06237335"            # RCT2100 mRNA LNP Phase 1/2

# Fundamental physical constants
BOLTZMANN_CONSTANT = 1.380649e-23                  # J / K
VACUUM_PERMITTIVITY = 8.8541878128e-12             # F / m
DEFAULT_TEMPERATURE_K = 310.15                     # Physiological body temperature (37 deg C)
WATER_VISCOSITY_PA_S = 0.6913e-3                   # Dynamic viscosity of water at 37 deg C (Pa*s)
RELATIVE_PERMITTIVITY_WATER = 74.3                 # Dielectric constant of water at 37 deg C
DEBYE_LENGTH_M = 0.8e-9                            # Physiological ionic strength (150 mM NaCl) -> kappa^-1 approx 0.8 nm
KAPPA_M_INV = 1.0 / DEBYE_LENGTH_M                 # Inverse Debye screening length (m^-1)

# Biological default parameters for CF mucus and LNPs
DEFAULT_LNP_DIAMETER_NM = 80.0                     # Nominal inhaled mRNA LNP hydrodynamic diameter (nm)
MUCIN_FIBER_RADIUS_NM = 3.5                        # Hydrodynamic radius of individual mucin polymer fiber (nm)
DEFAULT_CF_MUCUS_PORE_NM = 200.0                   # Average mesh pore diameter of CF respiratory mucus (nm)
MUCUS_SURFACE_POTENTIAL_V = -0.015                 # Effective negative surface potential of sialylated mucins (-15 mV)
CLOSEST_APPROACH_DISTANCE_M = 0.5e-9               # Surface-to-surface separation at collision (0.5 nm)
DEFAULT_MUCUS_THICKNESS_UM = 15.0                  # CF airway mucus plaque thickness (um)
DEFAULT_MCC_CLEARANCE_HALF_LIFE_S = 1200.0         # Mucociliary clearance half-life in obstructed airways (20 min)


@dataclass(frozen=True)
class LNPFormulation:
    """Biophysical formulation specification for an aerosolized mRNA lipid nanoparticle."""

    name: str
    description: str
    diameter_nm: float
    zeta_potential_mv: float
    peg_mol_pct: float
    peg_shedding_rate_per_hour: float
    ionizable_lipid_pka: float
    ionizable_lipid_name: str


@dataclass
class DiffusionProfileResult:
    """Quantitative diffusion, electrostatic energy, and transport kinetics result."""

    formulation: str
    diameter_nm: float
    zeta_potential_mv: float
    d_0_um2_s: float
    ogston_sieving_factor: float
    renkin_hindrance_factor: float
    d_steric_ratio: float
    electrostatic_energy_kb_t: float
    mobile_fraction: float
    d_eff_um2_s: float
    d_eff_over_d0: float
    mean_transit_time_seconds: float
    mucus_penetration_efficiency_pct: float
    shedding_half_life_hours: float
    uptake_competence_4h_pct: float
    protonation_extracellular_pct: float
    protonation_endosomal_pct: float
    endosomal_escape_probability_pct: float
    mucus_penetrating_particle_status: bool
    curative_delivery_competent: bool


# Canonical formulation library anchored in preclinical & clinical literature
CANONICAL_FORMULATIONS: Dict[str, LNPFormulation] = {
    "VX-522_Aerosol_Benchmark": LNPFormulation(
        name="VX-522_Aerosol_Benchmark",
        description="Clinical-stage inhaled CFTR mRNA LNP (RCT2100 / VX-522 proxy) with dense shedding PEG and tuned ionizable lipid",
        diameter_nm=80.0,
        zeta_potential_mv=-0.5,
        peg_mol_pct=3.0,
        peg_shedding_rate_per_hour=0.12,
        ionizable_lipid_pka=6.38,
        ionizable_lipid_name="Proprietary Ionizable Lipid (VX-522)",
    ),
    "Cationic_Unshielded_LNP": LNPFormulation(
        name="Cationic_Unshielded_LNP",
        description="Standard cationic transfection lipid (DOTAP/DODAP-like) without steric PEG shielding",
        diameter_nm=85.0,
        zeta_potential_mv=18.5,
        peg_mol_pct=0.5,
        peg_shedding_rate_per_hour=0.01,
        ionizable_lipid_pka=7.80,
        ionizable_lipid_name="DOTAP / Cationic",
    ),
    "Dense_NonShedding_PEG_LNP": LNPFormulation(
        name="Dense_NonShedding_PEG_LNP",
        description="Dense DSPE-PEG2000 non-shedding lipid nanoparticle (diffuses through mucus but blocked from cellular uptake)",
        diameter_nm=80.0,
        zeta_potential_mv=-1.2,
        peg_mol_pct=5.0,
        peg_shedding_rate_per_hour=0.005,
        ionizable_lipid_pka=6.40,
        ionizable_lipid_name="MC3 / Ionizable",
    ),
    "Oversized_Aerosol_LNP": LNPFormulation(
        name="Oversized_Aerosol_LNP",
        description="Oversized LNP formulation exceeding airway mucus pore dimensions (sterically filtered by Ogston sieving)",
        diameter_nm=220.0,
        zeta_potential_mv=-0.8,
        peg_mol_pct=3.0,
        peg_shedding_rate_per_hour=0.10,
        ionizable_lipid_pka=6.35,
        ionizable_lipid_name="SM-102 Analog",
    ),
    "Anionic_Charged_LNP": LNPFormulation(
        name="Anionic_Charged_LNP",
        description="Strongly anionic LNP experiencing electrostatic repulsion and hydration barrier",
        diameter_nm=80.0,
        zeta_potential_mv=-16.0,
        peg_mol_pct=1.5,
        peg_shedding_rate_per_hour=0.05,
        ionizable_lipid_pka=6.10,
        ionizable_lipid_name="DOPG / Anionic",
    ),
}


def stokes_einstein_diffusivity(
    diameter_nm: float,
    temperature_k: float = DEFAULT_TEMPERATURE_K,
    viscosity_pa_s: float = WATER_VISCOSITY_PA_S,
) -> float:
    """Calculate unhindered Brownian diffusion coefficient D_0 in m^2 / s."""
    radius_m = (diameter_nm * 1e-9) / 2.0
    d_0_m2_s = (BOLTZMANN_CONSTANT * temperature_k) / (6.0 * math.pi * viscosity_pa_s * radius_m)
    return d_0_m2_s


def ogston_sieving_factor(
    diameter_nm: float,
    pore_diameter_nm: float = DEFAULT_CF_MUCUS_PORE_NM,
    fiber_radius_nm: float = MUCIN_FIBER_RADIUS_NM,
) -> float:
    """Calculate Ogston steric mesh obstruction factor S_Ogston.

    S_Ogston = exp(-pi * (r_lnp + r_f)^2 / xi_mesh^2)
    where xi_mesh = pore_diameter_nm / 2.0.
    """
    r_lnp = diameter_nm / 2.0
    xi_mesh = pore_diameter_nm / 2.0
    exponent = -math.pi * ((r_lnp + fiber_radius_nm) ** 2) / (xi_mesh ** 2)
    return math.exp(exponent)


def renkin_hindrance_polynomial(
    diameter_nm: float,
    pore_diameter_nm: float = DEFAULT_CF_MUCUS_PORE_NM,
) -> float:
    """Calculate Renkin hydrodynamic pore hindrance polynomial F(lambda).

    lambda = d_lnp / d_pore
    F(lambda) = (1 - lambda)^2 * (1 - 2.104*lambda + 2.09*lambda^3 - 0.95*lambda^5)
    If lambda >= 1.0, F(lambda) = 0.0.
    """
    lambda_val = diameter_nm / pore_diameter_nm
    if lambda_val >= 1.0:
        return 0.0
    term1 = (1.0 - lambda_val) ** 2
    term2 = 1.0 - 2.104 * lambda_val + 2.09 * (lambda_val ** 3) - 0.95 * (lambda_val ** 5)
    return max(0.0, term1 * term2)


def dlvo_electrostatic_trap_energy(
    diameter_nm: float,
    zeta_potential_mv: float,
    mucus_potential_v: float = MUCUS_SURFACE_POTENTIAL_V,
    separation_m: float = CLOSEST_APPROACH_DISTANCE_M,
    temperature_k: float = DEFAULT_TEMPERATURE_K,
) -> float:
    """Calculate electrostatic double-layer interaction energy in units of k_B * T.

    U_elec = 2 * pi * epsilon_r * epsilon_0 * r_lnp * psi_mucus * zeta_lnp * exp(-kappa * h)
    Note: Both potentials in Volts. When zeta > 0 and psi < 0, U_elec is negative (attractive/adhesive).
    """
    radius_m = (diameter_nm * 1e-9) / 2.0
    zeta_v = zeta_potential_mv * 1e-3
    epsilon = RELATIVE_PERMITTIVITY_WATER * VACUUM_PERMITTIVITY
    u_elec_joules = (
        2.0
        * math.pi
        * epsilon
        * radius_m
        * mucus_potential_v
        * zeta_v
        * math.exp(-KAPPA_M_INV * separation_m)
    )
    kb_t_joules = BOLTZMANN_CONSTANT * temperature_k
    energy_kb_t = u_elec_joules / kb_t_joules
    return energy_kb_t


def calculate_mobile_fraction(energy_kb_t: float) -> float:
    """Calculate mobile non-adhered nanoparticle fraction based on electrostatic trap energy.

    If energy_kb_t >= 0.0, electrostatic repulsion/neutrality ensures zero adhesive trapping (f_mobile = 1.0).
    If energy_kb_t < 0.0, attractive adhesion occurs. When the trap depth is <= 1.0 k_B*T,
    thermal energy dominates and particles remain largely mobile (> 0.80).
    When trap depth exceeds 1.0 k_B*T, particles become immobilized exponentially.
    """
    if energy_kb_t >= 0.0:
        return 1.0
    trap_depth = -energy_kb_t
    if trap_depth <= 1.0:
        return 1.0 - 0.2 * trap_depth
    return max(0.0001, 0.80 * math.exp(-(trap_depth - 1.0) * 1.5))


def peg_shedding_kinetics(
    shedding_rate_per_hour: float,
    time_hours: float,
) -> Tuple[float, float]:
    """Calculate bound PEG fraction and cellular uptake competence at time t.

    bound_fraction = exp(-k_shed * t)
    uptake_competence = 1.0 - bound_fraction
    """
    bound_frac = math.exp(-shedding_rate_per_hour * time_hours)
    uptake_comp = 1.0 - bound_frac
    return bound_frac, uptake_comp


def ionizable_lipid_protonation(
    ph: float,
    pka: float,
) -> float:
    """Calculate fractional protonation alpha_prot of ionizable lipid via Henderson-Hasselbalch.

    alpha_prot = 1.0 / (1.0 + 10^(pH - pKa))
    """
    delta = ph - pka
    if delta > 10.0:
        return 0.0
    if delta < -10.0:
        return 1.0
    return 1.0 / (1.0 + (10.0 ** delta))


def endosomal_escape_probability(
    protonation_fraction: float,
    hill_coefficient: float = 2.5,
    half_activation_fraction: float = 0.50,
    max_escape_prob: float = 0.32,
) -> float:
    """Model cooperative endosomal escape probability as a function of lipid protonation.

    P_escape = max_escape_prob * (alpha^n) / (alpha_50^n + alpha^n)
    """
    if protonation_fraction <= 0.0:
        return 0.0
    term = (protonation_fraction / half_activation_fraction) ** hill_coefficient
    prob = max_escape_prob * (term / (1.0 + term))
    return prob


def evaluate_lnp_diffusion_profile(
    formulation: LNPFormulation,
    mucus_pore_diameter_nm: float = DEFAULT_CF_MUCUS_PORE_NM,
    mucus_thickness_um: float = DEFAULT_MUCUS_THICKNESS_UM,
    mcc_half_life_s: float = DEFAULT_MCC_CLEARANCE_HALF_LIFE_S,
    extracellular_ph: float = 7.2,
    endosomal_ph: float = 5.5,
) -> DiffusionProfileResult:
    """Evaluate full transport, steric sieving, charge neutralization, and uptake kinetics."""
    d_0_m2_s = stokes_einstein_diffusivity(formulation.diameter_nm)
    d_0_um2_s = d_0_m2_s * 1e12

    ogston = ogston_sieving_factor(formulation.diameter_nm, pore_diameter_nm=mucus_pore_diameter_nm)
    renkin = renkin_hindrance_polynomial(formulation.diameter_nm, pore_diameter_nm=mucus_pore_diameter_nm)
    d_steric_ratio = ogston * renkin

    u_elec = dlvo_electrostatic_trap_energy(
        formulation.diameter_nm,
        formulation.zeta_potential_mv,
    )
    f_mobile = calculate_mobile_fraction(u_elec)

    d_eff_um2_s = d_0_um2_s * d_steric_ratio * f_mobile
    d_eff_over_d0 = (d_eff_um2_s / d_0_um2_s) if d_0_um2_s > 0 else 0.0

    # Mean transit time through mucus plaque: tau = L^2 / (2 * D_eff)
    l_mucus_m = mucus_thickness_um * 1e-6
    d_eff_m2_s = d_eff_um2_s * 1e-12
    if d_eff_m2_s > 1e-18:
        tau_transit_s = (l_mucus_m ** 2) / (2.0 * d_eff_m2_s)
    else:
        tau_transit_s = 1e8  # Effectively immobilized

    # Penetration efficiency before clearance: P = exp(-tau_transit / tau_mcc)
    # Using clearance rate constant k_mcc = ln(2) / mcc_half_life_s
    k_mcc = math.log(2.0) / mcc_half_life_s
    penetration_prob = math.exp(-k_mcc * tau_transit_s)
    penetration_pct = penetration_prob * 100.0

    # Shedding kinetics
    shed_rate = formulation.peg_shedding_rate_per_hour
    shed_half_life_h = (math.log(2.0) / shed_rate) if shed_rate > 0 else 1000.0
    _, uptake_comp_4h = peg_shedding_kinetics(shed_rate, time_hours=4.0)
    uptake_comp_pct = uptake_comp_4h * 100.0

    # Endosomal escape kinetics
    prot_extra = ionizable_lipid_protonation(extracellular_ph, formulation.ionizable_lipid_pka)
    prot_endo = ionizable_lipid_protonation(endosomal_ph, formulation.ionizable_lipid_pka)
    p_escape = endosomal_escape_probability(prot_endo)
    p_escape_pct = p_escape * 100.0

    # Criteria for Mucus-Penetrating Particle (MPP) status
    # Per Suk et al. 2019: near-neutral charge (|zeta| <= 2.5 mV) and d_eff / d_0 >= 0.02
    is_mpp = abs(formulation.zeta_potential_mv) <= 2.5 and d_eff_over_d0 >= 0.02

    # Curative delivery competence:
    # Requires MPP status, penetration efficiency >= 10%, uptake competence >= 25% at 4h, and p_escape >= 10%
    is_competent = (
        is_mpp
        and penetration_pct >= 10.0
        and uptake_comp_pct >= 25.0
        and p_escape_pct >= 10.0
    )

    return DiffusionProfileResult(
        formulation=formulation.name,
        diameter_nm=formulation.diameter_nm,
        zeta_potential_mv=formulation.zeta_potential_mv,
        d_0_um2_s=round(d_0_um2_s, 3),
        ogston_sieving_factor=round(ogston, 4),
        renkin_hindrance_factor=round(renkin, 4),
        d_steric_ratio=round(d_steric_ratio, 4),
        electrostatic_energy_kb_t=round(u_elec, 3),
        mobile_fraction=round(f_mobile, 4),
        d_eff_um2_s=round(d_eff_um2_s, 4),
        d_eff_over_d0=round(d_eff_over_d0, 4),
        mean_transit_time_seconds=round(tau_transit_s, 1),
        mucus_penetration_efficiency_pct=round(penetration_pct, 2),
        shedding_half_life_hours=round(shed_half_life_h, 2),
        uptake_competence_4h_pct=round(uptake_comp_pct, 2),
        protonation_extracellular_pct=round(prot_extra * 100.0, 2),
        protonation_endosomal_pct=round(prot_endo * 100.0, 2),
        endosomal_escape_probability_pct=round(p_escape_pct, 2),
        mucus_penetrating_particle_status=is_mpp,
        curative_delivery_competent=is_competent,
    )


def run_benchmark_comparison() -> Dict[str, Any]:
    """Compare diffusion, steric sieving, charge trapping, and escape across all canonical formulations."""
    comparison: Dict[str, Any] = {}
    for name, formulation in CANONICAL_FORMULATIONS.items():
        res = evaluate_lnp_diffusion_profile(formulation)
        comparison[name] = asdict(res)
    return comparison


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Inhaled mRNA Lipid Nanoparticle Mucus Mesh Diffusion & Surface Charge Neutralization Model (T-0400)"
    )
    parser.add_argument(
        "--formulation",
        choices=list(CANONICAL_FORMULATIONS.keys()),
        default="VX-522_Aerosol_Benchmark",
        help="Target LNP formulation (default: VX-522_Aerosol_Benchmark)",
    )
    parser.add_argument(
        "--pore-size",
        type=float,
        default=DEFAULT_CF_MUCUS_PORE_NM,
        help=f"Airway mucus mesh pore diameter in nm (default: {DEFAULT_CF_MUCUS_PORE_NM})",
    )
    parser.add_argument(
        "--benchmark-all",
        action="store_true",
        help="Run comparative benchmark across all cataloged aerosol LNP formulations",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Emit results in machine-readable JSON format",
    )
    args = parser.parse_args()

    if args.benchmark_all:
        bench_data = run_benchmark_comparison()
        output = {
            "model": "vx522_inhaled_lnp_mucus_diffusion",
            "task": "T-0400",
            "tier": "T2",
            "nlm_attribution": NLM_ATTRIBUTION,
            "research_boundary": RESEARCH_BOUNDARY,
            "grounding": {
                "primary_literature": [PRIMARY_DOI_SUK, PRIMARY_DOI_KIM],
                "clinical_trials": [CLINICAL_TRIAL_REGISTRY],
            },
            "formulation_comparison": bench_data,
        }
        if args.json:
            print(json.dumps(output, indent=2))
        else:
            print("=" * 95)
            print("INHALED CFTR mRNA LNP MUCUS DIFFUSION & ENDOSOMAL ESCAPE BENCHMARK")
            print("=" * 95)
            for form_name, data in bench_data.items():
                print(
                    f"{form_name:28s} | d={data['diameter_nm']:3.0f}nm | zeta={data['zeta_potential_mv']:+5.1f}mV | "
                    f"D_eff/D0={data['d_eff_over_d0']:6.4f} | "
                    f"Penetration={data['mucus_penetration_efficiency_pct']:5.1f}% | "
                    f"Escape={data['endosomal_escape_probability_pct']:4.1f}% | "
                    f"[{'COMPETENT' if data['curative_delivery_competent'] else 'FAILED'}]"
                )
            print("=" * 95)
        return

    # Single formulation evaluation
    form = CANONICAL_FORMULATIONS[args.formulation]
    res = evaluate_lnp_diffusion_profile(form, mucus_pore_diameter_nm=args.pore_size)
    payload = {
        "model": "vx522_inhaled_lnp_mucus_diffusion",
        "task": "T-0400",
        "tier": "T2",
        "formulation": args.formulation,
        "nlm_attribution": NLM_ATTRIBUTION,
        "research_boundary": RESEARCH_BOUNDARY,
        "grounding": {
            "primary_literature": [PRIMARY_DOI_SUK, PRIMARY_DOI_KIM],
            "clinical_trials": [CLINICAL_TRIAL_REGISTRY],
        },
        "results": asdict(res),
    }

    if args.json:
        print(json.dumps(payload, indent=2))
    else:
        print("=" * 80)
        print(f"Inhaled mRNA LNP Diffusion & Transport Profile: {args.formulation}")
        print(f"Primary DOIs: {PRIMARY_DOI_SUK}, {PRIMARY_DOI_KIM}")
        print(f"Clinical Trial Registry: {CLINICAL_TRIAL_REGISTRY}")
        print("=" * 80)
        print(f"Hydrodynamic Diameter:         {res.diameter_nm:.1f} nm")
        print(f"Surface Zeta Potential:        {res.zeta_potential_mv:+.2f} mV")
        print(f"Unhindered Diffusivity (D0):   {res.d_0_um2_s:.3f} um^2/s")
        print(f"Steric Hindrance (D_steric/D0):{res.d_steric_ratio:.4f} (Ogston: {res.ogston_sieving_factor:.4f}, Renkin: {res.renkin_hindrance_factor:.4f})")
        print(f"Electrostatic Trap Energy:     {res.electrostatic_energy_kb_t:+.3f} k_B*T")
        print(f"Mobile Unadhered Fraction:     {res.mobile_fraction * 100.0:.1f}%")
        print(f"Effective Diffusivity (D_eff): {res.d_eff_um2_s:.4f} um^2/s ({res.d_eff_over_d0 * 100.0:.2f}% of D0)")
        print(f"Mean Mucus Transit Time:       {res.mean_transit_time_seconds:.1f} seconds")
        print(f"Penetration Before Clearance:  {res.mucus_penetration_efficiency_pct:.2f}%")
        print(f"PEG Shedding Half-Life:        {res.shedding_half_life_hours:.2f} hours")
        print(f"Uptake Competence (at 4h):     {res.uptake_competence_4h_pct:.1f}%")
        print(f"Lipid Protonation (pH 7.2):    {res.protonation_extracellular_pct:.1f}%")
        print(f"Lipid Protonation (pH 5.5):    {res.protonation_endosomal_pct:.1f}%")
        print(f"Endosomal Escape Probability:  {res.endosomal_escape_probability_pct:.1f}%")
        print(f"MPP Classification:            {'PASS' if res.mucus_penetrating_particle_status else 'FAIL'}")
        print(f"Curative Delivery Status:      {'COMPETENT' if res.curative_delivery_competent else 'SUB-THRESHOLD'}")
        print("=" * 80)


if __name__ == "__main__":
    main()
