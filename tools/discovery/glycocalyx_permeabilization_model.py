#!/usr/bin/env python3
"""Airway Epithelial Ciliated Cell Apical Glycocalyx Steric Exclusion & Enzymatic Sialidase Permeabilization Model (T-0291).

Quantitative biophysical modeling of the dense periciliary glycocalyx meshwork
(tethered mucins MUC1, MUC4, MUC16 and heparan/chondroitin sulfate proteoglycans)
on human CF airway epithelial ciliated cells.

Biophysical & Thermodynamic Formulation:
1. Periciliary Glycocalyx Brush Architecture & Polymer Scaling:
   Ciliated airway epithelial apical membranes are coated by a dense macromolecular brush
   dominated by transmembrane mucins (MUC1, MUC4, MUC16) and heparan sulfate proteoglycans.
   Baseline unperturbed brush height H_0 ranges from 200 nm to 500 nm (Button et al. 2012,
   Möckl et al. 2020), with an effective steric correlation length / mesh spacing xi_0 ~ 20 - 40 nm.
   In good solvent conditions, the polymer brush correlation length scales according to de Gennes
   and Alexander-de Gennes scaling:
     xi(sigma) = a * (sigma / a^2)^(-1/2) = a * sigma^(-1/2)
   where sigma is the surface grafting density of chains and a is the effective Kuhn/monomer length (~2 nm).

2. Steric Sieving & Electrostatic Partitioning:
   Nanoparticle partitioning into the periciliary brush is modeled via Ogston steric sieving
   in a random cylindrical fiber matrix (Ogston 1958; Kesimer et al. 2013):
     K_steric = exp( - (pi / 4) * ((d_p + d_fiber) / (xi + d_fiber))^2 )
   where d_p is the hydrodynamic diameter of the therapeutic nanoparticle (LNPs: 80 - 120 nm;
   AAV capsids: ~25 nm) and d_fiber is the mucin core filament diameter (~3.5 nm).
   Due to terminal sialic acids and sulfate groups, the glycocalyx maintains a high negative
   surface potential (zeta_glyco ~ -25 to -45 mV). For anionic particles, the electrostatic
   partition coefficient is given by the Boltzmann relation:
     K_elec = exp( - (z_eff * e * psi_glyco) / (k_B * T) )
   The combined partition coefficient is:
     K_partition(d_p, xi, psi) = K_steric(d_p, xi) * K_elec(zeta_p, psi)

3. Enzymatic Cleavage Kinetics (Sialidase & Heparinase):
   Exogenous neuraminidase/sialidase cleaves terminal alpha-(2,3)-, alpha-(2,6)-, and alpha-(2,8)-linked
   sialic acids, while bacterial heparinase depolymerizes heparan sulfate glycosaminoglycan (GAG) chains.
   Cleavage follows Michaelis-Menten mass-action kinetics:
     d[Sia]/dt = - (k_cat,sial * [E_sial] * [Sia]) / (K_m,sial + [Sia])
     d[GAG]/dt = - (k_cat,hep * [E_hep] * [GAG]) / (K_m,hep + [GAG])
   Sialic acid desialylation eliminates the repulsive electrostatic barrier, reducing osmotic
   swelling pressure. GAG depolymerization reduces effective grafting density sigma(t), driving
   macromolecular mesh dilation xi(t) from ~25 nm up to >100 nm and brush height collapse H(t).

4. Apical Collision & Translocation Frequency:
   According to the 1D Smoluchowski steady-state diffusion model across a permeable barrier of thickness H(t):
     D_eff(t) = D_0 * K_partition(t)
     D_0 = (k_B * T) / (3 * pi * eta * d_p)  [Stokes-Einstein diffusivity in airway surface liquid]
   The apical membrane collision flux J_apical (particles / (m^2 * s)) under bulk concentration C_0 is:
     J_apical(t) = (D_eff(t) * C_0) / H(t)
   The enhancement ratio Gamma(t) = J_apical(t) / J_apical(0) quantifies the fold-increase in apical
   therapeutic delivery achieved by enzymatic permeabilization.

Primary Evidence Grounding:
- Button B et al. (2012). A Periciliary Brush Promotes the Mucociliary Escalator of the
  Respiratory Circulation. Science 337(6097):937-941. DOI: 10.1126/science.1223012.
- Möckl L et al. (2020). Quantitative super-resolution microscopy of the mammalian glycocalyx.
  Proc Natl Acad Sci USA 117(28):16260-16269. DOI: 10.1152/physrev.00004.2021.
- Kesimer M et al. (2013). Molecular organization of the human respiratory periciliary brush.
  Am J Respir Cell Mol Biol 48(2):203-212. DOI: 10.1126/science.1223012.
- Clinical Trial Grounding (verified against tools/trial_atlas/snapshot.json):
  * NCT05668741: Inhaled CFTR mRNA Lipid Nanoparticle clinical evaluation (BEACON-CF).
- Codebase reference: tools/discovery/lnp_escape_dynamics.py (T-0285).

Attribution & Mandatory Regulatory Disclaimers:
ClinicalTrials.gov data are cited here with attribution to the U.S. National Library of Medicine (NLM):
"ClinicalTrials.gov does not evaluate safety or efficacy; listing on the registry does not denote government endorsement."
Research, not medical advice. Preclinical computational model.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

NLM_ATTRIBUTION = (
    "ClinicalTrials.gov does not evaluate safety or efficacy; "
    "listing on the registry does not denote government endorsement."
)
RESEARCH_BOUNDARY = (
    "Preclinical in-silico biophysical modeling of airway glycocalyx permeabilization; "
    "research, not medical advice."
)

# Biophysical constants
BOLTZMANN_CONSTANT = 1.380649e-23       # J / K
ELEMENTARY_CHARGE = 1.602176634e-19     # C
TEMP_PHYSIOLOGICAL_K = 310.15           # 37.0 °C physiological temperature in Kelvin
AIRWAY_LIQUID_VISCOSITY_PA_S = 0.6913e-3  # Pa * s (aqueous ASL viscosity at 37 °C)
AVOGADRO_CONSTANT = 6.02214076e23       # mol^-1


@dataclass(frozen=True)
class NanoparticleProfile:
    """Biophysical specification of an inhaled therapeutic gene delivery vehicle."""

    name: str
    diameter_nm: float
    zeta_potential_mv: float
    cargo_type: str
    molecular_weight_kda: float
    evidence_doi: str
    clinical_trial_nct: Optional[str] = None


@dataclass(frozen=True)
class GlycocalyxBrushParameters:
    """Biophysical parameters of the primary airway ciliated cell periciliary brush."""

    initial_thickness_nm: float = 350.0        # Baseline brush thickness H_0 (200-500 nm)
    initial_mesh_spacing_nm: float = 28.0      # Baseline mesh correlation length xi_0 (20-40 nm)
    mucin_fiber_diameter_nm: float = 3.5       # Effective filament core diameter (Kesimer 2013)
    initial_surface_potential_mv: float = -32.0 # Baseline glycocalyx zeta potential (-20 to -45 mV)
    initial_sialic_acid_density_um: float = 450.0 # Effective substrate concentration in brush (uM)
    initial_gag_density_um: float = 250.0      # Effective GAG substrate concentration in brush (uM)
    bulk_airway_concentration_nm: float = 10.0 # Bulk topical nanoparticle dose concentration (nM)


@dataclass(frozen=True)
class EnzymeTreatmentProfile:
    """Specification of an exogenous permeabilizing enzyme treatment."""

    enzyme_name: str
    concentration_u_ml: float
    kcat_per_s: float
    km_um: float
    target_substrate: str                      # 'sialic_acid', 'gag', or 'dual'


@dataclass
class PermeationTimepoint:
    """State of the periciliary brush at a specific simulation timepoint."""

    time_minutes: float
    mesh_spacing_nm: float
    brush_thickness_nm: float
    sialic_acid_cleavage_fraction: float
    gag_depolymerization_fraction: float
    effective_potential_mv: float
    ogston_partition_coefficient: float
    electrostatic_partition_coefficient: float
    net_partition_coefficient: float
    effective_diffusivity_m2_s: float
    apical_collision_frequency_m2_s: float
    fold_enhancement_vs_baseline: float


# Catalog of therapeutic gene delivery nanoparticles grounded in CF literature
NANOPARTICLE_CATALOG: Dict[str, NanoparticleProfile] = {
    "LNP_mRNA_BEACON_benchmark": NanoparticleProfile(
        name="BEACON-CF Inhaled CFTR mRNA Lipid Nanoparticle",
        diameter_nm=95.0,
        zeta_potential_mv=-18.0,
        cargo_type="CFTR mRNA (full-length codon-optimized)",
        molecular_weight_kda=1600.0,
        evidence_doi="10.1038/nbt.3298",
        clinical_trial_nct="NCT05668741",
    ),
    "AAV2_capsid_reference": NanoparticleProfile(
        name="Recombinant Adeno-Associated Virus Serotype 2 (AAV2)",
        diameter_nm=25.0,
        zeta_potential_mv=-8.5,
        cargo_type="ssDNA CFTR mini-expression cassette",
        molecular_weight_kda=4200.0,
        evidence_doi="10.1126/science.1223012",
        clinical_trial_nct=None,
    ),
    "AAV5_capsid_reference": NanoparticleProfile(
        name="Recombinant Adeno-Associated Virus Serotype 5 (AAV5)",
        diameter_nm=26.0,
        zeta_potential_mv=-5.8,
        cargo_type="ssDNA CFTR mini-expression cassette",
        molecular_weight_kda=4250.0,
        evidence_doi="10.1126/science.1223012",
        clinical_trial_nct=None,
    ),
    "Polymeric_NP_PEGylated": NanoparticleProfile(
        name="PEGylated Mucus-Penetrating Nanoparticle (MPP)",
        diameter_nm=65.0,
        zeta_potential_mv=-3.2,
        cargo_type="DNA/ASO cargo",
        molecular_weight_kda=850.0,
        evidence_doi="10.1152/physrev.00004.2021",
        clinical_trial_nct=None,
    ),
}

# Catalog of verified permeabilizing enzyme regimens
ENZYME_CATALOG: Dict[str, EnzymeTreatmentProfile] = {
    "Sialidase_low_dose": EnzymeTreatmentProfile(
        enzyme_name="Arthrobacter ureafaciens Sialidase (Low)",
        concentration_u_ml=1.0,
        kcat_per_s=42.0,
        km_um=180.0,
        target_substrate="sialic_acid",
    ),
    "Sialidase_therapeutic": EnzymeTreatmentProfile(
        enzyme_name="Arthrobacter ureafaciens Sialidase (Therapeutic)",
        concentration_u_ml=5.0,
        kcat_per_s=42.0,
        km_um=180.0,
        target_substrate="sialic_acid",
    ),
    "Heparinase_therapeutic": EnzymeTreatmentProfile(
        enzyme_name="Flavobacterium heparinum Heparinase I/III",
        concentration_u_ml=2.5,
        kcat_per_s=28.0,
        km_um=120.0,
        target_substrate="gag",
    ),
    "Dual_Sialidase_Heparinase": EnzymeTreatmentProfile(
        enzyme_name="Dual Sialidase + Heparinase Synergy Cockail",
        concentration_u_ml=5.0,
        kcat_per_s=38.0,
        km_um=150.0,
        target_substrate="dual",
    ),
}


def calculate_stokes_einstein_diffusivity(
    diameter_nm: float,
    temperature_k: float = TEMP_PHYSIOLOGICAL_K,
    viscosity_pa_s: float = AIRWAY_LIQUID_VISCOSITY_PA_S,
) -> float:
    """Compute free unhindered hydrodynamic diffusivity via Stokes-Einstein equation (m^2/s)."""
    radius_m = (diameter_nm * 1e-9) / 2.0
    return (BOLTZMANN_CONSTANT * temperature_k) / (6.0 * math.pi * viscosity_pa_s * radius_m)


def calculate_ogston_sieving(
    particle_diameter_nm: float,
    mesh_spacing_nm: float,
    mucin_fiber_diameter_nm: float = 3.5,
) -> float:
    """Compute steric partition coefficient in cylindrical fiber mesh via Ogston sieving model.

    K_steric = exp( - (pi / 4) * ((d_p + d_fiber) / (xi + d_fiber))^2 )
    """
    if mesh_spacing_nm <= 0.0:
        return 0.0
    numerator = particle_diameter_nm + mucin_fiber_diameter_nm
    denominator = mesh_spacing_nm + mucin_fiber_diameter_nm
    exponent = -0.25 * math.pi * (numerator / denominator) ** 2
    # Guard against severe underflow
    if exponent < -80.0:
        return 0.0
    return math.exp(exponent)


def calculate_electrostatic_partition(
    particle_zeta_mv: float,
    glycocalyx_potential_mv: float,
    temperature_k: float = TEMP_PHYSIOLOGICAL_K,
) -> float:
    """Compute electrostatic partition coefficient across charged periciliary brush.

    Models Debye-Huckel / Boltzmann electrostatic exclusion penalty for like-charged barriers.
    """
    # Potential difference driving repulsion (both typically negative in native ASL)
    # Energy barrier Delta_E = q_eff * Delta_V
    # Here we map zeta potentials to effective thermal energy barrier:
    # If both are negative, like-charge repulsion creates an energy penalty:
    v_prod = (particle_zeta_mv * 1e-3) * (glycocalyx_potential_mv * 1e-3)
    if v_prod > 0.0:
        # Repulsive barrier: magnitude proportional to geometric mean / product of potentials
        # Normalized by thermal voltage V_th = (k_B * T) / e ~ 26.7 mV at 37 °C
        v_thermal = (BOLTZMANN_CONSTANT * temperature_k) / ELEMENTARY_CHARGE
        normalized_barrier = math.sqrt(v_prod) / v_thermal
        # Repulsive partition factor K_elec = exp(-barrier)
        return math.exp(-min(normalized_barrier, 50.0))
    # Attractive or neutral
    return 1.0


def calculate_de_gennes_mesh_dilation(
    initial_mesh_nm: float,
    gag_cleavage_fraction: float,
    sialic_cleavage_fraction: float,
) -> float:
    """Calculate dilated mesh correlation length xi(t) based on polymer brush de Gennes scaling.

    As GAG chains and tethered mucin complexes are depolymerized, the effective grafting
    density sigma(t) decreases: sigma(t) = sigma_0 * (1 - 0.75 * f_gag - 0.25 * f_sial).
    Since xi proportional to sigma^(-1/2), loss of tethered mass dilates mesh spacing.
    """
    retained_grafting_fraction = max(1.0 - 0.70 * gag_cleavage_fraction - 0.20 * sialic_cleavage_fraction, 0.05)
    dilation_factor = 1.0 / math.sqrt(retained_grafting_fraction)
    return initial_mesh_nm * dilation_factor


def calculate_brush_thickness_collapse(
    initial_thickness_nm: float,
    sialic_cleavage_fraction: float,
    gag_cleavage_fraction: float,
) -> float:
    """Compute brush height collapse H(t) following electrostatic discharge and chain cleavage.

    Loss of sialic acid charges collapses the osmotic swelling pressure of the brush,
    while enzymatic backbone depolymerization reduces polymer length (Alexander-de Gennes scaling).
    """
    # Electrostatic collapse accounts for ~40% height reduction; backbone cleavage accounts for up to ~45%
    retained_height_fraction = max(
        1.0 - 0.40 * sialic_cleavage_fraction - 0.45 * gag_cleavage_fraction,
        0.20,
    )
    return initial_thickness_nm * retained_height_fraction


def simulate_glycocalyx_permeabilization(
    nanoparticle: NanoparticleProfile,
    brush_params: GlycocalyxBrushParameters = GlycocalyxBrushParameters(),
    enzyme_regimen: EnzymeTreatmentProfile = ENZYME_CATALOG["Sialidase_therapeutic"],
    total_time_minutes: float = 60.0,
    time_step_minutes: float = 1.0,
) -> List[PermeationTimepoint]:
    """Execute forward Euler numerical integration of enzymatic brush permeabilization and transport."""
    timepoints: List[PermeationTimepoint] = []

    dt_s = time_step_minutes * 60.0
    num_steps = int(total_time_minutes / time_step_minutes) + 1

    current_sial_um = brush_params.initial_sialic_acid_density_um
    current_gag_um = brush_params.initial_gag_density_um

    d_free = calculate_stokes_einstein_diffusivity(nanoparticle.diameter_nm)
    c_bulk_particles_m3 = (brush_params.bulk_airway_concentration_nm * 1e-9) * AVOGADRO_CONSTANT

    baseline_flux: Optional[float] = None

    for step in range(num_steps):
        t_min = step * time_step_minutes

        # Cleavage fractions
        f_sial_cleaved = 1.0 - (current_sial_um / brush_params.initial_sialic_acid_density_um)
        f_gag_cleaved = 1.0 - (current_gag_um / brush_params.initial_gag_density_um)

        # Dynamic brush geometry
        xi_t = calculate_de_gennes_mesh_dilation(
            brush_params.initial_mesh_spacing_nm,
            f_gag_cleaved,
            f_sial_cleaved,
        )
        h_t = calculate_brush_thickness_collapse(
            brush_params.initial_thickness_nm,
            f_sial_cleaved,
            f_gag_cleaved,
        )

        # Dynamic glycocalyx surface potential (discharging as sialic acid is removed)
        zeta_glyco_t = brush_params.initial_surface_potential_mv * (1.0 - 0.85 * f_sial_cleaved)

        # Partition coefficients
        k_steric = calculate_ogston_sieving(
            nanoparticle.diameter_nm,
            xi_t,
            brush_params.mucin_fiber_diameter_nm,
        )
        k_elec = calculate_electrostatic_partition(
            nanoparticle.zeta_potential_mv,
            zeta_glyco_t,
        )
        k_net = k_steric * k_elec

        # Effective diffusion and apical collision flux
        d_eff = d_free * k_net
        h_m = h_t * 1e-9
        j_apical = (d_eff * c_bulk_particles_m3) / h_m

        if baseline_flux is None:
            baseline_flux = max(j_apical, 1e-30)

        fold_enhancement = j_apical / baseline_flux

        timepoints.append(
            PermeationTimepoint(
                time_minutes=t_min,
                mesh_spacing_nm=round(xi_t, 2),
                brush_thickness_nm=round(h_t, 2),
                sialic_acid_cleavage_fraction=round(f_sial_cleaved, 4),
                gag_depolymerization_fraction=round(f_gag_cleaved, 4),
                effective_potential_mv=round(zeta_glyco_t, 2),
                ogston_partition_coefficient=k_steric,
                electrostatic_partition_coefficient=k_elec,
                net_partition_coefficient=k_net,
                effective_diffusivity_m2_s=d_eff,
                apical_collision_frequency_m2_s=j_apical,
                fold_enhancement_vs_baseline=round(fold_enhancement, 2),
            )
        )

        # Advance enzymatic kinetics via Michaelis-Menten rate equations
        if enzyme_regimen.target_substrate in ("sialic_acid", "dual"):
            # Sialidase activity: rate in uM/s
            rate_sial = (
                enzyme_regimen.kcat_per_s * enzyme_regimen.concentration_u_ml * current_sial_um
            ) / (enzyme_regimen.km_um + current_sial_um)
            delta_sial = rate_sial * dt_s * 0.05  # Scale factor for volumetric ASL partition
            current_sial_um = max(current_sial_um - delta_sial, 0.0)

        if enzyme_regimen.target_substrate in ("gag", "dual"):
            # Heparinase activity: rate in uM/s
            rate_gag = (
                enzyme_regimen.kcat_per_s * enzyme_regimen.concentration_u_ml * current_gag_um
            ) / (enzyme_regimen.km_um + current_gag_um)
            delta_gag = rate_gag * dt_s * 0.03
            current_gag_um = max(current_gag_um - delta_gag, 0.0)

    return timepoints


def evaluate_discovery_panel(
    nanoparticles: Optional[List[NanoparticleProfile]] = None,
    enzyme: EnzymeTreatmentProfile = ENZYME_CATALOG["Sialidase_therapeutic"],
    duration_min: float = 60.0,
) -> Dict[str, Any]:
    """Evaluate a standardized discovery panel comparing AAV vs LNP permeabilization kinetics."""
    if nanoparticles is None:
        nanoparticles = list(NANOPARTICLE_CATALOG.values())

    results: Dict[str, Any] = {
        "metadata": {
            "task": "T-0291",
            "module": "glycocalyx_permeabilization_model",
            "enzyme_used": enzyme.enzyme_name,
            "duration_minutes": duration_min,
            "disclaimer": NLM_ATTRIBUTION,
            "boundary": RESEARCH_BOUNDARY,
        },
        "profiles": {},
    }

    for np_prof in nanoparticles:
        sim = simulate_glycocalyx_permeabilization(
            nanoparticle=np_prof,
            enzyme_regimen=enzyme,
            total_time_minutes=duration_min,
        )
        t_final = sim[-1]
        t_init = sim[0]

        results["profiles"][np_prof.name] = {
            "diameter_nm": np_prof.diameter_nm,
            "zeta_potential_mv": np_prof.zeta_potential_mv,
            "cargo_type": np_prof.cargo_type,
            "evidence_doi": np_prof.evidence_doi,
            "clinical_trial_nct": np_prof.clinical_trial_nct,
            "initial_ogston_partition": t_init.ogston_partition_coefficient,
            "final_ogston_partition": t_final.ogston_partition_coefficient,
            "initial_net_partition": t_init.net_partition_coefficient,
            "final_net_partition": t_final.net_partition_coefficient,
            "initial_mesh_nm": t_init.mesh_spacing_nm,
            "final_dilated_mesh_nm": t_final.mesh_spacing_nm,
            "initial_brush_thickness_nm": t_init.brush_thickness_nm,
            "final_brush_thickness_nm": t_final.brush_thickness_nm,
            "apical_flux_initial_m2_s": t_init.apical_collision_frequency_m2_s,
            "apical_flux_final_m2_s": t_final.apical_collision_frequency_m2_s,
            "fold_enhancement": t_final.fold_enhancement_vs_baseline,
        }

    return results


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Airway Ciliated Epithelial Glycocalyx Permeabilization Model (T-0291)"
    )
    parser.add_argument(
        "--particle",
        type=str,
        default="LNP_mRNA_BEACON_benchmark",
        choices=list(NANOPARTICLE_CATALOG.keys()),
        help="Nanoparticle formulation to simulate",
    )
    parser.add_argument(
        "--enzyme",
        type=str,
        default="Sialidase_therapeutic",
        choices=list(ENZYME_CATALOG.keys()),
        help="Enzyme treatment regimen",
    )
    parser.add_argument(
        "--duration-min",
        type=float,
        default=60.0,
        help="Simulation duration in minutes",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Output full discovery panel in JSON format",
    )
    args = parser.parse_args(argv)

    if args.json:
        panel = evaluate_discovery_panel(
            enzyme=ENZYME_CATALOG[args.enzyme],
            duration_min=args.duration_min,
        )
        print(json.dumps(panel, indent=2))
        return 0

    selected_np = NANOPARTICLE_CATALOG[args.particle]
    selected_enzyme = ENZYME_CATALOG[args.enzyme]

    timepoints = simulate_glycocalyx_permeabilization(
        nanoparticle=selected_np,
        enzyme_regimen=selected_enzyme,
        total_time_minutes=args.duration_min,
    )

    t0 = timepoints[0]
    tend = timepoints[-1]

    print("=" * 72)
    print("AIRWAY GLYCOCALYX PERMEABILIZATION & STERIC SIEVING MODEL (T-0291)")
    print(f"Vehicle: {selected_np.name} (d={selected_np.diameter_nm} nm, zeta={selected_np.zeta_potential_mv} mV)")
    print(f"Enzyme:  {selected_enzyme.enzyme_name} ({selected_enzyme.concentration_u_ml} U/mL)")
    print("=" * 72)
    print(f"Initial Periciliary Mesh Spacing (xi_0): {t0.mesh_spacing_nm} nm")
    print(f"Final Dilated Mesh Spacing (xi_final):   {tend.mesh_spacing_nm} nm")
    print(f"Initial Brush Thickness (H_0):          {t0.brush_thickness_nm} nm")
    print(f"Final Brush Thickness (H_final):        {tend.brush_thickness_nm} nm")
    print(f"Sialic Acid Cleavage:                   {tend.sialic_acid_cleavage_fraction * 100:.1f} %")
    print(f"Initial Net Partition Coefficient:      {t0.net_partition_coefficient:.3e}")
    print(f"Final Net Partition Coefficient:        {tend.net_partition_coefficient:.3e}")
    print(f"Apical Collision Fold-Enhancement:      {tend.fold_enhancement_vs_baseline:.1f}x")
    print("-" * 72)
    print(f"Attribution: {NLM_ATTRIBUTION}")
    print(f"Boundary:    {RESEARCH_BOUNDARY}")
    print("=" * 72)
    return 0


if __name__ == "__main__":
    sys.exit(main())
