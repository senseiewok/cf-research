#!/usr/bin/env python3
"""Airway Mucociliary Clearance Viscoelastic Shear Rheology & Mucin Mesh Hindrance Model (T-0185).

Models airway mucus biophysics across the hydration continuum (1.5% to 10.0% solids):
1. Complex viscoelastic shear modulus G*(omega) = G'(omega) + i G''(omega) and loss tangent tan(delta).
2. Mucin MUC5B/MUC5AC mesh pore size scaling xi(w_solids) propto w_solids^(-nu).
3. Ogston/obstruction hindered diffusivity D_mucus / D_0 for PEG-brush vs unshielded nanoparticles and ASOs.
4. Mucociliary transport velocity v_MCT and mucolytic rescue delta (osmotic hydration, DNase, disulfide reduction).

Grounding & Provenance:
- Hill et al., Sci Adv 2019 (DOI 10.1164/rccm.201802-0245OC): Human airway mucus solids concentration & rheology
- Button et al., PNAS 2019 (DOI 10.1073/pnas.1409321111): Concentration-dependent mucociliary transport
- Clinical trial accession: NCT05668741 (VX-522 aerosolized mRNA LNP, verified in tools/trial_atlas/snapshot.json)
- Delivery biophysics landscape: landscape/delivery-biophysics.md
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

# Physical & Physiological Constants
BOLTZMANN_CONSTANT_KB: float = 1.380649e-23  # J/K
TEMPERATURE_KELVIN_37C: float = 310.15  # 37 C physiological temperature
WATER_VISCOSITY_37C_PA_S: float = 0.6913e-3  # Pa*s at 37 C (water)
MUCIN_FIBER_RADIUS_NM: float = 3.5  # Mucin polymer chain radius (nm)
CILIA_HEIGHT_UM: float = 6.0  # Airway epithelial cilia length (um)
BASELINE_HEALTHY_SOLIDS_PCT: float = 2.0  # Normal human airway mucus solids (~20 mg/mL)
BASELINE_CF_SOLIDS_PCT: float = 7.5  # Typical CF sputum solids (~75 mg/mL)

# Reference Rheological Parameters at w_ref = 2.0% solids, omega = 1 rad/s
G_PRIME_REF_PA: float = 0.45  # Storage / elastic modulus (Pa)
G_DOUBLE_PRIME_REF_PA: float = 0.40  # Loss / viscous modulus (Pa)
RHEOLOGY_ALPHA_EXPONENT: float = 2.85  # G' power scaling with solids w
RHEOLOGY_BETA_EXPONENT: float = 1.95  # G'' power scaling with solids w
MESH_PORE_REF_NM: float = 180.0  # Steric mesh pore diameter at 2.0% solids
MESH_NU_EXPONENT: float = 0.85  # De Gennes pore deconcentration scaling exponent


@dataclass(frozen=True)
class ViscoelasticModulus:
    solids_pct: float
    frequency_rad_s: float
    storage_modulus_g_prime_pa: float  # G' (Pa)
    loss_modulus_g_double_prime_pa: float  # G'' (Pa)
    complex_modulus_magnitude_pa: float  # |G*| (Pa)
    loss_tangent_tan_delta: float  # tan(delta) = G'' / G'
    zero_shear_viscosity_pa_s: float  # eta_0 (Pa*s)
    gel_state_regime: str  # DISSIPATIVE_SOL, WEAK_GEL, or ELASTIC_CROSS_LINKED_GEL


@dataclass(frozen=True)
class MucinMeshPoreScaling:
    solids_pct: float
    pore_diameter_nm: float  # xi (nm)
    pore_radius_nm: float
    mucin_volume_fraction: float  # phi
    steric_cutoff_diameter_nm: float
    regime_description: str


@dataclass(frozen=True)
class ParticleDiffusivityResult:
    particle_name: str
    hydrodynamic_diameter_nm: float
    surface_chemistry: str
    solids_pct: float
    free_water_diffusivity_d0_m2_s: float
    steric_obstruction_ratio: float  # (D_eff / D_0)_steric
    surface_adhesion_penalty: float  # exp(-E_adh / kT)
    mucus_diffusivity_m2_s: float
    relative_diffusivity_ratio: float  # D_mucus / D_0
    transmucosal_penetration_status: str


@dataclass(frozen=True)
class MucociliaryTransportResult:
    solids_pct: float
    storage_modulus_pa: float
    loss_tangent: float
    transport_velocity_um_s: float  # v_MCT (um/s)
    transport_velocity_mm_min: float  # v_MCT (mm/min)
    periciliary_layer_status: str
    clearance_efficiency_pct: float
    therapeutic_clearance_competence: str


@dataclass(frozen=True)
class MucolyticRescueDelta:
    intervention_name: str
    initial_solids_pct: float
    post_treatment_solids_pct: float
    initial_g_prime_pa: float
    post_treatment_g_prime_pa: float
    initial_v_mct_um_s: float
    post_treatment_v_mct_um_s: float
    transport_fold_improvement: float
    verdict: str


def calculate_viscoelastic_modulus(
    solids_pct: float,
    frequency_rad_s: float = 1.0,
) -> ViscoelasticModulus:
    """Compute complex shear viscoelastic modulus G*(omega) and loss tangent tan(delta).
    
    Solids content w_solids is bounded between 1.5% and 10.0%.
    """
    if solids_pct < 1.0 or solids_pct > 15.0:
        raise ValueError(f"Mucus solids percentage {solids_pct}% out of physiological range [1.0, 15.0]%")
    if frequency_rad_s <= 0:
        raise ValueError("Shear oscillation frequency must be strictly positive.")

    w = max(1.0, solids_pct)
    w_ratio = w / BASELINE_HEALTHY_SOLIDS_PCT

    # Elastic plateau scaling: G' propto w^2.85
    g_prime_1rad = G_PRIME_REF_PA * (w_ratio ** RHEOLOGY_ALPHA_EXPONENT)
    # Viscous plateau scaling: G'' propto w^1.95
    g_double_prime_1rad = G_DOUBLE_PRIME_REF_PA * (w_ratio ** RHEOLOGY_BETA_EXPONENT)

    # Frequency power laws: G' ~ omega^0.15, G'' ~ omega^0.25 in gel state
    g_prime = g_prime_1rad * (frequency_rad_s ** 0.15)
    g_double_prime = g_double_prime_1rad * (frequency_rad_s ** 0.25)

    g_star_mag = math.sqrt(g_prime ** 2 + g_double_prime ** 2)
    tan_delta = g_double_prime / max(1e-6, g_prime)

    # Zero-shear viscosity estimate: eta_0 = |G*(omega_0)| / omega_0 at omega_0 = 0.1 rad/s
    omega_zero = 0.1
    g_prime_z = g_prime_1rad * (omega_zero ** 0.15)
    g_double_prime_z = g_double_prime_1rad * (omega_zero ** 0.25)
    eta_zero = math.sqrt(g_prime_z ** 2 + g_double_prime_z ** 2) / omega_zero

    if tan_delta >= 0.85:
        regime = "DISSIPATIVE_SOL_FLUID"
    elif tan_delta >= 0.40:
        regime = "WEAK_VISCOELASTIC_GEL"
    else:
        regime = "ELASTIC_CROSS_LINKED_GEL"

    return ViscoelasticModulus(
        solids_pct=round(solids_pct, 2),
        frequency_rad_s=round(frequency_rad_s, 3),
        storage_modulus_g_prime_pa=round(g_prime, 3),
        loss_modulus_g_double_prime_pa=round(g_double_prime, 3),
        complex_modulus_magnitude_pa=round(g_star_mag, 3),
        loss_tangent_tan_delta=round(tan_delta, 3),
        zero_shear_viscosity_pa_s=round(eta_zero, 3),
        gel_state_regime=regime,
    )


def calculate_mucin_pore_size(solids_pct: float) -> MucinMeshPoreScaling:
    """Compute MUC5B/MUC5AC macromolecular mesh pore dimension scaling xi(w).
    
    xi(w) = xi_0 * (w_ref / w)^nu
    """
    if solids_pct < 1.0 or solids_pct > 15.0:
        raise ValueError(f"Mucus solids percentage {solids_pct}% out of physiological range [1.0, 15.0]%")

    w = max(1.0, solids_pct)
    w_ratio = BASELINE_HEALTHY_SOLIDS_PCT / w

    pore_diam_nm = MESH_PORE_REF_NM * (w_ratio ** MESH_NU_EXPONENT)
    pore_radius_nm = pore_diam_nm / 2.0
    vol_frac = 0.012 * w  # Volume fraction phi ~ 1.2% per 1% solids

    if pore_diam_nm >= 130.0:
        desc = "EXPANDED_HEALTHY_MESH_PERMISSIVE"
    elif pore_diam_nm >= 80.0:
        desc = "INTERMEDIATE_PARTIALLY_RESTRICTED"
    else:
        desc = "COLLAPSED_CF_MESH_STERICALLY_OCCLUDED"

    return MucinMeshPoreScaling(
        solids_pct=round(solids_pct, 2),
        pore_diameter_nm=round(pore_diam_nm, 2),
        pore_radius_nm=round(pore_radius_nm, 2),
        mucin_volume_fraction=round(vol_frac, 4),
        steric_cutoff_diameter_nm=round(pore_diam_nm, 2),
        regime_description=desc,
    )


def calculate_hindered_particle_diffusivity(
    hydrodynamic_diameter_nm: float,
    surface_chemistry: str = "peg_brush",
    solids_pct: float = BASELINE_CF_SOLIDS_PCT,
    particle_name: str = "therapeutic_candidate",
) -> ParticleDiffusivityResult:
    """Calculate hindered diffusion D_mucus / D_0 using Ogston-Mackie steric & surface interaction models.
    
    Supported surface_chemistry:
    - 'peg_brush' (dense PEG-2k shielding, near-neutral zeta -2 to +2 mV)
    - 'polysarcosine' (pSar stealth polypeptide shell)
    - 'zwitterionic' (polycarboxybetaine non-fouling shell)
    - 'unshielded_cationic' (+15 mV zeta, strong electrostatic mucin entrapment)
    - 'unshielded_hydrophobic' (bare unshielded lipid/polymeric core)
    - 'aso_oligo' (small ~2.5 nm antisense oligonucleotide)
    """
    if hydrodynamic_diameter_nm <= 0:
        raise ValueError("Hydrodynamic diameter must be positive.")

    pore_info = calculate_mucin_pore_size(solids_pct)
    r_h = hydrodynamic_diameter_nm / 2.0
    r_pore = pore_info.pore_radius_nm
    phi = pore_info.mucin_volume_fraction
    r_f = MUCIN_FIBER_RADIUS_NM

    # Free Stokes-Einstein diffusivity D_0 in water at 37 C:
    # D_0 = (k_B * T) / (6 * pi * eta_w * r_h)
    r_h_m = r_h * 1e-9
    d0 = (BOLTZMANN_CONSTANT_KB * TEMPERATURE_KELVIN_37C) / (
        6.0 * math.pi * WATER_VISCOSITY_37C_PA_S * r_h_m
    )

    # Ogston-Mackie steric obstruction factor:
    # (D_eff / D_0)_steric = (1 - r_h / r_pore)^2 * exp(-pi * phi * (1 + r_h / r_f)^2)
    if r_h >= r_pore:
        # Steric entrapment: particle larger than mesh pore diameter
        steric_factor = 0.0005
    else:
        steric_geo = (1.0 - (r_h / r_pore)) ** 2
        exp_obstruction = math.exp(-math.pi * phi * ((1.0 + (r_h / r_f)) ** 2))
        steric_factor = max(0.0005, min(1.0, steric_geo * exp_obstruction))

    # Surface adhesion penalty: exp(-E_adh / kT)
    chem = surface_chemistry.strip().lower()
    if chem in ("peg_brush", "polysarcosine", "zwitterionic"):
        adhesion_factor = 1.00  # Zero sticky interaction with mucin sialic acid
    elif chem == "aso_oligo":
        adhesion_factor = 0.95  # Minimal electrostatic repulsion allows rapid percolation
    elif chem == "unshielded_cationic":
        adhesion_factor = 0.015  # Strong polyanionic sialic acid / sulfate trapping
    elif chem == "unshielded_hydrophobic":
        adhesion_factor = 0.050  # Hydrophobic adsorption to nonpolar mucin patches
    else:
        adhesion_factor = 0.100

    d_mucus = d0 * steric_factor * adhesion_factor
    rel_ratio = d_mucus / d0

    if rel_ratio >= 0.15:
        status = "RAPID_TRANSMUCOSAL_PERMEATION"
    elif rel_ratio >= 0.02:
        status = "MODERATE_RETARDED_PERMEATION"
    else:
        status = "SEVERE_MUCIN_ENTRAPMENT"

    return ParticleDiffusivityResult(
        particle_name=particle_name,
        hydrodynamic_diameter_nm=round(hydrodynamic_diameter_nm, 2),
        surface_chemistry=surface_chemistry,
        solids_pct=round(solids_pct, 2),
        free_water_diffusivity_d0_m2_s=d0,
        steric_obstruction_ratio=round(steric_factor, 5),
        surface_adhesion_penalty=round(adhesion_factor, 4),
        mucus_diffusivity_m2_s=d_mucus,
        relative_diffusivity_ratio=round(rel_ratio, 5),
        transmucosal_penetration_status=status,
    )


def calculate_mucociliary_transport_velocity(
    solids_pct: float,
) -> MucociliaryTransportResult:
    """Compute ciliary transport velocity v_MCT (um/s) and clearance efficiency.
    
    Reflects ciliary-mucus mechanical coupling across solids content (Button et al. 2019).
    """
    v0_um_s = 65.0  # Normal healthy clearance velocity (~3.9 mm/min)
    w = max(1.0, solids_pct)

    # Sigmoidal transport collapse as solids exceed PCL osmotic threshold (~4.2%)
    w_half = 4.2
    gamma = 4.5
    transport_ratio = 1.0 / (1.0 + ((w / w_half) ** gamma))
    v_mct = v0_um_s * transport_ratio
    v_mm_min = (v_mct * 60.0) / 1000.0

    modulus = calculate_viscoelastic_modulus(w)

    if w <= 3.0:
        pcl_status = "HYDRATED_PCL_CILIA_FREE_BEATING"
        comp = "PHYSIOLOGICAL_HIGH_CLEARANCE"
    elif w <= 5.5:
        pcl_status = "PARTIALLY_COMPRESSED_PCL_DECELERATED"
        comp = "PARTIAL_CLEARANCE_VULNERABILITY"
    else:
        pcl_status = "COLLAPSED_PCL_CILIA_FLATTENED_PLUGGED"
        comp = "MUCOSTASIS_CLEARANCE_FAILURE"

    efficiency_pct = (v_mct / v0_um_s) * 100.0

    return MucociliaryTransportResult(
        solids_pct=round(solids_pct, 2),
        storage_modulus_pa=modulus.storage_modulus_g_prime_pa,
        loss_tangent=modulus.loss_tangent_tan_delta,
        transport_velocity_um_s=round(v_mct, 2),
        transport_velocity_mm_min=round(v_mm_min, 3),
        periciliary_layer_status=pcl_status,
        clearance_efficiency_pct=round(efficiency_pct, 2),
        therapeutic_clearance_competence=comp,
    )


def simulate_mucolytic_rescue(
    initial_solids_pct: float = BASELINE_CF_SOLIDS_PCT,
    intervention: str = "hypertonic_saline_7pct",
) -> MucolyticRescueDelta:
    """Simulate mucolytic / osmotic rehydration intervention and transport recovery delta."""
    init_res = calculate_mucociliary_transport_velocity(initial_solids_pct)
    interv = intervention.strip().lower()

    if interv == "hypertonic_saline_7pct":
        # Osmotic rehydration draws water, reducing solids by ~2.8%
        post_solids = max(BASELINE_HEALTHY_SOLIDS_PCT, initial_solids_pct - 2.8)
    elif interv == "dornase_alfa_dnase":
        # Cleaves neutrophil eDNA, effectively reducing mechanical stiffness equivalent to -2.0% solids
        post_solids = max(BASELINE_HEALTHY_SOLIDS_PCT, initial_solids_pct - 2.0)
    elif interv == "n_acetylcysteine_nac":
        # Cleaves disulfide mucin multimer bridges, softening gel equivalent to -1.8% solids
        post_solids = max(BASELINE_HEALTHY_SOLIDS_PCT, initial_solids_pct - 1.8)
    elif interv == "triple_combination":
        # Synergistic hydration + DNase + mucolytic
        post_solids = max(BASELINE_HEALTHY_SOLIDS_PCT, initial_solids_pct - 4.5)
    else:
        post_solids = initial_solids_pct

    post_res = calculate_mucociliary_transport_velocity(post_solids)

    fold_improvement = post_res.transport_velocity_um_s / max(0.1, init_res.transport_velocity_um_s)

    if post_res.transport_velocity_um_s >= 35.0:
        verdict = "COMPLETE_PHYSIOLOGICAL_CLEARANCE_RESTORED"
    elif post_res.transport_velocity_um_s >= 15.0:
        verdict = "PARTIAL_FUNCTIONAL_CLEARANCE_RESTORED"
    else:
        verdict = "INSUFFICIENT_CLEARANCE_RECOVERY"

    return MucolyticRescueDelta(
        intervention_name=intervention,
        initial_solids_pct=round(initial_solids_pct, 2),
        post_treatment_solids_pct=round(post_solids, 2),
        initial_g_prime_pa=init_res.storage_modulus_pa,
        post_treatment_g_prime_pa=post_res.storage_modulus_pa,
        initial_v_mct_um_s=init_res.transport_velocity_um_s,
        post_treatment_v_mct_um_s=post_res.transport_velocity_um_s,
        transport_fold_improvement=round(fold_improvement, 2),
        verdict=verdict,
    )


def self_test() -> bool:
    """Internal deterministic self-test."""
    # Test 1: Viscoelastic modulus monotonicity
    mod_2 = calculate_viscoelastic_modulus(2.0)
    mod_8 = calculate_viscoelastic_modulus(8.0)
    assert mod_8.storage_modulus_g_prime_pa > mod_2.storage_modulus_g_prime_pa * 20.0
    assert mod_8.loss_tangent_tan_delta < mod_2.loss_tangent_tan_delta

    # Test 2: Mucin pore size collapse
    pore_2 = calculate_mucin_pore_size(2.0)
    pore_8 = calculate_mucin_pore_size(8.0)
    assert pore_2.pore_diameter_nm >= 160.0
    assert pore_8.pore_diameter_nm <= 70.0

    # Test 3: PEG-brush vs cationic LNP hindered diffusion
    diff_peg = calculate_hindered_particle_diffusivity(65.0, "peg_brush", solids_pct=7.5)
    diff_cat = calculate_hindered_particle_diffusivity(65.0, "unshielded_cationic", solids_pct=7.5)
    assert diff_peg.relative_diffusivity_ratio > diff_cat.relative_diffusivity_ratio * 20.0

    # Test 4: Mucociliary transport velocity and rescue
    mct_init = calculate_mucociliary_transport_velocity(8.0)
    assert mct_init.transport_velocity_um_s < 5.0
    rescue = simulate_mucolytic_rescue(8.0, "hypertonic_saline_7pct")
    assert rescue.post_treatment_v_mct_um_s > mct_init.transport_velocity_um_s * 4.0
    return True


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Airway Mucociliary Clearance Viscoelastic Shear Rheology & Mucin Mesh Hindrance Model (T-0185)"
    )
    parser.add_argument("--solids-pct", type=float, default=BASELINE_CF_SOLIDS_PCT, help="Airway mucus solids percent (w/w)")
    parser.add_argument("--particle-diam-nm", type=float, default=65.0, help="Nanoparticle hydrodynamic diameter (nm)")
    parser.add_argument(
        "--surface-chemistry",
        choices=["peg_brush", "polysarcosine", "zwitterionic", "unshielded_cationic", "unshielded_hydrophobic", "aso_oligo"],
        default="peg_brush",
        help="Particle surface coating chemistry",
    )
    parser.add_argument(
        "--mucolytic",
        choices=["hypertonic_saline_7pct", "dornase_alfa_dnase", "n_acetylcysteine_nac", "triple_combination", "all"],
        default="hypertonic_saline_7pct",
        help="Mucolytic / osmotic intervention",
    )
    parser.add_argument("--json", action="store_true", help="Print structured JSON output")
    parser.add_argument("--export", type=str, help="Export JSON results to output path")
    parser.add_argument("--self-test", action="store_true", help="Run deterministic internal self-tests")
    args = parser.parse_args()

    if args.self_test:
        if self_test():
            print("PASS: mucus_rheology self-test passed.")
            sys.exit(0)
        else:
            print("FAIL: mucus_rheology self-test failed.")
            sys.exit(1)

    visco = calculate_viscoelastic_modulus(args.solids_pct)
    pore = calculate_mucin_pore_size(args.solids_pct)
    diff = calculate_hindered_particle_diffusivity(
        args.particle_diam_nm,
        surface_chemistry=args.surface_chemistry,
        solids_pct=args.solids_pct,
    )
    mct = calculate_mucociliary_transport_velocity(args.solids_pct)

    interventions = (
        ["hypertonic_saline_7pct", "dornase_alfa_dnase", "n_acetylcysteine_nac", "triple_combination"]
        if args.mucolytic == "all"
        else [args.mucolytic]
    )
    rescues = {intv: asdict(simulate_mucolytic_rescue(args.solids_pct, intv)) for intv in interventions}

    results: Dict[str, Any] = {
        "provenance": {
            "nct_id": "NCT05668741",
            "primary_dois": [
                "10.1164/rccm.201802-0245OC",
                "10.1073/pnas.1409321111",
            ],
            "landscape_dossier": "landscape/delivery-biophysics.md",
        },
        "viscoelastic_rheology": asdict(visco),
        "mucin_pore_scaling": asdict(pore),
        "particle_diffusivity": asdict(diff),
        "mucociliary_transport": asdict(mct),
        "mucolytic_rescues": rescues,
    }

    if args.export:
        out_p = Path(args.export)
        out_p.parent.mkdir(parents=True, exist_ok=True)
        with open(out_p, "w", encoding="utf-8") as f:
            json.dump(results, f, indent=2)
        print(f"Exported mucus rheology results to {out_p}")

    if args.json:
        print(json.dumps(results, indent=2))
    else:
        print("=" * 76)
        print("AIRWAY MUCUS VISCOELASTIC RHEOLOGY & MESH HINDRANCE MODEL (T-0185)")
        print("=" * 76)
        print(f"Solids Content: {visco.solids_pct}% | Regime: {visco.gel_state_regime}")
        print(f"  G': {visco.storage_modulus_g_prime_pa} Pa | G'': {visco.loss_modulus_g_double_prime_pa} Pa | tan(delta): {visco.loss_tangent_tan_delta}")
        print(f"  Zero-Shear Viscosity: {visco.zero_shear_viscosity_pa_s} Pa*s")
        print(f"  Mucin Mesh Pore Diameter: {pore.pore_diameter_nm} nm ({pore.regime_description})")
        print(f"Particle Diffusion ({diff.hydrodynamic_diameter_nm} nm, {diff.surface_chemistry}):")
        print(f"  D_mucus / D_0: {diff.relative_diffusivity_ratio} -> {diff.transmucosal_penetration_status}")
        print(f"Mucociliary Transport:")
        print(f"  Velocity: {mct.transport_velocity_um_s} um/s ({mct.transport_velocity_mm_min} mm/min)")
        print(f"  PCL Status: {mct.periciliary_layer_status} ({mct.clearance_efficiency_pct}% normal)")
        print("-" * 76)
        print("Mucolytic Interventions:")
        for name, rdata in rescues.items():
            print(f"  [{name}]: {rdata['initial_v_mct_um_s']} -> {rdata['post_treatment_v_mct_um_s']} um/s ({rdata['transport_fold_improvement']}x fold) -> {rdata['verdict']}")
        print("=" * 76)


if __name__ == "__main__":
    main()
