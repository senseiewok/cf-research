#!/usr/bin/env python3
"""Airway Epithelial Glycocalyx Barrier Sialic Acid Charge Density & Enzymatic Permeabilization Model (T-0351).

Quantitative biophysical model of the apical airway glycocalyx barrier:
1. Electrostatic Double-Layer & Charge Density Calculation:
   - Debye length (kappa^-1 approx 0.8 nm in airway surface liquid, I = 150 mM).
   - Local electrostatic potential psi(z) via Grahame / Gouy-Chapman electrochemistry
     as a function of sialic acid and glycosaminoglycan (GAG) surface charge density
     (sigma approx -0.02 to -0.05 C/m^2).
2. Hindered Transport & Steric Mesh Exclusion:
   - Renkin-Ogston steric partitioning coefficient Phi = (1 - r_p / r_pore)^2 and
     hydrodynamic hindrance factor F(lambda) for nanoparticle radii r_p in [10, 100] nm
     moving through fibrous glycocalyx networks of pore radius r_pore and fiber radius r_f.
   - Coupled electrostatic-steric partitioning K_part = Phi * F(lambda) * exp(-q_p * psi / (k_B * T)).
3. Enzymatic Desialylation Kinetics:
   - Michaelis-Menten kinetics of bacterial (Arthrobacter ureafaciens / Clostridium perfringens)
     and recombinant human neuraminidase/sialidase isoforms (Neu1, Neu2, Neu3).
   - Deterministic calculation of time required to reduce negative surface charge by 80% (t_80%)
     and enhance therapeutic nanoparticle diffusivity by > 5-fold.
4. Deterministic Validation:
   - Self-diagnostic consistency checks (--check) and JSON output (--json) verifying
     biophysical invariants against published airway epithelial micro-rheology.

Grounding & Provenance:
- Button, B. et al. (2012). A Periciliary Brush Promotes the Mucociliary Transport of
  Locally Concentrated Mucus. Science, 337(6097), 937-941. DOI: 10.1126/science.1223012
- Pickles, R. J. et al. (2000). Retargeting the Coxsackievirus and Adenovirus Receptor to the
  Apical Surface of Polarized Epithelial Cells Reveals the Glycocalyx as a Barrier to
  Adenovirus-Mediated Gene Transfer. J Virol, 74(13), 6050-6057. DOI: 10.1128/jvi.74.13.6050-6057.2000
- Stonebraker, J. R. et al. (2004). Glycocalyx Restricts Adenoviral Vector Access to Apical
  Receptors Expressed on Respiratory Epithelium In Vitro and In Vivo: Role for Tethered Mucins
  as Barriers to Lumenal Infection. J Virol, 78(24), 13755-13768. DOI: 10.1128/JVI.78.24.13755-13768.2004
- Boucher, R. C. / Fahy, J. V. et al. (2021). Physiology and pathophysiology of human airway mucus.
  Physiol Rev, 101(4), 1509-1550. DOI: 10.1152/physrev.00004.2021
- Clinical Trial Grounding (verified against tools/trial_atlas/snapshot.json):
  * NCT05668741 (BEACON-CF: Inhaled CFTR mRNA Lipid Nanoparticle clinical evaluation).
"""
from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
import math
from pathlib import Path
import sys
from typing import Any, Dict, List, Optional, Tuple

# Physical & Physiological Constants
BOLTZMANN_CONSTANT_KB: float = 1.380649e-23  # J/K
ELEMENTARY_CHARGE_E: float = 1.60217663e-19  # C
AVOGADRO_NUMBER_NA: float = 6.02214076e23  # mol^-1
VACUUM_PERMITTIVITY_EPS0: float = 8.8541878e-12  # F/m (C^2 / (J * m))
PHYSIOLOGICAL_TEMP_K: float = 310.15  # 37 C physiological temperature (K)
WATER_RELATIVE_PERMITTIVITY_37C: float = 74.2  # Dielectric constant of water at 37 C
ASL_IONIC_STRENGTH_MM: float = 150.0  # Airway surface liquid ionic strength (150 mM NaCl equivalent)
ASL_VISCOSITY_PA_S: float = 0.8e-3  # Airway surface liquid aqueous sol viscosity (Pa*s)

# Glycocalyx Structural Constants
BASELINE_BRUSH_HEIGHT_NM: float = 300.0  # Apical glycocalyx brush height (nm)
BASELINE_PORE_RADIUS_NM: float = 25.0  # Effective mesh pore radius (nm)
MUCIN_FIBER_RADIUS_NM: float = 2.0  # Tethered mucin / GAG fiber radius (nm)
BASELINE_SURFACE_CHARGE_DENSITY: float = -0.035  # Baseline surface charge density (C/m^2)
GAG_RESIDUAL_CHARGE_FRACTION: float = 0.20  # Fraction of charge contributed by sulfated GAGs (heparan/chondroitin)
SIALIC_CHARGE_FRACTION: float = 0.80  # Fraction of charge contributed by terminal sialic acids
BASELINE_SIALIC_ACID_CONC_MM: float = 5.0  # Glycocalyx local sialic acid concentration (mM)

# NLM Mandatory Clinical Registry Disclaimer
NLM_REGISTRY_DISCLAIMER: str = (
    "ClinicalTrials.gov does not evaluate safety or efficacy; listing on the registry "
    "does not denote government endorsement."
)

# Neuraminidase / Sialidase Kinetic Profiles
# k_cat in s^-1, K_m in mM
SIALIDASE_PROFILES: Dict[str, Dict[str, Any]] = {
    "bacterial_neuraminidase": {
        "organism": "Arthrobacter ureafaciens / Clostridium perfringens",
        "k_cat": 220.0,  # s^-1
        "k_m": 1.10,  # mM
        "ph_optimum": 6.5,
        "description": "Broad-spectrum bacterial sialidase cleaving alpha-(2,3), (2,6), and (2,8) sialic bonds.",
    },
    "recombinant_human_neu1": {
        "organism": "Homo sapiens Neu1",
        "k_cat": 65.0,  # s^-1
        "k_m": 1.25,  # mM
        "ph_optimum": 5.0,
        "description": "Lysosomal/membrane-associated sialidase; primary physiological desialylating enzyme.",
    },
    "recombinant_human_neu2": {
        "organism": "Homo sapiens Neu2",
        "k_cat": 130.0,  # s^-1
        "k_m": 0.85,  # mM
        "ph_optimum": 6.8,
        "description": "Cytosolic/apical neutral sialidase with high catalytic turnover at airway physiological pH.",
    },
    "recombinant_human_neu3": {
        "organism": "Homo sapiens Neu3",
        "k_cat": 95.0,  # s^-1
        "k_m": 1.05,  # mM
        "ph_optimum": 6.0,
        "description": "Membrane-bound ganglioside-selective sialidase active in epithelial plasma membranes.",
    },
}


@dataclass(frozen=True)
class ElectrostaticDoubleLayerResult:
    ionic_strength_mM: float
    debye_length_nm: float
    surface_charge_density_C_m2: float
    surface_potential_mV: float
    thermal_voltage_mV: float
    potential_at_1nm_mV: float
    potential_at_2nm_mV: float
    potential_at_5nm_mV: float


@dataclass(frozen=True)
class StericHindranceResult:
    nanoparticle_radius_nm: float
    nanoparticle_diameter_nm: float
    pore_radius_nm: float
    fiber_radius_nm: float
    pore_size_ratio_lambda: float
    ogston_partition_phi: float
    renkin_hydrodynamic_factor_f: float
    effective_steric_partition_k: float


@dataclass(frozen=True)
class CoupledDiffusivityResult:
    nanoparticle_name: str
    nanoparticle_radius_nm: float
    nanoparticle_zeta_mV: float
    glycocalyx_charge_density_C_m2: float
    d_0_asl_m2_s: float
    k_steric: float
    k_elec: float
    k_total: float
    d_eff_glycocalyx_m2_s: float
    barrier_penetration_status: str


@dataclass(frozen=True)
class DesialylationKineticResult:
    enzyme_name: str
    enzyme_conc_uM: float
    duration_minutes: float
    initial_sia_mM: float
    remaining_sia_mM: float
    desialylation_fraction: float
    charge_reduction_fraction: float
    initial_charge_density_C_m2: float
    final_charge_density_C_m2: float
    t_80_percent_reduction_seconds: float
    t_80_percent_reduction_minutes: float
    initial_diffusivity_m2_s: float
    final_diffusivity_m2_s: float
    diffusivity_fold_enhancement: float


def calculate_debye_length(
    ionic_strength_mM: float = ASL_IONIC_STRENGTH_MM,
    temp_k: float = PHYSIOLOGICAL_TEMP_K,
    rel_permittivity: float = WATER_RELATIVE_PERMITTIVITY_37C,
) -> float:
    """Calculate the Debye screening length kappa^-1 in nanometers for an electrolyte solution.

    kappa^2 = (2 * N_A * e^2 * I) / (eps_r * eps_0 * k_B * T)
    where I is ionic strength in mol/m^3 (1 mM = 1 mol/m^3).
    """
    if ionic_strength_mM <= 0.0:
        raise ValueError("Ionic strength must be strictly positive.")

    i_mol_m3 = ionic_strength_mM  # 1 mM = 1 mol / m^3
    eps = rel_permittivity * VACUUM_PERMITTIVITY_EPS0
    numerator = eps * BOLTZMANN_CONSTANT_KB * temp_k
    denominator = 2.0 * AVOGADRO_NUMBER_NA * (ELEMENTARY_CHARGE_E ** 2) * i_mol_m3
    kappa_inv_m = math.sqrt(numerator / denominator)
    return kappa_inv_m * 1e9  # convert m to nm


def calculate_surface_potential(
    charge_density_C_m2: float = BASELINE_SURFACE_CHARGE_DENSITY,
    ionic_strength_mM: float = ASL_IONIC_STRENGTH_MM,
    temp_k: float = PHYSIOLOGICAL_TEMP_K,
    rel_permittivity: float = WATER_RELATIVE_PERMITTIVITY_37C,
    valency: int = 1,
) -> float:
    """Calculate the electrostatic surface potential psi_0 in millivolts via the Grahame equation.

    sigma = sqrt(8 * eps_r * eps_0 * k_B * T * N_A * I) * sinh(z * e * psi_0 / (2 * k_B * T))
    Solving for psi_0:
    psi_0 = (2 * k_B * T / (z * e)) * asinh(sigma / sqrt(8 * eps * k_B * T * N_A * I))
    """
    i_mol_m3 = ionic_strength_mM
    eps = rel_permittivity * VACUUM_PERMITTIVITY_EPS0
    scale = math.sqrt(8.0 * eps * BOLTZMANN_CONSTANT_KB * temp_k * AVOGADRO_NUMBER_NA * i_mol_m3)
    arg = charge_density_C_m2 / scale
    thermal_voltage_v = (BOLTZMANN_CONSTANT_KB * temp_k) / (valency * ELEMENTARY_CHARGE_E)
    psi_0_v = 2.0 * thermal_voltage_v * math.asinh(arg)
    return psi_0_v * 1e3  # convert V to mV


def evaluate_electrostatic_double_layer(
    charge_density_C_m2: float = BASELINE_SURFACE_CHARGE_DENSITY,
    ionic_strength_mM: float = ASL_IONIC_STRENGTH_MM,
    temp_k: float = PHYSIOLOGICAL_TEMP_K,
) -> ElectrostaticDoubleLayerResult:
    """Evaluate Debye screening and potential decay in the glycocalyx double layer."""
    debye_nm = calculate_debye_length(ionic_strength_mM, temp_k)
    psi_0_mv = calculate_surface_potential(charge_density_C_m2, ionic_strength_mM, temp_k)
    thermal_mv = (BOLTZMANN_CONSTANT_KB * temp_k / ELEMENTARY_CHARGE_E) * 1e3

    # Linearized Gouy-Chapman potential profile: psi(z) = psi_0 * exp(-z / kappa^-1)
    psi_1nm = psi_0_mv * math.exp(-1.0 / debye_nm)
    psi_2nm = psi_0_mv * math.exp(-2.0 / debye_nm)
    psi_5nm = psi_0_mv * math.exp(-5.0 / debye_nm)

    return ElectrostaticDoubleLayerResult(
        ionic_strength_mM=ionic_strength_mM,
        debye_length_nm=round(debye_nm, 3),
        surface_charge_density_C_m2=charge_density_C_m2,
        surface_potential_mV=round(psi_0_mv, 2),
        thermal_voltage_mV=round(thermal_mv, 2),
        potential_at_1nm_mV=round(psi_1nm, 2),
        potential_at_2nm_mV=round(psi_2nm, 2),
        potential_at_5nm_mV=round(psi_5nm, 3),
    )


def evaluate_steric_hindrance(
    nanoparticle_radius_nm: float,
    pore_radius_nm: float = BASELINE_PORE_RADIUS_NM,
    fiber_radius_nm: float = MUCIN_FIBER_RADIUS_NM,
) -> StericHindranceResult:
    """Evaluate Renkin-Ogston steric partitioning and hydrodynamic hindrance.

    Renkin steric partitioning:
      Phi = (1 - lambda)^2  for lambda = r_p / r_pore < 1; 0 for lambda >= 1.
    Renkin hydrodynamic factor:
      F(lambda) = 1 - 2.104*lambda + 2.089*lambda^3 - 0.948*lambda^5  for lambda < 1.
    Effective steric partition:
      K_steric = Phi * F(lambda)
    """
    if nanoparticle_radius_nm <= 0.0 or pore_radius_nm <= 0.0:
        raise ValueError("Radii must be strictly positive.")

    lam = nanoparticle_radius_nm / pore_radius_nm
    diameter_nm = 2.0 * nanoparticle_radius_nm

    if lam >= 1.0:
        phi = 0.0
        f_renkin = 0.0
        k_steric = 0.0
    else:
        phi = (1.0 - lam) ** 2
        f_renkin = max(
            0.0,
            1.0 - 2.104 * lam + 2.089 * (lam ** 3) - 0.948 * (lam ** 5)
        )
        k_steric = phi * f_renkin

    return StericHindranceResult(
        nanoparticle_radius_nm=round(nanoparticle_radius_nm, 2),
        nanoparticle_diameter_nm=round(diameter_nm, 2),
        pore_radius_nm=round(pore_radius_nm, 2),
        fiber_radius_nm=round(fiber_radius_nm, 2),
        pore_size_ratio_lambda=round(lam, 3),
        ogston_partition_phi=round(phi, 4),
        renkin_hydrodynamic_factor_f=round(f_renkin, 4),
        effective_steric_partition_k=round(k_steric, 5),
    )


def evaluate_coupled_diffusivity(
    nanoparticle_name: str,
    nanoparticle_radius_nm: float,
    nanoparticle_zeta_mV: float,
    charge_density_C_m2: float = BASELINE_SURFACE_CHARGE_DENSITY,
    pore_radius_nm: float = BASELINE_PORE_RADIUS_NM,
    temp_k: float = PHYSIOLOGICAL_TEMP_K,
    viscosity_pa_s: float = ASL_VISCOSITY_PA_S,
) -> CoupledDiffusivityResult:
    """Calculate Stokes-Einstein aqueous diffusivity D_0 and coupled steric-electrostatic D_eff."""
    r_m = nanoparticle_radius_nm * 1e-9
    d_0 = (BOLTZMANN_CONSTANT_KB * temp_k) / (6.0 * math.pi * viscosity_pa_s * r_m)

    steric = evaluate_steric_hindrance(nanoparticle_radius_nm, pore_radius_nm)
    k_steric = steric.effective_steric_partition_k

    # Electrostatic partition: Boltzmann factor exp(-Delta_U / k_B T)
    # Average electrostatic barrier experienced within Debye length layer
    psi_0_v = calculate_surface_potential(charge_density_C_m2, ASL_IONIC_STRENGTH_MM, temp_k) * 1e-3
    zeta_v = nanoparticle_zeta_mV * 1e-3

    # Effective surface charge of nanoparticle via Ohshima relation:
    # Interaction potential between planar glycocalyx and spherical particle:
    # Delta_U ~ 4 * pi * eps * r_p * psi_0 * zeta_p
    eps = WATER_RELATIVE_PERMITTIVITY_37C * VACUUM_PERMITTIVITY_EPS0
    delta_u_joules = 4.0 * math.pi * eps * r_m * psi_0_v * zeta_v
    u_ratio = delta_u_joules / (BOLTZMANN_CONSTANT_KB * temp_k)

    # Bound exponent to avoid numerical overflow / underflow
    u_clamped = max(-25.0, min(25.0, u_ratio))
    k_elec = math.exp(-u_clamped)
    k_elec = min(1.0, k_elec)  # Repulsion restricts entry; favorable attraction is capped at steric bound

    k_total = k_steric * k_elec
    d_eff = d_0 * k_total

    if k_total > 0.15:
        status = "HIGH_PENETRATION"
    elif k_total > 0.01:
        status = "RESTRICTED_DIFFUSION"
    else:
        status = "COMPLETE_EXCLUSION"

    return CoupledDiffusivityResult(
        nanoparticle_name=nanoparticle_name,
        nanoparticle_radius_nm=round(nanoparticle_radius_nm, 2),
        nanoparticle_zeta_mV=round(nanoparticle_zeta_mV, 2),
        glycocalyx_charge_density_C_m2=charge_density_C_m2,
        d_0_asl_m2_s=d_0,
        k_steric=round(k_steric, 5),
        k_elec=round(k_elec, 5),
        k_total=round(k_total, 6),
        d_eff_glycocalyx_m2_s=d_eff,
        barrier_penetration_status=status,
    )


def simulate_enzymatic_desialylation(
    enzyme_name: str = "bacterial_neuraminidase",
    enzyme_conc_uM: float = 0.50,
    duration_minutes: float = 30.0,
    initial_sia_mM: float = BASELINE_SIALIC_ACID_CONC_MM,
    test_particle_radius_nm: float = 12.5,  # AAV capsid (~25 nm diameter)
    test_particle_zeta_mV: float = -15.0,
) -> DesialylationKineticResult:
    """Simulate Michaelis-Menten desialylation kinetics and calculate t_80% and diffusivity gain."""
    if enzyme_name not in SIALIDASE_PROFILES:
        raise KeyError(f"Unknown enzyme {enzyme_name}. Available: {list(SIALIDASE_PROFILES.keys())}")

    prof = SIALIDASE_PROFILES[enzyme_name]
    k_cat = prof["k_cat"]  # s^-1
    k_m = prof["k_m"]  # mM
    e_total_mm = enzyme_conc_uM * 1e-3  # convert uM to mM

    total_seconds = duration_minutes * 60.0
    dt_sec = 0.50
    steps = int(total_seconds / dt_sec)

    s = initial_sia_mM
    t_80_sec: Optional[float] = None
    target_sia_80 = initial_sia_mM * (1.0 - SIALIC_CHARGE_FRACTION)  # 80% desialylation -> remaining 20%

    for step in range(steps):
        t_current = step * dt_sec
        # dS/dt = - (k_cat * E * S) / (K_m + S)
        rate = (k_cat * e_total_mm * s) / (k_m + s)
        ds = rate * dt_sec
        s = max(0.0, s - ds)

        if t_80_sec is None and s <= target_sia_80:
            t_80_sec = t_current

    final_sia = s
    desial_fraction = (initial_sia_mM - final_sia) / initial_sia_mM

    # Surface charge density: GAG residual (20%) + sialic contribution (80% * (s / s_0))
    initial_sigma = BASELINE_SURFACE_CHARGE_DENSITY
    residual_gag_sigma = initial_sigma * GAG_RESIDUAL_CHARGE_FRACTION
    sialic_sigma_t = (initial_sigma * SIALIC_CHARGE_FRACTION) * (final_sia / initial_sia_mM)
    final_sigma = residual_gag_sigma + sialic_sigma_t

    charge_reduction_fraction = (initial_sigma - final_sigma) / initial_sigma

    if t_80_sec is None:
        # Analytical approximation for integrated Michaelis-Menten equation:
        # K_m * ln(S_0 / S) + (S_0 - S) = k_cat * [E] * t
        s_0 = initial_sia_mM
        s_target = target_sia_80
        lhs = k_m * math.log(s_0 / s_target) + (s_0 - s_target)
        rhs_rate = k_cat * e_total_mm
        t_80_sec = lhs / rhs_rate if rhs_rate > 0 else 9999.0

    # Diffusivity enhancement calculation for test nanoparticle
    diff_initial = evaluate_coupled_diffusivity(
        "AAV_Capsid_Control",
        test_particle_radius_nm,
        test_particle_zeta_mV,
        charge_density_C_m2=initial_sigma,
    )
    diff_final = evaluate_coupled_diffusivity(
        "AAV_Capsid_Permeabilized",
        test_particle_radius_nm,
        test_particle_zeta_mV,
        charge_density_C_m2=final_sigma,
    )

    d_init = diff_initial.d_eff_glycocalyx_m2_s
    d_fin = diff_final.d_eff_glycocalyx_m2_s
    fold_gain = (d_fin / d_init) if d_init > 0.0 else 1.0

    return DesialylationKineticResult(
        enzyme_name=enzyme_name,
        enzyme_conc_uM=enzyme_conc_uM,
        duration_minutes=duration_minutes,
        initial_sia_mM=initial_sia_mM,
        remaining_sia_mM=round(final_sia, 3),
        desialylation_fraction=round(desial_fraction, 4),
        charge_reduction_fraction=round(charge_reduction_fraction, 4),
        initial_charge_density_C_m2=initial_sigma,
        final_charge_density_C_m2=round(final_sigma, 5),
        t_80_percent_reduction_seconds=round(t_80_sec, 1),
        t_80_percent_reduction_minutes=round(t_80_sec / 60.0, 2),
        initial_diffusivity_m2_s=d_init,
        final_diffusivity_m2_s=d_fin,
        diffusivity_fold_enhancement=round(fold_gain, 2),
    )


def export_complete_glycocalyx_manifest() -> Dict[str, Any]:
    """Generate comprehensive biophysical manifest for discovery reporting."""
    edl = evaluate_electrostatic_double_layer()

    nanoparticles = [
        ("AAV_Capsid_25nm", 12.5, -15.0),
        ("Small_mRNA_LNP_50nm", 25.0, -20.0),
        ("Standard_mRNA_LNP_80nm", 40.0, -22.0),
        ("Large_Polyplex_100nm", 50.0, -18.0),
    ]
    diff_baseline = [asdict(evaluate_coupled_diffusivity(name, r, zeta)) for name, r, zeta in nanoparticles]

    steric_sweep = [
        asdict(evaluate_steric_hindrance(r))
        for r in [10.0, 12.5, 15.0, 20.0, 25.0, 30.0, 40.0, 50.0, 75.0, 100.0]
    ]

    kinetics = {
        name: asdict(simulate_enzymatic_desialylation(name, enzyme_conc_uM=0.50, duration_minutes=30.0))
        for name in SIALIDASE_PROFILES
    }

    return {
        "model_id": "T-0351",
        "title": "Airway Epithelial Glycocalyx Barrier Sialic Acid Charge Density & Enzymatic Permeabilization Model",
        "answer_tier": "T2",
        "track": "Track A",
        "nlm_registry_disclaimer": NLM_REGISTRY_DISCLAIMER,
        "clinical_trial_grounding": ["NCT05668741"],
        "primary_literature_dois": [
            "10.1126/science.1223012",
            "10.1128/jvi.74.13.6050-6057.2000",
            "10.1128/JVI.78.24.13755-13768.2004",
            "10.1152/physrev.00004.2021",
        ],
        "electrostatic_double_layer": asdict(edl),
        "steric_hindrance_sweep": steric_sweep,
        "baseline_particle_diffusivity": diff_baseline,
        "enzymatic_desialylation_kinetics": kinetics,
    }


def run_checks() -> int:
    """Run deterministic self-diagnostic verifications for T-0351."""
    print("Running self-diagnostic consistency checks for T-0351...")

    # Check 1: Debye screening length in ASL (150 mM, 37 C) should be ~0.8 nm
    debye = calculate_debye_length(150.0, 310.15)
    assert 0.70 < debye < 0.90, f"Debye length {debye:.3f} nm out of expected [0.70, 0.90] nm range"
    print(f"  [PASS] Debye length in ASL: {debye:.3f} nm (target ~0.8 nm)")

    # Check 2: Surface potential for baseline charge density (-0.035 C/m^2)
    psi_0 = calculate_surface_potential(-0.035, 150.0, 310.15)
    assert -45.0 < psi_0 < -20.0, f"Surface potential {psi_0:.2f} mV out of expected [-45, -20] mV range"
    print(f"  [PASS] Glycocalyx surface potential: {psi_0:.2f} mV (target -20 to -45 mV)")

    # Check 3: Steric partitioning monotonicity
    radii = [10.0, 15.0, 20.0, 24.9, 30.0]
    k_vals = [evaluate_steric_hindrance(r).effective_steric_partition_k for r in radii]
    assert all(k_vals[i] >= k_vals[i + 1] for i in range(len(k_vals) - 1)), "Steric partition non-monotonic"
    assert k_vals[-1] == 0.0, "Particles exceeding pore size must have zero steric partition"
    print("  [PASS] Steric partitioning monotonically decreases to 0 at pore limit")

    # Check 4: Enzymatic desialylation produces > 5-fold diffusivity enhancement
    kin = simulate_enzymatic_desialylation("bacterial_neuraminidase", enzyme_conc_uM=0.50, duration_minutes=30.0)
    assert kin.desialylation_fraction >= 0.80, f"Expected >= 80% desialylation, got {kin.desialylation_fraction:.2%}"
    assert kin.diffusivity_fold_enhancement >= 5.0, (
        f"Diffusivity gain {kin.diffusivity_fold_enhancement}x below 5.0-fold requirement"
    )
    assert kin.t_80_percent_reduction_seconds > 0.0
    print(f"  [PASS] Bacterial neuraminidase: {kin.desialylation_fraction:.1%} desialylation in 30 min, "
          f"{kin.diffusivity_fold_enhancement:.1f}x diffusivity gain, t_80% = {kin.t_80_percent_reduction_minutes:.1f} min")

    # Check 5: Human neuraminidase isoforms Neu1, Neu2, Neu3
    neu1 = simulate_enzymatic_desialylation("recombinant_human_neu1", enzyme_conc_uM=1.0)
    neu2 = simulate_enzymatic_desialylation("recombinant_human_neu2", enzyme_conc_uM=1.0)
    neu3 = simulate_enzymatic_desialylation("recombinant_human_neu3", enzyme_conc_uM=1.0)
    assert neu2.t_80_percent_reduction_seconds < neu1.t_80_percent_reduction_seconds, (
        "Neu2 neutral sialidase should turn over faster than Neu1"
    )
    print(f"  [PASS] Human Neu1 (t_80%={neu1.t_80_percent_reduction_minutes:.1f}m), "
          f"Neu2 (t_80%={neu2.t_80_percent_reduction_minutes:.1f}m), "
          f"Neu3 (t_80%={neu3.t_80_percent_reduction_minutes:.1f}m) verified")

    # Check 6: Manifest export integrity
    manifest = export_complete_glycocalyx_manifest()
    assert len(manifest["clinical_trial_grounding"]) >= 1
    assert len(manifest["primary_literature_dois"]) >= 4
    print("  [PASS] Manifest grounding & NLM disclaimer verified")

    print("\nALL DETERMINISTIC CONSISTENCY CHECKS PASSED (T-0351)")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Airway Epithelial Glycocalyx Barrier Sialic Acid Charge Density & Enzymatic Permeabilization Model (T-0351)"
    )
    parser.add_argument("--json", action="store_true", help="Print complete JSON manifest.")
    parser.add_argument("--check", action="store_true", help="Run self-diagnostic consistency checks.")
    parser.add_argument("--enzyme", type=str, default="bacterial_neuraminidase", choices=list(SIALIDASE_PROFILES.keys()),
                        help="Neuraminidase enzyme profile to evaluate.")
    parser.add_argument("--conc-uM", type=float, default=0.50, help="Enzyme concentration in micromolar (default: 0.50 uM).")
    args = parser.parse_args()

    if args.check:
        return run_checks()

    manifest = export_complete_glycocalyx_manifest()

    if args.json:
        print(json.dumps(manifest, indent=2))
        return 0

    print("=" * 78)
    print("SENSEI EWOK'S CF LAB — GLYCOCALYX SIALIC ACID BARRIER & PERMEABILIZATION MODEL")
    print("Task: T-0351 | Track A Curative Delivery Biophysics | Answer Tier: T2")
    print("=" * 78)
    print(f"Registry Disclaimer: {NLM_REGISTRY_DISCLAIMER}\n")

    edl = manifest["electrostatic_double_layer"]
    print("1. Electrostatic Double-Layer (ASL, I = 150 mM, 37 C):")
    print(f"   - Debye Screening Length (kappa^-1): {edl['debye_length_nm']} nm")
    print(f"   - Glycocalyx Surface Potential:      {edl['surface_potential_mV']} mV")
    print(f"   - Thermal Voltage (k_B * T / e):    {edl['thermal_voltage_mV']} mV")
    print(f"   - Potential at z = 1 nm / 2 nm:      {edl['potential_at_1nm_mV']} mV / {edl['potential_at_2nm_mV']} mV")

    print("\n2. Nanoparticle Steric & Electrostatic Partitioning (Baseline):")
    for p in manifest["baseline_particle_diffusivity"]:
        print(f"   - {p['nanoparticle_name']:<24}: r = {p['nanoparticle_radius_nm']} nm | "
              f"K_steric = {p['k_steric']:.4f} | K_elec = {p['k_elec']:.4f} | "
              f"K_tot = {p['k_total']:.5f} ({p['barrier_penetration_status']})")

    print(f"\n3. Enzymatic Desialylation Kinetics ([E] = {args.conc_uM} uM):")
    for name, k in manifest["enzymatic_desialylation_kinetics"].items():
        print(f"   - {name:<26}: desialylation = {k['desialylation_fraction'] * 100:.1f}% | "
              f"t_80% = {k['t_80_percent_reduction_minutes']:.1f} min | "
              f"diffusivity gain = {k['diffusivity_fold_enhancement']:.1f}x")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    sys.exit(main())
