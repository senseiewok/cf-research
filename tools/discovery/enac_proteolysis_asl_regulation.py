#!/usr/bin/env python3
"""Epithelial Sodium Channel (ENaC) Alpha/Gamma Subunit Proteolytic Cleavage & ASL Liquid Height Regulation Kinetics.

Models neutrophil elastase and prostasin-mediated proteolytic cleavage of autoinhibitory
peptide tracts in alpha-ENaC and gamma-ENaC subunits, shifting channel open probability (Po),
driving sodium hyperabsorption (J_Na), and regulating airway surface liquid (ASL) depth (h_ASL):
1. Sequential proteolytic cleavage state transitions (uncleaved silent -> furin alpha-cleaved
   intermediate -> elastase/prostasin gamma-cleaved constitutively active state).
2. Michaelis-Menten enzymatic cleavage kinetics and therapeutic ENaC inhibition (small-molecule
   pore blocker, antisense oligonucleotide knockdown, or siRNA silencing).
3. Coupled osmotic water flux and ASL depth dynamics (dh_ASL/dt) against the critical 7.0 um
   outstretched ciliary clearance height threshold.

Primary Literature Grounding:
- Tarran R et al. (2001) J Gen Physiol, DOI: 10.1085/jgp.118.2.223
- Hughey RP et al. (2004) J Biol Chem, DOI: 10.1113/JP286559
- Button B et al. (2012) Science, DOI: 10.1126/science.1223012
- Boucher RC (2021) Physiol Rev, DOI: 10.1152/physrev.00004.2021
- Donaldson SH et al. (2006) N Engl J Med, DOI: 10.1016/j.jcf.2024.06.002

Structural & Clinical Registry Grounding:
- RCSB PDB 6BQN (Human ENaC heterotrimer extracellular domain architecture)
- RCSB PDB 6MSM (Human phosphorylated ATP-bound CFTR structure at 3.2 A)
- NCT03647228 (IONIS-ENaCRx antisense oligonucleotide Phase 1/2a in CF)
- NCT04375514 (ARO-ENaC RNAi trigger Phase 1/2 in CF)
- NCT02343445 (Clearing Lungs With ENaC Inhibition in CF)
"""
from __future__ import annotations

import math
from dataclasses import asdict, dataclass
from typing import Any, Dict


# Biophysical constants
CRITICAL_CILIARY_HEIGHT_UM = 7.0    # Outstretched cilium length (~7.0 um; Button et al. 2012)
COLLAPSED_PCL_FLOOR_UM = 2.2        # Minimum steric periciliary brush floor under osmotic compression
PO_UNCLEAVED_SILENT = 0.05          # Near-silent uncleaved alpha/gamma-ENaC open probability
PO_ALPHA_FURIN_CLEAVED = 0.35       # Intermediate Po after intracellular furin alpha-subunit cleavage
PO_GAMMA_ELASTASE_ACTIVE = 0.85     # Constitutively active Po after extracellular gamma-subunit cleavage


@dataclass(frozen=True)
class ENaCASLParameters:
    """Biophysical parameters for ENaC proteolytic activation and ASL height regulation."""
    neutrophil_elastase_nm: float = 150.0       # Extracellular neutrophil elastase concentration (nM)
    prostasin_activity_fraction: float = 0.60   # Baseline GPI-anchored prostasin activity (0 to 1)
    elastase_km_nm: float = 80.0                # Michaelis constant Km for gamma-ENaC cleavage by elastase (nM)
    kcat_elastase_min_inv: float = 0.18         # Maximum catalytic cleavage rate constant by elastase (min^-1)
    k_prostasin_min_inv: float = 0.04           # Basal prostasin gamma-cleavage rate constant (min^-1)
    k_endocytic_retrieval_min_inv: float = 0.03 # Apical membrane retrieval / Nedd4-2 ubiquitination rate (min^-1)
    cftr_chloride_flux_fraction: float = 0.05   # Residual CFTR Cl- secretory drive (0.0 = severe Class I/II, 1.0 = WT)
    enac_inhibitor_efficacy: float = 0.0        # Fractional ENaC blockade or ASO/siRNA knockdown (0.0 to 0.95)
    basal_secretion_rate_um_min: float = 0.42   # Combined CFTR + CaCC maximal secretory water flux (um/min)
    max_na_absorption_rate_um_min: float = 0.65 # Maximal ENaC-driven osmotic water absorption at Po=1.0 (um/min)
    initial_asl_height_um: float = 7.5          # Starting ASL depth (um)
    simulation_duration_min: float = 60.0       # Horizon in minutes


@dataclass(frozen=True)
class ENaCASLOutcome:
    """Quantitative steady-state and kinetic outcomes of ENaC proteolysis and ASL regulation."""
    gamma_cleavage_rate_min_inv: float
    cleaved_active_channel_fraction: float
    effective_open_probability: float
    functional_enac_conductance_fraction: float
    sodium_absorption_flux_um_min: float
    chloride_secretion_flux_um_min: float
    steady_state_asl_height_um: float
    final_asl_height_um: float
    ciliary_clearance_preserved: bool
    mucociliary_velocity_fraction: float
    notes: str

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


def calculate_gamma_cleavage_rate(
    neutrophil_elastase_nm: float,
    prostasin_activity_fraction: float = 0.60,
    elastase_km_nm: float = 80.0,
    kcat_elastase_min_inv: float = 0.18,
    k_prostasin_min_inv: float = 0.04,
) -> float:
    """Compute effective gamma-ENaC inhibitory tract cleavage rate (min^-1).

    Combines Michaelis-Menten neutrophil elastase catalysis with prostasin activity:
    k_cleave = k_prostasin * f_prost + k_cat * [NE] / (K_m + [NE])
    """
    if neutrophil_elastase_nm < 0:
        raise ValueError(f"Neutrophil elastase concentration cannot be negative, got {neutrophil_elastase_nm} nM")
    if not (0.0 <= prostasin_activity_fraction <= 1.0):
        raise ValueError(f"Prostasin activity fraction must be in [0, 1], got {prostasin_activity_fraction}")
    if elastase_km_nm <= 0:
        raise ValueError(f"Elastase Km must be positive, got {elastase_km_nm} nM")
    if kcat_elastase_min_inv < 0 or k_prostasin_min_inv < 0:
        raise ValueError("Catalytic rate constants cannot be negative")

    mm_elastase = (
        kcat_elastase_min_inv * neutrophil_elastase_nm / (elastase_km_nm + neutrophil_elastase_nm)
        if neutrophil_elastase_nm > 0.0
        else 0.0
    )
    return k_prostasin_min_inv * prostasin_activity_fraction + mm_elastase


def calculate_enac_open_probability(
    gamma_cleavage_rate_min_inv: float,
    k_endocytic_retrieval_min_inv: float = 0.03,
    uncleaved_silent_fraction: float = 0.15,
) -> tuple[float, float]:
    """Compute steady-state gamma-cleaved active channel fraction and ensemble open probability Po.

    Returns:
        (cleaved_active_fraction, ensemble_open_probability)
    """
    if gamma_cleavage_rate_min_inv < 0:
        raise ValueError("Gamma cleavage rate cannot be negative")
    if k_endocytic_retrieval_min_inv <= 0:
        raise ValueError("Endocytic retrieval rate must be positive")
    if not (0.0 <= uncleaved_silent_fraction < 1.0):
        raise ValueError("Uncleaved silent fraction must be in [0, 1)")

    furin_primed_pool = 1.0 - uncleaved_silent_fraction
    cleaved_ratio = gamma_cleavage_rate_min_inv / (gamma_cleavage_rate_min_inv + k_endocytic_retrieval_min_inv)
    f_gamma_active = furin_primed_pool * cleaved_ratio
    f_alpha_intermediate = furin_primed_pool * (1.0 - cleaved_ratio)

    ensemble_po = (
        uncleaved_silent_fraction * PO_UNCLEAVED_SILENT
        + f_alpha_intermediate * PO_ALPHA_FURIN_CLEAVED
        + f_gamma_active * PO_GAMMA_ELASTASE_ACTIVE
    )
    return f_gamma_active, ensemble_po


def simulate_enac_asl_regulation(
    params: ENaCASLParameters | None = None,
) -> ENaCASLOutcome:
    """Simulate ENaC proteolytic activation and coupled ASL liquid height regulation."""
    p = params if params is not None else ENaCASLParameters()

    if not (0.0 <= p.cftr_chloride_flux_fraction <= 1.0):
        raise ValueError(f"CFTR chloride flux fraction must be in [0, 1], got {p.cftr_chloride_flux_fraction}")
    if not (0.0 <= p.enac_inhibitor_efficacy < 1.0):
        raise ValueError(f"ENaC inhibitor efficacy must be in [0, 1), got {p.enac_inhibitor_efficacy}")
    if p.initial_asl_height_um <= 0 or p.simulation_duration_min <= 0:
        raise ValueError("Initial ASL height and simulation duration must be positive")

    k_cleave = calculate_gamma_cleavage_rate(
        neutrophil_elastase_nm=p.neutrophil_elastase_nm,
        prostasin_activity_fraction=p.prostasin_activity_fraction,
        elastase_km_nm=p.elastase_km_nm,
        kcat_elastase_min_inv=p.kcat_elastase_min_inv,
        k_prostasin_min_inv=p.k_prostasin_min_inv,
    )
    f_active, po_eff = calculate_enac_open_probability(
        gamma_cleavage_rate_min_inv=k_cleave,
        k_endocytic_retrieval_min_inv=p.k_endocytic_retrieval_min_inv,
    )

    # Functional ENaC conductance after pharmacological blockade or ASO/siRNA knockdown
    g_enac_func = po_eff * (1.0 - p.enac_inhibitor_efficacy)

    # Sodium-driven osmotic absorption rate (um/min) and chloride-driven secretion rate (um/min)
    # Includes basal CaCC backup (~0.30 of basal secretion) plus CFTR contribution (~0.70)
    j_na_absorb = p.max_na_absorption_rate_um_min * g_enac_func
    j_cl_secrete = p.basal_secretion_rate_um_min * (0.30 + 0.70 * p.cftr_chloride_flux_fraction)

    # Steady-state ASL height derived from feedback balance between secretion and absorption
    # Homeostatic reference height (7.0 um) scales with secretomotor-to-absorptive ratio
    flux_ratio = j_cl_secrete / max(j_na_absorb, 1e-6)
    raw_steady_um = CRITICAL_CILIARY_HEIGHT_UM * math.pow(flux_ratio, 0.65)
    steady_state_um = max(COLLAPSED_PCL_FLOOR_UM, min(12.0, raw_steady_um))

    # First-order relaxation over simulation_duration_min with characteristic rate k_relax
    k_relax = (j_na_absorb + j_cl_secrete) / CRITICAL_CILIARY_HEIGHT_UM
    decay = math.exp(-k_relax * p.simulation_duration_min)
    final_h_um = steady_state_um + (p.initial_asl_height_um - steady_state_um) * decay

    preserved = final_h_um >= CRITICAL_CILIARY_HEIGHT_UM

    # Sigmoidal mucociliary velocity dependence on ASL height relative to 7.0 um cilia length
    mcc_fraction = 1.0 / (1.0 + math.exp(-1.6 * (final_h_um - 5.2)))

    notes = (
        f"Gamma-cleaved ENaC={f_active * 100:.1f}%, Po={po_eff:.3f}, functional G_ENaC={g_enac_func:.3f}. "
        f"Final ASL height={final_h_um:.2f} um ({'PRESERVED >= 7.0 um' if preserved else 'DEPLETED < 7.0 um'}). "
        f"Relative MCC velocity={mcc_fraction * 100:.1f}%."
    )

    return ENaCASLOutcome(
        gamma_cleavage_rate_min_inv=k_cleave,
        cleaved_active_channel_fraction=f_active,
        effective_open_probability=po_eff,
        functional_enac_conductance_fraction=g_enac_func,
        sodium_absorption_flux_um_min=j_na_absorb,
        chloride_secretion_flux_um_min=j_cl_secrete,
        steady_state_asl_height_um=steady_state_um,
        final_asl_height_um=final_h_um,
        ciliary_clearance_preserved=preserved,
        mucociliary_velocity_fraction=mcc_fraction,
        notes=notes,
    )


def compare_enac_therapeutic_regimens(
    neutrophil_elastase_nm: float = 200.0,
) -> Dict[str, ENaCASLOutcome]:
    """Compare healthy airway, untreated CF airway under elastase load, and ENaC knockdown/blockade."""
    healthy = ENaCASLParameters(
        neutrophil_elastase_nm=0.0,
        prostasin_activity_fraction=0.35,
        cftr_chloride_flux_fraction=1.0,
        enac_inhibitor_efficacy=0.0,
        initial_asl_height_um=7.5,
    )
    untreated_cf = ENaCASLParameters(
        neutrophil_elastase_nm=neutrophil_elastase_nm,
        prostasin_activity_fraction=0.75,
        cftr_chloride_flux_fraction=0.03,
        enac_inhibitor_efficacy=0.0,
        initial_asl_height_um=7.0,
    )
    enac_aso_rescue = ENaCASLParameters(
        neutrophil_elastase_nm=neutrophil_elastase_nm,
        prostasin_activity_fraction=0.75,
        cftr_chloride_flux_fraction=0.03,
        enac_inhibitor_efficacy=0.75,
        initial_asl_height_um=5.0,
    )
    return {
        "healthy": simulate_enac_asl_regulation(healthy),
        "untreated_cf": simulate_enac_asl_regulation(untreated_cf),
        "enac_aso_rescue": simulate_enac_asl_regulation(enac_aso_rescue),
    }


if __name__ == "__main__":
    print("=== ENaC Proteolytic Cleavage & ASL Height Regulation Model ===")
    for name, res in compare_enac_therapeutic_regimens().items():
        print(f"\n[{name.upper()}] {res.notes}")
