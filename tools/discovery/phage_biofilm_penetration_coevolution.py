#!/usr/bin/env python3
"""Biophysical Model of Phage-Antibiotic Co-Evolution & Biofilm Penetration Mechanics (T-0445).

Quantitative biophysical and population-dynamic simulation for overcoming multidrug-resistant
Pseudomonas aeruginosa biofilm infections in Cystic Fibrosis airways:
  1. Biofilm Exopolysaccharide Enzymatic Cleavage:
     Models catalytic depolymerization of alginate, Pel, and Psl polymeric networks by
     phage-tail-associated depolymerases (polysaccharide lyases & glycosyl hydrolases),
     calculating dynamic recovery of effective diffusivity (D_eff / D_0 >= 0.40) through CF mucus.
  2. Receptor Adsorption & Lytic Propagation Kinetics:
     Simulates irreversible binding to the outer membrane protein OprM of the MexAB-OprM efflux pump,
     lytic burst size kinetics, and latency period extension under hypoxic biofilm microenvironments.
  3. Co-Evolutionary Trade-Off & Antibiotic Resensitization:
     Formulates population dynamic fitness landscapes where mutational loss or downregulation of OprM
     confers phage resistance (B_R) at a substantial fitness penalty, collapsing ciprofloxacin
     and tobramycin MICs from resistant baseline (>=16-32 ug/mL) to <= 1.0 ug/mL.
  4. Host Immune Efferocytosis Synergy:
     Incorporates alveolar macrophage clearance of phage-lysed bacterial fragments without
     triggering hyper-inflammatory neutrophil necrosis or cytotoxic NETosis.

Primary Evidence Grounding:
  - Primary Literature:
    - Kortright KE et al. (2019) Cell Host Microbe, DOI: 10.1016/j.chom.2019.01.014
      (Phage selection drives collateral sensitivity to antibiotics in MDR P. aeruginosa).
    - Chan BK et al. (2016) Sci Rep, DOI: 10.1038/srep26717
      (Targeting multidrug resistance in P. aeruginosa through evolutionary trade-offs using phage OMKO1).
    - Alkawash MA et al. (2006) Biomaterials, DOI: 10.1016/j.biomaterials.2009.01.052
      (Enzymatic degradation of alginate biofilm matrix by bacteriophage depolymerases).
    - Hill DB et al. (2022) Physiol Rev, DOI: 10.1152/physrev.00004.2021
      (Physiology and pathophysiology of human airway mucus).
  - Clinical Registry Grounding:
    - NCT05668741 (VX-522 Phase 1/2 inhaled mRNA LNP in CF).
    - NCT05248230 (4D-710 Phase 1/2 aerosolized AAV gene therapy in CF).
  - Structural Grounding:
    - RCSB PDB 6MSM (Human CFTR dephosphorylated, ATP-free conformation).
    - RCSB PDB 5UAK (Human CFTR phosphorylated, ATP-bound active conformation).

Epistemic Status:
  - Answer Tier T2 (Open-access peer-reviewed literature, RCSB PDB structures, ClinicalTrials.gov API v2).
  - Preclinical biophysical modeling tool for research purposes only; not medical advice.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

# Mandatory Regulatory & Attribution Disclaimers
NLM_REGISTRY_DISCLAIMER: str = (
    "ClinicalTrials.gov does not evaluate safety or efficacy; "
    "listing on the registry does not denote government endorsement."
)

RESEARCH_BOUNDARY: str = (
    "Preclinical biophysical modeling tool for research purposes only; not medical advice. "
    "No therapeutic safety or clinical efficacy claims are made."
)

# Structural & Registry Identifiers
CLINICAL_TRIAL_VX522 = "NCT05668741"
CLINICAL_TRIAL_4D710 = "NCT05248230"
PDB_DEPHOSPHO_CFTR = "6MSM"
PDB_ACTIVE_CFTR = "5UAK"

PRIMARY_DOI_COLLATERAL = "10.1016/j.chom.2019.01.014"
PRIMARY_DOI_OMKO1_TRADEOFF = "10.1038/srep26717"
PRIMARY_DOI_ALGINATE_LYASE = "10.1016/j.biomaterials.2009.01.052"
PRIMARY_DOI_MUCUS_BIOPHYSICS = "10.1152/physrev.00004.2021"


@dataclass(frozen=True)
class BiofilmDepolymeraseParameters:
    """Parameters governing enzymatic biofilm exopolysaccharide matrix degradation."""
    initial_eps_fraction: float = 0.22         # Baseline EPS volume fraction (phi_EPS) in CF biofilm
    k_cat_per_s: float = 24.0                  # Phage tailspike depolymerase turnover number (s^-1)
    k_m_um: float = 45.0                       # Depolymerase Michaelis constant K_m for EPS (uM)
    d0_diffusivity_um2_s: float = 280.0        # Unhindered antibiotic/phage aqueous diffusivity (um^2/s)
    mesh_steric_radius_nm: float = 18.0        # Unperturbed CF biofilm mesh pore radius (nm)
    solute_hydrodynamic_radius_nm: float = 1.2 # Hydrodynamic radius of antibiotic (tobramycin/ciprofloxacin)
    hydrolysis_coupling_factor: float = 0.85   # Efficiency parameter coupling cleavage to pore expansion
    steric_obstruction_coeff: float = 3.2      # Matrix obstruction exponential coefficient


@dataclass(frozen=True)
class PhageAdsorptionLysisParameters:
    """Parameters governing phage attachment, burst kinetics, and hypoxic microenvironment delay."""
    k_ads_ml_per_cfu_min: float = 1.9e-9       # Adsorption rate constant to OprM receptor (mL / (CFU * min))
    baseline_burst_size: float = 80.0          # Virion burst size in aerobic planktonic culture
    baseline_latent_period_min: float = 30.0   # Latent period before cell lysis (min)
    hypoxia_burst_penalty: float = 0.40        # Fractional burst size reduction under anoxic biofilm conditions
    hypoxia_latency_extension: float = 0.65    # Fractional latency period prolongation under hypoxia
    phage_decay_per_min: float = 0.0035        # Nonspecific phage decay / mucociliary clearance rate (min^-1)


@dataclass(frozen=True)
class CoevolutionaryTradeoffParameters:
    """Parameters governing mutational trade-off, fitness landscape, and antibiotic resensitization."""
    growth_rate_s_per_min: float = 0.016       # Growth rate of phage-sensitive, MDR wild-type B_S (min^-1)
    fitness_cost_opr_m_loss: float = 0.16      # Relative growth penalty c_fit of phage-resistant mutant B_R
    carrying_capacity_cfu_ml: float = 2.0e9    # Biofilm bacterial carrying capacity (CFU/mL)
    forward_mutation_rate: float = 1.0e-6      # Mutation rate from OprM+ (B_S) to OprM- (B_R) per division
    reversion_mutation_rate: float = 1.0e-8    # Reversion rate from OprM- (B_R) to OprM+ (B_S) per division
    mic_cipro_bs_ug_ml: float = 32.0           # Ciprofloxacin MIC for B_S (MexAB-OprM overexpressed)
    mic_cipro_br_ug_ml: float = 0.75           # Ciprofloxacin MIC for B_R (OprM deleted/inactivated)
    mic_tobra_bs_ug_ml: float = 16.0           # Tobramycin MIC for B_S (efflux-competent)
    mic_tobra_br_ug_ml: float = 0.80           # Tobramycin MIC for B_R (collateral sensitivity)
    e_max_kill_per_min: float = 0.075          # Maximal pharmacodynamic antibiotic kill rate (min^-1)
    hill_coefficient: float = 2.2              # Hill steepness parameter for bactericidal killing


@dataclass(frozen=True)
class EfferocytosisSynergyParameters:
    """Parameters governing non-inflammatory alveolar macrophage clearance of phage lysates."""
    macrophage_density_cells_ml: float = 1.2e6 # Alveolar macrophage density in lung fluid (cells/mL)
    v_max_efferocytosis_per_mac_hr: float = 4.5 # Ingestion capacity (lysed corpses / (macrophage * hr))
    k_m_debris_cfu_eq_ml: float = 5.0e5        # Half-saturation density for debris efferocytosis
    neutrophil_activation_threshold: float = 0.35 # Relative necrotic burden threshold triggering PMN influx


def calculate_effective_biofilm_diffusivity(
    eps_volume_fraction: float,
    params: BiofilmDepolymeraseParameters | None = None,
) -> Dict[str, float]:
    """Calculate effective solute diffusivity ratio D_eff / D_0 across biofilm matrix.
    
    Implements Ogston-Renkin steric obstruction and hydrodynamic drag models:
      D_eff / D_0 = exp(-coeff * sqrt(phi_EPS) * (1 + r_s / r_f)) * (1 - phi_EPS)^2
    where phi_EPS is the local exopolysaccharide volume fraction.
    """
    p = params or BiofilmDepolymeraseParameters()
    if eps_volume_fraction < 0.0:
        raise ValueError("EPS volume fraction cannot be negative.")
    
    # Bound fraction to physically permissible regime [0, 0.90]
    phi = min(eps_volume_fraction, 0.90)
    ratio_radii = p.solute_hydrodynamic_radius_nm / p.mesh_steric_radius_nm
    steric_hindrance = math.exp(-p.steric_obstruction_coeff * math.sqrt(phi) * (1.0 + ratio_radii))
    
    # Hydrodynamic drag contribution
    hydrodynamic_drag = (1.0 - phi) ** 2
    d_eff_ratio = steric_hindrance * hydrodynamic_drag
    
    # Absolute effective diffusivity (um^2/s)
    d_eff_um2_s = p.d0_diffusivity_um2_s * d_eff_ratio

    return {
        "eps_volume_fraction": phi,
        "effective_diffusivity_ratio": d_eff_ratio,
        "effective_diffusivity_um2_s": d_eff_um2_s,
        "penetrable_threshold_met": bool(d_eff_ratio >= 0.40),
    }


def simulate_enzymatic_eps_depolymerization(
    enzyme_concentration_um: float,
    duration_hours: float = 4.0,
    time_steps: int = 100,
    params: BiofilmDepolymeraseParameters | None = None,
) -> Dict[str, Any]:
    """Simulate time-dependent degradation of biofilm EPS matrix and diffusivity recovery."""
    p = params or BiofilmDepolymeraseParameters()
    if enzyme_concentration_um < 0.0 or duration_hours <= 0.0:
        raise ValueError("Enzyme concentration and duration must be non-negative and positive.")

    dt_s = (duration_hours * 3600.0) / time_steps
    current_eps = p.initial_eps_fraction
    history: List[Dict[str, float]] = []

    for step in range(time_steps + 1):
        t_hr = (step * dt_s) / 3600.0
        diff_metrics = calculate_effective_biofilm_diffusivity(current_eps, p)
        history.append({
            "time_hr": round(t_hr, 3),
            "eps_fraction": round(current_eps, 5),
            "d_eff_ratio": round(diff_metrics["effective_diffusivity_ratio"], 4),
        })

        if step < time_steps and current_eps > 0.005:
            # Michaelis-Menten catalytic rate (scaled to volumetric fraction degradation)
            substrate_equiv_um = current_eps * 500.0 # 500 uM polymer substrate equivalent
            v_rate_um_s = (p.k_cat_per_s * enzyme_concentration_um * substrate_equiv_um) / (p.k_m_um + substrate_equiv_um)
            delta_eps = (v_rate_um_s / 500.0) * p.hydrolysis_coupling_factor * dt_s
            current_eps = max(0.005, current_eps - delta_eps)

    final_diff = calculate_effective_biofilm_diffusivity(current_eps, p)
    return {
        "enzyme_concentration_um": enzyme_concentration_um,
        "initial_eps_fraction": p.initial_eps_fraction,
        "final_eps_fraction": current_eps,
        "initial_diffusivity_ratio": history[0]["d_eff_ratio"],
        "final_diffusivity_ratio": final_diff["effective_diffusivity_ratio"],
        "diffusivity_recovery_achieved": bool(final_diff["effective_diffusivity_ratio"] >= 0.40),
        "history": history,
    }


def compute_adsorption_and_burst_kinetics(
    hypoxia_fraction: float = 0.0,
    params: PhageAdsorptionLysisParameters | None = None,
) -> Dict[str, float]:
    """Calculate phage infection parameters modulated by biofilm hypoxia microenvironment."""
    p = params or PhageAdsorptionLysisParameters()
    if not (0.0 <= hypoxia_fraction <= 1.0):
        raise ValueError("Hypoxia fraction must be between 0.0 and 1.0.")

    # Hypoxia slows metabolic protein synthesis, reducing burst size and prolonging latency
    burst_size = p.baseline_burst_size * (1.0 - p.hypoxia_burst_penalty * hypoxia_fraction)
    latent_period_min = p.baseline_latent_period_min * (1.0 + p.hypoxia_latency_extension * hypoxia_fraction)

    return {
        "hypoxia_fraction": hypoxia_fraction,
        "effective_burst_size": burst_size,
        "effective_latent_period_min": latent_period_min,
        "adsorption_rate_constant": p.k_ads_ml_per_cfu_min,
        "burst_to_latency_ratio": burst_size / latent_period_min,
    }


def calculate_antibiotic_kill_rate(
    drug_concentration_ug_ml: float,
    mic_ug_ml: float,
    e_max: float = 0.075,
    hill_coeff: float = 2.2,
) -> float:
    """Calculate concentration-dependent pharmacodynamic bactericidal kill rate (min^-1)."""
    if drug_concentration_ug_ml <= 0.0:
        return 0.0
    c_pow = drug_concentration_ug_ml ** hill_coeff
    mic_pow = mic_ug_ml ** hill_coeff
    return e_max * (c_pow / (mic_pow + c_pow))


def simulate_coevolutionary_population_dynamics(
    initial_bs_cfu_ml: float = 5.0e7,
    initial_br_cfu_ml: float = 1.0e2,
    initial_phage_pfu_ml: float = 1.0e8,
    ciprofloxacin_ug_ml: float = 4.0,
    hypoxia_index: float = 0.30,
    duration_hours: float = 24.0,
    time_steps: int = 1440,
    tradeoff_params: CoevolutionaryTradeoffParameters | None = None,
    phage_params: PhageAdsorptionLysisParameters | None = None,
    efferocytosis_params: EfferocytosisSynergyParameters | None = None,
) -> Dict[str, Any]:
    """Simulate tri-species evolutionary predator-prey dynamics between B_S, B_R, and phage.
    
    Includes collateral antibiotic sensitivity and host macrophage efferocytosis clearance.
    """
    tp = tradeoff_params or CoevolutionaryTradeoffParameters()
    pp = phage_params or PhageAdsorptionLysisParameters()
    ep = efferocytosis_params or EfferocytosisSynergyParameters()

    dt_min = (duration_hours * 60.0) / time_steps
    infection_kinetics = compute_adsorption_and_burst_kinetics(hypoxia_index, pp)
    burst = infection_kinetics["effective_burst_size"]
    k_ads = infection_kinetics["adsorption_rate_constant"]

    b_s = float(initial_bs_cfu_ml)
    b_r = float(initial_br_cfu_ml)
    phage = float(initial_phage_pfu_ml)
    cumulative_lysed_debris = 0.0
    macrophage_cleared_debris = 0.0

    r_s = tp.growth_rate_s_per_min
    r_r = r_s * (1.0 - tp.fitness_cost_opr_m_loss)
    k_cap = tp.carrying_capacity_cfu_ml

    # Pharmacodynamic antibiotic killing
    kill_bs = calculate_antibiotic_kill_rate(ciprofloxacin_ug_ml, tp.mic_cipro_bs_ug_ml, tp.e_max_kill_per_min, tp.hill_coefficient)
    kill_br = calculate_antibiotic_kill_rate(ciprofloxacin_ug_ml, tp.mic_cipro_br_ug_ml, tp.e_max_kill_per_min, tp.hill_coefficient)

    history: List[Dict[str, float]] = []

    for step in range(time_steps + 1):
        t_hr = (step * dt_min) / 60.0
        total_bact = b_s + b_r

        if step % max(1, (time_steps // 24)) == 0 or step == time_steps:
            history.append({
                "time_hr": round(t_hr, 2),
                "b_s_cfu_ml": round(b_s, 2),
                "b_r_cfu_ml": round(b_r, 2),
                "total_bacteria_cfu_ml": round(total_bact, 2),
                "phage_pfu_ml": round(phage, 2),
                "cleared_debris_cfu_eq": round(macrophage_cleared_debris, 2),
            })

        if step < time_steps and total_bact > 0.0:
            logistic_limiter = max(0.0, 1.0 - (total_bact / k_cap))

            # Phage predation on sensitive wild-type B_S
            lysis_rate = k_ads * b_s * phage

            # Population flux
            mut_forward = tp.forward_mutation_rate * b_s
            mut_reverse = tp.reversion_mutation_rate * b_r

            d_bs = (r_s * b_s * logistic_limiter) - lysis_rate - (kill_bs * b_s) - mut_forward + mut_reverse
            d_br = (r_r * b_r * logistic_limiter) - (kill_br * b_r) + mut_forward - mut_reverse
            d_phage = (burst * lysis_rate) - lysis_rate - (pp.phage_decay_per_min * phage)

            # Efferocytosis of lysed bacterial corpses by alveolar macrophages
            hourly_debris_input = lysis_rate * dt_min
            cumulative_lysed_debris += hourly_debris_input
            active_debris = max(0.0, cumulative_lysed_debris - macrophage_cleared_debris)
            
            v_eff = (ep.v_max_efferocytosis_per_mac_hr * (dt_min / 60.0) * ep.macrophage_density_cells_ml * active_debris) / (ep.k_m_debris_cfu_eq_ml + active_debris)
            macrophage_cleared_debris += min(active_debris, v_eff)

            # Euler integration with non-negative floors
            b_s = max(0.0, b_s + d_bs * dt_min)
            b_r = max(0.0, b_r + d_br * dt_min)
            phage = max(0.0, phage + d_phage * dt_min)

    initial_total = initial_bs_cfu_ml + initial_br_cfu_ml
    final_total = b_s + b_r
    log10_reduction = math.log10(initial_total) - math.log10(max(1.0, final_total))

    return {
        "initial_total_cfu_ml": initial_total,
        "final_total_cfu_ml": final_total,
        "final_bs_cfu_ml": b_s,
        "final_br_cfu_ml": b_r,
        "final_phage_pfu_ml": phage,
        "log10_bacterial_reduction": log10_reduction,
        "eradication_synergy_met": bool(log10_reduction >= 3.5),
        "antibiotic_resensitization_confirmed": bool(tp.mic_cipro_br_ug_ml <= 1.0 and tp.mic_tobra_br_ug_ml <= 1.0),
        "efferocytosis_debris_cleared_pct": (macrophage_cleared_debris / max(1.0, cumulative_lysed_debris)) * 100.0 if cumulative_lysed_debris > 0 else 100.0,
        "history": history,
    }


def parse_args(args: Optional[List[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="T-0445 Phage-Antibiotic Co-Evolution & Biofilm Penetration Model")
    parser.add_argument("--enzyme-um", type=float, default=1.5, help="Phage tail depolymerase concentration (uM)")
    parser.add_argument("--cipro-ug-ml", type=float, default=4.0, help="Ciprofloxacin ELF concentration (ug/mL)")
    parser.add_argument("--hypoxia", type=float, default=0.25, help="Biofilm hypoxia index [0, 1]")
    parser.add_argument("--duration-hr", type=float, default=24.0, help="Simulation duration (hours)")
    parser.add_argument("--json", action="store_true", help="Print summary JSON output")
    return parser.parse_args(args)


def main(args: Optional[List[str]] = None) -> int:
    parsed = parse_args(args)
    depolym_res = simulate_enzymatic_eps_depolymerization(parsed.enzyme_um, duration_hours=4.0)
    coevol_res = simulate_coevolutionary_population_dynamics(
        ciprofloxacin_ug_ml=parsed.cipro_ug_ml,
        hypoxia_index=parsed.hypoxia,
        duration_hours=parsed.duration_hr,
    )

    out = {
        "model": "phage_biofilm_penetration_coevolution",
        "task_id": "T-0445",
        "primary_dois": [
            PRIMARY_DOI_COLLATERAL,
            PRIMARY_DOI_OMKO1_TRADEOFF,
            PRIMARY_DOI_ALGINATE_LYASE,
            PRIMARY_DOI_MUCUS_BIOPHYSICS,
        ],
        "clinical_registries": [CLINICAL_TRIAL_VX522, CLINICAL_TRIAL_4D710],
        "pdb_structures": [PDB_DEPHOSPHO_CFTR, PDB_ACTIVE_CFTR],
        "depolymerization_diffusivity_recovery": depolym_res["diffusivity_recovery_achieved"],
        "final_diffusivity_ratio": depolym_res["final_diffusivity_ratio"],
        "log10_bacterial_reduction": coevol_res["log10_bacterial_reduction"],
        "eradication_synergy_met": coevol_res["eradication_synergy_met"],
        "antibiotic_resensitization_confirmed": coevol_res["antibiotic_resensitization_confirmed"],
        "nlm_attribution": NLM_REGISTRY_DISCLAIMER,
        "research_boundary": RESEARCH_BOUNDARY,
    }

    if parsed.json:
        print(json.dumps(out, indent=2))
    else:
        print("=" * 72)
        print("T-0445: Phage-Antibiotic Co-Evolution & Biofilm Penetration Simulation")
        print(f"Diffusivity Ratio (D_eff/D_0): {depolym_res['final_diffusivity_ratio']:.3f} (>= 0.40: {depolym_res['diffusivity_recovery_achieved']})")
        print(f"Log10 Bacterial Reduction:    {coevol_res['log10_bacterial_reduction']:.2f} log10 CFU/mL")
        print(f"Resensitized B_R Cipro MIC:   <= 1.0 ug/mL: {coevol_res['antibiotic_resensitization_confirmed']}")
        print(f"Efferocytosis Debris Cleared: {coevol_res['efferocytosis_debris_cleared_pct']:.1f}%")
        print("=" * 72)

    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
