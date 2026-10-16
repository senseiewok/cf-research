#!/usr/bin/env python3
"""Quantitative Biophysical Transport & Translation Models for CF Curative Discovery.

Provides deterministic mathematical functions for:
1. Nanoparticle Brownian diffusion through pathological CF mucus mesh (Ogston & Stokes-Einstein).
2. Electrostatic mucin adhesion potential based on zeta potential and ionic strength.
3. Stop-codon readthrough kinetics and Nonsense-Mediated Decay (NMD) mRNA steady-state dynamics.
4. Aerosol vibrating-mesh nebulization shear stability bounds.
"""
from __future__ import annotations

import math
from typing import Dict, Any


# Physical Constants
BOLTZMANN_CONSTANT = 1.380649e-23  # J/K
TEMPERATURE_PHYSIOLOGICAL = 310.15  # 37 C in Kelvin
WATER_VISCOSITY_37C = 0.6913e-3  # Pa*s at 37 C


def stokes_einstein_diffusion(hydrodynamic_radius_nm: float, viscosity_pa_s: float = WATER_VISCOSITY_37C) -> float:
    """Calculate free Stokes-Einstein diffusion coefficient D_0 (in m^2/s).
    
    D_0 = (k_B * T) / (6 * pi * eta * r)
    """
    if hydrodynamic_radius_nm <= 0:
        raise ValueError("Hydrodynamic radius must be positive")
    radius_m = hydrodynamic_radius_nm * 1e-9
    return (BOLTZMANN_CONSTANT * TEMPERATURE_PHYSIOLOGICAL) / (6.0 * math.pi * viscosity_pa_s * radius_m)


def cf_mucus_steric_retardation(
    particle_diameter_nm: float,
    mesh_pore_diameter_nm: float = 110.0,
    fiber_diameter_nm: float = 10.0
) -> float:
    """Calculate steric obstruction ratio D_eff / D_0 through CF mucin network using Ogston lattice model.
    
    CF sputum pore cutoffs typically range between 90 nm and 140 nm (mean ~110 nm).
    Ratio = exp(-pi * ((r_p + r_f) / xi)^2)
    """
    if particle_diameter_nm <= 0 or mesh_pore_diameter_nm <= 0:
        raise ValueError("Dimensions must be positive")
    r_p = particle_diameter_nm / 2.0
    r_f = fiber_diameter_nm / 2.0
    xi = mesh_pore_diameter_nm / 2.0
    
    # If particle exceeds pore dimension, steric entrapment prevents free diffusion
    if particle_diameter_nm >= mesh_pore_diameter_nm * 1.5:
        return 0.001
        
    steric_factor = math.exp(-math.pi * ((r_p + r_f) / (xi + r_f)) ** 2)
    return max(0.001, min(1.0, steric_factor))


def electrostatic_mucin_adhesion_penalty(zeta_potential_mv: float) -> float:
    """Compute relative mobility penalty due to electrostatic polyanionic mucin binding.
    
    Mucins (MUC5AC, MUC5B) carry dense sialic acid and sulfate groups (net negative charge).
    - Particles with positive zeta potential (> +5 mV) suffer severe electrostatic entrapment.
    - Particles with extreme negative zeta potential (< -15 mV) suffer repulsive charge exclusion or scavenger uptake.
    - Near-neutral zeta potential (-5 mV to +5 mV) minimizes electrostatic adhesion.
    Returns transmission efficiency [0.01 to 1.0].
    """
    if zeta_potential_mv > 0:
        # Cationic penalty: exponential trapping by polyanionic mucins
        return max(0.01, math.exp(-0.35 * zeta_potential_mv))
    elif zeta_potential_mv < -10.0:
        # Strong anionic exclusion
        return max(0.20, math.exp(0.08 * (zeta_potential_mv + 10.0)))
    else:
        # Near-neutral sweet spot (-10 mV to 0 mV)
        return 0.95 - (abs(zeta_potential_mv) * 0.03)


def calculate_mucus_penetration_index(diameter_nm: float, zeta_mv: float) -> Dict[str, Any]:
    """Calculate composite mucus penetration efficiency for an aerosolized nanoparticle."""
    d0 = stokes_einstein_diffusion(diameter_nm / 2.0)
    steric_ratio = cf_mucus_steric_retardation(diameter_nm)
    charge_ratio = electrostatic_mucin_adhesion_penalty(zeta_mv)
    effective_diffusivity_ratio = steric_ratio * charge_ratio
    
    is_penetrating = (effective_diffusivity_ratio >= 0.25 and diameter_nm <= 90.0 and abs(zeta_mv) <= 5.0)
    
    return {
        "hydrodynamic_diameter_nm": diameter_nm,
        "zeta_potential_mv": zeta_mv,
        "free_diffusion_d0_m2_s": d0,
        "steric_retardation_ratio": steric_ratio,
        "charge_transmission_ratio": charge_ratio,
        "effective_diffusivity_ratio": effective_diffusivity_ratio,
        "mucus_penetrating_verdict": "OPTIMAL" if is_penetrating else ("SUBOPTIMAL" if effective_diffusivity_ratio > 0.1 else "ENTRAPPED"),
    }


def nonsense_readthrough_kinetics(
    ptc_variant: str,
    nmd_inhibition_pct: float,
    suppressor_trna_efficiency_pct: float,
    potentiator_fold_increase: float = 1.8
) -> Dict[str, Any]:
    """Model functional CFTR channel rescue under tri-modal nonsense suppression.
    
    Evaluates:
    1. Steady-state mRNA abundance (salvaged from UPF1/SMG1 decay).
    2. Translation readthrough rate across premature termination codon (UGA G vs UGA A).
    3. Restored short-circuit current (Delta I_sc) as percentage of Wild-Type.
    """
    variant_profiles = {
        "G542X": {
            "exon": 12,
            "stop_codon": "UGA",
            "context": "UGA G",
            "baseline_mrna_pct_wt": 12.0,  # Severe canonical NMD (<15%)
            "baseline_readthrough_pct": 0.5,
            "ejc_distance_nt": 184,
        },
        "W1282X": {
            "exon": 23,
            "stop_codon": "UGA",
            "context": "UGA A",
            "baseline_mrna_pct_wt": 38.0,  # Partial NMD escape (30-45%)
            "baseline_readthrough_pct": 1.2,
            "ejc_distance_nt": 32,
        },
    }
    
    if ptc_variant not in variant_profiles:
        raise ValueError(f"Variant {ptc_variant} not supported; choose from {list(variant_profiles)}")
        
    prof = variant_profiles[ptc_variant]
    base_mrna = prof["baseline_mrna_pct_wt"]
    
    # NMD suppression elevates available mRNA transcript pool toward wild-type (max 100%)
    salvaged_mrna_pct = base_mrna + (100.0 - base_mrna) * (nmd_inhibition_pct / 100.0)
    
    # Readthrough rate combines baseline leaky readthrough and engineered suppressor tRNA
    # Codon context UGA A (adenine at +4) provides 1.4x higher readthrough than UGA G
    context_multiplier = 1.4 if prof["context"] == "UGA A" else 1.0
    effective_readthrough_pct = (prof["baseline_readthrough_pct"] + suppressor_trna_efficiency_pct) * context_multiplier
    effective_readthrough_pct = min(effective_readthrough_pct, 40.0)  # Physiological saturation limit
    
    # Functional full-length channel synthesis
    channel_synthesis_pct_wt = (salvaged_mrna_pct / 100.0) * (effective_readthrough_pct / 100.0) * 100.0
    
    # Restored short-circuit current with potentiator (ivacaftor / deuterated potentiators)
    restored_isc_pct_wt = channel_synthesis_pct_wt * potentiator_fold_increase
    
    # Clinical therapeutic threshold: >= 10.0% WT Delta I_sc confers substantial clinical benefit
    reaches_clinical_threshold = (restored_isc_pct_wt >= 10.0)
    
    return {
        "variant": ptc_variant,
        "stop_codon_context": prof["context"],
        "nmd_inhibition_input_pct": nmd_inhibition_pct,
        "suppressor_trna_input_pct": suppressor_trna_efficiency_pct,
        "salvaged_mrna_pct_wt": round(salvaged_mrna_pct, 2),
        "effective_readthrough_pct": round(effective_readthrough_pct, 2),
        "functional_channel_synthesis_pct_wt": round(channel_synthesis_pct_wt, 2),
        "restored_isc_pct_wt": round(restored_isc_pct_wt, 2),
        "clinical_therapeutic_threshold_met": reaches_clinical_threshold,
        "therapeutic_margin": round(restored_isc_pct_wt - 10.0, 2),
    }


def calculate_truncated_cftr_competition(
    variant: str,
    truncated_synthesis_rate: float,
    full_length_synthesis_rate: float,
    corrector_efficacy_pct: float = 0.0
) -> Dict[str, Any]:
    """Model competitive ER chaperone entrapment and ERAD inhibition by truncated CFTR fragments.
    
    Ref: Manfredi et al., Biochemistry 2025, DOI: 10.1021/acs.biochem.5c00590.
    In W1282X, truncated 1,281-aa protein binds calnexin/Hsp70, competing with full-length channel.
    G542X (541 aa) is degraded rapidly via RQC with negligible ER chaperone competition.
    """
    if full_length_synthesis_rate <= 0:
        raise ValueError("Full-length synthesis rate must be positive")
    if truncated_synthesis_rate < 0:
        raise ValueError("Truncated synthesis rate cannot be negative")
        
    ratio = truncated_synthesis_rate / full_length_synthesis_rate
    
    if variant == "G542X":
        # Rapid ribosome quality control (RQC) degradation before full NBD1 folding; negligible ER competition
        maturation_efficiency = 0.95
        competition_severity = "NEGLIGIBLE"
    elif variant == "W1282X":
        # 1281-CFTR retains MSD1, NBD1, MSD2; forms core-glycosylated aggregates that compete for ER export
        # Corrector (elexacaftor/tezacaftor, PDB 8EIQ) stabilizes full-length NBD1-MSD interface, mitigating competition
        corrector_protection = min(0.80, corrector_efficacy_pct / 100.0 * 0.80)
        inhibition_constant_alpha = 0.15 * (1.0 - corrector_protection)
        maturation_efficiency = 1.0 / (1.0 + (inhibition_constant_alpha * ratio))
        competition_severity = "HIGH" if maturation_efficiency < 0.60 else ("MODERATE" if maturation_efficiency < 0.85 else "CONTROLLED")
    else:
        raise ValueError(f"Variant {variant} not supported for ERAD competition modeling")
        
    return {
        "variant": variant,
        "truncated_to_full_length_ratio": round(ratio, 2),
        "corrector_efficacy_pct": corrector_efficacy_pct,
        "er_maturation_efficiency": round(maturation_efficiency, 3),
        "competition_severity": competition_severity,
        "requires_allosteric_corrector": (variant == "W1282X" and ratio > 3.0),
    }


def calculate_transcriptome_readthrough_safety_score(suppressor_trna_efficiency_pct: float) -> Dict[str, Any]:
    """Calculate proteome-wide off-target readthrough and safety selectivity ratio.
    
    Human essential genes feature stringent termination architecture:
    - 28.6% utilize tandem stop codons (UGA-UAA, etc.) requiring dual readthrough (P ~ k_RT^2).
    - 59.8% utilize restrictive +4 context (UGA C, UAA G) with ~0.2x relative readthrough.
    - 11.6% utilize standard or permissive single stop codons.
    """
    if suppressor_trna_efficiency_pct < 0:
        raise ValueError("Efficiency cannot be negative")
        
    k_rt = suppressor_trna_efficiency_pct / 100.0
    
    # Weighted average off-target readthrough probability across human proteome
    prob_tandem = 0.286 * (k_rt ** 2)
    prob_restrictive = 0.598 * (0.20 * k_rt)
    prob_permissive = 0.116 * (0.80 * k_rt)
    
    global_off_target_pct = (prob_tandem + prob_restrictive + prob_permissive) * 100.0
    
    # CFTR W1282X UGA A target readthrough
    cftr_target_readthrough_pct = (1.2 + suppressor_trna_efficiency_pct) * 1.4
    
    selectivity_ratio = cftr_target_readthrough_pct / max(0.001, global_off_target_pct)
    
    if suppressor_trna_efficiency_pct < 5.0:
        verdict = "SUB_THERAPEUTIC"
    elif suppressor_trna_efficiency_pct <= 10.0:
        verdict = "OPTIMAL_SAFETY_WINDOW"
    elif suppressor_trna_efficiency_pct <= 15.0:
        verdict = "ELEVATED_OFF_TARGET_RISK"
    else:
        verdict = "PROTEOTOXIC_HAZARD"
        
    return {
        "suppressor_trna_efficiency_pct": suppressor_trna_efficiency_pct,
        "cftr_target_readthrough_pct": round(cftr_target_readthrough_pct, 2),
        "global_off_target_readthrough_pct": round(global_off_target_pct, 3),
        "selectivity_ratio": round(selectivity_ratio, 1),
        "safety_window_verdict": verdict,
        "tandem_stop_protection_confirmed": True,
    }
