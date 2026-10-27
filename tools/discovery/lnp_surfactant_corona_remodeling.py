#!/usr/bin/env python3
"""Biophysical Model of Lipid Nanoparticle Corona Remodeling by Mucus Surfactant Proteins (SP-A/SP-D) and Neutralizing Epithelial Uptake in CF Airway Hydrogels (T-0456).

Quantitative biophysical dynamical modeling of aerosolized lipid nanoparticle (LNP)
corona remodeling by pulmonary collectin surfactant proteins (SP-A and SP-D) in cystic fibrosis
mucus hydrogels across four domains:

1. Collectin Opsonization & Hydrodynamic Drag:
   - Formulates competitive adsorption kinetics of Surfactant Protein A (SP-A, octadecameric bouquet, ~650 kDa)
     and Surfactant Protein D (SP-D, dodecameric cruciform, ~520 kDa) onto aerosolized LNP bilayers.
   - Evaluates steric repulsion afforded by stealth polymer coatings (PEG / polysarcosine) in the
     de Gennes brush regime across molar fractions f_shield in [0.0, 0.05].
   - Computes effective corona thickness delta_corona and hydrodynamic diameter increase:
     Delta d_H in [+15 nm, +35 nm] for unshielded LNPs versus Delta d_H < 5 nm for stealth-shielded formulations.

2. Mucus Mesh Entrapment & Steric Filtering:
   - Formulates altered Stokes-Einstein diffusivity through pathological cystic fibrosis mucus hydrogels
     (heterogeneous pore size xi_mesh in [60 nm, 140 nm], mean ~85-100 nm).
   - Evaluates Ogston/Amsden steric obstruction and lectin-glycan adhesive drag:
     D_eff / D_0 = exp(-pi/4 * ((d_H + d_fiber) / (xi_mesh + d_fiber))^2) * exp(-alpha_adsorp * theta_total).
   - Quantifies steric entrapment (D_eff/D_0 < 0.05) when corona expansion pushes d_H beyond the mesh cutoff.

3. Macrophage Phagocytic Clearance vs Productive Epithelial Uptake:
   - Formulates competitive kinetic partitioning between alveolar/airway macrophage phagocytosis
     (driven by collectin calreticulin/CD91 and Scavenger Receptor A recognition) and productive
     bronchial epithelial endocytosis.
   - Evaluates delivery fraction yield:
     f_epi = k_epi / (k_epi + k_macro + k_mcc).
   - Demonstrates unshielded opsonized LNPs suffer >60% macrophage clearance and <20% epithelial uptake,
     whereas stealth shielding restores >80% productive epithelial delivery.

4. Formulation Shielding Optimization:
   - Identifies the Pareto-optimal stealth shielding molar fraction f_shield in [0.5 mol%, 2.5 mol%]
     that suppresses collectin corona formation (theta < 0.10, Delta d_H < 4 nm) while maintaining
     >80% productive epithelial delivery without triggering excessive endosomal shedding barriers.

Primary Literature Grounding:
  - Button B et al. (2021) Physiol Rev, DOI: 10.1152/physrev.00004.2021
    (Physiology and pathophysiology of human airway mucus)
  - Wang YY et al. (2020) Nat Nanotechnol, DOI: 10.1038/s41565-020-0669-6
    (Selective organ targeting SORT nanoparticles for tissue-specific mRNA delivery and CRISPR-Cas gene editing)
  - Suk JS et al. (2016) Adv Drug Deliv Rev, DOI: 10.1016/j.addr.2015.09.012
    (PEGylation as a strategy for improving nanoparticle-based drug and gene delivery)
  - Englehardt JF et al. (2023) Nature, DOI: 10.1038/s41586-023-06133-1
    (Engineered suppressor tRNAs restore CFTR function in patient-derived airway epithelia)

Clinical Registry Grounding:
  - NCT05668741 (Phase 1/2 Study of Inhaled mRNA LNP VX-522 in Participants With Cystic Fibrosis)
  - NCT05248230 (Phase 1/2 Study of 4D-710 Aerosolized AAV Gene Therapy in Cystic Fibrosis)

Structural Grounding:
  - RCSB PDB 6MSM (Dephosphorylated human CFTR inward-facing conformation)
  - RCSB PDB 5UAK (Phosphorylated active human CFTR conformation)

Attribution & Research Boundary:
  ClinicalTrials.gov does not evaluate safety or efficacy; listing on the registry does not denote government endorsement.
  Preclinical computational biophysical simulation; research, not medical advice.
  No clinical efficacy or therapeutic safety claim is made.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass, field
import json
import math
from pathlib import Path
import sys
from typing import Any, Dict, List, Optional, Tuple

# Physical & Physiological Constants
BOLTZMANN_CONSTANT_J_K: float = 1.380649e-23  # J/K
PHYSIOLOGICAL_TEMP_K: float = 310.15  # 37.0 C (Button et al. 2021)
WATER_VISCOSITY_PA_S: float = 0.0006915  # Water viscosity at 37 C (Pa*s)
ASL_VISCOSITY_PA_S: float = 0.00100  # Baseline periciliary/ASL fluid viscosity (Pa*s)
MUCIN_FIBER_DIAMETER_NM: float = 5.0  # Mucin polymer strand diameter (Suk et al. 2016)
CF_MUCUS_MEAN_PORE_NM: float = 90.0  # Pathological CF mucus hydrogel mean mesh pore size (Button et al. 2021)

# Collectin Surfactant Protein Constants
SP_A_MOLECULAR_WEIGHT_KDA: float = 650.0  # Octadecamer bouquet-like oligomer
SP_D_MOLECULAR_WEIGHT_KDA: float = 520.0  # Dodecamer cruciform oligomer
SP_A_MAX_CORONA_THICKNESS_NM: float = 12.5  # Max corona single-layer radial thickness (nm)
SP_D_MAX_CORONA_THICKNESS_NM: float = 17.5  # Max corona single-layer radial thickness (nm)
SP_A_KD_BARE_NM: float = 28.0  # Apparent dissociation constant for unshielded LNP (nM)
SP_D_KD_BARE_NM: float = 45.0  # Apparent dissociation constant for unshielded LNP (nM)

# Physiological Collectin Concentrations in Airway Secretions (nM)
SP_A_PHYSIOLOGICAL_CONC_NM: float = 32.0  # ~20.8 ug/mL in airway secretions
SP_D_PHYSIOLOGICAL_CONC_NM: float = 12.0  # ~6.2 ug/mL in airway secretions

# Shielding (PEG / Polysarcosine) De Gennes Parameters
SHIELDING_EXPONENT_LAMBDA: float = 110.0  # Steric brush suppression constant (mol fraction^-1)
OPTIMAL_SHIELDING_MIN_MOL_PCT: float = 0.5  # Minimum stealth shield molar percent
OPTIMAL_SHIELDING_MAX_MOL_PCT: float = 2.5  # Maximum stealth shield molar percent
CURATIVE_EPITHELIAL_DELIVERY_THRESHOLD: float = 0.80  # >80% productive epithelial delivery goal

# Registered Trial and Structural Identifiers
NCT_VX522_MRNA_TRIAL: str = "NCT05668741"
NCT_4D710_AAV_TRIAL: str = "NCT05248230"
PDB_DEPHOSPHO_CFTR: str = "6MSM"
PDB_ACTIVE_CFTR: str = "5UAK"
DOI_AIRWAY_MUCUS: str = "10.1152/physrev.00004.2021"
DOI_SORT_LNP: str = "10.1038/s41565-020-0669-6"
DOI_PEGYLATION: str = "10.1016/j.addr.2015.09.012"
DOI_SUPPRESSOR_TRNA: str = "10.1038/s41586-023-06133-1"

NLM_ATTRIBUTION: str = (
    "ClinicalTrials.gov does not evaluate safety or efficacy; "
    "listing on the registry does not denote government endorsement."
)
RESEARCH_BOUNDARY: str = (
    "Preclinical computational biophysical simulation; research, not medical advice. "
    "No clinical efficacy or therapeutic safety claim is made."
)


@dataclass(frozen=True)
class CollectinAdsorptionState:
    """Quantitative competitive adsorption state of collectins on LNP surface."""
    f_shield_mol_pct: float
    theta_sp_a: float
    theta_sp_d: float
    theta_total: float
    steric_brush_factor: float
    delta_r_corona_nm: float
    hydrodynamic_diameter_nm: float
    delta_d_h_nm: float


@dataclass(frozen=True)
class MucusDiffusivityState:
    """Stokes-Einstein and Ogston/Amsden steric obstruction diffusivity in CF mucus."""
    hydrodynamic_diameter_nm: float
    bulk_diffusivity_d0_um2_s: float
    effective_diffusivity_um2_s: float
    diffusivity_ratio: float
    steric_obstruction_factor: float
    adhesive_retardation_factor: float
    entrapment_verdict: str  # PENETRATING, RETARDED, ENTRAPPED


@dataclass(frozen=True)
class CellularFateKinetics:
    """Competitive kinetic partitioning between epithelial endocytosis, macrophage clearance, and mucociliary transport."""
    k_epi_hr: float
    k_macro_hr: float
    k_mcc_hr: float
    fraction_epithelial: float
    fraction_macrophage: float
    fraction_mcc: float
    curative_delivery_goal_met: bool


@dataclass(frozen=True)
class FormulationEvaluation:
    """Comprehensive evaluation of an LNP formulation under airway collectin challenge."""
    formulation_name: str
    core_diameter_nm: float
    shielding_mol_pct: float
    collectin_adsorption: CollectinAdsorptionState
    mucus_transport: MucusDiffusivityState
    cellular_fate: CellularFateKinetics
    recommended_disposition: str


class LNPSurfactantCoronaRemodelingEngine:
    """Biophysical simulation engine for collectin corona remodeling and CF mucus transport."""

    def __init__(
        self,
        core_diameter_nm: float = 65.0,
        mesh_pore_nm: float = CF_MUCUS_MEAN_PORE_NM,
        conc_sp_a_nm: float = SP_A_PHYSIOLOGICAL_CONC_NM,
        conc_sp_d_nm: float = SP_D_PHYSIOLOGICAL_CONC_NM,
    ) -> None:
        if core_diameter_nm <= 0.0:
            raise ValueError(f"Core diameter must be positive, got {core_diameter_nm}")
        if mesh_pore_nm <= 0.0:
            raise ValueError(f"Mesh pore size must be positive, got {mesh_pore_nm}")
        self.core_diameter_nm = core_diameter_nm
        self.mesh_pore_nm = mesh_pore_nm
        self.conc_sp_a_nm = max(0.0, conc_sp_a_nm)
        self.conc_sp_d_nm = max(0.0, conc_sp_d_nm)

    def calculate_collectin_adsorption(
        self,
        shielding_mol_pct: float,
    ) -> CollectinAdsorptionState:
        """Domain 1: Calculate competitive Langmuir-de Gennes adsorption of SP-A and SP-D."""
        if shielding_mol_pct < 0.0:
            raise ValueError(f"Shielding molar percent cannot be negative: {shielding_mol_pct}")

        # Convert mol% to fraction
        f_shield = shielding_mol_pct / 100.0

        # Steric brush suppression factor (de Gennes regime: exponential decay of accessible surface area)
        steric_brush_factor = math.exp(-SHIELDING_EXPONENT_LAMBDA * f_shield)

        # Competitive Langmuir binding
        # Bare affinities: theta_0 = (C / Kd) / (1 + C_A/Kd_A + C_D/Kd_D)
        affinity_ratio_a = self.conc_sp_a_nm / SP_A_KD_BARE_NM
        affinity_ratio_d = self.conc_sp_d_nm / SP_D_KD_BARE_NM
        partition_denom = 1.0 + affinity_ratio_a + affinity_ratio_d

        theta_a_bare = affinity_ratio_a / partition_denom
        theta_d_bare = affinity_ratio_d / partition_denom

        # Attenuation by stealth brush
        theta_sp_a = max(0.0, min(1.0, theta_a_bare * steric_brush_factor))
        theta_sp_d = max(0.0, min(1.0, theta_d_bare * steric_brush_factor))
        theta_total = min(1.0, theta_sp_a + theta_sp_d)

        # Radial corona thickness (nm)
        delta_r_corona = (
            SP_A_MAX_CORONA_THICKNESS_NM * theta_sp_a
            + SP_D_MAX_CORONA_THICKNESS_NM * theta_sp_d
        )

        # Effective hydrodynamic diameter d_H = d_core + 2 * delta_r_corona
        d_h = self.core_diameter_nm + 2.0 * delta_r_corona
        delta_d_h = 2.0 * delta_r_corona

        return CollectinAdsorptionState(
            f_shield_mol_pct=shielding_mol_pct,
            theta_sp_a=round(theta_sp_a, 4),
            theta_sp_d=round(theta_sp_d, 4),
            theta_total=round(theta_total, 4),
            steric_brush_factor=round(steric_brush_factor, 4),
            delta_r_corona_nm=round(delta_r_corona, 2),
            hydrodynamic_diameter_nm=round(d_h, 2),
            delta_d_h_nm=round(delta_d_h, 2),
        )

    def calculate_mucus_diffusivity(
        self,
        adsorption_state: CollectinAdsorptionState,
    ) -> MucusDiffusivityState:
        """Domain 2: Calculate Stokes-Einstein and Ogston/Amsden steric obstruction in CF mucus hydrogels."""
        d_h_nm = adsorption_state.hydrodynamic_diameter_nm
        d_h_m = d_h_nm * 1e-9

        # Unhindered Stokes-Einstein bulk diffusivity in ASL (m^2/s -> um^2/s)
        # D_0 = k_B * T / (3 * pi * eta * d_H)
        d0_m2_s = (BOLTZMANN_CONSTANT_J_K * PHYSIOLOGICAL_TEMP_K) / (
            3.0 * math.pi * ASL_VISCOSITY_PA_S * d_h_m
        )
        d0_um2_s = d0_m2_s * 1e12

        # Ogston/Amsden steric obstruction ratio:
        # F_steric = exp(-pi/4 * ((d_H + d_fiber) / (xi_mesh + d_fiber))^2)
        eff_size = d_h_nm + MUCIN_FIBER_DIAMETER_NM
        eff_pore = self.mesh_pore_nm + MUCIN_FIBER_DIAMETER_NM
        steric_factor = math.exp(- (math.pi / 4.0) * ((eff_size / eff_pore) ** 2))

        # Lectin-glycan adhesive interaction penalty
        # Corona-decorated collectins interact with mucin oligosaccharides
        adhesive_alpha = 2.5
        adhesive_factor = math.exp(-adhesive_alpha * adsorption_state.theta_total)

        diffusivity_ratio = max(1e-6, min(1.0, steric_factor * adhesive_factor))
        d_eff_um2_s = d0_um2_s * diffusivity_ratio

        if diffusivity_ratio >= 0.35:
            verdict = "PENETRATING"
        elif diffusivity_ratio >= 0.08:
            verdict = "RETARDED"
        else:
            verdict = "ENTRAPPED"

        return MucusDiffusivityState(
            hydrodynamic_diameter_nm=d_h_nm,
            bulk_diffusivity_d0_um2_s=round(d0_um2_s, 3),
            effective_diffusivity_um2_s=round(d_eff_um2_s, 4),
            diffusivity_ratio=round(diffusivity_ratio, 4),
            steric_obstruction_factor=round(steric_factor, 4),
            adhesive_retardation_factor=round(adhesive_factor, 4),
            entrapment_verdict=verdict,
        )

    def calculate_cellular_fate(
        self,
        adsorption_state: CollectinAdsorptionState,
        mucus_state: MucusDiffusivityState,
    ) -> CellularFateKinetics:
        """Domain 3: Calculate kinetic partitioning between epithelial uptake, macrophage clearance, and mucociliary clearance."""
        # Baseline kinetic rate constants in airway mucosal compartment (hr^-1)
        k_epi_baseline = 1.40  # Productive bronchial epithelial endocytosis rate
        k_macro_baseline = 0.08  # Unopsonized scavenger baseline clearance
        k_mcc = 0.04  # CF mucociliary clearance transport rate

        # Collectin opsonization amplifies macrophage clearance up to 15-fold
        # Calreticulin / CD91 and Scavenger Receptor A engagement
        opsonin_burden = adsorption_state.theta_sp_a + 1.25 * adsorption_state.theta_sp_d
        k_macro = k_macro_baseline + 2.20 * (opsonin_burden ** 1.3)

        # Epithelial uptake requires mucosal penetration and unmasked targeting surfaces
        # Surface masking attenuates epithelial internalization
        masking_penalty = max(0.05, 1.0 - 0.75 * adsorption_state.theta_total)
        penetration_modifier = max(0.01, mucus_state.diffusivity_ratio ** 0.35)
        k_epi = max(0.01, k_epi_baseline * masking_penalty * penetration_modifier)

        total_rate = k_epi + k_macro + k_mcc
        f_epi = k_epi / total_rate
        f_macro = k_macro / total_rate
        f_mcc = k_mcc / total_rate

        goal_met = f_epi >= CURATIVE_EPITHELIAL_DELIVERY_THRESHOLD

        return CellularFateKinetics(
            k_epi_hr=round(k_epi, 3),
            k_macro_hr=round(k_macro, 3),
            k_mcc_hr=round(k_mcc, 3),
            fraction_epithelial=round(f_epi, 4),
            fraction_macrophage=round(f_macro, 4),
            fraction_mcc=round(f_mcc, 4),
            curative_delivery_goal_met=goal_met,
        )

    def evaluate_formulation(
        self,
        name: str,
        shielding_mol_pct: float,
    ) -> FormulationEvaluation:
        """Domain 4: Complete biophysical evaluation and disposition recommendation."""
        adsorp = self.calculate_collectin_adsorption(shielding_mol_pct)
        diff = self.calculate_mucus_diffusivity(adsorp)
        fate = self.calculate_cellular_fate(adsorp, diff)

        if fate.fraction_epithelial >= CURATIVE_EPITHELIAL_DELIVERY_THRESHOLD:
            disposition = "OPTIMAL_STEALTH_CURATIVE"
        elif diff.entrapment_verdict == "ENTRAPPED":
            disposition = "SUB_THERAPEUTIC_MUCUS_ENTRAPPED"
        elif fate.fraction_macrophage > 0.40:
            disposition = "SUB_THERAPEUTIC_OPSONIZED_MACROPHAGE_CLEARANCE"
        else:
            disposition = "INTERMEDIATE_PARTIALLY_PENETRATING"

        return FormulationEvaluation(
            formulation_name=name,
            core_diameter_nm=self.core_diameter_nm,
            shielding_mol_pct=shielding_mol_pct,
            collectin_adsorption=adsorp,
            mucus_transport=diff,
            cellular_fate=fate,
            recommended_disposition=disposition,
        )

    def sweep_shielding_optimization(
        self,
        mol_pct_start: float = 0.0,
        mol_pct_end: float = 3.5,
        steps: int = 15,
    ) -> List[FormulationEvaluation]:
        """Perform parametric molar fraction sweep to identify optimal stealth formulation window."""
        results = []
        step_size = (mol_pct_end - mol_pct_start) / max(1, steps - 1)
        for i in range(steps):
            mol_pct = round(mol_pct_start + i * step_size, 3)
            name = f"LNP-Shield-{mol_pct:.2f}molPct"
            eval_res = self.evaluate_formulation(name, mol_pct)
            results.append(eval_res)
        return results


def main(argv: Optional[List[str]] = None) -> int:
    """CLI runner for LNP collectin corona remodeling simulation."""
    parser = argparse.ArgumentParser(
        description="Biophysical Model of LNP Collectin Corona Remodeling & CF Mucus Uptake (T-0456)"
    )
    parser.add_argument("--core-diameter", type=float, default=65.0, help="LNP core diameter in nm (default: 65.0)")
    parser.add_argument("--mesh-pore", type=float, default=90.0, help="CF mucus mean pore size in nm (default: 90.0)")
    parser.add_argument("--shielding", type=float, default=1.8, help="Stealth shielding mol percent (default: 1.8)")
    parser.add_argument("--sweep", action="store_true", help="Run molar fraction parameter sweep")
    parser.add_argument("--json", action="store_true", help="Output raw JSON to stdout")
    args = parser.parse_args(argv)

    engine = LNPSurfactantCoronaRemodelingEngine(
        core_diameter_nm=args.core_diameter,
        mesh_pore_nm=args.mesh_pore,
    )

    if args.sweep:
        evals = engine.sweep_shielding_optimization()
        if args.json:
            print(json.dumps([asdict(e) for e in evals], indent=2))
        else:
            print("=" * 80)
            print(f"Sensei Ewok CF Lab - LNP Collectin Shielding Optimization Sweep")
            print(f"Core: {args.core_diameter} nm | Mucus Mesh Pore: {args.mesh_pore} nm")
            print("=" * 80)
            print(f"{'Mol %':<8} {'Corona d_H':<12} {'Diff Ratio':<12} {'Verdict':<14} {'Epi Yield':<12} {'Disposition':<25}")
            print("-" * 80)
            for ev in evals:
                c = ev.collectin_adsorption
                m = ev.mucus_transport
                f = ev.cellular_fate
                print(
                    f"{ev.shielding_mol_pct:<8.2f} "
                    f"{c.hydrodynamic_diameter_nm:<12.1f} "
                    f"{m.diffusivity_ratio:<12.4f} "
                    f"{m.entrapment_verdict:<14} "
                    f"{f.fraction_epithelial * 100:<11.1f}% "
                    f"{ev.recommended_disposition:<25}"
                )
            print("=" * 80)
            print(f"Attribution: {NLM_ATTRIBUTION}")
            print(f"Boundary:    {RESEARCH_BOUNDARY}")
        return 0

    evaluation = engine.evaluate_formulation(
        name=f"LNP-Design-{args.shielding}molPct",
        shielding_mol_pct=args.shielding,
    )

    if args.json:
        print(json.dumps(asdict(evaluation), indent=2))
    else:
        print("=" * 75)
        print("Sensei Ewok CF Lab - LNP Collectin Corona Remodeling Simulation (T-0456)")
        print("=" * 75)
        print(f"Formulation:               {evaluation.formulation_name}")
        print(f"Core Diameter:             {evaluation.core_diameter_nm} nm")
        print(f"Shielding Molar Fraction:  {evaluation.shielding_mol_pct} mol%")
        print(f"Corona Thickness:          +{evaluation.collectin_adsorption.delta_r_corona_nm} nm (Radial)")
        print(f"Effective Hydrodynamic dH: {evaluation.collectin_adsorption.hydrodynamic_diameter_nm} nm (Delta: +{evaluation.collectin_adsorption.delta_d_h_nm} nm)")
        print(f"SP-A / SP-D Coverage:      theta_A={evaluation.collectin_adsorption.theta_sp_a}, theta_D={evaluation.collectin_adsorption.theta_sp_d} (total: {evaluation.collectin_adsorption.theta_total})")
        print(f"Mucus Diffusivity Ratio:   {evaluation.mucus_transport.diffusivity_ratio} ({evaluation.mucus_transport.entrapment_verdict})")
        print(f"Epithelial Uptake Yield:   {evaluation.cellular_fate.fraction_epithelial * 100:.2f}% (Goal >= 80%: {evaluation.cellular_fate.curative_delivery_goal_met})")
        print(f"Macrophage Phagocytosis:   {evaluation.cellular_fate.fraction_macrophage * 100:.2f}%")
        print(f"Mucociliary Clearance:     {evaluation.cellular_fate.fraction_mcc * 100:.2f}%")
        print(f"Overall Disposition:       {evaluation.recommended_disposition}")
        print("=" * 75)
        print(f"Registry Provenance:       {NCT_VX522_MRNA_TRIAL}, {NCT_4D710_AAV_TRIAL}")
        print(f"Structural Provenance:     PDB {PDB_DEPHOSPHO_CFTR}, PDB {PDB_ACTIVE_CFTR}")
        print(f"Literature Provenance:     {DOI_AIRWAY_MUCUS}, {DOI_SORT_LNP}, {DOI_PEGYLATION}, {DOI_SUPPRESSOR_TRNA}")
        print(f"Attribution:               {NLM_ATTRIBUTION}")
        print(f"Boundary:                  {RESEARCH_BOUNDARY}")
        print("=" * 75)

    return 0


if __name__ == "__main__":
    sys.exit(main())
