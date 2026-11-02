#!/usr/bin/env python3
"""Biophysical Modeling of Macropinocytosis vs Clathrin-Mediated Endocytosis Flux in Airway Stem Cells (T-0561).

Quantitative biophysical model of coupled dual-pathway internalization flux,
endosomal sorting dynamics, and cytosolic mRNA transfection yield in airway
epithelial basal stem cells under inhaled lipid nanoparticle (LNP) delivery regimes.

Mathematical & Biophysical Formulation:
1. Coupled Dual-Pathway Cellular Internalization Flux:
   J_internalize = J_CME + J_macro
   where:
     J_CME   = V_max_CME * (C_apical / (K_m_CME + C_apical))   [Dynamin/clathrin-mediated endocytosis]
     J_macro = k_macro * C_apical * Phi_ruffle                  [Actin-driven macropinocytosis]
     C_apical   = Apical LNP concentration (nM)
     Phi_ruffle = Basal stem cell membrane ruffling index (dimensionless, 0.0 to 1.0)
     V_max_CME  = Maximum CME uptake velocity (particles / (cell * min))
     K_m_CME    = Half-saturation Michaelis constant for CME receptors (nM)
     k_macro    = Second-order macropinocytic rate constant (particles / (cell * min * nM))

2. Endosomal Sorting and Recycling Kinetics:
   dC_early/dt = J_internalize - (k_recycle + k_sort) * C_early
   where:
     f_recycle = k_recycle / (k_recycle + k_sort)              [Apical recycling fraction back to lumen]
     f_sort    = k_sort / (k_recycle + k_sort) = 1 - f_recycle  [Productive late endosomal maturation]

3. Stem Cell Transfection Yield Metric:
   Y_stem = int_0^t_delivery [ J_internalize(t) * eta_escape * (1 - f_recycle) ] dt
   For constant delivery window t_delivery:
     Y_stem = J_internalize * eta_escape * (1 - f_recycle) * t_delivery
   where:
     eta_escape = Fractional endosomal escape efficiency (typically 1.0 - 5.0%)
     t_delivery = Effective aerosol deposition / delivery duration (minutes)
     Y_stem     = Total cytosolic functional mRNA transcripts delivered per stem cell

Primary Evidence Grounding:
- DOI: 10.1038/s41565-020-0669-6 (Selective organ targeting (SORT) nanoparticles for tissue-specific mRNA delivery and CRISPR-Cas gene editing)
- DOI: 10.1038/s41587-023-01861-1 (Inhaled delivery of lipid nanoparticles to respiratory epithelia)
- DOI: 10.1152/physrev.00004.2021 (Physiology and pathophysiology of human airway mucus)

Clinical Trials:
- NCT05668741 (A Phase 1/2 Study of VX-522 in Participants With Cystic Fibrosis)

Mandatory Attribution & Disclaimers:
ClinicalTrials.gov data are cited here with attribution to the U.S. National Library of Medicine (NLM):
"ClinicalTrials.gov does not evaluate safety or efficacy; listing on the registry does not denote government endorsement."
Research boundary: Preclinical computational model; research, not medical advice. No clinical efficacy or therapeutic safety claim is made.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
import math
from pathlib import Path
import sys
from typing import Any, Dict, List, Optional, Tuple

NLM_REGISTRY_DISCLAIMER = (
    "ClinicalTrials.gov does not evaluate safety or efficacy; "
    "listing on the registry does not denote government endorsement."
)
RESEARCH_BOUNDARY = (
    "Preclinical in-silico biophysical modeling of basal stem cell endocytic internalization "
    "and cytosolic mRNA delivery yield; research, not medical advice. "
    "No clinical efficacy or therapeutic safety claim is made."
)

CLINICAL_TRIAL_VX522 = "NCT05668741"
PRIMARY_DOI_SORT_LNP = "10.1038/s41565-020-0669-6"
PRIMARY_DOI_INHALED_LNP = "10.1038/s41587-023-01861-1"
PRIMARY_DOI_MUCUS_PHYSIOLOGY = "10.1152/physrev.00004.2021"


@dataclass(frozen=True)
class EndocyticPathwayParameters:
    """Biophysical endocytic flux and sorting kinetics parameters."""
    profile_name: str = "Basal_Stem_SORT_Optimized"
    apical_lnp_concentration_nm: float = 30.0   # C_apical: Local apical LNP concentration (nM)
    v_max_cme_particles_min: float = 120.0       # V_max_CME: Saturated clathrin flux (particles/cell/min)
    km_cme_nm: float = 25.0                     # K_m_CME: CME half-saturation constant (nM)
    k_macro_per_nm_min: float = 3.5             # k_macro: Macropinocytic rate constant (particles/cell/min/nM)
    membrane_ruffling_index: float = 0.85        # Phi_ruffle: Dynamic actin ruffling factor (0.0 to 1.0)
    k_sort_per_min: float = 0.040               # k_sort: Rate of productive late endosome sorting (1/min)
    k_recycle_per_min: float = 0.020            # k_recycle: Rate of apical exocytic recycling (1/min)
    endosomal_escape_efficiency: float = 0.045  # eta_escape: Cytosolic escape fraction (e.g. 4.5% for SORT)
    delivery_duration_min: float = 120.0        # t_delivery: Exposure window duration (min)


@dataclass
class EndocyticPathwayProfile:
    """Quantitative internalization partitioning and stem cell transfection output."""
    profile_name: str
    apical_concentration_nm: float
    j_cme_particles_per_min: float
    j_macro_particles_per_min: float
    total_internalization_flux_per_min: float
    cme_flux_fraction: float
    macro_flux_fraction: float
    recycling_fraction: float
    endosomal_escape_efficiency: float
    cytosolic_yield_mrna_per_cell: float
    transfection_status: str


def compute_clathrin_mediated_flux(
    concentration_nm: float,
    v_max: float = 120.0,
    km_nm: float = 25.0,
) -> float:
    """Compute saturating clathrin-mediated endocytosis flux J_CME = V_max * C / (Km + C)."""
    c = max(concentration_nm, 0.0)
    if c == 0.0 or v_max == 0.0:
        return 0.0
    return v_max * (c / (km_nm + c))


def compute_macropinocytic_flux(
    concentration_nm: float,
    k_macro: float = 3.5,
    ruffling_index: float = 0.85,
) -> float:
    """Compute actin-driven macropinocytosis flux J_macro = k_macro * C * Phi_ruffle."""
    c = max(concentration_nm, 0.0)
    phi = min(max(ruffling_index, 0.0), 1.0)
    return k_macro * c * phi


def compute_recycling_fraction(k_sort: float, k_recycle: float) -> float:
    """Compute fraction of internalized endosomes recycled to lumen: f_recycle = k_recycle / (k_sort + k_recycle)."""
    total = k_sort + k_recycle
    if total <= 0.0:
        return 0.0
    return k_recycle / total


def simulate_endocytic_partitioning(params: EndocyticPathwayParameters) -> EndocyticPathwayProfile:
    """Simulate coupled CME vs macropinocytosis flux and stem cell cytosolic mRNA delivery yield."""
    c = params.apical_lnp_concentration_nm
    j_cme = compute_clathrin_mediated_flux(c, params.v_max_cme_particles_min, params.km_cme_nm)
    j_macro = compute_macropinocytic_flux(c, params.k_macro_per_nm_min, params.membrane_ruffling_index)
    j_total = j_cme + j_macro

    if j_total > 0.0:
        f_cme = j_cme / j_total
        f_macro = j_macro / j_total
    else:
        f_cme = 0.0
        f_macro = 0.0

    f_recycle = compute_recycling_fraction(params.k_sort_per_min, params.k_recycle_per_min)
    f_productive = max(1.0 - f_recycle, 0.0)

    # Stem cell transfection yield: Y_stem = J_total * eta_escape * (1 - f_recycle) * t_delivery
    y_stem = j_total * params.endosomal_escape_efficiency * f_productive * params.delivery_duration_min

    if y_stem >= 500.0:
        status = "ROBUST_STEM_CELL_TRANSFECTION"
    elif y_stem >= 100.0:
        status = "MODERATE_PARTIAL_TRANSFECTION"
    else:
        status = "SUB_THERAPEUTIC_DELIVERY_FAILURE"

    return EndocyticPathwayProfile(
        profile_name=params.profile_name,
        apical_concentration_nm=round(c, 2),
        j_cme_particles_per_min=round(j_cme, 2),
        j_macro_particles_per_min=round(j_macro, 2),
        total_internalization_flux_per_min=round(j_total, 2),
        cme_flux_fraction=round(f_cme, 3),
        macro_flux_fraction=round(f_macro, 3),
        recycling_fraction=round(f_recycle, 3),
        endosomal_escape_efficiency=round(params.endosomal_escape_efficiency, 4),
        cytosolic_yield_mrna_per_cell=round(y_stem, 1),
        transfection_status=status,
    )


def evaluate_endocytic_pathway_partitioning() -> Dict[str, Any]:
    """Evaluate endocytic partitioning across SORT-optimized, standard LNP, and quiescent stem cell states."""
    profiles = {
        "Basal_Stem_SORT_Optimized": EndocyticPathwayParameters(
            profile_name="Basal_Stem_SORT_Optimized",
            apical_lnp_concentration_nm=30.0,
            membrane_ruffling_index=0.85,
            endosomal_escape_efficiency=0.045,
            delivery_duration_min=120.0,
        ),
        "Basal_Stem_Standard_LNP": EndocyticPathwayParameters(
            profile_name="Basal_Stem_Standard_LNP",
            apical_lnp_concentration_nm=30.0,
            membrane_ruffling_index=0.45,
            endosomal_escape_efficiency=0.015,
            delivery_duration_min=120.0,
        ),
        "Quiescent_Stem_Untargeted": EndocyticPathwayParameters(
            profile_name="Quiescent_Stem_Untargeted",
            apical_lnp_concentration_nm=10.0,
            membrane_ruffling_index=0.15,
            endosomal_escape_efficiency=0.010,
            delivery_duration_min=60.0,
        ),
        "Clathrin_Inhibited_Rescue": EndocyticPathwayParameters(
            profile_name="Clathrin_Inhibited_Rescue",
            apical_lnp_concentration_nm=30.0,
            v_max_cme_particles_min=0.0,  # Dynamin / CME completely blocked
            membrane_ruffling_index=0.90,
            endosomal_escape_efficiency=0.045,
            delivery_duration_min=120.0,
        ),
    }

    results: Dict[str, Dict[str, Any]] = {}
    for name, p in profiles.items():
        res = simulate_endocytic_partitioning(p)
        results[name] = asdict(res)

    sort_opt = results["Basal_Stem_SORT_Optimized"]
    std_lnp = results["Basal_Stem_Standard_LNP"]
    quiescent = results["Quiescent_Stem_Untargeted"]
    cme_inhib = results["Clathrin_Inhibited_Rescue"]

    # Invariants verification
    sort_robust = (
        sort_opt["cytosolic_yield_mrna_per_cell"] >= 500.0
        and sort_opt["transfection_status"] == "ROBUST_STEM_CELL_TRANSFECTION"
    )
    std_moderate = (
        100.0 <= std_lnp["cytosolic_yield_mrna_per_cell"] < 500.0
        and std_lnp["transfection_status"] == "MODERATE_PARTIAL_TRANSFECTION"
    )
    quiescent_sub = (
        quiescent["cytosolic_yield_mrna_per_cell"] < 100.0
        and quiescent["transfection_status"] == "SUB_THERAPEUTIC_DELIVERY_FAILURE"
    )
    cme_inhib_macro_only = (
        cme_inhib["j_cme_particles_per_min"] == 0.0
        and cme_inhib["macro_flux_fraction"] == 1.0
        and cme_inhib["cytosolic_yield_mrna_per_cell"] >= 300.0
    )

    verified = sort_robust and std_moderate and quiescent_sub and cme_inhib_macro_only

    return {
        "task_id": "T-0561",
        "primary_evidence": {
            "dois": [
                PRIMARY_DOI_SORT_LNP,
                PRIMARY_DOI_INHALED_LNP,
                PRIMARY_DOI_MUCUS_PHYSIOLOGY,
            ],
            "clinical_trials": [CLINICAL_TRIAL_VX522],
        },
        "regulatory_disclaimer": NLM_REGISTRY_DISCLAIMER,
        "research_boundary": RESEARCH_BOUNDARY,
        "profiles": results,
        "sort_optimized_yield_mrna_per_cell": sort_opt["cytosolic_yield_mrna_per_cell"],
        "standard_lnp_yield_mrna_per_cell": std_lnp["cytosolic_yield_mrna_per_cell"],
        "endocytic_partitioning_verified": verified,
    }


def falsify_endocytic_pathway_partitioning() -> Dict[str, Any]:
    """AH-10 dynamic negative-control falsification check."""
    # Negative control: complete internalization shutdown (0 nM apical LNP + 0 membrane ruffling)
    neg_params = EndocyticPathwayParameters(
        profile_name="Poisoned_Negative_Control",
        apical_lnp_concentration_nm=0.0,
        membrane_ruffling_index=0.0,
    )
    neg_profile = simulate_endocytic_partitioning(neg_params)
    passed = (
        neg_profile.total_internalization_flux_per_min == 0.0
        and neg_profile.cytosolic_yield_mrna_per_cell == 0.0
        and neg_profile.transfection_status == "SUB_THERAPEUTIC_DELIVERY_FAILURE"
    )

    return {
        "rule": "AH-10",
        "task_id": "T-0561",
        "test_name": "falsify_endocytic_pathway_partitioning",
        "negative_control_passed": passed,
        "internalization_flux": neg_profile.total_internalization_flux_per_min,
        "cytosolic_yield": neg_profile.cytosolic_yield_mrna_per_cell,
        "verdict": "PASS" if passed else "FAIL",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Endocytic pathway partitioning simulation (T-0561)")
    parser.add_argument("--json", action="store_true", help="Emit raw JSON")
    parser.add_argument("--falsify", action="store_true", help="Run AH-10 negative-control falsification")
    args = parser.parse_args()

    if args.falsify:
        res = falsify_endocytic_pathway_partitioning()
        print(json.dumps(res, indent=2))
        return 0 if res["negative_control_passed"] else 1

    report = evaluate_endocytic_pathway_partitioning()
    if args.json:
        print(json.dumps(report, indent=2))
    else:
        print(f"Endocytic Pathway Partitioning (T-0561): {'PASS' if report['endocytic_partitioning_verified'] else 'FAIL'}")
        print(f"  SORT-Optimized Transfection Yield: {report['sort_optimized_yield_mrna_per_cell']} mRNA/cell")
        print(f"  Standard LNP Transfection Yield:   {report['standard_lnp_yield_mrna_per_cell']} mRNA/cell")
    return 0 if report["endocytic_partitioning_verified"] else 1


if __name__ == "__main__":
    sys.exit(main())
