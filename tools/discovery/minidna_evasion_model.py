#!/usr/bin/env python3
"""Nonviral Minimal DNA Vector (minidnacftr) CpG Depletion & TLR9 Evasion Model (T-0233).

Models CpG dinucleotide density (rho_CpG), endosomal Toll-Like Receptor 9 (TLR9)
pattern-recognition activation flux (J_TLR9), cytokine reactogenicity index
(R_cytokine), cytoplasmic/nuclear transport ODE dynamics, H3K9me3 heterochromatin
silencing, and Pol III anticodon-edited suppressor tRNA (ACE-tRNA) expression
persistence in human bronchial epithelia.

Primary Evidence Grounding:
  - Porter, Ko, Sorensen et al. (2026), Nucleic Acids Res, DOI: 10.1093/nar/gkag082 (PMID: 41641696)
  - Al-Nuaimi et al. (2023), Nat Commun, DOI: 10.1038/s41467-023-41392-8
  - ClinicalTrials.gov ID: NCT05668741 (inhaled nonviral CFTR delivery clinical surveillance)
  - NCBI RefSeq NM_000492.4 (CFTR mRNA reference)

Research boundary: Preclinical computational model; research, not medical advice.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

SCRIPT_PATH = Path(__file__).resolve()

NLM_ATTRIBUTION = (
    "This product uses publicly available data from the U.S. National Library of Medicine (NLM), "
    "National Institutes of Health, Department of Health and Human Services; NLM is not responsible "
    "for the product and does not endorse or recommend this or any other product."
)
RESEARCH_BOUNDARY = (
    "Preclinical in-silico minimal DNA vector & TLR9 evasion simulation; research, not medical advice. "
    "No clinical efficacy or safety claim is made."
)

# Default kinetic constants grounded in T-0233 spec
K_TLR9_TURNOVER_PER_HR = 0.42      # Catalytic TLR9 turnover rate per CpG motif (hr^-1)
KM_ENDO_NM = 35.0                  # Endosomal TLR9 Michaelis constant K_m (nM)
J_SAFE_FLUX_THRESHOLD = 4.5        # Clinical tolerability threshold J_safe (flux units)
K_ESC_ENDO_PER_HR = 0.085          # Endosomal escape rate constant k_esc (hr^-1)
K_DEG_ENDO_PER_HR = 0.140          # Endosomal lysosomal nuclease degradation k_deg,endo (hr^-1)
K_DEG_CYT_PER_HR = 0.032           # Cytosolic TREX1/DNase II degradation k_deg,cyt (hr^-1)
K_NUC_IMPORT_PER_HR = 0.048        # Nuclear pore import rate k_nuc (hr^-1, enhanced by DTS/Pol III TFIIIC)
K_SILENCE_BASE_PER_DAY = 0.0025    # Baseline episomal chromatin silencing rate (day^-1)
ALPHA_CPG_SILENCE_PER_DAY = 1.85   # H3K9me3 CpG-dependent heterochromatin silencing coefficient (bp * day^-1)
MU_DIV_AIRWAY_PER_DAY = 0.0055     # Mitotic dilution rate in differentiated airway epithelium (~125-day half-life)


@dataclass(frozen=True)
class VectorArchitectureSpec:
    """Biophysical specification of a nonviral DNA vector architecture."""

    vector_id: str
    label: str
    length_bp: int
    cpg_count: int
    ace_trna_copy_number: int
    nuclear_import_multiplier: float
    bacterial_backbone_present: bool
    citation_doi: str


VECTOR_CATALOG: dict[str, VectorArchitectureSpec] = {
    "bacterial_plasmid": VectorArchitectureSpec(
        vector_id="bacterial_plasmid",
        label="Conventional Bacterial Plasmid (pDNA 5.0 kb, ColE1 + KanR)",
        length_bp=5000,
        cpg_count=310,
        ace_trna_copy_number=1,
        nuclear_import_multiplier=0.55,
        bacterial_backbone_present=True,
        citation_doi="10.1093/nar/gkag082",
    ),
    "minicircle_v1": VectorArchitectureSpec(
        vector_id="minicircle_v1",
        label="First-Generation Minicircle DNA (mcDNA 1.8 kb)",
        length_bp=1800,
        cpg_count=28,
        ace_trna_copy_number=2,
        nuclear_import_multiplier=1.00,
        bacterial_backbone_present=False,
        citation_doi="10.1093/nar/gkag082",
    ),
    "nanoplasmid_r6k": VectorArchitectureSpec(
        vector_id="nanoplasmid_r6k",
        label="R6K-Origin Nanoplasmid (1.45 kb, RNA-OUT selection)",
        length_bp=1450,
        cpg_count=12,
        ace_trna_copy_number=4,
        nuclear_import_multiplier=1.18,
        bacterial_backbone_present=False,
        citation_doi="10.1093/nar/gkag082",
    ),
    "minidnacftr_4x": VectorArchitectureSpec(
        vector_id="minidnacftr_4x",
        label="CpG-Depleted minidnacftr Vector (1.2 kb, 4x Pol III ACE-tRNA)",
        length_bp=1200,
        cpg_count=3,
        ace_trna_copy_number=4,
        nuclear_import_multiplier=1.45,
        bacterial_backbone_present=False,
        citation_doi="10.1093/nar/gkag082",
    ),
    "dbdna_zero_cpg": VectorArchitectureSpec(
        vector_id="dbdna_zero_cpg",
        label="Synthetic Closed-Ended Doggybone DNA (dbDNA 0.95 kb, Zero-CpG 4x ACE-tRNA)",
        length_bp=950,
        cpg_count=0,
        ace_trna_copy_number=4,
        nuclear_import_multiplier=1.60,
        bacterial_backbone_present=False,
        citation_doi="10.1093/nar/gkag082",
    ),
}


def compute_cpg_density(cpg_count: int, length_bp: int) -> float:
    """Compute CpG dinucleotide frequency per base pair: rho_CpG = N_CpG / L_vector."""
    if length_bp <= 0:
        raise ValueError(f"length_bp must be positive, got {length_bp}")
    if cpg_count < 0:
        raise ValueError(f"cpg_count must be non-negative, got {cpg_count}")
    return cpg_count / float(length_bp)


def compute_tlr9_flux(
    cpg_count: int,
    endosomal_conc_nm: float,
    k_tlr9: float = K_TLR9_TURNOVER_PER_HR,
    km_nm: float = KM_ENDO_NM,
) -> float:
    """Compute endosomal TLR9 activation flux J_TLR9.

    J_TLR9 = k_TLR9 * N_CpG * ([V_endo] / (K_m + [V_endo]))
    """
    if cpg_count < 0:
        raise ValueError(f"cpg_count must be non-negative, got {cpg_count}")
    if endosomal_conc_nm < 0.0:
        raise ValueError(f"endosomal_conc_nm must be non-negative, got {endosomal_conc_nm}")
    if endosomal_conc_nm == 0.0 or cpg_count == 0:
        return 0.0
    occupancy = endosomal_conc_nm / (km_nm + endosomal_conc_nm)
    return k_tlr9 * float(cpg_count) * occupancy


def compute_cytokine_reactogenicity(
    j_tlr9: float,
    j_safe: float = J_SAFE_FLUX_THRESHOLD,
) -> float:
    """Compute dimensionless Cytokine Reactogenicity Index R_cytokine in [0, 1].

    R_cytokine = 1 - exp(-J_TLR9 / J_safe)
    """
    if j_tlr9 < 0.0:
        raise ValueError(f"j_tlr9 must be non-negative, got {j_tlr9}")
    if j_safe <= 0.0:
        raise ValueError(f"j_safe must be positive, got {j_safe}")
    return 1.0 - math.exp(-j_tlr9 / j_safe)


class MiniDNAEvasionModel:
    """Biophysical & kinetic simulator for nonviral minimal DNA ACE-tRNA vectors (T-0233)."""

    def __init__(self) -> None:
        self.vectors = VECTOR_CATALOG

    def evaluate_vector(
        self,
        vector_id: str,
        endosomal_dose_nm: float = 50.0,
        duration_days: int = 90,
        dt_hours: float = 0.5,
    ) -> dict[str, Any]:
        """Simulate TLR9 evasion, nuclear import, and 90-day episomal persistence for a DNA vector."""
        if vector_id not in self.vectors:
            raise ValueError(
                f"Unknown vector_id '{vector_id}'. Valid vectors: {sorted(self.vectors)}"
            )
        if endosomal_dose_nm < 0.0:
            raise ValueError(f"endosomal_dose_nm must be non-negative, got {endosomal_dose_nm}")
        if duration_days <= 0:
            raise ValueError(f"duration_days must be positive, got {duration_days}")
        if dt_hours <= 0.0 or dt_hours > 2.0:
            raise ValueError(f"dt_hours must be in (0, 2.0], got {dt_hours}")

        spec = self.vectors[vector_id]
        rho_cpg = compute_cpg_density(spec.cpg_count, spec.length_bp)
        j_tlr9_peak = compute_tlr9_flux(spec.cpg_count, endosomal_dose_nm)
        r_cytokine = compute_cytokine_reactogenicity(j_tlr9_peak)

        # Effective nuclear import and CpG-dependent H3K9me3 silencing rates (converted to hr^-1)
        k_nuc = K_NUC_IMPORT_PER_HR * spec.nuclear_import_multiplier
        k_silence_per_day = K_SILENCE_BASE_PER_DAY + ALPHA_CPG_SILENCE_PER_DAY * rho_cpg
        k_loss_nuc_per_hr = (k_silence_per_day + MU_DIV_AIRWAY_PER_DAY) / 24.0

        # Forward Euler integration of 3-compartment ODE:
        # d[V_endo]/dt = -(k_esc + k_deg,endo) [V_endo]
        # d[V_cyt]/dt  = k_esc [V_endo] - (k_deg,cyt + k_nuc) [V_cyt]
        # d[V_nuc]/dt  = k_nuc [V_cyt] - (k_silence + mu_div) [V_nuc]
        v_endo = float(endosomal_dose_nm)
        v_cyt = 0.0
        v_nuc = 0.0
        peak_v_nuc = 0.0

        total_steps = int(round((duration_days * 24.0) / dt_hours))
        day_snapshots: dict[int, float] = {0: 0.0}
        target_days = {1, 7, 14, 30, 60, 90, duration_days}

        for step in range(1, total_steps + 1):
            dv_endo = -(K_ESC_ENDO_PER_HR + K_DEG_ENDO_PER_HR) * v_endo
            dv_cyt = K_ESC_ENDO_PER_HR * v_endo - (K_DEG_CYT_PER_HR + k_nuc) * v_cyt
            dv_nuc = k_nuc * v_cyt - k_loss_nuc_per_hr * v_nuc

            v_endo = max(0.0, v_endo + dv_endo * dt_hours)
            v_cyt = max(0.0, v_cyt + dv_cyt * dt_hours)
            v_nuc = max(0.0, v_nuc + dv_nuc * dt_hours)
            if v_nuc > peak_v_nuc:
                peak_v_nuc = v_nuc

            elapsed_hours = step * dt_hours
            if abs(elapsed_hours % 24.0) < 1e-6:
                day_num = int(round(elapsed_hours / 24.0))
                if day_num in target_days:
                    day_snapshots[day_num] = round(v_nuc, 4)

        v_nuc_day30 = day_snapshots.get(30, v_nuc)
        v_nuc_day90 = day_snapshots.get(90, v_nuc)
        retention_day90_vs_peak_pct = (
            (v_nuc_day90 / peak_v_nuc) * 100.0 if peak_v_nuc > 0.0 else 0.0
        )

        # Steady-state Pol III ACE-tRNA expression and PTC readthrough rescue (% WT CFTR Cl- current)
        # Multicopy Pol III cassette amplifies active nuclear episomal concentration;
        # high cytokine reactogenicity (R_cytokine > 0.5) suppresses epithelial translation via PKR/IFN-beta.
        ifn_translation_penalty = 1.0 - 0.65 * (r_cytokine ** 1.5)
        effective_trna_drive = v_nuc_day30 * spec.ace_trna_copy_number * ifn_translation_penalty
        restored_cftr_current_day30_pct = min(
            38.0,
            42.0 * (effective_trna_drive / (12.0 + effective_trna_drive)),
        )

        if r_cytokine < 0.20:
            safety_tier = "CLINICALLY_TOLERABLE_LOW_REACTOGENICITY"
        elif r_cytokine < 0.55:
            safety_tier = "MODERATE_TLR9_INFLAMMATION_RISK"
        else:
            safety_tier = "SEVERE_TLR9_CYTOKINE_REACTIVITY"

        return {
            "vector_id": spec.vector_id,
            "label": spec.label,
            "length_bp": spec.length_bp,
            "cpg_count": spec.cpg_count,
            "rho_cpg_per_bp": round(rho_cpg, 6),
            "ace_trna_copy_number": spec.ace_trna_copy_number,
            "bacterial_backbone_present": spec.bacterial_backbone_present,
            "endosomal_dose_nm": round(endosomal_dose_nm, 2),
            "j_tlr9_flux": round(j_tlr9_peak, 4),
            "r_cytokine_index": round(r_cytokine, 4),
            "safety_tier": safety_tier,
            "k_silence_per_day": round(k_silence_per_day, 5),
            "episomal_half_life_days": round(
                math.log(2.0) / (k_silence_per_day + MU_DIV_AIRWAY_PER_DAY), 1
            ),
            "peak_nuclear_vector_nm": round(peak_v_nuc, 4),
            "nuclear_vector_day30_nm": round(v_nuc_day30, 4),
            "nuclear_vector_day90_nm": round(v_nuc_day90, 4),
            "retention_day90_vs_peak_pct": round(retention_day90_vs_peak_pct, 2),
            "restored_cftr_current_day30_pct_wt": round(restored_cftr_current_day30_pct, 2),
            "meets_curative_10pct_wt_day30": bool(
                restored_cftr_current_day30_pct >= 10.0 and r_cytokine < 0.25
            ),
            "day_snapshots_nm": day_snapshots,
            "citation_doi": spec.citation_doi,
            "clinical_surveillance_nct": "NCT05668741",
            "nlm_attribution": NLM_ATTRIBUTION,
            "research_boundary": RESEARCH_BOUNDARY,
        }

    def compare_all_vectors(self, endosomal_dose_nm: float = 50.0) -> dict[str, Any]:
        """Compare all 5 DNA vector architectures at a given endosomal dose."""
        reports = {
            vid: self.evaluate_vector(vid, endosomal_dose_nm=endosomal_dose_nm)
            for vid in self.vectors
        }
        plasmid_r = reports["bacterial_plasmid"]["r_cytokine_index"]
        minidna_r = reports["minidnacftr_4x"]["r_cytokine_index"]
        fold_reactogenicity_reduction = (
            plasmid_r / minidna_r if minidna_r > 0.0 else float("inf")
        )
        return {
            "endosomal_dose_nm": endosomal_dose_nm,
            "vectors": reports,
            "minidnacftr_vs_plasmid_reactogenicity_fold_reduction": round(
                fold_reactogenicity_reduction, 2
            ),
            "minidnacftr_vs_plasmid_day90_nuclear_ratio": round(
                reports["minidnacftr_4x"]["nuclear_vector_day90_nm"]
                / max(1e-6, reports["bacterial_plasmid"]["nuclear_vector_day90_nm"]),
                2,
            ),
            "nlm_attribution": NLM_ATTRIBUTION,
            "research_boundary": RESEARCH_BOUNDARY,
        }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Simulate nonviral minimal DNA vector CpG depletion and TLR9 evasion (T-0233)."
    )
    parser.add_argument(
        "--vector",
        default="minidnacftr_4x",
        choices=sorted(VECTOR_CATALOG.keys()),
        help="Vector architecture identifier",
    )
    parser.add_argument(
        "--dose-nm",
        type=float,
        default=50.0,
        help="Endosomal DNA vector concentration (nM)",
    )
    parser.add_argument(
        "--compare-all",
        action="store_true",
        help="Compare all 5 vector architectures across TLR9 reactogenicity and 90-day persistence",
    )
    parser.add_argument("--json", action="store_true", help="Emit JSON output")
    args = parser.parse_args(argv)

    model = MiniDNAEvasionModel()
    result = (
        model.compare_all_vectors(endosomal_dose_nm=args.dose_nm)
        if args.compare_all
        else model.evaluate_vector(args.vector, endosomal_dose_nm=args.dose_nm)
    )
    if args.json:
        print(json.dumps(result, indent=2))
    else:
        print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
