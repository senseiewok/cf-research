#!/usr/bin/env python3
"""Boltzmann Spliceosome-ASO Steric-Blocking Kinetics & NMD Rescue Model (T-0254).

Models competitive spliceosome recognition (U1 snRNP 5'SS and U2AF 3'SS assembly,
RCSB PDB 6FF4 and 7ABI) versus steric-blocking antisense oligonucleotides (2'-O-MOE-PS
and PMO morpholino) across Class V and deep-intronic CFTR splicing mutations:
  1. 3849+10kbC->T (c.3718-2477C>T): Intron 22 cryptic 5' donor -> 84-nt pseudoexon with in-frame UAA stop.
  2. 2789+5G->A (c.2657+5G>A): Intron 16 weakened canonical 5' donor -> Exon 16 skipping (frameshift).
  3. 1811+1.6kbA->G (c.1679+1.6kbA>G): Intron 11 strong cryptic 5' donor -> 49-nt pseudoexon (frameshift).

Grounding & Provenance (Answer Tier T2):
  - ClinicalTrials.gov: NCT06429176 (SPL84 inhaled ASO Phase 1/2 trial in 3849+10kbC->T CF)
  - ClinicalTrials.gov: NCT04684641 (CFTR modulator clinical benchmark in rare splice genotypes)
  - RCSB PDB: 6FF4 (Human pre-B spliceosome U1/U2 snRNP architecture), 7ABI (Human U1 snRNP 5'SS duplex),
    6MSM (Phosphorylated, ATP-bound human CFTR channel)
  - Primary Literature:
    * Oren et al., Am J Respir Crit Care Med 2021 (DOI: 10.1164/rccm.202105-1155OC)
    * Kim et al., Nat Commun 2022 (DOI: 10.1038/s41467-022-30495-1)
    * Michaels et al., Mol Ther 2020 (DOI: 10.1016/j.ymthe.2020.09.017)
    * Oren et al., J Cyst Fibros 2022 (DOI: 10.1016/j.jcf.2021.08.014)
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from dataclasses import asdict, dataclass
from typing import Any, Dict, List

GAS_CONSTANT_KCAL_MOL_K: float = 1.987204e-3
PHYSIOLOGICAL_TEMP_K: float = 310.15  # 37.0 C
RT_37C_KCAL_MOL: float = GAS_CONSTANT_KCAL_MOL_K * PHYSIOLOGICAL_TEMP_K  # ~0.6163 kcal/mol


@dataclass(frozen=True)
class SpliceKineticsAlleleSpec:
    allele_id: str
    legacy_name: str
    hgvs_cdna: str
    intron_location: str
    splice_defect_mechanism: str
    pseudoexon_length_nt: int
    contains_premature_stop: bool
    delta_g_u1_canonical_kcal: float
    delta_g_u1_cryptic_kcal: float
    aberrant_nmd_escape_fraction: float
    canonical_nmd_escape_fraction: float
    single_channel_conductance_scale: float
    eti_potentiation_multiplier: float
    clinical_trial_nct: str
    structural_pdbs: List[str]
    primary_dois: List[str]


ALLELE_KINETICS_SPECS: Dict[str, SpliceKineticsAlleleSpec] = {
    "3849+10kbC->T": SpliceKineticsAlleleSpec(
        allele_id="3849+10kbC->T",
        legacy_name="3849+10kbC->T",
        hgvs_cdna="c.3718-2477C>T",
        intron_location="Intron 22",
        splice_defect_mechanism="Cryptic 5' splice donor activation causing 84-nt pseudoexon inclusion with in-frame UAA stop codon",
        pseudoexon_length_nt=84,
        contains_premature_stop=True,
        delta_g_u1_canonical_kcal=-8.20,
        delta_g_u1_cryptic_kcal=-9.30,
        aberrant_nmd_escape_fraction=0.22,
        canonical_nmd_escape_fraction=1.00,
        single_channel_conductance_scale=0.68,
        eti_potentiation_multiplier=1.25,
        clinical_trial_nct="NCT06429176",
        structural_pdbs=["6FF4", "7ABI", "6MSM"],
        primary_dois=[
            "10.1164/rccm.202105-1155OC",
            "10.1038/s41467-022-30495-1",
            "10.1016/j.jcf.2021.08.014",
        ],
    ),
    "2789+5G->A": SpliceKineticsAlleleSpec(
        allele_id="2789+5G->A",
        legacy_name="2789+5G->A",
        hgvs_cdna="c.2657+5G>A",
        intron_location="Intron 16",
        splice_defect_mechanism="Weakened canonical Exon 16 5' donor (+5 G>A) competing with intronic splicing silencer (ISS) causing Exon 16 skipping",
        pseudoexon_length_nt=0,
        contains_premature_stop=True,
        delta_g_u1_canonical_kcal=-7.50,
        delta_g_u1_cryptic_kcal=-8.85,
        aberrant_nmd_escape_fraction=0.18,
        canonical_nmd_escape_fraction=1.00,
        single_channel_conductance_scale=0.62,
        eti_potentiation_multiplier=1.35,
        clinical_trial_nct="NCT04684641",
        structural_pdbs=["6FF4", "7ABI", "6MSM"],
        primary_dois=[
            "10.1038/s41467-022-30495-1",
            "10.1016/j.ymthe.2020.09.017",
        ],
    ),
    "1811+1.6kbA->G": SpliceKineticsAlleleSpec(
        allele_id="1811+1.6kbA->G",
        legacy_name="1811+1.6kbA->G",
        hgvs_cdna="c.1679+1.6kbA>G",
        intron_location="Intron 11",
        splice_defect_mechanism="High-affinity cryptic 5' donor creation driving 49-nt pseudoexon inclusion and frameshift termination",
        pseudoexon_length_nt=49,
        contains_premature_stop=True,
        delta_g_u1_canonical_kcal=-8.00,
        delta_g_u1_cryptic_kcal=-9.95,
        aberrant_nmd_escape_fraction=0.14,
        canonical_nmd_escape_fraction=1.00,
        single_channel_conductance_scale=0.65,
        eti_potentiation_multiplier=1.28,
        clinical_trial_nct="NCT06429176",
        structural_pdbs=["6FF4", "7ABI", "6MSM"],
        primary_dois=[
            "10.1016/j.ymthe.2020.09.017",
            "10.1164/rccm.202105-1155OC",
        ],
    ),
}

CHEMISTRY_KINETICS: Dict[str, Dict[str, float]] = {
    "2MOE": {
        "kd_nm": 38.0,
        "hill_n": 1.25,
        "steric_penalty_kcal": 3.40,
        "endosomal_escape_efficiency": 0.12,
    },
    "PMO": {
        "kd_nm": 64.0,
        "hill_n": 1.15,
        "steric_penalty_kcal": 2.85,
        "endosomal_escape_efficiency": 0.09,
    },
}


@dataclass(frozen=True)
class BoltzmannSpliceResult:
    allele_id: str
    chemistry: str
    aso_concentration_nm: float
    eti_co_administered: bool
    aso_target_occupancy: float
    effective_cryptic_delta_g_kcal: float
    canonical_splice_fraction: float
    aberrant_splice_fraction: float
    steady_state_mrna_pct_wt: float
    restored_isc_pct_wt: float
    curative_threshold_met: bool


def resolve_allele_spec(allele_id: str) -> SpliceKineticsAlleleSpec:
    normalized = allele_id.strip()
    if normalized in ALLELE_KINETICS_SPECS:
        return ALLELE_KINETICS_SPECS[normalized]
    aliases = {
        "3849+10kb C->T": "3849+10kbC->T",
        "3849": "3849+10kbC->T",
        "2789": "2789+5G->A",
        "1811": "1811+1.6kbA->G",
    }
    if normalized in aliases:
        return ALLELE_KINETICS_SPECS[aliases[normalized]]
    raise ValueError(f"Unsupported splice allele: {allele_id!r}. Valid: {sorted(ALLELE_KINETICS_SPECS)}")


def compute_boltzmann_splice_partition(
    allele_id: str,
    aso_concentration_nm: float,
    chemistry: str = "2MOE",
    eti_co_administered: bool = False,
) -> BoltzmannSpliceResult:
    """Compute Boltzmann spliceosome partition, NMD escape, and restored CFTR current."""
    if aso_concentration_nm < 0.0:
        raise ValueError("aso_concentration_nm must be non-negative")

    spec = resolve_allele_spec(allele_id)
    chem_key = chemistry.strip().upper()
    if chem_key not in CHEMISTRY_KINETICS:
        raise ValueError(f"Unsupported ASO chemistry: {chemistry!r}. Valid: {sorted(CHEMISTRY_KINETICS)}")

    chem = CHEMISTRY_KINETICS[chem_key]
    kd = chem["kd_nm"]
    hill_n = chem["hill_n"]

    if aso_concentration_nm == 0.0:
        theta_aso = 0.0
    else:
        c_n = math.pow(aso_concentration_nm, hill_n)
        theta_aso = c_n / (c_n + math.pow(kd, hill_n))

    # Steric ASO hybridization penalizes U1 snRNP / ISS binding at the cryptic/regulatory locus
    dg_cryptic_eff = spec.delta_g_u1_cryptic_kcal + theta_aso * chem["steric_penalty_kcal"]
    dg_canon = spec.delta_g_u1_canonical_kcal

    # Boltzmann weights: W_i = exp(-DeltaG_i / RT)
    # Use relative Delta-DeltaG to avoid floating-point overflow
    ddg = dg_cryptic_eff - dg_canon
    ratio_cryptic_to_canon = math.exp(-ddg / RT_37C_KCAL_MOL)

    f_canon = 1.0 / (1.0 + ratio_cryptic_to_canon)
    f_aberrant = 1.0 - f_canon

    # Steady-state mRNA pool accounting for NMD degradation of PTC-bearing aberrant isoform
    mrna_fraction = (
        f_canon * spec.canonical_nmd_escape_fraction
        + f_aberrant * spec.aberrant_nmd_escape_fraction
    )
    mrna_pct_wt = mrna_fraction * 100.0

    # Functional chloride transport comes from canonically spliced full-length CFTR transcripts
    isc_pct = f_canon * 100.0 * spec.single_channel_conductance_scale
    if eti_co_administered:
        isc_pct = min(95.0, isc_pct * spec.eti_potentiation_multiplier)

    return BoltzmannSpliceResult(
        allele_id=spec.allele_id,
        chemistry=chem_key,
        aso_concentration_nm=round(aso_concentration_nm, 2),
        eti_co_administered=eti_co_administered,
        aso_target_occupancy=round(theta_aso, 4),
        effective_cryptic_delta_g_kcal=round(dg_cryptic_eff, 3),
        canonical_splice_fraction=round(f_canon, 4),
        aberrant_splice_fraction=round(f_aberrant, 4),
        steady_state_mrna_pct_wt=round(mrna_pct_wt, 2),
        restored_isc_pct_wt=round(isc_pct, 2),
        curative_threshold_met=bool(isc_pct >= 10.0),
    )


def simulate_allele_dose_response(
    allele_id: str,
    chemistry: str = "2MOE",
    eti_co_administered: bool = False,
    concentrations_nm: List[float] | None = None,
) -> Dict[str, Any]:
    """Generate a full dose-response curve and compute EC50 for canonical exon restoration."""
    spec = resolve_allele_spec(allele_id)
    doses = concentrations_nm or [0.0, 5.0, 10.0, 25.0, 50.0, 100.0, 200.0, 500.0]
    points = [
        asdict(compute_boltzmann_splice_partition(spec.allele_id, d, chemistry, eti_co_administered))
        for d in doses
    ]
    baseline = points[0]["canonical_splice_fraction"]
    maximal = points[-1]["canonical_splice_fraction"]
    half_target = baseline + 0.5 * (maximal - baseline)

    # Estimate apparent EC50 by fine grid interpolation
    ec50_nm = CHEMISTRY_KINETICS[chemistry.strip().upper()]["kd_nm"]
    for test_c in range(1, 500):
        res = compute_boltzmann_splice_partition(spec.allele_id, float(test_c), chemistry, eti_co_administered)
        if res.canonical_splice_fraction >= half_target:
            ec50_nm = float(test_c)
            break

    return {
        "allele": asdict(spec),
        "chemistry": chemistry.strip().upper(),
        "eti_co_administered": eti_co_administered,
        "apparent_ec50_nm": ec50_nm,
        "dose_response": points,
    }


def main(argv: List[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Boltzmann Spliceosome-ASO Steric-Blocking Kinetics Model")
    parser.add_argument("--allele", default="3849+10kbC->T", help="CFTR splice allele")
    parser.add_argument("--chemistry", default="2MOE", choices=["2MOE", "PMO"])
    parser.add_argument("--eti", action="store_true", help="Include ETI modulator synergy")
    parser.add_argument("--json", action="store_true", help="Print JSON summary")
    args = parser.parse_args(argv)

    summary = simulate_allele_dose_response(args.allele, args.chemistry, args.eti)
    if args.json:
        print(json.dumps(summary, indent=2))
    else:
        print(f"Allele: {summary['allele']['allele_id']} ({summary['allele']['hgvs_cdna']})")
        print(f"Chemistry: {summary['chemistry']} | Apparent EC50: {summary['apparent_ec50_nm']} nM")
        for pt in summary["dose_response"]:
            print(
                f"  [{pt['aso_concentration_nm']:6.1f} nM] "
                f"theta={pt['aso_target_occupancy']:.3f} | "
                f"f_canon={pt['canonical_splice_fraction']*100:5.1f}% | "
                f"I_sc={pt['restored_isc_pct_wt']:5.1f}% WT"
            )
    return 0


if __name__ == "__main__":
    sys.exit(main())
