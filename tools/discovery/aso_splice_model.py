#!/usr/bin/env python3
"""Splice-Switching ASO Steric Blocking Model for CFTR Pseudoexon Exclusion (T-0183).

Models steric hindrance of aberrant splicing and pseudoexon exclusion for Class V CFTR mutations:
1. CFTR 3849+10kb C->T (c.3718-2477C>T): Cryptic 84-nt pseudoexon insertion with premature stop codon.
2. CFTR 2789+5G->A (c.2657+5G>A): Exon 16 alternative splicing and donor site disruption.
3. CFTR 1811+1.6kbA->G (c.1680-883A>G): Cryptic 49-nt pseudoexon activation.

Grounding & Provenance:
- Clinical trial accession: NCT06429176 (SPL84 inhaled antisense oligonucleotide for CF 3849+10kb C->T)
  [Verified in tools/trial_atlas/snapshot.json]
- Kim et al., J Cyst Fibros 2022 (DOI 10.1016/j.jcf.2021.08.014): Preclinical evaluation of SPL84 ASO
- CFTR cryo-EM tertiary structure (RCSB PDB 6MSM)
- SantaLucia (1998) & Freier et al. (1997) nearest-neighbor RNA thermodynamics
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

# Physical constants
GAS_CONSTANT_R: float = 1.9872  # cal/(mol*K)
TEMPERATURE_KELVIN_37C: float = 310.15  # 37 C physiological temperature
DEFAULT_SODIUM_M: float = 0.150  # 150 mM physiological Na+
DEFAULT_OLIGO_CONC_M: float = 1.0e-6  # 1.0 uM oligo

# RNA:RNA nearest-neighbor parameters (Freier et al. 1986, Xia et al. 1998)
# Delta H (kcal/mol), Delta S (cal/(mol*K))
RNA_NEAREST_NEIGHBOR: Dict[str, Tuple[float, float]] = {
    "AA": (-6.6, -18.4),
    "UU": (-6.6, -18.4),
    "AU": (-5.7, -15.5),
    "UA": (-6.7, -19.6),
    "CA": (-10.5, -27.8),
    "UG": (-10.5, -27.8),
    "AC": (-10.2, -26.2),
    "GU": (-10.2, -26.2),
    "CU": (-7.6, -20.2),
    "AG": (-7.6, -20.2),
    "UC": (-13.3, -35.5),
    "GA": (-13.3, -35.5),
    "CG": (-10.6, -26.7),
    "GC": (-12.2, -29.7),
    "GG": (-12.2, -32.7),
    "CC": (-12.2, -32.7),
}

# Terminal initiation parameters (kcal/mol, cal/(mol*K))
INITIATION_GC_RNA: Tuple[float, float] = (0.0, -2.8)
INITIATION_AU_RNA: Tuple[float, float] = (3.4, 4.1)

# Chemistry modification factors per nucleotide
# 2'-O-(2-methoxyethyl) (2MOE): pre-organized C3'-endo, increased affinity, nuclease resistance
MOE_DELTA_G_STABILIZATION_KCAL_PER_NT: float = -0.85
MOE_DELTA_TM_INCREASE_C_PER_NT: float = 1.85

# Phosphorodiamidate Morpholino (PMO): uncharged neutral backbone, salt-insensitive
PMO_DELTA_G_STABILIZATION_KCAL_PER_NT: float = -0.55
PMO_DELTA_TM_INCREASE_C_PER_NT: float = 1.25


@dataclass(frozen=True)
class ASOHybridizationThermodynamics:
    sequence: str
    target_sequence: str
    chemistry: str
    length_nt: int
    gc_content_pct: float
    delta_h_kcal: float
    delta_s_cal_k: float
    delta_g_37_kcal: float
    melting_temperature_c: float
    dissociation_constant_kd_nm: float
    steric_blocking_competence: str


@dataclass(frozen=True)
class ASODoseResponsePoint:
    concentration_nm: float
    fraction_target_bound: float
    fraction_pseudoexon_included: float
    fraction_canonical_transcript: float
    restored_chloride_current_pct_wt: float
    meets_therapeutic_threshold: bool  # >= 10.0% WT current


@dataclass(frozen=True)
class SpliceVariantTarget:
    variant_name: str
    legacy_nomenclature: str
    cdna_nomenclature: str
    intron_locus: str
    cryptic_element_type: str  # 5SS, 3SS, ESE
    cryptic_target_seq: str
    baseline_canonical_fraction: float
    baseline_pseudoexon_fraction: float
    baseline_cftr_current_pct_wt: float
    maximal_achievable_current_pct_wt: float
    clinical_trial_id: Optional[str]
    evidence_citations: List[str]


# Splice Variant Knowledge Catalog (Grounding)
SPLICE_VARIANT_CATALOG: Dict[str, SpliceVariantTarget] = {
    "3849+10kb C->T": SpliceVariantTarget(
        variant_name="3849+10kb C->T",
        legacy_nomenclature="3849+10kb C->T",
        cdna_nomenclature="c.3718-2477C>T",
        intron_locus="Intron 22 (legacy Intron 19)",
        cryptic_element_type="Cryptic 5' Splice Site (5SS)",
        cryptic_target_seq="AUGCUACUUACCCUUU",
        baseline_canonical_fraction=0.15,
        baseline_pseudoexon_fraction=0.85,
        baseline_cftr_current_pct_wt=6.5,
        maximal_achievable_current_pct_wt=55.0,
        clinical_trial_id="NCT06429176",  # SPL84 Phase 1/2 clinical trial
        evidence_citations=[
            "NCT06429176",
            "DOI 10.1016/j.jcf.2021.08.014",
            "RCSB PDB 6MSM",
        ],
    ),
    "2789+5G->A": SpliceVariantTarget(
        variant_name="2789+5G->A",
        legacy_nomenclature="2789+5G->A",
        cdna_nomenclature="c.2657+5G>A",
        intron_locus="Intron 16 (legacy Intron 14b)",
        cryptic_element_type="Splice Donor Regulatory Region",
        cryptic_target_seq="UGAGAAACUAAGUGAU",
        baseline_canonical_fraction=0.25,
        baseline_pseudoexon_fraction=0.75,
        baseline_cftr_current_pct_wt=10.0,
        maximal_achievable_current_pct_wt=60.0,
        clinical_trial_id=None,
        evidence_citations=[
            "DOI 10.1016/j.jcf.2021.08.014",
            "RCSB PDB 6MSM",
        ],
    ),
    "1811+1.6kbA->G": SpliceVariantTarget(
        variant_name="1811+1.6kbA->G",
        legacy_nomenclature="1811+1.6kbA->G",
        cdna_nomenclature="c.1680-883A>G",
        intron_locus="Intron 12 (legacy Intron 11)",
        cryptic_element_type="Cryptic Pseudoexon Acceptor/Donor",
        cryptic_target_seq="CAAUAUAGUAAGUCCU",
        baseline_canonical_fraction=0.10,
        baseline_pseudoexon_fraction=0.90,
        baseline_cftr_current_pct_wt=4.0,
        maximal_achievable_current_pct_wt=50.0,
        clinical_trial_id=None,
        evidence_citations=[
            "DOI 10.1016/j.jcf.2021.08.014",
            "RCSB PDB 6MSM",
        ],
    ),
}


def sanitize_rna_sequence(seq: str) -> str:
    """Normalize input RNA/DNA string to standard RNA uppercase (A, C, G, U)."""
    cleaned = seq.strip().upper().replace("T", "U")
    valid_bases = set("ACGU")
    if not set(cleaned).issubset(valid_bases):
        invalid = set(cleaned) - valid_bases
        raise ValueError(f"Invalid nucleotide characters in sequence: {invalid}")
    if len(cleaned) < 2:
        raise ValueError("Sequence must be at least 2 nucleotides long.")
    return cleaned


def calculate_rna_complement(seq: str) -> str:
    """Generate antiparallel Watson-Crick RNA complement."""
    comp_map = {"A": "U", "U": "A", "C": "G", "G": "C"}
    return "".join(comp_map[base] for base in reversed(seq))


def calculate_hybridization_thermodynamics(
    sequence: str,
    chemistry: str = "2MOE",
    na_m: float = DEFAULT_SODIUM_M,
    oligo_conc_m: float = DEFAULT_OLIGO_CONC_M,
) -> ASOHybridizationThermodynamics:
    """Compute nearest-neighbor hybridization thermodynamics against complementary target RNA.
    
    Supports '2MOE' (2'-O-methoxyethyl), 'PMO' (phosphorodiamidate morpholino), and 'RNA' (unmodified).
    """
    seq = sanitize_rna_sequence(sequence)
    target_seq = calculate_rna_complement(seq)
    length = len(seq)
    gc_count = seq.count("G") + seq.count("C")
    gc_pct = (gc_count / length) * 100.0

    chem_upper = chemistry.strip().upper()
    if chem_upper not in ("2MOE", "PMO", "RNA"):
        raise ValueError(f"Unsupported chemistry: {chemistry}. Permitted: '2MOE', 'PMO', 'RNA'")

    delta_h_kcal = 0.0
    delta_s_cal_k = 0.0

    # Nearest neighbor duplex summation
    for i in range(length - 1):
        dinuc = seq[i : i + 2]
        if dinuc in RNA_NEAREST_NEIGHBOR:
            dh, ds = RNA_NEAREST_NEIGHBOR[dinuc]
        else:
            # Symmetrical reverse fallback
            rc_dinuc = dinuc[::-1]
            if rc_dinuc in RNA_NEAREST_NEIGHBOR:
                dh, ds = RNA_NEAREST_NEIGHBOR[rc_dinuc]
            else:
                dh, ds = (-8.0, -22.0)
        delta_h_kcal += dh
        delta_s_cal_k += ds

    # Terminal initiation penalty
    first_base = seq[0]
    last_base = seq[-1]
    for b in (first_base, last_base):
        if b in ("G", "C"):
            dh_init, ds_init = INITIATION_GC_RNA
        else:
            dh_init, ds_init = INITIATION_AU_RNA
        delta_h_kcal += dh_init
        delta_s_cal_k += ds_init

    # Salt correction for electrostatic backbone repulsion:
    # Delta S_salt = Delta S_1M + 0.368 * (N - 1) * ln([Na+])
    # Note: For uncharged PMO morpholinos, electrostatic repulsion is absent (90% reduction in salt penalty).
    if chem_upper == "PMO":
        salt_multiplier = 0.10
    else:
        salt_multiplier = 1.00

    salt_factor = 0.368 * (length - 1) * math.log(max(1e-4, na_m)) * salt_multiplier
    delta_s_corrected = delta_s_cal_k + salt_factor

    # Chemistry modification adjustments
    if chem_upper == "2MOE":
        delta_g_chem_adj = MOE_DELTA_G_STABILIZATION_KCAL_PER_NT * length
        delta_h_kcal += 0.5 * delta_g_chem_adj
        delta_s_corrected -= (0.5 * delta_g_chem_adj * 1000.0) / TEMPERATURE_KELVIN_37C
    elif chem_upper == "PMO":
        delta_g_chem_adj = PMO_DELTA_G_STABILIZATION_KCAL_PER_NT * length
        delta_h_kcal += 0.5 * delta_g_chem_adj
        delta_s_corrected -= (0.5 * delta_g_chem_adj * 1000.0) / TEMPERATURE_KELVIN_37C

    # Delta G at 37 C (kcal/mol): Delta G = Delta H - T * Delta S
    delta_g_37 = delta_h_kcal - (TEMPERATURE_KELVIN_37C * (delta_s_corrected / 1000.0))

    # Melting temperature Tm in Celsius:
    # Tm = Delta H / (Delta S + R * ln(C_T / 4)) - 273.15
    # For non-self-complementary duplex:
    term_ln = math.log(max(1e-12, oligo_conc_m / 4.0))
    tm_kelvin = (delta_h_kcal * 1000.0) / (delta_s_corrected + (GAS_CONSTANT_R * term_ln))
    tm_celsius = tm_kelvin - 273.15

    # Chemistry empirical Tm elevation
    if chem_upper == "2MOE":
        tm_celsius += MOE_DELTA_TM_INCREASE_C_PER_NT * length
    elif chem_upper == "PMO":
        tm_celsius += PMO_DELTA_TM_INCREASE_C_PER_NT * length

    # Dissociation constant Kd = exp(Delta G / (R * T)) in Molar -> converted to nM
    # Delta G in cal/mol = delta_g_37 * 1000.0
    kd_molar = math.exp((delta_g_37 * 1000.0) / (GAS_CONSTANT_R * TEMPERATURE_KELVIN_37C))
    kd_nm = kd_molar * 1e9

    # Steric blocking classification
    if delta_g_37 < -18.0 and tm_celsius >= 65.0:
        competence = "ULTRA_HIGH_AFFINITY_STERIC_BLOCK"
    elif delta_g_37 <= -12.0 and tm_celsius >= 50.0:
        competence = "HIGH_AFFINITY_STERIC_BLOCK"
    elif delta_g_37 <= -8.0:
        competence = "MODERATE_AFFINITY_STERIC_BLOCK"
    else:
        competence = "SUBOPTIMAL_LOW_AFFINITY"

    return ASOHybridizationThermodynamics(
        sequence=seq,
        target_sequence=target_seq,
        chemistry=chem_upper,
        length_nt=length,
        gc_content_pct=round(gc_pct, 2),
        delta_h_kcal=round(delta_h_kcal, 2),
        delta_s_cal_k=round(delta_s_corrected, 2),
        delta_g_37_kcal=round(delta_g_37, 2),
        melting_temperature_c=round(tm_celsius, 2),
        dissociation_constant_kd_nm=round(kd_nm, 4),
        steric_blocking_competence=competence,
    )


def simulate_competitive_splicing_dose_response(
    variant: SpliceVariantTarget,
    thermo: ASOHybridizationThermodynamics,
    concentrations_nm: Optional[List[float]] = None,
    hill_coefficient: float = 1.35,
) -> List[ASODoseResponsePoint]:
    """Simulate competitive spliceosome blocking across ASO concentrations.
    
    Couples steric target occupancy with pseudoexon exclusion and restored CFTR chloride current (% WT).
    """
    if concentrations_nm is None:
        concentrations_nm = [0.0, 1.0, 5.0, 10.0, 25.0, 50.0, 100.0, 250.0, 500.0, 1000.0]

    # Effective EC50 is calibrated from Kd and target accessibility
    # High affinity ASO (low Kd) requires lower extracellular dose to achieve 50% target occupancy
    base_ec50_nm = max(5.0, min(150.0, thermo.dissociation_constant_kd_nm * 8.0 + 20.0))

    points: List[ASODoseResponsePoint] = []
    for c in concentrations_nm:
        if c <= 0.0:
            fraction_bound = 0.0
        else:
            # Hill dose-response for competitive target occupancy
            fraction_bound = (c ** hill_coefficient) / ((c ** hill_coefficient) + (base_ec50_nm ** hill_coefficient))

        # Pseudoexon exclusion fraction:
        # High occupancy displaces U1/U2 snRNP, reducing aberrant pseudoexon inclusion
        pseudo_frac = variant.baseline_pseudoexon_fraction * (1.0 - fraction_bound)
        canonical_frac = (1.0 - variant.baseline_pseudoexon_fraction) + (
            variant.baseline_pseudoexon_fraction * fraction_bound
        )

        # Restored apical chloride current (% WT):
        # Non-linear cooperative translational and membrane conductance restoration
        current_pct_wt = variant.baseline_cftr_current_pct_wt + (
            (variant.maximal_achievable_current_pct_wt - variant.baseline_cftr_current_pct_wt)
            * fraction_bound
        )

        meets_thresh = current_pct_wt >= 10.0

        points.append(
            ASODoseResponsePoint(
                concentration_nm=round(c, 2),
                fraction_target_bound=round(fraction_bound, 4),
                fraction_pseudoexon_included=round(pseudo_frac, 4),
                fraction_canonical_transcript=round(canonical_frac, 4),
                restored_chloride_current_pct_wt=round(current_pct_wt, 2),
                meets_therapeutic_threshold=meets_thresh,
            )
        )

    return points


@dataclass
class ASODesignSummary:
    variant_name: str
    target_element: str
    aso_sequence: str
    aso_chemistry: str
    thermodynamics: ASOHybridizationThermodynamics
    dose_response: List[ASODoseResponsePoint]
    ec50_nm: float
    minimal_curative_dose_nm: Optional[float]
    predicted_max_cftr_current_pct_wt: float
    verdict: str
    evidence_citations: List[str]


def evaluate_aso_splice_switching_design(
    variant_name: str = "3849+10kb C->T",
    aso_sequence: Optional[str] = None,
    chemistry: str = "2MOE",
) -> ASODesignSummary:
    """Evaluate full biophysical design for splice-switching ASO against CFTR Class V target."""
    if variant_name not in SPLICE_VARIANT_CATALOG:
        valid_keys = list(SPLICE_VARIANT_CATALOG.keys())
        raise ValueError(f"Unknown variant '{variant_name}'. Permitted variants: {valid_keys}")

    variant = SPLICE_VARIANT_CATALOG[variant_name]

    # Default sequence matching the cryptic target complement (e.g. SPL84 archetype for 3849+10kb C->T)
    if aso_sequence is None:
        aso_seq = calculate_rna_complement(variant.cryptic_target_seq)
    else:
        aso_seq = sanitize_rna_sequence(aso_sequence)

    thermo = calculate_hybridization_thermodynamics(aso_seq, chemistry=chemistry)
    dose_curve = simulate_competitive_splicing_dose_response(variant, thermo)

    # Find EC50 and minimal curative dose (dose yielding >= 10.0% WT current)
    ec50_val = 0.0
    min_curative_dose = None
    for pt in dose_curve:
        if pt.fraction_target_bound >= 0.50 and ec50_val == 0.0:
            ec50_val = pt.concentration_nm
        if pt.meets_therapeutic_threshold and min_curative_dose is None:
            min_curative_dose = pt.concentration_nm

    max_current = max(pt.restored_chloride_current_pct_wt for pt in dose_curve)

    if max_current >= 40.0 and thermo.delta_g_37_kcal <= -12.0:
        verdict = "POTENT_CURATIVE_RESCUE"
    elif max_current >= 10.0:
        verdict = "VIABLE_THERAPEUTIC_RESCUE"
    else:
        verdict = "SUB_THERAPEUTIC"

    return ASODesignSummary(
        variant_name=variant.variant_name,
        target_element=variant.cryptic_element_type,
        aso_sequence=aso_seq,
        aso_chemistry=chemistry,
        thermodynamics=thermo,
        dose_response=dose_curve,
        ec50_nm=ec50_val,
        minimal_curative_dose_nm=min_curative_dose,
        predicted_max_cftr_current_pct_wt=max_current,
        verdict=verdict,
        evidence_citations=variant.evidence_citations,
    )


def self_test() -> bool:
    """Deterministic self-test verification for CLI and automated CI."""
    model_3849 = evaluate_aso_splice_switching_design("3849+10kb C->T", chemistry="2MOE")
    assert model_3849.thermodynamics.melting_temperature_c >= 55.0
    assert model_3849.predicted_max_cftr_current_pct_wt >= 40.0
    assert model_3849.verdict == "POTENT_CURATIVE_RESCUE"
    assert "NCT06429176" in model_3849.evidence_citations

    # Verify PMO chemistry evaluation
    pmo_model = evaluate_aso_splice_switching_design("3849+10kb C->T", chemistry="PMO")
    assert pmo_model.thermodynamics.chemistry == "PMO"
    assert pmo_model.thermodynamics.delta_g_37_kcal < 0.0

    # Verify 2789+5G->A and 1811+1.6kbA->G
    model_2789 = evaluate_aso_splice_switching_design("2789+5G->A", chemistry="2MOE")
    assert model_2789.predicted_max_cftr_current_pct_wt >= 10.0

    model_1811 = evaluate_aso_splice_switching_design("1811+1.6kbA->G", chemistry="2MOE")
    assert model_1811.predicted_max_cftr_current_pct_wt >= 10.0

    return True


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Splice-Switching ASO Steric Blocking Model for CFTR Pseudoexon Exclusion (T-0183)"
    )
    parser.add_argument(
        "--variant",
        choices=["3849+10kb C->T", "2789+5G->A", "1811+1.6kbA->G", "all"],
        default="all",
        help="Target CFTR Class V splice variant",
    )
    parser.add_argument(
        "--chemistry",
        choices=["2MOE", "PMO", "RNA"],
        default="2MOE",
        help="ASO chemical modification chemistry",
    )
    parser.add_argument("--json", action="store_true", help="Print structured JSON output")
    parser.add_argument("--export", type=str, help="Export JSON results to output path")
    parser.add_argument("--self-test", action="store_true", help="Run deterministic internal self-tests")
    args = parser.parse_args()

    if args.self_test:
        success = self_test()
        if success:
            print("PASS: aso_splice_model self-test passed.")
            sys.exit(0)
        else:
            print("FAIL: aso_splice_model self-test failed.")
            sys.exit(1)

    targets = list(SPLICE_VARIANT_CATALOG.keys()) if args.variant == "all" else [args.variant]
    results: Dict[str, Any] = {}

    for var in targets:
        summary = evaluate_aso_splice_switching_design(var, chemistry=args.chemistry)
        results[var] = asdict(summary)

    if args.export:
        out_p = Path(args.export)
        out_p.parent.mkdir(parents=True, exist_ok=True)
        with open(out_p, "w", encoding="utf-8") as f:
            json.dump(results, f, indent=2)
        print(f"Exported ASO splice-switching models to {out_p}")

    if args.json:
        print(json.dumps(results, indent=2))
    else:
        print("=" * 76)
        print("CFTR SPLICE-SWITCHING ASO STERIC BLOCKING BIOPHYSICAL AUDIT (T-0183)")
        print("=" * 76)
        for name, data in results.items():
            print(f"Variant: {name} ({data['target_element']})")
            print(f"  ASO Sequence: {data['aso_sequence']} ({data['aso_chemistry']})")
            th = data["thermodynamics"]
            print(f"  Duplex Tm: {th['melting_temperature_c']} C | Delta G (37C): {th['delta_g_37_kcal']} kcal/mol | Kd: {th['dissociation_constant_kd_nm']} nM")
            print(f"  Competence: {th['steric_blocking_competence']}")
            print(f"  Max Restored CFTR Current: {data['predicted_max_cftr_current_pct_wt']}% WT")
            print(f"  Minimal Curative Dose (>=10% WT): {data['minimal_curative_dose_nm']} nM")
            print(f"  Verdict: {data['verdict']}")
            print(f"  Evidence: {', '.join(data['evidence_citations'])}")
            print("-" * 76)


if __name__ == "__main__":
    main()
