#!/usr/bin/env python3
"""Biophysical Feasibility Scorer for CFTR Prime Editing (T-0184).

Provides quantitative biophysical modeling and constraint scoring for:
1. pegRNA Primer Binding Site (PBS) nearest-neighbor thermodynamics and melting temperature (Tm 30-37C).
2. Reverse Transcriptase Template (RTT) length constraints and processivity limits (10-35 nt).
3. Engineered pegRNA (epegRNA) 3' structural motif stabilization (mpknot, tevpre).
4. PE3 / PE3b nicking sgRNA spacing and indel hazard mitigation.
5. Delivery vector packaging limits (Single AAV vs Dual AAV Split-Intein vs SORT-LNP mRNA).
6. Target CFTR loci: F508del (+CTT in-frame insertion) and W1282X (c.3846A>G nonsense reversal).

Scientific Grounding:
- Chen et al., Nat Biotechnol 2023 (DOI 10.1038/s41587-023-01851-6)
- Nelson et al., Nat Biotechnol 2022 (DOI 10.1038/s41587-021-01039-7)
- Anzalone et al., Nature 2019 (DOI 10.1038/s41586-019-1711-4)
- MTNA CFTR prime editing (DOI 10.1016/j.omtn.2022.04.025)
- CFTR Cryo-EM coordinates (RCSB PDB 6MSM)
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

# Physical & Thermodynamic Constants (SantaLucia 1998 unified parameters)
# Delta H (kcal/mol), Delta S (cal/(mol*K)) for DNA:DNA nearest-neighbor duplexe
NEAREST_NEIGHBOR_NN: Dict[str, Tuple[float, float]] = {
    "AA": (-7.6, -21.3),
    "TT": (-7.6, -21.3),
    "AT": (-7.2, -20.4),
    "TA": (-7.2, -21.3),
    "CA": (-8.5, -22.7),
    "TG": (-8.5, -22.7),
    "GT": (-8.4, -22.4),
    "AC": (-8.4, -22.4),
    "CT": (-7.8, -21.0),
    "AG": (-7.8, -21.0),
    "GA": (-8.2, -22.2),
    "TC": (-8.2, -22.2),
    "CG": (-10.6, -27.2),
    "GC": (-9.8, -24.4),
    "GG": (-8.0, -19.9),
    "CC": (-8.0, -19.9),
}

INITIATION_GC: Tuple[float, float] = (0.1, -2.8)
INITIATION_AT: Tuple[float, float] = (2.3, 4.1)
GAS_CONSTANT_R: float = 1.9872  # cal/(mol*K)
DEFAULT_SODIUM_M: float = 0.050  # 50 mM Na+
DEFAULT_OLIGO_CONC_M: float = 0.5e-6  # 0.5 uM primer

# Biophysical Design Thresholds for Prime Editing
PBS_TM_OPTIMAL_MIN_C: float = 30.0
PBS_TM_OPTIMAL_MAX_C: float = 37.0
PBS_LENGTH_MIN_NT: int = 10
PBS_LENGTH_MAX_NT: int = 17

RTT_LENGTH_MIN_NT: int = 10
RTT_LENGTH_OPTIMAL_MIN_NT: int = 12
RTT_LENGTH_OPTIMAL_MAX_NT: int = 22
RTT_LENGTH_MAX_NT: int = 35

NICK_DISTANCE_OPTIMAL_MIN_BP: int = 40
NICK_DISTANCE_OPTIMAL_MAX_BP: int = 90
NICK_DISTANCE_INDEL_HAZARD_BP: int = 35

AAV_PACKAGING_LIMIT_BP: int = 4700
LNP_MRNA_PACKAGING_LIMIT_NT: int = 10000

# Canonical Lengths for Prime Editor Architectures
PEMAX_CDS_BP: int = 6321  # SpCas9(H840A) + Linker + engineered M-MLV RT + bipartite NLS
EPEGRNA_AVERAGE_NT: int = 180  # Spacer (20) + Scaffold (76) + RTT (16) + PBS (13) + Linker (8) + Motif (47)
NICKING_SGRNA_NT: int = 96  # Spacer (20) + Scaffold (76)


def calculate_melting_temperature_nn(
    sequence: str,
    na_m: float = DEFAULT_SODIUM_M,
    oligo_m: float = DEFAULT_OLIGO_CONC_M,
) -> float:
    """Calculate melting temperature (Tm in Celsius) using SantaLucia (1998) nearest-neighbor thermodynamics.
    
    Includes sodium salt correction for entropy:
    Delta S_salt = Delta S_1M + 0.368 * (N - 1) * ln([Na+])
    """
    seq = sequence.upper().strip()
    if len(seq) < 2:
        raise ValueError("Sequence must be at least 2 nucleotides long.")
    valid_bases = set("ATGC")
    if not set(seq).issubset(valid_bases):
        raise ValueError(f"Invalid sequence characters. Only A, T, G, C permitted: {seq}")

    delta_h_kcal = 0.0
    delta_s_eu = 0.0

    # Nearest-neighbor duplex summation
    for i in range(len(seq) - 1):
        dinuc = seq[i : i + 2]
        if dinuc in NEAREST_NEIGHBOR_NN:
            dh, ds = NEAREST_NEIGHBOR_NN[dinuc]
            delta_h_kcal += dh
            delta_s_eu += ds
        else:
            # Fallback for symmetrical reverse complement
            rc_dinuc = dinuc[::-1]
            if rc_dinuc in NEAREST_NEIGHBOR_NN:
                dh, ds = NEAREST_NEIGHBOR_NN[rc_dinuc]
                delta_h_kcal += dh
                delta_s_eu += ds

    # Terminal initiation penalties
    first_base = seq[0]
    last_base = seq[-1]
    init_dh = 0.0
    init_ds = 0.0

    for base in (first_base, last_base):
        if base in ("G", "C"):
            init_dh += INITIATION_GC[0]
            init_ds += INITIATION_GC[1]
        else:
            init_dh += INITIATION_AT[0]
            init_ds += INITIATION_AT[1]

    delta_h_kcal += init_dh
    delta_s_eu += init_ds

    # Salt correction on entropy
    n_bases = len(seq)
    salt_correction = 0.368 * (n_bases - 1) * math.log(na_m)
    delta_s_corrected = delta_s_eu + salt_correction

    # Convert delta H to cal/mol
    delta_h_cal = delta_h_kcal * 1000.0

    # Tm equation for non-self-complementary oligo with excess primer
    # Tm = Delta H / (Delta S + R * ln(Ct / 4)) - 273.15
    denominator = delta_s_corrected + (GAS_CONSTANT_R * math.log(oligo_m / 4.0))
    if abs(denominator) < 1e-9:
        return 0.0

    tm_kelvin = delta_h_cal / denominator
    tm_celsius = tm_kelvin - 273.15
    return round(tm_celsius, 2)


def calculate_gc_content(sequence: str) -> float:
    """Compute GC content percentage (0.0 to 100.0)."""
    seq = sequence.upper().strip()
    if not seq:
        return 0.0
    gc_count = sum(1 for b in seq if b in ("G", "C"))
    return round((gc_count / len(seq)) * 100.0, 2)


@dataclass
class PBSScore:
    sequence: str
    length_nt: int
    melting_temperature_c: float
    gc_content_pct: float
    within_optimal_bounds: bool
    score: float
    details: str


@dataclass
class RTTScore:
    sequence: str
    length_nt: int
    nick_to_edit_distance_nt: int
    within_constraints: bool
    score: float
    details: str


@dataclass
class EpegrnaMotifScore:
    motif_name: str
    is_engineered: bool
    stability_retention_score: float
    linker_sequence: str
    details: str


@dataclass
class NickingSgRNAScore:
    nicking_sequence: str
    relative_distance_bp: int
    is_pe3b: bool
    indel_hazard_flag: bool
    score: float
    details: str


@dataclass
class DeliveryPackagingScore:
    system_name: str
    single_aav_feasible: bool
    dual_aav_split_intein_feasible: bool
    lnp_mrna_feasible: bool
    pemax_size_bp: int
    limiting_vector_kb: float
    recommendation: str


@dataclass
class PrimeEditingDesignResult:
    target_id: str
    cftr_mutation: str
    edit_type: str
    cryoem_pdb_context: str
    pbs: PBSScore
    rtt: RTTScore
    epegRNA: EpegrnaMotifScore
    nicking_guide: NickingSgRNAScore
    packaging: DeliveryPackagingScore
    composite_feasibility_score: float
    verdict: str
    evidence_citations: List[str] = field(default_factory=list)


def score_pbs(sequence: str) -> PBSScore:
    """Evaluate PBS thermodynamic melting temperature and length constraints."""
    seq = sequence.upper().strip()
    length_nt = len(seq)
    tm_c = calculate_melting_temperature_nn(seq)
    gc_pct = calculate_gc_content(seq)

    within_bounds = (PBS_TM_OPTIMAL_MIN_C <= tm_c <= PBS_TM_OPTIMAL_MAX_C) and (
        PBS_LENGTH_MIN_NT <= length_nt <= PBS_LENGTH_MAX_NT
    )

    # Score calculation
    # Penalty for deviating from 30-37 C
    if within_bounds:
        # Peak score near 33.5 C (center of optimal window)
        dist_center = abs(tm_c - 33.5)
        score = max(0.85, 1.0 - (dist_center * 0.03))
        details = f"Optimal PBS duplex stability: Tm={tm_c} C in [30.0, 37.0] C range; length={length_nt} nt."
    elif tm_c < PBS_TM_OPTIMAL_MIN_C:
        deficit = PBS_TM_OPTIMAL_MIN_C - tm_c
        score = max(0.2, 0.75 - (deficit * 0.08))
        details = f"Suboptimal low PBS Tm={tm_c} C (<30.0 C); reduced reverse transcription initiation kinetics."
    else:
        excess = tm_c - PBS_TM_OPTIMAL_MAX_C
        score = max(0.15, 0.70 - (excess * 0.07))
        details = f"High PBS Tm={tm_c} C (>37.0 C); secondary structure / flap-release inhibition risk."

    return PBSScore(
        sequence=seq,
        length_nt=length_nt,
        melting_temperature_c=tm_c,
        gc_content_pct=gc_pct,
        within_optimal_bounds=within_bounds,
        score=round(score, 3),
        details=details,
    )


def score_rtt(sequence: str, nick_to_edit_distance_nt: int) -> RTTScore:
    """Evaluate Reverse Transcriptase Template (RTT) length and distance constraints."""
    seq = sequence.upper().strip()
    length_nt = len(seq)

    within_constraints = (RTT_LENGTH_MIN_NT <= length_nt <= RTT_LENGTH_MAX_NT) and (
        nick_to_edit_distance_nt <= 25
    )

    if length_nt < RTT_LENGTH_MIN_NT:
        score = 0.1
        details = f"Inviable RTT length {length_nt} nt (< {RTT_LENGTH_MIN_NT} nt minimum homology threshold)."
    elif length_nt > RTT_LENGTH_MAX_NT:
        score = 0.25
        details = f"Exceeded processivity limit: RTT length {length_nt} nt (> {RTT_LENGTH_MAX_NT} nt); elevated synthesis error rate."
    elif RTT_LENGTH_OPTIMAL_MIN_NT <= length_nt <= RTT_LENGTH_OPTIMAL_MAX_NT:
        dist_penalty = max(0.0, (nick_to_edit_distance_nt - 5) * 0.02)
        score = max(0.7, 1.0 - dist_penalty)
        details = f"Optimal RTT length {length_nt} nt in [{RTT_LENGTH_OPTIMAL_MIN_NT}, {RTT_LENGTH_OPTIMAL_MAX_NT}] nt; nick-to-edit dist {nick_to_edit_distance_nt} nt."
    else:
        # Between 23 and 35 nt
        score = 0.75 - ((length_nt - RTT_LENGTH_OPTIMAL_MAX_NT) * 0.02)
        details = f"Acceptable RTT length {length_nt} nt; mild RT dissociation penalty above {RTT_LENGTH_OPTIMAL_MAX_NT} nt."

    return RTTScore(
        sequence=seq,
        length_nt=length_nt,
        nick_to_edit_distance_nt=nick_to_edit_distance_nt,
        within_constraints=within_constraints,
        score=round(score, 3),
        details=details,
    )


def score_epegRNA_motif(motif_name: str, linker_sequence: str = "CGCGGTTX") -> EpegrnaMotifScore:
    """Evaluate 3' structural RNA motif stabilization."""
    motif = motif_name.lower().strip()
    if motif in ("mpknot", "evopre", "tevpre"):
        return EpegrnaMotifScore(
            motif_name=motif,
            is_engineered=True,
            stability_retention_score=1.0,
            linker_sequence=linker_sequence,
            details=f"Engineered epegRNA with 3' {motif} structured motif resists cellular 3'->5' exonucleolytic decay (3.2x enhancement).",
        )
    elif motif in ("none", "linear", "unmodified"):
        return EpegrnaMotifScore(
            motif_name="unmodified",
            is_engineered=False,
            stability_retention_score=0.35,
            linker_sequence="",
            details="Linear unmodified pegRNA extension is vulnerable to intracellular exonucleolytic degradation.",
        )
    else:
        return EpegrnaMotifScore(
            motif_name=motif,
            is_engineered=True,
            stability_retention_score=0.75,
            linker_sequence=linker_sequence,
            details=f"Non-canonical structural motif '{motif}' provides moderate structural stabilization.",
        )


def score_nicking_guide(
    sequence: str,
    relative_distance_bp: int,
    is_pe3b: bool = False,
) -> NickingSgRNAScore:
    """Evaluate PE3 / PE3b nicking guide position and indel hazard."""
    seq = sequence.upper().strip()
    dist = abs(relative_distance_bp)

    indel_hazard = dist < NICK_DISTANCE_INDEL_HAZARD_BP and not is_pe3b

    if is_pe3b:
        score = 1.0
        details = f"PE3b edit-specific nicking design installed at {relative_distance_bp} bp; cuts only edited allele, reducing indels to <0.5%."
    elif indel_hazard:
        score = 0.2
        details = f"High indel hazard: nicking distance {dist} bp is under {NICK_DISTANCE_INDEL_HAZARD_BP} bp threshold; simultaneous double-strand break risk."
    elif NICK_DISTANCE_OPTIMAL_MIN_BP <= dist <= NICK_DISTANCE_OPTIMAL_MAX_BP:
        score = 0.90
        details = f"Optimal PE3 nicking interval: {dist} bp downstream/upstream promotes mismatch repair pathway preference."
    elif dist > 120:
        score = 0.50
        details = f"Distal nicking interval: {dist} bp exceeds 120 bp; reduced mismatch repair synergy."
    else:
        score = 0.70
        details = f"Acceptable nicking interval: {dist} bp."

    return NickingSgRNAScore(
        nicking_sequence=seq,
        relative_distance_bp=relative_distance_bp,
        is_pe3b=is_pe3b,
        indel_hazard_flag=indel_hazard,
        score=round(score, 3),
        details=details,
    )


def evaluate_packaging_limits(system_cds_bp: int = PEMAX_CDS_BP) -> DeliveryPackagingScore:
    """Evaluate delivery vector feasibility against strict packaging capacities."""
    single_aav_ok = system_cds_bp <= AAV_PACKAGING_LIMIT_BP
    # Dual AAV split-intein divides PEmax into ~3.3 kb and ~3.8 kb cassettes, both < 4.7 kb
    dual_aav_ok = (system_cds_bp / 2.0 + 800) <= AAV_PACKAGING_LIMIT_BP
    # LNP delivers mRNA up to ~10 kb
    lnp_ok = (system_cds_bp + EPEGRNA_AVERAGE_NT + NICKING_SGRNA_NT) <= LNP_MRNA_PACKAGING_LIMIT_NT

    if not single_aav_ok:
        rec = "Single AAV packaging exceeded (6.32 kb > 4.7 kb limit). Must use Dual AAV Split-Intein or SORT-LNP mRNA delivery."
    else:
        rec = "Single AAV delivery packaging feasible."

    return DeliveryPackagingScore(
        system_name="PEmax Prime Editor",
        single_aav_feasible=single_aav_ok,
        dual_aav_split_intein_feasible=dual_aav_ok,
        lnp_mrna_feasible=lnp_ok,
        pemax_size_bp=system_cds_bp,
        limiting_vector_kb=AAV_PACKAGING_LIMIT_BP / 1000.0,
        recommendation=rec,
    )


def build_cftr_f508del_prime_editing_design() -> PrimeEditingDesignResult:
    """Pre-computed, biophysically optimized prime editing strategy for CFTR F508del."""
    # F508del (c.1521_1523delCTT): restores CTT encoding Phe508 in NBD1 interface (PDB 6MSM)
    # PBS: 13 nt complementary to nicked non-target strand, engineered for optimal Tm ~32.9 C
    pbs_seq = "TATCATTGGTGTT"  # 13 nt, Tm = 32.94 C
    pbs_res = score_pbs(pbs_seq)

    # RTT: 15 nt incorporating +CTT in-frame insertion and flanking homology
    rtt_seq = "ACACCAAAGATGATA"  # encodes restored wild-type CTT codon context
    rtt_res = score_rtt(rtt_seq, nick_to_edit_distance_nt=4)

    epeg_res = score_epegRNA_motif("mpknot", linker_sequence="CGCGGTTX")
    nick_res = score_nicking_guide(
        sequence="GTTACTTACTGTCTTAACAG",
        relative_distance_bp=52,
        is_pe3b=True,
    )
    pack_res = evaluate_packaging_limits()

    # Weighted composite score
    composite = (
        (pbs_res.score * 0.30)
        + (rtt_res.score * 0.25)
        + (epeg_res.stability_retention_score * 0.20)
        + (nick_res.score * 0.15)
        + ((1.0 if pack_res.lnp_mrna_feasible else 0.0) * 0.10)
    ) * 100.0

    verdict = "OPTIMAL_FEASIBILITY" if composite >= 85.0 else "FEASIBLE"

    return PrimeEditingDesignResult(
        target_id="T-0184-F508DEL",
        cftr_mutation="F508del (c.1521_1523delCTT)",
        edit_type="In-frame CTT trinucleotide restoration (Phe508 insertion)",
        cryoem_pdb_context="Restores NBD1:ICL4 interdomain packing interface in Cryo-EM PDB 6MSM",
        pbs=pbs_res,
        rtt=rtt_res,
        epegRNA=epeg_res,
        nicking_guide=nick_res,
        packaging=pack_res,
        composite_feasibility_score=round(composite, 2),
        verdict=verdict,
        evidence_citations=[
            "DOI 10.1038/s41587-023-01851-6",
            "DOI 10.1016/j.omtn.2022.04.025",
            "RCSB PDB 6MSM",
        ],
    )


def build_cftr_w1282x_prime_editing_design() -> PrimeEditingDesignResult:
    """Pre-computed, biophysically optimized prime editing strategy for CFTR W1282X."""
    # W1282X (c.3846G>A): converts premature UGA termination codon back to UGG (Trp1282)
    # PBS: 13 nt engineered for optimal Tm ~34.7 C
    pbs_seq = "AGCTCTTGTTTAG"  # 13 nt, Tm = 34.65 C
    pbs_res = score_pbs(pbs_seq)

    # RTT: 16 nt template with A->G correction (converting TGA to TGG)
    rtt_seq = "AAGCCAGACTTTCTTG"  # template strand directing reverse transcription
    rtt_res = score_rtt(rtt_seq, nick_to_edit_distance_nt=3)

    epeg_res = score_epegRNA_motif("tevpre", linker_sequence="CGCGGTTX")
    nick_res = score_nicking_guide(
        sequence="TTCTTCCTCAGTGATGGTGC",
        relative_distance_bp=64,
        is_pe3b=False,
    )
    pack_res = evaluate_packaging_limits()

    composite = (
        (pbs_res.score * 0.30)
        + (rtt_res.score * 0.25)
        + (epeg_res.stability_retention_score * 0.20)
        + (nick_res.score * 0.15)
        + ((1.0 if pack_res.lnp_mrna_feasible else 0.0) * 0.10)
    ) * 100.0

    verdict = "OPTIMAL_FEASIBILITY" if composite >= 85.0 else "FEASIBLE"

    return PrimeEditingDesignResult(
        target_id="T-0184-W1282X",
        cftr_mutation="W1282X (c.3846G>A)",
        edit_type="Single-nucleotide transition A->G (restores Trp1282 codon TGG)",
        cryoem_pdb_context="Restores full-length NBD2 catalytic cassette in Cryo-EM PDB 6MSM",
        pbs=pbs_res,
        rtt=rtt_res,
        epegRNA=epeg_res,
        nicking_guide=nick_res,
        packaging=pack_res,
        composite_feasibility_score=round(composite, 2),
        verdict=verdict,
        evidence_citations=[
            "DOI 10.1038/s41587-023-01851-6",
            "DOI 10.1016/j.omtn.2022.04.025",
            "RCSB PDB 6MSM",
        ],
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Biophysical Feasibility Scorer for CFTR Prime Editing (T-0184)"
    )
    parser.add_argument(
        "--target",
        choices=["F508del", "W1282X", "all"],
        default="all",
        help="Target CFTR mutation to score",
    )
    parser.add_argument("--json", action="store_true", help="Print structured JSON output")
    parser.add_argument("--export", type=str, help="Export JSON results to path")
    args = parser.parse_args()

    results: Dict[str, Any] = {}

    if args.target in ("F508del", "all"):
        f508 = build_cftr_f508del_prime_editing_design()
        results["F508del"] = asdict(f508)

    if args.target in ("W1282X", "all"):
        w1282 = build_cftr_w1282x_prime_editing_design()
        results["W1282X"] = asdict(w1282)

    if args.export:
        out_path = Path(args.export)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with open(out_path, "w", encoding="utf-8") as f:
            json.dump(results, f, indent=2)
        print(f"Exported prime editing designs to {out_path}")

    if args.json:
        print(json.dumps(results, indent=2))
    else:
        print("========================================================================")
        print("CFTR PRIME EDITING BIOPHYSICAL FEASIBILITY AUDIT (T-0184)")
        print("========================================================================")
        for name, data in results.items():
            print(f"Mutation: {data['cftr_mutation']}")
            print(f"  Target: {data['target_id']} | Edit: {data['edit_type']}")
            print(f"  PDB Context: {data['cryoem_pdb_context']}")
            print(f"  PBS Tm: {data['pbs']['melting_temperature_c']} C (Optimal: {data['pbs']['within_optimal_bounds']})")
            print(f"  RTT Length: {data['rtt']['length_nt']} nt (Within constraints: {data['rtt']['within_constraints']})")
            print(f"  epegRNA Motif: {data['epegRNA']['motif_name']} (Engineered: {data['epegRNA']['is_engineered']})")
            print(f"  Packaging: Single AAV={data['packaging']['single_aav_feasible']} | Dual AAV={data['packaging']['dual_aav_split_intein_feasible']} | LNP={data['packaging']['lnp_mrna_feasible']}")
            print(f"  Composite Score: {data['composite_feasibility_score']} / 100.0 -> {data['verdict']}")
            print("------------------------------------------------------------------------")


if __name__ == "__main__":
    main()
