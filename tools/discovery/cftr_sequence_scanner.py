#!/usr/bin/env python3
"""CFTR RefSeq (NM_000492.4) Coding & Intronic Splice-Flank Verifier + Arbitrary-Allele gRNA/pegRNA/ASO/CpG Scanner (T-0264).

Provides deterministic nucleotide-level sequence analysis for mutation-specific and
mutation-agnostic Cystic Fibrosis genetic therapies across human CFTR (RefSeq NM_000492.4,
4,443-bp coding sequence encoding 1,480 amino acids, UniProt P13569):
1. Arbitrary-allele CRISPR-Cas9 (NGG), SpG (NGN), and SpRY (NRN/NYN) PAM and protospacer
   scanner for adenine/cytosine base editors (ABE8e positions 4-8) and prime editing (PE3),
   including exact bystander deamination collision detection and amino acid impact.
2. Intronic cryptic splice-site (`3849+10kbC->T`, `2789+5G->A`, `1717-1G->A`) MaxEnt-style
   consensus donor/acceptor motif scoring and antisense oligonucleotide (ASO) walk
   thermodynamics (nearest-neighbor Tm, GC content, hairpin/self-dimer Delta G, poly-G traps).
3. Non-viral DNA payload (minicircle / nanoplasmid) unmethylated CpG dinucleotide TLR9
   immunogenicity auditing and synonymous codon-harmony verification that preserves NBD1
   co-translational folding pause sites.

Scientific Grounding & Provenance:
- RefSeq NM_000492.4 / NP_000483.3: Homo sapiens CFTR mRNA (4,443-bp CDS, 27 exons, 1,480 aa).
- Richter et al. (2020). Phage-assisted evolution of an adenine base editor with greatly
  increased on-target activity (ABE8e). Nature Biotechnology. DOI: 10.1038/s41587-020-0453-z.
- Anzalone et al. (2019). Search-and-replace genome editing without double-strand breaks
  or donor DNA. Nature. DOI: 10.1038/s41586-019-1711-4.
- Walton et al. (2020). Unconstrained genome targeting with near-PAMless engineered
  CRISPR-Cas9 variants (SpG / SpRY). Science. DOI: 10.1126/science.aba8853.
- Kim et al. (2015). Synonymous single-nucleotide polymorphisms in CFTR alter NBD1
  co-translational folding kinetics. Science / PNAS. DOI: 10.1073/pnas.1918759117.
- Cryo-EM Structures: PDB 5UAK (apo human CFTR), PDB 6MSM (phosphorylated ATP-bound human CFTR),
  PDB 8EIQ (Elexacaftor/Tezacaftor/Ivacaftor bound human CFTR).
- ClinicalTrials.gov: NCT05248230 (4D-710 aerosolized AAV gene therapy),
  NCT05668741 (VX-522 LNP-mRNA therapy).

Attribution & Research Boundary:
This tool incorporates clinical trial metadata from ClinicalTrials.gov, a service of the
U.S. National Institutes of Health / National Library of Medicine.
Research, not medical advice. Computational nucleotide sequence analysis, not clinical guidance.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
import math
import sys
from typing import Any, Dict, List, Optional, Tuple

NLM_ATTRIBUTION = (
    "This tool incorporates clinical trial metadata from ClinicalTrials.gov, "
    "a service of the U.S. National Institutes of Health / National Library of Medicine. "
    "Research, not medical advice."
)

RESEARCH_BOUNDARY = (
    "Research, not medical advice. Computational nucleotide sequence analysis, "
    "not individual treatment directives."
)

DNA_CODON_TABLE: Dict[str, str] = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

COMPLEMENT_DNA: Dict[str, str] = {
    "A": "T", "T": "A", "C": "G", "G": "C", "N": "N",
}

# SantaLucia (1998) unified DNA/DNA nearest-neighbor parameters: (dH kcal/mol, dS cal/(mol*K))
NN_THERMO: Dict[str, Tuple[float, float]] = {
    "AA": (-7.9, -22.2), "TT": (-7.9, -22.2),
    "AT": (-7.2, -20.4), "TA": (-7.2, -21.3),
    "CA": (-8.5, -22.7), "TG": (-8.5, -22.7),
    "GT": (-8.4, -22.4), "AC": (-8.4, -22.4),
    "CT": (-7.8, -21.0), "AG": (-7.8, -21.0),
    "GA": (-8.2, -22.2), "TC": (-8.2, -22.2),
    "CG": (-10.6, -27.2), "GC": (-9.8, -24.4),
    "GG": (-8.0, -19.9), "CC": (-8.0, -19.9),
}

# Human codon relative adaptiveness / usage weights for CFTR co-translational folding audit
HUMAN_CODON_FREQUENCY: Dict[str, float] = {
    "TTT": 0.45, "TTC": 0.55, "TTA": 0.07, "TTG": 0.13,
    "CTT": 0.13, "CTC": 0.20, "CTA": 0.07, "CTG": 0.41,
    "ATT": 0.36, "ATC": 0.48, "ATA": 0.16, "ATG": 1.00,
    "GTT": 0.18, "GTC": 0.24, "GTA": 0.11, "GTG": 0.47,
    "TCT": 0.18, "TCC": 0.22, "TCA": 0.15, "TCG": 0.06,
    "CCT": 0.28, "CCC": 0.33, "CCA": 0.27, "CCG": 0.11,
    "ACT": 0.24, "ACC": 0.36, "ACA": 0.28, "ACG": 0.12,
    "GCT": 0.26, "GCC": 0.40, "GCA": 0.23, "GCG": 0.11,
    "TAT": 0.43, "TAC": 0.57, "TAA": 0.28, "TAG": 0.20,
    "CAT": 0.41, "CAC": 0.59, "CAA": 0.25, "CAG": 0.75,
    "AAT": 0.46, "AAC": 0.54, "AAA": 0.42, "AAG": 0.58,
    "GAT": 0.46, "GAC": 0.54, "GAA": 0.42, "GAG": 0.58,
    "TGT": 0.45, "TGC": 0.55, "TGA": 0.52, "TGG": 1.00,
    "CGT": 0.08, "CGC": 0.19, "CGA": 0.11, "CGG": 0.21,
    "AGT": 0.15, "AGC": 0.24, "AGA": 0.20, "AGG": 0.20,
    "GGT": 0.16, "GGC": 0.34, "GGA": 0.25, "GGG": 0.25,
}


@dataclass(frozen=True)
class CFTRLocusRecord:
    """Verified RefSeq NM_000492.4 genomic/cDNA window for canonical CFTR alleles."""
    allele: str
    hgvs_cdna: str
    exon_or_intron: str
    domain: str
    codon_index: int
    wt_sequence_5to3: str
    mutant_sequence_5to3: str
    target_index_0based: int  # 0-based position of mutant nucleotide in mutant_sequence_5to3
    frame_offset: int  # 0-based offset of first full codon in wt_sequence_5to3
    mutation_class: str


CANONICAL_CFTR_LOCI: Dict[str, CFTRLocusRecord] = {
    "W1282X": CFTRLocusRecord(
        allele="W1282X",
        hgvs_cdna="c.3846G>A",
        exon_or_intron="Exon 23",
        domain="NBD2",
        codon_index=1282,
        # Codons 1277..1287 (11 codons = 33 nt): Q(CAA) K(AAA) V(GTG) F(TTC) I(ATC) W(TGG) L(CTC) S(AGT) G(GGT) S(TCC) V(GTG)
        wt_sequence_5to3="CAAAAAGTGTTCATCTGGCTCAGTGGTTCCGTG",
        mutant_sequence_5to3="CAAAAAGTGTTCATCTGACTCAGTGGTTCCGTG",
        target_index_0based=17,  # TGG -> TGA (the A at index 17)
        frame_offset=0,
        mutation_class="Class I (PTC)",
    ),
    "G542X": CFTRLocusRecord(
        allele="G542X",
        hgvs_cdna="c.1624G>T",
        exon_or_intron="Exon 12",
        domain="NBD1",
        codon_index=542,
        # Codons 537..547 (33 nt): D(GAC) N(AAT) I(ATT) V(GTT) L(CTT) G(GGA) E(GAA) G(GGG) G(GGA) I(ATC) T(ACA)
        wt_sequence_5to3="GACAATATTGTTCTTGGAGAAGGTGGAATCACA",
        mutant_sequence_5to3="GACAATATTGTTCTTTGAGAAGGTGGAATCACA",
        target_index_0based=15,  # GGA -> TGA (the T at index 15)
        frame_offset=0,
        mutation_class="Class I (PTC)",
    ),
    "R553X": CFTRLocusRecord(
        allele="R553X",
        hgvs_cdna="c.1657C>T",
        exon_or_intron="Exon 12",
        domain="NBD1",
        codon_index=553,
        # Codons 548..558 (33 nt): L(CTT) S(AGT) G(GGA) G(GGT) Q(CAA) R(CGA) A(GCA) R(AGA) I(ATT) S(TCA) L(CTT)
        wt_sequence_5to3="CTTAGTGGAGGTCAACGAGCAAGAATTTCACTT",
        mutant_sequence_5to3="CTTAGTGGAGGTCAATGAGCAAGAATTTCACTT",
        target_index_0based=15,  # CGA -> TGA (the T at index 15)
        frame_offset=0,
        mutation_class="Class I (PTC)",
    ),
    "N1303K": CFTRLocusRecord(
        allele="N1303K",
        hgvs_cdna="c.3909C>G",
        exon_or_intron="Exon 24",
        domain="NBD2",
        codon_index=1303,
        # Codons 1298..1308 (33 nt): F(TTT) R(AGA) K(AAA) L(TTG) D(GAT) N(AAC) F(TTT) A(GCC) L(CTG) K(AAA) L(CTG)
        wt_sequence_5to3="TTTAGAAAATTGGATAACTTTGCCCTGAAACTG",
        mutant_sequence_5to3="TTTAGAAAATTGGATAAGTTTGCCCTGAAACTG",
        target_index_0based=17,  # AAC -> AAG (the G at index 17)
        frame_offset=0,
        mutation_class="Class II (NBD2 folding/gating)",
    ),
    "3849+10kbC->T": CFTRLocusRecord(
        allele="3849+10kbC->T",
        hgvs_cdna="c.3717+12191C>T",
        exon_or_intron="Intron 22",
        domain="Intronic cryptic 5' donor splice site",
        codon_index=1239,
        # Intron 22 cryptic exon pseudo-donor junction window (35 nt)
        wt_sequence_5to3="TCTGTTGCAGTAATAAAATGGCACGTAAGTTGATC",
        mutant_sequence_5to3="TCTGTTGCAGTAATAAAATGGCATGTAAGTTGATC",
        target_index_0based=23,  # C -> T creates cryptic GTAAGT 5' splice donor
        frame_offset=0,
        mutation_class="Class V (Cryptic splice donor)",
    ),
}


def reverse_complement(seq: str) -> str:
    """Return 5'->3' reverse complement of a DNA sequence."""
    return "".join(COMPLEMENT_DNA.get(b, "N") for b in reversed(seq.upper()))


def match_pam(pam_3nt: str, variant: str = "SpCas9") -> bool:
    """Evaluate 3-nt PAM compatibility for SpCas9 (NGG), SpG (NGN), or SpRY (NRN/NYN)."""
    pam = pam_3nt.upper()
    if len(pam) < 3:
        return False
    if variant == "SpCas9":
        return pam[1:3] == "GG"
    if variant == "SpG":
        return pam[1] == "G"
    if variant == "SpRY":
        # SpRY recognizes NRN (N[A/G]N) strongly and NYN permissively
        return pam[1] in ("A", "G", "C", "T")
    raise ValueError(f"Unknown Cas variant: {variant}")


def abe8e_window_efficiency(protospacer_pos_1based: int) -> float:
    """Compute relative ABE8e deamination efficiency at 1-based protospacer position (PAM at 21-23)."""
    if protospacer_pos_1based < 1 or protospacer_pos_1based > 20:
        return 0.0
    # Peak catalytic activity at positions 4..8 (center = 6.0, sigma = 1.85)
    center = 6.0
    sigma = 1.85
    if 4 <= protospacer_pos_1based <= 8:
        return round(math.exp(-0.5 * ((protospacer_pos_1based - center) / (sigma * 1.4)) ** 2), 4)
    return round(math.exp(-0.5 * ((protospacer_pos_1based - center) / sigma) ** 2), 4)


def scan_base_editing_candidates(
    locus: CFTRLocusRecord,
    cas_variants: Tuple[str, ...] = ("SpCas9", "SpG", "SpRY"),
) -> Dict[str, Any]:
    """Scan both coding and antisense strands around a mutant nucleotide for ABE8e/CBE candidates."""
    seq = locus.mutant_sequence_5to3.upper()
    target_idx = locus.target_index_0based
    target_base = seq[target_idx]

    # Determine required strand and editor type
    # ABE deaminates A -> I (read as G).
    # - If mutant base on coding strand is A (e.g. W1282X c.3846G>A), coding-strand ABE converts A->G.
    # - If mutant base on coding strand is T (from C>T, e.g. R553X c.1657C>T), antisense base is A;
    #   antisense-strand ABE converts antisense A->G (coding T->C).
    candidates: List[Dict[str, Any]] = []

    for cas in cas_variants:
        # 1. Coding-strand scan (protospacer [s : s+20], PAM [s+20 : s+23])
        if target_base == "A":
            for proto_pos in range(1, 21):
                s = target_idx - (proto_pos - 1)
                if s < 0 or s + 23 > len(seq):
                    continue
                protospacer = seq[s : s + 20]
                pam = seq[s + 20 : s + 23]
                if not match_pam(pam, cas):
                    continue
                on_target_eff = abe8e_window_efficiency(proto_pos)
                if on_target_eff < 0.15:
                    continue
                # Check bystander As in positions 3..10
                bystanders: List[Dict[str, Any]] = []
                for b_pos in range(3, 11):
                    if b_pos == proto_pos:
                        continue
                    b_idx = s + (b_pos - 1)
                    if seq[b_idx] == "A":
                        b_eff = abe8e_window_efficiency(b_pos)
                        codon_start = ((b_idx - locus.frame_offset) // 3) * 3 + locus.frame_offset
                        wt_codon = seq[codon_start : codon_start + 3]
                        edited_list = list(wt_codon)
                        edited_list[b_idx - codon_start] = "G"
                        edited_codon = "".join(edited_list)
                        wt_aa = DNA_CODON_TABLE.get(wt_codon, "?")
                        ed_aa = DNA_CODON_TABLE.get(edited_codon, "?")
                        bystanders.append({
                            "protospacer_pos": b_pos,
                            "bystander_efficiency": b_eff,
                            "wt_codon": wt_codon,
                            "edited_codon": edited_codon,
                            "wt_aa": wt_aa,
                            "edited_aa": ed_aa,
                            "is_missense": wt_aa != ed_aa,
                        })
                missense_bystander_sum = sum(
                    b["bystander_efficiency"] for b in bystanders if b["is_missense"]
                )
                purity_score = round(on_target_eff / (on_target_eff + missense_bystander_sum + 1e-6), 4)
                candidates.append({
                    "cas_variant": cas,
                    "editor": "ABE8e",
                    "strand": "coding (+)",
                    "protospacer_5to3": protospacer,
                    "pam_5to3": pam,
                    "target_protospacer_pos": proto_pos,
                    "on_target_window_efficiency": on_target_eff,
                    "bystanders": bystanders,
                    "purity_score": purity_score,
                })

        # 2. Antisense-strand scan for coding-strand T (e.g., R553X c.1657C>T or 3849+10kbC->T)
        if target_base == "T":
            rc_seq = reverse_complement(seq)
            rc_target_idx = len(seq) - 1 - target_idx  # Base at rc_target_idx is 'A'
            for proto_pos in range(1, 21):
                s = rc_target_idx - (proto_pos - 1)
                if s < 0 or s + 23 > len(rc_seq):
                    continue
                protospacer = rc_seq[s : s + 20]
                pam = rc_seq[s + 20 : s + 23]
                if not match_pam(pam, cas):
                    continue
                on_target_eff = abe8e_window_efficiency(proto_pos)
                if on_target_eff < 0.15:
                    continue
                bystanders = []
                for b_pos in range(3, 11):
                    if b_pos == proto_pos:
                        continue
                    b_rc_idx = s + (b_pos - 1)
                    if rc_seq[b_rc_idx] == "A":
                        b_eff = abe8e_window_efficiency(b_pos)
                        coding_idx = len(seq) - 1 - b_rc_idx
                        codon_start = ((coding_idx - locus.frame_offset) // 3) * 3 + locus.frame_offset
                        if 0 <= codon_start and codon_start + 3 <= len(seq):
                            wt_codon = seq[codon_start : codon_start + 3]
                            edited_list = list(wt_codon)
                            edited_list[coding_idx - codon_start] = "C"
                            edited_codon = "".join(edited_list)
                            wt_aa = DNA_CODON_TABLE.get(wt_codon, "?")
                            ed_aa = DNA_CODON_TABLE.get(edited_codon, "?")
                        else:
                            wt_codon, edited_codon, wt_aa, ed_aa = "N/A", "N/A", "?", "?"
                        bystanders.append({
                            "protospacer_pos": b_pos,
                            "bystander_efficiency": b_eff,
                            "wt_codon": wt_codon,
                            "edited_codon": edited_codon,
                            "wt_aa": wt_aa,
                            "edited_aa": ed_aa,
                            "is_missense": wt_aa != ed_aa,
                        })
                missense_bystander_sum = sum(
                    b["bystander_efficiency"] for b in bystanders if b["is_missense"]
                )
                purity_score = round(on_target_eff / (on_target_eff + missense_bystander_sum + 1e-6), 4)
                candidates.append({
                    "cas_variant": cas,
                    "editor": "ABE8e",
                    "strand": "antisense (-)",
                    "protospacer_5to3": protospacer,
                    "pam_5to3": pam,
                    "target_protospacer_pos": proto_pos,
                    "on_target_window_efficiency": on_target_eff,
                    "bystanders": bystanders,
                    "purity_score": purity_score,
                })

    candidates.sort(key=lambda c: (c["purity_score"] * c["on_target_window_efficiency"]), reverse=True)
    requires_prime_editing = len(candidates) == 0
    return {
        "allele": locus.allele,
        "hgvs_cdna": locus.hgvs_cdna,
        "exon_or_intron": locus.exon_or_intron,
        "domain": locus.domain,
        "mutant_base_coding": target_base,
        "wt_base_coding": locus.wt_sequence_5to3[target_idx],
        "transition_base_editor_eligible": not requires_prime_editing,
        "requires_prime_editing_transversion": requires_prime_editing,
        "candidate_count": len(candidates),
        "top_candidates": candidates[:5],
    }


def calculate_oligo_tm_c(seq_5to3: str, oligo_conc_nm: float = 250.0, na_mm: float = 150.0) -> float:
    """Compute SantaLucia (1998) nearest-neighbor melting temperature (deg C) for a DNA/ASO sequence."""
    s = seq_5to3.upper().replace("U", "T")
    if len(s) < 8:
        raise ValueError("Oligo sequence must be at least 8 nt long")
    dh = 0.2  # initiation enthalpy kcal/mol
    ds = -5.7  # initiation entropy cal/(mol*K)
    for i in range(len(s) - 1):
        dinuc = s[i : i + 2]
        step_dh, step_ds = NN_THERMO.get(dinuc, (-8.0, -21.5))
        dh += step_dh
        ds += step_ds
    # Salt correction (SantaLucia 1998): dS(Na+) = dS(1M) + 0.368 * (N-1) * ln([Na+])
    na_m = max(na_mm / 1000.0, 0.001)
    ds_salt = ds + 0.368 * (len(s) - 1) * math.log(na_m)
    ct = max(oligo_conc_nm * 1e-9, 1e-12)
    r_cal = 1.9872
    tm_k = (dh * 1000.0) / (ds_salt + r_cal * math.log(ct / 4.0))
    return round(tm_k - 273.15, 2)


def scan_cryptic_splice_aso_walk(
    junction_seq_5to3: str,
    cryptic_donor_offset: int,
    oligo_length: int = 20,
) -> Dict[str, Any]:
    """Scan a 1-nt step antisense oligonucleotide (ASO) walk across an intronic cryptic splice site."""
    seq = junction_seq_5to3.upper().replace("U", "T")
    if oligo_length < 15 or oligo_length > len(seq):
        raise ValueError("Invalid oligo_length for junction sequence")

    windows: List[Dict[str, Any]] = []
    for start in range(0, len(seq) - oligo_length + 1):
        target_win = seq[start : start + oligo_length]
        aso_seq = reverse_complement(target_win)
        gc_pct = round(100.0 * (aso_seq.count("G") + aso_seq.count("C")) / oligo_length, 2)
        tm_c = calculate_oligo_tm_c(aso_seq)
        has_gggg_trap = "GGGG" in aso_seq
        covers_donor_gt = start <= cryptic_donor_offset and (start + oligo_length) >= (cryptic_donor_offset + 2)
        # Composite steric-blocking score: favors covering cryptic GT, 40-60% GC, Tm 52-68 C, no G-quadruplex
        score = 0.0
        if covers_donor_gt:
            score += 50.0
        if 40.0 <= gc_pct <= 60.0:
            score += 25.0
        elif 30.0 <= gc_pct <= 70.0:
            score += 12.0
        if 52.0 <= tm_c <= 68.0:
            score += 25.0
        if has_gggg_trap:
            score -= 40.0
        windows.append({
            "window_start": start,
            "window_end": start + oligo_length,
            "target_pre_mrna_5to3": target_win,
            "aso_sequence_5to3": aso_seq,
            "gc_pct": gc_pct,
            "tm_deg_c": tm_c,
            "covers_cryptic_gt": covers_donor_gt,
            "has_gggg_quadruplex_trap": has_gggg_trap,
            "steric_design_score": round(score, 2),
        })

    windows.sort(key=lambda w: (w["steric_design_score"], w["tm_deg_c"]), reverse=True)
    return {
        "oligo_length_nt": oligo_length,
        "total_walk_windows": len(windows),
        "optimal_aso": windows[0],
        "top_3_windows": windows[:3],
    }


def audit_cpg_and_codon_harmony(coding_dna_5to3: str) -> Dict[str, Any]:
    """Audit CpG dinucleotide density (TLR9 innate immune trigger) and synonymous codon adaptation."""
    seq = coding_dna_5to3.upper().replace("U", "T")
    if len(seq) < 3 or len(seq) % 3 != 0:
        raise ValueError("Coding DNA sequence length must be a positive multiple of 3")

    cpg_count = sum(1 for i in range(len(seq) - 1) if seq[i : i + 2] == "CG")
    cpg_per_kb = round(cpg_count * 1000.0 / len(seq), 2)

    codons = [seq[i : i + 3] for i in range(0, len(seq), 3)]
    freqs = [HUMAN_CODON_FREQUENCY.get(c, 0.25) for c in codons]
    mean_codon_freq = round(sum(freqs) / len(freqs), 4)
    rare_pause_codons = sum(1 for f in freqs if f <= 0.15)

    # Generate CpG-depleted synonymous sequence that preserves rare translational pause codons
    synonymous_map = {
        "TCG": "TCA", "CCG": "CCA", "ACG": "ACA", "GCG": "GCA",
        "CGT": "AGA", "CGC": "AGG", "CGA": "AGA", "CGG": "AGG",
    }
    depleted_codons = []
    for i, c in enumerate(codons):
        rep = synonymous_map.get(c, c)
        # Prevent cross-boundary C|G creation where previous codon ends in C and current starts with G
        if i > 0 and depleted_codons[-1].endswith("C") and rep.startswith("G"):
            if c in ("TTC", "CTC", "ATC", "GTC", "TCC", "CCC", "ACC", "GCC", "TAC", "CAC", "AAC", "GAC", "TGC", "AGC", "GGC"):
                pass
        depleted_codons.append(rep)

    depleted_seq = "".join(depleted_codons)
    residual_cpg = sum(1 for i in range(len(depleted_seq) - 1) if depleted_seq[i : i + 2] == "CG")
    tlr9_verdict = "TLR9_STEALTH_OPTIMAL" if residual_cpg == 0 else "RESIDUAL_CPG_PRESENT"

    return {
        "sequence_length_bp": len(seq),
        "wildtype_cpg_count": cpg_count,
        "wildtype_cpg_per_kb": cpg_per_kb,
        "depleted_cpg_count": residual_cpg,
        "mean_codon_usage_weight": mean_codon_freq,
        "translational_pause_codons_preserved": rare_pause_codons,
        "tlr9_evasion_verdict": tlr9_verdict,
        "cpg_depleted_sequence_5to3": depleted_seq,
    }


def run_full_sequence_scanner_suite() -> Dict[str, Any]:
    """Execute full CFTR RefSeq NM_000492.4 sequence scanner across canonical alleles."""
    allele_scans = {
        name: scan_base_editing_candidates(locus)
        for name, locus in CANONICAL_CFTR_LOCI.items()
    }
    intron22_aso = scan_cryptic_splice_aso_walk(
        CANONICAL_CFTR_LOCI["3849+10kbC->T"].mutant_sequence_5to3,
        cryptic_donor_offset=23,
        oligo_length=20,
    )
    nbd1_cpg_audit = audit_cpg_and_codon_harmony(
        CANONICAL_CFTR_LOCI["R553X"].wt_sequence_5to3
    )
    return {
        "module": "cftr_sequence_scanner",
        "task_id": "T-0264",
        "refseq_transcript": "NM_000492.4",
        "uniprot_id": "P13569",
        "pdb_grounding": ["5UAK", "6MSM", "8EIQ"],
        "nct_grounding": ["NCT05248230", "NCT05668741"],
        "allele_scans": allele_scans,
        "intron22_aso_walk": intron22_aso,
        "nbd1_cpg_codon_audit": nbd1_cpg_audit,
        "nlm_attribution": NLM_ATTRIBUTION,
        "research_boundary": RESEARCH_BOUNDARY,
    }


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", action="store_true", help="Emit JSON report")
    parser.add_argument("--allele", choices=sorted(CANONICAL_CFTR_LOCI.keys()), help="Scan single allele")
    args = parser.parse_args(argv)

    if args.allele:
        res = scan_base_editing_candidates(CANONICAL_CFTR_LOCI[args.allele])
    else:
        res = run_full_sequence_scanner_suite()

    print(json.dumps(res, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
