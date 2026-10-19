#!/usr/bin/env python3
"""Base Editing (ABE8e / CBE) Bystander Deamination Window & PAM Accessibility Model for CFTR (T-0219).

Provides deterministic biophysical modeling for:
1. SpCas9 / SpRY / SpG / PAM variant accessibility (NGG, NRN, NYN, NGN).
2. Deamination activity window within the protospacer (Gaussian deamination probability profile).
3. On-target base editing efficiency across refractory CFTR mutations (W1282X, R553X, G542X, 2789+5G->A).
4. Off-target bystander edit kinetics (introducing missense perturbations into adjacent codons).
5. Cryo-EM PDB 6MSM structural stability impact (DeltaDelta G_fold) and bystander purity ratio (P).

Scientific Grounding & Provenance:
- Richter et al. (2020). Phage-assisted evolution of an adenine base editor with greatly increased on-target activity.
  Nature Biotechnology, 38(7):883-891. DOI: 10.1038/s41587-020-0453-z.
- Walton et al. (2020). Unconstrained genome targeting with near-PAMless engineered CRISPR-Cas9 variants.
  Science, 368(6488):290-296. DOI: 10.1126/science.aba8853.
- Koblan et al. (2018). Improving cytidine and adenine base editors by expression optimization and ancestral
  sequence reconstruction. Nature Biotechnology, 36(9):843-846. DOI: 10.1038/nbt.4172.
- RCSB PDB 6MSM: Human CFTR dephosphorylated, ATP-free conformation resolved by Cryo-EM at 3.2 A.
- RCSB PDB 8EIQ: Human CFTR phosphorylated, Trikafta-bound active conformation at 3.0 A.
- ClinicalTrials.gov NCT06237335: Clinical evaluation of investigational base and prime editing therapeutic agents
  in genetic respiratory disease.

NLM Attribution & Research Boundary:
This tool incorporates clinical trial metadata from ClinicalTrials.gov, a service of the
U.S. National Institutes of Health / National Library of Medicine.
Research, not medical advice. Computational biophysical model, not individual treatment directives.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass, field
import json
import math
from pathlib import Path
import re
import sys
from typing import Any, Dict, List, Optional, Tuple

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
CFTR_DATA_PATH = REPO_ROOT / "tools" / "cftr-model" / "web" / "cftr-data.js"

NLM_ATTRIBUTION = (
    "This tool incorporates clinical trial metadata from ClinicalTrials.gov, "
    "a service of the U.S. National Institutes of Health / National Library of Medicine. "
    "Research, not medical advice."
)

RESEARCH_BOUNDARY = (
    "Research, not medical advice. Computational biophysical model, "
    "not individual treatment directives."
)

# Standard Genetic Code (RNA Codon -> Amino Acid)
CODON_TABLE: Dict[str, str] = {
    "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
    "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
    "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M",
    "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
    "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
    "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "UAU": "Y", "UAC": "Y", "UAA": "*", "UAG": "*",
    "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "UGU": "C", "UGC": "C", "UGA": "*", "UGG": "W",
    "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

DNA_CODON_TABLE: Dict[str, str] = {
    codon.replace("U", "T"): aa for codon, aa in CODON_TABLE.items()
}

COMPLEMENT_MAP: Dict[str, str] = {
    "A": "T", "T": "A", "C": "G", "G": "C",
    "a": "t", "t": "a", "c": "g", "g": "c",
    "N": "N", "R": "Y", "Y": "R", "M": "K", "K": "M",
}


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    return "".join(COMPLEMENT_MAP.get(b, "N") for b in reversed(seq))


@dataclass(frozen=True)
class EditorProfile:
    """Biophysical deamination parameters for base editing deaminase architecture."""
    name: str
    editor_type: str  # "ABE" or "CBE"
    target_base: str  # "A" for ABE, "C" for CBE
    product_base: str  # "G" for ABE, "T" for CBE
    center_pos: float  # p_0 (typically 5.5)
    sigma: float  # window width parameter
    window_start: int  # 5' active boundary (1-indexed within protospacer)
    window_end: int  # 3' active boundary (1-indexed within protospacer)
    description: str


EDITORS: Dict[str, EditorProfile] = {
    "ABE8e": EditorProfile(
        name="ABE8e",
        editor_type="ABE",
        target_base="A",
        product_base="G",
        center_pos=5.5,
        sigma=1.8,
        window_start=3,
        window_end=9,
        description="Phage-assisted evolved TadA-8e adenine base editor with broadened activity window (Richter et al. 2020)",
    ),
    "ABE7.10": EditorProfile(
        name="ABE7.10",
        editor_type="ABE",
        target_base="A",
        product_base="G",
        center_pos=5.5,
        sigma=1.2,
        window_start=4,
        window_end=8,
        description="Canonical TadA-7.10 adenine base editor with narrow deamination window (Gaudelli et al. 2017)",
    ),
    "AncBE4max": EditorProfile(
        name="AncBE4max",
        editor_type="CBE",
        target_base="C",
        product_base="T",
        center_pos=5.5,
        sigma=1.3,
        window_start=4,
        window_end=8,
        description="Ancestral sequence reconstructed cytidine base editor (Koblan et al. 2018)",
    ),
}


@dataclass(frozen=True)
class CasVariantProfile:
    """CRISPR-Cas PAM recognition and affinity profile."""
    name: str
    canonical_pam: str
    pam_regex: str
    affinity_map: Dict[str, float]
    description: str


CAS_VARIANTS: Dict[str, CasVariantProfile] = {
    "SpCas9": CasVariantProfile(
        name="SpCas9",
        canonical_pam="NGG",
        pam_regex=r"^[ACGT]GG$",
        affinity_map={"NGG": 1.0, "NAG": 0.15},
        description="Wild-type Streptococcus pyogenes Cas9 with strict NGG PAM requirement",
    ),
    "SpRY": CasVariantProfile(
        name="SpRY",
        canonical_pam="NRN",
        pam_regex=r"^[ACGT][AG][ACGT]$|^[ACGT][CT][ACGT]$",
        affinity_map={"NRN": 0.85, "NYN": 0.65, "NNN": 0.50},
        description="Engineered near-PAMless SpCas9 recognizing NRN and NYN motifs (Walton et al. 2020)",
    ),
    "SpG": CasVariantProfile(
        name="SpG",
        canonical_pam="NGN",
        pam_regex=r"^[ACGT]G[ACGT]$",
        affinity_map={"NGN": 0.90, "NNN": 0.10},
        description="Engineered SpCas9 variant with expanded NGN PAM targeting (Walton et al. 2020)",
    ),
}


@dataclass
class BystanderMutation:
    """Biophysical and structural annotation of an off-target bystander deamination event."""
    position_in_protospacer: int  # 1 to 20
    original_base: str
    edited_base: str
    deamination_weight: float  # W_E(p)
    editing_efficiency: float  # W_E(p) * k_PAM * Lambda
    codon_residue_index: Optional[int]
    original_codon: Optional[str]
    edited_codon: Optional[str]
    original_aa: Optional[str]
    edited_aa: Optional[str]
    delta_delta_g_fold_kcal_per_mol: float
    classification: str  # "synonymous", "benign", "tolerated", "destabilizing"
    structural_rationale: str


@dataclass
class GuideRNAResult:
    """Complete evaluation result for a candidate base editing guide RNA."""
    grna_id: str
    allele_name: str
    strand: str  # "coding" or "non_coding"
    protospacer_sequence: str  # 20 nt
    pam_sequence: str  # 3 nt
    cas_variant: str
    editor_name: str
    target_base_position: int  # 1 to 20
    deamination_weight: float  # W_E(p_t)
    pam_affinity: float  # k_PAM
    chromatin_factor: float  # Lambda_chromatin
    on_target_efficiency: float  # Efficiency(A_t)
    bystander_count: int
    bystanders: List[BystanderMutation]
    bystander_purity_ratio: float  # P
    max_destabilizing_ddg: float
    overall_feasibility_score: float  # 0 to 100
    curative_verdict: str
    clinical_trial_reference: str
    cryoem_context: str


@dataclass
class AlleleGenomicProfile:
    """Genomic locus and biophysical parameters for priority CFTR mutations."""
    allele_name: str
    cdna_name: str
    mutation_type: str  # "nonsense", "canonical_splice"
    exon_intron: str
    wildtype_base: str
    mutant_base: str
    wt_codon: Optional[str]
    mutant_codon: Optional[str]
    residue_number: Optional[int]
    coding_sequence_context: str  # 5'->3' coding strand genomic context
    mutation_index_in_context: int  # 0-indexed position of mutant base in context
    target_strand_to_edit: str  # "coding" or "non_coding"
    target_base: str  # "A" for ABE
    desired_base: str  # "G" for ABE
    requires_transversion_editor: bool
    pdb_domain: str
    pdb_6msm_context: str


# Priority CFTR Allele Knowledge Catalog
CFTR_ALLELE_CATALOG: Dict[str, AlleleGenomicProfile] = {
    "W1282X": AlleleGenomicProfile(
        allele_name="W1282X",
        cdna_name="c.3846G>A",
        mutation_type="nonsense",
        exon_intron="Exon 23",
        wildtype_base="G",
        mutant_base="A",
        wt_codon="TGG",
        mutant_codon="TGA",
        residue_number=1282,
        # Flanking coding sequence around residue 1282 (aa 1272 to 1294)
        # Residues: 1272(Ala) 1273(Lys) 1274(Thr) 1275(Arg) 1276(Leu) 1277(Val) 1278(Gly) 1279(Leu) 1280(Ile) 1281(Glu)
        # 1282(Mutant:TGA) 1283(Gln) 1284(Ser) 1285(Thr) 1286(Arg) 1287(Gln) 1288(Arg) 1289(Leu) 1290(Phe) 1291(Gly) 1292(Leu) 1293(Met) 1294(Thr)
        coding_sequence_context="GCTAAGACTAGGTTGGTTGGTCTTATTGAATGACAAAGCACTAGGCAACGGTTGGGCTTGATGACT",
        mutation_index_in_context=28,  # 'A' at index 28 in context (...ATT GAA T-G-[A] CAA...)
        target_strand_to_edit="coding",
        target_base="A",
        desired_base="G",
        requires_transversion_editor=False,
        pdb_domain="NBD2 (ABC transporter 2, residues 1210-1443)",
        pdb_6msm_context="Cryo-EM PDB 6MSM: Codon 1282 C-alpha [-11.7, -24.6, -22.1] in NBD2 catalytic cassette; deletion traps CFTR in ERAD",
    ),
    "R553X": AlleleGenomicProfile(
        allele_name="R553X",
        cdna_name="c.1657C>T",
        mutation_type="nonsense",
        exon_intron="Exon 12",
        wildtype_base="C",
        mutant_base="T",
        wt_codon="CGA",
        mutant_codon="TGA",
        residue_number=553,
        # Flanking coding sequence around residue 553 (aa 545 to 566)
        # 545(Pro:CCA) 546(Glu:GAA) 547(Leu:CTT) 548(Arg:CGC) 549(Gly:GGT) 550(Ala:GCA) 551(Ile:ATC) 552(Lys:AAG)
        # 553(Mutant:TGA) 554(Gln:CAG) 555(Ala:GCT) 556(Arg:CGT) 557(Arg:CGT) 558(Phe:TTT) 559(Ser:TCT) 560(Tyr:TAT)
        # 561(Ala:GCA) 562(Gly:GGA) 563(Gln:CAA) 564(Leu:CTT) 565(Ser:AGT) 566(His:CAT)
        coding_sequence_context="CCAGAACTTCGCGGTGCAATCAAGTGACAGGCTCGTCGTTTTTCTTATGCAGGACAACTTAGTCAT",
        mutation_index_in_context=24,  # 'T' at index 24 in context ([T]-G-A)
        # On non-coding (antisense) strand, mutant base is 'A' (complement of T1657).
        # ABE deamination of this A restores G on non-coding strand, yielding wild-type CGA (Arg553).
        target_strand_to_edit="non_coding",
        target_base="A",
        desired_base="G",
        requires_transversion_editor=False,
        pdb_domain="NBD1 (ABC transporter 1, residues 423-646)",
        pdb_6msm_context="Cryo-EM PDB 6MSM: Codon 553 C-alpha [-7.3, 6.2, -25.2] in NBD1 alpha-subdomain; restores complete NBD1 fold",
    ),
    "G542X": AlleleGenomicProfile(
        allele_name="G542X",
        cdna_name="c.1624G>T",
        mutation_type="nonsense",
        exon_intron="Exon 12",
        wildtype_base="G",
        mutant_base="T",
        wt_codon="GGA",
        mutant_codon="TGA",
        residue_number=542,
        # Flanking coding sequence around residue 542 (aa 535 to 552)
        # 535(Ile:ATT) 536(Leu:CTC) 537(Thr:ACA) 538(Glu:GAA) 539(Asn:AAC) 540(Ile:ATC) 541(Gln:CAG)
        # 542(Mutant:TGA) 543(Leu:CTA) 544(Leu:CTT) 545(Pro:CCA) 546(Glu:GAA) 547(Leu:CTT) 548(Arg:CGC)
        # 549(Gly:GGT) 550(Ala:GCA) 551(Ile:ATC) 552(Lys:AAG)
        coding_sequence_context="ATTCTCACAGAAAACATCCAGTGACCTCTTTCCAGAACTTCGCGGTGCAATCAAG",
        mutation_index_in_context=21,  # 'T' at index 21 in context ([T]-G-A)
        target_strand_to_edit="coding",
        target_base="A",  # If ABE is applied to coding strand TGA, it converts A to G (TGG = Trp, missense!)
        desired_base="G",
        requires_transversion_editor=True,  # Reversal of c.1624G>T requires transversion editor or prime editing
        pdb_domain="NBD1 core beta-sheet turn",
        pdb_6msm_context="Cryo-EM PDB 6MSM: Codon 542 C-alpha [-10.6, 4.1, -13.6]; tight structural turn refractory to bulky aromatic missense",
    ),
    "2789+5G->A": AlleleGenomicProfile(
        allele_name="2789+5G->A",
        cdna_name="c.2657+5G>A",
        mutation_type="canonical_splice",
        exon_intron="Intron 16",
        wildtype_base="G",
        mutant_base="A",
        wt_codon=None,
        mutant_codon=None,
        residue_number=None,
        # Splice donor junction: Exon 16 (3' end: ...CAG) | Intron 16 (5' donor: gt aat t...)
        # WT donor: GTAAGT -> Mutant donor: GTAATT (position +5 is A)
        # Sequence: Exon 16 (...AATACTTCCTTATTGTTCAATTTCAG) | Intron 16 (GTAATTTTCTTGAGAAACTAAGTGATTTTTTAAA)
        coding_sequence_context="AATACTTCCTTATTGTTCAATTTCAGGTAATTTTCTTGAGAAACTAAGTGATTTTTTAAA",
        mutation_index_in_context=29,  # 'A' at index 29 (GTAA-[A]-T / +5 relative to exon boundary at index 24)
        target_strand_to_edit="coding",
        target_base="A",
        desired_base="G",
        requires_transversion_editor=False,
        pdb_domain="Intron 16 Splice Donor (Exon 16 exclusion rescue)",
        pdb_6msm_context="Cryo-EM PDB 6MSM: Restores canonical mRNA splicing to generate full-length 1480-aa CFTR in Cryo-EM 6MSM",
    ),
}


def compute_editing_window(position: int, editor: EditorProfile) -> float:
    """Calculate the asymmetric Gaussian deamination probability profile W_E(p).

    Equation:
      W_E(p) = exp( - (p - p_0)^2 / (2 * sigma^2) )
    where p is the 1-indexed position within the 20-nt protospacer (1 is 5' distal, 20 is adjacent to PAM).
    """
    if position < 1 or position > 20:
        return 0.0
    exponent = -((position - editor.center_pos) ** 2) / (2.0 * (editor.sigma ** 2))
    return math.exp(exponent)


def score_pam_affinity(pam_seq: str, cas_variant: str) -> float:
    """Evaluate CRISPR-Cas variant PAM binding kinetics factor k_PAM.

    Args:
      pam_seq: 3-nucleotide PAM sequence immediately 3' of the 20-nt protospacer.
      cas_variant: Name of Cas variant ("SpCas9", "SpRY", "SpG").
    """
    pam = pam_seq.upper().strip()
    if len(pam) != 3:
        return 0.0

    profile = CAS_VARIANTS.get(cas_variant)
    if not profile:
        return 0.0

    # SpCas9: NGG strict, NAG minor
    if cas_variant == "SpCas9":
        if pam[1:] == "GG":
            return profile.affinity_map.get("NGG", 1.0)
        elif pam[1:] == "AG":
            return profile.affinity_map.get("NAG", 0.15)
        return 0.0

    # SpRY: near-PAMless (NRN >= 0.85, NYN >= 0.65)
    elif cas_variant == "SpRY":
        middle_base = pam[1]
        if middle_base in ("A", "G"):  # R
            return profile.affinity_map.get("NRN", 0.85)
        elif middle_base in ("C", "T"):  # Y
            return profile.affinity_map.get("NYN", 0.65)
        return profile.affinity_map.get("NNN", 0.50)

    # SpG: NGN expanded PAM
    elif cas_variant == "SpG":
        if pam[1] == "G":
            return profile.affinity_map.get("NGN", 0.90)
        return profile.affinity_map.get("NNN", 0.10)

    return 0.0


def calculate_bystander_purity_ratio(on_target_efficiency: float, bystander_weights: List[float]) -> float:
    """Compute the bystander purity ratio P.

    Equation:
      P = Efficiency(A_t) / ( Efficiency(A_t) + sum_{b in B} W_E(p_b) )
    When no bystander adenines are present in the protospacer, P = 1.0 (100% purity).
    """
    bystander_sum = sum(bystander_weights)
    denominator = on_target_efficiency + bystander_sum
    if denominator <= 0.0:
        return 0.0
    return on_target_efficiency / denominator


def evaluate_structural_perturbation(
    allele: AlleleGenomicProfile,
    residue_num: Optional[int],
    orig_aa: Optional[str],
    edited_aa: Optional[str],
) -> Tuple[float, str, str]:
    """Annotate structural impact and DeltaDelta G_fold (kcal/mol) in Cryo-EM PDB 6MSM.

    Classifies missense outcomes:
      - synonymous: DeltaDelta G = 0.0 -> "synonymous"
      - benign: DeltaDelta G <= 0.5 kcal/mol -> "benign"
      - tolerated: 0.5 < DeltaDelta G <= 2.0 kcal/mol -> "tolerated"
      - destabilizing: DeltaDelta G > 2.0 kcal/mol -> "destabilizing"
    """
    if residue_num is None or orig_aa is None or edited_aa is None:
        # Non-coding (e.g. 2789+5G->A in Intron 16)
        return 0.0, "benign", "Intronic splice junction; zero coding missense perturbation."

    if orig_aa == edited_aa:
        return 0.0, "synonymous", "Synonymous codon transition; wild-type polypeptide sequence preserved."

    # Specific residue-level Cryo-EM biophysical rules
    # 1. Trp1282X (NBD2 catalytic core)
    if residue_num == 1282:
        if edited_aa == "W":
            return 0.0, "synonymous", "Restores native wild-type Trp1282; complete NBD2 cassette restored."
        elif edited_aa == "*":
            return 8.5, "destabilizing", "Premature stop codon persists (CFTR_1281 truncation)."
        elif edited_aa in ("R", "K"):
            return 2.8, "destabilizing", "Aromatic-to-basic substitution in NBD2 core disrupts ATP Site 2 interface."
        else:
            return 2.2, "destabilizing", f"Non-native missense p.Trp1282{edited_aa} destabilizes NBD2 fold."

    # 2. Glu1281 (adjacent to 1282)
    elif residue_num == 1281:
        if orig_aa == "E" and edited_aa == "G":
            return 1.4, "tolerated", "Acidic-to-neutral Gly substitution at NBD2 surface loop."
        elif orig_aa == "E" and edited_aa == "K":
            return 2.5, "destabilizing", "Charge reversal (Glu->Lys) disrupts local electrostatic network."
        return 1.2, "tolerated", f"Bystander substitution p.Glu1281{edited_aa} at NBD2 loop."

    # 3. Gln1283 (adjacent to 1282)
    elif residue_num == 1283:
        if orig_aa == "Q" and edited_aa == "R":
            return 0.8, "tolerated", "Polar-to-basic Gln->Arg substitution in solvent-exposed NBD2 loop."
        return 1.1, "tolerated", f"Bystander substitution p.Gln1283{edited_aa}."

    # 4. Arg553X (NBD1 alpha-subdomain)
    elif residue_num == 553:
        if edited_aa == "R":
            return 0.0, "synonymous", "Restores native wild-type Arg553; complete NBD1 alpha-subdomain restored."
        elif edited_aa == "*":
            return 9.2, "destabilizing", "Premature stop codon persists (CFTR_552 truncation)."
        elif edited_aa == "G":
            return 2.7, "destabilizing", "Loss of basic Arg553 sidechain destabilizes NBD1 core."
        return 2.1, "destabilizing", f"Non-native missense p.Arg553{edited_aa} in NBD1 alpha-subdomain."

    # 5. Lys552 (adjacent to 553)
    elif residue_num == 552:
        if orig_aa == "K" and edited_aa == "R":
            return 0.3, "benign", "Conservative Lys->Arg basic substitution on NBD1 surface."
        return 1.5, "tolerated", f"Bystander substitution p.Lys552{edited_aa}."

    # 6. Gly542X (NBD1 tight turn)
    elif residue_num == 542:
        if edited_aa == "G":
            return 0.0, "synonymous", "Restores native wild-type Gly542 in tight NBD1 turn."
        elif edited_aa == "W":
            return 3.6, "destabilizing", "Bulky aromatic Trp542 introduces severe steric clash in NBD1 core."
        elif edited_aa == "R":
            return 3.2, "destabilizing", "Bulky basic Arg542 disrupts NBD1 beta-turn."
        return 2.9, "destabilizing", f"Non-native missense p.Gly542{edited_aa}."

    # Generic empirical substitution model based on CFTR Cryo-EM secondary structure
    # Proline disruption or charge inversion in folded core
    if edited_aa == "P":
        return 2.6, "destabilizing", "Helix-breaking Proline substitution disrupts secondary structure."
    elif (orig_aa in "DE" and edited_aa in "KR") or (orig_aa in "KR" and edited_aa in "DE"):
        return 2.4, "destabilizing", "Charge reversal disrupts tertiary electrostatic salt bridges."
    elif (orig_aa in "ILVMF" and edited_aa in "KRDE") or (orig_aa in "KRDE" and edited_aa in "ILVMF"):
        return 2.1, "destabilizing", "Hydrophobic/hydrophilic transition disrupts core packing."
    else:
        return 0.9, "tolerated", f"Tolerated missense substitution p.{orig_aa}{residue_num}{edited_aa}."


class CFTRBaseEditingModel:
    """Deterministic in silico model for CFTR base editing gRNA evaluation."""

    def __init__(self, editors: Optional[Dict[str, EditorProfile]] = None, cas_variants: Optional[Dict[str, CasVariantProfile]] = None):
        self.editors = editors or EDITORS
        self.cas_variants = cas_variants or CAS_VARIANTS
        self.alleles = CFTR_ALLELE_CATALOG

    def design_candidate_grna(
        self,
        allele_name: str,
        target_pos_in_protospacer: int,
        editor_name: str = "ABE8e",
        cas_variant: str = "SpCas9",
        strand: str = "coding",
    ) -> Optional[GuideRNAResult]:
        """Construct and evaluate a specific guide RNA with the target base positioned at target_pos_in_protospacer.

        Args:
          allele_name: Allele identifier ("W1282X", "R553X", "G542X", "2789+5G->A").
          target_pos_in_protospacer: 1-indexed position (1 to 20) where the mutant nucleotide is placed.
          editor_name: "ABE8e", "ABE7.10", or "AncBE4max".
          cas_variant: "SpCas9", "SpRY", or "SpG".
          strand: "coding" or "non_coding".
        """
        allele = self.alleles.get(allele_name)
        editor = self.editors.get(editor_name)
        cas = self.cas_variants.get(cas_variant)

        if not allele or not editor or not cas:
            return None

        # Determine target strand sequence context
        if strand == "coding":
            ctx = allele.coding_sequence_context
            mut_idx = allele.mutation_index_in_context
        else:
            ctx = reverse_complement(allele.coding_sequence_context)
            mut_idx = len(ctx) - 1 - allele.mutation_index_in_context

        # We want the mutant base to be at target_pos_in_protospacer (1-indexed, so index target_pos - 1)
        # Protospacer starts at: mut_idx - (target_pos_in_protospacer - 1)
        protospacer_start = mut_idx - (target_pos_in_protospacer - 1)
        protospacer_end = protospacer_start + 20
        pam_start = protospacer_end
        pam_end = pam_start + 3

        if protospacer_start < 0 or pam_end > len(ctx):
            return None

        protospacer_seq = ctx[protospacer_start:protospacer_end]
        pam_seq = ctx[pam_start:pam_end]

        # Calculate PAM affinity
        k_pam = score_pam_affinity(pam_seq, cas_variant)
        if k_pam <= 0.0:
            return None

        # Calculate deamination weight for target
        w_target = compute_editing_window(target_pos_in_protospacer, editor)
        lambda_chromatin = 1.0
        on_target_eff = w_target * k_pam * lambda_chromatin

        # Analyze bystanders within protospacer (any target base other than the intended target)
        bystanders: List[BystanderMutation] = []
        target_char = editor.target_base

        for p_idx in range(len(protospacer_seq)):
            pos_1indexed = p_idx + 1
            base = protospacer_seq[p_idx]

            if pos_1indexed == target_pos_in_protospacer:
                continue

            if base == target_char:
                # Bystander deamination occurs
                w_bystander = compute_editing_window(pos_1indexed, editor)
                eff_bystander = w_bystander * k_pam * lambda_chromatin

                # Determine codon change if inside coding region
                if allele.mutation_type == "nonsense" and allele.residue_number is not None:
                    # Compute relative codon offset from mutant codon
                    # Mutant base is at mut_idx in ctx
                    base_ctx_idx = protospacer_start + p_idx
                    diff_nt = base_ctx_idx - mut_idx

                    # In coding context, each residue is 3 nt
                    # Mutant base in W1282X: c.3846G>A (3rd base of codon 1282) -> codon_offset 0, base_in_codon 2
                    # Mutant base in R553X: c.1657C>T (1st base of codon 553) -> codon_offset 0, base_in_codon 0
                    if allele.allele_name == "W1282X":
                        codon_res = 1282 + (diff_nt // 3)
                        orig_codon_str = "TGA" if codon_res == 1282 else ("CAA" if codon_res == 1283 else "GAA")
                    elif allele.allele_name == "R553X":
                        codon_res = 553 + (diff_nt // 3)
                        orig_codon_str = "TGA" if codon_res == 553 else ("CAG" if codon_res == 554 else "AAG")
                    else:
                        codon_res = 542 + (diff_nt // 3)
                        orig_codon_str = "TGA"

                    orig_aa = DNA_CODON_TABLE.get(orig_codon_str, "X")
                    # Simulate A->G or C->T in codon
                    edited_codon_str = orig_codon_str.replace("A", "G") if editor.editor_type == "ABE" else orig_codon_str.replace("C", "T")
                    edited_aa = DNA_CODON_TABLE.get(edited_codon_str, orig_aa)

                    ddg, classification, rationale = evaluate_structural_perturbation(
                        allele, codon_res, orig_aa, edited_aa
                    )
                else:
                    codon_res = None
                    orig_codon_str = None
                    edited_codon_str = None
                    orig_aa = None
                    edited_aa = None
                    ddg = 0.0
                    classification = "benign"
                    rationale = "Non-coding intronic base editing; zero polypeptide perturbation."

                bystanders.append(
                    BystanderMutation(
                        position_in_protospacer=pos_1indexed,
                        original_base=base,
                        edited_base=editor.product_base,
                        deamination_weight=round(w_bystander, 4),
                        editing_efficiency=round(eff_bystander, 4),
                        codon_residue_index=codon_res,
                        original_codon=orig_codon_str,
                        edited_codon=edited_codon_str,
                        original_aa=orig_aa,
                        edited_aa=edited_aa,
                        delta_delta_g_fold_kcal_per_mol=round(ddg, 2),
                        classification=classification,
                        structural_rationale=rationale,
                    )
                )

        bystander_weights = [b.deamination_weight for b in bystanders]
        purity_ratio = calculate_bystander_purity_ratio(on_target_eff, bystander_weights)

        max_ddg = max([b.delta_delta_g_fold_kcal_per_mol for b in bystanders], default=0.0)

        # Overall Feasibility Composite Score (0 to 100)
        # 40% on-target efficiency, 40% bystander purity, 20% structural tolerance
        eff_score = min(1.0, on_target_eff) * 40.0
        purity_score = purity_ratio * 40.0
        struct_penalty = min(20.0, max_ddg * 5.0)
        struct_score = max(0.0, 20.0 - struct_penalty)
        composite = eff_score + purity_score + struct_score

        if allele.requires_transversion_editor:
            verdict = "UNFEASIBLE_TRANSVERSION (Requires Dual Deaminase or Prime Editing)"
            composite = min(composite, 25.0)
        elif max_ddg > 2.0:
            verdict = "HIGH_BYSTANDER_HAZARD (Destabilizing Missense Edit)"
        elif composite >= 75.0:
            verdict = "HIGH_FIDELITY_RESCUE (Optimal Curative Candidate)"
        else:
            verdict = "MODERATE_FIDELITY"

        grna_id = f"gRNA-{allele_name}-{strand[:4]}-p{target_pos_in_protospacer}-{cas_variant}"

        return GuideRNAResult(
            grna_id=grna_id,
            allele_name=allele_name,
            strand=strand,
            protospacer_sequence=protospacer_seq,
            pam_sequence=pam_seq,
            cas_variant=cas_variant,
            editor_name=editor_name,
            target_base_position=target_pos_in_protospacer,
            deamination_weight=round(w_target, 4),
            pam_affinity=round(k_pam, 3),
            chromatin_factor=round(lambda_chromatin, 2),
            on_target_efficiency=round(on_target_eff, 4),
            bystander_count=len(bystanders),
            bystanders=bystanders,
            bystander_purity_ratio=round(purity_ratio, 4),
            max_destabilizing_ddg=round(max_ddg, 2),
            overall_feasibility_score=round(composite, 2),
            curative_verdict=verdict,
            clinical_trial_reference="ClinicalTrials.gov NCT06237335",
            cryoem_context=allele.pdb_6msm_context,
        )

    def scan_allele(
        self,
        allele_name: str,
        editor_name: str = "ABE8e",
        cas_variant: str = "SpCas9",
        min_efficiency: float = 0.05,
    ) -> List[GuideRNAResult]:
        """Scan all protospacer positions (positions 1-20) across both strands for candidate gRNAs."""
        allele = self.alleles.get(allele_name)
        if not allele:
            return []

        results: List[GuideRNAResult] = []
        strands = ["coding", "non_coding"] if allele.target_strand_to_edit == "both" else [allele.target_strand_to_edit]

        for s in strands:
            for pos in range(1, 21):
                res = self.design_candidate_grna(
                    allele_name=allele_name,
                    target_pos_in_protospacer=pos,
                    editor_name=editor_name,
                    cas_variant=cas_variant,
                    strand=s,
                )
                if res and res.on_target_efficiency >= min_efficiency:
                    results.append(res)

        # Sort by composite feasibility score descending
        results.sort(key=lambda x: x.overall_feasibility_score, reverse=True)
        return results

    def get_top_candidate(
        self,
        allele_name: str,
        editor_name: str = "ABE8e",
        cas_variant: str = "SpCas9",
    ) -> Optional[GuideRNAResult]:
        """Return the highest-scoring candidate gRNA for a given allele."""
        cands = self.scan_allele(allele_name, editor_name, cas_variant)
        return cands[0] if cands else None

    def evaluate_all_priority_alleles(self, editor_name: str = "ABE8e") -> Dict[str, Any]:
        """Generate comprehensive screening matrix across all priority CFTR refractory mutations."""
        out: Dict[str, Any] = {}
        for allele_name in self.alleles:
            # Try SpCas9 first; if no candidate meets threshold, try SpRY
            top = self.get_top_candidate(allele_name, editor_name, "SpCas9")
            if not top or top.on_target_efficiency < 0.20:
                top = self.get_top_candidate(allele_name, editor_name, "SpRY")

            if top:
                out[allele_name] = asdict(top)
            else:
                out[allele_name] = {
                    "allele_name": allele_name,
                    "curative_verdict": "NO_ELIGIBLE_PAM_WINDOW",
                    "requires_transversion_editor": self.alleles[allele_name].requires_transversion_editor,
                }
        return out


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Base Editing (ABE8e / CBE) Bystander Deamination Window & PAM Accessibility Model for CFTR (T-0219)"
    )
    parser.add_argument(
        "--allele",
        choices=["W1282X", "R553X", "G542X", "2789+5G->A", "all"],
        default="all",
        help="Target CFTR refractory allele",
    )
    parser.add_argument(
        "--editor",
        choices=["ABE8e", "ABE7.10", "AncBE4max"],
        default="ABE8e",
        help="Base editor architecture",
    )
    parser.add_argument(
        "--cas",
        choices=["SpCas9", "SpRY", "SpG", "all"],
        default="all",
        help="CRISPR-Cas PAM variant",
    )
    parser.add_argument("--json", action="store_true", help="Print structured JSON output")
    parser.add_argument("--export", type=str, help="Export JSON results to path")
    args = parser.parse_args()

    model = CFTRBaseEditingModel()

    target_alleles = [args.allele] if args.allele != "all" else list(CFTR_ALLELE_CATALOG.keys())
    cas_list = [args.cas] if args.cas != "all" else ["SpCas9", "SpRY", "SpG"]

    results: Dict[str, Any] = {}

    for a in target_alleles:
        results[a] = {}
        for c in cas_list:
            top = model.get_top_candidate(a, args.editor, c)
            if top:
                results[a][c] = asdict(top)

    if args.export:
        out_p = Path(args.export)
        out_p.parent.mkdir(parents=True, exist_ok=True)
        with open(out_p, "w", encoding="utf-8") as f:
            json.dump(results, f, indent=2)
        print(f"Exported base editing candidates to {out_p}")

    if args.json:
        print(json.dumps(results, indent=2))
    else:
        print("========================================================================")
        print("CFTR BASE EDITING (ABE8e / CBE) BIOPHYSICAL FEASIBILITY MODEL (T-0219)")
        print("========================================================================")
        print(f"Editor: {args.editor} | Active Window: pos 3-9 | PDB Context: 6MSM")
        print("------------------------------------------------------------------------")
        for a_name, cas_data in results.items():
            print(f"\nAllele: {a_name} ({CFTR_ALLELE_CATALOG[a_name].cdna_name}, {CFTR_ALLELE_CATALOG[a_name].exon_intron})")
            print(f"  PDB Context: {CFTR_ALLELE_CATALOG[a_name].pdb_6msm_context}")
            if not cas_data:
                print("  No viable gRNAs found within PAM constraints.")
                continue
            for c_name, g_info in cas_data.items():
                print(f"  [{c_name}] Guide: {g_info['grna_id']}")
                print(f"    Protospacer (20nt): 5'-{g_info['protospacer_sequence']}-3' | PAM: {g_info['pam_sequence']}")
                print(f"    Target Pos: pos {g_info['target_base_position']} (W_E={g_info['deamination_weight']}, Efficiency={g_info['on_target_efficiency']})")
                print(f"    Bystanders: {g_info['bystander_count']} | Purity Ratio (P): {g_info['bystander_purity_ratio']} | Max DDG: {g_info['max_destabilizing_ddg']} kcal/mol")
                print(f"    Composite Feasibility: {g_info['overall_feasibility_score']}/100 | Verdict: {g_info['curative_verdict']}")
        print("\n========================================================================")
        print("NLM Attribution & Research Boundary:")
        print(NLM_ATTRIBUTION)
        print(RESEARCH_BOUNDARY)
        print("========================================================================")


if __name__ == "__main__":
    main()
