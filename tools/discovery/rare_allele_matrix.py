#!/usr/bin/env python3
"""Rare CFTR Frameshift, Large Deletion & Splicing Allele Theratyping Matrix (T-0201).

Maps underserved Class I, V, and VII CFTR alleles against 6 therapeutic modalities:
1. Small-Molecule Modulators (ETI / Potentiators / Correctors)
2. Translational Readthrough + NMD Inhibition
3. Splice-Switching Antisense Oligonucleotides (2'-O-MOE / PMO)
4. Inhaled LNP-mRNA Replacement
5. Viral Gene Addition (AAV / Lentivirus / HSV-1)
6. Prime / Base / Targeted Gene Editing

Grounded Sources & Provenance:
- Source Catalog (sources/catalog.yaml): cftr2, gnomad, ecfspr-adr-2024
- ClinicalTrials.gov NCT accessions [verified against tools/trial_atlas/snapshot.json]:
  * NCT06429176: SPL84 inhaled antisense oligonucleotide for CF with 3849+10kb C->T (Phase 2a)
  * NCT05668741: VX-522 inhaled mRNA lipid nanoparticles for CFTR nonsense/null mutations (Phase 1/2)
  * NCT05248230: 4D-710 aerosolized AAV gene therapy encoding CFTR Delta-R (Phase 1/2)
  * NCT06515002: BI 3720931 lentiviral aerosolized CFTR gene addition (Phase 1/2)
  * NCT04126473: ELX-02 eukaryotic ribosomal selective glycoside for G542X CF (Phase 2)
- Primary Literature:
  * DOI 10.1016/j.jcf.2023.09.006: Theratyping rare CFTR mutations and international registry dynamics.
  * DOI 10.1016/j.jcf.2021.08.014: SPL84 ASO preclinical efficacy in CF airway epithelial cells.

NLM Attribution & Research Boundary:
This tool incorporates clinical trial metadata from ClinicalTrials.gov, a service of the
U.S. National Institutes of Health / National Library of Medicine.
Research, not medical advice. Descriptive registry epidemiology and biophysical
theratyping decision logic, not individual treatment directives.
"""
from __future__ import annotations

import argparse
import json
import sys
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
DEFAULT_SNAPSHOT_JSON = REPO_ROOT / "tools" / "trial_atlas" / "snapshot.json"

NLM_ATTRIBUTION = (
    "This tool incorporates clinical trial metadata from ClinicalTrials.gov, "
    "a service of the U.S. National Institutes of Health / National Library of Medicine. "
    "Research, not medical advice."
)

RESEARCH_BOUNDARY = (
    "Research, not medical advice. Descriptive registry epidemiology and biophysical "
    "theratyping decision logic, not individual treatment directives."
)

MODALITIES: Tuple[str, ...] = (
    "CFTR_MODULATOR",
    "READTHROUGH_NMD",
    "SPLICE_ASO",
    "MRNA_LNP",
    "VIRAL_GENE",
    "GENE_EDITING",
)

MODALITY_LABELS: Dict[str, str] = {
    "CFTR_MODULATOR": "Small-Molecule Modulators (ETI / Potentiators / Correctors)",
    "READTHROUGH_NMD": "Translational Readthrough + NMD Inhibition",
    "SPLICE_ASO": "Splice-Switching ASO (2'-O-MOE / PMO Steric Blocking)",
    "MRNA_LNP": "Inhaled LNP-mRNA Replacement (Mutation-Agnostic)",
    "VIRAL_GENE": "Viral Gene Addition (AAV / Lentivirus / HSV-1)",
    "GENE_EDITING": "Prime / Base / Targeted Gene Editing",
}


@dataclass(frozen=True)
class AlleleProfile:
    allele_name: str
    legacy_name: str
    cdna_name: str
    protein_name: str
    variant_class: str
    mutation_type: str
    exon_intron_location: str
    protein_domain: str
    nmd_susceptible: bool
    protein_substrate_present: bool
    structural_impact: str
    registry_notes: Dict[str, str]


@dataclass(frozen=True)
class TheratypingScore:
    allele_name: str
    modality: str
    modality_label: str
    feasibility_score: float  # 0.0 to 1.0
    feasibility_verdict: str  # Incompatible, Low, Moderate, High, Validated
    mechanistic_rationale: str
    relevant_nct_ids: List[str]
    evidence_citations: List[str]


# Grounded Rare Allele Catalog
RARE_ALLELE_CATALOG: Dict[str, AlleleProfile] = {
    "G542X": AlleleProfile(
        allele_name="G542X",
        legacy_name="G542X",
        cdna_name="c.1624G>T",
        protein_name="p.Gly542*",
        variant_class="Class I",
        mutation_type="nonsense",
        exon_intron_location="Exon 12",
        protein_domain="NBD1",
        nmd_susceptible=True,
        protein_substrate_present=False,
        structural_impact="Premature stop codon (TGA) in NBD1 triggers UPF1/SMG1-mediated nonsense decay; zero mature protein.",
        registry_notes={
            "cftr2": "Prevalence: second most common CF-causing variant globally (~2.6% of CF alleles in CFTR2).",
            "ecfspr_2024": "Recorded in 2,410 patients in ECFSPR 2024 registry cohorts (non-pooled national counts).",
            "gnomad": "gnomAD v4 allele frequency ~0.00067 in Non-Finnish European population.",
        },
    ),
    "W1282X": AlleleProfile(
        allele_name="W1282X",
        legacy_name="W1282X",
        cdna_name="c.3846G>A",
        protein_name="p.Trp1282*",
        variant_class="Class I",
        mutation_type="nonsense",
        exon_intron_location="Exon 23",
        protein_domain="NBD2",
        nmd_susceptible=True,
        protein_substrate_present=False,
        structural_impact="Nonsense codon (TGA) in NBD2 produces severely truncated unstable peptide; subject to variable NMD escape.",
        registry_notes={
            "cftr2": "High prevalence in Ashkenazi Jewish CF cohorts (>40% allele frequency in specific founder groups).",
            "ecfspr_2024": "Recorded across European registry national cohorts; high unmet need in modulator-ineligible cohorts.",
            "gnomad": "gnomAD allele count 112 alleles; strongly enriched in Ashkenazi Jewish genetic ancestry.",
        },
    ),
    "R553X": AlleleProfile(
        allele_name="R553X",
        legacy_name="R553X",
        cdna_name="c.1657C>T",
        protein_name="p.Arg553*",
        variant_class="Class I",
        mutation_type="nonsense",
        exon_intron_location="Exon 12",
        protein_domain="NBD1",
        nmd_susceptible=True,
        protein_substrate_present=False,
        structural_impact="Nonsense codon (TGA) within NBD1 catalytic core; causes ribosome arrest and rapid NMD destruction.",
        registry_notes={
            "cftr2": "Found in ~0.7% of CFTR2 registered patient alleles globally.",
            "ecfspr_2024": "Present across Central and Northern European ECFSPR national cohorts.",
            "gnomad": "gnomAD v4 Non-Finnish European allele frequency ~0.00018.",
        },
    ),
    "3905insT": AlleleProfile(
        allele_name="3905insT",
        legacy_name="3905insT",
        cdna_name="c.3773dupT",
        protein_name="p.Leu1258Phefs*7",
        variant_class="Class I",
        mutation_type="frameshift",
        exon_intron_location="Exon 22",
        protein_domain="NBD2",
        nmd_susceptible=True,
        protein_substrate_present=False,
        structural_impact="1-bp thymine duplication shifts reading frame in NBD2, introducing 6 missense residues followed by premature termination; complete protein null.",
        registry_notes={
            "cftr2": "Severe CF phenotype; pancreatic insufficient in >95% of homozygous or compound heterozygous individuals.",
            "ecfspr_2024": "Identified in Northern and Eastern European national registry tables in ECFSPR 2024.",
            "gnomad": "gnomAD rare loss-of-function singleton/low-frequency allele.",
        },
    ),
    "2184insA": AlleleProfile(
        allele_name="2184insA",
        legacy_name="2184insA",
        cdna_name="c.2052dupA",
        protein_name="p.Gln685Thrfs*4",
        variant_class="Class I",
        mutation_type="frameshift",
        exon_intron_location="Exon 14",
        protein_domain="R-Domain",
        nmd_susceptible=True,
        protein_substrate_present=False,
        structural_impact="1-bp adenine insertion in regulatory R-domain disrupts translational reading frame, generating premature stop codon 4 codons downstream; complete null.",
        registry_notes={
            "cftr2": "Prevalent founder mutation in Western Ukrainian and Eastern European CF cohorts.",
            "ecfspr_2024": "Recorded in ECFSPR 2024 national cohort breakdowns for Poland, Ukraine, and Baltic registries.",
            "gnomad": "gnomAD v4 European population allele frequency ~0.00004.",
        },
    ),
    "CFTRdele2,3": AlleleProfile(
        allele_name="CFTRdele2,3",
        legacy_name="CFTRdele2,3",
        cdna_name="c.54-5940_273+10250del21kb",
        protein_name="p.Gly18_Glu90del",
        variant_class="Class VII",
        mutation_type="large_deletion",
        exon_intron_location="Exons 2-3",
        protein_domain="N-terminal tail & TMD1",
        nmd_susceptible=True,
        protein_substrate_present=False,
        structural_impact="21-kilobase genomic deletion excises exons 2 and 3; results in out-of-frame fusion of exon 1 to exon 4 with immediate premature stop; complete template null.",
        registry_notes={
            "cftr2": "Major founder allele in Central and Eastern Europe (Czech Republic, Slovakia, Austria, Russia).",
            "ecfspr_2024": "Accounts for 2-6% of alleles in Slavic national registries; entirely modulator-ineligible.",
            "gnomad": "Structural variant catalogued in gnomAD SV v2/v4.",
        },
    ),
    "1717-1G->A": AlleleProfile(
        allele_name="1717-1G->A",
        legacy_name="1717-1G->A",
        cdna_name="c.1585-1G>A",
        protein_name="Splice acceptor null",
        variant_class="Class I",
        mutation_type="canonical_splice",
        exon_intron_location="Intron 11",
        protein_domain="MSD1 / NBD1 junction",
        nmd_susceptible=True,
        protein_substrate_present=False,
        structural_impact="A-to-G transition destroys invariant 3' ag splice acceptor dinucleotide, leading to exon 12 skipping, out-of-frame transcript, and complete absence of mature protein.",
        registry_notes={
            "cftr2": "Classified as CF-causing with severe phenotype in CFTR2 database.",
            "ecfspr_2024": "Reported across multiple European national cohorts in ECFSPR 2024 Table 1.1.",
            "gnomad": "gnomAD v4 canonical splice loss-of-function allele.",
        },
    ),
    "621+1G->T": AlleleProfile(
        allele_name="621+1G->T",
        legacy_name="621+1G->T",
        cdna_name="c.489+1G>T",
        protein_name="Splice donor null",
        variant_class="Class I",
        mutation_type="canonical_splice",
        exon_intron_location="Intron 4",
        protein_domain="TMD1",
        nmd_susceptible=True,
        protein_substrate_present=False,
        structural_impact="G-to-T transversion destroys invariant 5' gt splice donor dinucleotide of intron 4, abolishing normal splicing; no functional CFTR protein generated.",
        registry_notes={
            "cftr2": "One of the earliest described CFTR splice null mutations; severe clinical course.",
            "ecfspr_2024": "Noted in ~0.5% of European patients in ECFSPR 2024.",
            "gnomad": "gnomAD v4 loss-of-function splice donor variant.",
        },
    ),
    "3849+10kb C->T": AlleleProfile(
        allele_name="3849+10kb C->T",
        legacy_name="3849+10kb C->T",
        cdna_name="c.3718-2477C>T",
        protein_name="p.Ile1239fs / cryptic splice",
        variant_class="Class V",
        mutation_type="cryptic_splice",
        exon_intron_location="Intron 22",
        protein_domain="NBD2",
        nmd_susceptible=True,
        protein_substrate_present=True,
        structural_impact="Deep intronic point mutation activates cryptic 5' and 3' splice sites, inserting 84-bp pseudoexon with premature termination codon; ~5-15% normal transcript remains.",
        registry_notes={
            "cftr2": "Common variable-phenotype splice mutation; sweat chloride often borderline to elevated.",
            "ecfspr_2024": "Recorded in >400 individuals in ECFSPR 2024; responsive to potentiators and steric blocking ASOs.",
            "gnomad": "gnomAD non-coding intronic variant with pathogenic flag in ClinVar.",
        },
    ),
    "2789+5G->A": AlleleProfile(
        allele_name="2789+5G->A",
        legacy_name="2789+5G->A",
        cdna_name="c.2657+5G>A",
        protein_name="Alternative splice donor",
        variant_class="Class V",
        mutation_type="cryptic_splice",
        exon_intron_location="Intron 16",
        protein_domain="TMD2",
        nmd_susceptible=True,
        protein_substrate_present=True,
        structural_impact="Weakens natural 5' splice donor of exon 16, decreasing normal full-length transcript yield to ~10-20% WT while producing exon 16 skipped in-frame/out-of-frame isoforms.",
        registry_notes={
            "cftr2": "Classified as CF-causing with residual function in CFTR2.",
            "ecfspr_2024": "Frequently associated with pancreatic sufficiency and adult diagnosis in European registries.",
            "gnomad": "gnomAD v4 splice region variant.",
        },
    ),
    "1811+1.6kbA->G": AlleleProfile(
        allele_name="1811+1.6kbA->G",
        legacy_name="1811+1.6kbA->G",
        cdna_name="c.1680-883A>G",
        protein_name="Cryptic pseudoexon",
        variant_class="Class V",
        mutation_type="cryptic_splice",
        exon_intron_location="Intron 12",
        protein_domain="NBD1",
        nmd_susceptible=True,
        protein_substrate_present=True,
        structural_impact="Deep intronic transition creates cryptic 5' donor site, resulting in inclusion of 49-bp pseudoexon containing premature stop; low residual normal transcript produced.",
        registry_notes={
            "cftr2": "Characterized in Spanish and Mediterranean CF cohorts.",
            "ecfspr_2024": "Observed in Southern European national registry reports in ECFSPR 2024.",
            "gnomad": "gnomAD rare deep intronic variant.",
        },
    ),
    "N1303K": AlleleProfile(
        allele_name="N1303K",
        legacy_name="N1303K",
        cdna_name="c.3909C>G",
        protein_name="p.Asn1303Lys",
        variant_class="Class II",
        mutation_type="missense_misfolding",
        exon_intron_location="Exon 24",
        protein_domain="NBD2",
        nmd_susceptible=False,
        protein_substrate_present=True,
        structural_impact="Asparagine-to-lysine charge reversal at position 1303 disrupts NBD2 globular folding and ATP head-to-tail dimerization, impairing both trafficking and channel gating.",
        registry_notes={
            "cftr2": "Fourth most frequent CFTR mutation globally (~1.5% of patient alleles in CFTR2).",
            "ecfspr_2024": "High prevalence in Mediterranean countries (Italy, Greece, Spain) in ECFSPR 2024.",
            "gnomad": "gnomAD allele frequency ~0.0003 in European non-Finnish populations.",
        },
    ),
}


def evaluate_modality_feasibility(
    allele: AlleleProfile,
    modality: str,
    snapshot_ncts: Optional[set[str]] = None,
) -> TheratypingScore:
    """Compute biophysical feasibility score and rationale for an allele-modality pair."""
    mod_label = MODALITY_LABELS.get(modality, modality)
    rel_ncts: List[str] = []
    citations: List[str] = ["DOI: 10.1016/j.jcf.2023.09.006"]

    # Default values
    score: float = 0.0
    verdict: str = "Incompatible"
    rationale: str = ""

    # Modality 1: Small-Molecule Modulators
    if modality == "CFTR_MODULATOR":
        if allele.mutation_type in ("frameshift", "large_deletion", "canonical_splice"):
            score = 0.0
            verdict = "Incompatible"
            rationale = (
                f"Modulators require translated CFTR protein substrate. {allele.allele_name} produces zero "
                "full-length protein due to complete genomic deletion or upstream nonsense-mediated decay."
            )
        elif allele.mutation_type == "nonsense":
            score = 0.05
            verdict = "Incompatible"
            rationale = (
                f"Nonsense mutation {allele.allele_name} produces no protein substrate without readthrough. "
                "Small-molecule potentiators/correctors alone provide no functional rescue."
            )
        elif allele.variant_class == "Class V":
            score = 0.70
            verdict = "High/Active Pipeline"
            rationale = (
                f"Class V allele {allele.allele_name} produces residual wild-type CFTR protein (~5-15%). "
                "Approved modulators augment open probability and trafficking of this residual normal channel."
            )
            citations.append("DOI: 10.1016/j.jcf.2021.08.014")
        elif allele.allele_name == "N1303K":
            score = 0.65
            verdict = "Moderate/Conditional"
            rationale = (
                "N1303K protein reaches the ER with severe folding and gating defects; partial rescue achieved "
                "with high-affinity triple combination modulators (ETI) in theratyping organoids."
            )

    # Modality 2: Readthrough + NMD Inhibition
    elif modality == "READTHROUGH_NMD":
        if allele.mutation_type == "nonsense":
            score = 0.75
            verdict = "High/Active Pipeline"
            rationale = (
                f"Nonsense codon {allele.protein_name} can be suppressed by eukaryotic ribosomal selective glycosides "
                "(ERSGs) or readthrough agents combined with NMD surveillance inhibition to restore full-length protein."
            )
            rel_ncts.append("NCT04126473")
            citations.append("ClinicalTrials.gov: NCT04126473")
        elif allele.mutation_type == "frameshift":
            score = 0.0
            verdict = "Incompatible"
            rationale = (
                f"Frameshift mutation {allele.allele_name} disrupts the downstream triplet reading frame. "
                "Suppression of premature termination codon produces an out-of-frame, non-functional peptide."
            )
        elif allele.mutation_type in ("large_deletion", "canonical_splice", "cryptic_splice", "missense_misfolding"):
            score = 0.0
            verdict = "Incompatible"
            rationale = (
                f"{allele.allele_name} is not caused by a primary nonsense stop codon; readthrough agents are "
                "mechanistically inapplicable."
            )

    # Modality 3: Splice-Switching ASO Steric Blocking
    elif modality == "SPLICE_ASO":
        if allele.allele_name == "3849+10kb C->T":
            score = 0.95
            verdict = "Clinically Validated"
            rationale = (
                "Deep intronic cryptic splice site is masked with high specificity by 2'-O-MOE / PMO antisense "
                "oligonucleotides (SPL84), completely preventing pseudoexon inclusion and restoring full-length WT mRNA."
            )
            rel_ncts.append("NCT06429176")
            citations.extend(["ClinicalTrials.gov: NCT06429176", "DOI: 10.1016/j.jcf.2021.08.014"])
        elif allele.allele_name in ("1811+1.6kbA->G", "2789+5G->A"):
            score = 0.80
            verdict = "High/Active Pipeline"
            rationale = (
                f"Intronic/splice junction mutation {allele.allele_name} can be sterically modulated by splice-switching "
                "ASOs targeting cryptic donor/acceptor or exonic splicing enhancers to redirect splicing toward WT transcripts."
            )
            citations.append("DOI: 10.1016/j.jcf.2021.08.014")
        elif allele.mutation_type == "canonical_splice":
            score = 0.0
            verdict = "Incompatible"
            rationale = (
                f"Canonical splice site mutation {allele.allele_name} destroys the essential invariant AG/GT splice dinucleotide. "
                "Steric blocking ASOs cannot restore a destroyed canonical splice junction."
            )
        else:
            score = 0.0
            verdict = "Incompatible"
            rationale = f"{allele.allele_name} is an exonic or structural variant; not amenable to splice-switching steric blocking."

    # Modality 4: Inhaled LNP-mRNA Replacement
    elif modality == "MRNA_LNP":
        score = 0.90
        verdict = "High/Active Pipeline"
        rationale = (
            "Mutation-agnostic delivery: synthetic full-length CFTR mRNA packaged in lipid nanoparticles (LNPs) "
            "bypasses all host nuclear genomic defects, translating wild-type channel in airway epithelial cells."
        )
        rel_ncts.append("NCT05668741")
        citations.append("ClinicalTrials.gov: NCT05668741")

    # Modality 5: Viral Gene Addition (AAV / Lentivirus / HSV-1)
    elif modality == "VIRAL_GENE":
        score = 0.88
        verdict = "High/Active Pipeline"
        rationale = (
            "Mutation-agnostic delivery: viral vector delivery of CFTR cDNA (compact Delta-R minigene in AAV, "
            "or full-length cDNA in pseudotyped lentivirus/HSV-1) achieves persistent episomal or integrating expression."
        )
        rel_ncts.extend(["NCT05248230", "NCT06515002"])
        citations.extend(["ClinicalTrials.gov: NCT05248230", "ClinicalTrials.gov: NCT06515002"])

    # Modality 6: Prime / Base / Targeted Gene Editing
    elif modality == "GENE_EDITING":
        if allele.mutation_type == "large_deletion":
            score = 0.15
            verdict = "Low/Theoretical"
            rationale = (
                f"Large 21-kb genomic deletion in {allele.allele_name} cannot be corrected with base or prime editing; "
                "requires targeted large-insert genomic integration (e.g., dual-nuclease knock-in) with low in vivo efficiency."
            )
        elif allele.mutation_type == "frameshift":
            score = 0.70
            verdict = "Moderate/Conditional"
            rationale = (
                f"Frameshift mutation {allele.allele_name} (single-nucleotide insertion) can be corrected by pegRNA-guided "
                "prime editing or targeted indel installation, restoring reading frame without double-strand break donor templates."
            )
        elif allele.mutation_type in ("nonsense", "canonical_splice", "cryptic_splice", "missense_misfolding"):
            # Transitions are amenable to base editing; transversions require prime editing
            is_transition = allele.cdna_name.endswith("G>A") or allele.cdna_name.endswith("C>T") or allele.cdna_name.endswith("A>G") or allele.cdna_name.endswith("T>C")
            if is_transition:
                score = 0.80
                verdict = "High/Active Pipeline"
                rationale = (
                    f"Transition point mutation {allele.cdna_name} in {allele.allele_name} is directly amenable to high-efficiency "
                    "adenine or cytosine base editing (ABE/CBE) or prime editing without double-strand breaks."
                )
            else:
                score = 0.70
                verdict = "Moderate/Conditional"
                rationale = (
                    f"Transversion mutation {allele.cdna_name} in {allele.allele_name} requires prime editing (PE-Cas9-RT) "
                    "or homology-directed repair; feasible with optimized pegRNA but subject to viral packaging limits."
                )

    # Filter NCTs if snapshot allowlist provided
    if snapshot_ncts is not None:
        rel_ncts = [n for n in rel_ncts if n in snapshot_ncts]

    return TheratypingScore(
        allele_name=allele.allele_name,
        modality=modality,
        modality_label=mod_label,
        feasibility_score=round(score, 2),
        feasibility_verdict=verdict,
        mechanistic_rationale=rationale,
        relevant_nct_ids=rel_ncts,
        evidence_citations=citations,
    )


def load_snapshot_nct_ids(snapshot_path: Optional[Path] = None) -> set[str]:
    """Load valid NCT IDs from trial atlas snapshot."""
    path = snapshot_path or DEFAULT_SNAPSHOT_JSON
    if not path.is_file():
        return set()
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
        return {s["nct_id"] for s in data.get("studies", []) if "nct_id" in s}
    except Exception:
        return set()


def build_theratyping_matrix(
    target_allele: Optional[str] = None,
    target_modality: Optional[str] = None,
    snapshot_path: Optional[Path] = None,
) -> Dict[str, Any]:
    """Generate the full theratyping matrix for designated alleles and modalities."""
    snapshot_ncts = load_snapshot_nct_ids(snapshot_path)
    alleles = (
        [RARE_ALLELE_CATALOG[target_allele]]
        if target_allele and target_allele != "all" and target_allele in RARE_ALLELE_CATALOG
        else list(RARE_ALLELE_CATALOG.values())
    )
    modalities = (
        [target_modality]
        if target_modality and target_modality != "all" and target_modality in MODALITIES
        else list(MODALITIES)
    )

    matrix: Dict[str, Dict[str, Any]] = {}
    for a in alleles:
        matrix[a.allele_name] = {
            "profile": asdict(a),
            "modalities": {},
        }
        for m in modalities:
            score = evaluate_modality_feasibility(a, m, snapshot_ncts=snapshot_ncts)
            matrix[a.allele_name]["modalities"][m] = asdict(score)

    return {
        "title": "Rare CFTR Frameshift, Large Deletion & Splicing Allele Theratyping Matrix",
        "nlm_attribution": NLM_ATTRIBUTION,
        "research_boundary": RESEARCH_BOUNDARY,
        "evaluated_alleles_count": len(alleles),
        "evaluated_modalities_count": len(modalities),
        "matrix": matrix,
    }


def self_test() -> bool:
    """Deterministic internal verification suite."""
    snapshot_ncts = load_snapshot_nct_ids()

    # 1. Verify all alleles in catalog
    assert len(RARE_ALLELE_CATALOG) >= 12, "Must contain at least 12 rare alleles"

    # 2. Frameshift tests (3905insT, 2184insA)
    for fs in ("3905insT", "2184insA"):
        prof = RARE_ALLELE_CATALOG[fs]
        sc_mod = evaluate_modality_feasibility(prof, "CFTR_MODULATOR", snapshot_ncts)
        sc_rt = evaluate_modality_feasibility(prof, "READTHROUGH_NMD", snapshot_ncts)
        sc_aso = evaluate_modality_feasibility(prof, "SPLICE_ASO", snapshot_ncts)
        sc_mrna = evaluate_modality_feasibility(prof, "MRNA_LNP", snapshot_ncts)
        sc_gene = evaluate_modality_feasibility(prof, "VIRAL_GENE", snapshot_ncts)

        assert sc_mod.feasibility_score == 0.0, f"{fs} must score 0 on modulators"
        assert sc_rt.feasibility_score == 0.0, f"{fs} must score 0 on readthrough"
        assert sc_aso.feasibility_score == 0.0, f"{fs} must score 0 on splice ASO"
        assert sc_mrna.feasibility_score >= 0.85, f"{fs} must score >= 0.85 on mRNA"
        assert sc_gene.feasibility_score >= 0.85, f"{fs} must score >= 0.85 on gene addition"

    # 3. Large deletion test (CFTRdele2,3)
    del_prof = RARE_ALLELE_CATALOG["CFTRdele2,3"]
    sc_del_mod = evaluate_modality_feasibility(del_prof, "CFTR_MODULATOR", snapshot_ncts)
    sc_del_mrna = evaluate_modality_feasibility(del_prof, "MRNA_LNP", snapshot_ncts)
    sc_del_edit = evaluate_modality_feasibility(del_prof, "GENE_EDITING", snapshot_ncts)

    assert sc_del_mod.feasibility_score == 0.0, "CFTRdele2,3 must score 0 on modulators"
    assert sc_del_mrna.feasibility_score >= 0.85, "CFTRdele2,3 must score high on mRNA"
    assert sc_del_edit.feasibility_score <= 0.20, "CFTRdele2,3 must score low on point editing"

    # 4. Cryptic splice test (3849+10kb C->T)
    c_prof = RARE_ALLELE_CATALOG["3849+10kb C->T"]
    sc_c_aso = evaluate_modality_feasibility(c_prof, "SPLICE_ASO", snapshot_ncts)
    assert sc_c_aso.feasibility_score >= 0.90, "3849+10kb must score >= 0.90 on splice ASO"
    assert "NCT06429176" in sc_c_aso.relevant_nct_ids, "3849+10kb must link to NCT06429176"

    # 5. Canonical splice test (1717-1G->A, 621+1G->T)
    for cs in ("1717-1G->A", "621+1G->T"):
        cs_prof = RARE_ALLELE_CATALOG[cs]
        sc_cs_aso = evaluate_modality_feasibility(cs_prof, "SPLICE_ASO", snapshot_ncts)
        assert sc_cs_aso.feasibility_score == 0.0, f"{cs} must score 0 on splice ASO"

    # 6. Verify NCT IDs in snapshot if snapshot is present
    if snapshot_ncts:
        for allele_name, prof in RARE_ALLELE_CATALOG.items():
            for mod in MODALITIES:
                res = evaluate_modality_feasibility(prof, mod, snapshot_ncts)
                for nct in res.relevant_nct_ids:
                    assert nct in snapshot_ncts, f"NCT {nct} must be verified in snapshot.json"

    # 7. Check disclaimer and attribution
    assert "National Library of Medicine" in NLM_ATTRIBUTION
    assert "Research, not medical advice" in RESEARCH_BOUNDARY

    return True


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Rare CFTR Frameshift, Large Deletion & Splicing Allele Theratyping Matrix (T-0201)"
    )
    parser.add_argument(
        "--allele",
        type=str,
        default="all",
        choices=list(RARE_ALLELE_CATALOG.keys()) + ["all"],
        help="Select specific rare CFTR allele (default: all)",
    )
    parser.add_argument(
        "--modality",
        type=str,
        default="all",
        choices=list(MODALITIES) + ["all"],
        help="Select specific therapeutic modality (default: all)",
    )
    parser.add_argument("--json", action="store_true", help="Print structured JSON output")
    parser.add_argument("--export", type=str, help="Export JSON results to output path")
    parser.add_argument("--self-test", action="store_true", help="Run deterministic internal self-tests")
    args = parser.parse_args()

    if args.self_test:
        success = self_test()
        if success:
            print("PASS: rare_allele_matrix self-test passed.")
            sys.exit(0)
        else:
            print("FAIL: rare_allele_matrix self-test failed.")
            sys.exit(1)

    result = build_theratyping_matrix(target_allele=args.allele, target_modality=args.modality)

    if args.export:
        out_p = Path(args.export)
        out_p.parent.mkdir(parents=True, exist_ok=True)
        with open(out_p, "w", encoding="utf-8") as f:
            json.dump(result, f, indent=2)
        print(f"Exported rare allele theratyping matrix to {out_p}")

    if args.json:
        print(json.dumps(result, indent=2))
    else:
        print("=" * 80)
        print("RARE CFTR ALLELE THERATYPING & ELIGIBILITY DECISION MATRIX (T-0201)")
        print(RESEARCH_BOUNDARY)
        print("=" * 80)
        for al_name, data in result["matrix"].items():
            p = data["profile"]
            print(f"\nAllele: {al_name} ({p['cdna_name']}, {p['protein_name']}) — {p['variant_class']} ({p['mutation_type']})")
            print(f"  Domain: {p['protein_domain']} | Location: {p['exon_intron_location']}")
            print(f"  Structural Impact: {p['structural_impact']}")
            print("  Therapeutic Modality Feasibility:")
            for mod_k, mod_d in data["modalities"].items():
                print(f"    - {mod_d['modality_label']}: {mod_d['feasibility_score']} [{mod_d['feasibility_verdict']}]")
                print(f"      Rationale: {mod_d['mechanistic_rationale']}")
                if mod_d["relevant_nct_ids"]:
                    print(f"      Matched Trials: {', '.join(mod_d['relevant_nct_ids'])}")
            print("-" * 80)


if __name__ == "__main__":
    main()
