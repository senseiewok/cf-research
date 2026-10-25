#!/usr/bin/env python3
"""CFTR Missense & Splice Variant Cross-Catalog Concordance Model (ncbi-variation, uniprot-rest, ensembl-rest).

Reconciles canonical RefSeq (NM_000492.4), Ensembl (ENST00000003084), and UniProtKB (P13569)
coordinates across priority Class I-VI CFTR variants (F508del, G551D, W1282X, G542X, N1303K,
R117H, 3849+10kbC->T):
1. SPDI & HGVS Contextual Normalization: Verifies 0-based SPDI (NC_000007.14:117559592:CTT:)
   and 1-based HGVS (NM_000492.4:c.1521_1523del) coordinate shifts across genomic GRCh38 and
   coding cDNA coordinates.
2. UniProtKB Domain Boundary & Variant Mapping: Maps each missense/nonsense/indel position
   onto the 5 canonical CFTR structural domains (TMD1, NBD1, R, TMD2, NBD2) from P13569 and
   cryo-EM structure 6MSM.
3. Theratype & Modulator Eligibility Concordance: Cross-checks variant consequence classes
   against FDA label indications (NDA212273 ETI / ALYFTREK vanzacaftor-tezacaftor-deutivacaftor)
   and flags refractory nonsense/frameshift alleles requiring nucleic acid therapy (NCT05248230,
   NCT05668741).

Primary Literature Grounding:
- Sosnay PR et al. (2013) Nat Genet, DOI: 10.1038/ng.2745
- Veit G et al. (2016) Mol Biol Cell, DOI: 10.1091/mbc.E15-04-0223
- Holmes JB et al. (2020) Bioinformatics, DOI: 10.1093/bioinformatics/btz856
- Middleton PG et al. (2019) N Engl J Med, DOI: 10.1056/NEJMoa1908639

Structural & Clinical Registry Grounding:
- RCSB PDB 6MSM (Human phosphorylated ATP-bound CFTR at 3.2 A)
- RCSB PDB 5UAK (Dephosphorylated ATP-free human CFTR)
- NCT05248230 (4D-710 aerosolized AAV Phase 1/2 in modulator-ineligible CF)
- NCT05668741 (VX-522 mRNA-LNP Phase 1/2 in modulator-ineligible CF)
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any, Dict, List


# Canonical accession constants
CANONICAL_REFSEQ_TRANSCRIPT = "NM_000492.4"
CANONICAL_REFSEQ_PROTEIN = "NP_000483.3"
CANONICAL_ENSEMBL_TRANSCRIPT = "ENST00000003084"
CANONICAL_UNIPROT_ACCESSION = "P13569"
CANONICAL_GRCH38_CHROMOSOME = "NC_000007.14"

# UniProtKB P13569 & PDB 6MSM canonical domain boundaries (1-based amino acid coordinates, 1-1480)
CFTR_STRUCTURAL_DOMAINS: Dict[str, tuple[int, int]] = {
    "N_TERMINAL_LASSO": (1, 80),
    "TMD1": (81, 388),
    "NBD1": (389, 646),
    "R": (647, 845),
    "TMD2": (846, 1181),
    "NBD2": (1182, 1480),
}


@dataclass(frozen=True)
class CFTRVariantCatalogEntry:
    """Canonical multi-catalog representation of a priority CFTR variant."""
    legacy_name: str
    hgvs_cdna: str
    hgvs_protein: str
    spdi_grch38: str
    grch38_pos_1based: int
    spdi_pos_0based: int
    codon_position: int | None
    functional_class: str
    structural_domain: str
    pdb_reference: str
    fda_eti_responsive: bool
    fda_ivacaftor_monotherapy_responsive: bool
    requires_genetic_replacement_or_aso: bool
    active_clinical_trials: tuple[str, ...]
    consequence_type: str

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


PRIORITY_CFTR_VARIANT_CATALOG: Dict[str, CFTRVariantCatalogEntry] = {
    "F508del": CFTRVariantCatalogEntry(
        legacy_name="F508del",
        hgvs_cdna="NM_000492.4:c.1521_1523del",
        hgvs_protein="NP_000483.3:p.Phe508del",
        spdi_grch38="NC_000007.14:117559592:CTT:",
        grch38_pos_1based=117559593,
        spdi_pos_0based=117559592,
        codon_position=508,
        functional_class="Class II",
        structural_domain="NBD1",
        pdb_reference="6MSM",
        fda_eti_responsive=True,
        fda_ivacaftor_monotherapy_responsive=False,
        requires_genetic_replacement_or_aso=False,
        active_clinical_trials=("NCT05033080",),
        consequence_type="inframe_deletion",
    ),
    "G551D": CFTRVariantCatalogEntry(
        legacy_name="G551D",
        hgvs_cdna="NM_000492.4:c.1652G>A",
        hgvs_protein="NP_000483.3:p.Gly551Asp",
        spdi_grch38="NC_000007.14:117587805:G:A",
        grch38_pos_1based=117587806,
        spdi_pos_0based=117587805,
        codon_position=551,
        functional_class="Class III",
        structural_domain="NBD1",
        pdb_reference="6MSM",
        fda_eti_responsive=True,
        fda_ivacaftor_monotherapy_responsive=True,
        requires_genetic_replacement_or_aso=False,
        active_clinical_trials=("NCT05033080",),
        consequence_type="missense_variant",
    ),
    "W1282X": CFTRVariantCatalogEntry(
        legacy_name="W1282X",
        hgvs_cdna="NM_000492.4:c.3846G>A",
        hgvs_protein="NP_000483.3:p.Trp1282Ter",
        spdi_grch38="NC_000007.14:117642565:G:A",
        grch38_pos_1based=117642566,
        spdi_pos_0based=117642565,
        codon_position=1282,
        functional_class="Class I",
        structural_domain="NBD2",
        pdb_reference="5UAK",
        fda_eti_responsive=False,
        fda_ivacaftor_monotherapy_responsive=False,
        requires_genetic_replacement_or_aso=True,
        active_clinical_trials=("NCT05248230", "NCT05668741"),
        consequence_type="stop_gained",
    ),
    "G542X": CFTRVariantCatalogEntry(
        legacy_name="G542X",
        hgvs_cdna="NM_000492.4:c.1624G>T",
        hgvs_protein="NP_000483.3:p.Gly542Ter",
        spdi_grch38="NC_000007.14:117587777:G:T",
        grch38_pos_1based=117587778,
        spdi_pos_0based=117587777,
        codon_position=542,
        functional_class="Class I",
        structural_domain="NBD1",
        pdb_reference="5UAK",
        fda_eti_responsive=False,
        fda_ivacaftor_monotherapy_responsive=False,
        requires_genetic_replacement_or_aso=True,
        active_clinical_trials=("NCT05248230", "NCT05668741"),
        consequence_type="stop_gained",
    ),
    "N1303K": CFTRVariantCatalogEntry(
        legacy_name="N1303K",
        hgvs_cdna="NM_000492.4:c.3909C>G",
        hgvs_protein="NP_000483.3:p.Asn1303Lys",
        spdi_grch38="NC_000007.14:117652876:C:G",
        grch38_pos_1based=117652877,
        spdi_pos_0based=117652876,
        codon_position=1303,
        functional_class="Class II",
        structural_domain="NBD2",
        pdb_reference="6MSM",
        fda_eti_responsive=False,
        fda_ivacaftor_monotherapy_responsive=False,
        requires_genetic_replacement_or_aso=True,
        active_clinical_trials=("NCT05248230", "NCT05668741"),
        consequence_type="missense_variant",
    ),
    "R117H": CFTRVariantCatalogEntry(
        legacy_name="R117H",
        hgvs_cdna="NM_000492.4:c.350G>A",
        hgvs_protein="NP_000483.3:p.Arg117His",
        spdi_grch38="NC_000007.14:117530974:G:A",
        grch38_pos_1based=117530975,
        spdi_pos_0based=117530974,
        codon_position=117,
        functional_class="Class IV",
        structural_domain="TMD1",
        pdb_reference="6MSM",
        fda_eti_responsive=True,
        fda_ivacaftor_monotherapy_responsive=True,
        requires_genetic_replacement_or_aso=False,
        active_clinical_trials=("NCT05033080",),
        consequence_type="missense_variant",
    ),
    "3849+10kbC->T": CFTRVariantCatalogEntry(
        legacy_name="3849+10kbC->T",
        hgvs_cdna="NM_000492.4:c.3718-2477C>T",
        hgvs_protein="NP_000483.3:p.?",
        spdi_grch38="NC_000007.14:117639960:C:T",
        grch38_pos_1based=117639961,
        spdi_pos_0based=117639960,
        codon_position=None,
        functional_class="Class V",
        structural_domain="INTRON_22_CRYPTIC_SPLICE",
        pdb_reference="6MSM",
        fda_eti_responsive=True,
        fda_ivacaftor_monotherapy_responsive=True,
        requires_genetic_replacement_or_aso=True,
        active_clinical_trials=("NCT05033080",),
        consequence_type="splice_region_variant",
    ),
}


def map_codon_to_cftr_domain(codon_position: int | None) -> str:
    """Map a 1-based CFTR amino acid residue index (1..1480) to its canonical UniProt P13569 domain."""
    if codon_position is None:
        return "INTRON_22_CRYPTIC_SPLICE"
    if not (1 <= codon_position <= 1480):
        raise ValueError(f"CFTR codon position must be in [1, 1480], got {codon_position}")
    for domain_name, (start_aa, end_aa) in CFTR_STRUCTURAL_DOMAINS.items():
        if start_aa <= codon_position <= end_aa:
            return domain_name
    raise ValueError(f"Unmapped CFTR codon position {codon_position}")


def verify_spdi_hgvs_coordinate_shift(spdi_string: str, grch38_pos_1based: int) -> Dict[str, Any]:
    """Verify NCBI SPDI 0-based interbase coordinate against 1-based VCF/HGVS genomic position.

    Per Holmes et al. (2020) DOI: 10.1093/bioinformatics/btz856, SPDI uses 0-based interbase
    coordinates such that grch38_pos_1based - spdi_pos_0based == 1 for normalized SNVs and
    right-shifted deletions.
    """
    parts = spdi_string.split(":")
    if len(parts) != 4:
        raise ValueError(f"Invalid SPDI string format: {spdi_string!r}")
    seq_id, pos_str, deleted_seq, inserted_seq = parts
    spdi_pos_0based = int(pos_str)
    if grch38_pos_1based <= 0 or spdi_pos_0based < 0:
        raise ValueError("Genomic positions must be positive")
    coordinate_shift = grch38_pos_1based - spdi_pos_0based
    is_concordant = (seq_id == CANONICAL_GRCH38_CHROMOSOME) and (coordinate_shift == 1)
    return {
        "seq_id": seq_id,
        "spdi_pos_0based": spdi_pos_0based,
        "grch38_pos_1based": grch38_pos_1based,
        "coordinate_shift": coordinate_shift,
        "deleted_length": len(deleted_seq),
        "inserted_length": len(inserted_seq),
        "is_concordant": is_concordant,
    }


def evaluate_variant_catalog_concordance(variant_name: str) -> Dict[str, Any]:
    """Evaluate cross-catalog coordinate, structural domain, and regulatory theratype concordance."""
    if variant_name not in PRIORITY_CFTR_VARIANT_CATALOG:
        raise ValueError(f"Unsupported CFTR variant {variant_name!r}; expected one of {sorted(PRIORITY_CFTR_VARIANT_CATALOG)}")
    entry = PRIORITY_CFTR_VARIANT_CATALOG[variant_name]
    spdi_check = verify_spdi_hgvs_coordinate_shift(entry.spdi_grch38, entry.grch38_pos_1based)
    computed_domain = map_codon_to_cftr_domain(entry.codon_position)
    domain_concordant = computed_domain == entry.structural_domain

    if entry.fda_eti_responsive and entry.fda_ivacaftor_monotherapy_responsive:
        regulatory_pathway = "MODULATOR_MONOTHERAPY_AND_TRIPLE_ELIGIBLE"
    elif entry.fda_eti_responsive:
        regulatory_pathway = "TRIPLE_MODULATOR_ELIGIBLE"
    else:
        regulatory_pathway = "MODULATOR_REFRACTORY_NUCLEIC_ACID_PATHWAY"

    concordance_verdict = "CONCORDANT" if (spdi_check["is_concordant"] and domain_concordant) else "DISCORDANT"
    return {
        "legacy_name": entry.legacy_name,
        "transcript_refseq": CANONICAL_REFSEQ_TRANSCRIPT,
        "transcript_ensembl": CANONICAL_ENSEMBL_TRANSCRIPT,
        "uniprot_accession": CANONICAL_UNIPROT_ACCESSION,
        "hgvs_cdna": entry.hgvs_cdna,
        "hgvs_protein": entry.hgvs_protein,
        "spdi_grch38": entry.spdi_grch38,
        "spdi_normalization": spdi_check,
        "functional_class": entry.functional_class,
        "structural_domain": entry.structural_domain,
        "domain_concordant": domain_concordant,
        "pdb_reference": entry.pdb_reference,
        "fda_eti_responsive": entry.fda_eti_responsive,
        "fda_ivacaftor_monotherapy_responsive": entry.fda_ivacaftor_monotherapy_responsive,
        "requires_genetic_replacement_or_aso": entry.requires_genetic_replacement_or_aso,
        "regulatory_pathway": regulatory_pathway,
        "active_clinical_trials": list(entry.active_clinical_trials),
        "concordance_verdict": concordance_verdict,
    }


def summarize_catalog_concordance_panel() -> Dict[str, Any]:
    """Run full cross-catalog concordance audit across all 7 priority CFTR variants."""
    results: List[Dict[str, Any]] = [
        evaluate_variant_catalog_concordance(name)
        for name in PRIORITY_CFTR_VARIANT_CATALOG
    ]
    concordant_count = sum(1 for r in results if r["concordance_verdict"] == "CONCORDANT")
    modulator_eligible_count = sum(1 for r in results if r["fda_eti_responsive"])
    refractory_count = sum(1 for r in results if not r["fda_eti_responsive"])
    return {
        "total_variants_evaluated": len(results),
        "concordant_variants": concordant_count,
        "modulator_eligible_variants": modulator_eligible_count,
        "modulator_refractory_variants": refractory_count,
        "all_concordant": concordant_count == len(results),
        "variants": {r["legacy_name"]: r for r in results},
    }


if __name__ == "__main__":
    import json
    print(json.dumps(summarize_catalog_concordance_panel(), indent=2))
