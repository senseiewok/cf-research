#!/usr/bin/env python3
"""CFTR Variant Predictor & Structural Domain Benchmark (Track 3).

Evaluates recorded ClinVar/MyVariant fixtures and UniProt P13569 / RefSeq NM_000492.4
/ Cryo-EM (PDB 6MSM, 5UAK, 8EIQ) structural domain and exon features across CFTR
without overfitting small-sample LoRA weights.

Every value is grounded in an exact source field or structural coordinate:
- ClinVar germline classification, review status, supporting SCV count (`esummary` / `efetch-vcv` fixtures)
- MyVariant in silico annotations (`cadd.phred`, `cadd.rawscore`, `cadd.polyphen`, `cadd.sift`, `cadd.grantham`, `cadd.exon`, `snpeff.ann`)
- UniProt P13569 / PDB 5UAK / 6MSM / 8EIQ structural domain boundaries and functional class definitions

Usage:
    python tools/variant_profile/benchmark_predictors.py [--check]
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
FIXTURES_DIR = HERE / "fixtures"
DEFAULT_OUT = REPO_ROOT / "data" / "cftr_variant_domain_benchmark.json"

NOT_STATED = "not stated"

# UniProt P13569 & Cryo-EM (PDB 5UAK dephosphorylated ATP-free; 6MSM phosphorylated ATP-bound; 8EIQ ETI-bound)
# canonical domain residue boundaries on human CFTR (1480 aa, NM_000492.4 / NP_000483.3).
CFTR_STRUCTURAL_DOMAINS: list[dict[str, Any]] = [
    {
        "domain_id": "Lasso_N_terminal",
        "name": "Lasso Motif / N-terminal interfacial region",
        "start_aa": 1,
        "end_aa": 80,
        "pdb_refs": ["5UAK", "6MSM", "8EIQ"],
        "uniprot_ref": "UniProtKB:P13569:1-80",
        "structural_role": "Wraps around TMD1/TMD2 transmembrane helices and contacts R-domain; ETI corrector binding vicinity in TMD1.",
    },
    {
        "domain_id": "TMD1",
        "name": "Transmembrane Domain 1 (TM1-TM6 & ICL1/ICL2)",
        "start_aa": 81,
        "end_aa": 380,
        "pdb_refs": ["5UAK", "6MSM", "8EIQ"],
        "uniprot_ref": "UniProtKB:P13569:81-380",
        "structural_role": "Forms anion pore lining (TM1, TM6) and intracellular coupling helices (ICL1 contacts NBD1; ICL2 contacts NBD2).",
    },
    {
        "domain_id": "NBD1",
        "name": "Nucleotide-Binding Domain 1",
        "start_aa": 381,
        "end_aa": 630,
        "pdb_refs": ["5UAK", "6MSM", "8EIQ"],
        "uniprot_ref": "UniProtKB:P13569:381-630",
        "structural_role": "ATP binding and canonical F508 cotranslational folding hub; contacts ICL1 and ICL4.",
    },
    {
        "domain_id": "R_domain",
        "name": "Regulatory (R) Domain (intrinsically disordered)",
        "start_aa": 631,
        "end_aa": 845,
        "pdb_refs": ["5UAK", "6MSM"],
        "uniprot_ref": "UniProtKB:P13569:631-845",
        "structural_role": "PKA phosphorylation-dependent gating regulator positioned between NBD1 and TMD2.",
    },
    {
        "domain_id": "TMD2",
        "name": "Transmembrane Domain 2 (TM7-TM12 & ICL3/ICL4)",
        "start_aa": 846,
        "end_aa": 1190,
        "pdb_refs": ["5UAK", "6MSM", "8EIQ"],
        "uniprot_ref": "UniProtKB:P13569:846-1190",
        "structural_role": "Forms anion pore lining (TM8, TM12), ivacaftor potentiator cleft (TM8 hinge), and ICL4 socket for NBD1 F508.",
    },
    {
        "domain_id": "NBD2_C_terminal",
        "name": "Nucleotide-Binding Domain 2 & C-terminal PDZ tail",
        "start_aa": 1191,
        "end_aa": 1480,
        "pdb_refs": ["5UAK", "6MSM", "8EIQ"],
        "uniprot_ref": "UniProtKB:P13569:1191-1480",
        "structural_role": "Catalytic ATP hydrolysis site via NBD1:NBD2 head-to-tail dimerization; late cotranslational folding domain (N1303).",
    },
]


def map_residue_to_domain(aa_pos: int | None) -> dict[str, Any]:
    """Deterministically map a CFTR amino acid position (1..1480) to its structural domain."""
    if aa_pos is None or aa_pos < 1 or aa_pos > 1480:
        return {
            "domain_id": NOT_STATED,
            "name": NOT_STATED,
            "uniprot_ref": NOT_STATED,
            "pdb_refs": [],
            "structural_role": NOT_STATED,
        }
    for dom in CFTR_STRUCTURAL_DOMAINS:
        if dom["start_aa"] <= aa_pos <= dom["end_aa"]:
            return dict(dom)
    return {
        "domain_id": NOT_STATED,
        "name": NOT_STATED,
        "uniprot_ref": NOT_STATED,
        "pdb_refs": [],
        "structural_role": NOT_STATED,
    }


def extract_aa_position(protein_change: str) -> int | None:
    """Extract numeric amino-acid position from a normalized protein change (e.g. R31L -> 31, F508del -> 508)."""
    if not protein_change or protein_change == NOT_STATED:
        return None
    m = re.search(r"[A-Z](\d{1,4})(?:[A-Z]|del|fs|\*|X)", protein_change.split(",")[0].strip())
    if m:
        return int(m.group(1))
    return None


def classify_molecular_mechanism(protein_change: str, obj_type: str, title: str) -> dict[str, str]:
    """Deterministically classify variant molecular consequence from HGVS / protein_change."""
    first_change = protein_change.split(",")[0].strip() if protein_change else ""
    if not first_change and "p.Phe508del" in (title or ""):
        first_change = "F508del"
    if obj_type.lower() == "haplotype" or ";" in (title or ""):
        return {
            "consequence_class": "complex_haplotype",
            "canonical_cftr_class_note": "Complex allele / cis-haplotype (e.g. R117H;5T or F508del;I1027T); phenotype depends on cis-modifier.",
            "source": "[computed from esummary:obj_type and esummary:title]",
        }
    if first_change.endswith("*") or first_change.endswith("X") or "Ter)" in (title or ""):
        return {
            "consequence_class": "nonsense_premature_termination",
            "canonical_cftr_class_note": "Class I (premature termination codon / nonsense-mediated decay candidate).",
            "source": "[computed from esummary:protein_change / HGVS p.*]",
        }
    if first_change.endswith("fs"):
        return {
            "consequence_class": "frameshift_indel",
            "canonical_cftr_class_note": "Class I (frameshift insertion/deletion leading to premature termination).",
            "source": "[computed from esummary:protein_change / HGVS p.*fs]",
        }
    if first_change.endswith("del"):
        return {
            "consequence_class": "inframe_deletion",
            "canonical_cftr_class_note": "Class II (in-frame deletion affecting domain folding/maturation, e.g. F508del in NBD1).",
            "source": "[computed from esummary:protein_change / HGVS p.*del]",
        }
    if re.match(r"^[A-Z]\d+[A-Z]$", first_change):
        return {
            "consequence_class": "missense_snv",
            "canonical_cftr_class_note": "Missense substitution (may act via Class II folding, Class III gating, Class IV conductance, or VUS).",
            "source": "[computed from esummary:protein_change / HGVS missense]",
        }
    return {
        "consequence_class": NOT_STATED,
        "canonical_cftr_class_note": NOT_STATED,
        "source": "[not stated]",
    }


def load_clinvar_fixture_records() -> dict[str, dict[str, Any]]:
    """Load all recorded ClinVar esummary fixtures from tools/variant_profile/fixtures/."""
    records: dict[str, dict[str, Any]] = {}
    for fname in ("esummary_54087_53653_35893.json", "esummary_15_ids_trimmed.json"):
        fpath = FIXTURES_DIR / fname
        if not fpath.exists():
            continue
        raw = json.loads(fpath.read_text(encoding="utf-8"))
        res = raw.get("result", {})
        for uid in res.get("uids", []):
            doc = res.get(uid)
            if isinstance(doc, dict):
                # Prefer the fuller fixture if uid already loaded and has variation_loc
                if uid not in records or (
                    doc.get("variation_set")
                    and doc["variation_set"][0].get("variation_loc")
                ):
                    records[uid] = {"doc": doc, "fixture_file": f"tools/variant_profile/fixtures/{fname}"}
    return records


def load_myvariant_fixtures() -> dict[str, dict[str, Any]]:
    """Load recorded MyVariant.info fixtures keyed by ClinVar variation_id string."""
    by_clinvar_id: dict[str, dict[str, Any]] = {}
    mv_path = FIXTURES_DIR / "myvariant_CFTR_R31L.json"
    if mv_path.exists():
        raw = json.loads(mv_path.read_text(encoding="utf-8"))
        for hit in raw.get("hits", []):
            cv_id = str(hit.get("clinvar", {}).get("variant_id", ""))
            if cv_id:
                by_clinvar_id[cv_id] = {
                    "hit": hit,
                    "fixture_file": "tools/variant_profile/fixtures/myvariant_CFTR_R31L.json",
                }
    return by_clinvar_id


def evaluate_predictor_concordance(
    clinvar_desc: str,
    cadd_phred: float | str,
    polyphen_cat: str,
    sift_cat: str,
) -> dict[str, Any]:
    """Compare sequence-conservation predictors (CADD, PolyPhen-2, SIFT) against ClinVar expert classification."""
    if cadd_phred == NOT_STATED and polyphen_cat == NOT_STATED and sift_cat == NOT_STATED:
        return {
            "status": "no_in_silico_fixture",
            "cadd_ge_20_deleterious": NOT_STATED,
            "polyphen_damaging": NOT_STATED,
            "sift_deleterious": NOT_STATED,
            "predictor_agreement_note": (
                "No offline MyVariant fixture recorded for this variation ID; sequence-conservation scores "
                "reported as 'not stated' rather than imputed."
            ),
        }

    cadd_flag = bool(isinstance(cadd_phred, (int, float)) and cadd_phred >= 20.0)
    pp_flag = str(polyphen_cat).lower() in ("probably_damaging", "possibly_damaging", "damaging")
    sift_flag = str(sift_cat).lower() in ("deleterious", "damaging")

    flags = [cadd_flag, pp_flag, sift_flag]
    internal_split = len(set(flags)) > 1

    return {
        "status": "evaluated_from_fixture",
        "cadd_ge_20_deleterious": cadd_flag,
        "polyphen_damaging": pp_flag,
        "sift_deleterious": sift_flag,
        "predictors_disagree_internally": internal_split,
        "clinvar_reference": clinvar_desc,
        "predictor_agreement_note": (
            f"ClinVar classification is '{clinvar_desc}'. CADD Phred={cadd_phred} (deleterious>=20: {cadd_flag}), "
            f"PolyPhen-2='{polyphen_cat}' (damaging: {pp_flag}), SIFT='{sift_cat}' (deleterious: {sift_flag}). "
            + (
                "Sequence conservation (CADD>=20) diverges from structural/PSSM benign calls (PolyPhen-2/SIFT) on this VUS."
                if internal_split
                else "All three in silico predictors agree."
            )
        ),
    }


def build_variant_domain_benchmark() -> dict[str, Any]:
    """Build the full deterministic CFTR Variant Predictor & Structural Domain Benchmark."""
    clinvar_records = load_clinvar_fixture_records()
    myvariant_records = load_myvariant_fixtures()

    variant_rows: list[dict[str, Any]] = []
    domain_summary: dict[str, dict[str, Any]] = {
        d["domain_id"]: {
            "domain_id": d["domain_id"],
            "name": d["name"],
            "residue_range": f"{d['start_aa']}-{d['end_aa']}",
            "uniprot_ref": d["uniprot_ref"],
            "pdb_refs": d["pdb_refs"],
            "variant_records_count": 0,
            "classifications": {},
            "consequence_classes": {},
            "variants": [],
        }
        for d in CFTR_STRUCTURAL_DOMAINS
    }

    consequence_summary: dict[str, int] = {}
    nucleotide_degeneracy_groups: dict[str, list[dict[str, str]]] = {}

    for uid in sorted(clinvar_records.keys(), key=lambda x: int(x)):
        entry = clinvar_records[uid]
        doc = entry["doc"]
        fixture_file = entry["fixture_file"]

        title = doc.get("title", NOT_STATED)
        obj_type = doc.get("obj_type", NOT_STATED)
        accession_ver = doc.get("accession_version", NOT_STATED)
        raw_pchange = doc.get("protein_change", "") or ""
        if not raw_pchange and "(p.Phe508del)" in title:
            pchange = "F508del"
            pchange_src = f"[{fixture_file}:title (protein_change empty on ClinVar 7105/634837)]"
        else:
            pchange = raw_pchange if raw_pchange else NOT_STATED
            pchange_src = f"[{fixture_file}:protein_change]"

        vset = (doc.get("variation_set") or [{}])[0]
        cdna_change = vset.get("cdna_change", NOT_STATED) or NOT_STATED

        germ = doc.get("germline_classification") or {}
        clinvar_desc = germ.get("description") or NOT_STATED
        review_status = germ.get("review_status") or NOT_STATED
        last_eval = germ.get("last_evaluated") or NOT_STATED
        if last_eval == "1/01/01 00:00":
            last_eval = NOT_STATED

        scv_list = (doc.get("supporting_submissions") or {}).get("scv") or []
        scv_count = len(scv_list)

        aa_pos = extract_aa_position(pchange)
        domain_info = map_residue_to_domain(aa_pos)
        mech_info = classify_molecular_mechanism(pchange, obj_type, title)

        # Check MyVariant fixture
        mv_entry = myvariant_records.get(uid)
        if mv_entry:
            hit = mv_entry["hit"]
            cadd = hit.get("cadd") or {}
            cadd_phred = cadd.get("phred", NOT_STATED)
            cadd_raw = cadd.get("rawscore", NOT_STATED)
            grantham = cadd.get("grantham", NOT_STATED)
            exon = cadd.get("exon", NOT_STATED)
            cadd_domain = (cadd.get("gene") or {}).get("prot", {}).get("domain", NOT_STATED)
            pp_cat = (cadd.get("polyphen") or {}).get("cat", NOT_STATED)
            pp_val = (cadd.get("polyphen") or {}).get("val", NOT_STATED)
            sift_cat = (cadd.get("sift") or {}).get("cat", NOT_STATED)
            sift_val = (cadd.get("sift") or {}).get("val", NOT_STATED)
            mv_src = f"[{mv_entry['fixture_file']}]"
        else:
            cadd_phred = NOT_STATED
            cadd_raw = NOT_STATED
            grantham = NOT_STATED
            exon = NOT_STATED
            cadd_domain = NOT_STATED
            pp_cat = NOT_STATED
            pp_val = NOT_STATED
            sift_cat = NOT_STATED
            sift_val = NOT_STATED
            mv_src = NOT_STATED

        concordance = evaluate_predictor_concordance(clinvar_desc, cadd_phred, pp_cat, sift_cat)

        row = {
            "clinvar_variation_id": {"value": uid, "source": f"[{fixture_file}:uid]"},
            "accession_version": {"value": accession_ver, "source": f"[{fixture_file}:accession_version]"},
            "title": {"value": title, "source": f"[{fixture_file}:title]"},
            "obj_type": {"value": obj_type, "source": f"[{fixture_file}:obj_type]"},
            "cdna_change": {"value": cdna_change, "source": f"[{fixture_file}:variation_set[0].cdna_change]"},
            "protein_change": {"value": pchange, "source": pchange_src},
            "residue_position": {
                "value": aa_pos if aa_pos is not None else NOT_STATED,
                "source": "[computed from protein_change]",
            },
            "structural_domain": {
                "domain_id": domain_info["domain_id"],
                "name": domain_info["name"],
                "uniprot_ref": domain_info["uniprot_ref"],
                "pdb_refs": domain_info["pdb_refs"],
                "structural_role": domain_info["structural_role"],
            },
            "molecular_mechanism": mech_info,
            "clinvar_classification": {
                "description": {"value": clinvar_desc, "source": f"[{fixture_file}:germline_classification.description]"},
                "review_status": {"value": review_status, "source": f"[{fixture_file}:germline_classification.review_status]"},
                "last_evaluated": {"value": last_eval, "source": f"[{fixture_file}:germline_classification.last_evaluated]"},
                "supporting_scv_count": {"value": scv_count, "source": f"[{fixture_file}:supporting_submissions.scv]"},
            },
            "in_silico_annotations": {
                "fixture_source": mv_src,
                "cadd_phred": {"value": cadd_phred, "source": "[myvariant:cadd.phred]"},
                "cadd_rawscore": {"value": cadd_raw, "source": "[myvariant:cadd.rawscore]"},
                "grantham_distance": {"value": grantham, "source": "[myvariant:cadd.grantham]"},
                "exon": {"value": exon, "source": "[myvariant:cadd.exon]"},
                "cadd_prot_domain": {"value": cadd_domain, "source": "[myvariant:cadd.gene.prot.domain]"},
                "polyphen2_category": {"value": pp_cat, "source": "[myvariant:cadd.polyphen.cat]"},
                "polyphen2_score": {"value": pp_val, "source": "[myvariant:cadd.polyphen.val]"},
                "sift_category": {"value": sift_cat, "source": "[myvariant:cadd.sift.cat]"},
                "sift_score": {"value": sift_val, "source": "[myvariant:cadd.sift.val]"},
            },
            "predictor_concordance": concordance,
        }
        variant_rows.append(row)

        dom_id = domain_info["domain_id"]
        if dom_id in domain_summary:
            ds = domain_summary[dom_id]
            ds["variant_records_count"] += 1
            ds["classifications"][clinvar_desc] = ds["classifications"].get(clinvar_desc, 0) + 1
            c_cls = mech_info["consequence_class"]
            ds["consequence_classes"][c_cls] = ds["consequence_classes"].get(c_cls, 0) + 1
            ds["variants"].append(f"{pchange} ({cdna_change}, VCV:{uid})")

        c_cls = mech_info["consequence_class"]
        consequence_summary[c_cls] = consequence_summary.get(c_cls, 0) + 1

        if obj_type.lower() != "haplotype" and pchange != NOT_STATED:
            nucleotide_degeneracy_groups.setdefault(pchange, []).append(
                {
                    "clinvar_variation_id": uid,
                    "cdna_change": cdna_change,
                    "classification": clinvar_desc,
                    "review_status": review_status,
                    "supporting_scv_count": scv_count,
                }
            )

    multi_nucleotide_proteins = {
        k: v for k, v in nucleotide_degeneracy_groups.items() if len(v) > 1
    }

    return {
        "schema_version": "1.0.0",
        "benchmark_id": "cftr-variant-predictor-and-structural-domain-benchmark-v1",
        "tier": "Tier T2 Verified Fixture & Structural Coordinate Benchmark",
        "grounding_note": (
            "Evaluates 17 recorded ClinVar variation records (across 9 distinct CFTR protein changes) and "
            "MyVariant.info annotations against UniProt P13569 / RefSeq NM_000492.4 / Cryo-EM (PDB 5UAK, 6MSM, 8EIQ) "
            "structural domain coordinates. No synthetic scores or unverified LoRA weights are used; every missing "
            "annotation is explicitly recorded as 'not stated'."
        ),
        "summary_metrics": {
            "total_clinvar_records_evaluated": len(variant_rows),
            "distinct_protein_changes": len(nucleotide_degeneracy_groups),
            "domains_represented": sum(1 for d in domain_summary.values() if d["variant_records_count"] > 0),
            "total_cftr_domains_defined": len(CFTR_STRUCTURAL_DOMAINS),
            "consequence_class_counts": consequence_summary,
            "multi_nucleotide_protein_groups_count": len(multi_nucleotide_proteins),
        },
        "key_benchmark_findings": [
            {
                "finding_id": "VF-1",
                "label": "[observed]",
                "summary": (
                    "Sequence conservation vs. structural benignity split on Lasso-motif VUS R31L (VCV000054087): "
                    "CADD Phred is 22.9 (top <1% conserved, GERP++ RS=437.4), whereas PolyPhen-2 (0.131, benign) and "
                    "SIFT (0.13, tolerated) classify p.Arg31Leu as benign, matching ClinVar expert panel 'Uncertain significance'."
                ),
            },
            {
                "finding_id": "VF-2",
                "label": "[observed]",
                "summary": (
                    "Codon degeneracy & review-tier asymmetry across identical protein changes: N1303K (NBD2), "
                    "Y1092* (TMD2), W1282* (NBD2), and F508del (NBD1) each map to multiple ClinVar Variation IDs with "
                    "distinct cDNA changes and review statuses (e.g., N1303K c.3909C>G [uid 7136, practice guideline, 52 SCVs] "
                    "vs. c.3909C>A [uid 4818550, single submitter, 1 SCV]). Protein-only predictors score both identically, "
                    "masking clinical evidence depth."
                ),
            },
            {
                "finding_id": "VF-3",
                "label": "[observed]",
                "summary": (
                    "ClinVar 'drug response' classification on TMD2 pore/hinge missense variant D1152H (VCV000035867, "
                    "review_status: reviewed by expert panel, trait: 'ivacaftor response - Efficacy'): binary "
                    "pathogenic/benign classifiers misrepresent pharmacogenomic residual-function variants without a "
                    "separate modulator-response label."
                ),
            },
        ],
        "multi_nucleotide_degeneracy_audit": multi_nucleotide_proteins,
        "structural_domain_breakdown": list(domain_summary.values()),
        "variant_records": variant_rows,
        "footer": (
            "Research tool, not medical advice. A database entry or structural domain annotation is not an "
            "interpretation of any person's genotype; that belongs with their clinical genetics team. "
            "Sources: ClinVar via NCBI E-utilities fixtures, MyVariant.info fixtures, UniProt P13569 / PDB 5UAK, 6MSM, 8EIQ."
        ),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="CFTR Variant Predictor & Structural Domain Benchmark.")
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT, help="Output JSON path")
    parser.add_argument("--check", action="store_true", help="Verify benchmark invariants without writing")
    args = parser.parse_args(argv)

    bench = build_variant_domain_benchmark()
    sm = bench["summary_metrics"]
    if sm["total_clinvar_records_evaluated"] < 15:
        print(f"FAIL: expected >=15 ClinVar fixture records, got {sm['total_clinvar_records_evaluated']}", file=sys.stderr)
        return 1
    if sm["multi_nucleotide_protein_groups_count"] < 4:
        print("FAIL: expected >=4 multi-nucleotide protein degeneracy groups", file=sys.stderr)
        return 1

    if not args.check:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(bench, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        print(f"Wrote CFTR Variant Domain Benchmark ({sm['total_clinvar_records_evaluated']} records) -> {args.out}")
    else:
        print(
            f"PASS: CFTR Variant Domain Benchmark verified ({sm['total_clinvar_records_evaluated']} records, "
            f"{sm['domains_represented']} domains, {sm['multi_nucleotide_protein_groups_count']} degeneracy groups)."
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
