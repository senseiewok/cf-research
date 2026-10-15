"""Variant layer unifying ClinVar, MyVariant, and source catalog status.

Fulfills Step 2 of T-0102:
Gives one profile for a fixed public variant with every number tied to a source field.
Sources:
- ClinVar via NCBI E-utilities (catalog ids `ncbi-eutils`, `clinvar`)
- MyVariant.info API (catalog id `myvariant`)
- AlphaMissense (catalog id `alphamissense-predictions`, access: manual; not admitted for API)

Usage:
    python variant_layer.py CFTR R31L [--json]
"""

from __future__ import annotations

import argparse
import datetime
import json
import sys
from typing import Any

from . import myvariant, variant_profile

NOT_STATED = "not stated"
FOOTER = (
    "Research tool, not medical advice. A database entry is not an interpretation of any person's genotype; "
    "that belongs with their clinical genetics team. Sources: ClinVar via NCBI E-utilities, MyVariant.info API; "
    "retrieved {date}."
)


def run_variant_layer(gene: str, change: str, now: datetime.datetime | None = None) -> dict[str, Any]:
    """Run variant layer returning a single unified profile with every number tied to an exact source field."""
    now = now or datetime.datetime.now(datetime.timezone.utc)
    date_str = now.strftime("%Y-%m-%d (UTC)")

    # 1. Run ClinVar via NCBI E-utilities
    query_clinvar, profiles_clinvar = variant_profile.run(gene, change, want_submitters=True)
    first_profile = profiles_clinvar[0] if profiles_clinvar else {}

    # 2. Run MyVariant.info
    raw_mv = myvariant.query_myvariant(gene, change)
    hits = raw_mv.get("hits", [])
    first_hit = hits[0] if hits else {}
    parsed_mv = myvariant.parse_myvariant_record(first_hit)

    # 3. Cross-check between ClinVar NCBI and MyVariant ClinVar block
    ncbi_var_id = first_profile.get("clinvar_variation_id", NOT_STATED)
    mv_var_id = parsed_mv.get("clinvar_variant_id", (NOT_STATED,))[0]

    if ncbi_var_id != NOT_STATED and mv_var_id != NOT_STATED:
        if str(ncbi_var_id) == str(mv_var_id):
            cross_check = {
                "verdict": "PASS",
                "check": "ClinVar Variation ID Agreement",
                "detail": f"NCBI Variation ID {ncbi_var_id} matches MyVariant ClinVar ID {mv_var_id}",
            }
        else:
            cross_check = {
                "verdict": "CHECK",
                "check": "ClinVar Variation ID Disagreement",
                "detail": f"NCBI Variation ID {ncbi_var_id} differs from MyVariant ClinVar ID {mv_var_id}",
            }
    else:
        cross_check = {
            "verdict": "NOT STATED",
            "check": "ClinVar Variation ID Cross-check",
            "detail": f"One or both sources omitted Variation ID (NCBI={ncbi_var_id}, MyVariant={mv_var_id})",
        }

    # 4. AlphaMissense Catalog Status
    alphamissense_status = {
        "catalog_id": "alphamissense-predictions",
        "access": "manual",
        "admitted_for_api": False,
        "detail": (
            "AlphaMissense predictions are catalogued in sources/catalog.yaml with access: manual (bulk Zenodo download). "
            "No automated API route is admitted by the maintainer; unadmitted sources are not reached."
        ),
    }

    unified = {
        "variant": {"gene": gene, "protein_change": change},
        "query_metadata": {
            "retrieved_utc": now.strftime("%Y-%m-%dT%H:%M:%SZ"),
            "match_status": query_clinvar.get("match_status", NOT_STATED),
            "clinvar_records_count": len(profiles_clinvar),
            "myvariant_hits_count": raw_mv.get("total", 0),
        },
        "clinvar": first_profile,
        "in_silico_myvariant": parsed_mv,
        "alphamissense": alphamissense_status,
        "cross_checks": [cross_check],
        "footer": FOOTER.format(date=date_str),
    }
    return unified


def render_variant_layer_text(unified: dict[str, Any]) -> str:
    v = unified["variant"]
    qm = unified["query_metadata"]
    cv = unified["clinvar"]
    mv = unified["in_silico_myvariant"]
    am = unified["alphamissense"]
    cc = unified["cross_checks"][0]

    lines = [
        f"=== Unified Variant Layer Profile: {v['gene']} {v['protein_change']} ===",
        f"Match status: {qm['match_status']} (ClinVar matches: {qm['clinvar_records_count']}, MyVariant hits: {qm['myvariant_hits_count']})",
        "",
        "--- 1. Clinical Evidence (ClinVar via NCBI E-utilities) ---",
        f"  Variation ID: {cv.get('clinvar_variation_id', NOT_STATED)}    [esummary:uid]",
        f"  Accession: {cv.get('accession_version', NOT_STATED)}    [esummary:accession_version]",
        f"  Title: {cv.get('title', NOT_STATED)}    [esummary:title]",
        f"  Germline classification: {cv.get('germline_classification', NOT_STATED)}    [esummary:germline_classification.description]",
        f"  Review status: {cv.get('review_status', NOT_STATED)}    [esummary:germline_classification.review_status]",
        f"  Date last evaluated: {cv.get('date_last_evaluated', NOT_STATED)}    [esummary:germline_classification.last_evaluated]",
        f"  Number of submitters: {cv.get('number_of_submitters', NOT_STATED)}    [efetch-vcv:VariationArchive/@NumberOfSubmitters]",
        f"  Submitter agreement: {cv.get('submitters_disagree_flag', NOT_STATED)}    [computed from efetch-vcv per-submitter classifications]",
        "",
        "--- 2. In Silico & Functional Annotations (MyVariant.info) ---",
        f"  Genomic ID: {mv['variant_id_genomic'][0]}    {mv['variant_id_genomic'][1]}",
        f"  CADD Phred score: {mv['cadd_phred'][0]}    {mv['cadd_phred'][1]}",
        f"  CADD raw score: {mv['cadd_rawscore'][0]}    {mv['cadd_rawscore'][1]}",
        f"  CADD consequence: {mv['cadd_consequence'][0]}    {mv['cadd_consequence'][1]}",
        f"  PolyPhen-2 prediction: {mv['polyphen_category'][0]} {mv['polyphen_category'][1]} (score: {mv['polyphen_score'][0]} {mv['polyphen_score'][1]})",
        f"  SIFT prediction: {mv['sift_category'][0]} {mv['sift_category'][1]} (score: {mv['sift_score'][0]} {mv['sift_score'][1]})",
        f"  SnpEff effect: {mv['snpeff_effect'][0]} {mv['snpeff_effect'][1]} (impact: {mv['snpeff_putative_impact'][0]} {mv['snpeff_putative_impact'][1]})",
        f"  GRCh38 coordinates: chr7:{mv['hg38_start'][0]}-{mv['hg38_end'][0]}    {mv['hg38_start'][1]}",
        "",
        "--- 3. Structural Prediction Layer (AlphaMissense) ---",
        f"  Catalog ID: {am['catalog_id']} (access: {am['access']})",
        f"  Status: {am['detail']}",
        "",
        "--- 4. Cross-checks ---",
        f"  {cc['verdict']:<10} {cc['check']}: {cc['detail']}",
        "",
        unified["footer"],
    ]
    return "\n".join(lines)


def main(argv=None, out=sys.stdout, err=sys.stderr):
    ap = argparse.ArgumentParser(description="Unified variant layer profile (research tool, not medical advice).")
    ap.add_argument("gene", help="gene symbol, e.g. CFTR")
    ap.add_argument("protein_change", help="protein change, e.g. R31L, F508del, G542X")
    ap.add_argument("--json", action="store_true", help="output JSON")
    args = ap.parse_args(argv)

    try:
        profile = run_variant_layer(args.gene, args.protein_change)
        if args.json:
            print(json.dumps(profile, indent=2, ensure_ascii=False), file=out)
        else:
            print(render_variant_layer_text(profile), file=out)
        return 0
    except variant_profile.Refused as e:
        print(f"refused: {e}", file=err)
        return 2
    except variant_profile.ApiError as e:
        print(f"API error, not retried: {e}", file=err)
        return 3
    except Exception as e:
        print(f"unexpected error: {e}", file=err)
        return 3


if __name__ == "__main__":
    sys.exit(main())
