"""MyVariant.info variant annotation client (catalog entry `myvariant` in sources/catalog.yaml).

Queries MyVariant.info for in silico pathogenicity predictions (CADD, SIFT, PolyPhen-2)
and functional annotations (SnpEff) for a given gene and protein change.

Standard library only. Paced at <= 1.0 RPS.
"""

from __future__ import annotations

import argparse
import datetime
import json
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from typing import Any

from .variant_profile import (
    AA3,
    GENE_RE,
    ONE_LETTER_RE,
    THREE_LETTER_RE,
    ApiError,
    Refused,
    normalise_change,
    three_letter,
    validate,
)

BASE_URL = "https://myvariant.info/v1/"
CATALOG_ID = "myvariant"
MAX_RPS = 1.0
MIN_INTERVAL_S = 1.0 / MAX_RPS
TOOL_NAME = "senseiewok-myvariant"
USER_AGENT = TOOL_NAME + "/0.1 (+https://github.com/senseiewok/cf-research)"
TIMEOUT_S = 30
MAX_BYTES = 10_000_000
NOT_STATED = "not stated"
FOOTER = (
    "Research tool, not medical advice. A database entry is not an interpretation of any person's genotype; "
    "that belongs with their clinical genetics team. Source: MyVariant.info API, retrieved {date}."
)

_last_request = [None]


def _http_get(url: str) -> bytes:
    """The single place this module touches the network. Returns body bytes, paced at <= 1.0 RPS."""
    now = time.monotonic()
    if _last_request[0] is not None:
        elapsed = now - _last_request[0]
        if elapsed < MIN_INTERVAL_S:
            time.sleep(MIN_INTERVAL_S - elapsed)
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    try:
        with urllib.request.urlopen(req, timeout=TIMEOUT_S) as resp:
            body = resp.read(MAX_BYTES)
            _last_request[0] = time.monotonic()
            return body
    except urllib.error.HTTPError as e:
        _last_request[0] = time.monotonic()
        raise ApiError(f"HTTP {e.code}: {e.reason} for {url}") from e
    except Exception as e:
        _last_request[0] = time.monotonic()
        raise ApiError(f"request failed: {e} for {url}") from e


def query_myvariant(gene: str, change: str) -> dict[str, Any]:
    """Query MyVariant.info query API for gene and protein change."""
    gene, one = validate(gene, change)
    hgvs_p = three_letter(one)
    q = urllib.parse.quote(f"{gene} {hgvs_p}")
    url = f"{BASE_URL}query?q={q}&fields=cadd,clinvar,snpeff&size=5"
    raw = _http_get(url)
    try:
        return json.loads(raw.decode("utf-8"))
    except Exception as e:
        raise ApiError(f"failed to parse JSON from MyVariant: {e}") from e


def parse_myvariant_record(hit: dict[str, Any]) -> dict[str, Any]:
    """Extract structured fields from one MyVariant hit, tying every number to an exact source field."""
    record_id = hit.get("_id", NOT_STATED)
    cadd = hit.get("cadd", {}) if isinstance(hit.get("cadd"), dict) else {}
    clinvar = hit.get("clinvar", {}) if isinstance(hit.get("clinvar"), dict) else {}
    snpeff = hit.get("snpeff", {}) if isinstance(hit.get("snpeff"), dict) else {}
    snpeff_ann = snpeff.get("ann", {}) if isinstance(snpeff.get("ann"), dict) else {}

    polyphen = cadd.get("polyphen", {}) if isinstance(cadd.get("polyphen"), dict) else {}
    sift = cadd.get("sift", {}) if isinstance(cadd.get("sift"), dict) else {}
    hg38 = clinvar.get("hg38", {}) if isinstance(clinvar.get("hg38"), dict) else {}

    return {
        "variant_id_genomic": (record_id, "[myvariant:_id]"),
        "cadd_rawscore": (cadd.get("rawscore", NOT_STATED), "[myvariant:cadd.rawscore]"),
        "cadd_phred": (cadd.get("phred", NOT_STATED), "[myvariant:cadd.phred]"),
        "cadd_consequence": (cadd.get("consequence", NOT_STATED), "[myvariant:cadd.consequence]"),
        "polyphen_category": (polyphen.get("cat", NOT_STATED), "[myvariant:cadd.polyphen.cat]"),
        "polyphen_score": (polyphen.get("val", NOT_STATED), "[myvariant:cadd.polyphen.val]"),
        "sift_category": (sift.get("cat", NOT_STATED), "[myvariant:cadd.sift.cat]"),
        "sift_score": (sift.get("val", NOT_STATED), "[myvariant:cadd.sift.val]"),
        "clinvar_variant_id": (clinvar.get("variant_id", NOT_STATED), "[myvariant:clinvar.variant_id]"),
        "clinvar_rsid": (clinvar.get("rsid", NOT_STATED), "[myvariant:clinvar.rsid]"),
        "hg38_start": (hg38.get("start", NOT_STATED), "[myvariant:clinvar.hg38.start]"),
        "hg38_end": (hg38.get("end", NOT_STATED), "[myvariant:clinvar.hg38.end]"),
        "snpeff_effect": (snpeff_ann.get("effect", NOT_STATED), "[myvariant:snpeff.ann.effect]"),
        "snpeff_putative_impact": (snpeff_ann.get("putative_impact", NOT_STATED), "[myvariant:snpeff.ann.putative_impact]"),
        "snpeff_hgvs_c": (snpeff_ann.get("hgvs_c", NOT_STATED), "[myvariant:snpeff.ann.hgvs_c]"),
        "snpeff_hgvs_p": (snpeff_ann.get("hgvs_p", NOT_STATED), "[myvariant:snpeff.ann.hgvs_p]"),
    }


def render_myvariant_text(gene: str, change: str, parsed: dict[str, Any], now: datetime.datetime | None = None) -> str:
    now = now or datetime.datetime.now(datetime.timezone.utc)
    date_str = now.strftime("%Y-%m-%d (UTC)")
    lines = [
        f"MyVariant Annotation: {gene} {change}",
        f"Genomic ID: {parsed['variant_id_genomic'][0]}    {parsed['variant_id_genomic'][1]}",
        f"CADD Phred score: {parsed['cadd_phred'][0]}    {parsed['cadd_phred'][1]}",
        f"CADD raw score: {parsed['cadd_rawscore'][0]}    {parsed['cadd_rawscore'][1]}",
        f"CADD consequence: {parsed['cadd_consequence'][0]}    {parsed['cadd_consequence'][1]}",
        f"PolyPhen-2 prediction: {parsed['polyphen_category'][0]} {parsed['polyphen_category'][1]} (score: {parsed['polyphen_score'][0]} {parsed['polyphen_score'][1]})",
        f"SIFT prediction: {parsed['sift_category'][0]} {parsed['sift_category'][1]} (score: {parsed['sift_score'][0]} {parsed['sift_score'][1]})",
        f"SnpEff effect: {parsed['snpeff_effect'][0]} {parsed['snpeff_effect'][1]} (impact: {parsed['snpeff_putative_impact'][0]} {parsed['snpeff_putative_impact'][1]})",
        f"SnpEff HGVS c.: {parsed['snpeff_hgvs_c'][0]}    {parsed['snpeff_hgvs_c'][1]}",
        f"SnpEff HGVS p.: {parsed['snpeff_hgvs_p'][0]}    {parsed['snpeff_hgvs_p'][1]}",
        f"ClinVar variant ID: {parsed['clinvar_variant_id'][0]}    {parsed['clinvar_variant_id'][1]}",
        f"ClinVar rsID: {parsed['clinvar_rsid'][0]}    {parsed['clinvar_rsid'][1]}",
        f"GRCh38 coordinates: chr7:{parsed['hg38_start'][0]}-{parsed['hg38_end'][0]}    {parsed['hg38_start'][1]}",
        "",
        FOOTER.format(date=date_str),
    ]
    return "\n".join(lines)
