"""Print what ClinVar states about one protein variant in one gene, field by field, with the API field each value came from.

Variant-level only: it takes a gene symbol and a protein change (CFTR R31L, CFTR Arg31Leu) and nothing else, and refuses anything that looks like a
person (a name, a date of birth, a record number). It reads ClinVar through the documented NCBI E-utilities API: esearch, then esummary, then
efetch rettype=vcv for the per-submitter rows (skip that with --no-submitters). It never crawls web pages. A value the response does not
contain prints as "not stated"; nothing is filled in from elsewhere.

Usage: python variant_profile.py GENE PROTEIN_CHANGE [--json] [--no-submitters]
Exit 0: printed a profile. 1: no ClinVar record has that protein change. 2: refused input. 3: the API failed or answered unexpectedly.
Conduct (catalog entries `ncbi-eutils` and `clinvar` in sources/catalog.yaml): at most 1 request per second, one request per call and no retry,
a descriptive User-Agent and the `tool` parameter, no email address or other personal identifier, no API key, ClinVar credited as the source."""
import argparse
import datetime
import json
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"  # catalog `ncbi-eutils` base_url; the only host this tool calls
CATALOG_IDS = ("ncbi-eutils", "clinvar")
NCBI_PUBLISHED_RPS = 3.0  # NCBI's stated ceiling without an API key (catalog `ncbi-eutils` published_limit)
MAX_RPS = 1.0  # our cap, lower than the catalog's 2.5 and NCBI's 3
MIN_INTERVAL_S = 1.0 / MAX_RPS
TOOL_NAME = "senseiewok-variant-profile"  # sent as tool=; NCBI's email= parameter is deliberately not sent
USER_AGENT = TOOL_NAME + "/0.1 (+https://github.com/senseiewok/cf-research)"
TIMEOUT_S = 30
MAX_BYTES = 5_000_000
ESEARCH_RETMAX = 20
NOT_STATED = "not stated"
NO_CALL = "not provided"  # ClinVar's text for a submission that makes no classification; counted, but not a disagreement
EMPTY_DATE = "1/01/01 00:00"  # esummary's placeholder for an empty date
FOOTER = ("Research tool, not medical advice. A database entry is not an interpretation of any person's genotype; that belongs with their "
          "clinical genetics team. Source: ClinVar via NCBI E-utilities, retrieved {date}.")

SUBMITTER_FIELDS = ("scv", "submitter", "classification", "date_last_evaluated", "review_status")

GENE_RE =re.compile(r"^[A-Z0-9-]{2,12}$")
ONE_LETTER_RE = re.compile(r"^[A-Z][0-9]{1,5}([A-Z]|del|fs|X|\*)$")
AA3 = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I", "Leu": "L",
       "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V", "Ter": "X"}
_AA3_ALT = "|".join(AA3)
THREE_LETTER_RE = re.compile(r"^(" + _AA3_ALT + r")([0-9]{1,5})(" + _AA3_ALT + r"|del|fs|\*)$")


class Refused(ValueError):
    """Input this tool will not take (exit 2)."""


class ApiError(RuntimeError):
    """The API failed or answered with something this tool does not recognise (exit 3)."""


def normalise_change(change):
    """The one-letter form ClinVar's esummary uses in `protein_change` (R31L, F508del, G542X). Raises Refused for anything else."""
    if ONE_LETTER_RE.match(change):
        return change.replace("*", "X")
    m = THREE_LETTER_RE.match(change)
    if m:
        ref, pos, alt = m.groups()
        return AA3[ref] + pos + (AA3.get(alt) or alt.replace("*", "X"))
    raise Refused("protein change must look like R31L, F508del, G542X or Arg31Leu; this tool takes a variant, never personal details")


def validate(gene, change):
    """(gene, one-letter change) or Refused. Deliberately narrow: a name, a date of birth or a record number cannot pass either pattern."""
    if not GENE_RE.match(gene or ""):
        raise Refused("gene must be an upper-case symbol such as CFTR (2 to 12 of A-Z, 0-9, -); this tool takes a variant, never personal details")
    return gene, normalise_change(change or "")


# --- network: everything that leaves the machine goes through _http_get, one request per call, paced, never retried

_last_request = [None]


def _http_get(url):
    """The single place this tool touches the network. Returns the body bytes; raises ApiError on any failure (no retry)."""
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    try:
        with urllib.request.urlopen(req, timeout=TIMEOUT_S) as resp:  # noqa: S310 (https, fixed host)
            body = resp.read(MAX_BYTES + 1)
    except (urllib.error.URLError, OSError) as e:
        raise ApiError(f"request failed: {e}") from e
    if len(body) > MAX_BYTES:
        raise ApiError("response larger than the size cap")
    return body


def fetch(endpoint, params, clock=None, sleep=None):
    """GET BASE_URL+endpoint with `params` plus tool=, waiting so requests are at least MIN_INTERVAL_S apart."""
    clock, sleep = clock or time.monotonic, sleep or time.sleep
    if "email" in params or "api_key" in params:
        raise ValueError("this tool sends no email address and no API key")
    url = BASE_URL + endpoint + "?" + urllib.parse.urlencode(dict(params, tool=TOOL_NAME))
    if _last_request[0] is not None:
        wait = MIN_INTERVAL_S - (clock() - _last_request[0])
        if wait > 0:
            sleep(wait)
    try:
        return _http_get(url)
    finally:
        _last_request[0] = clock()


def _json(body):
    try:
        return json.loads(body.decode("utf-8"))
    except (UnicodeDecodeError, ValueError) as e:
        raise ApiError("the API did not answer with JSON") from e


def search_ids(gene, change):
    """ClinVar variation ids esearch returns for `GENE[gene] AND CHANGE` (it matches the change anywhere, so callers filter)."""
    data = _json(fetch("esearch.fcgi", {"db": "clinvar", "term": f"{gene}[gene] AND {change}", "retmode": "json", "retmax": ESEARCH_RETMAX}))
    try:
        return list(data["esearchresult"]["idlist"])
    except (KeyError, TypeError) as e:
        raise ApiError("esearch answer has no esearchresult.idlist") from e


def summaries(ids):
    """{id: esummary document} for `ids`, in one request."""
    data = _json(fetch("esummary.fcgi", {"db": "clinvar", "id": ",".join(ids), "retmode": "json"}))
    try:
        result = data["result"]
        return {uid: result[uid] for uid in result["uids"]}
    except (KeyError, TypeError) as e:
        raise ApiError("esummary answer has no result.uids") from e


def vcv_xml(variation_id):
    """The VCV XML record for one variation id. `is_variationid` is needed: without it the same request answered with an empty set."""
    return fetch("efetch.fcgi", {"db": "clinvar", "id": variation_id, "rettype": "vcv", "is_variationid": "true"})


# --- parsing: every value is a {"value", "source"} pair; missing or empty becomes NOT_STATED with the field that was looked for

def stated(value, source):
    if value is None or (isinstance(value, str) and value.strip() in ("", EMPTY_DATE)):
        return {"value": NOT_STATED, "source": source}
    return {"value": value, "source": source}


def _get(doc, path):
    cur = doc
    for key in path.split("."):
        if not isinstance(cur, dict) or key not in cur:
            return None
        cur = cur[key]
    return cur


def matches(summary, gene, change):
    """True when esummary's `protein_change` (comma-separated when there are several) holds `change` and `genes` holds `gene`."""
    changes = [c.strip() for c in str(summary.get("protein_change") or "").split(",")]
    genes = [g.get("symbol") for g in summary.get("genes") or [] if isinstance(g, dict)]
    return change in changes and gene in genes


def profile_from_summary(s):
    """The esummary part of a profile."""
    vs = s.get("variation_set") if isinstance(s.get("variation_set"), list) else []
    scv = _get(s, "supporting_submissions.scv")
    return {
        "variation_id": stated(s.get("uid"), "esummary:uid"),
        "accession": stated(s.get("accession_version"), "esummary:accession_version"),
        "title": stated(s.get("title"), "esummary:title"),
        "hgvs": [stated(v.get("variation_name"), f"esummary:variation_set[{i}].variation_name") for i, v in enumerate(vs)]
        + [stated(v.get("cdna_change"), f"esummary:variation_set[{i}].cdna_change") for i, v in enumerate(vs)]
        or [stated(None, "esummary:variation_set")],
        "protein_change": stated(s.get("protein_change"), "esummary:protein_change"),
        "germline_classification": stated(_get(s, "germline_classification.description"), "esummary:germline_classification.description"),
        "review_status": stated(_get(s, "germline_classification.review_status"), "esummary:germline_classification.review_status"),
        "last_evaluated": stated(_get(s, "germline_classification.last_evaluated"), "esummary:germline_classification.last_evaluated"),
        "conditions": [stated(t.get("trait_name"), f"esummary:germline_classification.trait_set[{i}].trait_name")
                       for i, t in enumerate(_get(s, "germline_classification.trait_set") or []) if isinstance(t, dict)]
        or [stated(None, "esummary:germline_classification.trait_set")],
        "supporting_scv_count": ({"value": len(scv), "source": "computed: len(esummary:supporting_submissions.scv)"}
                                 if isinstance(scv, list) else stated(None, "esummary:supporting_submissions.scv")),
    }


def parse_vcv(body):
    """(number_of_submissions, number_of_submitters, [submitter rows]) from a VCV XML record. Raises ApiError when it is not one."""
    if b"<!DOCTYPE" in body or b"<!ENTITY" in body:
        raise ApiError("VCV answer declares a DOCTYPE or ENTITY; refusing to parse it")
    try:
        root = ET.fromstring(body)
    except ET.ParseError as e:
        raise ApiError("VCV answer is not XML") from e
    va = root.find("VariationArchive")
    if va is None:
        raise ApiError("VCV answer has no VariationArchive (an empty set means the id was not read as a variation id)")
    src = "efetch-vcv:VariationArchive/@"
    rows = []
    for ca in va.iter("ClinicalAssertion"):
        acc = ca.find("ClinVarAccession")
        cl = ca.find("Classification")
        p = "efetch-vcv:ClinicalAssertion/"
        rows.append({  # keys are SUBMITTER_FIELDS
            "scv": stated(acc.get("Accession") if acc is not None else None, p + "ClinVarAccession/@Accession"),
            "submitter": stated(acc.get("SubmitterName") if acc is not None else None, p + "ClinVarAccession/@SubmitterName"),
            "classification": stated(cl.findtext("GermlineClassification") if cl is not None else None, p + "Classification/GermlineClassification"),
            "date_last_evaluated": stated(cl.get("DateLastEvaluated") if cl is not None else None, p + "Classification/@DateLastEvaluated"),
            "review_status": stated(cl.findtext("ReviewStatus") if cl is not None else None, p + "Classification/ReviewStatus"),
        })
    return stated(va.get("NumberOfSubmissions"), src + "NumberOfSubmissions"), stated(va.get("NumberOfSubmitters"), src + "NumberOfSubmitters"), rows


def tally(rows):
    """Counts per classification (letter case ignored, first spelling kept) and whether submitters disagree.

    They disagree when two or more different classifications remain after leaving out NO_CALL and NOT_STATED; a submitter that makes no call
    does not disagree with anyone."""
    counts, label = {}, {}
    for r in rows:
        text = r["classification"]["value"]
        key = text.casefold()
        label.setdefault(key, text)
        counts[key] = counts.get(key, 0) + 1
    calls = {k for k in counts if k not in (NO_CALL, NOT_STATED)}
    return {label[k]: n for k, n in counts.items()}, len(calls) > 1


def build_profile(summary, vcv_body=None):
    out = profile_from_summary(summary)
    if vcv_body is None:
        out["number_of_submissions"] = stated(None, "efetch-vcv not requested (--no-submitters)")
        out["submitters"] = None
        return out
    n_sub, n_submitters, rows = parse_vcv(vcv_body)
    counts, disagree = tally(rows)
    out.update({"number_of_submissions": n_sub, "number_of_submitters": n_submitters, "submitters": rows,
                "classification_counts": {"value": counts, "source": "computed: count of efetch-vcv:ClinicalAssertion/Classification/GermlineClassification, case ignored"},
                "submitters_disagree": {"value": disagree, "source": "computed: more than one classification other than 'not provided'"}})
    return out


def footer(now):
    return FOOTER.format(date=now.strftime("%Y-%m-%d") + " (UTC)")


def render_text(query, profiles, now):
    def line(label, f):
        return f"  {label}: {f['value']}    [{f['source']}]"
    lines = [f"Query: {query['gene']} {query['protein_change']} (searched as {query['searched']})", ""]
    for p in profiles:
        lines += [line("ClinVar variation ID", p["variation_id"]), line("Accession", p["accession"]), line("Title", p["title"])]
        lines += [line("HGVS", h) for h in p["hgvs"]]
        lines += [line("Protein change", p["protein_change"]), line("Germline classification", p["germline_classification"]),
                  line("Review status", p["review_status"]), line("Date last evaluated", p["last_evaluated"])]
        lines += [line("Condition", c) for c in p["conditions"]]
        lines += [line("Supporting SCV accessions", p["supporting_scv_count"]), line("Number of submissions", p["number_of_submissions"])]
        if p["submitters"] is None:
            lines.append("  Per-submitter rows: not requested")
        else:
            lines.append(line("Number of submitters", p["number_of_submitters"]))
            lines.append("  Per-submitter classifications:")
            if not p["submitters"]:
                lines.append(f"    {NOT_STATED}    [efetch-vcv:ClinicalAssertion]")
            for r in p["submitters"]:
                lines.append(f"    {r['scv']['value']} | {r['submitter']['value']} | {r['classification']['value']} | "
                             f"evaluated {r['date_last_evaluated']['value']} | {r['review_status']['value']}")
            if p["submitters"]:
                lines.append("    [fields: " + ", ".join(p["submitters"][0][k]["source"] for k in SUBMITTER_FIELDS) + "]")
            counts = ", ".join(f"{k}: {n}" for k, n in p["classification_counts"]["value"].items()) or NOT_STATED
            lines.append(f"  Classification counts: {counts}    [{p['classification_counts']['source']}]")
            flag = "submitters disagree" if p["submitters_disagree"]["value"] else "submitters do not disagree"
            lines.append(f"  {flag}    [{p['submitters_disagree']['source']}]")
        lines.append("")
    lines.append(footer(now))
    return "\n".join(lines)


def run(gene, change, want_submitters=True):
    """(query, profiles) for one validated variant. Raises Refused or ApiError."""
    gene, one = validate(gene, change)
    ids = search_ids(gene, one)
    sums = summaries(ids) if ids else {}
    chosen = [s for s in sums.values() if isinstance(s, dict) and matches(s, gene, one)]
    profiles = [build_profile(s, vcv_xml(s["uid"]) if want_submitters else None) for s in chosen]
    return {"gene": gene, "protein_change": change, "searched": f"{gene}[gene] AND {one}", "esearch_ids": ids}, profiles


def main(argv=None, out=sys.stdout, err=sys.stderr, now=None):
    ap = argparse.ArgumentParser(description="Print what ClinVar states about one protein variant (research tool, not medical advice).")
    ap.add_argument("gene", help="gene symbol, e.g. CFTR")
    ap.add_argument("protein_change", help="protein change, e.g. R31L, F508del, G542X or Arg31Leu")
    ap.add_argument("--json", action="store_true", help="print JSON")
    ap.add_argument("--no-submitters", action="store_true", help="skip the efetch VCV request (no per-submitter rows)")
    args = ap.parse_args(argv)
    now = now or datetime.datetime.now(datetime.timezone.utc)
    code, query, profiles, error = 0, None, [], None
    try:
        query, profiles = run(args.gene, args.protein_change, not args.no_submitters)
        if not profiles:
            code, error = 1, "no ClinVar record returned by esearch has esummary protein_change equal to the query"
    except Refused as e:
        code, error = 2, f"refused: {e}"
    except ApiError as e:
        code, error = 3, f"API problem, not retried: {e}"
    if args.json:
        doc = {"query": query, "profiles": profiles, "error": error, "retrieved_utc": now.strftime("%Y-%m-%dT%H:%M:%SZ"),
               "catalog_ids": list(CATALOG_IDS), "footer": footer(now)}
        print(json.dumps(doc, indent=2, ensure_ascii=False), file=out)
    else:
        if error:
            print(error, file=err)
        print(render_text(query, profiles, now) if query else footer(now), file=out)
    return code


if __name__ == "__main__":
    sys.exit(main())
