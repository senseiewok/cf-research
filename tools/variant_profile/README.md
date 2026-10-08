# tools/variant_profile

Prints what ClinVar states about one protein variant in one gene: the variation ID, HGVS names, germline classification, review status, date last evaluated, the number of submissions, and each submitter's classification and date, with a count per classification and a "submitters disagree" flag. Every value is printed next to the API field it came from. A value the API does not give prints as `not stated`.

Python 3.9+, standard library only. No API key, no `.env`.

## Run

```bash
python variant_profile.py CFTR R31L
python variant_profile.py CFTR Arg31Leu --json
python variant_profile.py CFTR F508del --no-submitters   # esummary only, no per-submitter rows
```

Exit 0: printed one or more profiles; read `Match status` (`single` or `ambiguous (N records match)`). 1: no record matches. 2: refused input. 3: the API failed or answered unexpectedly (never retried).

## How it reads ClinVar

Three or four requests to the documented NCBI E-utilities API, nothing else:

1. `esearch` on the HGVS protein name in the variant-name field: `CFTR[gene] AND "p.Phe508del"[varname]` (stops are written `Ter`, so G542X is searched as `p.Gly542Ter`). Only when that finds nothing, one second `esearch` with the plain term `CFTR[gene] AND F508del`, which matches any field. `retmax` is 100; when esearch reports more records than it returned, the output says the list may be incomplete.
2. `esummary` for the ids found. A record matches when its `genes` holds the gene and either `protein_change` holds the change (ClinVar writes stops as `*`, as in `W1282*`; the tool treats `*` and `X` as the same) or its `title` or a `variation_name` holds `(p.Phe508del)`. The title route is needed because the main F508del record (7105) has an empty `protein_change`. The output names the field each match was made on; `--json` also lists the esearch ids that did not match (`esearch_ids_not_matched`).
3. One `efetch rettype=vcv is_variationid=true` for all matches together, for the per-submitter rows (skipped above 10 matches, and said so).

When more than one record matches (two nucleotide changes with the same protein change, such as N1303K's c.3909C>G and c.3909C>A, or a haplotype record that contains the variant) every match is printed with its record type and the status says `ambiguous`. The tool never picks one.

ClinVar's aggregate classification is printed as stated even when it is not one of the five pathogenicity terms (D1152H's aggregate is `drug response`); a `Classification scale` line labels it, and it is never mapped onto the five.

Submitter classifications are counted ignoring letter case (one submitter writes "Uncertain Significance"); "not provided" is counted but is not a disagreement. Any other differing term, `drug response` included, sets the "submitters disagree" flag; that is deliberately conservative.

## Limits

From catalog entries `ncbi-eutils` and `clinvar` in [`sources/catalog.yaml`](../../sources/catalog.yaml):

- at most 1 request per second (`MAX_RPS`), below the catalog's 2.5 and NCBI's stated 3 without a key;
- one request per call, no retry; at most four requests per run (two esearch terms, esummary, efetch); a 30-second timeout and a 20 MB response cap;
- a descriptive User-Agent and the `tool` parameter; no email address or other personal identifier is sent (the catalog note mentions NCBI's email parameter; this tool leaves it out on purpose);
- only `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/`; it never crawls web pages (robots.txt there disallows crawlers).

## What it does not do

- It takes a variant, never a person. Input must be a gene symbol (`^[A-Z0-9-]{2,12}$`) and a protein change in one-letter form (`R31L`, `F508del`, `G542X`, `R553*`) or three-letter form (`Arg31Leu`); anything else, such as a name, a date of birth or a record number, is refused with exit 2.
- It does not interpret a classification, combine sources, or say anything about any person's genotype. CFTR2 remains the CF-specific authority.
- It does not read clinical-impact or oncogenicity classifications, allele frequencies, or citations.

## Source and rights

ClinVar via NCBI E-utilities (catalog ids `clinvar`, `ncbi-eutils`). ClinVar asks to be credited as the source, and its data are not for direct diagnostic use or medical decisions without review by a genetics professional. NCBI places no restrictions of its own, but submitters may hold rights in their submissions. Every output ends with the footer: "Research tool, not medical advice. A database entry is not an interpretation of any person's genotype; that belongs with their clinical genetics team. Source: ClinVar via NCBI E-utilities, retrieved <UTC date>."

## Tests

```bash
python -m unittest discover -s tools/variant_profile -v   # from the repo root; no network
```

Tested against the CFTR R31L fixture recorded 2026-10-07, and fixtures for F508del, R117H, N1303K, W1282X, Y1092X, D1152H and R104fs recorded 2026-10-08 (see [`fixtures/README.md`](fixtures/README.md)).
