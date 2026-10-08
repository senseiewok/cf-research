# Fixtures

Answers from NCBI E-utilities, retrieved 2026-10-07 (UTC). The three small files are saved byte for byte; the large VCV record is trimmed (see below). Each request sent the `tool=senseiewok-variant-profile` parameter and the User-Agent `senseiewok-variant-profile/0.1 (+https://github.com/senseiewok/cf-research)`, no email address and no API key, about one second apart. Source: ClinVar via NCBI E-utilities.

| File | Request | HTTP status |
| --- | --- | --- |
| `esearch_CFTR_R31L.json` | `esearch.fcgi?db=clinvar&term=CFTR[gene] AND R31L&retmode=json&retmax=20` | 200 (3 ids: 54087, 53653, 35893) |
| `esummary_54087_53653_35893.json` | `esummary.fcgi?db=clinvar&id=54087,53653,35893&retmode=json` | 200 |
| `efetch_vcv_54087.xml` | `efetch.fcgi?db=clinvar&id=54087&rettype=vcv` | 200, but an empty `<set/>`: without `is_variationid` the id was not read as a variation id |
| `efetch_vcv_54087_is_variationid.xml` | `efetch.fcgi?db=clinvar&id=54087&rettype=vcv&is_variationid=true` | 200 (VCV000054087 version 29, 13 submissions). **Trimmed after saving** to the record's attributes and, per submission, the accession, submitter name, review status, classification and date last evaluated, so the fixture holds no free text, comments or contact details |

The esearch term matches `R31L` in any field, so it also returned R31C (35893) and R104fs (53653); the tool keeps only records whose esummary `protein_change` equals the query.

These are ClinVar data as of the retrieval date and will drift; the tests check the parser, not today's classification. ClinVar asks to be credited as the source; submitters may hold their own rights in their submissions.

## Follow-up fixtures, retrieved 2026-10-08 (UTC)

Recorded after a 10-variant audit found six variants the first version could not resolve. Same User-Agent and `tool` parameter, no email address, no API key, at least one second apart. Five live requests, all HTTP 200.

| File | Request | Kept |
| --- | --- | --- |
| `esearch_varname_CFTR_F508del.json` | `esearch.fcgi?db=clinvar&term=CFTR[gene] AND "p.Phe508del"[varname]&retmode=json&retmax=100` | byte for byte (3 ids: 4072070, 634837, 7105) |
| `esearch_varname_CFTR_six_variants.json` | the same shape with `("p.Arg117His"[varname] OR "p.Trp1282Ter"[varname] OR "p.Tyr1092Ter"[varname] OR "p.Asn1303Lys"[varname] OR "p.Asp1152His"[varname] OR "p.Arg31Leu"[varname])` | byte for byte (12 ids). One OR query to save requests; tests feed it to every single-variant query, so the tool's filter has to split it |
| `esearch_varname_CFTR_R104fs.json` | the same shape with `"p.Arg104fs"[varname]` | byte for byte (1 id: 53653) |
| `esummary_15_ids_trimmed.json` | `esummary.fcgi?db=clinvar&id=<the 15 ids above>&retmode=json` | **trimmed** to uid, obj_type, accession, title, variation names, supporting submissions, germline classification (description, dates, review status, trait names), gene symbols and protein_change; locations and allele frequencies dropped |
| `efetch_vcv_35867_trimmed.xml` | `efetch.fcgi?db=clinvar&id=35867&rettype=vcv&is_variationid=true` | **trimmed** like the R31L record: record attributes, and per submission the accession, submitter name, review status, classification and date last evaluated. No free text, comments or contact details. All 41 submitter names are organisations |

What they showed: ClinVar writes stops as `*` in `protein_change` (`W1282*`, `Y1092*`); the main F508del record 7105 (101 supporting SCVs) has an empty `protein_change`; R117H and W1282X each have a haplotype record carrying the same protein change; N1303K has two nucleotide changes (7136, c.3909C>G; 4818550, c.3909C>A); D1152H's aggregate germline description is `drug response`.
