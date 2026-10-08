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
