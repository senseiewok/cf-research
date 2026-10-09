# API candidates for the source catalog: terms read 2026-10-07, re-read 2026-10-08 and 2026-10-09

Status: proposed. Six candidate entries were added to `sources/catalog.yaml`. All six are `access: manual`, and none has a `base_url` or `url`, so no tool can reach any of them until a human admits one. The decisions are the human's. The notes below are the lab's reading of the terms as quoted, not legal advice.

Provenance: read by `claude-opus-5-5` on 2026-10-07. Pages were fetched directly with curl, using the user agent `SenseiEwokCFLab/0.1 (reading terms pages for our source catalog)`, at about one request per second (10 seconds on Zenodo, per its Crawl-delay). That came to 27 requests, robots.txt first on each host. A blind review then re-read some pages on 2026-10-08 (DailyMed download terms, Ensembl robots.txt, the AlphaMissense v3 record) and a council re-read the terms pages on 2026-10-09; those requests are not in the 27. No data API was called, no key was requested and no account was made. Details and quotes are in each entry's `notes`.

## Decisions asked of the maintainer

1. DailyMed: accept NLM's download conditions and admit it for label text? Copyright of company-submitted label text is not settled, so the lab would read and link, not copy.
2. Semantic Scholar: accept Ai2's licence, which binds the organisation, carries indemnity and Washington law, lets Ai2 name the organisation in marketing, and sends query text to Ai2? If not, leave it manual.
3. GTEx, Ensembl, Unpaywall: read their terms in a browser, or leave them manual.
4. AlphaMissense: confirm v3 only (CC BY 4.0), cite the Science paper, and download by hand only if a file is needed.

## What each entry says

| Candidate (catalog id) | Catalog access now | Maintainer must do | What the terms read say about a public, non-commercial repo | Recommendation |
| --- | --- | --- | --- | --- |
| GTEx Portal API (`gtex-portal-api`) | manual | Read the GTEx terms and API docs in a browser; robots.txt is `Disallow: /` for all but Googlebot and Bingbot | Not read: robots.txt disallows our agent on the whole host | Keep manual. Link GTEx for CFTR expression by hand |
| DailyMed web services (`dailymed-web-services`) | manual | Nothing to sign up for (none stated); accept NLM's download conditions (acknowledge NLM, no implied endorsement, flag currency if redistributed) | NLM's conditions are stated for NLM data. Copyright of company-submitted label text is not established; NLM's web policies say reproducing copyrighted items beyond fair use needs the holders' written permission | Decision 1. If admitted, read label text and link it; keep approvals on `openfda-label`, because DailyMed labels may differ from FDA-approved labeling |
| Ensembl REST incl. VEP (`ensembl-rest`) | manual | Read the current Ensembl disclaimer and legal terms (the URL redirects to a different host, which we did not follow) | Not stated: terms not read | Our inference, not checked: it may suit `api` once the terms are read; the published rate-limit page is a 2014 example |
| AlphaMissense predictions (`alphamissense-predictions`) | manual | Nothing to sign up for (none stated on the record page); pick the version deliberately | v3 states CC BY 4.0 (attribution, cite the Science paper). v1 states CC BY-NC-SA 4.0 (non-commercial only, share-alike). Zenodo's terms and the v3 README were not read | Use v3 only; download by hand only if a file is needed (v1 files total 6.1 GB; v3 sizes not checked) |
| Semantic Scholar Academic Graph API (`semantic-scholar-api`) | manual | Decision 2. Optional: request a key by form (sent by email) | The licence read grants use to access and display data through the API; data licences vary (CC BY-NC, ODC-BY named) and were not read | Admit keyless only if citation-graph lookups are needed and you accept the licence; our inference, not checked: Europe PMC and OpenAlex may already cover it |
| Unpaywall API (`unpaywall-api`) | manual | Read the API docs in a browser (pages need JavaScript) | Not read: the pages need JavaScript | Keep manual; read the terms and any email requirement first |

## What admission takes

Both DailyMed and Semantic Scholar already carry `terms_url`, so admission would be a small edit: set `access: api`, add a `base_url`, and choose a `max_rps`. Adding `base_url` alone lets the evidence tool's gate pass but triggers no call, because cf-skills has no provider or command for these sources yet. Writing one is a separate, reviewed change. A `url` with a `filename` would let `fetch_sources.py` download an entry, so none is added here.

## What was not covered

- Not read: the GTEx terms, the Ensembl legal terms, the Unpaywall documentation and terms, Zenodo's own terms of use, the AlphaMissense v1 README text, the DailyMed About page on 2026-10-09, and any AlphaMissense hosting outside Zenodo.
- DailyMed: no copyright statement for company-submitted label text was found. NLM's web policies say content contributed by companies may be protected by copyright and that NLM cannot guarantee the copyright status of any item.
- `sources/manual-downloads.md` is generated from the catalog and was regenerated in this change.
- The public cf-skills catalog copy will differ from this catalog once this merges, and no check enforces the sync. Run `sync_skill_catalog.py` after merge.

## Dependencies and network calls

No dependency was added. Admitting a source would add a new external network call only when a cf-skills provider for it exists, which is a separate, reviewed change.
