# API candidates for the source catalog: terms read 2026-10-07

Status: proposed. Six candidate entries were added to `sources/catalog.yaml` without a `base_url` or `url`, so no tool can reach them until a human admits one. The decisions are the human's.
Provenance: read by `claude-opus-5-5` on 2026-10-07. Pages were fetched directly with curl, using the user agent `SenseiEwokCFLab/0.1 (reading terms pages for our source catalog)`, at about one request per second (10 seconds on Zenodo, per its Crawl-delay). That came to 27 requests, robots.txt first on each host. No data API was called, no key was requested and no account was made. Details and quotes are in each entry's `notes`.

| Candidate (catalog id) | access set | Maintainer must do | Public non-commercial repo? | Recommendation |
| --- | --- | --- | --- | --- |
| GTEx Portal API (`gtex-portal-api`) | manual | Read the GTEx terms and API docs in a browser; robots.txt is `Disallow: /` for all but Googlebot and Bingbot | Could not read | Keep manual. Link GTEx for CFTR expression by hand; reconsider only after a human reads the terms |
| DailyMed web services (`dailymed-web-services`) | api | Nothing to sign up for (none stated); accept NLM's download conditions (acknowledge NLM, no implied endorsement, flag currency if redistributed) | Yes, per the NLM terms read | Admit for label text (add `base_url`) if label sections are needed; keep approvals on `openfda-label`, because DailyMed labels may differ from FDA-approved labeling |
| Ensembl REST incl. VEP (`ensembl-rest`) | manual | Read the current Ensembl disclaimer and legal terms (the URL redirects to a different host, which we did not follow) | Not stated: terms not read | Strong candidate: re-read the terms, then move to `api` with `max_rps: 1.0`; the published rate-limit page is a 2014 example |
| AlphaMissense predictions (`alphamissense-predictions`) | manual | None to download from Zenodo; pick the version deliberately | v3 (CC BY 4.0): yes with attribution and citation. v1 (CC BY-NC-SA 4.0): non-commercial only, share-alike | Use v3 only, cite the Science paper, download by hand only if a file is needed (files are GB-sized) |
| Semantic Scholar Academic Graph API (`semantic-scholar-api`) | api | Optional: request a key by form (sent by email). Using the API at all accepts Ai2's licence for the organisation (indemnity, Washington law) | Yes, with "Semantic Scholar" attribution; data licences vary (CC BY-NC, ODC-BY named) | Admit keyless at 1 request/second if citation-graph lookups are needed and the maintainer accepts the licence; our inference, not checked: Europe PMC and OpenAlex may already cover it |
| Unpaywall API (`unpaywall-api`) | manual | Read the API docs in a browser (pages need JavaScript) | Could not read | Keep manual; a human reads the terms and any email requirement first |

## What was not covered

- Not read: the GTEx terms, the Ensembl legal terms, the Unpaywall documentation and terms, Zenodo's own terms of use, and any AlphaMissense hosting outside Zenodo.
- DailyMed: no copyright statement for company-submitted label text was found; NLM's copyright page covers only government works.
- `sources/manual-downloads.md` is generated from the catalog and now lacks the four new manual entries. Regenerate it with `tools/sources/make_manual_downloads.py`, and sync the cf-skills catalog copy with `sync_skill_catalog.py`, when these entries are accepted.

## Dependencies and network calls

No dependency was added. Admitting DailyMed or Semantic Scholar (adding `base_url`) would add a new external network call from `cf-evidence-loop`, which is a separate, reviewed change.
