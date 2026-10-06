# Registry data access

What we can and cannot retrieve programmatically from the CF registry publishers, tested rather than assumed. Relevant to any project that reads published aggregate CF statistics.

Checked 2026-10-04 from a datacentre IP with a default Python HTTP client. Results may differ from a residential IP or a browser.

## What responded

| Host | Status | Evidence | Label |
| --- | --- | --- | --- |
| `pr.ecfs.eu` | 200. Annual-reports index lists 32 report PDFs | https://pr.ecfs.eu/annual-reports/ (2026-10-04) | `[verified]` |
| `www.ecfs.eu` | 200, `application/pdf` | 2021 annual report PDF fetched directly (2026-10-04) | `[verified]` |
| `www.cff.org`, `cff.org` | **403 Access Denied** from the CDN edge, on the site root, on `/media/<id>/download` report paths, **and on `/robots.txt`** | https://www.cff.org/ (2026-10-04) | `[verified]` |
| `databases.lovd.nl` | **402 Payment Required**, with an explicit notice that scraping is not permitted and permission must be requested via their contact form | https://databases.lovd.nl/shared/genes/CFTR (2026-10-04) | `[verified]` |
| `cftr2.org` | 200 | https://cftr2.org/ (2026-10-04) | `[verified]` |
| `eutils.ncbi.nlm.nih.gov`, `pmc.ncbi.nlm.nih.gov` | 200 | E-utilities esearch, PMC article page (2026-10-04) | `[verified]` |
| `clinicaltrials.gov` API v2 | 200, `application/json` | `/api/v2/studies?query.cond=cystic+fibrosis` (2026-10-04) | `[verified]` |

## Crawl policy

- ECFSPR publishes `User-agent: * / Disallow:` — an empty disallow list, which permits automated retrieval of the whole site. Source: https://pr.ecfs.eu/robots.txt (2026-10-04). `[verified]`
- CFF's `robots.txt` is itself behind the 403, so **we have no statement of their crawl policy and no basis for inferring permission**. The only signal available is refusal. `[verified]`
- LOVD states the position in the response body: no scraping, contact them for access. `[verified]`

## Position on browser automation

**We do not use Playwright, Selenium, a spoofed user-agent, a mirror, a cache, or an archive to retrieve CFF material.** `[hypothesis — this is our policy, not an external fact]`

The reasoning, so a stranger can check it:

1. A 403 served to an automated client is the publisher declining automated access. Driving a headless browser so the same request presents as a human one evades that decision. The tool is irrelevant; the intent is what is being judged.
2. CFF's `robots.txt` is unreachable, so we cannot even claim to be acting within a published policy. Acting anyway means asserting a permission nobody granted.
3. A person downloading a public report in their own browser is a person reading a public page. That is not the same act as a script populating a pipeline at machine speed, and the difference is not a technicality we get to argue away afterwards.
4. The cost of being wrong is asymmetric. A volunteer project with no institutional backing trades entirely on being trustworthy. Being correct about a terms-of-service question is worth more to us than a fortnight of calendar time.

It is plausible that the 403 is a datacentre-IP reputation rule rather than a deliberate stance against research reuse — CFF publishes these reports for the community, and blocking researchers would be an odd goal. `[unverified]` The way to resolve that is to ask them, not to find out which client gets through.

## Sanctioned routes for CFF material

1. A human downloads the annual data reports in a browser and places them in a local folder; pipelines read local files and record a SHA-256 per input so a third party can reproduce the extraction from the same bytes.
2. Ask the CF Foundation directly whether tabular reuse and programmatic retrieval are acceptable, and what citation string they want. If they say yes, automate freely and cite this note as the date the question was asked.

## Email draft: permissions request

Send to the CF Foundation registry/communications contact, and an equivalent to ECFSPR. Adjust before sending; this is a draft, not a sent message.

> Subject: Permission question — reuse of published Patient Registry Annual Data Report figures
>
> Hello,
>
> I run a small independent, unfunded, open-source project that is building a machine-readable table of aggregate statistics that are already published in CF registry annual data reports, so that each figure can be traced back to the exact report and page it came from. No patient-level data is involved and none is being requested. The project has no institutional affiliation and does not represent itself as endorsed by any CF organisation.
>
> Three questions:
>
> 1. Is it acceptable to extract numeric values from the published annual data reports and redistribute them in a table that cites the report, year, page, and URL for every value?
> 2. What citation or attribution string would you like used?
> 3. Automated requests to www.cff.org are refused at the CDN (403), so the reports are currently being downloaded by hand. Is there an endpoint or arrangement you would prefer for programmatic retrieval, or would you rather we continue downloading manually?
>
> If the answer to (1) is no, we will publish only the extraction script and no values, so that users run it against their own copies of the reports.
>
> Thank you for the registry and for publishing the reports openly.

## Where to download the reports

**No report PDFs are committed to this repo.** They are third-party copyrighted publications; we link to them, record a SHA-256 of each file we used, and leave distribution to the publishers. Put downloaded files in a local folder outside the repo, or in an ignored path.

Suggested local folder: `~/Documents/cf-registry-sources/` (Windows: `%USERPROFILE%\Documents\cf-registry-sources\`).

Keep the publisher's original filename. The manifest records the URL, the download date, and the SHA-256, so the filename is a convenience rather than the identifier.

### ECFSPR, fetchable automatically, no hand-downloading needed

Their robots.txt permits automated retrieval, so a script can pull these. URLs confirmed from the annual-reports index 2026-10-04. `[verified]`

| Data year | File |
| --- | --- |
| 2024 | https://pr.ecfs.eu/wp-content/uploads/2026/05/Annual-Report_2024_vs1.0_ECFSPR_20260518.pdf |
| 2023 | https://pr.ecfs.eu/wp-content/uploads/2025/08/Annual-Report_2023_vs1.1_ECFSPR_20250715.pdf |
| 2022 | https://pr.ecfs.eu/wp-content/uploads/2025/06/Annual-Report_2022_vs1.0_ECFSPR_20250130.pdf |
| 2021 | https://pr.ecfs.eu/wp-content/uploads/2025/06/Annual-Report_2021_09Jun2023_ECFSPR_final.pdf |
| 2020 | https://pr.ecfs.eu/wp-content/uploads/2025/06/ECFSPR_Report_2020_v1.0-07Jun2022_website.pdf |

Companion short reports, same index: 2024 Highlights at `/wp-content/uploads/2025/12/ECFS-Patient-Registry-2024-Highlights-Report-v1.0-ECFSPR-20251218.pdf` and 2023 Highlights at `/wp-content/uploads/2025/06/Highlights_Report_2023_vs1.1_ECFSPR_20250508.pdf`. The index lists 32 PDFs going back to 2003, so earlier years are available if the series needs extending. `[verified]`

### CFF, download by hand in a browser

Automated requests are refused (see above), so these must be fetched manually.

**Start here:** https://www.cff.org/medical-professionals/patient-registry carries a `Download (PDF)` link for every year. As of the 2026-10-04 search index it listed the 2025 Highlights Report, the 2024 and 2023 Annual Data Reports with their Highlights, the 2022 Annual Data Report and Highlights handout, the 2021 and 2020 Annual Data Reports, and the Annual Data Report Technical Supplement. `[verified]`

Grab, at minimum:

- The five most recent **Annual Data Reports**, the long ones published each autumn.
- The **Annual Data Report Technical Supplement**, which explains inclusion criteria, limitations, and changes in reporting methods over time. It is the single most important file for the `comparability_notes` field; do not skip it.

Two direct URLs seen in the search index on 2026-10-04, offered as a shortcut rather than a substitute for the landing page: the 2022 Annual Data Report at https://www.cff.org/media/31216/download and the Technical Supplement at https://www.cff.org/media/24916/download. `[verified, URL and document title confirmed from the search index, not opened from this machine]`

One older file is mislabelled in a way that will bite an automated pipeline: a PDF served at a path containing `2019-Patient-Registry-Annual-Data-Report.pdf` carries the title *Patient Registry Annual Data Report 2020*. **Trust the title page, never the filename, when assigning a data year.** `[verified, filename/title mismatch observed in the search index, 2026-10-04]`

### Citation and contacts

- CFF's own suggested citation form, from the Technical Supplement: *Cystic Fibrosis Foundation Patient Registry. Annual Data Report. Bethesda, Maryland. (c)[year] Cystic Fibrosis Foundation.* `[verified, Technical Supplement, 2026-10-04]`
- Registry questions go to the CFF Patient Registry help desk. Individual-level data requests, out of scope for us, go to its separate data-requests desk. The addresses are omitted from this note. `[verified, Technical Supplement and a published study data-availability statement, 2026-10-04]`
- Send the permissions email above to the registry contact, not to a general inbox.

## Open questions

- Whether any openly maintained cross-year, cross-registry table of published CF indicators already exists. Searched 2026-10-04 without finding one; absence of evidence only. `[unverified]`
- Whether ECFSPR's open crawl policy extends to redistributing extracted values, which robots.txt does not speak to. Ask in the same email. `[unverified]`
