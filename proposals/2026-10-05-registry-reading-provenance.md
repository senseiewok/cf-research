# Provenance of the registry reading skills (2026-10-05)

Status: record, written 2026-10-05, for `cffpr-report-reading` and `cf-genetics-reading` in this repo. It says which model read what, what a script verified, where the models disagreed and what was kept. The two skills are reading aids written in the lab's own words; they contain no tables or charts from the reports.

## Sources read

Sixteen registry PDFs saved by hand into the git-ignored `sources/downloads/` (CFF annual reports 2020 to 2024, the Technical Supplement, five Highlights handouts, the UK 2024, Australia 2025, Ireland 2024, Canada 2024 and ECFSPR 2024 reports), and eleven CFTR-France variant workbooks (ten distinct versions) with the text of the site's home and download pages. Text came out with `pdftotext` page by page. None of it is in a repo.

## The pipeline

| Step | Who | Result |
| --- | --- | --- |
| Pick the pages and make plain text | A script (`extract.py`, outside the repos) | Page files per report; CFTR-France counts per version from a script that reads the workbooks and prints aggregates only |
| Extract facts, each with exact quotations, doubts and a list of what not to claim | Four read-only Fable reviewers, one per area | 50 (CFF structure), 45 (CFF genetics), 42 (other registries) and 23 (CFTR-France) items |
| Keep only the verbatim ones | A script (`verify_facts.py`): whitespace and quote marks normalised, quotation at least 25 characters, found on the cited page, a claim about a difference needs two reports | 42, 44, 35 and 17 kept on the sets as I entered them (the last two as shortened subsets of the reviewers' lists; the others dropped for a quotation not on the page or too short) |
| Check each claim against its quotations | The local worker (Qwen3.8 27B), a different model family, one packet per fact | On the two sets whose packets held the full sentences: 25 of 42 and 35 of 44 supported, the rest "overstated" or "unsupported". The other two sets were reviewed from shorthand labels and the verdicts were discarded (ledger E72) |
| Write the skills | The controlling agent, using only claims the quotations support, narrowed where the reviewer disagreed | Drafts |
| Check the final wording against all the evidence | A fresh Fable reader per skill, each given one file, and a local Qwen pass | 14 points on the CFF skill and 23 on the genetics skill from Fable; the Qwen result is recorded below |
| Resolve each point | The controlling agent, against the page text | See below |

## What the final review found

On the CFF skill, 13 of the 14 points were year ranges or hedges stretched beyond the quotations, and all were narrowed (for example "through 2023" became "in the 2023 summary table"). Two were checked on the pages and kept: the Foundation's figure-permissions notice and the Supplement's copyright year come from the catalog entry read on the title and credits pages, and the two-page offset between PDF and printed numbers was confirmed on the page footers.

On the genetics skill, 11 of the 23 points were rated high. Four were real errors in the draft and were corrected: a ratio (127 to 564 is 4.4 times, not "roughly five"), a claimed page and workbook agreement on the 1,047 and 167 split (it is on the page only), the CFF genotype base, and the UK being credited with "minimal function" wording. Seven were quotations trimmed from the evidence file, which were confirmed on the pages (the Canadian July 2024 expansion, the ECFSPR figure of F508del and two class I variants, the 2025 EMA authorisation, the three countries with national criteria, the UK age 6 header, free-text variant handling, and the use of "allele frequency"). The remaining points were narrowed or labelled "our inference".

## Limits

- The reviewers were told to use read-only local tools; the transcripts show no web tool and two directory listings in a shell (ledger E66).
- Quotations were verified on the extracted text, not on the rendered PDF, and chart text was refused by rule; only sentences and table rows were used.
- A different model family judged the claims against the quotations, but the blind readers of the final text were Fable, the same family as the extractors. The local pass is the independent check; its result and any changes it caused are noted under "Local pass" below.
- Nothing here has been read by a registry scientist (board T-0011). Cross-registry comparisons are described as limits, and no figure from one registry is set against another's.
- The skills record what the reports said on 2026-10-05; they are not a statement about CF genetics and must not be used to interpret a person's genotype.

## Local pass and the cloud fallback

- CFF skill: the local worker in thinking mode returned 3 low points (the Supplement's copyright year, the page offset, and the mental health guideline wording). The first two had been checked already; the third was narrowed.
- Genetics skill: whole-note thinking runs stopped at the output limit twice and a fast run returned no problems, so under the maintainer's rule to fall back to cloud, an Opus reader (a different model from the Fable extractors and readers) checked the corrected text. It confirmed every number and calculation and raised 1 high and 9 low wording points, all resolved: the INSERM and CHU attribution came from the catalog and not the pages and was reworded; the chapter on people ineligible by genotype is described as restricted to ages 12 and over without a transplant only for 2022; the "non-CF" subclass is described as growing, not new; a claim of what the CFF pages do not list was removed.
- The same note split in two was then run locally in thinking mode (89 and 73 seconds, one low point each, both already settled). The working method is recorded in `model-qwen3-8-27b`.
