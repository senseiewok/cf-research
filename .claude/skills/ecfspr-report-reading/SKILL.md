---
name: ecfspr-report-reading
description: How to read an ECFS Patient Registry (ECFSPR) Annual Data Report before extracting numbers - where things are, which population a figure uses, how the layout changed between report years 2020 to 2024, and what the report itself cautions. Use before extracting or quoting any ECFSPR figure.
license: CC0-1.0
compatibility: Works with any Agent-Skills-spec-compatible tool (Claude Code, GitHub Copilot, Qwen3-Coder, Cursor, etc.)
---

# Reading ECFSPR annual reports

Use this before pulling a number out of an ECFSPR Annual Data Report (2020 to 2024). It records what was observed in those five reports on 2026-10-04, not a general rule about registries. Check `sources/catalog.yaml` in this repo for how each report may be obtained.

Licence note: the CC0 dedication in the header covers the wording of this note only. The short quotations are the registry's, attributed to it and kept brief; the reports themselves are not re-hosted. Before adding more text from a report, run `python tools/sources/check_source_overlap.py` on the draft (see `tools/sources/README.md`).

Before a claim taken from a report goes into a note, write it with its exact quote in a claims file and run `tools/claims/check_claims.py` (see `tools/claims/README.md`); the answer tiers T0 to T3 proposed in `proposals/2026-10-08-antihallucination-strategy.md` say what a claim must rest on before it is used.

## Where things are

Each report has a contents page, a "Data report" with numbered sections, then appendices. Numbers and titles move between years, so read the contents of the report you are using instead of assuming a section number.

- Sections 7 and 8 are "Respiratory complications and therapies" and "Gastro-intestinal complications and therapies" in 2020 and 2021, and "Complications" and "Therapies" from 2022.
- Late sections differ. 2021 numbers "12. Longitudinal data", 2022 numbers "12. Data Quality", and 2023 and 2024 number "13. Data quality". The 2020 and 2021 contents also list "Data Quality" and "Data quality" as unnumbered entries after the numbered sections, so look for an unnumbered entry as well as a numbered one.
- The 2023 and 2024 contents add "10. Pregnancy". Transplantation is "10. Transplantation" in 2020 to 2022 and "11. Lung and other organ transplants" in 2023 and 2024.
- The opening table is titled "Summary of data report" in 2020 to 2023; in 2024 it is an overview of 2024 set beside the earlier years 2014 and 2019.
- The contents list an appendix of variables and definitions in 2020 to 2023 (Appendix 3 in 2020 and 2021, Appendix 4 in 2022 and 2023). The 2024 contents list none.

## Page numbering

In the 2020 report the summary heading is printed as page 7 in the contents but is on PDF page 11. In 2021 to 2024 the printed and PDF page of that heading match. This was checked only at that heading, so say which numbering you cite and re-check before relying on an offset elsewhere.

## Which population a figure uses

The summary table has three different footnotes, and each applies to different rows. Read the footnote attached to the row you are quoting.

- Each report from 2020 to 2024 has three. The first says who was seen during the year, the second who was alive at the end of the year, and the third again who was seen during the year, followed by a United Kingdom note and a total of its own.
- 2021: the footnotes limit the rows to people seen by clinical staff during the year, people alive on 31 December 2021, and again people seen during the year.
- 2020: the same three limits, worded in terms of patients (seen during the year for the first and third, alive on 31 December 2020 for the second).
- The United Kingdom note says that all individuals with a confirmed diagnosis of CF were included in 2021 to 2024; the 2020 wording says all patients.
- In 2024 the footnotes give totals for three years (2024, 2019 and 2014); the first states that total, and the second limits the row to people alive at the end of each year.

Copy the printed denominator verbatim beside every value. Never apply one row's footnote to its neighbour.

## Definitions and cross-references

- All five introductions use chronic Pseudomonas aeruginosa infection as an example of a definition the registry applies, and define it by the modified Leeds criteria, an elevated anti-Pseudomonas antibody level, or both (2024 wording).
- They point to the full list of variables and definitions in different places: Appendix 3 on page 152 (2020 and 2021), Appendix 4 (2022 and 2023), and the registry's website (2024).
- The 2021 pointer to page 152 does not match the 2021 contents, which put Appendix 3 on printed page 171. Check a cross-reference against that report's contents before trusting it.
- Every report says it used common reference populations for percent-predicted values; the 2023 text adds that another reference population would give slightly different values.

## What the reports caution

- 2020: the report advises comparing only countries with similar coverage.
- 2020: countries with fewer than five patients in an age group are left out of both the graphs and the tables.
- 2024: the tables and graphs were made from the database as it stood on 20 November 2025.

## Text quality and comparing years

- The 2022 PDF prints some letter pairs as odd characters, for example "Gene/cs" and "Nutri/on" in its contents and "popula?ons" in its text. Copy such text as printed, normalise it with a deterministic rule, and have a human confirm the label. Never let a model repair it from memory.
- Do not say a definition or a layout changed between years from a partial quote. Read the full sentence in each report. A truncated quote can make two identical definitions look different.
- Do not guess why a section moved or was renamed.
