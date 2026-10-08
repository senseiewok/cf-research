---
name: cffpr-report-reading
description: How to read a Cystic Fibrosis Foundation Patient Registry (CFFPR) Annual Data Report before extracting numbers - who is counted, how a value is made, which definitions and methods changed between report years 2020 and 2024, and what the reports themselves warn about. Use before extracting, comparing or quoting any CFFPR figure.
license: CC0-1.0
compatibility: Works with any Agent-Skills-spec-compatible tool (Claude Code, GitHub Copilot, Qwen3-Coder, Cursor, etc.)
---

# Reading CFFPR annual reports

Use this before pulling a number out of a CFFPR Annual Data Report (2020 to 2024), its Technical Supplement or a Highlights handout. It records what was read in those reports on 2026-10-05, in the lab's own words. It is not a general rule about registries. Check `sources/catalog.yaml` for how each report may be obtained and reused: the Foundation asks readers to contact the Registry team before reusing its charts or data, so this lab cites and links and does not reproduce charts or tables. The CC0 dedication in the header covers this note's wording only; the short quotations are the Foundation's.

How it was made: a script chose the pages, a model extracted facts with exact quotations, a script kept only the facts whose quotations are verbatim on the cited page, a second model from a different family checked that each claim says no more than its quotations, and the controlling agent kept the narrower claim where they disagreed. See `proposals/2026-10-05-registry-reading-provenance.md`. Anything marked "our inference" is not from the report.

Before a claim taken from a report goes into a note, write it with its exact quote in a claims file and run `tools/claims/check_claims.py` (see `tools/claims/README.md`); the answer tiers T0 to T3 proposed in `proposals/2026-10-08-antihallucination-strategy.md` say what a claim must rest on before it is used.

## Who is counted

- A report covers people with a CF diagnosis who consented to the Registry and were seen at a CF Care Center during the calendar year, including people born, diagnosed or who died that year (all five reports; the Supplement words it as seen at a center or born, diagnosed or died in the year).
- People who have had a lung transplant appear only in the chapters on demographics, diagnosis, the CFTR gene, transplantation and survival. By the reports' own account every other chapter excludes them. Read the chapter's own note before using a denominator.
- People recorded as CRMS/CFSPID or with a CFTR-related disorder are left out of every figure except the one on new diagnoses in the reporting year. Whether someone is entered as CF or as CRMS is a clinical judgement, and non-CF entries are lost to follow-up far more often than people with CF.
- Lost to follow-up means seen in the previous reporting year, not in this one, and not known to have died. Our inference: a death that was never reported is therefore counted as lost.
- Transplant recipients are often followed at transplant centers, so some drop out of the Registry. The 2023 and 2024 mortality chapters, which split deaths by transplant status, note that missing data may affect the findings for this reason.
- New diagnoses for the latest year are not final: centers have until February to enter the previous year's data, and the reports expect more diagnoses to be added later.

## How a value is made

- Encounter data are rolled up to one value per person per year. The 2020 report says continuous measures are averaged and yes/no indicators become ever/never; the 2021 and 2022 reports add that annual lung function, weight and height are the mean of each quarter's highest value.
- Data are reprocessed every year. Compare a new figure with the adjusted history in the same report, not with a number printed in an older one (Supplement).
- Center-level plots show each center's own percentage or median, so they look tighter than a patient-level plot of the same measure (Supplement).
- The Supplement reports an audit of Registry records against medical records: 96.5 percent of clinic visits and 89.7 percent of hospitalizations were found, with similar methods used in later years. It was written around the 2016 report and carries a 2018 copyright, so it may not describe later changes.

## What changed between report years

| Topic | What the reports show |
| --- | --- |
| Lung-function equations | The 2022 report used the 2012 GLI equations. From the 2023 report every percent-predicted value, historical ones included, is recalculated with the 2022 GLI equations; the 2024 report calls them race-neutral. Our inference: do not compare a 2023 or 2024 value with a figure printed in an earlier report. |
| Effect of that change | The 2023 report says the change raised the population median because most participants are white, and that the population pattern may not reflect the change for an individual. The impairment cut-points (below 40 severe, 40 to 69 moderate, 70 to 89 mild, 90 and above normal) were not changed. |
| Pulmonary exacerbation days | Mean days of treatment in the 2023 summary table; medians in 2024. |
| Median age at diagnosis | Months in the 2022 summary table, days in 2024. |
| Adult BMI row | Ages 20 to 40 in the 2021 summary table, ages 20 and older in 2022. |
| Race | Reports before 2021 let a person count in more than one race category; the 2021 and 2024 summary tables footnote this. |
| Visits row | The 2022 summary table has a "phone, phone with video or other" row; the 2023 table has a "telehealth visits" row and a footnote that outpatient visits are encounters with Clinic as the location. |
| Conditional survival | The figure stops at attained age 40 in the 2023 report and extends to age 50 in 2024. |
| Diabetes screening | In 2024 the records include people who used continuous glucose monitoring, in every diabetes report of that edition. |
| Pregnancy | The 2021 report counts women pregnant in the year; the 2024 chart also shows pregnancies reported the year before, because a pregnancy can span two calendar years. Do not read the 2024 counts as distinct pregnancies. |
| New content | 2022 introduced CFTR modulator status and advanced lung disease as new categories; 2024 introduced a chapter on people aged 40 and older and a pregnancy chart. |
| Transplant in complications chapters | Described as removed in 2020 and as censored at the year of transplant in 2024. |

## Definitions the reports state

- A pulmonary exacerbation is a course of IV antibiotics, in hospital or at home. Days-of-treatment measures use only people with at least one.
- A therapy percentage counts anyone on it at any clinical visit in the year, after removing people noted as intolerant or allergic.
- Microbiology denominators differ by organism: S. aureus among people with a bacterial culture, mycobacteria among people with a mycobacterial culture. MSSA and MRSA are not mutually exclusive, so their percentages do not add up.
- Advanced lung disease is flagged by FEV1 below 40 percent, supplemental oxygen or other clinical signs, and stays flagged if lung function improves. The 2024 wording adds a pre-modulator criterion that 2023 does not have.
- Modulator eligibility is judged by age and CFTR genotype. The yearly modulator chart counts only each person's latest medication form and omits people on a modulator who are not eligible. In 2024 people eligible only for the newly approved VTD (Alyftrek) are counted as ineligible unless the text says otherwise.
- The "eligible but not prescribed" group is people aged 12 and older who were eligible in 2020 and had no modulator prescription in a window of years. The window was 2018 to 2022 in the 2022 report and 2022 to 2024 in the 2024 report.
- Median predicted survival is the age that half of today's newborns are expected to exceed if current age-specific death rates hold; it is shown in five-year birth cohorts because single years are unstable. Median conditional survival is for people who have already reached an age. Median age at death describes that year's deaths only and must not be used to predict survival or to compare periods. The annual mortality rate is deaths in the calendar year divided by people in the Registry that year.
- The guidelines behind the mental health screening figures recommend annual screening from age 12, and the 2020 report says caregiver screening is not in the Registry.

## What the reports warn about

- Since 2020, the reports ask readers to take prevalence and incidence in light of less in-person care. They say fewer cultures per person may lower the chance of finding an organism, which would make microbiology after 2019 partly a sampling effect, and that fewer spirometry measurements, some from home devices, and fewer height measurements in children may affect the precision of annual FEV1.
- Telehealth could be recorded as a location only from 28 May 2020. Earlier telehealth visits stayed as "other", so telehealth in March to June 2020 is understated.
- By-age figures mix survival bias at older ages (older people in the Registry have survived and are likely healthier) with a different diagnostic mix among very young people in recent years. The 2024 report puts a favourable lipid pattern in older adults down to survivor bias, not improvement with age.
- The Supplement says median predicted survival is a population figure that does not account for genotype and is not for decisions about one person.

## Pages and text

- In the 2023 and 2024 reports the PDF page runs two ahead of the printed page (file page 22 prints 20; file page 60 prints 58). Earlier years were not checked. Say which numbering you cite.
- Chart text extracts with numbers apart from their labels. Take a figure from a sentence or a table with its row label, never from a chart's loose numbers.

## Genetics

Genotype and variant content (how variants are grouped, the retired class system, who is eligible by genotype) is in `cf-genetics-reading`.
