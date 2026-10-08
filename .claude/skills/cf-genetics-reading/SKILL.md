---
name: cf-genetics-reading
description: How to read genotype and CFTR-variant content in CF registry reports (CFF, ECFSPR, UK, Australia, Ireland, Canada) and in the public CFTR-France variant lists before extracting or comparing a number - what the base population is, what "unknown" means, how variants are grouped, why eligibility by genotype is not one population, and how the CFTR-France classification scheme changed between 2019 and 2026. Use before extracting, comparing or quoting any genotype, variant, class or eligibility figure.
license: CC0-1.0
compatibility: Works with any Agent-Skills-spec-compatible tool (Claude Code, GitHub Copilot, Qwen3-Coder, Cursor, etc.)
---

# Reading genetics in CF registry reports and variant lists

Use this before taking a genotype, variant, class or "eligible by genotype" number from a registry report or from a CFTR-France list. It records what six registries' reports and ten saved CFTR-France versions showed when read on 2026-10-05, in the lab's own words. It is a reading aid, not a statement about CF genetics. Check `sources/catalog.yaml` for how each source may be obtained and reused: most reports may be read, cited and linked but not reproduced, so this note carries no tables. The CC0 dedication in the header covers this note's wording only; the short quotations belong to the publishers.

Method: a script chose the pages, a model extracted facts with exact quotations, a script kept only the verbatim ones, and a model from a different family checked the final wording against the quotations. See `proposals/2026-10-05-registry-reading-provenance.md`. "Our inference" marks a step the sources do not state. Related: `cffpr-report-reading`, `ecfspr-report-reading`.

Before a claim taken from a report goes into a note, write it with its exact quote in a claims file and run `tools/claims/check_claims.py` (see `tools/claims/README.md`); the answer tiers T0 to T3 proposed in `proposals/2026-10-08-antihallucination-strategy.md` say what a claim must rest on before it is used.

This is a reading aid for lab tools. It does not interpret anyone's genotype, and a class or an eligibility figure here is never a statement about a person.

## First, the unit: variants, alleles, people

- A list of variants (CFTR-France) is not a list of people. A registry table of "the most common variants" counts people with at least one copy (UK, Canada, CFF), so a person with two different variants appears twice and the groups overlap. An allele frequency (ECFSPR, and Australia's per-allele figures) counts chromosomes. These do not convert into each other without the original data.
- The same variant can show two figures for the same country and year. Ireland's G551D appears as 8.3 percent in the ECFSPR 2024 report and as 15.4 percent in the Irish 2024 report. The Irish text does not say which base it uses; our inference is that one is per allele and the other per person. Never put the two side by side as a trend.
- Terms differ. ECFSPR and Australia write "variant" (ECFSPR adds that these are also called mutations), Canada and Ireland write "mutation", and the UK uses both. The reports give different totals: over 1,400 mutations that can cause CF (UK glossary), about 2,000 CFTR mutations identified (Ireland), more than 2,000 (Canada). The CFF text says more than 2,000 in the 2020 to 2022 reports and more than 4,000 in 2023 and 2024. None of them gives a reason for the difference; do not supply one.

## The base population differs

| Report | The genotype figures are based on |
| --- | --- |
| ECFSPR 2024 | People who had DNA analysis, in the table of people with two known variants |
| UK 2024 | Everyone registered who has a known genotype |
| Australia 2025 | Everyone in the registry that year, including lung-transplant recipients, who are kept in this section and left out of the clinical-outcome sections |
| Ireland 2024 | People alive and not lost to follow-up at the end of 2024 |
| Canada 2024 | All people reported on in 2024, including those with no mutation recorded |
| CFF 2020 to 2024 | Individuals genotyped (99.4 percent of the Registry each year); the F508del prevalence footnote says so in the 2021 to 2024 reports |

Two further population notes: the UK's "genotyped" share in its 2024 summary is now calculated from all registered patients, where earlier reports used only people with an annual review, so it may differ from older publications. The UK's charts of variant combinations by devolved nation place people by the location of their CF center, which does not account for people who travel between nations for care.

## How complete is genotyping, and what "unknown" means

- Genotyping is almost universal: 99.4 percent of the CFF Registry in each year 2020 to 2024, more than 99 percent in ECFSPR, 99.5 percent with at least one allele known in Australia, and 98.9 percent with at least one mutation recorded in Canada.
- Completeness is not the same as identification. In ECFSPR 2024, 97 percent of those tested had at least two variants identified, 1.9 percent of variants stayed unidentified, and 3.0 percent of people had at least one unidentified variant. ECFSPR warns that the share of unidentified variants differs greatly between countries partly because testing differs, from kits for a limited number of common variants to whole-gene analysis.
- ECFSPR separates "Not done" (no DNA analysis) from "Unknown" (analysis done, but one or both variants not found). In its genotype grouping it folds "Unknown" into "other", on the reasoning that F508del is in every kit and would have been found.
- Ireland says its "Unknown" may reflect a genetic report missing from the patient file, not a failed test. Canada gives only "at least one mutation recorded" in its text, not a figure for both.
- ECFSPR codes variants against a list of the 1,600 most common ones (from the CFTR1 database), takes others as free text, cleans misspellings, and reports the original ("legacy") name where one exists. It also captures complex alleles (two or more variants on one allele) but does not present them.

## How variants are grouped, and why the groupings are not severity

- CFF: every report from 2020 to 2024 says the five-class system was used in previous reports and calls it an oversimplification, because many variants cause more than one defect in CFTR function. The reports instead list individual variants and compare sweat chloride between genotypes with little or no CFTR function and those with residual function. The 2020 to 2022 reports also mention classifying a variant by its response to a given modulator ("theratyping") as an approach under evaluation. Do not describe the five-class system as the scheme in use in these reports.
- CFF Technical Supplement: describes the five-class system as used in the 2016 report, says a proposed sixth class was not used, and says some variants have no known class. It collapses two alleles into three groups, in this order: class IV-V if at least one allele is class 4 or 5; class I-III if both are class 1 to 3; unknown if one allele is class 1 to 3 and the other has no class, or if neither has one. It notes a modest rise of class IV-V genotypes among children aged 5 and younger, born in the era of universal newborn screening, and puts the larger share at older ages down to survivor bias. It does not say whether these rules apply to a 2020 to 2024 figure.
- ECFSPR: its genotype figure shows the prevalence of F508del and two class I variants by country, with the remaining people as "other", and says plainly that the grouping does not reflect severity. It describes class I as no functional protein and notes that in 2025 the EMA authorised two modulator combinations for people with at least one non-class I variant.
- Canada: a six-class functional classification, with variants of unknown impact left unclassified. It cautions that people with the same class or the same mutation do not all have the same clinical picture.
- UK: sorts people aged 6 and over into genotype groups by F508del status and by lists of responsive variants (FDA lists and a French compassionate-use list). People in group 4 are described as likely not eligible for modulators.
- Do not use "minimal function" or "residual function" as registry genotype groups. Ireland's modulator section uses "minimal function" and "residual" and the UK's uses "residual function", in their descriptions of modulator use, not as groups in a genetics section.

## Eligibility by genotype is not one population

- Each report anchors eligibility to something different: ECFSPR to the EMA criteria for 2024, in all countries but three (Israel, the Russian Federation and Switzerland), which use their national criteria; the UK to 178 FDA-listed variants plus 7 in the French programme, deliberately ignoring the December 2024 expansion to 272 named variants because it came at the very end of the review year; Australia to the criteria as of 30 June 2025; Canada to Health Canada criteria as of 31 December 2024; Ireland to age and genotype without naming a list or a date.
- Canada notes that a July 2024 expansion of 152 further mutations means people newly eligible may not yet have had access, and leaves out of its modulator table the 7.6 percent not eligible by age and genotype.
- CFF counts eligibility by age and genotype. Its 2023 and 2024 reports say 92.3 percent of the Registry is eligible by genotype for at least one modulator; the 2024 report counts people eligible only for the newly approved vanzacaftor combination as ineligible in its charts. Its chapter on people not eligible by genotype is restricted in 2022 to people aged 12 and older with no lung transplant; the 2023 and 2024 chapters count people ineligible by age or genotype and report data for adolescents and adults. All three note that Black and Hispanic individuals are a larger share of that group than of the whole Registry.
- Ireland's 2024 report gives 95.9 percent eligible by age and genotype, with 6.7 percent of children and 2.5 percent of adults ineligible. Australia assesses eligibility using age at the start of the year while reporting age groups at year end.
- So never put one registry's "eligible" share next to another's, and do not compute across them (the lab's rule is that cross-registry arithmetic waits for a registry scientist).

## What the CFF reports say about specific variants

- F508del prevalence is the share of genotyped people with at least one copy, not of everyone with CF. It fell slowly across the five reports, from 85.8 percent in 2020 to 85.0 percent in 2024. The "most common variants" tables count people with one or two copies.
- Layout change: the 2020 to 2022 reports print the 25 most common variants for the whole Registry; the 2023 and 2024 reports print the 10 most common by race and ethnicity and say frequency varies across groups. F508del is lower among Black (58.6 percent in 2023, 57.8 in 2024) and Hispanic (67.9 and 67.7) individuals. 3120+1G>A is the second most common among Black individuals (about 22 percent), and 3876delA the third among Hispanic individuals. Do not attribute these within-group figures to the whole Registry, and do not guess why the table changed.
- R117H needs care. Each report contrasts the share among people genotyped in 1993 (under 1 percent) with those genotyped in the report year and suggests newborn screening algorithms as one possible explanation. About a tenth of people with R117H had a sweat chloride below 30 (10.5 percent in 2020, 11.1 in 2022, 10.1 in 2024), which the 2022 and 2024 reports say adds to diagnostic complexity. More than half of the people with R117H have no recorded poly-T status. From 2024 the Registry merges R117H and any poly-T result (and 5T with any poly-TG result) into single complex alleles, which lowered the stand-alone percentages for R117H and 5T.
- Sweat chloride by genotype uses each person's highest recorded test. In 2020 to 2022 the footnote warns that for some people (132 in 2020, 191 in 2021 and 2022) that value may be post-modulator; from 2023 post-modulator values are excluded. The reports' explanation of the median trend changes from year to year, so quote it per year.

## CFTR-France public variant lists

What they are: public downloads of CFTR-France's list of variants and their classification (the site is at cftr.chu-montpellier.fr). The site's own home page says its variants come from 10 French laboratories with a focus on phenotypes. The workbooks hold variant names, positions, segment, DNA and protein type, class, subclass and a CFTR2 column. They hold no counts of people, no frequencies and no genotypes. Eleven dated versions are listed on the download page, from July 2019 to September 2026; ten were saved and read as counts only, the 06/10/2025 version was not saved, and two files for 03/07/2019 have identical counts.

The scheme changed:

- 2019 to 2024: three classes, which are disease-causing, unclassified and non disease-causing, with subclasses that include VUS1 to VUS5 (the counts do not show which class each sits under). The column is named `classe` in 2019 and `class` from October 2020. The October 2020 file calls the middle class `VUS`; every other file calls it `unclassified`.
- September 2026: five classes. `likely pathogenic` (58) and `likely benign` (38) are new, and the VUS tiers are gone. New subclass labels `likely CF` (22) and `likely CFTR-RD` (28) appear, and `non-CF` grows to 42 from a single variant in each of the 2023 and 2024 files.
- Columns arrived in steps: protein name and protein type in December 2022, DNA type in September 2023, the CFTR2 column in September 2024, and in September 2026 "Other names" plus hg38 coordinates (earlier versions have hg19 only, transcript NM_000492.3 throughout). Header spelling changed too (`type_segment` became `type segment` in 2020), so read by the header of the version in hand.

What follows for anyone counting:

- A count of one label is not comparable across versions until the scheme is mapped. A script counting `unclassified` finds none in the October 2020 file. `unclassified` falls from 504 (September 2024) to 437 (September 2026) while the list grows from 1,120 to 1,214 and 96 variants gain the two new classes; the aggregates cannot show whether variants were resolved or relabelled. Do not describe the fall as progress.
- The VUS3 subclass reads 5, 236, 121, absent and 5 across the October 2021, March 2022, December 2022, September 2023 and September 2024 files, and `undefined` moves the other way. Treat any trend in either as unreliable until a person understands it. A file with the same row count can still differ: the 2019 versions of July and September both have 934 rows, yet their class counts differ by one (our inference: one variant was reclassified).
- The CFTR2 column is YES or NO as exported. Its meaning is not defined on the pages read, and the split flips between 2024 (693 NO, 427 YES) and 2026 (774 YES, 440 NO). Do not call it a CFTR2 classification or compare the two years.
- `unclassified` and `undefined` are not benign. In 2026 the benign-leaning labels are 89 `non disease-causing` and 38 `likely benign`, 127 together; adding the 437 unclassified would make 564, more than four times as many, and that would be wrong.
- A row is not a patient. The download page says 167 of the 1,214 variants have no individual in CFTR-France and 1,047 have at least one. The home page gives a third figure, 22,840 variants (1,047 different) among 5,963 individuals, which is a different measure again.
- The 2026 download page and workbook agree on the total, the likely pathogenic and likely CF and CFTR-RD counts, the 437 unclassified and the 89 non disease-causing, with one gap: the page says 41 VUS non-CF, the workbook has 42 `non-CF`. The 167 and 1,047 split appears only on the page. The page's "89 disease causing (CF or CFTR-RD)" has no workbook label; it equals 592 minus the 503 variants with a CF, CFTR-RD or varying-consequence subclass (our arithmetic).
- Load-time checks that pass for the 2026 workbook: segment, class, subclass and CFTR2 counts each add up to 1,214 rows. A future download that fails them has changed its rules.

How to cite a number from these lists: say that they are variants listed in CFTR-France as of the file date, under the named scheme, name the transcript and genome build, and cite Claustres et al., Human Mutation 2017 (doi 10.1002/humu.23276), which the site asks for. Never write "CF variants". The site's agreement says the data are for educational purposes, to be read in each patient's clinical context, with no warranty of accuracy or completeness; keep that with any reuse, never present a class as a diagnosis, and never touch any individual-level database the site keeps behind a login (the catalog entry records it as out of scope).

## Do not claim

- That a count from a variant list is a number of patients, carriers or alleles, or says how common a variant is anywhere.
- That unclassified, VUS, undefined or an unknown tier means benign.
- That one registry's genotype share, eligibility share or variant frequency is comparable with another's without the base population, the date and the list behind each.
- That a genotype grouping ranks severity.
- That a reason exists for a change (variant totals, table layouts, class labels) when the report gives none.
- Anything about the CFTR-France 06/10/2025 version, which was not saved.
