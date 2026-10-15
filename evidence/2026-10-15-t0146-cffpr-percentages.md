---
task_id: T-0146
title: US CFFPR 2024 Modulator Ineligibility Demographics
date: 2026-10-15
answer_tier: T2
review_tier: full
catalog_sources:
  - cffpr-adr-2024
claims_file: docs/cffpr-waiting-claims.json
status: verified
---

# Task Evidence Packet: T-0146

## 1. Curative Mission & Biophysical Relevance
- **Track**: Track A (Curative Biophysics & Target Discovery) / Track B (Clinical Translation & Public Platform)
- **Biophysical Context**: In Cystic Fibrosis, CFTR modulators (elexacaftor/tezacaftor/ivacaftor, ivacaftor, tezacaftor/ivacaftor, lumacaftor/ivacaftor) restore trafficking and gating for approved genotypes. Individuals ineligible for these therapies carry premature termination codons (PTC / nonsense mutations like W1282X, G542X), splice-site defects, or ultra-rare structural variants.
- **Mission Impact**: Quantifying the racial and ethnic composition of the modulator-ineligible population highlights cohorts that do not benefit from current small-molecule modulators, underscoring the biophysical urgency of genetic and alternative conductance therapies (PTC readthrough, NMD suppression, AAV/LNP delivery, TMEM16A/SLC26A9 activation).

## 2. Primary Sources & Exact Locators
- **Source ID**: `cffpr-adr-2024`
- **Full Citation**: Cystic Fibrosis Foundation Patient Registry. 2024 Annual Data Report. Bethesda, Maryland. ©2025 Cystic Fibrosis Foundation.
- **Report Scope**: 120 pages, calendar year 2024.
- **Exact Locators**:
  - **Table**: "Demographics of CFTRm Ineligible Individuals Aged 12 Years and Older, 2024", Page 63.
  - **Baseline Demographic Comparison**: Section "National Registry Demographics", Page 8 (and Page 63 prose).

## 3. Verbatim Excerpts & Quantitative Extractions
> "Data were available for 1,226 adolescents and adults in 2024... with 14.9 percent Black or African American race and 23.7 percent Hispanic ethnicity compared to 3.4 percent and 10.3 percent in the overall CF population, respectively." (CFFPR 2024 Annual Data Report, Page 63)

### Extracted Demographic Metrics:
| Metric Description | Exact Value | Source Page | Claim Identifier |
| :--- | :---: | :---: | :--- |
| Ineligible Cohort: Black or African American | 14.9% | Page 63 | `cffpr-ineligible-black-14-9` |
| Ineligible Cohort: Hispanic | 23.7% | Page 63 | `cffpr-ineligible-hispanic-23-7` |
| Registry Baseline: Black or African American | 3.4% | Page 8, 63 | `cffpr-registry-baseline-black-3-4` |
| Registry Baseline: Hispanic | 10.3% | Page 8, 63 | `cffpr-registry-baseline-hispanic-10-3` |

## 4. Denominators, Exclusions, & Cohort Boundaries ($N$)
- **Exact Denominator ($N$)**: **$N = 1,226$**.
- **Inclusion Criteria**:
  - Individuals with CF aged 12 years and older registered in the CFFPR in 2024.
  - Not eligible for any approved CFTR modulator based on age and CFTR genotype.
  - Documented health insurance encounter data available in 2024.
- **Exclusions**:
  - Children with CF under 12 years of age (excluded from this specific table).
  - Eligible individuals who were not prescribed a modulator due to provider choice, patient preference, or adverse effects ($N = 3,037$, distinct cohort).
  - Individuals with missing encounter insurance data.

## 5. Invariants Against Misuse (The Guardrails)
1. **No Total Registry Equivalence**: Do NOT report 14.9% or 23.7% as percentages of the *total* US CF registry. They represent percentages of the *modulator-ineligible sub-cohort aged 12+ with insurance data* ($N=1,226$).
2. **The Non-Pooling Imperative**: Do NOT sum or average these figures with other national registries (e.g. UKCFR, CCFR, ACFDR, ECFSPR). Different national registries employ distinct encounter filters and demographic categorization schemes.
3. **Temporal Invariant (ALYFTREK Approval Boundary)**: These figures predate FDA approval and commercial availability of ALYFTREK (vanzacaftor/tezacaftor/deutivacaftor, approved late 2024/early 2025). During calendar year 2024, VTD was not an approved commercial therapy in the registry data tables.

## 6. Limitations & Explicit "Not Stated" Register
- The source does **not state** the exact variant breakdown (e.g. W1282X vs G542X vs non-F508del missense) specifically within the 14.9% Black or 23.7% Hispanic subgroups in Table 1.3 on page 63.
- The source does **not state** geographic or socioeconomic distributions specifically cross-tabulated with race/ethnicity in this ineligible cohort.

## 7. Deterministic Verification Command
```powershell
python ../cf-research/tools/claims/check_claims.py docs/cffpr-waiting-claims.json
```
```bash
pytest ../cf-research/tools/claims/test_check_claims.py
```
