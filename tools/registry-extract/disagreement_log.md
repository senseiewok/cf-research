# ECFSPR Extraction Pipeline Disagreement Log

**Task Anchor**: `T-0009` ("ECFSPR extraction pipeline v0.1 (5 report years, ~15 indicators)")  
**Audit Protocol**: 100-cell random sample across 5 report years (2020, 2021, 2022, 2023, 2024)  
**Acceptance Criterion**: ">=99% extraction accuracy on a 100-cell random audit (Wilson 95% CI reported); 100% rows carry url+page; disagreement log committed."  
**Audit Outcome**:
- **Sample Size**: 100 discrete data cells
- **Matches**: 100 / 100
- **Accuracy**: 100.0% (threshold $\ge 99.0\%$)
- **Wilson 95% Score CI**: `[0.9630, 1.0000]` ($96.30\% - 100.00\%$)
- **URL & Page Coverage**: 100.0% (185 / 185 rows carry valid canonical `url` and `page`)

---

## 1. Audit Methodology & Cell Sampling Protocol

The audit engine ([`audit_accuracy.py`](file:///tools/registry-extract/audit_accuracy.py)) extracted all discrete quantitative and categorical cells from [`ecfspr_extracted_dataset.json`](file:///tools/registry-extract/ecfspr_extracted_dataset.json). Each cell represents an atomic value (`count_n`, `percentage`, `value_point`, `value_low`, `value_high`, `median_iqr.median`, `median_iqr.q1`, `median_iqr.q3`, or `denominator_n`).

A pseudo-random number generator initialized with deterministic seed `42` selected 100 cells without replacement. Each cell was verified against:
1. The exact character sequence in `raw_text_quote`.
2. The primary source PDF document in `sources/downloads/` matching `source_sha256`.
3. The printed page number and official download URL from `sources/manifest.json`.

All 100 audited cells matched the ground truth without divergence. The detailed record of every audited cell is committed in [`audit_100_cells.json`](file:///tools/registry-extract/audit_100_cells.json).

---

## 2. Documented Discrepancies, Ambiguities & Edge Cases

During the implementation and verification of the 5-year extraction pipeline, several formatting shifts, typographic inconsistencies, and reporting differences were identified and resolved.

### Case 1: Parametric vs Non-Parametric Reporting Shift (2020 vs 2021–2024)
- **Observation**:
  - In the 2020 report (page 7 / PDF index 10), continuous clinical variables (`age_at_followup`, `age_at_death`, and `age_at_diagnosis`) were reported as separate parametric (mean) and non-parametric (median) rows. For example, age at follow-up reported `mean 21.5, 22.1, 21.8` and `median 18.8, 19.8, 19.4`. Quartiles ($Q_1$ and $Q_3$) were omitted from the summary table.
  - Starting in 2021 and continuing through 2024, the ECFSPR standardized to reporting non-parametric medians with interquartile ranges: `median (25th pctl-75th pctl)` (e.g., `19.4 (10.2-31.8)` in 2022).
- **Resolution**:
  - For 2020, `value_point` captures the scalar median (and dedicated rows capture the mean as `age_at_followup_mean` / `age_at_death_mean`), while `median_iqr` is set to `null` because quartiles were not stated in the source table.
  - For 2021–2024, `median_iqr` is populated with `{"median": ..., "q1": ..., "q3": ...}`, and `value_low` / `value_high` capture the quartile bounds.

### Case 2: Inferred vs Verbatim Denominators (`denominator_n`)
- **Observation**:
  - In 2020, `Patients >= 18 years (patients alive on 31/12/2020)` reported `12788 (51.9%)` for females and `27568 (53.1%)` for total. An initial draft populated `denominator_n` with `24799` and `52246` (the registered counts). However, the cell quote only stated `"27568 (53.1%)"`, and the denominator for patients alive on 31/12 was actually 48,795 (footnote **).
- **Resolution**:
  - In strict compliance with `cf-research-context` Answer Tier T2 invariants, `denominator_n` is set to `null` unless the denominator is explicitly stated in the cell or footnote. Inferred denominators without verbatim cell grounding are prohibited.

### Case 3: Character Ligatures and Special Typographic Symbols in PDFs
- **Observation**:
  - ECFSPR PDFs contain custom OpenType glyphs:
    - Footnote markers: 2020, 2021, and 2023 use private-use character `\uf02a` [U+F02A] or asterisk operator `\u2217` (`*`), whereas 2022 uses standard ASCII `*`.
    - Mathematical operator: `\u2265` (`≥`) is used for ">= 18 years".
    - Ligatures: `\ufb01` (`fi`) in words like "Fibrosis" in PDF stream text.
- **Resolution**:
  - Text normalization handles ligatures and Unicode operators gracefully. Quotes preserve verbatim character streams from the PDF extractor while normalized indicator IDs use standard ASCII `pwcf_adult_ge18`.

### Case 4: 2024 Overview Table Layout Reorganization
- **Observation**:
  - In the 2024 report (page 9 / PDF index 8), Table 1 was restructured from a sex-stratified single-year summary into `Overview for 2024 and previous years (2014 and 2019)`, presenting the total cohort for 2024 alongside historical reference columns for 2019 and 2014.
- **Resolution**:
  - The pipeline extracts the primary 2024 total cohort indicators for all 9 summary outcomes, plus the published historical benchmarks (2019 and 2014) as dedicated comparative rows with appropriate notes.

### Case 5: UK Registry Inclusion in Mortality Analysis
- **Observation**:
  - All 5 reports include footnote `***` explaining that for mortality evaluations, the United Kingdom contributed all individuals with a confirmed diagnosis of CF rather than only those seen during the year:
    - 2020: UK $N = 10,837$; total mortality cohort $N = 50,026$.
    - 2021: UK $N = 10,907$; total mortality cohort $N = 51,901$.
    - 2022: UK $N = 11,148$; total mortality cohort $N = 51,919$.
    - 2023: UK $N = 11,318$; total mortality cohort $N = 53,926$.
    - 2024: UK $N = 11,381$; total mortality cohort $N = 55,199$.
- **Resolution**:
  - The pipeline extracts both `pwcf_mortality_cohort` and `pwcf_uk_included` as explicit standardized indicators, ensuring complete epidemiological transparency.

### Case 6: 5-Tier Denominator Cascade (`denominator_scope`) & Country Stratification (`stratum_country`)
- **Observation**:
  - Across a single ECFSPR report year (e.g., 2024), five distinct population denominators govern different analyses:
    1. `all_registered`: `57,288` (Table 1 & Table 1.1 total registered PwCF)
    2. `mortality_cohort_uk_all`: `55,199` (Table 1 footnote `***` — mortality denominator including all `11,381` UK confirmed diagnoses)
    3. `seen_during_year`: `54,241` (Table 1 footnote `*` & Table 1.1 total seen during 2024)
    4. `alive_end_of_year`: `53,988` (Table 1 footnote `**` — denominator for adult $\ge 18$ proportion `60%`)
    5. `seen_never_transplanted`: `51,294` (Page 13 text — clinical outcome denominator from Chapter 4 onwards)
  - Conflating these 5 tiers is a classic source of epidemiological claim misinterpretation.
- **Resolution**:
  - Added explicit `denominator_scope` and `stratum_country` schema dimensions (`row_schema.json`), extracting all 5 denominator tiers and 16 European country cohorts from Table 1.1 (page 12) with verbatim quotes.

---

## 3. Disagreement Status Summary

| Audit Item | Expected Value | Extracted Value | Disagreement? | Disposition |
| :--- | :--- | :--- | :---: | :--- |
| Random 100-cell sample | Exact PDF values | Exact PDF values | **None (0/100)** | Verified 100.0% accurate |
| URL & Page coverage | 100% rows | 185 / 185 rows | **None** | Full compliance |
| Schema validation | 0 errors | 0 errors | **None** | Passes `row_schema.json` |
