# CF Registry Extract Row Schema (`cf-registry-extract`)

**Status**: Frozen  
**Task Anchor**: `T-0008` (Freeze the cf-registry-extract row schema including `value_low`/`value_high`, `median_iqr`, and `source_sha256`)  
**Version**: 1.0.0  
**Date**: 2026-10-16  
**License**: CC0-1.0  
**Epistemic Standard**: Answer Tier T2 minimum per `cf-research-context`  

---

## 1. Overview & Architectural Rationale

This document specifies and freezes the standardized record row schema for epidemiological data extracted from Cystic Fibrosis patient registry annual reports (ECFSPR, CFFPR, UKCFR, ACFDR, CCFR, and CFRI).

### Why Freeze Now?
Per project scoping memo and task `T-0008`, everything downstream—automated table extractors, local-worker verification prompts, comparison harnesses, and golden regression fixtures—encodes this record schema. Freezing the contract now:
1. **Prevents Multiplied Downstream Cost**: Changes to field names, value shapes, or provenance keys after building extraction tools invalidate all existing golden corpora and prompt templates.
2. **Standardizes Parametric & Non-Parametric Metrics**: Accommodates both discrete categorical indicators (e.g., patient counts `count_n` and proportions `percentage`) and continuous, non-normally distributed metrics (e.g., `median_iqr` with lower/upper quartile bounds `value_low` / `value_high`).
3. **Enforces Cryptographic Lineage**: Embeds the exact `source_sha256` hash of the underlying source PDF, ensuring extracted data can always be deterministically audited against the original document file in `sources/downloads/`.
4. **Preserves Verbatim Quote Primacy**: Every row records `raw_text_quote` capturing the exact character sequence from the table cell, enabling automated quote matching under `check_claims.py` and `check_numbers.py`.

---

## 2. Schema Specification

The formal JSON Schema is defined in [`row_schema.json`](file:///tools/registry-extract/row_schema.json). Every extraction row represents an atomic stratum observation from a report table.

### 2.1 Field Reference

| Field Name | Type | Required? | Description & Constraints |
| :--- | :--- | :---: | :--- |
| `row_id` | `string` | **Yes** | Deterministic row identifier (e.g. `ECFSPR-2022-P08-R01`). Pattern: `^[A-Z0-9_-]+$`. |
| `registry_id` | `string` | **Yes** | Canonical registry key matching `sources/catalog.yaml`. Enum: `["ecfspr", "cffpr", "ukcfr", "acfdr", "ccfr", "cfri"]`. |
| `report_title` | `string` | **Yes** | Official report publication title. |
| `report_year` | `integer` | **Yes** | Reporting calendar/cohort year (e.g., `2022`). Range: 1990–2035. |
| `page_number` | `integer` | **Yes** | 1-indexed printed page number in the report PDF. |
| `page` | `integer` | No | Alias for `page_number` (1-indexed page in PDF). |
| `url` | `string` | No | Canonical download or landing URL matching `sources/manifest.json`. |
| `table_id` | `string` | **Yes** | Identifier or title of the table (e.g. `Table 1.1`, `Summary of data report`). |
| `indicator_id` | `string` | **Yes** | Normalized snake_case indicator key (e.g. `pwcf_registered`, `age_at_followup`, `pwcf_adult_ge18`). Pattern: `^[a-z0-9_]+$`. |
| `indicator_name` | `string` | **Yes** | Full verbatim title of the indicator from the table row header. |
| `stratum_sex` | `string` | **Yes** | Biological sex stratum. Enum: `["females", "males", "total", "not_applicable"]`. |
| `count_n` | `integer` \| `null` | No | Absolute patient count ($N$), or `null` if the indicator reports only continuous metrics. |
| `percentage` | `number` \| `null` | No | Observed percentage ($0.0 \le \% \le 100.0$), or `null` if not applicable. |
| `denominator_n` | `integer` \| `null` | No | Denominator cohort count for percentages, if reported or footnoted. |
| `value_point` | `number` \| `null` | No | Point estimate (mean, rate, median) when scalar analysis is needed. |
| `value_low` | `number` \| `null` | No | Lower bound of the interval (e.g. 25th percentile, minimum, or 95% CI lower). |
| `value_high` | `number` \| `null` | No | Upper bound of the interval (e.g. 75th percentile, maximum, or 95% CI upper). |
| `median_iqr` | `object` \| `null` | No | Structured non-parametric distribution containing `{"median": float, "q1": float, "q3": float}`. |
| `unit` | `string` | **Yes** | Measurement unit. Enum: `["count_and_percentage", "count", "percentage", "years", "months", "mmol_l", "score", "ratio", "other"]`. |
| `source_sha256` | `string` | **Yes** | 64-character lowercase hex SHA-256 digest of the source PDF. Pattern: `^[0-9a-f]{64}$`. |
| `raw_text_quote` | `string` | **Yes** | Verbatim text string from the cell (e.g. `"25964 (47.6)"`, `"19.4 (10.2-31.8)"`). |
| `extraction_kind` | `string` | **Yes** | Epistemic kind per `cf-research-context`: `"observed"`, `"computed"`, or `"inferred"`. |
| `notes` | `string` \| `null` | No | Associated report footnotes, methodological caveats, or sample restrictions. |

---

## 3. Representation Rules for Value Fields

1. **Counts & Percentages**:
   - For cells formatted as `n (%)` (e.g., `"25964 (47.6)"`):
     - `count_n`: `25964`
     - `percentage`: `47.6`
     - `unit`: `"count_and_percentage"`
     - `raw_text_quote`: `"25964 (47.6)"`
2. **Medians and Interquartile Ranges (IQR)**:
   - For cells formatted as `median (25th pctl-75th pctl)` (e.g., `"19.4 (10.2-31.8)"`):
     - `value_point`: `19.4`
     - `value_low`: `10.2`
     - `value_high`: `31.8`
     - `median_iqr`: `{"median": 19.4, "q1": 10.2, "q3": 31.8}`
     - `unit`: `"years"` (or `"months"`)
     - `raw_text_quote`: `"19.4 (10.2-31.8)"`
3. **Footnote Invariant**:
   - Subpopulation restrictions (e.g., "seen during the year", "alive on 31/12/2022") must be recorded verbatim in `notes` and reflected in `denominator_n` where explicitly stated.

---

## 4. Gold Sample Verification (ECFSPR 2022 p.8)

A hand-extracted reference dataset comprising **27 rows** across all 9 summary indicators and 3 sex strata (females, males, total) from page 8 of the *European Cystic Fibrosis Society Patient Registry 2022 Annual Data Report* is provided in [`gold_sample_ecfspr_2022_p8.json`](file:///tools/registry-extract/gold_sample_ecfspr_2022_p8.json).

- **Source Document**: `sources/downloads/ECFSPR_Annual_Report_2022_v1.0.pdf`
- **Source SHA-256**: `22868219bcaeb5e4b2090643b0d5b1e8be91fddb1cf37ba74bf1bfaeec11506e`
- **Validation**: Verified against `row_schema.json` via `tools/registry-extract/test_registry_extract.py`.
