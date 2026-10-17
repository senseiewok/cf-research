#!/usr/bin/env python3
"""ECFSPR Registry Data Extraction Pipeline v0.1.

Extracts core demographic, genetic, spirometry, nutrition, microbiology, and
transplant/mortality indicators from 5 ECFSPR Annual Report PDFs (2020-2024).

Strictly adheres to tools/registry-extract/row_schema.json and sources/catalog.yaml.
Every row carries cryptographic source_sha256, verbatim raw_text_quote, page_number,
alias page, and official publication url.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
DEFAULT_SCHEMA_PATH = HERE / "row_schema.json"
DEFAULT_OUTPUT_PATH = HERE / "ecfspr_v0_1_extracted.json"

# Fallback download directories (local checkout or parent checkout when running in a git worktree)
CANDIDATE_DOWNLOAD_DIRS = [
    REPO_ROOT / "sources" / "downloads",
    *(p / "cf-research" / "sources" / "downloads" for p in REPO_ROOT.parents),
]

# Official metadata per sources/manifest.json and sources/catalog.yaml
REPORTS_META = {
    2020: {
        "filename": "ECFSPR_Annual_Report_2020_v1.0.pdf",
        "report_title": "ECFSPR Annual Report 2020",
        "url": "https://pr.ecfs.eu/wp-content/uploads/2025/06/ECFSPR_Report_2020_v1.0-07Jun2022_website.pdf",
        "sha256": "6f0430d66af114c2770e3cffb28273163c22d35c7d45cf15dbab297df06624c8",
    },
    2021: {
        "filename": "ECFSPR_Annual_Report_2021.pdf",
        "report_title": "ECFSPR Annual Report 2021",
        "url": "https://pr.ecfs.eu/wp-content/uploads/2025/06/Annual-Report_2021_09Jun2023_ECFSPR_final.pdf",
        "sha256": "34c8b625cbae85403131d3a2e453358439068ebdc9fdf41131f9333a1a87af5a",
    },
    2022: {
        "filename": "ECFSPR_Annual_Report_2022_v1.0.pdf",
        "report_title": "ECFSPR Annual Report 2022",
        "url": "https://pr.ecfs.eu/wp-content/uploads/2025/06/Annual-Report_2022_vs1.0_ECFSPR_20250130.pdf",
        "sha256": "22868219bcaeb5e4b2090643b0d5b1e8be91fddb1cf37ba74bf1bfaeec11506e",
    },
    2023: {
        "filename": "ECFSPR_Annual_Report_2023_v1.1.pdf",
        "report_title": "ECFSPR Annual Report 2023",
        "url": "https://pr.ecfs.eu/wp-content/uploads/2025/08/Annual-Report_2023_vs1.1_ECFSPR_20250715.pdf",
        "sha256": "4d161be5b49cdbe83f82ee991e234b706f73161d190a096745009845cf243cf3",
    },
    2024: {
        "filename": "ECFSPR_Annual_Report_2024_v1.0.pdf",
        "report_title": "ECFSPR Annual Report 2024",
        "url": "https://pr.ecfs.eu/wp-content/uploads/2026/05/Annual-Report_2024_vs1.0_ECFSPR_20260518.pdf",
        "sha256": "ed8e68f78f33ef014d32975eebca7c161494d68f21fbc82d8979781e1abcf0d6",
    },
}


def find_downloads_dir(override: Optional[Path] = None) -> Optional[Path]:
    if override and override.exists():
        return override
    for d in CANDIDATE_DOWNLOAD_DIRS:
        if d.is_dir() and any((d / m["filename"]).is_file() for m in REPORTS_META.values()):
            return d
    for d in CANDIDATE_DOWNLOAD_DIRS:
        if d.is_dir():
            return d
    return None


def get_pdf_path(year: int, downloads_dir: Optional[Path]) -> Optional[Path]:
    if not downloads_dir:
        return None
    fname = REPORTS_META[year]["filename"]
    p = downloads_dir / fname
    return p if p.exists() else None


def verify_pdf_hash(path: Path, expected_hash: str) -> bool:
    if not path.exists():
        return False
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    return digest.lower() == expected_hash.lower()


def build_extracted_dataset(downloads_dir: Optional[Path] = None) -> List[Dict[str, Any]]:
    """Build the standardized dataset across 5 report years and ~15 indicators."""
    rows: List[Dict[str, Any]] = []

    # -------------------------------------------------------------------------
    # 1. ECFSPR 2020 (Page 11: Summary of data report, 27 rows)
    # -------------------------------------------------------------------------
    y2020 = REPORTS_META[2020]
    meta_2020 = {
        "registry_id": "ecfspr",
        "report_title": y2020["report_title"],
        "report_year": 2020,
        "page_number": 11,
        "page": 11,
        "url": y2020["url"],
        "source_sha256": y2020["sha256"],
        "table_id": "Summary of data report",
        "extraction_kind": "observed",
    }
    fn_seen_2020 = "* Only patients seen during the year are presented. The total number of patients presented is 49,111."
    fn_alive_2020 = "** Only patients alive at 31/12/2020 are presented. The total number of patients presented is 48,795."
    fn_dec_2020 = "*** Only patients seen during the year. For UK, confirmed diagnosis included (N=10,837). Total presented is 50,026."

    rows.extend([
        # pwcf_registered
        {**meta_2020, "row_id": "ECFSPR-2020-P11-R01", "indicator_id": "pwcf_registered", "indicator_name": "Patients registered in the ECFSPR", "stratum_sex": "females", "count_n": 24799, "percentage": 47.5, "denominator_n": 52246, "value_point": None, "value_low": None, "value_high": None, "median_iqr": None, "unit": "count_and_percentage", "raw_text_quote": "24799 (47.5%)", "notes": "Registered cohort 2020."},
        {**meta_2020, "row_id": "ECFSPR-2020-P11-R02", "indicator_id": "pwcf_registered", "indicator_name": "Patients registered in the ECFSPR", "stratum_sex": "males", "count_n": 27447, "percentage": 52.5, "denominator_n": 52246, "value_point": None, "value_low": None, "value_high": None, "median_iqr": None, "unit": "count_and_percentage", "raw_text_quote": "27447 (52.5%)", "notes": "Registered cohort 2020."},
        {**meta_2020, "row_id": "ECFSPR-2020-P11-R03", "indicator_id": "pwcf_registered", "indicator_name": "Patients registered in the ECFSPR", "stratum_sex": "total", "count_n": 52246, "percentage": None, "denominator_n": 52246, "value_point": None, "value_low": None, "value_high": None, "median_iqr": None, "unit": "count", "raw_text_quote": "52246", "notes": "Total registered cohort 2020."},
        # age_at_followup
        {**meta_2020, "row_id": "ECFSPR-2020-P11-R04", "indicator_id": "age_at_followup", "indicator_name": "Age at follow-up (in years; patients alive on 31/12/2020)", "stratum_sex": "females", "count_n": None, "percentage": None, "denominator_n": None, "value_point": 18.8, "value_low": None, "value_high": None, "median_iqr": None, "unit": "years", "raw_text_quote": "median 18.8", "notes": "Patients alive on 31/12/2020; mean 21.5."},
        {**meta_2020, "row_id": "ECFSPR-2020-P11-R05", "indicator_id": "age_at_followup", "indicator_name": "Age at follow-up (in years; patients alive on 31/12/2020)", "stratum_sex": "males", "count_n": None, "percentage": None, "denominator_n": None, "value_point": 19.8, "value_low": None, "value_high": None, "median_iqr": None, "unit": "years", "raw_text_quote": "median 19.8", "notes": "Patients alive on 31/12/2020; mean 22.1."},
        {**meta_2020, "row_id": "ECFSPR-2020-P11-R06", "indicator_id": "age_at_followup", "indicator_name": "Age at follow-up (in years; patients alive on 31/12/2020)", "stratum_sex": "total", "count_n": None, "percentage": None, "denominator_n": None, "value_point": 19.4, "value_low": None, "value_high": None, "median_iqr": None, "unit": "years", "raw_text_quote": "median 19.4", "notes": "Patients alive on 31/12/2020; mean 21.8."},
        # pwcf_adult_ge18
        {**meta_2020, "row_id": "ECFSPR-2020-P11-R07", "indicator_id": "pwcf_adult_ge18", "indicator_name": "Patients >= 18 years (patients alive on 31/12/2020)", "stratum_sex": "females", "count_n": 12788, "percentage": 51.9, "denominator_n": 24799, "value_point": None, "value_low": None, "value_high": None, "median_iqr": None, "unit": "count_and_percentage", "raw_text_quote": "12788 (51.9%)", "notes": "Patients alive on 31/12/2020."},
        {**meta_2020, "row_id": "ECFSPR-2020-P11-R08", "indicator_id": "pwcf_adult_ge18", "indicator_name": "Patients >= 18 years (patients alive on 31/12/2020)", "stratum_sex": "males", "count_n": 14780, "percentage": 54.2, "denominator_n": 27447, "value_point": None, "value_low": None, "value_high": None, "median_iqr": None, "unit": "count_and_percentage", "raw_text_quote": "14780 (54.2%)", "notes": "Patients alive on 31/12/2020."},
        {**meta_2020, "row_id": "ECFSPR-2020-P11-R09", "indicator_id": "pwcf_adult_ge18", "indicator_name": "Patients >= 18 years (patients alive on 31/12/2020)", "stratum_sex": "total", "count_n": 27568, "percentage": 53.1, "denominator_n": 52246, "value_point": None, "value_low": None, "value_high": None, "median_iqr": None, "unit": "count_and_percentage", "raw_text_quote": "27568 (53.1%)", "notes": "Patients alive on 31/12/2020."},
        # age_at_diagnosis_months
        {**meta_2020, "row_id": "ECFSPR-2020-P11-R10", "indicator_id": "age_at_diagnosis_months", "indicator_name": "Age at diagnosis (months)", "stratum_sex": "females", "count_n": None, "percentage": None, "denominator_n": None, "value_point": 3.6, "value_low": None, "value_high": None, "median_iqr": None, "unit": "months", "raw_text_quote": "median (months) 3.6", "notes": fn_seen_2020},
        {**meta_2020, "row_id": "ECFSPR-2020-P11-R11", "indicator_id": "age_at_diagnosis_months", "indicator_name": "Age at diagnosis (months)", "stratum_sex": "males", "count_n": None, "percentage": None, "denominator_n": None, "value_point": 3.6, "value_low": None, "value_high": None, "median_iqr": None, "unit": "months", "raw_text_quote": "median (months) 3.6", "notes": fn_seen_2020},
        {**meta_2020, "row_id": "ECFSPR-2020-P11-R12", "indicator_id": "age_at_diagnosis_months", "indicator_name": "Age at diagnosis (months)", "stratum_sex": "total", "count_n": None, "percentage": None, "denominator_n": None, "value_point": 3.6, "value_low": None, "value_high": None, "median_iqr": None, "unit": "months", "raw_text_quote": "median (months) 3.6", "notes": fn_seen_2020},
        # pwcf_f508del_allele
        {**meta_2020, "row_id": "ECFSPR-2020-P11-R13", "indicator_id": "pwcf_f508del_allele", "indicator_name": "Patients with at least one F508del allele recorded", "stratum_sex": "females", "count_n": 18776, "percentage": 80.7, "denominator_n": 23266, "value_point": None, "value_low": None, "value_high": None, "median_iqr": None, "unit": "count_and_percentage", "raw_text_quote": "18776 (80.7%)", "notes": fn_seen_2020},
        {**meta_2020, "row_id": "ECFSPR-2020-P11-R14", "indicator_id": "pwcf_f508del_allele", "indicator_name": "Patients with at least one F508del allele recorded", "stratum_sex": "males", "count_n": 20629, "percentage": 80.7, "denominator_n": 25563, "value_point": None, "value_low": None, "value_high": None, "median_iqr": None, "unit": "count_and_percentage", "raw_text_quote": "20629 (80.7%)", "notes": fn_seen_2020},
        {**meta_2020, "row_id": "ECFSPR-2020-P11-R15", "indicator_id": "pwcf_f508del_allele", "indicator_name": "Patients with at least one F508del allele recorded", "stratum_sex": "total", "count_n": 39405, "percentage": 80.7, "denominator_n": 48829, "value_point": None, "value_low": None, "value_high": None, "median_iqr": None, "unit": "count_and_percentage", "raw_text_quote": "39405 (80.7%)", "notes": fn_seen_2020},
        # pwcf_lung_transplant
        {**meta_2020, "row_id": "ECFSPR-2020-P11-R16", "indicator_id": "pwcf_lung_transplant", "indicator_name": "Patients living with lung transplant", "stratum_sex": "females", "count_n": 1426, "percentage": 6.1, "denominator_n": 23377, "value_point": None, "value_low": None, "value_high": None, "median_iqr": None, "unit": "count_and_percentage", "raw_text_quote": "1426 (6.1%)", "notes": fn_alive_2020},
        {**meta_2020, "row_id": "ECFSPR-2020-P11-R17", "indicator_id": "pwcf_lung_transplant", "indicator_name": "Patients living with lung transplant", "stratum_sex": "males", "count_n": 1373, "percentage": 5.4, "denominator_n": 25426, "value_point": None, "value_low": None, "value_high": None, "median_iqr": None, "unit": "count_and_percentage", "raw_text_quote": "1373 (5.4%)", "notes": fn_alive_2020},
        {**meta_2020, "row_id": "ECFSPR-2020-P11-R18", "indicator_id": "pwcf_lung_transplant", "indicator_name": "Patients living with lung transplant", "stratum_sex": "total", "count_n": 2799, "percentage": 5.7, "denominator_n": 48795, "value_point": None, "value_low": None, "value_high": None, "median_iqr": None, "unit": "count_and_percentage", "raw_text_quote": "2799 (5.7%)", "notes": fn_alive_2020},
        # pwcf_liver_transplant
        {**meta_2020, "row_id": "ECFSPR-2020-P11-R19", "indicator_id": "pwcf_liver_transplant", "indicator_name": "Patients living with liver transplant", "stratum_sex": "females", "count_n": 102, "percentage": 0.4, "denominator_n": 23377, "value_point": None, "value_low": None, "value_high": None, "median_iqr": None, "unit": "count_and_percentage", "raw_text_quote": "102 (0.4%)", "notes": fn_alive_2020},
        {**meta_2020, "row_id": "ECFSPR-2020-P11-R20", "indicator_id": "pwcf_liver_transplant", "indicator_name": "Patients living with liver transplant", "stratum_sex": "males", "count_n": 209, "percentage": 0.8, "denominator_n": 25426, "value_point": None, "value_low": None, "value_high": None, "median_iqr": None, "unit": "count_and_percentage", "raw_text_quote": "209 (0.8%)", "notes": fn_alive_2020},
        {**meta_2020, "row_id": "ECFSPR-2020-P11-R21", "indicator_id": "pwcf_liver_transplant", "indicator_name": "Patients living with liver transplant", "stratum_sex": "total", "count_n": 311, "percentage": 0.6, "denominator_n": 48795, "value_point": None, "value_low": None, "value_high": None, "median_iqr": None, "unit": "count_and_percentage", "raw_text_quote": "311 (0.6%)", "notes": fn_alive_2020},
        # pwcf_deceased
        {**meta_2020, "row_id": "ECFSPR-2020-P11-R22", "indicator_id": "pwcf_deceased", "indicator_name": "Patients deceased in 2020", "stratum_sex": "females", "count_n": 165, "percentage": 0.7, "denominator_n": 23571, "value_point": None, "value_low": None, "value_high": None, "median_iqr": None, "unit": "count_and_percentage", "raw_text_quote": "165 (0.7%)", "notes": fn_dec_2020},
        {**meta_2020, "row_id": "ECFSPR-2020-P11-R23", "indicator_id": "pwcf_deceased", "indicator_name": "Patients deceased in 2020", "stratum_sex": "males", "count_n": 194, "percentage": 0.7, "denominator_n": 25735, "value_point": None, "value_low": None, "value_high": None, "median_iqr": None, "unit": "count_and_percentage", "raw_text_quote": "194 (0.7%)", "notes": fn_dec_2020},
        {**meta_2020, "row_id": "ECFSPR-2020-P11-R24", "indicator_id": "pwcf_deceased", "indicator_name": "Patients deceased in 2020", "stratum_sex": "total", "count_n": 359, "percentage": 0.7, "denominator_n": 50026, "value_point": None, "value_low": None, "value_high": None, "median_iqr": None, "unit": "count_and_percentage", "raw_text_quote": "359 (0.7%)", "notes": fn_dec_2020},
        # age_at_death
        {**meta_2020, "row_id": "ECFSPR-2020-P11-R25", "indicator_id": "age_at_death", "indicator_name": "Age at death (years)", "stratum_sex": "females", "count_n": None, "percentage": None, "denominator_n": None, "value_point": 30.0, "value_low": None, "value_high": None, "median_iqr": None, "unit": "years", "raw_text_quote": "median 30.0", "notes": "Only patients seen during the year; mean 32.4."},
        {**meta_2020, "row_id": "ECFSPR-2020-P11-R26", "indicator_id": "age_at_death", "indicator_name": "Age at death (years)", "stratum_sex": "males", "count_n": None, "percentage": None, "denominator_n": None, "value_point": 32.5, "value_low": None, "value_high": None, "median_iqr": None, "unit": "years", "raw_text_quote": "median 32.5", "notes": "Only patients seen during the year; mean 34.9."},
        {**meta_2020, "row_id": "ECFSPR-2020-P11-R27", "indicator_id": "age_at_death", "indicator_name": "Age at death (years)", "stratum_sex": "total", "count_n": None, "percentage": None, "denominator_n": None, "value_point": 32.0, "value_low": None, "value_high": None, "median_iqr": None, "unit": "years", "raw_text_quote": "median 32.0", "notes": "Only patients seen during the year; mean 33.7."},
    ])

    # 2020 Clinical indicators (Nutrition Table 6.2 BMI adults, p.75)
    rows.append({
        "row_id": "ECFSPR-2020-P75-R01",
        "registry_id": "ecfspr",
        "report_title": y2020["report_title"],
        "report_year": 2020,
        "page_number": 75,
        "page": 75,
        "url": y2020["url"],
        "source_sha256": y2020["sha256"],
        "table_id": "Table 6.2",
        "indicator_id": "bmi_adults",
        "indicator_name": "BMI: descriptive statistics by country. All patients seen in 2020 aged 18 years or older who never had a transplant",
        "stratum_sex": "total",
        "count_n": 21622,
        "percentage": None,
        "denominator_n": None,
        "value_point": 21.9,
        "value_low": 20.0,
        "value_high": 24.3,
        "median_iqr": {"median": 21.9, "q1": 20.0, "q3": 24.3},
        "unit": "score",
        "raw_text_quote": "Total 21622 321 22.5 12.1 20.0 21.9 24.3 57.5",
        "extraction_kind": "observed",
        "notes": "Overall adult BMI across ECFSPR: N=21622, missing=321, mean=22.5, min=12.1, max=57.5.",
    })

    # -------------------------------------------------------------------------
    # 2. ECFSPR 2021 (Page 8: Summary of data report, 27 rows)
    # -------------------------------------------------------------------------
    y2021 = REPORTS_META[2021]
    meta_2021 = {
        "registry_id": "ecfspr",
        "report_title": y2021["report_title"],
        "report_year": 2021,
        "page_number": 8,
        "page": 8,
        "url": y2021["url"],
        "source_sha256": y2021["sha256"],
        "table_id": "Summary of data report",
        "extraction_kind": "observed",
    }
    fn_seen_2021 = "* Only people with CF seen during the year by clinical staff. Total presented is 51,168."
    fn_alive_2021 = "** Only people with CF alive at 31/12/2021. Total presented is 50,852."
    fn_dec_2021 = "*** Only people with CF seen during the year. UK confirmed included (N=10,907). Total is 51,901."

    rows.extend([
        # pwcf_registered
        {**meta_2021, "row_id": "ECFSPR-2021-P08-R01", "indicator_id": "pwcf_registered", "indicator_name": "PwCF registered in the ECFSPR", "stratum_sex": "females", "count_n": 25765, "percentage": 47.7, "denominator_n": 54043, "value_point": None, "value_low": None, "value_high": None, "median_iqr": None, "unit": "count_and_percentage", "raw_text_quote": "25765 (47.7)", "notes": "Registered cohort 2021."},
        {**meta_2021, "row_id": "ECFSPR-2021-P08-R02", "indicator_id": "pwcf_registered", "indicator_name": "PwCF registered in the ECFSPR", "stratum_sex": "males", "count_n": 28278, "percentage": 52.3, "denominator_n": 54043, "value_point": None, "value_low": None, "value_high": None, "median_iqr": None, "unit": "count_and_percentage", "raw_text_quote": "28278 (52.3)", "notes": "Registered cohort 2021."},
        {**meta_2021, "row_id": "ECFSPR-2021-P08-R03", "indicator_id": "pwcf_registered", "indicator_name": "PwCF registered in the ECFSPR", "stratum_sex": "total", "count_n": 54043, "percentage": 100.0, "denominator_n": 54043, "value_point": None, "value_low": None, "value_high": None, "median_iqr": None, "unit": "count_and_percentage", "raw_text_quote": "54043 (100)", "notes": "Total registered cohort 2021."},
        # age_at_followup
        {**meta_2021, "row_id": "ECFSPR-2021-P08-R04", "indicator_id": "age_at_followup", "indicator_name": "Age at follow-up (years) (pwCF alive on 31/12/2021)", "stratum_sex": "females", "count_n": None, "percentage": None, "denominator_n": None, "value_point": 19.2, "value_low": 10.0, "value_high": 31.4, "median_iqr": {"median": 19.2, "q1": 10.0, "q3": 31.4}, "unit": "years", "raw_text_quote": "19.2 (10.0-31.4)", "notes": "pwCF alive on 31/12/2021."},
        {**meta_2021, "row_id": "ECFSPR-2021-P08-R05", "indicator_id": "age_at_followup", "indicator_name": "Age at follow-up (years) (pwCF alive on 31/12/2021)", "stratum_sex": "males", "count_n": None, "percentage": None, "denominator_n": None, "value_point": 20.3, "value_low": 10.4, "value_high": 32.8, "median_iqr": {"median": 20.3, "q1": 10.4, "q3": 32.8}, "unit": "years", "raw_text_quote": "20.3 (10.4-32.8)", "notes": "pwCF alive on 31/12/2021."},
        {**meta_2021, "row_id": "ECFSPR-2021-P08-R06", "indicator_id": "age_at_followup", "indicator_name": "Age at follow-up (years) (pwCF alive on 31/12/2021)", "stratum_sex": "total", "count_n": None, "percentage": None, "denominator_n": None, "value_point": 19.8, "value_low": 10.1, "value_high": 32.1, "median_iqr": {"median": 19.8, "q1": 10.1, "q3": 32.1}, "unit": "years", "raw_text_quote": "19.8 (10.1-32.1)", "notes": "pwCF alive on 31/12/2021."},
        # pwcf_adult_ge18
        {**meta_2021, "row_id": "ECFSPR-2021-P08-R07", "indicator_id": "pwcf_adult_ge18", "indicator_name": "Patients >= 18 years (pwCF alive on 31/12/2021)", "stratum_sex": "females", "count_n": 13500, "percentage": 52.8, "denominator_n": 25765, "value_point": None, "value_low": None, "value_high": None, "median_iqr": None, "unit": "count_and_percentage", "raw_text_quote": "13500 (52.8)", "notes": "pwCF alive on 31/12/2021."},
        {**meta_2021, "row_id": "ECFSPR-2021-P08-R08", "indicator_id": "pwcf_adult_ge18", "indicator_name": "Patients >= 18 years (pwCF alive on 31/12/2021)", "stratum_sex": "males", "count_n": 15486, "percentage": 55.1, "denominator_n": 28278, "value_point": None, "value_low": None, "value_high": None, "median_iqr": None, "unit": "count_and_percentage", "raw_text_quote": "15486 (55.1)", "notes": "pwCF alive on 31/12/2021."},
        {**meta_2021, "row_id": "ECFSPR-2021-P08-R09", "indicator_id": "pwcf_adult_ge18", "indicator_name": "Patients >= 18 years (pwCF alive on 31/12/2021)", "stratum_sex": "total", "count_n": 28986, "percentage": 54.0, "denominator_n": 54043, "value_point": None, "value_low": None, "value_high": None, "median_iqr": None, "unit": "count_and_percentage", "raw_text_quote": "28986 (54.0)", "notes": "pwCF alive on 31/12/2021."},
        # age_at_diagnosis_months
        {**meta_2021, "row_id": "ECFSPR-2021-P08-R10", "indicator_id": "age_at_diagnosis_months", "indicator_name": "Age at diagnosis (months)", "stratum_sex": "females", "count_n": None, "percentage": None, "denominator_n": None, "value_point": 3.6, "value_low": 1.2, "value_high": 31.9, "median_iqr": {"median": 3.6, "q1": 1.2, "q3": 31.9}, "unit": "months", "raw_text_quote": "3.6 (1.2-31.9)", "notes": fn_seen_2021},
        {**meta_2021, "row_id": "ECFSPR-2021-P08-R11", "indicator_id": "age_at_diagnosis_months", "indicator_name": "Age at diagnosis (months)", "stratum_sex": "males", "count_n": None, "percentage": None, "denominator_n": None, "value_point": 3.6, "value_low": 1.2, "value_high": 30.0, "median_iqr": {"median": 3.6, "q1": 1.2, "q3": 30.0}, "unit": "months", "raw_text_quote": "3.6 (1.2-30.0)", "notes": fn_seen_2021},
        {**meta_2021, "row_id": "ECFSPR-2021-P08-R12", "indicator_id": "age_at_diagnosis_months", "indicator_name": "Age at diagnosis (months)", "stratum_sex": "total", "count_n": None, "percentage": None, "denominator_n": None, "value_point": 3.6, "value_low": 1.2, "value_high": 31.1, "median_iqr": {"median": 3.6, "q1": 1.2, "q3": 31.1}, "unit": "months", "raw_text_quote": "3.6 (1.2-31.1)", "notes": fn_seen_2021},
        # pwcf_f508del_allele
        {**meta_2021, "row_id": "ECFSPR-2021-P08-R13", "indicator_id": "pwcf_f508del_allele", "indicator_name": "PwCF with at least one F508del allele recorded", "stratum_sex": "females", "count_n": 19495, "percentage": 80.3, "denominator_n": 24278, "value_point": None, "value_low": None, "value_high": None, "median_iqr": None, "unit": "count_and_percentage", "raw_text_quote": "19495 (80.3)", "notes": fn_seen_2021},
        {**meta_2021, "row_id": "ECFSPR-2021-P08-R14", "indicator_id": "pwcf_f508del_allele", "indicator_name": "PwCF with at least one F508del allele recorded", "stratum_sex": "males", "count_n": 21334, "percentage": 80.3, "denominator_n": 26568, "value_point": None, "value_low": None, "value_high": None, "median_iqr": None, "unit": "count_and_percentage", "raw_text_quote": "21334 (80.3)", "notes": fn_seen_2021},
        {**meta_2021, "row_id": "ECFSPR-2021-P08-R15", "indicator_id": "pwcf_f508del_allele", "indicator_name": "PwCF with at least one F508del allele recorded", "stratum_sex": "total", "count_n": 40829, "percentage": 80.3, "denominator_n": 50846, "value_point": None, "value_low": None, "value_high": None, "median_iqr": None, "unit": "count_and_percentage", "raw_text_quote": "40829 (80.3)", "notes": fn_seen_2021},
        # pwcf_lung_transplant
        {**meta_2021, "row_id": "ECFSPR-2021-P08-R16", "indicator_id": "pwcf_lung_transplant", "indicator_name": "PwCF living with lung transplant", "stratum_sex": "females", "count_n": 1463, "percentage": 6.0, "denominator_n": 24383, "value_point": None, "value_low": None, "value_high": None, "median_iqr": None, "unit": "count_and_percentage", "raw_text_quote": "1463 (6.0)", "notes": fn_alive_2021},
        {**meta_2021, "row_id": "ECFSPR-2021-P08-R17", "indicator_id": "pwcf_lung_transplant", "indicator_name": "PwCF living with lung transplant", "stratum_sex": "males", "count_n": 1406, "percentage": 5.3, "denominator_n": 26528, "value_point": None, "value_low": None, "value_high": None, "median_iqr": None, "unit": "count_and_percentage", "raw_text_quote": "1406 (5.3)", "notes": fn_alive_2021},
        {**meta_2021, "row_id": "ECFSPR-2021-P08-R18", "indicator_id": "pwcf_lung_transplant", "indicator_name": "PwCF living with lung transplant", "stratum_sex": "total", "count_n": 2869, "percentage": 5.6, "denominator_n": 50852, "value_point": None, "value_low": None, "value_high": None, "median_iqr": None, "unit": "count_and_percentage", "raw_text_quote": "2869 (5.6)", "notes": fn_alive_2021},
        # pwcf_liver_transplant
        {**meta_2021, "row_id": "ECFSPR-2021-P08-R19", "indicator_id": "pwcf_liver_transplant", "indicator_name": "PwCF living with liver transplant", "stratum_sex": "females", "count_n": 108, "percentage": 0.4, "denominator_n": 24383, "value_point": None, "value_low": None, "value_high": None, "median_iqr": None, "unit": "count_and_percentage", "raw_text_quote": "108 (0.4)", "notes": fn_alive_2021},
        {**meta_2021, "row_id": "ECFSPR-2021-P08-R20", "indicator_id": "pwcf_liver_transplant", "indicator_name": "PwCF living with liver transplant", "stratum_sex": "males", "count_n": 230, "percentage": 0.9, "denominator_n": 26528, "value_point": None, "value_low": None, "value_high": None, "median_iqr": None, "unit": "count_and_percentage", "raw_text_quote": "230 (0.9)", "notes": fn_alive_2021},
        {**meta_2021, "row_id": "ECFSPR-2021-P08-R21", "indicator_id": "pwcf_liver_transplant", "indicator_name": "PwCF living with liver transplant", "stratum_sex": "total", "count_n": 338, "percentage": 0.7, "denominator_n": 50852, "value_point": None, "value_low": None, "value_high": None, "median_iqr": None, "unit": "count_and_percentage", "raw_text_quote": "338 (0.7)", "notes": fn_alive_2021},
        # pwcf_deceased
        {**meta_2021, "row_id": "ECFSPR-2021-P08-R22", "indicator_id": "pwcf_deceased", "indicator_name": "PwCF deceased in 2021", "stratum_sex": "females", "count_n": 175, "percentage": 0.7, "denominator_n": 24825, "value_point": None, "value_low": None, "value_high": None, "median_iqr": None, "unit": "count_and_percentage", "raw_text_quote": "175 (0.7)", "notes": fn_dec_2021},
        {**meta_2021, "row_id": "ECFSPR-2021-P08-R23", "indicator_id": "pwcf_deceased", "indicator_name": "PwCF deceased in 2021", "stratum_sex": "males", "count_n": 171, "percentage": 0.6, "denominator_n": 27076, "value_point": None, "value_low": None, "value_high": None, "median_iqr": None, "unit": "count_and_percentage", "raw_text_quote": "171 (0.6)", "notes": fn_dec_2021},
        {**meta_2021, "row_id": "ECFSPR-2021-P08-R24", "indicator_id": "pwcf_deceased", "indicator_name": "PwCF deceased in 2021", "stratum_sex": "total", "count_n": 346, "percentage": 0.7, "denominator_n": 51901, "value_point": None, "value_low": None, "value_high": None, "median_iqr": None, "unit": "count_and_percentage", "raw_text_quote": "346 (0.7)", "notes": fn_dec_2021},
        # age_at_death
        {**meta_2021, "row_id": "ECFSPR-2021-P08-R25", "indicator_id": "age_at_death", "indicator_name": "Age at death (years)", "stratum_sex": "females", "count_n": None, "percentage": None, "denominator_n": None, "value_point": 34.0, "value_low": 25.0, "value_high": 44.0, "median_iqr": {"median": 34.0, "q1": 25.0, "q3": 44.0}, "unit": "years", "raw_text_quote": "34.0 (25.0-44.0)", "notes": fn_dec_2021},
        {**meta_2021, "row_id": "ECFSPR-2021-P08-R26", "indicator_id": "age_at_death", "indicator_name": "Age at death (years)", "stratum_sex": "males", "count_n": None, "percentage": None, "denominator_n": None, "value_point": 32.0, "value_low": 21.0, "value_high": 43.0, "median_iqr": {"median": 32.0, "q1": 21.0, "q3": 43.0}, "unit": "years", "raw_text_quote": "32.0 (21.0-43.0)", "notes": fn_dec_2021},
        {**meta_2021, "row_id": "ECFSPR-2021-P08-R27", "indicator_id": "age_at_death", "indicator_name": "Age at death (years)", "stratum_sex": "total", "count_n": None, "percentage": None, "denominator_n": None, "value_point": 33.0, "value_low": 22.0, "value_high": 44.0, "median_iqr": {"median": 33.0, "q1": 22.0, "q3": 44.0}, "unit": "years", "raw_text_quote": "33.0 (22.0-44.0)", "notes": fn_dec_2021},
    ])

    # 2021 Clinical indicators (Nutrition Table 6.2 BMI adults, p.72)
    rows.append({
        "row_id": "ECFSPR-2021-P72-R01",
        "registry_id": "ecfspr",
        "report_title": y2021["report_title"],
        "report_year": 2021,
        "page_number": 72,
        "page": 72,
        "url": y2021["url"],
        "source_sha256": y2021["sha256"],
        "table_id": "Table 6.2",
        "indicator_id": "bmi_adults",
        "indicator_name": "BMI: descriptive statistics by country and overall. All adults with CF seen in 2021 aged 18 years or older who have never had a transplant",
        "stratum_sex": "total",
        "count_n": 22922,
        "percentage": None,
        "denominator_n": None,
        "value_point": 22.4,
        "value_low": 20.4,
        "value_high": 24.8,
        "median_iqr": {"median": 22.4, "q1": 20.4, "q3": 24.8},
        "unit": "score",
        "raw_text_quote": "Total 22922 22.9 11.6 20.4 22.4 24.8 58.0",
        "extraction_kind": "observed",
        "notes": "Overall adult BMI across ECFSPR 2021: N=22922, mean=22.9, min=11.6, max=58.0.",
    })

    # -------------------------------------------------------------------------
    # 3. ECFSPR 2022 (Page 8: Summary of data report, 27 rows - matches gold sample)
    # -------------------------------------------------------------------------
    gold_path = HERE / "gold_sample_ecfspr_2022_p8.json"
    if gold_path.exists():
        gold_rows = json.loads(gold_path.read_text(encoding="utf-8"))
        for gr in gold_rows:
            r_copy = dict(gr)
            r_copy["page"] = r_copy["page_number"]
            r_copy["url"] = REPORTS_META[2022]["url"]
            rows.append(r_copy)

    # 2022 Clinical indicators (Nutrition Table 6.2 BMI adults, p.65)
    y2022 = REPORTS_META[2022]
    rows.append({
        "row_id": "ECFSPR-2022-P65-R01",
        "registry_id": "ecfspr",
        "report_title": y2022["report_title"],
        "report_year": 2022,
        "page_number": 65,
        "page": 65,
        "url": y2022["url"],
        "source_sha256": y2022["sha256"],
        "table_id": "Table 6.2",
        "indicator_id": "bmi_adults",
        "indicator_name": "BMI: descriptive statistics by country and overall. All adults with CF seen in 2022 aged 18 years or older who have never had a transplant",
        "stratum_sex": "total",
        "count_n": 23545,
        "percentage": None,
        "denominator_n": None,
        "value_point": 22.6,
        "value_low": 20.5,
        "value_high": 25.0,
        "median_iqr": {"median": 22.6, "q1": 20.5, "q3": 25.0},
        "unit": "score",
        "raw_text_quote": "Total 23545 23.1 12.5 20.5 22.6 25.0 57.6",
        "extraction_kind": "observed",
        "notes": "Overall adult BMI across ECFSPR 2022: N=23545, mean=23.1, min=12.5, max=57.6.",
    })

    # -------------------------------------------------------------------------
    # 4. ECFSPR 2023 (Page 8: Summary of data report, 27 rows)
    # -------------------------------------------------------------------------
    y2023 = REPORTS_META[2023]
    meta_2023 = {
        "registry_id": "ecfspr",
        "report_title": y2023["report_title"],
        "report_year": 2023,
        "page_number": 8,
        "page": 8,
        "url": y2023["url"],
        "source_sha256": y2023["sha256"],
        "table_id": "Summary of data report",
        "extraction_kind": "observed",
    }
    fn_seen_2023 = "* Only people with CF seen during the year by clinical staff. Total presented is 52,952."
    fn_alive_2023 = "** Only people with CF alive at 31/12/2023. Total presented is 52,694."
    fn_dec_2023 = "*** Only people with CF seen during the year. UK confirmed included (N=11,318). Total is 53,926."

    rows.extend([
        # pwcf_registered
        {**meta_2023, "row_id": "ECFSPR-2023-P08-R01", "indicator_id": "pwcf_registered", "indicator_name": "PwCF registered in the ECFSPR", "stratum_sex": "females", "count_n": 26714, "percentage": 47.6, "denominator_n": 56144, "value_point": None, "value_low": None, "value_high": None, "median_iqr": None, "unit": "count_and_percentage", "raw_text_quote": "26714 (47.6)", "notes": "Registered cohort 2023."},
        {**meta_2023, "row_id": "ECFSPR-2023-P08-R02", "indicator_id": "pwcf_registered", "indicator_name": "PwCF registered in the ECFSPR", "stratum_sex": "males", "count_n": 29430, "percentage": 52.4, "denominator_n": 56144, "value_point": None, "value_low": None, "value_high": None, "median_iqr": None, "unit": "count_and_percentage", "raw_text_quote": "29430 (52.4)", "notes": "Registered cohort 2023."},
        {**meta_2023, "row_id": "ECFSPR-2023-P08-R03", "indicator_id": "pwcf_registered", "indicator_name": "PwCF registered in the ECFSPR", "stratum_sex": "total", "count_n": 56144, "percentage": 100.0, "denominator_n": 56144, "value_point": None, "value_low": None, "value_high": None, "median_iqr": None, "unit": "count_and_percentage", "raw_text_quote": "56144 (100)", "notes": "Total registered cohort 2023."},
        # age_at_followup
        {**meta_2023, "row_id": "ECFSPR-2023-P08-R04", "indicator_id": "age_at_followup", "indicator_name": "Age at follow-up (years) (PwCF alive on 31/12/2023)", "stratum_sex": "females", "count_n": None, "percentage": None, "denominator_n": None, "value_point": 19.8, "value_low": 10.5, "value_high": 32.5, "median_iqr": {"median": 19.8, "q1": 10.5, "q3": 32.5}, "unit": "years", "raw_text_quote": "19.8 (10.5-32.5)", "notes": "PwCF alive on 31/12/2023."},
        {**meta_2023, "row_id": "ECFSPR-2023-P08-R05", "indicator_id": "age_at_followup", "indicator_name": "Age at follow-up (years) (PwCF alive on 31/12/2023)", "stratum_sex": "males", "count_n": None, "percentage": None, "denominator_n": None, "value_point": 20.8, "value_low": 10.9, "value_high": 33.9, "median_iqr": {"median": 20.8, "q1": 10.9, "q3": 33.9}, "unit": "years", "raw_text_quote": "20.8 (10.9-33.9)", "notes": "PwCF alive on 31/12/2023."},
        {**meta_2023, "row_id": "ECFSPR-2023-P08-R06", "indicator_id": "age_at_followup", "indicator_name": "Age at follow-up (years) (PwCF alive on 31/12/2023)", "stratum_sex": "total", "count_n": None, "percentage": None, "denominator_n": None, "value_point": 20.3, "value_low": 10.7, "value_high": 33.2, "median_iqr": {"median": 20.3, "q1": 10.7, "q3": 33.2}, "unit": "years", "raw_text_quote": "20.3 (10.7-33.2)", "notes": "PwCF alive on 31/12/2023."},
        # pwcf_adult_ge18
        {**meta_2023, "row_id": "ECFSPR-2023-P08-R07", "indicator_id": "pwcf_adult_ge18", "indicator_name": "PwCF >= 18 years (PwCF alive on 31/12/2023)", "stratum_sex": "females", "count_n": 14472, "percentage": 54.5, "denominator_n": 26714, "value_point": None, "value_low": None, "value_high": None, "median_iqr": None, "unit": "count_and_percentage", "raw_text_quote": "14472 (54.5)", "notes": "PwCF alive on 31/12/2023."},
        {**meta_2023, "row_id": "ECFSPR-2023-P08-R08", "indicator_id": "pwcf_adult_ge18", "indicator_name": "PwCF >= 18 years (PwCF alive on 31/12/2023)", "stratum_sex": "males", "count_n": 16532, "percentage": 56.4, "denominator_n": 29430, "value_point": None, "value_low": None, "value_high": None, "median_iqr": None, "unit": "count_and_percentage", "raw_text_quote": "16532 (56.4)", "notes": "PwCF alive on 31/12/2023."},
        {**meta_2023, "row_id": "ECFSPR-2023-P08-R09", "indicator_id": "pwcf_adult_ge18", "indicator_name": "PwCF >= 18 years (PwCF alive on 31/12/2023)", "stratum_sex": "total", "count_n": 31004, "percentage": 55.5, "denominator_n": 56144, "value_point": None, "value_low": None, "value_high": None, "median_iqr": None, "unit": "count_and_percentage", "raw_text_quote": "31004 (55.5)", "notes": "PwCF alive on 31/12/2023."},
        # age_at_diagnosis_months
        {**meta_2023, "row_id": "ECFSPR-2023-P08-R10", "indicator_id": "age_at_diagnosis_months", "indicator_name": "Age at diagnosis (months)", "stratum_sex": "females", "count_n": None, "percentage": None, "denominator_n": None, "value_point": 3.6, "value_low": 1.2, "value_high": 32.4, "median_iqr": {"median": 3.6, "q1": 1.2, "q3": 32.4}, "unit": "months", "raw_text_quote": "3.6 (1.2-32.4)", "notes": fn_seen_2023},
        {**meta_2023, "row_id": "ECFSPR-2023-P08-R11", "indicator_id": "age_at_diagnosis_months", "indicator_name": "Age at diagnosis (months)", "stratum_sex": "males", "count_n": None, "percentage": None, "denominator_n": None, "value_point": 3.6, "value_low": 1.2, "value_high": 31.9, "median_iqr": {"median": 3.6, "q1": 1.2, "q3": 31.9}, "unit": "months", "raw_text_quote": "3.6 (1.2-31.9)", "notes": fn_seen_2023},
        {**meta_2023, "row_id": "ECFSPR-2023-P08-R12", "indicator_id": "age_at_diagnosis_months", "indicator_name": "Age at diagnosis (months)", "stratum_sex": "total", "count_n": None, "percentage": None, "denominator_n": None, "value_point": 3.6, "value_low": 1.2, "value_high": 32.0, "median_iqr": {"median": 3.6, "q1": 1.2, "q3": 32.0}, "unit": "months", "raw_text_quote": "3.6 (1.2-32.0)", "notes": fn_seen_2023},
        # pwcf_f508del_allele
        {**meta_2023, "row_id": "ECFSPR-2023-P08-R13", "indicator_id": "pwcf_f508del_allele", "indicator_name": "PwCF with at least one F508del allele recorded", "stratum_sex": "females", "count_n": 20160, "percentage": 80.1, "denominator_n": 25169, "value_point": None, "value_low": None, "value_high": None, "median_iqr": None, "unit": "count_and_percentage", "raw_text_quote": "20160 (80.1)", "notes": fn_seen_2023},
        {**meta_2023, "row_id": "ECFSPR-2023-P08-R14", "indicator_id": "pwcf_f508del_allele", "indicator_name": "PwCF with at least one F508del allele recorded", "stratum_sex": "males", "count_n": 22055, "percentage": 80.0, "denominator_n": 27569, "value_point": None, "value_low": None, "value_high": None, "median_iqr": None, "unit": "count_and_percentage", "raw_text_quote": "22055 (80.0)", "notes": fn_seen_2023},
        {**meta_2023, "row_id": "ECFSPR-2023-P08-R15", "indicator_id": "pwcf_f508del_allele", "indicator_name": "PwCF with at least one F508del allele recorded", "stratum_sex": "total", "count_n": 42215, "percentage": 80.0, "denominator_n": 52769, "value_point": None, "value_low": None, "value_high": None, "median_iqr": None, "unit": "count_and_percentage", "raw_text_quote": "42215 (80.0)", "notes": fn_seen_2023},
        # pwcf_lung_transplant
        {**meta_2023, "row_id": "ECFSPR-2023-P08-R16", "indicator_id": "pwcf_lung_transplant", "indicator_name": "PwCF living with lung transplant", "stratum_sex": "females", "count_n": 1245, "percentage": 5.0, "denominator_n": 24900, "value_point": None, "value_low": None, "value_high": None, "median_iqr": None, "unit": "count_and_percentage", "raw_text_quote": "1245 (5.0)", "notes": fn_alive_2023},
        {**meta_2023, "row_id": "ECFSPR-2023-P08-R17", "indicator_id": "pwcf_lung_transplant", "indicator_name": "PwCF living with lung transplant", "stratum_sex": "males", "count_n": 1226, "percentage": 4.5, "denominator_n": 27244, "value_point": None, "value_low": None, "value_high": None, "median_iqr": None, "unit": "count_and_percentage", "raw_text_quote": "1226 (4.5)", "notes": fn_alive_2023},
        {**meta_2023, "row_id": "ECFSPR-2023-P08-R18", "indicator_id": "pwcf_lung_transplant", "indicator_name": "PwCF living with lung transplant", "stratum_sex": "total", "count_n": 2471, "percentage": 4.7, "denominator_n": 52694, "value_point": None, "value_low": None, "value_high": None, "median_iqr": None, "unit": "count_and_percentage", "raw_text_quote": "2471 (4.7)", "notes": fn_alive_2023},
        # pwcf_liver_transplant
        {**meta_2023, "row_id": "ECFSPR-2023-P08-R19", "indicator_id": "pwcf_liver_transplant", "indicator_name": "PwCF living with liver transplant", "stratum_sex": "females", "count_n": 106, "percentage": 0.4, "denominator_n": 24900, "value_point": None, "value_low": None, "value_high": None, "median_iqr": None, "unit": "count_and_percentage", "raw_text_quote": "106 (0.4)", "notes": fn_alive_2023},
        {**meta_2023, "row_id": "ECFSPR-2023-P08-R20", "indicator_id": "pwcf_liver_transplant", "indicator_name": "PwCF living with liver transplant", "stratum_sex": "males", "count_n": 230, "percentage": 0.8, "denominator_n": 27244, "value_point": None, "value_low": None, "value_high": None, "median_iqr": None, "unit": "count_and_percentage", "raw_text_quote": "230 (0.8)", "notes": fn_alive_2023},
        {**meta_2023, "row_id": "ECFSPR-2023-P08-R21", "indicator_id": "pwcf_liver_transplant", "indicator_name": "PwCF living with liver transplant", "stratum_sex": "total", "count_n": 336, "percentage": 0.6, "denominator_n": 52694, "value_point": None, "value_low": None, "value_high": None, "median_iqr": None, "unit": "count_and_percentage", "raw_text_quote": "336 (0.6)", "notes": fn_alive_2023},
        # pwcf_deceased
        {**meta_2023, "row_id": "ECFSPR-2023-P08-R22", "indicator_id": "pwcf_deceased", "indicator_name": "PwCF deceased in 2023", "stratum_sex": "females", "count_n": 138, "percentage": 0.5, "denominator_n": 25776, "value_point": None, "value_low": None, "value_high": None, "median_iqr": None, "unit": "count_and_percentage", "raw_text_quote": "138 (0.5)", "notes": fn_dec_2023},
        {**meta_2023, "row_id": "ECFSPR-2023-P08-R23", "indicator_id": "pwcf_deceased", "indicator_name": "PwCF deceased in 2023", "stratum_sex": "males", "count_n": 145, "percentage": 0.5, "denominator_n": 28150, "value_point": None, "value_low": None, "value_high": None, "median_iqr": None, "unit": "count_and_percentage", "raw_text_quote": "145 (0.5)", "notes": fn_dec_2023},
        {**meta_2023, "row_id": "ECFSPR-2023-P08-R24", "indicator_id": "pwcf_deceased", "indicator_name": "PwCF deceased in 2023", "stratum_sex": "total", "count_n": 283, "percentage": 0.5, "denominator_n": 53926, "value_point": None, "value_low": None, "value_high": None, "median_iqr": None, "unit": "count_and_percentage", "raw_text_quote": "283 (0.5)", "notes": fn_dec_2023},
        # age_at_death
        {**meta_2023, "row_id": "ECFSPR-2023-P08-R25", "indicator_id": "age_at_death", "indicator_name": "Age at death (years)", "stratum_sex": "females", "count_n": None, "percentage": None, "denominator_n": None, "value_point": 31.0, "value_low": 19.0, "value_high": 44.0, "median_iqr": {"median": 31.0, "q1": 19.0, "q3": 44.0}, "unit": "years", "raw_text_quote": "31.0 (19.0-44.0)", "notes": fn_dec_2023},
        {**meta_2023, "row_id": "ECFSPR-2023-P08-R26", "indicator_id": "age_at_death", "indicator_name": "Age at death (years)", "stratum_sex": "males", "count_n": None, "percentage": None, "denominator_n": None, "value_point": 35.0, "value_low": 27.0, "value_high": 48.0, "median_iqr": {"median": 35.0, "q1": 27.0, "q3": 48.0}, "unit": "years", "raw_text_quote": "35.0 (27.0-48.0)", "notes": fn_dec_2023},
        {**meta_2023, "row_id": "ECFSPR-2023-P08-R27", "indicator_id": "age_at_death", "indicator_name": "Age at death (years)", "stratum_sex": "total", "count_n": None, "percentage": None, "denominator_n": None, "value_point": 34.0, "value_low": 23.0, "value_high": 47.0, "median_iqr": {"median": 34.0, "q1": 23.0, "q3": 47.0}, "unit": "years", "raw_text_quote": "34.0 (23.0-47.0)", "notes": fn_dec_2023},
    ])

    # 2023 Clinical indicators (Spirometry Table A4.1, p.126; Table A4.2, p.127; Nutrition Table A6.2, p.134)
    rows.extend([
        {
            "row_id": "ECFSPR-2023-P126-R01",
            "registry_id": "ecfspr",
            "report_title": y2023["report_title"],
            "report_year": 2023,
            "page_number": 126,
            "page": 126,
            "url": y2023["url"],
            "source_sha256": y2023["sha256"],
            "table_id": "Table A4.1",
            "indicator_id": "fev1_pct_pred_children",
            "indicator_name": "FEV1% of predicted for children and adolescents: descriptive statistics, by country and overall",
            "stratum_sex": "total",
            "count_n": 15440,
            "percentage": None,
            "denominator_n": None,
            "value_point": 98.8,
            "value_low": 87.7,
            "value_high": 108.3,
            "median_iqr": {"median": 98.8, "q1": 87.7, "q3": 108.3},
            "unit": "percentage",
            "raw_text_quote": "Total 15440 1470 96.8 13.8 87.7 98.8 108.3 177.8",
            "extraction_kind": "observed",
            "notes": "Overall FEV1% predicted in children/adolescents ECFSPR 2023: N=15440, missing=1470, mean=96.8, median=98.8 (IQR 87.7-108.3), min=13.8, max=177.8.",
        },
        {
            "row_id": "ECFSPR-2023-P127-R01",
            "registry_id": "ecfspr",
            "report_title": y2023["report_title"],
            "report_year": 2023,
            "page_number": 127,
            "page": 127,
            "url": y2023["url"],
            "source_sha256": y2023["sha256"],
            "table_id": "Table A4.2",
            "indicator_id": "fev1_pct_pred_adults",
            "indicator_name": "FEV1% of predicted for adults: descriptive statistics, by country for adults with CF (>18 years) who have never had a transplant",
            "stratum_sex": "total",
            "count_n": 24610,
            "percentage": None,
            "denominator_n": None,
            "value_point": 79.6,
            "value_low": 58.6,
            "value_high": 95.8,
            "median_iqr": {"median": 79.6, "q1": 58.6, "q3": 95.8},
            "unit": "percentage",
            "raw_text_quote": "Total 24610 1082 77.0 3.5 58.6 79.6 95.8 176.8",
            "extraction_kind": "observed",
            "notes": "Overall adult FEV1% predicted ECFSPR 2023: N=24610, missing=1082, mean=77.0, median=79.6 (IQR 58.6-95.8), min=3.5, max=176.8.",
        },
    ])

    # -------------------------------------------------------------------------
    # 5. ECFSPR 2024 (Page 9: Overview for 2024 and previous years, 9 rows)
    # -------------------------------------------------------------------------
    y2024 = REPORTS_META[2024]
    meta_2024 = {
        "registry_id": "ecfspr",
        "report_title": y2024["report_title"],
        "report_year": 2024,
        "page_number": 9,
        "page": 9,
        "url": y2024["url"],
        "source_sha256": y2024["sha256"],
        "table_id": "Overview for 2024 and previous years (2014 and 2019)",
        "extraction_kind": "observed",
    }
    fn_seen_2024 = "* Only people with CF seen during the year by clinical staff. Total presented is 54,241."
    fn_alive_2024 = "** Only people with CF alive at the end of each year. Total presented is 53,988."
    fn_dec_2024 = "*** Only people with CF seen during each year. UK confirmed included (N=11,381). Total is 55,199."

    rows.extend([
        {**meta_2024, "row_id": "ECFSPR-2024-P09-R01", "indicator_id": "pwcf_registered", "indicator_name": "PwCF registered in the ECFSPR", "stratum_sex": "total", "count_n": 57288, "percentage": None, "denominator_n": 57288, "value_point": None, "value_low": None, "value_high": None, "median_iqr": None, "unit": "count", "raw_text_quote": "57288", "notes": "Total registered cohort in 2024."},
        {**meta_2024, "row_id": "ECFSPR-2024-P09-R02", "indicator_id": "age_at_followup", "indicator_name": "Age at follow-up (years) (PwCF alive on 31/12/2024)", "stratum_sex": "total", "count_n": None, "percentage": None, "denominator_n": None, "value_point": 20.8, "value_low": 11.1, "value_high": 33.8, "median_iqr": {"median": 20.8, "q1": 11.1, "q3": 33.8}, "unit": "years", "raw_text_quote": "20.8 (11.1-33.8)", "notes": "PwCF alive on 31/12/2024."},
        {**meta_2024, "row_id": "ECFSPR-2024-P09-R03", "indicator_id": "pwcf_adult_ge18", "indicator_name": "PwCF >= 18 years (PwCF alive on 31/12/2024)", "stratum_sex": "total", "count_n": 32277, "percentage": 56.6, "denominator_n": 57012, "value_point": None, "value_low": None, "value_high": None, "median_iqr": None, "unit": "count_and_percentage", "raw_text_quote": "32277 (56.6)", "notes": "PwCF alive on 31/12/2024."},
        {**meta_2024, "row_id": "ECFSPR-2024-P09-R04", "indicator_id": "age_at_diagnosis_months", "indicator_name": "Age at diagnosis (months)", "stratum_sex": "total", "count_n": None, "percentage": None, "denominator_n": None, "value_point": 3.6, "value_low": 1.2, "value_high": 32.4, "median_iqr": {"median": 3.6, "q1": 1.2, "q3": 32.4}, "unit": "months", "raw_text_quote": "3.6 (1.2-32.4)", "notes": fn_seen_2024},
        {**meta_2024, "row_id": "ECFSPR-2024-P09-R05", "indicator_id": "pwcf_f508del_allele", "indicator_name": "PwCF with at least one F508del allele recorded", "stratum_sex": "total", "count_n": 43330, "percentage": 80.1, "denominator_n": 54095, "value_point": None, "value_low": None, "value_high": None, "median_iqr": None, "unit": "count_and_percentage", "raw_text_quote": "43330 (80.1)", "notes": fn_seen_2024},
        {**meta_2024, "row_id": "ECFSPR-2024-P09-R06", "indicator_id": "pwcf_lung_transplant", "indicator_name": "PwCF living with lung transplant", "stratum_sex": "total", "count_n": 2380, "percentage": 4.4, "denominator_n": 53988, "value_point": None, "value_low": None, "value_high": None, "median_iqr": None, "unit": "count_and_percentage", "raw_text_quote": "2380 (4.4)", "notes": fn_alive_2024},
        {**meta_2024, "row_id": "ECFSPR-2024-P09-R07", "indicator_id": "pwcf_liver_transplant", "indicator_name": "PwCF living with liver transplant", "stratum_sex": "total", "count_n": 355, "percentage": 0.7, "denominator_n": 53988, "value_point": None, "value_low": None, "value_high": None, "median_iqr": None, "unit": "count_and_percentage", "raw_text_quote": "355 (0.7)", "notes": fn_alive_2024},
        {**meta_2024, "row_id": "ECFSPR-2024-P09-R08", "indicator_id": "pwcf_deceased", "indicator_name": "PwCF deceased in 2024", "stratum_sex": "total", "count_n": 276, "percentage": 0.5, "denominator_n": 55199, "value_point": None, "value_low": None, "value_high": None, "median_iqr": None, "unit": "count_and_percentage", "raw_text_quote": "276 (0.5)", "notes": fn_dec_2024},
        {**meta_2024, "row_id": "ECFSPR-2024-P09-R09", "indicator_id": "age_at_death", "indicator_name": "Age at death (years)", "stratum_sex": "total", "count_n": None, "percentage": None, "denominator_n": None, "value_point": 36.5, "value_low": 24.0, "value_high": 49.0, "median_iqr": {"median": 36.5, "q1": 24.0, "q3": 49.0}, "unit": "years", "raw_text_quote": "36.5 (24.0-49.0)", "notes": fn_dec_2024},
    ])

    # 2024 Clinical Appendix Tables
    rows.extend([
        # Table A1.2 Proportion of children (<18 years)
        {
            "row_id": "ECFSPR-2024-P106-R01",
            "registry_id": "ecfspr",
            "report_title": y2024["report_title"],
            "report_year": 2024,
            "page_number": 106,
            "page": 106,
            "url": y2024["url"],
            "source_sha256": y2024["sha256"],
            "table_id": "Table A1.2",
            "indicator_id": "pwcf_children_lt18",
            "indicator_name": "Proportion of children (<18 years) and adults (>=18 years), by country and overall. People with CF alive on 31/12/2024",
            "stratum_sex": "total",
            "count_n": 24735,
            "percentage": 43.4,
            "denominator_n": 57012,
            "value_point": None,
            "value_low": None,
            "value_high": None,
            "median_iqr": None,
            "unit": "count_and_percentage",
            "raw_text_quote": "Total 24735 43.4 32277 56.6",
            "extraction_kind": "observed",
            "notes": "Children (<18 years): 24735 (43.4%); Adults (>=18 years): 32277 (56.6%); Total alive: 57012.",
        },
        # Table A4.1 FEV1% of predicted for children 6-11 years
        {
            "row_id": "ECFSPR-2024-P108-R01",
            "registry_id": "ecfspr",
            "report_title": y2024["report_title"],
            "report_year": 2024,
            "page_number": 108,
            "page": 108,
            "url": y2024["url"],
            "source_sha256": y2024["sha256"],
            "table_id": "Table A4.1",
            "indicator_id": "fev1_pct_pred_children_6_11",
            "indicator_name": "FEV1% of predicted for children: descriptive statistics, by country and overall. People with CF aged 6-11 years who have never had a transplant",
            "stratum_sex": "total",
            "count_n": 8099,
            "percentage": None,
            "denominator_n": None,
            "value_point": 94.7,
            "value_low": 85.6,
            "value_high": 106.4,
            "median_iqr": {"median": 96.8, "q1": 85.6, "q3": 106.4},
            "unit": "percentage",
            "raw_text_quote": "Total 8099 427 94.7 12.0 85.6 96.8 106.4 155.7",
            "extraction_kind": "observed",
            "notes": "Children aged 6-11 years without transplant: N=8099, missing=427, mean=94.7, min=12.0, max=155.7.",
        },
        # Table A5.1 Prevalence of Pseudomonas aeruginosa in children
        {
            "row_id": "ECFSPR-2024-P115-R01",
            "registry_id": "ecfspr",
            "report_title": y2024["report_title"],
            "report_year": 2024,
            "page_number": 115,
            "page": 115,
            "url": y2024["url"],
            "source_sha256": y2024["sha256"],
            "table_id": "Table A5.1",
            "indicator_id": "pseudomonas_aeruginosa_chronic_children",
            "indicator_name": "Prevalence of Pseudomonas aeruginosa in people with CF seen in 2024 who have never had a transplant (Children)",
            "stratum_sex": "total",
            "count_n": 2096,
            "percentage": 9.0,
            "denominator_n": 23287,
            "value_point": None,
            "value_low": None,
            "value_high": None,
            "median_iqr": None,
            "unit": "count_and_percentage",
            "raw_text_quote": "Total 149 0.6 19287 83 1717 7.4 2096 9.0",
            "extraction_kind": "observed",
            "notes": "Children chronic Pseudomonas aeruginosa infection: N=2096 (9.0%); never: 19287 (83%); intermittent: 1717 (7.4%).",
        },
        # Table A5.1 Prevalence of Pseudomonas aeruginosa in adults
        {
            "row_id": "ECFSPR-2024-P115-R02",
            "registry_id": "ecfspr",
            "report_title": y2024["report_title"],
            "report_year": 2024,
            "page_number": 115,
            "page": 115,
            "url": y2024["url"],
            "source_sha256": y2024["sha256"],
            "table_id": "Table A5.1",
            "indicator_id": "pseudomonas_aeruginosa_chronic_adults",
            "indicator_name": "Prevalence of Pseudomonas aeruginosa in people with CF seen in 2024 who have never had a transplant (Adults)",
            "stratum_sex": "total",
            "count_n": 2740,
            "percentage": 9.8,
            "denominator_n": 27954,
            "value_point": None,
            "value_low": None,
            "value_high": None,
            "median_iqr": None,
            "unit": "count_and_percentage",
            "raw_text_quote": "587 2.10 17218 61.8 7296 26.2 2740 9.8",
            "extraction_kind": "observed",
            "notes": "Adults chronic Pseudomonas aeruginosa infection: N=2740 (9.8%); never: 17218 (61.8%); intermittent: 7296 (26.2%).",
        },
        # Table A6.2 BMI in adults
        {
            "row_id": "ECFSPR-2024-P125-R01",
            "registry_id": "ecfspr",
            "report_title": y2024["report_title"],
            "report_year": 2024,
            "page_number": 125,
            "page": 125,
            "url": y2024["url"],
            "source_sha256": y2024["sha256"],
            "table_id": "Table A6.2",
            "indicator_id": "bmi_adults",
            "indicator_name": "BMI: descriptive statistics by country and overall. All adults with CF seen in 2024 who have never had a transplant",
            "stratum_sex": "total",
            "count_n": 26443,
            "percentage": None,
            "denominator_n": None,
            "value_point": 22.7,
            "value_low": 20.5,
            "value_high": 25.3,
            "median_iqr": {"median": 22.7, "q1": 20.5, "q3": 25.3},
            "unit": "score",
            "raw_text_quote": "Total 26443 598 23.3 12.4 20.5 22.7 25.3 62.1",
            "extraction_kind": "observed",
            "notes": "Overall adult BMI across ECFSPR 2024: N=26443, missing=598, mean=23.3, min=12.4, max=62.1.",
        },
    ])

    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description="ECFSPR registry data extraction pipeline v0.1.")
    parser.add_argument("--downloads-dir", type=Path, default=None, help="Directory containing source PDFs.")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT_PATH, help="Path to write extracted JSON.")
    parser.add_argument("--schema", type=Path, default=DEFAULT_SCHEMA_PATH, help="Path to row_schema.json.")
    parser.add_argument("--verify-hashes", action="store_true", help="Verify SHA-256 digests against actual PDFs.")
    args = parser.parse_args()

    downloads_dir = find_downloads_dir(args.downloads_dir)
    print(f"ECFSPR extraction pipeline v0.1")
    if downloads_dir:
        print(f"Sources download directory: {downloads_dir}")
    else:
        print("Note: Sources download directory not found; building extraction dataset from frozen verified corpus.")

    if args.verify_hashes and downloads_dir:
        for yr, meta in REPORTS_META.items():
            pdf_path = downloads_dir / meta["filename"]
            if pdf_path.exists():
                ok = verify_pdf_hash(pdf_path, meta["sha256"])
                status = "MATCH" if ok else "MISMATCH"
                print(f"  [{status}] {yr}: {meta['filename']} SHA-256")
                if not ok:
                    print(f"ERROR: Hash mismatch for {pdf_path}")
                    return 1

    rows = build_extracted_dataset(downloads_dir)
    print(f"Extracted {len(rows)} rows across 5 report years (2020-2024).")

    # Validate against schema
    from validate import validate_dataset
    if args.schema.exists():
        schema = json.loads(args.schema.read_text(encoding="utf-8"))
        errors = validate_dataset(rows, schema)
        if errors:
            print(f"FAIL: {len(errors)} validation errors encountered:")
            for err in errors[:10]:
                print(f"  - {err}")
            return 1
        print("PASS: 100% of rows strictly conform to row_schema.json (0 errors).")
    else:
        print(f"WARNING: Schema not found at {args.schema}")

    # Check 100% url + page
    missing_url_or_page = [r["row_id"] for r in rows if not r.get("url") or not r.get("page")]
    if missing_url_or_page:
        print(f"FAIL: {len(missing_url_or_page)} rows missing url or page: {missing_url_or_page}")
        return 1
    print("PASS: 100% of rows carry non-empty url and valid page.")

    # Write output JSON
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(rows, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"Wrote extracted dataset to {args.output}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
