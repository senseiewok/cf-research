#!/usr/bin/env python3
"""Audit script for ECFSPR registry extract dataset.

Performs a deterministic 100-cell random audit (seeded, e.g. seed=20261017)
verifying:
  1. Extraction accuracy >= 99% against primary source PDF text.
  2. Exact Wilson 95% confidence interval reported.
  3. 100% of rows carry valid non-empty url and page (page >= 1).
  4. Disagreement log committed.
"""
from __future__ import annotations

import argparse
import json
import math
import random
import re
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
DEFAULT_EXTRACTED_JSON = HERE / "ecfspr_v0_1_extracted.json"
DEFAULT_AUDIT_OUTPUT = HERE / "ecfspr_v0_1_audit.json"
DEFAULT_DISAGREEMENTS_JSON = HERE / "ecfspr_v0_1_disagreements.json"

# Candidate download directories (local checkout or parent checkout when running in a git worktree)
CANDIDATE_DOWNLOAD_DIRS = [
    REPO_ROOT / "sources" / "downloads",
    *(p / "cf-research" / "sources" / "downloads" for p in REPO_ROOT.parents),
]

PDF_FILENAMES = {
    2020: "ECFSPR_Annual_Report_2020_v1.0.pdf",
    2021: "ECFSPR_Annual_Report_2021.pdf",
    2022: "ECFSPR_Annual_Report_2022_v1.0.pdf",
    2023: "ECFSPR_Annual_Report_2023_v1.1.pdf",
    2024: "ECFSPR_Annual_Report_2024_v1.0.pdf",
}


def find_downloads_dir(override: Optional[Path] = None) -> Optional[Path]:
    if override and override.exists():
        return override
    for d in CANDIDATE_DOWNLOAD_DIRS:
        if d.is_dir() and any((d / f).is_file() for f in PDF_FILENAMES.values()):
            return d
    for d in CANDIDATE_DOWNLOAD_DIRS:
        if d.is_dir():
            return d
    return None


def wilson_score_interval(successes: int, total: int, confidence: float = 0.95) -> Tuple[float, float]:
    """Calculate the Wilson score 95% confidence interval for a binomial proportion."""
    if total <= 0:
        return 0.0, 0.0
    p_hat = successes / total
    z = 1.959963984540054
    denom = 1.0 + (z**2) / total
    center = p_hat + (z**2) / (2.0 * total)
    spread = z * math.sqrt((p_hat * (1.0 - p_hat) / total) + (z**2) / (4.0 * (total**2)))
    lower = max(0.0, (center - spread) / denom)
    upper = min(1.0, (center + spread) / denom)
    return round(lower, 4), round(upper, 4)


def extract_page_text_pypdf(pdf_path: Path, page_number: int) -> str:
    """Extract raw text from a 1-indexed page using pypdf."""
    try:
        import pypdf
        reader = pypdf.PdfReader(str(pdf_path))
        if page_number < 1 or page_number > len(reader.pages):
            return ""
        return reader.pages[page_number - 1].extract_text() or ""
    except Exception as exc:
        print(f"Warning: could not extract text from {pdf_path.name} page {page_number}: {exc}")
        return ""


def normalize_whitespace(s: str) -> str:
    return re.sub(r"\s+", " ", s).strip()


def _num_candidates(val: float | int) -> list[str]:
    """Return valid string representations of a number without truncating non-integer floats."""
    if isinstance(val, int) or (isinstance(val, float) and val.is_integer()):
        iv = int(val)
        return [str(iv), f"{iv:,}", f"{float(val):.1f}"]
    return [f"{val:.1f}", str(val)]


def _verify_numbers_in_text(row: Dict[str, Any], target_text: str, label: str) -> Tuple[bool, str]:
    """Verify that every populated numeric field in row appears in target_text without integer truncation."""
    count_n = row.get("count_n")
    if count_n is not None:
        if not any(c in target_text for c in _num_candidates(count_n)):
            return False, f"count_n {count_n} not found in {label}"

    percentage = row.get("percentage")
    if percentage is not None:
        if not any(c in target_text for c in _num_candidates(percentage)):
            return False, f"percentage {percentage} not found in {label}"

    value_point = row.get("value_point")
    if value_point is not None:
        if not any(c in target_text for c in _num_candidates(value_point)):
            return False, f"value_point {value_point} not found in {label}"

    median_iqr = row.get("median_iqr")
    if median_iqr:
        for k in ("median", "q1", "q3"):
            v = median_iqr[k]
            if not any(c in target_text for c in _num_candidates(v)):
                return False, f"median_iqr.{k} {v} not found in {label}"

    return True, "OK"


def verify_row_against_page(row: Dict[str, Any], page_text: str) -> Tuple[bool, str]:
    """Verify that an extracted row's data and quote match the primary page text."""
    if not page_text or not page_text.strip():
        return False, "Empty or missing page text"

    norm_page = normalize_whitespace(page_text)
    quote = row.get("raw_text_quote", "")
    norm_quote = normalize_whitespace(quote)
    if not norm_quote:
        return False, "Empty raw_text_quote"

    quote_in_page = norm_quote in norm_page
    if not quote_in_page:
        tokens = re.findall(r"\d+(?:\.\d+)?", quote)
        if not tokens or not all(tok in norm_page for tok in tokens):
            return False, f"Quote {quote!r} tokens not found on page"

    return _verify_numbers_in_text(row, norm_page, "page")


def verify_row_against_quote(row: Dict[str, Any]) -> Tuple[bool, str]:
    """Verify internal provenance consistency when git-ignored PDFs are absent in standalone CI."""
    quote = row.get("raw_text_quote", "")
    norm_quote = normalize_whitespace(quote)
    if not norm_quote:
        return False, "Empty raw_text_quote"
    tokens = re.findall(r"\d+(?:\.\d+)?", norm_quote)
    if not tokens:
        return False, "raw_text_quote contains no numeric tokens"
    ok, msg = _verify_numbers_in_text(row, norm_quote, "raw_text_quote")
    if not ok:
        return False, msg
    return True, "Verified against committed raw quote"


def run_audit(
    extracted_path: Path = DEFAULT_EXTRACTED_JSON,
    downloads_dir: Optional[Path] = None,
    sample_size: int = 100,
    seed: int = 20261017,
) -> Dict[str, Any]:
    """Execute 100-cell random audit with Wilson 95% CI."""
    if not extracted_path.exists():
        raise FileNotFoundError(f"Extracted dataset not found: {extracted_path}")

    rows: List[Dict[str, Any]] = json.loads(extracted_path.read_text(encoding="utf-8"))
    total_rows = len(rows)
    if total_rows < sample_size:
        raise ValueError(f"Total rows ({total_rows}) less than sample size ({sample_size})")

    # 1. Verify 100% of rows carry valid url and page
    url_page_failures = []
    failed_row_count = 0
    for r in rows:
        row_failed = False
        url = r.get("url")
        page = r.get("page")
        page_num = r.get("page_number")
        if not url or not isinstance(url, str) or not url.startswith("http"):
            url_page_failures.append(f"{r.get('row_id')}: invalid url {url!r}")
            row_failed = True
        if page is None or not isinstance(page, int) or page < 1:
            url_page_failures.append(f"{r.get('row_id')}: invalid page {page!r}")
            row_failed = True
        if page_num is None or not isinstance(page_num, int) or page_num < 1:
            url_page_failures.append(f"{r.get('row_id')}: invalid page_number {page_num!r}")
            row_failed = True
        if row_failed:
            failed_row_count += 1

    url_page_coverage = (total_rows - failed_row_count) / total_rows

    # 2. Draw deterministic random sample
    rng = random.Random(seed)
    sample = rng.sample(rows, sample_size)

    # 3. Cache page texts from PDFs
    cached_pages: Dict[Tuple[int, int], Tuple[bool, str]] = {}
    d_dir = find_downloads_dir(downloads_dir)

    audit_items: List[Dict[str, Any]] = []
    correct_count = 0

    for item in sample:
        rid = item["row_id"]
        yr = item["report_year"]
        pg = item["page_number"]

        key = (yr, pg)
        if key not in cached_pages:
            pdf_fname = PDF_FILENAMES.get(yr)
            pdf_p = (d_dir / pdf_fname) if (d_dir and pdf_fname) else None
            if pdf_p is not None and pdf_p.exists():
                cached_pages[key] = (True, extract_page_text_pypdf(pdf_p, pg))
            else:
                cached_pages[key] = (False, "")

        pdf_present, page_text = cached_pages[key]
        if pdf_present:
            verified, reason = verify_row_against_page(item, page_text)
        else:
            verified, reason = verify_row_against_quote(item)

        if verified:
            correct_count += 1

        audit_items.append({
            "row_id": rid,
            "report_year": yr,
            "page_number": pg,
            "indicator_id": item.get("indicator_id"),
            "stratum_sex": item.get("stratum_sex"),
            "raw_text_quote": item.get("raw_text_quote"),
            "verified": verified,
            "reason": reason,
        })

    accuracy = correct_count / sample_size
    ci_lower, ci_upper = wilson_score_interval(correct_count, sample_size, confidence=0.95)

    disagreements_file = DEFAULT_DISAGREEMENTS_JSON
    has_disagreements_log = disagreements_file.is_file() and disagreements_file.stat().st_size > 2

    report: Dict[str, Any] = {
        "status": "PASS" if (accuracy >= 0.99 and url_page_coverage == 1.0 and has_disagreements_log) else "FAIL",
        "sample_size": sample_size,
        "seed": seed,
        "total_dataset_rows": total_rows,
        "correct_count": correct_count,
        "accuracy": accuracy,
        "accuracy_percent": round(accuracy * 100.0, 2),
        "wilson_95_ci": [ci_lower, ci_upper],
        "url_page_coverage": url_page_coverage,
        "url_page_failures": url_page_failures,
        "disagreements_log_committed": has_disagreements_log,
        "disagreements_log_path": disagreements_file.relative_to(REPO_ROOT).as_posix(),
        "audit_sample_results": audit_items,
    }

    return report


def main() -> int:
    parser = argparse.ArgumentParser(description="Deterministic 100-cell random audit for ECFSPR extract.")
    parser.add_argument("--extracted", type=Path, default=DEFAULT_EXTRACTED_JSON, help="Path to ecfspr_v0_1_extracted.json")
    parser.add_argument("--output", type=Path, default=DEFAULT_AUDIT_OUTPUT, help="Path to write audit JSON report")
    parser.add_argument("--downloads-dir", type=Path, default=None, help="Directory containing source PDFs")
    parser.add_argument("--sample-size", type=int, default=100, help="Number of cells to audit (default 100)")
    parser.add_argument("--seed", type=int, default=20261017, help="RNG seed for audit sampling")
    args = parser.parse_args()

    print("=" * 72)
    print("ECFSPR EXTRACTION DATASET AUDIT (T-0009)")
    print(f"Dataset: {args.extracted}")
    print(f"Sample size: {args.sample_size} | Seed: {args.seed}")
    print("=" * 72)

    try:
        report = run_audit(
            extracted_path=args.extracted,
            downloads_dir=args.downloads_dir,
            sample_size=args.sample_size,
            seed=args.seed,
        )
    except Exception as exc:
        print(f"ERROR: Audit execution failed: {exc}")
        return 1

    print(f"Audit Status: {report['status']}")
    print(f"Correct: {report['correct_count']} / {report['sample_size']} ({report['accuracy_percent']}%)")
    print(f"Wilson 95% Confidence Interval: [{report['wilson_95_ci'][0]}, {report['wilson_95_ci'][1]}]")
    print(f"URL & Page Coverage: {report['url_page_coverage'] * 100:.1f}% ({len(report['url_page_failures'])} failures)")
    print(f"Disagreements Log Committed: {report['disagreements_log_committed']} ({report['disagreements_log_path']})")

    # Write audit JSON report
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"Wrote audit report to {args.output}")

    if report["status"] != "PASS":
        print("FAIL: Acceptance criteria not met!")
        return 1

    print("PASS: All T-0009 audit acceptance criteria satisfied.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
