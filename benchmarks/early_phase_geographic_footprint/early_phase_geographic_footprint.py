#!/usr/bin/env python3
"""Early-Phase Multi-Country Clinical Trial Site-Network Asymmetry Benchmark (`PRED-0018`).

Tests the pre-registered multi-country site-network hypothesis
(`benchmarks/early_phase_geographic_footprint/done.md`) across the locked 71 early-phase
non-healthy-volunteer-only PRISMA-940 Cystic Fibrosis interventional trials from `PRED-0015`
through `PRED-0017` (`benchmarks/delivery_route_compartmentalization/results_2026-11-24.json`).

Invariants (`cf-research-context`, `discovery-practice`, `security-browsing`, `AH-1`..`AH-4`):
- 100% pure Python 3.10+ standard library + local `tools/trial_atlas/` modules; zero network/LLM calls.
- Strict `AH-1` verbatim NLM attribution, `AH-2` PRISMA-940 cohort boundary, and zero PII/contact fields.
- Omitted/`locations: []`/all-`WITHDRAWN` trials retained in full arm denominators with separate counts.
- Sawtooth-safe upward/downward `delta` bracket scan + multi-peak `p2` golden-section polish for Chan-Zhang 95% CIs.
- CodeQL-safe identifiers (`RESEARCH_BOUNDARY`, `NLM_ATTRIBUTION`, `quarantine_notice`).
"""
from __future__ import annotations

import argparse
import copy
import datetime as dt
from fractions import Fraction
import functools
import hashlib
import json
import math
from pathlib import Path
import statistics
import sys
import tempfile
import unicodedata
from typing import Any

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
PRED0015_DIR = REPO_ROOT / "benchmarks" / "delivery_route_compartmentalization"
TRIAL_ATLAS_DIR = REPO_ROOT / "tools" / "trial_atlas"
for p_dir in (str(PRED0015_DIR), str(TRIAL_ATLAS_DIR)):
    if p_dir not in sys.path:
        sys.path.append(p_dir)

import delivery_route_compartmentalization as p15_mod  # noqa: E402
import scope as scope_mod  # noqa: E402

BENCHMARK_VERSION = "2026-11-25.v1"
EXPECTED_PRED0015_RESULTS_SHA256 = "c4288a7c40783c746e1e4c103e9b1b8867818cb97cc371d19348715d9251dcb1"
EXPECTED_SNAPSHOT_SHA256 = "de0f330c4441d77b84dc5f547e4dfcf99195c1b9f223e4dbeebd628054925b2d"

DEFAULT_DONE_MD = HERE / "done.md"
DEFAULT_PRED0015_RESULTS_JSON = PRED0015_DIR / "results_2026-11-24.json"
DEFAULT_SNAPSHOT_JSON = TRIAL_ATLAS_DIR / "snapshot.json"
DEFAULT_FIXTURE_JSON = HERE / "early_phase_locations_2026-11-25.json"
DEFAULT_RESULTS_JSON = HERE / "results_2026-11-25.json"

NLM_ATTRIBUTION = (
    "ClinicalTrials.gov does not evaluate safety or efficacy; "
    "listing on the registry does not denote government endorsement."
)
RESEARCH_BOUNDARY = (
    "Independent open-source research benchmark only (Tier T1/T2 registry metadata); "
    "not clinical guidance and not an assessment of personal trial eligibility, site availability, or treatment efficacy."
)

ALLOWED_LOCATION_STATUSES: frozenset[str | None] = frozenset(
    {
        "RECRUITING", "NOT_YET_RECRUITING", "ACTIVE_NOT_RECRUITING", "COMPLETED",
        "ENROLLING_BY_INVITATION", "SUSPENDED", "TERMINATED", "WITHDRAWN", "AVAILABLE", "UNKNOWN", None,
    }
)
FORBIDDEN_PII_KEYS: frozenset[str] = frozenset(
    {
        "centralContacts", "overallOfficials", "contacts", "email", "phone",
        "phoneExt", "name", "role", "facility", "city", "state", "zip", "geoPoint",
    }
)
ALLOWED_STUDY_TOP_KEYS: frozenset[str] = frozenset(
    {"nct_id", "arm", "dual_route_modulator_co_admin", "contacts_locations_module_present", "locations_key_present", "locations"}
)
ALLOWED_LOCATION_KEYS: frozenset[str] = frozenset({"country", "status"})

EU27_SOVEREIGN_STATES: frozenset[str] = frozenset(
    {
        "Austria", "Belgium", "Bulgaria", "Croatia", "Cyprus", "Czechia", "Denmark",
        "Estonia", "Finland", "France", "Germany", "Greece", "Hungary", "Ireland",
        "Italy", "Latvia", "Lithuania", "Luxembourg", "Malta", "Netherlands",
        "Poland", "Portugal", "Romania", "Slovakia", "Slovenia", "Spain", "Sweden",
    }
)
NON_EU_ECFSPR_STATES: tuple[str, ...] = ("United Kingdom", "Switzerland", "Israel", "Russia")

_CANONICAL_ISO_SOVEREIGNS: tuple[str, ...] = (
    "Afghanistan", "Albania", "Algeria", "Andorra", "Angola", "Antigua and Barbuda", "Argentina", "Armenia",
    "Australia", "Austria", "Azerbaijan", "Bahamas", "Bahrain", "Bangladesh", "Barbados", "Belarus", "Belgium",
    "Belize", "Benin", "Bhutan", "Bolivia", "Bosnia and Herzegovina", "Botswana", "Brazil", "Brunei", "Bulgaria",
    "Burkina Faso", "Burundi", "Cabo Verde", "Cambodia", "Cameroon", "Canada", "Central African Republic", "Chad",
    "Chile", "China", "Colombia", "Comoros", "Congo", "Costa Rica", "Côte d'Ivoire", "Croatia", "Cuba", "Cyprus",
    "Czechia", "Democratic Republic of the Congo", "Denmark", "Djibouti", "Dominica", "Dominican Republic",
    "Ecuador", "Egypt", "El Salvador", "Equatorial Guinea", "Eritrea", "Estonia", "Eswatini", "Ethiopia", "Fiji",
    "Finland", "France", "Gabon", "Gambia", "Georgia", "Germany", "Ghana", "Greece", "Grenada", "Guatemala",
    "Guinea", "Guinea-Bissau", "Guyana", "Haiti", "Honduras", "Hungary", "Iceland", "India", "Indonesia", "Iran",
    "Iraq", "Ireland", "Israel", "Italy", "Jamaica", "Japan", "Jordan", "Kazakhstan", "Kenya", "Kiribati", "Kosovo",
    "Kuwait", "Kyrgyzstan", "Laos", "Latvia", "Lebanon", "Lesotho", "Liberia", "Libya", "Liechtenstein", "Lithuania",
    "Luxembourg", "Madagascar", "Malawi", "Malaysia", "Maldives", "Mali", "Malta", "Marshall Islands", "Mauritania",
    "Mauritius", "Mexico", "Micronesia", "Moldova", "Monaco", "Mongolia", "Montenegro", "Morocco", "Mozambique",
    "Myanmar", "Namibia", "Nauru", "Nepal", "Netherlands", "New Zealand", "Nicaragua", "Niger", "Nigeria",
    "North Korea", "North Macedonia", "Norway", "Oman", "Pakistan", "Palau", "Palestine", "Panama",
    "Papua New Guinea", "Paraguay", "Peru", "Philippines", "Poland", "Portugal", "Qatar", "Romania", "Russia",
    "Rwanda", "Saint Kitts and Nevis", "Saint Lucia", "Saint Vincent and the Grenadines", "Samoa", "San Marino",
    "Sao Tome and Principe", "Saudi Arabia", "Senegal", "Serbia", "Seychelles", "Sierra Leone", "Singapore",
    "Slovakia", "Slovenia", "Solomon Islands", "Somalia", "South Africa", "South Korea", "South Sudan", "Spain",
    "Sri Lanka", "Sudan", "Suriname", "Sweden", "Switzerland", "Syria", "Taiwan", "Tajikistan", "Tanzania",
    "Thailand", "Timor-Leste", "Togo", "Tonga", "Trinidad and Tobago", "Tunisia", "Türkiye", "Turkmenistan",
    "Tuvalu", "Uganda", "Ukraine", "United Arab Emirates", "United Kingdom", "United States", "Uruguay",
    "Uzbekistan", "Vanuatu", "Vatican City", "Venezuela", "Vietnam", "Yemen", "Zambia", "Zimbabwe",
)

_TERRITORY_AND_RENAME_ALIASES: dict[str, str] = {
    "Puerto Rico": "United States", "Guam": "United States", "U.S. Virgin Islands": "United States",
    "Virgin Islands (U.S.)": "United States", "American Samoa": "United States",
    "Northern Mariana Islands": "United States", "United States of America": "United States",
    "England": "United Kingdom", "Scotland": "United Kingdom", "Wales": "United Kingdom",
    "Northern Ireland": "United Kingdom", "Jersey": "United Kingdom", "Guernsey": "United Kingdom",
    "Isle of Man": "United Kingdom", "Gibraltar": "United Kingdom",
    "Reunion": "France", "Réunion": "France", "Guadeloupe": "France", "Martinique": "France",
    "French Guiana": "France", "New Caledonia": "France", "French Polynesia": "France",
    "Hong Kong": "China", "Macau": "China", "Macao": "China",
    "Faroe Islands": "Denmark", "Greenland": "Denmark",
    "Aruba": "Netherlands", "Curacao": "Netherlands", "Curaçao": "Netherlands",
    "Czech Republic": "Czechia", "Slovak Republic": "Slovakia", "Republic of Ireland": "Ireland",
    "The Netherlands": "Netherlands", "Russian Federation": "Russia",
    "Turkey": "Türkiye", "Turkiye": "Türkiye",
    "Korea, Republic of": "South Korea", "Republic of Korea": "South Korea",
    "Macedonia, The Former Yugoslav Republic of": "North Macedonia",
    "Moldova, Republic of": "Moldova", "Iran, Islamic Republic of": "Iran",
    "Syrian Arab Republic": "Syria", "Venezuela, Bolivarian Republic of": "Venezuela",
    "Bolivia, Plurinational State of": "Bolivia", "Viet Nam": "Vietnam",
}


def _norm_country_key(raw: str) -> str:
    return " ".join(unicodedata.normalize("NFC", raw).split()).casefold()


SOVEREIGN_COUNTRY_MAP: dict[str, str] = {
    **{_norm_country_key(c): c for c in _CANONICAL_ISO_SOVEREIGNS},
    **{_norm_country_key(k): v for k, v in _TERRITORY_AND_RENAME_ALIASES.items()},
}

Z_95 = 1.959963984540054
ALPHA_ONE_SIDED_THRESHOLD = Fraction(1, 100)

_CONFUSABLE_SKELETON_MAP = str.maketrans(
    {
        "\u0430": "a", "\u0435": "e", "\u043e": "o", "\u0440": "p",
        "\u0441": "c", "\u0443": "y", "\u0445": "x", "\u0501": "d",
        "\u0456": "i", "\u0455": "s", "\u0442": "t", "\u043d": "n",
        "\u0410": "A", "\u0412": "B", "\u0421": "C", "\u0500": "D",
        "\u0415": "E", "\u041d": "H", "\u0406": "I", "\u041a": "K",
        "\u041c": "M", "\u041e": "O", "\u0420": "P", "\u0405": "S",
        "\u0422": "T", "\u0425": "X", "\u0423": "Y",
        "\u2215": "/", "\u2044": "/", "\uff0f": "/",
        "\uff1c": "<", "\u27e8": "<", "\u3008": "<", "\u2329": "<", "\u2039": "<", "\ufe64": "<", "\u226e": "<",
        "\uff1e": ">", "\u27e9": ">", "\u3009": ">", "\u232a": ">", "\u203a": ">", "\ufe65": ">", "\u226f": ">",
    }
)


class EarlyPhaseGeographicBenchmarkError(RuntimeError):
    """Raised when a pre-registered schema, PII, provenance, positive-control, or falsification check fails."""


def _is_default_ignorable_or_control(ch: str) -> bool:
    cp = ord(ch)
    cat = unicodedata.category(ch)
    if cat in ("Cf", "Cc") and ch not in ("\n", "\t"):
        return True
    return (
        cp in (0x00AD, 0x034F, 0x061C, 0x115F, 0x1160, 0x17B4, 0x17B5, 0x3164, 0xFFA0)
        or (0x180B <= cp <= 0x180F)
        or (0xFE00 <= cp <= 0xFE0F)
        or (0xFFF0 <= cp <= 0xFFF8)
        or (0x1BCA0 <= cp <= 0x1BCA3)
        or (0x1D173 <= cp <= 0x1D17A)
        or (0xE0000 <= cp <= 0xE0FFF)
    )


def _map_angle_or_char(ch: str) -> str:
    return unicodedata.normalize("NFC", ch).translate(_CONFUSABLE_SKELETON_MAP)


def _validate_2x2_counts(n1: int, k1: int, n2: int, k2: int) -> None:
    for name, val in (("n1", n1), ("k1", k1), ("n2", n2), ("k2", k2)):
        if isinstance(val, bool) or not isinstance(val, int):
            raise ValueError(f"Invalid integer type for {name}: {val!r}")
    if n1 <= 0 or n2 <= 0 or not (0 <= k1 <= n1) or not (0 <= k2 <= n2):
        raise ValueError(f"Invalid 2x2 counts: n1={n1}, k1={k1}, n2={n2}, k2={k2}")


def sha256_lf_file(path: Path) -> str:
    txt = path.read_text(encoding="utf-8").replace("\r\n", "\n")
    return hashlib.sha256(txt.encode("utf-8")).hexdigest()


def _char_skeleton(ch: str) -> str:
    if _is_default_ignorable_or_control(ch):
        return ""
    norm = unicodedata.normalize("NFKD", _map_angle_or_char(ch))
    return "".join(
        sub.lower()
        for sub in norm
        if unicodedata.category(sub) != "Mn" and not _is_default_ignorable_or_control(sub) and not sub.isspace()
    )


def sanitize_registry_text(value: Any) -> str:
    """Strict linear O(N) escape of format/bidi/default-ignorable characters and closing-tag lookalikes."""
    if value is None:
        return ""
    raw = str(value)
    n = len(raw)
    mapped_chars = [_map_angle_or_char(ch) for ch in raw]
    skel_pieces = [_char_skeleton(ch) for ch in raw]

    s_pos = [0] * (n + 1)
    offset = 0
    for i, piece in enumerate(skel_pieces):
        s_pos[i] = offset
        offset += len(piece)
    s_pos[n] = offset
    skel_str = "".join(skel_pieces)
    m = len(skel_str)

    next_data = [m + 1] * (m + 1)
    next_untrusted = [m + 1] * (m + 1)
    last_d = m + 1
    last_u = m + 1
    for p in range(m - 1, -1, -1):
        if skel_str.startswith("data", p):
            last_d = p
        if skel_str.startswith("untrusted", p):
            last_u = p
        next_data[p] = last_d
        next_untrusted[p] = last_u

    next_gt = [-1] * (n + 1)
    run_end_gt = [-1] * (n + 1)
    last_gt = -1
    for i in range(n - 1, -1, -1):
        if mapped_chars[i] == ">":
            last_gt = i
            run_end_gt[i] = run_end_gt[i + 1] if (i + 1 < n and mapped_chars[i + 1] == ">") else i
        next_gt[i] = last_gt

    tag_safe_parts: list[str] = []
    idx = 0
    while idx < n:
        if mapped_chars[idx] == "<" and next_gt[idx] != -1:
            end_idx = run_end_gt[next_gt[idx]]
            a = s_pos[idx]
            b = s_pos[end_idx + 1]
            if (
                b - a >= 6
                and skel_str[a] == "<"
                and skel_str[a + 1] == "/"
                and skel_str[b - 1] == ">"
                and (next_data[a + 2] + 4 <= b - 1 or next_untrusted[a + 2] + 9 <= b - 1)
            ):
                tag_safe_parts.append("[ESCAPED_CLOSING_TAG]")
                idx = end_idx + 1
                continue
        tag_safe_parts.append(raw[idx])
        idx += 1

    inter = "".join(tag_safe_parts)
    out_chars: list[str] = []
    for ch in inter:
        if _is_default_ignorable_or_control(ch):
            out_chars.append(f"<U+{ord(ch):04X}>")
        else:
            out_chars.append(ch)
    return "".join(out_chars)


def assert_zero_at_sign_anywhere(obj: Any, path_label: str = "root") -> None:
    """Recursively verify zero '@' email markers and zero unescaped bidi/control chars in any string."""
    if isinstance(obj, str):
        if "@" in obj:
            raise EarlyPhaseGeographicBenchmarkError(f"{path_label}: '@' email marker detected")
        if any(_is_default_ignorable_or_control(ch) for ch in obj):
            raise EarlyPhaseGeographicBenchmarkError(f"{path_label}: raw bidi/control/default-ignorable character detected")
        if sanitize_registry_text(obj) != obj:
            raise EarlyPhaseGeographicBenchmarkError(f"{path_label}: unescaped closing-tag lookalike detected")
    elif isinstance(obj, dict):
        for k, v in obj.items():
            assert_zero_at_sign_anywhere(k, f"{path_label}.<key>")
            assert_zero_at_sign_anywhere(v, f"{path_label}.{k}")
    elif isinstance(obj, list):
        for i, item in enumerate(obj):
            assert_zero_at_sign_anywhere(item, f"{path_label}[{i}]")


def assert_zero_pii_in_study_record(record: Any, nct_id: str = "synthetic") -> None:
    """Fail closed if any contact/PII/address key, '@', or raw control char is present."""
    if not isinstance(record, dict):
        raise EarlyPhaseGeographicBenchmarkError(f"{nct_id}: study record must be a dict")
    bad_top = set(record.keys()) & FORBIDDEN_PII_KEYS
    if bad_top:
        raise EarlyPhaseGeographicBenchmarkError(f"{nct_id}: forbidden PII/contact keys at study level: {sorted(bad_top)}")
    extra_top = set(record.keys()) - ALLOWED_STUDY_TOP_KEYS
    if extra_top:
        raise EarlyPhaseGeographicBenchmarkError(f"{nct_id}: unexpected study keys outside allowlist: {sorted(extra_top)}")
    assert_zero_at_sign_anywhere(record, path_label=nct_id)
    locs = record.get("locations")
    if locs is not None:
        if not isinstance(locs, list):
            raise EarlyPhaseGeographicBenchmarkError(f"{nct_id}: locations must be list or None")
        for idx, loc in enumerate(locs):
            if not isinstance(loc, dict):
                raise EarlyPhaseGeographicBenchmarkError(f"{nct_id}: location[{idx}] must be a dict")
            bad_loc = set(loc.keys()) & FORBIDDEN_PII_KEYS
            if bad_loc:
                raise EarlyPhaseGeographicBenchmarkError(
                    f"{nct_id}: forbidden PII/facility keys in location[{idx}]: {sorted(bad_loc)}"
                )
            extra_loc = set(loc.keys()) - ALLOWED_LOCATION_KEYS
            if extra_loc:
                raise EarlyPhaseGeographicBenchmarkError(
                    f"{nct_id}: unexpected keys in location[{idx}]: {sorted(extra_loc)}"
                )


def canonicalize_sovereign_country(raw_country: Any, nct_id: str = "synthetic") -> str:
    """Normalize and map a facility country string to its canonical ISO-3166-1 sovereign parent, failing closed."""
    if not isinstance(raw_country, str):
        raise EarlyPhaseGeographicBenchmarkError(f"{nct_id}: facility country must be a non-empty string, got {raw_country!r}")
    if any(_is_default_ignorable_or_control(ch) for ch in raw_country):
        raise EarlyPhaseGeographicBenchmarkError(f"{nct_id}: raw bidi/control character in country {raw_country!r}")
    if sanitize_registry_text(raw_country) != raw_country or "[ESCAPED_CLOSING_TAG]" in raw_country or "<U+" in raw_country:
        raise EarlyPhaseGeographicBenchmarkError(f"{nct_id}: quarantine marker or tag lookalike in country {raw_country!r}")
    key = _norm_country_key(raw_country)
    if not key:
        raise EarlyPhaseGeographicBenchmarkError(f"{nct_id}: blank facility country string {raw_country!r}")
    canonical = SOVEREIGN_COUNTRY_MAP.get(key)
    if canonical is None:
        raise EarlyPhaseGeographicBenchmarkError(
            f"{nct_id}: unmapped country string {raw_country!r} outside locked ISO-3166-1 sovereign allowlist"
        )
    return canonical


def validate_location_status(raw_status: Any, nct_id: str = "synthetic") -> str | None:
    """Validate facility `status` against NLM `RecruitmentStatus` allowlist."""
    if raw_status is not None and not isinstance(raw_status, str):
        raise EarlyPhaseGeographicBenchmarkError(f"{nct_id}: non-string facility status {raw_status!r}")
    if raw_status not in ALLOWED_LOCATION_STATUSES:
        raise EarlyPhaseGeographicBenchmarkError(
            f"{nct_id}: unexpected NLM facility RecruitmentStatus {raw_status!r}"
        )
    return raw_status


def evaluate_locations_list(
    locations: list[tuple[Any, Any]] | list[dict[str, Any]] | None,
    nct_id: str = "synthetic",
) -> dict[str, Any]:
    """Evaluate canonicalized sovereign countries on a list of `{country, status}` dicts or `(country, status)` tuples."""
    if locations is None:
        locations_seq: list[dict[str, Any]] = []
    else:
        locations_seq = []
        for idx, item in enumerate(locations):
            if isinstance(item, tuple) and len(item) == 2:
                locations_seq.append({"country": item[0], "status": item[1]})
            elif isinstance(item, dict):
                locations_seq.append(item)
            else:
                raise EarlyPhaseGeographicBenchmarkError(f"{nct_id}: invalid location entry at {idx}: {item!r}")

    non_withdrawn_countries: set[str] = set()
    all_unfiltered_countries: set[str] = set()
    non_withdrawn_facility_count = 0
    withdrawn_facility_count = 0
    null_status_facility_count = 0

    for loc in locations_seq:
        sov = canonicalize_sovereign_country(loc.get("country"), nct_id=nct_id)
        st = validate_location_status(loc.get("status"), nct_id=nct_id)
        all_unfiltered_countries.add(sov)
        if st == "WITHDRAWN":
            withdrawn_facility_count += 1
        else:
            non_withdrawn_countries.add(sov)
            non_withdrawn_facility_count += 1
            if st is None:
                null_status_facility_count += 1

    sorted_non_withdrawn = sorted(non_withdrawn_countries)
    sorted_unfiltered = sorted(all_unfiltered_countries)
    distinct_n = len(sorted_non_withdrawn)
    return {
        "total_facility_count": len(locations_seq),
        "non_withdrawn_facility_count": non_withdrawn_facility_count,
        "withdrawn_facility_count": withdrawn_facility_count,
        "null_status_non_withdrawn_facility_count": null_status_facility_count,
        "sovereign_countries_non_withdrawn": sorted_non_withdrawn,
        "sovereign_countries_unfiltered": sorted_unfiltered,
        "distinct_sovereign_countries": distinct_n,
        "distinct_sovereign_countries_unfiltered": len(sorted_unfiltered),
        "is_multi_country": distinct_n >= 2,
        "is_multi_country_unfiltered": len(sorted_unfiltered) >= 2,
        "has_at_least_one_non_withdrawn_country": distinct_n >= 1,
    }


def parse_study_locations_record(
    record: dict[str, Any],
    nct_id: str = "synthetic",
    expected_arm: str | None = None,
    expected_dual_route: bool | None = None,
) -> dict[str, Any]:
    """Validate PII guard, study identity/arm parity, flag consistency, and sovereign locations."""
    assert_zero_pii_in_study_record(record, nct_id=nct_id)
    rec_nct = record.get("nct_id")
    if rec_nct is not None and rec_nct != nct_id:
        raise EarlyPhaseGeographicBenchmarkError(f"{nct_id}: record nct_id mismatch {rec_nct!r}")
    if expected_arm is not None and record.get("arm") != expected_arm:
        raise EarlyPhaseGeographicBenchmarkError(f"{nct_id}: record arm {record.get('arm')!r} != {expected_arm!r}")
    if "arm" in record and record.get("arm") not in {"systemic_cftr", "inhaled_cftr"}:
        raise EarlyPhaseGeographicBenchmarkError(f"{nct_id}: invalid arm {record.get('arm')!r}")
    if "dual_route_modulator_co_admin" in record and not isinstance(record.get("dual_route_modulator_co_admin"), bool):
        raise EarlyPhaseGeographicBenchmarkError(f"{nct_id}: dual_route_modulator_co_admin must be bool")
    if expected_dual_route is not None and record.get("dual_route_modulator_co_admin") is not expected_dual_route:
        raise EarlyPhaseGeographicBenchmarkError(f"{nct_id}: dual_route_modulator_co_admin mismatch")

    raw_locs = record.get("locations")
    mod_present = record.get("contacts_locations_module_present")
    if mod_present is not None and not isinstance(mod_present, bool):
        raise EarlyPhaseGeographicBenchmarkError(f"{nct_id}: contacts_locations_module_present must be bool")
    locs_key_present = record.get("locations_key_present")
    if locs_key_present is not None and not isinstance(locs_key_present, bool):
        raise EarlyPhaseGeographicBenchmarkError(f"{nct_id}: locations_key_present must be bool")

    if mod_present is False and (locs_key_present is True or raw_locs is not None):
        raise EarlyPhaseGeographicBenchmarkError(f"{nct_id}: module absent but locations_key_present/locations populated")
    if locs_key_present is False and raw_locs is not None:
        raise EarlyPhaseGeographicBenchmarkError(f"{nct_id}: locations_key_present=False but locations is not None")
    if locs_key_present is True and raw_locs is None:
        raise EarlyPhaseGeographicBenchmarkError(f"{nct_id}: locations_key_present=True but locations is None")

    loc_eval = evaluate_locations_list(raw_locs, nct_id=nct_id)
    omitted_or_empty = raw_locs is None or len(raw_locs) == 0
    all_withdrawn = (not omitted_or_empty) and (loc_eval["distinct_sovereign_countries"] == 0)

    return {
        "nct_id": nct_id,
        "contacts_locations_module_present": bool(mod_present) if mod_present is not None else (raw_locs is not None),
        "locations_key_present": bool(locs_key_present) if locs_key_present is not None else (raw_locs is not None),
        "omitted_or_empty_locations": omitted_or_empty,
        "all_withdrawn_locations": all_withdrawn,
        **loc_eval,
    }


def fisher_hypergeometric_exact(n1: int, k1: int, n2: int, k2: int) -> dict[str, Any]:
    _validate_2x2_counts(n1, k1, n2, k2)
    total_n = n1 + n2
    total_k = k1 + k2
    denom = math.comb(total_n, total_k)
    j_min = max(0, total_k - n2)
    j_max = min(n1, total_k)
    pmf = {j: Fraction(math.comb(n1, j) * math.comb(n2, total_k - j), denom) for j in range(j_min, j_max + 1)}
    support_sum = sum(pmf.values(), Fraction(0, 1))
    if support_sum != Fraction(1, 1):
        raise EarlyPhaseGeographicBenchmarkError(f"Hypergeometric PMF sum != 1: {support_sum}")
    p_one = sum((pmf[j] for j in range(k1, j_max + 1)), Fraction(0, 1))
    p_obs = pmf[k1]
    p_two = sum((prob for prob in pmf.values() if prob <= p_obs), Fraction(0, 1))
    return {
        "p_one_sided_exact": f"{p_one.numerator}/{p_one.denominator}",
        "p_one_sided_float": float(p_one),
        "p_one_sided_fraction": p_one,
        "p_two_sided_exact": f"{p_two.numerator}/{p_two.denominator}",
        "p_two_sided_float": float(p_two),
        "support_sum_exact": "1/1",
    }


def _wilson_single_ci(k: int, n: int, z: float = Z_95) -> tuple[float, float]:
    z2 = z * z
    p_hat = k / n
    denom = 2.0 * (n + z2)
    disc = math.sqrt(z2 + 4.0 * k * (1.0 - p_hat))
    low = 0.0 if k == 0 else (2.0 * k + z2 - z * disc) / denom
    high = 1.0 if k == n else (2.0 * k + z2 + z * disc) / denom
    return low, high


def newcombe_wilson_rd_ci95(n1: int, k1: int, n2: int, k2: int, z: float = Z_95) -> tuple[float, float]:
    _validate_2x2_counts(n1, k1, n2, k2)
    if isinstance(z, bool) or not isinstance(z, (int, float)) or not math.isfinite(z) or z <= 0.0:
        raise ValueError(f"Invalid z-score parameter: {z!r}")
    p1 = k1 / n1
    p2 = k2 / n2
    rd = p1 - p2
    l1, u1 = _wilson_single_ci(k1, n1, z)
    l2, u2 = _wilson_single_ci(k2, n2, z)
    low = rd - math.sqrt((p1 - l1) ** 2 + (u2 - p2) ** 2)
    high = rd + math.sqrt((u1 - p1) ** 2 + (p2 - l2) ** 2)
    return max(-1.0, low), min(1.0, high)


def _mn_restricted_p2(x1: int, n1: int, x2: int, n2: int, delta: float) -> float:
    low_b = max(0.0, -delta)
    high_b = min(1.0, 1.0 - delta)
    if high_b <= low_b + 1e-15:
        return 0.5 * (low_b + high_b)
    total_n = float(n1 + n2)
    total_x = float(x1 + x2)
    l3 = total_n
    l2 = (n1 + 2 * n2) * delta - total_n - total_x
    l1 = (n2 * delta - total_n - 2 * x2) * delta + total_x
    l0 = x2 * delta * (1.0 - delta)
    q = (l2 * l2) / (9.0 * l3 * l3) - l1 / (3.0 * l3)
    r = (l2 ** 3) / (27.0 * (l3 ** 3)) - (l1 * l2) / (6.0 * l3 * l3) + l0 / (2.0 * l3)
    if q <= 1e-30:
        cand = -l2 / (3.0 * l3)
    else:
        u = math.copysign(math.sqrt(q), r if r != 0.0 else 1.0)
        cos_arg = max(-1.0, min(1.0, r / (u ** 3)))
        w = (math.pi + math.acos(cos_arg)) / 3.0
        cand = 2.0 * u * math.cos(w) - l2 / (3.0 * l3)
    return max(low_b, min(high_b, cand))


def _mn_score_z(x1: int, n1: int, x2: int, n2: int, delta: float) -> float:
    num = (x1 / n1) - (x2 / n2) - delta
    if abs(num) < 1e-15:
        return 0.0
    p2_tilde = _mn_restricted_p2(x1, n1, x2, n2, delta)
    p1_tilde = min(1.0, max(0.0, p2_tilde + delta))
    total_n = n1 + n2
    var = (p1_tilde * (1.0 - p1_tilde) / n1 + p2_tilde * (1.0 - p2_tilde) / n2) * (total_n / (total_n - 1.0))
    if var <= 1e-24:
        return math.copysign(1e12, num)
    return num / math.sqrt(var)


def _binom_pmf_vec(n: int, combs: list[int], p: float) -> list[float]:
    if p <= 0.0:
        return [1.0] + [0.0] * n
    if p >= 1.0:
        return [0.0] * n + [1.0]
    q = 1.0 - p
    return [combs[i] * (p ** i) * (q ** (n - i)) for i in range(n + 1)]


def chan_zhang_unconditional_rd_ci95(
    n1: int, k1: int, n2: int, k2: int, alpha_two_sided: float = 0.05
) -> tuple[float, float]:
    """Conservative unconditional Chan-Zhang (1999) 95% CI with multi-peak `p2` polish."""
    _validate_2x2_counts(n1, k1, n2, k2)
    if (
        isinstance(alpha_two_sided, bool)
        or not isinstance(alpha_two_sided, (int, float))
        or not math.isfinite(alpha_two_sided)
        or not (0.0 < alpha_two_sided < 1.0)
    ):
        raise ValueError(f"Invalid alpha_two_sided parameter: {alpha_two_sided!r}")
    return _chan_zhang_cached(int(n1), int(k1), int(n2), int(k2), float(alpha_two_sided))


@functools.lru_cache(maxsize=128)
def _chan_zhang_cached(
    n1: int, k1: int, n2: int, k2: int, alpha_two_sided: float
) -> tuple[float, float]:
    alpha_tail = 0.5 * alpha_two_sided
    comb1 = [math.comb(n1, i) for i in range(n1 + 1)]
    comb2 = [math.comb(n2, j) for j in range(n2 + 1)]
    inv_phi = (math.sqrt(5.0) - 1.0) / 2.0

    def tail_supremum(delta: float, upper: bool) -> float:
        low_p2 = max(0.0, -delta)
        high_p2 = min(1.0, 1.0 - delta)
        z_obs = _mn_score_z(k1, n1, k2, n2, delta)
        tie_tol = 1e-11
        x1_slices: list[tuple[int, int, int]] = []
        pair_total = 0
        for x2 in range(n2 + 1):
            if upper:
                thresh = z_obs - tie_tol
                lo_x, hi_x = 0, n1 + 1
                while lo_x < hi_x:
                    mid_x = (lo_x + hi_x) // 2
                    if _mn_score_z(mid_x, n1, x2, n2, delta) >= thresh:
                        hi_x = mid_x
                    else:
                        lo_x = mid_x + 1
                if lo_x <= n1:
                    x1_slices.append((x2, lo_x, n1 + 1))
                    pair_total += (n1 + 1) - lo_x
            else:
                thresh = z_obs + tie_tol
                lo_x, hi_x = 0, n1 + 1
                while lo_x < hi_x:
                    mid_x = (lo_x + hi_x) // 2
                    if _mn_score_z(mid_x, n1, x2, n2, delta) <= thresh:
                        lo_x = mid_x + 1
                    else:
                        hi_x = mid_x
                if lo_x > 0:
                    x1_slices.append((x2, 0, lo_x))
                    pair_total += lo_x
        if pair_total == 0:
            return 0.0
        if pair_total == (n1 + 1) * (n2 + 1):
            return 1.0

        def prob_at(p2: float) -> float:
            p1 = min(1.0, max(0.0, p2 + delta))
            b1 = _binom_pmf_vec(n1, comb1, p1)
            b2 = _binom_pmf_vec(n2, comb2, p2)
            return sum(b2[x2] * sum(b1[s_lo:s_hi]) for x2, s_lo, s_hi in x1_slices)

        if high_p2 <= low_p2 + 1e-14:
            return prob_at(0.5 * (low_p2 + high_p2))

        grid_n = 240
        step = (high_p2 - low_p2) / grid_n
        vals = [prob_at(low_p2 + step * idx) for idx in range(grid_n + 1)]
        best_val = max(vals)

        peak_indices: list[int] = []
        for idx in range(grid_n + 1):
            if vals[idx] <= 0.0:
                continue
            left_ok = idx == 0 or vals[idx] >= vals[idx - 1]
            right_ok = idx == grid_n or vals[idx] > vals[idx + 1]
            if left_ok and right_ok:
                peak_indices.append(idx)
        peak_indices.sort(key=lambda i: vals[i], reverse=True)
        peak_indices = peak_indices[:6]

        for idx in peak_indices:
            center = low_p2 + step * idx
            a = max(low_p2, center - step)
            b = min(high_p2, center + step)
            for _ in range(36):
                c = b - inv_phi * (b - a)
                d = a + inv_phi * (b - a)
                fc = prob_at(c)
                fd = prob_at(d)
                if fc > best_val:
                    best_val = fc
                if fd > best_val:
                    best_val = fd
                if fc < fd:
                    a = c
                else:
                    b = d
        return min(1.0, max(0.0, best_val))

    rd_obs = (k1 / n1) - (k2 / n2)
    delta_grid = 250

    if k1 == 0 and k2 == n2:
        low_ci = -1.0
    else:
        d_start = -1.0 + 1e-12
        if tail_supremum(d_start, True) > alpha_tail:
            low_ci = -1.0
        elif rd_obs <= d_start + 1e-12:
            low_ci = -1.0
        else:
            step_d = (rd_obs - d_start) / delta_grid
            lo_b = d_start
            hi_b = rd_obs
            for m in range(1, delta_grid + 1):
                cand_d = rd_obs if m == delta_grid else (d_start + step_d * m)
                if tail_supremum(cand_d, True) > alpha_tail:
                    lo_b = d_start + step_d * (m - 1)
                    hi_b = cand_d
                    break
            for _ in range(48):
                mid = 0.5 * (lo_b + hi_b)
                if tail_supremum(mid, True) > alpha_tail:
                    hi_b = mid
                else:
                    lo_b = mid
            low_ci = hi_b

    if k1 == n1 and k2 == 0:
        high_ci = 1.0
    else:
        d_end = 1.0 - 1e-12
        if tail_supremum(d_end, False) > alpha_tail:
            high_ci = 1.0
        elif d_end <= rd_obs + 1e-12:
            high_ci = 1.0
        else:
            step_d = (d_end - rd_obs) / delta_grid
            lo_b = rd_obs
            hi_b = d_end
            for m in range(1, delta_grid + 1):
                cand_d = rd_obs if m == delta_grid else (d_end - step_d * m)
                if tail_supremum(cand_d, False) > alpha_tail:
                    lo_b = cand_d
                    hi_b = d_end - step_d * (m - 1)
                    break
            for _ in range(48):
                mid = 0.5 * (lo_b + hi_b)
                if tail_supremum(mid, False) > alpha_tail:
                    lo_b = mid
                else:
                    hi_b = mid
            high_ci = lo_b

    return low_ci, high_ci


def run_positive_controls(
    pred0015_doc: dict[str, Any],
    snapshot_path: Path = DEFAULT_SNAPSHOT_JSON,
) -> dict[str, Any]:
    """Run all four pre-registered positive controls before reading empirical contactsLocationsModule."""
    cal0 = fisher_hypergeometric_exact(12, 9, 10, 0)
    if cal0["p_one_sided_fraction"] != Fraction(11, 24871):
        raise EarlyPhaseGeographicBenchmarkError(f"Calibration (12,9,10,0) Fisher mismatch: {cal0['p_one_sided_exact']}")
    cz0 = chan_zhang_unconditional_rd_ci95(12, 9, 10, 0)
    nw0 = newcombe_wilson_rd_ci95(12, 9, 10, 0)
    if cz0[0] <= 0.35 or nw0[0] <= 0.35:
        raise EarlyPhaseGeographicBenchmarkError(f"Calibration (12,9,10,0) CI95 <= 0.35: cz={cz0}, nw={nw0}")

    cal_pos = fisher_hypergeometric_exact(54, 35, 15, 3)
    if cal_pos["p_one_sided_fraction"] >= ALPHA_ONE_SIDED_THRESHOLD:
        raise EarlyPhaseGeographicBenchmarkError(
            f"Calibration (54,35,15,3) Fisher p >= 0.01: {cal_pos['p_one_sided_exact']}"
        )
    cz_pos = chan_zhang_unconditional_rd_ci95(54, 35, 15, 3)
    nw_pos = newcombe_wilson_rd_ci95(54, 35, 15, 3)
    if cz_pos[0] <= 0.07 or nw_pos[0] <= 0.15:
        raise EarlyPhaseGeographicBenchmarkError(
            f"Calibration (54,35,15,3) CI95 bounds failed: cz={cz_pos}, nw={nw_pos}"
        )

    snap_doc = json.loads(snapshot_path.read_text(encoding="utf-8"))
    raw_studies = snap_doc.get("studies") or []
    prisma_studies: list[dict[str, Any]] = []
    for st in raw_studies:
        rec_st = p15_mod.snap.study_record(p15_mod._to_raw_study(st))
        rule, _ = scope_mod.classify(rec_st)
        if rule is None:
            p15_mod.assert_flattened_study_schema(st)
            prisma_studies.append(st)
    if len(prisma_studies) != 940:
        raise EarlyPhaseGeographicBenchmarkError(f"Snapshot PRISMA parity failed: prisma={len(prisma_studies)}")

    rederived_sys: list[str] = []
    rederived_inh: list[str] = []
    rederived_dual: list[str] = []
    for st in prisma_studies:
        if not p15_mod.is_early_phase(st.get("phases")):
            continue
        if p15_mod.is_healthy_volunteer_only(st.get("brief_title", ""), st.get("official_title", "")):
            continue
        reg = p15_mod.classify_study_regimen(st)
        nid = str(st["nct_id"])
        if reg["arm"] == "systemic_cftr":
            rederived_sys.append(nid)
        elif reg["arm"] == "inhaled_cftr":
            rederived_inh.append(nid)
            if reg["dual_route_modulator_co_admin"]:
                rederived_dual.append(nid)

    pe15 = pred0015_doc["primary_estimand_early_phase"]
    if (
        sorted(rederived_sys) != sorted(pe15["systemic_ncts"])
        or sorted(rederived_inh) != sorted(pe15["inhaled_ncts"])
        or sorted(rederived_dual) != sorted(pe15["dual_route_inhaled_ncts"])
    ):
        raise EarlyPhaseGeographicBenchmarkError("Locked 71-NCT partition parity mismatch against PRED-0015")

    single_country_fixtures = [
        [("United States", None), ("Puerto Rico", "RECRUITING")],
        [("Czech Republic", None), ("Czechia", "COMPLETED")],
        [("England", None), ("Scotland", None)],
        [("United States", "RECRUITING"), ("Canada", "WITHDRAWN")],
    ]
    for fx in single_country_fixtures:
        ev = evaluate_locations_list(fx)
        if ev["distinct_sovereign_countries"] != 1 or ev["is_multi_country"]:
            raise EarlyPhaseGeographicBenchmarkError(f"Single-country/WITHDRAWN control failed on {fx}: {ev}")

    ev_multi = evaluate_locations_list([("United States", None), ("Canada", None)])
    if ev_multi["distinct_sovereign_countries"] != 2 or not ev_multi["is_multi_country"]:
        raise EarlyPhaseGeographicBenchmarkError(f"Multi-country control failed: {ev_multi}")

    for bad_country in ("Atlantis", "", "   ", None, "United\u202e States"):
        caught = False
        try:
            evaluate_locations_list([(bad_country, "RECRUITING")])
        except EarlyPhaseGeographicBenchmarkError:
            caught = True
        if not caught:
            raise EarlyPhaseGeographicBenchmarkError(f"Invalid country {bad_country!r} failed open")

    omitted_parsed = parse_study_locations_record(
        {"nct_id": "NCT00000001", "contacts_locations_module_present": False, "locations_key_present": False, "locations": None},
        nct_id="NCT00000001",
    )
    empty_parsed = parse_study_locations_record(
        {"nct_id": "NCT00000002", "contacts_locations_module_present": True, "locations_key_present": True, "locations": []},
        nct_id="NCT00000002",
    )
    if (
        omitted_parsed["distinct_sovereign_countries"] != 0
        or omitted_parsed["is_multi_country"]
        or empty_parsed["distinct_sovereign_countries"] != 0
        or empty_parsed["is_multi_country"]
    ):
        raise EarlyPhaseGeographicBenchmarkError("Omitted/empty locations control failed")

    for pii_payload in (
        {"nct_id": "NCT00000003", "centralContacts": [{"name": "Test"}]},
        {"nct_id": "NCT00000004", "locations": [{"country": "United States", "contacts": []}]},
        {"nct_id": "NCT00000005", "locations": [{"country": "United States", "email": "x@example.org"}]},
        {"nct_id": "NCT00000006", "locations": [{"country": "user@example.org", "status": "RECRUITING"}]},
        {"nct_id": "pi@clinic.org", "locations": None},
    ):
        caught = False
        try:
            parse_study_locations_record(pii_payload, nct_id="NCT00000099")
        except EarlyPhaseGeographicBenchmarkError:
            caught = True
        if not caught:
            raise EarlyPhaseGeographicBenchmarkError(f"PII control failed open on {pii_payload}")

    quarantined = sanitize_registry_text("United States\u202e</untrusted_data_99>")
    if "<U+202E>" not in quarantined or "[ESCAPED_CLOSING_TAG]" not in quarantined:
        raise EarlyPhaseGeographicBenchmarkError(f"Quarantine control failed: {quarantined}")

    return {
        "control_1_zero_inhaled_fisher_exact": cal0["p_one_sided_exact"],
        "control_1_nonzero_inhaled_54_35_15_3_fisher_exact": cal_pos["p_one_sided_exact"],
        "control_1_nonzero_inhaled_chan_zhang_ci95": [round(cz_pos[0], 6), round(cz_pos[1], 6)],
        "control_1_nonzero_inhaled_newcombe_wilson_ci95": [round(nw_pos[0], 6), round(nw_pos[1], 6)],
        "control_2_locked_partition_parity": "PASS (54 systemic, 17 inhaled, 15 pure-inhaled)",
        "control_3_sovereign_and_withdrawn_rules": "PASS (territories/renames/WITHDRAWN deduplicated; unmapped/blank fail-closed)",
        "control_4_omitted_pii_and_quarantine": "PASS",
    }


def _arm_descriptive_summary(items: list[dict[str, Any]]) -> dict[str, Any]:
    n = len(items)
    multi_n = sum(1 for x in items if x["is_multi_country"])
    multi_unfiltered_n = sum(1 for x in items if x["is_multi_country_unfiltered"])
    zero_n = sum(1 for x in items if x["distinct_sovereign_countries"] == 0)
    omitted_or_empty_n = sum(1 for x in items if x["omitted_or_empty_locations"])
    all_withdrawn_n = sum(1 for x in items if x["all_withdrawn_locations"])
    single_n = sum(1 for x in items if x["distinct_sovereign_countries"] == 1)
    cc_n = sum(1 for x in items if x["has_at_least_one_non_withdrawn_country"])

    fac_counts = [int(x["non_withdrawn_facility_count"]) for x in items]
    fac_counts_pos = [c for c in fac_counts if c >= 1]
    country_counts = [int(x["distinct_sovereign_countries"]) for x in items]

    us_n = sum(1 for x in items if "United States" in x["sovereign_countries_non_withdrawn"])
    ca_n = sum(1 for x in items if "Canada" in x["sovereign_countries_non_withdrawn"])
    uk_n = sum(1 for x in items if "United Kingdom" in x["sovereign_countries_non_withdrawn"])
    au_n = sum(1 for x in items if "Australia" in x["sovereign_countries_non_withdrawn"])
    eu27_n = sum(1 for x in items if bool(set(x["sovereign_countries_non_withdrawn"]) & EU27_SOVEREIGN_STATES))
    non_eu_ecfspr = {
        st_name: sum(1 for x in items if st_name in x["sovereign_countries_non_withdrawn"])
        for st_name in NON_EU_ECFSPR_STATES
    }
    distinct_union = sorted({c for x in items for c in x["sovereign_countries_non_withdrawn"]})

    return {
        "n": n,
        "multi_country_count": multi_n,
        "multi_country_proportion": round(multi_n / n, 6) if n else 0.0,
        "unfiltered_including_withdrawn_multi_country_count": multi_unfiltered_n,
        "unfiltered_including_withdrawn_multi_country_proportion": round(multi_unfiltered_n / n, 6) if n else 0.0,
        "country_count_breakdown": {
            "zero_countries_total": zero_n,
            "zero_countries_omitted_or_empty_locations": omitted_or_empty_n,
            "zero_countries_all_withdrawn": all_withdrawn_n,
            "single_country_count": single_n,
            "multi_country_ge_2_count": multi_n,
            "max_sovereign_countries_in_single_study": max(country_counts) if country_counts else 0,
        },
        "complete_case_ge_1_country": {
            "n_complete": cc_n,
            "multi_country_count": multi_n,
            "multi_country_proportion": round(multi_n / cc_n, 6) if cc_n else None,
        },
        "non_withdrawn_facility_summary": {
            "median_facilities_all_studies": float(statistics.median(fac_counts)) if fac_counts else 0.0,
            "median_facilities_complete_case": float(statistics.median(fac_counts_pos)) if fac_counts_pos else None,
            "min_facilities": min(fac_counts) if fac_counts else 0,
            "max_facilities": max(fac_counts) if fac_counts else 0,
            "total_non_withdrawn_facilities_sum": sum(fac_counts),
        },
        "jurisdiction_site_presence_non_withdrawn": {
            "United States": {"study_count": us_n, "study_proportion": round(us_n / n, 6) if n else 0.0},
            "Canada": {"study_count": ca_n, "study_proportion": round(ca_n / n, 6) if n else 0.0},
            "United Kingdom": {"study_count": uk_n, "study_proportion": round(uk_n / n, 6) if n else 0.0},
            "Australia": {"study_count": au_n, "study_proportion": round(au_n / n, 6) if n else 0.0},
            "EU-27": {"study_count": eu27_n, "study_proportion": round(eu27_n / n, 6) if n else 0.0},
            "non_eu_ecfspr_states_separate": non_eu_ecfspr,
        },
        "distinct_sovereign_countries_across_arm": distinct_union,
    }


def _complete_case_2x2_summary(sys_items: list[dict[str, Any]], inh_items: list[dict[str, Any]]) -> dict[str, Any]:
    n1 = sum(1 for x in sys_items if x["has_at_least_one_non_withdrawn_country"])
    k1 = sum(1 for x in sys_items if x["is_multi_country"])
    n2 = sum(1 for x in inh_items if x["has_at_least_one_non_withdrawn_country"])
    k2 = sum(1 for x in inh_items if x["is_multi_country"])
    if n1 == 0 or n2 == 0:
        return {"n_sys_cc": n1, "k_sys": k1, "n_inh_cc": n2, "k_inh": k2, "evaluatable": False}
    f_cc = fisher_hypergeometric_exact(n1, k1, n2, k2)
    cz_cc = chan_zhang_unconditional_rd_ci95(n1, k1, n2, k2)
    nw_cc = newcombe_wilson_rd_ci95(n1, k1, n2, k2)
    return {
        "evaluatable": True,
        "n_sys_cc": n1,
        "k_sys": k1,
        "p_sys_cc": round(k1 / n1, 6),
        "n_inh_cc": n2,
        "k_inh": k2,
        "p_inh_cc": round(k2 / n2, 6),
        "risk_difference_cc": round((k1 / n1) - (k2 / n2), 6),
        "p_fisher_one_sided_exact": f_cc["p_one_sided_exact"],
        "p_fisher_one_sided": round(f_cc["p_one_sided_float"], 8),
        "chan_zhang_unconditional_ci95": [round(cz_cc[0], 6), round(cz_cc[1], 6)],
        "newcombe_wilson_ci95": [round(nw_cc[0], 6), round(nw_cc[1], 6)],
    }


def evaluate_benchmark(
    done_path: Path = DEFAULT_DONE_MD,
    pred0015_results_path: Path = DEFAULT_PRED0015_RESULTS_JSON,
    snapshot_path: Path = DEFAULT_SNAPSHOT_JSON,
    fixture_path: Path = DEFAULT_FIXTURE_JSON,
    enforce_canonical_hashes: bool = True,
) -> dict[str, Any]:
    if not done_path.is_file():
        raise EarlyPhaseGeographicBenchmarkError(f"Missing done.md: {done_path}")
    if not pred0015_results_path.is_file():
        raise EarlyPhaseGeographicBenchmarkError(f"Missing PRED-0015 results: {pred0015_results_path}")
    if not fixture_path.is_file():
        raise EarlyPhaseGeographicBenchmarkError(f"Missing frozen locations fixture: {fixture_path}")

    done_sha = sha256_lf_file(done_path)
    p15_sha = sha256_lf_file(pred0015_results_path)
    snap_sha = sha256_lf_file(snapshot_path)
    fixture_sha = sha256_lf_file(fixture_path)

    if enforce_canonical_hashes:
        if p15_sha != EXPECTED_PRED0015_RESULTS_SHA256:
            raise EarlyPhaseGeographicBenchmarkError(f"PRED-0015 SHA-256 mismatch: {p15_sha}")
        if snap_sha != EXPECTED_SNAPSHOT_SHA256:
            raise EarlyPhaseGeographicBenchmarkError(f"snapshot.json SHA-256 mismatch: {snap_sha}")

    pred0015_doc = json.loads(pred0015_results_path.read_text(encoding="utf-8"))
    controls = run_positive_controls(pred0015_doc, snapshot_path=snapshot_path)

    fixture_doc = json.loads(fixture_path.read_text(encoding="utf-8"))
    if not isinstance(fixture_doc, dict):
        raise EarlyPhaseGeographicBenchmarkError("Fixture root must be a JSON object")
    assert_zero_at_sign_anywhere(fixture_doc, path_label="fixture_doc")
    if fixture_doc.get("nlm_attribution") != NLM_ATTRIBUTION:
        raise EarlyPhaseGeographicBenchmarkError("Fixture missing verbatim NLM attribution (AH-1)")
    if fixture_doc.get("research_boundary") != RESEARCH_BOUNDARY:
        raise EarlyPhaseGeographicBenchmarkError("Fixture missing verbatim RESEARCH_BOUNDARY")
    if fixture_doc.get("study_count") != 71 or isinstance(fixture_doc.get("study_count"), bool):
        raise EarlyPhaseGeographicBenchmarkError(f"Fixture study_count {fixture_doc.get('study_count')!r} != 71")

    eval_local_str = str(fixture_doc.get("evaluated_on_local") or "")
    try:
        eval_dt = dt.date.fromisoformat(eval_local_str)
    except ValueError as exc:
        raise EarlyPhaseGeographicBenchmarkError(f"Invalid evaluated_on_local date {eval_local_str!r}") from exc
    if eval_dt > dt.date.today():
        raise EarlyPhaseGeographicBenchmarkError(f"evaluated_on_local {eval_local_str} is in the future")

    pe15 = pred0015_doc["primary_estimand_early_phase"]
    sys_ncts = list(pe15["systemic_ncts"])
    inh_all_ncts = list(pe15["inhaled_ncts"])
    dual_ncts = set(pe15["dual_route_inhaled_ncts"])
    expected_all = set(sys_ncts) | set(inh_all_ncts)

    studies_map = fixture_doc.get("studies")
    if not isinstance(studies_map, dict) or set(studies_map.keys()) != expected_all:
        raise EarlyPhaseGeographicBenchmarkError(
            f"Fixture NCT keys ({len(studies_map) if isinstance(studies_map, dict) else 'invalid'}) != locked 71 NCTs ({len(expected_all)})"
        )

    sys_items = [
        parse_study_locations_record(studies_map[n], nct_id=n, expected_arm="systemic_cftr", expected_dual_route=False)
        for n in sorted(sys_ncts)
    ]
    inh_all_items = [
        parse_study_locations_record(
            studies_map[n], nct_id=n, expected_arm="inhaled_cftr", expected_dual_route=(n in dual_ncts)
        )
        for n in sorted(inh_all_ncts)
    ]
    inh_pure_items = [x for x in inh_all_items if x["nct_id"] not in dual_ncts]

    n_sys, n_inh_pure, n_inh_all = len(sys_items), len(inh_pure_items), len(inh_all_items)
    if (n_sys, n_inh_pure, n_inh_all) != (54, 15, 17):
        raise EarlyPhaseGeographicBenchmarkError(
            f"Cohort margin mismatch: (n_sys={n_sys}, n_inh_pure={n_inh_pure}, n_inh_all={n_inh_all}) != (54, 15, 17)"
        )

    sys_matched = [x["nct_id"] for x in sys_items if x["is_multi_country"]]
    inh_pure_matched = [x["nct_id"] for x in inh_pure_items if x["is_multi_country"]]
    inh_all_matched = [x["nct_id"] for x in inh_all_items if x["is_multi_country"]]

    k_sys, k_inh_pure, k_inh_all = len(sys_matched), len(inh_pure_matched), len(inh_all_matched)
    p_sys, p_inh_pure, p_inh_all = k_sys / n_sys, k_inh_pure / n_inh_pure, k_inh_all / n_inh_all
    rd_pure, rd_all = p_sys - p_inh_pure, p_sys - p_inh_all

    fisher_pure = fisher_hypergeometric_exact(n_sys, k_sys, n_inh_pure, k_inh_pure)
    cz_pure = chan_zhang_unconditional_rd_ci95(n_sys, k_sys, n_inh_pure, k_inh_pure)
    nw_pure = newcombe_wilson_rd_ci95(n_sys, k_sys, n_inh_pure, k_inh_pure)

    fisher_all = fisher_hypergeometric_exact(n_sys, k_sys, n_inh_all, k_inh_all)
    cz_all = chan_zhang_unconditional_rd_ci95(n_sys, k_sys, n_inh_all, k_inh_all)
    nw_all = newcombe_wilson_rd_ci95(n_sys, k_sys, n_inh_all, k_inh_all)

    crit_1_margins = (n_sys == 54) and (n_inh_pure == 15) and (n_inh_all == 17)
    crit_2_directional_pure = p_sys > p_inh_pure
    crit_3_fisher_pure = fisher_pure["p_one_sided_fraction"] < ALPHA_ONE_SIDED_THRESHOLD
    crit_4_ci_bounds_pure = (cz_pure[0] > 0.0) and (nw_pure[0] > 0.0)
    crit_5_all17_preserved = (p_sys > p_inh_all) and (fisher_all["p_one_sided_fraction"] < ALPHA_ONE_SIDED_THRESHOLD)

    is_supported = (
        crit_1_margins and crit_2_directional_pure and crit_3_fisher_pure and crit_4_ci_bounds_pure and crit_5_all17_preserved
    )
    outcome = "supported" if is_supported else "refuted"

    return {
        "schema_version": "pred0018-early-phase-geographic-footprint-v1",
        "benchmark_version": BENCHMARK_VERSION,
        "prediction_id": "PRED-0018",
        "evaluated_on_local": eval_local_str,
        "done_sha256": done_sha,
        "pred0015_results_sha256": p15_sha,
        "snapshot_sha256": snap_sha,
        "fixture_sha256": fixture_sha,
        "nlm_attribution": NLM_ATTRIBUTION,
        "research_boundary": RESEARCH_BOUNDARY,
        "positive_controls": controls,
        "dataset": {
            "catalog_source_id": "clinicaltrials-gov",
            "prisma_included_interventional_studies": 940,
            "locked_early_phase_total": 71,
            "systemic_cftr_count": n_sys,
            "inhaled_cftr_pure_count": n_inh_pure,
            "inhaled_cftr_all_count": n_inh_all,
            "dual_route_inhaled_ncts": sorted(dual_ncts),
        },
        "primary_estimand_pure_route_multi_country": {
            "n_sys": n_sys, "k_sys": k_sys, "p_sys": round(p_sys, 6),
            "n_inh_pure": n_inh_pure, "k_inh_pure": k_inh_pure, "p_inh_pure": round(p_inh_pure, 6),
            "risk_difference": round(rd_pure, 6),
            "p_fisher_one_sided_exact": fisher_pure["p_one_sided_exact"],
            "p_fisher_one_sided": round(fisher_pure["p_one_sided_float"], 8),
            "p_fisher_two_sided_exact": fisher_pure["p_two_sided_exact"],
            "p_fisher_two_sided": round(fisher_pure["p_two_sided_float"], 8),
            "chan_zhang_unconditional_ci95": [round(cz_pure[0], 6), round(cz_pure[1], 6)],
            "newcombe_wilson_ci95": [round(nw_pure[0], 6), round(nw_pure[1], 6)],
            "systemic_matched_ncts": sys_matched,
            "inhaled_pure_matched_ncts": inh_pure_matched,
        },
        "secondary_gate_all_17_inhaled_including_dual_route": {
            "n_sys": n_sys, "k_sys": k_sys, "p_sys": round(p_sys, 6),
            "n_inh_all": n_inh_all, "k_inh_all": k_inh_all, "p_inh_all": round(p_inh_all, 6),
            "risk_difference": round(rd_all, 6),
            "p_fisher_one_sided_exact": fisher_all["p_one_sided_exact"],
            "p_fisher_one_sided": round(fisher_all["p_one_sided_float"], 8),
            "chan_zhang_unconditional_ci95": [round(cz_all[0], 6), round(cz_all[1], 6)],
            "newcombe_wilson_ci95": [round(nw_all[0], 6), round(nw_all[1], 6)],
            "inhaled_all_matched_ncts": inh_all_matched,
        },
        "support_criteria_evaluation": {
            "crit_1_exact_54_15_17_margins": crit_1_margins,
            "crit_2_directional_pure": crit_2_directional_pure,
            "crit_3_fisher_pure_lt_0_01": crit_3_fisher_pure,
            "crit_4_chan_zhang_and_newcombe_ci95_low_gt_0": crit_4_ci_bounds_pure,
            "crit_5_all_17_inhaled_preserved_lt_0_01": crit_5_all17_preserved,
            "outcome": outcome,
        },
        "descriptive_cross_scale_non_gated": {
            "systemic_cftr": _arm_descriptive_summary(sys_items),
            "inhaled_cftr_pure": _arm_descriptive_summary(inh_pure_items),
            "inhaled_cftr_all": _arm_descriptive_summary(inh_all_items),
            "complete_case_sensitivity_pure": _complete_case_2x2_summary(sys_items, inh_pure_items),
            "complete_case_sensitivity_all_17": _complete_case_2x2_summary(sys_items, inh_all_items),
            "engine_a_registry_context_non_pooled": {
                "cffpr_2024_and_2025_us": (
                    "CFFPR 2024: 8.0% ineligible by genotype (29,314 age-and-genotype eligible prescribing denominator); "
                    "2025 Highlights: 2,253 modulator-ineligible of 34,526 in the Registry"
                ),
                "ccfr_2024_canada": "CCFR 2024: 7.6% ineligible by age and genotype",
                "ukcfr_2024_uk": "UKCFR 2024: 12.4% never used modulator among 11,381 registered",
                "acfdr_2025_australia": (
                    "ACFDR 2025: 5.8% adult (of 2,326 non-transplant adults) and 11.5% pediatric "
                    "(170 of 1,473 non-transplant pediatric participants) PBS-ineligible"
                ),
                "ecfspr_2024_eu27_separation": (
                    "EU-27 sovereign site presence is reported separately from non-EU ECFSPR states "
                    "(United Kingdom, Switzerland, Israel, Russia)"
                ),
                "non_pooling_note": (
                    "National registry percentages use distinct regulator lists, transplant filters, and reporting years "
                    "and are juxtaposed strictly as non-pooled epidemiological context without cross-dividing into trial counts."
                ),
            },
        },
    }


def _set_study_locations(st: dict[str, Any], locs: list[dict[str, Any]]) -> None:
    st["contacts_locations_module_present"] = True
    st["locations_key_present"] = True
    st["locations"] = locs


def run_falsification_suite(
    done_path: Path = DEFAULT_DONE_MD,
    pred0015_results_path: Path = DEFAULT_PRED0015_RESULTS_JSON,
    snapshot_path: Path = DEFAULT_SNAPSHOT_JSON,
    fixture_path: Path = DEFAULT_FIXTURE_JSON,
) -> dict[str, Any]:
    """Run 7 planted-defect negative controls verifying deterministic falsification."""
    base_fixture = json.loads(fixture_path.read_text(encoding="utf-8"))
    checks_passed = 0

    with tempfile.TemporaryDirectory() as tmp:
        tmp_dir = Path(tmp)

        eq_doc = copy.deepcopy(base_fixture)
        sys_idx, inh_idx = 0, 0
        for nid in sorted(eq_doc["studies"]):
            st = eq_doc["studies"][nid]
            if st["arm"] == "systemic_cftr":
                active = sys_idx < 27
                sys_idx += 1
            else:
                active = inh_idx < 8
                inh_idx += 1
            _set_study_locations(
                st,
                [{"country": "United States", "status": "COMPLETED"}, {"country": "Canada", "status": "COMPLETED"}]
                if active
                else [{"country": "United States", "status": "COMPLETED"}],
            )
        eq_path = tmp_dir / "eq.json"
        eq_path.write_text(json.dumps(eq_doc), encoding="utf-8")
        if evaluate_benchmark(done_path, pred0015_results_path, snapshot_path, eq_path)["support_criteria_evaluation"]["outcome"] != "refuted":
            raise EarlyPhaseGeographicBenchmarkError("Falsify 1 failed: equalized proportions did not refute")
        checks_passed += 1

        terr_doc = copy.deepcopy(base_fixture)
        for st in terr_doc["studies"].values():
            _set_study_locations(
                st,
                [{"country": "United States", "status": "RECRUITING"}, {"country": "Puerto Rico", "status": "RECRUITING"}],
            )
        terr_path = tmp_dir / "terr.json"
        terr_path.write_text(json.dumps(terr_doc), encoding="utf-8")
        terr_res = evaluate_benchmark(done_path, pred0015_results_path, snapshot_path, terr_path)
        if terr_res["primary_estimand_pure_route_multi_country"]["k_sys"] != 0 or terr_res["support_criteria_evaluation"]["outcome"] != "refuted":
            raise EarlyPhaseGeographicBenchmarkError("Falsify 2 failed: US + Puerto Rico counted as multi-country")
        checks_passed += 1

        wd_doc = copy.deepcopy(base_fixture)
        for st in wd_doc["studies"].values():
            _set_study_locations(
                st,
                [{"country": "United States", "status": "RECRUITING"}, {"country": "Canada", "status": "WITHDRAWN"}],
            )
        wd_path = tmp_dir / "wd.json"
        wd_path.write_text(json.dumps(wd_doc), encoding="utf-8")
        wd_res = evaluate_benchmark(done_path, pred0015_results_path, snapshot_path, wd_path)
        if wd_res["primary_estimand_pure_route_multi_country"]["k_sys"] != 0 or wd_res["support_criteria_evaluation"]["outcome"] != "refuted":
            raise EarlyPhaseGeographicBenchmarkError("Falsify 3 failed: WITHDRAWN foreign site counted as active multi-country")
        checks_passed += 1

        drop_doc = copy.deepcopy(base_fixture)
        first_key = sorted(drop_doc["studies"])[0]
        del drop_doc["studies"][first_key]
        drop_doc["study_count"] = 70
        drop_path = tmp_dir / "drop.json"
        drop_path.write_text(json.dumps(drop_doc), encoding="utf-8")
        caught = False
        try:
            evaluate_benchmark(done_path, pred0015_results_path, snapshot_path, drop_path)
        except EarlyPhaseGeographicBenchmarkError:
            caught = True
        if not caught:
            raise EarlyPhaseGeographicBenchmarkError("Falsify 4 failed: 70/71 incomplete fixture did not fail closed")
        checks_passed += 1

        bad_c_doc = copy.deepcopy(base_fixture)
        _set_study_locations(bad_c_doc["studies"][first_key], [{"country": "Atlantis", "status": "RECRUITING"}])
        bc_path = tmp_dir / "bad_country.json"
        bc_path.write_text(json.dumps(bad_c_doc), encoding="utf-8")
        caught = False
        try:
            evaluate_benchmark(done_path, pred0015_results_path, snapshot_path, bc_path)
        except EarlyPhaseGeographicBenchmarkError:
            caught = True
        if not caught:
            raise EarlyPhaseGeographicBenchmarkError("Falsify 5 failed: unmapped country did not fail closed")
        checks_passed += 1

        pii_doc = copy.deepcopy(base_fixture)
        pii_doc["studies"][first_key]["centralContacts"] = [{"email": "investigator@example.org"}]
        pii_path = tmp_dir / "pii.json"
        pii_path.write_text(json.dumps(pii_doc), encoding="utf-8")
        caught = False
        try:
            evaluate_benchmark(done_path, pred0015_results_path, snapshot_path, pii_path)
        except EarlyPhaseGeographicBenchmarkError:
            caught = True
        if not caught:
            raise EarlyPhaseGeographicBenchmarkError("Falsify 6 failed: planted PII centralContacts did not fail closed")
        checks_passed += 1

        no_nlm_doc = copy.deepcopy(base_fixture)
        no_nlm_doc["nlm_attribution"] = "ClinicalTrials implies endorsement."
        nn_path = tmp_dir / "no_nlm.json"
        nn_path.write_text(json.dumps(no_nlm_doc), encoding="utf-8")
        caught = False
        try:
            evaluate_benchmark(done_path, pred0015_results_path, snapshot_path, nn_path)
        except EarlyPhaseGeographicBenchmarkError:
            caught = True
        if not caught:
            raise EarlyPhaseGeographicBenchmarkError("Falsify 7a failed: tampered NLM attribution did not fail closed")

        border_doc = copy.deepcopy(base_fixture)
        s_count = 0
        for nid in sorted(border_doc["studies"]):
            st = border_doc["studies"][nid]
            if st["arm"] == "systemic_cftr" and s_count < 16:
                _set_study_locations(
                    st, [{"country": "United States", "status": "COMPLETED"}, {"country": "United Kingdom", "status": "COMPLETED"}]
                )
                s_count += 1
            else:
                _set_study_locations(st, [{"country": "United States", "status": "COMPLETED"}])
        bd_path = tmp_dir / "border.json"
        bd_path.write_text(json.dumps(border_doc), encoding="utf-8")
        if evaluate_benchmark(done_path, pred0015_results_path, snapshot_path, bd_path)["support_criteria_evaluation"]["outcome"] != "refuted":
            raise EarlyPhaseGeographicBenchmarkError("Falsify 7b failed: p >= 0.01 borderline table did not refute")
        checks_passed += 1

    return {"status": "PASS", "controls_passed": f"{checks_passed}/7"}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--done", type=Path, default=DEFAULT_DONE_MD)
    ap.add_argument("--pred0015-results", type=Path, default=DEFAULT_PRED0015_RESULTS_JSON)
    ap.add_argument("--snapshot", type=Path, default=DEFAULT_SNAPSHOT_JSON)
    ap.add_argument("--fixture", type=Path, default=DEFAULT_FIXTURE_JSON)
    ap.add_argument("--results", type=Path, default=DEFAULT_RESULTS_JSON)
    ap.add_argument("--write-results", type=Path, default=None)
    ap.add_argument("--controls-only", action="store_true")
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--falsify", action="store_true")
    ap.add_argument("--json", action="store_true")
    args = ap.parse_args(argv)

    if args.controls_only:
        p15_doc = json.loads(args.pred0015_results.read_text(encoding="utf-8"))
        ctrl = run_positive_controls(p15_doc, snapshot_path=args.snapshot)
        print(f"PASS positive_controls: {json.dumps(ctrl)}")
        return 0

    res = evaluate_benchmark(args.done, args.pred0015_results, args.snapshot, args.fixture)
    if args.write_results is not None:
        args.write_results.write_text(json.dumps(res, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    if args.check:
        if not args.results.is_file():
            raise EarlyPhaseGeographicBenchmarkError(f"Missing saved results file for --check: {args.results}")
        saved = json.loads(args.results.read_text(encoding="utf-8"))
        if saved != res:
            raise EarlyPhaseGeographicBenchmarkError(f"Saved {args.results.name} differs from live evaluation")

    falsify_rep = None
    if args.falsify:
        falsify_rep = run_falsification_suite(args.done, args.pred0015_results, args.snapshot, args.fixture)

    if args.json:
        print(json.dumps({"evaluation": res, "falsify": falsify_rep}, indent=2))
    else:
        pe = res["primary_estimand_pure_route_multi_country"]
        sg = res["secondary_gate_all_17_inhaled_including_dual_route"]
        out = res["support_criteria_evaluation"]["outcome"]
        msg = (
            f"PASS early_phase_geographic_footprint: outcome={out} "
            f"sys={pe['k_sys']}/{pe['n_sys']} ({pe['p_sys']:.4f}) "
            f"inh_pure={pe['k_inh_pure']}/{pe['n_inh_pure']} ({pe['p_inh_pure']:.4f}) "
            f"inh_all={sg['k_inh_all']}/{sg['n_inh_all']} ({sg['p_inh_all']:.4f}) "
            f"p_pure_1s={pe['p_fisher_one_sided_exact']} ({pe['p_fisher_one_sided']:.8f}) "
            f"cz_ci95={pe['chan_zhang_unconditional_ci95']} nw_ci95={pe['newcombe_wilson_ci95']}"
        )
        if falsify_rep:
            msg += f" falsify={falsify_rep['controls_passed']}"
        print(msg)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
