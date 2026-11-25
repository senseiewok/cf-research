#!/usr/bin/env python3
"""Independent zero-`re`, zero-`math`, zero-`float` audit of `PRED-0018`.

Enforces Rule 2 of `discovery-practice` and `done.md` Section 4:
- Imports ZERO `re`, ZERO `math`, ZERO `importlib`, ZERO `tools/`, and ZERO benchmark helper modules.
- Uses ZERO IEEE `float()` calls and ZERO `float` literals anywhere in the module (`28`-digit `Decimal` bracket scan refined at `50`-digit `Decimal` precision on final bisection steps).
- Re-checks `done_sha256`, `pred0015_results_sha256`, `snapshot_sha256`, `fixture_sha256`, verbatim NLM attribution & research boundary, zero PII/contact/`@`/raw-control markers, `ledger.json` `checked_on`, and `54/17/15` NCT set equality.
- Independently canonicalizes raw `locations[].country` strings, filters `status == "WITHDRAWN"`, and audits primary, secondary, complete-case sensitivity, facility medians/sums, country breakdowns, and jurisdiction counts.
"""
from __future__ import annotations

import argparse
import ast
import datetime as dt
from decimal import Decimal, getcontext, localcontext
from fractions import Fraction
import functools
import hashlib
import json
from pathlib import Path
import sys
import unicodedata
from typing import Any

getcontext().prec = 50

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
DEFAULT_DONE = HERE / "done.md"
DEFAULT_LEDGER = HERE / "ledger.json"
DEFAULT_PRED0015_RESULTS = REPO_ROOT / "benchmarks" / "delivery_route_compartmentalization" / "results_2026-11-24.json"
DEFAULT_SNAPSHOT = REPO_ROOT / "tools" / "trial_atlas" / "snapshot.json"
DEFAULT_FIXTURE = HERE / "early_phase_locations_2026-11-25.json"
DEFAULT_RESULTS = HERE / "results_2026-11-25.json"

EXPECTED_PRED0015_RESULTS_SHA256 = "c4288a7c40783c746e1e4c103e9b1b8867818cb97cc371d19348715d9251dcb1"
EXPECTED_SNAPSHOT_SHA256 = "de0f330c4441d77b84dc5f547e4dfcf99195c1b9f223e4dbeebd628054925b2d"

EXPECTED_NLM_ATTRIBUTION = (
    "ClinicalTrials.gov does not evaluate safety or efficacy; "
    "listing on the registry does not denote government endorsement."
)
EXPECTED_RESEARCH_BOUNDARY = (
    "Independent open-source research benchmark only (Tier T1/T2 registry metadata); "
    "not clinical guidance and not an assessment of personal trial eligibility, site availability, or treatment efficacy."
)

VALID_STATUSES = frozenset(
    {
        "RECRUITING", "NOT_YET_RECRUITING", "ACTIVE_NOT_RECRUITING", "COMPLETED",
        "ENROLLING_BY_INVITATION", "SUSPENDED", "TERMINATED", "WITHDRAWN", "AVAILABLE", "UNKNOWN", None,
    }
)
FORBIDDEN_PII = frozenset(
    {
        "centralContacts", "overallOfficials", "contacts", "email", "phone",
        "phoneExt", "name", "role", "facility", "city", "state", "zip", "geoPoint",
    }
)
ALLOWED_STUDY_KEYS = frozenset(
    {"nct_id", "arm", "dual_route_modulator_co_admin", "contacts_locations_module_present", "locations_key_present", "locations"}
)
ALLOWED_LOC_KEYS = frozenset({"country", "status"})

EU27_SET = frozenset(
    {
        "Austria", "Belgium", "Bulgaria", "Croatia", "Cyprus", "Czechia", "Denmark",
        "Estonia", "Finland", "France", "Germany", "Greece", "Hungary", "Ireland",
        "Italy", "Latvia", "Lithuania", "Luxembourg", "Malta", "Netherlands",
        "Poland", "Portugal", "Romania", "Slovakia", "Slovenia", "Spain", "Sweden",
    }
)

_SOVEREIGNS_AUDIT: tuple[str, ...] = (
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

_ALIASES_AUDIT: dict[str, str] = {
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


def _audit_country_lookup() -> dict[str, str]:
    table: dict[str, str] = {}
    for s in _SOVEREIGNS_AUDIT:
        table[" ".join(unicodedata.normalize("NFC", s).split()).casefold()] = s
    for raw_k, target_v in _ALIASES_AUDIT.items():
        table[" ".join(unicodedata.normalize("NFC", raw_k).split()).casefold()] = target_v
    return table


AUDIT_COUNTRY_TABLE = _audit_country_lookup()


def _is_ignorable_ctrl(ch: str) -> bool:
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


_AUDIT_TAG_CONFUSABLES = str.maketrans(
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


def _has_closing_data_tag_lookalike(text: str) -> bool:
    pieces: list[str] = []
    mapped = [unicodedata.normalize("NFC", ch).translate(_AUDIT_TAG_CONFUSABLES) for ch in text]
    for ch, mch in zip(text, mapped):
        if _is_ignorable_ctrl(ch):
            pieces.append("")
            continue
        norm = unicodedata.normalize("NFKD", mch)
        pieces.append(
            "".join(
                sub.lower()
                for sub in norm
                if unicodedata.category(sub) != "Mn" and not _is_ignorable_ctrl(sub) and not sub.isspace()
            )
        )
    n = len(text)
    s_pos = [0] * (n + 1)
    off = 0
    for i, pc in enumerate(pieces):
        s_pos[i] = off
        off += len(pc)
    s_pos[n] = off
    skel = "".join(pieces)
    m = len(skel)
    if m < 6:
        return False
    next_d = [m + 1] * (m + 1)
    next_u = [m + 1] * (m + 1)
    ld, lu = m + 1, m + 1
    for p in range(m - 1, -1, -1):
        if skel.startswith("data", p):
            ld = p
        if skel.startswith("untrusted", p):
            lu = p
        next_d[p] = ld
        next_u[p] = lu
    next_gt = [-1] * (n + 1)
    run_end = [-1] * (n + 1)
    lg = -1
    for i in range(n - 1, -1, -1):
        if mapped[i] == ">":
            lg = i
            run_end[i] = run_end[i + 1] if (i + 1 < n and mapped[i + 1] == ">") else i
        next_gt[i] = lg
    for idx in range(n):
        if mapped[idx] == "<" and next_gt[idx] != -1:
            end_idx = run_end[next_gt[idx]]
            a, b = s_pos[idx], s_pos[end_idx + 1]
            if (
                b - a >= 6
                and skel[a] == "<"
                and skel[a + 1] == "/"
                and skel[b - 1] == ">"
                and (next_d[a + 2] + 4 <= b - 1 or next_u[a + 2] + 9 <= b - 1)
            ):
                return True
    return False


def _scan_strings_for_pii_or_ctrl(obj: Any, label: str, mismatches: list[str]) -> None:
    if isinstance(obj, str):
        if "@" in obj:
            mismatches.append(f"{label}: '@' email marker detected")
        if any(_is_ignorable_ctrl(ch) for ch in obj):
            mismatches.append(f"{label}: raw bidi/control character detected")
        if "<U+" in obj or "[ESCAPED_CLOSING_TAG]" in obj:
            mismatches.append(f"{label}: quarantine escape token detected")
        if _has_closing_data_tag_lookalike(obj):
            mismatches.append(f"{label}: unescaped closing-tag lookalike detected")
    elif isinstance(obj, dict):
        for k, v in obj.items():
            _scan_strings_for_pii_or_ctrl(k, f"{label}.<key>", mismatches)
            _scan_strings_for_pii_or_ctrl(v, f"{label}.{k}", mismatches)
    elif isinstance(obj, list):
        for idx, elem in enumerate(obj):
            _scan_strings_for_pii_or_ctrl(elem, f"{label}[{idx}]", mismatches)


def _validate_2x2_counts(n1: int, k1: int, n2: int, k2: int) -> None:
    for name, val in (("n1", n1), ("k1", k1), ("n2", n2), ("k2", k2)):
        if isinstance(val, bool) or not isinstance(val, int):
            raise ValueError(f"Invalid integer type for {name}: {val!r}")
    if n1 <= 0 or n2 <= 0 or not (0 <= k1 <= n1) or not (0 <= k2 <= n2):
        raise ValueError(f"Invalid 2x2 counts: n1={n1}, k1={k1}, n2={n2}, k2={k2}")


def _sha256_lf(path: Path) -> str:
    txt = path.read_text(encoding="utf-8").replace("\r\n", "\n")
    return hashlib.sha256(txt.encode("utf-8")).hexdigest()


def verify_ast_independence(
    self_path: Path = Path(__file__),
    source_override: str | None = None,
) -> None:
    """Module-wide AST guard forbidding `re`, `math`, `importlib`, `tools`, `benchmarks`, dynamic imports, `AsyncFunctionDef`, `float(...)` calls, and `float` literals."""
    src = source_override if source_override is not None else self_path.read_text(encoding="utf-8")
    tree = ast.parse(src, filename=str(self_path))
    banned = {
        "re", "math", "importlib", "builtins", "tools", "benchmarks",
        "early_phase_geographic_footprint", "acquire_locations_fixture",
        "delivery_route_compartmentalization", "scope", "snapshot",
    }
    dyn_builtin = "_" + "_import_" + "_"
    dyn_mod = "import_" + "module"
    builtins_dunder = "_" + "_builtins_" + "_"
    for node in ast.walk(tree):
        if isinstance(node, ast.AsyncFunctionDef):
            raise RuntimeError(f"Forbidden AsyncFunctionDef in independent_audit.py: {node.name}")
        if isinstance(node, ast.Constant):
            if isinstance(node.value, float):
                raise RuntimeError(f"Forbidden IEEE float literal in independent_audit.py: {node.value!r}")
            if isinstance(node.value, str) and node.value in {dyn_builtin, dyn_mod, builtins_dunder}:
                if node.value != "_":
                    raise RuntimeError(f"Forbidden dynamic-import string token: {node.value!r}")
        if isinstance(node, ast.Name) and node.id in {dyn_builtin, builtins_dunder}:
            raise RuntimeError(f"Forbidden dynamic/builtins identifier: {node.id}")
        if isinstance(node, ast.Import):
            for alias in node.names:
                if set(alias.name.split(".")) & banned:
                    raise RuntimeError(f"Forbidden import in independent_audit.py: {alias.name}")
        elif isinstance(node, ast.ImportFrom):
            if (node.level or 0) > 0:
                raise RuntimeError("Forbidden relative import in independent_audit.py")
            mod_parts = set(node.module.split(".")) if node.module else set()
            alias_parts = {alias.name.split(".")[0] for alias in node.names}
            if (mod_parts | alias_parts) & banned:
                raise RuntimeError(f"Forbidden from-import in independent_audit.py: {node.module}")
        elif isinstance(node, ast.Call):
            if isinstance(node.func, ast.Name) and node.func.id in {dyn_builtin, "float", "eval", "exec", "getattr"}:
                raise RuntimeError(f"Forbidden call in independent_audit.py: {node.func.id}()")
            if isinstance(node.func, ast.Attribute) and node.func.attr in {dyn_mod, dyn_builtin}:
                raise RuntimeError(f"Forbidden attribute call in independent_audit.py: {node.func.attr}()")


def _binom_fraction(n: int, k: int) -> Fraction:
    if k < 0 or k > n:
        return Fraction(0, 1)
    m = min(k, n - k)
    acc = Fraction(1, 1)
    for r in range(1, m + 1):
        acc *= Fraction(n - m + r, r)
    return acc


def _fisher_tails_multiplicative(n1: int, k1: int, n2: int, k2: int) -> tuple[Fraction, Fraction]:
    _validate_2x2_counts(n1, k1, n2, k2)
    total_n = n1 + n2
    total_k = k1 + k2
    denom = _binom_fraction(total_n, total_k)
    j_min = max(0, total_k - n2)
    j_max = min(n1, total_k)
    pmf = {j: (_binom_fraction(n1, j) * _binom_fraction(n2, total_k - j)) / denom for j in range(j_min, j_max + 1)}
    if sum(pmf.values(), Fraction(0, 1)) != Fraction(1, 1):
        raise ValueError("Hypergeometric mass sum != 1")
    p_one = sum((pmf[j] for j in range(k1, j_max + 1)), Fraction(0, 1))
    p_obs = pmf[k1]
    p_two = sum((prob for prob in pmf.values() if prob <= p_obs), Fraction(0, 1))
    return p_one, p_two


def _median_decimal(vals: list[int]) -> Decimal | None:
    if not vals:
        return None
    s = sorted(vals)
    m = len(s)
    mid = m // 2
    if m % 2 == 1:
        return Decimal(s[mid])
    return (Decimal(s[mid - 1]) + Decimal(s[mid])) / Decimal(2)


def _newcombe_wilson_decimal(n1: int, k1: int, n2: int, k2: int) -> tuple[Decimal, Decimal]:
    _validate_2x2_counts(n1, k1, n2, k2)
    z = Decimal("1.959963984540054")
    z2 = z * z

    def wilson(k: int, n: int) -> tuple[Decimal, Decimal]:
        nd = Decimal(n)
        kd = Decimal(k)
        p = kd / nd
        denom = Decimal(2) * (nd + z2)
        disc = (z2 + Decimal(4) * kd * (Decimal(1) - p)).sqrt()
        lo = Decimal(0) if k == 0 else (Decimal(2) * kd + z2 - z * disc) / denom
        hi = Decimal(1) if k == n else (Decimal(2) * kd + z2 + z * disc) / denom
        return lo, hi

    p1 = Decimal(k1) / Decimal(n1)
    p2 = Decimal(k2) / Decimal(n2)
    rd = p1 - p2
    l1, u1 = wilson(k1, n1)
    l2, u2 = wilson(k2, n2)
    low = max(Decimal(-1), rd - ((p1 - l1) ** 2 + (u2 - p2) ** 2).sqrt())
    high = min(Decimal(1), rd + ((u1 - p1) ** 2 + (p2 - l2) ** 2).sqrt())
    return low, high


def _newton_bisect_mn_p2_decimal(x1: int, n1: int, x2: int, n2: int, delta: Decimal, prec_digits: int = 28) -> Decimal:
    """Safeguarded Newton-bisection on constrained binomial log-likelihood derivative (zero cubic expansion)."""
    zero = Decimal(0)
    one = Decimal(1)
    half = Decimal("0.5")
    rel_eps = Decimal("1e-32") if prec_digits >= 50 else Decimal("1e-12")
    floor_p = Decimal("1e-46") if prec_digits >= 50 else Decimal("1e-24")

    lo = max(zero, -delta)
    hi = min(one, one - delta)
    if hi <= lo + Decimal(4) * floor_p:
        return half * (lo + hi)

    dx1, dy1 = Decimal(x1), Decimal(n1 - x1)
    dx2, dy2 = Decimal(x2), Decimal(n2 - x2)

    if delta >= zero and x2 == 0:
        if delta == zero and x1 == 0:
            return zero
        if delta > zero and (dx1 / delta - dy1 / (one - delta) - dy2) <= zero:
            return lo
    elif delta < zero and x1 == 0:
        neg_d = -delta
        if (-dy1 + dx2 / neg_d - dy2 / (one - neg_d)) <= zero:
            return lo

    if delta <= zero and x2 == n2:
        if delta == zero and x1 == n1:
            return one
        if delta < zero and (dx1 / (one + delta) - dy1 / (-delta) + dx2) >= zero:
            return hi
    elif delta > zero and x1 == n1:
        omd = one - delta
        if (dx1 + dx2 / omd - dy2 / delta) >= zero:
            return hi

    margin = max((hi - lo) * rel_eps, floor_p)
    a = lo + margin
    b = hi - margin
    if b <= a:
        return half * (lo + hi)

    def grad_ell(p2: Decimal) -> Decimal:
        p2_c = min(one - floor_p, max(floor_p, p2))
        p1_c = min(one - floor_p, max(floor_p, p2 + delta))
        t1 = (dx1 / p1_c if x1 > 0 else zero) - (dy1 / (one - p1_c) if x1 < n1 else zero)
        t2 = (dx2 / p2_c if x2 > 0 else zero) - (dy2 / (one - p2_c) if x2 < n2 else zero)
        return t1 + t2

    def dgrad_ell(p2: Decimal) -> Decimal:
        p2_c = min(one - floor_p, max(floor_p, p2))
        p1_c = min(one - floor_p, max(floor_p, p2 + delta))
        d1 = -(dx1 / (p1_c * p1_c) if x1 > 0 else zero) - (dy1 / ((one - p1_c) * (one - p1_c)) if x1 < n1 else zero)
        d2 = -(dx2 / (p2_c * p2_c) if x2 > 0 else zero) - (dy2 / ((one - p2_c) * (one - p2_c)) if x2 < n2 else zero)
        return d1 + d2

    if grad_ell(a) <= zero:
        return lo
    if grad_ell(b) >= zero:
        return hi

    x = half * (a + b)
    tol = Decimal("1e-42") if prec_digits >= 50 else Decimal("1e-18")
    max_iters = 64 if prec_digits >= 50 else 36
    for _ in range(max_iters):
        fx = grad_ell(x)
        if abs(fx) <= tol or (b - a) <= tol:
            return x
        if fx > zero:
            a = x
        else:
            b = x
        dfx = dgrad_ell(x)
        if dfx < -floor_p:
            cand = x - fx / dfx
            if a < cand < b:
                x = cand
                continue
        x = half * (a + b)
    return half * (a + b)


def _mn_z_decimal(x1: int, n1: int, x2: int, n2: int, delta: Decimal, prec_digits: int = 28) -> Decimal:
    zero = Decimal(0)
    one = Decimal(1)
    num_tol = Decimal("1e-46") if prec_digits >= 50 else Decimal("1e-24")
    num = (Decimal(x1) / Decimal(n1)) - (Decimal(x2) / Decimal(n2)) - delta
    if abs(num) < num_tol:
        return zero
    p2_t = _newton_bisect_mn_p2_decimal(x1, n1, x2, n2, delta, prec_digits=prec_digits)
    p1_t = min(one, max(zero, p2_t + delta))
    tot = Decimal(n1 + n2)
    var = (p1_t * (one - p1_t) / Decimal(n1) + p2_t * (one - p2_t) / Decimal(n2)) * (tot / (tot - one))
    var_floor = Decimal("1e-90") if prec_digits >= 50 else Decimal("1e-36")
    if var <= var_floor:
        return Decimal("1e18") if num > zero else Decimal("-1e18")
    return num / var.sqrt()


def _chan_zhang_ci95_decimal(n1: int, k1: int, n2: int, k2: int) -> tuple[Decimal, Decimal]:
    """100% pure `Decimal` (zero `float` pre-screen; 28-digit scan refined at 50-digit precision) Chan-Zhang 95% CI."""
    _validate_2x2_counts(n1, k1, n2, k2)
    return _chan_zhang_ci95_decimal_inner(int(n1), int(k1), int(n2), int(k2))


@functools.lru_cache(maxsize=64)
def _chan_zhang_ci95_decimal_inner(n1: int, k1: int, n2: int, k2: int) -> tuple[Decimal, Decimal]:
    zero = Decimal(0)
    one = Decimal(1)
    half = Decimal("0.5")
    third = one / Decimal(3)
    alpha_tail = Decimal("0.025")
    tie_tol = Decimal("1e-11")

    c1 = [Decimal(int(_binom_fraction(n1, i))) for i in range(n1 + 1)]
    c2 = [Decimal(int(_binom_fraction(n2, j))) for j in range(n2 + 1)]

    def tail_sup(delta: Decimal, upper: bool, prec_digits: int = 28) -> Decimal:
        with localcontext() as ctx:
            ctx.prec = prec_digits
            lo_p2 = max(zero, -delta)
            hi_p2 = min(one, one - delta)
            z_obs = _mn_z_decimal(k1, n1, k2, n2, delta, prec_digits=prec_digits)
            x1_slices: list[tuple[int, int, int]] = []
            pair_total = 0
            for x2 in range(n2 + 1):
                if upper:
                    thresh = z_obs - tie_tol
                    lo_x, hi_x = 0, n1 + 1
                    while lo_x < hi_x:
                        mid_x = (lo_x + hi_x) // 2
                        if _mn_z_decimal(mid_x, n1, x2, n2, delta, prec_digits=prec_digits) >= thresh:
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
                        if _mn_z_decimal(mid_x, n1, x2, n2, delta, prec_digits=prec_digits) <= thresh:
                            lo_x = mid_x + 1
                        else:
                            hi_x = mid_x
                    if lo_x > 0:
                        x1_slices.append((x2, 0, lo_x))
                        pair_total += lo_x
            if pair_total == 0:
                return zero
            if pair_total == (n1 + 1) * (n2 + 1):
                return one

            def _binom_vec_dec(n_sz: int, combs_d: list[Decimal], p_val: Decimal) -> list[Decimal]:
                if p_val <= zero:
                    return [one] + [zero] * n_sz
                if p_val >= one:
                    return [zero] * n_sz + [one]
                q_val = one - p_val
                p_pows = [one] * (n_sz + 1)
                q_pows = [one] * (n_sz + 1)
                for idx_p in range(1, n_sz + 1):
                    p_pows[idx_p] = p_pows[idx_p - 1] * p_val
                    q_pows[idx_p] = q_pows[idx_p - 1] * q_val
                return [combs_d[idx_p] * p_pows[idx_p] * q_pows[n_sz - idx_p] for idx_p in range(n_sz + 1)]

            def prob_dec(p2: Decimal) -> Decimal:
                p1 = min(one, max(zero, p2 + delta))
                b1 = _binom_vec_dec(n1, c1, p1)
                b2 = _binom_vec_dec(n2, c2, p2)
                acc = zero
                for x2, s_lo, s_hi in x1_slices:
                    acc += b2[x2] * sum(b1[s_lo:s_hi], zero)
                return acc

            if hi_p2 <= lo_p2 + Decimal("1e-20"):
                return prob_dec(half * (lo_p2 + hi_p2))

            grid_m = 160
            span = hi_p2 - lo_p2
            step = span / Decimal(grid_m)
            vals = [prob_dec(lo_p2 + step * Decimal(idx)) for idx in range(grid_m + 1)]
            best_v = max(vals)

            peak_indices: list[int] = []
            for idx in range(grid_m + 1):
                if vals[idx] <= zero:
                    continue
                l_ok = idx == 0 or vals[idx] >= vals[idx - 1]
                r_ok = idx == grid_m or vals[idx] > vals[idx + 1]
                if l_ok and r_ok:
                    peak_indices.append(idx)
            peak_indices.sort(key=lambda i: vals[i], reverse=True)
            peak_indices = peak_indices[:6]

            polish_steps = 40 if prec_digits >= 50 else 24
            for idx in peak_indices:
                center = lo_p2 + step * Decimal(idx)
                left = max(lo_p2, center - step)
                right = min(hi_p2, center + step)
                for _ in range(polish_steps):
                    m1 = left + (right - left) * third
                    m2 = right - (right - left) * third
                    v1 = prob_dec(m1)
                    v2 = prob_dec(m2)
                    if v1 > best_v:
                        best_v = v1
                    if v2 > best_v:
                        best_v = v2
                    if v1 < v2:
                        left = m1
                    else:
                        right = m2
            return min(one, max(zero, best_v))

    rd_obs = (Decimal(k1) / Decimal(n1)) - (Decimal(k2) / Decimal(n2))
    delta_grid = 250

    if k1 == 0 and k2 == n2:
        low_ci = Decimal(-1)
    else:
        d_start = Decimal("-1") + Decimal("1e-12")
        if tail_sup(d_start, True, prec_digits=28) > alpha_tail:
            low_ci = Decimal(-1)
        elif rd_obs <= d_start + Decimal("1e-12"):
            low_ci = Decimal(-1)
        else:
            step_d = (rd_obs - d_start) / Decimal(delta_grid)
            lo_b = d_start
            hi_b = rd_obs
            for m in range(1, delta_grid + 1):
                cand_d = rd_obs if m == delta_grid else (d_start + step_d * Decimal(m))
                if tail_sup(cand_d, True, prec_digits=28) > alpha_tail:
                    lo_b = d_start + step_d * Decimal(m - 1)
                    hi_b = cand_d
                    break
            for it in range(44):
                mid = half * (lo_b + hi_b)
                prec_use = 50 if it >= 28 else 28
                if tail_sup(mid, True, prec_digits=prec_use) > alpha_tail:
                    hi_b = mid
                else:
                    lo_b = mid
            low_ci = hi_b

    if k1 == n1 and k2 == 0:
        high_ci = one
    else:
        d_end = one - Decimal("1e-12")
        if tail_sup(d_end, False, prec_digits=28) > alpha_tail:
            high_ci = one
        elif d_end <= rd_obs + Decimal("1e-12"):
            high_ci = one
        else:
            step_d = (d_end - rd_obs) / Decimal(delta_grid)
            lo_b = rd_obs
            hi_b = d_end
            for m in range(1, delta_grid + 1):
                cand_d = rd_obs if m == delta_grid else (d_end - step_d * Decimal(m))
                if tail_sup(cand_d, False, prec_digits=28) > alpha_tail:
                    lo_b = cand_d
                    hi_b = d_end - step_d * Decimal(m - 1)
                    break
            for it in range(44):
                mid = half * (lo_b + hi_b)
                prec_use = 50 if it >= 28 else 28
                if tail_sup(mid, False, prec_digits=prec_use) > alpha_tail:
                    lo_b = mid
                else:
                    hi_b = mid
            high_ci = lo_b

    return low_ci, high_ci


def _check_ci_pair_dec(
    label: str,
    got_pair: Any,
    exp_pair: tuple[Decimal, Decimal],
    tol: Decimal,
    mismatches: list[str],
) -> None:
    if not isinstance(got_pair, list) or len(got_pair) != 2:
        mismatches.append(f"{label}: missing CI pair")
        return
    for idx, bound_name in enumerate(("low", "high")):
        if isinstance(got_pair[idx], bool) or not isinstance(got_pair[idx], (int, float)):
            mismatches.append(f"{label}.{bound_name}: non-numeric CI bound")
            continue
        diff = abs(Decimal(str(got_pair[idx])) - exp_pair[idx])
        if diff > tol:
            mismatches.append(f"{label}.{bound_name} drift {diff} > {tol}: got={got_pair[idx]} exp={exp_pair[idx]}")


def run_independent_audit(
    done_path: Path = DEFAULT_DONE,
    pred0015_results_path: Path = DEFAULT_PRED0015_RESULTS,
    fixture_path: Path = DEFAULT_FIXTURE,
    results_path: Path = DEFAULT_RESULTS,
    snapshot_path: Path = DEFAULT_SNAPSHOT,
    ledger_path: Path | None = DEFAULT_LEDGER,
) -> dict[str, Any]:
    verify_ast_independence()
    done_sha = _sha256_lf(done_path)
    p15_sha = _sha256_lf(pred0015_results_path)
    snap_sha = _sha256_lf(snapshot_path)
    fixture_sha = _sha256_lf(fixture_path)

    p15_doc = json.loads(pred0015_results_path.read_text(encoding="utf-8"))
    pe15 = p15_doc["primary_estimand_early_phase"]
    sys_set = set(pe15["systemic_ncts"])
    inh_set = set(pe15["inhaled_ncts"])
    dual_set = set(pe15["dual_route_inhaled_ncts"])
    inh_pure_set = inh_set - dual_set

    fix_doc = json.loads(fixture_path.read_text(encoding="utf-8"))
    if not isinstance(fix_doc, dict):
        return {"passed": False, "mismatches": ["Fixture root is not a dict"]}
    studies = fix_doc.get("studies")
    if not isinstance(studies, dict):
        return {"passed": False, "mismatches": ["Fixture studies is not a dict"]}

    mismatches: list[str] = []
    _scan_strings_for_pii_or_ctrl(fix_doc, "fixture_doc", mismatches)
    if fix_doc.get("nlm_attribution") != EXPECTED_NLM_ATTRIBUTION:
        mismatches.append("Fixture missing verbatim NLM attribution (AH-1)")
    if fix_doc.get("research_boundary") != EXPECTED_RESEARCH_BOUNDARY:
        mismatches.append("Fixture missing verbatim RESEARCH_BOUNDARY")
    if fix_doc.get("study_count") != 71 or isinstance(fix_doc.get("study_count"), bool):
        mismatches.append(f"Fixture study_count {fix_doc.get('study_count')!r} != 71")
    if p15_sha != EXPECTED_PRED0015_RESULTS_SHA256:
        mismatches.append(f"Canonical PRED-0015 SHA-256 mismatch: {p15_sha}")
    if snap_sha != EXPECTED_SNAPSHOT_SHA256:
        mismatches.append(f"Canonical snapshot.json SHA-256 mismatch: {snap_sha}")
    if set(studies.keys()) != (sys_set | inh_set):
        mismatches.append("Fixture NCT set != locked 71 NCTs")

    valid_sys_n, valid_inh_pure_n, valid_inh_all_n = 0, 0, 0
    k_sys, k_inh_pure, k_inh_all = 0, 0, 0
    cc_counts = {"sys": 0, "inh_pure": 0, "inh_all": 0}
    unfiltered_multi = {"sys": 0, "inh_pure": 0, "inh_all": 0}
    zero_total = {"sys": 0, "inh_pure": 0, "inh_all": 0}
    zero_omitted = {"sys": 0, "inh_pure": 0, "inh_all": 0}
    zero_withdrawn = {"sys": 0, "inh_pure": 0, "inh_all": 0}
    single_counts = {"sys": 0, "inh_pure": 0, "inh_all": 0}
    max_countries = {"sys": 0, "inh_pure": 0, "inh_all": 0}
    fac_lists: dict[str, list[int]] = {"sys": [], "inh_pure": [], "inh_all": []}
    us_counts = {"sys": 0, "inh_pure": 0, "inh_all": 0}
    ca_counts = {"sys": 0, "inh_pure": 0, "inh_all": 0}
    uk_counts = {"sys": 0, "inh_pure": 0, "inh_all": 0}
    au_counts = {"sys": 0, "inh_pure": 0, "inh_all": 0}
    eu27_counts = {"sys": 0, "inh_pure": 0, "inh_all": 0}
    sys_matched: list[str] = []
    inh_pure_matched: list[str] = []
    inh_all_matched: list[str] = []

    for nid in sorted(studies):
        if nid not in (sys_set | inh_set):
            continue
        st = studies[nid]
        if not isinstance(st, dict):
            mismatches.append(f"{nid}: study record is not a dict")
            continue
        if (set(st.keys()) & FORBIDDEN_PII) or (set(st.keys()) - ALLOWED_STUDY_KEYS):
            mismatches.append(f"{nid}: forbidden PII or unexpected study keys: {sorted(st.keys())}")
            continue
        exp_arm = "systemic_cftr" if nid in sys_set else "inhaled_cftr"
        exp_dual = nid in dual_set
        if st.get("nct_id") != nid:
            mismatches.append(f"{nid}: nct_id mismatch {st.get('nct_id')!r}")
            continue
        if st.get("arm") != exp_arm:
            mismatches.append(f"{nid}: arm mismatch {st.get('arm')!r} != {exp_arm!r}")
            continue
        if not isinstance(st.get("dual_route_modulator_co_admin"), bool) or st.get("dual_route_modulator_co_admin") is not exp_dual:
            mismatches.append(f"{nid}: dual_route_modulator_co_admin mismatch")
            continue

        mod_pres = st.get("contacts_locations_module_present")
        key_pres = st.get("locations_key_present")
        raw_locs = st.get("locations")
        if not isinstance(mod_pres, bool) or not isinstance(key_pres, bool):
            mismatches.append(f"{nid}: contacts_locations_module_present/locations_key_present must be bool")
            continue
        if (not mod_pres and (key_pres or raw_locs is not None)) or (not key_pres and raw_locs is not None) or (key_pres and raw_locs is None):
            mismatches.append(f"{nid}: inconsistent module/key presence flags vs locations")
            continue
        if raw_locs is not None and not isinstance(raw_locs, list):
            mismatches.append(f"{nid}: locations must be list or None")
            continue

        nw_set: set[str] = set()
        all_set: set[str] = set()
        nw_fac_n = 0
        loc_err = False
        for idx, loc in enumerate(raw_locs or []):
            if not isinstance(loc, dict):
                mismatches.append(f"{nid}: location[{idx}] not a dict")
                loc_err = True
                break
            if (set(loc.keys()) & FORBIDDEN_PII) or (set(loc.keys()) - ALLOWED_LOC_KEYS):
                mismatches.append(f"{nid}: forbidden PII or unexpected location keys at {idx}")
                loc_err = True
                break
            raw_c = loc.get("country")
            raw_s = loc.get("status")
            if not isinstance(raw_c, str) or "@" in raw_c or any(_is_ignorable_ctrl(ch) for ch in raw_c):
                mismatches.append(f"{nid}: invalid, PII, or raw-control country at {idx}")
                loc_err = True
                break
            if raw_s not in VALID_STATUSES:
                mismatches.append(f"{nid}: invalid facility status {raw_s!r}")
                loc_err = True
                break
            norm_k = " ".join(unicodedata.normalize("NFC", raw_c).split()).casefold()
            sov = AUDIT_COUNTRY_TABLE.get(norm_k)
            if sov is None:
                mismatches.append(f"{nid}: unmapped country {raw_c!r}")
                loc_err = True
                break
            all_set.add(sov)
            if raw_s != "WITHDRAWN":
                nw_set.add(sov)
                nw_fac_n += 1
        if loc_err:
            continue

        c_cnt = len(nw_set)
        omitted_empty = raw_locs is None or len(raw_locs) == 0
        hit = c_cnt >= 2
        hit_unf = len(all_set) >= 2
        has_cc = c_cnt >= 1

        def _bump(bucket: str) -> None:
            fac_lists[bucket].append(nw_fac_n)
            if c_cnt > max_countries[bucket]:
                max_countries[bucket] = c_cnt
            if c_cnt == 0:
                zero_total[bucket] += 1
                if omitted_empty:
                    zero_omitted[bucket] += 1
                else:
                    zero_withdrawn[bucket] += 1
            elif c_cnt == 1:
                single_counts[bucket] += 1
            if has_cc:
                cc_counts[bucket] += 1
            if hit_unf:
                unfiltered_multi[bucket] += 1
            if "United States" in nw_set:
                us_counts[bucket] += 1
            if "Canada" in nw_set:
                ca_counts[bucket] += 1
            if "United Kingdom" in nw_set:
                uk_counts[bucket] += 1
            if "Australia" in nw_set:
                au_counts[bucket] += 1
            if nw_set & EU27_SET:
                eu27_counts[bucket] += 1

        if nid in sys_set:
            valid_sys_n += 1
            _bump("sys")
            if hit:
                k_sys += 1
                sys_matched.append(nid)
        else:
            valid_inh_all_n += 1
            _bump("inh_all")
            if hit:
                k_inh_all += 1
                inh_all_matched.append(nid)
            if nid in inh_pure_set:
                valid_inh_pure_n += 1
                _bump("inh_pure")
                if hit:
                    k_inh_pure += 1
                    inh_pure_matched.append(nid)

    if valid_sys_n == 0 or valid_inh_pure_n == 0 or valid_inh_all_n == 0:
        mismatches.append(f"Zero valid resolved studies in one or more arms: {(valid_sys_n, valid_inh_pure_n, valid_inh_all_n)}")
        return {"passed": False, "mismatches": mismatches}

    n_sys, n_inh_pure, n_inh_all = valid_sys_n, valid_inh_pure_n, valid_inh_all_n

    p_pure_1s, p_pure_2s = _fisher_tails_multiplicative(n_sys, k_sys, n_inh_pure, k_inh_pure)
    p_all_1s, _ = _fisher_tails_multiplicative(n_sys, k_sys, n_inh_all, k_inh_all)
    nw_pure_lo, nw_pure_hi = _newcombe_wilson_decimal(n_sys, k_sys, n_inh_pure, k_inh_pure)
    cz_pure_lo, cz_pure_hi = _chan_zhang_ci95_decimal(n_sys, k_sys, n_inh_pure, k_inh_pure)
    nw_all_lo, nw_all_hi = _newcombe_wilson_decimal(n_sys, k_sys, n_inh_all, k_inh_all)
    cz_all_lo, cz_all_hi = _chan_zhang_ci95_decimal(n_sys, k_sys, n_inh_all, k_inh_all)

    crit_1 = (valid_sys_n, valid_inh_pure_n, valid_inh_all_n) == (54, 15, 17)
    crit_2 = Fraction(k_sys, n_sys) > Fraction(k_inh_pure, n_inh_pure)
    crit_3 = p_pure_1s < Fraction(1, 100)
    crit_4 = (cz_pure_lo > Decimal(0)) and (nw_pure_lo > Decimal(0))
    crit_5 = (Fraction(k_sys, n_sys) > Fraction(k_inh_all, n_inh_all)) and (p_all_1s < Fraction(1, 100))
    expected_outcome = "supported" if (crit_1 and crit_2 and crit_3 and crit_4 and crit_5) else "refuted"

    res_doc = json.loads(results_path.read_text(encoding="utf-8"))
    prim_res = res_doc.get("primary_estimand_pure_route_multi_country") or {}
    sec_res = res_doc.get("secondary_gate_all_17_inhaled_including_dual_route") or {}
    sc_res = res_doc.get("support_criteria_evaluation") or {}
    desc_res = res_doc.get("descriptive_cross_scale_non_gated") or {}

    if res_doc.get("nlm_attribution") != EXPECTED_NLM_ATTRIBUTION:
        mismatches.append("results_doc missing verbatim NLM attribution")
    if res_doc.get("research_boundary") != EXPECTED_RESEARCH_BOUNDARY:
        mismatches.append("results_doc missing verbatim RESEARCH_BOUNDARY")
    if res_doc.get("done_sha256") != done_sha:
        mismatches.append(f"done_sha256 mismatch: {res_doc.get('done_sha256')} != {done_sha}")
    if res_doc.get("pred0015_results_sha256") != p15_sha:
        mismatches.append("pred0015_results_sha256 mismatch")
    if res_doc.get("snapshot_sha256") != snap_sha:
        mismatches.append("snapshot_sha256 mismatch")
    if res_doc.get("fixture_sha256") != fixture_sha:
        mismatches.append("fixture_sha256 mismatch")
    if prim_res.get("n_sys") != n_sys or prim_res.get("n_inh_pure") != n_inh_pure:
        mismatches.append(f"Primary denominators mismatch: {prim_res.get('n_sys')}/{prim_res.get('n_inh_pure')} != {n_sys}/{n_inh_pure}")
    if prim_res.get("k_sys") != k_sys or prim_res.get("k_inh_pure") != k_inh_pure:
        mismatches.append(f"Primary counts mismatch: {prim_res.get('k_sys')}/{prim_res.get('k_inh_pure')} != {k_sys}/{k_inh_pure}")
    if sec_res.get("n_sys") != n_sys or sec_res.get("k_sys") != k_sys or sec_res.get("n_inh_all") != n_inh_all:
        mismatches.append("Secondary denominators/k_sys mismatch")
    if sec_res.get("k_inh_all") != k_inh_all:
        mismatches.append(f"All-17 inhaled count mismatch: {sec_res.get('k_inh_all')} != {k_inh_all}")
    if prim_res.get("systemic_matched_ncts") != sys_matched:
        mismatches.append("systemic_matched_ncts mismatch")
    if prim_res.get("inhaled_pure_matched_ncts") != inh_pure_matched:
        mismatches.append("inhaled_pure_matched_ncts mismatch")
    if sec_res.get("inhaled_all_matched_ncts") != inh_all_matched:
        mismatches.append("inhaled_all_matched_ncts mismatch")
    if prim_res.get("p_fisher_one_sided_exact") != f"{p_pure_1s.numerator}/{p_pure_1s.denominator}":
        mismatches.append("p_fisher_one_sided_exact pure mismatch")
    if prim_res.get("p_fisher_two_sided_exact") != f"{p_pure_2s.numerator}/{p_pure_2s.denominator}":
        mismatches.append("p_fisher_two_sided_exact pure mismatch")
    if sec_res.get("p_fisher_one_sided_exact") != f"{p_all_1s.numerator}/{p_all_1s.denominator}":
        mismatches.append("p_fisher_one_sided_exact all-17 mismatch")

    tol = Decimal("0.000001")
    tol_p = Decimal("0.00000001")
    p_sys_dec = Decimal(k_sys) / Decimal(n_sys)
    p_inh_pure_dec = Decimal(k_inh_pure) / Decimal(n_inh_pure)
    p_inh_all_dec = Decimal(k_inh_all) / Decimal(n_inh_all)
    rd_pure_dec = p_sys_dec - p_inh_pure_dec
    rd_all_dec = p_sys_dec - p_inh_all_dec
    p_pure_1s_dec = Decimal(p_pure_1s.numerator) / Decimal(p_pure_1s.denominator)
    p_pure_2s_dec = Decimal(p_pure_2s.numerator) / Decimal(p_pure_2s.denominator)
    p_all_1s_dec = Decimal(p_all_1s.numerator) / Decimal(p_all_1s.denominator)

    for field_label, got_val, exp_dec, bound_tol in (
        ("prim.p_sys", prim_res.get("p_sys"), p_sys_dec, tol),
        ("prim.p_inh_pure", prim_res.get("p_inh_pure"), p_inh_pure_dec, tol),
        ("prim.risk_difference", prim_res.get("risk_difference"), rd_pure_dec, tol),
        ("prim.p_fisher_one_sided", prim_res.get("p_fisher_one_sided"), p_pure_1s_dec, tol_p),
        ("prim.p_fisher_two_sided", prim_res.get("p_fisher_two_sided"), p_pure_2s_dec, tol_p),
        ("sec.p_sys", sec_res.get("p_sys"), p_sys_dec, tol),
        ("sec.p_inh_all", sec_res.get("p_inh_all"), p_inh_all_dec, tol),
        ("sec.risk_difference", sec_res.get("risk_difference"), rd_all_dec, tol),
        ("sec.p_fisher_one_sided", sec_res.get("p_fisher_one_sided"), p_all_1s_dec, tol_p),
    ):
        if isinstance(got_val, bool) or not isinstance(got_val, (int, float)):
            mismatches.append(f"{field_label}: missing or non-numeric ({got_val!r})")
        elif abs(Decimal(str(got_val)) - exp_dec) > bound_tol:
            mismatches.append(f"{field_label}: {got_val} != {exp_dec}")

    _check_ci_pair_dec("chan_zhang_pure", prim_res.get("chan_zhang_unconditional_ci95"), (cz_pure_lo, cz_pure_hi), tol, mismatches)
    _check_ci_pair_dec("newcombe_pure", prim_res.get("newcombe_wilson_ci95"), (nw_pure_lo, nw_pure_hi), tol, mismatches)
    _check_ci_pair_dec("chan_zhang_all17", sec_res.get("chan_zhang_unconditional_ci95"), (cz_all_lo, cz_all_hi), tol, mismatches)
    _check_ci_pair_dec("newcombe_all17", sec_res.get("newcombe_wilson_ci95"), (nw_all_lo, nw_all_hi), tol, mismatches)

    for crit_key, exp_bool in (
        ("crit_1_exact_54_15_17_margins", crit_1),
        ("crit_2_directional_pure", crit_2),
        ("crit_3_fisher_pure_lt_0_01", crit_3),
        ("crit_4_chan_zhang_and_newcombe_ci95_low_gt_0", crit_4),
        ("crit_5_all_17_inhaled_preserved_lt_0_01", crit_5),
    ):
        if sc_res.get(crit_key) is not exp_bool:
            mismatches.append(f"{crit_key} mismatch: {sc_res.get(crit_key)!r} != {exp_bool!r}")

    if sc_res.get("outcome") != expected_outcome:
        mismatches.append(f"outcome mismatch: {sc_res.get('outcome')} != {expected_outcome}")

    for arm_key, short_key, exp_n, exp_k in (
        ("systemic_cftr", "sys", n_sys, k_sys),
        ("inhaled_cftr_pure", "inh_pure", n_inh_pure, k_inh_pure),
        ("inhaled_cftr_all", "inh_all", n_inh_all, k_inh_all),
    ):
        arm_blk = desc_res.get(arm_key) or {}
        if arm_blk.get("n") != exp_n or arm_blk.get("multi_country_count") != exp_k:
            mismatches.append(f"descriptive.{arm_key} n/multi_country_count mismatch")
        if arm_blk.get("unfiltered_including_withdrawn_multi_country_count") != unfiltered_multi[short_key]:
            mismatches.append(f"descriptive.{arm_key}.unfiltered_multi mismatch")
        ccb = arm_blk.get("country_count_breakdown") or {}
        if (
            ccb.get("zero_countries_total") != zero_total[short_key]
            or ccb.get("zero_countries_omitted_or_empty_locations") != zero_omitted[short_key]
            or ccb.get("zero_countries_all_withdrawn") != zero_withdrawn[short_key]
            or ccb.get("single_country_count") != single_counts[short_key]
            or ccb.get("multi_country_ge_2_count") != exp_k
            or ccb.get("max_sovereign_countries_in_single_study") != max_countries[short_key]
        ):
            mismatches.append(f"descriptive.{arm_key}.country_count_breakdown mismatch")

        cc_blk = arm_blk.get("complete_case_ge_1_country") or {}
        if cc_blk.get("n_complete") != cc_counts[short_key] or cc_blk.get("multi_country_count") != exp_k:
            mismatches.append(f"descriptive.{arm_key}.n_complete mismatch")

        fac_blk = arm_blk.get("non_withdrawn_facility_summary") or {}
        f_all = fac_lists[short_key]
        f_pos = [v for v in f_all if v >= 1]
        med_all = _median_decimal(f_all)
        med_pos = _median_decimal(f_pos)
        got_med_all = fac_blk.get("median_facilities_all_studies")
        got_med_pos = fac_blk.get("median_facilities_complete_case")
        if med_all is not None and (isinstance(got_med_all, bool) or not isinstance(got_med_all, (int, float)) or Decimal(str(got_med_all)) != med_all):
            mismatches.append(f"descriptive.{arm_key}.median_facilities_all_studies mismatch")
        if med_pos is None:
            if got_med_pos is not None:
                mismatches.append(f"descriptive.{arm_key}.median_facilities_complete_case expected None")
        elif isinstance(got_med_pos, bool) or not isinstance(got_med_pos, (int, float)) or Decimal(str(got_med_pos)) != med_pos:
            mismatches.append(f"descriptive.{arm_key}.median_facilities_complete_case mismatch")
        if (
            fac_blk.get("min_facilities") != min(f_all)
            or fac_blk.get("max_facilities") != max(f_all)
            or fac_blk.get("total_non_withdrawn_facilities_sum") != sum(f_all)
        ):
            mismatches.append(f"descriptive.{arm_key}.facility min/max/sum mismatch")

        jur_blk = arm_blk.get("jurisdiction_site_presence_non_withdrawn") or {}
        for jur_name, jur_exp in (
            ("United States", us_counts[short_key]),
            ("Canada", ca_counts[short_key]),
            ("United Kingdom", uk_counts[short_key]),
            ("Australia", au_counts[short_key]),
            ("EU-27", eu27_counts[short_key]),
        ):
            if (jur_blk.get(jur_name) or {}).get("study_count") != jur_exp:
                mismatches.append(f"descriptive.{arm_key}.jurisdiction.{jur_name} mismatch")

    for cc_key, n_inh_cc_exp, k_inh_cc_exp in (
        ("complete_case_sensitivity_pure", cc_counts["inh_pure"], k_inh_pure),
        ("complete_case_sensitivity_all_17", cc_counts["inh_all"], k_inh_all),
    ):
        cc_res = desc_res.get(cc_key) or {}
        n_sys_cc_exp = cc_counts["sys"]
        if n_sys_cc_exp > 0 and n_inh_cc_exp > 0:
            p_cc_1s, _ = _fisher_tails_multiplicative(n_sys_cc_exp, k_sys, n_inh_cc_exp, k_inh_cc_exp)
            nw_cc = _newcombe_wilson_decimal(n_sys_cc_exp, k_sys, n_inh_cc_exp, k_inh_cc_exp)
            cz_cc = _chan_zhang_ci95_decimal(n_sys_cc_exp, k_sys, n_inh_cc_exp, k_inh_cc_exp)
            p_cc_dec = Decimal(p_cc_1s.numerator) / Decimal(p_cc_1s.denominator)
            if (
                cc_res.get("n_sys_cc") != n_sys_cc_exp
                or cc_res.get("k_sys") != k_sys
                or cc_res.get("n_inh_cc") != n_inh_cc_exp
                or cc_res.get("k_inh") != k_inh_cc_exp
                or cc_res.get("p_fisher_one_sided_exact") != f"{p_cc_1s.numerator}/{p_cc_1s.denominator}"
            ):
                mismatches.append(f"descriptive.{cc_key} counts/exact-p mismatch")
            got_cc_p = cc_res.get("p_fisher_one_sided")
            if isinstance(got_cc_p, bool) or not isinstance(got_cc_p, (int, float)) or abs(Decimal(str(got_cc_p)) - p_cc_dec) > tol_p:
                mismatches.append(f"descriptive.{cc_key}.p_fisher_one_sided mismatch")
            _check_ci_pair_dec(f"{cc_key}.chan_zhang", cc_res.get("chan_zhang_unconditional_ci95"), cz_cc, tol, mismatches)
            _check_ci_pair_dec(f"{cc_key}.newcombe", cc_res.get("newcombe_wilson_ci95"), nw_cc, tol, mismatches)

    if ledger_path is not None and ledger_path.is_file():
        led_doc = json.loads(ledger_path.read_text(encoding="utf-8"))
        if led_doc.get("done_sha256") != done_sha:
            mismatches.append("ledger.json done_sha256 mismatch")
        today_iso = dt.date.today().isoformat()
        for cl in led_doc.get("claims") or []:
            c_on = str(cl.get("checked_on") or "")
            if c_on != str(res_doc.get("evaluated_on_local")) or c_on > today_iso:
                mismatches.append(f"ledger.json {cl.get('id')} checked_on {c_on!r} mismatch vs {res_doc.get('evaluated_on_local')!r}")

    return {
        "passed": len(mismatches) == 0,
        "mismatches": mismatches,
        "n_sys": n_sys,
        "k_sys": k_sys,
        "n_inh_pure": n_inh_pure,
        "k_inh_pure": k_inh_pure,
        "n_inh_all": n_inh_all,
        "k_inh_all": k_inh_all,
        "outcome": expected_outcome,
        "p_fisher_pure_one_sided_exact": f"{p_pure_1s.numerator}/{p_pure_1s.denominator}",
        "p_fisher_all_one_sided_exact": f"{p_all_1s.numerator}/{p_all_1s.denominator}",
        "chan_zhang_pure_ci95": [str(cz_pure_lo.quantize(Decimal("0.000001"))), str(cz_pure_hi.quantize(Decimal("0.000001")))],
        "newcombe_wilson_pure_ci95": [str(nw_pure_lo.quantize(Decimal("0.000001"))), str(nw_pure_hi.quantize(Decimal("0.000001")))],
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--done", type=Path, default=DEFAULT_DONE)
    ap.add_argument("--pred0015-results", type=Path, default=DEFAULT_PRED0015_RESULTS)
    ap.add_argument("--snapshot", type=Path, default=DEFAULT_SNAPSHOT)
    ap.add_argument("--fixture", type=Path, default=DEFAULT_FIXTURE)
    ap.add_argument("--results", type=Path, default=DEFAULT_RESULTS)
    ap.add_argument("--ledger", type=Path, default=DEFAULT_LEDGER)
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--json", action="store_true")
    args = ap.parse_args(argv)

    rep = run_independent_audit(args.done, args.pred0015_results, args.fixture, args.results, args.snapshot, args.ledger)
    if args.json:
        print(json.dumps(rep, indent=2))
    if not rep["passed"]:
        for m in rep["mismatches"]:
            print(f"FAIL independent_audit: {m}", file=sys.stderr)
        return 1
    if not args.json:
        print(
            f"PASS independent_audit: outcome={rep['outcome']} sys={rep['k_sys']}/{rep['n_sys']} "
            f"inh_pure={rep['k_inh_pure']}/{rep['n_inh_pure']} inh_all={rep['k_inh_all']}/{rep['n_inh_all']} "
            f"p_pure_1s={rep['p_fisher_pure_one_sided_exact']} "
            f"cz_ci95={rep['chan_zhang_pure_ci95']} nw_ci95={rep['newcombe_wilson_pure_ci95']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
