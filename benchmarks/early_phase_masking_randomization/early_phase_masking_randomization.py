#!/usr/bin/env python3
"""Early-Phase Randomization & Masking Asymmetry Benchmark (`PRED-0017`).

Tests the pre-registered design-asymmetry hypothesis (`benchmarks/early_phase_masking_randomization/done.md`)
across the locked 71 early-phase non-healthy-volunteer-only PRISMA-940 Cystic Fibrosis interventional trials
from `PRED-0015` and `PRED-0016` (`benchmarks/delivery_route_compartmentalization/results_2026-11-24.json`):

- Primary comparison: `systemic_cftr` (`N_sys == 54`) vs pure airway-restricted `inhaled_cftr_pure`
  (`N_inh_pure == 15`, excluding dual-route background systemic modulator trials `NCT05668741` and `NCT06237335`),
  with all-17 `inhaled_cftr` (`N_inh_all == 17`) evaluated alongside as a mandatory secondary gate.
- Primary indicator (`is_randomized_and_masked`): `True` iff `designInfo.allocation == "RANDOMIZED"`
  AND `designInfo.maskingInfo.masking in {"SINGLE", "DOUBLE", "TRIPLE", "QUADRUPLE"}` in frozen
  ClinicalTrials.gov API v2 `early_phase_design_eligibility_2026-11-24.json`.
- Non-gated descriptive cross-scale context: component-wise `is_randomized` and `is_masked`, complete-case
  sensitivity, phase-stratified (`EARLY_PHASE1`/`PHASE1` vs `PHASE1/PHASE2` + `PHASE2`) proportions,
  median `enrollmentInfo.count`, and `eligibilityModule.stdAges` pediatric `"CHILD"` inclusion juxtaposed
  with `ACFDR 2025` (`11.5%` pediatric [`170/1,473`] vs `5.8%` adult [`2,326` adult cohort] non-transplant
  PBS modulator ineligibility) under the Non-Pooling Invariant.

Invariants (`cf-research-context`, `discovery-practice`, `security-browsing`, `AH-1`..`AH-4`):
- 100% pure Python 3.10+ standard library + local `tools/trial_atlas/` modules; zero network/LLM calls.
- Strict `AH-1` verbatim NLM attribution and `AH-2` PRISMA-940 cohort boundary.
- Omitted/`NA` `allocation` and omitted/`NONE` `masking` retained in denominators with separate missingness counts.
- Sawtooth-safe upward/downward `delta` bracket scan + multi-peak `p2` golden-section polish for Chan-Zhang 95% CIs.
- CodeQL-safe identifiers (`RESEARCH_BOUNDARY`, `NLM_ATTRIBUTION`, `quarantine_notice`).
"""
from __future__ import annotations

import argparse
import copy
from fractions import Fraction
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

BENCHMARK_VERSION = "2026-11-24.v1"
EXPECTED_PRED0015_RESULTS_SHA256 = "c4288a7c40783c746e1e4c103e9b1b8867818cb97cc371d19348715d9251dcb1"
EXPECTED_SNAPSHOT_SHA256 = "de0f330c4441d77b84dc5f547e4dfcf99195c1b9f223e4dbeebd628054925b2d"

DEFAULT_DONE_MD = HERE / "done.md"
DEFAULT_PRED0015_RESULTS_JSON = PRED0015_DIR / "results_2026-11-24.json"
DEFAULT_SNAPSHOT_JSON = TRIAL_ATLAS_DIR / "snapshot.json"
DEFAULT_FIXTURE_JSON = HERE / "early_phase_design_eligibility_2026-11-24.json"
DEFAULT_RESULTS_JSON = HERE / "results_2026-11-24.json"

NLM_ATTRIBUTION = (
    "ClinicalTrials.gov does not evaluate safety or efficacy; "
    "listing on the registry does not denote government endorsement."
)
RESEARCH_BOUNDARY = (
    "Independent open-source research benchmark only (Tier T1/T2 registry metadata); "
    "not clinical guidance and not an assessment of personal trial eligibility or treatment efficacy."
)

ALLOWED_ALLOCATIONS: frozenset[str | None] = frozenset({"RANDOMIZED", "NON_RANDOMIZED", "NA", None})
ALLOWED_MASKINGS: frozenset[str | None] = frozenset({"NONE", "SINGLE", "DOUBLE", "TRIPLE", "QUADRUPLE", None})
MASKED_VALUES: frozenset[str] = frozenset({"SINGLE", "DOUBLE", "TRIPLE", "QUADRUPLE"})
ALLOWED_STD_AGES: frozenset[str] = frozenset({"CHILD", "ADULT", "OLDER_ADULT"})
ALLOWED_ENROLLMENT_TYPES: frozenset[str | None] = frozenset({"ACTUAL", "ESTIMATED", None})
ALLOWED_PHASES: frozenset[str] = frozenset({"EARLY_PHASE1", "PHASE1", "PHASE2", "PHASE3", "PHASE4", "NA"})

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
        "\uff1c": "<", "\u27e8": "<", "\u3008": "<", "\u2039": "<", "\ufe64": "<", "\u226e": "<",
        "\uff1e": ">", "\u27e9": ">", "\u3009": ">", "\u203a": ">", "\ufe65": ">", "\u226f": ">",
    }
)


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


def _validate_2x2_counts(n1: int, k1: int, n2: int, k2: int) -> None:
    for name, val in (("n1", n1), ("k1", k1), ("n2", n2), ("k2", k2)):
        if isinstance(val, bool) or not isinstance(val, int):
            raise ValueError(f"Invalid integer type for {name}: {val!r}")
    if n1 <= 0 or n2 <= 0 or not (0 <= k1 <= n1) or not (0 <= k2 <= n2):
        raise ValueError(f"Invalid 2x2 counts: n1={n1}, k1={k1}, n2={n2}, k2={k2}")


class EarlyPhaseMaskingBenchmarkError(RuntimeError):
    """Raised when a pre-registered schema, provenance, positive-control, or falsification check fails."""


def sha256_lf_file(path: Path) -> str:
    txt = path.read_text(encoding="utf-8").replace("\r\n", "\n")
    return hashlib.sha256(txt.encode("utf-8")).hexdigest()


def _is_closing_data_tag_span(span: str) -> bool:
    skeleton_chars: list[str] = []
    for ch in span:
        if _is_default_ignorable_or_control(ch):
            continue
        norm = unicodedata.normalize("NFKD", ch.translate(_CONFUSABLE_SKELETON_MAP))
        for sub in norm:
            if unicodedata.category(sub) == "Mn" or _is_default_ignorable_or_control(sub) or sub.isspace():
                continue
            skeleton_chars.append(sub.lower())
    folded = "".join(skeleton_chars)
    return folded.startswith("</") and folded.endswith(">") and ("data" in folded or "untrusted" in folded)


def sanitize_registry_text(value: Any) -> str:
    """Linear O(N) escape of format/bidi/default-ignorable characters and closing-tag lookalikes."""
    if value is None:
        return ""
    raw = str(value)
    n = len(raw)
    tag_safe_parts: list[str] = []
    idx = 0
    while idx < n:
        ch = raw[idx]
        mapped = ch.translate(_CONFUSABLE_SKELETON_MAP)
        if mapped == "<":
            end_idx = -1
            scan = idx + 1
            while scan < n:
                s_mapped = raw[scan].translate(_CONFUSABLE_SKELETON_MAP)
                if s_mapped == "<":
                    break
                if s_mapped == ">":
                    end_idx = scan
                    break
                scan += 1
            if end_idx != -1 and _is_closing_data_tag_span(raw[idx : end_idx + 1]):
                tag_safe_parts.append("[ESCAPED_CLOSING_TAG]")
                idx = end_idx + 1
                continue
        tag_safe_parts.append(ch)
        idx += 1

    inter = "".join(tag_safe_parts)
    out_chars: list[str] = []
    for ch in inter:
        if _is_default_ignorable_or_control(ch):
            out_chars.append(f"<U+{ord(ch):04X}>")
        else:
            out_chars.append(ch)
    return "".join(out_chars)


def classify_design_enums(allocation: Any, masking: Any, nct_id: str = "synthetic") -> dict[str, bool]:
    """Validate NLM `allocation` and `masking` enums fail-closed and compute design indicators."""
    if allocation is not None and not isinstance(allocation, str):
        raise EarlyPhaseMaskingBenchmarkError(f"{nct_id}: non-string allocation {allocation!r}")
    if masking is not None and not isinstance(masking, str):
        raise EarlyPhaseMaskingBenchmarkError(f"{nct_id}: non-string masking {masking!r}")
    if allocation not in ALLOWED_ALLOCATIONS:
        raise EarlyPhaseMaskingBenchmarkError(
            f"{nct_id}: unexpected NLM allocation enum {allocation!r} outside {sorted(x for x in ALLOWED_ALLOCATIONS if x)}"
        )
    if masking not in ALLOWED_MASKINGS:
        raise EarlyPhaseMaskingBenchmarkError(
            f"{nct_id}: unexpected NLM masking enum {masking!r} outside {sorted(x for x in ALLOWED_MASKINGS if x)}"
        )
    is_rand = allocation == "RANDOMIZED"
    is_mask = masking in MASKED_VALUES
    return {
        "is_randomized": is_rand,
        "is_masked": is_mask,
        "is_randomized_and_masked": is_rand and is_mask,
        "is_complete_case_alloc_and_mask": (allocation in {"RANDOMIZED", "NON_RANDOMIZED"}) and (masking is not None),
    }


def parse_study_design_eligibility(record: dict[str, Any], nct_id: str = "synthetic") -> dict[str, Any]:
    """Validate and normalize a frozen ClinicalTrials.gov study design/eligibility record."""
    if not isinstance(record, dict):
        raise EarlyPhaseMaskingBenchmarkError(f"{nct_id}: study record must be a dict")

    design_mod = record.get("design_module")
    if design_mod is None:
        design_mod = {}
    if not isinstance(design_mod, dict):
        raise EarlyPhaseMaskingBenchmarkError(f"{nct_id}: design_module must be dict or None")

    raw_design_info = design_mod.get("designInfo")
    if raw_design_info is not None and not isinstance(raw_design_info, dict):
        raise EarlyPhaseMaskingBenchmarkError(f"{nct_id}: designInfo must be dict or None")
    design_info_present = raw_design_info is not None
    design_info = {} if raw_design_info is None else raw_design_info

    allocation_key_present = "allocation" in design_info and design_info.get("allocation") is not None
    allocation = design_info.get("allocation")

    raw_masking_info = design_info.get("maskingInfo")
    if raw_masking_info is not None and not isinstance(raw_masking_info, dict):
        raise EarlyPhaseMaskingBenchmarkError(f"{nct_id}: maskingInfo must be dict or None")
    masking_info_present = raw_masking_info is not None
    masking_info = {} if raw_masking_info is None else raw_masking_info
    masking_key_present = "masking" in masking_info and masking_info.get("masking") is not None
    masking = masking_info.get("masking")

    cls = classify_design_enums(allocation, masking, nct_id=nct_id)

    intervention_model = design_info.get("interventionModel")
    if intervention_model is not None:
        intervention_model = sanitize_registry_text(intervention_model)

    raw_enroll_info = design_mod.get("enrollmentInfo")
    if raw_enroll_info is not None and not isinstance(raw_enroll_info, dict):
        raise EarlyPhaseMaskingBenchmarkError(f"{nct_id}: enrollmentInfo must be dict or None")
    enroll_info = {} if raw_enroll_info is None else raw_enroll_info
    enroll_count = enroll_info.get("count")
    if enroll_count is not None and (isinstance(enroll_count, bool) or not isinstance(enroll_count, int) or enroll_count < 0):
        raise EarlyPhaseMaskingBenchmarkError(f"{nct_id}: invalid enrollmentInfo.count {enroll_count!r}")
    enroll_type = enroll_info.get("type")
    if enroll_type not in ALLOWED_ENROLLMENT_TYPES:
        raise EarlyPhaseMaskingBenchmarkError(f"{nct_id}: invalid enrollmentInfo.type {enroll_type!r}")

    elig_mod = record.get("eligibility_module")
    if elig_mod is None:
        elig_mod = {}
    if not isinstance(elig_mod, dict):
        raise EarlyPhaseMaskingBenchmarkError(f"{nct_id}: eligibility_module must be dict or None")

    std_ages_raw = elig_mod.get("stdAges")
    if std_ages_raw is None:
        std_ages: list[str] = []
    elif isinstance(std_ages_raw, list) and all(isinstance(x, str) and x in ALLOWED_STD_AGES for x in std_ages_raw):
        std_ages = list(std_ages_raw)
    else:
        raise EarlyPhaseMaskingBenchmarkError(f"{nct_id}: invalid eligibilityModule.stdAges {std_ages_raw!r}")

    hv_raw = elig_mod.get("healthyVolunteers")
    if hv_raw is not None and not isinstance(hv_raw, bool):
        raise EarlyPhaseMaskingBenchmarkError(f"{nct_id}: invalid eligibilityModule.healthyVolunteers {hv_raw!r}")

    min_age = sanitize_registry_text(elig_mod.get("minimumAge")) if elig_mod.get("minimumAge") is not None else None
    max_age = sanitize_registry_text(elig_mod.get("maximumAge")) if elig_mod.get("maximumAge") is not None else None
    sex = sanitize_registry_text(elig_mod.get("sex")) if elig_mod.get("sex") is not None else None

    phases_raw = design_mod.get("phases")
    if phases_raw is None:
        phases: list[str] = []
    elif isinstance(phases_raw, list) and all(isinstance(x, str) and x in ALLOWED_PHASES for x in phases_raw):
        phases = list(phases_raw)
    else:
        raise EarlyPhaseMaskingBenchmarkError(f"{nct_id}: invalid designModule.phases {phases_raw!r}")

    return {
        "nct_id": nct_id,
        "design_info_present": design_info_present,
        "allocation_key_present": allocation_key_present,
        "allocation": allocation,
        "masking_info_present": masking_info_present,
        "masking_key_present": masking_key_present,
        "masking": masking,
        "intervention_model": intervention_model,
        "phases": phases,
        "enrollment_count": enroll_count,
        "enrollment_type": enroll_type,
        "std_ages": std_ages,
        "includes_child_std_age": "CHILD" in std_ages,
        "healthy_volunteers": hv_raw,
        "minimum_age": min_age,
        "maximum_age": max_age,
        "sex": sex,
        **cls,
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
        raise EarlyPhaseMaskingBenchmarkError(f"Hypergeometric PMF sum != 1: {support_sum}")
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
    """Conservative unconditional Chan-Zhang (1999) 95% CI with multi-peak `p2` polish
    and upward/downward `delta` bracket scanning (`inf{delta: P_U* > alpha/2}`,
    `sup{delta: P_L* > alpha/2}`) safe under non-zero `k2 > 0` sawtooth discontinuities.
    """
    _validate_2x2_counts(n1, k1, n2, k2)
    if (
        isinstance(alpha_two_sided, bool)
        or not isinstance(alpha_two_sided, (int, float))
        or not math.isfinite(alpha_two_sided)
        or not (0.0 < alpha_two_sided < 1.0)
    ):
        raise ValueError(f"Invalid alpha_two_sided parameter: {alpha_two_sided!r}")
    alpha_tail = 0.5 * alpha_two_sided
    comb1 = [math.comb(n1, i) for i in range(n1 + 1)]
    comb2 = [math.comb(n2, j) for j in range(n2 + 1)]
    inv_phi = (math.sqrt(5.0) - 1.0) / 2.0

    def tail_supremum(delta: float, upper: bool) -> float:
        low_p2 = max(0.0, -delta)
        high_p2 = min(1.0, 1.0 - delta)
        z_obs = _mn_score_z(k1, n1, k2, n2, delta)
        tie_tol = 1e-11
        active_pairs: list[tuple[int, int]] = []
        for x1 in range(n1 + 1):
            p1_hat = x1 / n1
            for x2 in range(n2 + 1):
                num_xy = p1_hat - (x2 / n2) - delta
                if upper and z_obs > tie_tol and num_xy < 0.0:
                    continue
                if (not upper) and z_obs < -tie_tol and num_xy > 0.0:
                    continue
                z_xy = _mn_score_z(x1, n1, x2, n2, delta)
                if (upper and z_xy >= z_obs - tie_tol) or ((not upper) and z_xy <= z_obs + tie_tol):
                    active_pairs.append((x1, x2))
        if not active_pairs:
            return 0.0
        if len(active_pairs) == (n1 + 1) * (n2 + 1):
            return 1.0

        def prob_at(p2: float) -> float:
            p1 = min(1.0, max(0.0, p2 + delta))
            b1 = _binom_pmf_vec(n1, comb1, p1)
            b2 = _binom_pmf_vec(n2, comb2, p2)
            return sum(b1[x1] * b2[x2] for x1, x2 in active_pairs)

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
    """Run all four pre-registered positive controls before reading empirical designModule/eligibilityModule."""
    # Control 1a: k_inh == 0 calibration (12, 9, 10, 0)
    cal0 = fisher_hypergeometric_exact(12, 9, 10, 0)
    if cal0["p_one_sided_fraction"] != Fraction(11, 24871):
        raise EarlyPhaseMaskingBenchmarkError(f"Calibration (12,9,10,0) Fisher mismatch: {cal0['p_one_sided_exact']}")
    cz0 = chan_zhang_unconditional_rd_ci95(12, 9, 10, 0)
    nw0 = newcombe_wilson_rd_ci95(12, 9, 10, 0)
    if cz0[0] <= 0.35 or nw0[0] <= 0.35:
        raise EarlyPhaseMaskingBenchmarkError(f"Calibration (12,9,10,0) CI95 <= 0.35: cz={cz0}, nw={nw0}")

    # Control 1b: Non-zero k_inh > 0 calibration (54, 35, 15, 3)
    cal_pos = fisher_hypergeometric_exact(54, 35, 15, 3)
    if cal_pos["p_one_sided_fraction"] >= ALPHA_ONE_SIDED_THRESHOLD:
        raise EarlyPhaseMaskingBenchmarkError(
            f"Calibration (54,35,15,3) Fisher p >= 0.01: {cal_pos['p_one_sided_exact']}"
        )
    cz_pos = chan_zhang_unconditional_rd_ci95(54, 35, 15, 3)
    nw_pos = newcombe_wilson_rd_ci95(54, 35, 15, 3)
    if cz_pos[0] <= 0.07 or nw_pos[0] <= 0.15:
        raise EarlyPhaseMaskingBenchmarkError(
            f"Calibration (54,35,15,3) CI95 bounds failed: cz={cz_pos}, nw={nw_pos}"
        )

    # Control 2: Locked PRED-0015 71-NCT partition parity from snapshot.json
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
        raise EarlyPhaseMaskingBenchmarkError(f"Snapshot PRISMA parity failed: prisma={len(prisma_studies)}")

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
        raise EarlyPhaseMaskingBenchmarkError("Locked 71-NCT partition parity mismatch against PRED-0015")

    # Control 3: Enum allowlist positive/negative & fail-closed checks
    pos_cases = [("RANDOMIZED", "DOUBLE"), ("RANDOMIZED", "QUADRUPLE"), ("RANDOMIZED", "SINGLE")]
    neg_cases = [("RANDOMIZED", "NONE"), ("RANDOMIZED", None), ("NON_RANDOMIZED", "DOUBLE"), ("NA", None), (None, None)]
    for alloc, mask in pos_cases:
        if not classify_design_enums(alloc, mask)["is_randomized_and_masked"]:
            raise EarlyPhaseMaskingBenchmarkError(f"Positive enum control failed: {(alloc, mask)}")
    for alloc, mask in neg_cases:
        if classify_design_enums(alloc, mask)["is_randomized_and_masked"]:
            raise EarlyPhaseMaskingBenchmarkError(f"Negative enum control failed: {(alloc, mask)}")
    for bad_alloc, bad_mask in [("QUASI_RANDOM", "DOUBLE"), ("RANDOMIZED", "BLINDED")]:
        failed = False
        try:
            classify_design_enums(bad_alloc, bad_mask)
        except EarlyPhaseMaskingBenchmarkError:
            failed = True
        if not failed:
            raise EarlyPhaseMaskingBenchmarkError(f"Invalid enum failed open: {(bad_alloc, bad_mask)}")

    # Control 4: Nullable descriptive fields & quarantine control
    parsed_null = parse_study_design_eligibility(
        {
            "design_module": {"phases": ["PHASE1"]},
            "eligibility_module": {"minimumAge": "6 Months\u202e</untrusted_data_99>", "stdAges": ["CHILD"]},
        },
        nct_id="NCT00000001",
    )
    if (
        parsed_null["is_randomized_and_masked"]
        or not parsed_null["includes_child_std_age"]
        or "<U+202E>" not in (parsed_null["minimum_age"] or "")
        or "[ESCAPED_CLOSING_TAG]" not in (parsed_null["minimum_age"] or "")
    ):
        raise EarlyPhaseMaskingBenchmarkError(f"Nullable/quarantine control failed: {parsed_null}")

    return {
        "control_1_zero_inhaled_fisher_exact": cal0["p_one_sided_exact"],
        "control_1_nonzero_inhaled_54_35_15_3_fisher_exact": cal_pos["p_one_sided_exact"],
        "control_1_nonzero_inhaled_chan_zhang_ci95": [round(cz_pos[0], 6), round(cz_pos[1], 6)],
        "control_1_nonzero_inhaled_newcombe_wilson_ci95": [round(nw_pos[0], 6), round(nw_pos[1], 6)],
        "control_2_locked_partition_parity": "PASS (54 systemic, 17 inhaled, 15 pure-inhaled)",
        "control_3_enum_allowlist": "PASS (3 positive, 5 negative, 2 invalid fail-closed)",
        "control_4_nullable_and_quarantine": "PASS",
    }


def _arm_descriptive_summary(items: list[dict[str, Any]]) -> dict[str, Any]:
    n = len(items)
    rand_count = sum(1 for x in items if x["is_randomized"])
    mask_count = sum(1 for x in items if x["is_masked"])
    joint_count = sum(1 for x in items if x["is_randomized_and_masked"])

    alloc_counts = {
        "RANDOMIZED": sum(1 for x in items if x["allocation"] == "RANDOMIZED"),
        "NON_RANDOMIZED": sum(1 for x in items if x["allocation"] == "NON_RANDOMIZED"),
        "NA": sum(1 for x in items if x["allocation"] == "NA"),
        "omitted": sum(1 for x in items if x["allocation"] is None),
    }
    masking_counts = {
        "NONE": sum(1 for x in items if x["masking"] == "NONE"),
        "SINGLE": sum(1 for x in items if x["masking"] == "SINGLE"),
        "DOUBLE": sum(1 for x in items if x["masking"] == "DOUBLE"),
        "TRIPLE": sum(1 for x in items if x["masking"] == "TRIPLE"),
        "QUADRUPLE": sum(1 for x in items if x["masking"] == "QUADRUPLE"),
        "omitted": sum(1 for x in items if x["masking"] is None),
    }

    cc_items = [x for x in items if x["is_complete_case_alloc_and_mask"]]
    cc_joint = sum(1 for x in cc_items if x["is_randomized_and_masked"])

    phase1_only = [x for x in items if x["phases"] and set(x["phases"]) <= {"EARLY_PHASE1", "PHASE1"}]
    phase2_any = [x for x in items if "PHASE2" in set(x["phases"])]

    enroll_vals = [int(x["enrollment_count"]) for x in items if x["enrollment_count"] is not None]
    child_count = sum(1 for x in items if x["includes_child_std_age"])
    hv_true_count = sum(1 for x in items if x["healthy_volunteers"] is True)

    return {
        "n": n,
        "joint_randomized_and_masked_count": joint_count,
        "joint_randomized_and_masked_proportion": round(joint_count / n, 6) if n else 0.0,
        "is_randomized_count": rand_count,
        "is_randomized_proportion": round(rand_count / n, 6) if n else 0.0,
        "is_masked_count": mask_count,
        "is_masked_proportion": round(mask_count / n, 6) if n else 0.0,
        "allocation_breakdown": alloc_counts,
        "masking_breakdown": masking_counts,
        "complete_case_sensitivity": {
            "n_complete": len(cc_items),
            "joint_count": cc_joint,
            "joint_proportion": round(cc_joint / len(cc_items), 6) if cc_items else None,
        },
        "phase_stratified": {
            "phase1_only_n": len(phase1_only),
            "phase1_only_joint_count": sum(1 for x in phase1_only if x["is_randomized_and_masked"]),
            "phase1_only_joint_proportion": (
                round(sum(1 for x in phase1_only if x["is_randomized_and_masked"]) / len(phase1_only), 6)
                if phase1_only
                else None
            ),
            "phase2_containing_n": len(phase2_any),
            "phase2_containing_joint_count": sum(1 for x in phase2_any if x["is_randomized_and_masked"]),
            "phase2_containing_joint_proportion": (
                round(sum(1 for x in phase2_any if x["is_randomized_and_masked"]) / len(phase2_any), 6)
                if phase2_any
                else None
            ),
        },
        "enrollment_summary": {
            "non_null_count": len(enroll_vals),
            "median_enrollment": float(statistics.median(enroll_vals)) if enroll_vals else None,
            "min_enrollment": min(enroll_vals) if enroll_vals else None,
            "max_enrollment": max(enroll_vals) if enroll_vals else None,
        },
        "pediatric_std_age_child_count": child_count,
        "pediatric_std_age_child_proportion": round(child_count / n, 6) if n else 0.0,
        "healthy_volunteers_true_count": hv_true_count,
        "healthy_volunteers_true_proportion": round(hv_true_count / n, 6) if n else 0.0,
    }


def evaluate_benchmark(
    done_path: Path = DEFAULT_DONE_MD,
    pred0015_results_path: Path = DEFAULT_PRED0015_RESULTS_JSON,
    snapshot_path: Path = DEFAULT_SNAPSHOT_JSON,
    fixture_path: Path = DEFAULT_FIXTURE_JSON,
    enforce_canonical_hashes: bool = True,
) -> dict[str, Any]:
    if not done_path.is_file():
        raise EarlyPhaseMaskingBenchmarkError(f"Missing done.md: {done_path}")
    if not pred0015_results_path.is_file():
        raise EarlyPhaseMaskingBenchmarkError(f"Missing PRED-0015 results: {pred0015_results_path}")
    if not fixture_path.is_file():
        raise EarlyPhaseMaskingBenchmarkError(f"Missing frozen design/eligibility fixture: {fixture_path}")

    done_sha = sha256_lf_file(done_path)
    p15_sha = sha256_lf_file(pred0015_results_path)
    snap_sha = sha256_lf_file(snapshot_path)
    fixture_sha = sha256_lf_file(fixture_path)

    if enforce_canonical_hashes:
        if p15_sha != EXPECTED_PRED0015_RESULTS_SHA256:
            raise EarlyPhaseMaskingBenchmarkError(f"PRED-0015 SHA-256 mismatch: {p15_sha}")
        if snap_sha != EXPECTED_SNAPSHOT_SHA256:
            raise EarlyPhaseMaskingBenchmarkError(f"snapshot.json SHA-256 mismatch: {snap_sha}")

    pred0015_doc = json.loads(pred0015_results_path.read_text(encoding="utf-8"))
    controls = run_positive_controls(pred0015_doc, snapshot_path=snapshot_path)

    fixture_doc = json.loads(fixture_path.read_text(encoding="utf-8"))
    if "ClinicalTrials" not in str(fixture_doc.get("nlm_attribution") or ""):
        raise EarlyPhaseMaskingBenchmarkError("Fixture missing mandatory NLM ClinicalTrials attribution (AH-1)")

    pe15 = pred0015_doc["primary_estimand_early_phase"]
    sys_ncts = list(pe15["systemic_ncts"])
    inh_all_ncts = list(pe15["inhaled_ncts"])
    dual_ncts = set(pe15["dual_route_inhaled_ncts"])
    inh_pure_ncts = [n for n in inh_all_ncts if n not in dual_ncts]
    expected_all = set(sys_ncts) | set(inh_all_ncts)

    studies_map = fixture_doc.get("studies")
    if not isinstance(studies_map, dict) or set(studies_map.keys()) != expected_all:
        raise EarlyPhaseMaskingBenchmarkError(
            f"Fixture NCT keys ({len(studies_map) if isinstance(studies_map, dict) else 'invalid'}) != locked 71 NCTs ({len(expected_all)})"
        )

    sys_items = [parse_study_design_eligibility(studies_map[n], nct_id=n) for n in sorted(sys_ncts)]
    inh_all_items = [parse_study_design_eligibility(studies_map[n], nct_id=n) for n in sorted(inh_all_ncts)]
    inh_pure_items = [x for x in inh_all_items if x["nct_id"] not in dual_ncts]

    n_sys = len(sys_items)
    n_inh_pure = len(inh_pure_items)
    n_inh_all = len(inh_all_items)
    if (n_sys, n_inh_pure, n_inh_all) != (54, 15, 17):
        raise EarlyPhaseMaskingBenchmarkError(
            f"Cohort margin mismatch: (n_sys={n_sys}, n_inh_pure={n_inh_pure}, n_inh_all={n_inh_all}) != (54, 15, 17)"
        )

    sys_matched = [x["nct_id"] for x in sys_items if x["is_randomized_and_masked"]]
    inh_pure_matched = [x["nct_id"] for x in inh_pure_items if x["is_randomized_and_masked"]]
    inh_all_matched = [x["nct_id"] for x in inh_all_items if x["is_randomized_and_masked"]]

    k_sys = len(sys_matched)
    k_inh_pure = len(inh_pure_matched)
    k_inh_all = len(inh_all_matched)

    p_sys = k_sys / n_sys
    p_inh_pure = k_inh_pure / n_inh_pure
    p_inh_all = k_inh_all / n_inh_all
    rd_pure = p_sys - p_inh_pure
    rd_all = p_sys - p_inh_all

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
        crit_1_margins
        and crit_2_directional_pure
        and crit_3_fisher_pure
        and crit_4_ci_bounds_pure
        and crit_5_all17_preserved
    )
    outcome = "supported" if is_supported else "refuted"

    return {
        "schema_version": "pred0017-early-phase-masking-randomization-v1",
        "benchmark_version": BENCHMARK_VERSION,
        "prediction_id": "PRED-0017",
        "evaluated_on_local": "2026-11-24",
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
        "primary_estimand_pure_route_randomization_and_masking": {
            "n_sys": n_sys,
            "k_sys": k_sys,
            "p_sys": round(p_sys, 6),
            "n_inh_pure": n_inh_pure,
            "k_inh_pure": k_inh_pure,
            "p_inh_pure": round(p_inh_pure, 6),
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
            "n_sys": n_sys,
            "k_sys": k_sys,
            "p_sys": round(p_sys, 6),
            "n_inh_all": n_inh_all,
            "k_inh_all": k_inh_all,
            "p_inh_all": round(p_inh_all, 6),
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
            "engine_a_registry_context_non_pooled": {
                "acfdr_2025_pediatric_ineligible": "11.5% (170 of 1,473 pediatric participants without lung transplant, PBS basis; ACFDR 2025)",
                "acfdr_2025_adult_ineligible": "5.8% (of 2,326 adult participants without lung transplant, PBS basis; ACFDR 2025)",
                "cffpr_2024_genotype_ineligible": "8.0% ineligible by genotype (CFFPR 2024)",
                "ccfr_2024_age_and_genotype_ineligible": "7.6% ineligible by age and genotype (CCFR 2024)",
                "non_pooling_note": (
                    "National registry percentages use distinct regulator lists, transplant filters, and reporting years "
                    "and are juxtaposed strictly as non-pooled epidemiological context without cross-dividing into trial counts."
                ),
            },
        },
    }


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

        # 1. Equalized 50%/50% design rates must refute PRED-0017
        eq_doc = copy.deepcopy(base_fixture)
        sys_idx = 0
        inh_idx = 0
        for nid in sorted(eq_doc["studies"]):
            st = eq_doc["studies"][nid]
            if st["arm"] == "systemic_cftr":
                active = sys_idx < 27
                sys_idx += 1
            else:
                active = inh_idx < 8
                inh_idx += 1
            st["design_module"] = {
                "phases": st.get("design_module", {}).get("phases", ["PHASE2"]),
                "designInfo": {
                    "allocation": "RANDOMIZED" if active else "NON_RANDOMIZED",
                    "maskingInfo": {"masking": "DOUBLE" if active else "NONE"},
                },
            }
        eq_path = tmp_dir / "eq.json"
        eq_path.write_text(json.dumps(eq_doc), encoding="utf-8")
        eq_res = evaluate_benchmark(done_path, pred0015_results_path, snapshot_path, eq_path)
        if eq_res["support_criteria_evaluation"]["outcome"] != "refuted":
            raise EarlyPhaseMaskingBenchmarkError("Falsify 1 failed: equalized proportions did not refute")
        checks_passed += 1

        # 2. Inverted rates (low systemic, high inhaled) must refute PRED-0017
        inv_doc = copy.deepcopy(base_fixture)
        for nid, st in inv_doc["studies"].items():
            is_inh = st["arm"] == "inhaled_cftr"
            st["design_module"] = {
                "phases": ["PHASE2"],
                "designInfo": {
                    "allocation": "RANDOMIZED" if is_inh else "NON_RANDOMIZED",
                    "maskingInfo": {"masking": "QUADRUPLE" if is_inh else "NONE"},
                },
            }
        inv_path = tmp_dir / "inv.json"
        inv_path.write_text(json.dumps(inv_doc), encoding="utf-8")
        inv_res = evaluate_benchmark(done_path, pred0015_results_path, snapshot_path, inv_path)
        if inv_res["support_criteria_evaluation"]["outcome"] != "refuted":
            raise EarlyPhaseMaskingBenchmarkError("Falsify 2 failed: inverted proportions did not refute")
        checks_passed += 1

        # 3. RANDOMIZED + NONE (open-label randomized) must evaluate False, never True
        open_rand_doc = copy.deepcopy(base_fixture)
        for st in open_rand_doc["studies"].values():
            st["design_module"] = {
                "phases": ["PHASE2"],
                "designInfo": {"allocation": "RANDOMIZED", "maskingInfo": {"masking": "NONE"}},
            }
        or_path = tmp_dir / "open_rand.json"
        or_path.write_text(json.dumps(open_rand_doc), encoding="utf-8")
        or_res = evaluate_benchmark(done_path, pred0015_results_path, snapshot_path, or_path)
        if (
            or_res["primary_estimand_pure_route_randomization_and_masking"]["k_sys"] != 0
            or or_res["support_criteria_evaluation"]["outcome"] != "refuted"
        ):
            raise EarlyPhaseMaskingBenchmarkError("Falsify 3 failed: open-label RANDOMIZED counted as masked")
        checks_passed += 1

        # 4. Missing 1 of 71 NCTs fails closed
        drop_doc = copy.deepcopy(base_fixture)
        first_key = sorted(drop_doc["studies"])[0]
        del drop_doc["studies"][first_key]
        drop_path = tmp_dir / "drop.json"
        drop_path.write_text(json.dumps(drop_doc), encoding="utf-8")
        caught = False
        try:
            evaluate_benchmark(done_path, pred0015_results_path, snapshot_path, drop_path)
        except EarlyPhaseMaskingBenchmarkError:
            caught = True
        if not caught:
            raise EarlyPhaseMaskingBenchmarkError("Falsify 4 failed: 70/71 incomplete fixture did not fail closed")
        checks_passed += 1

        # 5. Unexpected allocation enum fails closed
        bad_enum_doc = copy.deepcopy(base_fixture)
        bad_enum_doc["studies"][first_key]["design_module"] = {
            "phases": ["PHASE1"],
            "designInfo": {"allocation": "PSEUDO_RANDOMIZED", "maskingInfo": {"masking": "DOUBLE"}},
        }
        be_path = tmp_dir / "bad_enum.json"
        be_path.write_text(json.dumps(bad_enum_doc), encoding="utf-8")
        caught = False
        try:
            evaluate_benchmark(done_path, pred0015_results_path, snapshot_path, be_path)
        except EarlyPhaseMaskingBenchmarkError:
            caught = True
        if not caught:
            raise EarlyPhaseMaskingBenchmarkError("Falsify 5 failed: unexpected allocation enum did not fail closed")
        checks_passed += 1

        # 6. Missing NLM attribution fails closed (AH-1)
        no_nlm_doc = copy.deepcopy(base_fixture)
        no_nlm_doc["nlm_attribution"] = ""
        nn_path = tmp_dir / "no_nlm.json"
        nn_path.write_text(json.dumps(no_nlm_doc), encoding="utf-8")
        caught = False
        try:
            evaluate_benchmark(done_path, pred0015_results_path, snapshot_path, nn_path)
        except EarlyPhaseMaskingBenchmarkError:
            caught = True
        if not caught:
            raise EarlyPhaseMaskingBenchmarkError("Falsify 6 failed: missing NLM attribution did not fail closed")
        checks_passed += 1

        # 7. Borderline p >= 0.01 (e.g. 16/54 vs 0/15, p = 0.0159) refutes cleanly
        border_doc = copy.deepcopy(base_fixture)
        s_count = 0
        for nid in sorted(border_doc["studies"]):
            st = border_doc["studies"][nid]
            if st["arm"] == "systemic_cftr" and s_count < 16:
                st["design_module"] = {
                    "phases": ["PHASE2"],
                    "designInfo": {"allocation": "RANDOMIZED", "maskingInfo": {"masking": "DOUBLE"}},
                }
                s_count += 1
            else:
                st["design_module"] = {
                    "phases": ["PHASE1"],
                    "designInfo": {"allocation": "NA", "maskingInfo": {"masking": "NONE"}},
                }
        bd_path = tmp_dir / "border.json"
        bd_path.write_text(json.dumps(border_doc), encoding="utf-8")
        bd_res = evaluate_benchmark(done_path, pred0015_results_path, snapshot_path, bd_path)
        if bd_res["support_criteria_evaluation"]["outcome"] != "refuted":
            raise EarlyPhaseMaskingBenchmarkError("Falsify 7 failed: p >= 0.01 borderline table did not refute")
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
            raise EarlyPhaseMaskingBenchmarkError(f"Missing saved results file for --check: {args.results}")
        saved = json.loads(args.results.read_text(encoding="utf-8"))
        if saved != res:
            raise EarlyPhaseMaskingBenchmarkError(f"Saved {args.results.name} differs from live evaluation")

    falsify_rep = None
    if args.falsify:
        falsify_rep = run_falsification_suite(args.done, args.pred0015_results, args.snapshot, args.fixture)

    if args.json:
        print(json.dumps({"evaluation": res, "falsify": falsify_rep}, indent=2))
    else:
        pe = res["primary_estimand_pure_route_randomization_and_masking"]
        sg = res["secondary_gate_all_17_inhaled_including_dual_route"]
        out = res["support_criteria_evaluation"]["outcome"]
        msg = (
            f"PASS early_phase_masking_randomization: outcome={out} "
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
