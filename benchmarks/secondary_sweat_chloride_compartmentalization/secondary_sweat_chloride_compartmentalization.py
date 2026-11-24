#!/usr/bin/env python3
"""Full-Protocol (`primaryOutcomes U secondaryOutcomes`) Sweat-Chloride Compartmentalization Benchmark (`PRED-0016`).

Tests the pre-registered compartmentalization hypothesis (`benchmarks/secondary_sweat_chloride_compartmentalization/done.md`)
across the locked 71 early-phase non-healthy-volunteer PRISMA-940 Cystic Fibrosis interventional trials
from `PRED-0015` (`benchmarks/delivery_route_compartmentalization/results_2026-11-24.json`):

- Primary comparison: `systemic_cftr` (`N_sys == 54`) vs pure airway-restricted `inhaled_cftr_pure`
  (`N_inh_pure == 15`, excluding dual-route background systemic modulator trials `NCT05668741` and `NCT06237335`),
  with all-17 `inhaled_cftr` (`N_inh_all == 17`) evaluated alongside as a mandatory secondary gate.
- Outcome indicator (`has_registered_sweat_chloride`): Evaluated over `primaryOutcomes[].measure U secondaryOutcomes[].measure`
  from frozen ClinicalTrials.gov API v2 `early_phase_outcomes_2026-11-24.json` (avoiding primary-promotion censoring
  of systemic trials that already list primary sweat chloride, while reporting secondary-only and `measure + description`
  as non-gated descriptive context).

Invariants (`cf-research-context`, `discovery-practice`, `security-browsing`, `AH-1`..`AH-4`):
- 100% pure Python 3.10+ standard library + local `tools/trial_atlas/` modules; zero network/LLM calls.
- Strict `AH-1` verbatim NLM attribution and `AH-2` PRISMA-940 vs raw-1053 cohort boundary.
- Omitted `secondaryOutcomes` mapped to `[]` (retained in arm denominator) with per-arm missingness reported.
- CodeQL-safe identifiers (`RESEARCH_BOUNDARY`, `NLM_ATTRIBUTION`).
"""
from __future__ import annotations

import argparse
import copy
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import re
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
        sys.path.insert(0, p_dir)

import delivery_route_compartmentalization as p15_mod  # noqa: E402
import scope as scope_mod  # noqa: E402

BENCHMARK_VERSION = "2026-11-24.v1"
EXPECTED_PRED0015_RESULTS_SHA256 = "c4288a7c40783c746e1e4c103e9b1b8867818cb97cc371d19348715d9251dcb1"
EXPECTED_SNAPSHOT_SHA256 = "de0f330c4441d77b84dc5f547e4dfcf99195c1b9f223e4dbeebd628054925b2d"

DEFAULT_DONE_MD = HERE / "done.md"
DEFAULT_PRED0015_RESULTS_JSON = PRED0015_DIR / "results_2026-11-24.json"
DEFAULT_SNAPSHOT_JSON = TRIAL_ATLAS_DIR / "snapshot.json"
DEFAULT_OUTCOMES_JSON = HERE / "early_phase_outcomes_2026-11-24.json"
DEFAULT_RESULTS_JSON = HERE / "results_2026-11-24.json"

NLM_ATTRIBUTION = (
    "ClinicalTrials.gov does not evaluate safety or efficacy; "
    "listing on the registry does not denote government endorsement."
)
RESEARCH_BOUNDARY = (
    "Independent open-source research benchmark only (Tier T1/T2 registry metadata); "
    "not clinical guidance and not an assessment of personal trial eligibility or treatment efficacy."
)

UNICODE_DASH_MAP = str.maketrans({"\u2212": "-", "\u2013": "-", "\u2014": "-"})
_INTRA_CLAUSE_SEP = r"[\s\-\[\]\(\)]+"
SWEAT_CHLORIDE_FULL_RE = re.compile(
    rf"\bsweat{_INTRA_CLAUSE_SEP}(?:chloride|cl)\b"
    rf"|\bchloride{_INTRA_CLAUSE_SEP}concentration{_INTRA_CLAUSE_SEP}in{_INTRA_CLAUSE_SEP}sweat\b"
    rf"|\bchloride{_INTRA_CLAUSE_SEP}in{_INTRA_CLAUSE_SEP}sweat\b"
    rf"|\bsweat{_INTRA_CLAUSE_SEP}electrolytes?\b",
    re.I,
)
DIAGNOSTIC_EXCLUSION_RE = re.compile(
    rf"\b(?:diagnosis|diagnostic|inclusion|eligibility|screening{_INTRA_CLAUSE_SEP}criteria)\b",
    re.I,
)
TROJAN_BIDI_CODEPOINTS: frozenset[int] = frozenset(
    {0x200B, 0x200C, 0x200D, 0x200E, 0x200F, 0x202A, 0x202B, 0x202C, 0x202D, 0x202E, 0x2060, 0x2066, 0x2067, 0x2068, 0x2069, 0xFEFF}
)
CONFUSABLE_ASCII_MAP = str.maketrans(
    {
        "\u0430": "a", "\u0435": "e", "\u043e": "o", "\u0440": "p", "\u0441": "c", "\u0443": "y",
        "\u0445": "x", "\u0501": "d", "\u0456": "i", "\u0455": "s", "\u0442": "t", "\u043d": "n",
        "\u0410": "A", "\u0412": "B", "\u0421": "C", "\u0500": "D", "\u0415": "E", "\u041d": "H",
        "\u0406": "I", "\u041a": "K", "\u041c": "M", "\u041e": "O", "\u0420": "P", "\u0405": "S",
        "\u0422": "T", "\u0425": "X", "\u0423": "Y",
        "\u03b1": "a", "\u03b5": "e", "\u03b9": "i", "\u03bf": "o", "\u03c1": "p", "\u03c4": "t",
        "\u03c5": "u", "\u03bd": "n", "\u03ba": "k", "\u03c7": "x", "\u03c2": "s", "\u03c3": "s",
        "\u0391": "A", "\u0392": "B", "\u0395": "E", "\u0397": "H", "\u0399": "I", "\u039a": "K",
        "\u039c": "M", "\u039d": "N", "\u039f": "O", "\u03a1": "P", "\u03a4": "T", "\u03a7": "X", "\u03a5": "Y",
        "\uFF1C": "<", "\uFF1E": ">", "\uFF0F": "/", "\u2044": "/", "\u2215": "/",
        "\u2329": "<", "\u232A": ">", "\u27E8": "<", "\u27E9": ">",
        "\u2039": "<", "\u203A": ">", "\u3008": "<", "\u3009": ">",
        "\uFE64": "<", "\uFE65": ">", "\u02C2": "<", "\u02C3": ">",
    }
)
UNTRUSTED_CLOSE_LOOKALIKE_RE = re.compile(r"<\s*/\s*[^\s>]*data|/\s*untrusted", re.I)
REQUIRED_OUTCOME_KEYS: frozenset[str] = frozenset({"measure", "description", "time_frame"})
Z_95 = 1.959963984540054


class FullProtocolCompartmentalizationError(RuntimeError):
    """Raised when schema, provenance, positive-control, or falsification checks fail."""


def sha256_lf_text(path: Path) -> str:
    text = path.read_text(encoding="utf-8").replace("\r\n", "\n")
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def sanitize_untrusted_text(raw: Any) -> str:
    if raw is None:
        return ""
    s = str(raw)
    stripped_for_probe = "".join(
        ch
        for ch in s
        if ord(ch) not in TROJAN_BIDI_CODEPOINTS and unicodedata.category(ch) != "Cf"
    )
    folded = unicodedata.normalize("NFKD", stripped_for_probe).translate(CONFUSABLE_ASCII_MAP)
    is_close_lookalike = bool(UNTRUSTED_CLOSE_LOOKALIKE_RE.search(folded))
    out_chars: list[str] = []
    for ch in s:
        cp = ord(ch)
        if cp in TROJAN_BIDI_CODEPOINTS or unicodedata.category(ch) == "Cf":
            out_chars.append(f"<U+{cp:04X}>")
            continue
        if is_close_lookalike:
            mapped = unicodedata.normalize("NFKD", ch).translate(CONFUSABLE_ASCII_MAP)
            if "<" in mapped:
                out_chars.append("&lt;")
                continue
            if ">" in mapped:
                out_chars.append("&gt;")
                continue
            if "/" in mapped:
                out_chars.append("&#47;")
                continue
        out_chars.append(ch)
    return "".join(out_chars)


def matches_sweat_chloride_measure(text: str) -> bool:
    norm = (text or "").translate(UNICODE_DASH_MAP).strip()
    if not norm:
        return False
    return bool(SWEAT_CHLORIDE_FULL_RE.search(norm) and not DIAGNOSTIC_EXCLUSION_RE.search(norm))


def normalize_outcomes_array(raw_items: Any, nct_id: str, field_name: str) -> list[dict[str, str]]:
    """Map omitted/None outcomes list to `[]`; reject non-list or malformed entries."""
    if raw_items is None:
        return []
    if not isinstance(raw_items, list):
        raise FullProtocolCompartmentalizationError(
            f"{nct_id}: {field_name} must be a list or omitted/None, got {type(raw_items).__name__}"
        )
    out: list[dict[str, str]] = []
    for idx, entry in enumerate(raw_items):
        if not isinstance(entry, dict) or not REQUIRED_OUTCOME_KEYS.issubset(entry.keys()):
            raise FullProtocolCompartmentalizationError(
                f"{nct_id}: {field_name}[{idx}] missing required keys {sorted(REQUIRED_OUTCOME_KEYS)}"
            )
        out.append(
            {
                "measure": sanitize_untrusted_text(entry.get("measure") or ""),
                "description": sanitize_untrusted_text(entry.get("description") or ""),
                "time_frame": sanitize_untrusted_text(entry.get("time_frame") or ""),
            }
        )
    return out


def fisher_hypergeometric_exact(n1: int, k1: int, n2: int, k2: int) -> dict[str, Any]:
    if n1 <= 0 or n2 <= 0 or not (0 <= k1 <= n1) or not (0 <= k2 <= n2):
        raise ValueError(f"Invalid 2x2 counts: n1={n1}, k1={k1}, n2={n2}, k2={k2}")
    total_n = n1 + n2
    total_k = k1 + k2
    denom = math.comb(total_n, total_k)
    j_min = max(0, total_k - n2)
    j_max = min(n1, total_k)
    pmf = {j: Fraction(math.comb(n1, j) * math.comb(n2, total_k - j), denom) for j in range(j_min, j_max + 1)}
    support_sum = sum(pmf.values(), Fraction(0, 1))
    if support_sum != Fraction(1, 1):
        raise FullProtocolCompartmentalizationError(f"Hypergeometric PMF sum != 1: {support_sum}")
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


def chan_zhang_unconditional_rd_ci95(
    n1: int, k1: int, n2: int, k2: int, alpha_two_sided: float = 0.05
) -> tuple[float, float]:
    alpha_tail = 0.5 * alpha_two_sided
    comb1 = [math.comb(n1, i) for i in range(n1 + 1)]
    comb2 = [math.comb(n2, j) for j in range(n2 + 1)]

    def binom_pmf(n: int, combs: list[int], p: float) -> list[float]:
        if p <= 0.0:
            return [1.0] + [0.0] * n
        if p >= 1.0:
            return [0.0] * n + [1.0]
        q = 1.0 - p
        return [combs[i] * (p ** i) * (q ** (n - i)) for i in range(n + 1)]

    def tail_supremum(delta: float, upper: bool) -> float:
        low_p2 = max(0.0, -delta)
        high_p2 = min(1.0, 1.0 - delta)
        z_obs = _mn_score_z(k1, n1, k2, n2, delta)
        active_pairs: list[tuple[int, int]] = []
        for x1 in range(n1 + 1):
            for x2 in range(n2 + 1):
                z_xy = _mn_score_z(x1, n1, x2, n2, delta)
                if (upper and z_xy >= z_obs - 1e-11) or ((not upper) and z_xy <= z_obs + 1e-11):
                    active_pairs.append((x1, x2))
        if not active_pairs:
            return 0.0
        if len(active_pairs) == (n1 + 1) * (n2 + 1):
            return 1.0

        def prob_at(p2: float) -> float:
            p1 = min(1.0, max(0.0, p2 + delta))
            b1 = binom_pmf(n1, comb1, p1)
            b2 = binom_pmf(n2, comb2, p2)
            return sum(b1[x1] * b2[x2] for x1, x2 in active_pairs)

        if high_p2 <= low_p2 + 1e-14:
            return prob_at(0.5 * (low_p2 + high_p2))

        grid_n = 240
        best_p2 = low_p2
        best_val = -1.0
        for idx in range(grid_n + 1):
            p2 = low_p2 + (high_p2 - low_p2) * (idx / grid_n)
            val = prob_at(p2)
            if val > best_val:
                best_val = val
                best_p2 = p2
        step = (high_p2 - low_p2) / grid_n
        a = max(low_p2, best_p2 - step)
        b = min(high_p2, best_p2 + step)
        inv_phi = (math.sqrt(5.0) - 1.0) / 2.0
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
    if k1 == 0 and k2 == n2:
        low_ci = -1.0
    else:
        lo, hi = -1.0 + 1e-12, rd_obs
        if tail_supremum(lo, True) > alpha_tail:
            low_ci = -1.0
        else:
            for _ in range(52):
                mid = 0.5 * (lo + hi)
                if tail_supremum(mid, True) > alpha_tail:
                    hi = mid
                else:
                    lo = mid
            low_ci = hi

    if k1 == n1 and k2 == 0:
        high_ci = 1.0
    else:
        lo, hi = rd_obs, 1.0 - 1e-12
        if tail_supremum(hi, False) > alpha_tail:
            high_ci = 1.0
        else:
            for _ in range(52):
                mid = 0.5 * (lo + hi)
                if tail_supremum(mid, False) > alpha_tail:
                    lo = mid
                else:
                    hi = mid
            high_ci = lo
    return low_ci, high_ci


def run_positive_controls(
    pred0015_doc: dict[str, Any],
    snapshot_path: Path = DEFAULT_SNAPSHOT_JSON,
) -> dict[str, Any]:
    """Run all four pre-registered positive controls before reading empirical secondaryOutcomes."""
    # Control 1: Hypergeometric & CI calibration
    cal = fisher_hypergeometric_exact(12, 9, 10, 0)
    if cal["p_one_sided_fraction"] != Fraction(11, 24871):
        raise FullProtocolCompartmentalizationError(f"Calibration Fisher mismatch: {cal['p_one_sided_exact']}")
    cz_cal = chan_zhang_unconditional_rd_ci95(12, 9, 10, 0)
    nw_cal = newcombe_wilson_rd_ci95(12, 9, 10, 0)
    if cz_cal[0] <= 0.35 or nw_cal[0] <= 0.35:
        raise FullProtocolCompartmentalizationError(f"Calibration CI95 <= 0.35: cz={cz_cal}, nw={nw_cal}")

    b17 = fisher_hypergeometric_exact(54, 17, 15, 0)["p_one_sided_fraction"]
    b16 = fisher_hypergeometric_exact(54, 16, 15, 0)["p_one_sided_fraction"]
    if not (b17 < Fraction(1, 100) < b16):
        raise FullProtocolCompartmentalizationError(f"(54,15) alpha=0.01 boundary failed: p17={b17}, p16={b16}")

    # Control 2: Locked PRED-0015 71-NCT partition & primary baseline parity from snapshot.json
    snap_doc = json.loads(snapshot_path.read_text(encoding="utf-8"))
    snap_meta = snap_doc.get("metadata") or {}
    raw_studies = snap_doc.get("studies") or []
    prisma_studies: list[dict[str, Any]] = []
    for st in raw_studies:
        rec_st = p15_mod.snap.study_record(p15_mod._to_raw_study(st))
        rule, _ = scope_mod.classify(rec_st)
        if rule is None:
            p15_mod.assert_flattened_study_schema(st)
            prisma_studies.append(st)
    if len(raw_studies) != 1053 or len(prisma_studies) != 940:
        raise FullProtocolCompartmentalizationError(
            f"Snapshot PRISMA parity failed: raw={len(raw_studies)}, prisma={len(prisma_studies)}"
        )

    rederived_sys: list[str] = []
    rederived_inh: list[str] = []
    rederived_dual: list[str] = []
    rederived_sys_prim_cl: list[str] = []
    rederived_inh_prim_cl: list[str] = []
    snap_primary_by_nct: dict[str, list[str]] = {}

    for st in prisma_studies:
        if not p15_mod.is_early_phase(st.get("phases")):
            continue
        if p15_mod.is_healthy_volunteer_only(st.get("brief_title", ""), st.get("official_title", "")):
            continue
        reg = p15_mod.classify_study_regimen(st)
        prim_out = p15_mod.classify_primary_outcomes(st)
        nid = str(st["nct_id"])
        if reg["arm"] == "systemic_cftr":
            rederived_sys.append(nid)
            snap_primary_by_nct[nid] = [str(o.get("measure") or "").strip() for o in st["primary_outcomes"]]
            if prim_out["has_primary_sweat_chloride"]:
                rederived_sys_prim_cl.append(nid)
        elif reg["arm"] == "inhaled_cftr":
            rederived_inh.append(nid)
            snap_primary_by_nct[nid] = [str(o.get("measure") or "").strip() for o in st["primary_outcomes"]]
            if reg["dual_route_modulator_co_admin"]:
                rederived_dual.append(nid)
            if prim_out["has_primary_sweat_chloride"]:
                rederived_inh_prim_cl.append(nid)

    pe15 = pred0015_doc["primary_estimand_early_phase"]
    if (
        sorted(rederived_sys) != sorted(pe15["systemic_ncts"])
        or sorted(rederived_inh) != sorted(pe15["inhaled_ncts"])
        or sorted(rederived_dual) != sorted(pe15["dual_route_inhaled_ncts"])
        or sorted(rederived_sys_prim_cl) != sorted(pe15["systemic_sweat_cl_ncts"])
        or rederived_inh_prim_cl != []
    ):
        raise FullProtocolCompartmentalizationError("Rederived 71-NCT partition diverged from locked PRED-0015 results")

    # Control 3: Lexical positive & negative synthetic outcome controls
    pos_phrases = (
        "Change from baseline in sweat chloride concentration",
        "Sweat [Cl\u2212] at Day 28",
        "Sweat electrolytes (chloride)",
    )
    neg_phrases = (
        "Diagnostic sweat chloride >= 60 mmol/L inclusion",
        "Sweat sodium concentration",
        "Sweat rate evaporimetry",
        "Percent predicted FEV1",
    )
    for phr in pos_phrases:
        if not matches_sweat_chloride_measure(phr):
            raise FullProtocolCompartmentalizationError(f"Lexical positive control failed: {phr!r}")
    for phr in neg_phrases:
        if matches_sweat_chloride_measure(phr):
            raise FullProtocolCompartmentalizationError(f"Lexical negative control failed: {phr!r}")

    # Control 4: Omitted secondaryOutcomes & untrusted_data quarantine control
    if normalize_outcomes_array(None, "NCT00000000", "secondaryOutcomes") != []:
        raise FullProtocolCompartmentalizationError("Omitted secondaryOutcomes did not map to []")
    escaped = sanitize_untrusted_text("test\u202E </untrusted_data>")
    if "<U+202E>" not in escaped or "</untrusted_data>" in escaped:
        raise FullProtocolCompartmentalizationError(f"Quarantine sanitizer failed: {escaped!r}")
    for adv_ctrl in ("\u27E8/untrusted_data\u27E9", "\uFF1C/untrusted_data\uFF1E", "</un\u200Btrusted_da\u200Bta>"):
        esc_adv = sanitize_untrusted_text(adv_ctrl)
        if "&lt;" not in esc_adv or "&#47;" not in esc_adv or "&gt;" not in esc_adv:
            raise FullProtocolCompartmentalizationError(f"Quarantine homoglyph/Cf sanitizer failed: {esc_adv!r}")

    return {
        "calibration_12_9_10_0": {
            "p_fisher_1s_exact": cal["p_one_sided_exact"],
            "chan_zhang_ci95_low": round(cz_cal[0], 6),
            "newcombe_wilson_ci95_low": round(nw_cal[0], 6),
        },
        "boundary_54_15_alpha_0_01": {
            "p_fisher_1s_k17_exact": f"{b17.numerator}/{b17.denominator}",
            "p_fisher_1s_k16_exact": f"{b16.numerator}/{b16.denominator}",
        },
        "locked_71_nct_parity_verified": True,
        "lexical_and_quarantine_controls_verified": True,
        "_snap_data_timestamp": str(snap_meta.get("data_timestamp") or ""),
        "_snap_primary_by_nct": snap_primary_by_nct,
    }


def evaluate_benchmark(
    done_path: Path = DEFAULT_DONE_MD,
    pred0015_results_path: Path = DEFAULT_PRED0015_RESULTS_JSON,
    snapshot_path: Path = DEFAULT_SNAPSHOT_JSON,
    outcomes_path: Path = DEFAULT_OUTCOMES_JSON,
    enforce_expected_hashes: bool = True,
) -> dict[str, Any]:
    done_sha = sha256_lf_text(done_path)
    p15_sha = sha256_lf_text(pred0015_results_path)
    snap_sha = sha256_lf_text(snapshot_path)
    outcomes_sha = sha256_lf_text(outcomes_path)

    if enforce_expected_hashes:
        if p15_sha != EXPECTED_PRED0015_RESULTS_SHA256:
            raise FullProtocolCompartmentalizationError(
                f"PRED-0015 results SHA-256 mismatch: {p15_sha} != {EXPECTED_PRED0015_RESULTS_SHA256}"
            )
        if snap_sha != EXPECTED_SNAPSHOT_SHA256:
            raise FullProtocolCompartmentalizationError(
                f"snapshot.json SHA-256 mismatch: {snap_sha} != {EXPECTED_SNAPSHOT_SHA256}"
            )

    p15_doc = json.loads(pred0015_results_path.read_text(encoding="utf-8"))
    pos_ctrl = run_positive_controls(p15_doc, snapshot_path=snapshot_path)
    snap_data_ts = pos_ctrl.pop("_snap_data_timestamp")
    snap_primary_by_nct: dict[str, list[str]] = pos_ctrl.pop("_snap_primary_by_nct")

    outcomes_doc = json.loads(outcomes_path.read_text(encoding="utf-8"))
    if outcomes_doc.get("nlm_attribution") != NLM_ATTRIBUTION or "ClinicalTrials" not in NLM_ATTRIBUTION:
        raise FullProtocolCompartmentalizationError("Missing verbatim NLM attribution in outcomes fixture")
    studies_map = outcomes_doc.get("studies")
    if not isinstance(studies_map, dict):
        raise FullProtocolCompartmentalizationError("outcomes fixture missing 'studies' dict")

    pe15 = p15_doc["primary_estimand_early_phase"]
    sys_ncts: list[str] = sorted(pe15["systemic_ncts"])
    inh_ncts: list[str] = sorted(pe15["inhaled_ncts"])
    dual_ncts: set[str] = set(pe15["dual_route_inhaled_ncts"])
    expected_71 = set(sys_ncts) | set(inh_ncts)
    if len(expected_71) != 71 or set(studies_map.keys()) != expected_71:
        raise FullProtocolCompartmentalizationError(
            f"outcomes fixture study keys ({len(studies_map)}) != locked 71 NCTs"
        )

    sys_any_ncts: list[str] = []
    sys_sec_only_ncts: list[str] = []
    sys_desc_any_ncts: list[str] = []
    sys_empty_sec = 0

    inh_all_any_ncts: list[str] = []
    inh_pure_any_ncts: list[str] = []
    inh_sec_only_ncts: list[str] = []
    inh_desc_any_ncts: list[str] = []
    inh_all_empty_sec = 0
    inh_pure_empty_sec = 0

    primary_drift_ncts: list[dict[str, Any]] = []

    for nct_id in sorted(expected_71):
        rec = studies_map[nct_id]
        if not isinstance(rec, dict) or rec.get("nct_id") != nct_id:
            raise FullProtocolCompartmentalizationError(f"Malformed study entry for {nct_id}")
        prim_list = normalize_outcomes_array(rec.get("primary_outcomes"), nct_id, "primary_outcomes")
        sec_list = normalize_outcomes_array(rec.get("secondary_outcomes"), nct_id, "secondary_outcomes")

        api_prim_measures = [o["measure"].strip() for o in prim_list]
        if api_prim_measures != snap_primary_by_nct.get(nct_id, []):
            primary_drift_ncts.append(
                {
                    "nct_id": nct_id,
                    "snapshot_primary_count": len(snap_primary_by_nct.get(nct_id, [])),
                    "api_primary_count": len(api_prim_measures),
                }
            )

        prim_match = any(matches_sweat_chloride_measure(o["measure"]) for o in prim_list)
        sec_match = any(matches_sweat_chloride_measure(o["measure"]) for o in sec_list)
        any_match = prim_match or sec_match
        desc_match = any(
            matches_sweat_chloride_measure(o["measure"])
            or matches_sweat_chloride_measure(f"{o['measure']} {o['description']}")
            for o in (prim_list + sec_list)
        )

        if nct_id in sys_ncts:
            if len(sec_list) == 0:
                sys_empty_sec += 1
            if any_match:
                sys_any_ncts.append(nct_id)
            if sec_match and not prim_match:
                sys_sec_only_ncts.append(nct_id)
            if desc_match:
                sys_desc_any_ncts.append(nct_id)
        else:
            is_dual = nct_id in dual_ncts
            if len(sec_list) == 0:
                inh_all_empty_sec += 1
                if not is_dual:
                    inh_pure_empty_sec += 1
            if any_match:
                inh_all_any_ncts.append(nct_id)
                if not is_dual:
                    inh_pure_any_ncts.append(nct_id)
            if sec_match and not prim_match:
                inh_sec_only_ncts.append(nct_id)
            if desc_match:
                inh_desc_any_ncts.append(nct_id)

    n_sys = len(sys_ncts)
    k_sys = len(sys_any_ncts)
    n_inh_pure = len(inh_ncts) - len(dual_ncts)
    k_inh_pure = len(inh_pure_any_ncts)
    n_inh_all = len(inh_ncts)
    k_inh_all = len(inh_all_any_ncts)

    fisher_pure = fisher_hypergeometric_exact(n_sys, k_sys, n_inh_pure, k_inh_pure)
    cz_pure = chan_zhang_unconditional_rd_ci95(n_sys, k_sys, n_inh_pure, k_inh_pure)
    nw_pure = newcombe_wilson_rd_ci95(n_sys, k_sys, n_inh_pure, k_inh_pure)

    fisher_all = fisher_hypergeometric_exact(n_sys, k_sys, n_inh_all, k_inh_all)
    cz_all = chan_zhang_unconditional_rd_ci95(n_sys, k_sys, n_inh_all, k_inh_all)
    nw_all = newcombe_wilson_rd_ci95(n_sys, k_sys, n_inh_all, k_inh_all)

    p_sys = k_sys / n_sys
    p_inh_pure = k_inh_pure / n_inh_pure
    p_inh_all = k_inh_all / n_inh_all

    c1_sizes = (n_sys == 54 and n_inh_pure == 15 and n_inh_all == 17)
    c2_dir = p_sys > p_inh_pure
    c3_fisher_pure = fisher_pure["p_one_sided_fraction"] < Fraction(1, 100)
    c4_cis_pure = (cz_pure[0] > 0.0 and nw_pure[0] > 0.0)
    c5_all17 = (p_sys > p_inh_all and fisher_all["p_one_sided_fraction"] < Fraction(1, 100))

    supported = bool(c1_sizes and c2_dir and c3_fisher_pure and c4_cis_pure and c5_all17)
    outcome = "supported" if supported else "refuted"
    ledger_standing = "open" if supported else "refuted"

    return {
        "benchmark_version": BENCHMARK_VERSION,
        "prediction_id": "PRED-0016",
        "evaluated_on_local": "2026-11-24",
        "done_sha256": done_sha,
        "pred0015_results_sha256": p15_sha,
        "snapshot_sha256": snap_sha,
        "outcomes_sha256": outcomes_sha,
        "nlm_attribution": NLM_ATTRIBUTION,
        "research_boundary": RESEARCH_BOUNDARY,
        "dataset": {
            "catalog_source_id": "clinicaltrials-gov",
            "snapshot_path": "tools/trial_atlas/snapshot.json",
            "snapshot_data_timestamp": snap_data_ts,
            "outcomes_fixture_path": "benchmarks/secondary_sweat_chloride_compartmentalization/early_phase_outcomes_2026-11-24.json",
            "outcomes_retrieved_at_utc": str(outcomes_doc.get("retrieved_at_utc") or ""),
            "raw_interventional_count": 1053,
            "prisma_included_count": 940,
            "locked_early_phase_total": 71,
            "primary_outcomes_drift_vs_snapshot_count": len(primary_drift_ncts),
            "primary_outcomes_drift_ncts": primary_drift_ncts,
            "empty_secondary_outcomes": {
                "systemic_empty_count": sys_empty_sec,
                "systemic_empty_rate": round(sys_empty_sec / n_sys, 6),
                "inhaled_pure_empty_count": inh_pure_empty_sec,
                "inhaled_pure_empty_rate": round(inh_pure_empty_sec / n_inh_pure, 6),
                "inhaled_all_empty_count": inh_all_empty_sec,
                "inhaled_all_empty_rate": round(inh_all_empty_sec / n_inh_all, 6),
            },
        },
        "positive_controls": pos_ctrl,
        "primary_estimand_pure_route_full_protocol": {
            "n_sys": n_sys,
            "k_sys": k_sys,
            "p_sys": round(p_sys, 6),
            "p_sys_exact": f"{k_sys}/{n_sys}",
            "n_inh_pure": n_inh_pure,
            "k_inh_pure": k_inh_pure,
            "p_inh_pure": round(p_inh_pure, 6),
            "p_inh_pure_exact": f"{k_inh_pure}/{n_inh_pure}",
            "risk_difference": round(p_sys - p_inh_pure, 6),
            "p_fisher_one_sided_exact": fisher_pure["p_one_sided_exact"],
            "p_fisher_one_sided": round(fisher_pure["p_one_sided_float"], 8),
            "p_fisher_two_sided_exact": fisher_pure["p_two_sided_exact"],
            "p_fisher_two_sided": round(fisher_pure["p_two_sided_float"], 8),
            "chan_zhang_unconditional_ci95": [round(cz_pure[0], 6), round(cz_pure[1], 6)],
            "newcombe_wilson_ci95": [round(nw_pure[0], 6), round(nw_pure[1], 6)],
            "criteria": {
                "exact_cohort_sizes_54_15_17": c1_sizes,
                "directional_p_sys_gt_p_inh_pure": c2_dir,
                "fisher_pure_p_lt_0_01": c3_fisher_pure,
                "unconditional_and_nw_ci95_gt_0": c4_cis_pure,
                "all_17_inhaled_fisher_p_lt_0_01": c5_all17,
            },
            "outcome": outcome,
            "ledger_standing": ledger_standing,
            "systemic_sweat_cl_any_ncts": sys_any_ncts,
            "systemic_sweat_cl_secondary_added_ncts": sys_sec_only_ncts,
            "inhaled_pure_sweat_cl_any_ncts": inh_pure_any_ncts,
        },
        "sensitivity_all_17_inhaled_including_dual_route": {
            "n_sys": n_sys,
            "k_sys": k_sys,
            "n_inh_all": n_inh_all,
            "k_inh_all": k_inh_all,
            "p_inh_all": round(p_inh_all, 6),
            "p_inh_all_exact": f"{k_inh_all}/{n_inh_all}",
            "risk_difference": round(p_sys - p_inh_all, 6),
            "p_fisher_one_sided_exact": fisher_all["p_one_sided_exact"],
            "p_fisher_one_sided": round(fisher_all["p_one_sided_float"], 8),
            "p_fisher_two_sided_exact": fisher_all["p_two_sided_exact"],
            "chan_zhang_unconditional_ci95": [round(cz_all[0], 6), round(cz_all[1], 6)],
            "newcombe_wilson_ci95": [round(nw_all[0], 6), round(nw_all[1], 6)],
            "inhaled_all_sweat_cl_any_ncts": inh_all_any_ncts,
            "dual_route_inhaled_ncts": sorted(dual_ncts),
        },
        "descriptive_context_non_gated": {
            "systemic_primary_only_count": len(pe15["systemic_sweat_cl_ncts"]),
            "systemic_secondary_only_added_count": len(sys_sec_only_ncts),
            "inhaled_all_secondary_only_added_count": len(inh_sec_only_ncts),
            "inhaled_all_secondary_only_added_ncts": inh_sec_only_ncts,
            "measure_plus_description_systemic_count": len(sys_desc_any_ncts),
            "measure_plus_description_inhaled_all_count": len(inh_desc_any_ncts),
            "measure_plus_description_inhaled_all_ncts": inh_desc_any_ncts,
        },
    }


def run_falsification_suite(
    done_path: Path = DEFAULT_DONE_MD,
    pred0015_results_path: Path = DEFAULT_PRED0015_RESULTS_JSON,
    snapshot_path: Path = DEFAULT_SNAPSHOT_JSON,
    outcomes_path: Path = DEFAULT_OUTCOMES_JSON,
) -> list[str]:
    """Run 6 deterministic negative controls verifying fail-closed behavior."""
    passed_controls: list[str] = []
    base_outcomes = json.loads(outcomes_path.read_text(encoding="utf-8"))

    # 1. Diagnostic eligibility sweat-chloride wording must NOT match
    if matches_sweat_chloride_measure("Sweat chloride >= 60 mmol/L diagnostic eligibility criteria"):
        raise FullProtocolCompartmentalizationError("Falsify 1 failed: diagnostic wording matched")
    passed_controls.append("diagnostic_eligibility_sweat_chloride_rejected")

    # 2. Sweat sodium / sweat rate without chloride/electrolyte must NOT match
    if matches_sweat_chloride_measure("Mean change in sweat sodium and sweat rate"):
        raise FullProtocolCompartmentalizationError("Falsify 2 failed: sweat sodium/rate matched")
    passed_controls.append("sweat_sodium_and_rate_rejected")

    # 3. Malformed non-list secondary_outcomes must fail closed
    try:
        normalize_outcomes_array("not-a-list", "NCT05248230", "secondary_outcomes")
        raise FullProtocolCompartmentalizationError("Falsify 3 failed: string secondary_outcomes accepted")
    except FullProtocolCompartmentalizationError as exc:
        if "Falsify 3 failed" in str(exc):
            raise
    passed_controls.append("malformed_secondary_outcomes_fails_closed")

    # 4 & 5. Dropping 1 study (70/71) or omitting verbatim NLM attribution must fail closed
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_root = Path(tmp_dir)
        tampered_70 = copy.deepcopy(base_outcomes)
        first_key = sorted(tampered_70["studies"].keys())[0]
        del tampered_70["studies"][first_key]
        tmp_path = tmp_root / "falsify_70.json"
        tmp_path.write_text(json.dumps(tampered_70), encoding="utf-8")
        try:
            evaluate_benchmark(
                done_path=done_path,
                pred0015_results_path=pred0015_results_path,
                snapshot_path=snapshot_path,
                outcomes_path=tmp_path,
                enforce_expected_hashes=False,
            )
            raise FullProtocolCompartmentalizationError("Falsify 4 failed: 70/71 fixture accepted")
        except FullProtocolCompartmentalizationError as exc:
            if "Falsify 4 failed" in str(exc):
                raise
        passed_controls.append("truncated_70_of_71_fixture_rejected")

        tampered_nlm = copy.deepcopy(base_outcomes)
        tampered_nlm["nlm_attribution"] = "missing"
        tmp_nlm = tmp_root / "falsify_nlm.json"
        tmp_nlm.write_text(json.dumps(tampered_nlm), encoding="utf-8")
        try:
            evaluate_benchmark(
                done_path=done_path,
                pred0015_results_path=pred0015_results_path,
                snapshot_path=snapshot_path,
                outcomes_path=tmp_nlm,
                enforce_expected_hashes=False,
            )
            raise FullProtocolCompartmentalizationError("Falsify 5 failed: missing NLM attribution accepted")
        except FullProtocolCompartmentalizationError as exc:
            if "Falsify 5 failed" in str(exc):
                raise
        passed_controls.append("missing_nlm_attribution_rejected")

    # 6. Null equal-proportion fixture (10/54 vs 3/15) must score as refuted
    f_null = fisher_hypergeometric_exact(54, 10, 15, 3)
    cz_null = chan_zhang_unconditional_rd_ci95(54, 10, 15, 3)
    if f_null["p_one_sided_fraction"] < Fraction(1, 100) or cz_null[0] > 0.0:
        raise FullProtocolCompartmentalizationError("Falsify 6 failed: null equal-rate fixture passed")
    passed_controls.append("null_equal_proportion_refuted")

    # 7. Homoglyph and zero-width Cf-split </untrusted_data> closing-tag lookalikes must be neutralized
    for adv in (
        "\u27E8/untrusted_data\u27E9",
        "\uFF1C/untrusted_data\uFF1E",
        "</un\u200Btrusted_da\u200Bta>",
        "\u202E</untrusted_data>",
    ):
        cleaned = sanitize_untrusted_text(adv)
        if "&lt;" not in cleaned or "&#47;" not in cleaned or "&gt;" not in cleaned:
            raise FullProtocolCompartmentalizationError(f"Falsify 7 failed: lookalike unescaped: {cleaned!r}")
    passed_controls.append("homoglyph_and_cf_split_untrusted_close_neutralized")

    return passed_controls


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--done", type=Path, default=DEFAULT_DONE_MD)
    ap.add_argument("--pred0015-results", type=Path, default=DEFAULT_PRED0015_RESULTS_JSON)
    ap.add_argument("--snapshot", type=Path, default=DEFAULT_SNAPSHOT_JSON)
    ap.add_argument("--outcomes", type=Path, default=DEFAULT_OUTCOMES_JSON)
    ap.add_argument("--results", type=Path, default=DEFAULT_RESULTS_JSON)
    ap.add_argument("--write-results", action="store_true")
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--falsify", action="store_true")
    ap.add_argument("--json", action="store_true")
    args = ap.parse_args(argv)

    computed = evaluate_benchmark(
        done_path=args.done,
        pred0015_results_path=args.pred0015_results,
        snapshot_path=args.snapshot,
        outcomes_path=args.outcomes,
    )
    if args.write_results:
        args.results.write_text(json.dumps(computed, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    if args.check:
        if not args.results.is_file():
            print(f"FAIL: missing {args.results}", file=sys.stderr)
            return 1
        on_disk = json.loads(args.results.read_text(encoding="utf-8"))
        if on_disk != computed:
            print("FAIL: results_2026-11-24.json out of sync with deterministic recomputation", file=sys.stderr)
            return 1

    falsify_list: list[str] = []
    if args.falsify:
        falsify_list = run_falsification_suite(
            done_path=args.done,
            pred0015_results_path=args.pred0015_results,
            snapshot_path=args.snapshot,
            outcomes_path=args.outcomes,
        )

    pe = computed["primary_estimand_pure_route_full_protocol"]
    sa = computed["sensitivity_all_17_inhaled_including_dual_route"]
    if args.json:
        print(json.dumps({"results": computed, "falsify_controls": falsify_list}, indent=2))
    else:
        msg = (
            f"PASS PRED-0016 ({pe['outcome']}): sys={pe['k_sys']}/{pe['n_sys']} ({pe['p_sys']:.4f}) "
            f"vs inh_pure={pe['k_inh_pure']}/{pe['n_inh_pure']} ({pe['p_inh_pure']:.4f}) "
            f"[inh_all={sa['k_inh_all']}/{sa['n_inh_all']}], "
            f"p_fisher_pure_1s={pe['p_fisher_one_sided_exact']} ({pe['p_fisher_one_sided']:.6f}), "
            f"cz_ci95={pe['chan_zhang_unconditional_ci95']}, nw_ci95={pe['newcombe_wilson_ci95']}"
        )
        if falsify_list:
            msg += f"; falsify={len(falsify_list)}/7 PASS"
        print(msg)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
