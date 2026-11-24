#!/usr/bin/env python3
"""Independent regex-free, math-free, zero-shared-import audit of `PRED-0016`.

Enforces Rule 2 of `discovery-practice` and `done.md` Section 4:
- Imports ZERO `re`, ZERO `math`, ZERO `tools/`, and ZERO benchmark helper modules.
- Re-checks `done_sha256`, `pred0015_results_sha256`, `outcomes_sha256`, and `54/17/15` NCT set equality.
- Classifies `primary_outcomes[].measure` and `secondary_outcomes[].measure` via character-by-character
  alphanumeric + bracket token scanning.
- Computes hypergeometric Fisher p-values via multiplicative `fractions.Fraction` products (`tolerance = 0`).
- Recomputes BOTH Newcombe-Wilson Method-10 AND Chan-Zhang / Miettinen-Nurminen unconditional profile-score
  95% confidence bounds in 50-digit `decimal.Decimal` arithmetic (`tolerance <= 1e-6`).
"""
from __future__ import annotations

import argparse
from decimal import Decimal, getcontext
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys
from typing import Any

getcontext().prec = 50

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
DEFAULT_DONE = HERE / "done.md"
DEFAULT_PRED0015_RESULTS = REPO_ROOT / "benchmarks" / "delivery_route_compartmentalization" / "results_2026-11-24.json"
DEFAULT_OUTCOMES = HERE / "early_phase_outcomes_2026-11-24.json"
DEFAULT_RESULTS = HERE / "results_2026-11-24.json"

REQUIRED_OUTCOME_KEYS = frozenset({"measure", "description", "time_frame"})
DIAG_WORDS = frozenset({"diagnosis", "diagnostic", "inclusion", "eligibility"})


def _sha256_lf(path: Path) -> str:
    txt = path.read_text(encoding="utf-8").replace("\r\n", "\n")
    return hashlib.sha256(txt.encode("utf-8")).hexdigest()


def _clause_tokens(text: str) -> tuple[str, ...]:
    """Character-level tokenizer without `re`, treating `isalnum() or '_'` as word chars
    and inserting `""` clause breaks at punctuation outside whitespace, `-`, `[`, `]`, `(`, `)`.
    """
    norm_chars: list[str] = []
    for ch in (text or "").lower():
        if ch in ("\u2212", "\u2013", "\u2014"):
            norm_chars.append("-")
        else:
            norm_chars.append(ch)
    norm = "".join(norm_chars)

    buf: list[str] = []
    toks: list[str] = []
    for ch in norm:
        if ch.isalnum() or ch == "_":
            buf.append(ch)
        else:
            if buf:
                toks.append("".join(buf))
                buf.clear()
            if not (ch.isspace() or ch in "-[]()"):
                if toks and toks[-1] != "":
                    toks.append("")
    if buf:
        toks.append("".join(buf))
    return tuple(toks)


def _has_seq(toks: tuple[str, ...], seq: tuple[str, ...]) -> bool:
    m = len(seq)
    if m == 0 or len(toks) < m:
        return False
    for idx in range(len(toks) - m + 1):
        if toks[idx : idx + m] == seq:
            return True
    return False


def _matches_sweat_chloride_independent(measure: str) -> bool:
    toks = _clause_tokens(measure)
    if not toks:
        return False
    is_sweat_cl = (
        _has_seq(toks, ("sweat", "chloride"))
        or _has_seq(toks, ("sweat", "cl"))
        or _has_seq(toks, ("chloride", "concentration", "in", "sweat"))
        or _has_seq(toks, ("chloride", "in", "sweat"))
        or _has_seq(toks, ("sweat", "electrolyte"))
        or _has_seq(toks, ("sweat", "electrolytes"))
    )
    if not is_sweat_cl:
        return False
    is_diag = any(t in DIAG_WORDS for t in toks) or _has_seq(toks, ("screening", "criteria"))
    return not is_diag


def _binom_fraction(n: int, k: int) -> Fraction:
    if k < 0 or k > n:
        return Fraction(0, 1)
    m = min(k, n - k)
    acc = Fraction(1, 1)
    for r in range(1, m + 1):
        acc *= Fraction(n - m + r, r)
    return acc


def _fisher_tails_multiplicative(n1: int, k1: int, n2: int, k2: int) -> tuple[Fraction, Fraction]:
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


def _newcombe_wilson_decimal(n1: int, k1: int, n2: int, k2: int) -> tuple[Decimal, Decimal]:
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


def _bisect_mn_restricted_p2_decimal(x1: int, n1: int, x2: int, n2: int, delta: Decimal) -> Decimal:
    """50-digit Decimal monotone polynomial bisection for Miettinen-Nurminen restricted MLE."""
    zero = Decimal(0)
    one = Decimal(1)
    half = Decimal("0.5")
    eps = Decimal("1e-24")
    lo = max(zero, -delta)
    hi = min(one, one - delta)
    if hi <= lo + eps:
        return half * (lo + hi)

    dx1 = Decimal(x1)
    dn1_x1 = Decimal(n1 - x1)
    dx2 = Decimal(x2)
    dn2_x2 = Decimal(n2 - x2)

    def poly_score(p2: Decimal) -> Decimal:
        p1 = p2 + delta
        q1 = one - p1
        q2 = one - p2
        return dx1 * q1 * p2 * q2 - dn1_x1 * p1 * p2 * q2 + dx2 * p1 * q1 * q2 - dn2_x2 * p1 * q1 * p2

    if poly_score(lo + eps) <= zero:
        return lo
    if poly_score(hi - eps) >= zero:
        return hi
    a, b = lo, hi
    for _ in range(90):
        mid = half * (a + b)
        if poly_score(mid) > zero:
            a = mid
        else:
            b = mid
    return half * (a + b)


def _mn_z_decimal(x1: int, n1: int, x2: int, n2: int, delta: Decimal) -> Decimal:
    zero = Decimal(0)
    one = Decimal(1)
    num = (Decimal(x1) / Decimal(n1)) - (Decimal(x2) / Decimal(n2)) - delta
    if abs(num) < Decimal("1e-24"):
        return zero
    p2_t = _bisect_mn_restricted_p2_decimal(x1, n1, x2, n2, delta)
    p1_t = min(one, max(zero, p2_t + delta))
    tot = Decimal(n1 + n2)
    var = (p1_t * (one - p1_t) / Decimal(n1) + p2_t * (one - p2_t) / Decimal(n2)) * (tot / (tot - one))
    if var <= Decimal("1e-36"):
        return Decimal("1e18") if num > zero else Decimal("-1e18")
    return num / var.sqrt()


def _mn_z_float_screen(x1: int, n1: int, x2: int, n2: int, delta: float) -> float:
    num = (x1 / n1) - (x2 / n2) - delta
    if abs(num) < 1e-15:
        return 0.0
    lo = max(0.0, -delta)
    hi = min(1.0, 1.0 - delta)
    if hi <= lo + 1e-15:
        p2_t = 0.5 * (lo + hi)
    else:
        def pf(p2: float) -> float:
            p1 = p2 + delta
            return (
                x1 * (1.0 - p1) * p2 * (1.0 - p2)
                - (n1 - x1) * p1 * p2 * (1.0 - p2)
                + x2 * p1 * (1.0 - p1) * (1.0 - p2)
                - (n2 - x2) * p1 * (1.0 - p1) * p2
            )
        if pf(lo + 1e-15) <= 0.0:
            p2_t = lo
        elif pf(hi - 1e-15) >= 0.0:
            p2_t = hi
        else:
            a, b = lo, hi
            for _ in range(52):
                mid = 0.5 * (a + b)
                if pf(mid) > 0.0:
                    a = mid
                else:
                    b = mid
            p2_t = 0.5 * (a + b)
    p1_t = min(1.0, max(0.0, p2_t + delta))
    tot = n1 + n2
    var = (p1_t * (1.0 - p1_t) / n1 + p2_t * (1.0 - p2_t) / n2) * (tot / (tot - 1.0))
    if var <= 1e-24:
        return 1e12 if num > 0.0 else -1e12
    return num / (var ** 0.5)


def _chan_zhang_ci95_decimal(n1: int, k1: int, n2: int, k2: int) -> tuple[Decimal, Decimal]:
    """50-digit Decimal Chan-Zhang unconditional 95% CI."""
    zero = Decimal(0)
    one = Decimal(1)
    half = Decimal("0.5")
    third = one / Decimal(3)
    alpha_tail = Decimal("0.025")
    tie_tol = Decimal("1e-11")

    c1 = [Decimal(int(_binom_fraction(n1, i))) for i in range(n1 + 1)]
    c2 = [Decimal(int(_binom_fraction(n2, j))) for j in range(n2 + 1)]
    c1_f = [float(v) for v in c1]
    c2_f = [float(v) for v in c2]

    def tail_sup(delta: Decimal, upper: bool, full_decimal_prob: bool = False) -> Decimal:
        lo_p2 = max(zero, -delta)
        hi_p2 = min(one, one - delta)
        z_obs = _mn_z_decimal(k1, n1, k2, n2, delta)
        df = float(delta)
        zf_obs = float(z_obs)
        pairs: list[tuple[int, int]] = []
        for x1 in range(n1 + 1):
            for x2 in range(n2 + 1):
                zf_xy = _mn_z_float_screen(x1, n1, x2, n2, df)
                if abs(zf_xy - zf_obs) <= 1e-5 or full_decimal_prob:
                    z_xy = _mn_z_decimal(x1, n1, x2, n2, delta)
                    if (upper and z_xy >= z_obs - tie_tol) or ((not upper) and z_xy <= z_obs + tie_tol):
                        pairs.append((x1, x2))
                else:
                    if (upper and zf_xy >= zf_obs) or ((not upper) and zf_xy <= zf_obs):
                        pairs.append((x1, x2))
        if not pairs:
            return zero
        if len(pairs) == (n1 + 1) * (n2 + 1):
            return one

        def prob_dec(p2: Decimal) -> Decimal:
            p1 = min(one, max(zero, p2 + delta))
            if p1 <= zero:
                b1 = [one] + [zero] * n1
            elif p1 >= one:
                b1 = [zero] * n1 + [one]
            else:
                q1 = one - p1
                b1 = [c1[i] * (p1 ** i) * (q1 ** (n1 - i)) for i in range(n1 + 1)]
            if p2 <= zero:
                b2 = [one] + [zero] * n2
            elif p2 >= one:
                b2 = [zero] * n2 + [one]
            else:
                q2 = one - p2
                b2 = [c2[j] * (p2 ** j) * (q2 ** (n2 - j)) for j in range(n2 + 1)]
            acc = zero
            for x1, x2 in pairs:
                acc += b1[x1] * b2[x2]
            return acc

        def prob_flt(p2: float) -> float:
            p1 = min(1.0, max(0.0, p2 + df))
            b1 = [1.0] + [0.0] * n1 if p1 <= 0.0 else ([0.0] * n1 + [1.0] if p1 >= 1.0 else [c1_f[i] * (p1 ** i) * ((1.0 - p1) ** (n1 - i)) for i in range(n1 + 1)])
            b2 = [1.0] + [0.0] * n2 if p2 <= 0.0 else ([0.0] * n2 + [1.0] if p2 >= 1.0 else [c2_f[j] * (p2 ** j) * ((1.0 - p2) ** (n2 - j)) for j in range(n2 + 1)])
            return sum(b1[x1] * b2[x2] for x1, x2 in pairs)

        if hi_p2 <= lo_p2 + Decimal("1e-20"):
            return prob_dec(half * (lo_p2 + hi_p2))

        lo_f = float(lo_p2)
        hi_f = float(hi_p2)
        grid_m = 200
        best_idx = 0
        best_vf = -1.0
        for idx in range(grid_m + 1):
            cand_f = lo_f + (hi_f - lo_f) * (idx / grid_m)
            vf = prob_flt(cand_f)
            if vf > best_vf:
                best_vf = vf
                best_idx = idx
        span = hi_p2 - lo_p2
        step = span / Decimal(grid_m)
        best_p2 = lo_p2 + span * Decimal(best_idx) / Decimal(grid_m)
        best_v = prob_dec(best_p2)
        left = max(lo_p2, best_p2 - step)
        right = min(hi_p2, best_p2 + step)
        for _ in range(32):
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
    if k1 == 0 and k2 == n2:
        low_ci = Decimal(-1)
    else:
        lo = Decimal("-1") + Decimal("1e-12")
        hi = rd_obs
        if tail_sup(lo, True) > alpha_tail:
            low_ci = Decimal(-1)
        else:
            for _ in range(52):
                mid = half * (lo + hi)
                if tail_sup(mid, True) > alpha_tail:
                    hi = mid
                else:
                    lo = mid
            low_ci = hi

    if k1 == n1 and k2 == 0:
        high_ci = one
    else:
        lo = rd_obs
        hi = one - Decimal("1e-12")
        if tail_sup(hi, False) > alpha_tail:
            high_ci = one
        else:
            for _ in range(52):
                mid = half * (lo + hi)
                if tail_sup(mid, False) > alpha_tail:
                    lo = mid
                else:
                    hi = mid
            high_ci = lo
    return low_ci, high_ci


def run_independent_audit(
    done_path: Path = DEFAULT_DONE,
    pred0015_results_path: Path = DEFAULT_PRED0015_RESULTS,
    outcomes_path: Path = DEFAULT_OUTCOMES,
    results_path: Path = DEFAULT_RESULTS,
) -> dict[str, Any]:
    done_sha = _sha256_lf(done_path)
    p15_sha = _sha256_lf(pred0015_results_path)
    outcomes_sha = _sha256_lf(outcomes_path)

    p15_doc = json.loads(pred0015_results_path.read_text(encoding="utf-8"))
    pe15 = p15_doc["primary_estimand_early_phase"]
    sys_set = set(pe15["systemic_ncts"])
    inh_set = set(pe15["inhaled_ncts"])
    dual_set = set(pe15["dual_route_inhaled_ncts"])
    inh_pure_set = inh_set - dual_set

    outcomes_doc = json.loads(outcomes_path.read_text(encoding="utf-8"))
    studies_map = outcomes_doc.get("studies") or {}

    mismatches: list[str] = []
    if len(sys_set) != 54 or len(inh_set) != 17 or len(inh_pure_set) != 15:
        mismatches.append(f"PRED-0015 set sizes != 54/17/15: {len(sys_set)}, {len(inh_set)}, {len(inh_pure_set)}")
    if set(studies_map.keys()) != (sys_set | inh_set):
        mismatches.append("outcomes_2026-11-24.json NCT keys != locked 71 NCT union")

    k_sys = k_inh_pure = k_inh_all = 0
    sys_empty_sec = inh_pure_empty_sec = inh_all_empty_sec = 0
    sys_sec_only_count = inh_sec_only_count = 0
    sys_desc_count = inh_desc_count = 0
    sys_matched: list[str] = []
    inh_pure_matched: list[str] = []
    inh_all_matched: list[str] = []

    for nct_id in sorted(sys_set | inh_set):
        rec = studies_map.get(nct_id) or {}
        prim_arr = rec.get("primary_outcomes")
        sec_arr = rec.get("secondary_outcomes")
        if prim_arr is None:
            prim_arr = []
        if sec_arr is None:
            sec_arr = []
        if not isinstance(prim_arr, list) or not isinstance(sec_arr, list):
            mismatches.append(f"Non-list outcomes on {nct_id}")
            continue
        prim_hit = False
        sec_hit = False
        desc_hit = False
        for o in prim_arr:
            if not isinstance(o, dict) or not REQUIRED_OUTCOME_KEYS.issubset(o.keys()):
                mismatches.append(f"Malformed primary outcome dict on {nct_id}")
                continue
            m_str = str(o.get("measure") or "")
            d_str = str(o.get("description") or "")
            if _matches_sweat_chloride_independent(m_str):
                prim_hit = True
                desc_hit = True
            elif _matches_sweat_chloride_independent(f"{m_str} {d_str}"):
                desc_hit = True
        for o in sec_arr:
            if not isinstance(o, dict) or not REQUIRED_OUTCOME_KEYS.issubset(o.keys()):
                mismatches.append(f"Malformed secondary outcome dict on {nct_id}")
                continue
            m_str = str(o.get("measure") or "")
            d_str = str(o.get("description") or "")
            if _matches_sweat_chloride_independent(m_str):
                sec_hit = True
                desc_hit = True
            elif _matches_sweat_chloride_independent(f"{m_str} {d_str}"):
                desc_hit = True

        hit = prim_hit or sec_hit
        if nct_id in sys_set:
            if len(sec_arr) == 0:
                sys_empty_sec += 1
            if hit:
                k_sys += 1
                sys_matched.append(nct_id)
            if sec_hit and not prim_hit:
                sys_sec_only_count += 1
            if desc_hit:
                sys_desc_count += 1
        else:
            if len(sec_arr) == 0:
                inh_all_empty_sec += 1
                if nct_id in inh_pure_set:
                    inh_pure_empty_sec += 1
            if hit:
                k_inh_all += 1
                inh_all_matched.append(nct_id)
                if nct_id in inh_pure_set:
                    k_inh_pure += 1
                    inh_pure_matched.append(nct_id)
            if sec_hit and not prim_hit:
                inh_sec_only_count += 1
            if desc_hit:
                inh_desc_count += 1

    n_sys = len(sys_set)
    n_inh_pure = len(inh_pure_set)
    n_inh_all = len(inh_set)

    p_pure_1s, p_pure_2s = _fisher_tails_multiplicative(n_sys, k_sys, n_inh_pure, k_inh_pure)
    p_all_1s, p_all_2s = _fisher_tails_multiplicative(n_sys, k_sys, n_inh_all, k_inh_all)
    nw_pure_lo, nw_pure_hi = _newcombe_wilson_decimal(n_sys, k_sys, n_inh_pure, k_inh_pure)
    cz_pure_lo, cz_pure_hi = _chan_zhang_ci95_decimal(n_sys, k_sys, n_inh_pure, k_inh_pure)
    nw_all_lo, nw_all_hi = _newcombe_wilson_decimal(n_sys, k_sys, n_inh_all, k_inh_all)
    cz_all_lo, cz_all_hi = _chan_zhang_ci95_decimal(n_sys, k_sys, n_inh_all, k_inh_all)

    res_doc = json.loads(results_path.read_text(encoding="utf-8"))
    prim_res = res_doc.get("primary_estimand_pure_route_full_protocol") or {}
    sens_res = res_doc.get("sensitivity_all_17_inhaled_including_dual_route") or {}
    desc_res = res_doc.get("descriptive_context_non_gated") or {}
    empty_sec_res = ((res_doc.get("dataset") or {}).get("empty_secondary_outcomes")) or {}

    if res_doc.get("done_sha256") != done_sha:
        mismatches.append(f"done_sha256 mismatch: {res_doc.get('done_sha256')} != {done_sha}")
    if res_doc.get("pred0015_results_sha256") != p15_sha:
        mismatches.append("pred0015_results_sha256 mismatch")
    if res_doc.get("outcomes_sha256") != outcomes_sha:
        mismatches.append("outcomes_sha256 mismatch")

    if (n_sys, k_sys, n_inh_pure, k_inh_pure) != (
        prim_res.get("n_sys"),
        prim_res.get("k_sys"),
        prim_res.get("n_inh_pure"),
        prim_res.get("k_inh_pure"),
    ):
        mismatches.append(
            f"Pure-route 2x2 mismatch: {(n_sys, k_sys, n_inh_pure, k_inh_pure)} vs "
            f"{(prim_res.get('n_sys'), prim_res.get('k_sys'), prim_res.get('n_inh_pure'), prim_res.get('k_inh_pure'))}"
        )
    if (n_inh_all, k_inh_all) != (sens_res.get("n_inh_all"), sens_res.get("k_inh_all")):
        mismatches.append(f"All-17 inhaled mismatch: {(n_inh_all, k_inh_all)} vs {(sens_res.get('n_inh_all'), sens_res.get('k_inh_all'))}")

    if sys_matched != prim_res.get("systemic_sweat_cl_any_ncts"):
        mismatches.append("systemic_sweat_cl_any_ncts list mismatch")
    if inh_pure_matched != prim_res.get("inhaled_pure_sweat_cl_any_ncts"):
        mismatches.append("inhaled_pure_sweat_cl_any_ncts list mismatch")
    if inh_all_matched != sens_res.get("inhaled_all_sweat_cl_any_ncts"):
        mismatches.append("inhaled_all_sweat_cl_any_ncts list mismatch")

    if p_pure_1s != Fraction(str(prim_res.get("p_fisher_one_sided_exact", "0/1"))):
        mismatches.append(f"Pure Fisher 1s mismatch: {p_pure_1s} vs {prim_res.get('p_fisher_one_sided_exact')}")
    if p_pure_2s != Fraction(str(prim_res.get("p_fisher_two_sided_exact", "0/1"))):
        mismatches.append(f"Pure Fisher 2s mismatch: {p_pure_2s} vs {prim_res.get('p_fisher_two_sided_exact')}")
    if p_all_1s != Fraction(str(sens_res.get("p_fisher_one_sided_exact", "0/1"))):
        mismatches.append(f"All-17 Fisher 1s mismatch: {p_all_1s} vs {sens_res.get('p_fisher_one_sided_exact')}")
    if p_all_2s != Fraction(str(sens_res.get("p_fisher_two_sided_exact", "0/1"))):
        mismatches.append(f"All-17 Fisher 2s mismatch: {p_all_2s} vs {sens_res.get('p_fisher_two_sided_exact')}")

    tol = Decimal("0.000001")
    rd_pure_dec = Decimal(k_sys) / Decimal(n_sys) - Decimal(k_inh_pure) / Decimal(n_inh_pure)
    rd_all_dec = Decimal(k_sys) / Decimal(n_sys) - Decimal(k_inh_all) / Decimal(n_inh_all)
    if prim_res.get("risk_difference") is None or abs(rd_pure_dec - Decimal(str(prim_res.get("risk_difference")))) > tol:
        mismatches.append(f"Pure risk_difference mismatch: {rd_pure_dec:.7f} vs {prim_res.get('risk_difference')}")
    if sens_res.get("risk_difference") is None or abs(rd_all_dec - Decimal(str(sens_res.get("risk_difference")))) > tol:
        mismatches.append(f"All-17 risk_difference mismatch: {rd_all_dec:.7f} vs {sens_res.get('risk_difference')}")

    for label, d_lo, d_hi, ref_pair in (
        ("pure NW", nw_pure_lo, nw_pure_hi, prim_res.get("newcombe_wilson_ci95")),
        ("pure CZ", cz_pure_lo, cz_pure_hi, prim_res.get("chan_zhang_unconditional_ci95")),
        ("all17 NW", nw_all_lo, nw_all_hi, sens_res.get("newcombe_wilson_ci95")),
        ("all17 CZ", cz_all_lo, cz_all_hi, sens_res.get("chan_zhang_unconditional_ci95")),
    ):
        if (
            not isinstance(ref_pair, list)
            or len(ref_pair) != 2
            or abs(d_lo - Decimal(str(ref_pair[0]))) > tol
            or abs(d_hi - Decimal(str(ref_pair[1]))) > tol
        ):
            mismatches.append(f"{label} CI95 mismatch (> 1e-6): ({d_lo:.7f}, {d_hi:.7f}) vs {ref_pair}")

    expected_criteria = {
        "exact_cohort_sizes_54_15_17": (n_sys == 54 and n_inh_pure == 15 and n_inh_all == 17),
        "directional_p_sys_gt_p_inh_pure": (Fraction(k_sys, n_sys) > Fraction(k_inh_pure, n_inh_pure)),
        "fisher_pure_p_lt_0_01": (p_pure_1s < Fraction(1, 100)),
        "unconditional_and_nw_ci95_gt_0": (cz_pure_lo > Decimal(0) and nw_pure_lo > Decimal(0)),
        "all_17_inhaled_fisher_p_lt_0_01": (
            Fraction(k_sys, n_sys) > Fraction(k_inh_all, n_inh_all) and p_all_1s < Fraction(1, 100)
        ),
    }
    expected_supported = all(expected_criteria.values())
    expected_outcome = "supported" if expected_supported else "refuted"
    expected_standing = "open" if expected_supported else "refuted"
    if prim_res.get("criteria") != expected_criteria:
        mismatches.append(f"criteria mismatch: {prim_res.get('criteria')} != {expected_criteria}")
    if prim_res.get("outcome") != expected_outcome:
        mismatches.append(f"outcome mismatch: {prim_res.get('outcome')} != {expected_outcome}")
    if prim_res.get("ledger_standing") != expected_standing:
        mismatches.append(f"ledger_standing mismatch: {prim_res.get('ledger_standing')} != {expected_standing}")

    if (
        empty_sec_res.get("systemic_empty_count"),
        empty_sec_res.get("inhaled_pure_empty_count"),
        empty_sec_res.get("inhaled_all_empty_count"),
    ) != (sys_empty_sec, inh_pure_empty_sec, inh_all_empty_sec):
        mismatches.append("empty secondary_outcomes provenance count mismatch")

    if (
        desc_res.get("systemic_secondary_only_added_count"),
        desc_res.get("inhaled_all_secondary_only_added_count"),
        desc_res.get("measure_plus_description_systemic_count"),
        desc_res.get("measure_plus_description_inhaled_all_count"),
    ) != (sys_sec_only_count, inh_sec_only_count, sys_desc_count, inh_desc_count):
        mismatches.append(
            "descriptive_context_non_gated count mismatch: "
            f"{(sys_sec_only_count, inh_sec_only_count, sys_desc_count, inh_desc_count)} vs "
            f"{(desc_res.get('systemic_secondary_only_added_count'), desc_res.get('inhaled_all_secondary_only_added_count'), desc_res.get('measure_plus_description_systemic_count'), desc_res.get('measure_plus_description_inhaled_all_count'))}"
        )

    return {
        "passed": len(mismatches) == 0,
        "mismatches": mismatches,
        "done_sha256": done_sha,
        "pred0015_results_sha256": p15_sha,
        "outcomes_sha256": outcomes_sha,
        "n_sys": n_sys,
        "k_sys": k_sys,
        "n_inh_pure": n_inh_pure,
        "k_inh_pure": k_inh_pure,
        "n_inh_all": n_inh_all,
        "k_inh_all": k_inh_all,
        "outcome": expected_outcome,
        "p_fisher_pure_one_sided_exact": f"{p_pure_1s.numerator}/{p_pure_1s.denominator}",
        "p_fisher_pure_one_sided": round(float(p_pure_1s), 8),
        "p_fisher_all_one_sided_exact": f"{p_all_1s.numerator}/{p_all_1s.denominator}",
        "chan_zhang_pure_ci95": [round(float(cz_pure_lo), 6), round(float(cz_pure_hi), 6)],
        "newcombe_wilson_pure_ci95": [round(float(nw_pure_lo), 6), round(float(nw_pure_hi), 6)],
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--done", type=Path, default=DEFAULT_DONE)
    ap.add_argument("--pred0015-results", type=Path, default=DEFAULT_PRED0015_RESULTS)
    ap.add_argument("--outcomes", type=Path, default=DEFAULT_OUTCOMES)
    ap.add_argument("--results", type=Path, default=DEFAULT_RESULTS)
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--json", action="store_true")
    args = ap.parse_args(argv)

    rep = run_independent_audit(args.done, args.pred0015_results, args.outcomes, args.results)
    if args.json:
        print(json.dumps(rep, indent=2))
    if not rep["passed"]:
        for m in rep["mismatches"]:
            print(f"FAIL independent_audit: {m}", file=sys.stderr)
        return 1
    if not args.json:
        print(
            f"PASS independent_audit: sys={rep['k_sys']}/{rep['n_sys']} "
            f"inh_pure={rep['k_inh_pure']}/{rep['n_inh_pure']} inh_all={rep['k_inh_all']}/{rep['n_inh_all']} "
            f"p_pure_1s={rep['p_fisher_pure_one_sided_exact']} "
            f"cz_ci95={rep['chan_zhang_pure_ci95']} nw_ci95={rep['newcombe_wilson_pure_ci95']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
