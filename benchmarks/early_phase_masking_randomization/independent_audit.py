#!/usr/bin/env python3
"""Independent zero-`re`, zero-`math`, zero-`float`-prescreen audit of `PRED-0017`.

Enforces Rule 2 of `discovery-practice` and `done.md` Section 4:
- Imports ZERO `re`, ZERO `math`, ZERO `tools/`, and ZERO benchmark helper modules.
- Uses ZERO IEEE `float` pre-screening inside its 50-digit `decimal.Decimal` Chan-Zhang verifier.
- Re-checks `done_sha256`, `pred0015_results_sha256`, `fixture_sha256`, and `54/17/15` NCT set equality.
- Evaluates raw `designInfo.allocation` and `designInfo.maskingInfo.masking` enums directly from JSON.
- Computes hypergeometric Fisher p-values via multiplicative `fractions.Fraction` products (`tolerance = 0`).
- Recomputes BOTH Newcombe-Wilson Method-10 AND conservative sawtooth-safe multi-peak Chan-Zhang 95% CIs
  entirely in `decimal.Decimal` arithmetic (`tolerance <= 1e-6`).
"""
from __future__ import annotations

import argparse
import ast
from decimal import Decimal, getcontext, localcontext
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
DEFAULT_FIXTURE = HERE / "early_phase_design_eligibility_2026-11-24.json"
DEFAULT_RESULTS = HERE / "results_2026-11-24.json"

VALID_ALLOC = frozenset({"RANDOMIZED", "NON_RANDOMIZED", "NA", None})
VALID_MASK = frozenset({"NONE", "SINGLE", "DOUBLE", "TRIPLE", "QUADRUPLE", None})
MASKED_SET = frozenset({"SINGLE", "DOUBLE", "TRIPLE", "QUADRUPLE"})


def _sha256_lf(path: Path) -> str:
    txt = path.read_text(encoding="utf-8").replace("\r\n", "\n")
    return hashlib.sha256(txt.encode("utf-8")).hexdigest()


def verify_ast_independence(self_path: Path = Path(__file__)) -> None:
    """AST guard proving zero imports of `re`, `math`, `tools`, or sibling benchmark modules."""
    tree = ast.parse(self_path.read_text(encoding="utf-8"), filename=str(self_path))
    banned = {"re", "math", "early_phase_masking_randomization", "delivery_route_compartmentalization", "scope", "snapshot"}
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                root_mod = alias.name.split(".")[0]
                if root_mod in banned or root_mod == "tools":
                    raise RuntimeError(f"Forbidden import in independent_audit.py: {alias.name}")
        elif isinstance(node, ast.ImportFrom) and node.module:
            root_mod = node.module.split(".")[0]
            if root_mod in banned or root_mod == "tools":
                raise RuntimeError(f"Forbidden from-import in independent_audit.py: {node.module}")


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


def _newton_bisect_mn_p2_decimal(x1: int, n1: int, x2: int, n2: int, delta: Decimal) -> Decimal:
    """Safeguarded Newton-bisection on exact constrained binomial log-likelihood derivative (zero cubic expansion)."""
    zero = Decimal(0)
    one = Decimal(1)
    half = Decimal("0.5")
    eps = Decimal("1e-15")

    lo = max(zero, -delta)
    hi = min(one, one - delta)
    if hi <= lo + Decimal("1e-14"):
        return half * (lo + hi)

    dx1 = Decimal(x1)
    dy1 = Decimal(n1 - x1)
    dx2 = Decimal(x2)
    dy2 = Decimal(n2 - x2)

    # Check exact analytical boundary limits on (lo, hi) before interior root-finding:
    # As p2 -> lo+: either p2 -> 0+ (delta >= 0) or p1 -> 0+ (delta < 0).
    # Only finite <= 0 when the vanishing denominator has zero numerator count:
    if delta >= zero and x2 == 0:
        # lim_{p2->0+} d ell / d p2 = x1/delta - (n1-x1)/(1-delta) - n2 (when delta > 0)
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

    a = lo + (hi - lo) * eps
    b = hi - (hi - lo) * eps

    def grad_ell(p2: Decimal) -> Decimal:
        p1 = p2 + delta
        t1 = (dx1 / p1 if x1 > 0 else zero) - (dy1 / (one - p1) if x1 < n1 else zero)
        t2 = (dx2 / p2 if x2 > 0 else zero) - (dy2 / (one - p2) if x2 < n2 else zero)
        return t1 + t2

    def dgrad_ell(p2: Decimal) -> Decimal:
        p1 = p2 + delta
        d1 = -(dx1 / (p1 * p1) if x1 > 0 else zero) - (dy1 / ((one - p1) * (one - p1)) if x1 < n1 else zero)
        d2 = -(dx2 / (p2 * p2) if x2 > 0 else zero) - (dy2 / ((one - p2) * (one - p2)) if x2 < n2 else zero)
        return d1 + d2

    if grad_ell(a) <= zero:
        return lo
    if grad_ell(b) >= zero:
        return hi

    x = half * (a + b)
    tol = Decimal("1e-18")
    for _ in range(36):
        fx = grad_ell(x)
        if abs(fx) <= tol or (b - a) <= tol:
            return x
        if fx > zero:
            a = x
        else:
            b = x
        dfx = dgrad_ell(x)
        if dfx < -eps:
            cand = x - fx / dfx
            if a < cand < b:
                x = cand
                continue
        x = half * (a + b)
    return half * (a + b)


def _mn_z_decimal(x1: int, n1: int, x2: int, n2: int, delta: Decimal) -> Decimal:
    zero = Decimal(0)
    one = Decimal(1)
    num = (Decimal(x1) / Decimal(n1)) - (Decimal(x2) / Decimal(n2)) - delta
    if abs(num) < Decimal("1e-24"):
        return zero
    p2_t = _newton_bisect_mn_p2_decimal(x1, n1, x2, n2, delta)
    p1_t = min(one, max(zero, p2_t + delta))
    tot = Decimal(n1 + n2)
    var = (p1_t * (one - p1_t) / Decimal(n1) + p2_t * (one - p2_t) / Decimal(n2)) * (tot / (tot - one))
    if var <= Decimal("1e-36"):
        return Decimal("1e18") if num > zero else Decimal("-1e18")
    return num / var.sqrt()


def _chan_zhang_ci95_decimal(n1: int, k1: int, n2: int, k2: int) -> tuple[Decimal, Decimal]:
    """100% pure `Decimal` (zero `float` pre-screen) sawtooth-safe multi-peak Chan-Zhang 95% CI."""
    zero = Decimal(0)
    one = Decimal(1)
    half = Decimal("0.5")
    third = one / Decimal(3)
    alpha_tail = Decimal("0.025")
    tie_tol = Decimal("1e-11")

    c1 = [Decimal(int(_binom_fraction(n1, i))) for i in range(n1 + 1)]
    c2 = [Decimal(int(_binom_fraction(n2, j))) for j in range(n2 + 1)]
    p1_hats = [Decimal(x1) / Decimal(n1) for x1 in range(n1 + 1)]
    p2_hats = [Decimal(x2) / Decimal(n2) for x2 in range(n2 + 1)]

    def tail_sup(delta: Decimal, upper: bool, prec_digits: int = 28) -> Decimal:
        with localcontext() as ctx:
            ctx.prec = prec_digits
            lo_p2 = max(zero, -delta)
            hi_p2 = min(one, one - delta)
            z_obs = _mn_z_decimal(k1, n1, k2, n2, delta)
            pairs: list[tuple[int, int]] = []
            for x1 in range(n1 + 1):
                p1_h = p1_hats[x1]
                for x2 in range(n2 + 1):
                    num_xy = p1_h - p2_hats[x2] - delta
                    if upper and z_obs > tie_tol and num_xy < zero:
                        continue
                    if (not upper) and z_obs < -tie_tol and num_xy > zero:
                        continue
                    z_xy = _mn_z_decimal(x1, n1, x2, n2, delta)
                    if (upper and z_xy >= z_obs - tie_tol) or ((not upper) and z_xy <= z_obs + tie_tol):
                        pairs.append((x1, x2))
            if not pairs:
                return zero
            if len(pairs) == (n1 + 1) * (n2 + 1):
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
                for x1, x2 in pairs:
                    acc += b1[x1] * b2[x2]
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

            for idx in peak_indices:
                center = lo_p2 + step * Decimal(idx)
                left = max(lo_p2, center - step)
                right = min(hi_p2, center + step)
                for _ in range(24):
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
            for it in range(42):
                mid = half * (lo_b + hi_b)
                prec_use = 45 if it >= 28 else 28
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
            for it in range(42):
                mid = half * (lo_b + hi_b)
                prec_use = 45 if it >= 28 else 28
                if tail_sup(mid, False, prec_digits=prec_use) > alpha_tail:
                    lo_b = mid
                else:
                    hi_b = mid
            high_ci = lo_b

    return low_ci, high_ci


def run_independent_audit(
    done_path: Path = DEFAULT_DONE,
    pred0015_results_path: Path = DEFAULT_PRED0015_RESULTS,
    fixture_path: Path = DEFAULT_FIXTURE,
    results_path: Path = DEFAULT_RESULTS,
) -> dict[str, Any]:
    verify_ast_independence()
    done_sha = _sha256_lf(done_path)
    p15_sha = _sha256_lf(pred0015_results_path)
    fixture_sha = _sha256_lf(fixture_path)

    p15_doc = json.loads(pred0015_results_path.read_text(encoding="utf-8"))
    pe15 = p15_doc["primary_estimand_early_phase"]
    sys_set = set(pe15["systemic_ncts"])
    inh_set = set(pe15["inhaled_ncts"])
    dual_set = set(pe15["dual_route_inhaled_ncts"])
    inh_pure_set = inh_set - dual_set

    fix_doc = json.loads(fixture_path.read_text(encoding="utf-8"))
    studies = fix_doc.get("studies") or {}

    mismatches: list[str] = []
    if set(studies.keys()) != (sys_set | inh_set):
        mismatches.append("Fixture NCT set != locked 71 NCTs")

    k_sys = 0
    k_inh_pure = 0
    k_inh_all = 0
    sys_matched: list[str] = []
    inh_pure_matched: list[str] = []
    inh_all_matched: list[str] = []

    for nid in sorted(studies):
        st = studies[nid]
        dmod = st.get("design_module") or {}
        dinfo = dmod.get("designInfo") or {}
        alloc = dinfo.get("allocation")
        minfo = dinfo.get("maskingInfo") or {}
        mask = minfo.get("masking")
        if alloc not in VALID_ALLOC:
            mismatches.append(f"{nid}: invalid allocation {alloc!r}")
            continue
        if mask not in VALID_MASK:
            mismatches.append(f"{nid}: invalid masking {mask!r}")
            continue
        hit = (alloc == "RANDOMIZED") and (mask in MASKED_SET)
        if nid in sys_set:
            if hit:
                k_sys += 1
                sys_matched.append(nid)
        else:
            if hit:
                k_inh_all += 1
                inh_all_matched.append(nid)
                if nid in inh_pure_set:
                    k_inh_pure += 1
                    inh_pure_matched.append(nid)

    n_sys = len(sys_set)
    n_inh_pure = len(inh_pure_set)
    n_inh_all = len(inh_set)

    p_pure_1s, p_pure_2s = _fisher_tails_multiplicative(n_sys, k_sys, n_inh_pure, k_inh_pure)
    p_all_1s, _ = _fisher_tails_multiplicative(n_sys, k_sys, n_inh_all, k_inh_all)
    nw_pure_lo, nw_pure_hi = _newcombe_wilson_decimal(n_sys, k_sys, n_inh_pure, k_inh_pure)
    cz_pure_lo, cz_pure_hi = _chan_zhang_ci95_decimal(n_sys, k_sys, n_inh_pure, k_inh_pure)
    nw_all_lo, nw_all_hi = _newcombe_wilson_decimal(n_sys, k_sys, n_inh_all, k_inh_all)
    cz_all_lo, cz_all_hi = _chan_zhang_ci95_decimal(n_sys, k_sys, n_inh_all, k_inh_all)

    expected_supported = (
        (n_sys, n_inh_pure, n_inh_all) == (54, 15, 17)
        and Fraction(k_sys, n_sys) > Fraction(k_inh_pure, n_inh_pure)
        and p_pure_1s < Fraction(1, 100)
        and cz_pure_lo > Decimal(0)
        and nw_pure_lo > Decimal(0)
        and Fraction(k_sys, n_sys) > Fraction(k_inh_all, n_inh_all)
        and p_all_1s < Fraction(1, 100)
    )
    expected_outcome = "supported" if expected_supported else "refuted"

    res_doc = json.loads(results_path.read_text(encoding="utf-8"))
    prim_res = res_doc.get("primary_estimand_pure_route_randomization_and_masking") or {}
    sec_res = res_doc.get("secondary_gate_all_17_inhaled_including_dual_route") or {}
    sc_res = res_doc.get("support_criteria_evaluation") or {}

    if res_doc.get("done_sha256") != done_sha:
        mismatches.append(f"done_sha256 mismatch: {res_doc.get('done_sha256')} != {done_sha}")
    if res_doc.get("pred0015_results_sha256") != p15_sha:
        mismatches.append("pred0015_results_sha256 mismatch")
    if res_doc.get("fixture_sha256") != fixture_sha:
        mismatches.append("fixture_sha256 mismatch")
    if prim_res.get("k_sys") != k_sys or prim_res.get("k_inh_pure") != k_inh_pure:
        mismatches.append(f"Primary counts mismatch: {prim_res.get('k_sys')}/{prim_res.get('k_inh_pure')} != {k_sys}/{k_inh_pure}")
    if sec_res.get("k_inh_all") != k_inh_all:
        mismatches.append(f"All-17 inhaled count mismatch: {sec_res.get('k_inh_all')} != {k_inh_all}")
    if prim_res.get("systemic_matched_ncts") != sys_matched:
        mismatches.append("systemic_matched_ncts mismatch")
    if prim_res.get("inhaled_pure_matched_ncts") != inh_pure_matched:
        mismatches.append("inhaled_pure_matched_ncts mismatch")
    if prim_res.get("p_fisher_one_sided_exact") != f"{p_pure_1s.numerator}/{p_pure_1s.denominator}":
        mismatches.append("p_fisher_one_sided_exact pure mismatch")
    if prim_res.get("p_fisher_two_sided_exact") != f"{p_pure_2s.numerator}/{p_pure_2s.denominator}":
        mismatches.append("p_fisher_two_sided_exact pure mismatch")
    if sec_res.get("p_fisher_one_sided_exact") != f"{p_all_1s.numerator}/{p_all_1s.denominator}":
        mismatches.append("p_fisher_one_sided_exact all-17 mismatch")

    tol = Decimal("0.000001")
    for label, got_pair, exp_pair in (
        ("chan_zhang_pure", prim_res.get("chan_zhang_unconditional_ci95"), (cz_pure_lo, cz_pure_hi)),
        ("newcombe_pure", prim_res.get("newcombe_wilson_ci95"), (nw_pure_lo, nw_pure_hi)),
        ("chan_zhang_all17", sec_res.get("chan_zhang_unconditional_ci95"), (cz_all_lo, cz_all_hi)),
        ("newcombe_all17", sec_res.get("newcombe_wilson_ci95"), (nw_all_lo, nw_all_hi)),
    ):
        if not isinstance(got_pair, list) or len(got_pair) != 2:
            mismatches.append(f"{label}: missing CI pair")
            continue
        for idx, bound_name in enumerate(("low", "high")):
            diff = abs(Decimal(str(got_pair[idx])) - exp_pair[idx])
            if diff > tol:
                mismatches.append(f"{label}.{bound_name} drift {diff} > {tol}: got={got_pair[idx]} exp={exp_pair[idx]}")

    if sc_res.get("outcome") != expected_outcome:
        mismatches.append(f"outcome mismatch: {sc_res.get('outcome')} != {expected_outcome}")

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
    ap.add_argument("--fixture", type=Path, default=DEFAULT_FIXTURE)
    ap.add_argument("--results", type=Path, default=DEFAULT_RESULTS)
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--json", action="store_true")
    args = ap.parse_args(argv)

    rep = run_independent_audit(args.done, args.pred0015_results, args.fixture, args.results)
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
