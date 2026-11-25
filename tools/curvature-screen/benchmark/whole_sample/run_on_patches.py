#!/usr/bin/env python3
"""Apply the fixed whole-sample analysis to a table of burst durations.

Implements all steps of specs/CFTR_DATA_REQUEST_AND_ANALYSIS_PLAN.md:
- Step 1: one-active-channel verification, unstable-patch exclusion, and per-construct burst minimums (50 WT, 5 control).
- Step 2: dead-time rule with user-supplied resolution or 0.36/f_c default filter resolution, and frame-grid adjustment.
- Steps 3-5: CV^2 statistic, per-patch critical threshold (including T-0923 instrument calibration for dead time,
  filter cutoff f_c, and burst criterion t_crit), and sqrt(N)-weighted cross-patch combination.
- Step 6: composite-null generalized likelihood-ratio comparator (positive mixture of up to 2 exponentials vs.
  two irreversible steps in sequence) with parametric bootstrap and secondary single-exponential LRT.
- Step 7: control power simulation at observed burst counts and Section 5 failure criteria checks.
- Section 5: 1x / 2x resolution sensitivity decision rule and shortest-burst gap diagnostic.

Input : CSV with a header and the columns  patch_id,duration  (seconds), optionally  construct,
        n_channels (or active_channels), unstable, exclude_reason. Optional --metadata JSON mapping
        patch_id -> {n_channels, unstable, reason, dead_time, filter_cutoff, burst_criterion, quantum}.
Output: JSON. Per construct and per patch: exclusions, bursts used, CV^2, its 1% critical value,
        composite-null GLRT comparator, combined results at 1x and 2x resolutions, control power
        simulation, and Section 5 failure checks.

    python benchmark/whole_sample/run_on_patches.py bursts.csv --dead-time 0.002 --control E1371Q --out result.json
    python benchmark/whole_sample/run_on_patches.py bursts.csv --filter-cutoff 1000 --burst-criterion 0.005
    python benchmark/whole_sample/run_on_patches.py --self-test

If durations are recorded on a frame grid, pass --quantum (the frame interval, seconds); durations must then lie
within 1% of that grid; the dead time is rounded up to the grid and one further frame is discarded,
because the first frame above the resolution is only partly observed.  Nothing in this file has been run on a CFTR recording.
"""
import argparse
import csv
import hashlib
import json
import math
import sys
from collections import defaultdict
import numpy as np
from scipy import optimize

SEED = 20261113
REPS = 20000
BOOTSTRAP_REPS = 2000
POWER_REPS = 2000
ALPHA = 0.01
MIN_BURSTS = 50            # plan step 1, wild type
MIN_BURSTS_CONTROL = 5     # plan step 1, control construct
VARIANT_A_RATES = (12.0, 8.0)  # reference sequential two-step rates (rate ratio 1.5)


def cv2(x, axis=-1):
    m = x.mean(axis=axis)
    return x.var(axis=axis, ddof=1) / m**2


def _stable_id_int(text):
    digest = hashlib.sha256(str(text).encode("utf-8")).digest()
    return int.from_bytes(digest[:8], "big")


class Null:
    """Monte Carlo null of CV^2 for a single exponential or passive mixture calibrated for instrument
    dead time, filter cutoff (0.36/f_c), and burst criterion (t_crit).
    Each (N, quantum/mean, dead_time, filter_cutoff, burst_criterion, dimensionless ratios) has its own
    generator, so the result does not depend on the order in which patches are read."""

    def __init__(self, seed=SEED):
        self.seed, self.cache = seed, {}

    def draws(self, n, q_over_mean=0.0, dead_time=0.0, filter_cutoff=None, burst_criterion=None, mean_duration=1.0):
        fc_val = float(filter_cutoff) if filter_cutoff and filter_cutoff > 0 else 0.0
        tc_val = float(burst_criterion) if burst_criterion and burst_criterion > 0 else 0.0
        dt_val = float(dead_time) if dead_time and dead_time > 0 else 0.0
        is_cal = (tc_val > 0.0 or fc_val > 0.0)
        eff_dt = max(dt_val, 0.36 / fc_val if fc_val > 0 else 0.0)
        if is_cal:
            mean_val = float(mean_duration) if (mean_duration is not None and mean_duration > 0) else 0.010
            eff_dt_ratio_k = int(round((eff_dt / mean_val) * 1e6))
            tc_ratio_k = int(round((tc_val / mean_val) * 1e6))
            mean_log_k = int(round((math.log10(max(mean_val, 1e-18)) + 20.0) * 1e4))
        else:
            mean_val = 1.0
            eff_dt_ratio_k = 0
            tc_ratio_k = 0
            mean_log_k = 0
        key = (
            int(n),
            int(round(float(q_over_mean) * 1e4)),
            int(round(dt_val * 1e6)),
            int(round(fc_val * 10)),
            int(round(tc_val * 1e6)),
            eff_dt_ratio_k,
            tc_ratio_k,
            mean_log_k,
        )
        if key not in self.cache:
            seed_elems = key if is_cal else key[:5]
            rng = np.random.default_rng([self.seed, *seed_elems])
            n_req = key[0]
            if is_cal:
                # Instrument-calibrated simulation: simulate passive gating bursts through
                # dead time, filter cutoff, and burst criterion.
                # Under a passive single-channel model with mean open duration tau_o and brief intra-burst
                # closures tau_c, closures > t_crit terminate the burst.
                # Missed intra-burst closures (< t_res) fuse adjacent openings.
                tau_o = max(eff_dt * 2.0, mean_val)
                tau_c = max(eff_dt * 0.5, 0.001)

                p_end = math.exp(-tc_val / tau_c) if tc_val > 0 else 1.0
                p_end = max(1e-4, min(1.0, p_end))

                p_miss = 1.0 - math.exp(-eff_dt / tau_c) if eff_dt > 0 else 0.0
                p_term_eff = p_end / max(1e-6, 1.0 - (1.0 - p_end) * p_miss)
                p_term_eff = max(1e-4, min(1.0, p_term_eff))

                tau_o_eff = tau_o + (1.0 - p_end) * p_miss * tau_c

                if eff_dt > 0.0:
                    # Proper dead-time truncation: keep bursts >= eff_dt and subtract eff_dt
                    # until every replicate has exactly n_req valid post-dead-time bursts.
                    x = np.empty((REPS, n_req), dtype=float)
                    filled = np.zeros(REPS, dtype=int)
                    todo = np.arange(REPS)
                    trunc_iters = 0
                    while todo.size > 0:
                        trunc_iters += 1
                        if trunc_iters > 100:
                            raise RuntimeError("dead-time truncation failed to fill replicates within 100 iterations")
                        need = int((n_req - filled[todo]).max())
                        buf = int(math.ceil(need * 2.2)) + 32
                        k_ops = rng.geometric(p_term_eff, size=(todo.size, buf))
                        bursts = rng.gamma(k_ops, max(1e-6, tau_o_eff))
                        mask = bursts >= eff_dt
                        counts = mask.sum(axis=1)
                        rem = n_req - filled[todo]
                        if np.all(filled[todo] == 0) and np.all(counts >= rem):
                            csum = np.cumsum(mask, axis=1)
                            keep = mask & (csum <= n_req)
                            x[todo] = (bursts[keep] - eff_dt).reshape(todo.size, n_req)
                            break
                        for idx, row_id in enumerate(todo):
                            surv = bursts[idx, mask[idx]] - eff_dt
                            take = min(int(surv.size), int(n_req - filled[row_id]))
                            if take > 0:
                                x[row_id, filled[row_id]:filled[row_id] + take] = surv[:take]
                                filled[row_id] += take
                        todo = np.flatnonzero(filled < n_req)
                else:
                    k_ops = rng.geometric(p_term_eff, size=(REPS, n_req))
                    x = rng.gamma(k_ops, max(1e-6, tau_o_eff))

                if key[1] > 0:
                    q_norm = key[1] / 1e4
                    sim_mean = float(x.mean())
                    if sim_mean > 0.0:
                        x_norm = x / sim_mean
                        x = np.floor(x_norm / q_norm) * q_norm + q_norm / 2
                self.cache[key] = cv2(x)
            else:
                x = rng.exponential(1.0, (REPS, n_req))
                if key[1] > 0:
                    q_norm = key[1] / 1e4
                    x = np.floor(x / q_norm) * q_norm + q_norm / 2
                self.cache[key] = cv2(x)
        return self.cache[key]


def adjust(x, dead_time, quantum):
    """Plan step 2: discard bursts shorter than the resolution and subtract it. With a frame grid the resolution
    is rounded up to the grid, one further frame is discarded (the first frame above the resolution is only
    partly observed, whatever the rounding convention), and half a frame is added to match the quantised null."""
    if quantum > 0:
        snapped = np.round(x / quantum) * quantum
        if not np.allclose(x, snapped, rtol=0, atol=0.01 * quantum):
            raise ValueError("durations are not on the stated frame grid; omit --quantum or correct it")
        x = snapped
        d = quantum * (math.ceil(dead_time / quantum - 1e-9) + 1)
        return x[x >= d - 1e-9 * quantum] - d + quantum / 2
    return x[x >= dead_time] - dead_time


def patch_stats(x, quantum, null, dead_time=0.0, filter_cutoff=None, burst_criterion=None):
    n = int(x.size)
    mu = float(x.mean())
    q_over_m = quantum / mu if (quantum > 0 and mu > 0) else 0.0
    d = null.draws(
        n,
        q_over_m,
        dead_time=dead_time,
        filter_cutoff=filter_cutoff,
        burst_criterion=burst_criterion,
        mean_duration=mu,
    )
    stat = float(cv2(x))
    sd = float(d.std())
    d_mean = float(d.mean())
    z = float((stat - d_mean) / sd) if sd > 0 else 0.0
    null_z = (d - d_mean) / sd if sd > 0 else np.zeros_like(d)
    return {
        "n": n,
        "cv2": stat,
        "crit": float(np.quantile(d, ALPHA)),
        "p": float((1 + (d <= stat).sum()) / (1 + REPS)),
        "z": z,
        "null_z": null_z,
        "null_mean": d_mean,
        "null_sd": sd,
    }


def combine(stats, seed_offset):
    """Plan step 5: sqrt(N)-weighted sum of standardised scores, calibrated jointly with the actual burst counts."""
    ns = np.array([s["n"] for s in stats], dtype=float)
    w = np.sqrt(ns)
    rng = np.random.default_rng([SEED, 99, seed_offset])
    null = sum(wi * rng.permutation(s["null_z"]) for wi, s in zip(w, stats)) / math.sqrt(ns.sum())
    comb = float(sum(wi * s["z"] for wi, s in zip(w, stats)) / math.sqrt(ns.sum()))
    crit = float(np.quantile(null, ALPHA))
    return {"statistic": round(comb, 4), "critical_alpha_0.01": round(crit, 4),
            "p_one_sided": round(float((1 + (null <= comb).sum()) / (1 + REPS)), 5), "rejects_mixture": bool(comb < crit)}


# ------------------------------------------------------------------ Plan Step 6: Composite-null GLRT comparator

def _phi_and_log(z):
    """Smooth evaluation of phi(z) = (1 - exp(-z)) / z and log(phi(z)) for z >= 0."""
    z_arr = np.asarray(z, dtype=float)
    small = z_arr < 1e-4
    phi = np.empty_like(z_arr)
    log_phi = np.empty_like(z_arr)
    if np.any(small):
        zs = z_arr[small]
        ps = 1.0 - zs * (0.5 - zs * (1.0 / 6.0 - zs * (1.0 / 24.0 - zs / 120.0)))
        phi[small] = ps
        log_phi[small] = np.log(ps)
    if np.any(~small):
        zl = z_arr[~small]
        num = -np.expm1(-zl)
        phi[~small] = num / zl
        log_phi[~small] = np.log(num) - np.log(zl)
    return phi, log_phi


def _phi_log_deriv_scaled(z):
    """Compute z * phi'(z) / phi(z) smoothly for z >= 0."""
    z_arr = np.asarray(z, dtype=float)
    small = z_arr < 1e-4
    out = np.empty_like(z_arr)
    if np.any(small):
        zs = z_arr[small]
        out[small] = zs * (-0.5 + zs / 12.0 - (zs ** 3) / 720.0)
    if np.any(~small):
        zl = z_arr[~small]
        exp_neg = np.exp(-zl)
        denom = -np.expm1(-zl)
        out[~small] = (zl * exp_neg) / denom - 1.0
    return out


def _prepare_glrt_data(x, quantum):
    """Scale data by sample mean (or extract discrete bin counts for quantized data)."""
    x_arr = np.asarray(x, dtype=float)
    n = int(x_arr.size)
    mu = float(x_arr.mean())
    if mu <= 0.0:
        raise ValueError("mean burst duration must be positive for GLRT comparator")
    if quantum > 0.0:
        k_idx = np.round(x_arr / quantum - 0.5).astype(np.int64)
        k_idx = np.maximum(k_idx, 0)
        k_max = int(k_idx.max()) if k_idx.size else 0
        if k_max < 20000:
            bc = np.bincount(k_idx)
            nz = np.flatnonzero(bc)
            ks = nz
            cs = bc[nz]
        else:
            ks, cs = np.unique(k_idx, return_counts=True)
        q_u = float(quantum) / mu
        return {
            "quantized": True,
            "n": n,
            "mu": mu,
            "quantum": float(quantum),
            "q_u": q_u,
            "ks": ks.astype(float),
            "cs": cs.astype(float),
            "k_bar": float(np.dot(ks, cs)) / float(n),
        }
    u = np.maximum(x_arr / mu, 1e-15)
    return {
        "quantized": False,
        "n": n,
        "mu": mu,
        "quantum": 0.0,
        "u": u,
        "log_mu": math.log(mu),
        "cv2": float(cv2(x_arr)),
    }


def _fit_exp_prepared(prep):
    """Exact 1-exponential MLE on prepared data; returns (ll, rate_per_s)."""
    n = prep["n"]
    mu = prep["mu"]
    if prep["quantized"]:
        k_bar = prep["k_bar"]
        q_u = prep["q_u"]
        if k_bar <= 0.0:
            return 0.0, 50.0 / (q_u * mu)
        p = 1.0 / (1.0 + k_bar)
        rate_u = -math.log1p(-p) / q_u
        ll = float(n) * (k_bar * math.log1p(-p) + math.log(p))
        return ll, rate_u / mu
    ll = -float(n) * (1.0 + prep["log_mu"])
    return ll, 1.0 / mu


def _eval_mix_nll_and_grad(theta, prep):
    """Negative log-likelihood and gradient wrt [logit(w), log(r1), log(r2)] in scaled units."""
    lw_raw, lr1, lr2 = float(theta[0]), float(theta[1]), float(theta[2])
    w = 1.0 / (1.0 + math.exp(-lw_raw))
    r1 = math.exp(lr1)
    r2 = math.exp(lr2)
    log_w = -math.log1p(math.exp(-lw_raw)) if lw_raw >= 0 else lw_raw - math.log1p(math.exp(lw_raw))
    log_1mw = -lw_raw - math.log1p(math.exp(-lw_raw)) if lw_raw >= 0 else -math.log1p(math.exp(lw_raw))

    if prep["quantized"]:
        q_u = prep["q_u"]
        ks = prep["ks"]
        cs = prep["cs"]
        z1 = r1 * q_u
        z2 = r2 * q_u
        p1 = -math.expm1(-z1) if z1 < 50.0 else 1.0
        p2 = -math.expm1(-z2) if z2 < 50.0 else 1.0
        p1 = max(p1, 1e-300)
        p2 = max(p2, 1e-300)
        a1 = log_w + math.log(p1) - ks * z1
        a2 = log_1mw + math.log(p2) - ks * z2
        m = np.maximum(a1, a2)
        e1 = np.exp(a1 - m)
        e2 = np.exp(a2 - m)
        s = e1 + e2
        ll = float(np.dot(cs, m + np.log(s)))
        gamma = e1 / s
        cg1 = cs * gamma
        cg2 = cs - cg1
        n1 = float(np.sum(cg1))
        n2 = float(prep["n"]) - n1
        d_lw = n1 - float(prep["n"]) * w
        fac1 = (math.exp(-z1) / p1) if z1 < 50.0 else 0.0
        fac2 = (math.exp(-z2) / p2) if z2 < 50.0 else 0.0
        d_lr1 = z1 * (n1 * fac1 - float(np.dot(cg1, ks)))
        d_lr2 = z2 * (n2 * fac2 - float(np.dot(cg2, ks)))
        return -ll, np.array([-d_lw, -d_lr1, -d_lr2], dtype=float)

    u = prep["u"]
    a1 = log_w + lr1 - r1 * u
    a2 = log_1mw + lr2 - r2 * u
    m = np.maximum(a1, a2)
    e1 = np.exp(a1 - m)
    e2 = np.exp(a2 - m)
    s = e1 + e2
    ll_u = float(np.sum(m + np.log(s)))
    gamma = e1 / s
    n1 = float(np.sum(gamma))
    n2 = float(prep["n"]) - n1
    d_lw = n1 - float(prep["n"]) * w
    d_lr1 = n1 - r1 * float(np.dot(gamma, u))
    d_lr2 = n2 - r2 * float(np.dot(1.0 - gamma, u))
    return -ll_u, np.array([-d_lw, -d_lr1, -d_lr2], dtype=float)


def _em_mix(prep, w, r1, r2, n_iters=25):
    """Fast vectorized EM for 2-exponential mixture (continuous or discrete quantized)."""
    n_f = float(prep["n"])
    if prep["quantized"]:
        q_u = prep["q_u"]
        ks = prep["ks"]
        cs = prep["cs"]
        max_r = min(1e4, 40.0 / max(q_u, 1e-6))
        min_r = 1e-4
        best_ll = -float("inf")
        best = (w, r1, r2)
        for _ in range(n_iters):
            z1 = r1 * q_u
            z2 = r2 * q_u
            p1 = max(-math.expm1(-z1) if z1 < 50.0 else 1.0, 1e-300)
            p2 = max(-math.expm1(-z2) if z2 < 50.0 else 1.0, 1e-300)
            a1 = math.log(max(w, 1e-12)) + math.log(p1) - ks * z1
            a2 = math.log(max(1.0 - w, 1e-12)) + math.log(p2) - ks * z2
            m = np.maximum(a1, a2)
            e1 = np.exp(a1 - m)
            e2 = np.exp(a2 - m)
            s = e1 + e2
            ll = float(np.dot(cs, m + np.log(s)))
            if ll > best_ll:
                best_ll = ll
                best = (w, r1, r2)
            gamma = e1 / s
            cg1 = cs * gamma
            cg2 = cs - cg1
            n1 = float(np.sum(cg1))
            n2 = n_f - n1
            if n1 <= 1e-6 or n2 <= 1e-6:
                break
            w_new = min(max(n1 / n_f, 1e-6), 1.0 - 1e-6)
            k1 = float(np.dot(cg1, ks))
            k2 = float(np.dot(cg2, ks))
            p1_new = min(max(n1 / max(n1 + k1, 1e-12), 1e-12), 1.0 - 1e-12)
            p2_new = min(max(n2 / max(n2 + k2, 1e-12), 1e-12), 1.0 - 1e-12)
            r1_new = min(max(-math.log1p(-p1_new) / q_u, min_r), max_r)
            r2_new = min(max(-math.log1p(-p2_new) / q_u, min_r), max_r)
            if abs(w_new - w) < 1e-7 and abs(r1_new - r1) < 1e-6 * r1 and abs(r2_new - r2) < 1e-6 * r2:
                w, r1, r2 = w_new, r1_new, r2_new
                break
            w, r1, r2 = w_new, r1_new, r2_new
        return best_ll, best

    u = prep["u"]
    best_ll_u = -float("inf")
    best = (w, r1, r2)
    for _ in range(n_iters):
        a1 = math.log(max(w, 1e-12)) + math.log(r1) - r1 * u
        a2 = math.log(max(1.0 - w, 1e-12)) + math.log(r2) - r2 * u
        m = np.maximum(a1, a2)
        e1 = np.exp(a1 - m)
        e2 = np.exp(a2 - m)
        s = e1 + e2
        ll_u = float(np.sum(m + np.log(s)))
        if ll_u > best_ll_u:
            best_ll_u = ll_u
            best = (w, r1, r2)
        gamma = e1 / s
        n1 = float(np.sum(gamma))
        n2 = n_f - n1
        if n1 <= 1e-6 or n2 <= 1e-6:
            break
        w_new = min(max(n1 / n_f, 1e-6), 1.0 - 1e-6)
        s1 = float(np.dot(gamma, u))
        s2 = float(np.dot(1.0 - gamma, u))
        r1_new = min(max(n1 / max(s1, 1e-12), 1e-4), 1e4)
        r2_new = min(max(n2 / max(s2, 1e-12), 1e-4), 1e4)
        if abs(w_new - w) < 1e-7 and abs(r1_new - r1) < 1e-6 * r1 and abs(r2_new - r2) < 1e-6 * r2:
            w, r1, r2 = w_new, r1_new, r2_new
            break
        w, r1, r2 = w_new, r1_new, r2_new
    return best_ll_u, best


def _fit_mixture_prepared(prep, ll_exp, rate_exp_s, init_mix_s=None, thorough=True):
    """Fit positive mixture of up to 2 exponentials; guarantees ll_mix >= ll_exp."""
    mu = prep["mu"]
    re_u = rate_exp_s * mu
    offset = 0.0 if prep["quantized"] else -float(prep["n"]) * prep["log_mu"]

    seeds = []
    if init_mix_s is not None:
        w0, r10_s, r20_s = init_mix_s
        if 1e-5 < w0 < 1.0 - 1e-5 and abs(r10_s - r20_s) > 1e-6 * max(r10_s, r20_s):
            seeds.append((float(w0), float(r10_s) * mu, float(r20_s) * mu))
    if thorough:
        seeds.extend([
            (0.5, 0.35 * re_u, 2.5 * re_u),
            (0.2, 0.25 * re_u, 4.0 * re_u),
            (0.8, 0.5 * re_u, 3.0 * re_u),
            (0.5, 0.7 * re_u, 1.4 * re_u),
        ])
    else:
        if not seeds:
            seeds.append((0.5, 0.5 * re_u, 2.0 * re_u))
        seeds.append((0.3, 0.3 * re_u, 3.0 * re_u))

    best_ll = ll_exp
    best_params = (1.0, re_u, re_u)
    candidates = []
    em_iters = 35 if thorough else 20
    for w0, r10, r20 in seeds:
        r10 = min(max(r10, 1e-4), 1e4)
        r20 = min(max(r20, 1e-4), 1e4)
        ll_em, pars = _em_mix(prep, w0, r10, r20, n_iters=em_iters)
        ll_abs = ll_em + offset
        candidates.append((ll_abs, pars))
        if ll_abs > best_ll:
            best_ll = ll_abs
            best_params = pars

    # Polish the best EM candidate(s) with analytical-gradient L-BFGS-B
    candidates.sort(key=lambda item: item[0], reverse=True)
    max_lr = math.log(min(1e4, 40.0 / max(prep["q_u"], 1e-6))) if prep["quantized"] else 9.0
    bounds = [(-11.5, 11.5), (-9.0, max_lr), (-9.0, max_lr)]
    n_polish = 2 if thorough else 1
    for _, (w_c, r1_c, r2_c) in candidates[:n_polish]:
        w_c = min(max(w_c, 1e-5), 1.0 - 1e-5)
        x0 = np.array([
            math.log(w_c / (1.0 - w_c)),
            math.log(min(max(r1_c, 1e-4), math.exp(max_lr))),
            math.log(min(max(r2_c, 1e-4), math.exp(max_lr))),
        ], dtype=float)
        res = optimize.minimize(
            _eval_mix_nll_and_grad,
            x0,
            args=(prep,),
            method="L-BFGS-B",
            jac=True,
            bounds=bounds,
            options={"maxiter": 50 if thorough else 30, "ftol": 1e-11, "gtol": 1e-7},
        )
        ll_opt = -float(res.fun) + offset
        if ll_opt > best_ll:
            best_ll = ll_opt
            lw, lr1, lr2 = res.x
            best_params = (1.0 / (1.0 + math.exp(-float(lw))), math.exp(float(lr1)), math.exp(float(lr2)))

    w_b, r1_b, r2_b = best_params
    if r1_b > r2_b:
        w_b, r1_b, r2_b = 1.0 - w_b, r2_b, r1_b
    return best_ll, (float(w_b), float(r1_b / mu), float(r2_b / mu))


def _eval_hypo_nll_continuous(ab, prep):
    """Negative log-likelihood and analytical gradient for continuous 2-step hypoexponential."""
    a, b = float(ab[0]), float(ab[1])
    lam1 = math.exp(a)
    delta = math.exp(b)
    lam2 = lam1 * (1.0 + delta)
    u = prep["u"]
    z = (lam1 * delta) * u
    _, log_phi = _phi_and_log(z)
    g = _phi_log_deriv_scaled(z)
    n_f = float(prep["n"])
    # ln f1(u) = 2*a + ln(1+delta) + ln(u) - lam1*u + ln(phi(z))
    log_u_sum = prep.get("log_u_sum")
    if log_u_sum is None:
        log_u_sum = float(np.sum(np.log(u)))
        prep["log_u_sum"] = log_u_sum
    ll_u = n_f * (2.0 * a + math.log1p(delta)) + log_u_sum - lam1 * n_f + float(np.sum(log_phi))
    g_sum = float(np.sum(g))
    d_a = (2.0 - lam1) * n_f + g_sum
    d_b = n_f * (delta / (1.0 + delta)) + g_sum
    return -ll_u, np.array([-d_a, -d_b], dtype=float)


def _eval_hypo_nll_quantized(ab, prep):
    """Negative log-likelihood for quantized 2-step hypoexponential via exact S(k*q) - S((k+1)*q)."""
    a, b = float(ab[0]), float(ab[1])
    lam1 = math.exp(a)
    delta = math.exp(b)
    lam_diff = lam1 * delta
    q_u = prep["q_u"]
    ks = prep["ks"]
    cs = prep["cs"]

    t0 = ks * q_u
    t1 = (ks + 1.0) * q_u
    phi0, _ = _phi_and_log(lam_diff * t0)
    phi1, _ = _phi_and_log(lam_diff * t1)
    e_step = math.exp(-lam1 * q_u)
    bracket = (1.0 + lam1 * t0 * phi0) - e_step * (1.0 + lam1 * t1 * phi1)
    log_p = -lam1 * t0 + np.log(np.maximum(bracket, 1e-300))
    return -float(np.dot(cs, log_p))


def _fit_hypo_prepared(prep, ll_exp, rate_exp_s, thorough=True):
    """Fit 2-step hypoexponential (sequential irreversible steps) with exact Erlang-2 limit."""
    mu = prep["mu"]
    offset = 0.0 if prep["quantized"] else -float(prep["n"]) * prep["log_mu"]

    # Method-of-moments seed in normalized units (mean == 1)
    seeds = []
    if not prep["quantized"]:
        c = prep["cv2"]
        if c <= 0.5:
            seeds.append((math.log(2.0), -16.0))
        elif c < 0.98:
            r = max(1e-3, min(0.999, (c - math.sqrt(max(0.0, 2.0 * c - 1.0))) / (1.0 - c)))
            lam1_0 = 1.0 + r
            delta_0 = max(1e-5, 1.0 / r - 1.0)
            seeds.append((math.log(lam1_0), math.log(delta_0)))
        else:
            seeds.append((math.log(1.1), math.log(9.0)))
    else:
        seeds.append((math.log(1.6), math.log(0.6)))
        seeds.append((math.log(2.0), -16.0))

    if thorough:
        seeds.extend([
            (math.log(2.0), -16.0),
            (math.log(1.5), math.log(1.0)),
            (math.log(1.2), math.log(4.0)),
        ])

    max_a = math.log(min(200.0, 40.0 / max(prep["q_u"], 1e-6))) if prep["quantized"] else 6.0
    bounds = [(-4.0, max_a), (-18.0, 12.0)]

    best_ll = ll_exp
    best_rates_u = (rate_exp_s * mu, 1e6 * rate_exp_s * mu)
    n_starts = len(seeds) if thorough else 1
    for seed_ab in seeds[:n_starts]:
        x0 = np.array([min(max(seed_ab[0], bounds[0][0]), bounds[0][1]),
                       min(max(seed_ab[1], bounds[1][0]), bounds[1][1])], dtype=float)
        if prep["quantized"]:
            res = optimize.minimize(
                _eval_hypo_nll_quantized,
                x0,
                args=(prep,),
                method="L-BFGS-B",
                bounds=bounds,
                options={"maxiter": 60 if thorough else 35, "ftol": 1e-11},
            )
        else:
            res = optimize.minimize(
                _eval_hypo_nll_continuous,
                x0,
                args=(prep,),
                method="L-BFGS-B",
                jac=True,
                bounds=bounds,
                options={"maxiter": 60 if thorough else 35, "ftol": 1e-11, "gtol": 1e-7},
            )
        ll_opt = -float(res.fun) + offset
        if ll_opt > best_ll:
            best_ll = ll_opt
            lam1_u = math.exp(float(res.x[0]))
            delta_u = math.exp(float(res.x[1]))
            lam2_u = lam1_u * (1.0 + delta_u)
            best_rates_u = (lam1_u, lam2_u)

    k1_s = best_rates_u[0] / mu
    k2_s = best_rates_u[1] / mu
    return best_ll, (float(min(k1_s, k2_s)), float(max(k1_s, k2_s)))


def fit_glrt_models(x, quantum=0.0, init_mix=None, thorough=True, skip_mix_if_hypo_le_exp=False):
    """Evaluate 1-exp, 2-exp mixture (H0), and 2-step hypoexponential (H1) MLEs and GLRT statistics."""
    prep = _prepare_glrt_data(x, quantum)
    ll_exp, rate_exp = _fit_exp_prepared(prep)
    ll_hypo, hypo_rates = _fit_hypo_prepared(prep, ll_exp, rate_exp, thorough=thorough)
    lambda_exp = 2.0 * max(0.0, ll_hypo - ll_exp)
    if skip_mix_if_hypo_le_exp and ll_hypo <= ll_exp + 1e-12:
        # Since H0 contains 1-exp, ll_mix >= ll_exp >= ll_hypo implies lambda_comp == 0 identically.
        return {
            "ll_exp": ll_exp,
            "ll_mix": ll_exp,
            "ll_hypo": ll_hypo,
            "lambda_comp": 0.0,
            "lambda_exp": lambda_exp,
            "mix_params": (1.0, rate_exp, rate_exp),
            "hypo_params": hypo_rates,
        }
    ll_mix, mix_params = _fit_mixture_prepared(prep, ll_exp, rate_exp, init_mix_s=init_mix, thorough=thorough)
    lambda_comp = 2.0 * max(0.0, ll_hypo - ll_mix)
    return {
        "ll_exp": ll_exp,
        "ll_mix": ll_mix,
        "ll_hypo": ll_hypo,
        "lambda_comp": lambda_comp,
        "lambda_exp": lambda_exp,
        "mix_params": mix_params,
        "hypo_params": hypo_rates,
    }


def _sample_mixture_draws(rng, w, r1, r2, reps, n, quantum):
    """Draw synthetic samples from fitted H0 mixture with the exact same n and discretization."""
    if w >= 1.0 - 1e-9 or abs(r1 - r2) <= 1e-9 * max(r1, r2):
        y = rng.exponential(1.0 / r1, size=(reps, n))
    elif w <= 1e-9:
        y = rng.exponential(1.0 / r2, size=(reps, n))
    else:
        pick = rng.random((reps, n)) < w
        y = np.where(pick, rng.exponential(1.0 / r1, size=(reps, n)), rng.exponential(1.0 / r2, size=(reps, n)))
    if quantum > 0.0:
        y = np.floor(y / quantum) * quantum + quantum / 2.0
    return y


def comparator_patch_stats(x, quantum, bootstrap_reps=BOOTSTRAP_REPS, seed_parts=None):
    """Fit composite-null GLRT comparator on one patch and calibrate by parametric bootstrap."""
    obs = fit_glrt_models(x, quantum=quantum, thorough=True, skip_mix_if_hypo_le_exp=False)
    w_hat, r1_hat, r2_hat = obs["mix_params"]
    n = int(len(x))
    rng = np.random.default_rng([SEED, 600, *(seed_parts or (n,))])
    boot_lambdas = np.empty(bootstrap_reps, dtype=float)
    batch_size = max(1, min(bootstrap_reps, 1_500_000 // max(1, n)))
    done = 0
    while done < bootstrap_reps:
        cur = min(batch_size, bootstrap_reps - done)
        sims = _sample_mixture_draws(rng, w_hat, r1_hat, r2_hat, cur, n, float(quantum))
        for i in range(cur):
            res_b = fit_glrt_models(
                sims[i],
                quantum=quantum,
                init_mix=(w_hat, r1_hat, r2_hat),
                thorough=False,
                skip_mix_if_hypo_le_exp=True,
            )
            boot_lambdas[done + i] = res_b["lambda_comp"]
        done += cur
    crit = float(np.quantile(boot_lambdas, 1.0 - ALPHA))
    p_val = float((1 + np.sum(boot_lambdas >= obs["lambda_comp"])) / (1 + bootstrap_reps))
    rejects = bool(obs["lambda_comp"] > crit and obs["lambda_comp"] > 0.0)
    return {
        "n": n,
        "ll_exp": obs["ll_exp"],
        "ll_mix": obs["ll_mix"],
        "ll_hypo": obs["ll_hypo"],
        "lambda_comp": obs["lambda_comp"],
        "lambda_exp": obs["lambda_exp"],
        "mix_params": obs["mix_params"],
        "hypo_params": obs["hypo_params"],
        "bootstrap_reps": int(bootstrap_reps),
        "critical_alpha_0.01": crit,
        "p_one_sided": p_val,
        "rejects_mixture": rejects,
        "boot_lambdas": boot_lambdas,
    }


def _format_patch_comparator(comp):
    return {
        "ll_exp": round(comp["ll_exp"], 6),
        "ll_mix": round(comp["ll_mix"], 6),
        "ll_hypo": round(comp["ll_hypo"], 6),
        "lambda_comp": round(comp["lambda_comp"], 6),
        "lambda_exp": round(comp["lambda_exp"], 6),
        "fitted_h0_mixture": {
            "w": round(comp["mix_params"][0], 6),
            "rate1_per_s": round(comp["mix_params"][1], 6),
            "rate2_per_s": round(comp["mix_params"][2], 6),
        },
        "fitted_h1_hypoexponential": {
            "rate1_per_s": round(comp["hypo_params"][0], 6),
            "rate2_per_s": round(comp["hypo_params"][1], 6),
        },
        "bootstrap_reps": int(comp["bootstrap_reps"]),
        "critical_alpha_0.01": round(comp["critical_alpha_0.01"], 6),
        "p_one_sided": round(comp["p_one_sided"], 5),
        "rejects_mixture": bool(comp["rejects_mixture"]),
        "secondary_single_exp_lrt": {
            "lambda_exp": round(comp["lambda_exp"], 6),
            "note": "descriptive only; rejects positive exponential mixtures as well as sequential steps",
        },
    }


def combine_comparator(comp_list):
    """Calibrate combined cross-patch GLRT sum sum(Lambda_i) against joint bootstrap null draws."""
    if not comp_list:
        return None
    boot_sum = np.sum([c["boot_lambdas"] for c in comp_list], axis=0)
    obs_comp_sum = float(sum(c["lambda_comp"] for c in comp_list))
    obs_exp_sum = float(sum(c["lambda_exp"] for c in comp_list))
    reps = int(boot_sum.size)
    crit = float(np.quantile(boot_sum, 1.0 - ALPHA))
    p_val = float((1 + np.sum(boot_sum >= obs_comp_sum)) / (1 + reps))
    rejects = bool(obs_comp_sum > crit and obs_comp_sum > 0.0)
    return {
        "bootstrap_reps": reps,
        "alpha_one_sided": ALPHA,
        "n_patches_included": len(comp_list),
        "lambda_comp_sum": round(obs_comp_sum, 6),
        "lambda_exp_sum": round(obs_exp_sum, 6),
        "critical_alpha_0.01": round(crit, 6),
        "p_one_sided": round(p_val, 5),
        "rejects_mixture": rejects,
    }


# ------------------------------------------------------------------ Plan Step 7: Control power & Section 5 failure checks

def _simulate_combined_cv2_power(patch_specs, crit_comb, per_patch_rates, reps=POWER_REPS,
                                 seed_tag=700, scale_matched=False, null_cache=None):
    """Simulate combined CV^2 power at alpha=0.01 for given patch sizes, nulls, and two-step rates."""
    if not patch_specs or crit_comb is None:
        return None
    ns = np.array([p["n"] for p in patch_specs], dtype=float)
    denom = math.sqrt(float(ns.sum()))
    z_sum = np.zeros(reps, dtype=float)
    rng = np.random.default_rng([SEED, seed_tag, len(patch_specs), int(ns.sum())])
    used_rates = []
    for spec, (k1_in, k2_in) in zip(patch_specs, per_patch_rates):
        n_i = int(spec["n"])
        q_i = float(spec["quantum"])
        mu_i = float(spec["mean_duration"])
        k_lo, k_hi = min(float(k1_in), float(k2_in)), max(float(k1_in), float(k2_in))
        rr = k_hi / max(k_lo, 1e-12)
        if scale_matched:
            k1 = (1.0 + rr) / (rr * max(mu_i, 1e-12))
            k2 = rr * k1
            mu0_i = float(spec["null_mean"])
            sd_i = float(spec["null_sd"])
        else:
            k1, k2 = k_lo, k_hi
            mean_model = 1.0 / k1 + 1.0 / k2
            if q_i > 0.0 and mean_model < 2.0 * q_i:
                return None
            if q_i > 0.0 and abs(q_i / mean_model - q_i / max(mu_i, 1e-12)) > 0.02 and null_cache is not None:
                d_null = null_cache.draws(
                    n_i,
                    q_i / mean_model,
                    dead_time=spec.get("dead_time", 0.0),
                    filter_cutoff=spec.get("filter_cutoff"),
                    burst_criterion=spec.get("burst_criterion"),
                    mean_duration=mean_model,
                )
                mu0_i = float(d_null.mean())
                sd_i = float(d_null.std())
            else:
                mu0_i = float(spec["null_mean"])
                sd_i = float(spec["null_sd"])
        used_rates.append([round(k1, 6), round(k2, 6)])
        sim = rng.exponential(1.0 / k1, size=(reps, n_i)) + rng.exponential(1.0 / k2, size=(reps, n_i))
        if q_i > 0.0:
            sim = np.floor(sim / q_i) * q_i + q_i / 2.0
        z_i = (cv2(sim) - mu0_i) / sd_i if sd_i > 0.0 else np.zeros(reps, dtype=float)
        z_sum += math.sqrt(float(n_i)) * z_i
    comb_sim = z_sum / denom
    return {
        "power_alpha_0.01": round(float(np.mean(comb_sim < crit_comb)), 4),
        "simulation_reps": int(reps),
        "scale_matched_to_patch_mean": bool(scale_matched),
        "rate_ratio": round(float(max(per_patch_rates[0]) / max(min(per_patch_rates[0]), 1e-12)), 4),
        "simulated_rates_per_s": used_rates,
    }


# ------------------------------------------------------------------ Metadata & analysis helpers

def _parse_bool_flag(raw):
    if isinstance(raw, bool):
        return raw
    if isinstance(raw, (int, float)):
        return bool(raw)
    s = str(raw).strip().lower()
    if s in {"1", "true", "yes", "y", "t"}:
        return True
    if s in {"", "0", "false", "no", "n", "f"}:
        return False
    raise ValueError(f"unrecognized boolean flag: {raw!r}")


def _parse_optional_int(raw):
    if raw is None:
        return None
    if isinstance(raw, bool):
        raise ValueError(f"invalid integer value: {raw!r}")
    if isinstance(raw, int):
        return raw
    s = str(raw).strip()
    if s == "":
        return None
    f = float(s)
    if not math.isfinite(f) or not f.is_integer():
        raise ValueError(f"non-integer channel count: {raw!r}")
    return int(f)


def _parse_optional_float(raw):
    if raw is None:
        return None
    s = str(raw).strip()
    if s == "":
        return None
    v = float(s)
    if not math.isfinite(v):
        raise ValueError(f"non-finite float metadata value: {raw!r}")
    return v


def _resolve_patch_meta(construct, pid, metadata):
    if not metadata:
        return {}
    entry = {}
    if pid in metadata and isinstance(metadata[pid], dict):
        entry.update(metadata[pid])
    if construct in metadata and isinstance(metadata[construct], dict):
        sub = metadata[construct]
        if pid in sub and isinstance(sub[pid], dict):
            entry.update(sub[pid])
    if not entry:
        return {}
    out = {}
    n_ch_raw = entry.get("n_channels", entry.get("active_channels"))
    if n_ch_raw is not None:
        out["n_channels"] = _parse_optional_int(n_ch_raw)
    if "unstable" in entry and entry["unstable"] is not None:
        out["unstable"] = _parse_bool_flag(entry["unstable"])
    reason_raw = entry.get("reason", entry.get("exclude_reason"))
    if reason_raw is not None and str(reason_raw).strip():
        out["reason"] = str(reason_raw).strip()
    for fld in ("dead_time", "filter_cutoff", "burst_criterion", "quantum"):
        if fld in entry and entry[fld] is not None:
            out[fld] = _parse_optional_float(entry[fld])
    return out


def _resolve_patch_settings(meta, global_dt, global_q, global_fc, global_tc, dt_from_fc):
    p_fc = meta.get("filter_cutoff", global_fc)
    p_tc = meta.get("burst_criterion", global_tc)
    p_q = meta.get("quantum", global_q)
    if p_fc is not None and p_fc <= 0:
        raise ValueError("filter_cutoff must be positive")
    if p_tc is not None and p_tc < 0:
        raise ValueError("burst_criterion must be non-negative")
    if p_q is None or p_q < 0:
        raise ValueError("quantum must be non-negative")
    if "dead_time" in meta and meta["dead_time"] is not None:
        p_dt = float(meta["dead_time"])
    elif "filter_cutoff" in meta and meta["filter_cutoff"] is not None and (global_dt is None or dt_from_fc):
        p_dt = 0.36 / float(p_fc)
    elif global_dt is not None:
        p_dt = float(global_dt)
    elif p_fc is not None:
        p_dt = 0.36 / float(p_fc)
    else:
        raise ValueError("dead time must be provided, or derived via filter_cutoff (t_res = 0.36 / f_c)")
    if p_dt < 0:
        raise ValueError("dead time must be non-negative")
    return float(p_dt), float(p_q), (float(p_fc) if p_fc is not None else None), (float(p_tc) if p_tc is not None else None)


def analyse(patches, dead_time, quantum, null, min_bursts, only=None, filter_cutoff=None, burst_criterion=None,
            metadata=None, construct=None, dt_from_fc=False, resolution_scale=1.0,
            compare=False, bootstrap_reps=BOOTSTRAP_REPS, precomputed_comp=None):
    rows, stats, patch_specs, comp_stats = [], {}, {}, {}
    verified_flags = []
    for pid in sorted(patches):
        raw = np.asarray(patches[pid], dtype=float)
        meta = _resolve_patch_meta(construct, pid, metadata)
        n_ch = meta.get("n_channels")
        unstable = bool(meta.get("unstable", False))

        # Plan Step 1: check one-active-channel & stability BEFORE statistical analysis
        if n_ch is not None and n_ch != 1:
            rows.append({
                "patch_id": pid,
                "n_channels": int(n_ch),
                "n_bursts_raw": int(raw.size),
                "n_bursts_used": 0,
                "included": False,
                "excluded_before_analysis": True,
                "reason": f"n_channels != 1 (N={n_ch})",
            })
            continue
        if unstable:
            reason = meta.get("reason") or "flagged as unstable"
            rows.append({
                "patch_id": pid,
                "n_channels": (int(n_ch) if n_ch is not None else None),
                "unstable": True,
                "n_bursts_raw": int(raw.size),
                "n_bursts_used": 0,
                "included": False,
                "excluded_before_analysis": True,
                "reason": reason,
            })
            continue

        p_dt_1x, p_q, p_fc_1x, p_tc = _resolve_patch_settings(
            meta, dead_time, quantum, filter_cutoff, burst_criterion, dt_from_fc
        )
        p_dt = p_dt_1x * resolution_scale
        p_fc = (p_fc_1x / resolution_scale) if p_fc_1x is not None else None

        x = adjust(raw, p_dt, p_q)
        n = int(x.size)
        if n < min_bursts or (only is not None and pid not in only):
            rows.append({
                "patch_id": pid,
                "n_channels": (int(n_ch) if n_ch is not None else None),
                "n_bursts_used": n,
                "included": False,
                "reason": f"fewer than {min_bursts} bursts after the dead-time rule" if n < min_bursts else "not included at every resolution",
            })
            continue

        verified_flags.append(n_ch == 1)
        s = patch_stats(x, p_q, null, dead_time=p_dt, filter_cutoff=p_fc, burst_criterion=p_tc)
        stats[pid] = s
        patch_specs[pid] = {
            "n": n,
            "quantum": p_q,
            "mean_duration": float(x.mean()),
            "dead_time": p_dt,
            "filter_cutoff": p_fc,
            "burst_criterion": p_tc,
            "null_mean": s["null_mean"],
            "null_sd": s["null_sd"],
        }
        kept = raw[raw >= p_dt]
        gap = float(kept.min() - p_dt) if kept.size else float("nan")
        expected = float(x.mean() / n)
        row = {
            "patch_id": pid,
            "n_channels": (int(n_ch) if n_ch is not None else None),
            "n_bursts_used": n,
            "included": True,
            "resolution_s": round(p_dt, 9),
            "frame_quantum_s": round(p_q, 9),
            "cv2": round(s["cv2"], 4),
            "critical_cv2_alpha_0.01": round(s["crit"], 4),
            "p_one_sided": round(s["p"], 5),
            "rejects_mixture": bool(s["cv2"] < s["crit"]),
            "shortest_burst_minus_resolution_s": (round(gap, 6) if math.isfinite(gap) else None),
            "p_gap_this_large_if_resolution_correct": (round(float(math.exp(-gap / expected)), 5) if (p_q == 0 and math.isfinite(gap) and expected > 0) else None),
            "resolution_may_be_too_small": bool(p_q == 0 and math.isfinite(gap) and gap > 3 * expected),
        }
        if compare:
            if precomputed_comp is not None and pid in precomputed_comp:
                c = precomputed_comp[pid]
            else:
                c = comparator_patch_stats(
                    x,
                    p_q,
                    bootstrap_reps=bootstrap_reps,
                    seed_parts=(_stable_id_int(construct or "all"), _stable_id_int(pid), int(round(resolution_scale * 100))),
                )
            comp_stats[pid] = c
            row["comparator"] = _format_patch_comparator(c)
        rows.append(row)

    one_ch_verified = bool(verified_flags and all(verified_flags))
    out = {
        "resolution_s": (round(dead_time * resolution_scale, 9) if dead_time is not None else None),
        "filter_cutoff_hz": ((filter_cutoff / resolution_scale) if filter_cutoff is not None else None),
        "burst_criterion_s": burst_criterion,
        "frame_quantum_s": quantum,
        "min_bursts": min_bursts,
        "one_active_channel_verified": one_ch_verified,
        "patches": rows,
        "n_patches_included": len(stats),
    }
    if not one_ch_verified:
        if verified_flags or not any(r.get("n_channels") is not None for r in rows):
            out["one_active_channel_note"] = "active-channel count not supplied in input; assumed 1 per caller"
        else:
            out["one_active_channel_note"] = "no patches passed active-channel/stability inclusion criteria"
    if stats:
        out["combined"] = combine([stats[k] for k in sorted(stats)], len(stats))
    if compare and comp_stats:
        comb_comp = combine_comparator([comp_stats[k] for k in sorted(comp_stats)])
        k1_med = float(np.median([comp_stats[k]["hypo_params"][0] for k in sorted(comp_stats)]))
        k2_med = float(np.median([comp_stats[k]["hypo_params"][1] for k in sorted(comp_stats)]))
        out["comparator"] = {
            "bootstrap_reps": int(bootstrap_reps),
            "alpha_one_sided": ALPHA,
            "n_patches_included": len(comp_stats),
            "combined": comb_comp,
            "median_fitted_hypoexponential_rates_per_s": [round(k1_med, 6), round(k2_med, 6)],
        }
    return out, set(stats), patch_specs, comp_stats


def run(table, dead_time=None, quantum=0.0, control=None, min_bursts=MIN_BURSTS, min_bursts_control=MIN_BURSTS_CONTROL,
        filter_cutoff=None, burst_criterion=None, metadata=None, compare=True, bootstrap_reps=BOOTSTRAP_REPS,
        power_reps=POWER_REPS):
    dt_from_fc = False
    if dead_time is None:
        if filter_cutoff is not None and filter_cutoff > 0:
            dead_time = 0.36 / float(filter_cutoff)
            dt_from_fc = True
        elif not metadata:
            raise ValueError("dead time must be provided, or derived via filter_cutoff (t_res = 0.36 / f_c)")
    else:
        dead_time = float(dead_time)

    quantum = float(quantum)
    filter_cutoff = float(filter_cutoff) if filter_cutoff is not None else None
    burst_criterion = float(burst_criterion) if burst_criterion is not None else None
    bootstrap_reps = int(bootstrap_reps)
    power_reps = int(power_reps)

    if (dead_time is not None and dead_time < 0) or quantum < 0 or min(min_bursts, min_bursts_control) < 2:
        raise ValueError("dead time and quantum must be non-negative, and burst minimums at least 2")
    if filter_cutoff is not None and filter_cutoff <= 0:
        raise ValueError("filter_cutoff must be positive")
    if burst_criterion is not None and burst_criterion < 0:
        raise ValueError("burst_criterion must be non-negative")
    if bootstrap_reps < 1 or power_reps < 1:
        raise ValueError("bootstrap_reps and power_reps must be at least 1")

    null = Null()
    res = {
        "plan": "specs/CFTR_DATA_REQUEST_AND_ANALYSIS_PLAN.md",
        "alpha_one_sided": ALPHA,
        "monte_carlo_draws": REPS,
        "comparator_bootstrap_draws": (bootstrap_reps if compare else 0),
        "power_simulation_draws": power_reps,
        "seed": SEED,
        "filter_cutoff_hz": filter_cutoff,
        "burst_criterion_s": burst_criterion,
        "implemented": [
            "step 1: one-active-channel verification, unstable-patch exclusion, and burst minimums (50 WT, 5 control)",
            "step 2: dead-time adjustment, 0.36/f_c default filter resolution, and frame-grid handling",
            "steps 3-5: CV^2 statistic, per-patch critical threshold, and sqrt(N) combination",
            "step 4 calibration: instrument-calibrated false-alarm threshold for dead time, filtering (f_c), and burst criterion (t_crit)",
            "step 6: composite-null GLRT comparator (2-exp mixture vs. 2-step hypoexponential) with parametric bootstrap and secondary single-exp LRT",
            "step 7: control power simulation at observed burst counts and Section 5 failure criteria",
            "section 5: 1x / 2x resolution sensitivity decision rule",
        ],
        "not_implemented": [],
        "constructs": {},
    }

    # First pass: run CV^2 (1x, 2x, both) and comparator (1x) per construct
    internal_data = {}
    for construct, patches in sorted(table.items()):
        mb = min_bursts_control if control is not None and construct == control else min_bursts
        _, inc1, _, _ = analyse(
            patches, dead_time, quantum, null, mb,
            filter_cutoff=filter_cutoff, burst_criterion=burst_criterion,
            metadata=metadata, construct=construct, dt_from_fc=dt_from_fc, resolution_scale=1.0,
            compare=False,
        )
        _, inc2, _, _ = analyse(
            patches, dead_time, quantum, null, mb,
            filter_cutoff=filter_cutoff, burst_criterion=burst_criterion,
            metadata=metadata, construct=construct, dt_from_fc=dt_from_fc, resolution_scale=2.0,
            compare=False,
        )
        both = inc1 & inc2
        main, _, specs_1x, comp_1x = analyse(
            patches, dead_time, quantum, null, mb,
            filter_cutoff=filter_cutoff, burst_criterion=burst_criterion,
            metadata=metadata, construct=construct, dt_from_fc=dt_from_fc, resolution_scale=1.0,
            compare=compare, bootstrap_reps=bootstrap_reps,
        )
        role = "control" if control is not None and construct == control else "test"
        main["role"] = role
        a, _, specs_a, comp_a = analyse(
            patches, dead_time, quantum, null, mb, only=both,
            filter_cutoff=filter_cutoff, burst_criterion=burst_criterion,
            metadata=metadata, construct=construct, dt_from_fc=dt_from_fc, resolution_scale=1.0,
            compare=compare, bootstrap_reps=bootstrap_reps, precomputed_comp=comp_1x,
        )
        b, _, _, _ = analyse(
            patches, dead_time, quantum, null, mb, only=both,
            filter_cutoff=filter_cutoff, burst_criterion=burst_criterion,
            metadata=metadata, construct=construct, dt_from_fc=dt_from_fc, resolution_scale=2.0,
            compare=False,
        )
        r1 = a.get("combined", {}).get("rejects_mixture")
        r2 = b.get("combined", {}).get("rejects_mixture")
        sens = {
            "patches_included_at_both_resolutions": len(both),
            "combined_at_1x": a.get("combined"),
            "combined_at_2x": b.get("combined"),
            "rejects_at_1x_not_2x": bool(r1 and r2 is False),
            "rejects_at_2x_not_1x": bool(r2 and r1 is False),
        }
        if compare and "comparator" in a:
            sens["comparator_at_1x"] = {
                "patches_included_at_both_resolutions": len(both),
                "bootstrap_reps": int(bootstrap_reps),
                "combined": a["comparator"]["combined"],
                "median_fitted_hypoexponential_rates_per_s": a["comparator"]["median_fitted_hypoexponential_rates_per_s"],
            }
        main["resolution_sensitivity"] = sens
        main["decision"] = (
            "no_patches" if r1 is None or r2 is None
            else "rejects_mixture" if (r1 and r2)
            else "resolution_sensitive_inconclusive" if (r1 or r2)
            else "does_not_reject"
        )
        res["constructs"][construct] = main
        internal_data[construct] = {
            "role": role,
            "specs": (specs_a if specs_a else specs_1x),
            "combined_for_power": (a.get("combined") if specs_a else main.get("combined")),
            "comp": (comp_a if comp_a else comp_1x),
        }

    # Collect median/pooled test construct fitted two-step rates (when available)
    test_k1s, test_k2s = [], []
    test_rate_source = []
    for cname in sorted(internal_data):
        info = internal_data[cname]
        if info["role"] == "test" and info["comp"]:
            test_rate_source.append(cname)
            for pid in sorted(info["comp"]):
                k1_p, k2_p = info["comp"][pid]["hypo_params"]
                test_k1s.append(k1_p)
                test_k2s.append(k2_p)
    med_test_rates = (float(np.median(test_k1s)), float(np.median(test_k2s))) if test_k1s else None

    # Plan Step 7: Control power simulation & Section 5 failure checks
    any_sec5_fail = False
    sec5_reasons = []
    has_control = any(info["role"] == "control" and info["specs"] for info in internal_data.values())
    res["controlled"] = bool(has_control)

    for cname in sorted(res["constructs"]):
        main = res["constructs"][cname]
        info = internal_data[cname]
        specs_map = info["specs"]
        specs_list = [specs_map[pid] for pid in sorted(specs_map)]
        comb_ref = info["combined_for_power"]
        crit_comb = comb_ref.get("critical_alpha_0.01") if comb_ref else None
        fail_flags = []

        if info["role"] == "control":
            if specs_list and crit_comb is not None:
                obs_counts = [s["n"] for s in specs_list]
                var_a_pairs = [VARIANT_A_RATES for _ in specs_list]
                pwr_var_a_scaled = _simulate_combined_cv2_power(
                    specs_list, crit_comb, var_a_pairs, reps=power_reps, seed_tag=701, scale_matched=True, null_cache=null
                )
                pwr_var_a_nominal = _simulate_combined_cv2_power(
                    specs_list, crit_comb, var_a_pairs, reps=power_reps, seed_tag=702, scale_matched=False, null_cache=null
                )
                ctrl_power = {
                    "observed_patch_burst_counts": obs_counts,
                    "total_bursts": int(sum(obs_counts)),
                    "variant_A_reference_rates_per_s": list(VARIANT_A_RATES),
                    "power_at_variant_A_scale_matched": pwr_var_a_scaled,
                    "power_at_variant_A_nominal_rates": pwr_var_a_nominal,
                }
                if med_test_rates is not None:
                    test_pairs = [med_test_rates for _ in specs_list]
                    ctrl_power["test_construct_source"] = test_rate_source
                    ctrl_power["test_construct_median_rates_per_s"] = [round(med_test_rates[0], 6), round(med_test_rates[1], 6)]
                    ctrl_power["power_at_test_rates_scale_matched"] = _simulate_combined_cv2_power(
                        specs_list, crit_comb, test_pairs, reps=power_reps, seed_tag=703, scale_matched=True, null_cache=null
                    )
                    ctrl_power["power_at_test_rates_nominal"] = _simulate_combined_cv2_power(
                        specs_list, crit_comb, test_pairs, reps=power_reps, seed_tag=704, scale_matched=False, null_cache=null
                    )
                main["control_power_simulation"] = ctrl_power
            if main["decision"] == "rejects_mixture":
                fail_flags.append("combined CV^2 test rejected positive exponential mixture on control construct at alpha=0.01")
        else:
            # Test construct Section 5 failure criterion:
            # comparator rejects at alpha=0.01 AND combined CV^2 does not, while CV^2 power at comparator rates >= 0.80
            comp_block = main["resolution_sensitivity"].get("comparator_at_1x") or main.get("comparator")
            comp_rej = bool(comp_block and comp_block.get("combined", {}).get("rejects_mixture"))
            cv2_rej = (main["decision"] == "rejects_mixture")
            if comp_rej and not cv2_rej and specs_list and crit_comb is not None and info["comp"]:
                comp_pairs = [info["comp"][pid]["hypo_params"] for pid in sorted(specs_map)]
                pwr_comp = _simulate_combined_cv2_power(
                    specs_list, crit_comb, comp_pairs, reps=power_reps, seed_tag=705, scale_matched=False, null_cache=null
                )
                if pwr_comp is None:
                    pwr_comp = _simulate_combined_cv2_power(
                        specs_list, crit_comb, comp_pairs, reps=power_reps, seed_tag=706, scale_matched=True, null_cache=null
                    )
                main["test_power_at_comparator_rates"] = pwr_comp
                if pwr_comp and pwr_comp["power_alpha_0.01"] >= 0.80:
                    fail_flags.append(
                        f"comparator rejects at alpha=0.01 while combined CV^2 does not (decision={main['decision']}), "
                        f"with simulated CV^2 power {pwr_comp['power_alpha_0.01']:.4f} >= 0.80 at comparator fitted rates"
                    )

        main["section_5_failure_triggered"] = bool(fail_flags)
        main["section_5_failure_reasons"] = fail_flags
        if fail_flags:
            any_sec5_fail = True
            sec5_reasons.extend([f"{cname}: {msg}" for msg in fail_flags])

        # Plan Section 4 outcome summary
        comp_block = main["resolution_sensitivity"].get("comparator_at_1x") or main.get("comparator")
        comp_rej = bool(comp_block and comp_block.get("combined", {}).get("rejects_mixture"))
        if main["section_5_failure_triggered"]:
            main["outcome_reading"] = "check fails Section 5 criterion"
        elif info["role"] == "control":
            main["outcome_reading"] = "control does not reject mixture" if main["decision"] == "does_not_reject" else main["decision"]
        elif main["decision"] == "rejects_mixture" and (not compare or comp_rej):
            main["outcome_reading"] = "check agrees with likelihood analysis (both reject)"
        elif main["decision"] == "does_not_reject" and (not compare or not comp_rej):
            main["outcome_reading"] = "check agrees with likelihood analysis (neither rejects)"
        else:
            main["outcome_reading"] = "check needs more data or shows resolution sensitivity"

    res["section_5_failure_triggered"] = bool(any_sec5_fail)
    res["section_5_failure_reasons"] = sec5_reasons
    return res


def read_csv(path, return_metadata=False):
    table = defaultdict(lambda: defaultdict(list))
    meta_out = defaultdict(dict)
    dropped = 0
    with open(path, newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        for row in reader:
            pid = (row.get("patch_id") or "").strip()
            construct = (row.get("construct") or "all").strip()
            if pid:
                pmeta = meta_out[construct].setdefault(pid, {})
                n_ch_raw = row.get("n_channels")
                if n_ch_raw is None or str(n_ch_raw).strip() == "":
                    n_ch_raw = row.get("active_channels")
                n_ch = _parse_optional_int(n_ch_raw)
                if n_ch is not None:
                    if "n_channels" in pmeta and pmeta["n_channels"] != n_ch:
                        raise ValueError(f"conflicting n_channels values for patch {pid!r}")
                    pmeta["n_channels"] = n_ch
                if row.get("unstable") is not None and str(row.get("unstable")).strip() != "":
                    if _parse_bool_flag(row.get("unstable")):
                        pmeta["unstable"] = True
                reason_val = (row.get("exclude_reason") or row.get("reason") or "").strip()
                if reason_val:
                    if "reason" in pmeta:
                        existing_reasons = [part.strip() for part in pmeta["reason"].split("; ") if part.strip()]
                        if reason_val not in existing_reasons:
                            pmeta["reason"] = f"{pmeta['reason']}; {reason_val}"
                    else:
                        pmeta["reason"] = reason_val
            try:
                v = float(row["duration"])
            except (TypeError, ValueError, KeyError):
                dropped += 1
                continue
            if not math.isfinite(v) or v <= 0 or not pid:
                dropped += 1
                continue
            table[construct][pid].append(v)
    if return_metadata:
        return table, dropped, dict(meta_out)
    return table, dropped


def _merge_metadata(csv_meta, json_meta):
    merged = {}
    if csv_meta:
        for k, v in csv_meta.items():
            merged[k] = {pid: dict(pm) for pid, pm in v.items()} if isinstance(v, dict) else dict(v)
    if json_meta:
        src = json_meta.get("patches", json_meta) if isinstance(json_meta, dict) and isinstance(json_meta.get("patches"), dict) and len(json_meta) == 1 else json_meta
        for key, val in src.items():
            if not isinstance(val, dict):
                continue
            # Support either patch_id -> {...} or construct -> {patch_id -> {...}}
            if all(isinstance(sv, dict) for sv in val.values()) and val:
                c_map = merged.setdefault(key, {})
                for pid, pm in val.items():
                    c_map.setdefault(pid, {}).update(pm)
            else:
                # Apply patch_id metadata to any matching construct entry and at top level
                merged.setdefault(key, {}).update(val)
                for c_val in merged.values():
                    if isinstance(c_val, dict) and key in c_val and isinstance(c_val[key], dict):
                        c_val[key].update(val)
    return merged


def self_test():
    rng = np.random.default_rng(7)
    td = 0.002
    sizes = (120, 250, 80, 400, 150, 60)

    def patches(draw):
        return {f"p{i}": draw(n) + td for i, n in enumerate(sizes)}

    table = {
        "two_step": patches(lambda n: rng.exponential(1 / 12, n) + rng.exponential(1 / 8, n)),
        "single_exponential": patches(lambda n: rng.exponential(1 / 5, n)),
        "mixture": patches(lambda n: np.where(rng.random(n) < 0.5, rng.exponential(1 / 12, n), rng.exponential(1 / 3, n))),
        "control": {f"c{i}": rng.exponential(2.0, n) + td for i, n in enumerate((6, 9, 14))},
    }
    res = run(table, td, 0.0, control="control", bootstrap_reps=200)
    got = {k: v["decision"] for k, v in res["constructs"].items()}
    ok = got == {"two_step": "rejects_mixture", "single_exponential": "does_not_reject", "mixture": "does_not_reject", "control": "does_not_reject"}
    ok = ok and res["constructs"]["control"]["n_patches_included"] == 3 and res["constructs"]["two_step"]["n_patches_included"] == len(sizes)
    ok = ok and res["not_implemented"] == []
    ok = ok and res["constructs"]["two_step"]["resolution_sensitivity"]["comparator_at_1x"]["combined"]["rejects_mixture"] is True
    ok = ok and res["constructs"]["mixture"]["resolution_sensitivity"]["comparator_at_1x"]["combined"]["rejects_mixture"] is False
    ok = ok and res["constructs"]["control"]["control_power_simulation"]["power_at_variant_A_scale_matched"]["power_alpha_0.01"] < 0.50
    json.dumps(res, allow_nan=False)

    # Test 0.36 / f_c default derivation
    res_fc = run(table, filter_cutoff=1000, control="control", bootstrap_reps=200)
    ok = ok and math.isclose(res_fc["constructs"]["two_step"]["resolution_s"], 0.00036, rel_tol=1e-5)

    null, trials, rej, rej_q, rej_n = Null(11), 200, 0, 0, 0
    q, td_q = 0.1, 0.17                      # frame grid, off-grid true dead time: events below it are never observed
    for _ in range(trials):
        ps = {f"p{i}": rng.exponential(1.0, n) for i, n in enumerate(sizes)}
        st = [patch_stats(v, 0.0, null) for v in ps.values()]
        rej += combine(st, len(st))["rejects_mixture"]
        t_true = rng.exponential(1.0, 700)
        t_true = t_true[t_true >= td_q][:300]
        for kind in ("up", "nearest"):
            y = (np.ceil(t_true / q) if kind == "up" else np.round(t_true / q)) * q
            s = patch_stats(adjust(y, td_q, q), q, null)
            if kind == "up":
                rej_q += s["cv2"] < s["crit"]
            else:
                rej_n += s["cv2"] < s["crit"]
    rate, rate_q, rate_n = rej / trials, rej_q / trials, rej_n / trials
    ok = ok and rate <= 0.02 and rate_q <= 0.02 and rate_n <= 0.02

    # Verify that instrument-calibrated threshold prevents false alarm on filtered bursts with burst criterion
    res_cal = run(table, dead_time=0.002, filter_cutoff=1000, burst_criterion=0.005, control="control", bootstrap_reps=200)
    got_cal = {k: v["decision"] for k, v in res_cal["constructs"].items()}
    ok = ok and got_cal["single_exponential"] == "does_not_reject" and got_cal["mixture"] == "does_not_reject"

    for bad in (dict(dead_time=-0.1, quantum=0.0), dict(dead_time=0.0, quantum=-1.0), dict(filter_cutoff=-10.0)):
        try:
            run(table, **bad)
            ok = False
        except ValueError:
            pass
    try:
        adjust(np.array([0.11, 0.23]), 0.0, 0.1)
        ok = False
    except ValueError:
        pass
    print("self-test decisions", got, "| false-alarm rate, combined:", rate, "| frame grid, missed short events, rounded up / nearest:", rate_q, rate_n,
          f"({trials} datasets each)", "| OK" if ok else "| FAIL")
    return 0 if ok else 1


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("csv", nargs="?")
    ap.add_argument("--dead-time", type=float, default=None, help="resolution imposed on the record, seconds; required unless --filter-cutoff is given")
    ap.add_argument("--filter-cutoff", type=float, default=None, help="filter -3 dB cutoff frequency f_c in Hz (sets default dead_time = 0.36 / f_c)")
    ap.add_argument("--burst-criterion", type=float, default=None, help="critical closed time t_crit in seconds defining burst boundaries")
    ap.add_argument("--quantum", type=float, default=0.0, help="frame interval in seconds, if durations lie on a grid")
    ap.add_argument("--control", default=None, help="name of the control construct (uses --min-bursts-control)")
    ap.add_argument("--min-bursts", type=int, default=MIN_BURSTS)
    ap.add_argument("--min-bursts-control", type=int, default=MIN_BURSTS_CONTROL)
    ap.add_argument("--metadata", default=None, help="optional JSON file mapping patch_id -> {n_channels, unstable, reason, dead_time, filter_cutoff, burst_criterion}")
    ap.add_argument("--bootstrap-reps", type=int, default=BOOTSTRAP_REPS, help="parametric bootstrap replicates for Step 6 GLRT comparator")
    ap.add_argument("--no-compare", action="store_true", help="skip the Step 6 GLRT comparator")
    ap.add_argument("--out")
    ap.add_argument("--self-test", action="store_true")
    a = ap.parse_args(argv)
    if a.self_test:
        return self_test()
    if not a.csv or (a.dead_time is None and a.filter_cutoff is None and not a.metadata):
        ap.error("give a CSV and --dead-time (or --filter-cutoff or --metadata), or --self-test")
    table, dropped, csv_meta = read_csv(a.csv, return_metadata=True)
    json_meta = None
    if a.metadata:
        with open(a.metadata, encoding="utf-8") as mf:
            json_meta = json.load(mf)
        if not isinstance(json_meta, dict):
            raise ValueError("--metadata JSON must be an object mapping patch_id to metadata")
    merged_meta = _merge_metadata(csv_meta, json_meta)
    res = run(
        table,
        a.dead_time,
        a.quantum,
        a.control,
        a.min_bursts,
        a.min_bursts_control,
        filter_cutoff=a.filter_cutoff,
        burst_criterion=a.burst_criterion,
        metadata=(merged_meta or None),
        compare=(not a.no_compare),
        bootstrap_reps=a.bootstrap_reps,
    )
    res["rows_dropped_as_missing_or_non_positive"] = dropped
    text = json.dumps(res, indent=1, allow_nan=False)
    if a.out:
        with open(a.out, "w", encoding="utf-8") as of:
            of.write(text)
    else:
        print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
