#!/usr/bin/env python3
"""Whole-sample variability test for event lists (burst or dwell durations).

Claim tested: "these durations come from a positive mixture of exponentials" (what any reversible gating scheme gives).
Fact used:    every positive mixture of exponentials has squared coefficient of variation CV^2 >= 1, with equality only
              for a single exponential (law of total variance). A sequential multi-step dwell has CV^2 < 1.
Test:         one-sided, reject when CV^2 is below a critical value calibrated by Monte Carlo at the boundary null
              (single exponential) for the same sample size and the same dead-time / frame-quantization handling.
              Like the curvature test it is one-way: CV^2 < 1 rules the mixture out; CV^2 >= 1 shows nothing.
The statistic itself is classical (the "randomness parameter" of single-molecule kinetics). What is done here is the
calibration, the dead-time handling, the per-patch combination, and the comparison with the package's 3-bin test.

Part A is simulation only (DERIVED, NOT MEASURED). Part B uses the cached F1-ATPase dwell times (Zenodo 5567245).
No CFTR recording is analysed. Usage: python benchmark/whole_sample/cv2_burst_test.py
"""
import hashlib, json, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BENCH = os.path.dirname(HERE)
OUT = os.path.join(HERE, "cv2_burst_test.json")
SEED = 20261113
REPS = 20000
rng = np.random.default_rng(SEED)


def cv2(x, axis=-1):
    m = x.mean(axis=axis)
    return x.var(axis=axis, ddof=1) / m**2


def hypo(a, b, size):
    """Two sequential steps with rates a and b (per second)."""
    return rng.exponential(1 / a, size) + rng.exponential(1 / b, size)


def crit(n, alpha, reps=REPS):
    return float(np.quantile(cv2(rng.exponential(1.0, (reps, n))), alpha))


def r(v, k=4):
    return round(float(v), k)


out = {"seed": SEED, "reps": REPS, "label": "Part A: DERIVED, NOT MEASURED (simulation). Part B: public F1-ATPase dwell times."}

# ------------------------------------------------------------------ Part A: simulation
variants = {"A_12_8_per_s": (12.0, 8.0), "B_5.988_4_per_s": (5.988, 4.0)}
out["population_cv2"] = {k: r((1 / a**2 + 1 / b**2) / (1 / a + 1 / b) ** 2) for k, (a, b) in variants.items()}
power = {}
for n in (30, 50, 75, 100, 200, 300, 1000):
    row = {}
    for alpha in (0.05, 0.01):
        c = crit(n, alpha)
        row[f"critical_cv2_alpha_{alpha}"] = r(c)
        for k, (a, b) in variants.items():
            row[f"power_{k}_alpha_{alpha}"] = r((cv2(hypo(a, b, (REPS, n))) < c).mean())
        mix = np.where(rng.random((REPS, n)) < 0.5, rng.exponential(1 / 12, (REPS, n)), rng.exponential(1 / 3, (REPS, n)))
        row[f"false_alarm_two_rate_mixture_alpha_{alpha}"] = r((cv2(mix) < c).mean())
        row[f"false_alarm_single_exponential_alpha_{alpha}"] = r((cv2(rng.exponential(1.0, (REPS, n))) < c).mean())
    power[f"N_{n}"] = row
out["partA_power"] = power

# power depends on how similar the two step rates are: one dominant slow step looks almost exponential (CV^2 -> 1)
ratio = {}
for rr in (1.0, 1.5, 3.0, 5.0, 10.0):
    slow = 4.0
    fast = slow * rr
    row = {"rates_per_s": [fast, slow], "population_cv2": r((1 / fast**2 + 1 / slow**2) / (1 / fast + 1 / slow) ** 2)}
    for n_ in (100, 300, 1000):
        row[f"power_N_{n_}_alpha_0.01"] = r((cv2(hypo(fast, slow, (REPS, n_))) < crit(n_, 0.01)).mean())
    ratio[f"rate_ratio_{rr:g}"] = row
out["partA_power_vs_rate_ratio"] = dict(note="two sequential steps; ratio = faster rate / slower rate; one-sided alpha 0.01", **ratio)

# level near the boundary of the null: two-component mixtures close to a single exponential
near = {}
for w_, ratio_ in ((0.5, 1.2), (0.5, 2.0), (0.1, 3.0), (0.9, 3.0), (0.02, 10.0)):
    row = {"weight_of_fast_component": w_, "rate_ratio": ratio_}
    for n_ in (100, 300):
        pick = rng.random((REPS, n_)) < w_
        x_ = np.where(pick, rng.exponential(1 / ratio_, (REPS, n_)), rng.exponential(1.0, (REPS, n_)))
        row[f"false_alarm_N_{n_}_alpha_0.01"] = r((cv2(x_) < crit(n_, 0.01)).mean())
    near[f"mixture_w{w_:g}_ratio{ratio_:g}"] = row
out["partA_false_alarm_near_boundary_mixtures"] = dict(
    note="the single exponential appears least favourable; this is shown by simulation, not proven for finite N", **near)

dead = {}
n, c = 300, crit(300, 0.01)
a, b = variants["A_12_8_per_s"]
for td in (0.005, 0.020, 0.050):
    pw, fa = [], []
    for _ in range(4000):
        x = hypo(a, b, 3 * n)
        x = x[x > td][:n] - td
        pw.append(cv2(x) < c)
        e = rng.exponential(1 / 6.0, 3 * n)
        e = e[e > td][:n] - td
        fa.append(cv2(e) < c)
    dead[f"dead_time_{int(td * 1000)}_ms"] = {"power_variant_A_alpha_0.01": r(np.mean(pw)), "false_alarm_single_exponential_alpha_0.01": r(np.mean(fa)), "runs": 4000}
out["partA_dead_time_N300"] = dict(rule="discard events shorter than the dead time and subtract it; an exponential stays exponential", **dead)

P, m = 10, 100
null_patch = cv2(rng.exponential(1.0, (REPS, m)))
mu0, sd0 = null_patch.mean(), null_patch.std()
c_pool = crit(P * m, 0.01)
z_null = ((cv2(rng.exponential(1.0, (REPS, P, m))) - mu0) / sd0).sum(axis=1) / np.sqrt(P)
c_comb = float(np.quantile(z_null, 0.01))
pool, comb = [], []
for _ in range(4000):
    scale = np.exp(rng.normal(0, 0.5, P))[:, None]
    x = hypo(a, b, (P, m)) * scale
    pool.append(cv2(x.ravel()) < c_pool)
    comb.append(((cv2(x) - mu0) / sd0).sum() / np.sqrt(P) < c_comb)
out["partA_heterogeneous_patches"] = {
    "design": "10 patches x 100 bursts, variant A, each patch's time scale drawn log-normally (sigma 0.5)",
    "power_pooling_raw_durations_alpha_0.01": r(np.mean(pool)),
    "power_per_patch_scores_combined_alpha_0.01": r(np.mean(comb)), "runs": 4000,
    "why": "pooling patches of different speed is itself a mixture and pushes CV^2 up toward 1"}

ref = json.load(open(os.path.join(BENCH, "out", "single_channel_sample_size.json"), encoding="utf-8"))
cal = {row["N"]: row for row in ref["calibration"]["variant_A_binned_vs_hypoexponential_LRT"]}
out["partA_reference_3bin_test"] = {
    "source": "benchmark/out/single_channel_sample_size.json calibration.variant_A_binned_vs_hypoexponential_LRT",
    "power_3bin_D2_alpha_0.01": {f"N_{k}": cal[k]["d2_power_z_gt_2p326"] for k in (300, 1000)},
    "power_likelihood_ratio_test_alpha_0.01": {f"N_{k}": cal[k]["lrt_power_at_null99"] for k in (300, 1000)}}

# ------------------------------------------------------------------ Part B: real dwell times (F1-ATPase, Zenodo 5567245)
npz_path = os.path.join(BENCH, "out", "zenodo5567245_f1_atpase.npz")
z = np.load(npz_path, allow_pickle=False)
conditions = {"dwell_2d": "Fig2D wild type, 3 mM ATP (saturating; hydrolysis-limited dwell)",
              "dwell_3a": "Fig3A wild type, 1 mM ATP-gamma-S (slowly hydrolysed analog)",
              "dwell_2c": "Fig2C wild type, 300 nM ATP (binding-limited dwell)"}
real = {}
for key, label in conditions.items():
    x = np.asarray(z[key], dtype=float)
    n = x.size
    u = np.unique(x)
    q = float(np.min(np.diff(u)))          # frame quantum, in the units of the cached array
    xmin = float(x.min())
    raw = float(cv2(x))
    adj = x - xmin + q / 2                   # dead-time handling: time since the shortest recordable dwell
    stat = float(cv2(adj))
    null_raw = cv2(rng.exponential(1.0, (REPS, n)))
    mu = adj.mean()
    sim = np.floor(rng.exponential(mu, (REPS, n)) / q) * q + q / 2   # same quantization as the data
    null_adj = cv2(sim)
    real[key] = {"condition": label, "n_dwells": int(n), "frame_quantum": r(q, 6), "shortest_dwell": r(xmin, 6),
                 "cv2_raw": r(raw), "p_one_sided_raw_vs_exponential": r((null_raw <= raw).mean(), 5),
                 "cv2_dead_time_adjusted": r(stat), "p_one_sided_adjusted_vs_quantized_exponential": r((null_adj <= stat).mean(), 5),
                 "p_upper_tail_adjusted": r((null_adj >= stat).mean(), 5), "critical_cv2_alpha_0.01_adjusted": r(np.quantile(null_adj, 0.01)), "reject_mixture_at_0.01": bool(stat < np.quantile(null_adj, 0.01))}
out["partB_f1_atpase"] = {
    "cache": "benchmark/out/zenodo5567245_f1_atpase.npz", "cache_sha256": hashlib.sha256(open(npz_path, "rb").read()).hexdigest(),
    "note": "rotary F1-ATPase 120-degree step dwells, not CFTR; a step-finding algorithm that misses short dwells lowers CV^2, which the adjusted statistic only partly addresses",
    "conditions": real}
out["not_run"] = {"BK_channel_dwells": "the package cache holds only binned dwell histograms truncated at 400 bins; truncation lowers CV^2 even for a mixture, so the test is not valid on that cache and was not run",
                  "CFTR": "no public event-level CFTR recording was found; nothing is analysed"}

json.dump(out, open(OUT, "w", encoding="utf-8"), indent=1)
for k, v in power.items():
    print(k, {kk.replace("power_", "").replace("_alpha_0.01", ""): vv for kk, vv in v.items() if kk.endswith("alpha_0.01")})
for k, v in ratio.items():
    print(k, v)
for k, v in near.items():
    print(k, v)
print("dead time", {k: v for k, v in dead.items()})
print("patches", out["partA_heterogeneous_patches"]["power_pooling_raw_durations_alpha_0.01"], out["partA_heterogeneous_patches"]["power_per_patch_scores_combined_alpha_0.01"])
for k, v in real.items():
    print(k, v)
print("3-bin reference", out["partA_reference_3bin_test"])
