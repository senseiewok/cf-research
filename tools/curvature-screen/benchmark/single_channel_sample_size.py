#!/usr/bin/env python3
"""Single-channel burst-duration test: reconcile sample-size claims, calibrate power, compare to a likelihood test.

Everything here is computed (nothing is hard-coded). Run:  python single_channel_sample_size.py [--quick]
Writes out/single_channel_sample_size.json.

Observable (one definition used throughout)
  Burst durations t_j (one per opening burst). Three equal-width, disjoint bins
  [0,w), [w,2w), [2w,3w). With total event count N (Poisson/multinomial; events outside the 3 bins
  still count in N), the statistic is
      D2 = ln n1 - 2 ln n2 + ln n3 .
  Delta-method variance (zero-sum contrast, exact for Poisson counts):
      N * Var(D2) = sum_i 1/p_i  *  (c_i^2)   with c = (1,-2,1)
  and 'SNR = 5' means  E[D2] / SD(D2) = 5.  SNR is NOT a detection probability; see the calibration section.
Models (survival functions, closed form)
  hypoexponential (two sequential steps with rates a, b):  S(t) = (a e^{-bt} - b e^{-at}) / (a - b)
  Variant A  a=12, b=8 per s  (time constants 83 and 125 ms)  [the paper's T36 numbers]
  Variant B  time constants 167 and 250 ms (a=5.99, b=4.0 per s)  [the Loop summary numbers]
Three different calculations were previously quoted as one:
  (1) disjoint 50-ms bins, variant A          -> 2,580 events
  (2) nested empirical survival S(40,80,120ms), variant A -> 4,588 events (different observable, positively correlated)
  (3) disjoint 50-ms bins, variant B          -> 7,873 events
"""
import json, math, os, sys, time
import numpy as np
from scipy import optimize, stats

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'out', 'single_channel_sample_size.json')
QUICK = '--quick' in sys.argv
rng = np.random.default_rng(20261113)
C = np.array([1.0, -2.0, 1.0])


def S(t, a, b):
    return (a * np.exp(-b * t) - b * np.exp(-a * t)) / (a - b)


def d2(x):
    return math.log(x[0]) - 2 * math.log(x[1]) + math.log(x[2])


def disjoint_bins(a, b, w):
    p = [S(w * i, a, b) - S(w * (i + 1), a, b) for i in range(3)]
    var = sum(c * c / q for c, q in zip(C, p))       # N Var(D2)
    snr1 = d2(p) / math.sqrt(var)                      # SNR for N=1 (scales as sqrt N)
    return p, d2(p), var, snr1


def nested_survival(a, b, ts):
    s = [S(t, a, b) for t in ts]
    var = sum(C[i] * C[j] / (s[i] * s[j]) * (s[max(i, j)] * (1 - s[min(i, j)])) for i in range(3) for j in range(3))
    return s, d2(s), var, d2(s) / math.sqrt(var)


def n_for_snr(snr1, target):
    return math.ceil((target / abs(snr1)) ** 2)


def rates_from_tau_ms(t1, t2):
    return 1000.0 / t1, 1000.0 / t2


res = {}
# ---------------------------------------------------------------- reconcile
aA, bA = 12.0, 8.0
aB, bB = rates_from_tau_ms(167.0, 250.0)
rows = {}
for label, (a, b), kind, w in [('A_disjoint50ms', (aA, bA), 'disjoint', 0.05), ('A_nested40_80_120', (aA, bA), 'nested', None),
                               ('B_disjoint50ms', (aB, bB), 'disjoint', 0.05), ('B_nested40_80_120', (aB, bB), 'nested', None)]:
    if kind == 'disjoint':
        p, d, v, s1 = disjoint_bins(a, b, w)
    else:
        p, d, v, s1 = nested_survival(a, b, (0.04, 0.08, 0.12))
    rows[label] = dict(rates_per_s=[round(a, 4), round(b, 4)], observable=kind, bin_or_grid=w if w else [0.04, 0.08, 0.12],
                       cell_probs=[float(x) for x in p], D2=float(d), N_var_D2=float(v), snr_per_sqrtN=float(s1),
                       N_for_SNR5=n_for_snr(s1, 5.0), N_for_SNR3=n_for_snr(s1, 3.0))
res['reconciliation'] = rows
print({k: v['N_for_SNR5'] for k, v in rows.items()})

# ---------------------------------------------------------------- event rate -> recording time
rate_rows = []
for key in ('A_disjoint50ms', 'B_disjoint50ms'):
    for burst_rate in (0.5, 1.0, 2.0, 5.0):          # opening bursts per second in an active single-channel patch
        n = rows[key]['N_for_SNR5']
        rate_rows.append(dict(variant=key, bursts_per_s=burst_rate, N=n, minutes=round(n / burst_rate / 60.0, 1)))
res['recording_time_vs_burst_rate'] = rate_rows
res['recording_time_note'] = ('Time = N / (bursts per second). The old 13-21 min figure corresponds to a burst rate of '
                              f"{rows['A_disjoint50ms']['N_for_SNR5'] / (13 * 60):.1f}-{rows['A_disjoint50ms']['N_for_SNR5'] / (21 * 60):.1f} "
                              'bursts/s, i.e. a very active channel; patches with 0.5-2 bursts/s need 20-90 min under variant A.')

# ---------------------------------------------------------------- calibrated test (type-I error and power by Monte Carlo)
def sample_hypo(n, a, b, size_rng):
    return size_rng.exponential(1 / a, n) + size_rng.exponential(1 / b, n)


def d2_stat(x, w):
    e = np.array([0, w, 2 * w, 3 * w])
    h = np.histogram(x, bins=e)[0].astype(float)
    if np.any(h == 0):
        return np.nan
    return math.log(h[0]) - 2 * math.log(h[1]) + math.log(h[2]), h


def se_plug(h):
    return math.sqrt(np.sum(C ** 2 / h))                     # plug-in Poisson SE of D2 from counts


def nll_exp(x):
    m = x.mean()
    return len(x) * (math.log(m) + 1.0)


def nll_hypo(x):
    """Fit hypoexponential (rates a >= b); nests exponential limit a->inf (boundary)."""
    n = len(x)

    def f(p):
        a, b = np.exp(p)
        if abs(a - b) < 1e-6 * a:
            a = a * (1 + 1e-6)
        pdf = a * b / (a - b) * (np.exp(-b * x) - np.exp(-a * x))
        pdf = np.maximum(pdf, 1e-300)
        return -np.sum(np.log(pdf))
    m = x.mean(); v = x.var()
    best = None
    for fa in (3.0, 1.5, 0.5, 8.0):
        a0, b0 = fa / m, 1.0 / (m * (1 - 1 / (fa + 1) + 1e-3))
        r = optimize.minimize(f, np.log([max(a0, b0 * 1.05), b0]), method='Nelder-Mead', options=dict(xatol=1e-6, fatol=1e-9, maxiter=2000))
        if best is None or r.fun < best.fun:
            best = r
    return min(best.fun, nll_exp(x))


def lrt(x):
    return max(2 * (nll_exp(x) - nll_hypo(x)), 0.0)


def run_power(a, b, w, Ns, nrep, null_rate_ms, dead_ms=0.0):
    out = []
    mean_t = 1 / a + 1 / b
    for N in Ns:
        # null: single exponential with the SAME mean (matched first moment)
        z0, z1, l0, l1 = [], [], [], []
        for _ in range(nrep):
            for kind in (0, 1):
                x = rng.exponential(mean_t, N) if kind == 0 else sample_hypo(N, a, b, rng)
                if dead_ms > 0:
                    x = x[x > dead_ms * 1e-3]
                r = d2_stat(x, w)
                zz = np.nan if r is np.nan or isinstance(r, float) else r[0] / se_plug(r[1])
                ll = lrt(x)
                (z0 if kind == 0 else z1).append(zz); (l0 if kind == 0 else l1).append(ll)
        z0, z1, l0, l1 = map(np.array, (z0, z1, l0, l1))
        # one-sided D2 < 0 is the serial/log-concave direction for hypoexponential bins -> test on -z
        zt = -z0[~np.isnan(z0)]; za = -z1[~np.isnan(z1)]
        out.append(dict(N=N, mc_reps=nrep,
                        d2_null_mean_z=float(np.mean(zt)), d2_null_sd_z=float(np.std(zt, ddof=1)),
                        d2_alt_mean_z=float(np.mean(za)),
                        d2_typeI_z_gt_1p645=float(np.mean(zt > 1.645)), d2_power_z_gt_1p645=float(np.mean(za > 1.645)),
                        d2_typeI_z_gt_2p326=float(np.mean(zt > 2.326)), d2_power_z_gt_2p326=float(np.mean(za > 2.326)),
                        lrt_crit_95=float(np.quantile(l0, 0.95)), lrt_crit_99=float(np.quantile(l0, 0.99)),
                        lrt_typeI_at_own_95=float(np.mean(l0 > np.quantile(l0, 0.95))),
                        lrt_power_at_null95=float(np.mean(l1 > np.quantile(l0, 0.95))),
                        lrt_power_at_null99=float(np.mean(l1 > np.quantile(l0, 0.99)))))
        print('power', a, b, N, out[-1]['d2_power_z_gt_2p326'], out[-1]['lrt_power_at_null99'])
    return out


nrep = 150 if QUICK else 600
Ns = [300, 1000, 2580, 7873] if not QUICK else [300, 1000]
res['calibration'] = dict(
    note=('Null = exponential burst durations with matched mean. alpha is calibrated by Monte Carlo (LRT) or the plug-in '
          'normal approximation (D2), whose empirical type-I error is reported. SNR=5 means E[z]=5; it is not a detection probability.'),
    variant_A_binned_vs_hypoexponential_LRT=run_power(aA, bA, 0.05, Ns, nrep, None),
    variant_B_binned_vs_hypoexponential_LRT=run_power(aB, bB, 0.05, [1000, 7873] if not QUICK else [1000], nrep, None))

# analytic z=3.17 gives 80% power at one-sided alpha=0.01 (z_a + z_b)
zreq = stats.norm.ppf(0.99) + stats.norm.ppf(0.80)
res['detection_probability_conversion'] = dict(
    z_required_for_power0p8_alpha0p01=float(zreq),
    N_variant_A=int(math.ceil(rows['A_disjoint50ms']['N_for_SNR5'] * (zreq / 5) ** 2)),
    N_variant_B=int(math.ceil(rows['B_disjoint50ms']['N_for_SNR5'] * (zreq / 5) ** 2)),
    note='Asymptotic; confirm with the Monte-Carlo power columns above.')

# ---------------------------------------------------------------- missed short events (dead time / filter) sensitivity
dead_rows = []
for dead in (0.0, 1.0, 3.0, 5.0):
    dead_rows.append(dict(dead_ms=dead, A_power_z_gt_2p326_N2580=run_power(aA, bA, 0.05, [2580], 60 if not QUICK else 20, None, dead_ms=dead)[0]['d2_power_z_gt_2p326']))
res['dead_time_sensitivity_N2580_variantA'] = dead_rows

res['limitations'] = [
    'Burst durations are assumed independent with known sequential two-step kinetics; real records contain flickery closures, subconductance states, multiple channels and missed events (see dead_time_sensitivity).',
    'A likelihood-based hypoexponential LRT is the established comparator; it is expected to be at least as powerful when the model is correct. The binned statistic is a model-light screen, not a replacement.',
    'No experimental single-channel recording is analysed in this script.']
os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, 'w', encoding='utf-8') as f:
    json.dump(res, f, indent=1)
print('wrote', OUT)
