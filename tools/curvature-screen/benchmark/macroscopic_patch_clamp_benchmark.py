#!/usr/bin/env python3
"""Macroscopic patch-clamp current-noise autocorrelation Wronskian W_C: generating code (nothing hard-coded).

  python macroscopic_patch_clamp_benchmark.py [--quick]      -> out/macroscopic_patch_clamp_benchmark.json

Computes, from rate matrices:
  (1) Reversible-chain null: random reversible Markov chains, W_C >= 0 at equally spaced lags (a theorem; checked numerically).
  (2) A driven 3-state cycle C -> O1 -> O2 -> C built from the *listed* rates: eigenvalues are computed (they are
      complex for these rates), spectral weights, cycle affinity, fraction of lag triples with W_C < 0.
  (3) A real-spectrum driven cycle (parameters searched and recorded) for comparison.
  (4) Time to a 5-sigma W_C detection by Monte Carlo in the Gaussian (many-channel) limit with white measurement noise,
      for n = 1, 10, 100, 1000 channels, plus the scaling check z ~ sqrt(T).
Assumptions are stated in the JSON. W_C < 0 is a sufficient signature of broken detailed balance, not a necessary one.
"""
import json, math, os, sys, time
import numpy as np
from scipy.linalg import expm

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'out', 'macroscopic_patch_clamp_benchmark.json')
QUICK = '--quick' in sys.argv
rng = np.random.default_rng(20261113)


def gen_from_rates(k):
    """k: dict (i,j)->rate for i != j; returns generator Q (rows sum to zero, Q[i,j] = rate i->j)."""
    m = 1 + max(max(i, j) for i, j in k)
    Q = np.zeros((m, m))
    for (i, j), r in k.items():
        Q[i, j] = r
    Q[np.diag_indices(m)] = -Q.sum(1)
    return Q


def stationary(Q):
    w, v = np.linalg.eig(Q.T)
    p = np.real(v[:, np.argmin(abs(w))])
    return p / p.sum()


def autocov(Q, g, taus):
    p = stationary(Q); mu = p @ g; d = g - mu
    return np.array([(p * d) @ expm(Q * t) @ d for t in taus])


def lag_triples(n_t1=15, n_d=15, lo=0.02, hi=0.40):
    return [(t1, d) for t1 in np.linspace(lo, hi, n_t1) for d in np.linspace(lo, hi, n_d)]


def frac_negative(Q, g, tol=1e-14):
    neg = tot = 0
    for t1, d in lag_triples():
        c = autocov(Q, g, [t1, t1 + d, t1 + 2 * d])
        W = c[0] * c[2] - c[1] ** 2
        tot += 1; neg += W < -tol * max(abs(c[1]) ** 2, 1e-300)
    return neg / tot


def spectral(Q, g):
    """Return eigenvalues (nonzero), and the weights a_k of C(tau) = sum a_k exp(lambda_k tau) (complex)."""
    w, V = np.linalg.eig(Q)
    Vi = np.linalg.inv(V)
    p = stationary(Q); mu = p @ g; d = g - mu
    a = (p * d) @ V * (Vi @ d)
    order = np.argsort(abs(w))
    keep = [i for i in order if abs(w[i]) > 1e-9]
    return w[keep], a[keep]


def random_reversible(m):
    pi = rng.dirichlet(np.ones(m))
    S = rng.uniform(0.2, 5.0, (m, m)); S = (S + S.T) / 2
    Q = np.zeros((m, m))
    for i in range(m):
        for j in range(m):
            if i != j:
                Q[i, j] = S[i, j] / pi[i] ** 0.5 * pi[j] ** 0.5       # detailed balance: pi_i Q_ij = pi_j Q_ji
    Q[np.diag_indices(m)] = -Q.sum(1)
    return Q, pi


res = {'assumptions': [
    'Stationary, ergodic, finite-state Markov channel; n independent identical channels; additive white measurement noise; current = sum of conductances.',
    'W_C(t1,d) = C(t1) C(t1+2d) - C(t1+d)^2 on equally spaced lags; reversible => C completely monotone => W_C >= 0.',
    'Filter (anti-aliasing) effects are not modelled here; lags must lie well beyond the filter time constant.',
    'No experimental recording is analysed here. This is a simulation benchmark of the statistic.']}

# ---------------------------------------------------------------- (1) reversible null
t0 = time.time()
n_chains = 300 if QUICK else 2000
viol = tot = 0; min_rel = np.inf
for _ in range(n_chains):
    m = int(rng.choice([3, 4]))
    Q, pi = random_reversible(m)
    assert np.allclose(pi[:, None] * Q - (pi[:, None] * Q).T, 0, atol=1e-9), 'detailed balance construction failed'
    g = rng.uniform(0, 1, m)
    trip = lag_triples(6, 6) if QUICK else lag_triples(9, 9)
    for t1, d in trip:
        c = autocov(Q, g, [t1, t1 + d, t1 + 2 * d]); W = c[0] * c[2] - c[1] ** 2
        rel = W / max(c[1] ** 2, 1e-300)
        min_rel = min(min_rel, rel); tot += 1; viol += rel < -1e-9
res['reversible_null'] = dict(n_chains=n_chains, states=[3, 4], triples_per_chain=len(trip), evaluations=tot,
                              violations_W_C_lt_minus1e9_rel=int(viol), min_W_over_C2sq=float(min_rel),
                              note='min W/C(t2)^2 ~ 0 up to float error; sign constraint holds.')
print('reversible', res['reversible_null'], round(time.time() - t0, 1), 's')

# ---------------------------------------------------------------- (2),(3) driven cycles
g = np.array([0.0, 1.0, 1.0])


def describe(name, rates):
    Q = gen_from_rates(rates)
    ev, a = spectral(Q, g)
    A = math.log(rates[(0, 1)] * rates[(1, 2)] * rates[(2, 0)] / (rates[(1, 0)] * rates[(2, 1)] * rates[(0, 2)]))
    f = frac_negative(Q, g)
    tri = (0.05, 0.20, 0.35)
    c = autocov(Q, g, tri); Wopt = c[0] * c[2] - c[1] ** 2
    return dict(name=name, rates_per_s={f'{i}->{j}': r for (i, j), r in rates.items()},
                eigenvalues_nonzero=[[float(z.real), float(z.imag)] for z in ev],
                eigenvalues_strictly_real=bool(np.all(abs(ev.imag) < 1e-9)),
                spectral_weights=[[float(z.real), float(z.imag)] for z in a],
                cycle_affinity_kBT=float(A), fraction_lag_triples_W_C_lt_0=float(f),
                W_C_at_50_200_350ms_per_unit_current_sq=float(Wopt),
                stationary=[float(x) for x in stationary(Q)])


listed = {(0, 1): 2.5, (1, 0): 0.05, (1, 2): 4.0, (2, 1): 0.05, (2, 0): 3.5, (0, 2): 0.05}
res['driven_cycle_listed_rates'] = describe('listed rates (k12=2.5,k23=4.0,k31=3.5; reverse 0.05 /s)', listed)
print('listed eig', res['driven_cycle_listed_rates']['eigenvalues_nonzero'], 'A', res['driven_cycle_listed_rates']['cycle_affinity_kBT'])

# real-spectrum driven cycle: first candidate from a seeded search with A >= 6 and strictly real spectrum, max W_C<0 fraction
best = None
for _ in range(1500 if QUICK else 6000):
    r = np.exp(rng.uniform(math.log(0.02), math.log(30.0), 6))
    k = {(0, 1): r[0], (1, 0): r[1], (1, 2): r[2], (2, 1): r[3], (2, 0): r[4], (0, 2): r[5]}
    Q = gen_from_rates(k)
    if np.all(abs(np.linalg.eigvals(Q).imag) < 1e-9):
        A = math.log(r[0] * r[2] * r[4] / (r[1] * r[3] * r[5]))
        if A >= 6.0:
            # keep rates in a CFTR-like gating range (slowest relaxation 0.5-10 /s, fastest <= 40 /s)
            ev = np.sort(-np.linalg.eigvals(Q).real)[1:]
            if 0.5 <= ev[0] <= 10 and ev[1] <= 40:
                f = frac_negative(Q, g)
                if best is None or f > best[0]:
                    best = (f, k)
res['driven_cycle_real_spectrum'] = describe('seeded search: strictly real spectrum, A>=6, slow rate 0.5-10/s, fast <=40/s', best[1])
res['real_spectrum_search_seed'] = 20261113

# affinity sweep on the listed topology (scale the forward rates)
sweep = []
for scale in (1.0, 2.0, 5.0, 10.0, 30.0, 100.0):
    k = dict(listed)
    for key in ((0, 1), (1, 2), (2, 0)):
        k[key] = listed[key] * scale / 1.0
    d = describe(f'forward x{scale}', k)
    sweep.append(dict(forward_scale=scale, affinity_kBT=d['cycle_affinity_kBT'], strictly_real=d['eigenvalues_strictly_real'],
                      fraction_W_C_lt_0=d['fraction_lag_triples_W_C_lt_0']))
res['affinity_sweep_listed_topology'] = sweep

# ---------------------------------------------------------------- (4) detection time, Gaussian many-channel limit
dt = 0.002; lag_idx = (25, 100, 175)          # 50, 200, 350 ms
K = 16


def circulant_sampler(Cfun_vals, L):
    """Gaussian stationary process of length L with covariance Cfun_vals[0..L-1] (circulant embedding)."""
    lam = np.fft.rfft(np.concatenate([Cfun_vals[:L // 2 + 1], Cfun_vals[1:L // 2][::-1]]))
    lam = np.maximum(lam.real, 0)
    def draw():
        z = rng.standard_normal(L)
        return np.fft.irfft(np.fft.rfft(z) * np.sqrt(lam), n=L) * 1.0
    return draw


def simulate_z(Q, n, sigma_w, T, reps):
    L = int(T / K / dt)
    L -= L % 2
    taus = np.arange(L // 2 + 1) * dt
    ev_, a_ = spectral(Q, g)
    Cs = np.real(np.exp(np.outer(taus, ev_)) @ a_) * n
    Cs_full = Cs.copy()
    Cs_full[0] += sigma_w ** 2                             # white measurement noise at lag 0 only
    draw = circulant_sampler(Cs_full, L)
    zs = []
    for _ in range(reps):
        blk = np.empty((K, 3))
        for k in range(K):
            x = draw(); x -= x.mean()
            f = np.fft.rfft(x, n=2 * L)
            ac = np.fft.irfft(f * np.conj(f))[:max(lag_idx) + 1] / L
            blk[k] = ac[list(lag_idx)]
        vals = []
        for a in range(K):
            for b in range(a + 1, K):
                vals.append(0.5 * (blk[a, 0] * blk[b, 2] + blk[b, 0] * blk[a, 2]) - blk[a, 1] * blk[b, 1])
        u = np.mean(vals)
        jk = []
        for j in range(K):
            idx = [i for i in range(K) if i != j]
            s = blk[idx]
            v = [0.5 * (s[a, 0] * s[b, 2] + s[b, 0] * s[a, 2]) - s[a, 1] * s[b, 1] for a in range(K - 1) for b in range(a + 1, K - 1)]
            jk.append(np.mean(v))
        jk = np.array(jk)
        se = math.sqrt((K - 1) / K * np.sum((jk - jk.mean()) ** 2))
        zs.append(u / se)
    return np.array(zs)


timing = {}
for label, key in (('listed_complex_spectrum', 'driven_cycle_listed_rates'), ('real_spectrum', 'driven_cycle_real_spectrum')):
    rates = {tuple(int(s) for s in kk.split('->')): v for kk, v in res[key]['rates_per_s'].items()}
    Q = gen_from_rates(rates)
    per_n = []
    for n in ((1, 10, 100, 1000) if not QUICK else (10, 100)):
        T0 = 300.0
        z = simulate_z(Q, n, 1.0, T0, 12 if QUICK else 40)
        per_n.append(dict(n_channels=n, sigma_w_over_i=1.0, T_s=T0, mean_z=float(z.mean()), sd_z=float(z.std(ddof=1)),
                          frac_z_lt_minus2=float(np.mean(z < -2)), frac_z_gt_2=float(np.mean(z > 2))))
        print(label, per_n[-1])
    timing[label] = dict(per_n=per_n)
    # time to 5 sigma by sqrt(T) scaling from n=100, then a confirmation run at the predicted time
    base = next(x for x in per_n if x['n_channels'] == 100)
    z0 = abs(base['mean_z'])
    T5 = base['T_s'] * (5.0 / max(z0, 1e-6)) ** 2
    timing[label]['T5sigma_s_by_sqrtT_scaling_n100'] = float(T5)
    timing[label]['T5sigma_minutes'] = float(T5 / 60)
    Tcur = T5; hist = []
    for it in range(3 if not QUICK else 0):
        if Tcur > 6 * 3600:
            break
        zc = simulate_z(Q, 100, 1.0, Tcur, 30)
        hist.append(dict(T_s=float(Tcur), mean_z=float(zc.mean()), sd_z=float(zc.std(ddof=1))))
        Tcur = Tcur * (5.0 / max(abs(zc.mean()), 1e-6)) ** 2
    timing[label]['refinement_runs'] = hist
    timing[label]['T_mean_z_5_after_refinement_s'] = float(Tcur)
    timing[label]['T_mean_z_5_after_refinement_minutes'] = float(Tcur / 60)
res['detection_time_gaussian_limit'] = dict(
    method='Circulant-embedding Gaussian process with exact C(tau) from the rate matrix + white noise sigma_w; K=16 block U-statistic; jackknife SE',
    caveat='Gaussian (many-channel) limit; finite-n non-Gaussianity, filtering, baseline drift and rundown are not modelled. z is signed: W_C < 0 gives negative z.',
    results=timing)

with open(OUT, 'w', encoding='utf-8') as f:
    json.dump(res, f, indent=1)
print('wrote', OUT, round(time.time() - t0, 1), 's')
