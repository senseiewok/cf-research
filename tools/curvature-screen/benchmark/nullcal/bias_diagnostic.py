#!/usr/bin/env python3
"""Diagnostic only (NOT the package estimator): is the negative null offset of z caused by per-block mean subtraction?

Runs the W_C = 0 boundary null (n = 100, T = 300 s) through
  (a) the verbatim package simulate_z (wc_pipeline.simulate_z), and
  (b) a copy differing in exactly one line: the block mean is NOT subtracted (the true mean, 0, is known in simulation).
Same seeds for (a) and (b), so both see identical Gaussian draws.  Writes bias_diagnostic.json.
Also gives the first-order prediction of the offset: E[ac_k] ~ (1-k/L) C_k - S/L (S = sum_j C_j over all lags, two-sided),
so E[W_hat] ~ W - (S/L)(C1 + C3 - 2 C2) (the U-statistic removes the within-block product bias).
"""
import inspect, json, math, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import wc_pipeline as W
import nullcal

src = inspect.getsource(W.simulate_z)
assert src.count('x -= x.mean()') == 1
ns = {}
exec(src.replace('x = draw(); x -= x.mean()', 'x = draw()').replace('def simulate_z(', 'def simulate_z_nomean('), W.__dict__, ns)
nomean = ns['simulate_z_nomean']

Q = W.gen_from_rates(nullcal.rev_lumpable_boundary())
out = {}
for T in (300.0, 1800.0):
    reps = 3000 if T == 300.0 else 600
    W.set_seed(np.random.SeedSequence([20261201, int(T)])); za = W.simulate_z(Q, 100, 1.0, T, reps)
    W.set_seed(np.random.SeedSequence([20261201, int(T)])); zb = nomean(Q, 100, 1.0, T, reps)
    # first-order prediction of the offset in W and in z
    L = int(T / W.K / W.dt); L -= L % 2
    ev, a = W.spectral(Q, W.g)
    taus = np.arange(L) * W.dt
    C = np.real(np.exp(np.outer(taus, ev)) @ a) * 100; C[0] += 1.0
    S = C[0] + 2 * C[1:].sum()
    c = C[list(W.lag_idx)]
    out[f'T{int(T)}'] = dict(reps=reps, package_mean_z=float(za.mean()), package_sd_z=float(za.std(ddof=1)),
                            package_frac_lt_m2=float(np.mean(za < -2)),
                            no_mean_subtraction_mean_z=float(zb.mean()), no_mean_subtraction_sd_z=float(zb.std(ddof=1)),
                            no_mean_subtraction_frac_lt_m2=float(np.mean(zb < -2)),
                            predicted_W_offset=float(-(S / L) * (c[0] + c[2] - 2 * c[1])))
    print(T, out[f'T{int(T)}'])
json.dump(out, open(os.path.join(HERE, 'bias_diagnostic.json'), 'w', encoding='utf-8'), indent=1)
