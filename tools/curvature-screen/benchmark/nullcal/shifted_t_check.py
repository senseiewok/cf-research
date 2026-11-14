"""Compare null tail fractions with a t_15 (and N(0,1)) shifted by the observed null mean. Reads z_*.npy; writes shifted_t_check.json."""
import json, os, numpy as np
from scipy import stats
H = os.path.dirname(os.path.abspath(__file__)); out = {}
for f in ('z_rev_lumpable_boundary_WC0_n100_T300.npy', 'z_rev_lumpable_boundary_WC0_n100_T1800.npy',
          'z_rev_lumpable_boundary_WC0_n10_T300.npy', 'z_rev_lumpable_boundary_WC0_n10_T1800.npy'):
    z = np.load(os.path.join(H, f)); m = float(z.mean()); d = {'N': len(z), 'mean': m}
    for thr in (-3, -2, 2, 3):
        obs = float(np.mean(z < thr) if thr < 0 else np.mean(z > thr))
        cdf = (lambda F: F.cdf(thr - m) if thr < 0 else F.sf(thr - m))
        d[f'thr_{thr}'] = dict(observed=obs, shifted_normal=float(cdf(stats.norm)), shifted_t15=float(cdf(stats.t(15))))
    d['ks_vs_shifted_t15_p'] = float(stats.kstest(z - m, stats.t(15).cdf).pvalue)
    d['ks_vs_shifted_normal_p'] = float(stats.kstest(z - m, stats.norm.cdf).pvalue)
    out[f] = d; print(f, json.dumps(d))
json.dump(out, open(os.path.join(H, 'shifted_t_check.json'), 'w', encoding='utf-8'), indent=1)
