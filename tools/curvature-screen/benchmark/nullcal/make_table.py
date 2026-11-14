#!/usr/bin/env python3
"""Render nullcal_results.json (+ bias_diagnostic.json) into nullcal_table.md (UTF-8)."""
import json, os
HERE = os.path.dirname(os.path.abspath(__file__))
r = json.load(open(os.path.join(HERE, 'nullcal_results.json'), encoding='utf-8'))
b = json.load(open(os.path.join(HERE, 'bias_diagnostic.json'), encoding='utf-8'))
L = ['# Null calibration of the W_C jackknife z (package Part B pipeline)', '',
     f"Reproduction: whole package JSON identical = {r['reproduction'].get('whole_json_identical')}; "
     f"listed n=100 T=300 mean_z = {r['reproduction']['scratch_copy_full_run_listed_n100_T300']['mean_z']:.4f}", '',
     '## All runs (JSON `runs[i]`)', '',
     '| i | config | n | T (s) | reps | W_C pop/ch^2 | mean z (SE) | SD z | z<-2 (SE) | z<-3 (SE) | z>+2 (SE) | z>+3 (SE) | q05 | q01 |',
     '|---|---|---|---|---|---|---|---|---|---|---|---|---|---|']
f = lambda d: f"{d['frac']:.4f} ({d['se']:.4f})"
for i, s in enumerate(r['runs']):
    L.append(f"| {i} | {s['config']} | {s['n_channels']} | {s['T_s']:.0f} | {s['N']} | {s['W_C_population_per_channel_sq']:.2e} | "
             f"{s['mean_z']:.3f} ({s['se_mean']:.3f}) | {s['sd_z']:.3f} | {f(s['frac_z_lt_m2'])} | {f(s['frac_z_lt_m3'])} | "
             f"{f(s['frac_z_gt_p2'])} | {f(s['frac_z_gt_p3'])} | {s['q05']:.3f} | {s['q01']:.3f} |")
ref = r['reference_tails']
L += ['', 'Reference tails: ' + '; '.join(f"{k}: " + ', '.join(f"{kk}={vv:.4f}" for kk, vv in v.items()) for k, v in ref.items()), '',
      '## Power at calibrated thresholds (JSON `power_calibrated[i]`)', '',
      '| i | T (s) | n | alpha | least-fav. null | threshold z (95% CI) | power (SE) | power range over thr CI | power @ N(0,1) thr | power @ t15 thr | power @ z<-3 |',
      '|---|---|---|---|---|---|---|---|---|---|---|']
for i, p in enumerate(r['power_calibrated']):
    c = p['threshold_bootstrap_95ci']
    L.append(f"| {i} | {p['T_s']:.0f} | {p['n_channels']} | {p['alpha_one_sided']} | {p['least_favourable_null']} | {p['threshold_z']:.3f} ({c[0]:.3f}, {c[1]:.3f}) | "
             f"{p['power']:.3f} ({p['power_binom_se']:.3f}) | {p['power_range_over_threshold_ci'][0]:.3f}-{p['power_range_over_threshold_ci'][1]:.3f} | "
             f"{p['power_at_nominal_normal']:.3f} | {p['power_at_nominal_t15']:.3f} | {p['power_at_minus3']:.3f} |")
x = r['xpcs_short_lag_multiplicity']
L += ['', '## XPCS short-lag multiplicity (JSON `xpcs_short_lag_multiplicity`)', '',
      '| q (1/A) | z | p one-sided (t11) | p two-sided (t11) | Bonf 1s | Holm 1s |', '|---|---|---|---|---|---|']
for k in range(10):
    L.append(f"| {x['q_angstrom_inv'][k]} | {x['z'][k]} | {x['p_one_sided_lower'][k]:.4f} | {x['p_two_sided'][k]:.4f} | "
             f"{x['bonferroni_reject_one_sided'][k]} | {x['holm_reject_one_sided'][k]} |")
L += ['', f"t11 critical (Bonferroni 0.005): one-sided {x['t11_crit_one_sided']:.3f}, two-sided {x['t11_crit_two_sided']:.3f}. "
      f"Sign test {x['sign_test']['n_negative']}/10 negative: p one-sided {x['sign_test']['p_one_sided']:.4f}, two-sided {x['sign_test']['p_two_sided']:.4f}.", '',
      '## Bias diagnostic (bias_diagnostic.json; boundary null n=100)', '']
for k, v in b.items():
    L.append(f"- {k}: package mean z {v['package_mean_z']:.3f} (z<-2 {v['package_frac_lt_m2']:.3f}); without block-mean subtraction "
             f"mean z {v['no_mean_subtraction_mean_z']:.3f} (z<-2 {v['no_mean_subtraction_frac_lt_m2']:.3f}); reps {v['reps']}")
open(os.path.join(HERE, 'nullcal_table.md'), 'w', encoding='utf-8').write('\n'.join(L) + '\n')
print('\n'.join(L))
s = json.load(open(os.path.join(HERE, 'shifted_t_check.json'), encoding='utf-8'))
M = ['', '## Shape of the null z (shifted_t_check.json; boundary null): observed vs N(mean,1) and mean + t15', '']
for k, v in s.items():
    M.append(f"- {k}: N={v['N']}, mean={v['mean']:.3f}; " + '; '.join(
        f"{t}: obs {v[t]['observed']:.4f} / shifted N {v[t]['shifted_normal']:.4f} / shifted t15 {v[t]['shifted_t15']:.4f}"
        for t in ('thr_-3', 'thr_-2', 'thr_2', 'thr_3')) + f"; KS p (shifted t15) {v['ks_vs_shifted_t15_p']:.3f}, (shifted N) {v['ks_vs_shifted_normal_p']:.2g}")
open(os.path.join(HERE, 'nullcal_table.md'), 'a', encoding='utf-8').write('\n'.join(M) + '\n')
