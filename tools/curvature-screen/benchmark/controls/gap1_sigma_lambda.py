"""GAP 1: recompute sigma_lambda and sigma_lambda,diff for the Hellenkamp 2018 smFRET mixtures with the exact
windows recorded in the ledger (anchor = peak_bin, offset 50 bins, M = 30 bins, D = 60 / 75 bins, 16 ps bins,
background bins 2200..3000, B_eff = 800), using verbatim copies of wstats / wdiff_stats and the exact
sigma formulas of positive_control_hellenkamp2018_block() (validate_fluorescence_benchmark.py lines 507-524).
Deterministic (no randomness). Reads the cache read-only.
Run: python -X utf8 -B gap1_sigma_lambda.py
"""
import json, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from refstats import wstats, wdiff_stats, r, sha256, HELLENKAMP_NPZ, LEDGER_JSON, REF_SCRIPT

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'gap1_sigma_lambda.json')
QUOTED = {  # numbers quoted in the briefing / FINDINGS / earlier draft, for comparison only
    '1a_1b_Mix.ht3': dict(sigma_lambda_diff_quoted=[0.2966, 0.2583], sigma_lambda_unsplit_quoted=[0.3245, 0.2848],
                          burst_mix_sigma_lambda_diff_quoted_FINDINGS=[0.2738, 0.2584]),
    '2a_2b_Mix.ht3': dict(sigma_lambda_diff_quoted=[0.3204, 0.2907], sigma_lambda_diff_earlier_draft=[0.3016, 0.2827],
                          sigma_lambda_unsplit_quoted=[0.3566, 0.3151]),
}
npz = np.load(HELLENKAMP_NPZ, allow_pickle=False)
ledger = json.load(open(LEDGER_JSON, encoding='utf-8'))['positive_control_hellenkamp2018']['files']
bg_slice, B_eff = slice(2200, 3000), 800


def one(hA, hB, pk, D_bins, off=50, M=30):
    """Exact copy of the per-window computation in positive_control_hellenkamp2018_block (lines 509-524)."""
    bA, bB = float(hA[bg_slice].mean()), float(hB[bg_slice].mean())
    w4 = [slice(pk + off + D_bins * k, pk + off + D_bins * k + M) for k in range(4)]
    s3 = wstats(hA, hB, w4[:3], bA=bA, bB=bB, B_eff=B_eff)
    sd = wdiff_stats(hA, hB, w4)
    D_ns = D_bins * 0.016
    yA2 = float(hA[w4[1]].sum() - M * bA)
    yB2 = float(hB[w4[1]].sum() - M * bB)
    dA2, dB2 = float(sd['dA'][1]), float(sd['dB'][1])
    sig_lam = (1.0 / D_ns) * np.sqrt(max(0.0, s3['Ws'] / max(yA2 * yB2, 1.0))) if yA2 * yB2 > 0 else 0.0
    sig_lam_d = (1.0 / D_ns) * np.sqrt(max(0.0, sd['W'] / max(dA2 * dB2, 1.0))) if dA2 * dB2 > 0 else 0.0
    return dict(windows=[[s.start, s.stop] for s in w4], W_split=s3['Ws'], z_split_exact=s3['zs_exact'],
                W_rel=s3['W_rel'], sigma_lambda_ns_inv=sig_lam,
                W_diff_split=sd['W'], z_diff_split_exact=sd['z_exact'], W_diff_rel=sd['W_rel'],
                sigma_lambda_diff_ns_inv=sig_lam_d, yA2=yA2, yB2=yB2, dA2=dA2, dB2=dB2)


out = dict(script='benchmark/controls/' + os.path.basename(__file__), reference_script_sha256=sha256(REF_SCRIPT),
           cache_npz=HELLENKAMP_NPZ, cache_sha256=sha256(HELLENKAMP_NPZ),
           formula=dict(sigma_lambda='(1/(D*0.016 ns)) * sqrt(max(0, W_split / max(yA2*yB2, 1)))',
                        sigma_lambda_diff='(1/(D*0.016 ns)) * sqrt(max(0, W_diff_split / max(dA2*dB2, 1)))',
                        source='validate_fluorescence_benchmark.py lines 517-518'),
           files={})
for tag, fname in (('1a_1b', '1a_1b_Mix.ht3'), ('2a_2b', '2a_2b_Mix.ht3')):
    pk = int(npz[f'{tag}_peak_bin'])
    fo = dict(anchor_peak_bin=pk, sub_ensembles={})
    for sub in ('full_stream', 'burst_mix_1ms', 'low_FRET_1ms', 'high_FRET_1ms', 'donly_1ms'):
        hA = npz[f'{tag}_{sub}_hA'].astype(float); hB = npz[f'{tag}_{sub}_hB'].astype(float)
        so = {}
        for D_bins in (60, 75):
            res = one(hA, hB, pk, D_bins)
            key = f'D{D_bins}_bins_{D_bins * 16}ps'
            led = ledger[fname]['sub_ensembles'][sub]['lifetime_matched_windows'][key]
            res_r = {k: (r(v, 4) if isinstance(v, float) else v) for k, v in res.items()}
            res_r['ledger_sigma_lambda_ns_inv'] = led['sigma_lambda_ns_inv']
            res_r['ledger_sigma_lambda_diff_ns_inv'] = led['sigma_lambda_diff_ns_inv']
            res_r['matches_ledger'] = bool(r(res['sigma_lambda_ns_inv'], 4) == led['sigma_lambda_ns_inv']
                                           and r(res['sigma_lambda_diff_ns_inv'], 4) == led['sigma_lambda_diff_ns_inv']
                                           and r(res['z_split_exact'], 3) == led['z_split_exact']
                                           and r(res['z_diff_split_exact'], 3) == led['z_diff_split_exact'])
            so[key] = res_r
        fo['sub_ensembles'][sub] = so
    fo['quoted_in_briefing'] = QUOTED[fname]
    fs = fo['sub_ensembles']['full_stream']
    fo['full_stream_recomputed'] = dict(
        sigma_lambda_D60_D75=[fs['D60_bins_960ps']['sigma_lambda_ns_inv'], fs['D75_bins_1200ps']['sigma_lambda_ns_inv']],
        sigma_lambda_diff_D60_D75=[fs['D60_bins_960ps']['sigma_lambda_diff_ns_inv'], fs['D75_bins_1200ps']['sigma_lambda_diff_ns_inv']])
    out['files'][fname] = fo
    print(fname, fo['full_stream_recomputed'])

# Which quoted sigma_lambda,diff values appear ANYWHERE among the ledger-formula outputs above (any sub-ensemble,
# either sigma, 4-decimal match)? This is a search for the source, not a re-tuning of windows.
all_vals = {}
for fname, fo in out['files'].items():
    for sub, so in fo['sub_ensembles'].items():
        for key, res in so.items():
            for s in ('sigma_lambda_ns_inv', 'sigma_lambda_diff_ns_inv'):
                all_vals[f'{fname}/{sub}/{key}/{s}'] = res[s]
hits = {}
for q in (0.2966, 0.2583, 0.3204, 0.2907, 0.3016, 0.2827, 0.2738, 0.2584):
    hits[str(q)] = [k for k, v in all_vals.items() if v is not None and abs(v - q) < 0.00005]
out['quoted_value_search_in_ledger_windows'] = hits
json.dump(out, open(OUT, 'w', encoding='utf-8'), indent=1)
print('hits', hits)
