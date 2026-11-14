"""GAP 3: photon-matched specificity test for the tissue FLIM superpixel map (Zenodo 7656540, cached 8x8 grid).

The cache holds per-tile checkerboard-split decay histograms (grid8_chk_hA/hB: 64 x 272; grid16: 256 x 272), so
the test can be done. Estimators: verbatim copies of wstats / wdiff_stats; per-tile background = mean of bins
0..5 of each half, B_eff = 6; a truncated statistic is set to 0 (ledger rule, validate_fluorescence_benchmark.py
lines 881-889). Windows: the three superpixel windows of the ledger (off10_M10_D12, off10_M12_D15, off16_M12_D16),
anchor = whole-image peak_bin.

(1) Ledger reproduction: dark tiles (N < 1,500), stained tiles (N >= 20,000), counts above z thresholds.
(2) Downward: every stained tile (N >= 20,000) is thinned WITHOUT replacement (multivariate hypergeometric over
    the 2 x 272 (half, bin) cells, so the checkerboard split is thinned with the photons) to N_target total photons,
    N_target in {1,500, 5,000} (requested) and {median dark-tile N} (pre-declared extra level, the dark tiles'
    own typical count). R = 1,000 repeats per tile and level, fixed seed. The SAME thinned sample is evaluated
    on all three windows. Detection rates at z > 2 and z > 3 are reported per window for z_split and z_diff,split,
    pooled over tiles x repeats, with the per-tile spread and the fraction of truncated evaluations (z := 0).
(3) Upward: the dark tiles are pooled into one histogram (sum of hA, sum of hB) and evaluated; also the 13 tiles
    with N < 2,500 (the ledger's W_diff,split low-count set).
Seed 20261113. Run: python -X utf8 -B gap3_tissue_thinning.py
"""
import json, os, sys, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from refstats import wstats, wdiff_stats, r, sha256, TISSUE_NPZ, LEDGER_JSON, REF_SCRIPT

SEED = 20261113
R = 1000
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'gap3_tissue_thinning.json')
npz = np.load(TISSUE_NPZ, allow_pickle=False)
pk = int(npz['peak_bin'])
B_eff = 6
WINS = [('off10_M10_D12', 10, 10, 12), ('off10_M12_D15', 10, 12, 15), ('off16_M12_D16', 16, 12, 16)]
HA = npz['grid8_chk_hA'].astype(np.int64); HB = npz['grid8_chk_hB'].astype(np.int64)
NBIN = HA.shape[1]
N_tile = (HA + HB).sum(axis=1)
dark = np.where(N_tile < 1500)[0]; low2500 = np.where(N_tile < 2500)[0]
stained = np.where(N_tile >= 20000)[0]
ledger = json.load(open(LEDGER_JSON, encoding='utf-8'))['dataset3_zenodo7656540_leica_stellaris8_tissue_flim'][
    'lc_fce_dynamical_microscope_spatial_demonstration']['grid_8x8_superpixels']['windows']


def tile_z(ha, hb, off, M, D):
    ha = ha.astype(float); hb = hb.astype(float)
    w4 = [slice(pk + off + k * D, pk + off + k * D + M) for k in range(4)]
    ba, bb = float(ha[:6].mean()), float(hb[:6].mean())
    s3 = wstats(ha, hb, w4[:3], bA=ba, bB=bb, B_eff=B_eff)
    sd = wdiff_stats(ha, hb, w4)
    zs = s3['zs_exact'] if not s3['truncated'] else 0.0
    zd = sd['z_exact'] if not sd['truncated'] else 0.0
    return zs, zd, bool(sd['truncated']), s3['W_rel'], sd['W_rel']


out = dict(script='benchmark/controls/' + os.path.basename(__file__), reference_script_sha256=sha256(REF_SCRIPT), cache_npz=TISSUE_NPZ,
           cache_sha256=sha256(TISSUE_NPZ), seed=SEED, repeats_per_tile_and_level=R, grid='8x8 (64 tiles, 64x64 px each)',
           split='checkerboard (cached grid8_chk_hA / grid8_chk_hB)', anchor_peak_bin=pk, B_eff=B_eff,
           n_dark_lt_1500=int(len(dark)), dark_tile_N=[int(N_tile[i]) for i in dark],
           n_lt_2500=int(len(low2500)), n_stained_ge_20000=int(len(stained)),
           stained_tile_N_min_median_max=[int(N_tile[stained].min()), int(np.median(N_tile[stained])), int(N_tile[stained].max())])

# (1) ledger reproduction at full counts
repro = {}
for lab, off, M, D in WINS:
    zs = np.array([tile_z(HA[i], HB[i], off, M, D)[0] for i in range(64)])
    zd = np.array([tile_z(HA[i], HB[i], off, M, D)[1] for i in range(64)])
    led = ledger[lab]
    mine = dict(dark_fp_z_split_gt2=int((zs[dark] > 2).sum()), dark_fp_z_diff_gt2=int((zd[dark] > 2).sum()),
                dark_max_z_split=r(zs[dark].max(), 2), dark_max_z_diff=r(zd[dark].max(), 2),
                lt2500_fp_z_diff_gt2=int((zd[low2500] > 2).sum()),
                stained_n_z_split_gt3=int((zs[stained] > 3).sum()), stained_n_z_split_gt2=int((zs[stained] > 2).sum()),
                stained_n_z_diff_gt2=int((zd[stained] > 2).sum()), stained_n_z_diff_gt3=int((zd[stained] > 3).sum()),
                dark_z_split=[r(x, 2) for x in zs[dark]], dark_z_diff=[r(x, 2) for x in zd[dark]])
    ledger_vals = dict(dark_fp_z_split_gt2=led['dark_lumen_bg_lt_1500']['false_positives_z_split_gt_2'],
                       dark_fp_z_diff_gt2=led['dark_lumen_bg_lt_1500']['false_positives_z_diff_split_gt_2'],
                       lt2500_fp_z_diff_gt2=led['low_count_bg_lt_2500']['false_positives_z_diff_split_gt_2'],
                       stained_n_z_split_gt3=led['stained_tissue_fg_ge_20000']['n_z_split_gt_3'],
                       stained_n_z_diff_gt2=led['stained_tissue_fg_ge_20000']['n_z_diff_split_gt_2'])
    mine['matches_ledger'] = all(mine[k] == v for k, v in ledger_vals.items())
    mine['ledger'] = ledger_vals
    repro[lab] = mine
    print('repro', lab, mine['matches_ledger'], {k: mine[k] for k in ledger_vals})
out['full_count_reproduction'] = repro

# (2) downward thinning of stained tiles
g = np.random.default_rng(SEED)
levels = [1500, 5000, int(np.median(N_tile[dark]))]
out['thinning_levels'] = dict(requested=[1500, 5000], extra_median_dark_tile_N=levels[2])
thin = {}
t0 = time.time()
for Nt in levels:
    ZS = np.zeros((len(WINS), len(stained), R)); ZD = np.zeros_like(ZS); TR = np.zeros_like(ZS, dtype=bool)
    for ti, i in enumerate(stained):
        colors = np.concatenate([HA[i], HB[i]])
        for rep in range(R):
            s = g.multivariate_hypergeometric(colors, Nt)
            ha, hb = s[:NBIN], s[NBIN:]
            for wi, (lab, off, M, D) in enumerate(WINS):
                zs, zd, tr, _, _ = tile_z(ha, hb, off, M, D)
                ZS[wi, ti, rep] = zs; ZD[wi, ti, rep] = zd; TR[wi, ti, rep] = tr
    lev = {}
    for wi, (lab, off, M, D) in enumerate(WINS):
        zs, zd = ZS[wi], ZD[wi]
        per_tile_split2 = (zs > 2).mean(axis=1); per_tile_diff2 = (zd > 2).mean(axis=1)
        lev[lab] = dict(
            n_evaluations=int(zs.size),
            rate_z_split_gt2=r((zs > 2).mean(), 4), rate_z_split_gt3=r((zs > 3).mean(), 4),
            rate_z_diff_gt2=r((zd > 2).mean(), 4), rate_z_diff_gt3=r((zd > 3).mean(), 4),
            rate_z_split_lt_m2=r((zs < -2).mean(), 4), rate_z_diff_lt_m2=r((zd < -2).mean(), 4),
            mc_se_max=r(0.5 / np.sqrt(zs.size), 4),
            median_z_split=r(np.median(zs), 3), median_z_diff=r(np.median(zd), 3),
            per_tile_rate_z_split_gt2_min_median_max=[r(per_tile_split2.min(), 3), r(np.median(per_tile_split2), 3), r(per_tile_split2.max(), 3)],
            per_tile_rate_z_diff_gt2_min_median_max=[r(per_tile_diff2.min(), 3), r(np.median(per_tile_diff2), 3), r(per_tile_diff2.max(), 3)],
            fraction_truncated_diff_set_to_zero=r(TR[wi].mean(), 4),
            prob_0_of_n_dark_exceed_z_split_gt2_if_dark_were_like_thinned_stained=r((1 - (zs > 2).mean()) ** len(dark), 4),
            prob_0_of_n_dark_exceed_z_diff_gt2_if_dark_were_like_thinned_stained=r((1 - (zd > 2).mean()) ** len(dark), 4))
        print(f'N={Nt} {lab}: split>2 {lev[lab]["rate_z_split_gt2"]} >3 {lev[lab]["rate_z_split_gt3"]} | '
              f'diff>2 {lev[lab]["rate_z_diff_gt2"]} >3 {lev[lab]["rate_z_diff_gt3"]} | trunc {lev[lab]["fraction_truncated_diff_set_to_zero"]}'
              f' ({time.time() - t0:.0f}s)')
    thin[f'N_{Nt}'] = lev
out['stained_tiles_thinned'] = thin

# (3) upward: pooled dark tiles and pooled N < 2,500 tiles
pooled = {}
for name, idx in (('dark_lt_1500_pooled', dark), ('lt_2500_pooled', low2500)):
    ha, hb = HA[idx].sum(axis=0), HB[idx].sum(axis=0)
    pooled[name] = dict(n_tiles=int(len(idx)), N_photons=int((ha + hb).sum()), windows={})
    for lab, off, M, D in WINS:
        zs, zd, tr, wr, wdr = tile_z(ha, hb, off, M, D)
        pooled[name]['windows'][lab] = dict(z_split_exact=r(zs, 2), z_diff_split_exact=r(zd, 2), W_rel=r(wr, 4),
                                            W_diff_rel=r(wdr, 4), diff_truncated=tr)
    print(name, pooled[name])
out['pooled_low_count_tiles'] = pooled

# (4) thinned stained tiles at the POOLED dark-tile count, for a like-for-like comparison with (3)
Np = pooled['dark_lt_1500_pooled']['N_photons']
ZS = np.zeros((len(WINS), len(stained), R)); ZD = np.zeros_like(ZS)
for ti, i in enumerate(stained):
    colors = np.concatenate([HA[i], HB[i]])
    for rep in range(R):
        s = g.multivariate_hypergeometric(colors, Np)
        for wi, (lab, off, M, D) in enumerate(WINS):
            zs, zd, _, _, _ = tile_z(s[:NBIN], s[NBIN:], off, M, D)
            ZS[wi, ti, rep] = zs; ZD[wi, ti, rep] = zd
out['stained_tiles_thinned_to_pooled_dark_N'] = dict(N_target=Np, windows={
    lab: dict(rate_z_split_gt2=r((ZS[wi] > 2).mean(), 4), rate_z_split_gt3=r((ZS[wi] > 3).mean(), 4),
              rate_z_diff_gt2=r((ZD[wi] > 2).mean(), 4), rate_z_diff_gt3=r((ZD[wi] > 3).mean(), 4),
              median_z_split=r(np.median(ZS[wi]), 3), median_z_diff=r(np.median(ZD[wi]), 3),
              fraction_thinned_z_split_le_pooled_dark=r((ZS[wi] <= pooled['dark_lt_1500_pooled']['windows'][lab]['z_split_exact']).mean(), 4),
              fraction_thinned_z_diff_le_pooled_dark=r((ZD[wi] <= pooled['dark_lt_1500_pooled']['windows'][lab]['z_diff_split_exact']).mean(), 4))
    for wi, (lab, off, M, D) in enumerate(WINS)})
print('thinned to pooled dark N', out['stained_tiles_thinned_to_pooled_dark_N'])
json.dump(out, open(OUT, 'w', encoding='utf-8'), indent=1)
