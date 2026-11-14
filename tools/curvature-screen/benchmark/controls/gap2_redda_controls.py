"""GAP 2: matched controls for the RED_DA (sensitized acceptor) onset-window headline.

Every channel in the cache that has no FRET feeding is run through the SAME estimators (verbatim copies of
wstats / wdiff_stats) with the SAME anchor / offset / M / D grid that the ledger used for RED_DA
(validate_fluorescence_benchmark.py lines 559-577: anchor 804, windows anchor+off+D*k .. +M, k = 0..3).

Channels (per file, 1a_1b_Mix and 2a_2b_Mix):
  RED_DA            red detectors, green excitation (the headline channel); bg bins 100..700, B_eff 600 (ledger)
  RED_AA            red detectors, red PIE excitation (direct acceptor, no FRET feeding); same red bg convention.
                    Bin 804 is meaningless for RED_AA, so two anchors are reported:
                      'peak_minus_lag' = red_aa_peak - (red_da_peak - 804)  (same position relative to own peak as RED_DA)
                      'own_peak'       = red_aa_peak                         (same position relative to own peak as the donor)
  DONOR_full_stream donor detectors, full stream (no feeding; decays from its own excitation); anchor 804
  DONOR_donly_1ms   donor detectors, donor-only bursts (no acceptor at all); anchor 804
  DONOR_low_FRET_1ms / DONOR_high_FRET_1ms   donor detectors of FRET bursts (quenched, but no feeding); anchor 804
  Donor background: ledger donor convention bins 2200..3000, B_eff 800 (primary); 100..700 / 600 (sensitivity).
  NOT cached: the red channel of donor-only bursts, any IRF/scatter measurement. W_diff,split needs no background.

Tables:
  A. ledger grid exactly as recorded for RED_DA (off0 M15 D20, off0 M15 D15, off5 M15 D20, off10 M15 D20,
     off50 M15 D60, off50 M30 D60), every channel.
  B. ONE (M, D) pair for both files, M = 15, D = 20, at offsets 0, 5, 10, 50 from the anchor, every channel.
  C. matched time since each channel's OWN peak: anchor = own peak, offsets -30..+80 bins step 5, M = 15, D = 20.
Relative statistics: W_rel = W_split / ybar2^2, W_diff_rel = W_diff,split / dbar2^2 (count-free).
z_scaled_to_RED_DA_counts = z * sqrt(N4_RED_DA / N4_channel) (N4 = raw counts in the four windows, both halves):
an approximate count normalisation (z ~ sqrt(N) at fixed shape), labelled as such.
Deterministic: no random numbers. Run: python -X utf8 -B gap2_redda_controls.py
"""
import json, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from refstats import wstats, wdiff_stats, r, sha256, HELLENKAMP_NPZ, LEDGER_JSON, REF_SCRIPT

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'gap2_redda_controls.json')
npz = np.load(HELLENKAMP_NPZ, allow_pickle=False)
ledger = json.load(open(LEDGER_JSON, encoding='utf-8'))['positive_control_hellenkamp2018']['files']
ANC = 804
LEDGER_GRID = [('feeding_zone_off0_M15_D20', 0, 15, 20), ('feeding_zone_off0_M15_D15', 0, 15, 15),
               ('feeding_zone_off5_M15_D20', 5, 15, 20), ('feeding_zone_off10_M15_D20', 10, 15, 20),
               ('late_tail_off50_M15_D60', 50, 15, 60), ('late_tail_off50_M30_D60', 50, 30, 60)]
ONE_MD = (15, 20)
UNIFORM_OFFSETS = (0, 5, 10, 50)
SCAN_OFFSETS = list(range(-30, 81, 5))


def calc(hA, hB, anchor, off, M, D, bg, B_eff):
    bA, bB = float(hA[bg].mean()), float(hB[bg].mean())
    w4 = [slice(anchor + off + D * k, anchor + off + D * k + M) for k in range(4)]
    s3 = wstats(hA, hB, w4[:3], bA=bA, bB=bB, B_eff=B_eff)
    sd = wdiff_stats(hA, hB, w4)
    n4 = float(sum(hA[s].sum() + hB[s].sum() for s in w4))
    return dict(anchor_bin=int(anchor), offset_bins=off, M_bins=M, D_bins=D,
                windows=[[s.start, s.stop] for s in w4], N4_counts=int(n4),
                W_split=r(s3['Ws'], 1), z_split_exact=r(s3['zs_exact'], 2), W_rel=r(s3['W_rel'], 4),
                W_diff_split=r(sd['W'], 1), z_diff_split_exact=r(sd['z_exact'], 2), W_diff_rel=r(sd['W_rel'], 4),
                floor_binding_diff=sd['floor_binding'], truncated=bool(s3['truncated'] or sd['truncated']))


def channels(tag):
    rA = npz[f'{tag}_red_acceptor_hA'].astype(float); rB = npz[f'{tag}_red_acceptor_hB'].astype(float)
    pk_da, pk_aa = int(npz[f'{tag}_red_da_peak_bin']), int(npz[f'{tag}_red_aa_peak_bin'])
    pk_d = int(npz[f'{tag}_peak_bin'])
    red_bg, red_B = slice(100, 700), 600
    don_bg, don_B = slice(2200, 3000), 800
    ch = {
        'RED_DA': dict(hA=rA, hB=rB, anchor=ANC, own_peak=pk_da, bg=red_bg, B_eff=red_B, feeding=True,
                       note='sensitized acceptor (headline channel)'),
        'RED_AA_anchor_peak_minus_lag': dict(hA=rA, hB=rB, anchor=pk_aa - (pk_da - ANC), own_peak=pk_aa, bg=red_bg,
                                             B_eff=red_B, feeding=False,
                                             note=f'direct acceptor; anchor = {pk_aa} - {pk_da - ANC} (RED_DA lag before own peak)'),
        'RED_AA_anchor_own_peak': dict(hA=rA, hB=rB, anchor=pk_aa, own_peak=pk_aa, bg=red_bg, B_eff=red_B,
                                       feeding=False, note='direct acceptor; anchor = own peak'),
    }
    for sub in ('full_stream', 'donly_1ms', 'low_FRET_1ms', 'high_FRET_1ms'):
        hA = npz[f'{tag}_{sub}_hA'].astype(float); hB = npz[f'{tag}_{sub}_hB'].astype(float)
        own = int(np.argmax(hA + hB))
        ch[f'DONOR_{sub}'] = dict(hA=hA, hB=hB, anchor=ANC, own_peak=own, bg=don_bg, B_eff=don_B, feeding=False,
                                  note=f'donor detectors, {sub}; file donor peak_bin {pk_d}; raw argmax {own}')
        ch[f'DONOR_{sub}_bg100_700'] = dict(hA=hA, hB=hB, anchor=ANC, own_peak=own, bg=red_bg, B_eff=red_B,
                                            feeding=False, note='background sensitivity: bins 100..700, B_eff 600')
    return ch, pk_da, pk_aa, pk_d


out = dict(script='benchmark/controls/' + os.path.basename(__file__), reference_script_sha256=sha256(REF_SCRIPT),
           cache_npz=HELLENKAMP_NPZ, cache_sha256=sha256(HELLENKAMP_NPZ),
           not_in_cache=['red channel of donor-only bursts (only the full-stream red histogram is cached)',
                         'any IRF / scatter / buffer-only measurement'],
           one_M_D_pair=dict(M=ONE_MD[0], D=ONE_MD[1]), files={})
for tag, fname in (('1a_1b', '1a_1b_Mix.ht3'), ('2a_2b', '2a_2b_Mix.ht3')):
    ch, pk_da, pk_aa, pk_d = channels(tag)
    fo = dict(red_da_peak=pk_da, red_aa_peak=pk_aa, donor_peak_bin=pk_d,
              channel_notes={k: v['note'] for k, v in ch.items()},
              channel_own_peaks={k: v['own_peak'] for k, v in ch.items()},
              A_ledger_grid={}, B_uniform_M15_D20={}, C_scan_since_own_peak_M15_D20={})
    led = ledger[fname]['red_acceptor_modalities']['red_sensitized_acceptor_RED_DA']['windows']
    # reproduction check of the ledger RED_DA numbers
    chk = {}
    for lab, off, M, D in LEDGER_GRID:
        c = calc(ch['RED_DA']['hA'], ch['RED_DA']['hB'], ANC, off, M, D, slice(100, 700), 600)
        chk[lab] = dict(recomputed=[c['z_split_exact'], c['z_diff_split_exact'], c['W_rel'], c['W_diff_rel']],
                        ledger=[led[lab]['z_split_exact'], led[lab]['z_diff_split_exact'], led[lab]['W_rel'], led[lab]['W_diff_rel']])
        chk[lab]['match'] = chk[lab]['recomputed'] == chk[lab]['ledger']
    fo['ledger_reproduction_RED_DA'] = chk
    for name, c in ch.items():
        fo['A_ledger_grid'][name] = {lab: calc(c['hA'], c['hB'], c['anchor'], off, M, D, c['bg'], c['B_eff'])
                                     for lab, off, M, D in LEDGER_GRID}
        fo['B_uniform_M15_D20'][name] = {f'off{off}': calc(c['hA'], c['hB'], c['anchor'], off, *ONE_MD, c['bg'], c['B_eff'])
                                         for off in UNIFORM_OFFSETS}
        fo['C_scan_since_own_peak_M15_D20'][name] = {
            f'off{off:+d}': calc(c['hA'], c['hB'], c['own_peak'], off, *ONE_MD, c['bg'], c['B_eff']) for off in SCAN_OFFSETS}
    # approximate count normalisation relative to RED_DA (same window labels)
    for tbl in ('A_ledger_grid', 'B_uniform_M15_D20', 'C_scan_since_own_peak_M15_D20'):
        ref = fo[tbl]['RED_DA']
        for name, rows in fo[tbl].items():
            for lab, row in rows.items():
                f = np.sqrt(ref[lab]['N4_counts'] / max(row['N4_counts'], 1))
                row['z_split_scaled_to_RED_DA_counts'] = r(row['z_split_exact'] * f, 2)
                row['z_diff_split_scaled_to_RED_DA_counts'] = r(row['z_diff_split_exact'] * f, 2)
    # last offset (since own peak) at which W_rel / W_diff_rel is negative, and first sign change to positive
    summ = {}
    for name, rows in fo['C_scan_since_own_peak_M15_D20'].items():
        offs = SCAN_OFFSETS
        wr = [rows[f'off{o:+d}']['W_rel'] for o in offs]
        wd = [rows[f'off{o:+d}']['W_diff_rel'] for o in offs]
        zs = [rows[f'off{o:+d}']['z_split_exact'] for o in offs]
        zd = [rows[f'off{o:+d}']['z_diff_split_exact'] for o in offs]
        summ[name] = dict(
            last_offset_W_split_significantly_negative_z_lt_m2=max([o for o, z in zip(offs, zs) if z is not None and z < -2], default=None),
            last_offset_W_diff_significantly_negative_z_lt_m2=max([o for o, z in zip(offs, zd) if z is not None and z < -2], default=None),
            W_rel_by_offset=dict(zip([f'{o:+d}' for o in offs], wr)),
            W_diff_rel_by_offset=dict(zip([f'{o:+d}' for o in offs], wd)))
    fo['C_summary_negativity_persistence'] = summ
    # D. matched lag before own peak: anchor = own_peak - (red_da_peak - 804) for every channel (for RED_DA this is
    #    exactly the ledger anchor 804), offsets 0/5/10/50, M = 15, D = 20.
    lag = pk_da - ANC
    fo['D_matched_lag_before_own_peak_M15_D20'] = dict(lag_bins=lag, rows={
        name: {f'off{off}': calc(c['hA'], c['hB'], c['own_peak'] - lag, off, *ONE_MD, c['bg'], c['B_eff'])
               for off in UNIFORM_OFFSETS} for name, c in ch.items()})
    ref = fo['D_matched_lag_before_own_peak_M15_D20']['rows']['RED_DA']
    for name, rows in fo['D_matched_lag_before_own_peak_M15_D20']['rows'].items():
        for lab, row in rows.items():
            f = np.sqrt(ref[lab]['N4_counts'] / max(row['N4_counts'], 1))
            row['z_split_scaled_to_RED_DA_counts'] = r(row['z_split_exact'] * f, 2)
            row['z_diff_split_scaled_to_RED_DA_counts'] = r(row['z_diff_split_exact'] * f, 2)
    # E. descriptive rise times (background-subtracted, both halves pooled): bins from 10 % to 90 % of own peak height
    rise = {}
    for name, c in ch.items():
        if name.endswith('bg100_700'):
            continue
        h = c['hA'] + c['hB']; b = float(h[c['bg']].mean()); pk = c['own_peak']; amp = h[pk] - b
        i10 = next(i for i in range(pk, 0, -1) if h[i] - b < 0.1 * amp)
        i90 = next(i for i in range(pk, 0, -1) if h[i] - b < 0.9 * amp)
        rise[name] = dict(own_peak=pk, bin_10pct=int(i10), bin_90pct=int(i90), rise_10_90_bins=int(i90 - i10),
                          rise_10_90_ps=float((i90 - i10) * 16))
    fo['E_rise_time_10_90'] = rise
    out['files'][fname] = fo

json.dump(out, open(OUT, 'w', encoding='utf-8'), indent=1)

# compact console view
for fname, fo in out['files'].items():
    print('==', fname, 'peaks DA/AA/donor', fo['red_da_peak'], fo['red_aa_peak'], fo['donor_peak_bin'],
          'ledger repro', all(v['match'] for v in fo['ledger_reproduction_RED_DA'].values()))
    for tbl in ('A_ledger_grid', 'B_uniform_M15_D20', 'D_matched_lag_before_own_peak_M15_D20'):
        print(' --', tbl)
        tb = fo[tbl]['rows'] if tbl.startswith('D_') else fo[tbl]
        for name, rows in tb.items():
            if name.endswith('bg100_700'):
                continue
            print(f'  {name:32s}', ' | '.join(f"{lab.replace('feeding_zone_','').replace('late_tail_','')}: "
                                              f"{row['z_split_exact']:+.1f}/{row['z_diff_split_exact']:+.1f} "
                                              f"rel {row['W_rel']:+.3f}/{row['W_diff_rel'] if row['W_diff_rel'] is not None else float('nan'):+.3f}"
                                              for lab, row in rows.items()))
    print(' -- C last offset since own peak with z<-2 (split, diff)')
    for name, s in fo['C_summary_negativity_persistence'].items():
        if name.endswith('bg100_700'):
            continue
        print(f'  {name:32s}', s['last_offset_W_split_significantly_negative_z_lt_m2'], s['last_offset_W_diff_significantly_negative_z_lt_m2'])
    print(' -- E rise 10-90 %', {k: (v['own_peak'], v['rise_10_90_bins']) for k, v in fo['E_rise_time_10_90'].items()})
