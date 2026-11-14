#!/usr/bin/env python3
"""Null calibration and power of the noisy W_C jackknife z test (package Part B pipeline).

Run:  python -X utf8 -B nullcal.py [--quick]
Reads (read-only) the package script and its JSON; writes nullcal_results.json and nullcal_table.md next to this file.
Estimator functions come from wc_pipeline.py (verbatim copies, checked against the package source by AST below).
Seeds: chain construction 20261113 + offsets; simulation chunks SeedSequence([20261199, config_id, chunk_id]).
"""
import ast, hashlib, json, math, os, sys, time
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import wc_pipeline as W

PKG = os.path.dirname(HERE)
PKG_SCRIPT = os.path.join(PKG, 'macroscopic_patch_clamp_benchmark.py')
PKG_JSON = os.path.join(PKG, 'out', 'macroscopic_patch_clamp_benchmark.json')
REPRO_JSON = os.path.join(HERE, 'repro_macroscopic_patch_clamp_benchmark.json')   # from running pkg_repro.py (copy, OUT redirected)
OUT_JSON = os.path.join(HERE, 'nullcal_results.json')
OUT_MD = os.path.join(HERE, 'nullcal_table.md')
QUICK = '--quick' in sys.argv
REUSE = '--reuse' in sys.argv     # skip simulation; reuse saved z_*.npy (same seeds => same numbers)
CHUNK = 250
BASE_SEED = 20261199

LISTED = {(0, 1): 2.5, (1, 0): 0.05, (1, 2): 4.0, (2, 1): 0.05, (2, 0): 3.5, (0, 2): 0.05}


def sha(p):
    return hashlib.sha256(open(p, 'rb').read()).hexdigest()


def check_verbatim():
    """Compare AST dumps of each copied function with the package's definitions."""
    pkg = {n.name: ast.dump(n) for n in ast.parse(open(PKG_SCRIPT, encoding='utf-8').read()).body if isinstance(n, ast.FunctionDef)}
    mine = {n.name: ast.dump(n) for n in ast.parse(open(os.path.join(HERE, 'wc_pipeline.py'), encoding='utf-8').read()).body
            if isinstance(n, ast.FunctionDef)}
    names = ['gen_from_rates', 'stationary', 'autocov', 'spectral', 'random_reversible', 'circulant_sampler', 'simulate_z']
    out = {nm: (nm in pkg and pkg[nm] == mine.get(nm)) for nm in names}
    src = open(PKG_SCRIPT, encoding='utf-8').read()
    out['globals_match'] = all(s in src for s in ["g = np.array([0.0, 1.0, 1.0])", "dt = 0.002; lag_idx = (25, 100, 175)", "K = 16"])
    return out


# ------------------------------------------------------------------ chains
def rev_same_pi_from_listed():
    """Detailed-balance chain with the listed cycle's stationary distribution and the same per-edge total rate."""
    pi = W.stationary(W.gen_from_rates(LISTED))
    k = {}
    for i, j in ((0, 1), (1, 2), (2, 0)):
        s = LISTED[(i, j)] + LISTED[(j, i)]
        k[(i, j)] = s * pi[j] / (pi[i] + pi[j]); k[(j, i)] = s * pi[i] / (pi[i] + pi[j])
    return k


def rev_geo_from_listed():
    """Keep the listed forward rates; set all three reverse rates to (k01 k12 k20)^(1/3) so the cycle affinity is 0."""
    r = (LISTED[(0, 1)] * LISTED[(1, 2)] * LISTED[(2, 0)]) ** (1 / 3)
    return {(0, 1): 2.5, (1, 2): 4.0, (2, 0): 3.5, (1, 0): r, (2, 1): r, (0, 2): r}


def rev_lumpable_boundary():
    """Reversible chain whose open states O1, O2 have equal exit rate to C: current is an exact 2-state process,
    C(tau) is a single exponential and W_C = 0 exactly (the least-favourable reversible null)."""
    return {(0, 1): 2.0, (1, 0): 3.0, (0, 2): 1.0, (2, 0): 3.0, (1, 2): 1.0, (2, 1): 2.0}


def random_rev_chains(nchains):
    W.set_seed(20261113 + 7)
    out = []
    for _ in range(nchains):
        Q, pi = W.random_reversible(3)
        out.append({(i, j): float(Q[i, j]) for i in range(3) for j in range(3) if i != j})
    return out


def describe_chain(k):
    Q = W.gen_from_rates(k)
    pi = W.stationary(Q)
    A = math.log(k[(0, 1)] * k[(1, 2)] * k[(2, 0)] / (k[(1, 0)] * k[(2, 1)] * k[(0, 2)]))
    db = float(np.max(np.abs(pi[:, None] * Q - (pi[:, None] * Q).T)))
    ev, a = W.spectral(Q, W.g)
    taus = np.array(W.lag_idx) * W.dt
    c = np.real(np.exp(np.outer(taus, ev)) @ a)           # same C(tau) construction as simulate_z (per channel)
    Wp = c[0] * c[2] - c[1] ** 2
    return dict(rates_per_s={f'{i}->{j}': round(v, 6) for (i, j), v in k.items()}, stationary=[float(x) for x in pi],
                cycle_affinity_kBT=float(A), max_detailed_balance_residual=db,
                eigenvalues_nonzero=[[float(z.real), float(z.imag)] for z in ev],
                spectral_weights=[[float(z.real), float(z.imag)] for z in a],
                C_at_lags_per_channel=[float(x) for x in c], W_C_population_per_channel_sq=float(Wp),
                W_C_rel=float(Wp / c[1] ** 2))


# ------------------------------------------------------------------ simulation
def _work(args):
    cfg_id, chunk_id, rates, n, T, reps = args
    W.set_seed(np.random.SeedSequence([BASE_SEED, cfg_id, chunk_id]))
    return cfg_id, W.simulate_z(W.gen_from_rates({tuple(map(int, kk.split(','))): v for kk, v in rates.items()}), n, 1.0, T, reps)


def binom(k, N):
    p = k / N
    return dict(count=int(k), frac=p, se=math.sqrt(max(p * (1 - p), 1.0 / N ** 2) / N))


def summarise(z):
    N = len(z)
    d = dict(N=N, mean_z=float(z.mean()), sd_z=float(z.std(ddof=1)), se_mean=float(z.std(ddof=1) / math.sqrt(N)),
             median_z=float(np.median(z)),
             q01=float(np.quantile(z, 0.01)), q05=float(np.quantile(z, 0.05)), q95=float(np.quantile(z, 0.95)), q99=float(np.quantile(z, 0.99)))
    for lab, m in (('lt_m2', z < -2), ('lt_m3', z < -3), ('gt_p2', z > 2), ('gt_p3', z > 3)):
        d['frac_z_' + lab] = binom(int(m.sum()), N)
    return d


def reference_tails(df):
    return {name: {'P(z<-2)': float(f.cdf(-2)), 'P(z<-3)': float(f.cdf(-3)), 'P(z>2)': float(f.sf(2)), 'P(z>3)': float(f.sf(3)),
                   'q01': float(f.ppf(0.01)), 'q05': float(f.ppf(0.05)), 'sd': float(f.std())}
            for name, f in (('normal', stats.norm), (f't_{df}', stats.t(df)))}


def tail_comparison(s, df):
    """z-scores of observed tail fractions against N(0,1) and t_df expectations."""
    ref = reference_tails(df)
    out = {}
    for key, rk in (('frac_z_lt_m2', 'P(z<-2)'), ('frac_z_lt_m3', 'P(z<-3)'), ('frac_z_gt_p2', 'P(z>2)'), ('frac_z_gt_p3', 'P(z>3)')):
        obs = s[key]['count']; N = s['N']
        row = {}
        for nm, r in ref.items():
            p = r[rk]; row[nm] = dict(expected=p, expected_count=p * N, binom_two_sided_p=float(stats.binomtest(obs, N, p).pvalue))
        out[key] = row
    return out


def main():
    t0 = time.time()
    res = dict(script='nullcal.py', python=sys.version.split()[0], numpy=np.__version__, scipy=__import__('scipy').__version__,
               package_script=PKG_SCRIPT, package_script_sha256=sha(PKG_SCRIPT), package_json_sha256=sha(PKG_JSON),
               verbatim_function_check=check_verbatim(), quick=QUICK, chunk_reps=CHUNK, seeds=dict(chains=20261113 + 7, sim_base=BASE_SEED))
    assert all(res['verbatim_function_check'].values()), res['verbatim_function_check']

    # ---------------- reproduction check
    pkg = json.load(open(PKG_JSON, encoding='utf-8'))
    pk100 = next(x for x in pkg['detection_time_gaussian_limit']['results']['listed_complex_spectrum']['per_n'] if x['n_channels'] == 100)
    rep = dict(package_listed_n100_T300=pk100)
    if os.path.exists(REPRO_JSON):
        rj = json.load(open(REPRO_JSON, encoding='utf-8'))
        r100 = next(x for x in rj['detection_time_gaussian_limit']['results']['listed_complex_spectrum']['per_n'] if x['n_channels'] == 100)
        rep['scratch_copy_full_run_listed_n100_T300'] = r100
        rep['exact_match_mean_z'] = r100['mean_z'] == pk100['mean_z']
        rep['whole_json_identical'] = rj == pkg
    res['reproduction'] = rep

    # ---------------- chains
    nulls = {'rev_from_listed_same_pi': rev_same_pi_from_listed(), 'rev_from_listed_geo_reverse': rev_geo_from_listed(),
             'rev_lumpable_boundary_WC0': rev_lumpable_boundary()}
    for i, k in enumerate(random_rev_chains(6)):
        nulls[f'rev_random_{i + 1}'] = k
    chains = dict(nulls); chains['driven_listed'] = dict(LISTED)
    res['chains'] = {nm: describe_chain(k) for nm, k in chains.items()}
    for nm in nulls:
        assert res['chains'][nm]['max_detailed_balance_residual'] < 1e-9 and abs(res['chains'][nm]['cycle_affinity_kBT']) < 1e-9

    # ---------------- design (replicates)
    if REUSE:      # recompute only the post-processing from the saved z_*.npy and the previous JSON's run summaries
        prev = json.load(open(OUT_JSON, encoding='utf-8'))
        res['runs'] = runs = prev['runs']; res['workers'] = prev['workers']; res['reference_tails'] = prev['reference_tails']
        res['simulation_elapsed_s_previous'] = prev['elapsed_s']
    else:
        sc =0.1 if QUICK else 1.0
        key_nulls = ['rev_lumpable_boundary_WC0', 'rev_from_listed_same_pi', 'rev_from_listed_geo_reverse']
        design = []
        for T in (300.0, 1800.0):
            for n in (100, 10):
                for nm in key_nulls:
                    base = {(300.0, 100): 10000, (1800.0, 100): 4000, (300.0, 10): 4000, (1800.0, 10): 2000}[(T, n)]
                    if nm == 'rev_lumpable_boundary_WC0':
                        base = {(300.0, 100): 40000, (1800.0, 100): 10000, (300.0, 10): 10000, (1800.0, 10): 4000}[(T, n)]
                    design.append((nm, n, T, base))
                if n == 100:
                    for i in range(6):
                        design.append((f'rev_random_{i + 1}', n, T, 2000 if T == 300.0 else 1000))
                design.append(('driven_listed', n, T, 4000 if T == 300.0 else 2000))
        design = [(nm, n, T, max(CHUNK, int(r * sc) // CHUNK * CHUNK)) for nm, n, T, r in design]
        tasks = []
        for cid, (nm, n, T, reps) in enumerate(design):
            rates = {f'{i},{j}': v for (i, j), v in chains[nm].items()}
            for ch in range(reps // CHUNK):
                tasks.append((cid, ch, rates, n, T, CHUNK))
        tasks.sort(key=lambda a: -a[4])          # long records first
        zs = {cid: [] for cid in range(len(design))}
        nw = max(1, (os.cpu_count() or 2) - 2)
        with ProcessPoolExecutor(max_workers=nw) as ex:
            for cid, z in ex.map(_work, tasks, chunksize=1):
                zs[cid].append(z)
        res['workers'] = nw
    
        runs = []
        for cid, (nm, n, T, reps) in enumerate(design):
            z = np.concatenate(zs[cid])
            s = summarise(z)
            s.update(config=nm, n_channels=n, T_s=T, sigma_w_over_i=1.0, K=W.K, lags_s=[x * W.dt for x in W.lag_idx], dt_s=W.dt,
                     W_C_population_per_channel_sq=res['chains'][nm]['W_C_population_per_channel_sq'])
            if nm != 'driven_listed':
                s['tail_vs_reference'] = tail_comparison(s, W.K - 1)
            runs.append(s)
            np.save(os.path.join(HERE, f'z_{nm}_n{n}_T{int(T)}.npy'), z)
        res['runs'] = runs
        res['reference_tails'] = reference_tails(W.K - 1)

    # ---------------- calibrated thresholds and power
    rng = np.random.default_rng(BASE_SEED + 1)
    power = []
    for T in (300.0, 1800.0):
        for n in (100, 10):
            nul = {r['config']: np.load(os.path.join(HERE, f"z_{r['config']}_n{n}_T{int(T)}.npy"))
                   for r in runs if r['n_channels'] == n and r['T_s'] == T and r['config'] != 'driven_listed'}
            zd = np.load(os.path.join(HERE, f'z_driven_listed_n{n}_T{int(T)}.npy'))
            for alpha in (0.05, 0.01):
                # least favourable null: the reversible configuration with the MOST NEGATIVE lower-tail quantile
                # (threshold = min over adequately resolved nulls, so the FPR is <= alpha for every one of them)
                qs_all = {nm: float(np.quantile(z, alpha)) for nm, z in nul.items()}
                big = max(len(z) for z in nul.values())
                qs = {nm: q for nm, q in qs_all.items() if len(nul[nm]) >= 0.2 * big}   # adequately resolved nulls only
                lf = min(qs, key=qs.get)
                zn = nul[lf]
                thr = qs[lf]
                boot = np.array([np.quantile(rng.choice(zn, len(zn)), alpha) for _ in range(400)])
                pw = float(np.mean(zd < thr))
                pw_lo = float(np.mean(zd < np.quantile(boot, 0.025))); pw_hi = float(np.mean(zd < np.quantile(boot, 0.975)))
                power.append(dict(T_s=T, n_channels=n, alpha_one_sided=alpha, least_favourable_null=lf, null_quantiles_by_config=qs_all, configs_eligible_for_threshold=sorted(qs),
                                  null_N=int(len(zn)), threshold_z=thr, threshold_bootstrap_95ci=[float(np.quantile(boot, 0.025)), float(np.quantile(boot, 0.975))],
                                  driven_N=int(len(zd)), power=pw, power_binom_se=math.sqrt(pw * (1 - pw) / len(zd)),
                                  power_range_over_threshold_ci=[pw_lo, pw_hi],
                                  power_at_nominal_normal=float(np.mean(zd < stats.norm.ppf(alpha))),
                                  power_at_nominal_t15=float(np.mean(zd < stats.t(W.K - 1).ppf(alpha))),
                                  power_at_minus3=float(np.mean(zd < -3))))
    res['power_calibrated'] = power

    # ---------------- XPCS multiplicity (arithmetic only)
    zx = [1.01, -2.12, -0.34, -0.18, -1.77, -2.77, -2.12, -1.95, -1.39, -0.94]
    qx = [0.0074, 0.0102, 0.0129, 0.0157, 0.0185, 0.0213, 0.024, 0.0268, 0.0296, 0.0324]
    t11 = stats.t(11)
    p1 = [float(t11.cdf(z)) for z in zx]                    # one-sided, H1: W_diff < 0
    p2 = [float(2 * t11.sf(abs(z))) for z in zx]

    def holm(p, a=0.05):
        order = np.argsort(p); rej = [False] * len(p)
        for r, i in enumerate(order):
            if p[i] <= a / (len(p) - r):
                rej[i] = True
            else:
                break
        return rej
    nneg = sum(z < 0 for z in zx)
    res['xpcs_short_lag_multiplicity'] = dict(
        z=zx, q_angstrom_inv=qx, df=11, p_one_sided_lower=p1, p_two_sided=p2,
        min_p_one_sided=min(p1), min_p_two_sided=min(p2),
        bonferroni_alpha_per_test=0.005, t11_crit_one_sided=float(t11.ppf(0.005)), t11_crit_two_sided=float(t11.ppf(0.0025)),
        bonferroni_reject_one_sided=[p <= 0.005 for p in p1], bonferroni_reject_two_sided=[p <= 0.005 for p in p2],
        holm_reject_one_sided=holm(p1), holm_reject_two_sided=holm(p2),
        sign_test=dict(n_negative=nneg, n=10, p_one_sided=float(stats.binomtest(nneg, 10, 0.5, alternative='greater').pvalue),
                       p_two_sided=float(stats.binomtest(nneg, 10, 0.5).pvalue),
                       assumption='independent rings AND null median of z exactly 0; neither is established'),
        stouffer_if_independent_illustration=dict(z=float(sum(zx) / math.sqrt(10)), p_one_sided=float(stats.norm.cdf(sum(zx) / math.sqrt(10)))),
        other_regimes_tested=['two_pop_intermediate_t0_5s_M30s_D60s', 'two_pop_macro_lag_t0_5s_M250s_D500s', 'short_lag_compressed_exp_t0_5s_M20s_D20s'],
        statistics_per_regime=['W_C_U', 'W_C_diff_U'])
    res['elapsed_s'] = time.time() - t0
    with open(OUT_JSON, 'w', encoding='utf-8') as f:
        json.dump(res, f, indent=1)
    print('wrote', OUT_JSON, round(res['elapsed_s'], 1), 's')


if __name__ == '__main__':
    main()
