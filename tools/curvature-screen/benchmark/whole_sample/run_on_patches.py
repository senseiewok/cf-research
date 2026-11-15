#!/usr/bin/env python3
"""Apply the fixed whole-sample analysis to a table of burst durations.

Implements, from specs/CFTR_DATA_REQUEST_AND_ANALYSIS_PLAN.md: the dead-time rule of step 2 with the resolution
supplied by the user or derived from the filter cut-off (0.36/f_c default), steps 3 to 5, the instrument-calibrated
false-alarm threshold for dead time, filtering, and burst criterion (T-0923), the per-construct burst minimums of step 1,
and the 1x / 2x resolution runs.
Not implemented: the one-active-channel check (from metadata), the control power simulation, and the
likelihood-ratio comparator (step 6).

Input : CSV with a header and the columns  patch_id,duration  (seconds) and optionally  construct.
Output: JSON. Per construct and per patch: bursts used, CV^2, its 1% critical value and the decision; the combined
        result over patches; the same at twice the resolution on the same patches; and a diagnostic for a
        resolution that looks too small.

    python benchmark/whole_sample/run_on_patches.py bursts.csv --dead-time 0.002 --control E1371Q --out result.json
    python benchmark/whole_sample/run_on_patches.py bursts.csv --filter-cutoff 1000 --burst-criterion 0.005
    python benchmark/whole_sample/run_on_patches.py --self-test

If durations are recorded on a frame grid, pass --quantum (the frame interval, seconds); durations must then lie
within 1% of that grid; the dead time is rounded up to the grid and one further frame is discarded,
because the first frame above the dead time is only partly observed.  Nothing in this file has been run on a CFTR recording.
"""
import argparse, csv, json, math, sys
from collections import defaultdict
import numpy as np

SEED = 20261113
REPS = 20000
ALPHA = 0.01
MIN_BURSTS = 50            # plan step 1, wild type
MIN_BURSTS_CONTROL = 5     # plan step 1, control construct


def cv2(x, axis=-1):
    m = x.mean(axis=axis)
    return x.var(axis=axis, ddof=1) / m**2


class Null:
    """Monte Carlo null of CV^2 for a single exponential or passive mixture calibrated for instrument
    dead time, filter cutoff (0.36/f_c), and burst criterion (t_crit).
    Each (N, quantum/mean, dead_time, filter_cutoff, burst_criterion) has its own generator, so the
    result does not depend on the order in which patches are read."""

    def __init__(self, seed=SEED):
        self.seed, self.cache = seed, {}

    def draws(self, n, q_over_mean=0.0, dead_time=0.0, filter_cutoff=None, burst_criterion=None, mean_duration=1.0):
        fc_val = float(filter_cutoff) if filter_cutoff and filter_cutoff > 0 else 0.0
        tc_val = float(burst_criterion) if burst_criterion and burst_criterion > 0 else 0.0
        dt_val = float(dead_time) if dead_time and dead_time > 0 else 0.0
        key = (
            int(n),
            int(round(float(q_over_mean) * 1e4)),
            int(round(dt_val * 1e6)),
            int(round(fc_val * 10)),
            int(round(tc_val * 1e6)),
        )
        if key not in self.cache:
            rng = np.random.default_rng([self.seed, key[0], key[1], key[2], key[3], key[4]])
            if tc_val > 0 or fc_val > 0:
                # Instrument-calibrated simulation: simulate passive gating bursts through
                # dead time, filter cutoff, and burst criterion.
                # Under a passive single-channel model with mean open duration tau_o and brief intra-burst
                # closures tau_c, closures > t_crit terminate the burst.
                # Missed intra-burst closures (< t_res) fuse adjacent openings.
                eff_dt = max(dt_val, 0.36 / fc_val if fc_val > 0 else 0.0)
                tau_o = max(eff_dt * 2.0, mean_duration if mean_duration > 0 else 0.010)
                tau_c = max(eff_dt * 0.5, 0.001)

                p_end = math.exp(-tc_val / tau_c) if tc_val > 0 else 1.0
                p_end = max(1e-4, min(1.0, p_end))

                p_miss = 1.0 - math.exp(-eff_dt / tau_c) if eff_dt > 0 else 0.0
                p_term_eff = p_end / max(1e-6, 1.0 - (1.0 - p_end) * p_miss)
                p_term_eff = max(1e-4, min(1.0, p_term_eff))

                tau_o_eff = tau_o + (1.0 - p_end) * p_miss * tau_c
                buf = int(key[0] * 1.5) + 20
                k_ops = rng.geometric(p_term_eff, size=(REPS, buf))
                bursts = rng.gamma(k_ops, max(1e-6, tau_o_eff))

                valid = np.maximum(bursts - eff_dt, 0.0) if eff_dt > 0 else bursts
                x = valid[:, :key[0]]
                if key[1] > 0:
                    q = key[1] / 1e4
                    x = np.floor(x / q) * q + q / 2
                self.cache[key] = cv2(x)
            else:
                x = rng.exponential(1.0, (REPS, key[0]))
                if key[1] > 0:
                    q = key[1] / 1e4
                    x = np.floor(x / q) * q + q / 2
                self.cache[key] = cv2(x)
        return self.cache[key]


def adjust(x, dead_time, quantum):
    """Plan step 2: discard bursts shorter than the resolution and subtract it. With a frame grid the resolution
    is rounded up to the grid, one further frame is discarded (the first frame above the resolution is only
    partly observed, whatever the rounding convention), and half a frame is added to match the quantised null."""
    if quantum > 0:
        snapped = np.round(x / quantum) * quantum
        if not np.allclose(x, snapped, rtol=0, atol=0.01 * quantum):
            raise ValueError("durations are not on the stated frame grid; omit --quantum or correct it")
        x = snapped
        d = quantum * (math.ceil(dead_time / quantum - 1e-9) + 1)
        return x[x >= d - 1e-9 * quantum] - d + quantum / 2
    return x[x >= dead_time] - dead_time


def patch_stats(x, quantum, null, dead_time=0.0, filter_cutoff=None, burst_criterion=None):
    n = int(x.size)
    q_over_m = quantum / x.mean() if quantum > 0 else 0.0
    d = null.draws(
        n,
        q_over_m,
        dead_time=dead_time,
        filter_cutoff=filter_cutoff,
        burst_criterion=burst_criterion,
        mean_duration=float(x.mean()),
    )
    stat = float(cv2(x))
    return {"n": n, "cv2": stat, "crit": float(np.quantile(d, ALPHA)), "p": float((1 + (d <= stat).sum()) / (1 + REPS)),
            "z": float((stat - d.mean()) / d.std()), "null_z": (d - d.mean()) / d.std()}


def combine(stats, seed_offset):
    """Plan step 5: sqrt(N)-weighted sum of standardised scores, calibrated jointly with the actual burst counts."""
    ns = np.array([s["n"] for s in stats], dtype=float)
    w = np.sqrt(ns)
    rng = np.random.default_rng([SEED, 99, seed_offset])
    null = sum(wi * rng.permutation(s["null_z"]) for wi, s in zip(w, stats)) / math.sqrt(ns.sum())
    comb = float(sum(wi * s["z"] for wi, s in zip(w, stats)) / math.sqrt(ns.sum()))
    crit = float(np.quantile(null, ALPHA))
    return {"statistic": round(comb, 4), "critical_alpha_0.01": round(crit, 4),
            "p_one_sided": round(float((1 + (null <= comb).sum()) / (1 + REPS)), 5), "rejects_mixture": bool(comb < crit)}


def analyse(patches, dead_time, quantum, null, min_bursts, only=None, filter_cutoff=None, burst_criterion=None):
    rows, stats = [], {}
    for pid in sorted(patches):
        raw = np.asarray(patches[pid], dtype=float)
        x = adjust(raw, dead_time, quantum)
        n = int(x.size)
        if n < min_bursts or (only is not None and pid not in only):
            rows.append({"patch_id": pid, "n_bursts_used": n, "included": False,
                         "reason": f"fewer than {min_bursts} bursts after the dead-time rule" if n < min_bursts else "not included at every resolution"})
            continue
        s = patch_stats(x, quantum, null, dead_time=dead_time, filter_cutoff=filter_cutoff, burst_criterion=burst_criterion)
        stats[pid] = s
        kept = raw[raw >= dead_time]
        gap = float(kept.min() - dead_time) if kept.size else float("nan")
        expected = float(x.mean() / n)        # under an exponential the shortest adjusted burst is about mean / N
        rows.append({"patch_id": pid, "n_bursts_used": n, "included": True, "cv2": round(s["cv2"], 4), "critical_cv2_alpha_0.01": round(s["crit"], 4),
                     "p_one_sided": round(s["p"], 5), "rejects_mixture": bool(s["cv2"] < s["crit"]),
                     "shortest_burst_minus_resolution_s": round(gap, 6),
                     # chance of a gap at least this large if the resolution is right (about exp(-gap * N / mean));
                     # the flag has a false-flag rate of about 5% per patch and is not computed on a frame grid
                     "p_gap_this_large_if_resolution_correct": (round(float(math.exp(-gap / expected)), 5) if quantum == 0 else None),
                     "resolution_may_be_too_small": bool(quantum == 0 and gap > 3 * expected)})
    out = {
        "resolution_s": dead_time,
        "filter_cutoff_hz": filter_cutoff,
        "burst_criterion_s": burst_criterion,
        "frame_quantum_s": quantum,
        "min_bursts": min_bursts,
        "patches": rows,
        "n_patches_included": len(stats),
    }
    if stats:
        out["combined"] = combine([stats[k] for k in sorted(stats)], len(stats))
    return out, set(stats)


def run(table, dead_time=None, quantum=0.0, control=None, min_bursts=MIN_BURSTS, min_bursts_control=MIN_BURSTS_CONTROL,
        filter_cutoff=None, burst_criterion=None):
    if dead_time is None:
        if filter_cutoff is not None and filter_cutoff > 0:
            dead_time = 0.36 / float(filter_cutoff)
        else:
            raise ValueError("dead time must be provided, or derived via filter_cutoff (t_res = 0.36 / f_c)")
    dead_time = float(dead_time)
    quantum = float(quantum)
    filter_cutoff = float(filter_cutoff) if filter_cutoff is not None else None
    burst_criterion = float(burst_criterion) if burst_criterion is not None else None

    if dead_time < 0 or quantum < 0 or min(min_bursts, min_bursts_control) < 2:
        raise ValueError("dead time and quantum must be non-negative, and burst minimums at least 2")
    if filter_cutoff is not None and filter_cutoff <= 0:
        raise ValueError("filter_cutoff must be positive")
    if burst_criterion is not None and burst_criterion < 0:
        raise ValueError("burst_criterion must be non-negative")

    null = Null()
    res = {
        "plan": "specs/CFTR_DATA_REQUEST_AND_ANALYSIS_PLAN.md",
        "alpha_one_sided": ALPHA,
        "monte_carlo_draws": REPS,
        "seed": SEED,
        "filter_cutoff_hz": filter_cutoff,
        "burst_criterion_s": burst_criterion,
        "implemented": [
            "step 1: burst minimums (50 WT, 5 control)",
            "step 2: dead-time adjustment and 0.36/f_c default filter resolution",
            "steps 3-5: CV^2 statistic, per-patch critical threshold, and sqrt(N) combination",
            "step 4 calibration: instrument-calibrated false-alarm threshold for dead time, filtering (f_c), and burst criterion (t_crit)",
            "section 5: 1x / 2x resolution sensitivity decision rule",
        ],
        "not_implemented": ["one-active-channel check (from metadata)", "control power simulation", "likelihood-ratio comparator (step 6)"],
        "constructs": {},
    }
    for construct, patches in sorted(table.items()):
        mb = min_bursts_control if control is not None and construct == control else min_bursts
        fc2 = (filter_cutoff / 2.0) if filter_cutoff is not None else None
        _, inc1 = analyse(patches, dead_time, quantum, null, mb, filter_cutoff=filter_cutoff, burst_criterion=burst_criterion)
        _, inc2 = analyse(patches, 2 * dead_time, quantum, null, mb, filter_cutoff=fc2, burst_criterion=burst_criterion)
        both = inc1 & inc2
        main, _ = analyse(patches, dead_time, quantum, null, mb, filter_cutoff=filter_cutoff, burst_criterion=burst_criterion)
        main["role"] = "control" if control is not None and construct == control else "test"
        a, _ = analyse(patches, dead_time, quantum, null, mb, only=both, filter_cutoff=filter_cutoff, burst_criterion=burst_criterion)
        b, _ = analyse(patches, 2 * dead_time, quantum, null, mb, only=both, filter_cutoff=fc2, burst_criterion=burst_criterion)
        r1, r2 = a.get("combined", {}).get("rejects_mixture"), b.get("combined", {}).get("rejects_mixture")
        main["resolution_sensitivity"] = {
            "patches_included_at_both_resolutions": len(both),
            "combined_at_1x": a.get("combined"), "combined_at_2x": b.get("combined"),
            "rejects_at_1x_not_2x": bool(r1 and r2 is False), "rejects_at_2x_not_1x": bool(r2 and r1 is False)}
        # the decision: reject only if the combined test rejects at both resolutions on the patches included at both
        main["decision"] = ("no_patches" if r1 is None or r2 is None else "rejects_mixture" if (r1 and r2)
                            else "resolution_sensitive_inconclusive" if (r1 or r2) else "does_not_reject")
        res["constructs"][construct] = main
    return res


def read_csv(path):
    table, dropped = defaultdict(lambda: defaultdict(list)), 0
    with open(path, newline="", encoding="utf-8-sig") as f:
        for row in csv.DictReader(f):
            try:
                v = float(row["duration"])
            except (TypeError, ValueError, KeyError):
                dropped += 1
                continue
            if not math.isfinite(v) or v <= 0:
                dropped += 1
                continue
            table[(row.get("construct") or "all").strip()][row["patch_id"].strip()].append(v)
    return table, dropped


def self_test():
    rng = np.random.default_rng(7)
    td = 0.002
    sizes = (120, 250, 80, 400, 150, 60)

    def patches(draw):
        return {f"p{i}": draw(n) + td for i, n in enumerate(sizes)}

    table = {
        "two_step": patches(lambda n: rng.exponential(1 / 12, n) + rng.exponential(1 / 8, n)),
        "single_exponential": patches(lambda n: rng.exponential(1 / 5, n)),
        "mixture": patches(lambda n: np.where(rng.random(n) < 0.5, rng.exponential(1 / 12, n), rng.exponential(1 / 3, n))),
        "control": {f"c{i}": rng.exponential(2.0, n) + td for i, n in enumerate((6, 9, 14))},
    }
    res = run(table, td, 0.0, control="control")
    got = {k: v["decision"] for k, v in res["constructs"].items()}
    ok = got == {"two_step": "rejects_mixture", "single_exponential": "does_not_reject", "mixture": "does_not_reject", "control": "does_not_reject"}
    ok = ok and res["constructs"]["control"]["n_patches_included"] == 3 and res["constructs"]["two_step"]["n_patches_included"] == len(sizes)
    json.dumps(res, allow_nan=False)

    # Test 0.36 / f_c default derivation
    res_fc = run(table, filter_cutoff=1000, control="control")
    ok = ok and math.isclose(res_fc["constructs"]["two_step"]["resolution_s"], 0.00036, rel_tol=1e-5)

    null, trials, rej, rej_q, rej_n = Null(11), 2000, 0, 0, 0
    q, td_q = 0.1, 0.17                      # frame grid, off-grid true dead time: events below it are never observed
    for _ in range(trials):
        ps = {f"p{i}": rng.exponential(1.0, n) for i, n in enumerate(sizes)}
        st = [patch_stats(v, 0.0, null) for v in ps.values()]
        rej += combine(st, len(st))["rejects_mixture"]
        t_true = rng.exponential(1.0, 700)
        t_true = t_true[t_true >= td_q][:300]
        for kind in ("up", "nearest"):
            y = (np.ceil(t_true / q) if kind == "up" else np.round(t_true / q)) * q
            s = patch_stats(adjust(y, td_q, q), q, null)
            if kind == "up":
                rej_q += s["cv2"] < s["crit"]
            else:
                rej_n += s["cv2"] < s["crit"]
    rate, rate_q, rate_n = rej / trials, rej_q / trials, rej_n / trials
    ok = ok and rate <= 0.02 and rate_q <= 0.02 and rate_n <= 0.02

    # Verify that instrument-calibrated threshold prevents false alarm on filtered bursts with burst criterion
    res_cal = run(table, dead_time=0.002, filter_cutoff=1000, burst_criterion=0.005, control="control")
    got_cal = {k: v["decision"] for k, v in res_cal["constructs"].items()}
    ok = ok and got_cal["single_exponential"] == "does_not_reject" and got_cal["mixture"] == "does_not_reject"

    for bad in (dict(dead_time=-0.1, quantum=0.0), dict(dead_time=0.0, quantum=-1.0), dict(filter_cutoff=-10.0)):
        try:
            run(table, **bad)
            ok = False
        except ValueError:
            pass
    try:
        adjust(np.array([0.11, 0.23]), 0.0, 0.1)
        ok = False
    except ValueError:
        pass
    print("self-test decisions", got, "| false-alarm rate, combined:", rate, "| frame grid, missed short events, rounded up / nearest:", rate_q, rate_n,
          f"({trials} datasets each)", "| OK" if ok else "| FAIL")
    return 0 if ok else 1


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("csv", nargs="?")
    ap.add_argument("--dead-time", type=float, default=None, help="resolution imposed on the record, seconds; required unless --filter-cutoff is given")
    ap.add_argument("--filter-cutoff", type=float, default=None, help="filter -3 dB cutoff frequency f_c in Hz (sets default dead_time = 0.36 / f_c)")
    ap.add_argument("--burst-criterion", type=float, default=None, help="critical closed time t_crit in seconds defining burst boundaries")
    ap.add_argument("--quantum", type=float, default=0.0, help="frame interval in seconds, if durations lie on a grid")
    ap.add_argument("--control", default=None, help="name of the control construct (uses --min-bursts-control)")
    ap.add_argument("--min-bursts", type=int, default=MIN_BURSTS)
    ap.add_argument("--min-bursts-control", type=int, default=MIN_BURSTS_CONTROL)
    ap.add_argument("--out")
    ap.add_argument("--self-test", action="store_true")
    a = ap.parse_args(argv)
    if a.self_test:
        return self_test()
    if not a.csv or (a.dead_time is None and a.filter_cutoff is None):
        ap.error("give a CSV and --dead-time (or --filter-cutoff), or --self-test")
    table, dropped = read_csv(a.csv)
    res = run(table, a.dead_time, a.quantum, a.control, a.min_bursts, a.min_bursts_control,
              filter_cutoff=a.filter_cutoff, burst_criterion=a.burst_criterion)
    res["rows_dropped_as_missing_or_non_positive"] = dropped
    text = json.dumps(res, indent=1, allow_nan=False)
    if a.out:
        open(a.out, "w", encoding="utf-8").write(text)
    else:
        print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
