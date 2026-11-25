#!/usr/bin/env python3
"""Full-pipeline validation on cached F1-ATPase 120-degree step dwells (Zenodo 5567245).

Source dataset: benchmark/out/zenodo5567245_f1_atpase.npz (Kobayashi, Ueno, Li et al., CC0 1.0).
Scope note: rotary F1-ATPase 120-degree step dwells, NOT CFTR. No CFTR recording is analyzed.

Usage:
    python benchmark/whole_sample/run_f1_atpase_validation.py
    python benchmark/whole_sample/run_f1_atpase_validation.py --check
"""
import argparse
import hashlib
import json
import os
import sys
import time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BENCH = os.path.dirname(HERE)
NPZ_PATH = os.path.join(BENCH, "out", "zenodo5567245_f1_atpase.npz")
OUT_PATH = os.path.join(HERE, "f1_atpase_full_pipeline.json")

if HERE not in sys.path:
    sys.path.insert(0, HERE)

import run_on_patches  # noqa: E402


CONDITIONS = {
    "dwell_2d": {
        "label": "Fig2D wild type, 3 mM ATP (saturating; hydrolysis-limited dwell)",
        "role": "test",
    },
    "dwell_3a": {
        "label": "Fig3A wild type, 1 mM ATP-gamma-S (slowly hydrolysed analog)",
        "role": "test",
    },
    "dwell_2c": {
        "label": "Fig2C wild type, 300 nM ATP (binding-limited dwell, reference control)",
        "role": "control",
    },
}


def _condition_grid(arr):
    u = np.unique(arr)
    q = float(np.min(np.diff(u)))
    xmin = float(arr.min())
    # Matching cv2_burst_test.py: adj = x - xmin + q/2 corresponds to adjust(x, xmin - q, q)
    # because adjust rounds dead_time up to the grid and discards one additional partial frame.
    dead_time = max(0.0, xmin - q)
    if abs(dead_time) < 1e-12:
        dead_time = 0.0
    return q, xmin, dead_time


def build_validation_payload(bootstrap_reps=run_on_patches.BOOTSTRAP_REPS):
    if not os.path.isfile(NPZ_PATH):
        raise FileNotFoundError(f"Cached F1-ATPase dataset not found: {NPZ_PATH}")
    with open(NPZ_PATH, "rb") as f:
        cache_sha256 = hashlib.sha256(f.read()).hexdigest()

    z = np.load(NPZ_PATH, allow_pickle=False)
    arrays = {}
    grids = {}
    for key in CONDITIONS:
        arr = np.asarray(z[key], dtype=float)
        q, xmin, dt = _condition_grid(arr)
        arrays[key] = arr
        adj = run_on_patches.adjust(arr, dt, q)
        grids[key] = {
            "condition": CONDITIONS[key]["label"],
            "role": CONDITIONS[key]["role"],
            "n_raw": int(arr.size),
            "n_after_adjust": int(adj.size),
            "frame_quantum_s": round(q, 9),
            "shortest_dwell_s": round(xmin, 9),
            "effective_dead_time_arg_s": round(dt, 9),
            "mean_adjusted_s": round(float(adj.mean()), 6),
        }

    # 1. Clean per-condition full-pipeline runs (each condition has its own frame quantum and dead time)
    per_condition = {}
    for key in ("dwell_2d", "dwell_3a", "dwell_2c"):
        q, _, dt = _condition_grid(arrays[key])
        ctrl = key if CONDITIONS[key]["role"] == "control" else None
        per_condition[key] = run_on_patches.run(
            {key: {key: arrays[key]}},
            dead_time=dt,
            quantum=q,
            control=ctrl,
            bootstrap_reps=bootstrap_reps,
        )

    # 2. Paired runs treating 300 nM ATP (dwell_2c) as the binding-limited reference control
    paired_runs = {}
    q_2c, _, dt_2c = _condition_grid(arrays["dwell_2c"])
    for test_key in ("dwell_2d", "dwell_3a"):
        q_t, _, dt_t = _condition_grid(arrays[test_key])
        pair_name = f"{test_key}_vs_dwell_2c_control"
        table = {
            test_key: {test_key: arrays[test_key]},
            "dwell_2c": {"dwell_2c": arrays["dwell_2c"]},
        }
        metadata = {
            test_key: {test_key: {"dead_time": dt_t, "quantum": q_t}},
            "dwell_2c": {"dwell_2c": {"dead_time": dt_2c, "quantum": q_2c}},
        }
        paired_runs[pair_name] = run_on_patches.run(
            table,
            dead_time=None,
            quantum=0.0,
            control="dwell_2c",
            metadata=metadata,
            bootstrap_reps=bootstrap_reps,
        )

    return {
        "schema_version": "1.0",
        "dataset": {
            "source": "Zenodo 5567245 (Kobayashi, Ueno, Li et al., CC0 1.0)",
            "cache_path": "benchmark/out/zenodo5567245_f1_atpase.npz",
            "cache_sha256": cache_sha256,
            "system_note": "rotary F1-ATPase 120-degree step dwells, NOT CFTR",
            "cftr_boundary": "No CFTR recording is analyzed; every CFTR result in this package is a simulation.",
        },
        "dead_time_alignment_note": (
            "cv2_burst_test.py sets adj = x - xmin + q/2. In run_on_patches.adjust(x, dead_time, q), "
            "dead_time = max(0, xmin - q) reproduces the exact same adjusted dwells and burst counts "
            "(1599 for dwell_2d where xmin=2*q, 192 for dwell_3a where xmin=q, 443 for dwell_2c where xmin=q)."
        ),
        "bootstrap_reps": int(bootstrap_reps),
        "conditions_metadata": grids,
        "per_condition_runs": per_condition,
        "paired_control_runs": paired_runs,
    }


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bootstrap-reps", type=int, default=run_on_patches.BOOTSTRAP_REPS)
    ap.add_argument("--out", default=OUT_PATH)
    ap.add_argument("--check", action="store_true", help="verify committed JSON matches deterministic recomputation")
    args = ap.parse_args(argv)

    t0 = time.perf_counter()
    payload = build_validation_payload(bootstrap_reps=args.bootstrap_reps)
    elapsed = time.perf_counter() - t0
    rendered = json.dumps(payload, indent=1, allow_nan=False) + "\n"

    if args.check:
        if not os.path.isfile(args.out):
            print(f"FAIL: missing output file {args.out}", file=sys.stderr)
            return 1
        with open(args.out, encoding="utf-8") as f:
            existing = f.read()
        if existing == rendered:
            print(f"PASS: {os.path.relpath(args.out, BENCH)} matches deterministic recomputation ({elapsed:.2f}s)")
            return 0
        try:
            if json.loads(existing) == payload:
                print(f"PASS (JSON equal): {os.path.relpath(args.out, BENCH)} ({elapsed:.2f}s)")
                return 0
        except (json.JSONDecodeError, TypeError, ValueError):
            pass
        print(f"FAIL: {args.out} differs from recomputed output", file=sys.stderr)
        return 1

    with open(args.out, "w", encoding="utf-8", newline="\n") as f:
        f.write(rendered)
    print(f"Wrote {args.out} in {elapsed:.2f}s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
