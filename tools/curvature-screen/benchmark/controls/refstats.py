"""Verbatim copies of the estimator functions from
benchmark/validate_fluorescence_benchmark.py
(lines 109-111 constants, 96-97 r(), 122-164 wstats, 193-212 wdiff_stats, 439-483 _h3_diff_split_stats).
The module cannot be imported directly: it calls sys.exit() at import time when the Zenodo 14848839 raw data
are absent, and the raw data are not in the package. No statistics were rewritten; only copied.
"""
import hashlib
import numpy as np

import os
PKG = os.path.dirname(os.path.dirname(os.path.abspath(__file__))).replace('\\', '/')
REF_SCRIPT = PKG + '/validate_fluorescence_benchmark.py'
HELLENKAMP_NPZ = PKG + '/out/hellenkamp2018_zen1249497_histograms.npz'
TISSUE_NPZ = PKG + '/out/zenodo7656540_leica_stellaris8_tissue_flim.npz'
LEDGER_JSON = PKG + '/out/validated_benchmark_suite.json'

BG = 300
K3 = np.array([[0.0, 0.0, 0.5], [0.0, -1.0, 0.0], [0.5, 0.0, 0.0]])
Q6 = np.block([[np.zeros((3, 3)), K3], [K3, np.zeros((3, 3))]])


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def r(x, n=4):
    return None if x is None or not np.isfinite(x) else round(float(x), n)


def wstats(hA, hB, w, bA=None, bB=None, B_eff=BG, shared=False):
    """All estimators + plug-in and exact variances. Returns dict (floats)."""
    M = w[0].stop - w[0].start
    n = len(hA)
    truncated = any(s.stop > n for s in w)
    if bA is None:
        bA = hA[:BG].mean()
    if bB is None:
        bB = hB[:BG].mean()
    YA = np.array([hA[s].sum() for s in w], float)
    YB = np.array([hB[s].sum() for s in w], float)
    yA, yB = YA - M * bA, YB - M * bB
    cA, cB = (M ** 2 / B_eff) * max(bA, 0.0), (M ** 2 / B_eff) * max(bB, 0.0)
    SA = np.diag(np.maximum(YA, 0.0)) + cA * np.ones((3, 3))
    SB = np.diag(np.maximum(YB, 0.0)) + cB * np.ones((3, 3))
    C = (cA if shared else 0.0) * np.ones((3, 3))
    S = np.block([[SA, C], [C.T, SB]])
    z = np.concatenate([yA, yB])
    Ws = float(yA @ K3 @ yB)
    QSQ = Q6 @ S @ Q6
    v_plug = float(z @ QSQ @ z)
    tr_s = 0.5 * float(np.trace(QSQ @ S))
    v_exact = max(v_plug - tr_s, tr_s)
    se_plug = np.sqrt(max(v_plug, 1e-12))
    se_exact = np.sqrt(max(v_exact, 1e-12))
    yp = 0.5 * (yA + yB)
    Sp = 0.25 * (SA + SB + C + C.T)
    Wu = float(yp @ K3 @ yp)
    Wu_bc = Wu - float(np.trace(K3 @ Sp))
    vu_plug = 4.0 * float(yp @ K3 @ Sp @ K3 @ yp)
    tr_u = 2.0 * float(np.trace(K3 @ Sp @ K3 @ Sp))
    vu_exact = max(vu_plug - tr_u, tr_u)
    se_u = np.sqrt(max(vu_plug, 1e-12))
    se_u_exact = np.sqrt(max(vu_exact, 1e-12))
    Wun = float(yA[0] * yB[2] - yA[1] * yB[1])
    L = float(np.log(yp[0]) - 2 * np.log(yp[1]) + np.log(yp[2])) if np.all(yp > 0) else float('nan')
    evaluable = bool(np.all(yp > 2.0 * np.sqrt(np.diag(Sp)))) and not truncated
    y2sq = yp[1] ** 2 if yp[1] != 0 else float('nan')
    return dict(Ws=Ws, Wu=Wu, Wu_bc=Wu_bc, Wun=Wun, L=L,
                se_plug=se_plug, se_exact=se_exact, se_u=se_u, se_u_exact=se_u_exact,
                zs_plug=Ws / se_plug, zs_exact=Ws / se_exact, zu=Wu / se_u, zu_bc=Wu_bc / se_u_exact,
                W_rel=Ws / y2sq, truncated=truncated, evaluable=evaluable, yp=yp.tolist())


def wdiff_stats(hA, hB, w4):
    YA = np.array([hA[s].sum() for s in w4], float)
    YB = np.array([hB[s].sum() for s in w4], float)
    dA, dB = YA[:-1] - YA[1:], YB[:-1] - YB[1:]

    def tri(Y):
        return np.array([[Y[0] + Y[1], -Y[1], 0.0], [-Y[1], Y[1] + Y[2], -Y[2]], [0.0, -Y[2], Y[2] + Y[3]]])
    SA, SB = tri(YA), tri(YB)
    W = float(dA @ K3 @ dB)
    v_plug = float(dB @ K3 @ SA @ K3 @ dB + dA @ K3 @ SB @ K3 @ dA)
    tr_ab = float(np.trace(K3 @ SA @ K3 @ SB))
    v_exact = max(v_plug - tr_ab, tr_ab)
    se_plug, se_exact = np.sqrt(max(v_plug, 1e-12)), np.sqrt(max(v_exact, 1e-12))
    d_bar = 0.5 * (dA + dB)
    d2sq = d_bar[1] ** 2 if d_bar[1] != 0 else float('nan')
    return dict(W=W, W_rel=float(W / d2sq) if np.isfinite(d2sq) else None,
                se_plug=float(se_plug), se_exact=float(se_exact), z_plug=float(W / se_plug),
                z_exact=float(W / se_exact), floor_binding=bool(v_plug - tr_ab < tr_ab),
                truncated=any(s.stop > len(hA) for s in w4) or bool(np.min(YA) < 5 or np.min(YB) < 5),
                dA=dA.tolist(), dB=dB.tolist())
