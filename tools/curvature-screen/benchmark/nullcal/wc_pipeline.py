"""Verbatim copies of the estimator functions from the package script
benchmark/macroscopic_patch_clamp_benchmark.py
(SHA-256 088fa230283829078ec06ec63371360465fe06068fc903f80ce6a4f3a9113ea9).

The package script has no __main__ guard (importing it would rerun everything and rewrite its JSON), so the
functions below are copied text-for-text; nullcal.py checks at run time, via the AST, that each copied function
is identical to the one in the package file (read-only).  Module globals used by the copied functions
(rng, g, dt, lag_idx, K) carry the package's values; rng is re-seeded per work chunk by set_seed().
"""
import math
import numpy as np
from scipy.linalg import expm

rng = np.random.default_rng(20261113)
g = np.array([0.0, 1.0, 1.0])
dt = 0.002; lag_idx = (25, 100, 175)          # 50, 200, 350 ms
K = 16


def set_seed(seed):
    global rng
    rng = np.random.default_rng(seed)


# ---------------- verbatim copies below ----------------
def gen_from_rates(k):
    """k: dict (i,j)->rate for i != j; returns generator Q (rows sum to zero, Q[i,j] = rate i->j)."""
    m = 1 + max(max(i, j) for i, j in k)
    Q = np.zeros((m, m))
    for (i, j), r in k.items():
        Q[i, j] = r
    Q[np.diag_indices(m)] = -Q.sum(1)
    return Q


def stationary(Q):
    w, v = np.linalg.eig(Q.T)
    p = np.real(v[:, np.argmin(abs(w))])
    return p / p.sum()


def autocov(Q, g, taus):
    p = stationary(Q); mu = p @ g; d = g - mu
    return np.array([(p * d) @ expm(Q * t) @ d for t in taus])


def spectral(Q, g):
    """Return eigenvalues (nonzero), and the weights a_k of C(tau) = sum a_k exp(lambda_k tau) (complex)."""
    w, V = np.linalg.eig(Q)
    Vi = np.linalg.inv(V)
    p = stationary(Q); mu = p @ g; d = g - mu
    a = (p * d) @ V * (Vi @ d)
    order = np.argsort(abs(w))
    keep = [i for i in order if abs(w[i]) > 1e-9]
    return w[keep], a[keep]


def random_reversible(m):
    pi = rng.dirichlet(np.ones(m))
    S = rng.uniform(0.2, 5.0, (m, m)); S = (S + S.T) / 2
    Q = np.zeros((m, m))
    for i in range(m):
        for j in range(m):
            if i != j:
                Q[i, j] = S[i, j] / pi[i] ** 0.5 * pi[j] ** 0.5       # detailed balance: pi_i Q_ij = pi_j Q_ji
    Q[np.diag_indices(m)] = -Q.sum(1)
    return Q, pi


def circulant_sampler(Cfun_vals, L):
    """Gaussian stationary process of length L with covariance Cfun_vals[0..L-1] (circulant embedding)."""
    lam = np.fft.rfft(np.concatenate([Cfun_vals[:L // 2 + 1], Cfun_vals[1:L // 2][::-1]]))
    lam = np.maximum(lam.real, 0)
    def draw():
        z = rng.standard_normal(L)
        return np.fft.irfft(np.fft.rfft(z) * np.sqrt(lam), n=L) * 1.0
    return draw


def simulate_z(Q, n, sigma_w, T, reps):
    L = int(T / K / dt)
    L -= L % 2
    taus = np.arange(L // 2 + 1) * dt
    ev_, a_ = spectral(Q, g)
    Cs = np.real(np.exp(np.outer(taus, ev_)) @ a_) * n
    Cs_full = Cs.copy()
    Cs_full[0] += sigma_w ** 2                             # white measurement noise at lag 0 only
    draw = circulant_sampler(Cs_full, L)
    zs = []
    for _ in range(reps):
        blk = np.empty((K, 3))
        for k in range(K):
            x = draw(); x -= x.mean()
            f = np.fft.rfft(x, n=2 * L)
            ac = np.fft.irfft(f * np.conj(f))[:max(lag_idx) + 1] / L
            blk[k] = ac[list(lag_idx)]
        vals = []
        for a in range(K):
            for b in range(a + 1, K):
                vals.append(0.5 * (blk[a, 0] * blk[b, 2] + blk[b, 0] * blk[a, 2]) - blk[a, 1] * blk[b, 1])
        u = np.mean(vals)
        jk = []
        for j in range(K):
            idx = [i for i in range(K) if i != j]
            s = blk[idx]
            v = [0.5 * (s[a, 0] * s[b, 2] + s[b, 0] * s[a, 2]) - s[a, 1] * s[b, 1] for a in range(K - 1) for b in range(a + 1, K - 1)]
            jk.append(np.mean(v))
        jk = np.array(jk)
        se = math.sqrt((K - 1) / K * np.sum((jk - jk.mean()) ** 2))
        zs.append(u / se)
    return np.array(zs)
