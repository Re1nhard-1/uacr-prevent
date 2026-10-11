import math
import numpy as np
from scipy.special import expit, logit
from scipy.stats import t as tdist, norm, beta as betad

def cr1_var(y, cl):
    y = np.asarray(y, float)
    n = len(y)
    m = y.mean()
    ids, inv = np.unique(cl, return_inverse=True)
    G = len(ids)
    tot = np.bincount(inv, weights=y - m, minlength=G)
    return (G / (G - 1) * np.sum(tot ** 2) / n ** 2, G)

def kg_interval(x, n, deff, df):
    neff = n / max(deff, 1.0)
    neff *= (norm.ppf(0.975) / tdist.ppf(0.975, df)) ** 2
    xe = x / n * neff
    lo = 0.0 if xe <= 0 else betad.ppf(0.025, xe, neff - xe + 1)
    hi = 1.0 if xe >= neff else betad.ppf(0.975, xe + 1, neff - xe)
    return (float(lo), float(hi))

def wilson(x, n, z=1.959963984540054):
    p = x / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (c - h, c + h)

def prop_ci(y, cl):
    y = np.asarray(y, float)
    n = len(y)
    x = float(y.sum())
    th = x / n
    V, G = cr1_var(y, cl)
    df = G - 1
    base = th * (1 - th) / n if 0 < th < 1 else 0.0
    deff = V / base if base > 0 else float('nan')
    Vs = max(V, base)
    out = {'n': n, 'events': int(round(x)), 'estimate': th, 'deff': deff, 'G': G, 'df': df, 'wilson': wilson(x, n)}
    if x < 10 or x > n - 10 or (not 0 < th < 1):
        lo, hi = kg_interval(x, n, deff if deff == deff else 1.0, df)
        out.update(lower=lo, upper=hi, method='Korn-Graubard (CR1 deff)')
    else:
        se = math.sqrt(Vs) / (th * (1 - th))
        q = tdist.ppf(0.975, df)
        out.update(lower=float(expit(logit(th) - q * se)), upper=float(expit(logit(th) + q * se)), method='CR1 logit-t, deff floored at 1')
    return out

def paired_ci(d, cl):
    d = np.asarray(d, float)
    n = len(d)
    m = d.mean()
    V, G = cr1_var(d, cl)
    Vs = max(V, d.var(ddof=1) / n)
    q = tdist.ppf(0.975, G - 1)
    return {'n': n, 'estimate': m, 'lower': m - q * math.sqrt(Vs), 'upper': m + q * math.sqrt(Vs), 'G': G, 'method': 'CR1 floored, Wald t(G-1)'}

def tango(b, c, n, z=1.959963984540054):

    def T(dl):
        A = 2 * n
        B = -b - c + (2 * n - b + c) * dl
        C = -c * dl * (1 - dl)
        q = (math.sqrt(max(B * B - 4 * A * C, 0)) - B) / (2 * A)
        den = n * (2 * q + dl * (1 - dl))
        return (b - c - n * dl) / math.sqrt(den) if den > 0 else math.inf if b - c - n * dl > 0 else -math.inf
    est = (b - c) / n

    def root(lo, hi, target):
        for _ in range(200):
            mid = (lo + hi) / 2
            if (T(mid) - target) * (T(lo) - target) <= 0:
                hi = mid
            else:
                lo = mid
        return (lo + hi) / 2
    lo = root(-1 + 1e-12, est, z) if T(-1 + 1e-12) > z else -1.0
    hi = root(est, 1 - 1e-12, -z) if T(1 - 1e-12) < -z else 1.0
    return (lo, hi)

def rcs_basis(x, knots):
    x = np.asarray(x, float)
    k1, k2, k3 = knots
    p = lambda u: np.maximum(u, 0) ** 3
    nl = (p(x - k1) - p(x - k2) * (k3 - k1) / (k3 - k2) + p(x - k3) * (k2 - k1) / (k3 - k2)) / (k3 - k1) ** 2
    return np.column_stack([x, nl])

def logit_fit(X, y, w=None, iters=50):
    X = np.asarray(X, float)
    y = np.asarray(y, float)
    w = np.ones(len(y)) if w is None else np.asarray(w, float)
    b = np.zeros(X.shape[1])
    for _ in range(iters):
        p = expit(X @ b)
        W = w * p * (1 - p)
        H = X.T @ (X * W[:, None]) + 1e-10 * np.eye(X.shape[1])
        g = X.T @ (w * (y - p))
        step = np.linalg.solve(H, g)
        b += step
        if np.max(np.abs(step)) < 1e-10:
            break
    return b

def icc_a1(x1, x2):
    x = np.column_stack([x1, x2])
    n, k = x.shape
    gm = x.mean()
    msr = k * np.sum((x.mean(1) - gm) ** 2) / (n - 1)
    msc = n * np.sum((x.mean(0) - gm) ** 2) / (k - 1)
    sse = np.sum((x - x.mean(1, keepdims=True) - x.mean(0, keepdims=True) + gm) ** 2)
    mse = sse / ((n - 1) * (k - 1))
    return float((msr - mse) / (msr + (k - 1) * mse + k * (msc - mse) / n))

def boot_weights(cl, B, seed):
    ids, inv = np.unique(cl, return_inverse=True)
    G = len(ids)
    rng = np.random.default_rng(seed)
    counts = rng.multinomial(G, np.full(G, 1 / G), size=B)
    return counts[:, inv].astype(float)

def pct_ci(v):
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    return [float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))]

def batch_mcse(v, nb=4):
    v = np.asarray(v, float)
    parts = np.array_split(v, nb)
    lo = [np.percentile(p, 2.5) for p in parts]
    hi = [np.percentile(p, 97.5) for p in parts]
    return (float(np.std(lo, ddof=1) / math.sqrt(nb)), float(np.std(hi, ddof=1) / math.sqrt(nb)))

def strat_taylor(y, strata, psu):
    y = np.asarray(y, float)
    n = len(y)
    th = y.mean()
    z = (y - th) / n
    var = 0.0
    npsu = 0
    nst = 0
    for h in np.unique(strata):
        m = strata == h
        tots = np.array([z[m & (psu == p)].sum() for p in np.unique(psu[m])])
        if len(tots) < 2:
            continue
        var += len(tots) / (len(tots) - 1) * np.sum((tots - tots.mean()) ** 2)
        npsu += len(tots)
        nst += 1
    return (th, var, npsu - nst)

def nchs_flag(x, n, lo, hi, df):
    width = hi - lo
    flags = []
    if n < 30:
        flags.append('n<30')
    if df < 8:
        flags.append('df<8')
    if width >= 0.3:
        flags.append('absolute CI width >=0.30')
    p = x / n if n else 0
    if width > 0.05 and p > 0 and (width / p > 1.3):
        flags.append('relative CI width >130%')
    return flags
