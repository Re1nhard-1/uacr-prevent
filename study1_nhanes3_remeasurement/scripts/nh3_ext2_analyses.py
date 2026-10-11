"""Post hoc extensions, set 2 (2026-10-10), Study 1 NHANES III: E3 antihypertensive drug-initiation decision under the 2025 AHA/ACC"""
from pathlib import Path
import hashlib
import json
import math
import numpy as np
import pandas as pd
from scipy.special import expit, logit
from scipy.stats import norm, t as tdist
import nh3_build
import nh3_models as M
import nh3_stats as S
import nh3_testretest_stage12 as s12
from nh3_stage3 import fixed, lp_of, visit_inputs, SEED, B_BOOT, CUTS
ROOT = Path(__file__).resolve().parents[1]
ST = ROOT / 'studies/nhanes3_testretest'
R3 = json.loads((ST / 'analysis_v1/stage3/stage3_results.json').read_text(encoding='utf-8'))
R4 = json.loads((ST / 'analysis_v1/stage4/stage4_results.json').read_text(encoding='utf-8'))
LK = json.loads((ST / 'analysis_v1_posthoc/reclassification_lookup.json').read_text(encoding='utf-8'))
PROTO = ST / 'extensions2_protocol_20261010.md'
OUT = ST / 'analysis_v1_posthoc/ext2_results.json'
LT = logit(np.array(CUTS) / 100)

def dbp_reads():
    pos1 = {'PEP6G3': (1399, 1401), 'PEP6H3': (1409, 1411), 'PEP6I3': (1419, 1421)}
    pos2 = {'PER6G3': (159, 161), 'PER6H3': (169, 171), 'PER6I3': (179, 181)}
    return (s12.read('exam', pos1), s12.read('examse', pos2))

def dbp_value(rec, p, mode):
    v = {k: s12.num(rec[f'{p}6{k}3']) for k in 'GHI'}
    v = {k: x if x is not None and x > 0 else None for k, x in v.items()}
    if mode == 'all3':
        a = [x for x in v.values() if x is not None]
        return sum(a) / len(a) if a else None
    a = [x for x in (v['H'], v['I']) if x is not None]
    if a:
        return sum(a) / len(a)
    a = [x for x in v.values() if x is not None]
    return sum(a) / len(a) if a else None

def stage(sbp, dbp):
    return np.where((sbp >= 140) | (dbp >= 90), 2, np.where((sbp >= 130) | (dbp >= 80), 1, 0))

def decision(sbp, dbp, highrisk):
    st = stage(sbp, dbp)
    return np.where(st == 2, 3, np.where(st == 1, np.where(highrisk, 2, 1), 0))

def predict_pi(lp1, model):
    knots, beta, ind = (np.asarray(model['knots']), np.asarray(model['beta']), model['threshold_indicator'])
    dist = np.abs(lp1[:, None] - LT[None, :])
    d = dist.min(1)
    near = dist.argmin(1)
    X = np.column_stack([np.ones(len(d)), S.rcs_basis(d, knots)])
    if ind:
        X = np.column_stack([X, near == 1, near == 2]).astype(float)
    return expit(X @ beta)

def prop(y, cl):
    return S.prop_ci(np.asarray(y, float), cl)

def wild_ci(y, g, B=9999, seed=20261010):
    rng = np.random.default_rng(seed)
    ids, inv = np.unique(g, return_inverse=True)
    G = len(ids)
    ybar = y.mean()
    u = y - ybar
    V, _ = S.cr1_var(y, g)
    se = math.sqrt(V)
    webb = np.array([-math.sqrt(1.5), -1, -math.sqrt(0.5), math.sqrt(0.5), 1, math.sqrt(1.5)])
    ts = []
    for _ in range(B):
        wg = rng.choice(webb, G)
        ys = ybar + wg[inv] * u
        Vs, _ = S.cr1_var(ys, g)
        ts.append((ys.mean() - ybar) / math.sqrt(Vs))
    q = np.quantile(ts, [0.025, 0.975])
    return {'lower': float(ybar - q[1] * se), 'upper': float(ybar - q[0] * se), 'replicates': B, 'weights': 'Webb six-point'}

def e3(B, F, cl):
    ex, es = dbp_reads()
    out = {}

    def run(mode='r23', ckd_mode='egfr_or_uacr', dm_mode='visit', keep=None):
        sb1 = B.sbp1.to_numpy(float) if mode == 'r23' else B.sbp1_avail.to_numpy(float)
        sb2 = B.sbp2.to_numpy(float) if mode == 'r23' else B.sbp2_avail.to_numpy(float)
        db1 = np.array([dbp_value(ex[s], 'PEP', mode) if s in ex else None for s in B.SEQN], dtype=float)
        db2 = np.array([dbp_value(es[s], 'PER', mode) if s in es else None for s in B.SEQN], dtype=float)
        I1, I2 = (visit_inputs(B, 1), visit_inputs(B, 2))
        if dm_mode == 'confirmed':
            h1, h2 = (B.hba1c1.to_numpy(float), B.hba1c2.to_numpy(float))
            g1, g2 = (np.where(np.isfinite(h1), h1, h2), np.where(np.isfinite(h2), h2, h1))
            conf = np.where((B.fixed_dm.to_numpy() == 1) | (np.nan_to_num(g1) >= 6.5) & (np.nan_to_num(g2) >= 6.5), 1.0, 0.0)
            I1 = {**I1, 'dm': conf}
            I2 = {**I2, 'dm': conf}
        if mode != 'r23':
            I1 = {**I1, 'sbp': M.clamp(sb1, 90, 200)}
            I2 = {**I2, 'sbp': M.clamp(sb2, 90, 200)}

        def flags(I, egfr_raw, uacr_raw):
            ckd = egfr_raw < 60
            if ckd_mode == 'egfr_or_uacr':
                ckd = ckd | (np.nan_to_num(uacr_raw, nan=0.0) >= 30)
            cvd = M.risk_from_lp(lp_of(F, I, 'base', 'total_cvd'))
            return ((I['dm'] == 1) | ckd | (cvd >= 7.5), cvd)
        e1, e2 = (B.egfr1.to_numpy(float), B.egfr2.to_numpy(float))
        u1, u2 = (B.uacr1.to_numpy(float), B.uacr2.to_numpy(float))
        hr1, cvd1 = flags(I1, e1, u1)
        hr2, cvd2 = flags(I2, e2, u2)
        if mode == 'r23' and dm_mode == 'visit':
            P = pd.read_csv(ST / 'analysis_v1/stage3/stage3_participants.csv').set_index('SEQN').loc[B.SEQN]
            assert np.max(np.abs(cvd1 - P.riskT1.to_numpy())) < 1e-09 and np.max(np.abs(cvd2 - P.riskT2.to_numpy())) < 1e-09, 'total-CVD risk not reproduced'
        untreated = (B.bptreat.to_numpy() == 0) & np.isfinite(db1) & np.isfinite(db2)
        if keep is not None:
            untreated &= keep
        idx = np.where(untreated)[0]
        D1, D2 = (decision(sb1, db1, hr1), decision(sb2, db2, hr2))
        rx1, rx2 = (D1 >= 2, D2 >= 2)
        excl = np.where((B.bptreat.to_numpy() == 0) & ~(np.isfinite(db1) & np.isfinite(db2)))[0]
        res = {'excluded_seqn': [int(B.SEQN.iloc[i]) for i in excl], 'untreated_missing_uacr_any_visit': int((~(np.isfinite(u1) & np.isfinite(u2)))[idx].sum()), 'n_untreated_eligible': int(len(idx)), 'excluded_missing_dbp': int(((B.bptreat.to_numpy() == 0) & ~(np.isfinite(db1) & np.isfinite(db2))).sum()), 'theta_rx': prop((rx1 != rx2)[idx], cl[idx]), 'drug_now_v1': prop(rx1[idx], cl[idx]), 'drug_now_v2': prop(rx2[idx], cl[idx]), 'htn_label_flip': prop(((stage(sb1, db1) >= 1) != (stage(sb2, db2) >= 1))[idx], cl[idx]), 'stage2_flip': prop(((stage(sb1, db1) == 2) != (stage(sb2, db2) == 2))[idx], cl[idx]), 'transitions': pd.crosstab(pd.Series(D1[idx], name='v1'), pd.Series(D2[idx], name='v2')).reindex(index=range(4), columns=range(4), fill_value=0).values.tolist()}
        a, b_ = (rx1[idx], rx2[idx])
        res['p_not_rx2_given_rx1'] = {'events': int((a & ~b_).sum()), 'n': int(a.sum()), 'estimate': float((a & ~b_).sum() / max(a.sum(), 1))}
        res['p_rx2_given_not_rx1'] = {'events': int((~a & b_).sum()), 'n': int((~a).sum()), 'estimate': float((~a & b_).sum() / max((~a).sum(), 1))}
        st1, st2 = (stage(sb1, db1), stage(sb2, db2))
        risk_only = (st1 == 1) & (st2 == 1) & (rx1 != rx2)
        res['flips_stage1_both_highrisk_change'] = int(risk_only[idx].sum())
        dmck_same = ((I1['dm'] == 1) == (I2['dm'] == 1)) & ((e1 < 60) == (e2 < 60)) & ((np.nan_to_num(u1) >= 30) == (np.nan_to_num(u2) >= 30))
        res['flips_stage1_both_cvd75_only'] = int((risk_only & dmck_same & ((cvd1 >= 7.5) != (cvd2 >= 7.5)))[idx].sum())
        res['flips_total'] = int((rx1 != rx2)[idx].sum())
        if mode == 'r23' and ckd_mode == 'egfr_or_uacr' and (dm_mode == 'visit') and (keep is None):

            def hyb(Ibp, Inon, sbp, dbp, egfr_raw, uacr_raw):
                I = {**Inon, 'sbp': Ibp['sbp']}
                hr, _ = flags(I, egfr_raw, uacr_raw)
                return decision(sbp, dbp, hr) >= 2
            bp_only_r1 = hyb(I2, I1, sb2, db2, e1, u1)
            non_only_r1 = hyb(I1, I2, sb1, db1, e2, u2)
            c = {k_: v_.astype(float) for k_, v_ in {'bp_r1': bp_only_r1 != rx1, 'non_r1': non_only_r1 != rx1, 'bp_r2': non_only_r1 != rx2, 'non_r2': bp_only_r1 != rx2, 'all': rx1 != rx2}.items()}
            res['single_group'] = {'bp_only_sym': float((c['bp_r1'] + 0.0 + c['bp_r2'])[idx].mean() / 2), 'nonbp_only_sym': float((c['non_r1'] + 0.0 + c['non_r2'])[idx].mean() / 2), 'bp_only_ref1': prop(c['bp_r1'][idx], cl[idx]), 'nonbp_only_ref1': prop(c['non_r1'][idx], cl[idx])}
            psi_bp = 0.5 * (0.5 * (c['bp_r1'] + 0.0) + 0.5 * (c['all'] - c['non_r1'] + 0.0)) + 0.5 * (0.5 * (c['bp_r2'] + 0.0) + 0.5 * (c['all'] - c['non_r2'] + 0.0))
            psi_non = c['all'] + 0.0 - psi_bp
            res['shapley'] = {'BP': S.paired_ci(psi_bp[idx], cl[idx]), 'nonBP': S.paired_ci(psi_non[idx], cl[idx])}
            geo = lambda x, y: np.sqrt(x * y)
            h1, h2 = (B.hba1c1.to_numpy(float), B.hba1c2.to_numpy(float))
            hbar = np.where(np.isfinite(h1) & np.isfinite(h2), (h1 + h2) / 2, np.where(np.isfinite(h1), h1, h2))
            dmA = np.where((B.fixed_dm.to_numpy() == 1) | (np.nan_to_num(hbar) >= 6.5), 1.0, 0.0)
            scA = geo(B.scr1.to_numpy(float), B.scr2.to_numpy(float))
            egA = np.asarray(M.ckd_epi_2021(scA, B.age, B.female), float)
            uA = np.where(np.isfinite(u1) & np.isfinite(u2), geo(np.maximum(u1, 1e-09), np.maximum(u2, 1e-09)), np.where(np.isfinite(u1), u1, u2))
            sbA, dbA = ((sb1 + sb2) / 2, (db1 + db2) / 2)
            IA = {'sbp': M.clamp(sbA, 90, 200), 'tc': M.clamp(geo(B.tc1.to_numpy(float), B.tc2.to_numpy(float)), 130, 320), 'hdl': M.clamp(geo(B.hdl1.to_numpy(float), B.hdl2.to_numpy(float)), 20, 100), 'egfr': M.clamp(egA, 15, 140), 'dm': dmA, 'uacr': uA}
            hrA, _ = flags(IA, egA, uA)
            rxA = decision(sbA, dbA, hrA) >= 2
            res['vs_two_visit_mean'] = {'v1_differs': prop((rx1 != rxA)[idx], cl[idx]), 'v2_differs': prop((rx2 != rxA)[idx], cl[idx]), 'mean_single_visit_differs': float(((rx1 != rxA) + 0.0 + (rx2 != rxA))[idx].mean() / 2), 'mean_single_visit_differs_ci': S.paired_ci((((rx1 != rxA) + 0.0 + (rx2 != rxA)) / 2)[idx], cl[idx]), 'drug_now_two_visit': prop(rxA[idx], cl[idx])}
            tr = (B.bptreat.to_numpy() == 1) & np.isfinite(db1) & np.isfinite(db2)
            ti = np.where(tr)[0]
            res['treated_above_goal_flip'] = prop(((stage(sb1, db1) >= 1) != (stage(sb2, db2) >= 1))[ti], cl[ti])
            exm = s12.read('exam', {'PEPTECH': (1498, 1502)}, keep=set(B.SEQN))
            e1c = np.array([exm[q]['PEPTECH'].strip() for q in B.SEQN])
            r_ = S.prop_ci((rx1 != rx2)[idx].astype(float), e1c[idx])
            res['theta_rx_visit1_physician'] = {k: r_[k] for k in ['estimate', 'lower', 'upper', 'deff', 'G']}
            res['theta_rx_visit1_physician_wild'] = wild_ci((rx1 != rx2)[idx].astype(float), e1c[idx])
            res['flips_with_stage_change'] = int(((stage(sb1, db1) != stage(sb2, db2)) & (rx1 != rx2))[idx].sum())
            res['_idx'] = idx
            res['_flip'] = rx1 != rx2
        return res
    out['primary'] = run()
    out['sens_all3_readings'] = run(mode='all3')
    out['sens_ckd_egfr_only'] = run(ckd_mode='egfr')
    out['sens_confirmed_dm_fixed'] = run(dm_mode='confirmed')
    out['sens_exclude_referred'] = run(keep=B.peplevel.to_numpy() == 3)
    return out

def npmle(m, s_m, w=None, K=400, iters=3000, init=None, tol=1e-09, info=None):
    w = np.ones(len(m)) if w is None else np.asarray(w, float)
    grid = np.linspace(m.min() - 1, m.max() + 1, K)
    L = norm.pdf((m[:, None] - grid[None, :]) / s_m) / s_m
    p = np.full(K, 1 / K) if init is None else init.copy()
    for _ in range(iters):
        f = L @ p
        post = L * p / np.maximum(f[:, None], 1e-300)
        pn = (w[:, None] * post).sum(0) / w.sum()
        dmax = float(np.max(np.abs(pn - p)))
        p = pn
        if dmax < tol:
            if info is not None:
                info.update(converged=True, iterations=_ + 1, last_max_change=dmax)
            break
    else:
        if info is not None:
            info.update(converged=False, iterations=iters, last_max_change=dmax)
    if info is not None:
        info['loglik'] = float(np.sum(w * np.log(L @ p)))
    return (grid, p)

def npmle_sq(m, s_m, w=None, K=400, tol=1e-09, maxit=60000):
    w = np.ones(len(m)) if w is None else np.asarray(w, float)
    W = w / w.sum()
    grid = np.linspace(m.min() - 1, m.max() + 1, K)
    L = norm.pdf((m[:, None] - grid[None, :]) / s_m) / s_m
    F = lambda q: (W[:, None] * (L * q / np.maximum((L @ q)[:, None], 1e-300))).sum(0)
    ll = lambda q: float(np.sum(w * np.log(np.maximum(L @ q, 1e-300))))
    p = np.full(K, 1 / K)
    n_eval = 0
    conv = False
    while n_eval < maxit:
        p1 = F(p)
        p2 = F(p1)
        n_eval += 2
        r = p1 - p
        v = p2 - p1 - r
        nv = np.linalg.norm(v)
        if nv == 0:
            p = p2
            conv = True
            break
        a = min(-1.0, -np.linalg.norm(r) / nv)
        pn = np.maximum(p - 2 * a * r + a * a * v, 0)
        pn = F(pn / pn.sum())
        n_eval += 1
        if ll(pn) < ll(p2):
            pn = p2
        dmax = float(np.max(np.abs(pn - p)))
        p = pn
        if dmax < tol:
            conv = True
            break
    return (grid, p, {'converged': conv, 'em_evaluations': n_eval, 'loglik': ll(p)})

def gsmooth(m, s_m, K=400, deg=7, c=1.0):
    from scipy.optimize import minimize
    from scipy.special import logsumexp
    grid = np.linspace(m.min() - 1, m.max() + 1, K)
    z = (grid - grid.mean()) / grid.std()
    Q = np.column_stack([z ** k for k in range(1, deg + 1)])
    logL = norm.logpdf((m[:, None] - grid[None, :]) / s_m) - math.log(s_m)

    def f(a):
        eta = Q @ a
        lg = eta - logsumexp(eta)
        ll = logsumexp(logL + lg[None, :], axis=1)
        post = np.exp(logL + lg[None, :] - ll[:, None])
        g_ = np.exp(lg)
        grad = -(Q.T @ (post.sum(0) - len(m) * g_)) + 2 * c * a
        return (-ll.sum() + c * np.sum(a ** 2), grad)
    r = minimize(f, np.zeros(deg), jac=True, method='L-BFGS-B', options={'maxiter': 5000, 'gtol': 1e-08})
    eta = Q @ r.x
    return (grid, np.exp(eta - logsumexp(eta)), bool(r.success))

def cat_bounds(c):
    lo = np.concatenate([[-np.inf], LT])[c]
    hi = np.concatenate([LT, [np.inf]])[c]
    return (lo, hi)

def kappa_theta(grid, p, delta, sigma, err=None):
    cg = np.searchsorted(LT, grid, side='right')
    lo, hi = cat_bounds(cg)
    edges = np.concatenate([[-np.inf], LT, [np.inf]])
    if err is None:
        kap = 0.0
        for s in (-1, 1):
            a = s * delta / 2
            stay = norm.cdf((hi - grid - a) / sigma) - norm.cdf((lo - grid - a) / sigma)
            kap += 0.5 * np.sum(p * (1 - stay))
        P1 = np.stack([norm.cdf((edges[j + 1] - grid + delta / 2) / sigma) - norm.cdf((edges[j] - grid + delta / 2) / sigma) for j in range(4)], 1)
        P2 = np.stack([norm.cdf((edges[j + 1] - grid - delta / 2) / sigma) - norm.cdf((edges[j] - grid - delta / 2) / sigma) for j in range(4)], 1)
        th = np.sum(p * (1 - (P1 * P2).sum(1)))
        return (float(kap), float(th))
    e = np.asarray(err)
    kap = 0.0
    P = []
    for s in (-1, 1):
        x = grid[:, None] + s * delta / 2 + e[None, :]
        cx = np.searchsorted(LT, x, side='right')
        kap += 0.5 * np.sum(p * (cx != cg[:, None]).mean(1))
        P.append(np.stack([(cx == j).mean(1) for j in range(4)], 1))
    th = np.sum(p * (1 - (P[0] * P[1]).sum(1)))
    return (float(kap), float(th))

def kappa_r(grid, p, delta, sigma, risks):
    out = {}
    for r in risks:
        l1 = logit(r / 100)
        post = p * norm.pdf((l1 - grid + delta / 2) / sigma)
        post /= post.sum()
        out[str(r)] = float(np.sum(post * (np.searchsorted(LT, grid, side='right') != np.searchsorted(LT, l1, side='right'))))
    return out

def e4(P, cl):
    L1, L2 = (P.lpB_m0.to_numpy(float), P.lpB_m15.to_numpy(float))
    ch = P.changeB.to_numpy(float)
    assert np.max(np.abs(100 * expit(L1) - P.riskB_m0)) < 1e-09 and np.max(np.abs(100 * expit(L2) - P.riskB_m15)) < 1e-09
    d = L2 - L1
    m = (L1 + L2) / 2

    def fit(w=None, init=None, risks=None, err_mode=None, age=None):
        w_ = np.ones(len(d)) if w is None else w
        delta = np.sum(w_ * d) / np.sum(w_)
        sigma = math.sqrt(np.sum(w_ * (d - delta) ** 2) / (np.sum(w_) - 1) / 2)
        th_obs = np.sum(w_ * ch) / np.sum(w_)
        if age is not None:
            kap = th = 0.0
            for lo_, hi_ in [(30, 39), (40, 59), (60, 79)]:
                sel = (age >= lo_) & (age <= hi_)
                wb = w_ * sel
                dl = np.sum(wb * d) / np.sum(wb)
                sg = math.sqrt(np.sum(wb * (d - dl) ** 2) / (np.sum(wb) - 1) / 2)
                g, p = npmle(m[sel], sg / math.sqrt(2), wb[sel])
                k_, t_ = kappa_theta(g, p, dl, sg)
                kap += k_ * np.sum(wb) / np.sum(w_)
                th += t_ * np.sum(wb) / np.sum(w_)
            return {'kappa_model': kap, 'theta_model': th, 'theta_obs': th_obs, 'kappa': kap * th_obs / th}
        info = {}
        g, p = npmle(m, sigma / math.sqrt(2), w_, init=init, iters=100000 if init is None else 400, info=info)
        err = None
        if err_mode == 'empirical':
            e = (d - delta) / math.sqrt(2)
            err = np.concatenate([e, -e])
        k_, t_ = kappa_theta(g, p, delta, sigma, err)
        res = {'delta': delta, 'sigma': sigma, 'kappa_model': k_, 'theta_model': t_, 'theta_obs': th_obs, 'kappa': k_ * th_obs / t_, '_p': p, 'npmle': info}
        if risks is not None:
            res['kappa_r'] = kappa_r(g, p, delta, sigma, risks)
        return res
    risks = [e_['risk_pct'] for e_ in LK['lookup']]
    base = fit(risks=risks)
    W = S.boot_weights(cl, B_BOOT, SEED)[:1000]
    boots = [fit(W[b], init=base['_p'], risks=risks) for b in range(len(W))]
    out = {'primary': {k: v for k, v in base.items() if not k.startswith('_')}, 'grid_points': 400, 'bootstrap_em': 'warm start from the full-sample solution, at most 400 EM iterations per replicate', 'ci': {'kappa': S.pct_ci([b['kappa'] for b in boots]), 'kappa_model': S.pct_ci([b['kappa_model'] for b in boots]), 'ratio_kappa_theta': S.pct_ci([b['kappa_model'] / b['theta_model'] for b in boots]), 'kappa_r': {r: S.pct_ci([b['kappa_r'][r] for b in boots]) for r in base['kappa_r']}}, 'lower_bound_theta_half': float(ch.mean() / 2), 'sens_empirical_errors': {k: v for k, v in fit(err_mode='empirical').items() if not k.startswith('_')}, 'sens_sigma_by_age': fit(age=P.age.to_numpy()), 'sens_smooth_g': (lambda gs: {'converged': gs[2], **dict(zip(['kappa_model', 'theta_model'], kappa_theta(gs[0], gs[1], base['delta'], base['sigma']))), 'kappa_r': kappa_r(gs[0], gs[1], base['delta'], base['sigma'], risks)})(gsmooth(m, base['sigma'] / math.sqrt(2))), 'n_boot': len(W)}
    return out

def e5(P, B, cl):
    ex = s12.read('exam', {'PEPTECH': (1498, 1502)}, keep=set(B.SEQN))
    es = s12.read('examse', {'PERTECH': (258, 262)}, keep=set(B.SEQN))
    e1 = np.array([ex[s]['PEPTECH'].strip() for s in P.index])
    e2 = np.array([es[s]['PERTECH'].strip() for s in P.index])
    y = P.changeB.to_numpy(float)
    out = {'n_examiners_v1': int(len(set(e1))), 'n_examiners_v2': int(len(set(e2)))}
    for nm, g in [('pseudo_psu', cl), ('visit1_physician', e1), ('visit2_physician', e2), ('physician_pair', np.char.add(np.char.add(e1, '_'), e2))]:
        r = S.prop_ci(y, g)
        out[nm] = {k: r[k] for k in ['estimate', 'lower', 'upper', 'deff', 'G', 'method']}
    rng = np.random.default_rng(20261010)
    ids, inv = np.unique(e1, return_inverse=True)
    G = len(ids)
    ybar = y.mean()
    u = y - ybar
    V, _ = S.cr1_var(y, e1)
    se = math.sqrt(V)
    webb = np.array([-math.sqrt(1.5), -1, -math.sqrt(0.5), math.sqrt(0.5), 1, math.sqrt(1.5)])
    ts = []
    for _ in range(9999):
        wg = rng.choice(webb, G)
        ys = ybar + wg[inv] * u
        Vs, _ = S.cr1_var(ys, e1)
        ts.append((ys.mean() - ybar) / math.sqrt(Vs))
    q = np.quantile(ts, [0.025, 0.975])
    out['visit1_physician_wild_bootstrap'] = {'lower': float(ybar - q[1] * se), 'upper': float(ybar - q[0] * se), 'replicates': 9999, 'weights': 'Webb six-point'}
    same = (e1 == e2) & (e1 != '88888')
    pi = predict_pi(P.lpB_m0.to_numpy(float), LK['models']['base'])
    for nm, sel in [('same_physician', same), ('different_physician', ~same)]:
        out[nm] = {'n': int(sel.sum()), 'changes': int(y[sel].sum()), 'estimate': float(y[sel].mean()), 'expected_changes_pi_clin': float(pi[sel].sum())}
    out['pi_clin_mean_equals_theta'] = float(abs(pi.mean() - y.mean()))
    return out

def e6(P, cl):
    W = S.boot_weights(cl, B_BOOT, SEED)
    psi = {g: P[f'psiB_{g}'].to_numpy(float) for g in ['BP', 'LIPID', 'KIDNEY', 'GLYC']}
    bw = lambda v: W @ v / W.sum(1)
    ci_bp = S.pct_ci(100 * bw(psi['BP']))
    fz = R3['K2_base']['Psi']['BP']['ci_pp']
    assert abs(ci_bp[0] - fz[0]) < 1e-09 and abs(ci_bp[1] - fz[1]) < 1e-09, ('frozen Psi_BP CI not reproduced', ci_bp, fz)
    diff = 100 * (psi['BP'] - psi['LIPID'])
    return {'psi_bp_minus_lipid_pp': float(diff.mean()), 'ci_pp': S.pct_ci(W @ diff / W.sum(1)), 'p_bp_gt_lipid': float(np.mean(bw(psi['BP']) > bw(psi['LIPID']))), 'frozen_psi_bp_ci_reproduced': True}

def e7(P, B, T, cl):
    F = fixed(B)
    L1, L2 = (P.lpB_m0.to_numpy(float), P.lpB_m15.to_numpy(float))
    ch = P.changeB.to_numpy(float)
    d = L2 - L1
    rd = s12.read('exam', {'DMARETHN': s12.POS['exam']['DMARETHN']}, keep=set(T.SEQN))
    race_t = np.array([s12.code(rd[s]['DMARETHN']) for s in T.SEQN])
    race_s = P.dmarethn.to_numpy()

    def cell(age, race):
        a = np.digitize(age, [40, 60])
        r = np.where(race == 2, 1, np.where(race == 3, 2, 0))
        return a * 3 + r
    cs, ct = (cell(P.age.to_numpy(), race_s), cell(T.age.to_numpy(), race_t))
    FT = {k: T[k].to_numpy(float) for k in ['female', 'age', 'smoking', 'bptreat', 'statin']}
    IT = {'sbp': M.clamp(T.sbp1, 90, 200), 'tc': M.clamp(T.tc1, 130, 320), 'hdl': M.clamp(T.hdl1, 20, 100), 'egfr': M.clamp(T.egfr1, 15, 140), 'dm': T.dm1.to_numpy(float), 'uacr': None}
    LTB = lp_of(FT, IT, 'base')
    wT = T.w.to_numpy(float)
    stT = T.stratum.to_numpy()
    psT = (T.stratum * 10 + T.psu).to_numpy()

    def dist(lp):
        dd = np.abs(lp[:, None] - LT[None, :])
        return (dd.min(1), dd.argmin(1))
    ds, ns = dist(L1)
    dt, nt = dist(LTB)

    def sig(w):
        sa = math.sqrt(np.sum(w * (d - np.sum(w * d) / np.sum(w)) ** 2) / np.sum(w))
        sc = np.array([math.sqrt(np.sum(w * (cs == c) * (d - np.sum(w * (cs == c) * d) / np.sum(w * (cs == c))) ** 2) / np.sum(w * (cs == c))) if np.sum(w * (cs == c)) > 0 else np.nan for c in range(9)])
        return (sa, np.where(np.isfinite(sc) & (sc > 0), sc, SC0[0] if SC0 else sa))

    def design(dv, nr, knots):
        X = np.column_stack([np.ones(len(dv)), S.rcs_basis(dv, knots)])
        return np.column_stack([X, nr == 1, nr == 2]).astype(float)
    SC0 = []
    sa, sc = sig(np.ones(len(d)))
    SC0.append(sc.copy())
    assert np.all(np.bincount(cs, minlength=9) >= 10), np.bincount(cs, minlength=9).tolist()
    dss = ds * sa / sc[cs]
    knots = np.percentile(dss, [10, 50, 90])
    beta = S.logit_fit(design(dss, ns, knots), ch)
    th_het = float(np.sum(wT * expit(design(dt * sa / sc[ct], nt, knots) @ beta)) / np.sum(wT))
    W = S.boot_weights(cl, B_BOOT, SEED)
    rng = np.random.default_rng(SEED + 1)
    vals, diffs, v0s = ([], [], [])
    kf, bf = (np.asarray(LK['models']['base']['knots']), np.asarray(LK['models']['base']['beta']))
    Xs0 = design(ds, ns, kf)
    Xt0 = design(dt, nt, kf)
    st_ids = np.unique(stT)
    for b in range(2000):
        wb = W[b]
        sab, scb = sig(wb)
        beta_b = S.logit_fit(design(ds * sab / scb[cs], ns, knots), ch, wb)
        beta0_b = S.logit_fit(Xs0, ch, wb)
        wt = wT.copy()
        for h in st_ids:
            mk = stT == h
            ps_ = np.unique(psT[mk])
            nh = len(ps_)
            if nh < 2:
                continue
            draw = rng.choice(ps_, nh - 1, replace=True)
            cnt = {p_: int((draw == p_).sum()) for p_ in ps_}
            wt[mk] = wT[mk] * nh / (nh - 1) * np.array([cnt[p_] for p_ in psT[mk]])
        v_het = np.sum(wt * expit(design(dt * sab / scb[ct], nt, knots) @ beta_b)) / np.sum(wt)
        v0 = np.sum(wt * expit(Xt0 @ beta0_b)) / np.sum(wt)
        vals.append(v_het)
        diffs.append(v_het - v0)
        v0s.append(v0)
    th0 = float(np.sum(wT * expit(Xt0 @ bf)) / np.sum(wT))
    assert abs(th0 - R3['K4_base']['theta_std']) < 1e-09, ('frozen theta_std not reproduced', th0)
    ci0 = S.pct_ci(v0s)
    assert max(abs(ci0[0] - R3['K4_base']['theta_std_ci'][0]), abs(ci0[1] - R3['K4_base']['theta_std_ci'][1])) < 1e-09, ('frozen K4 CI not reproduced', ci0)
    return {'theta_std_het': th_het, 'ci': S.pct_ci(vals), 'theta_std_frozen': th0, 'frozen_ci_reproduced': ci0, 'difference_pp': 100 * (th_het - th0), 'difference_ci_pp': [100 * x for x in S.pct_ci(diffs)], 'sigma_all': sa, 'sigma_cells': sc.tolist(), 'cell_n_sample': np.bincount(cs, minlength=9).tolist(), 'cell_weight_target': [float(wT[ct == c].sum() / wT.sum()) for c in range(9)], 'cells': 'age 30-39/40-59/60-79 x (NHW+other, NHB, MA)', 'replicates': 2000}

def taylor_prop(y, w, strata, psu):
    y, w = (np.asarray(y, float), np.asarray(w, float))
    W_ = w.sum()
    p = np.sum(w * y) / W_
    z = w * (y - p) / W_
    V = 0.0
    nps = 0
    hs = np.unique(strata)
    for h in hs:
        mk = strata == h
        ps_ = np.unique(psu[mk])
        nps += len(ps_)
        tot = np.array([z[mk & (psu == q)].sum() for q in ps_])
        if len(ps_) > 1:
            V += len(ps_) / (len(ps_) - 1) * np.sum((tot - tot.mean()) ** 2)
    df = nps - len(hs)
    se = math.sqrt(V) / (p * (1 - p))
    q = tdist.ppf(0.975, df)
    return {'estimate': float(p), 'lower': float(expit(logit(p) - q * se)), 'upper': float(expit(logit(p) + q * se)), 'df': int(df), 'events': int(y.sum()), 'n': int(len(y))}

def e8():
    S2P = pd.read_csv(ROOT / 'studies/uacr_prevent/final_and_external_v1/primary_participant_results.csv')
    ev = json.loads((ST / 'manuscript_evidence_v1/study2_evidence.json').read_text(encoding='utf-8'))

    def risk(u):
        return M.risk_from_lp(M.lp(S2P.sex.to_numpy(float), S2P.age.to_numpy(float), S2P.tc.to_numpy(float), S2P.hdl.to_numpy(float), S2P.sbp.to_numpy(float), S2P.dm.to_numpy(float), S2P.smoking.to_numpy(float), S2P.egfr.to_numpy(float), S2P.bptreat.to_numpy(float), S2P.statin.to_numpy(float), uacr=u, model='uacr', outcome='ascvd'))
    r1, r2 = (risk(S2P.URDACT.to_numpy(float)), risk(S2P.URDACT2.to_numpy(float)))
    rep = float(max(np.max(np.abs(r1 - S2P.risk_first)), np.max(np.abs(r2 - S2P.risk_second))))
    if rep > 1e-09:
        return {'run': False, 'reason': f'locked risks not reproduced (max diff {rep:.2e} pp)'}
    w, st_, ps_ = (S2P.WTMEC2YR.to_numpy(float), S2P.SDMVSTRA.to_numpy(), S2P.SDMVPSU.to_numpy())
    c1, c2 = (M.category(r1), M.category(r2))
    prim = taylor_prop(c1 != c2, w, st_, ps_)
    assert abs(prim['estimate'] - ev['primary']['any_category_change']['estimate']) < 1e-12 and abs(prim['lower'] - ev['primary']['any_category_change']['lower95']) < 1e-09, 'Study 2 primary not reproduced'
    lr = np.log(S2P.URDACT2.to_numpy(float) / S2P.URDACT.to_numpy(float))
    gmr = float(np.exp(np.sum(w * lr) / np.sum(w)))
    r2s = risk(S2P.URDACT2.to_numpy(float) / gmr)
    c2s = M.category(r2s)
    return {'run': True, 'max_reproduction_diff_pp': rep, 'primary_reproduced': prim, 'weighted_gmr_first_morning_over_random': gmr, 'sd_log_ratio_weighted': float(math.sqrt(np.sum(w * (lr - np.log(gmr)) ** 2) / np.sum(w))), 'shift_removed_any': taylor_prop(c1 != c2s, w, st_, ps_), 'shift_removed_down': taylor_prop(c2s < c1, w, st_, ps_), 'shift_removed_up': taylor_prop(c2s > c1, w, st_, ps_), 'original_down': taylor_prop(c2 < c1, w, st_, ps_), 'original_up': taylor_prop(c2 > c1, w, st_, ps_)}

def main():
    E, B, T = nh3_build.build()
    P = pd.read_csv(ST / 'analysis_v1/stage3/stage3_participants.csv').set_index('SEQN').loc[B.SEQN]
    assert len(P) == 759 and int(P.changeB.sum()) == R3['primary']['events']
    cl = B.cluster.to_numpy()
    F = fixed(B)
    res = {'label': __doc__.split('\n')[0], 'protocol': str(PROTO.relative_to(ROOT)), 'protocol_sha256': hashlib.sha256(PROTO.read_bytes()).hexdigest()}
    print('E3...', flush=True)
    r3 = e3(B, F, cl)
    for k in r3:
        r3[k] = {kk: vv for kk, vv in r3[k].items() if not kk.startswith('_')}
    res['E3'] = r3
    print('E5...', flush=True)
    res['E5'] = e5(P, B, cl)
    print('E6...', flush=True)
    res['E6'] = e6(P, cl)
    print('E4...', flush=True)
    res['E4'] = e4(P, cl)
    print('E7...', flush=True)
    res['E7'] = e7(P, B, T, cl)
    print('E8...', flush=True)
    res['E8'] = e8()
    OUT.write_text(json.dumps(res, indent=1, default=lambda o: o.item() if hasattr(o, 'item') else o.tolist() if hasattr(o, 'tolist') else str(o)) + '\n', encoding='utf-8')
    print('saved', OUT)
if __name__ == '__main__':
    main()
