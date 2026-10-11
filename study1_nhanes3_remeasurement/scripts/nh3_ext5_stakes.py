from pathlib import Path
import hashlib
import json
import math
import numpy as np
import pandas as pd
from scipy.special import expit, logit
from scipy.stats import norm
import nh3_build
import nh3_stats as S
import nh3_models as M
from nh3_stage3 import fixed, visit_inputs
from nh3_ext2_analyses import npmle, gsmooth, npmle_sq
from nh3_ext2021_burden import zones
from nh3_ext3_statin import risk30
ROOT = Path(__file__).resolve().parents[1]
ST = ROOT / 'studies/nhanes3_testretest'
PH = ST / 'analysis_v1_posthoc'
PROTO = ST / 'extensions5_protocol_20261010.md'
OUT = PH / 'ext5_stakes_results.json'
SEED_BOOT, N_BOOT, MC, MC_B, SEED_MC = (20261009, 1000, 200000, 20000, 20261011)

def fit_g(L1, L2, w=None, init=None):
    w_ = np.ones(len(L1)) if w is None else np.asarray(w, float)
    d = L2 - L1
    m = (L1 + L2) / 2
    delta = np.sum(w_ * d) / np.sum(w_)
    sigma = math.sqrt(np.sum(w_ * (d - delta) ** 2) / (np.sum(w_) - 1) / 2)
    info = {}
    g, p = npmle(m, sigma / math.sqrt(2), w_, init=init, iters=100000 if init is None else 400, info=info)
    return (g, p, delta, sigma, info)

def fit_g_sq(L1, L2, w=None, extra_var=0.0, maxit=60000):
    w_ = np.ones(len(L1)) if w is None else np.asarray(w, float)
    d = L2 - L1
    m = (L1 + L2) / 2
    delta = np.sum(w_ * d) / np.sum(w_)
    sigma = math.sqrt(np.sum(w_ * (d - delta) ** 2) / (np.sum(w_) - 1) / 2)
    g, p, info = npmle_sq(m, math.sqrt(sigma ** 2 / 2 + extra_var * sigma ** 2), w_, maxit=maxit)
    return (g, p, delta, sigma, info)

def metrics(g, p, delta, sigma, pt, zone=None, mc=MC, seed=SEED_MC, err=None, sd_single=None, sd_two=None):
    t = logit(pt)
    r = expit(g)
    usual = g >= t
    stake = np.abs(r - pt) / (1 - pt)
    achievable = float(np.sum(p * np.maximum(r - pt, 0)) / (1 - pt))
    if err is None:
        s1 = sigma if sd_single is None else sd_single
        s2 = sigma / math.sqrt(2) if sd_two is None else sd_two
        ps = 0.5 * (1 - norm.cdf((t - g + delta / 2) / s1) + (1 - norm.cdf((t - g - delta / 2) / s1)))
        p2 = 1 - norm.cdf((t - g) / s2)
    else:
        e = np.asarray(err)
        rng = np.random.default_rng(seed)
        ps = 0.5 * ((g[:, None] - delta / 2 + e[None, :] >= t).mean(1) + (g[:, None] + delta / 2 + e[None, :] >= t).mean(1))
        pair = (rng.choice(e, 20000) + rng.choice(e, 20000)) / 2
        p2 = (g[:, None] + pair[None, :] >= t).mean(1)
    dis1 = np.where(usual, 1 - ps, ps)
    dis2 = np.where(usual, 1 - p2, p2)
    reg1, reg2 = (float(np.sum(p * stake * dis1)), float(np.sum(p * stake * dis2)))
    kap = float(np.sum(p * dis1))
    out = {'pt': pt, 'achievable_nb_per_1000': 1000 * achievable, 'kappa_single_vs_usual': kap, 'regret_single_per_1000': 1000 * reg1, 'loss_fraction_single': reg1 / achievable if achievable > 0 else None, 'mean_abs_risk_gap_pp_among_disagreements': float(100 * np.sum(p * dis1 * np.abs(r - pt)) / kap) if kap > 0 else None, 'regret_two_day_per_1000': 1000 * reg2, 'loss_fraction_two_day': reg2 / achievable if achievable > 0 else None}
    nb_all = float(np.sum(p * (r - pt)) / (1 - pt))
    ach_best = achievable - max(0.0, nb_all)
    out.update({'nb_treat_all_per_1000': 1000 * nb_all, 'loss_fraction_single_vs_best_default': reg1 / ach_best if ach_best > 0 else None, 'extra_measurements_per_net_tp_two_day': 1 / (reg1 - reg2) if reg1 > reg2 else None})
    if zone is not None and err is None:
        rng = np.random.default_rng(seed)
        i = rng.choice(len(g), mc, p=p / p.sum())
        T = g[i]
        s = rng.choice([-1.0, 1.0], mc)
        L1 = T + s * delta / 2 + rng.normal(0, sigma, mc)
        L2 = T - s * delta / 2 + rng.normal(0, sigma, mc)
        r1 = 100 * expit(L1)
        inz = (r1 >= zone[0]) & (r1 <= zone[1])
        dec = np.where(inz, (L1 + L2) / 2 >= t, L1 >= t)
        u = T >= t
        regc = float(np.mean(stake[i] * (dec != u)))
        reg1mc = float(np.mean(stake[i] * ((L1 >= t) != u)))
        out.update({'targeted_zone_risk_pct': list(zone), 'targeted_fraction_remeasured': float(inz.mean()), 'regret_targeted_per_1000': 1000 * regc, 'loss_fraction_targeted': regc / achievable if achievable > 0 else None, 'targeted_share_of_single_regret_recovered': (reg1mc - regc) / reg1mc if reg1mc > 0 else None, 'mc_check_regret_single_per_1000': 1000 * reg1mc, 'extra_measurements_per_net_tp_targeted': float(inz.mean()) / (reg1mc - regc) if reg1mc > regc else None})
    return out

def empirical(L1, L2, pt):
    t = logit(pt)
    m = (L1 + L2) / 2
    rm = expit(m)
    stake = np.abs(rm - pt) / (1 - pt)
    um = m >= t
    ach = float(np.mean(np.maximum(rm - pt, 0)) / (1 - pt))
    dis = [(L >= t) != um for L in (L1, L2)]
    reg = float(np.mean([np.mean(stake * d_) for d_ in dis]))
    nd = sum((d_.sum() for d_ in dis))
    return {'reference': 'two-visit mean (understates regret)', 'achievable_nb_per_1000': 1000 * ach, 'regret_single_per_1000': 1000 * reg, 'loss_fraction_single': reg / ach, 'disagreement_single_vs_mean': float(np.mean([d_.mean() for d_ in dis])), 'mean_abs_risk_gap_pp_among_disagreements': float(100 * sum(((np.abs(rm - pt) * d_).sum() for d_ in dis)) / nd) if nd else None}

def reading(x):
    return 'minor (<5%)' if x < 0.05 else 'modest (5-<15%)' if x < 0.15 else 'substantial (>=15%)'
CI_KEYS = ['loss_fraction_single', 'loss_fraction_two_day', 'regret_single_per_1000', 'achievable_nb_per_1000', 'kappa_single_vs_usual', 'loss_fraction_targeted', 'targeted_share_of_single_regret_recovered', 'targeted_fraction_remeasured', 'mean_abs_risk_gap_pp_among_disagreements', 'regret_two_day_per_1000', 'loss_fraction_single_vs_best_default']
_WK = {}

def _init_worker(L1, L2, W, specs):
    _WK.update(L1=L1, L2=L2, W=W, specs=specs)

def _boot_uniform(b):
    gb, pb, db, sb, info = fit_g_sq(_WK['L1'], _WK['L2'], _WK['W'][b], maxit=20000)
    return (gb, pb, db, sb, info['converged'])

def _cached(name, key, fn):
    import pickle
    d = PH / '_cache_e11'
    d.mkdir(exist_ok=True)
    f = d / f'{name}_{key[:16]}.pkl'
    if f.exists():
        return pickle.loads(f.read_bytes())
    v = fn()
    f.write_bytes(pickle.dumps(v))
    return v

def by_band(L1, L2, age, bands, pt):
    reg = ach = 0.0
    parts = {}
    for lo, hi in bands:
        s = (age >= lo) & (age <= hi)
        g, p, d, sg, info = fit_g_sq(L1[s], L2[s])
        m = metrics(g, p, d, sg, pt)
        wgt = s.sum() / len(L1)
        reg += wgt * m['regret_single_per_1000']
        ach += wgt * m['achievable_nb_per_1000']
        parts[f'{lo}-{hi}'] = {'n': int(s.sum()), 'sigma': sg, 'loss_fraction_single': m['loss_fraction_single'], 'converged': info['converged']}
    return {'loss_fraction_single': reg / ach, 'regret_single_per_1000': reg, 'achievable_nb_per_1000': ach, 'bands': parts}

def inflated(L1, L2, pt, k):
    g, p, d, sg, info = fit_g_sq(L1, L2, extra_var=k - 1)
    m = metrics(g, p, d, sg, pt, sd_single=math.sqrt(k) * sg, sd_two=math.sqrt(sg ** 2 / 2 + (k - 1) * sg ** 2))
    return {k_: m[k_] for k_ in ['loss_fraction_single', 'loss_fraction_two_day', 'regret_single_per_1000', 'achievable_nb_per_1000', 'kappa_single_vs_usual']} | {'converged': info['converged']}

def analyse(L1, L2, cl_full, sel, pts_zones, label, sens=False, age=None, bands=None, pool_n=12):
    cl = cl_full[sel]
    g, p, delta, sigma, info = fit_g(L1, L2)
    res = {'label': label, 'n': len(L1), 'delta': delta, 'sigma': sigma, 'npmle': info, 'thresholds': {}}
    W = S.boot_weights(cl_full, 4000, SEED_BOOT)[:N_BOOT][:, sel]
    key = hashlib.sha256(L1.tobytes() + L2.tobytes() + W.tobytes()).hexdigest()
    fits = _cached('warm', key, lambda: [fit_g(L1, L2, W[b], init=p) for b in range(N_BOOT)])
    gq, pq, dq, sq, infq = fit_g_sq(L1, L2)

    def _run_pool():
        from multiprocessing import Pool
        with Pool(pool_n, initializer=_init_worker, initargs=(L1, L2, W, pts_zones)) as pool:
            return pool.map(_boot_uniform, range(N_BOOT), chunksize=8)
    ufits = _cached('uniform', key, _run_pool)
    uni = [([metrics(gb, pb, db, sb, pt, zone, mc=MC_B, seed=SEED_MC + b) for pt, zone in pts_zones], c) for b, (gb, pb, db, sb, c) in enumerate(ufits)]
    res['post_audit'] = {'squarem_full_sample': infq, 'uniform_bootstrap_converged_at_tol': int(sum((c for _, c in uni))), 'n_boot': N_BOOT, 'uniform_refit_cap': '20,000 EM evaluations (functional stable to <0.01 pp in checks)'}
    for j, (pt, zone) in enumerate(pts_zones):
        t = logit(pt)
        obs = ((L1 >= t) != (L2 >= t)).astype(float)
        mm = metrics(g, p, delta, sigma, pt, zone)
        bs = [metrics(gb, pb, db, sb, pt, zone, mc=MC_B, seed=SEED_MC + b) for b, (gb, pb, db, sb, _) in enumerate(fits)]
        ci = {k: S.pct_ci([x[k] for x in bs]) for k in CI_KEYS if k in mm and mm[k] is not None}
        ciu = {k: S.pct_ci([u[0][j][k] for u in uni if u[0][j].get(k) is not None]) for k in CI_KEYS if k in mm and mm[k] is not None}
        mq = metrics(gq, pq, dq, sq, pt, zone)
        entry = {**mm, 'ci': ci, 'observed_flip_between_visits': S.prop_ci(obs, cl)['estimate'], 'empirical': empirical(L1, L2, pt), 'reading_point': reading(mm['loss_fraction_single']), 'reading_ci': [reading(v) for v in ci['loss_fraction_single']], 'ci_uniform_start': ciu, 'reading_ci_uniform_start': [reading(v) for v in ciu['loss_fraction_single']], 'squarem_check': {'loss_fraction_single': mq['loss_fraction_single'], 'abs_diff_vs_plain_em': abs(mq['loss_fraction_single'] - mm['loss_fraction_single'])}, 'note': 'kappa_single_vs_usual = single-visit decision disagreement with the usual-level decision (not Cohen kappa)'}
        if sens:
            d = L2 - L1
            e = (d - delta) / math.sqrt(2)
            err = np.concatenate([e, -e])
            entry['sens_empirical_errors'] = metrics(g, p, delta, sigma, pt, err=err)
            gs = gsmooth((L1 + L2) / 2, sigma / math.sqrt(2))
            entry['sens_smooth_g'] = {'converged': gs[2], **metrics(gs[0], gs[1], delta, sigma, pt, zone)}
        if bands is not None:
            entry['post_audit_sens_sigma_by_age'] = by_band(L1, L2, age, bands, pt)
        entry['post_audit_sens_variance_inflation'] = {str(k): inflated(L1, L2, pt, k) for k in (1.5, 2.0)}
        res['thresholds'][f'{100 * pt:g}'] = entry
    return res

def main():
    ptxt = PROTO.read_bytes()
    out = {'label': 'E11: net benefit lost by single-visit risk-threshold decisions (decision-curve logic), Study 1 NHANES III', 'protocol': str(PROTO.relative_to(ROOT)), 'protocol_sha256': hashlib.sha256(ptxt).hexdigest(), 'protocol_sha256_before_addendum': hashlib.sha256(ptxt.split(b'\n## Addendum')[0]).hexdigest(), 'units': 'net benefit in true-positive-equivalent units over the risk horizon; loss fraction = regret / achievable net benefit (relative risk reduction cancels)'}
    P = pd.read_csv(ST / 'analysis_v1/stage3/stage3_participants.csv')
    z25 = zones('0.25')
    CV = json.loads((PH / 'cvd75_curve.json').read_text(encoding='utf-8'))
    E10 = json.loads((PH / 'ext4_30yr_results.json').read_text(encoding='utf-8'))
    BANDS3 = [(30, 39), (40, 59), (60, 79)]
    s1 = (P.statin == 0).to_numpy()
    agep = P.age.to_numpy(float)
    L1, L2 = (P.lpB_m0.to_numpy(float)[s1], P.lpB_m15.to_numpy(float)[s1])
    assert np.max(np.abs(100 * expit(L1) - P.riskB_m0.to_numpy()[s1])) < 1e-09
    out['ascvd10'] = analyse(L1, L2, P.cluster.to_numpy(), s1, [(0.05, z25[1]), (0.03, z25[0]), (0.1, z25[2])], 'base PREVENT 10-year ASCVD, statin-naive', sens=True, age=agep[s1], bands=BANDS3)
    sub = s1 & (agep >= 30) & (agep <= 59)
    gS, pS, dS, sS, iS = fit_g_sq(P.lpB_m0.to_numpy(float)[sub], P.lpB_m15.to_numpy(float)[sub])
    out['ascvd10']['post_audit_sens_5pct_ages_30_59'] = {'n': int(sub.sum()), 'converged': iS['converged'], **metrics(gS, pS, dS, sS, 0.05)}
    s2 = (P.bptreat == 0).to_numpy()
    T1, T2 = (logit(P.riskT1.to_numpy(float)[s2] / 100), logit(P.riskT2.to_numpy(float)[s2] / 100))
    out['cvd10_bp'] = analyse(T1, T2, P.cluster.to_numpy(), s2, [(0.075, tuple(CV['zones_risk_pct']['0.25']))], 'PREVENT 10-year total CVD, not on BP medication (all BP levels)', age=agep[s2], bands=BANDS3)
    E, B, T = nh3_build.build()
    F = fixed(B)
    age = B.age.to_numpy(float)
    q1, q2 = (risk30(F, visit_inputs(B, 1)), risk30(F, visit_inputs(B, 2)))
    a1 = (age >= 30) & (age <= 59) & (B.statin.to_numpy() != 1)
    assert a1.sum() == E10['part_A']['n_A1']
    out['ascvd30'] = analyse(logit(q1[a1] / 100), logit(q2[a1] / 100), B.cluster.to_numpy(), a1, [(0.1, tuple(E10['part_A']['pi30']['zones_risk30_pct']['0.25']))], 'base PREVENT 30-year ASCVD, ages 30-59 statin-naive', age=age[a1], bands=[(30, 39), (40, 59)])
    out['limits'] = 'assumes calibration and that the guideline threshold is the benefit-harm exchange rate; risk-threshold components only (BP stage, LDL-C 160/190 and HbA1c 6.5% components not covered); usual level = average over the about 2-week process (variance-inflation sensitivity for a longer-term level); loss fractions are specific to this population; 1988-1994 measurement; model-based'
    OUT.write_text(json.dumps(out, indent=1, default=lambda o: o.item() if hasattr(o, 'item') else o.tolist() if hasattr(o, 'tolist') else str(o)) + '\n', encoding='utf-8')
    print('saved', OUT)
if __name__ == '__main__':
    main()
