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
from nh3_stage3 import fixed, lp_of, visit_inputs, SEED, B_BOOT, CUTS
from nh3_reclassification_lookup import refit, design, dist_near, LT, R3
from nh3_ext2021_burden import zones as frozen_zones
from nh3_ext3_statin import risk30
import nh3_ext5_stakes as E11
ROOT = Path(__file__).resolve().parents[1]
ST = ROOT / 'studies/nhanes3_testretest'
PH = ST / 'analysis_v1_posthoc'
PROTO = ST / 'extensions6_protocol_20261011.md'
OUT = PH / 'ext6_review_results.json'
FINE = np.linspace(0, 1.5, 1501)
MIDS = [(LT[0] + LT[1]) / 2, (LT[1] + LT[2]) / 2]
EFLM_CVI = {'tc': 5.4, 'hdl': 5.7, 'scr': 4.4, 'h': 1.2}
CVA_OBS = {'tc': 1.33, 'hdl': 2.49, 'scr': 2.7, 'h': 2.0}

def pct(v):
    return S.pct_ci(v) if len(v) else None

def dstar_matrix(pf, lev):
    below = pf < lev
    exists = ~below[0]
    first = np.argmax(below, axis=0)
    never = ~below.any(axis=0)
    idx = np.where(never, len(FINE) - 1, first - 1)
    ds = FINE[np.clip(idx, 0, len(FINE) - 1)].astype(float)
    ds[~exists] = np.nan
    return ds

def truncate(ti, lo, hi):
    if ti == 0:
        hi = np.minimum(hi, MIDS[0])
    if ti == 1:
        lo = np.maximum(lo, MIDS[0])
        hi = np.minimum(hi, MIDS[1])
    if ti == 2:
        lo = np.maximum(lo, MIDS[1])
    return (lo, hi)

def zone_intervals(P):
    cl = P.cluster.to_numpy()
    lp1 = P.lpB_m0.to_numpy(float)
    ch = P.changeB.to_numpy(float)
    r1 = P.riskB_m0.to_numpy(float)
    W = S.boot_weights(cl, B_BOOT, SEED)
    knots, beta, boots, ind, worst = refit(lp1, ch, W, R3['K3_base']['clin'])
    frozen = R3['K3_base']['clin']['by_threshold']
    out = {'reproduction_max_abs_diff_vs_frozen': worst, 'thresholds': {}}
    allb = np.vstack([beta[None, :], boots])
    inz = {lev: np.zeros((len(allb), len(r1)), bool) for lev in (0.25, 0.1)}
    for ti, t in enumerate(CUTS):
        Xf = design(FINE, np.full(len(FINE), ti), knots, ind)
        pf = expit(Xf @ allb.T)
        e = {'peak_pi_at_threshold': float(pf[0, 0]), 'peak_ci': pct(pf[0, 1:])}
        assert abs(e['peak_pi_at_threshold'] - frozen[str(t)]['pi_at_d']['0.0']) < 1e-09
        for lev in (0.25, 0.1):
            ds = dstar_matrix(pf, lev)
            lo, hi = (LT[ti] - ds, LT[ti] + ds)
            unt = [float(100 * expit(lo[0])), float(100 * expit(hi[0]))]
            fz = frozen[str(t)]['instability_zones_risk_pct'][str(lev)]
            assert fz is not None and np.allclose(unt, fz, atol=1e-09), (t, lev, unt, fz)
            lo, hi = truncate(ti, lo, hi)
            rlo, rhi = (100 * expit(lo), 100 * expit(hi))
            ok = np.isfinite(ds)
            for b in range(len(allb)):
                if ok[b]:
                    inz[lev][b] |= (r1 >= rlo[b]) & (r1 <= rhi[b])
            e[str(lev)] = {'point': [float(rlo[0]), float(rhi[0])], 'untruncated_point': unt, 'lower_limit_ci': pct(rlo[1:]), 'upper_limit_ci': pct(rhi[1:]), 'width_pp_ci': pct((rhi - rlo)[1:]), 'width_pp': float(rhi[0] - rlo[0]), 'dstar_logit': float(ds[0]), 'dstar_logit_ci': pct(ds[1:]), 'width_logit': float(2 * ds[0]), 'width_logit_ci': pct(2 * ds[1:]), 'lower_limit_at_truncation_share': float(np.mean(np.isclose(lo[1:][ok[1:]], MIDS[0] if ti == 1 else -np.inf))) if ti == 1 else None, 'empty_replicates': int((~ok[1:]).sum()), 'replicates': int(len(allb) - 1), 'note': 'limit, width and d* intervals are percentile intervals among replicates with a non-empty zone; limits truncated at the logit midpoints'}
        out['thresholds'][str(t)] = e
    for lev, key in [(0.25, '0.25'), (0.1, '0.1')]:
        fz = frozen_zones(key)
        for ti, t in enumerate(CUTS):
            assert np.allclose(out['thresholds'][str(t)][str(lev)]['point'], fz[ti], atol=1e-09), (t, lev)
    for lev in (0.25, 0.1):
        fixed_share = float(inz[lev][0].mean())
        rep = [float(np.sum(W[b] * inz[lev][b + 1]) / np.sum(W[b])) for b in range(B_BOOT)]
        samp = [float(np.sum(W[b] * inz[lev][0]) / np.sum(W[b])) for b in range(B_BOOT)]
        out[f'share_in_union_zone_{lev}'] = {'fixed_zone_share': fixed_share, 'ci_sampling_only': pct(samp), 'ci_zone_and_sampling': pct(rep), 'n': int(len(r1)), 'in_zone_n': int(inz[lev][0].sum())}
    CV = json.loads((PH / 'cvd75_curve.json').read_text(encoding='utf-8'))
    rT1, rT2 = (P.riskT1.to_numpy(float), P.riskT2.to_numpy(float))
    x = ((rT1 >= 7.5) != (rT2 >= 7.5)).astype(float)
    d = np.abs(logit(rT1 / 100) - logit(0.075))
    kn = np.percentile(d, [10, 50, 90])
    X = lambda dv: np.column_stack([np.ones(len(dv)), S.rcs_basis(dv, kn)])
    b0 = S.logit_fit(X(d), x)
    assert np.allclose(b0, CV['beta'], atol=1e-09) and np.allclose(kn, CV['knots'])
    bb = np.array([S.logit_fit(X(d), x, W[b]) for b in range(B_BOOT)])
    allc = np.vstack([b0[None, :], bb])
    pf = expit(X(FINE) @ allc.T)
    c = {'events': int(x.sum()), 'peak_pi_at_threshold': float(pf[0, 0]), 'peak_ci': pct(pf[0, 1:])}
    for lev in (0.25, 0.1):
        ds = dstar_matrix(pf, lev)
        rlo, rhi = (100 * expit(logit(0.075) - ds), 100 * expit(logit(0.075) + ds))
        fz = CV['zones_risk_pct'].get(str(lev))
        assert fz is not None and np.allclose([rlo[0], rhi[0]], fz, atol=1e-09), (lev, rlo[0], rhi[0], fz)
        ok = np.isfinite(ds)
        c[str(lev)] = {'point': [float(rlo[0]), float(rhi[0])], 'lower_limit_ci': pct(rlo[1:]), 'upper_limit_ci': pct(rhi[1:]), 'width_pp': float(rhi[0] - rlo[0]), 'width_pp_ci': pct((rhi - rlo)[1:]), 'empty_replicates': int((~ok[1:]).sum()), 'replicates': int(B_BOOT)}
    out['cvd75_curve'] = c
    return out

def zone_regret(g, p, delta, sigma, pt, zone, mc, seed):
    t = logit(pt)
    r = expit(g)
    stake = np.abs(r - pt) / (1 - pt)
    rng = np.random.default_rng(seed)
    i = rng.choice(len(g), mc, p=p / p.sum())
    T = g[i]
    s = rng.choice([-1.0, 1.0], mc)
    L1 = T + s * delta / 2 + rng.normal(0, sigma, mc)
    r1 = 100 * expit(L1)
    inz = (r1 >= zone[0]) & (r1 <= zone[1])
    u = T >= t
    dis = (L1 >= t) != u
    reg = stake[i] * dis
    tot = reg.sum()
    return {'share_of_single_regret_from_zone': float(reg[inz].sum() / tot) if tot > 0 else None, 'fraction_in_zone': float(inz.mean()), 'regret_per_1000_in_zone': float(1000 * reg[inz].mean()) if inz.any() else None, 'regret_per_1000_outside_zone': float(1000 * reg[~inz].mean()), 'disagreement_in_zone': float(dis[inz].mean()) if inz.any() else None, 'disagreement_outside_zone': float(dis[~inz].mean()), 'regret_single_per_1000_mc': float(1000 * reg.mean())}

def exchange_rate(g, p, delta, sigma, pt, pe):
    t = logit(pt)
    r = expit(g)
    usual = (g >= t).astype(float)
    ps = 0.5 * (1 - norm.cdf((t - g + delta / 2) / sigma) + (1 - norm.cdf((t - g - delta / 2) / sigma)))
    w = (r - pe) / (1 - pe)
    nb_u = float(np.sum(p * w * usual))
    nb_s = float(np.sum(p * w * ps))
    nb_opt = float(np.sum(p * np.maximum(w, 0)))
    return {'decision_threshold': pt, 'exchange_rate': pe, 'nb_usual_rule_per_1000': 1000 * nb_u, 'nb_single_visit_per_1000': 1000 * nb_s, 'loss_per_1000': 1000 * (nb_u - nb_s), 'loss_fraction_of_usual_rule_nb': (nb_u - nb_s) / nb_u if nb_u > 0 else None, 'achievable_optimal_at_exchange_rate_per_1000': 1000 * nb_opt, 'loss_fraction_of_achievable_at_exchange_rate': (nb_u - nb_s) / nb_opt if nb_opt > 0 else None}

def stakes(P, B):
    z25 = frozen_zones('0.25')
    CV = json.loads((PH / 'cvd75_curve.json').read_text(encoding='utf-8'))
    E10 = json.loads((PH / 'ext4_30yr_results.json').read_text(encoding='utf-8'))
    agep = P.age.to_numpy(float)
    cl_full = P.cluster.to_numpy()
    s1 = (P.statin == 0).to_numpy()
    s2 = (P.bptreat == 0).to_numpy()
    F = fixed(B)
    age = B.age.to_numpy(float)
    q1, q2 = (risk30(F, visit_inputs(B, 1)), risk30(F, visit_inputs(B, 2)))
    a1 = (age >= 30) & (age <= 59) & (B.statin.to_numpy() != 1)
    assert a1.sum() == E10['part_A']['n_A1']
    pops = {'ascvd10': (P.lpB_m0.to_numpy(float)[s1], P.lpB_m15.to_numpy(float)[s1], cl_full, s1, [(0.05, z25[1], [0.03, 0.075]), (0.03, z25[0], [0.02, 0.05]), (0.1, z25[2], [0.075, 0.15])]), 'cvd10_bp': (logit(P.riskT1.to_numpy(float)[s2] / 100), logit(P.riskT2.to_numpy(float)[s2] / 100), cl_full, s2, [(0.075, tuple(CV['zones_risk_pct']['0.25']), [0.05, 0.1])]), 'ascvd30': (logit(q1[a1] / 100), logit(q2[a1] / 100), B.cluster.to_numpy(), a1, [(0.1, tuple(E10['part_A']['pi30']['zones_risk30_pct']['0.25']), [0.075, 0.15])])}
    X5 = json.loads((PH / 'ext5_stakes_results.json').read_text(encoding='utf-8'))
    out = {}
    for name, (L1, L2, clf, sel, specs) in pops.items():
        gkey = hashlib.sha256(b'e12fitg' + L1.tobytes() + L2.tobytes()).hexdigest()
        g, p, delta, sigma, info = E11._cached('e12_fitg', gkey, lambda: E11.fit_g(L1, L2))
        assert abs(sigma - X5[name]['sigma']) < 1e-12 and abs(delta - X5[name]['delta']) < 1e-12
        W = S.boot_weights(clf, 4000, E11.SEED_BOOT)[:E11.N_BOOT][:, sel]
        key = hashlib.sha256(L1.tobytes() + L2.tobytes() + W.tobytes()).hexdigest()
        f = PH / '_cache_e11' / f'uniform_{key[:16]}.pkl'
        assert f.exists(), f'E11 uniform-start cache missing for {name}'
        ufits = E11._cached('uniform', key, lambda: (_ for _ in ()).throw(RuntimeError('cache missing')))
        res = {'n': int(len(L1)), 'sigma': sigma, 'delta': delta, 'thresholds': {}, 'ci_rule': 'percentile intervals over the 1,000 uniform-start replicates; a ratio interval is reported only when the point and every replicate are defined (audit item 2)'}

        def ci_of(vals, point):
            und = int(sum((v is None for v in vals)))
            return (None if point is None or und else pct([v for v in vals]), und)
        for pt, zone, pes in specs:
            zr = zone_regret(g, p, delta, sigma, pt, zone, E11.MC, E11.SEED_MC)
            e5 = X5[name]['thresholds'][f'{100 * pt:g}']
            assert abs(zr['regret_single_per_1000_mc'] - e5['mc_check_regret_single_per_1000']) < 1e-09 and abs(zr['fraction_in_zone'] - e5['targeted_fraction_remeasured']) < 1e-12
            bz = [zone_regret(gb, pb, db, sb, pt, zone, E11.MC_B, E11.SEED_MC + b) for b, (gb, pb, db, sb, _) in enumerate(ufits)]
            cis = {k: ci_of([x[k] for x in bz], zr[k]) for k in zr}
            ent = {'zone_risk_pct': list(zone), **zr, 'ci_uniform_start': {k: v[0] for k, v in cis.items() if v[0] is not None}, 'undefined_replicates': {k: v[1] for k, v in cis.items() if v[1]}, 'e11_regret_single_per_1000': e5['regret_single_per_1000'], 'e11_achievable_nb_per_1000': e5['achievable_nb_per_1000'], 'e11_loss_fraction_single': e5['loss_fraction_single'], 'exchange_rate_sensitivity': []}
            for pe in pes:
                xr = exchange_rate(g, p, delta, sigma, pt, pe)
                bx = [exchange_rate(gb, pb, db, sb, pt, pe) for gb, pb, db, sb, _ in ufits]
                cx = {k: ci_of([x[k] for x in bx], xr[k]) for k in ['loss_per_1000', 'loss_fraction_of_usual_rule_nb', 'loss_fraction_of_achievable_at_exchange_rate']}
                xr['ci_uniform_start'] = {k: v[0] for k, v in cx.items() if v[0] is not None}
                xr['undefined_replicates'] = {k: v[1] for k, v in cx.items() if v[1]}
                xr['replicates_with_nonpositive_rule_nb'] = int(sum((x['nb_usual_rule_per_1000'] <= 0 for x in bx)))
                ent['exchange_rate_sensitivity'].append(xr)
            x0 = exchange_rate(g, p, delta, sigma, pt, pt)
            assert abs(x0['loss_per_1000'] - e5['regret_single_per_1000']) < 1e-09
            ent['exchange_rate_equal_to_threshold_check'] = x0
            res['thresholds'][f'{100 * pt:g}'] = ent
        out[name] = res
        print('stakes', name, 'done', flush=True)
    return out

def subsample(P, E, B):
    assert np.array_equal(P.SEQN.to_numpy(), B.SEQN.to_numpy())
    cl = P.cluster.to_numpy()
    ch = P.changeB.to_numpy(float)
    rr = B.rr_days.to_numpy(float)
    ok = np.isfinite(rr)
    out = {'interval_missing_n': int((~ok).sum()), 'interval_summary': {'min': float(np.nanmin(rr)), 'median': float(np.nanmedian(rr)), 'max': float(np.nanmax(rr))}}
    cuts = np.percentile(rr[ok], [100 / 3, 200 / 3])
    out['tertile_cutpoints_days'] = cuts.tolist()
    grp = np.digitize(rr, cuts, right=True)
    tert = []
    for k in range(3):
        s = ok & (grp == k)
        e = S.prop_ci(ch[s], cl[s])
        e.update(days_min=float(rr[s].min()), days_max=float(rr[s].max()), days_median=float(np.median(rr[s])))
        tert.append(e)
    out['theta_by_interval_tertile'] = tert
    if (~ok).any():
        out['theta_interval_missing'] = S.prop_ci(ch[~ok], cl[~ok])
    X = np.column_stack([np.ones(ok.sum()), rr[ok] / 10])
    y = ch[ok]
    b = S.logit_fit(X, y)
    pr = expit(X @ b)
    Wd = pr * (1 - pr)
    H = X.T @ (X * Wd[:, None])
    Hi = np.linalg.inv(H)
    sc = X * (y - pr)[:, None]
    ids, inv = np.unique(cl[ok], return_inverse=True)
    G = len(ids)
    Sg = np.zeros((G, X.shape[1]))
    np.add.at(Sg, inv, sc)
    meat = Sg.T @ Sg * G / (G - 1)
    V = Hi @ meat @ Hi
    se = math.sqrt(V[1, 1])
    q = tdist.ppf(0.975, G - 1)
    out['trend_per_10_days'] = {'log_or': float(b[1]), 'se_cr1': se, 'or': float(math.exp(b[1])), 'or_ci': [float(math.exp(b[1] - q * se)), float(math.exp(b[1] + q * se))], 'G': int(G), 'method': 'logistic regression of the change indicator on interval/10 days; CR1 cluster-robust (pseudo-PSU) variance; t(G-1)'}
    assert len(E) == 979 and int(E.complete.sum()) == 759
    comp = E.complete.to_numpy(bool)
    rows = {}

    def bin_(name, v):
        v = np.asarray(v, float)
        a, b_ = (v[comp], v[~comp])
        a, b_ = (a[np.isfinite(a)], b_[np.isfinite(b_)])
        pa, pb = (a.mean(), b_.mean())
        sd = math.sqrt((pa * (1 - pa) + pb * (1 - pb)) / 2)
        rows[name] = {'complete_pct': 100 * pa, 'incomplete_pct': 100 * pb, 'complete_n': int(len(a)), 'incomplete_n': int(len(b_)), 'std_diff': (pa - pb) / sd if sd > 0 else None}
    ages = E.age.to_numpy(float)
    a, b_ = (ages[comp], ages[~comp])
    rows['age_years'] = {'complete_mean': a.mean(), 'complete_sd': a.std(ddof=1), 'incomplete_mean': b_.mean(), 'incomplete_sd': b_.std(ddof=1), 'std_diff': (a.mean() - b_.mean()) / math.sqrt((a.var(ddof=1) + b_.var(ddof=1)) / 2)}
    bin_('women', E.female)
    dm = E.dmarethn.to_numpy(float)
    for code, lab in [(1, 'non_hispanic_white'), (2, 'non_hispanic_black'), (3, 'mexican_american'), (4, 'other')]:
        bin_(lab, dm == code)
    bin_('current_smoking', E.smoking)
    bin_('bp_treatment', E.bptreat)
    bin_('self_reported_diabetes', E.srd)
    st = E.statin.to_numpy(dtype=object)
    known = np.array([x is not None and x == x for x in st])
    stv = np.where(known, np.array([float(x) if x is not None and x == x else np.nan for x in st]), np.nan)
    bin_('statin_use_where_known', stv)
    rows['statin_unknown_n'] = {'complete': int((~known & comp).sum()), 'incomplete': int((~known & ~comp).sum())}
    out['complete_vs_incomplete'] = {'n_complete': int(comp.sum()), 'n_incomplete': int((~comp).sum()), 'rows': rows}
    out['published_values_note'] = 'Lacher 2005 and Selvin 2013 derived their within-person variation estimates from the same NHANES III second-examination subsample; agreement with them checks data processing, not representativeness.'
    return out

def scenarios(B, P):
    n = len(B)
    F = fixed(B)
    sbp1, sbp2 = (B.sbp1.to_numpy(float), B.sbp2.to_numpy(float))
    tc1, tc2, hd1, hd2 = (B[c].to_numpy(float) for c in ['tc1', 'tc2', 'hdl1', 'hdl2'])
    sc1, sc2 = (B.scr1.to_numpy(float), B.scr2.to_numpy(float))
    h1, h2 = (B.hba1c1.to_numpy(float), B.hba1c2.to_numpy(float))
    fx = B.fixed_dm.to_numpy(float)
    geo = lambda a, b: np.sqrt(a * b)
    c = {'sbp': (sbp1 + sbp2) / 2, 'tc': geo(tc1, tc2), 'hdl': geo(hd1, hd2), 'scr': geo(sc1, sc2), 'h': np.where(np.isfinite(h1) & np.isfinite(h2), (h1 + h2) / 2, np.where(np.isfinite(h1), h1, h2))}
    sd_obs = {'sbp': np.std(sbp2 - sbp1, ddof=1) / math.sqrt(2), 'tc': np.std(np.log(tc2 / tc1), ddof=1) / math.sqrt(2), 'hdl': np.std(np.log(hd2 / hd1), ddof=1) / math.sqrt(2), 'scr': np.std(np.log(sc2 / sc1), ddof=1) / math.sqrt(2), 'h': np.nanstd(h2 - h1, ddof=1) / math.sqrt(2)}
    RS = json.loads((PH / 'repeat_strategies.json').read_text(encoding='utf-8'))
    assert all((abs(sd_obs[k] - RS['sd_used'][k]) < 1e-12 for k in sd_obs))
    from nh3_figures import eflm_cvi
    ef = eflm_cvi()
    assert {'tc': ef['tc'], 'hdl': ef['hdl'], 'scr': ef['scr_cal'], 'h': ef['hba1c_log']} == EFLM_CVI, ef
    hmean = float(np.nanmean(c['h']))
    cv_mod = {k: math.sqrt(EFLM_CVI[k] ** 2 + CVA_OBS[k] ** 2) / 100 for k in EFLM_CVI}
    sd_mod = {k: math.sqrt(math.log(1 + cv_mod[k] ** 2)) for k in ['tc', 'hdl', 'scr']}
    sd_mod['h'] = cv_mod['h'] * hmean
    LAB = ['tc', 'hdl', 'scr', 'h']

    def category(zbar, sd):
        sb = c['sbp'] + zbar['sbp'] * sd['sbp']
        t_ = c['tc'] * np.exp(zbar['tc'] * sd['tc'])
        hd = c['hdl'] * np.exp(zbar['hdl'] * sd['hdl'])
        sc = c['scr'] * np.exp(zbar['scr'] * sd['scr'])
        hv = c['h'] + zbar['h'] * sd['h']
        dm = np.where((fx == 1) | (np.nan_to_num(hv, nan=0) >= 6.5), 1.0, 0.0)
        I = {'sbp': M.clamp(sb, 90, 200), 'tc': M.clamp(t_, 130, 320), 'hdl': M.clamp(hd, 20, 100), 'egfr': M.clamp(M.ckd_epi_2021(np.maximum(sc, 1e-06), F['age'], F['female']), 15, 140), 'dm': dm, 'uacr': None}
        r = M.risk_from_lp(lp_of(F, I, 'base'))
        return (M.category(r), r)
    R_SIM = 200
    grid = []
    for li, lab in enumerate(['observed', 'modern']):
        for mi, mult in enumerate([0.75, 1.0, 1.25, 1.5]):
            sd = dict(sd_obs) if lab == 'observed' else {**sd_mod, 'sbp': sd_obs['sbp']}
            sd['sbp'] = sd_obs['sbp'] * mult
            rng = np.random.default_rng(SEED + 200 + 4 * li + mi)
            dis = np.zeros(n)
            ra_all = np.empty((R_SIM, n))
            ch_all = np.empty((R_SIM, n))
            for r_ in range(R_SIM):
                za = {k: rng.standard_normal(n) for k in ['sbp'] + LAB}
                zb = {k: rng.standard_normal(n) for k in ['sbp'] + LAB}
                ca, ra = category(za, sd)
                cb, _ = category(zb, sd)
                ch_all[r_] = ca != cb
                ra_all[r_] = ra
                dis += ca != cb
            theta = float(dis.mean() / R_SIM)
            lp = logit(ra_all.ravel() / 100)
            ch = ch_all.ravel()
            d, near = dist_near(lp)
            knots = np.percentile(d, [10, 50, 90])
            X = design(d, near, knots, True)
            beta = S.logit_fit(X, ch)
            e = {'laboratory': lab, 'sbp_multiplier': mult, 'sd_used': {k: float(v) for k, v in sd.items()}, 'cvw_pct_implied': {k: 100 * math.sqrt(math.exp(sd[k] ** 2) - 1) for k in ['tc', 'hdl', 'scr']} | {'h_sd_units': sd['h'], 'sbp_mm_hg': sd['sbp']}, 'theta_sim': theta, 'rows_for_curve': int(len(ch)), 'thresholds': {}}
            inz = {0.25: np.zeros(n, bool), 0.1: np.zeros(n, bool)}
            r1 = P.riskB_m0.to_numpy(float)
            inz_sim = {0.25: np.zeros(ra_all.shape, bool), 0.1: np.zeros(ra_all.shape, bool)}
            for ti, t in enumerate(CUTS):
                pf = expit(design(FINE, np.full(len(FINE), ti), knots, True) @ beta)[:, None]
                et = {'peak': float(pf[0, 0])}
                for lev in (0.25, 0.1):
                    ds = dstar_matrix(pf, lev)[0]
                    if not np.isfinite(ds):
                        et[str(lev)] = None
                        continue
                    lo, hi = truncate(ti, np.array([LT[ti] - ds]), np.array([LT[ti] + ds]))
                    rlo, rhi = (float(100 * expit(lo[0])), float(100 * expit(hi[0])))
                    et[str(lev)] = [rlo, rhi]
                    inz[lev] |= (r1 >= rlo) & (r1 <= rhi)
                    inz_sim[lev] |= (ra_all >= rlo) & (ra_all <= rhi)
                e['thresholds'][str(t)] = et
            e['share_simulated_first_assessment_in_union_zone'] = {str(lev): float(inz_sim[lev].mean()) for lev in inz_sim}
            e['share_observed_first_visit_risk_in_union_zone'] = {str(lev): float(inz[lev].mean()) for lev in inz}
            grid.append(e)
            print('scenario', lab, mult, 'theta_sim', round(100 * theta, 2), flush=True)
    base = [g for g in grid if g['laboratory'] == 'observed' and g['sbp_multiplier'] == 1.0][0]
    assert abs(base['theta_sim'] - R3['S12']['A_base']['theta_sim']) < 0.01, 'model A reference not reproduced'
    return {'model_A_frozen_theta_sim': R3['S12']['A_base']['theta_sim'], 'observed_theta': R3['primary']['estimate'], 'sd_observed': {k: float(v) for k, v in sd_obs.items()}, 'modern_laboratory_cvw_pct': {k: 100 * v for k, v in cv_mod.items()}, 'eflm_cvi_pct': EFLM_CVI, 'cva_pct': CVA_OBS, 'hba1c_mean_used': hmean, 'R_SIM': R_SIM, 'zones_from_simulated_first_assessment': True, 'grid': grid, 'cva_provenance': 'frozen NHANES III QC transcription: TC 1.33 (Q17 pool; Q20 pool 1.25), HDL 2.49 (AQ13 pool), creatinine 2.7 (1.7 mg/dL pool), HbA1c 2.0 (documentation inter-assay note; manual pools 1.57-1.73)', 'note': "model-based; independent errors around each participant's two-visit mean; modern laboratory level = sqrt(EFLM CVI^2 + observed NHANES III CVa^2) (research-protocol biological variation plus analytical imprecision; routine pre-analytical variation not included); SBP multipliers span plausible device and protocol effects and are not estimates of current practice; zones fitted to the simulated first assessment with the prespecified spline form (no intervals: scenario values)"}

def main():
    ptxt = PROTO.read_bytes()
    out = {'label': 'E12: analyses prompted by the independent referee review of revision 7 draft 1 (zone intervals; per-zone and exchange-rate stakes; subsample checks; error-variance scenarios)', 'protocol': str(PROTO.relative_to(ROOT)), 'protocol_sha256': hashlib.sha256(ptxt).hexdigest(), 'protocol_sha256_before_addendum': hashlib.sha256(ptxt.split(b'\n## Addendum')[0]).hexdigest()}
    P = pd.read_csv(ST / 'analysis_v1/stage3/stage3_participants.csv')
    E, B, T = nh3_build.build()
    out['A_zone_intervals'] = zone_intervals(P)
    print('A done', flush=True)
    out['C_subsample'] = subsample(P, E, B)
    print('C done', flush=True)
    out['D_error_scenarios'] = scenarios(B, P)
    print('D done', flush=True)
    out['B_stakes'] = stakes(P, B)
    print('B done', flush=True)
    OUT.write_text(json.dumps(out, indent=1, default=lambda o: o.item() if hasattr(o, 'item') else o.tolist() if hasattr(o, 'tolist') else str(o)) + '\n', encoding='utf-8')
    print('saved', OUT)
if __name__ == '__main__':
    main()
