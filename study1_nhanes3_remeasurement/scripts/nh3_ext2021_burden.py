"""E2 (post hoc, exploratory, 2026-10-10; extensions_protocol_20261010.md): transport of the Study 1 threshold-distance curves to US adults"""
from pathlib import Path
import json
import math
import numpy as np
import pandas as pd
import pyreadstat
from scipy.special import expit, logit
from scipy.stats import t as tdist
import nh3_models as M
import nh3_stats as S
from nh3_reclassification_lookup import refit, design, dist_near, CUTS, B_BOOT, SEED
ROOT = Path(__file__).resolve().parents[1]
ST = ROOT / 'studies/nhanes3_testretest'
D = ST / 'data_2021_2023/raw'
OUT = ST / 'analysis_v1_posthoc/ext2021_burden.json'
R3 = json.loads((ST / 'analysis_v1/stage3/stage3_results.json').read_text(encoding='utf-8'))
CV = json.loads((ST / 'analysis_v1_posthoc/cvd75_curve.json').read_text(encoding='utf-8'))

def xpt(name, cols):
    df, _ = pyreadstat.read_xport(str(D / f'{name}.xpt'), usecols=['SEQN'] + cols, encoding='latin1')
    return df

def zones(level):
    K3 = R3['K3_base']['clin']['by_threshold']
    lt = logit(np.array(CUTS) / 100)
    b35, b510 = (100 * expit((lt[0] + lt[1]) / 2), 100 * expit((lt[1] + lt[2]) / 2))
    z = {t: list(K3[t]['instability_zones_risk_pct'][level]) for t in ['3.0', '5.0', '10.0']}
    z['3.0'][1] = min(z['3.0'][1], b35)
    z['5.0'][0] = max(z['5.0'][0], b35)
    z['5.0'][1] = min(z['5.0'][1], b510)
    z['10.0'][0] = max(z['10.0'][0], b510)
    return [tuple(v) for v in z.values()]

class Svy:

    def __init__(self, w, st, ps):
        self.w, self.st, self.ps = (np.asarray(w, float), np.asarray(st), np.asarray(st) * 10 + np.asarray(ps))
        self.df = len(np.unique(self.ps)) - len(np.unique(self.st))

    def _var_total(self, z):
        v = 0.0
        for h in np.unique(self.st):
            m = self.st == h
            ids = np.unique(self.ps[m])
            tots = np.array([z[m & (self.ps == p)].sum() for p in ids])
            if len(ids) > 1:
                v += len(ids) / (len(ids) - 1) * np.sum((tots - tots.mean()) ** 2)
        return v

    def mean(self, y, dom=None, logit_ci=True):
        y = np.asarray(y, float)
        dom = np.ones(len(y), bool) if dom is None else np.asarray(dom, bool)
        wd = self.w * dom
        th = np.sum(wd * y) / np.sum(wd)
        z = wd * (y - th) / np.sum(wd)
        se = math.sqrt(self._var_total(z))
        q = tdist.ppf(0.975, self.df)
        if logit_ci and 0 < th < 1:
            s = se / (th * (1 - th))
            lo, hi = (expit(logit(th) - q * s), expit(logit(th) + q * s))
        else:
            lo, hi = (th - q * se, th + q * se)
        return {'estimate': float(th), 'lower': float(lo), 'upper': float(hi), 'se': se, 'n': int(dom.sum()), 'events': int(np.sum(dom * (y > 0.5))) if set(np.unique(y)) <= {0.0, 1.0} else None}

    def total(self, y, dom=None):
        y = np.asarray(y, float)
        dom = np.ones(len(y), bool) if dom is None else np.asarray(dom, bool)
        z = self.w * dom * y
        T = float(z.sum())
        se = math.sqrt(self._var_total(z))
        q = tdist.ppf(0.975, self.df)
        return {'millions': T / 1000000.0, 'lower': (T - q * se) / 1000000.0, 'upper': (T + q * se) / 1000000.0}

def main():
    demo = xpt('DEMO_L', ['RIDSTATR', 'RIAGENDR', 'RIDAGEYR', 'RIDEXPRG', 'WTMEC2YR', 'SDMVSTRA', 'SDMVPSU'])
    d = demo
    for name, cols in [('MCQ_L', ['MCQ160B', 'MCQ160C', 'MCQ160D', 'MCQ160E', 'MCQ160F']), ('BPXO_L', ['BPXOSY1', 'BPXOSY2', 'BPXOSY3', 'BPXODI1', 'BPXODI2', 'BPXODI3']), ('TCHOL_L', ['LBXTC', 'WTPH2YR']), ('HDL_L', ['LBDHDD']), ('BIOPRO_L', ['LBXSCR']), ('GHB_L', ['LBXGH']), ('ALB_CR_L', ['URDACT']), ('BPQ_L', ['BPQ020', 'BPQ150', 'BPQ101D']), ('DIQ_L', ['DIQ010', 'DIQ050', 'DIQ070']), ('SMQ_L', ['SMQ020', 'SMQ040'])]:
        d = d.merge(xpt(name, cols), on='SEQN', how='left')
    flow = [('DEMO_L records', len(d))]
    d = d[(d.RIDSTATR == 2) & (d.WTPH2YR > 0)]
    flow.append(('Examined with positive phlebotomy weight', len(d)))
    d = d[(d.RIDAGEYR >= 30) & (d.RIDAGEYR <= 79)]
    flow.append(('Age 30-79 years', len(d)))
    d = d[(d[['MCQ160B', 'MCQ160C', 'MCQ160D', 'MCQ160E', 'MCQ160F']] == 2).all(axis=1)]
    flow.append(('Explicitly no CHF, CHD, angina, MI or stroke', len(d)))
    d = d[d.RIDEXPRG != 1]
    flow.append(('Not pregnant', len(d)))
    d = d.copy()
    w_elig = float(d.WTPH2YR.sum())
    sy = d[['BPXOSY2', 'BPXOSY3']].mean(axis=1, skipna=False).fillna(d[['BPXOSY1', 'BPXOSY2', 'BPXOSY3']].mean(axis=1))
    dy = d[['BPXODI2', 'BPXODI3']].mean(axis=1, skipna=False).fillna(d[['BPXODI1', 'BPXODI2', 'BPXODI3']].mean(axis=1))
    d['sbp'], d['dbp'] = (sy, dy)
    qpos = (d.DIQ010 == 1) | (d.DIQ050 == 1) | (d.DIQ070 == 1)
    d['dm'] = np.where(qpos | (d.LBXGH >= 6.5), 1.0, np.where(d.LBXGH.notna() & d.DIQ010.isin([2, 3]), 0.0, np.nan))
    d['smoking'] = np.where((d.SMQ020 == 1) & d.SMQ040.isin([1, 2]), 1.0, np.where((d.SMQ020 == 2) | (d.SMQ020 == 1) & (d.SMQ040 == 3), 0.0, np.nan))
    d['bptreat'] = np.where(d.BPQ150 == 1, 1.0, np.where((d.BPQ150 == 2) | (d.BPQ020 == 2), 0.0, np.nan))
    d['statin'] = np.where(d.BPQ101D == 1, 1.0, np.where(d.BPQ101D == 2, 0.0, np.nan))
    d['female'] = (d.RIAGENDR == 2).astype(float)
    need = ['sbp', 'LBXTC', 'LBDHDD', 'LBXSCR', 'dm', 'smoking', 'bptreat', 'statin']
    d = d[d[need].notna().all(axis=1)].copy()
    flow.append(('Complete base PREVENT inputs', len(d)))
    eg = M.ckd_epi_2021(d.LBXSCR.to_numpy(float), d.RIDAGEYR.to_numpy(float), d.female.to_numpy(float))
    args = dict(female=d.female.to_numpy(float), age=d.RIDAGEYR.to_numpy(float), tc=M.clamp(d.LBXTC, 130, 320), hdl=M.clamp(d.LBDHDD, 20, 100), sbp=M.clamp(d.sbp, 90, 200), dm=d.dm.to_numpy(float), smoking=d.smoking.to_numpy(float), egfr=M.clamp(eg, 15, 140), bptreat=d.bptreat.to_numpy(float), statin=d.statin.to_numpy(float))
    lpB = M.lp(args['female'], args['age'], args['tc'], args['hdl'], args['sbp'], args['dm'], args['smoking'], args['egfr'], args['bptreat'], args['statin'], model='base')
    lpT = M.lp(args['female'], args['age'], args['tc'], args['hdl'], args['sbp'], args['dm'], args['smoking'], args['egfr'], args['bptreat'], args['statin'], model='base', outcome='total_cvd')
    rB, rT = (M.risk_from_lp(lpB), M.risk_from_lp(lpT))
    sv = Svy(d.WTPH2YR, d.SDMVSTRA, d.SDMVPSU)
    cat = M.category(rB)
    scale = w_elig / float(d.WTPH2YR.sum())
    res = {'label': __doc__.split('\n')[0], 'flow': flow, 'df': sv.df, 'weight': 'WTPH2YR (phlebotomy)', 'population_millions': float(d.WTPH2YR.sum() / 1000000.0), 'eligible_population_millions': w_elig / 1000000.0, 'completeness_scale': scale, 'note_millions': 'millions are complete-case weighted totals; *_scaled multiply by the eligible/complete weight ratio, assuming missingness unrelated to the outcome', 'oscillometric_sbp_mean': sv.mean(d.sbp.to_numpy(float), logit_ci=False), 'median_risk_unweighted': float(np.median(rB)), 'categories': {lab: {**sv.mean(cat == k), **{'millions': sv.total(cat == k)}} for k, lab in enumerate(['<3', '3-<5', '5-<10', '>=10'])}}
    inz = lambda r, Z: np.any([(r >= a) & (r <= b) for a, b in Z], axis=0)
    for lev in ('0.25', '0.1'):
        m = inz(rB, zones(lev)).astype(float)
        res[f'zone_{lev}'] = {**sv.mean(m), 'millions': sv.total(m)}
    P = pd.read_csv(ST / 'analysis_v1/stage3/stage3_participants.csv')
    W = S.boot_weights(P.cluster.to_numpy(), B_BOOT, SEED)
    chB = (np.searchsorted(CUTS, P.riskB_m15, side='right') != np.searchsorted(CUTS, P.riskB_m0, side='right')).astype(float)
    knots, beta, boots, ind, _ = refit(P.lpB_m0.to_numpy(float), chB, W, R3['K3_base']['clin'])
    dd, near = dist_near(lpB)
    X = design(dd, near, knots, ind)
    pi = expit(X @ beta)
    th = sv.mean(pi, logit_ci=False)
    wn = d.WTPH2YR.to_numpy(float) / d.WTPH2YR.sum()
    curve = np.array([wn @ expit(X @ b) for b in boots])
    sd_tot = math.sqrt(np.var(curve, ddof=1) + th['se'] ** 2)
    res['transported_theta'] = {'estimate': th['estimate'], 'lower': th['estimate'] - 1.96 * sd_tot, 'upper': th['estimate'] + 1.96 * sd_tot, 'curve_percentile': S.pct_ci(curve), 'survey_se': th['se'], 'millions': float(th['estimate'] * d.WTPH2YR.sum() / 1000000.0)}
    res['nhanes3_standardized_theta'] = R3['K4_base']['theta_std']
    res['e14_transport_2009_2010'] = R3['E14']['transported_theta']
    s1 = (d.bptreat == 0) & ((d.sbp >= 130) & (d.sbp < 140) & (d.dbp < 90) | (d.dbp >= 80) & (d.dbp < 90) & (d.sbp < 140))
    s1 = s1.to_numpy() & np.isfinite(d.dbp.to_numpy(float))
    dbp_only = s1 & ~((d.sbp >= 130) & (d.sbp < 140)).to_numpy()
    xc = P.cvd_cross75.to_numpy(float)
    dc = np.abs(logit(P.riskT1.to_numpy(float) / 100) - logit(0.075))
    Xc = np.column_stack([np.ones(len(dc)), S.rcs_basis(dc, np.array(CV['knots']))])
    bc0 = S.logit_fit(Xc, xc)
    assert np.allclose(bc0, CV['beta'], atol=1e-09)
    bcb = np.array([S.logit_fit(Xc, xc, W[b]) for b in range(W.shape[0])])
    z10 = CV['zones_risk_pct']['0.1']
    kn, bc = (np.array(CV['knots']), np.array(CV['beta']))
    dT = np.abs(logit(rT / 100) - logit(0.075))
    piT = expit(np.column_stack([np.ones(len(dT)), S.rcs_basis(dT, kn)]) @ bc)
    res['stage1_untreated'] = {'n': int(s1.sum()), 'millions': sv.total(s1.astype(float)), 'cvd_ge_7_5': {**sv.mean((rT >= 7.5).astype(float), s1), 'millions': sv.total((rT >= 7.5).astype(float), s1)}, 'cvd_in_10pct_zone': {**sv.mean(((rT >= z10[0]) & (rT <= z10[1])).astype(float), s1), 'millions': sv.total(((rT >= z10[0]) & (rT <= z10[1])).astype(float), s1)}, 'expected_crossing_7_5': sv.mean(piT, s1, logit_ci=False), 'zone_10pct_risk_pct': z10, 'dbp_only_n': int(dbp_only.sum())}
    XT = np.column_stack([np.ones(len(dT)), S.rcs_basis(dT, kn)])
    ws1 = d.WTPH2YR.to_numpy(float) * s1
    ws1 = ws1 / ws1.sum()
    cc = np.array([ws1 @ expit(XT @ b) for b in bcb])
    ec = res['stage1_untreated']['expected_crossing_7_5']
    sdc = math.sqrt(np.var(cc, ddof=1) + ec['se'] ** 2)
    ec.update(lower=ec['estimate'] - 1.96 * sdc, upper=ec['estimate'] + 1.96 * sdc, note='curve bootstrap and survey variance combined')
    res['zone_intervals_note'] = 'zone limits fixed at the Study 1 point estimates; zone intervals reflect survey sampling only'
    for key in ['zone_0.25', 'zone_0.1']:
        res[key]['millions_scaled'] = {k: v * scale for k, v in res[key]['millions'].items()}
    res['transported_theta']['millions_scaled'] = res['transported_theta']['millions'] * scale
    for key in ['cvd_ge_7_5', 'cvd_in_10pct_zone']:
        res['stage1_untreated'][key]['millions_scaled'] = {k: v * scale for k, v in res['stage1_untreated'][key]['millions'].items()}
    res['stage1_untreated']['millions_scaled'] = {k: v * scale for k, v in res['stage1_untreated']['millions'].items()}
    print('eligible millions', round(w_elig / 1000000.0, 1), 'scale', round(scale, 3), 'zone25 scaled', round(res['zone_0.25']['millions_scaled']['millions'], 1), 'zone10 scaled', round(res['zone_0.1']['millions_scaled']['millions'], 1), 'theta scaled', round(res['transported_theta']['millions_scaled'], 1), 'stage1 scaled', round(res['stage1_untreated']['millions_scaled']['millions'], 1), 'stage1 zone scaled', round(res['stage1_untreated']['cvd_in_10pct_zone']['millions_scaled']['millions'], 2))
    res['assumptions'] = 'NHANES III (1988-1994) within-person variation and distance-only dependence transported; oscillometric BP; statin proxy BPQ101D'
    OUT.write_text(json.dumps(res, indent=1, default=lambda o: o.item() if hasattr(o, 'item') else str(o)) + '\n', encoding='utf-8')
    print('flow', flow)
    print('pop millions', round(res['population_millions'], 1), 'df', sv.df)
    for k, v in res['categories'].items():
        print(k, round(100 * v['estimate'], 1), round(v['millions']['millions'], 1))
    for lev in ('0.25', '0.1'):
        z = res[f'zone_{lev}']
        print('zone', lev, round(100 * z['estimate'], 1), (round(100 * z['lower'], 1), round(100 * z['upper'], 1)), round(z['millions']['millions'], 1), (round(z['millions']['lower'], 1), round(z['millions']['upper'], 1)))
    t = res['transported_theta']
    print('transported', round(100 * t['estimate'], 2), round(100 * t['lower'], 2), round(100 * t['upper'], 2), 'millions', round(t['millions'], 1))
    s = res['stage1_untreated']
    print('stage1', s['n'], round(s['millions']['millions'], 1), 'cvd>=7.5', round(100 * s['cvd_ge_7_5']['estimate'], 1), round(s['cvd_ge_7_5']['millions']['millions'], 2), 'zone10', round(100 * s['cvd_in_10pct_zone']['estimate'], 1), round(s['cvd_in_10pct_zone']['millions']['millions'], 2), (round(s['cvd_in_10pct_zone']['millions']['lower'], 2), round(s['cvd_in_10pct_zone']['millions']['upper'], 2)), 'exp cross', round(100 * s['expected_crossing_7_5']['estimate'], 2))
if __name__ == '__main__':
    main()
