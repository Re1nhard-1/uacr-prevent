from pathlib import Path
import json
import math
import numpy as np
import pandas as pd
import nh3_build
import nh3_models as M
from nh3_testretest_stage12 import num, code, sbp23, dm_gate, egfr, cal, uacr
ROOT = Path(__file__).resolve().parents[1]
STUDY = ROOT / 'studies/nhanes3_testretest'
OUT = STUDY / 'manuscript_evidence_v1'
BOUNDS = {'sbp': (90, 200), 'tc': (130, 320), 'hdl': (20, 100), 'egfr': (15, 140)}

def wq(x, w, q):
    x, w = (np.asarray(x, float), np.asarray(w, float))
    o = np.argsort(x)
    x, w = (x[o], w[o])
    c = (np.cumsum(w) - 0.5 * w) / w.sum()
    return float(np.interp(q, c, x))

def summ(x, w=None, kind='mean'):
    x = np.asarray(x, float)
    m = np.isfinite(x)
    x = x[m]
    w = np.ones_like(x) if w is None else np.asarray(w, float)[m]
    if kind == 'pct':
        p = float(np.dot(w, x) / w.sum())
        return {'n': int(m.sum()), 'count': int(x.sum()), 'pct': 100 * p, 'p': p}
    mu = float(np.dot(w, x) / w.sum())
    sd = float(math.sqrt(np.dot(w, (x - mu) ** 2) / w.sum()))
    return {'n': int(m.sum()), 'mean': mu, 'sd': sd, 'median': wq(x, w, 0.5), 'q1': wq(x, w, 0.25), 'q3': wq(x, w, 0.75)}

def smd(a, b, kind):
    if a is None or b is None:
        return None
    if kind == 'pct':
        p1, p2 = (a['p'], b['p'])
        den = math.sqrt((p1 * (1 - p1) + p2 * (1 - p2)) / 2)
        return None if den == 0 else (p1 - p2) / den
    den = math.sqrt((a['sd'] ** 2 + b['sd'] ** 2) / 2)
    return None if den == 0 else (a['mean'] - b['mean']) / den

def risk(df):
    I = {k: M.clamp(df[k], *BOUNDS[k]) for k in BOUNDS}
    L = M.lp(df.female.to_numpy(float), df.age.to_numpy(float), I['tc'], I['hdl'], I['sbp'], df.dm.to_numpy(float), df.smoking.to_numpy(float), I['egfr'], df.bptreat.to_numpy(float), df.statin.to_numpy(float), model='base')
    return M.risk_from_lp(L)

def main():
    OUT.mkdir(exist_ok=True)
    raw, _ = nh3_build.load()
    E, B, T = nh3_build.build()
    exam, lab, adult = (raw['exam'], raw['lab'], raw['adult'])
    N = E[~E.complete].copy()
    assert len(N) == 220
    for col, f in {'sbp': lambda s: sbp23(exam[s], 'PEP6'), 'tc': lambda s: num(lab[s]['TCP']) if s in lab else None, 'hdl': lambda s: num(lab[s]['HDP']) if s in lab else None, 'scr_raw': lambda s: num(lab[s]['CEP']) if s in lab else None, 'hba1c': lambda s: num(lab[s]['GHP']) if s in lab else None, 'uacr': lambda s: uacr(lab[s]['UBP'], lab[s]['URP']) if s in lab else None, 'dm': lambda s: dm_gate(adult[s], lab[s]['PHPINSU'], lab[s]['GHP']) if s in lab else None}.items():
        N[col] = [f(s) for s in N.SEQN]
    N[['sbp', 'tc', 'hdl', 'scr_raw', 'hba1c', 'uacr', 'dm']] = N[['sbp', 'tc', 'hdl', 'scr_raw', 'hba1c', 'uacr', 'dm']].astype(float)
    N['egfr'] = [egfr(cal(c), a, f == 1) if np.isfinite(c) else np.nan for c, a, f in zip(N.scr_raw, N.age, N.female)]
    S = B.rename(columns={'sbp1': 'sbp', 'tc1': 'tc', 'hdl1': 'hdl', 'egfr1': 'egfr', 'hba1c1': 'hba1c', 'uacr1': 'uacr', 'dm1': 'dm'}).copy()
    Tt = T.rename(columns={'sbp1': 'sbp', 'tc1': 'tc', 'hdl1': 'hdl', 'egfr1': 'egfr', 'uacr1': 'uacr', 'dm1': 'dm'}).copy()
    Tt['dmarethn'] = [code(exam[s]['DMARETHN']) for s in Tt.SEQN]
    Tt['hba1c'] = [num(lab[s]['GHP']) for s in Tt.SEQN]
    Tt['hba1c'] = Tt.hba1c.astype(float)
    S['risk'] = risk(S)
    Tt['risk'] = risk(Tt)
    ok = N[['sbp', 'tc', 'hdl', 'egfr', 'dm', 'smoking', 'bptreat', 'statin']].notna().all(1)
    N['risk'] = np.nan
    N.loc[ok, 'risk'] = risk(N[ok].astype({'statin': float}))
    for d in (S, N, Tt):
        d['nonhdl'] = d.tc - d.hdl
        for k, (lo, hi) in {'a30': (30, 39), 'a40': (40, 59), 'a60': (60, 79)}.items():
            d[k] = d.age.between(lo, hi).astype(float)
        for k, v in {'nhw': 1, 'nhb': 2, 'mex': 3, 'oth': 4}.items():
            d[k] = (d.dmarethn == v).astype(float)
        c = np.where(d.risk.notna(), M.category(d.risk.fillna(0)), -1)
        for j, k in enumerate(['c0', 'c1', 'c2', 'c3']):
            d[k] = np.where(c < 0, np.nan, (c == j).astype(float))
        d['phase2'] = (d.phase == 2).astype(float) if 'phase' in d else np.nan
    rows = [('Age, y', 'age', 'mean'), ('Age 30-39 y', 'a30', 'pct'), ('Age 40-59 y', 'a40', 'pct'), ('Age 60-79 y', 'a60', 'pct'), ('Female', 'female', 'pct'), ('Non-Hispanic White', 'nhw', 'pct'), ('Non-Hispanic Black', 'nhb', 'pct'), ('Mexican American', 'mex', 'pct'), ('Other race or ethnicity', 'oth', 'pct'), ('Current smoking', 'smoking', 'pct'), ('Antihypertensive medication', 'bptreat', 'pct'), ('Statin', 'statin', 'pct'), ('Diabetes (visit-1 rule)', 'dm', 'pct'), ('Systolic BP, mm Hg', 'sbp', 'mean'), ('Total cholesterol, mg/dL', 'tc', 'mean'), ('HDL cholesterol, mg/dL', 'hdl', 'mean'), ('eGFR, mL/min/1.73 m2', 'egfr', 'mean'), ('HbA1c, %', 'hba1c', 'mean'), ('UACR, mg/g', 'uacr', 'median'), ('PREVENT 10-year ASCVD risk, %', 'risk', 'median'), ('Risk <3%', 'c0', 'pct'), ('Risk 3-<5%', 'c1', 'pct'), ('Risk 5-<10%', 'c2', 'pct'), ('Risk >=10%', 'c3', 'pct')]
    table = []
    for label, col, kind in rows:
        k2 = 'pct' if kind == 'pct' else 'mean'
        a = summ(S[col], kind=k2)
        b = summ(N[col], kind=k2)
        t = summ(Tt[col], Tt.w, kind=k2)
        if kind == 'median':
            lg = lambda d: np.log(d[col].astype(float).where(d[col] > 0))
            la, lb, lt = (summ(lg(S)), summ(lg(N)), summ(lg(Tt), Tt.w))
            sa, sb = (smd(la, lt, 'mean'), smd(la, lb, 'mean'))
        else:
            sa, sb = (smd(a, t, k2), smd(a, b, k2))
        table.append({'label': label, 'var': col, 'kind': kind, 'analysis_759': a, 'noncompleters_220': b, 'target_9478_weighted': t, 'smd_analysis_vs_target': sa, 'smd_analysis_vs_noncompleters': sb, 'smd_scale': 'log' if kind == 'median' else 'original'})
    s7 = json.loads((STUDY / 'analysis_v1/stage3/stage3_results.json').read_text(encoding='utf-8'))['S7']
    B['nonhdl1'], B['nonhdl2'] = (B.tc1 - B.hdl1, B.tc2 - B.hdl2)
    inputs = []
    for label, a, b, key in [('Systolic BP, mm Hg', 'sbp1', 'sbp2', 'sbp'), ('Total cholesterol, mg/dL', 'tc1', 'tc2', 'tc'), ('HDL cholesterol, mg/dL', 'hdl1', 'hdl2', 'hdl'), ('Non-HDL cholesterol, mg/dL', 'nonhdl1', 'nonhdl2', 'nonhdl'), ('Serum creatinine (calibrated), mg/dL', 'scr1', 'scr2', 'scr_cal'), ('eGFR, mL/min/1.73 m2', 'egfr1', 'egfr2', 'egfr'), ('HbA1c, %', 'hba1c1', 'hba1c2', 'hba1c_log'), ('Urine albumin, mg/L', 'ualb1', 'ualb2', 'ualb'), ('Urine creatinine, mg/dL', 'ucr1', 'ucr2', 'ucr'), ('UACR, mg/g', 'uacr1', 'uacr2', 'uacr')]:
        sub = B if key not in ('ualb', 'ucr', 'uacr') else B[B.in_uacr]
        if key == 'hba1c_log':
            sub = B[B.hba1c1.notna() & B.hba1c2.notna()]
        x1, x2 = (summ(sub[a]), summ(sub[b]))
        pos = (sub[a] > 0) & (sub[b] > 0)
        gmr = float(np.exp(np.mean(np.log(sub[b][pos].astype(float)) - np.log(sub[a][pos].astype(float))))) if key != 'sbp' else None
        inputs.append({'label': label, 'n': int(len(sub)), 'v1': x1, 'v2': x2, 's7': s7[key], 'geometric_mean_ratio_v2_v1': gmr, 'sd_w_linear': s7['sbp']['sd_w'] if key == 'sbp' else s7['hba1c']['sd_w'] if key == 'hba1c_log' else None, 'icc_linear': s7['hba1c']['icc'] if key == 'hba1c_log' else None})
    dm = {'v1': int(B.dm1.sum()), 'v2': int(B.dm2.sum()), 'discordant': int((B.dm1 != B.dm2).sum())}
    assert dm['discordant'] == s7['diabetes_discordant']
    other = {'interval_days': {'median': float(B.rr_days.median()), 'q1': float(B.rr_days.quantile(0.25)), 'q3': float(B.rr_days.quantile(0.75)), 'min': float(B.rr_days.min()), 'max': float(B.rr_days.max())}, 'phase': {'1': int((B.phase == 1).sum()), '2': int((B.phase == 2).sum())}, 'same_arm': int(B.same_arm.sum()), 'same_examiner': int(B.same_examiner.sum()), 'n_clusters': int(B.cluster.nunique()), 'n_strata': int(B.stratum.nunique()), 'noncompleter_risk_n': int(ok.sum()), 'target_groups': T.group.value_counts().to_dict(), 'target_w_sum': float(T.w.sum())}
    res = {'note': 'Descriptive; unweighted for the 759 and 220; target weighted by WTPFEX6. Standardized differences use pooled SDs (binary: sqrt of mean Bernoulli variance). Visit-1 diabetes: protocol rule for the 759 (identical to the gate visit-1 rule for all 759, P3), gate visit-1 rule for others. Risk uses the primary clamp; non-completer risk only where all visit-1 inputs are present.', 'table1': table, 'inputs_by_visit': inputs, 'diabetes': dm, 'other': other}
    (OUT / 'table1.json').write_text(json.dumps(res, indent=1, default=float) + '\n', encoding='utf-8')
    for r in table:
        f = lambda s: (f"{s['pct']:.1f}" if 'pct' in s else f"{s['mean']:.1f} ({s['sd']:.1f}) med {s['median']:.2f}") + f" n={s['n']}"
        print(f"{r['label']:<36} {f(r['analysis_759']):<30} {f(r['noncompleters_220']):<30} {f(r['target_9478_weighted']):<30} smdT={(r['smd_analysis_vs_target'] if r['smd_analysis_vs_target'] is None else round(r['smd_analysis_vs_target'], 3))} smdN={(r['smd_analysis_vs_noncompleters'] if r['smd_analysis_vs_noncompleters'] is None else round(r['smd_analysis_vs_noncompleters'], 3))}")
    print(json.dumps(other, default=float))
    print(dm)
if __name__ == '__main__':
    main()
