"""E1 (post hoc, 2026-10-10; extensions_protocol_20261010.md): laboratory-only remeasurement of PREVENT inputs in the NHANES 2001-2002"""
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
ROOT = Path(__file__).resolve().parents[1]
ST = ROOT / 'studies/nhanes3_testretest'
G01 = ROOT / 'studies/uacr_prevent/external_2001_gate/data'
R01 = ST / 'data_2001_repeat/raw'
OUT = ST / 'analysis_v1_posthoc/ext2001_labs.json'

def xpt(p, cols):
    df, _ = pyreadstat.read_xport(str(p), usecols=['SEQN'] + cols, encoding='latin1')
    return df

def svy_prop(y, st, ps):
    y = np.asarray(y, float)
    n = len(y)
    x = float(y.sum())
    th = x / n
    _, var, df = S.strat_taylor(y, st, ps)
    out = {'n': n, 'events': int(round(x)), 'estimate': th, 'df': int(df)}
    base = th * (1 - th) / n if 0 < th < 1 else 0.0
    if x < 10 or x > n - 10:
        lo, hi = S.kg_interval(x, n, var / base if base > 0 else 1.0, df)
        out.update(lower=lo, upper=hi, method='Korn-Graubard')
    else:
        se = math.sqrt(var) / (th * (1 - th))
        q = tdist.ppf(0.975, df)
        out.update(lower=float(expit(logit(th) - q * se)), upper=float(expit(logit(th) + q * se)), method='stratified Taylor, logit t')
    return out

def lp(d, sfx, model):
    return M.lp(d.female.to_numpy(float), d.age.to_numpy(float), M.clamp(d['tc' + sfx], 130, 320), M.clamp(d['hdl' + sfx], 20, 100), M.clamp(d.sbp, 90, 200), d['dm' + sfx].to_numpy(float), d.smoking.to_numpy(float), M.clamp(d['egfr' + sfx], 15, 140), d.bptreat.to_numpy(float), d.statin.to_numpy(float), uacr=M.clamp(d['uacr' + sfx], 0.1, 25000) if model == 'uacr' else None, model=model)

def main():
    X = pd.read_csv(ROOT / 'studies/uacr_prevent/final_and_external_v1/external_bridged_participant_results.csv')
    assert len(X) == 278
    demo = xpt(G01 / 'DEMO_B.xpt', ['RIAGENDR'])
    v1 = xpt(G01 / 'L13_B.xpt', ['LBXTC', 'LBDHDL']).merge(xpt(G01 / 'L40_B.xpt', ['LBDSCR']), on='SEQN', how='outer').merge(xpt(G01 / 'L10_B.xpt', ['LBXGH']), on='SEQN', how='outer')
    diq = xpt(G01 / 'DIQ_B.xpt', ['DIQ010', 'DIQ050', 'DIQ070'])
    v2 = xpt(R01 / 'L13_2_B.xpt', ['LB2DAY', 'LB2TC', 'LB2HDL']).rename(columns={'LB2DAY': 'lab_days'}).merge(xpt(R01 / 'L40_2_B.xpt', ['LB2SCR']), on='SEQN', how='outer').merge(xpt(R01 / 'L10_2_B.xpt', ['LB2GH']), on='SEQN', how='outer')
    d = X.merge(demo, on='SEQN', how='left').merge(v1, on='SEQN', how='left').merge(diq, on='SEQN', how='left').merge(v2, on='SEQN', how='left')
    d['female'] = (d.RIAGENDR == 2).astype(float)
    assert np.array_equal(d.female.to_numpy(), d.sex.to_numpy(float)), 'sex coding'
    assert np.allclose(d.tc, d.LBXTC) and np.allclose(d.hdl, d.LBDHDL)
    eg1 = M.ckd_epi_2021(d.LBDSCR.to_numpy(float), d.age.to_numpy(float), d.female.to_numpy(float))
    assert np.allclose(eg1, d.egfr, atol=1e-08)
    qpos = ((d.DIQ010 == 1) | (d.DIQ050 == 1) | (d.DIQ070 == 1)).to_numpy()
    dm1 = (qpos | (np.nan_to_num(d.LBXGH.to_numpy(float), nan=0) >= 6.5)).astype(float)
    assert np.array_equal(dm1, d.dm.to_numpy(float)), 'diabetes rule'
    d['tc1'], d['hdl1'], d['egfr1'], d['dm1'], d['uacr1'] = (d.tc, d.hdl, d.egfr, dm1, d.uacr_first_bridged)
    r1u = M.risk_from_lp(lp(d, '1', 'uacr'))
    assert np.max(np.abs(r1u - d.risk_first)) < 1e-09, np.max(np.abs(r1u - d.risk_first))
    flow = {'study2_external': len(d)}
    ok = d[['LB2TC', 'LB2HDL', 'LB2SCR', 'LB2GH']].notna().all(axis=1)
    flow['all_second_exam_labs'] = int(ok.sum())
    e = d[ok].copy()
    e['tc2'], e['hdl2'] = (e.LB2TC, e.LB2HDL)
    e['egfr2'] = M.ckd_epi_2021(e.LB2SCR.to_numpy(float), e.age.to_numpy(float), e.female.to_numpy(float))
    e['dm2'] = (qpos[ok.to_numpy()] | (e.LB2GH.to_numpy(float) >= 6.5)).astype(float)
    e['uacr2'] = e.uacr_repeat_bridged
    st, ps = (e.SDMVSTRA.to_numpy(), (e.SDMVSTRA * 10 + e.SDMVPSU).to_numpy())
    cat = lambda v: M.category(M.risk_from_lp(v))
    res = {'label': __doc__.split('\n')[0], 'flow': flow, 'interval_days': {k: float(v) for k, v in e.lab_days.describe(percentiles=[0.25, 0.5, 0.75]).items() if k in ('min', '25%', '50%', '75%', 'max')}}
    for model, mix in [('base', None), ('uacr', None)]:
        a, b = (lp(e, '1', model), lp(e, '2', model))
        ca, cb = (cat(a), cat(b))
        ch = (ca != cb).astype(float)
        res[model] = {'theta_lab': svy_prop(ch, st, ps), 'up': int((cb > ca).sum()), 'down': int((cb < ca).sum()), 'mean_abs_delta_pp': float(np.mean(np.abs(M.risk_from_lp(b) - M.risk_from_lp(a)))), 'transitions_visit1_rows': pd.crosstab(ca, cb).to_dict(orient='index')}
    clamp2 = ((e.tc2 < 130) | (e.tc2 > 320) | (e.hdl2 < 20) | (e.hdl2 > 100) | (e.egfr2 < 15) | (e.egfr2 > 140)).to_numpy()
    clamp1 = ((e.tc1 < 130) | (e.tc1 > 320) | (e.hdl1 < 20) | (e.hdl1 > 100) | (e.egfr1 < 15) | (e.egfr1 > 140) | (e.sbp < 90) | (e.sbp > 200)).to_numpy()

    def raw_lp(d_, sfx):
        return M.lp(d_.female.to_numpy(float), d_.age.to_numpy(float), d_['tc' + sfx].to_numpy(float), d_['hdl' + sfx].to_numpy(float), d_.sbp.to_numpy(float), d_['dm' + sfx].to_numpy(float), d_.smoking.to_numpy(float), d_['egfr' + sfx].to_numpy(float), d_.bptreat.to_numpy(float), d_.statin.to_numpy(float), model='base')
    ch_b = (cat(lp(e, '1', 'base')) != cat(lp(e, '2', 'base'))).astype(float)
    keep = ~(clamp1 | clamp2)
    res['clamping'] = {'visit1_clamped_n': int(clamp1.sum()), 'visit2_clamped_n': int(clamp2.sum()), 'visit2_tc_below_130': int((e.tc2 < 130).sum()), 'visit2_hdl_above_100': int((e.hdl2 > 100).sum()), 'base_excluding_clamped': {'n': int(keep.sum()), 'events': int(ch_b[keep].sum()), 'estimate': float(ch_b[keep].mean())}, 'base_without_clamping': {'n': len(e), 'events': int((cat(raw_lp(e, '1')) != cat(raw_lp(e, '2'))).sum()), 'estimate': float((cat(raw_lp(e, '1')) != cat(raw_lp(e, '2'))).mean())}}
    lg = lambda a, b: float(100 * np.std(np.log(b / a), ddof=1) / math.sqrt(2))
    res['cvw_pct'] = {'tc': lg(e.tc1, e.tc2), 'hdl': lg(e.hdl1, e.hdl2), 'creatinine': lg(e.LBDSCR, e.LB2SCR)}
    res['hba1c_sdw'] = float(np.std(e.LB2GH - e.LBXGH, ddof=1) / math.sqrt(2))
    res['diabetes_discordant'] = int((e.dm1 != e.dm2).sum())
    P = pd.read_csv(ST / 'analysis_v1/stage3/stage3_participants.csv')
    c3 = lambda r: np.searchsorted([3.0, 5.0, 10.0], r, side='right')
    lab1 = (c3(P.riskB_m14) != c3(P.riskB_m0)).astype(float)
    U = P[P.in_uacr.astype(bool)]
    labu = (c3(U.riskU_m30) != c3(U.riskU_m0)).astype(float)
    res['study1_comparator'] = {'base_lab_only_ref1': S.prop_ci(lab1, P.cluster.to_numpy()), 'uacr_lab_urine_ref1': S.prop_ci(labu, U.cluster.to_numpy())}
    S7 = json.loads((ST / 'analysis_v1/stage3/stage3_results.json').read_text(encoding='utf-8'))['S7']
    res['study1_cvw_pct'] = {'tc': S7['tc']['cvw_total_pct'], 'hdl': S7['hdl']['cvw_total_pct'], 'creatinine_calibrated': S7['scr_cal']['cvw_total_pct']}
    OUT.write_text(json.dumps(res, indent=1, default=lambda o: o.item() if hasattr(o, 'item') else str(o)) + '\n', encoding='utf-8')
    for k in ['flow', 'interval_days', 'cvw_pct', 'study1_cvw_pct', 'hba1c_sdw', 'diabetes_discordant']:
        print(k, res[k])
    for m in ['base', 'uacr']:
        print(m, {kk: round(100 * vv, 2) if isinstance(vv, float) else vv for kk, vv in res[m]['theta_lab'].items()}, 'up', res[m]['up'], 'down', res[m]['down'])
    print('study1', {k: (round(100 * v['estimate'], 2), round(100 * v['lower'], 2), round(100 * v['upper'], 2)) for k, v in res['study1_comparator'].items()})
if __name__ == '__main__':
    main()
