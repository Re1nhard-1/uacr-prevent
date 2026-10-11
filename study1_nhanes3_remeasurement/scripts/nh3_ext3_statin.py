from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
import nh3_build
import nh3_models as M
import nh3_stats as S
import nh3_testretest_stage12 as s12
from nh3_stage3 import fixed, lp_of, visit_inputs
from nh3_ext2_analyses import wild_ci
ROOT = Path(__file__).resolve().parents[1]
ST = ROOT / 'studies/nhanes3_testretest'
PROTO = ST / 'extensions3_protocol_20261010.md'
OUT = ST / 'analysis_v1_posthoc/ext3_statin_results.json'
LEVEL = {'N': 0, 'L': 0, 'R': 1, 'C': 1, 'S': 2, 'I': 2, 'H': 3}
REC = ['I', 'S', 'H']
REC3 = {'v1': ['I', 'C', 'S', 'H'], 'v2': ['I', 'R', 'C', 'S', 'H']}
LET = {'v1': 'NLCSHI', 'v2': 'NLRCSHI'}
T30 = M._DAT['base_30yr']
assert len(T30) == 24 and T30.beta_coefficients.iloc[1] == 'Age squared' and (T30.beta_coefficients.iloc[-1] == 'Constant')

def sampson(tc, hdl, tg):
    nonhdl = tc - hdl
    return tc / 0.948 - hdl / 0.971 - (tg / 8.56 + tg * nonhdl / 2140 - tg ** 2 / 16100) - 9.44

def tiers(ldl, risk, dm, egfr, age, rule='v1', r30=None):
    a4075 = (age >= 40) & (age <= 75)
    if rule == 'v1':
        ind = (ldl >= 190) | a4075 & ((dm == 1) | (egfr >= 15) & (egfr < 60))
        low2 = np.zeros(len(ldl), bool)
    else:
        a3059 = (age >= 30) & (age <= 59)
        assert np.isfinite(r30[a3059]).all(), '30-year risk missing at age 30-59'
        ind = (ldl >= 190) | a4075 & (dm == 1) | a4075 & (egfr < 60) & (ldl >= 70)
        low2 = a3059 & ((ldl >= 160) | (np.where(a3059, r30, 0.0) >= 10))
    cat = M.category(risk)
    out = np.empty(len(ldl), dtype=object)
    for i in range(len(ldl)):
        if ind[i]:
            out[i] = 'I'
        elif ldl[i] < 70:
            out[i] = 'N'
        elif cat[i] == 0 and low2[i]:
            out[i] = 'R'
        else:
            out[i] = 'LCSH'[cat[i]]
    return out

def risk30(F, I):
    X = M.features(F['female'], F['age'], I['tc'], I['hdl'], I['sbp'], I['dm'], F['smoking'], I['egfr'], F['bptreat'], F['statin'])
    a = (np.asarray(F['age'], float) - 55) / 10
    X = np.column_stack([X[:, :1], a ** 2, X[:, 1:]])
    f = np.asarray(F['female'], int) == 1
    return M.risk_from_lp(np.where(f, X @ np.asarray(T30['female_ascvd'], float), X @ np.asarray(T30['male_ascvd'], float)))

def check30(F, I, r, idx):
    worst = 0.0
    for i in idx:
        s = M.prevent_base(sex=float(F['female'][i]), age=float(F['age'][i]), tc=float(I['tc'][i]), hdl=float(I['hdl'][i]), sbp=float(I['sbp'][i]), dm=float(I['dm'][i]), smoking=float(F['smoking'][i]), bmi=float('nan'), egfr=float(I['egfr'][i]), bptreat=float(F['bptreat'][i]), statin=float(F['statin'][i]))['prevent_base_30yr_ASCVD']
        worst = max(worst, abs(s - r[i])) if np.isfinite(s) and np.isfinite(r[i]) else np.inf
    return worst

def prop(y, cl):
    r = S.prop_ci(np.asarray(y, float), cl)
    return {k: r[k] for k in ['n', 'events', 'estimate', 'lower', 'upper', 'deff', 'method']}

def main():
    E, B, T = nh3_build.build()
    F = fixed(B)
    cl = B.cluster.to_numpy()
    age = B.age.to_numpy(float)
    keep = set(B.SEQN)
    l1 = s12.read('lab', {'TGP': (1606, 1609), 'LCP': (1615, 1617), 'PHPFAST': (1263, 1267)}, keep=keep)
    l2 = s12.read('labse', {'TGR': (358, 361), 'LCR': (367, 369), 'PHRFAST': (28, 32)}, keep=keep)
    g = lambda d, s, k: s12.num(d[s][k]) if s in d else None
    tg1 = np.array([g(l1, s, 'TGP') for s in B.SEQN], dtype=float)
    tg2 = np.array([g(l2, s, 'TGR') for s in B.SEQN], dtype=float)
    lc1 = np.array([g(l1, s, 'LCP') for s in B.SEQN], dtype=float)
    lc2 = np.array([g(l2, s, 'LCR') for s in B.SEQN], dtype=float)
    tc1, tc2, hd1, hd2 = (B[c].to_numpy(float) for c in ['tc1', 'tc2', 'hdl1', 'hdl2'])
    I1, I2 = (visit_inputs(B, 1), visit_inputs(B, 2))
    r1, r2 = (M.risk_from_lp(lp_of(F, I1, 'base')), M.risk_from_lp(lp_of(F, I2, 'base')))
    Pp = pd.read_csv(ST / 'analysis_v1/stage3/stage3_participants.csv').set_index('SEQN').loc[B.SEQN]
    assert np.max(np.abs(r1 - Pp.riskB_m0.to_numpy())) < 1e-09 and np.max(np.abs(r2 - Pp.riskB_m15.to_numpy())) < 1e-09, 'visit risks not reproduced'
    e1, e2 = (B.egfr1.to_numpy(float), B.egfr2.to_numpy(float))
    statin = B.statin.to_numpy() == 1
    tg_ok = np.isfinite(tg1) & np.isfinite(tg2) & (tg1 <= 800) & (tg2 <= 800)
    base_pop = ~statin & tg_ok
    tg_miss = ~(np.isfinite(tg1) & np.isfinite(tg2))
    flow = {'analysis_set': int(len(B)), 'statin_users_excluded': int(statin.sum()), 'tg_missing_excluded': int((~statin & tg_miss).sum()), 'tg_over_800_excluded': int((~statin & ~tg_miss & ~tg_ok).sum()), 'analysed': int(base_pop.sum()), 'note': 'TG >800 exclusion is in the protocol; the missing-TG exclusion is a protocol departure (addendum item 1)'}
    exm = s12.read('exam', {'PEPTECH': (1498, 1502)}, keep=keep)
    phys = np.array([exm[s]['PEPTECH'].strip() for s in B.SEQN])
    phys_single = np.where(phys == '88888', np.array([f'missing_{s}' for s in B.SEQN]), phys)
    a3059 = (age >= 30) & (age <= 59)
    r30_1, r30_2 = (risk30(F, I1), risk30(F, I2))
    rec = lambda t: np.isin(t, REC)

    def run(ldl1, ldl2, dm1, dm2, pop, label, extra=True, rule='v1', R=(None, None), R30=(None, None)):
        Ra, Rb = (r1, r2) if R[0] is None else R
        idx = np.where(pop & np.isfinite(ldl1) & np.isfinite(ldl2))[0]
        t1, t2 = (tiers(ldl1, Ra, dm1, e1, age, rule, R30[0]), tiers(ldl2, Rb, dm2, e2, age, rule, R30[1]))
        rec3 = lambda t: np.isin(t, REC3[rule])
        lv1, lv2 = (np.array([LEVEL[x] for x in t1]), np.array([LEVEL[x] for x in t2]))
        ch_rec = rec(t1) != rec(t2)
        res = {'label': label, 'rule': rule, 'n': int(len(idx)), 'primary_recommended_change': prop(ch_rec[idx], cl[idx]), 'recommended_v1': prop(rec(t1)[idx], cl[idx]), 'recommended_v2': prop(rec(t2)[idx], cl[idx]), 'recommended_or_considered_change' if rule == 'v1' else 'recommended_or_reasonable_change': prop((rec3(t1) != rec3(t2))[idx], cl[idx]), 'action_level_change': prop((lv1 != lv2)[idx], cl[idx]), 'high_intensity_change': prop(((t1 == 'H') != (t2 == 'H'))[idx], cl[idx]), 'tier_counts_v1': {k: int((t1[idx] == k).sum()) for k in LET[rule]}, 'tier_counts_v2': {k: int((t2[idx] == k).sum()) for k in LET[rule]}}
        if extra:
            catch = M.category(Ra) != M.category(Rb)
            ldlx = np.digitize(ldl1, [70, 190]) != np.digitize(ldl2, [70, 190])
            dmx = dm1 != dm2
            egx = (e1 < 60) != (e2 < 60)
            path = list('LCSH') + (['R'] if rule == 'v2' else [])

            def mech(sel, names, flags):
                pat = pd.Series([' + '.join((n for n, f in zip(names, fl) if f)) or 'none' for fl in zip(*[f[idx][sel] for f in flags])]).value_counts()
                alone = catch.copy()
                for f in flags[1:]:
                    alone &= ~f
                return {'changes': int(sel.sum()), **{f'{n}_flag': int(f[idx][sel].sum()) for n, f in zip(names, flags)}, 'prevent_category_alone': int(alone[idx][sel].sum()), 'overlap_pattern': {k: int(v) for k, v in pat.items()}, 'prevent_category_changed_risk_pathway_both_visits': int((catch & np.isin(t1, path) & np.isin(t2, path))[idx][sel].sum())}
            m = mech(ch_rec[idx], ['PREVENT', 'LDL', 'diabetes', 'eGFR'], [catch, ldlx, dmx, egx])
            res['mechanisms_among_changes'] = {'changes': m['changes'], 'prevent_category_changed': m['PREVENT_flag'], 'ldl_crossed_70_or_190': m['LDL_flag'], 'diabetes_flag_changed': m['diabetes_flag'], 'egfr60_flag_changed': m['eGFR_flag'], 'prevent_category_alone': m['prevent_category_alone'], 'overlap_pattern': m['overlap_pattern'], 'prevent_category_changed_risk_pathway_both_visits': m['prevent_category_changed_risk_pathway_both_visits']}
            if rule == 'v2':
                lowx = a3059 & ((M.category(Ra) == 0) | (M.category(Rb) == 0))
                l160 = lowx & ((ldl1 >= 160) != (ldl2 >= 160))
                r30x = lowx & ((R30[0] >= 10) != (R30[1] >= 10))
                res['mechanisms_among_reasonable_changes'] = mech((rec3(t1) != rec3(t2))[idx], ['PREVENT', 'LDL70_190', 'diabetes', 'eGFR', 'LDL160', '30-year'], [catch, ldlx, dmx, egx, l160, r30x])
            tr = pd.crosstab(pd.Series(t1[idx], name='v1'), pd.Series(t2[idx], name='v2')).reindex(index=list(LET[rule]), columns=list(LET[rule]), fill_value=0)
            res[f'transitions_{LET[rule]}'] = tr.values.tolist()
            r_ = S.prop_ci(ch_rec[idx].astype(float), phys[idx])
            res['primary_visit1_physician'] = {k: r_[k] for k in ['lower', 'upper', 'deff', 'G']}
            res['primary_visit1_physician_wild'] = wild_ci(ch_rec[idx].astype(float), phys[idx])
            res['physician_note'] = 'code 88888 (missing physician) forms one cluster, as in E3/E5; singleton version below'
            r_ = S.prop_ci(ch_rec[idx].astype(float), phys_single[idx])
            res['primary_visit1_physician_missing_singletons'] = {k: r_[k] for k in ['lower', 'upper', 'deff', 'G']}
        return (res, idx, t1, t2)
    ptxt = PROTO.read_bytes()
    out = {'label': 'E9: statin-recommendation tier under the 2026 ACC/AHA dyslipidemia guideline at two visits (NHANES III Study 1)', 'protocol': str(PROTO.relative_to(ROOT)), 'protocol_sha256': hashlib.sha256(ptxt).hexdigest(), 'protocol_sha256_before_addendum': hashlib.sha256(ptxt.split(b'\n## Addendum')[0]).hexdigest(), 'protocol_sha256_before_second_addendum': hashlib.sha256(ptxt.split(b'\n## Second addendum')[0]).hexdigest(), 'flow': flow}
    dm1, dm2 = (I1['dm'], I2['dm'])
    ls1, ls2 = (sampson(tc1, hd1, tg1), sampson(tc2, hd2, tg2))
    geo = lambda a, b: np.sqrt(a * b)
    h1, h2 = (B.hba1c1.to_numpy(float), B.hba1c2.to_numpy(float))
    hbar = np.where(np.isfinite(h1) & np.isfinite(h2), (h1 + h2) / 2, np.where(np.isfinite(h1), h1, h2))
    dmA = np.where((B.fixed_dm.to_numpy() == 1) | (np.nan_to_num(hbar) >= 6.5), 1.0, 0.0)
    egA = np.asarray(M.ckd_epi_2021(geo(B.scr1.to_numpy(float), B.scr2.to_numpy(float)), B.age, B.female), float)
    tcA, hdA, tgA = (geo(tc1, tc2), geo(hd1, hd2), geo(tg1, tg2))
    IA = {'sbp': M.clamp((B.sbp1.to_numpy(float) + B.sbp2.to_numpy(float)) / 2, 90, 200), 'tc': M.clamp(tcA, 130, 320), 'hdl': M.clamp(hdA, 20, 100), 'egfr': M.clamp(egA, 15, 140), 'dm': dmA, 'uacr': None}
    rA, r30A = (M.risk_from_lp(lp_of(F, IA, 'base')), risk30(F, IA))
    g1, g2 = (np.where(np.isfinite(h1), h1, h2), np.where(np.isfinite(h2), h2, h1))
    conf = np.where((B.fixed_dm.to_numpy() == 1) | (np.nan_to_num(g1) >= 6.5) & (np.nan_to_num(g2) >= 6.5), 1.0, 0.0)
    I1c, I2c = ({**I1, 'dm': conf}, {**I2, 'dm': conf})
    rc1, rc2 = (M.risk_from_lp(lp_of(F, I1c, 'base')), M.risk_from_lp(lp_of(F, I2c, 'base')))
    fh1 = np.array([g(l1, s, 'PHPFAST') for s in B.SEQN], dtype=float)
    fh2 = np.array([g(l2, s, 'PHRFAST') for s in B.SEQN], dtype=float)
    fin = np.isfinite(fh1) & np.isfinite(fh2)
    same = fin & ((fh1 >= 9) == (fh2 >= 9))

    def analyse(rule):
        R30 = (r30_1, r30_2) if rule == 'v2' else (None, None)
        prim, idx, t1, t2 = run(ls1, ls2, dm1, dm2, base_pop, 'Sampson LDL-C (primary)', rule=rule, R30=R30)
        tA = tiers(sampson(tcA, hdA, tgA), rA, dmA, egA, age, rule, r30A if rule == 'v2' else None)
        d1, d2 = (rec(t1) != rec(tA), rec(t2) != rec(tA))
        prim['vs_two_visit_mean'] = {'v1_differs': prop(d1[idx], cl[idx]), 'v2_differs': prop(d2[idx], cl[idx]), 'mean_single_visit_differs': S.paired_ci((d1 + 0.0 + d2)[idx] / 2, cl[idx]), 'recommended_two_visit': prop(rec(tA)[idx], cl[idx])}
        sens = {}
        fw_ok = np.isfinite(lc1) & np.isfinite(lc2)
        sens['sens_nhanes_friedewald_reported'] = run(lc1, lc2, dm1, dm2, ~statin & fw_ok, 'NHANES-reported Friedewald LDL-C at both visits', extra=False, rule=rule, R30=R30)[0]
        fwc1, fwc2 = (tc1 - hd1 - tg1 / 5, tc2 - hd2 - tg2 / 5)
        sens['sens_friedewald_computed'] = run(fwc1, fwc2, dm1, dm2, ~statin & np.isfinite(tg1) & np.isfinite(tg2) & (tg1 < 400) & (tg2 < 400), 'Friedewald computed, TG <400 at both visits', extra=False, rule=rule, R30=R30)[0]
        R30c = (risk30(F, I1c), risk30(F, I2c)) if rule == 'v2' else (None, None)
        t1c, t2c = (tiers(ls1, rc1, conf, e1, age, rule, R30c[0]), tiers(ls2, rc2, conf, e2, age, rule, R30c[1]))
        sens['sens_confirmed_dm_fixed'] = {'n': int(len(idx)), 'primary_recommended_change': prop((rec(t1c) != rec(t2c))[idx], cl[idx])}
        chp = rec(t1) != rec(t2)
        js = idx[same[idx]]
        fast = {'fasting_hours_missing_either_visit': int((~fin)[idx].sum()), 'fasting_status_differs': int((fin & ~same)[idx].sum()), 'changes_with_fasting_status_differing': int((chp & fin & ~same)[idx].sum()), 'primary_change_same_fasting_status': prop(chp[js], cl[js])}
        return (prim, sens, fast, idx, t1, t2)
    prim, sens, fast, idx, t1, t2 = analyse('v1')
    out['primary'] = prim
    out.update(sens)
    out['audit_fasting_check'] = fast
    out['ldl_sampson_summary'] = {'median_v1': float(np.nanmedian(ls1[idx])), 'median_v2': float(np.nanmedian(ls2[idx])), 'cvw_pct': float(100 * np.nanstd(np.log(np.maximum(ls2[idx], 1) / np.maximum(ls1[idx], 1)), ddof=1) / np.sqrt(2)), 'in_70_189_v1': int(((ls1[idx] >= 70) & (ls1[idx] < 190)).sum())}
    j = idx[a3059[idx]]
    chk = max(check30(F, I1, r30_1, idx), check30(F, I2, r30_2, idx), check30(F, IA, r30A, idx))
    assert chk < 1e-09, ('30-year risk mismatch', chk)
    vp, vs, vf, _, u1, u2 = analyse('v2')
    vp['risk30'] = {'table_vs_pyprevent_max_abs_diff_pp': chk, 'aged_30_59': int(len(j)), 'median_v1': float(np.median(r30_1[j])), 'ge10_v1': int((r30_1[j] >= 10).sum()), 'ge10_v2': int((r30_2[j] >= 10).sum()), 'R_tier_v1_by_reason': {'ldl_ge160': int(((u1 == 'R') & (ls1 >= 160))[idx].sum()), 'risk30_ge10_only': int(((u1 == 'R') & (ls1 < 160))[idx].sum())}}
    vp['tier_change_vs_rule_v1'] = {'visit1_tier_differs': int((u1 != t1)[idx].sum()), 'visit2_tier_differs': int((u2 != t2)[idx].sum()), 'visit1_recommended_differs': int((rec(u1) != rec(t1))[idx].sum()), 'visit2_recommended_differs': int((rec(u2) != rec(t2))[idx].sum())}
    out['verified_rule'] = {'note': 'rule v2 from the protocol second addendum (primary guideline text); reference E9 result', 'primary': vp, **vs, 'audit_fasting_check': vf}
    OUT.write_text(json.dumps(out, indent=1, default=lambda o: o.item() if hasattr(o, 'item') else o.tolist() if hasattr(o, 'tolist') else str(o)) + '\n', encoding='utf-8')
    print('saved', OUT)
if __name__ == '__main__':
    main()
