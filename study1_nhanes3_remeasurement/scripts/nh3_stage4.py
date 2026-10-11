from datetime import datetime, timezone
from pathlib import Path
import json
import math
import shutil
import numpy as np
import pandas as pd
import nh3_build
import nh3_models as M
import nh3_stats as S
from nh3_stage3 import visit_inputs, fixed, lp_of, hybrid, scalar_check, shapley, PRIMARY, GROUPS, GKEYS, SEED
ROOT = Path(__file__).resolve().parents[1]
A1 = ROOT / 'studies/nhanes3_testretest/analysis_v1'
OUT = A1 / 'stage4'

def save(name, x):
    (OUT / name).write_text(json.dumps(x, indent=2, default=lambda o: o.item() if hasattr(o, 'item') else str(o)) + '\n', encoding='utf-8')

def change(F, I1, I2, model='base', outcome='ascvd', rounded=False):
    oi = 0 if outcome == 'ascvd' else 1
    l1, l2 = (lp_of(F, I1, model, outcome), lp_of(F, I2, model, outcome))
    w = max(scalar_check(F, I1, model, oi, l1, range(len(l1))), scalar_check(F, I2, model, oi, l2, range(len(l2))))
    assert w < (0.0001 if outcome == 'total_cvd' else 1e-10), w
    r1, r2 = (M.risk_from_lp(l1), M.risk_from_lp(l2))
    if rounded:
        r1, r2 = (np.round(r1, 1), np.round(r2, 1))
    if outcome == 'total_cvd':
        return ((r1 >= 7.5) != (r2 >= 7.5)).astype(float)
    return (M.category(r1) != M.category(r2)).astype(float)

def sub(d, idx):
    return {k: None if v is None else v[idx] for k, v in d.items()}

def main():
    if OUT.exists():
        raise FileExistsError('Keep existing execution')
    st3 = json.loads((A1 / 'stage3/stage3_results.json').read_text(encoding='utf-8'))
    OUT.mkdir()
    shutil.copy2(__file__, OUT / 'executed_nh3_stage4.py')
    E, B, T = nh3_build.build()
    n = len(B)
    cl = B.cluster.to_numpy()
    F = fixed(B)
    U = np.where(B.in_uacr.to_numpy())[0]
    I1, I2 = (visit_inputs(B, 1), visit_inputs(B, 2))
    theta0 = st3['primary']['estimate']
    thetaU0 = st3['K1']['theta_U']['estimate']
    res = {}

    def run(name, I1x, I2x, idx=None, Fx=None, also_u=False, rounded=False, key=False):
        idx = np.arange(n) if idx is None else np.asarray(idx)
        Fx = F if Fx is None else Fx
        ch = change(sub(Fx, idx), sub(I1x, idx), sub(I2x, idx), rounded=rounded)
        out = {'theta': S.prop_ci(ch, cl[idx]), 'delta_vs_primary_pp': 100 * (ch.mean() - theta0), 'key': key}
        if key:
            out['abstract_flag_ge_1.5pp'] = abs(ch.mean() - theta0) >= 0.015
        if also_u:
            iu = np.intersect1d(idx, U)
            chu = change(sub(Fx, iu), sub(I1x, iu), sub(I2x, iu), model='uacr')
            out['theta_U'] = S.prop_ci(chu, cl[iu])
            out['delta_U_vs_primary_pp'] = 100 * (chu.mean() - thetaU0)
            cv = change(sub(Fx, idx), sub(I1x, idx), sub(I2x, idx), outcome='total_cvd')
            out['S10_cross75'] = S.prop_ci(cv, cl[idx])
        res[name] = out
        return ch
    raw = lambda v: {k: B[f'{k}{v}'].to_numpy(float) for k in ['sbp', 'tc', 'hdl', 'egfr']}

    def outdom(v, b):
        return np.any([~((raw(v)[k] >= b[k][0]) & (raw(v)[k] <= b[k][1])) for k in b], 0)
    B180 = {**PRIMARY, 'sbp': (90.0, 180.0)}
    run('SA1a_exclude_out_of_domain', I1, I2, np.where(~outdom(1, PRIMARY) & ~outdom(2, PRIMARY))[0], also_u=True, key=True)
    run('SA1b_clamp_sbp_90_180', visit_inputs(B, 1, B180), visit_inputs(B, 2, B180), also_u=True)
    run('SA1c_in_domain_v1', I1, I2, np.where(~outdom(1, PRIMARY))[0], also_u=True)
    run('SA1d_paper1_exclusion', visit_inputs(B, 1, B180), visit_inputs(B, 2, B180), np.where(~outdom(1, B180) & ~outdom(2, B180))[0], also_u=True)
    s2 = lambda a, b: ({**I1, 'sbp': M.clamp(B[a], 90, 200)}, {**I2, 'sbp': M.clamp(B[b], 90, 200)})
    run('SA2a_single_reading2', *s2('sbp1_r2', 'sbp2_r2'), key=True)
    run('SA2b_mean_available_1_3', *s2('sbp1_avail', 'sbp2_avail'))
    run('SA2c_exclude_different_arm', I1, I2, np.where(B.same_arm.to_numpy())[0])
    h1, h2 = (B.hba1c1.to_numpy(float), B.hba1c2.to_numpy(float))
    g1, g2 = (np.where(np.isfinite(h1), h1, h2), np.where(np.isfinite(h2), h2, h1))
    conf = np.where((B.fixed_dm.to_numpy() == 1) | (np.nan_to_num(g1) >= 6.5) & (np.nan_to_num(g2) >= 6.5), 1.0, 0.0)
    run('SA3a_confirmed_diabetes_fixed', {**I1, 'dm': conf}, {**I2, 'dm': conf}, also_u=True, key=True)
    res['SA3a_confirmed_diabetes_fixed']['interpretive_e_abstract_both'] = abs(res['SA3a_confirmed_diabetes_fixed']['theta']['estimate'] - theta0) >= 0.015
    run('SA3b_exclude_hba1c_method_switch', I1, I2, np.where(~B.hba1c_method_switch.to_numpy(bool))[0], also_u=True)
    run('SA3c_gate_rule', {**I1, 'dm': B.dm1_gate.to_numpy(float)}, {**I2, 'dm': B.dm2_gate.to_numpy(float)}, also_u=True)
    res['SA3c_gate_rule']['discordant'] = int((B.dm1_gate != B.dm2_gate).sum())
    eg = lambda scr: M.clamp(M.ckd_epi_2021(scr, B.age, B.female), 15, 140)
    run('SA4a_raw_creatinine', {**I1, 'egfr': eg(B.scr1_raw)}, {**I2, 'egfr': eg(B.scr2_raw)}, also_u=True, key=True)
    run('SA4b_coresh_minus_0.23', {**I1, 'egfr': eg(np.maximum(B.scr1_raw - 0.23, 1e-06))}, {**I2, 'egfr': eg(np.maximum(B.scr2_raw - 0.23, 1e-06))}, also_u=True)
    lev = B.peplevel.to_numpy()
    keep5 = np.where(lev == 3)[0]
    drop5 = np.where(lev != 3)[0]
    run('SA5_referral_guard', I1, I2, keep5, also_u=True, key=True)
    dsbp = B.sbp2.to_numpy() - B.sbp1.to_numpy()
    res['SA5_referral_guard']['mean_delta_sbp'] = {'excluded': S.paired_ci(dsbp[drop5], cl[drop5]), 'kept': S.paired_ci(dsbp[keep5], cl[keep5])}
    ex = res['SA5_referral_guard']
    ex['interpretive_d_treatment_feedback'] = bool(abs(ex['theta']['estimate'] - theta0) >= 0.015 and ex['mean_delta_sbp']['excluded']['upper'] < 0 and (ex['mean_delta_sbp']['excluded']['estimate'] < ex['mean_delta_sbp']['kept']['estimate']))

    def extra_player(name, I1x, I2x, F1, F2, player):
        groups = GROUPS[:4] + [player]
        K = 5
        L = []
        for m in range(1 << K):
            H = hybrid(I1x, I2x, m & 15, GROUPS[:4])
            Fm = F2 if m >> 4 & 1 else F1
            if player == 'AGE':
                scr = B.scr2.to_numpy() if m >> 2 & 1 else B.scr1.to_numpy()
                H['egfr'] = M.clamp(M.ckd_epi_2021(scr, Fm['age'], Fm['female']), 15, 140)
            l_ = lp_of(Fm, H, 'base')
            assert scalar_check(Fm, H, 'base', 0, l_, range(n)) < 1e-10
            L.append(l_)
        C = M.category(M.risk_from_lp(np.column_stack(L)))
        D1 = (C != C[:, [0]]).astype(float)
        D2 = np.column_stack([C[:, 31 ^ m] != C[:, 31] for m in range(32)]).astype(float)
        psi = (shapley(D1, K) + shapley(D2, K)) / 2
        ch = (C[:, -1] != C[:, 0]).astype(float)
        assert np.max(np.abs(psi.sum(1) - ch)) < 1e-12
        res[name] = {'theta': S.prop_ci(ch, cl), 'delta_vs_primary_pp': 100 * (ch.mean() - theta0), 'Psi_pp': {g: float(100 * psi[:, k].mean()) for k, g in enumerate(groups)}}
    age2 = np.minimum(np.floor((B.age_months + B.rr_days / 30.4375) / 12), 79).to_numpy(float)
    F2a = {**F, 'age': age2}
    I2a = {**I2, 'egfr': M.clamp(M.ckd_epi_2021(B.scr2, age2, B.female), 15, 140)}
    extra_player('SA6_visit_specific_age', I1, I2a, F, F2a, 'AGE')
    res['SA6_visit_specific_age']['reaching_80_clamped'] = int((np.floor((B.age_months + B.rr_days / 30.4375) / 12) >= 80).sum())
    sm1 = np.where(B.smk_v1_5day.notna(), B.smk_v1_5day, B.smoking).astype(float)
    sm2 = np.where(B.smk_v2_5day.notna(), B.smk_v2_5day, B.smoking).astype(float)
    extra_player('SA7_visit_specific_smoking', I1, I2, {**F, 'smoking': sm1}, {**F, 'smoking': sm2}, 'SMOKING')
    fl = ((np.abs(dsbp) > 40) | (np.abs(B.tc2 - B.tc1) > 80) | (np.abs(B.hdl2 - B.hdl1) > 25) | ~(B.scr2 / B.scr1).between(0.67, 1.5) | (np.abs(B.hba1c2 - B.hba1c1) > 1.5) | B.in_uacr & ~(B.uacr2 / B.uacr1).between(0.1, 10)).to_numpy()
    run('SA8_exclude_flagged', I1, I2, np.where(~fl)[0])
    res['SA8_exclude_flagged']['n_flagged'] = int(fl.sum())
    for nm, m in [('rose', B.rose), ('angina_users', B.angina_user), ('preg_unknown', B.preg_unknown)]:
        run(f'SA9_exclude_{nm}', I1, I2, np.where(~m.to_numpy(bool))[0])
    allx = (B.rose | B.angina_user | B.preg_unknown).to_numpy(bool)
    run('SA9_exclude_all_three', I1, I2, np.where(~allx)[0])
    Ec = E.copy()
    comp = Ec.complete.to_numpy(float)
    ak = np.percentile(Ec.age, [10, 50, 90])

    def design(df):
        return np.column_stack([np.ones(len(df)), S.rcs_basis(df.age.to_numpy(float), ak), df.female, df.dmarethn == 2, df.dmarethn == 3, df.dmarethn == 4, df.phase == 2, df.smoking, df.bptreat, df.srd == 1]).astype(float)
    Xc = design(Ec)
    clE = Ec.cluster.to_numpy()
    chmap = dict(zip(B.SEQN, (M.category(M.risk_from_lp(lp_of(F, I1))) != M.category(M.risk_from_lp(lp_of(F, I2)))).astype(float)))
    ych = np.array([chmap.get(s, np.nan) for s in Ec.SEQN])

    def ipw(wrow=None):
        b = S.logit_fit(Xc, comp, wrow)
        p = 1 / (1 + np.exp(-(Xc @ b)))
        w0 = np.ones(len(comp)) if wrow is None else wrow
        pm = np.sum(w0 * comp) / np.sum(w0)
        sw = np.where(comp == 1, pm / p, 0.0)
        lo, hi = np.percentile(sw[comp == 1], [1, 99])
        sw = np.where(comp == 1, np.clip(sw, lo, hi), 0.0)
        m = comp == 1
        return (float(np.sum((w0 * sw)[m] * ych[m]) / np.sum((w0 * sw)[m])), sw)
    est, sw = ipw()
    Wb = S.boot_weights(clE, 4000, SEED + 10)
    bvals = [ipw(Wb[b])[0] for b in range(4000)]
    std_diff = {}
    for col in ['age', 'female', 'smoking', 'bptreat']:
        a_, b_ = (Ec.loc[Ec.complete, col].astype(float), Ec.loc[~Ec.complete, col].astype(float))
        std_diff[col] = float((a_.mean() - b_.mean()) / math.sqrt((a_.var() + b_.var()) / 2))
    res['SA10_completeness_weighting'] = {'theta_weighted': est, 'ci': S.pct_ci(bvals), 'delta_vs_primary_pp': 100 * (est - theta0), 'key': True, 'abstract_flag_ge_1.5pp': abs(est - theta0) >= 0.015, 'n_eligible': int(len(Ec)), 'n_complete': int(comp.sum()), 'std_diff_complete_vs_not': std_diff, 'weight_range': [float(sw[comp == 1].min()), float(sw[comp == 1].max())]}
    run('SA11_lipid_lowering_proxy', I1, I2, Fx={**F, 'statin': B.lipid_proxy.to_numpy(float)})
    run('SA12_same_examiner', I1, I2, np.where(B.same_examiner.to_numpy(bool))[0])
    run('SA13_risk_rounded_0.1', I1, I2, rounded=True)
    ua_ok = np.where(B.in_uacr.to_numpy() & ~B.ualb1_bdl.to_numpy(bool) & ~B.ualb2_bdl.to_numpy(bool))[0]

    def uacr_shapley(name, idx, I1x, I2x):
        Fx, a, b = (sub(F, idx), sub(I1x, idx), sub(I2x, idx))
        L = np.column_stack([lp_of(Fx, hybrid(a, b, m, GROUPS), 'uacr') for m in range(32)])
        C = M.category(M.risk_from_lp(L))
        D1 = (C != C[:, [0]]).astype(float)
        D2 = np.column_stack([C[:, 31 ^ m] != C[:, 31] for m in range(32)]).astype(float)
        psi = (shapley(D1, 5) + shapley(D2, 5)) / 2
        ch = (C[:, -1] != C[:, 0]).astype(float)
        assert scalar_check(Fx, a, 'uacr', 0, L[:, 0], range(len(idx))) < 1e-10 and scalar_check(Fx, b, 'uacr', 0, L[:, -1], range(len(idx))) < 1e-10
        res[name] = {'theta_U': S.prop_ci(ch, cl[idx]), 'Psi_URINE_pp': float(100 * psi[:, 4].mean()), 'Psi_pp': {g: float(100 * psi[:, k].mean()) for k, g in enumerate(GROUPS)}}
    uacr_shapley('SA14a_exclude_albumin_bdl', ua_ok, I1, I2)
    k79 = np.where(B.uacr1_keep79.notna().to_numpy() & B.uacr2_keep79.notna().to_numpy())[0]
    uacr_shapley('SA14b_keep_ucr_bdl', k79, {**I1, 'uacr': B.uacr1_keep79.to_numpy(float)}, {**I2, 'uacr': B.uacr2_keep79.to_numpy(float)})
    psiU0 = st3['K2_uacr']['Psi']['URINE']['pp']
    pa = res['SA14a_exclude_albumin_bdl']['Psi_URINE_pp']
    res['SA14a_exclude_albumin_bdl']['interpretive_f_sensitive'] = bool(abs(pa - psiU0) >= 0.5 and (abs(pa - psiU0) / abs(psiU0) >= 0.5 if psiU0 else True))
    ch0 = (M.category(M.risk_from_lp(lp_of(F, I1))) != M.category(M.risk_from_lp(lp_of(F, I2)))).astype(float)
    assert abs(ch0.mean() - theta0) < 1e-12
    strata = B.stratum.to_numpy().copy()
    psu = B.cluster.to_numpy()
    th, var, df = S.strat_taylor(ch0, strata, psu)
    from scipy.stats import t as tdist
    se = math.sqrt(var)
    q = tdist.ppf(0.975, df)
    Wb2 = S.boot_weights(cl, 4000, SEED)
    V49, G49 = S.cr1_var(ch0, strata)
    se49 = math.sqrt(max(V49, th * (1 - th) / n)) / (th * (1 - th))
    q49 = tdist.ppf(0.975, G49 - 1)
    from scipy.special import expit, logit
    res['SA15_uncertainty'] = {'wilson': S.wilson(ch0.sum(), n), 'stratified_taylor': {'lower': th - q * se, 'upper': th + q * se, 'df': df}, 'cluster_bootstrap_percentile': S.pct_ci(Wb2 @ ch0 / Wb2.sum(1)), 'cr1_pseudo_stratum': {'lower': float(expit(logit(th) - q49 * se49)), 'upper': float(expit(logit(th) + q49 * se49)), 'G': G49}}
    w_primary = st3['primary']['upper'] - st3['primary']['lower']
    w49 = res['SA15_uncertainty']['cr1_pseudo_stratum']['upper'] - res['SA15_uncertainty']['cr1_pseudo_stratum']['lower']
    res['SA15_uncertainty']['stratum_ci_more_than_20pct_wider'] = w49 > 1.2 * w_primary
    res['completed_utc'] = datetime.now(timezone.utc).isoformat()
    save('stage4_results.json', res)
    print(json.dumps({k: round(100 * v['theta']['estimate'], 2) if 'theta' in v else v.get('theta_weighted') for k, v in res.items() if isinstance(v, dict) and ('theta' in v or 'theta_weighted' in v)}, indent=1))
if __name__ == '__main__':
    main()
