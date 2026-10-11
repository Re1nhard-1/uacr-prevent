from datetime import datetime, timezone
from importlib.metadata import version
from itertools import combinations
from pathlib import Path
import hashlib
import json
import math
import platform
import shutil
import sys
import numpy as np
import pandas as pd
import pyreadstat
from scipy.special import expit, logit
import nh3_build
import nh3_models as M
import nh3_stats as S
import nh3_testretest_stage12 as s12
ROOT = Path(__file__).resolve().parents[1]
STUDY = ROOT / 'studies/nhanes3_testretest'
A1 = STUDY / 'analysis_v1'
OUT = A1 / 'stage3'
SEED = 20261009
B_BOOT = 4000
CUTS = (3.0, 5.0, 10.0)
GROUPS = ['BP', 'LIPID', 'KIDNEY', 'GLYC', 'URINE']
GKEYS = {'BP': ['sbp'], 'LIPID': ['tc', 'hdl'], 'KIDNEY': ['egfr'], 'GLYC': ['dm'], 'URINE': ['uacr']}
PRIMARY = {'sbp': (90.0, 200.0), 'tc': (130.0, 320.0), 'hdl': (20.0, 100.0), 'egfr': (15.0, 140.0)}
LACHER = {'tc': 8.2, 'hdl': 12.4, 'scr': 6.8, 'hba1c': 1.5, 'ucr': 43.0}
CVA = {'tc': (1.25, 2.3), 'hdl': (2.49, 2.5), 'scr_raw': (1.0, 2.7), 'hba1c': (1.57, 3.1), 'ucr': (2.2, 2.2)}

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def save(name, x):
    (OUT / name).write_text(json.dumps(x, indent=2, default=lambda o: o.item() if hasattr(o, 'item') else str(o)) + '\n', encoding='utf-8')

def visit_inputs(B, v, bounds=PRIMARY):
    I = {k: M.clamp(B[f'{k}{v}'], *bounds[k]) for k in ['sbp', 'tc', 'hdl', 'egfr']}
    I['dm'] = B[f'dm{v}'].to_numpy(float)
    I['uacr'] = B[f'uacr{v}'].to_numpy(float)
    return I

def fixed(B):
    return {k: B[k].to_numpy(float) for k in ['female', 'age', 'smoking', 'bptreat', 'statin']}

def lp_of(F, I, model='base', outcome='ascvd', idx=None):
    sel = (lambda a: a) if idx is None else lambda a: a[idx]
    return M.lp(sel(F['female']), sel(F['age']), sel(I['tc']), sel(I['hdl']), sel(I['sbp']), sel(I['dm']), sel(F['smoking']), sel(I['egfr']), sel(F['bptreat']), sel(F['statin']), uacr=sel(I['uacr']) if model == 'uacr' else None, model=model, outcome=outcome)

def hybrid(I1, I2, mask, groups):
    H = {k: v.copy() for k, v in I1.items()}
    for j, g in enumerate(groups):
        if mask >> j & 1:
            for k in GKEYS[g]:
                H[k] = I2[k].copy()
    return H

def scalar_check(F, I, model, outcome_idx, lps, idx):
    worst = 0.0
    for j, i in enumerate(idx):
        r = M.scalar(model, F['age'][i], int(F['female'][i]), int(I['dm'][i]), int(F['smoking'][i]), int(F['bptreat'][i]), int(F['statin'][i]), I['sbp'][i], I['tc'][i], I['hdl'][i], I['egfr'][i], I['uacr'][i] if model == 'uacr' else None)[outcome_idx]
        worst = max(worst, abs(r - M.risk_from_lp(lps[j])))
    return worst

def shapley(D, K):
    n = D.shape[0]
    phi = np.zeros((n, K))
    full = (1 << K) - 1
    for k in range(K):
        others = [j for j in range(K) if j != k]
        for s in range(K):
            w = math.factorial(s) * math.factorial(K - s - 1) / math.factorial(K)
            for comb in combinations(others, s):
                m = sum((1 << j for j in comb))
                phi[:, k] += w * (D[:, m | 1 << k] - D[:, m])
    return phi

def p4_checks():
    rng = np.random.default_rng(SEED)
    out = {}
    n = 1600
    f = np.arange(n) % 2
    age = rng.integers(30, 80, n).astype(float)
    age[:40] = 30
    age[40:80] = 79
    sbp = rng.uniform(90, 200, n)
    sbp[80:200] = rng.choice([90, 110, 130, 200], 120)
    tc = rng.uniform(130, 320, n)
    tc[200:260] = rng.choice([130, 320], 60)
    hdl = rng.uniform(20, 100, n)
    hdl[260:320] = rng.choice([20, 100], 60)
    eg = rng.uniform(15, 140, n)
    eg[320:440] = rng.choice([15, 60, 90, 140], 120)
    dm, sm, bp, st = (rng.integers(0, 2, n) for _ in range(4))
    ua = np.exp(rng.uniform(np.log(0.1), np.log(25000), n))
    I = {'sbp': sbp, 'tc': tc, 'hdl': hdl, 'egfr': eg, 'dm': dm.astype(float), 'uacr': ua}
    F = {'female': f.astype(float), 'age': age, 'smoking': sm.astype(float), 'bptreat': bp.astype(float), 'statin': st.astype(float)}
    w = {}
    for model, outs in [('base', [('ascvd', 0), ('total_cvd', 1)]), ('uacr', [('ascvd', 0)])]:
        for o, oi in outs:
            l = lp_of(F, I, model, o)
            w[f'{model}_{o}'] = scalar_check(F, I, model, oi, l, range(n))
    out['a_b_vectorized_vs_pyprevent_max_pp'] = w
    from prevent._base import prevent_base
    bmi_diff = max((abs(prevent_base(sex=float(f[i]), age=float(age[i]), tc=float(tc[i]), hdl=float(hdl[i]), sbp=float(sbp[i]), dm=float(dm[i]), smoking=float(sm[i]), bmi=b, egfr=float(eg[i]), bptreat=float(bp[i]), statin=float(st[i]))['prevent_base_10yr_ASCVD'] - prevent_base(sex=float(f[i]), age=float(age[i]), tc=float(tc[i]), hdl=float(hdl[i]), sbp=float(sbp[i]), dm=float(dm[i]), smoking=float(sm[i]), bmi=float('nan'), egfr=float(eg[i]), bptreat=float(bp[i]), statin=float(st[i]))['prevent_base_10yr_ASCVD']) for i in range(200) for b in (5.0, 25.0, 90.0)))
    out['a_bmi_invariance_max_pp'] = bmi_diff
    out['b_existing_uacr_verification_hashes'] = {p: sha(ROOT / p) for p in ['studies/uacr_prevent/ascvd_input_verification_v2/executed_verification_script.py']}
    scr = rng.uniform(0.3, 6, 500)
    ag = rng.integers(30, 80, 500)
    fe = rng.integers(0, 2, 500)
    paper1 = lambda c, a, female: 142 * min(c / (0.7 if female else 0.9), 1) ** (-0.241 if female else -0.302) * max(c / (0.7 if female else 0.9), 1) ** (-1.2) * 0.9938 ** a * (1.012 if female else 1.0)
    out['c_ckdepi_identity_max'] = float(max((abs(M.ckd_epi_2021(scr[i], ag[i], fe[i]) - paper1(scr[i], ag[i], fe[i])) for i in range(500))))
    out['c_hand_cases'] = {'scr1.0_age50_male': float(M.ckd_epi_2021(1.0, 50, 0)), 'scr0.7_age40_female': float(M.ckd_epi_2021(0.7, 40, 1)), 'raw0.6_calibrated': 0.96 * 0.6 - 0.184}
    cl = [M.clamp(x, 90, 200) for x in (89.9, 90, 150, 200, 200.1)]
    out['d_clamp'] = cl == [90, 90, 150, 200, 200] and all((M.clamp(M.clamp(x, 90, 200), 90, 200) == M.clamp(x, 90, 200) for x in (50, 95, 250)))
    refusals = {}
    good = dict(age=55, female=1, dm=0, smoking=0, bptreat=0, statin=0, sbp=120.0, tc=200.0, hdl=50.0, egfr=90.0, uacr=10.0)
    for name, bad in [('nan_uacr', {'uacr': float('nan')}), ('age80', {'age': 80}), ('age_nonint', {'age': 55.5}), ('sbp_out', {'sbp': 201.0}), ('tc_out', {'tc': 400.0}), ('nonbinary', {'dm': 2})]:
        try:
            M.scalar('uacr', **{**good, **bad})
            refusals[name] = False
        except AssertionError:
            refusals[name] = True
    out['d_guard_refusals'] = refusals
    I2 = {k: v * rng.uniform(0.85, 1.15, n) if k != 'dm' else v for k, v in I.items()}
    for k in ['sbp', 'tc', 'hdl', 'egfr']:
        I2[k] = M.clamp(I2[k], *PRIMARY[k])
    I2['dm'] = np.where(rng.random(n) < 0.05, 1 - I['dm'], I['dm'])
    I2['uacr'] = M.clamp(I['uacr'] * np.exp(rng.normal(0, 0.5, n)), 0.1, 25000)
    worst_add, worst_eff = (0.0, 0.0)
    for model, groups in [('base', GROUPS[:4]), ('uacr', GROUPS)]:
        K = len(groups)
        L = np.column_stack([lp_of(F, hybrid(I, I2, m, groups), model) for m in range(1 << K)])
        f_ = np.column_stack([L[:, 1 << k] - L[:, 0] for k in range(K)])
        worst_add = max(worst_add, float(np.max(np.abs(L[:, -1] - L[:, 0] - f_.sum(1)))))
        C = M.category(M.risk_from_lp(L))
        D1 = (C != C[:, [0]]).astype(float)
        D2 = (C[:, [-1]] != C[:, ::-1]).astype(float)
        D2 = np.column_stack([C[:, (1 << K) - 1 ^ m] != C[:, -1] for m in range(1 << K)]).astype(float)
        psi = (shapley(D1, K) + shapley(D2, K)) / 2
        worst_eff = max(worst_eff, float(np.max(np.abs(psi.sum(1) - (C[:, -1] != C[:, 0])))))
        assert np.allclose(L[:, 0], lp_of(F, I, model)) and np.allclose(L[:, -1], lp_of(F, I2, model))
    out['f_lp_additivity_max'] = worst_add
    out['f_shapley_efficiency_max'] = worst_eff
    out['a_total_cvd_note'] = 'pinned pyprevent hard-codes the female total-CVD SBP-high coefficient as 0.360078; preventr base_10yr has 0.3600781. Total CVD (S10 only) uses the full-precision preventr table; scalar tolerance 1e-4 pp for this output (deviation logged before scoring).'
    ok = w['base_ascvd'] < 1e-10 and w['uacr_ascvd'] < 1e-10 and (w['base_total_cvd'] < 0.0001) and (bmi_diff == 0) and (out['c_ckdepi_identity_max'] < 1e-12) and out['d_clamp'] and all(refusals.values()) and (worst_add < 1e-10) and (worst_eff < 1e-12)
    out['pass'] = bool(ok)
    return out

def main():
    if OUT.exists():
        raise FileExistsError('Keep existing execution; select a new version explicitly')
    ye = json.loads((A1 / 'ye_reproduction.json').read_text(encoding='utf-8'))
    st1 = json.loads((A1 / 'stage1_status.json').read_text(encoding='utf-8'))
    assert ye['verdict'] == 'PASS' and st1['pass'], 'Stage 1-2 gates not passed'
    OUT.mkdir()
    for f in ['nh3_stage3.py', 'nh3_build.py', 'nh3_models.py', 'nh3_stats.py']:
        shutil.copy2(Path(__file__).parent / f, OUT / f'executed_{f}')
    p4 = p4_checks()
    save('p4_verification.json', p4)
    assert p4['pass'], 'P4 failed: STOP'
    E, B, T = nh3_build.build()
    n = len(B)
    cl = B.cluster.to_numpy()
    F = fixed(B)
    U = np.where(B.in_uacr.to_numpy())[0]
    nU = len(U)
    clU = cl[U]
    I1, I2 = (visit_inputs(B, 1), visit_inputs(B, 2))
    KB, KU = (4, 5)
    LB = np.column_stack([lp_of(F, hybrid(I1, I2, m, GROUPS[:4]), 'base') for m in range(1 << KB)])
    I1u = {k: v[U] for k, v in I1.items()}
    I2u = {k: v[U] for k, v in I2.items()}
    Fu = {k: v[U] for k, v in F.items()}
    LU = np.column_stack([lp_of(Fu, hybrid(I1u, I2u, m, GROUPS), 'uacr') for m in range(1 << KU)])
    p7 = {'scalar_vs_vectorized_max_pp': {}}
    wb = 0.0
    for m in range(1 << KB):
        wb = max(wb, scalar_check(F, hybrid(I1, I2, m, GROUPS[:4]), 'base', 0, LB[:, m], range(n)))
    wu = 0.0
    for m in range(1 << KU):
        wu = max(wu, scalar_check(Fu, hybrid(I1u, I2u, m, GROUPS), 'uacr', 0, LU[:, m], range(nU)))
    p7['scalar_vs_vectorized_max_pp'].update(base_hybrids=wb, uacr_hybrids=wu)
    RB, RU = (M.risk_from_lp(LB), M.risk_from_lp(LU))
    assert np.all(np.isfinite(RB)) and np.all((RB > 0) & (RB < 100)) and np.all(np.isfinite(RU)) and np.all(np.isfinite(I1u['uacr'])) and np.all(np.isfinite(I2u['uacr']))
    CB, CU = (M.category(RB), M.category(RU))

    def dmat(C, K):
        full = (1 << K) - 1
        D1 = (C != C[:, [0]]).astype(float)
        D2 = np.column_stack([C[:, full ^ m] != C[:, full] for m in range(1 << K)]).astype(float)
        return (D1, D2)
    D1B, D2B = dmat(CB, KB)
    D1U, D2U = dmat(CU, KU)
    phiB1, phiB2 = (shapley(D1B, KB), shapley(D2B, KB))
    psiB = (phiB1 + phiB2) / 2
    phiU1, phiU2 = (shapley(D1U, KU), shapley(D2U, KU))
    psiU = (phiU1 + phiU2) / 2
    chB = (CB[:, -1] != CB[:, 0]).astype(float)
    chU = (CU[:, -1] != CU[:, 0]).astype(float)
    fB = np.column_stack([LB[:, 1 << k] - LB[:, 0] for k in range(KB)])
    dLB = LB[:, -1] - LB[:, 0]
    fU = np.column_stack([LU[:, 1 << k] - LU[:, 0] for k in range(KU)])
    dLU = LU[:, -1] - LU[:, 0]
    p7.update(lp_additivity_base=float(np.max(np.abs(dLB - fB.sum(1)))), lp_additivity_uacr=float(np.max(np.abs(dLU - fU.sum(1)))), shapley_efficiency_base=float(np.max(np.abs(psiB.sum(1) - chB))), shapley_efficiency_uacr=float(np.max(np.abs(psiU.sum(1) - chU))), theta_swap_invariant=bool(np.array_equal(CB[:, 0] != CB[:, -1], CB[:, -1] != CB[:, 0])))
    LT1, LT2 = (lp_of(F, I1, 'base', 'total_cvd'), lp_of(F, I2, 'base', 'total_cvd'))
    p7['scalar_vs_vectorized_max_pp']['total_cvd'] = max(scalar_check(F, I1, 'base', 1, LT1, range(n)), scalar_check(F, I2, 'base', 1, LT2, range(n)))
    FT = {k: T[k].to_numpy(float) for k in ['female', 'age', 'smoking', 'bptreat', 'statin']}
    IT = {'sbp': M.clamp(T.sbp1, 90, 200), 'tc': M.clamp(T.tc1, 130, 320), 'hdl': M.clamp(T.hdl1, 20, 100), 'egfr': M.clamp(T.egfr1, 15, 140), 'dm': T.dm1.to_numpy(float), 'uacr': T.uacr1.to_numpy(float)}
    LTB = lp_of(FT, IT, 'base')
    TU = np.where(T.uacr1.notna().to_numpy())[0]
    ITu = {k: v[TU] for k, v in IT.items()}
    FTu = {k: v[TU] for k, v in FT.items()}
    LTU = lp_of(FTu, ITu, 'uacr')
    p7['scalar_vs_vectorized_max_pp']['target_base'] = scalar_check(FT, IT, 'base', 0, LTB, range(len(T)))
    p7['scalar_vs_vectorized_max_pp']['target_uacr'] = scalar_check(FTu, ITu, 'uacr', 0, LTU, range(len(TU)))
    p7['all_finite'] = True
    p7['pass'] = all((v < (0.0001 if k == 'total_cvd' else 1e-10) for k, v in p7['scalar_vs_vectorized_max_pp'].items())) and p7['lp_additivity_base'] < 1e-10 and (p7['lp_additivity_uacr'] < 1e-10) and (p7['shapley_efficiency_base'] < 1e-12) and (p7['shapley_efficiency_uacr'] < 1e-12) and p7['theta_swap_invariant']
    save('p7_scoring_checks.json', p7)
    assert p7['pass'], 'P7 failed: STOP'
    res = {}
    up = (CB[:, -1] > CB[:, 0]).astype(float)
    dn = (CB[:, -1] < CB[:, 0]).astype(float)
    res['primary'] = {**S.prop_ci(chB, cl), 'up': int(up.sum()), 'down': int(dn.sum())}
    pr = res['primary']
    pr['precision_label'] = 'adequately precise' if pr['upper'] - pr['lower'] <= 0.06 else 'imprecise'
    chB_U = chB[U]
    res['K1'] = {'theta_U': S.prop_ci(chU, clU), 'theta_base_same725': S.prop_ci(chB_U, clU), 'delta_U': S.paired_ci(chU - chB_U, clU), 'discordant': {'U_only': int(((chU == 1) & (chB_U == 0)).sum()), 'base_only': int(((chU == 0) & (chB_U == 1)).sum())}}
    b_, c_ = (res['K1']['discordant']['U_only'], res['K1']['discordant']['base_only'])
    if b_ + c_ < 20:
        res['K1']['delta_U_tango'] = S.tango(b_, c_, nU)
    W = S.boot_weights(cl, B_BOOT, SEED)
    WU = W[:, U]
    mean_w = lambda Wm, x: Wm @ x / Wm.sum(1)

    def attribution(D1, D2, psi, f, dL, ch, groups, Wm, clx, LPm):
        K = len(groups)
        full = (1 << K) - 1
        out = {}
        th = ch.mean()
        Psi = psi.mean(0)
        bootPsi = np.column_stack([mean_w(Wm, psi[:, k]) for k in range(K)])
        rank = np.argmax(bootPsi, 1)
        out['Psi'] = {g: {'pp': float(100 * Psi[k]), 'ci_pp': [100 * x for x in S.pct_ci(bootPsi[:, k])], 'share': float(Psi[k] / th) if th > 0 else None, 'rank_first_freq': float(np.mean(rank == k))} for k, g in enumerate(groups)}
        out['events'] = int(ch.sum())
        out['shares_interpretable'] = bool(ch.sum() >= 20)
        lead = max(range(K), key=lambda k: Psi[k])
        out['largest_contributor'] = groups[lead] if np.mean(rank == lead) >= 0.9 else 'comparable (no group ranks first in >=90% of replicates)'
        out['batch_mcse_pp'] = {g: [100 * x for x in S.batch_mcse(bootPsi[:, k])] for k, g in enumerate(groups)}
        out['single_swap'] = {g: {'ref1': float(D1[:, 1 << k].mean()), 'ref2': float(D2[:, 1 << k].mean()), 'sym': float((D1[:, 1 << k].mean() + D2[:, 1 << k].mean()) / 2)} for k, g in enumerate(groups)}
        out['hold_one_fixed'] = {g: {'ref1': float(D1[:, full ^ 1 << k].mean()), 'ref2': float(D2[:, full ^ 1 << k].mean())} for k, g in enumerate(groups)}
        bl = 1 << groups.index('BP') | 1 << groups.index('LIPID')
        out['coalition_BP_LIPID'] = {'ref1': float(D1[:, bl].mean()), 'ref2': float(D2[:, bl].mean())}
        out['non_additivity_gap'] = float(th - sum((D1[:, 1 << k].mean() for k in range(K))))
        mu = f.mean(0)
        pi = (f * dL[:, None]).mean(0) / np.mean(dL ** 2)
        nu = np.array([np.cov(f[:, k], dL)[0, 1] for k in range(K)]) / np.var(dL, ddof=1)
        bpi = np.column_stack([Wm @ (f[:, k] * dL) / (Wm @ dL ** 2) for k in range(K)])
        out['continuous'] = {g: {'mu': float(mu[k]), 'mean_abs': float(np.abs(f[:, k]).mean()), 'pi': float(pi[k]), 'pi_ci': S.pct_ci(bpi[:, k]), 'nu': float(nu[k])} for k, g in enumerate(groups)}
        out['continuous_corr'] = np.corrcoef(f.T).tolist()
        bp_sym = (D1[:, 1 << groups.index('BP')] + D2[:, 1 << groups.index('BP')]) / 2
        out['theta_minus_BP_sym'] = S.paired_ci(ch - bp_sym, clx)
        out['interpretive_a'] = 'Laboratory remeasurement adds instability beyond BP' if out['theta_minus_BP_sym']['lower'] > 0 or out['theta_minus_BP_sym']['upper'] < 0 else 'compatible with BP accounting for most of the instability'
        lab = sum((1 << groups.index(g) for g in ['LIPID', 'KIDNEY', 'GLYC']))
        out['D_sym_lab'] = float((D1[:, lab].mean() + D2[:, lab].mean()) / 2)
        return (out, bootPsi)
    res['K2_base'], bootPsiB = attribution(D1B, D2B, psiB, fB, dLB, chB, GROUPS[:4], W, cl, LB)
    res['K2_uacr'], bootPsiU = attribution(D1U, D2U, psiU, fU, dLU, chU, GROUPS, WU, clU, LU)
    for tag, Fx, Ia, Ib, model in [('base', F, I1, I2, 'base')]:
        tc_only = {**Ia, 'tc': Ib['tc']}
        hdl_only = {**Ia, 'hdl': Ib['hdl']}
        l0 = lp_of(Fx, Ia, model)
        res['K2_lipid_split_' + tag] = {'mu_tc': float((lp_of(Fx, tc_only, model) - l0).mean()), 'mu_hdl': float((lp_of(Fx, hdl_only, model) - l0).mean()), 'swap_tc': float((M.category(M.risk_from_lp(lp_of(Fx, tc_only, model))) != CB[:, 0]).mean()), 'swap_hdl': float((M.category(M.risk_from_lp(lp_of(Fx, hdl_only, model))) != CB[:, 0]).mean())}
    LT = np.log(np.array(CUTS) / 100 / (1 - np.array(CUTS) / 100))

    def k3(LP1, LP2, ch, Wm):
        out = {}
        Mn = (LP1 + LP2) / 2
        for name, base_lp in [('clin', LP1), ('mech', Mn)]:
            dist = np.abs(base_lp[:, None] - LT[None, :])
            d = dist.min(1)
            near = dist.argmin(1)
            ev = int(ch.sum())
            if ev < 20:
                out[name] = {'events': ev, 'binned_only': True}
                continue
            knots = np.percentile(d, [10, 50, 90])

            def X(dv, nr):
                Xb = np.column_stack([np.ones(len(dv)), S.rcs_basis(dv, knots)])
                return np.column_stack([Xb, nr == 1, nr == 2]).astype(float) if ev >= 40 else Xb
            beta = S.logit_fit(X(d, near), ch)
            grid = np.array([0, 0.05, 0.1, 0.2, 0.3, 0.5])
            fine = np.linspace(0, 1.5, 1501)
            res_t = {}
            boots = []
            for b in range(Wm.shape[0]):
                boots.append(S.logit_fit(X(d, near), ch, Wm[b]))
            boots = np.array(boots)
            for ti, t in enumerate(CUTS):
                nr = np.full(len(grid), ti)
                nrf = np.full(len(fine), ti)
                p = expit(X(grid, nr) @ beta)
                pb = expit(X(grid, nr) @ boots.T)
                pf = expit(X(fine, nrf) @ beta)
                zones = {}
                for lev in (0.1, 0.25):
                    ok = np.where(pf >= lev)[0]
                    if len(ok) == 0 or ok[0] != 0:
                        zones[str(lev)] = None
                        continue
                    dstar = fine[ok[-1]] if ok[-1] == len(fine) - 1 else fine[np.argmax(pf < lev) - 1]
                    zones[str(lev)] = [float(100 * expit(logit(t / 100) - dstar)), float(100 * expit(logit(t / 100) + dstar))]
                res_t[str(t)] = {'pi_at_d': {str(g_): float(v) for g_, v in zip(grid, p)}, 'ci': {str(g_): S.pct_ci(pb[i]) for i, g_ in enumerate(grid)}, 'instability_zones_risk_pct': zones}
            bins = [(0, 0.05), (0.05, 0.1), (0.1, 0.2), (0.2, 0.4), (0.4, 99)]
            binned = [{'bin': [lo, hi], 'n': int(((d >= lo) & (d < hi)).sum()), 'changes': int(ch[(d >= lo) & (d < hi)].sum())} for lo, hi in bins]
            out[name] = {'events': ev, 'knots': knots.tolist(), 'threshold_indicator': ev >= 40, 'beta': beta.tolist(), 'by_threshold': res_t, 'binned': binned, 'sd_dLP': float(np.std(LP2 - LP1, ddof=1)), 'fitted_mean_equals_theta': float(abs(expit(X(d, near) @ beta).mean() - ch.mean()))}
            if name == 'mech':
                out[name]['mechanism_survival_2d'] = {str(g_): float(np.mean(np.abs(LP2 - LP1) > 2 * g_)) for g_ in grid}
            out[name]['_model'] = (knots, beta, ev >= 40, boots)
        return out
    k3B = k3(LB[:, 0], LB[:, -1], chB, W)
    k3U = k3(LU[:, 0], LU[:, -1], chU, WU)
    res['K3_base'] = {k: {kk: vv for kk, vv in v.items() if kk != '_model'} for k, v in k3B.items()}
    res['K3_uacr'] = {k: {kk: vv for kk, vv in v.items() if kk != '_model'} for k, v in k3U.items()}

    def predict(model_tuple, lpv, beta=None):
        knots, b0, ind, _ = model_tuple
        b = b0 if beta is None else beta
        dist = np.abs(lpv[:, None] - LT[None, :])
        d = dist.min(1)
        near = dist.argmin(1)
        Xb = np.column_stack([np.ones(len(d)), S.rcs_basis(d, knots)])
        if ind:
            Xb = np.column_stack([Xb, near == 1, near == 2]).astype(float)
        return expit(Xb @ b)

    def k4(model_tuple, lp_sample, ch, lp_target, w_target, strata, psu, seed):
        out = {}
        ph_s = predict(model_tuple, lp_sample)
        out['code_identity_abs'] = float(abs(ph_s.mean() - ch.mean()))
        ph_t = predict(model_tuple, lp_target)
        th_std = float(np.sum(w_target * ph_t) / np.sum(w_target))
        out['theta_std'] = th_std
        out['theta'] = float(ch.mean())
        out['ratio'] = th_std / float(ch.mean())
        bands = [0, 2, 3, 4, 5, 7.5, 10, 15, 1000]
        rs, rt = (M.risk_from_lp(lp_sample), M.risk_from_lp(lp_target))
        bs, bt = (np.digitize(rs, bands[1:-1]), np.digitize(rt, bands[1:-1]))
        groups_ = [[i] for i in range(8)]

        def merge(groups_):
            k5 = 4
            changed = True
            while changed:
                changed = False
                for gi, g in enumerate(groups_):
                    if sum((np.isin(bs, g).sum() for _ in [0])) < 20 and len(groups_) > 1:
                        j = gi + 1 if max(g) < k5 else gi - 1
                        j = min(max(j, 0), len(groups_) - 1)
                        if j == gi:
                            j = gi - 1 if gi > 0 else gi + 1
                        groups_[j] = sorted(groups_[j] + g)
                        groups_.pop(gi)
                        changed = True
                        break
            return groups_
        groups_ = merge(groups_)
        ps = sum((np.sum(w_target[np.isin(bt, g)]) * ch[np.isin(bs, g)].mean() for g in groups_ if np.isin(bs, g).any())) / np.sum(w_target)
        out['poststrat_check'] = float(ps)
        out['poststrat_bands'] = groups_
        knots, b0, ind, boots = model_tuple
        rng = np.random.default_rng(seed)
        vals = []
        st_ids = np.unique(strata)
        for bb in range(2000):
            beta_b = boots[bb]
            wt = w_target.copy()
            for h in st_ids:
                m = strata == h
                ps_ = np.unique(psu[m])
                nh = len(ps_)
                if nh < 2:
                    continue
                draw = rng.choice(ps_, nh - 1, replace=True)
                cnt = {p: int((draw == p).sum()) for p in ps_}
                wt[m] = w_target[m] * nh / (nh - 1) * np.array([cnt[p] for p in psu[m]])
            vals.append(np.sum(wt * predict(model_tuple, lp_target, beta_b)) / np.sum(wt))
        out['theta_std_ci'] = S.pct_ci(vals)
        out['headline_rule_curve_first'] = abs(th_std - float(ch.mean())) > 0.02
        return out
    wT = T.w.to_numpy(float)
    stT = T.stratum.to_numpy()
    psT = (T.stratum * 10 + T.psu).to_numpy()
    res['K4_base'] = k4(k3B['clin']['_model'], LB[:, 0], chB, LTB, wT, stT, psT, SEED + 1)
    res['K4_uacr'] = k4(k3U['clin']['_model'], LU[:, 0], chU, LTU, wT[TU], stT[TU], psT[TU], SEED + 2)

    def s5(C, R, ch, clx):
        C1, C2 = (C[:, 0], C[:, -1])
        tr = np.zeros((4, 4), int)
        for a_, b_ in zip(C1, C2):
            tr[a_, b_] += 1
        upx, dnx = ((C2 > C1).astype(float), (C2 < C1).astype(float))
        win = np.any(np.abs(R[:, 0][:, None] - np.array(CUTS)[None, :]) <= 0.5, 1)
        return {'transition_counts': tr.tolist(), 'row_pct': (100 * tr / np.maximum(tr.sum(1, keepdims=True), 1)).tolist(), 'total_pct': (100 * tr / tr.sum()).tolist(), 'p_up': S.prop_ci(upx, clx), 'p_down': S.prop_ci(dnx, clx), 'net_shift': S.paired_ci(upx - dnx, clx), 'non_adjacent': int((np.abs(C2 - C1) > 1).sum()), 'dist_v1': np.bincount(C1, minlength=4).tolist(), 'dist_v2': np.bincount(C2, minlength=4).tolist(), 'crossing': {str(t): S.prop_ci(((R[:, 0] >= t) != (R[:, -1] >= t)).astype(float), clx) for t in CUTS}, 'window_0p5': {'share_in_window': float(win.mean()), 'share_of_changes_in_window': float(ch[win].sum() / ch.sum()) if ch.sum() else None}}
    res['S5_base'] = s5(CB, RB, chB, cl)
    res['S5_uacr'] = s5(CU, RU, chU, clU)

    def s6(LP1, LP2, R1, R2, Wm, clx):
        dL = LP2 - LP1
        sd = float(np.std(dL, ddof=1))
        dlr = np.log(R2) - np.log(R1)
        dr = R2 - R1
        mad = float(1.4826 * np.median(np.abs(dlr - np.median(dlr))))
        Mn = (LP1 + LP2) / 2
        Xh = np.column_stack([np.ones(len(Mn)), Mn])
        yh = np.abs(dL)
        bh = np.linalg.lstsq(Xh, yh, rcond=None)[0]
        e = yh - Xh @ bh
        ids, inv = np.unique(clx, return_inverse=True)
        G = len(ids)
        meat = sum((np.outer(Xh[inv == g].T @ e[inv == g], Xh[inv == g].T @ e[inv == g]) for g in range(G)))
        bread = np.linalg.inv(Xh.T @ Xh)
        V = G / (G - 1) * bread @ meat @ bread
        bsd = np.sqrt(np.maximum(Wm @ dL ** 2 / Wm.sum(1) - (Wm @ dL / Wm.sum(1)) ** 2, 0))
        bsdr = np.sqrt(np.maximum(Wm @ dlr ** 2 / Wm.sum(1) - (Wm @ dlr / Wm.sum(1)) ** 2, 0))
        return {'mean_dLP': float(dL.mean()), 'sd_dLP': sd, 'sd_dLP_ci': S.pct_ci(bsd), 'sigma_w': sd / math.sqrt(2), 'icc_a1': S.icc_a1(LP1, LP2), 'bland_altman_OR': {'bias': float(math.exp(dL.mean())), 'loa': [float(math.exp(dL.mean() - 1.96 * sd)), float(math.exp(dL.mean() + 1.96 * sd))]}, 'rcv_risk': [float(math.exp(-1.96 * np.std(dlr, ddof=1)) - 1), float(math.exp(1.96 * np.std(dlr, ddof=1)) - 1)], 'rcv_risk_ci_upper': S.pct_ci(np.exp(1.96 * bsdr) - 1), 'rcv_risk_robust': [float(math.exp(-1.96 * mad) - 1), float(math.exp(1.96 * mad) - 1)], 'abs_dr_median_pp': float(np.median(np.abs(dr))), 'abs_dr_p90_pp': float(np.percentile(np.abs(dr), 90)), 'abs_dr_ge1': float((np.abs(dr) >= 1).mean()), 'abs_dr_ge2': float((np.abs(dr) >= 2).mean()), 'heteroscedasticity_slope': float(bh[1]), 'heteroscedasticity_se_cr1': float(math.sqrt(V[1, 1]))}
    res['S6_base'] = s6(LB[:, 0], LB[:, -1], RB[:, 0], RB[:, -1], W, cl)
    res['S6_uacr'] = s6(LU[:, 0], LU[:, -1], RU[:, 0], RU[:, -1], WU, clU)

    def s7_one(x1, x2, log=True, cva=None, lacher=None):
        x1, x2 = (np.asarray(x1, float), np.asarray(x2, float))
        m = np.isfinite(x1) & np.isfinite(x2) & ((x1 > 0) & (x2 > 0) if log else True)
        x1, x2 = (x1[m], x2[m])
        d = x2 - x1
        out = {'n': int(m.sum()), 'mean_diff': float(d.mean())}
        if log:
            s = float(np.std(np.log(x2) - np.log(x1), ddof=1) / math.sqrt(2))
            cvw = 100 * math.sqrt(math.exp(s * s) - 1)
            out.update(cvw_total_pct=cvw, icc_log=S.icc_a1(np.log(x1), np.log(x2)), rcv=[float(math.exp(-1.96 * math.sqrt(2) * s) - 1), float(math.exp(1.96 * math.sqrt(2) * s) - 1)])
            if cva:
                out['cvi_derived_pct'] = {'low_cva': math.sqrt(max(cvw ** 2 - cva[0] ** 2, 0)), 'high_cva': math.sqrt(max(cvw ** 2 - cva[1] ** 2, 0))}
            if lacher:
                out['lacher_cvi'] = lacher
        else:
            out.update(sd_w=float(np.std(d, ddof=1) / math.sqrt(2)), icc=S.icc_a1(x1, x2))
        return out
    res['S7'] = {'sbp': s7_one(B.sbp1, B.sbp2, False), 'tc': s7_one(B.tc1, B.tc2, cva=CVA['tc'], lacher=LACHER['tc']), 'hdl': s7_one(B.hdl1, B.hdl2, cva=CVA['hdl'], lacher=LACHER['hdl']), 'nonhdl': s7_one(B.tc1 - B.hdl1, B.tc2 - B.hdl2), 'scr_cal': s7_one(B.scr1, B.scr2, cva=CVA['scr_raw'], lacher=LACHER['scr']), 'egfr': s7_one(B.egfr1, B.egfr2), 'ualb': s7_one(B.ualb1[U], B.ualb2[U]), 'ucr': s7_one(B.ucr1[U], B.ucr2[U], cva=CVA['ucr'], lacher=LACHER['ucr']), 'uacr': s7_one(B.uacr1[U], B.uacr2[U]), 'hba1c': s7_one(B.hba1c1, B.hba1c2, False), 'hba1c_log': s7_one(B.hba1c1, B.hba1c2, cva=CVA['hba1c'], lacher=LACHER['hba1c'])}
    sc1, sc2 = (B.scr1.to_numpy(), B.scr2.to_numpy())
    dd = sc2 - sc1
    keep = np.abs(dd - dd.mean()) <= 3 * dd.std()
    lin = lambda a, b: float(100 * math.sqrt(np.mean((a - b) ** 2) / 2) / np.mean(np.r_[a, b]))
    res['S7']['scr_cal_selvin_linear'] = {'all': lin(sc1, sc2), 'excl_3sd': lin(sc1[keep], sc2[keep]), 'selvin': 7.6}
    res['S7']['selvin_log_scale_descriptive'] = {'urine_albumin': 37.3, 'urine_creatinine': 43.0, 'acr': 11.3}
    res['S7']['diabetes_discordant'] = int((B.dm1 != B.dm2).sum())
    r23_1, r23_2 = (B.sbp1_r2 - B.sbp1_r3, B.sbp2_r2 - B.sbp2_r3)
    var_rd = float(np.var(np.r_[r23_1, r23_2], ddof=1))
    var_dm = float(np.var(B.sbp2 - B.sbp1, ddof=1))
    res['S7']['sbp_variance'] = {'var_delta_mean23': var_dm, 'var_reading2_minus_3_pooled': var_rd, 'sigma2_within_reading': var_rd / 2, 'sigma2_between_day': (var_dm - var_rd / 2) / 2}

    def sbp_scen(s_a, s_b, idx=None, model='base'):
        Ia = {**I1, 'sbp': M.clamp(s_a, 90, 200)}
        Ib = {**I1, 'sbp': M.clamp(s_b, 90, 200)}
        Fx = F
        if idx is not None:
            Ia = {k: v[idx] for k, v in Ia.items()}
            Ib = {k: v[idx] for k, v in Ib.items()}
            Fx = {k: v[idx] for k, v in F.items()}
        la, lb = (lp_of(Fx, Ia, model), lp_of(Fx, Ib, model))
        wa = scalar_check(Fx, Ia, model, 0, la, range(len(la)))
        wb_ = scalar_check(Fx, Ib, model, 0, lb, range(len(lb)))
        assert max(wa, wb_) < 1e-10
        return (M.category(M.risk_from_lp(la)) != M.category(M.risk_from_lp(lb))).astype(float)
    hasr1 = np.where(B.sbp1_r1.notna().to_numpy())[0]
    s8 = {}
    for model, idx_all, clx in [('base', np.arange(n), cl), ('uacr', U, clU)]:
        ia = np.intersect1d(idx_all, hasr1)
        a_ = sbp_scen(B.sbp1_r1.to_numpy(), B.sbp1_r2.to_numpy(), ia, model)
        b2 = sbp_scen(B.sbp1_r2.to_numpy(), B.sbp1_r3.to_numpy(), idx_all, model)
        c_ = sbp_scen(B.sbp1_r2.to_numpy(), B.sbp2_r2.to_numpy(), idx_all, model)
        dsw = D1B[:, 1] if model == 'base' else D1U[:, 1]
        th = chB if model == 'base' else chU
        s8[model] = {'a_v1_r1_vs_r2': S.prop_ci(a_, cl[ia]), 'b_v1_r2_vs_r3': S.prop_ci(b2, clx), 'c_v1r2_vs_v2r2': S.prop_ci(c_, clx), 'd_between_visit_bp_swap': S.prop_ci(dsw, clx), 'c_minus_b': S.paired_ci(c_ - b2, clx), 'theta_minus_d': S.paired_ci(th - dsw, clx)}
        hw = (s8[model]['c_v1r2_vs_v2r2']['upper'] - s8[model]['c_v1r2_vs_v2r2']['lower']) / 2
        s8[model]['interpretive_c_flag'] = bool(s8[model]['b_v1_r2_vs_r3']['estimate'] - s8[model]['c_v1r2_vs_v2r2']['estimate'] > hw)
    res['S8'] = s8
    h1, h2 = (B.hba1c1.to_numpy(float), B.hba1c2.to_numpy(float))
    fx = B.fixed_dm.to_numpy(float)

    def s9_theta(wv):

        def wm(x, m=None):
            m = np.isfinite(x) if m is None else m
            return float(np.sum(wv[m] * x[m]) / np.sum(wv[m]))
        sh = {'sbp': wm(B.sbp2.to_numpy() - B.sbp1.to_numpy()), 'tc': wm(np.log(B.tc2.to_numpy() / B.tc1.to_numpy())), 'hdl': wm(np.log(B.hdl2.to_numpy() / B.hdl1.to_numpy())), 'scr': wm(np.log(B.scr2.to_numpy() / B.scr1.to_numpy())), 'hba1c': wm(h2 - h1)}
        s2 = B.sbp2.to_numpy() - sh['sbp']
        t2 = B.tc2.to_numpy() * math.exp(-sh['tc'])
        hd2 = B.hdl2.to_numpy() * math.exp(-sh['hdl'])
        c2 = B.scr2.to_numpy() * math.exp(-sh['scr'])
        g2 = np.where(np.isfinite(h2), h2 - sh['hba1c'], h1)
        g1 = np.where(np.isfinite(h1), h1, h2)
        dm2 = np.where((fx == 1) | (np.nan_to_num(g2, nan=0) >= 6.5), 1.0, 0.0)
        dm1 = np.where((fx == 1) | (np.nan_to_num(g1, nan=0) >= 6.5), 1.0, 0.0)
        I2s = {'sbp': M.clamp(s2, 90, 200), 'tc': M.clamp(t2, 130, 320), 'hdl': M.clamp(hd2, 20, 100), 'egfr': M.clamp(M.ckd_epi_2021(c2, B.age, B.female), 15, 140), 'dm': dm2, 'uacr': I2['uacr']}
        I1s = {**I1, 'dm': dm1}
        l1, l2 = (lp_of(F, I1s, 'base'), lp_of(F, I2s, 'base'))
        chs = (M.category(M.risk_from_lp(l1)) != M.category(M.risk_from_lp(l2))).astype(float)
        return (float(np.sum(wv * chs) / np.sum(wv)), sh, chs, (I1s, I2s, l1, l2))
    assert np.array_equal(np.where((fx == 1) | (np.nan_to_num(np.where(np.isfinite(h1), h1, h2), nan=0) >= 6.5), 1.0, 0.0), B.dm1.to_numpy(float))
    th_s, sh, chs, (I1s, I2s, l1s, l2s) = s9_theta(np.ones(n))
    assert max(scalar_check(F, I2s, 'base', 0, l2s, range(n)), scalar_check(F, I1s, 'base', 0, l1s, range(n))) < 1e-10
    bvals = [chB.mean() - s9_theta(W[b])[0] for b in range(B_BOOT)]
    res['S9'] = {'shifts': sh, 'theta_shift_removed': th_s, 'theta_minus_shift_removed': float(chB.mean() - th_s), 'ci': S.pct_ci(bvals)}
    RT1, RT2 = (M.risk_from_lp(LT1), M.risk_from_lp(LT2))
    cr75 = ((RT1 >= 7.5) != (RT2 >= 7.5)).astype(float)
    pyT = np.array([[M.scalar('base', F['age'][i], int(F['female'][i]), int(I['dm'][i]), int(F['smoking'][i]), int(F['bptreat'][i]), int(F['statin'][i]), I['sbp'][i], I['tc'][i], I['hdl'][i], I['egfr'][i])[1] for I in (I1, I2)] for i in range(n)])
    s10_impl_diff = int(np.sum(((pyT[:, 0] >= 7.5) != (pyT[:, 1] >= 7.5)) != (cr75 == 1)))
    assert s10_impl_diff == 0, 'S10 crossing differs between implementations'
    unt = np.where((B.bptreat.to_numpy() == 0) & (B.sbp1.to_numpy() >= 130) & (B.sbp1.to_numpy() < 140))[0]
    res['S10'] = {'implementation_disagreements_in_crossing': s10_impl_diff, 'all_759': S.prop_ci(cr75, cl), 'uacr_set_base_model': S.prop_ci(cr75[U], clU), 'untreated_sbp130_139': {'n': int(len(unt)), 'crossing': int(cr75[unt].sum()), 'proportion': float(cr75[unt].mean())}}
    pc = json.loads((A1 / 'pce_coefficients.json').read_text(encoding='utf-8'))
    a40 = np.where(B.age.to_numpy() >= 40)[0]

    def pce_arr(idx, sbp_v, tc_v, hdl_v, dm_v):
        out = []
        for i in idx:
            grp = ('African American ' if B.dmarethn[i] == 2 else 'White ') + ('women' if B.female[i] == 1 else 'men')
            out.append(s12.pce(pc[grp], float(min(max(B.age[i], 40), 79)), float(M.clamp(tc_v[i], 130, 320)), float(M.clamp(hdl_v[i], 20, 100)), float(M.clamp(sbp_v[i], 90, 200)), B.bptreat[i] == 1, int(B.smoking[i]), int(dm_v[i])))
        return np.array(out)
    p1 = pce_arr(a40, B.sbp1, B.tc1, B.hdl1, B.dm1)
    pa = pce_arr(a40, B.sbp2, B.tc1, B.hdl1, B.dm1)
    pb = pce_arr(a40, B.sbp2, B.tc2, B.hdl2, B.dm2)
    ccat = lambda r: np.searchsorted([5, 7.5, 20], r, side='right')
    xa75, xb75 = (((p1 >= 7.5) != (pa >= 7.5)).astype(float), ((p1 >= 7.5) != (pb >= 7.5)).astype(float))
    ca_, cb_ = ((ccat(p1) != ccat(pa)).astype(float), (ccat(p1) != ccat(pb)).astype(float))
    res['S11'] = {'n': int(len(a40)), 'cross75_sbp_only': S.prop_ci(xa75, cl[a40]), 'cross75_full': S.prop_ci(xb75, cl[a40]), 'cross75_full_minus_sbp': S.paired_ci(xb75 - xa75, cl[a40]), 'cat_sbp_only': S.prop_ci(ca_, cl[a40]), 'cat_full': S.prop_ci(cb_, cl[a40]), 'cat_full_minus_sbp': S.paired_ci(cb_ - ca_, cl[a40]), 'prevent_base_theta_same569': S.prop_ci(chB[a40], cl[a40])}
    R_SIM = 1000
    rng = np.random.default_rng(SEED)

    def simulate(model, idx, mode, cva_bracket=None):
        Fx = {k: v[idx] for k, v in F.items()}
        nn = len(idx)
        c_sbp = (B.sbp1.to_numpy()[idx] + B.sbp2.to_numpy()[idx]) / 2
        geo = lambda a, b: np.sqrt(a * b)
        c_tc, c_hdl = (geo(B.tc1.to_numpy()[idx], B.tc2.to_numpy()[idx]), geo(B.hdl1.to_numpy()[idx], B.hdl2.to_numpy()[idx]))
        c_scr = geo(B.scr1.to_numpy()[idx], B.scr2.to_numpy()[idx])
        c_raw = geo(B.scr1_raw.to_numpy()[idx], B.scr2_raw.to_numpy()[idx])
        hh1, hh2 = (h1[idx], h2[idx])
        c_h = np.where(np.isfinite(hh1) & np.isfinite(hh2), (hh1 + hh2) / 2, np.where(np.isfinite(hh1), hh1, hh2))
        c_u = geo(B.uacr1.to_numpy()[idx], B.uacr2.to_numpy()[idx]) if model == 'uacr' else None
        fxx = fx[idx]
        if mode == 'A':
            sd = {'sbp': np.std(B.sbp2 - B.sbp1, ddof=1) / math.sqrt(2), 'tc': np.std(np.log(B.tc2 / B.tc1), ddof=1) / math.sqrt(2), 'hdl': np.std(np.log(B.hdl2 / B.hdl1), ddof=1) / math.sqrt(2), 'scr': np.std(np.log(B.scr2 / B.scr1), ddof=1) / math.sqrt(2), 'h': np.nanstd(h2 - h1, ddof=1) / math.sqrt(2)}
            if model == 'uacr':
                sd['u'] = np.std(np.log(B.uacr2[U] / B.uacr1[U]), ddof=1) / math.sqrt(2)
        probs = np.zeros(nn)
        for r in range(R_SIM):
            cats = []
            for v in range(2):
                z = rng.standard_normal((6, nn))
                if mode == 'A':
                    sb = c_sbp + z[0] * sd['sbp']
                    t_ = c_tc * np.exp(z[1] * sd['tc'])
                    hd = c_hdl * np.exp(z[2] * sd['hdl'])
                    sc = c_scr * np.exp(z[3] * sd['scr'])
                    hv = c_h + z[4] * sd['h']
                    uu = c_u * np.exp(z[5] * sd['u']) if model == 'uacr' else None
                else:
                    lo, hi = cva_bracket
                    sig = lambda key: math.sqrt(math.log(1 + (CVA[key][0 if lo else 1] / 100) ** 2))
                    sb = c_sbp
                    t_ = np.round(c_tc * np.exp(z[1] * sig('tc')))
                    hd = np.round(c_hdl * np.exp(z[2] * sig('hdl')))
                    sc = 0.96 * np.round(c_raw * np.exp(z[3] * sig('scr_raw')), 1) - 0.184
                    hv = np.round(c_h * np.exp(z[4] * sig('hba1c')), 1)
                    uu = None
                dmv = np.where((fxx == 1) | (np.nan_to_num(hv, nan=0) >= 6.5), 1.0, 0.0)
                Iv = {'sbp': M.clamp(sb, 90, 200), 'tc': M.clamp(t_, 130, 320), 'hdl': M.clamp(hd, 20, 100), 'egfr': M.clamp(M.ckd_epi_2021(np.maximum(sc, 1e-06), Fx['age'], Fx['female']), 15, 140), 'dm': dmv, 'uacr': M.clamp(uu, 0.1, 25000) if uu is not None else None}
                cats.append(M.category(M.risk_from_lp(lp_of(Fx, Iv, model))))
            probs += cats[0] != cats[1]
        probs /= R_SIM
        mcse = math.sqrt(np.sum(probs * (1 - probs)) / R_SIM) / nn
        return (probs, mcse)
    pA, mA = simulate('base', np.arange(n), 'A')
    pAU, mAU = simulate('uacr', U, 'A')
    pBl, mBl = simulate('base', np.arange(n), 'B', (True, False))
    pBh, mBh = simulate('base', np.arange(n), 'B', (False, True))
    lab = sum((1 << GROUPS.index(g) for g in ['LIPID', 'KIDNEY', 'GLYC']))
    dsym = (D1B[:, lab] + D2B[:, lab]) / 2
    ratio_boot = lambda p: [float(np.sum(W[b] * p) / np.sum(W[b] * dsym)) if np.sum(W[b] * dsym) > 0 else np.nan for b in range(B_BOOT)]
    res['S12'] = {'A_base': {'theta_sim': float(pA.mean()), 'ratio_theta_over_sim': float(chB.mean() / pA.mean()), 'mcse': mA, 'ratio_ci': S.pct_ci([np.sum(W[b] * chB) / np.sum(W[b] * pA) for b in range(B_BOOT)])}, 'A_uacr': {'theta_sim': float(pAU.mean()), 'ratio_theta_over_sim': float(chU.mean() / pAU.mean()), 'mcse': mAU}, 'B_base_low_cva': {'theta_sim': float(pBl.mean()), 'ratio_to_Dsym_lab': float(pBl.mean() / dsym.mean()) if dsym.mean() > 0 else None, 'ratio_ci': S.pct_ci(ratio_boot(pBl)), 'mcse': mBl}, 'B_base_high_cva': {'theta_sim': float(pBh.mean()), 'ratio_to_Dsym_lab': float(pBh.mean() / dsym.mean()) if dsym.mean() > 0 else None, 'ratio_ci': S.pct_ci(ratio_boot(pBh)), 'mcse': mBh}, 'D_sym_lab_observed': float(dsym.mean()), 'B_uacr': 'not computed: no NHANES III urine albumin analytical CV (prefreeze resolution)', 'note': 'model-based; simulation SDs fixed at full-sample values (bootstrap varies people, not SDs)'}

    def interp_b(ci):
        if ci[1] < 0.5:
            return 'mainly within-person biological or pre-analytical'
        if ci[0] <= 1 <= ci[1]:
            return 'compatible with analytical imprecision alone'
        return 'partly analytical'
    wl, wh = (interp_b(res['S12']['B_base_low_cva']['ratio_ci']), interp_b(res['S12']['B_base_high_cva']['ratio_ci']))
    res['S12']['interpretive_b'] = wl if wl == wh else 'depends on the assumed analytical CV'
    Mn = (LB[:, 0] + LB[:, -1]) / 2
    Cm = M.category(M.risk_from_lp(Mn))
    res['S13'] = {'by_visit1_category': {str(k): S.prop_ci(chB[CB[:, 0] == k], cl[CB[:, 0] == k]) for k in range(4) if (CB[:, 0] == k).sum() > 0}, 'by_visitmean_category': {str(k): S.prop_ci(chB[Cm == k], cl[Cm == k]) for k in range(4) if (Cm == k).sum() > 0}, 'rtm_slope_on_lp1': float(np.polyfit(LB[:, 0] - LB[:, 0].mean(), dLB, 1)[0]), 'rtm_slope_on_mean': float(np.polyfit(Mn - Mn.mean(), dLB, 1)[0])}
    mid = (RB[:, 0] >= 5) & (RB[:, 0] < 10)
    res['S13']['ye_style_mid'] = {'n': int(mid.sum()), 'changing_category': float(chB[mid].mean()), 'crossing_5_or_10': float((((RB[:, 0] >= 5) != (RB[:, -1] >= 5)) | ((RB[:, 0] >= 10) != (RB[:, -1] >= 10)))[mid].mean())}
    subs = {}
    tert = pd.qcut(B.rr_days, 3, labels=False, duplicates='drop').to_numpy()
    for name, lab_arr in [('sex', np.where(B.female == 1, 'female', 'male')), ('age_band', pd.cut(B.age, [29, 39, 59, 79], labels=['30-39', '40-59', '60-79']).astype(str)), ('race', B.dmarethn.map({1: 'NHW', 2: 'NHB', 3: 'MexAm', 4: 'other'}).to_numpy()), ('bptreat', B.bptreat.astype(int).astype(str).to_numpy()), ('interval_tertile', tert.astype(str)), ('phase', B.phase.astype(str).to_numpy()), ('same_examiner', B.same_examiner.astype(str).to_numpy())]:
        subs[name] = {}
        for lv in np.unique(lab_arr):
            m = lab_arr == lv
            ci = S.prop_ci(chB[m], cl[m])
            ci['nchs_flags'] = S.nchs_flag(ci['events'], ci['n'], ci['lower'], ci['upper'], ci['df'])
            subs[name][str(lv)] = ci
    res['S13']['subgroups_descriptive'] = subs
    P1 = pd.read_csv(ROOT / 'studies/uacr_prevent/final_and_external_v1/primary_participant_results.csv')
    demo, _ = pyreadstat.read_xport(ROOT / 'studies/uacr_prevent/data/DEMO_F.xpt', usecols=['SEQN', 'WTMEC2YR'])
    demo['SEQN'] = demo.SEQN.astype(int)
    if 'WTMEC2YR' not in P1.columns:
        P1 = P1.merge(demo, on='SEQN', how='left')
    FP = {'female': P1.sex.to_numpy(float), 'age': P1.age.to_numpy(float), 'smoking': P1.smoking.to_numpy(float), 'bptreat': P1.bptreat.to_numpy(float), 'statin': P1.statin.to_numpy(float)}
    IP = {'sbp': M.clamp(P1.sbp, 90, 200), 'tc': M.clamp(P1.tc, 130, 320), 'hdl': M.clamp(P1.hdl, 20, 100), 'egfr': M.clamp(P1.egfr, 15, 140), 'dm': P1.dm.to_numpy(float), 'uacr': None}
    LP1p = lp_of(FP, IP, 'base')
    ph = predict(k3B['clin']['_model'], LP1p)
    wp = P1.WTMEC2YR.to_numpy(float)
    res['E14'] = {'n': int(len(P1)), 'transported_theta': float(np.sum(wp * ph) / np.sum(wp)), 'note': 'exploratory; assumes 1988-94 within-person variability applies to 2009-10; Paper 1 SBP = mean of available readings (mismatch)'}
    ex = B[['SEQN', 'cluster', 'stratum', 'psu', 'female', 'age', 'smoking', 'bptreat', 'statin', 'dmarethn', 'in_uacr']].copy()
    for v, I in ((1, I1), (2, I2)):
        for k in ['sbp', 'tc', 'hdl', 'egfr', 'dm', 'uacr']:
            ex[f'{k}{v}_used'] = I[k]
        for k in ['sbp', 'tc', 'hdl', 'egfr']:
            ex[f'{k}{v}_clamped'] = B[f'{k}{v}'].to_numpy() != I[k]
    for m in range(1 << KB):
        ex[f'riskB_m{m}'] = RB[:, m]
        ex[f'lpB_m{m}'] = LB[:, m]
    ux = np.full((n, 1 << KU), np.nan)
    ux[U] = RU
    for m in range(1 << KU):
        ex[f'riskU_m{m}'] = ux[:, m]
    for k, g in enumerate(GROUPS[:4]):
        ex[f'psiB_{g}'] = psiB[:, k]
        ex[f'fB_{g}'] = fB[:, k]
    pu = np.full((n, KU), np.nan)
    pu[U] = psiU
    for k, g in enumerate(GROUPS):
        ex[f'psiU_{g}'] = pu[:, k]
    ex['changeB'] = chB
    ex['cvd_cross75'] = cr75
    ex['riskT1'] = RT1
    ex['riskT2'] = RT2
    s11c = np.full((n, 6), np.nan)
    s11c[a40] = np.column_stack([p1, pa, pb, xa75, xb75, cb_])
    for j, nm in enumerate(['pce_v1', 'pce_sbp_only', 'pce_full', 'pce_x75_sbp', 'pce_x75_full', 'pce_cat_full']):
        ex[nm] = s11c[:, j]
    ex.to_csv(OUT / 'stage3_participants.csv', index=False)
    coef = {m: {c: np.asarray(TABLE[c], float).tolist() for c in ['female_ascvd', 'male_ascvd', 'female_total_cvd', 'male_total_cvd']} for m, TABLE in M.TABLE.items()}
    (OUT / 'preventr_coefficients.json').write_text(json.dumps(coef) + '\n', encoding='utf-8')
    res['environment'] = {'python': platform.python_version(), **{p: version(p) for p in ['numpy', 'pandas', 'scipy', 'pyreadstat', 'rdata']}}
    res['completed_utc'] = datetime.now(timezone.utc).isoformat()
    save('stage3_results.json', res)
    pri = res['primary']
    print(json.dumps({'theta': [pri['events'], pri['n'], round(100 * pri['estimate'], 2), round(100 * pri['lower'], 2), round(100 * pri['upper'], 2), pri['precision_label']], 'up_down': [pri['up'], pri['down']], 'theta_U': round(100 * res['K1']['theta_U']['estimate'], 2), 'delta_U_pp': [round(100 * res['K1']['delta_U'][k], 2) for k in ('estimate', 'lower', 'upper')], 'Psi_base': {g: round(v['pp'], 2) for g, v in res['K2_base']['Psi'].items()}, 'largest': res['K2_base']['largest_contributor'], 'theta_std': round(100 * res['K4_base']['theta_std'], 2), 'rcv_risk': res['S6_base']['rcv_risk']}, indent=1))
if __name__ == '__main__':
    main()
