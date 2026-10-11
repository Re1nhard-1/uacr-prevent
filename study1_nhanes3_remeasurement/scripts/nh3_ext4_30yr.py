from pathlib import Path
import hashlib
import json
import math
import numpy as np
import pandas as pd
from scipy.special import expit, logit
import nh3_build
import nh3_models as M
import nh3_stats as S
import nh3_testretest_stage12 as s12
from nh3_stage3 import fixed, lp_of, visit_inputs
from nh3_ext3_statin import sampson, tiers, risk30, check30, prop, REC3
from nh3_ext2021_burden import xpt, Svy, D as D21
ROOT = Path(__file__).resolve().parents[1]
ST = ROOT / 'studies/nhanes3_testretest'
PROTO = ST / 'extensions4_protocol_20261010.md'
OUT = ST / 'analysis_v1_posthoc/ext4_30yr_results.json'
SEED, B_BOOT = (20261009, 4000)
LOOK30 = [5, 7, 8, 9, 10, 11, 12, 13, 15, 20]

def curve(dist, y, W, grid_dist, cut_logit=None):
    knots = np.percentile(dist, [10, 50, 90])
    X = lambda dv: np.column_stack([np.ones(len(dv)), S.rcs_basis(dv, knots)])
    beta = S.logit_fit(X(dist), y)
    boots = np.array([S.logit_fit(X(dist), y, W[b]) for b in range(W.shape[0])])
    p, pb = (expit(X(grid_dist) @ beta), expit(X(grid_dist) @ boots.T))
    fine = np.linspace(0, 2.0, 2001)
    pf = expit(X(fine) @ beta)
    zone = {}
    for lev in (0.1, 0.25):
        ok = np.where(pf >= lev)[0]
        if len(ok) and ok[0] == 0:
            zone[str(lev)] = float(fine[np.argmax(pf < lev) - 1]) if np.any(pf < lev) else float(fine[-1])
    p0 = expit(X(np.array([0.0])) @ boots.T)[0]
    return {'knots': knots, 'beta': beta, 'boots': boots, 'X': X, 'p': p, 'lower': np.percentile(pb, 2.5, axis=1), 'upper': np.percentile(pb, 97.5, axis=1), 'zone_dist': zone, 'pi_at_threshold': float(expit(X(np.array([0.0])) @ beta)[0]), 'pi_at_threshold_ci': S.pct_ci(p0)}

def main():
    ptxt = PROTO.read_bytes()
    out = {'label': 'E10: 30-year-risk statin pathway, remeasurement stability (Study 1) and US reach (NHANES 2021-2023, exploratory)', 'protocol': str(PROTO.relative_to(ROOT)), 'protocol_sha256': hashlib.sha256(ptxt).hexdigest(), 'protocol_sha256_before_addendum': hashlib.sha256(ptxt.split(b'\n## Addendum')[0]).hexdigest()}
    E, B, T = nh3_build.build()
    F = fixed(B)
    cl = B.cluster.to_numpy()
    age = B.age.to_numpy(float)
    keep = set(B.SEQN)
    I1, I2 = (visit_inputs(B, 1), visit_inputs(B, 2))
    r1, r2 = (M.risk_from_lp(lp_of(F, I1, 'base')), M.risk_from_lp(lp_of(F, I2, 'base')))
    Pp = pd.read_csv(ST / 'analysis_v1/stage3/stage3_participants.csv').set_index('SEQN').loc[B.SEQN]
    assert np.max(np.abs(r1 - Pp.riskB_m0.to_numpy())) < 1e-09 and np.max(np.abs(r2 - Pp.riskB_m15.to_numpy())) < 1e-09, 'visit risks not reproduced'
    q1, q2 = (risk30(F, I1), risk30(F, I2))
    statin = B.statin.to_numpy() == 1
    A1 = (age >= 30) & (age <= 59) & ~statin
    idx = np.where(A1)[0]
    chk = max(check30(F, I1, q1, idx), check30(F, I2, q2, idx))
    assert chk < 1e-09, ('30-year check', chk)
    flip = ((q1 >= 10) != (q2 >= 10)).astype(float)
    A = {'n_A1': int(len(idx)), 'statin_users_30_59_excluded': int(((age >= 30) & (age <= 59) & statin).sum()), 'risk30_check_max_abs_diff_pp': chk, 'risk30_v1_median': float(np.median(q1[idx])), 'risk30_v1_iqr': [float(np.percentile(q1[idx], 25)), float(np.percentile(q1[idx], 75))], 'ge10_v1': int((q1[idx] >= 10).sum()), 'ge10_v2': int((q2[idx] >= 10).sum()), 'up': int(((q1 < 10) & (q2 >= 10))[idx].sum()), 'down': int(((q1 >= 10) & (q2 < 10))[idx].sum()), 'theta30': prop(flip[idx], cl[idx])}
    exm = s12.read('exam', {'PEPTECH': (1498, 1502)}, keep=keep)
    phys = np.array([exm[s]['PEPTECH'].strip() for s in B.SEQN])
    phys_single = np.where(phys == '88888', np.array([f'missing_{s}' for s in B.SEQN]), phys)
    for nm, g in [('theta30_visit1_physician', phys), ('theta30_visit1_physician_missing_singletons', phys_single)]:
        r_ = S.prop_ci(flip[idx], g[idx])
        A[nm] = {k: r_[k] for k in ['lower', 'upper', 'deff', 'G']}
    W = S.boot_weights(cl, B_BOOT, SEED)[:, idx]
    dA = np.abs(logit(q1[idx] / 100) - logit(0.1))
    grid = np.round(np.arange(2.0, 40.0001, 0.05), 2)
    dg = np.abs(logit(grid / 100) - logit(0.1))
    C = curve(dA, flip[idx], W, dg)
    look = []
    for r in LOOK30:
        i = int(np.argmin(np.abs(grid - r)))
        look.append({'risk30_pct': r, 'p': float(C['p'][i]), 'lower': float(C['lower'][i]), 'upper': float(C['upper'][i])})
    zone30 = {k: [float(100 * expit(logit(0.1) - v)), float(100 * expit(logit(0.1) + v))] for k, v in C['zone_dist'].items()}
    A['pi30'] = {'knots': C['knots'].tolist(), 'beta': C['beta'].tolist(), 'pi_at_threshold': C['pi_at_threshold'], 'pi_at_threshold_ci': C['pi_at_threshold_ci'], 'zones_risk30_pct': zone30, 'lookup': look, 'events': int(flip[idx].sum()), 'in_zone_A1': {k: int(((q1 >= a) & (q1 <= b))[idx].sum()) for k, (a, b) in zone30.items()}, 'note': 'chance that remeasurement changes the 30-year >=10% classification, by first-visit 30-year risk; conditions on visit 1 (includes regression to the mean)', 'grid': {'risk30_pct': grid.tolist(), 'p': C['p'].round(5).tolist(), 'lower': C['lower'].round(5).tolist(), 'upper': C['upper'].round(5).tolist()}, 'observed_bins': [{'risk30_lo': lo, 'risk30_hi': hi, 'n': int(((q1 >= lo) & (q1 < hi))[idx].sum()), 'changes': int((flip.astype(bool) & (q1 >= lo) & (q1 < hi))[idx].sum())} for lo, hi in [(0, 5), (5, 7), (7, 8.5), (8.5, 10), (10, 11.5), (11.5, 13), (13, 16), (16, 25), (25, 100)]]}
    l1 = s12.read('lab', {'TGP': (1606, 1609)}, keep=keep)
    l2 = s12.read('labse', {'TGR': (358, 361)}, keep=keep)
    g_ = lambda d, s, k: s12.num(d[s][k]) if s in d else None
    tg1 = np.array([g_(l1, s, 'TGP') for s in B.SEQN], dtype=float)
    tg2 = np.array([g_(l2, s, 'TGR') for s in B.SEQN], dtype=float)
    tc1, tc2, hd1, hd2 = (B[c].to_numpy(float) for c in ['tc1', 'tc2', 'hdl1', 'hdl2'])
    ls1, ls2 = (sampson(tc1, hd1, tg1), sampson(tc2, hd2, tg2))
    tg_ok = np.isfinite(tg1) & np.isfinite(tg2) & (tg1 <= 800) & (tg2 <= 800)
    jl = np.where(A1 & tg_ok)[0]
    fl = ((ls1 >= 160) != (ls2 >= 160)).astype(float)
    A['ldl160'] = {'n': int(len(jl)), 'ge160_v1': int((ls1[jl] >= 160).sum()), 'flip': prop(fl[jl], cl[jl])}
    zone_ldl = None
    if fl[jl].sum() >= 15:
        Wl = S.boot_weights(cl, B_BOOT, SEED)[:, jl]
        dl = np.abs(np.log(ls1[jl]) - np.log(160.0))
        gl = np.arange(100.0, 240.01, 1.0)
        dgl = np.abs(np.log(gl) - np.log(160.0))
        CL = curve(dl, fl[jl], Wl, dgl)
        zone_ldl = {k: [float(160 * math.exp(-v)), float(160 * math.exp(v))] for k, v in CL['zone_dist'].items()}
        lk = [{'ldl': v, 'p': float(CL['p'][int(v - 100)]), 'lower': float(CL['lower'][int(v - 100)]), 'upper': float(CL['upper'][int(v - 100)])} for v in [130, 140, 150, 155, 160, 165, 170, 180]]
        A['ldl160']['curve'] = {'knots': CL['knots'].tolist(), 'beta': CL['beta'].tolist(), 'pi_at_threshold': CL['pi_at_threshold'], 'pi_at_threshold_ci': CL['pi_at_threshold_ci'], 'zones_ldl_mgdl': zone_ldl, 'lookup': lk, 'grid': {'ldl_mgdl': gl.tolist(), 'p': CL['p'].round(5).tolist(), 'lower': CL['lower'].round(5).tolist(), 'upper': CL['upper'].round(5).tolist()}, 'observed_bins': [{'ldl_lo': lo, 'ldl_hi': hi, 'n': int(((ls1 >= lo) & (ls1 < hi))[jl].sum()), 'changes': int((fl.astype(bool) & (ls1 >= lo) & (ls1 < hi))[jl].sum())} for lo, hi in [(0, 130), (130, 145), (145, 155), (155, 165), (165, 175), (175, 190), (190, 400)]]}
    else:
        A['ldl160']['curve'] = 'not estimated (<15 crossings; protocol)'
    e1, e2 = (B.egfr1.to_numpy(float), B.egfr2.to_numpy(float))
    dm1, dm2 = (I1['dm'], I2['dm'])
    t1, t2 = (tiers(ls1, r1, dm1, e1, age, 'v2', q1), tiers(ls2, r2, dm2, e2, age, 'v2', q2))
    rr = lambda t: np.isin(t, REC3['v2'])
    dom = A1 & tg_ok & (t1 != 'I') & (ls1 >= 70) & (ls1 < 160) & (r1 < 3)
    jd = np.where(dom)[0]
    el = q1 >= 10
    A['pathway_domain_v1'] = {'n': int(len(jd)), 'eligible_30yr_v1': int(el[jd].sum()), 'risk30_flag_differs_v2': prop(flip[jd], cl[jd]), 'cor1_or_2a_differs_v2': prop((rr(t1) != rr(t2))[jd].astype(float), cl[jd])}
    dom2 = A1 & tg_ok & (t2 != 'I') & (ls2 >= 70) & (ls2 < 160) & (r2 < 3)
    elig2 = dom2 & (q2 >= 10)
    for nm, sub in [('among_eligible_v1', dom & el), ('among_not_eligible_v1', dom & ~el)]:
        js = np.where(sub)[0]
        A['pathway_domain_v1'][nm] = {'n': int(len(js)), 'risk30_flag_differs_v2': prop(flip[js], cl[js]), 'cor1_or_2a_differs_v2': prop((rr(t1) != rr(t2))[js].astype(float), cl[js]), 'tier_v2_counts': {k: int((t2[js] == k).sum()) for k in 'NLRCSHI'}}
        if nm == 'among_eligible_v1':
            A['pathway_domain_v1'][nm]['not_pathway_eligible_v2'] = prop((~elig2)[js].astype(float), cl[js])
            A['pathway_domain_v1'][nm]['not_pathway_eligible_v2_reasons'] = {'risk30_below_10': int((q2 < 10)[js].sum()), 'risk30_ge10_but_left_domain': {'indication': int(((q2 >= 10) & (t2 == 'I'))[js].sum()), 'ldl_ge160': int(((q2 >= 10) & (t2 != 'I') & (ls2 >= 160))[js].sum()), 'ldl_lt70': int(((q2 >= 10) & (t2 != 'I') & (ls2 < 70))[js].sum()), 'tenyear_ge3': int(((q2 >= 10) & (t2 != 'I') & (ls2 >= 70) & (ls2 < 160) & (r2 >= 3))[js].sum())}}
        else:
            A['pathway_domain_v1'][nm]['became_pathway_eligible_v2'] = prop(elig2[js].astype(float), cl[js])
    A['pathway_domain_v1']['note'] = 'conditions on visit-1 status, so includes regression to the mean; risk30_flag_differs_v2 = the 30-year >=10% criterion only; not_pathway_eligible_v2 = out of the pathway domain or 30-year <10% at visit 2; cor1_or_2a = any rule-v2 recommendation (I, R, C, S, H)'
    geo = lambda a, b: np.sqrt(a * b)
    h1, h2 = (B.hba1c1.to_numpy(float), B.hba1c2.to_numpy(float))
    hbar = np.where(np.isfinite(h1) & np.isfinite(h2), (h1 + h2) / 2, np.where(np.isfinite(h1), h1, h2))
    dmA = np.where((B.fixed_dm.to_numpy() == 1) | (np.nan_to_num(hbar) >= 6.5), 1.0, 0.0)
    egA = np.asarray(M.ckd_epi_2021(geo(B.scr1.to_numpy(float), B.scr2.to_numpy(float)), B.age, B.female), float)
    IA = {'sbp': M.clamp((B.sbp1.to_numpy(float) + B.sbp2.to_numpy(float)) / 2, 90, 200), 'tc': M.clamp(geo(tc1, tc2), 130, 320), 'hdl': M.clamp(geo(hd1, hd2), 20, 100), 'egfr': M.clamp(egA, 15, 140), 'dm': dmA, 'uacr': None}
    qA = risk30(F, IA)
    assert check30(F, IA, qA, idx) < 1e-09
    dv1, dv2 = ((q1 >= 10) != (qA >= 10), (q2 >= 10) != (qA >= 10))
    A['vs_two_visit_mean'] = {'mean_single_visit_differs': S.paired_ci((dv1 + 0.0 + dv2)[idx] / 2, cl[idx]), 'ge10_two_visit_mean': int((qA[idx] >= 10).sum())}
    out['part_A'] = A
    demo = xpt('DEMO_L', ['RIDSTATR', 'RIAGENDR', 'RIDAGEYR', 'RIDEXPRG', 'SDMVSTRA', 'SDMVPSU'])
    d = demo
    for name, cols in [('MCQ_L', ['MCQ160B', 'MCQ160C', 'MCQ160D', 'MCQ160E', 'MCQ160F']), ('BPXO_L', ['BPXOSY1', 'BPXOSY2', 'BPXOSY3']), ('TCHOL_L', ['LBXTC']), ('HDL_L', ['LBDHDD']), ('BIOPRO_L', ['LBXSCR']), ('GHB_L', ['LBXGH']), ('BPQ_L', ['BPQ020', 'BPQ150', 'BPQ101D']), ('DIQ_L', ['DIQ010', 'DIQ050', 'DIQ070']), ('SMQ_L', ['SMQ020', 'SMQ040'])]:
        d = d.merge(xpt(name, cols), on='SEQN', how='left')
    import pyreadstat
    tr, _ = pyreadstat.read_xport(str(D21.parent / 'raw_trigly/TRIGLY_L.xpt'), usecols=['SEQN', 'WTSAF2YR', 'LBDLDLN', 'LBXTLG'], encoding='latin1')
    d = d.merge(tr, on='SEQN', how='left')
    flow = [('DEMO_L records', len(d))]
    d = d[(d.RIDSTATR == 2) & (d.WTSAF2YR > 0)]
    flow.append(('Examined, fasting subsample with positive weight', len(d)))
    d = d[(d.RIDAGEYR >= 30) & (d.RIDAGEYR <= 59)]
    flow.append(('Age 30-59 years', len(d)))
    d = d[(d[['MCQ160B', 'MCQ160C', 'MCQ160D', 'MCQ160E', 'MCQ160F']] == 2).all(axis=1)]
    flow.append(('Explicitly no CHF, CHD, angina, MI or stroke', len(d)))
    d = d[d.RIDEXPRG != 1]
    flow.append(('Not pregnant', len(d)))
    d = d[d.BPQ101D.isin([2])].copy()
    flow.append(('Not taking cholesterol medication (BPQ101D = 2)', len(d)))
    w_elig = float(d.WTSAF2YR.sum())
    d['sbp'] = d[['BPXOSY2', 'BPXOSY3']].mean(axis=1, skipna=False).fillna(d[['BPXOSY1', 'BPXOSY2', 'BPXOSY3']].mean(axis=1))
    qpos = (d.DIQ010 == 1) | (d.DIQ050 == 1) | (d.DIQ070 == 1)
    d['dm'] = np.where(qpos | (d.LBXGH >= 6.5), 1.0, np.where(d.LBXGH.notna() & d.DIQ010.isin([2, 3]), 0.0, np.nan))
    d['smoking'] = np.where((d.SMQ020 == 1) & d.SMQ040.isin([1, 2]), 1.0, np.where((d.SMQ020 == 2) | (d.SMQ020 == 1) & (d.SMQ040 == 3), 0.0, np.nan))
    d['bptreat'] = np.where(d.BPQ150 == 1, 1.0, np.where((d.BPQ150 == 2) | (d.BPQ020 == 2), 0.0, np.nan))
    d['female'] = (d.RIAGENDR == 2).astype(float)
    need = ['sbp', 'LBXTC', 'LBDHDD', 'LBXSCR', 'dm', 'smoking', 'bptreat', 'LBDLDLN']
    d = d[d[need].notna().all(axis=1)].copy()
    flow.append(('Complete base PREVENT inputs and LDL-C (NIH equation 2)', len(d)))
    eg = np.asarray(M.ckd_epi_2021(d.LBXSCR.to_numpy(float), d.RIDAGEYR.to_numpy(float), d.female.to_numpy(float)), float)
    Fb = {'female': d.female.to_numpy(float), 'age': d.RIDAGEYR.to_numpy(float), 'smoking': d.smoking.to_numpy(float), 'bptreat': d.bptreat.to_numpy(float), 'statin': np.zeros(len(d))}
    Ib = {'tc': M.clamp(d.LBXTC, 130, 320), 'hdl': M.clamp(d.LBDHDD, 20, 100), 'sbp': M.clamp(d.sbp, 90, 200), 'dm': d.dm.to_numpy(float), 'egfr': M.clamp(eg, 15, 140), 'uacr': None}
    rb = M.risk_from_lp(lp_of(Fb, Ib, 'base'))
    qb = risk30(Fb, Ib)
    chkB = check30(Fb, Ib, qb, np.arange(len(d)))
    assert chkB < 1e-09, ('30-year check B', chkB)
    ldl, ageb, dmb = (d.LBDLDLN.to_numpy(float), d.RIDAGEYR.to_numpy(float), d.dm.to_numpy(float))
    a4075 = (ageb >= 40) & (ageb <= 75)
    ind = (ldl >= 190) | a4075 & (dmb == 1) | a4075 & (eg < 60) & (ldl >= 70)
    B1 = ~ind & (ldl >= 70) & (ldl < 160) & (rb < 3)
    sv = Svy(d.WTSAF2YR, d.SDMVSTRA, d.SDMVPSU)
    scale = w_elig / float(d.WTSAF2YR.sum())
    el = (qb >= 10).astype(float)
    Bres = {'flow': flow, 'df': sv.df, 'weight': 'WTSAF2YR (fasting subsample)', 'completeness_scale': scale, 'risk30_check_max_abs_diff_pp': chkB, 'population_millions': float(d.WTSAF2YR.sum() / 1000000.0), 'B1_n': int(B1.sum()), 'B1_share_of_population': sv.mean(B1.astype(float)), 'B1_millions': sv.total(B1.astype(float)), 'B1_pathway_eligible': {**sv.mean(el, B1), 'millions': sv.total(el, B1)}, 'B1_risk30_median_unweighted': float(np.median(qb[B1]))}
    inz = lambda lo, hi: ((qb >= lo) & (qb <= hi)).astype(float)
    fine = np.linspace(0, 2.0, 2001)
    Xf = C['X'](fine)
    for k, (lo, hi) in zone30.items():
        m = inz(lo, hi)
        Bres[f'B1_in_pi30_zone_{k}'] = {**sv.mean(m, B1), 'millions': sv.total(m, B1)}
        Bres[f'eligible_in_pi30_zone_{k}'] = {**sv.mean(m, B1 & (el == 1)), 'millions': sv.total(m, B1 & (el == 1))}
        lev = float(k)
        dom_e = B1 & (el == 1)
        we = d.WTSAF2YR.to_numpy(float) * dom_e
        we = we / we.sum()
        shares = []
        for b in C['boots']:
            pf = expit(Xf @ b)
            if pf[0] < lev:
                shares.append(0.0)
                continue
            dz = fine[np.argmax(pf < lev) - 1] if np.any(pf < lev) else fine[-1]
            shares.append(float(we @ inz(100 * expit(logit(0.1) - dz), 100 * expit(logit(0.1) + dz))))
        shares = np.array(shares)
        est = Bres[f'eligible_in_pi30_zone_{k}']['estimate']
        sdz = math.sqrt(np.var(shares, ddof=1) + Bres[f'eligible_in_pi30_zone_{k}']['se'] ** 2)
        Bres[f'eligible_in_pi30_zone_{k}']['with_curve_uncertainty'] = {'lower': max(0.0, est - 1.96 * sdz), 'upper': min(1.0, est + 1.96 * sdz), 'replicate_percentile': S.pct_ci(shares), 'note': 'zone limits recomputed in each of the 4,000 curve replicates; replicate variance plus survey variance'}
    Bres['zone_intervals_note'] = 'zone limits fixed at the Part A point estimates; plain intervals reflect survey sampling only (as E2); with_curve_uncertainty adds the curve replicate variance. The zones are where remeasurement changes the 30-year >=10% classification; eligibility can also change through the 10-year 3% cut, LDL-C 70/160 and the indications'
    bins = np.arange(0, 40.01, 0.5)
    wfull = d.WTSAF2YR.to_numpy(float)
    Bres['B1_risk30_histogram'] = {'bin_edges': bins.tolist(), 'weighted_millions_all_B1': [float(wfull[B1 & (qb >= a) & (qb < b)].sum() / 1000000.0) for a, b in zip(bins[:-1], bins[1:])], 'over_40_millions': float(wfull[B1 & (qb >= 40)].sum() / 1000000.0)}
    dB = np.abs(logit(qb / 100) - logit(0.1))
    XB = C['X'](dB)
    piB = expit(XB @ C['beta'])
    th = sv.mean(piB, B1, logit_ci=False)
    wB = d.WTSAF2YR.to_numpy(float) * B1
    wB = wB / wB.sum()
    cb = np.array([wB @ expit(XB @ b) for b in C['boots']])
    sd = math.sqrt(np.var(cb, ddof=1) + th['se'] ** 2)
    Bres['B1_transported_expected_flip'] = {'estimate': th['estimate'], 'lower': th['estimate'] - 1.96 * sd, 'upper': th['estimate'] + 1.96 * sd, 'curve_percentile': S.pct_ci(cb), 'survey_se': th['se'], 'note': 'curve bootstrap and survey variance combined'}
    if zone_ldl and '0.25' in zone_ldl:
        lo, hi = zone_ldl['0.25']
        m = ((ldl >= lo) & (ldl <= hi)).astype(float)
        Bres['B1_ldl_in_ldl160_zone_0.25'] = {**sv.mean(m, B1), 'millions': sv.total(m, B1)}
    for k, v in list(Bres.items()):
        if isinstance(v, dict) and 'millions' in v and isinstance(v['millions'], dict):
            v['millions_scaled'] = {kk: vv * scale for kk, vv in v['millions'].items()}
    Bres['B1_millions_scaled'] = {kk: vv * scale for kk, vv in Bres['B1_millions'].items()}
    Bres['assumptions'] = 'PH-P, exploratory transport. NHANES III (1988-1994) within-person variation and distance-only dependence transported; mercury (1988-1994) versus oscillometric (2021-2023) BP; fasting LDL-C by NIH equation 2; statin proxy BPQ101D (any cholesterol medicine); risk enhancers, CAC, HIV and LDL-C goals not modelled; projection, not an observed remeasurement'
    Bres['millions_note'] = 'millions are complete-case weighted totals; *_scaled multiply by the eligible/complete fasting-weight ratio, assuming missingness unrelated to the outcome'
    out['part_B'] = Bres
    OUT.write_text(json.dumps(out, indent=1, default=lambda o: o.item() if hasattr(o, 'item') else o.tolist() if hasattr(o, 'tolist') else str(o)) + '\n', encoding='utf-8')
    print('saved', OUT)
if __name__ == '__main__':
    main()
