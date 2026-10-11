"""Post hoc, model-based (2026-10-10): how much would averaging repeat measurements, or repeating them only near category"""
from pathlib import Path
import json
import math
import numpy as np
from scipy.special import expit, logit
import nh3_build
import nh3_models as M
import nh3_stats as S
from nh3_stage3 import fixed, lp_of, SEED, B_BOOT, CUTS
ROOT = Path(__file__).resolve().parents[1]
ST = ROOT / 'studies/nhanes3_testretest'
R3 = json.loads((ST / 'analysis_v1/stage3/stage3_results.json').read_text(encoding='utf-8'))
OUT = ST / 'analysis_v1_posthoc/repeat_strategies.json'
R_SIM = 1000

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

def main():
    E, B, T = nh3_build.build()
    n = len(B)
    cl = B.cluster.to_numpy()
    F = fixed(B)
    assert n == R3['primary']['n']
    sbp1, sbp2 = (B.sbp1.to_numpy(float), B.sbp2.to_numpy(float))
    tc1, tc2, hd1, hd2 = (B[c].to_numpy(float) for c in ['tc1', 'tc2', 'hdl1', 'hdl2'])
    sc1, sc2 = (B.scr1.to_numpy(float), B.scr2.to_numpy(float))
    h1, h2 = (B.hba1c1.to_numpy(float), B.hba1c2.to_numpy(float))
    fx = B.fixed_dm.to_numpy(float)
    geo = lambda a, b: np.sqrt(a * b)
    c = {'sbp': (sbp1 + sbp2) / 2, 'tc': geo(tc1, tc2), 'hdl': geo(hd1, hd2), 'scr': geo(sc1, sc2), 'h': np.where(np.isfinite(h1) & np.isfinite(h2), (h1 + h2) / 2, np.where(np.isfinite(h1), h1, h2))}
    sd = {'sbp': np.std(sbp2 - sbp1, ddof=1) / math.sqrt(2), 'tc': np.std(np.log(tc2 / tc1), ddof=1) / math.sqrt(2), 'hdl': np.std(np.log(hd2 / hd1), ddof=1) / math.sqrt(2), 'scr': np.std(np.log(sc2 / sc1), ddof=1) / math.sqrt(2), 'h': np.nanstd(h2 - h1, ddof=1) / math.sqrt(2)}
    LAB = ['tc', 'hdl', 'scr', 'h']
    Z25, Z10 = (zones('0.25'), zones('0.1'))
    inzone = lambda r, Z: np.any([(r >= a) & (r <= b) for a, b in Z], axis=0)

    def category(zbar):
        sb = c['sbp'] + zbar['sbp'] * sd['sbp']
        t_ = c['tc'] * np.exp(zbar['tc'] * sd['tc'])
        hd = c['hdl'] * np.exp(zbar['hdl'] * sd['hdl'])
        sc = c['scr'] * np.exp(zbar['scr'] * sd['scr'])
        hv = c['h'] + zbar['h'] * sd['h']
        dm = np.where((fx == 1) | (np.nan_to_num(hv, nan=0) >= 6.5), 1.0, 0.0)
        I = {'sbp': M.clamp(sb, 90, 200), 'tc': M.clamp(t_, 130, 320), 'hdl': M.clamp(hd, 20, 100), 'egfr': M.clamp(M.ckd_epi_2021(np.maximum(sc, 1e-06), F['age'], F['female']), 15, 140), 'dm': dm, 'uacr': None}
        r = M.risk_from_lp(lp_of(F, I, 'base'))
        return (M.category(r), r)

    def assess(rng, k_bp, k_lab, target=None, target_inputs=None):
        if target is None:
            zbar = {'sbp': rng.standard_normal((k_bp, n)).mean(0)}
            for k in LAB:
                zbar[k] = rng.standard_normal((k_lab, n)).mean(0)
            return (category(zbar)[0], np.zeros(n, bool))
        z1 = {k: rng.standard_normal(n) for k in ['sbp'] + LAB}
        cat1, r1 = category(z1)
        rep = inzone(r1, target)
        z2 = dict(z1)
        for k in target_inputs:
            z2[k] = np.where(rep, (z1[k] + rng.standard_normal(n)) / 2, z1[k])
        return (category(z2)[0], rep)
    strategies = [('single', 'One measurement of each input (model A)', dict(k_bp=1, k_lab=1)), ('bp2', 'BP averaged over 2 days; laboratory inputs once', dict(k_bp=2, k_lab=1)), ('bp3', 'BP averaged over 3 days; laboratory inputs once', dict(k_bp=3, k_lab=1)), ('lab2', 'BP once; laboratory inputs averaged over 2 draws', dict(k_bp=1, k_lab=2)), ('all2', 'Everything measured on 2 days and averaged', dict(k_bp=2, k_lab=2)), ('all3', 'Everything measured on 3 days and averaged', dict(k_bp=3, k_lab=3)), ('t25_all', 'Second day of all inputs only if first risk in the >=25% zone', dict(k_bp=1, k_lab=1, target=Z25, target_inputs=['sbp'] + LAB)), ('t10_all', 'Second day of all inputs only if first risk in the >=10% zone', dict(k_bp=1, k_lab=1, target=Z10, target_inputs=['sbp'] + LAB)), ('t10_bp', 'Second BP day only if first risk in the >=10% zone', dict(k_bp=1, k_lab=1, target=Z10, target_inputs=['sbp']))]
    res = {'label': __doc__.split('\n')[0], 'sd_used': sd, 'zones_25': Z25, 'zones_10': Z10, 'R_SIM': R_SIM, 'strategies': {}}
    W = S.boot_weights(cl, B_BOOT, SEED)
    probs = {}
    for i, (key, desc, kw) in enumerate(strategies):
        rng = np.random.default_rng(SEED + 100 + i)
        pr = np.zeros(n)
        wl = np.zeros(n)
        for _ in range(R_SIM):
            ca, ra = assess(rng, **kw)
            cb, rb = assess(rng, **kw)
            pr += ca != cb
            wl += (ra.astype(float) + rb.astype(float)) / 2
        pr /= R_SIM
        wl /= R_SIM
        probs[key] = pr
        res['strategies'][key] = {'description': desc, 'theta_sim': float(pr.mean()), 'mcse': float(math.sqrt(np.sum(pr * (1 - pr)) / R_SIM) / n), 'repeat_workload': float(wl.mean()) if 'target' in kw else 0.0 if kw['k_bp'] == kw['k_lab'] == 1 else 1.0, 'extra_bp_days_per_person': float(wl.mean()) if 'target' in kw else float(kw['k_bp'] - 1), 'extra_lab_draws_per_person': (float(wl.mean()) if kw.get('target_inputs') and 'tc' in kw['target_inputs'] else 0.0) if 'target' in kw else float(kw['k_lab'] - 1)}
        print(f"{key:8s} theta_sim {100 * pr.mean():5.2f}%  workload {res['strategies'][key]['repeat_workload']:.3f}", flush=True)
    ref = probs['single']
    for key, pr in probs.items():
        rel = [float(np.sum(W[b] * pr) / np.sum(W[b] * ref)) for b in range(B_BOOT)]
        res['strategies'][key]['relative_to_single'] = float(pr.mean() / ref.mean())
        res['strategies'][key]['relative_ci'] = S.pct_ci(rel)
    assert abs(res['strategies']['single']['theta_sim'] - R3['S12']['A_base']['theta_sim']) < 0.01, 'model A reference not reproduced'
    res['model_A_frozen_theta_sim'] = R3['S12']['A_base']['theta_sim']
    res['observed_theta'] = R3['primary']['estimate']
    OUT.write_text(json.dumps(res, indent=1, default=float) + '\n', encoding='utf-8')
    print('saved', OUT)
if __name__ == '__main__':
    main()
