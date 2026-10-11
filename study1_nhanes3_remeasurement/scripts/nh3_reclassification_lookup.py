from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.special import expit, logit
import nh3_stats as S
ROOT = Path(__file__).resolve().parents[1]
ST = ROOT / 'studies/nhanes3_testretest'
R3 = json.loads((ST / 'analysis_v1/stage3/stage3_results.json').read_text(encoding='utf-8'))
OUT = ST / 'analysis_v1_posthoc'
CUTS = (3.0, 5.0, 10.0)
LT = logit(np.array(CUTS) / 100)
SEED = 20261009
B_BOOT = 4000
GRID_D = np.array([0, 0.05, 0.1, 0.2, 0.3, 0.5])

def design(d, near, knots, ind):
    Xb = np.column_stack([np.ones(len(d)), S.rcs_basis(d, knots)])
    return np.column_stack([Xb, near == 1, near == 2]).astype(float) if ind else Xb

def dist_near(lp):
    dist = np.abs(np.asarray(lp)[:, None] - LT[None, :])
    return (dist.min(1), dist.argmin(1))

def refit(lp1, ch, Wm, frozen):
    d, near = dist_near(lp1)
    ev = int(ch.sum())
    ind = ev >= 40
    knots = np.percentile(d, [10, 50, 90])
    X = design(d, near, knots, ind)
    beta = S.logit_fit(X, ch)
    boots = np.array([S.logit_fit(X, ch, Wm[b]) for b in range(Wm.shape[0])])
    assert ind == frozen['threshold_indicator'] and np.allclose(knots, frozen['knots'], atol=1e-12) and np.allclose(beta, frozen['beta'], atol=1e-09), 'coefficients differ'
    worst = 0.0
    for ti, t in enumerate(CUTS):
        Xg = design(GRID_D, np.full(len(GRID_D), ti), knots, ind)
        p, pb = (expit(Xg @ beta), expit(Xg @ boots.T))
        fr = frozen['by_threshold'][str(t)]
        for i, g in enumerate(GRID_D):
            worst = max(worst, abs(p[i] - fr['pi_at_d'][str(g)]), *[abs(a - b) for a, b in zip(S.pct_ci(pb[i]), fr['ci'][str(g)])])
    assert worst < 1e-09, worst
    return (knots, beta, boots, ind, worst)

def curve(risk_pct, knots, beta, boots, ind):
    d, near = dist_near(logit(np.asarray(risk_pct) / 100))
    X = design(d, near, knots, ind)
    p, pb = (expit(X @ beta), expit(X @ boots.T))
    lo, hi = (np.percentile(pb, 2.5, axis=1), np.percentile(pb, 97.5, axis=1))
    return (p, lo, hi, d, near)

def main():
    P = pd.read_csv(ST / 'analysis_v1/stage3/stage3_participants.csv')
    cl = P.cluster.to_numpy()
    U = np.where(P.in_uacr.to_numpy().astype(bool))[0]
    cat = lambda r: np.searchsorted(np.array(CUTS), r, side='right')
    lpB1 = P.lpB_m0.to_numpy(float)
    chB = (cat(P.riskB_m15.to_numpy(float)) != cat(P.riskB_m0.to_numpy(float))).astype(float)
    assert np.array_equal(chB, P.changeB.to_numpy(float)) and int(chB.sum()) == R3['primary']['events']
    rU1, rU2 = (P.riskU_m0.to_numpy(float)[U], P.riskU_m31.to_numpy(float)[U])
    lpU1 = logit(rU1 / 100)
    chU = (cat(rU2) != cat(rU1)).astype(float)
    assert int(chU.sum()) == R3['K1']['theta_U']['events']
    W = S.boot_weights(cl, B_BOOT, SEED)
    WU = W[:, U]
    mods = {'base': refit(lpB1, chB, W, R3['K3_base']['clin']), 'uacr': refit(lpU1, chU, WU, R3['K3_uacr']['clin'])}
    rcv = {'base': R3['S6_base']['rcv_risk'], 'uacr': R3['S6_uacr']['rcv_risk']}
    grid = np.round(np.arange(0.5, 30.0001, 0.05), 2)
    rows = {'risk_pct': grid}
    for m, (knots, beta, boots, ind, _) in mods.items():
        p, lo, hi, d, near = curve(grid, knots, beta, boots, ind)
        rows.update({f'{m}_p': p, f'{m}_lower': lo, f'{m}_upper': hi})
    rows['nearest_threshold_pct'] = np.array(CUTS)[dist_near(logit(grid / 100))[1]]
    rows['logit_distance'] = dist_near(logit(grid / 100))[0]
    G = pd.DataFrame(rows)
    G.to_csv(OUT / 'reclassification_grid.csv', index=False, float_format='%.6f')
    sel = [1, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10, 12, 15, 20]
    table = []
    for r in sel:
        e = {'risk_pct': r, 'nearest_threshold_pct': float(np.array(CUTS)[dist_near(logit([r / 100]))[1][0]])}
        for m, (knots, beta, boots, ind, _) in mods.items():
            p, lo, hi, _, _ = curve([r], knots, beta, boots, ind)
            e[m] = {'p': float(p[0]), 'lower': float(lo[0]), 'upper': float(hi[0])}
        e['repeat_range_base_pct'] = [r * (1 + rcv['base'][0]), r * (1 + rcv['base'][1])]
        table.append(e)
    bands = [(0, 2), (2, 3), (3, 4), (4, 5), (5, 7.1), (7.1, 10), (10, 15), (15, 100)]
    r1 = P.riskB_m0.to_numpy(float)
    obs = [{'band': f'{a}-<{b}', 'n': int(((r1 >= a) & (r1 < b)).sum()), 'changed': int(chB[(r1 >= a) & (r1 < b)].sum())} for a, b in bands]
    res = {'label': 'Post hoc display of the prespecified pi_clin model (K3, Study 1); same 4,000 cluster-bootstrap replicates (seed 20261009)', 'reproduction_max_abs_diff_vs_frozen': {m: v[4] for m, v in mods.items()}, 'models': {m: {'knots': v[0].tolist(), 'beta': v[1].tolist(), 'threshold_indicator': v[3]} for m, v in mods.items()}, 'rcv_risk': rcv, 'logit_cuts': LT.tolist(), 'logit_midpoints_risk_pct': [100 * expit((LT[0] + LT[1]) / 2), 100 * expit((LT[1] + LT[2]) / 2)], 'lookup': table, 'observed_by_band_base': obs, 'interpretation': 'pi_clin conditions on the first-visit risk and therefore includes regression to the mean; it is the probability that a second standardized measurement of all routine inputs about 2 weeks later gives a different category, in the NHANES III analysis sample. The RCV range ignores regression to the mean.'}
    (OUT / 'reclassification_lookup.json').write_text(json.dumps(res, indent=1) + '\n', encoding='utf-8')
    for e in table:
        print(f"{e['risk_pct']:>5}%  near {e['nearest_threshold_pct']:>4}  base {100 * e['base']['p']:5.1f} ({100 * e['base']['lower']:4.1f}-{100 * e['base']['upper']:4.1f})  uacr {100 * e['uacr']['p']:5.1f} ({100 * e['uacr']['lower']:4.1f}-{100 * e['uacr']['upper']:4.1f})  repeat {e['repeat_range_base_pct'][0]:.2f}-{e['repeat_range_base_pct'][1]:.2f}")
    print('reproduction', res['reproduction_max_abs_diff_vs_frozen'])
    print(obs)
if __name__ == '__main__':
    main()
