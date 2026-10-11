"""Post hoc (2026-10-10): chance that remeasurement moves PREVENT 10-year total CVD risk across 7.5%, the threshold the 2025"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.special import expit, logit
import nh3_stats as S
ROOT = Path(__file__).resolve().parents[1]
ST = ROOT / 'studies/nhanes3_testretest'
R3 = json.loads((ST / 'analysis_v1/stage3/stage3_results.json').read_text(encoding='utf-8'))
OUT = ST / 'analysis_v1_posthoc/cvd75_curve.json'

def main():
    P = pd.read_csv(ST / 'analysis_v1/stage3/stage3_participants.csv')
    r1, r2 = (P.riskT1.to_numpy(float), P.riskT2.to_numpy(float))
    x = ((r1 >= 7.5) != (r2 >= 7.5)).astype(float)
    assert np.array_equal(x, P.cvd_cross75.to_numpy(float)) and int(x.sum()) == R3['S10']['all_759']['events']
    d = np.abs(logit(r1 / 100) - logit(0.075))
    knots = np.percentile(d, [10, 50, 90])
    X = lambda dv: np.column_stack([np.ones(len(dv)), S.rcs_basis(dv, knots)])
    beta = S.logit_fit(X(d), x)
    W = S.boot_weights(P.cluster.to_numpy(), 4000, 20261009)
    boots = np.array([S.logit_fit(X(d), x, W[b]) for b in range(W.shape[0])])
    grid = np.round(np.arange(1.0, 30.0001, 0.05), 2)
    dg = np.abs(logit(grid / 100) - logit(0.075))
    p, pb = (expit(X(dg) @ beta), expit(X(dg) @ boots.T))
    lo, hi = (np.percentile(pb, 2.5, axis=1), np.percentile(pb, 97.5, axis=1))
    sel = [3, 4, 5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 10, 12, 15]
    look = []
    for r in sel:
        i = int(np.argmin(np.abs(grid - r)))
        look.append({'risk_pct': r, 'p': float(p[i]), 'lower': float(lo[i]), 'upper': float(hi[i])})
    fine = np.linspace(0, 1.5, 1501)
    pf = expit(X(fine) @ beta)
    zone = {}
    for lev in (0.1, 0.25):
        ok = np.where(pf >= lev)[0]
        if len(ok) and ok[0] == 0:
            ds = fine[np.argmax(pf < lev) - 1]
            zone[str(lev)] = [float(100 * expit(logit(0.075) - ds)), float(100 * expit(logit(0.075) + ds))]
    unt = (P.bptreat == 0) & (P.sbp1_used >= 130) & (P.sbp1_used < 140)
    res = {'label': __doc__.split('\n')[0], 'n': int(len(x)), 'events': int(x.sum()), 'knots': knots.tolist(), 'beta': beta.tolist(), 'pi_at_threshold': float(expit(X(np.array([0.0])) @ beta)[0]), 'pi_at_threshold_ci': S.pct_ci(expit(X(np.array([0.0])) @ boots.T)[0]), 'zones_risk_pct': zone, 'lookup': look, 'grid': {'risk_pct': grid.tolist(), 'p': p.round(5).tolist(), 'lower': lo.round(5).tolist(), 'upper': hi.round(5).tolist()}, 'untreated_sbp130_139': {'n': int(unt.sum()), 'crossing': int(x[unt].sum()), 'risk_median': float(np.median(r1[unt])), 'in_10pct_zone': int(((r1 >= zone['0.1'][0]) & (r1 <= zone['0.1'][1]) & unt).sum()) if '0.1' in zone else None}, 'note': 'post hoc; 23 crossings, so the curve is imprecise; total CVD risk computed with the full-precision preventr table (S10)'}
    OUT.write_text(json.dumps(res, indent=1) + '\n', encoding='utf-8')
    print({k: res[k] for k in ['events', 'pi_at_threshold', 'pi_at_threshold_ci', 'zones_risk_pct', 'untreated_sbp130_139']})
    for e in look:
        print(f"{e['risk_pct']:>5}%  {100 * e['p']:5.1f} ({100 * e['lower']:4.1f}-{100 * e['upper']:4.1f})")
if __name__ == '__main__':
    main()
