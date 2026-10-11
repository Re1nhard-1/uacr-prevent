"""Descriptive numbers for the 2026-statin-guideline research letter, aggregated from audited result files (E9 rule v2, E10) plus participant"""
from pathlib import Path
import json
import numpy as np
import nh3_build
import nh3_testretest_stage12 as s12
from nh3_ext3_statin import sampson
ROOT = Path(__file__).resolve().parents[1]
PH = ROOT / 'studies/nhanes3_testretest/analysis_v1_posthoc'
OUT = PH / 'letter_numbers.json'

def main():
    E9 = json.loads((PH / 'ext3_statin_results.json').read_text(encoding='utf-8'))
    V = E9['verified_rule']['primary']
    tr = np.array(V['transitions_NLRCSHI'])
    lab = list('NLRCSHI')
    grp = {'none': ['N', 'L'], 'cor2a_only': ['R', 'C'], 'cor1': ['S', 'H', 'I']}
    gi = {g: [lab.index(x) for x in v] for g, v in grp.items()}
    tab = {g1: {g2: int(tr[np.ix_(i1, i2)].sum()) for g2, i2 in gi.items()} for g1, i1 in gi.items()}
    n = int(tr.sum())
    assert n == V['n']
    off = n - sum((tab[g][g] for g in grp))
    cor1_change = sum((tab[a][b] for a in grp for b in grp if (a == 'cor1') != (b == 'cor1')))
    any_change = sum((tab[a][b] for a in grp for b in grp if (a == 'none') != (b == 'none')))
    assert cor1_change == V['primary_recommended_change']['events'] and any_change == V['recommended_or_reasonable_change']['events']
    E, B, T = nh3_build.build()
    keep = set(B.SEQN)
    l1 = s12.read('lab', {'TGP': (1606, 1609)}, keep=keep)
    l2 = s12.read('labse', {'TGR': (358, 361)}, keep=keep)
    g_ = lambda d, s, k: s12.num(d[s][k]) if s in d else None
    tg1 = np.array([g_(l1, s, 'TGP') for s in B.SEQN], dtype=float)
    tg2 = np.array([g_(l2, s, 'TGR') for s in B.SEQN], dtype=float)
    sel = (B.statin.to_numpy() != 1) & np.isfinite(tg1) & np.isfinite(tg2) & (tg1 <= 800) & (tg2 <= 800)
    assert sel.sum() == n
    Bs = B[sel]
    age = Bs.age.to_numpy(float)
    days = Bs.rr_days.to_numpy(float)
    ldl1 = sampson(Bs.tc1.to_numpy(float), Bs.hdl1.to_numpy(float), tg1[sel])
    race = {1: 'non-Hispanic White', 2: 'non-Hispanic Black', 3: 'Mexican American', 4: 'other'}
    desc = {'n': n, 'age_mean': float(age.mean()), 'age_sd': float(age.std(ddof=1)), 'age_range': [float(age.min()), float(age.max())], 'female_n': int(Bs.female.sum()), 'female_pct': float(100 * Bs.female.mean()), 'race_ethnicity_n': {race.get(int(k), str(k)): int(v) for k, v in Bs.dmarethn.value_counts().sort_index().items()}, 'interval_days_median': float(np.nanmedian(days)), 'interval_days_iqr': [float(np.nanpercentile(days, 25)), float(np.nanpercentile(days, 75))], 'interval_days_range': [float(np.nanmin(days)), float(np.nanmax(days))], 'interval_missing': int(np.isnan(days).sum()), 'aged_30_59_n': int(((age >= 30) & (age <= 59)).sum()), 'ldl_v1_median': float(np.median(ldl1))}
    import nh3_models as M
    import nh3_stats as S
    from nh3_stage3 import fixed, lp_of, visit_inputs
    from nh3_ext3_statin import tiers, risk30
    F = fixed(B)
    I1, I2 = (visit_inputs(B, 1), visit_inputs(B, 2))
    r1, r2 = (M.risk_from_lp(lp_of(F, I1, 'base')), M.risk_from_lp(lp_of(F, I2, 'base')))
    q1, q2 = (risk30(F, I1), risk30(F, I2))
    L1 = sampson(B.tc1.to_numpy(float), B.hdl1.to_numpy(float), tg1)
    L2 = sampson(B.tc2.to_numpy(float), B.hdl2.to_numpy(float), tg2)
    ag = B.age.to_numpy(float)
    t1 = tiers(L1, r1, I1['dm'], B.egfr1.to_numpy(float), ag, 'v2', q1)
    t2 = tiers(L2, r2, I2['dm'], B.egfr2.to_numpy(float), ag, 'v2', q2)
    lev = lambda t: np.select([np.isin(t, grp['cor1']), np.isin(t, grp['cor2a_only'])], [2, 1], 0)
    v1, v2 = (lev(t1)[sel], lev(t2)[sel])
    assert all((int(((v1 == a) & (v2 == b)).sum()) == tab[ga][gb] for a, ga in enumerate(['none', 'cor2a_only', 'cor1']) for b, gb in enumerate(['none', 'cor2a_only', 'cor1'])))
    r_ = S.prop_ci((v1 != v2).astype(float), B.cluster.to_numpy()[sel])
    three_level = {k: r_[k] for k in ['n', 'events', 'estimate', 'lower', 'upper', 'method']}
    out = {'label': __doc__.split('\n')[0], 'three_level_change': three_level, 'sources': ['ext3_statin_results.json (verified_rule)', 'ext4_30yr_results.json'], 'recommendation_transitions': {'groups': {k: ''.join(v) for k, v in grp.items()}, 'rows_visit1_cols_visit2': tab, 'n': n, 'cor1_status_changed': cor1_change, 'any_recommendation_changed': any_change, 'off_diagonal_3level': off, 'visit1_totals': {g: sum(tab[g].values()) for g in grp}, 'visit2_totals': {g: sum((tab[a][g] for a in grp)) for g in grp}}, 'descriptives_e9_set': desc}
    OUT.write_text(json.dumps(out, indent=1) + '\n', encoding='utf-8')
    print(json.dumps(out, indent=1))
if __name__ == '__main__':
    main()
