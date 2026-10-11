from pathlib import Path
import json
import math
import numpy as np
import pandas as pd
import nh3_build
import nh3_stats as S
import nh3_models as M
ROOT = Path(__file__).resolve().parents[1]
ST = ROOT / 'studies/nhanes3_testretest'
OUT = ST / 'analysis_v1_posthoc'
GROUPS = ['BP', 'LIPID', 'KIDNEY', 'GLYC', 'URINE']

def pce_cat(r):
    return np.searchsorted([5.0, 7.5, 20.0], r, side='right')

def cvw_log(x1, x2):
    d = np.log(np.asarray(x2, float)) - np.log(np.asarray(x1, float))
    s = np.std(d, ddof=1) / math.sqrt(2)
    return 100 * math.sqrt(math.exp(s * s) - 1)

def main():
    r3 = json.loads((ST / 'analysis_v1/stage3/stage3_results.json').read_text(encoding='utf-8'))
    P = pd.read_csv(ST / 'analysis_v1/stage3/stage3_participants.csv')
    res = {'label': 'POST HOC 2026-10-09 (after manuscript review wf_47ea3196-9ae); descriptive or reporting-only; frozen outputs unchanged'}
    s11 = P[P.pce_v1.notna()].copy()
    n = len(s11)
    assert n == r3['S11']['n']
    full, sbp = (s11.pce_x75_full.to_numpy(int), s11.pce_x75_sbp.to_numpy(int))
    assert full.sum() == r3['S11']['cross75_full']['events'] and sbp.sum() == r3['S11']['cross75_sbp_only']['events']
    b, c = (int(((full == 1) & (sbp == 0)).sum()), int(((full == 0) & (sbp == 1)).sum()))
    lo, hi = S.tango(b, c, n)
    catf = s11.pce_cat_full.to_numpy(int)
    cats = (pce_cat(s11.pce_v1.to_numpy()) != pce_cat(s11.pce_sbp_only.to_numpy())).astype(int)
    assert cats.sum() == r3['S11']['cat_sbp_only']['events'] and catf.sum() == r3['S11']['cat_full']['events']
    bc, cc = (int(((catf == 1) & (cats == 0)).sum()), int(((catf == 0) & (cats == 1)).sum()))
    res['tango_S11_cross75'] = {'n': n, 'full_only': b, 'sbp_only_only': c, 'discordant': b + c, 'estimate': (b - c) / n, 'tango95': [lo, hi], 'cr1_95': [r3['S11']['cross75_full_minus_sbp']['lower'], r3['S11']['cross75_full_minus_sbp']['upper']], 'wider': 'CR1' if r3['S11']['cross75_full_minus_sbp']['upper'] - r3['S11']['cross75_full_minus_sbp']['lower'] > hi - lo else 'Tango'}
    res['discordant_S11_category'] = {'full_only': bc, 'sbp_only_only': cc, 'discordant': bc + cc, 'tango_required': bc + cc < 20}
    l0, h0 = S.tango(10, 10, 500)
    assert abs(l0 + h0) < 1e-09 and l0 < 0 < h0
    E, B, T = nh3_build.build()
    raw = cvw_log(B.scr1_raw, B.scr2_raw)
    cal = cvw_log(B.scr1, B.scr2)
    assert abs(cal - r3['S7']['scr_cal']['cvw_total_pct']) < 1e-09
    res['creatinine_scale'] = {'cvw_raw_pct': raw, 'cvw_calibrated_pct': cal, 'ratio_cal_to_raw': cal / raw, 'derived_cvi_raw_pct': {'cva_1.0': math.sqrt(max(raw ** 2 - 1.0 ** 2, 0)), 'cva_2.7': math.sqrt(max(raw ** 2 - 2.7 ** 2, 0))}, 'lacher_cvi_pct': 6.8, 'mean_raw_v1': float(B.scr1_raw.mean()), 'mean_cal_v1': float(B.scr1.mean())}
    dlp = P.lpB_m15 - P.lpB_m0
    assert abs(dlp.std(ddof=1) - r3['S6_base']['sd_dLP']) < 1e-09
    band = pd.cut(P.age, [29, 39, 59, 79], labels=['30-39', '40-59', '60-79'])
    race = P.dmarethn.map({1: 'Non-Hispanic White', 2: 'Non-Hispanic Black', 3: 'Mexican American', 4: 'Other'})
    sub = {}
    for name, g in [('age_band', band), ('bptreat', P.bptreat.map({0: 'No', 1: 'Yes'})), ('race', race), ('sex', P.female.map({0: 'Men', 1: 'Women'}))]:
        sub[name] = {str(k): {'n': int(v.size), 'sd_dLP': float(v.std(ddof=1))} for k, v in dlp.groupby(g, observed=True)}
    res['sd_dLP_by_subgroup'] = {'overall': float(dlp.std(ddof=1)), 'groups': sub}
    tables = {}
    for model, K, pref, mask in [('base', 4, 'riskB_m', np.ones(len(P), bool)), ('uacr', 5, 'riskU_m', P.in_uacr.to_numpy(bool))]:
        Q = P[mask]
        full_ = (1 << K) - 1
        C = np.column_stack([M.category(Q[f'{pref}{m}'].to_numpy()) for m in range(1 << K)])
        rows = []
        for m in range(1, 1 << K):
            names = [GROUPS[j] for j in range(K) if m >> j & 1]
            rows.append({'groups_remeasured': names, 'ref_visit1': float((C[:, m] != C[:, 0]).mean()), 'ref_visit2': float((C[:, full_ ^ m] != C[:, full_]).mean())})
        tables[model] = {'n': int(mask.sum()), 'rows': rows}
    assert abs(tables['base']['rows'][-1]['ref_visit1'] - r3['primary']['estimate']) < 1e-12
    assert abs(tables['base']['rows'][0]['ref_visit1'] - r3['K2_base']['single_swap']['BP']['ref1']) < 1e-12
    res['k2_subset_tables'] = tables
    (OUT / 'review_round1.json').write_text(json.dumps(res, indent=1) + '\n', encoding='utf-8')
    print(json.dumps({k: v for k, v in res.items() if k != 'k2_subset_tables'}, indent=1))
if __name__ == '__main__':
    main()
