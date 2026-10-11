from pathlib import Path
import json
import math
import numpy as np
import pandas as pd
import nh3_build
import nh3_models as M
import nh3_stats as S
from nh3_stage3 import visit_inputs, fixed, lp_of, B_BOOT, SEED
ROOT = Path(__file__).resolve().parents[1]
ST = ROOT / 'studies/nhanes3_testretest'
OUT = ST / 'analysis_v1_corrections'

def main():
    OUT.mkdir(exist_ok=True)
    r3 = json.loads((ST / 'analysis_v1/stage3/stage3_results.json').read_text(encoding='utf-8'))
    E, B, T = nh3_build.build()
    n = len(B)
    cl = B.cluster.to_numpy()
    F = fixed(B)
    I1, I2 = (visit_inputs(B, 1), visit_inputs(B, 2))
    L1, L2 = (lp_of(F, I1, 'base'), lp_of(F, I2, 'base'))
    chB = (M.category(M.risk_from_lp(L1)) != M.category(M.risk_from_lp(L2))).astype(float)
    assert abs(chB.mean() - r3['primary']['estimate']) < 1e-12
    W = S.boot_weights(cl, B_BOOT, SEED)
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
        return (float(np.sum(wv * chs) / np.sum(wv)), chs)
    th_s, chs = s9_theta(np.ones(n))
    assert abs(th_s - r3['S9']['theta_shift_removed']) < 1e-12
    old = [chB.mean() - s9_theta(W[b])[0] for b in range(B_BOOT)]
    assert np.allclose(S.pct_ci(old), r3['S9']['ci']), (S.pct_ci(old), r3['S9']['ci'])
    new = [float(np.sum(W[b] * chB) / np.sum(W[b])) - s9_theta(W[b])[0] for b in range(B_BOOT)]
    res = {'label': 'Correction of a documented coding error (2026-10-09, after results); original stage3_results.json S9 interval retained for the record', 'theta': float(chB.mean()), 'theta_shift_removed': th_s, 'difference': float(chB.mean() - th_s), 'ci_corrected_bootstrap': S.pct_ci(new), 'ci_frozen_erroneous': r3['S9']['ci'], 'cr1_wald_paired': S.paired_ci(chB - chs, cl), 'participants_differing_between_scenarios': int((chB != chs).sum()), 'changed_only_original': int(((chB == 1) & (chs == 0)).sum()), 'changed_only_shift_removed': int(((chB == 0) & (chs == 1)).sum())}
    Pp = pd.read_csv(ST / 'analysis_v1/stage3/stage3_participants.csv')
    C0 = M.category(Pp.riskB_m0.to_numpy())
    Cf = M.category(Pp.riskB_m15.to_numpy())
    sw = np.column_stack([M.category(Pp[f'riskB_m{1 << k}'].to_numpy()) != C0 for k in range(4)])
    full = Cf != C0
    k_sw = sw.sum(1)
    res['single_swap_overlap'] = {'sum_of_single_swap_events': int(sw.sum()), 'theta_events': int(full.sum()), 'participants_changing_under_any_single_swap': int((k_sw > 0).sum()), 'participants_changing_under_two_or_more_single_swaps': int((k_sw >= 2).sum()), 'changed_under_a_single_swap_but_not_full': int(((k_sw > 0) & ~full).sum()), 'changed_under_full_but_no_single_swap': int((full & (k_sw == 0)).sum())}
    (OUT / 's9_and_swap_counts.json').write_text(json.dumps(res, indent=1, default=float) + '\n', encoding='utf-8')
    print(json.dumps(res, indent=1, default=float))
if __name__ == '__main__':
    main()
