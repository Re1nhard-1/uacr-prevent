from pathlib import Path
import json
import numpy as np
import nh3_build
import nh3_models as M
from nh3_stage3 import visit_inputs, fixed, lp_of, hybrid, shapley, GROUPS, PRIMARY
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'studies/nhanes3_testretest/analysis_v1_posthoc'

def attribution(B, scr_cols):
    B = B.copy()
    for v, col in zip((1, 2), scr_cols):
        B[f'egfr{v}'] = M.ckd_epi_2021(B[col], B.age, B.female)
    F = fixed(B)
    I1, I2 = (visit_inputs(B, 1), visit_inputs(B, 2))
    L = np.column_stack([lp_of(F, hybrid(I1, I2, m, GROUPS[:4]), 'base') for m in range(16)])
    C = M.category(M.risk_from_lp(L))
    D1 = (C != C[:, [0]]).astype(float)
    D2 = np.column_stack([C[:, 15 ^ m] != C[:, 15] for m in range(16)]).astype(float)
    psi = (shapley(D1, 4) + shapley(D2, 4)) / 2
    ch = (C[:, 15] != C[:, 0]).astype(float)
    g1 = M.clamp(B.egfr1, *PRIMARY['egfr'])
    g2 = M.clamp(B.egfr2, *PRIMARY['egfr'])
    lo = (g1 < 60) | (g2 < 60)
    return {'theta': float(ch.mean()), 'events': int(ch.sum()), 'Psi_pp': {g: float(100 * psi[:, k].mean()) for k, g in enumerate(GROUPS[:4])}, 'kidney_single_swap_ref1': float(D1[:, 1 << 2].mean()), 'egfr_v1_median': float(np.median(B.egfr1)), 'egfr_lt60_either_visit': int(lo.sum()), 'changers_with_egfr_lt60_either_visit': int((ch.astype(bool) & lo).sum()), 'mean_abs_kidney_lp_term': float(np.mean(np.abs(L[:, 1 << 2] - L[:, 0])))}

def main():
    OUT.mkdir(exist_ok=True)
    E, B, T = nh3_build.build()
    cal = attribution(B, ['scr1', 'scr2'])
    raw = attribution(B, ['scr1_raw', 'scr2_raw'])
    r4 = json.loads((ROOT / 'studies/nhanes3_testretest/analysis_v1/stage4/stage4_results.json').read_text(encoding='utf-8'))
    r3 = json.loads((ROOT / 'studies/nhanes3_testretest/analysis_v1/stage3/stage3_results.json').read_text(encoding='utf-8'))
    assert abs(cal['theta'] - r3['primary']['estimate']) < 1e-12 and abs(raw['theta'] - r4['SA4a_raw_creatinine']['theta']['estimate']) < 1e-12
    assert all((abs(cal['Psi_pp'][g] - r3['K2_base']['Psi'][g]['pp']) < 1e-09 for g in GROUPS[:4]))
    res = {'label': 'POST HOC 2026-10-09: decomposition of SA4a (uncalibrated creatinine) versus primary; descriptive, no intervals', 'calibrated_primary': cal, 'uncalibrated_SA4a': raw, 'checks': 'theta reproduces primary and SA4a exactly; primary Psi reproduces Stage 3 K2_base'}
    (OUT / 'sa4a_decomposition.json').write_text(json.dumps(res, indent=1) + '\n', encoding='utf-8')
    print(json.dumps(res, indent=1))
if __name__ == '__main__':
    main()
