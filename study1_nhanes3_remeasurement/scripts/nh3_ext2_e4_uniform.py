"""Post-audit addition to E4 (2026-10-11; third addendum of extensions2_protocol_20261010.md, written before computing): the kappa interval"""
from pathlib import Path
import hashlib
import json
import math
import numpy as np
import pandas as pd
import nh3_stats as S
from nh3_ext2_analyses import npmle_sq, kappa_theta, kappa_r, LK
ROOT = Path(__file__).resolve().parents[1]
ST = ROOT / 'studies/nhanes3_testretest'
PH = ST / 'analysis_v1_posthoc'
PROTO = ST / 'extensions2_protocol_20261010.md'
OUT = PH / 'ext2_e4_uniform_ci.json'
SEED_BOOT, N_BOOT, MAXIT = (20261009, 1000, 20000)
_W = {}

def _init(m, d, ch, W, risks):
    _W.update(m=m, d=d, ch=ch, W=W, risks=risks)

def _fit(w, maxit):
    m, d, ch = (_W['m'], _W['d'], _W['ch'])
    delta = np.sum(w * d) / np.sum(w)
    sigma = math.sqrt(np.sum(w * (d - delta) ** 2) / (np.sum(w) - 1) / 2)
    th_obs = np.sum(w * ch) / np.sum(w)
    g, p, info = npmle_sq(m, sigma / math.sqrt(2), w, maxit=maxit)
    k_, t_ = kappa_theta(g, p, delta, sigma)
    return {'kappa_model': k_, 'theta_model': t_, 'theta_obs': float(th_obs), 'kappa': k_ * th_obs / t_, 'kappa_r': kappa_r(g, p, delta, sigma, _W['risks']), 'converged': info['converged'], 'em_evaluations': info['em_evaluations']}

def _rep(b):
    return _fit(_W['W'][b], MAXIT)

def main():
    ptxt = PROTO.read_bytes()
    P = pd.read_csv(ST / 'analysis_v1/stage3/stage3_participants.csv')
    L1, L2, ch = (P.lpB_m0.to_numpy(float), P.lpB_m15.to_numpy(float), P.changeB.to_numpy(float))
    d, m = (L2 - L1, (L1 + L2) / 2)
    risks = [e['risk_pct'] for e in LK['lookup']]
    X2 = json.loads((PH / 'ext2_results.json').read_text(encoding='utf-8'))
    E4 = X2['E4']
    _init(m, d, ch, None, risks)
    full = _fit(np.ones(len(d)), 60000)
    chk = abs(full['kappa'] - E4['primary']['kappa'])
    W = S.boot_weights(P.cluster.to_numpy(), 4000, SEED_BOOT)[:N_BOOT]
    from multiprocessing import Pool
    with Pool(12, initializer=_init, initargs=(m, d, ch, W, risks)) as pool:
        reps = pool.map(_rep, range(N_BOOT), chunksize=8)
    out = {'label': __doc__.split('\n')[0], 'protocol': str(PROTO.relative_to(ROOT)), 'protocol_sha256': hashlib.sha256(ptxt).hexdigest(), 'protocol_sha256_before_third_addendum': hashlib.sha256(ptxt.split(b'\n## Third addendum')[0]).hexdigest(), 'n': int(len(d)), 'n_boot': N_BOOT, 'refit_cap_em_evaluations': MAXIT, 'replicates_converged_at_tol': int(sum((r['converged'] for r in reps))), 'full_sample_squarem': {k: full[k] for k in ['kappa_model', 'theta_model', 'kappa', 'converged', 'em_evaluations']}, 'e4_point_estimate': {k: E4['primary'][k] for k in ['kappa_model', 'theta_model', 'kappa']}, 'abs_diff_kappa_squarem_vs_plain_em': chk, 'ci_warm_start_e4': E4['ci'], 'ci_uniform_start': {'kappa': S.pct_ci([r['kappa'] for r in reps]), 'kappa_model': S.pct_ci([r['kappa_model'] for r in reps]), 'ratio_kappa_theta': S.pct_ci([r['kappa_model'] / r['theta_model'] for r in reps]), 'kappa_r': {k: S.pct_ci([r['kappa_r'][k] for r in reps]) for k in reps[0]['kappa_r']}}}
    assert chk < 0.0002, ('SQUAREM full-sample fit differs from E4 plain EM', chk)
    OUT.write_text(json.dumps(out, indent=1, default=lambda o: o.item() if hasattr(o, 'item') else str(o)) + '\n', encoding='utf-8')
    print('saved', OUT)
    print('kappa', E4['primary']['kappa'], 'warm CI', E4['ci']['kappa'], 'uniform CI', out['ci_uniform_start']['kappa'], 'chk', chk)
if __name__ == '__main__':
    main()
