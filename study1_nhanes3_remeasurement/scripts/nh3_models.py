from pathlib import Path
import math
import sys
import warnings
import numpy as np
import rdata
from scipy.special import expit
ROOT = Path(__file__).resolve().parents[1]
SOURCES = ROOT / 'studies/uacr_prevent/method_sources'
sys.path.insert(0, str(SOURCES / 'pyprevent'))
from prevent._base import prevent_base
from prevent._uacr import prevent_uacr
with warnings.catch_warnings():
    warnings.filterwarnings('ignore')
    _DAT = rdata.read_rda(SOURCES / 'preventr/R/sysdata.rda')
TABLE = {'base': _DAT['base_10yr'], 'uacr': _DAT['uacr_10yr']}
assert len(TABLE['base']) == 23 and len(TABLE['uacr']) == 25
DOMAIN = {'sbp': (90.0, 200.0), 'tc': (130.0, 320.0), 'hdl': (20.0, 100.0), 'egfr': (15.0, 140.0), 'uacr': (0.1, 25000.0)}

def clamp(x, lo, hi):
    return np.minimum(np.maximum(np.asarray(x, float), lo), hi)

def features(female, age, tc, hdl, sbp, dm, smoking, egfr, bptreat, statin, uacr=None):
    age = (np.asarray(age, float) - 55) / 10
    tc, hdl, sbp, egfr = (np.asarray(v, float) for v in (tc, hdl, sbp, egfr))
    nonhdl = (tc - hdl) * 0.02586 - 3.5
    hd = (hdl * 0.02586 - 1.3) / 0.3
    sbpl = (np.minimum(sbp, 110) - 110) / 20
    sbph = (np.maximum(sbp, 110) - 130) / 20
    gl = (np.minimum(egfr, 60) - 60) / -15
    gh = (np.maximum(egfr, 60) - 90) / -15
    dm, sm, bp, st = (np.asarray(v, float) for v in (dm, smoking, bptreat, statin))
    z = np.zeros_like(age)
    one = np.ones_like(age)
    cols = [age, nonhdl, hd, sbpl, sbph, dm, sm, z, z, gl, gh, bp, st, bp * sbph, st * nonhdl, age * nonhdl, age * hd, age * sbph, age * dm, age * sm, z, age * gl]
    if uacr is not None:
        cols += [np.log(np.asarray(uacr, float)), z]
    return np.column_stack(cols + [one])

def lp(female, *args, model='base', outcome='ascvd', **kw):
    X = features(female, *args, **kw)
    t = TABLE[model]
    bf = np.asarray(t[f'female_{outcome}'], float)
    bm = np.asarray(t[f'male_{outcome}'], float)
    f = np.asarray(female, int) == 1
    return np.where(f, X @ bf, X @ bm)

def risk_from_lp(v):
    return 100 * expit(v)

def guard(age, female, dm, smoking, bptreat, statin, sbp, tc, hdl, egfr, uacr=None, sbp_hi=200.0):
    assert float(age) == int(age) and 30 <= age <= 79, ('age', age)
    for k, v in dict(female=female, dm=dm, smoking=smoking, bptreat=bptreat, statin=statin).items():
        assert v in (0, 1), (k, v)
    for k, v, (lo, hi) in [('sbp', sbp, (90.0, sbp_hi)), ('tc', tc, DOMAIN['tc']), ('hdl', hdl, DOMAIN['hdl']), ('egfr', egfr, DOMAIN['egfr'])]:
        assert math.isfinite(v) and lo <= v <= hi, (k, v)
    if uacr is not None:
        assert math.isfinite(uacr) and 0.1 <= uacr <= 25000, ('uacr', uacr)

def scalar(model, age, female, dm, smoking, bptreat, statin, sbp, tc, hdl, egfr, uacr=None, sbp_hi=200.0):
    guard(age, female, dm, smoking, bptreat, statin, sbp, tc, hdl, egfr, uacr if model == 'uacr' else None, sbp_hi)
    args = dict(sex=float(female), age=float(age), tc=float(tc), hdl=float(hdl), sbp=float(sbp), dm=float(dm), smoking=float(smoking), bmi=float('nan'), egfr=float(egfr), bptreat=float(bptreat), statin=float(statin))
    if model == 'base':
        r = prevent_base(**args)
        return (r['prevent_base_10yr_ASCVD'], r['prevent_base_10yr_CVD'])
    assert uacr is not None and math.isfinite(uacr)
    r = prevent_uacr(**args, uacr=float(uacr))
    return (r['prevent_uacr_10yr_ASCVD'], r['prevent_uacr_10yr_CVD'])

def ckd_epi_2021(scr, age, female):
    scr, age, female = (np.asarray(v, float) for v in (scr, age, female))
    k = np.where(female == 1, 0.7, 0.9)
    a = np.where(female == 1, -0.241, -0.302)
    r = scr / k
    return 142 * np.minimum(r, 1) ** a * np.maximum(r, 1) ** (-1.2) * 0.9938 ** age * np.where(female == 1, 1.012, 1.0)

def category(r, cuts=(3.0, 5.0, 10.0)):
    return np.searchsorted(np.asarray(cuts), r, side='right')
