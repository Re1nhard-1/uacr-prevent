"""Verify ASCVD-only inputs with unmodified sources, including absent unused BMI."""
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path
import json
import platform
import shutil

import numpy as np
from scipy.special import expit

from verify_prevent_uacr import coefficients, hash_file, prevent_uacr, SOURCES

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'studies/uacr_prevent/ascvd_input_verification_v2'


def from_r_ascvd_coefficients(x, table):
    """Independent pinned R-coefficient check; omit proven-zero BMI terms, no imputation."""
    if x['sex'] not in (0,1):
        raise ValueError('Observed sex0/1 required')
    required = ['age','tc','hdl','sbp','dm','smoking','egfr','bptreat','statin','uacr']
    if not all(np.isfinite(x[k]) for k in required) or not .1<=x['uacr']<=25000:
        raise ValueError('Observed finite ASCVD inputs and UACR within0.1–25000 required')
    assert len(table)==25
    for i in [7,8,20]:
        assert 'BMI' in str(table['beta_coefficients'].iloc[i])
    column = 'female_ascvd' if x['sex']==1 else 'male_ascvd'
    beta = np.asarray(table[column],dtype=float)
    assert np.isfinite(beta).all() and np.all(beta[[7,8,20]]==0)
    age=(x['age']-55)/10
    nonhdl=(x['tc']-x['hdl'])*.02586-3.5
    hdl=(x['hdl']*.02586-1.3)/.3
    sbpl=(min(x['sbp'],110)-110)/20
    sbph=(max(x['sbp'],110)-130)/20
    gl=(min(x['egfr'],60)-60)/-15
    gh=(max(x['egfr'],60)-90)/-15
    dm,sm,bp,st=(x[k] for k in ['dm','smoking','bptreat','statin'])
    terms=np.array([age,nonhdl,hdl,sbpl,sbph,dm,sm,np.nan,np.nan,gl,gh,bp,st,
                    bp*sbph,st*nonhdl,age*nonhdl,age*hdl,age*sbph,age*dm,age*sm,
                    np.nan,age*gl,np.log(x['uacr']),0.,1.])
    active=beta!=0
    assert np.isfinite(terms[active]).all(), 'A required nonzero-coefficient term is nonfinite'
    result=float(100*expit(terms[active]@beta[active]))
    assert np.isfinite(result) and 0<result<100
    return result


def verify_finite_agreement(actual, reference, tolerance=1e-10):
    errors=np.abs(np.asarray(actual,dtype=float)-np.asarray(reference,dtype=float))
    if errors.size==0 or not np.isfinite(errors).all():
        raise ValueError('All comparison errors must be finite and nonempty')
    if float(errors.max())>=tolerance:
        raise ValueError('Source/reference disagreement exceeds tolerance')
    return float(errors.max())


def main():
    if OUT.exists():
        raise FileExistsError('Refusing to overwrite completed verification')
    for item in json.loads((SOURCES/'manifest.json').read_text()):
        assert hash_file(SOURCES/item['file'])==item['sha256']
    table=coefficients()
    seed=20261005
    rng=np.random.default_rng(seed)
    actual, reference, bmi_invariance=[],[],[]
    cases=[]
    for i in range(360):
        x={'sex':i%2,'age':float(rng.choice([30,55,79])),
           'tc':float(rng.uniform(130,320)),'hdl':float(rng.uniform(20,100)),
           'sbp':float(rng.choice([90,110,180,200])),
           'dm':int(rng.integers(2)),'smoking':int(rng.integers(2)),
           'egfr':float(rng.choice([15,60,140])),
           'bptreat':int(rng.integers(2)),'statin':int(rng.integers(2)),
           'uacr':float(rng.choice([.1,1,30,300,25000]))}
        checked=from_r_ascvd_coefficients(x,table)
        values=[]
        for bmi in [np.nan,5.,25.,90.]:
            score=prevent_uacr(**x,bmi=bmi)['prevent_uacr_10yr_ASCVD']
            actual.append(score); reference.append(checked); values.append(score)
        bmi_invariance.append(float(np.max(values)-np.min(values)))
        cases.append({'inputs':x,'reference_pct':checked,'python_absent_bmi_pct':values[0],
                      'max_bmi_variant_difference_pp':bmi_invariance[-1]})
    maximum=verify_finite_agreement(actual,reference)
    assert max(bmi_invariance)==0
    rejected=0
    for key in ['age','tc','hdl','sbp','dm','egfr','uacr']:
        invalid=dict(cases[0]['inputs']); invalid[key]=np.nan
        try:
            from_r_ascvd_coefficients(invalid,table)
        except ValueError:
            rejected+=1
        else:
            raise AssertionError(f'Missing required {key} was accepted')
    # A NaN can be hidden by Python max depending on order; both orders must fail.
    for values in [[0.,np.nan],[np.nan,0.]]:
        try:
            verify_finite_agreement(values,[0.,0.])
        except ValueError:
            rejected+=1
        else:
            raise AssertionError('Nonfinite comparison silently accepted')
    for uacr in [.099,25000.01]:
        invalid=dict(cases[0]['inputs']); invalid['uacr']=uacr
        try:
            from_r_ascvd_coefficients(invalid,table)
        except ValueError:
            rejected+=1
        else:
            raise AssertionError('Out-of-study-domain UACR accepted')
    report={'status':'PASS_NUMERICAL_ASCVD_INPUT_VERIFICATION_ONLY',
            'executed_utc':datetime.now(timezone.utc).isoformat(),'seed':seed,
            'synthetic_profiles':len(cases),'source_function_evaluations':len(actual),
            'maximum_reference_difference_pp':maximum,
            'maximum_bmi_variant_difference_pp':max(bmi_invariance),
            'missing_or_nonfinite_cases_rejected':rejected,
            'source_manifest_sha256':hash_file(SOURCES/'manifest.json'),
            'script_sha256':hash_file(__file__),
            'command':"& '.venv/Scripts/python.exe' -X utf8 scripts/verify_prevent_ascvd_inputs.py",
            'environment':{'python':platform.python_version(),**{p:version(p) for p in ['numpy','pandas','scipy','rdata']}},
            'scope':'Observed UACR0.1–25000,10yearASCVD,bothsexes,SBP90–200. BMI is not imputed; proven-zero terms are omitted from the independent coefficient check.',
            'supersedes':'ascvd_input_verification/verification.json retained; v2 enforces the declared UACR domain and tests its rejection boundary.',
            'limits':['No participant risks computed.','No R runtime executed.','No clinical calibration/outcome validation.','HF/CVD/30year outputs are outside verified scope.'],
            'cases':cases}
    OUT.mkdir()
    (OUT/'verification.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    shutil.copy2(__file__,OUT/'executed_verification_script.py')
    print(json.dumps({k:v for k,v in report.items() if k!='cases'},indent=2))


if __name__=='__main__':
    main()
