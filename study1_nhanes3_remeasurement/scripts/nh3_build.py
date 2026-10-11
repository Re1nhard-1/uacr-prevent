import math
import numpy as np
import pandas as pd
import pyreadstat
import nh3_testretest_stage12 as s12
from nh3_testretest_stage12 import num, code, sbp23, smoking, bp_treated, srd, dm_gate, dm_protocol, egfr, cal, uacr
POS = {k: dict(v) for k, v in s12.POS.items()}
POS['exam'].update({'MAPA1': (5072, 5074)})
POS['examse'].update({'MARA1': (3308, 3310)})
POS['adult'].update({'HAE9D': (1628, 1628)})

def load():
    raw = {f: s12.read(f, POS[f]) for f in ['exam', 'lab', 'adult', 'examse', 'labse']}
    rx, _ = pyreadstat.read_xport(s12.DATA / 'pupremed.xpt')
    rx['SEQN'] = rx.SEQN.astype(int)
    return (raw, rx)

def rx_sets(rx):
    statin_yes = set(rx.loc[rx.HQRXCODE.isin(s12.STATIN), 'SEQN'])
    unknown = set(rx.groupby('SEQN').HQKX1B.apply(lambda s: (s == 8).all()).loc[lambda s: s].index) - statin_yes

    def has_class(codes):
        m = np.zeros(len(rx), bool)
        for c in ['HQKX9D', 'HQKX9E', 'HQKX9F']:
            if c not in rx.columns:
                continue
            v = rx[c]
            m |= v.astype(str).str.strip().str.lstrip('0').isin({x.lstrip('0') for x in codes}) if v.dtype == object else v.isin([int(x) for x in codes])
        return set(rx.loc[m, 'SEQN'])
    return (statin_yes, unknown, has_class(['0503', '0508']), has_class(['0912']))

def sbp_single(r, p, k):
    return num(r[p + k + '1'])

def sbp_avail(r, p):
    v = [num(r[p + k + '1']) for k in 'GHI']
    v = [x for x in v if x is not None]
    return sum(v) / len(v) if v else None

def cigs(raw):
    v = num(raw)
    if v is None or raw.strip() == '999':
        return None
    return 1 if v > 0 else 0

def build():
    raw, rx = load()
    exam, lab, adult, examse, labse = (raw[k] for k in ['exam', 'lab', 'adult', 'examse', 'labse'])
    statin_yes, unknown, angina_users, lipid_class = rx_sets(rx)
    mec2 = sorted((s for s, r in examse.items() if r['PERTECH'].strip()))
    home2 = {s for s, r in examse.items() if r['HXRRDAYS'].strip()}
    age_m = lambda s: num(exam[s]['MXPAXTMR'])
    elig = [s for s in mec2 if s in exam and code(exam[s]['DMPSTAT']) == 2 and (age_m(s) is not None) and (360 <= age_m(s) <= 959) and (s in adult) and (not any((code(adult[s][k]) == 1 for k in ('HAF10', 'HAC1D', 'HAC1C')))) and (not (code(exam[s]['MAPF12R']) == 1 or code(examse[s]['MARF12R']) == 1)) and all((code(adult[s][k]) == 2 for k in ('HAF10', 'HAC1D', 'HAC1C')))]
    assert len(elig) == 979
    rows = []
    for s in elig:
        e1, e2, a = (exam[s], examse[s], adult[s])
        l1, l2 = (lab.get(s), labse.get(s))
        complete = l1 is not None and l2 is not None and (code(e1['HSSEX']) in (1, 2)) and (sbp23(e1, 'PEP6') is not None) and (sbp23(e2, 'PER6') is not None) and (num(l1['TCP']) is not None) and (num(l2['TCR']) is not None) and (num(l1['HDP']) is not None) and (num(l2['HDR']) is not None) and (num(l1['CEP']) is not None) and (num(l2['CER']) is not None) and (dm_gate(a, l1['PHPINSU'], l1['GHP']) is not None) and (dm_gate(a, l2['PHRINSU'], l2['GHR']) is not None) and (smoking(a) is not None) and (bp_treated(a) is not None) and (s not in unknown)
        age = math.floor(age_m(s) / 12)
        r = {'SEQN': s, 'complete': complete, 'stratum': code(e1['SDPSTRA6']), 'psu': code(e1['SDPPSU6']), 'phase': code(e1['SDPPHASE']), 'female': int(code(e1['HSSEX']) == 2), 'age': age, 'age_months': age_m(s), 'dmarethn': code(e1['DMARETHN']), 'smoking': smoking(a), 'bptreat': bp_treated(a), 'statin': None if s in unknown else int(s in statin_yes), 'srd': srd(a), 'rr_days': num(e2['MXRRDAYS']), 'peplevel': code(e1['PEPLEVEL']), 'same_arm': e1['PEP6B'] == e2['PER6B'], 'same_examiner': e1['PEPTECH'].strip() == e2['PERTECH'].strip() and e1['PEPTECH'].strip() != '88888', 'examiner_both_bba': e1['PEPTECH'].strip() == '88888' and e2['PERTECH'].strip() == '88888', 'rose': code(a['HAF1']) == 1 and (code(a['HAF2']) == 1 or code(a['HAF3']) == 1) and (code(a['HAF4']) == 1) and (code(a['HAF5']) == 1) and (code(a['HAF6']) == 1), 'angina_user': s in angina_users, 'lipid_proxy': code(a['HAE9D']) == 1 or s in lipid_class, 'preg_unknown': code(e1['HSSEX']) == 2 and 30 <= age <= 59 and (code(e1['MAPF12R']) not in (1, 2) or code(e2['MARF12R']) not in (1, 2)), 'smk_v1_5day': cigs(e1['MAPA1']), 'smk_v2_5day': cigs(e2['MARA1'])}
        if complete:
            dmp = dm_protocol(a, l1, l2)
            fixed_dm = srd(a) == 1 or code(a['HAD6']) == 1 or code(a['HAD10']) == 1 or (code(l1['PHPINSU']) == 1) or (code(l2['PHRINSU']) == 1)
            meth = lambda x: None if not x.strip() or x.strip() == '8' else 'diamat' if x.strip() in '123' else 'affinity'
            r.update({'sbp1': sbp23(e1, 'PEP6'), 'sbp2': sbp23(e2, 'PER6'), 'sbp1_r1': sbp_single(e1, 'PEP6', 'G'), 'sbp1_r2': sbp_single(e1, 'PEP6', 'H'), 'sbp1_r3': sbp_single(e1, 'PEP6', 'I'), 'sbp2_r1': sbp_single(e2, 'PER6', 'G'), 'sbp2_r2': sbp_single(e2, 'PER6', 'H'), 'sbp2_r3': sbp_single(e2, 'PER6', 'I'), 'sbp1_avail': sbp_avail(e1, 'PEP6'), 'sbp2_avail': sbp_avail(e2, 'PER6'), 'tc1': num(l1['TCP']), 'tc2': num(l2['TCR']), 'hdl1': num(l1['HDP']), 'hdl2': num(l2['HDR']), 'scr1_raw': num(l1['CEP']), 'scr2_raw': num(l2['CER']), 'hba1c1': num(l1['GHP']), 'hba1c2': num(l2['GHR']), 'dm1': dmp[0], 'dm2': dmp[1], 'dm1_gate': dm_gate(a, l1['PHPINSU'], l1['GHP']), 'dm2_gate': dm_gate(a, l2['PHRINSU'], l2['GHR']), 'fixed_dm': int(fixed_dm), 'ins_any': int(code(l1['PHPINSU']) == 1 or code(l2['PHRINSU']) == 1), 'hba1c_method_switch': bool(meth(l1['GHPMETH']) and meth(l2['GHRMETH']) and (meth(l1['GHPMETH']) != meth(l2['GHRMETH']))), 'ualb1': num(l1['UBP']), 'ucr1': num(l1['URP']), 'ualb2': num(l2['UBR']), 'ucr2': num(l2['URR']), 'uacr1': uacr(l1['UBP'], l1['URP']), 'uacr2': uacr(l2['UBR'], l2['URR']), 'ucr1_bdl': l1['URP'].strip() == '007.9', 'ucr2_bdl': l2['URR'].strip() == '007.9', 'ualb1_bdl': l1['UBP'].strip() == '0000.4', 'ualb2_bdl': l2['UBR'].strip() == '000.4'})

            def keep79(alb, cr):
                a_, c_ = (num(alb), num(cr))
                return None if a_ is None or c_ is None else 100 * a_ / c_
            r['uacr1_keep79'], r['uacr2_keep79'] = (keep79(l1['UBP'], l1['URP']), keep79(l2['UBR'], l2['URR']))
        rows.append(r)
    E = pd.DataFrame(rows)
    E['cluster'] = E.stratum * 10 + E.psu
    B = E[E.complete].copy().reset_index(drop=True)
    assert len(B) == 759
    for v in (1, 2):
        B[f'scr{v}'] = cal(B[f'scr{v}_raw'])
        B[f'egfr{v}'] = [egfr(c, a, f == 1) for c, a, f in zip(B[f'scr{v}'], B.age, B.female)]
    B['in_uacr'] = B.uacr1.notna() & B.uacr2.notna()
    assert int(B.in_uacr.sum()) == 725
    mec_set = set(mec2)
    T = []
    for s, e1 in exam.items():
        a, l1, ag = (adult.get(s), lab.get(s), num(e1['MXPAXTMR']))
        if code(e1['DMPSTAT']) != 2 or ag is None or (not 360 <= ag <= 959) or (a is None) or (l1 is None):
            continue
        if any((code(a[k]) == 1 for k in ('HAF10', 'HAC1D', 'HAC1C'))) or code(e1['MAPF12R']) == 1:
            continue
        if not all((code(a[k]) == 2 for k in ('HAF10', 'HAC1D', 'HAC1C'))):
            continue
        w = num(e1['WTPFEX6'])
        if not (code(e1['HSSEX']) in (1, 2) and sbp23(e1, 'PEP6') is not None and (num(l1['TCP']) is not None) and (num(l1['HDP']) is not None) and (num(l1['CEP']) is not None) and (dm_gate(a, l1['PHPINSU'], l1['GHP']) is not None) and (smoking(a) is not None) and (bp_treated(a) is not None) and (s not in unknown) and (w is not None) and (w > 0)):
            continue
        age = math.floor(ag / 12)
        fem = code(e1['HSSEX']) == 2
        T.append({'SEQN': s, 'group': 'mec' if s in mec_set else 'home' if s in home2 else 'none', 'w': w, 'stratum': code(e1['SDPSTRA6']), 'psu': code(e1['SDPPSU6']), 'female': int(fem), 'age': age, 'sbp1': sbp23(e1, 'PEP6'), 'tc1': num(l1['TCP']), 'hdl1': num(l1['HDP']), 'egfr1': egfr(cal(num(l1['CEP'])), age, fem), 'dm1': dm_gate(a, l1['PHPINSU'], l1['GHP']), 'smoking': smoking(a), 'bptreat': bp_treated(a), 'statin': int(s in statin_yes), 'uacr1': uacr(l1['UBP'], l1['URP'])})
    T = pd.DataFrame(T)
    assert len(T) == 9478 and int(T.uacr1.notna().sum()) == 9316
    return (E, B, T)
