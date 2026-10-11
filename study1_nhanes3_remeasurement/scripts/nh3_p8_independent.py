import csv
import json
import math
from itertools import combinations
from pathlib import Path
D = Path(__file__).resolve().parents[1] / 'studies/nhanes3_testretest/analysis_v1/stage3'
GROUPS = ['BP', 'LIPID', 'KIDNEY', 'GLYC', 'URINE']
KEYS = {'BP': ['sbp'], 'LIPID': ['tc', 'hdl'], 'KIDNEY': ['egfr'], 'GLYC': ['dm'], 'URINE': ['uacr']}

def cat(r):
    return sum((r >= c for c in (3.0, 5.0, 10.0)))

def lp(coef, female, age, tc, hdl, sbp, dm, sm, egfr, bp, st, uacr=None):
    a = (age - 55) / 10
    nonhdl = (tc - hdl) * 0.02586 - 3.5
    hd = (hdl * 0.02586 - 1.3) / 0.3
    sl = (min(sbp, 110) - 110) / 20
    sh = (max(sbp, 110) - 130) / 20
    gl = (min(egfr, 60) - 60) / -15
    gh = (max(egfr, 60) - 90) / -15
    x = [a, nonhdl, hd, sl, sh, dm, sm, 0, 0, gl, gh, bp, st, bp * sh, st * nonhdl, a * nonhdl, a * hd, a * sh, a * dm, a * sm, 0, a * gl]
    if uacr is not None:
        x += [math.log(uacr), 0]
    x.append(1)
    b = coef['female' if female == 1 else 'male']
    return math.fsum((xi * bi for xi, bi in zip(x, b)))

def risk(v):
    return 100 / (1 + math.exp(-v))

def shapley_mean(D, K, rows):
    tot = [0.0] * K
    for d in rows:
        for k in range(K):
            others = [j for j in range(K) if j != k]
            for s in range(K):
                w = math.factorial(s) * math.factorial(K - s - 1) / math.factorial(K)
                for comb in combinations(others, s):
                    m = sum((1 << j for j in comb))
                    tot[k] += w * (d[m | 1 << k] - d[m])
    return [t / len(rows) for t in tot]

def main():
    rows = list(csv.DictReader(open(D / 'stage3_participants.csv', encoding='utf-8')))
    co = json.loads((D / 'preventr_coefficients.json').read_text(encoding='utf-8'))
    res = json.loads((D / 'stage3_results.json').read_text(encoding='utf-8'))
    f = lambda r, k: float(r[k])
    worst = {'base': 0.0, 'uacr': 0.0, 'total_cvd': 0.0}
    D1B, D2B, D1U, D2U, ch, chU, clu, cluU = ([], [], [], [], [], [], [], [])
    for r in rows:
        fx = dict(female=f(r, 'female'), age=f(r, 'age'), sm=f(r, 'smoking'), bp=f(r, 'bptreat'), st=f(r, 'statin'))
        v = {i: {k: f(r, f'{k}{i}_used') for k in ['sbp', 'tc', 'hdl', 'egfr', 'dm']} for i in (1, 2)}
        inU = r['in_uacr'] == 'True'
        if inU:
            for i in (1, 2):
                v[i]['uacr'] = f(r, f'uacr{i}_used')
        for model, K, pref in [('base', 4, 'riskB_m'), ('uacr', 5, 'riskU_m')]:
            if model == 'uacr' and (not inU):
                continue
            cats = []
            for m in range(1 << K):
                h = {}
                for j, g in enumerate(GROUPS[:K]):
                    for k in KEYS[g]:
                        h[k] = v[2 if m >> j & 1 else 1][k]
                coef = {'female': co[model]['female_ascvd'], 'male': co[model]['male_ascvd']}
                rr = risk(lp(coef, fx['female'], fx['age'], h['tc'], h['hdl'], h['sbp'], h['dm'], fx['sm'], h['egfr'], fx['bp'], fx['st'], h.get('uacr') if model == 'uacr' else None))
                worst[model] = max(worst[model], abs(rr - f(r, f'{pref}{m}')))
                cats.append(cat(f(r, f'{pref}{m}')))
            full = (1 << K) - 1
            d1 = [float(c != cats[0]) for c in cats]
            d2 = [float(cats[full ^ m] != cats[full]) for m in range(1 << K)]
            if model == 'base':
                D1B.append(d1)
                D2B.append(d2)
                ch.append(d1[full])
                clu.append(r['cluster'])
            else:
                D1U.append(d1)
                D2U.append(d2)
                chU.append(d1[full])
                cluU.append(r['cluster'])
        coefT = {'female': co['base']['female_total_cvd'], 'male': co['base']['male_total_cvd']}
        for i in (1, 2):
            rr = risk(lp(coefT, fx['female'], fx['age'], v[i]['tc'], v[i]['hdl'], v[i]['sbp'], v[i]['dm'], fx['sm'], v[i]['egfr'], fx['bp'], fx['st']))
            worst['total_cvd'] = max(worst['total_cvd'], abs(rr - f(r, f'riskT{i}')))
    n = len(rows)
    theta = sum(ch) / n
    out = {'risk_max_abs_diff_pp': worst, 'n': n, 'theta': theta, 'theta_matches': abs(theta - res['primary']['estimate']) < 1e-12}
    tr = [[0] * 4 for _ in range(4)]
    for r in rows:
        tr[cat(f(r, 'riskB_m0'))][cat(f(r, 'riskB_m15'))] += 1
    out['transitions_match'] = tr == res['S5_base']['transition_counts']
    sw = {g: sum((d[1 << k] for d in D1B)) / n for k, g in enumerate(GROUPS[:4])}
    out['single_swap_match'] = all((abs(sw[g] - res['K2_base']['single_swap'][g]['ref1']) < 1e-12 for g in sw))
    hf = {g: sum((d[15 ^ 1 << k] for d in D1B)) / n for k, g in enumerate(GROUPS[:4])}
    out['hold_one_fixed_match'] = all((abs(hf[g] - res['K2_base']['hold_one_fixed'][g]['ref1']) < 1e-12 for g in hf))
    psi = [(a + b) / 2 for a, b in zip(shapley_mean(D1B, 4, D1B), shapley_mean(D2B, 4, D2B))]
    out['shapley_base_match'] = all((abs(100 * psi[k] - res['K2_base']['Psi'][g]['pp']) < 1e-09 for k, g in enumerate(GROUPS[:4])))
    psiU = [(a + b) / 2 for a, b in zip(shapley_mean(D1U, 5, D1U), shapley_mean(D2U, 5, D2U))]
    out['shapley_uacr_match'] = all((abs(100 * psiU[k] - res['K2_uacr']['Psi'][g]['pp']) < 1e-09 for k, g in enumerate(GROUPS)))
    out['shapley_sums_to_theta'] = abs(sum(psi) - theta) < 1e-12
    out['sanity_no_swap_zero'] = sum((d[0] for d in D1B)) == 0
    x75 = [float((f(r, 'riskT1') >= 7.5) != (f(r, 'riskT2') >= 7.5)) for r in rows]
    out['s10_match'] = all((abs(a - f(r, 'cvd_cross75')) == 0 for a, r in zip(x75, rows)))
    s11 = [r for r in rows if r['pce_v1'] not in ('', 'nan')]
    out['s11_match'] = all((float((f(r, 'pce_v1') >= 7.5) != (f(r, 'pce_full') >= 7.5)) == f(r, 'pce_x75_full') for r in s11)) and len(s11) == res['S11']['n']
    groups = {}
    for y, c in zip(ch, clu):
        groups.setdefault(c, []).append(y - theta)
    G = len(groups)
    V = G / (G - 1) * math.fsum((math.fsum(v) ** 2 for v in groups.values())) / n ** 2
    Vres = res['primary']['deff'] * theta * (1 - theta) / n
    out['cr1_variance'] = V
    out['cr1_variance_match'] = abs(V - Vres) < 1e-10
    out['pass'] = worst['base'] < 1e-10 and worst['uacr'] < 1e-10 and (worst['total_cvd'] < 1e-10) and all((v for k, v in out.items() if k.endswith('match') or k.startswith('sanity') or k == 'shapley_sums_to_theta'))
    (D.parent / 'p8_independent_recomputation.json').write_text(json.dumps(out, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(out, indent=1))
if __name__ == '__main__':
    main()
