from pathlib import Path
import json
import math
import re
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
sys.path.insert(0, str(Path(__file__).resolve().parent / 'vendor'))
from figure_tools import publication_style, export_figure
from scipy.special import expit, logit
import nh3_stats as S
ROOT = Path(__file__).resolve().parents[1]
A = ROOT / 'studies/nhanes3_testretest/analysis_v1'
OUT = ROOT / 'studies/nhanes3_testretest/figures_v1'
R3 = json.loads((A / 'stage3/stage3_results.json').read_text(encoding='utf-8'))
R4 = json.loads((A / 'stage4/stage4_results.json').read_text(encoding='utf-8'))
P = pd.read_csv(A / 'stage3/stage3_participants.csv')
RRP = json.loads((ROOT / 'studies/nhanes3_testretest/analysis_v1_posthoc/review_round1.json').read_text(encoding='utf-8'))

def _lg(p):
    return math.log(p / (1 - p))
B35, B510 = (100 * expit((_lg(0.03) + _lg(0.05)) / 2), 100 * expit((_lg(0.05) + _lg(0.1)) / 2))

def zone_box(K3):
    out = []
    for lvl in ('0.1', '0.25'):
        z = {t: list(K3['by_threshold'][t]['instability_zones_risk_pct'][lvl]) for t in ['3.0', '5.0', '10.0']}
        z['3.0'][1] = min(z['3.0'][1], B35)
        z['5.0'][0] = max(z['5.0'][0], B35)
        z['5.0'][1] = min(z['5.0'][1], B510)
        z['10.0'][0] = max(z['10.0'][0], B510)
        out.append(f'≥{int(float(lvl) * 100)}%:\n' + ', '.join((f'{a:.1f}–{b:.1f}' for a, b in (z[t] for t in ['3.0', '5.0', '10.0']))))
    return 'Visit-1 risk (%) where change\nprobability ' + '\n'.join(out)
INK, TEAL, ORANGE, GREY, GRID = ('#203443', '#187D84', '#CB6C38', '#B9C3CB', '#E3E8EC')
GC = {'BP': '#3766A0', 'LIPID': '#B0548C', 'KIDNEY': '#4C9A6A', 'GLYC': '#C9A227', 'URINE': '#7A6A58'}
GL = {'BP': 'Blood pressure', 'LIPID': 'Lipids', 'KIDNEY': 'Creatinine (eGFR)', 'GLYC': 'HbA1c (diabetes)', 'URINE': 'Urine ACR'}
SRC = [A / 'stage3/stage3_results.json', A / 'stage4/stage4_results.json', A / 'stage3/stage3_participants.csv', ROOT / 'studies/nhanes3_testretest/manuscript_evidence_v1/study2_evidence.json', ROOT / 'studies/uacr_prevent/final_and_external_v1/primary_participant_results.csv', ROOT / 'studies/nhanes3_testretest/analysis_v1_posthoc/repeat_strategies.json', Path(__file__)]
mm = lambda w, h: (w / 25.4, h / 25.4)

def lab(ax, s, x=-0.14, y=1.04):
    ax.text(x, y, s, transform=ax.transAxes, fontweight='bold', fontsize=10, color=INK, va='bottom')

def pc(e):
    return (100 * e['estimate'], 100 * e['lower'], 100 * e['upper'])

def fig1():
    fig = plt.figure(figsize=mm(180, 112))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 180)
    ax.set_ylim(0, 112)
    ax.axis('off')

    def box(x, y, w, h, t, fc='white', ec=INK, fs=7.5, bold=False):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.4,rounding_size=1.2', fc=fc, ec=ec, lw=0.8))
        ax.text(x + w / 2, y + h / 2, t, ha='center', va='center', fontsize=fs, color=INK, fontweight='bold' if bold else 'normal', linespacing=1.25)

    def arrow(x1, y1, x2, y2):
        ax.annotate('', (x2, y2), (x1, y1), arrowprops=dict(arrowstyle='-|>', lw=0.8, color=INK, mutation_scale=8))
    steps = [('NHANES III second-examination records', 2596), ('MEC second exam (737 home second exams excluded)', 1859), ('Linked to a MEC first exam', 1859), ('Aged 30–79 years at visit 1', 1106), ('No self-reported MI, stroke or heart failure', 990), ('Not pregnant at either visit', 988), ('History determinable (all three answered No)', 979)]
    y = 103
    for i, (t, nn) in enumerate(steps):
        box(8, y - 4.5, 74, 8.5, f'{t}\nn = {nn:,}', fs=6.8)
        if i:
            arrow(45, y + 6.2, 45, y + 4.3)
        y -= 10.8
    box(8, y - 6, 74, 10, 'Complete base PREVENT inputs at both visits\n(analysis set)\nn = 759', fc='#EEF3F6', bold=True, fs=7)
    arrow(45, y + 6.2, 45, y + 4.3)
    box(8, y - 18, 74, 8, 'Also urine ACR at both visits (UACR set)\nn = 725', fc='#EEF3F6', fs=7)
    arrow(45, y - 6.4, 45, y - 9.7)
    box(96, 85, 78, 21, 'Inputs remeasured at visit 2 (median 16 days, 3–48)\nSBP (mean of physician readings 2–3), total and HDL\ncholesterol, creatinine (eGFR), HbA1c, urine ACR\nHeld fixed: age, sex, smoking, antihypertensive\nand statin use (single household interview)', fs=6.6)
    box(96, 65, 78, 13, 'Standardization target: 9,478 eligible MEC participants\nwith complete visit-1 inputs (9,316 with urine ACR),\nweighted by examination weights', fs=6.6)
    box(96, 45, 78, 13, '220 eligible without complete inputs at both visits\n(missing at visit 1/visit 2: SBP 41/104, total cholesterol\n39/93, HDL 45/100, creatinine 52/99; overlapping)', fs=6.6)
    box(96, 25, 78, 13, 'Values outside SBP 90–200 mm Hg, TC 130–320 mg/dL,\nHDL 20–100 mg/dL or eGFR 15–140 set to the nearest\nbound: 49 of 759 participants at either visit', fs=6.6)
    box(96, 2, 78, 16, 'Pipeline screen against Ye et al. 2014 (PCE, SBP only):\nreconstructed cohort 1,424 → 1,190 → 1,119 (Ye: 956)\n→ 1,078 → 910 (Ye: 808); all four criteria passed\nbefore any PREVENT value was computed', fs=6.6)
    for y0, y1 in [(92.2, 95.5), (70.6, 71.5), (32.8, 51.5), (26.4, 31.5)]:
        ax.plot([83.5, 95], [y0, y1], color=GREY, lw=0.7, ls=(0, (2, 2)))
    return (fig, [{'id': 'flow', 'question': 'How was the analysis set derived?', 'units': 'participants', 'sample_unit': 'participant', 'transform': 'none', 'estimate': 'counts', 'interval': 'none; counts', 'n': 759}])

def scatter_panel(ax, r1, r2, title, xl='10-year ASCVD risk, visit 1 (%)', yl='10-year ASCVD risk, visit 2 (%)', second='visit 2'):
    c1 = np.searchsorted([3, 5, 10], r1, side='right')
    c2 = np.searchsorted([3, 5, 10], r2, side='right')
    col = np.where(c2 > c1, ORANGE, np.where(c2 < c1, TEAL, GREY))
    order = np.argsort(c1 != c2)
    ax.scatter(r1[order], r2[order], s=7, c=col[order], lw=0, alpha=0.85, rasterized=True)
    lim = [0.05, 60]
    ax.plot(lim, lim, color='#72808C', lw=0.6)
    for t in (3, 5, 10):
        ax.axvline(t, color=GRID, lw=0.8, zorder=0)
        ax.axhline(t, color=GRID, lw=0.8, zorder=0)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    tk = [0.1, 1, 3, 5, 10, 50]
    ax.set_xticks(tk)
    ax.set_yticks(tk)
    ax.set_xticklabels(['0.1', '1', '3', '5', '10', '50'])
    ax.set_yticklabels(['0.1', '1', '3', '5', '10', '50'])
    ax.minorticks_off()
    ax.set_xlabel(xl)
    ax.set_ylabel(yl)
    ax.set_title(title, fontsize=8, loc='left')
    for l, c in [('Same category', GREY), (f'Higher with {second}', ORANGE), (f'Lower with {second}', TEAL)]:
        ax.scatter([], [], s=12, c=c, label=l)
    ax.legend(loc='upper left', fontsize=6.8, handletextpad=0.2)

def trans_panel(ax, tr, title, xl='Visit-2 category (%)', yl='Visit-1 category (%)'):
    tr = np.array(tr)
    labs = ['<3', '3–<5', '5–<10', '≥10']
    for i in range(4):
        for j in range(4):
            fc = '#F2F4F6' if i == j else '#FBE7DC' if j > i else '#DDEFF0'
            ax.add_patch(plt.Rectangle((j, 3 - i), 1, 1, fc=fc, ec='white', lw=1.5))
            ax.text(j + 0.5, 3 - i + 0.5, f'{tr[i, j]}', ha='center', va='center', fontsize=7.5, color=INK)
    ax.set_xlim(0, 4)
    ax.set_ylim(0, 4)
    ax.set_xticks(np.arange(4) + 0.5)
    ax.set_xticklabels(labs)
    ax.set_yticks(np.arange(4) + 0.5)
    ax.set_yticklabels(labs[::-1])
    ax.set_xlabel(xl)
    ax.set_ylabel(yl)
    ax.set_title(title, fontsize=8, loc='left')
    for s in ax.spines.values():
        s.set_visible(False)
    ax.tick_params(length=0)

def fig2():
    fig = plt.figure(figsize=mm(180, 150))
    a = fig.add_axes([0.08, 0.57, 0.4, 0.37])
    b = fig.add_axes([0.62, 0.6, 0.3, 0.31])
    c = fig.add_axes([0.08, 0.08, 0.4, 0.37])
    d = fig.add_axes([0.62, 0.11, 0.3, 0.31])
    rb1, rb2 = (P.riskB_m0.to_numpy(), P.riskB_m15.to_numpy())
    u = P.in_uacr.to_numpy(bool)
    ru1, ru2 = (P.riskU_m0.to_numpy()[u], P.riskU_m31.to_numpy()[u])
    e, eu = (R3['primary'], R3['K1']['theta_U'])
    scatter_panel(a, rb1, rb2, f"Base PREVENT-ASCVD (n = 759)\nchanged category: {e['events']} ({100 * e['estimate']:.1f}%)")
    trans_panel(b, R3['S5_base']['transition_counts'], 'Base model: visit 1 × visit 2 (counts)')
    scatter_panel(c, ru1, ru2, f"UACR-enhanced PREVENT-ASCVD (n = 725)\nchanged category: {eu['events']} ({100 * eu['estimate']:.1f}%)")
    trans_panel(d, R3['S5_uacr']['transition_counts'], 'UACR model: visit 1 × visit 2 (counts)')
    for ax, l in zip([a, b, c, d], 'ABCD'):
        lab(ax, l, -0.16 if l in 'AC' else -0.22)
    return (fig, [{'id': k, 'question': q, 'units': 'percent risk / counts', 'sample_unit': 'participant', 'transform': 'log10 axes' if k in 'AC' else 'none', 'estimate': 'raw paired observations' if k in 'AC' else 'counts', 'interval': 'none; raw observations', 'n': n_} for k, q, n_ in [('A', 'Base risk at both visits', 759), ('B', 'Base category transitions', 759), ('C', 'UACR-model risk at both visits', 725), ('D', 'UACR-model transitions', 725)]])

def fig3():
    fig = plt.figure(figsize=mm(180, 128))
    a = fig.add_axes([0.17, 0.58, 0.3, 0.34])
    b = fig.add_axes([0.67, 0.58, 0.3, 0.34])
    c = fig.add_axes([0.17, 0.09, 0.3, 0.34])
    d = fig.add_axes([0.67, 0.09, 0.3, 0.34])
    for ax, key, groups, title in [(a, 'K2_base', ['BP', 'LIPID', 'KIDNEY', 'GLYC'], 'Base model: Shapley split of θ = 10.4%'), (b, 'K2_uacr', ['BP', 'LIPID', 'KIDNEY', 'GLYC', 'URINE'], 'UACR model: Shapley split of $\\theta_U$ = 11.6%')]:
        Ps = R3[key]['Psi']
        for i, g in enumerate(groups):
            v = Ps[g]
            ax.barh(i, v['pp'], color=GC[g], height=0.55)
            ax.plot(v['ci_pp'], [i, i], color=INK, lw=1)
            ax.text(max(v['ci_pp'][1], v['pp']) + 0.25, i, f"{v['pp']:.1f}", va='center', fontsize=7, color=INK)
        ax.set_yticks(range(len(groups)))
        ax.set_yticklabels([GL[g] for g in groups])
        ax.invert_yaxis()
        ax.axvline(0, color='#72808C', lw=0.6)
        ax.set_xlim(-0.5, 8)
        ax.set_xlabel('Contribution to $\\theta_U$, pp (95% CI)' if key == 'K2_uacr' else 'Contribution to θ, pp (95% CI)')
        ax.set_title(title, fontsize=8, loc='left')
    a.text(0.98, 0.04, 'BP ranked first in 92% of\nbootstrap replicates', transform=a.transAxes, ha='right', va='bottom', fontsize=6.8, color=INK)
    b.text(0.98, 0.5, 'No group ranked first in ≥90%\n(urine 61%, BP 38%)', transform=b.transAxes, ha='right', va='center', fontsize=6.8, color=INK)
    groups = ['BP', 'LIPID', 'KIDNEY', 'GLYC', 'URINE']
    for i, g in enumerate(groups):
        for key, off, mk in [('K2_base', -0.15, 'o'), ('K2_uacr', 0.15, 's')]:
            if g not in R3[key]['single_swap']:
                continue
            sw = 100 * R3[key]['single_swap'][g]['ref1']
            hf = 100 * R3[key]['hold_one_fixed'][g]['ref1']
            th = 100 * (R3['primary']['estimate'] if key == 'K2_base' else R3['K1']['theta_U']['estimate'])
            c.plot([sw], [i + off], mk, color=GC[g], ms=5, mfc=GC[g] if key == 'K2_base' else 'white')
            c.plot([th - hf], [i + off], mk, color=GC[g], ms=5, alpha=0.45, mfc=GC[g] if key == 'K2_base' else 'white')
    c.set_yticks(range(5))
    c.set_yticklabels([GL[g] for g in groups])
    c.invert_yaxis()
    c.set_xlim(-1, 8.6)
    c.set_ylim(4.6, -0.75)
    c.axvline(0, color='#72808C', lw=0.6)
    c.set_xlabel('Percentage of participants (%)')
    c.set_title('Single-group swap (dark) and θ minus\nhold-one-fixed (light); visit 1 as reference', fontsize=8, loc='left')
    for row, x, l in [(0, 100 * 14 / 808, 'Ye (PCE)'), (0, 5.8, 'Ogunwale (PCE)'), (1, 6.3, 'Cole (simulation)'), (4, 5.34, 'Study 2 (specimen)')]:
        c.plot([x, x], [row - 0.38, row + 0.38], color='#8A96A0', lw=0.9, ls=(0, (2, 2)))
        c.text(x, row - 0.42, l, fontsize=6.5, ha='center', va='bottom', color='#5A6670')
    c.plot([], [], 'o', color=INK, label='Base model')
    c.plot([], [], 's', color=INK, mfc='white', label='UACR model')
    c.plot([], [], color='#8A96A0', lw=0.9, ls=(0, (2, 2)), label='Benchmark (other studies or Study 2)')
    c.legend(fontsize=6.6, loc='center right', bbox_to_anchor=(1.0, 0.4), framealpha=1, facecolor='white')
    for i, g in enumerate(['BP', 'LIPID', 'KIDNEY', 'GLYC']):
        v = R3['K2_base']['continuous'][g]
        d.plot([100 * v['pi']], [i - 0.12], 'o', color=GC[g], ms=5)
        d.plot([100 * x for x in v['pi_ci']], [i - 0.12, i - 0.12], color=GC[g], lw=1)
    for i, g in enumerate(groups):
        v = R3['K2_uacr']['continuous'][g]
        d.plot([100 * v['pi']], [i + 0.12], 's', color=GC[g], ms=5, mfc='white')
        d.plot([100 * x for x in v['pi_ci']], [i + 0.12, i + 0.12], color=GC[g], lw=1)
    d.set_yticks(range(5))
    d.set_yticklabels([GL[g] for g in groups])
    d.invert_yaxis()
    d.axvline(0, color='#72808C', lw=0.6)
    d.set_xlim(-5, 65)
    d.set_xlabel('Share of mean squared logit change (%)')
    d.set_title('Continuous LP partition (95% CI)', fontsize=8, loc='left')
    for ax, l in zip([a, b, c, d], 'ABCD'):
        lab(ax, l, -0.52)
    return (fig, [{'id': 'A', 'question': 'Which inputs drive base-model category change?', 'units': 'percentage points', 'sample_unit': 'participant', 'transform': 'none', 'estimate': 'symmetrized Shapley contribution', 'interval': '95% cluster-bootstrap percentile (98 pseudo-PSUs, B=4000)', 'n': 759}, {'id': 'B', 'question': 'Same for the UACR model', 'units': 'percentage points', 'sample_unit': 'participant', 'transform': 'none', 'estimate': 'symmetrized Shapley contribution', 'interval': '95% cluster-bootstrap percentile (97 clusters, B=4000)', 'n': 725}, {'id': 'C', 'question': 'Single-group swap and hold-one-fixed proportions', 'units': 'percent', 'sample_unit': 'participant', 'transform': 'none', 'estimate': 'proportions', 'interval': 'none; point estimates', 'n': 759}, {'id': 'D', 'question': 'Continuous LP partition', 'units': 'percent share', 'sample_unit': 'participant', 'transform': 'logit scale', 'estimate': 'E[f_k dLP]/E[dLP^2]', 'interval': '95% cluster-bootstrap percentile', 'n': 759}])

def curve(ax, key, name, title, color_by_t=True):
    k = R3[key][name]
    knots = np.array(k['knots'])
    beta = np.array(k['beta'])
    xs = np.linspace(0, 0.6, 301)
    cols = {'3.0': '#3766A0', '5.0': INK, '10.0': '#B0548C'}
    for ti, t in enumerate(['3.0', '5.0', '10.0']):
        dx = (ti - 1) * 0.007
        Xb = np.column_stack([np.ones(len(xs)), S.rcs_basis(xs, knots), np.full(len(xs), ti == 1), np.full(len(xs), ti == 2)]).astype(float)
        ax.plot(xs, 100 * expit(Xb @ beta), color=cols[t], lw=1.2, label=f'Nearest threshold {t[:-2]}%')
        g = k['by_threshold'][t]
        for d_, ci in g['ci'].items():
            if float(d_) <= 0.6:
                ax.plot([float(d_) + dx] * 2, [100 * ci[0], 100 * ci[1]], color=cols[t], lw=0.7, alpha=0.6)
    mids = [0.025, 0.075, 0.15, 0.3, 0.55]
    for bn, m_ in zip(k['binned'], mids):
        p = bn['changes'] / bn['n']
        ax.plot(m_, 100 * p, 'D', color='white', mec=INK, ms=4.5, zorder=5)
        ax.text(m_ + 0.012, 100 * p + 2.5, f"{bn['changes']}/{bn['n']}", fontsize=6.5, ha='left', color='#5A6670')
    ax.plot([], [], 'D', color='white', mec=INK, ms=4.5, label='Observed, binned (all thresholds)')
    if 'mechanism_survival_2d' in k:
        ms_ = k['mechanism_survival_2d']
        ax.plot([float(x) for x in ms_], [100 * v for v in ms_.values()], ':', color=ORANGE, lw=1.2, label='P(|ΔLP| > 2d), closed form')
    ax.axhline(10, color=GRID, lw=0.8, zorder=0)
    ax.axhline(25, color=GRID, lw=0.8, zorder=0)
    for yy in (10, 25):
        ax.text(0.598, yy + 1, f'{yy}%', ha='right', va='bottom', fontsize=6.5, color='#8A96A0')
    ax.set_xlim(0, 0.6)
    ax.set_ylim(0, 100)
    ax.set_xlabel('Distance from nearest threshold (logit units)')
    ax.set_ylabel('Probability of category change (%)')
    ax.set_title(title, fontsize=8, loc='left')
    ax.legend(fontsize=6.6, loc='upper right')

def fig4():
    fig = plt.figure(figsize=mm(180, 140))
    axs = [fig.add_axes([0.08, 0.57, 0.38, 0.36]), fig.add_axes([0.58, 0.57, 0.38, 0.36]), fig.add_axes([0.08, 0.08, 0.38, 0.36]), fig.add_axes([0.58, 0.08, 0.38, 0.36])]
    curve(axs[0], 'K3_base', 'clin', 'Base model, conditional on visit-1 risk ($\\pi_{clin}$)')
    curve(axs[1], 'K3_base', 'mech', 'Base model, conditional on visit-mean risk ($\\pi_{mech}$)')
    curve(axs[2], 'K3_uacr', 'clin', 'UACR model, $\\pi_{clin}$')
    curve(axs[3], 'K3_uacr', 'mech', 'UACR model, $\\pi_{mech}$')
    for ax, key in [(axs[0], 'K3_base'), (axs[2], 'K3_uacr')]:
        ax.text(0.585, 47, zone_box(R3[key]['clin']), ha='right', va='center', fontsize=6.5, color=INK, linespacing=1.25, bbox=dict(boxstyle='round,pad=0.35', fc='white', ec=GRID, lw=0.6))
    z = R3['K3_base']['clin']['by_threshold']
    for ax, l in zip(axs, 'ABCD'):
        lab(ax, l, -0.16)
    return (fig, [{'id': l, 'question': 'Change probability by distance to the nearest category threshold', 'units': 'percent', 'sample_unit': 'participant', 'transform': 'logit-scale distance', 'estimate': 'logistic regression with restricted cubic spline (3 knots) and nearest-threshold indicator', 'interval': 'pointwise 95% cluster-bootstrap percentile at d = 0, 0.05, 0.1, 0.2, 0.3, 0.5', 'n': 759 if l in 'AB' else 725} for l in 'ABCD'])

def eflm_cvi():
    facts = json.loads((ROOT / 'studies/nhanes3_testretest/prefreeze_facts_v1.json').read_text(encoding='utf-8'))['cva']['items']
    pick = {'tc': 'Cholesterol (serum', 'hdl': 'HDL cholesterol', 'scr_cal': 'Creatinine (serum', 'hba1c_log': 'HbA1c'}
    out = {}
    for k, key in pick.items():
        it = [x for x in facts if str(x.get('label', '')).startswith('EFLM') and key in x['label']]
        assert len(it) == 1, k
        out[k] = float(re.search('CVI (\\d+(?:\\.\\d+)?)%', it[0]['value']).group(1))
    assert out == {'tc': 5.4, 'hdl': 5.7, 'scr_cal': 4.4, 'hba1c_log': 1.2}, out
    return out

def fig5():
    fig = plt.figure(figsize=mm(180, 150))
    a = fig.add_axes([0.2, 0.58, 0.28, 0.35])
    b = fig.add_axes([0.69, 0.58, 0.28, 0.35])
    c = fig.add_axes([0.2, 0.08, 0.28, 0.31])
    dd = fig.add_axes([0.69, 0.08, 0.28, 0.31])
    s7 = R3['S7']
    rows = [('Total cholesterol', 'tc'), ('HDL cholesterol', 'hdl'), ('Creatinine (calibrated)', 'scr_cal'), ('HbA1c', 'hba1c_log'), ('Urine creatinine', 'ucr')]
    eflm = eflm_cvi()
    for i, (nm, k) in enumerate(rows):
        v = s7[k]
        a.plot(v['cvw_total_pct'], i + 0.09, 'o', color=INK, ms=4.5)
        a.plot([v['cvi_derived_pct']['high_cva'], v['cvi_derived_pct']['low_cva']], [i + 0.27] * 2, color=TEAL, lw=2.2)
        a.plot(v['lacher_cvi'], i - 0.09, '^', color=ORANGE, ms=4.5)
        if k in eflm:
            a.plot(eflm[k], i - 0.27, 'D', color=GREY, ms=3.8)
        if k == 'scr_cal':
            a.plot(RRP['creatinine_scale']['cvw_raw_pct'], i + 0.09, 'o', color=INK, mfc='white', ms=4.5)
    a.set_yticks(range(len(rows)))
    a.set_yticklabels([r[0] for r in rows])
    a.invert_yaxis()
    a.set_xscale('log')
    a.set_xlim(1, 100)
    a.set_xticks([1, 2, 5, 10, 20, 50, 100])
    a.set_xticklabels(['1', '2', '5', '10', '20', '50', '100'])
    a.minorticks_off()
    a.set_xlabel('Within-person CV (%)')
    a.set_title('Analyte variation: replication', fontsize=8, loc='left')
    a.plot([], [], 'o', color=INK, label='Total CVw (calibrated creatinine)')
    a.plot([], [], color=TEAL, lw=2.2, label='CVw net of CVa (model-based)')
    a.plot([], [], '^', color=ORANGE, label='Lacher 2005 CVI')
    a.plot([], [], 'D', color=GREY, ms=3.8, label='EFLM meta-analysis CVI')
    a.plot([], [], 'o', color=INK, mfc='white', ms=4.5, label='CVw, uncalibrated creatinine')
    a.legend(fontsize=6.6, loc='upper center', bbox_to_anchor=(0.4, -0.17), ncol=2, columnspacing=1.6, handletextpad=0.5, handlelength=1.6, labelspacing=0.3)
    s8 = R3['S8']['base']
    ladder = [('Within visit: reading 1 vs 2', s8['a_v1_r1_vs_r2']), ('Within visit: reading 2 vs 3', s8['b_v1_r2_vs_r3']), ('Between days: single reading 2', s8['c_v1r2_vs_v2r2']), ('Between days: mean of readings 2–3', s8['d_between_visit_bp_swap']), ('All inputs remeasured (θ)', R3['primary'])]
    for i, (nm, e) in enumerate(ladder):
        m, lo, hi = pc(e)
        col = INK if i == 4 else GC['BP']
        b.plot([lo, hi], [i, i], color=col, lw=1.2)
        b.plot(m, i, 'o', color=col, ms=5)
    x2 = 100 * S2['bp_benchmark']['reading1_vs_reading2']['bp']['estimate']
    for x in (x2, 5.8):
        b.axvline(x, color=GREY, lw=0.7, ls=(0, (2, 2)))
    b.text(x2 - 0.15, -0.85, f'Study 2\nreading\nswap {x2:.1f}%', fontsize=6.5, ha='right', va='center', color='#5A6670')
    b.text(5.95, 5.3, 'Ogunwale (ARIC)\n5.8%', fontsize=6.5, ha='left', va='center', color='#5A6670')
    b.set_yticks(range(5))
    b.set_yticklabels([l[0] for l in ladder])
    b.invert_yaxis()
    b.set_xlim(0, 15)
    b.set_ylim(5.9, -1.5)
    b.set_xlabel('Changed base category, % (95% CI)')
    b.set_title('Blood pressure measurement ladder', fontsize=8, loc='left')
    s12 = R3['S12']
    th = 100 * R3['primary']['estimate']
    items = [('Observed θ', th, None), ('Simulated, observed\nvariation (A)', 100 * s12['A_base']['theta_sim'], None), ('Observed lab-only\nremeasurement', 100 * s12['D_sym_lab_observed'], None), ('Simulated analytical\nonly (B, low CVa)', 100 * s12['B_base_low_cva']['theta_sim'], None), ('Simulated analytical\nonly (B, high CVa)', 100 * s12['B_base_high_cva']['theta_sim'], None)]
    for i, (nm, v, _) in enumerate(items):
        c.barh(i, v, color=[INK, '#72808C', GC['LIPID'], '#D8C3D1', '#D8C3D1'][i], height=0.55)
        c.text(v + 0.2, i, f'{v:.1f}', va='center', fontsize=7)
    c.set_yticks(range(5))
    c.set_yticklabels([x[0] for x in items])
    c.invert_yaxis()
    c.set_xlim(0, 14)
    c.set_xlabel('Changed category (%)')
    c.set_title('Model-based benchmarks', fontsize=8, loc='left')
    sv = R3['S7']['sbp_variance']
    comps = [('Within visit,\nbetween readings', math.sqrt(sv['sigma2_within_reading'])), ('Between days', math.sqrt(sv['sigma2_between_day']))]
    for i, (nm, v) in enumerate(comps):
        dd.barh(i, v, color=GC['BP'] if i else '#A9BCD8', height=0.5)
        dd.text(v + 0.1, i, f'{v:.1f}', va='center', fontsize=7)
    dd.set_yticks(range(2))
    dd.set_yticklabels([x[0] for x in comps])
    dd.invert_yaxis()
    dd.set_xlim(0, 9)
    dd.set_xlabel('Systolic BP standard deviation (mmHg)')
    dd.set_title('SBP variation components', fontsize=8, loc='left')
    for ax, l in zip([a, b, c, dd], 'ABCD'):
        lab(ax, l, -0.62)
    return (fig, [{'id': 'A', 'question': 'Do analyte CVs replicate Lacher 2005?', 'units': 'percent', 'sample_unit': 'participant', 'transform': 'log scale', 'estimate': 'total CVw = sqrt(exp(s^2)-1); derived CVI; Lacher 2005 and EFLM CVI (frozen transcriptions)', 'interval': 'none; point estimates with CVa bracket', 'n': 759}, {'id': 'B', 'question': 'BP measurement ladder', 'units': 'percent', 'sample_unit': 'participant', 'transform': 'none', 'estimate': 'proportion changing category', 'interval': 'CR1 logit-t 95% CI', 'n': 759}, {'id': 'C', 'question': 'Model-based simulation benchmarks', 'units': 'percent', 'sample_unit': 'participant', 'transform': 'none', 'estimate': 'simulated or observed proportions', 'interval': 'none; model-based point estimates', 'n': 759}, {'id': 'D', 'question': 'SBP variance components', 'units': 'mmHg', 'sample_unit': 'participant', 'transform': 'square root of variance components', 'estimate': 'within-reading and between-day SD of a single reading', 'interval': 'none; point estimates', 'n': 759}])

def fig_s1():
    fig = plt.figure(figsize=mm(180, 200))
    a = fig.add_axes([0.38, 0.05, 0.25, 0.9])
    b = fig.add_axes([0.77, 0.05, 0.2, 0.9])
    T = [('Primary (n = 759)', R3['primary'])]
    lab_ = {'SA1a_exclude_out_of_domain': 'SA1a Exclude out-of-domain inputs', 'SA1b_clamp_sbp_90_180': 'SA1b Clamp SBP at 90–180', 'SA1c_in_domain_v1': 'SA1c In domain at visit 1', 'SA1d_paper1_exclusion': 'SA1d Exclude, SBP 90–180 (Study 2 pilot)', 'SA2a_single_reading2': 'SA2a Single reading (reading 2)', 'SA2b_mean_available_1_3': 'SA2b Mean of available readings 1–3', 'SA2c_exclude_different_arm': 'SA2c Exclude different arm', 'SA3a_confirmed_diabetes_fixed': 'SA3a Confirmed diabetes held fixed', 'SA3b_exclude_hba1c_method_switch': 'SA3b Exclude HbA1c method switch', 'SA3c_gate_rule': 'SA3c Count-gate diabetes rule', 'SA4a_raw_creatinine': 'SA4a Uncalibrated creatinine', 'SA4b_coresh_minus_0.23': 'SA4b Coresh correction (−0.23 mg/dL)', 'SA5_referral_guard': 'SA5 Exclude referral Level I/II or unknown', 'SA6_visit_specific_age': 'SA6 Visit-specific age', 'SA7_visit_specific_smoking': 'SA7 Visit-specific smoking', 'SA8_exclude_flagged': 'SA8 Exclude implausible-change flags', 'SA9_exclude_rose': 'SA9 Exclude Rose exertional chest pain', 'SA9_exclude_angina_users': 'SA9 Exclude antianginal users', 'SA9_exclude_preg_unknown': 'SA9 Exclude unknown pregnancy status', 'SA9_exclude_all_three': 'SA9 All three exclusions', 'SA11_lipid_lowering_proxy': 'SA11 Lipid-lowering proxy for statin', 'SA12_same_examiner': 'SA12 Same physician examiner only', 'SA13_risk_rounded_0.1': 'SA13 Risk rounded to 0.1 pp'}
    for k, l in lab_.items():
        e = R4[k]['theta']
        T.append((f"{l} (n = {e['n']})", e))
        if k == 'SA9_exclude_all_three':
            w = R4['SA10_completeness_weighting']
            T.append(('SA10 Completeness weighting (n = 759)', {'estimate': w['theta_weighted'], 'lower': w['ci'][0], 'upper': w['ci'][1]}))
    s15 = R4['SA15_uncertainty']
    th = R3['primary']['estimate']
    for l, lo, hi in [('SA15a Wilson interval', *s15['wilson']), ('SA15b Stratified Taylor (df 49)', s15['stratified_taylor']['lower'], s15['stratified_taylor']['upper']), ('SA15c Cluster-bootstrap percentile', *s15['cluster_bootstrap_percentile']), ('SA15d CR1 on pseudo-stratum (G = 49)', s15['cr1_pseudo_stratum']['lower'], s15['cr1_pseudo_stratum']['upper'])]:
        T.append((l, {'estimate': th, 'lower': lo, 'upper': hi}))
    for i, (l, e) in enumerate(T):
        m, lo, hi = pc(e)
        col = INK if i == 0 else ORANGE if R4.get('SA4a_raw_creatinine', {}).get('theta') is e else GC['BP']
        a.plot([lo, hi], [i, i], color=col, lw=1)
        a.plot(m, i, 'o' if i else 'D', color=col, ms=3.6)
    a.axvline(100 * th, color=GREY, lw=0.7, ls=(0, (2, 2)))
    a.axvspan(100 * th - 1.5, 100 * th + 1.5, color=GRID, alpha=0.6, lw=0)
    a.set_yticks(range(len(T)))
    a.set_yticklabels([t[0] for t in T], fontsize=6.6)
    a.set_ylim(len(T) - 0.4, -0.6)
    a.set_xlim(5, 17)
    a.set_xlabel('Changed category, % (95% CI)')
    a.set_title('A  Base PREVENT-ASCVD ($\\theta$)', fontsize=8, loc='left')
    U = [('Primary (n = 725)', R3['K1']['theta_U'])]
    for k in ['SA1a_exclude_out_of_domain', 'SA1b_clamp_sbp_90_180', 'SA1c_in_domain_v1', 'SA1d_paper1_exclusion', 'SA3a_confirmed_diabetes_fixed', 'SA3b_exclude_hba1c_method_switch', 'SA3c_gate_rule', 'SA4a_raw_creatinine', 'SA4b_coresh_minus_0.23', 'SA5_referral_guard', 'SA14a_exclude_albumin_bdl', 'SA14b_keep_ucr_bdl']:
        e = R4[k]['theta_U']
        U.append((k.split('_')[0] + f" (n = {e['n']})", e))
    for i, (l, e) in enumerate(U):
        m, lo, hi = pc(e)
        col = INK if i == 0 else GC['URINE']
        b.plot([lo, hi], [i, i], color=col, lw=1)
        b.plot(m, i, 'o' if i else 'D', color=col, ms=3.6)
    thU = 100 * R3['K1']['theta_U']['estimate']
    b.axvline(thU, color=GREY, lw=0.7, ls=(0, (2, 2)))
    b.axvspan(thU - 1.5, thU + 1.5, color=GRID, alpha=0.6, lw=0)
    b.set_yticks(range(len(U)))
    b.set_yticklabels([u[0] for u in U], fontsize=6.6)
    b.set_ylim(len(T) - 0.4, -0.6)
    b.set_xlim(7, 17)
    b.set_xlabel('Changed category, % (95% CI)')
    b.set_title('B  UACR-enhanced ($\\theta_U$)', fontsize=8, loc='left')
    return (fig, [{'id': 'A', 'question': 'Sensitivity of theta to SA1-SA15', 'units': 'percent', 'sample_unit': 'participant', 'transform': 'none', 'estimate': 'proportion changing base category', 'interval': 'CR1 logit-t 95% CI (SA10: cluster bootstrap with refitted weights; SA15: alternative intervals for the primary estimate)', 'n': 759}, {'id': 'B', 'question': 'Sensitivity of theta_U', 'units': 'percent', 'sample_unit': 'participant', 'transform': 'none', 'estimate': 'proportion changing UACR-model category', 'interval': 'CR1 logit-t 95% CI', 'n': 725}])
S2 = json.loads((ROOT / 'studies/nhanes3_testretest/manuscript_evidence_v1/study2_evidence.json').read_text(encoding='utf-8'))
P2 = pd.read_csv(ROOT / 'studies/uacr_prevent/final_and_external_v1/primary_participant_results.csv')

def fig6():
    fig = plt.figure(figsize=mm(180, 195))
    a = fig.add_axes([0.08, 0.64, 0.38, 0.3])
    b = fig.add_axes([0.62, 0.665, 0.3, 0.245])
    c = fig.add_axes([0.42, 0.06, 0.55, 0.48])
    E1 = json.loads((ROOT / 'studies/nhanes3_testretest/analysis_v1_posthoc/ext2001_labs.json').read_text(encoding='utf-8'))
    pa = S2['primary']['any_category_change']
    assert len(P2) == pa['n'] and int((P2.category_first != P2.category_second).sum()) == pa['events']
    scatter_panel(a, P2.risk_first.to_numpy(), P2.risk_second.to_numpy(), f"Study 2: NHANES 2009–2010 (n = {pa['n']:,})\nchanged category: {pa['events']} (weighted {100 * pa['estimate']:.1f}%)", xl='10-year ASCVD risk, random urine (%)', yl='Risk, first-morning urine (%)', second='first-morning urine')
    trans_panel(b, S2['transition_counts'], 'Random × first-morning urine (counts)', xl='First-morning urine category (%)', yl='Random urine category (%)')
    s8 = R3['S8']['base']
    s8u = R3['S8']['uacr']
    b12 = S2['bp_benchmark']['reading1_vs_reading2']
    ex = S2['external']['any_category_change']
    pc2 = lambda e: (100 * e['estimate'], 100 * e['lower95'], 100 * e['upper95'])
    rows = [('Within one visit', None), ('SBP reading 1 vs 2, PREVENT-UACR', ('B', pc2(b12['bp']))), ('SBP reading 1 vs 2, base PREVENT', ('A', pc(s8['a_v1_r1_vs_r2']))), ('SBP reading 2 vs 3, base PREVENT', ('A', pc(s8['b_v1_r2_vs_r3']))), ('Urine specimen', None), ('Random → first-morning UACR', ('B', pc2(pa))), ('Random → random UACR, 2001–2002 repeat', ('C', pc2(ex))), ('Random → random UACR (urine only, UACR model)', ('A', (100 * R3['K2_uacr']['single_swap']['URINE']['ref1'], None, None))), ('Another day, one input group', None), ('Single SBP reading', ('A', pc(s8['c_v1r2_vs_v2r2']))), ('Mean of SBP readings 2–3', ('A', pc(s8['d_between_visit_bp_swap']))), ('Lipids only', ('A', (100 * R3['K2_base']['single_swap']['LIPID']['ref1'], None, None))), ('Another day, laboratory inputs (BP fixed)', None), ('Laboratory inputs, base PREVENT', ('A', pc(E1['study1_comparator']['base_lab_only_ref1']))), ('Laboratory inputs, base PREVENT, 2001–2002', ('C', pc(E1['base']['theta_lab']))), ('Laboratory inputs + UACR, PREVENT-UACR', ('A', pc(E1['study1_comparator']['uacr_lab_urine_ref1']))), ('Laboratory inputs + UACR, PREVENT-UACR, 2001–2002', ('C', pc(E1['uacr']['theta_lab']))), ('Another day, all inputs', None), ('All routinely measured inputs, base PREVENT', ('A', pc(R3['primary']))), ('All inputs including UACR, PREVENT-UACR', ('A', pc(R3['K1']['theta_U'])))]
    mk = {'A': ('o', '#3766A0', 'Study 1 (NHANES III)'), 'B': ('s', ORANGE, 'Study 2 (NHANES 2009–2010)'), 'C': ('^', '#72808C', 'NHANES 2001–2002 repeat')}
    labels = []
    for i, (lbl, v) in enumerate(rows):
        if v is None:
            c.text(-0.02, i, lbl, transform=c.get_yaxis_transform(), ha='right', va='center', fontsize=7, fontweight='bold', color=INK)
            labels.append('')
            continue
        labels.append(lbl)
        key, (m, lo, hi) = v
        mark, col, _ = mk[key]
        if lo is not None:
            c.plot([lo, hi], [i, i], color=col, lw=1.2)
        c.plot(m, i, mark, color=col, ms=5, mfc=col if lo is not None else 'white')
        c.text((hi if lo is not None else m) + 0.3, i, f'{m:.1f}', va='center', ha='left', fontsize=6.8, color=INK)
    c.set_yticks(range(len(rows)))
    c.set_yticklabels(labels, fontsize=7)
    c.set_ylim(len(rows) - 0.5, -0.6)
    c.set_xlim(0, 16)
    c.tick_params(axis='y', length=0)
    for k_, (mark, col, l) in mk.items():
        c.plot([], [], mark, color=col, ms=5, label=l)
    c.plot([], [], 'o', color=INK, mfc='white', ms=5, label='No interval (attribution)')
    c.legend(fontsize=6.6, loc='upper right', bbox_to_anchor=(1.0, 1.0), framealpha=1, facecolor='white')
    c.set_xlabel('Changed PREVENT-ASCVD category, % (95% CI)')
    c.set_title('Category change by source of measurement variation', fontsize=8, loc='left')
    lab(a, 'A', -0.16)
    lab(b, 'B', -0.22)
    lab(c, 'C', -0.62)
    return (fig, [{'id': 'A', 'question': 'Study 2 paired risks', 'units': 'percent risk', 'sample_unit': 'participant', 'transform': 'log10 axes', 'estimate': 'raw paired observations', 'interval': 'none', 'n': pa['n']}, {'id': 'B', 'question': 'Study 2 transitions', 'units': 'counts', 'sample_unit': 'participant', 'transform': 'none', 'estimate': 'unweighted counts', 'interval': 'none', 'n': pa['n']}, {'id': 'C', 'question': 'Category change by source of variation across studies', 'units': 'percent', 'sample_unit': 'participant', 'transform': 'none', 'estimate': 'proportions changing category (NHANES 2009-2010 survey-weighted; others unweighted)', 'interval': 'Study 1 CR1 logit-t; Study 2 Taylor-logit; 2001-2002 cluster-adjusted', 'n': 759}])
RS = json.loads((ROOT / 'studies/nhanes3_testretest/analysis_v1_posthoc/repeat_strategies.json').read_text(encoding='utf-8'))

def fig7():
    fig = plt.figure(figsize=mm(180, 92))
    a = fig.add_axes([0.25, 0.16, 0.25, 0.74])
    b = fig.add_axes([0.62, 0.16, 0.35, 0.74])
    st = RS['strategies']
    rows = [('single', 'One measurement'), ('bp2', 'BP on 2 days'), ('bp3', 'BP on 3 days'), ('lab2', 'Laboratory inputs on 2 days'), ('all2', 'All inputs on 2 days'), ('all3', 'All inputs on 3 days'), ('t10_bp', 'Second BP day if in ≥10% zone'), ('t25_all', 'Second day, all inputs,\nif in ≥25% zone'), ('t10_all', 'Second day, all inputs,\nif in ≥10% zone')]
    what = {'single': 'none', 'bp2': 'bp', 'bp3': 'bp', 'lab2': 'lab', 'all2': 'all', 'all3': 'all', 't10_bp': 'bp', 't25_all': 'all', 't10_all': 'all'}
    MK = {'none': 'D', 'bp': '^', 'lab': 's', 'all': 'o'}
    COL = {'none': INK, 'bp': GC['BP'], 'lab': GC['LIPID'], 'all': TEAL}
    targeted = lambda k: k.startswith('t')
    for i, (k, l) in enumerate(rows):
        v = 100 * st[k]['theta_sim']
        col = COL[what[k]]
        a.barh(i, v, color=col if not targeted(k) else 'white', edgecolor=col, height=0.6, lw=1.2)
        a.text(v + 0.25, i, f'{v:.1f}', va='center', fontsize=6.8, color=INK)
    a.axvline(100 * RS['observed_theta'], color=GREY, lw=0.9, ls=(0, (2, 2)))
    a.text(100 * RS['observed_theta'] + 0.2, len(rows) - 0.35, f"observed θ {100 * RS['observed_theta']:.1f}%", fontsize=6.5, color='#5A6670')
    a.set_yticks(range(len(rows)))
    a.set_yticklabels([l for _, l in rows], fontsize=6.8)
    a.invert_yaxis()
    a.set_xlim(0, 14)
    a.set_ylim(len(rows) - 0.2, -0.6)
    a.set_xlabel('Simulated category disagreement (%)')
    a.set_title('Disagreement between two assessments', fontsize=8, loc='left')
    for k, _ in rows[1:]:
        x = st[k]['extra_bp_days_per_person'] if what[k] != 'lab' else st[k]['extra_lab_draws_per_person']
        y = 100 * (1 - st[k]['relative_to_single'])
        lo, hi = (100 * (1 - st[k]['relative_ci'][1]), 100 * (1 - st[k]['relative_ci'][0]))
        col = COL[what[k]]
        b.plot([x, x], [lo, hi], color=col, lw=1.1)
        b.plot(x, y, MK[what[k]], color=col, ms=6, mfc='white' if targeted(k) else col, mew=1.3)
    lbl = {'bp2': ('BP, 2 days', (6, -4)), 'bp3': ('BP, 3 days', (-7, -3)), 'lab2': ('Laboratory, 2 days', (6, 0)), 'all2': ('All, 2 days', (6, 2)), 'all3': ('All, 3 days', (-6, 5)), 't10_bp': ('BP if ≥10% zone', (5, -13)), 't25_all': ('All if ≥25% zone (17%)', (6, -12)), 't10_all': ('All if ≥10% zone\n(37%)', (0.001, 19))}
    for k, (t, off) in lbl.items():
        x = st[k]['extra_bp_days_per_person'] if what[k] != 'lab' else st[k]['extra_lab_draws_per_person']
        b.annotate(t, (x, 100 * (1 - st[k]['relative_to_single'])), xytext=off, textcoords='offset points', fontsize=6.6, color=INK, ha='right' if off[0] < 0 else 'center' if off[0] < 1 else 'left', va='center')
    b.set_xlim(-0.05, 2.25)
    b.set_ylim(0, 50)
    b.axhline(0, color='#72808C', lw=0.6)
    b.set_xlabel('Additional measurement days per person')
    b.set_ylabel('Reduction in disagreement (%)')
    b.set_title('Reduction versus extra measurement', fontsize=8, loc='left')
    for k_, l_ in [('all', 'All inputs repeated'), ('bp', 'BP only'), ('lab', 'Laboratory inputs only')]:
        b.plot([], [], MK[k_], color=COL[k_], ms=5, label=l_)
    b.plot([], [], 'o', color=INK, mfc='white', ms=5, label='Targeted (open symbols)')
    b.legend(fontsize=6.5, loc='lower right', ncol=2, columnspacing=1.0, handletextpad=0.4, framealpha=1, facecolor='white')
    lab(a, 'A', -0.93)
    lab(b, 'B', -0.12)
    return (fig, [{'id': 'A', 'question': 'Simulated test-retest category disagreement under repeat-measurement strategies', 'units': 'percent', 'sample_unit': 'participant', 'transform': 'none', 'estimate': 'model-based simulation (1,000 replicates per participant; extends prespecified model A)', 'interval': 'none', 'n': 759}, {'id': 'B', 'question': 'Relative reduction versus additional measurement days', 'units': 'percent', 'sample_unit': 'participant', 'transform': '1 - ratio to single', 'estimate': 'ratio of mean simulated disagreement', 'interval': 'cluster bootstrap of participants (4,000), simulation SDs fixed', 'n': 759}])
LKP = ROOT / 'studies/nhanes3_testretest/analysis_v1_posthoc'
SRC_CI = SRC + [LKP / 'reclassification_lookup.json', LKP / 'reclassification_grid.csv']

def fig_ci():
    import build_nh3_manuscript as M
    Wd, Ht = (180, 104)
    fig = plt.figure(figsize=mm(Wd, Ht))
    bg = fig.add_axes([0, 0, 1, 1])
    bg.set_xlim(0, Wd)
    bg.set_ylim(0, Ht)
    bg.axis('off')

    def card(x, y, w, h, ec, fc='white'):
        bg.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0,rounding_size=2.2', fc=fc, ec=ec, lw=0.9))

    def T(x, y, s, **k):
        bg.text(x, y, s, **{'fontsize': 7, 'color': INK, 'va': 'top', 'linespacing': 1.25, **k})
    pr, pa, Psi = (M.pr, M.pa, M.K2['Psi'])
    T(4, Ht - 3.5, 'How stable is a PREVENT risk category when the same person is measured again?', fontsize=10.5, fontweight='bold')
    T(4, Ht - 9.6, 'PREVENT 10-year ASCVD risk categories <3%, 3–<5%, 5–<10% and ≥10%; two NHANES repeat-measurement studies', fontsize=7.2, color='#4A5966')
    card(4, 38, 55, 50, TEAL)
    T(7, 85, 'Study 1 · NHANES III, 1988–1994', fontsize=7.6, fontweight='bold', color=TEAL)
    T(7, 80.3, f"{pr['n']} adults; blood pressure, cholesterol,\ncreatinine and HbA1c measured again\na median {M.iv['median']:.0f} days later (range {M.iv['min']:.0f}–{M.iv['max']:.0f})", fontsize=6.9)
    T(7, 67.2, f"{M.P(pr['estimate'])}%", fontsize=18, fontweight='bold', color=TEAL)
    T(32, 67.6, f"changed category\n(base PREVENT)\n95% CI {M.cin(pr)}\n{pr['up']} up, {pr['down']} down", fontsize=6.9)
    T(7, 53.6, 'Contribution, percentage points', fontsize=6.6, color='#4A5966')
    x0, sc = (7, 49 / (100 * pr['estimate']))
    for g in ['BP', 'LIPID', 'GLYC', 'KIDNEY']:
        w = Psi[g]['pp'] * sc
        bg.add_patch(plt.Rectangle((x0, 47.2), w, 3.0, fc=GC[g], ec='white', lw=0.5))
        x0 += w
    for g, nm, lx, ly in [('BP', 'Blood pressure', 7, 45.6), ('LIPID', 'Lipids', 33, 45.6), ('GLYC', 'HbA1c', 7, 41.6), ('KIDNEY', 'Creatinine', 33, 41.6)]:
        bg.add_patch(plt.Rectangle((lx, ly - 2.3), 2, 2, fc=GC[g], ec='none'))
        T(lx + 2.8, ly + 0.1, f"{nm} {M.F(Psi[g]['pp'])}", fontsize=6.6)
    card(4, 5, 55, 30, GC['URINE'])
    T(7, 32, 'Study 2 · NHANES 2009–2010', fontsize=7.6, fontweight='bold', color=GC['URINE'])
    T(7, 27.3, f"{pa['n']:,} adults; only the urine specimen\nchanged (random → first-morning UACR)", fontsize=6.9)
    T(7, 18.6, f"{M.P(pa['estimate'])}%", fontsize=18, fontweight='bold', color=GC['URINE'])
    T(27, 17.8, f'changed category\n(PREVENT-UACR)\n95% CI {M.cin(pa)}', fontsize=6.9)
    g = pd.read_csv(LKP / 'reclassification_grid.csv')
    g = g[(g.risk_pct >= 1) & (g.risk_pct <= 25)].reset_index(drop=True)
    T(66, 85, 'Modelled chance of a different category\non remeasurement (Study 1)', fontsize=7.6, fontweight='bold')
    ax = fig.add_axes([73 / Wd, 17 / Ht, 55 / Wd, 58 / Ht])
    for t, (a, b) in M.ZB25.items():
        ax.axvspan(a, b, color='#F6DCC8', lw=0, zorder=0)
    ax.fill_between(g.risk_pct, 100 * g.base_lower, 100 * g.base_upper, color=TEAL, alpha=0.18, lw=0, zorder=1)
    ax.plot(g.risk_pct, 100 * g.base_p, color=TEAL, lw=1.6, zorder=2)
    for t, v in zip([3, 5, 10], M.pi0):
        ax.axvline(t, color='#72808C', lw=0.7, ls=(0, (3, 2)), zorder=1)
        up = 100 * g.base_upper[(g.risk_pct - t).abs().idxmin()]
        ax.text(t, up + 2.5, f'{100 * v:.0f}%', ha='center', fontsize=7, fontweight='bold', color=INK)
    ax.axhline(25, color=ORANGE, lw=0.7, ls=(0, (1, 1.5)))
    ax.set_xscale('log')
    ax.set_xlim(1, 25)
    ax.set_ylim(0, 80)
    ax.set_xticks([1, 2, 3, 5, 10, 20])
    ax.set_xticklabels(['1', '2', '3', '5', '10', '20'])
    ax.minorticks_off()
    ax.set_yticks([0, 25, 50])
    ax.set_yticklabels(['0%', '25%', '50%'])
    ax.set_xlabel('First PREVENT-ASCVD 10-year risk (%)', fontsize=7)
    bg.add_patch(plt.Rectangle((73, 4.6), 3, 2.4, fc='#F6DCC8', ec='none'))
    T(77, 7.6, 'Chance ≥25% (instability zone)', fontsize=6.6)
    bg.add_patch(plt.Rectangle((117, 4.6), 3, 2.4, fc=TEAL, alpha=0.18, ec='none'))
    T(121, 7.6, '95% CI', fontsize=6.6)
    card(134, 5, 42, 83, ORANGE)
    T(137, 85, 'What it implies', fontsize=7.6, fontweight='bold', color=ORANGE)
    lk = {r['risk_pct']: r['base']['p'] for r in M.LK['lookup']}
    T(137, 79.5, f'Near 3%, 5% or 10%, one\nset of measurements gives\nan uncertain category; far\nfrom them, the modelled\nchance was {M.P(lk[1])}% at 1%\nrisk and {M.P(lk[15])}% at 15%.', fontsize=6.9)
    x3 = M.x3
    T(137, 59, 'Blood pressure treatment\n(post hoc)', fontsize=6.9, fontweight='bold', color=ORANGE)
    T(137, 52, f"{M.P(x3['theta_rx']['estimate'])}%", fontsize=14, fontweight='bold', color=ORANGE)
    T(137, 46.3, f"of {x3['n_untreated_eligible']} untreated adults had\na different 2025-guideline\nrecommendation to start\nmedication at the second\nvisit (95% CI {M.cin(x3['theta_rx'])})", fontsize=6.9)
    T(137, 29.5, f"Blood pressure accounted\nfor {M.F(100 * x3['shapley']['BP']['estimate'])} of the {M.P(x3['theta_rx']['estimate'])}\npercentage points; {x3['flips_with_stage_change']} of {x3['flips_total']}\nchanges involved a change\nin blood pressure stage", fontsize=6.9)
    return (fig, [{'id': 'centre', 'question': 'Modelled probability of a different base PREVENT-ASCVD category on remeasurement by first-visit risk', 'units': 'percent', 'sample_unit': 'participant', 'transform': 'log x axis', 'estimate': 'prespecified pi_clin model (K3) evaluated on a risk grid (post hoc display)', 'interval': 'pointwise 95% cluster bootstrap (4,000)', 'n': 759}])

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    titles = {1: 'Study flow', 2: 'Paired PREVENT risks and category transitions', 3: 'Input attribution', 4: 'Threshold-distance curves', 5: 'Analyte variation and benchmarks'}
    titles['s1'] = 'Sensitivity analyses'
    titles[6] = 'Urine specimen choice and sources of category change'
    titles[7] = 'Repeat-measurement strategies (model-based)'
    titles['ci'] = 'Central illustration'
    only = sys.argv[1:]
    for i, f in list(enumerate([fig1, fig2, fig3, fig4, fig5, fig6, fig7], 1)) + [('s1', fig_s1), ('ci', fig_ci)]:
        if only and str(i) not in only:
            continue
        with publication_style(font_size=7.5):
            fig, panels = f()
            export_figure(fig, OUT / f'figure_{i}', sources=SRC_CI if i == 'ci' else SRC, panels=panels, title=titles[i], dpi=300, overwrite=True)
        plt.close(fig)
        print('figure', i, 'ok')
if __name__ == '__main__':
    main()
