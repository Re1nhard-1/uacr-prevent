from pathlib import Path
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import json
from scipy.special import expit, logit
import nh3_figures as NF
import build_nh3_manuscript as M
import nh3_fig_e11 as E11F
from figure_tools import publication_style, export_figure
ROOT = NF.ROOT
R3, P = (NF.R3, NF.P)
OUT = ROOT / 'studies/nhanes3_testretest/figures_v7'
INK, TEAL, ORANGE, GREY, GRID, GC, GL, BLUE = (NF.INK, NF.TEAL, NF.ORANGE, NF.GREY, NF.GRID, NF.GC, NF.GL, '#3766A0')
mm, lab = (NF.mm, NF.lab)

def pi_risk_panel(ax):
    g = pd.read_csv(NF.LKP / 'reclassification_grid.csv')
    g = g[(g.risk_pct >= 1) & (g.risk_pct <= 25)].reset_index(drop=True)
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
    ax.set_xlabel('First-visit PREVENT-ASCVD 10-year risk (%)')
    ax.set_ylabel('Probability of a different category\nat the second visit (%)')
    ax.set_title('Conditional on first-visit risk ($\\pi_{clin}$), risk scale', fontsize=8, loc='left')
    ax.text(0.985, 0.985, 'Shaded: zones with\nprobability ≥25%', transform=ax.transAxes, ha='right', va='top', fontsize=6.5, color=INK, bbox=dict(boxstyle='round,pad=0.3', fc='white', ec=GRID, lw=0.6))

def shapley_panel(ax, key, groups, title):
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

def lp_partition_panel(ax):
    groups = ['BP', 'LIPID', 'KIDNEY', 'GLYC', 'URINE']
    for i, g in enumerate(['BP', 'LIPID', 'KIDNEY', 'GLYC']):
        v = R3['K2_base']['continuous'][g]
        ax.plot([100 * v['pi']], [i - 0.12], 'o', color=GC[g], ms=5)
        ax.plot([100 * x for x in v['pi_ci']], [i - 0.12, i - 0.12], color=GC[g], lw=1)
    for i, g in enumerate(groups):
        v = R3['K2_uacr']['continuous'][g]
        ax.plot([100 * v['pi']], [i + 0.12], 's', color=GC[g], ms=5, mfc='white')
        ax.plot([100 * x for x in v['pi_ci']], [i + 0.12, i + 0.12], color=GC[g], lw=1)
    ax.set_yticks(range(5))
    ax.set_yticklabels([GL[g] for g in groups])
    ax.invert_yaxis()
    ax.axvline(0, color='#72808C', lw=0.6)
    ax.set_xlim(-5, 65)
    ax.plot([], [], 'o', color=INK, label='Base model')
    ax.plot([], [], 's', color=INK, mfc='white', label='UACR model')
    ax.legend(fontsize=6.5, loc='lower right')
    ax.set_xlabel('Share of mean squared change in the linear predictor (%)')
    ax.set_title('Continuous scale: partition of the change', fontsize=8, loc='left')

def benchmarks_panel(ax):
    s12 = R3['S12']
    th = 100 * R3['primary']['estimate']
    items = [('Observed θ', th), ('Simulated, observed\nvariation (model A)', 100 * s12['A_base']['theta_sim']), ('Observed laboratory-only\nremeasurement', 100 * s12['D_sym_lab_observed']), ('Simulated analytical\nerror only (B, low CVa)', 100 * s12['B_base_low_cva']['theta_sim']), ('Simulated analytical\nerror only (B, high CVa)', 100 * s12['B_base_high_cva']['theta_sim'])]
    for i, (nm, v) in enumerate(items):
        ax.barh(i, v, color=[INK, '#72808C', GC['LIPID'], '#D8C3D1', '#D8C3D1'][i], height=0.55)
        ax.text(v + 0.2, i, f'{v:.1f}', va='center', fontsize=7)
    ax.set_yticks(range(5))
    ax.set_yticklabels([x[0] for x in items])
    ax.invert_yaxis()
    ax.set_xlim(0, 14)
    ax.set_xlabel('Changed base category (%)')
    ax.set_title('Model-based benchmarks', fontsize=8, loc='left')

def e11_loss_panel(ax, ylabels=True):
    E, LAB, y = (E11F.E, E11F.LAB, E11F.y)
    for j, (key, col, mk, lb) in enumerate([('loss_fraction_single', INK, 'o', 'Single visit'), ('loss_fraction_two_day', BLUE, 's', 'Two-day average'), ('loss_fraction_targeted', TEAL, 'D', 'Second day only in the ≥25% zone')]):
        dy = (1 - j) * 0.22
        for yi, e in zip(y, E):
            ci = e['ci_uniform_start'][key]
            ax.plot([100 * ci[0], 100 * ci[1]], [yi + dy] * 2, color=col, lw=1, alpha=0.8)
            ax.plot(100 * e[key], yi + dy, mk, color=col, ms=4.5, zorder=5)
        ax.plot([], [], mk, color=col, ms=4.5, label=lb)
    ax.set_yticks(y)
    ax.set_yticklabels(LAB if ylabels else [], fontsize=7)
    ax.set_xlim(0, 15)
    ax.set_ylim(-0.6, len(y) - 0.3)
    ax.set_xlabel('Expected net benefit forfeited (% of achievable)')
    ax.legend(fontsize=6.5, loc='upper right', frameon=True)
    ax.set_title('Net benefit forfeited, by measurement strategy', fontsize=8, loc='left')

def e11_regret_panel(ax):
    E, y = (E11F.E, E11F.y)
    for yi, e in zip(y, E):
        ci = e['ci_uniform_start']['regret_single_per_1000']
        ax.barh(yi, e['regret_single_per_1000'], color=INK, height=0.5)
        ax.plot(ci, [yi] * 2, color=ORANGE, lw=1.2)
        ax.text(max(ci[1], e['regret_single_per_1000']) + 0.03, yi, f"of {e['achievable_nb_per_1000']:.1f}", va='center', fontsize=6.5, color='#5A6670')
    ax.set_yticks(y)
    ax.set_yticklabels([])
    ax.set_xlim(0, 1.25)
    ax.set_xlabel('Net true positives lost per 1,000\n(single visit; achievable printed beside bar)')
    ax.set_title('Absolute scale', fontsize=8, loc='left')

def e11_sens_panel(ax):
    E, y = (E11F.E, E11F.y)
    short = ['3%', '5%', '7.5% CVD', '30-y 10%', '10%']
    for yi, e in zip(y, E):
        sb = e.get('post_audit_sens_sigma_by_age', {}).get('loss_fraction_single')
        vi = e['post_audit_sens_variance_inflation']
        ax.plot(100 * e['loss_fraction_single'], yi, 'o', color=INK, ms=4, zorder=5)
        if sb is not None:
            ax.plot(100 * sb, yi, '^', color=BLUE, ms=4, zorder=5)
        ax.plot(100 * vi['1.5']['loss_fraction_single'], yi, '>', color=ORANGE, ms=4, zorder=5)
        ax.plot(100 * vi['2.0']['loss_fraction_single'], yi, '>', mfc='white', color=ORANGE, ms=4, zorder=5)
    for mk, col, mfc, lb in [('o', INK, INK, 'Main'), ('^', BLUE, BLUE, 'Error SD by age'), ('>', ORANGE, ORANGE, 'Variance ×1.5'), ('>', ORANGE, 'white', 'Variance ×2')]:
        ax.plot([], [], mk, color=col, mfc=mfc, ms=4, label=lb)
    ax.set_yticks(y)
    ax.set_yticklabels(short, fontsize=7)
    ax.set_xlim(0, 20)
    ax.set_xlabel('Forfeited at a single visit (%)')
    ax.legend(fontsize=6.5, loc='upper right', frameon=True, handletextpad=0.3, borderpad=0.3)
    ax.set_title('Sensitivity to the error model', fontsize=8, loc='left')

def repeat_panel(ax):
    st = NF.RS['strategies']
    rows = ['bp2', 'bp3', 'lab2', 'all2', 'all3', 't10_bp', 't25_all', 't10_all']
    what = {'bp2': 'bp', 'bp3': 'bp', 'lab2': 'lab', 'all2': 'all', 'all3': 'all', 't10_bp': 'bp', 't25_all': 'all', 't10_all': 'all'}
    MK = {'bp': '^', 'lab': 's', 'all': 'o'}
    COL = {'bp': GC['BP'], 'lab': GC['LIPID'], 'all': TEAL}
    xs = lambda k: st[k]['extra_bp_days_per_person'] if what[k] != 'lab' else st[k]['extra_lab_draws_per_person']
    for k in rows:
        x = xs(k)
        yv = 100 * (1 - st[k]['relative_to_single'])
        lo, hi = (100 * (1 - st[k]['relative_ci'][1]), 100 * (1 - st[k]['relative_ci'][0]))
        col = COL[what[k]]
        ax.plot([x, x], [lo, hi], color=col, lw=1.1)
        ax.plot(x, yv, MK[what[k]], color=col, ms=6, mfc='white' if k.startswith('t') else col, mew=1.3)
    lbl = {'bp2': ('BP, 2 days', (6, 0)), 'bp3': ('BP, 3 days', (-7, -3)), 'lab2': ('Laboratory, 2 days', (6, 0)), 'all2': ('All, 2 days', (6, 2)), 'all3': ('All, 3 days', (-6, 5)), 't10_bp': ('BP if ≥10% zone', (0.001, -13)), 't25_all': ('All if ≥25% zone (17%)', (6, -12)), 't10_all': ('All if ≥10% zone\n(37%)', (0.001, 19))}
    for k, (t, off) in lbl.items():
        ax.annotate(t, (xs(k), 100 * (1 - st[k]['relative_to_single'])), xytext=off, textcoords='offset points', fontsize=6.6, color=INK, ha='right' if off[0] < 0 else 'center' if off[0] < 1 else 'left', va='center')
    ax.set_xlim(-0.05, 2.25)
    ax.set_ylim(0, 57)
    ax.axhline(0, color='#72808C', lw=0.6)
    ax.set_xlabel('Additional measurement days per person')
    ax.set_ylabel('Reduction in category disagreement (%)')
    ax.set_title('Repeat strategies: disagreement vs workload', fontsize=8, loc='left')
    for k_, l_ in [('all', 'All inputs repeated'), ('bp', 'BP only'), ('lab', 'Laboratory inputs only')]:
        ax.plot([], [], MK[k_], color=COL[k_], ms=5, label=l_)
    ax.plot([], [], 'o', color=INK, mfc='white', ms=5, label='Targeted (open symbols)')
    ax.legend(fontsize=6.5, loc='upper left', ncol=1, handletextpad=0.4, borderpad=0.3, labelspacing=0.3, framealpha=1, facecolor='white')

def ladder_by_design():
    E1 = json.loads((NF.LKP / 'ext2001_labs.json').read_text(encoding='utf-8'))
    S2 = NF.S2
    s8 = R3['S8']['base']
    pa = S2['primary']['any_category_change']
    b12 = S2['bp_benchmark']['reading1_vs_reading2']
    ex = S2['external']['any_category_change']
    pc = NF.pc
    pc2 = lambda e: (100 * e['estimate'], 100 * e['lower95'], 100 * e['upper95'])
    A = [('SBP reading 1 vs 2, within one visit', pc(s8['a_v1_r1_vs_r2'])), ('SBP reading 2 vs 3, within one visit', pc(s8['b_v1_r2_vs_r3'])), ('Single SBP reading, other day', pc(s8['c_v1r2_vs_v2r2'])), ('Mean of SBP readings 2–3, other day', pc(s8['d_between_visit_bp_swap'])), ('Lipids only, other day (attribution)', (100 * R3['K2_base']['single_swap']['LIPID']['ref1'], None, None)), ('Urine only, other day (UACR model; attribution)', (100 * R3['K2_uacr']['single_swap']['URINE']['ref1'], None, None)), ('Laboratory inputs, other day (BP fixed)', pc(E1['study1_comparator']['base_lab_only_ref1'])), ('Laboratory inputs + UACR, other day (UACR model)', pc(E1['study1_comparator']['uacr_lab_urine_ref1'])), ('All routinely measured inputs (θ)', pc(R3['primary'])), ('All inputs including UACR (UACR model)', pc(R3['K1']['theta_U']))]
    B = [(f"SBP reading 1 → 2, same {b12['n']:,} participants", pc2(b12['bp'])), ('Random → first-morning urine, same participants', pc2(b12['urine_same_participants'])), (f"Random → first-morning urine, all {pa['n']:,} participants", pc2(pa))]
    C = [('Random → random urine, two examinations', pc2(ex)), ('Laboratory inputs, other examination (BP fixed)', pc(E1['base']['theta_lab'])), ('Laboratory inputs + UACR (UACR model)', pc(E1['uacr']['theta_lab']))]
    fig = plt.figure(figsize=mm(180, 150))
    axA = fig.add_axes([0.5, 0.49, 0.46, 0.44])
    axB = fig.add_axes([0.5, 0.27, 0.46, 0.13])
    axC = fig.add_axes([0.5, 0.07, 0.46, 0.13])
    specs = [(axA, A, 'o', BLUE, f"Study 1: NHANES III (n = {R3['primary']['n']}; unweighted)"), (axB, B, 's', ORANGE, 'Study 2: NHANES 2009–2010 (weighted; UACR model)'), (axC, C, '^', '#72808C', f"NHANES 2001–2002 repeat (post hoc; n = {E1['base']['theta_lab']['n']})")]
    for ax, rows, mk, col, title in specs:
        for i, (lbl_, (m, lo, hi)) in enumerate(rows):
            if lo is not None:
                ax.plot([lo, hi], [i, i], color=col, lw=1.2)
            ax.plot(m, i, mk, color=col, ms=5, mfc=col if lo is not None else 'white', mew=1.2)
            ax.text((hi if lo is not None else m) + 0.3, i, f'{m:.1f}', va='center', ha='left', fontsize=6.8, color=INK)
        ax.set_yticks(range(len(rows)))
        ax.set_yticklabels([r[0] for r in rows], fontsize=7)
        ax.set_ylim(len(rows) - 0.5, -0.6)
        ax.set_xlim(0, 16)
        ax.tick_params(axis='y', length=0)
        ax.set_title(title, fontsize=8, loc='left')
        if ax is not axC:
            ax.set_xticklabels([])
    axC.set_xlabel('Changed PREVENT-ASCVD category, % (95% CI)')
    axA.plot([], [], 'o', color=INK, mfc='white', ms=5, label='Open circle: attribution\nswap (no interval)')
    axA.legend(fontsize=6.6, loc='lower left', frameon=False)
    for ax, l in zip([axA, axB, axC], 'ABC'):
        lab(ax, l, -1.05)
    return (fig, [{'id': 'A', 'question': 'Study 1 category change by source of variation', 'units': 'percent', 'sample_unit': 'participant', 'transform': 'none', 'estimate': 'proportions (unweighted)', 'interval': 'CR1 logit-t; none for attribution swaps', 'n': int(R3['primary']['n'])}, {'id': 'B', 'question': 'Study 2 specimen vs SBP reading swap', 'units': 'percent', 'sample_unit': 'participant', 'transform': 'none', 'estimate': 'survey-weighted proportions', 'interval': 'Taylor-logit', 'n': int(pa['n'])}, {'id': 'C', 'question': 'NHANES 2001-2002 repeat-sample proportions', 'units': 'percent', 'sample_unit': 'participant', 'transform': 'none', 'estimate': 'proportions', 'interval': 'cluster-adjusted', 'n': int(E1['base']['theta_lab']['n'])}])

def fig1():
    fig = plt.figure(figsize=mm(180, 150))
    a = fig.add_axes([0.08, 0.57, 0.4, 0.37])
    b = fig.add_axes([0.62, 0.6, 0.3, 0.31])
    c = fig.add_axes([0.09, 0.08, 0.37, 0.37])
    d = fig.add_axes([0.58, 0.08, 0.38, 0.37])
    e = R3['primary']
    NF.scatter_panel(a, P.riskB_m0.to_numpy(), P.riskB_m15.to_numpy(), f"Base PREVENT-ASCVD (n = {e['n']})\nchanged category: {e['events']} ({100 * e['estimate']:.1f}%)")
    NF.trans_panel(b, R3['S5_base']['transition_counts'], 'Visit 1 × visit 2 category (counts)')
    pi_risk_panel(c)
    NF.curve(d, 'K3_base', 'mech', 'Conditional on visit-mean risk ($\\pi_{mech}$), distance scale')
    fpp = lambda x: 100 * expit(logit(0.05) + x) - 5
    finv = lambda pp: logit((np.asarray(pp) + 5) / 100) - logit(0.05)
    sec = d.secondary_xaxis('top', functions=(fpp, finv))
    sec.set_xticks([0, 1, 2, 3])
    sec.tick_params(labelsize=7)
    sec.set_xlabel('')
    d.set_title('Visit-mean risk ($\\pi_{mech}$); top axis: pp above 5%', fontsize=8, loc='left', pad=20)
    for ax, l in zip([a, b, c, d], 'ABCD'):
        lab(ax, l, -0.16 if l in 'AC' else -0.22, 1.04 if l != 'D' else 1.15)
    n = int(e['n'])
    return (fig, [{'id': 'A', 'question': 'Base risk at both visits', 'units': 'percent risk', 'sample_unit': 'participant', 'transform': 'log10 axes', 'estimate': 'raw paired observations', 'interval': 'none', 'n': n}, {'id': 'B', 'question': 'Category transitions', 'units': 'counts', 'sample_unit': 'participant', 'transform': 'none', 'estimate': 'counts', 'interval': 'none', 'n': n}, {'id': 'C', 'question': 'Probability of a different category by first-visit risk, with >=25% zones', 'units': 'percent', 'sample_unit': 'participant', 'transform': 'log x; logit distance model', 'estimate': 'prespecified pi_clin (RCS logistic) on a risk grid (post hoc display)', 'interval': 'pointwise 95% cluster bootstrap (4,000)', 'n': n}, {'id': 'D', 'question': 'Probability by distance from the nearest threshold, conditional on visit-mean risk', 'units': 'percent', 'sample_unit': 'participant', 'transform': 'logit distance', 'estimate': 'pi_mech RCS logistic with threshold indicator', 'interval': 'pointwise 95% cluster bootstrap (4,000)', 'n': n}])

def fig2():
    fig = plt.figure(figsize=mm(180, 128))
    a = fig.add_axes([0.17, 0.58, 0.3, 0.34])
    b = fig.add_axes([0.67, 0.58, 0.3, 0.34])
    c = fig.add_axes([0.17, 0.09, 0.3, 0.34])
    d = fig.add_axes([0.67, 0.09, 0.3, 0.34])
    shapley_panel(a, 'K2_base', ['BP', 'LIPID', 'KIDNEY', 'GLYC'], 'Base model: Shapley split of θ = 10.4%')
    shapley_panel(b, 'K2_uacr', ['BP', 'LIPID', 'KIDNEY', 'GLYC', 'URINE'], 'UACR model: Shapley split of $\\theta_U$ = 11.6%')
    a.text(0.98, 0.04, 'BP ranked first in 92% of\nbootstrap replicates;\nBP − lipids 1.4 pp (−0.6 to 3.4)', transform=a.transAxes, ha='right', va='bottom', fontsize=6.8, color=INK)
    b.text(0.98, 0.5, 'No group ranked first in ≥90%\n(urine 61%, BP 38%)', transform=b.transAxes, ha='right', va='center', fontsize=6.8, color=INK)
    lp_partition_panel(c)
    benchmarks_panel(d)
    for ax, l in zip([a, b, c, d], 'ABCD'):
        lab(ax, l, -0.52)
    return (fig, [{'id': 'A', 'question': 'Which inputs drive base-model category change?', 'units': 'percentage points', 'sample_unit': 'participant', 'transform': 'none', 'estimate': 'symmetrized Shapley', 'interval': '95% cluster bootstrap (4,000)', 'n': 759}, {'id': 'B', 'question': 'Same for the UACR model', 'units': 'percentage points', 'sample_unit': 'participant', 'transform': 'none', 'estimate': 'symmetrized Shapley', 'interval': '95% cluster bootstrap (4,000)', 'n': 725}, {'id': 'C', 'question': 'Continuous-scale partition', 'units': 'percent share', 'sample_unit': 'participant', 'transform': 'logit scale', 'estimate': 'E[f_k dLP]/E[dLP^2]', 'interval': '95% cluster bootstrap', 'n': 759}, {'id': 'D', 'question': 'Model-based simulation benchmarks', 'units': 'percent', 'sample_unit': 'participant', 'transform': 'none', 'estimate': 'simulated or observed proportions', 'interval': 'none', 'n': 759}])

def fig4():
    fig = plt.figure(figsize=mm(180, 150))
    a = fig.add_axes([0.21, 0.58, 0.33, 0.35])
    b = fig.add_axes([0.65, 0.58, 0.32, 0.35])
    c = fig.add_axes([0.1, 0.08, 0.32, 0.35])
    d = fig.add_axes([0.56, 0.08, 0.41, 0.35])
    e11_loss_panel(a)
    e11_regret_panel(b)
    e11_sens_panel(c)
    repeat_panel(d)
    for ax, l, x in [(a, 'A', -0.62), (b, 'B', -0.08), (c, 'C', -0.26), (d, 'D', -0.16)]:
        lab(ax, l, x)
    return (fig, [{'id': 'A', 'question': 'Expected net benefit forfeited by single-visit, two-day and targeted decisions', 'units': 'percent of achievable', 'sample_unit': 'participant (model-based)', 'transform': 'none', 'estimate': 'decision-curve regret / achievable net benefit', 'interval': '95% percentile, 1,000 uniform-start bootstrap refits', 'n': 749}, {'id': 'B', 'question': 'Regret per 1,000 against achievable net benefit', 'units': 'net true positives per 1,000', 'sample_unit': 'participant', 'transform': 'none', 'estimate': 'model-based', 'interval': 'as A', 'n': 749}, {'id': 'C', 'question': 'Sensitivity of the single-visit loss to the error model', 'units': 'percent', 'sample_unit': 'participant', 'transform': 'none', 'estimate': 'age-band SD; variance x1.5, x2', 'interval': 'none', 'n': 749}, {'id': 'D', 'question': 'Reduction in simulated disagreement versus extra measurement days', 'units': 'percent', 'sample_unit': 'participant', 'transform': '1 - ratio to single', 'estimate': 'simulation (1,000 replicates per participant)', 'interval': 'cluster bootstrap (4,000), simulation SDs fixed', 'n': 759}])

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    src = NF.SRC_CI + [NF.LKP / 'ext5_stakes_results.json', NF.LKP / 'repeat_strategies.json', Path(__file__)]
    titles = {1: 'Reclassification and where it happens', 2: 'Sources of category change', 3: 'Category change by source of measurement variation, by design', 4: 'Decisions, net benefit at stake and repeat strategies'}
    only = sys.argv[1:]
    for i, f in [(1, fig1), (2, fig2), (3, ladder_by_design), (4, fig4)]:
        if only and str(i) not in only:
            continue
        with publication_style(font_size=7.5):
            fig, panels = f()
            export_figure(fig, OUT / f'fig{i}', sources=src, panels=panels, title=titles[i], dpi=300, overwrite=True)
        plt.close(fig)
        print('fig', i, 'ok')
if __name__ == '__main__':
    main()
