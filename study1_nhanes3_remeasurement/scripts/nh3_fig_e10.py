from pathlib import Path
import json
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, str(Path(__file__).resolve().parent / 'vendor'))
from figure_tools import publication_style, export_figure
ROOT = Path(__file__).resolve().parents[1]
SRCJ = ROOT / 'studies/nhanes3_testretest/analysis_v1_posthoc/ext4_30yr_results.json'
OUT = ROOT / 'studies/nhanes3_testretest/figures_v1'
R = json.loads(SRCJ.read_text(encoding='utf-8'))
A, B = (R['part_A'], R['part_B'])
INK, TEAL, ORANGE, GREY, GRID = ('#203443', '#187D84', '#CB6C38', '#B9C3CB', '#E3E8EC')
Z25, Z10 = ('#F2C9AE', '#FBE9DD')
mm = lambda w, h: (w / 25.4, h / 25.4)

def lab(ax, s, x=-0.16, y=1.04):
    ax.text(x, y, s, transform=ax.transAxes, fontweight='bold', fontsize=10, color=INK, va='bottom')

def curve_panel(ax, g, xkey, bins, lokey, hikey, zones, thr, xlim, xlabel, title, binmid):
    x = np.array(g[xkey])
    m = (x >= xlim[0]) & (x <= xlim[1])
    for lev, col in (('0.1', Z10), ('0.25', Z25)):
        if lev in zones:
            ax.axvspan(*zones[lev], color=col, lw=0, zorder=0)
    ax.fill_between(x[m], 100 * np.array(g['lower'])[m], 100 * np.array(g['upper'])[m], color=TEAL, alpha=0.18, lw=0)
    ax.plot(x[m], 100 * np.array(g['p'])[m], color=TEAL, lw=1.4)
    ax.axvline(thr, color=INK, lw=0.8, ls='--')
    for b in bins:
        if b['n'] == 0:
            continue
        xm = binmid(b[lokey], b[hikey])
        if not xlim[0] <= xm <= xlim[1]:
            continue
        p = 100 * b['changes'] / b['n']
        ax.plot(xm, p, 'D', color='white', mec=INK, ms=4, zorder=5)
        ax.text(xm, p + 6.5, f"{b['changes']}/{b['n']}", fontsize=6.5, ha='center', va='bottom', color='#5A6670', zorder=6, bbox=dict(boxstyle='round,pad=0.12', fc='white', ec='none', alpha=0.85))
    ax.set_xlim(*xlim)
    ax.set_ylim(0, 100)
    ax.set_xlabel(xlabel)
    ax.set_ylabel('Chance of a different classification\non remeasurement (%)')
    ax.set_title(title, fontsize=8, loc='left')

def fig():
    f = plt.figure(figsize=mm(180, 138))
    ax = [f.add_axes([0.09, 0.58, 0.38, 0.34]), f.add_axes([0.59, 0.58, 0.38, 0.34]), f.add_axes([0.27, 0.13, 0.21, 0.29]), f.add_axes([0.59, 0.09, 0.38, 0.31])]
    p30 = A['pi30']
    zl = p30['zones_risk30_pct']
    curve_panel(ax[0], p30['grid'], 'risk30_pct', p30['observed_bins'], 'risk30_lo', 'risk30_hi', zl, 10, (2, 25), 'First-visit 30-year ASCVD risk (%)', f"A  30-year risk ≥10% (Study 1, age 30–59, n={A['n_A1']})", lambda lo, hi: (lo + min(hi, 30)) / 2 if hi < 100 else 27)
    ax[0].text(24.6, 92, f"Shaded: chance ≥25% ({zl['0.25'][0]:.1f}–{zl['0.25'][1]:.1f}%)\nand ≥10% ({zl['0.1'][0]:.1f}–{zl['0.1'][1]:.1f}%)", ha='right', va='top', fontsize=6.5, color=INK, bbox=dict(boxstyle='round,pad=0.3', fc='white', ec=GRID, lw=0.6))
    L = A['ldl160']['curve']
    zL = L['zones_ldl_mgdl']
    curve_panel(ax[1], L['grid'], 'ldl_mgdl', L['observed_bins'], 'ldl_lo', 'ldl_hi', zL, 160, (110, 220), 'First-visit LDL cholesterol (mg/dL)', f"B  LDL cholesterol ≥160 mg/dL (Study 1, n={A['ldl160']['n']})", lambda lo, hi: 120 if lo == 0 else 205 if hi >= 400 else (lo + hi) / 2)
    ax[1].text(217, 92, f"Shaded: chance ≥25% ({zL['0.25'][0]:.0f}–{zL['0.25'][1]:.0f})\nand ≥10% ({zL['0.1'][0]:.0f}–{zL['0.1'][1]:.0f} mg/dL)", ha='right', va='top', fontsize=6.5, color=INK, bbox=dict(boxstyle='round,pad=0.3', fc='white', ec=GRID, lw=0.6))
    E = A['pathway_domain_v1']['among_eligible_v1']
    rs = E['not_pathway_eligible_v2_reasons']
    left = rs['risk30_ge10_but_left_domain']
    n = E['n']
    still = n - E['not_pathway_eligible_v2']['events']
    rows = [('Still eligible via 30-year risk', still, TEAL), ('30-year risk fell below 10%', rs['risk30_below_10'], ORANGE), ('LDL-C rose to ≥160 mg/dL', left['ldl_ge160'], '#B0548C'), ('10-year risk rose to ≥3%', left['tenyear_ge3'], '#3766A0'), ('Met a risk-independent indication', left['indication'], '#7A6A58')]
    y = np.arange(len(rows))[::-1]
    ax[2].barh(y, [r[1] for r in rows], color=[r[2] for r in rows], height=0.62)
    for yi, r in zip(y, rows):
        ax[2].text(r[1] + 0.8, yi, f'{r[1]}', va='center', fontsize=7, color=INK)
    ax[2].set_yticks(y)
    ax[2].set_yticklabels([r[0] for r in rows], fontsize=7)
    ax[2].set_xlim(0, max((r[1] for r in rows)) * 1.25)
    ax[2].set_xlabel('Participants at visit 2 (n)')
    f.text(0.02, 0.455, f'C  Eligible via 30-year risk at visit 1 (Study 1, n={n}): status at visit 2', fontsize=8, ha='left', va='bottom', color=INK)
    f.text(0.02, 0.015, f"Not eligible via this route at visit 2: {E['not_pathway_eligible_v2']['events']}/{n} ({100 * E['not_pathway_eligible_v2']['estimate']:.0f}%);\nlost any statin recommendation: {E['cor1_or_2a_differs_v2']['events']}/{n}. Conditions on visit 1.", fontsize=6.5, ha='left', va='bottom', color='#5A6670')
    H = B['B1_risk30_histogram']
    e = np.array(H['bin_edges'])
    w = np.array(H['weighted_millions_all_B1'])
    mid = (e[:-1] + e[1:]) / 2
    col = [TEAL if a >= 10 else GREY for a in e[:-1]]
    ax[3].axvspan(*zl['0.25'], color=Z25, lw=0, zorder=0)
    ax[3].bar(mid, w, width=0.5, color=col, lw=0)
    ax[3].axvline(10, color=INK, lw=0.8, ls='--')
    ax[3].set_xlim(0, 30)
    ax[3].set_xlabel('30-year ASCVD risk (%)')
    ax[3].set_ylabel('US adults (millions per 0.5% band)')
    el = B['B1_pathway_eligible']
    z = B['eligible_in_pi30_zone_0.25']
    ax[3].set_title('D  Pathway domain, US adults 30–59\n    (NHANES 2021–2023, survey-weighted)', fontsize=8, loc='left')
    ax[3].text(29.6, ax[3].get_ylim()[1] * 0.97, f"Eligible via 30-year risk\n≥10% (teal): {el['millions']['millions']:.1f} million\n({el['millions']['lower']:.1f}–{el['millions']['upper']:.1f})\nIn the ≥25% zone (shaded):\n{100 * z['estimate']:.0f}% of them ({100 * z['with_curve_uncertainty']['lower']:.0f}–{100 * z['with_curve_uncertainty']['upper']:.0f}%),\n{z['millions']['millions']:.1f} million\nExploratory: assumes\n1988–1994 variation", ha='right', va='top', fontsize=6.5, color=INK, linespacing=1.2, bbox=dict(boxstyle='round,pad=0.3', fc='white', ec=GRID, lw=0.6))
    panels = [{'id': 'A', 'question': 'Chance that remeasurement changes the 30-year >=10% classification, by first-visit 30-year risk', 'units': 'percent', 'sample_unit': 'participant', 'transform': 'logit distance from 10%', 'estimate': 'logistic regression, restricted cubic spline (3 knots)', 'interval': 'pointwise 95% cluster-bootstrap percentile (4,000 replicates)', 'n': A['n_A1']}, {'id': 'B', 'question': 'Chance that remeasurement changes the LDL-C >=160 classification', 'units': 'percent', 'sample_unit': 'participant', 'transform': 'log distance from 160 mg/dL', 'estimate': 'logistic regression, restricted cubic spline (3 knots)', 'interval': 'pointwise 95% cluster-bootstrap percentile (4,000 replicates)', 'n': A['ldl160']['n']}, {'id': 'C', 'question': 'Visit-2 status of visit-1 pathway-eligible participants', 'units': 'count', 'sample_unit': 'participant', 'transform': 'none', 'estimate': 'counts (conditions on visit 1; includes regression to the mean)', 'interval': 'none shown', 'n': n}, {'id': 'D', 'question': 'Survey-weighted distribution of 30-year risk in the pathway domain, NHANES 2021-2023', 'units': 'millions', 'sample_unit': 'participant (survey-weighted, fasting subsample)', 'transform': 'none', 'estimate': 'weighted totals, complete case', 'interval': 'text: Taylor survey interval; zone share with curve uncertainty', 'n': B['B1_n']}]
    return (f, panels)

def main():
    with publication_style(font_size=7.5):
        f, panels = fig()
        export_figure(f, OUT / 'figure_e10', sources=[SRCJ, Path(__file__)], panels=panels, title='30-year-risk statin pathway: remeasurement stability and US reach', dpi=300, overwrite=True)
    plt.close(f)
    print('ok')
if __name__ == '__main__':
    main()
