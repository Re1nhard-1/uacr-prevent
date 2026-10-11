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
SRCJ = ROOT / 'studies/nhanes3_testretest/analysis_v1_posthoc/ext5_stakes_results.json'
OUT = ROOT / 'studies/nhanes3_testretest/figures_v1'
R = json.loads(SRCJ.read_text(encoding='utf-8'))
INK, TEAL, ORANGE, GREY, GRID, BLUE = ('#203443', '#187D84', '#CB6C38', '#B9C3CB', '#E3E8EC', '#3766A0')
mm = lambda w, h: (w / 25.4, h / 25.4)
ROWS = [('ascvd10', '3', '10-y ASCVD 3%'), ('ascvd10', '5', '10-y ASCVD 5%\n(statin, Class 1)'), ('cvd10_bp', '7.5', '10-y total CVD 7.5%\n(BP drug threshold)'), ('ascvd30', '10', '30-y ASCVD 10%\n(statin, Class 2a; 30–59 y)'), ('ascvd10', '10', '10-y ASCVD 10%\n(high intensity)')]
E = [R[k]['thresholds'][t] for k, t, _ in ROWS]
LAB = [l for _, _, l in ROWS]
y = np.arange(len(ROWS))[::-1]

def fig():
    f = plt.figure(figsize=mm(180, 82))
    a = f.add_axes([0.21, 0.17, 0.3, 0.72])
    b = f.add_axes([0.6, 0.17, 0.17, 0.72])
    c = f.add_axes([0.82, 0.17, 0.16, 0.72])
    for j, (key, col, mk, lab) in enumerate([('loss_fraction_single', INK, 'o', 'Single visit'), ('loss_fraction_two_day', BLUE, 's', 'Two-day average'), ('loss_fraction_targeted', TEAL, 'D', 'Second day only in the ≥25% zone')]):
        dy = (1 - j) * 0.22
        for yi, e in zip(y, E):
            v = 100 * e[key]
            ci = e['ci_uniform_start'][key]
            a.plot([100 * ci[0], 100 * ci[1]], [yi + dy] * 2, color=col, lw=1, alpha=0.8)
            a.plot(v, yi + dy, mk, color=col, ms=4.5, zorder=5)
        a.plot([], [], mk, color=col, ms=4.5, label=lab)
    a.axvspan(0, 5, color='#EEF3F5', lw=0, zorder=0)
    a.axvspan(5, 15, color='#FBE9DD', lw=0, zorder=0)
    a.text(2.5, -0.5, 'minor', ha='center', fontsize=6.5, color='#5A6670')
    a.text(10, -0.5, 'modest', ha='center', fontsize=6.5, color='#5A6670')
    a.set_yticks(y)
    a.set_yticklabels(LAB, fontsize=7)
    a.set_xlim(0, 15)
    a.set_ylim(-0.6, len(y) - 0.3)
    a.set_xlabel('Expected net benefit forfeited (% of achievable)')
    a.legend(fontsize=6.5, loc='upper right', frameon=True)
    a.set_title('A  Net benefit forfeited, by measurement strategy', fontsize=8, loc='left')
    for yi, e in zip(y, E):
        ci = e['ci_uniform_start']['regret_single_per_1000']
        b.barh(yi, e['regret_single_per_1000'], color=INK, height=0.5)
        b.plot(ci, [yi] * 2, color=ORANGE, lw=1.2)
        b.text(max(ci[1], e['regret_single_per_1000']) + 0.03, yi, f"of {e['achievable_nb_per_1000']:.1f}", va='center', fontsize=6.5, color='#5A6670')
    b.set_yticks(y)
    b.set_yticklabels([])
    b.set_xlim(0, 1.25)
    b.set_xlabel('Net true positives lost\nper 1,000 (single visit)')
    b.set_title('B  Absolute scale', fontsize=8, loc='left')
    for yi, e in zip(y, E):
        sb = e.get('post_audit_sens_sigma_by_age', {}).get('loss_fraction_single')
        vi = e['post_audit_sens_variance_inflation']
        c.plot(100 * e['loss_fraction_single'], yi, 'o', color=INK, ms=4, zorder=5)
        if sb is not None:
            c.plot(100 * sb, yi, '^', color=BLUE, ms=4, zorder=5)
        c.plot(100 * vi['1.5']['loss_fraction_single'], yi, '>', color=ORANGE, ms=4, zorder=5)
        c.plot(100 * vi['2.0']['loss_fraction_single'], yi, '>', mfc='white', color=ORANGE, ms=4, zorder=5)
    for mk, col, mfc, lab in [('o', INK, INK, 'Main'), ('^', BLUE, BLUE, 'Error SD by age'), ('>', ORANGE, ORANGE, 'Variance ×1.5'), ('>', ORANGE, 'white', 'Variance ×2')]:
        c.plot([], [], mk, color=col, mfc=mfc, ms=4, label=lab)
    c.set_yticks(y)
    c.set_yticklabels([])
    c.set_xlim(0, 20)
    c.set_xlabel('Forfeited (%)')
    c.legend(fontsize=6.5, loc='upper right', frameon=True, handletextpad=0.3, borderpad=0.3)
    c.set_title('C  Sensitivity', fontsize=8, loc='left')
    panels = [{'id': 'A', 'question': 'Expected net benefit forfeited by single-visit, two-day and targeted decisions relative to usual-level decisions', 'units': 'percent of achievable net benefit', 'sample_unit': 'participant (model-based over the fitted usual-level distribution)', 'transform': 'none', 'estimate': 'decision-curve regret / achievable net benefit', 'interval': '95% percentile from 1,000 cluster bootstrap replicates refitted from a uniform start', 'n': 749}, {'id': 'B', 'question': 'Regret per 1,000 people (single visit) against achievable net benefit (populations: 749 statin-naive; 623 not on BP medication; 499 aged 30-59 statin-naive)', 'units': 'net true positives per 1,000', 'sample_unit': 'participant', 'transform': 'none', 'estimate': 'model-based', 'interval': 'as A', 'n': 749}, {'id': 'C', 'question': 'Sensitivity of the single-visit loss fraction to the error model', 'units': 'percent', 'sample_unit': 'participant', 'transform': 'none', 'estimate': 'age-band error SD; error variance inflated 1.5x and 2x', 'interval': 'none', 'n': 749}]
    return (f, panels)

def main():
    with publication_style(font_size=7.5):
        f, panels = fig()
        export_figure(f, OUT / 'figure_e11', sources=[SRCJ, Path(__file__)], panels=panels, title='Net benefit at stake in measurement-driven decision changes (model-based)', dpi=300, overwrite=True)
    plt.close(f)
    print('ok')
if __name__ == '__main__':
    main()
