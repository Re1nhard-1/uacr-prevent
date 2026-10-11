from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import nh3_figures as NF
import nh3_figures_v7 as F7
import build_nh3_manuscript as M
from figure_tools import publication_style, export_figure
OUT = NF.ROOT / 'studies/nhanes3_testretest/submission_ajpc/figures'

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    pr = M.pr
    with publication_style(font_size=7.5):
        fig = plt.figure(figsize=F7.mm(180, 110))
        a = fig.add_axes([0.08, 0.13, 0.37, 0.62])
        b = fig.add_axes([0.66, 0.13, 0.32, 0.62])
        F7.pi_risk_panel(a)
        F7.e11_loss_panel(b)
        a.set_title('Where the category changes', fontsize=8, loc='left')
        b.set_title('What is at stake (model-based)', fontsize=8, loc='left')
        fig.text(0.5, 0.955, f"Remeasuring PREVENT inputs 2 weeks later changed the risk category of {100 * pr['estimate']:.1f}% of adults,", ha='center', va='top', fontsize=9, fontweight='bold', color=NF.INK)
        fig.text(0.5, 0.905, 'almost only near the 3%, 5% and 10% thresholds, where little net benefit is at stake', ha='center', va='top', fontsize=9, fontweight='bold', color=NF.INK)
        fig.text(0.5, 0.855, 'NHANES III reexamination, n = 759, median 16 days, same protocol; 2025 BP and 2026 lipid guideline thresholds', ha='center', va='top', fontsize=7, color='#5A6670')
        for ax, l, x in [(a, 'A', -0.17), (b, 'B', -0.7)]:
            F7.lab(ax, l, x)
        panels = [{'id': 'A', 'question': 'Probability of a different category by first-visit risk, with >=25% zones', 'units': 'percent', 'sample_unit': 'participant', 'transform': 'log x; logit distance model', 'estimate': 'prespecified pi_clin', 'interval': 'pointwise 95% cluster bootstrap (4,000)', 'n': int(pr['n'])}, {'id': 'B', 'question': 'Expected net benefit forfeited by single-visit, two-day and targeted decisions', 'units': 'percent of achievable', 'sample_unit': 'participant (model-based)', 'transform': 'none', 'estimate': 'decision-curve regret / achievable net benefit', 'interval': '95% percentile, 1,000 uniform-start bootstrap refits', 'n': 749}]
        src = NF.SRC_CI + [NF.LKP / 'ext5_stakes_results.json', Path(__file__)]
        export_figure(fig, OUT / 'Central_Illustration', sources=src, panels=panels, title='Central Illustration', dpi=400, formats=('png', 'pdf', 'svg'), overwrite=True)
    plt.close(fig)
    print('saved', OUT / 'Central_Illustration.png')
if __name__ == '__main__':
    main()
