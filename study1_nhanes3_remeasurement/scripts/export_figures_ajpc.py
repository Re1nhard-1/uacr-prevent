from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import nh3_figures as NF
import nh3_figures_v7 as F7
from figure_tools import publication_style, export_figure
OUT = NF.ROOT / 'studies/nhanes3_testretest/submission_ajpc/figures'

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    src = NF.SRC_CI + [NF.LKP / 'ext5_stakes_results.json', NF.LKP / 'repeat_strategies.json', Path(__file__), Path(F7.__file__)]
    titles = {1: 'Reclassification and where it happens', 2: 'Sources of category change', 3: 'Category change by source of measurement variation, by design', 4: 'Decisions, net benefit at stake and repeat strategies'}
    for i, f in [(1, F7.fig1), (2, F7.fig2), (3, F7.ladder_by_design), (4, F7.fig4)]:
        with publication_style(font_size=7.5):
            fig, panels = f()
            export_figure(fig, OUT / f'Figure_{i}', sources=src, panels=panels, title=titles[i], dpi=400, formats=('png', 'pdf', 'svg'), overwrite=True)
        plt.close(fig)
        print('Figure', i, 'ok')
if __name__ == '__main__':
    main()
