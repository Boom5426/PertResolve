"""
Supplementary Figure 1 (rankability verdict sensitivity to criterion choice).

Grouped bar chart: un-rankable fraction for 7 datasets under 5 rankability criteria.
Data: results/rankability_sensitivity.csv (5 criteria x 7 datasets = 35 rows).

Note on the ratio criterion: a variant is un-rankable when D_self/D_null >= 0.8
(within-replicate distance is a large fraction of the variant-to-WT distance, i.e.
signal does not clearly exceed replicate noise). An earlier version of the
`ratio < 0.8` column stored the complementary (rankable) fraction; it has been
sign-corrected so all five criteria report the un-rankable fraction on the same
convention. The remote generator should carry the same sign fix.

Run from: AllelePerturb/scripts/figures/
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
CSV = os.path.join(ROOT, 'results', 'rankability_sensitivity.csv')

# criterion label -> display label (order preserved)
CRIT_ORDER = ['S > W (current)', 'D_self_hi < D_null', 'S > 1.5×W', 'S > 2×W', 'ratio < 0.8']
CRIT_DISPLAY = ['S > W\n(default)', 'CI(S) > 0', 'S > 1.5W', 'S > 2W', 'D/D_null\n< 0.8']

# dataset -> colour (allele genes share the manuscript GENE_COLORS; atlases distinct)
DS_ORDER = ['TP53', 'KRAS', 'GATA1', 'JAK1', 'Replogle', 'Norman', 'Adamson']
DS_COLORS = {
    'TP53': '#5185C0', 'KRAS': '#E99D4E', 'GATA1': '#8281B9', 'JAK1': '#55966B',
    'Replogle': '#C96144', 'Norman': '#D6B26B', 'Adamson': '#6BAAA7',
}


def main():
    df = pd.read_csv(CSV)
    # pivot to dataset x criterion (%)
    pivot = (df.pivot(index='dataset', columns='criterion', values='unrankable_frac') * 100.0)
    pivot = pivot.reindex(index=DS_ORDER, columns=CRIT_ORDER)

    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'font.size': 7, 'axes.labelsize': 8, 'xtick.labelsize': 7, 'ytick.labelsize': 7,
        'legend.fontsize': 6.5, 'figure.dpi': 300, 'savefig.dpi': 300,
        'axes.spines.top': False, 'axes.spines.right': False,
        'axes.linewidth': 0.6, 'xtick.major.width': 0.5, 'ytick.major.width': 0.5,
    })

    n_crit = len(CRIT_ORDER)
    n_ds = len(DS_ORDER)
    fig, ax = plt.subplots(figsize=(7.0, 3.2))
    group_w = 0.8
    bar_w = group_w / n_ds
    x = np.arange(n_crit)
    for i, ds in enumerate(DS_ORDER):
        vals = pivot.loc[ds].values
        offs = (i - (n_ds - 1) / 2.0) * bar_w
        ax.bar(x + offs, vals, bar_w, label=ds, color=DS_COLORS[ds],
               edgecolor='white', linewidth=0.2)

    ax.axhline(50, ls=':', lw=0.6, color='0.5', zorder=0)
    ax.set_ylabel('Un-rankable (%)')
    ax.set_ylim(0, 105)
    ax.set_yticks([0, 20, 40, 60, 80, 100])
    ax.set_xticks(x)
    ax.set_xticklabels(CRIT_DISPLAY)
    ax.set_xlim(-0.5, n_crit - 0.5)
    ax.legend(ncol=4, loc='upper center', bbox_to_anchor=(0.5, 1.20),
              frameon=False, columnspacing=1.0, handlelength=1.1)
    fig.tight_layout()

    outs = [
        os.path.join(ROOT, 'manuscript', 'latex', 'figures', 'ED_fig1.pdf'),
        os.path.join(ROOT, 'manuscript', 'figures', 'ED_fig1.pdf'),
        os.path.join(ROOT, 'manuscript', 'figures', 'extended_data', 'ED_fig1_rankability_sensitivity.pdf'),
    ]
    for o in outs:
        os.makedirs(os.path.dirname(o), exist_ok=True)
        fig.savefig(o, bbox_inches='tight')
        print('wrote', os.path.relpath(o, ROOT))
    plt.close(fig)


if __name__ == '__main__':
    main()
