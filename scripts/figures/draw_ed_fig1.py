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

Style: this panel is drawn through manuscript/figures/nm_style.py, the same shared
house style as the six main figures. It previously set its own rcParams with a
font stack of Arial/Helvetica/DejaVu Sans and, since neither Arial nor Helvetica is
installed, fell through to DejaVu Sans, so Extended Data Fig. 1 was the only
display item in the paper set in a different typeface. Its three external atlases
also carried invented hues that collided with the gene palette.

Run from: AllelePerturb/scripts/figures/
"""
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, os.path.join(ROOT, 'manuscript', 'figures'))
import nm_style as S  # noqa: E402

CSV = os.path.join(ROOT, 'results', 'rankability_sensitivity.csv')

# criterion label -> display label (order preserved)
CRIT_ORDER = ['S > W (current)', 'D_self_hi < D_null', 'S > 1.5×W', 'S > 2×W', 'ratio < 0.8']
CRIT_DISPLAY = ['S > W\n(default)', 'CI(S) > 0', 'S > 1.5W', 'S > 2W',
                r'$D_\mathrm{self}/D_\mathrm{null}$' '\n< 0.8']

# dataset -> colour: the shared house mapping (gene hues for the four allele genes,
# neutral slate ramp for the external atlases), identical to Fig. 6.
DS_ORDER = ['TP53', 'KRAS', 'GATA1', 'JAK1', 'Replogle', 'Norman', 'Adamson']
DS_COLORS = {d: S.DATASET_COLORS[d] for d in DS_ORDER}


def main():
    df = pd.read_csv(CSV)
    # pivot to dataset x criterion (%)
    pivot = (df.pivot(index='dataset', columns='criterion', values='unrankable_frac') * 100.0)
    pivot = pivot.reindex(index=DS_ORDER, columns=CRIT_ORDER)

    S.apply_rcparams()
    plt.rcParams.update({'figure.dpi': 300, 'savefig.dpi': 300})

    n_crit = len(CRIT_ORDER)
    n_ds = len(DS_ORDER)
    fig, ax = plt.subplots(figsize=(178.0 * S.MM, 81.0 * S.MM))
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
