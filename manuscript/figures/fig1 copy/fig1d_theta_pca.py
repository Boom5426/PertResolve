"""Figure 1d (inset) - Variants placed in a shared protein-feature space.

Data-direct companion to the (AI-drawn) Figure 1d schematic. It shows the 470
real variants embedded by PCA of the 6-dimensional biophysical theta vector,
coloured by gene, to make concrete that variants from different genes share one
feature space (which is what enables held-out variant prediction). No result is
implied; this is a representation cartoon.

Source (committed): data/allele_perturb_bench.csv, theta columns
  d_hydro, d_vol, d_charge, fold_core, cat_switch, is_hotspot
(standardised again before PCA so all 6 dims contribute comparably).

Run:  python fig1d_theta_pca.py
Out:  fig1d_theta_pca.pdf (+ .png preview)
"""
from __future__ import annotations

import os
import sys

import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import nm_style as S

THETA_COLS = ["d_hydro", "d_vol", "d_charge", "fold_core", "cat_switch", "is_hotspot"]


def main() -> None:
    S.apply_rcparams()
    df = S.load_bench(exclude_wt=True)

    X = StandardScaler().fit_transform(df[THETA_COLS].to_numpy(float))
    pca = PCA(n_components=2, random_state=0).fit(X)
    Z = pca.transform(X)
    ev = pca.explained_variance_ratio_ * 100

    fig, ax = S.panel(46, 46)
    for gene in S.GENE_ORDER:
        m = (df.gene == gene).to_numpy()
        ax.scatter(Z[m, 0], Z[m, 1], s=5, color=S.GENE_COLORS[gene], alpha=0.75,
                   linewidths=0, label=gene)  # vector: ~470 points, no rasterization

    ax.set_xlabel(f"PC1 ({ev[0]:.0f}%)", labelpad=1)
    ax.set_ylabel(f"PC2 ({ev[1]:.0f}%)", labelpad=1)
    ax.set_xticks([])
    ax.set_yticks([])
    S.despine(ax, keep=("left", "bottom"))
    ax.axhline(0, color=S.LIGHT_GREY, lw=0.4, zorder=0)
    ax.axvline(0, color=S.LIGHT_GREY, lw=0.4, zorder=0)
    leg = ax.legend(loc="upper left", handletextpad=0.2, borderpad=0.2,
                    labelspacing=0.25, markerscale=1.4, bbox_to_anchor=(-0.02, 1.03))
    for t, g in zip(leg.get_texts(), S.GENE_ORDER):
        t.set_color(S.GENE_COLORS[g])

    S.save(fig, "fig1d_theta_pca")
    print("explained variance %:", np.round(ev, 1), " n variants:", len(df))


if __name__ == "__main__":
    main()
