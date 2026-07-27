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

    # aspect is equal because PC1 and PC2 carry comparable variance (26% / 21%),
    # so a unit of PC1 must be a unit of PC2 on the page; limits are padded by
    # 4% of the data range, which also removes the empty margin of the old
    # square panel.
    pad = 0.04
    xr = float(np.ptp(Z[:, 0]))
    yr = float(np.ptp(Z[:, 1]))
    xlim = (Z[:, 0].min() - pad * xr, Z[:, 0].max() + pad * xr)
    ylim = (Z[:, 1].min() - pad * yr, Z[:, 1].max() + pad * yr)
    w_mm = 57.4
    h_mm = w_mm * (ylim[1] - ylim[0]) / (xlim[1] - xlim[0])

    fig, ax = S.panel(w_mm, h_mm)
    # draw the most numerous gene first so the smaller sets are not buried; this
    # changes only paint order, not a single plotted coordinate
    order = sorted(S.GENE_ORDER, key=lambda g: -(df.gene == g).sum())
    for gene in order:
        m = (df.gene == gene).to_numpy()
        ax.scatter(Z[m, 0], Z[m, 1], s=3.5, color=S.GENE_COLORS[gene], alpha=0.7,
                   linewidths=0, label=gene)  # vector: 470 points, no rasterization

    ax.set_xlabel(f"PC1 ({ev[0]:.0f}%)", labelpad=1)
    ax.set_ylabel(f"PC2 ({ev[1]:.0f}%)", labelpad=1)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect("equal", adjustable="box")
    S.despine(ax, keep=("left", "bottom"))
    ax.axhline(0, color=S.LIGHT_GREY, lw=0.4, zorder=0)
    ax.axvline(0, color=S.LIGHT_GREY, lw=0.4, zorder=0)
    # No gene legend here. Gene colours are direct-labelled by the track titles
    # of panel b and the tick labels of panel c; Figure 1 carries ONE gene key,
    # not one per panel.
    ax.text(0.0, 1.015, f"{len(df)} variants, coloured by gene",
            transform=ax.transAxes, ha="left", va="bottom", fontsize=5.5,
            color=S.GREY)

    S.save(fig, "fig1d_theta_pca")
    print("explained variance %:", np.round(ev, 1), " n variants:", len(df))


if __name__ == "__main__":
    main()
