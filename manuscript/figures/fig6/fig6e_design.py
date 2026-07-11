"""Figure 6e - Effect-size-conditioned design landscape.

Data-direct. Source: results/second_probe_rankability_table.csv via
D.native_rankability(gene), one row per allele-gene perturbation at its deepest
split-half bin. Plots variant effect size (x, log) against cells per variant
(y, log), coloured by gene, with fill encoding rankability: filled = rankable,
open (white face) = un-rankable.

Message: rankability is set by *effect size*, not by a universal cell cutoff.
High-effect JAK1 variants are rankable even at low depth, while low-effect
TP53/KRAS variants stay un-rankable even at 300 cells; GATA1 is intermediate.
No cell-number cutoff line is drawn and no contour is fabricated.

Run:  python fig6e_design.py  ->  fig6e_design.pdf (+ .png)
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib.ticker as mticker
from matplotlib.lines import Line2D

import nm_style as S
import fig6_data as D

# faint gene-cluster label anchors in data coordinates (x=effect size, y=cells);
# placed in open regions near each cluster, verified against the plotted ranges.
GENE_LABEL_POS = {
    "TP53":  (0.78, 175),
    "KRAS":  (0.70, 320),
    "GATA1": (1.65, 78),
    "JAK1":  (13.6, 210),
}
GENE_LABEL_HA = {"TP53": "left", "KRAS": "left", "GATA1": "left", "JAK1": "right"}


def main() -> None:
    S.apply_rcparams()
    fig, ax = S.panel(58, 48)

    for gene in S.GENE_ORDER:
        n = D.native_rankability(gene)
        es = n["effect_size"].to_numpy(float)
        nc = n["n_cells"].to_numpy(float)
        rank = n["rankable"].to_numpy(bool)
        good = np.isfinite(es) & np.isfinite(nc) & (es > 0) & (nc > 0)
        es, nc, rank = es[good], nc[good], rank[good]
        col = S.GENE_COLORS[gene]

        # un-rankable: hollow (white face, gene-colour edge)
        m = ~rank
        if m.any():
            ax.scatter(es[m], nc[m], s=9, facecolors="white", edgecolors=col,
                       linewidths=0.5, alpha=0.85, zorder=3)
        # rankable: filled gene colour
        m = rank
        if m.any():
            ax.scatter(es[m], nc[m], s=9, facecolors=col, edgecolors=col,
                       linewidths=0.5, alpha=0.85, zorder=4)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.6, 16)
    ax.set_ylim(55, 340)
    ax.set_xticks([1, 3, 10])
    ax.set_xticklabels(["1", "3", "10"])
    ax.set_yticks([100, 200, 300])
    ax.set_yticklabels(["100", "200", "300"])
    # log scale would otherwise stamp minor-tick labels (e.g. "6x10^1"); suppress
    ax.yaxis.set_minor_formatter(mticker.NullFormatter())
    ax.xaxis.set_minor_formatter(mticker.NullFormatter())
    ax.set_xlabel("Variant effect size")
    ax.set_ylabel("Cells per variant")
    S.despine(ax)
    ax.tick_params(length=2.2)

    # faint gene-cluster labels
    for gene, (lx, ly) in GENE_LABEL_POS.items():
        ax.text(lx, ly, gene, color=S.GENE_COLORS[gene], fontsize=6,
                ha=GENE_LABEL_HA[gene], va="center", alpha=0.9, zorder=5)

    # fill legend: filled=rankable, open=un-rankable (grey, encoding fill only)
    fill_handles = [
        Line2D([0], [0], marker="o", linestyle="none", markersize=3.2,
               markerfacecolor=S.GREY, markeredgecolor=S.GREY,
               markeredgewidth=0.5, label="rankable"),
        Line2D([0], [0], marker="o", linestyle="none", markersize=3.2,
               markerfacecolor="white", markeredgecolor=S.GREY,
               markeredgewidth=0.5, label="un-rankable"),
    ]
    leg = ax.legend(handles=fill_handles, loc="lower left",
                    bbox_to_anchor=(0.0, 0.0),
                    handletextpad=0.3, labelspacing=0.3, borderpad=0.3,
                    labelcolor=S.GREY)

    S.save(fig, "fig6e_design")
    print("drawn: effect size vs cells/variant, fill=rankable; genes", S.GENE_ORDER)


if __name__ == "__main__":
    main()
