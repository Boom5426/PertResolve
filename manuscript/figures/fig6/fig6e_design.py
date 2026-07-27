"""Figure 6e - Effect-size-conditioned design landscape.

Data-direct. Source: results/second_probe_rankability_table.csv via
D.native_rankability(gene), one row per allele-gene perturbation at its deepest
split-half bin. Plots variant effect size (x, log) against the number of cells
scored at that bin (y, log), coloured by gene, with fill encoding rankability:
filled = rankable, open = un-rankable.

Message: rankability is set by *effect size*, not by a universal cell cutoff.
High-effect JAK1 variants are rankable even at low depth, while low-effect
TP53/KRAS variants stay un-rankable at the deepest rung; GATA1 is intermediate.
No cell-number cutoff line is drawn and no contour is fabricated.

Readability rework (2026-07-27), no plotted value moved:
  * The deepest split-half rung for the allele genes is n_work = 150, i.e.
    n_cells = 300, so 341 of the 464 perturbations (73%) sit exactly on y = 300.
    That censoring is now stated on the y-axis label and marked with an explicit
    hairline cap guide, instead of appearing as an unexplained dense line.
  * Un-rankable points are drawn small and semi-transparent so the cap band reads
    as a density band; the rankable points (the panel's actual message) are drawn
    larger and opaque on top. No jitter is applied: every x and y is exact.
  * Gene names moved out of the point cloud into a header band above the cap,
    each over its gene's median effect size with a hairline leader, so the
    overlapping TP53 / KRAS labels no longer collide with the data or each other.

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

CAP = 300.0          # deepest split-half rung for the allele genes (n_work = 150)
Y_LAB = 392.0        # header band: gene name baseline
LEAD_TOP, LEAD_BOT = 358.0, 316.0   # hairline leader from the name down to the cap

# horizontal anchoring of each gene name in the header band. TP53 and KRAS have
# overlapping effect-size ranges (medians 1.41 and 1.11), so they are anchored on
# opposite sides of their own leader; GATA1 and JAK1 are anchored on the side that
# keeps them inside the axes.
GENE_LABEL_HA = {"KRAS": "right", "TP53": "left", "GATA1": "left", "JAK1": "right"}


def main() -> None:
    S.apply_rcparams()
    fig, ax = S.panel(58, 50)

    medians, counts = {}, {}
    for gene in S.GENE_ORDER:
        n = D.native_rankability(gene)
        es = n["effect_size"].to_numpy(float)
        nc = n["n_cells"].to_numpy(float)
        rank = n["rankable"].to_numpy(bool)
        good = np.isfinite(es) & np.isfinite(nc) & (es > 0) & (nc > 0)
        es, nc, rank = es[good], nc[good], rank[good]
        col = S.GENE_COLORS[gene]
        medians[gene] = float(np.median(es))
        counts[gene] = (int(rank.sum()), int(rank.size))

        # un-rankable: small, open, semi-transparent -> the cap line reads as a
        # density band rather than an opaque bar
        m = ~rank
        if m.any():
            ax.scatter(es[m], nc[m], s=3.6, facecolors="none", edgecolors=col,
                       linewidths=0.3, alpha=0.5, zorder=3)
        # rankable: larger, opaque, white-rimmed, on top
        m = rank
        if m.any():
            ax.scatter(es[m], nc[m], s=13, facecolors=col, edgecolors="white",
                       linewidths=0.4, alpha=1.0, zorder=5)

    # explicit censoring guide: the analysis cap, not a trend
    ax.axhline(CAP, ls=(0, (1.2, 1.6)), lw=0.4, color=S.GREY, zorder=1)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.62, 17.0)
    ax.set_ylim(55, 480)
    ax.set_xticks([1, 3, 10])
    ax.set_xticklabels(["1", "3", "10"])
    ax.set_yticks([100, 200, 300])
    ax.set_yticklabels(["100", "200", "300"])
    # log scale would otherwise stamp minor-tick labels (e.g. "6x10^1"); suppress
    ax.yaxis.set_minor_formatter(mticker.NullFormatter())
    ax.xaxis.set_minor_formatter(mticker.NullFormatter())
    # minor ticks would otherwise stamp marks beside the empty header band
    ax.yaxis.set_minor_locator(mticker.NullLocator())
    ax.set_xlabel("Variant effect size")
    ax.set_ylabel("Cells scored per variant\n(deepest split-half rung, cap 300)")
    S.despine(ax)
    # the left spine must not run up through the empty header band
    ax.spines["left"].set_bounds(55, CAP)
    ax.tick_params(length=2.2)

    # direct gene labels in the header band, each over its own median effect size
    for gene in S.GENE_ORDER:
        col = S.GENE_COLORS[gene]
        mx = medians[gene]
        ax.text(mx, Y_LAB, gene, color=col, fontsize=6,
                ha=GENE_LABEL_HA[gene], va="center", zorder=6)
        ax.plot([mx, mx], [LEAD_BOT, LEAD_TOP], color=col, lw=0.4, alpha=0.9,
                solid_capstyle="butt", zorder=5)

    # fill key: the panel's only legend (the gene names above are direct labels).
    # Placed in the verified-empty lower-left block: no perturbation has an effect
    # size below 1.72 together with fewer than 190 scored cells.
    fill_handles = [
        Line2D([0], [0], marker="o", linestyle="none", markersize=3.0,
               markerfacecolor=S.GREY, markeredgecolor="white",
               markeredgewidth=0.4, label="rankable"),
        Line2D([0], [0], marker="o", linestyle="none", markersize=2.4,
               markerfacecolor="none", markeredgecolor=S.GREY,
               markeredgewidth=0.4, label="un-rankable"),
    ]
    ax.legend(handles=fill_handles, loc="lower left", bbox_to_anchor=(-0.01, -0.01),
              handletextpad=0.3, labelspacing=0.3, borderpad=0.2,
              # 6.0 pt here, not 5.5: this panel is placed at ~0.91 composite scale
              labelcolor=S.GREY, fontsize=6.0, frameon=False)

    S.save(fig, "fig6e_design")
    for g in S.GENE_ORDER:
        print("%-6s rankable %d/%d  median effect size %.2f"
              % (g, counts[g][0], counts[g][1], medians[g]))


if __name__ == "__main__":
    main()
