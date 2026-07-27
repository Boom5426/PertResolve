"""Figure 5g - The benchmark-validity transition.

Data-direct. Source (committed): results/_remote/unified/controlled_recovery.csv
(loaded via remote_data.controlled()). Each row is a gene x sequencing-depth cell
with the oracle ceiling PDS attainable at that depth and two recovery probabilities
from the controlled-perturbation simulation:
  - P_correct_order : probability the benchmark recovers the full true ranking of
    the graded synthetic predictors ("recover full order").
  - P_winner        : probability it merely selects the single best model
    ("select best model").
Both are plotted against the oracle ceiling (x). The message: below the measurement
window (ceiling near the 0.5 floor, i.e. TP53/KRAS) both probabilities are low and a
leaderboard is unstable; once the ceiling clears ~0.52-0.55 (GATA1) both jump toward
1, and by 0.79+ (JAK1) recovery is essentially certain. Selecting a winner is always
easier than recovering the full order, so its curve leads.

Legend economy: the gene colour key that used to sit inside this panel is removed.
The TP53/KRAS/GATA1/JAK1 colour mapping is already established by direct labels in
panels b, c and f, and the four gene groups are direct-labelled here too, so the
panel now carries exactly one key (the metric key) and the figure carries no
repeated gene legend. Marker shape encodes the metric; fill follows the same rule as
panel b, filled = the primary/limiting quantity, open = the secondary one.

Axis limits and ticks are identical to panel f so the two panels share one
oracle-ceiling axis, but the title and the x-axis label make the unit of analysis
explicit, because g is NOT panel f re-plotted: f is one point per DATASET at a
matched shallow depth (nine benchmarks, one metric), whereas g is one point per
GENE x SEQUENCING DEPTH for the four allele genes (15 cells, depth 25 to 250) and
carries two metrics. Where the two appear to overlap they disagree by construction
(GATA1 resolution 0.512 in f versus P_correct_order 0.989 at depth 50 here).

Run:  python fig5g_validity.py  ->  fig5g_validity.pdf (+ .png)
Vector check:  pdfimages -list fig5g_validity.pdf | tail -n +3 | wc -l  # -> 0
"""
from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import nm_style as S
import remote_data as D

# marker style per metric (shape encodes the metric; colour encodes the gene)
ORDER_C = S.INK           # trend colour for "recover full order"
WINNER_C = S.GREY         # trend colour for "select best model" (house grey)


def _trend(ax, x, y, color):
    """Light monotone-ish guide: sort by x, connect, no markers."""
    o = np.argsort(x)
    ax.plot(np.asarray(x)[o], np.asarray(y)[o], "-", color=color,
            lw=0.6, alpha=0.35, zorder=1)


def main() -> None:
    S.apply_rcparams()
    df = D.controlled().copy()

    # height net of the panel title, so the placed height is unchanged (57.9 mm)
    fig, ax = S.panel(96.9, 61.2)

    # Same two regime bands as panel f, on the same y-axis, in the same row: they
    # are labelled once (in f) and not repeated here, which is the point.
    ax.axhspan(0.90, 1.06, color="#F4F4F4", alpha=1.0, zorder=0)
    ax.axhspan(-0.05, 0.50, color="#EDEDED", alpha=1.0, zorder=0)

    x = df.ceiling_pds.to_numpy()
    # light trend lines (one per metric)
    _trend(ax, x, df.P_winner.to_numpy(), WINNER_C)
    _trend(ax, x, df.P_correct_order.to_numpy(), ORDER_C)

    # points: colour = gene, shape = metric
    for _, r in df.iterrows():
        gc = S.GENE_COLORS[r.gene]
        # "select best model" - open diamond (upper series)
        ax.scatter([r.ceiling_pds], [r.P_winner], s=17, marker="D",
                   facecolor="none", edgecolor=gc, linewidths=0.7, zorder=3)
        # "recover full order" - filled circle
        ax.scatter([r.ceiling_pds], [r.P_correct_order], s=15, marker="o",
                   facecolor=gc, edgecolor="white", linewidths=0.35, zorder=4)

    # regime annotation on the crowded low-ceiling cluster
    ax.annotate("below the measurement window:\nranking unstable",
                xy=(0.512, 0.24), xytext=(0.585, 0.16), fontsize=5.5,
                color=S.GREY, ha="left", va="center", linespacing=1.2,
                arrowprops=dict(arrowstyle="-", lw=0.4, color=S.GREY,
                                shrinkA=1, shrinkB=2))

    # ---- direct gene labels (no gene legend; see module docstring) ----------
    # Each label is placed next to its own group; TP53 and KRAS share a very narrow
    # ceiling range, so both use a short leader line into their own points.
    from matplotlib.lines import Line2D
    lead = dict(arrowstyle="-", lw=0.4, shrinkA=1, shrinkB=2)
    ax.annotate("TP53", xy=(0.4955, 0.135), xytext=(0.464, 0.44), fontsize=6,
                color=S.GENE_COLORS["TP53"], ha="left", va="center",
                arrowprops=dict(color=S.GENE_COLORS["TP53"], **lead))
    ax.annotate("KRAS", xy=(0.5085, 0.600), xytext=(0.532, 0.70), fontsize=6,
                color=S.GENE_COLORS["KRAS"], ha="left", va="center",
                arrowprops=dict(color=S.GENE_COLORS["KRAS"], **lead))
    ax.text(0.556, 0.905, "GATA1", fontsize=6, color=S.GENE_COLORS["GATA1"],
            ha="left", va="center")
    ax.text(0.840, 0.885, "JAK1", fontsize=6, color=S.GENE_COLORS["JAK1"],
            ha="center", va="center")

    # --- metric legend: shape encodes the metric (gene colour neutralised) ---
    leg_h = [
        Line2D([0], [0], marker="D", markerfacecolor="none", markeredgecolor=S.INK,
               markeredgewidth=0.7, markersize=3.6, linestyle="none",
               label="select best model"),
        Line2D([0], [0], marker="o", markerfacecolor=S.INK, markeredgecolor="white",
               markeredgewidth=0.35, markersize=3.6, linestyle="none",
               label="recover full order"),
    ]
    leg = ax.legend(handles=leg_h, loc="lower right", handletextpad=0.4,
                    labelspacing=0.35, borderpad=0.3, fontsize=5.5,
                    labelcolor=S.GREY)
    leg.get_frame().set_linewidth(0.0)

    ax.set_title("Four allele genes, across sequencing depth", fontsize=6.5,
                 loc="left", pad=3)
    ax.set_xlabel("Oracle ceiling per gene and depth")
    ax.set_ylabel("P(recovery): order or winner")
    ax.set_xlim(0.46, 0.98)
    ax.set_ylim(-0.05, 1.06)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_xticks([0.5, 0.6, 0.7, 0.8, 0.9])
    S.despine(ax)
    ax.tick_params(length=2.2)

    S.save(fig, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             "fig5g_validity"))
    print(df[["gene", "depth_m", "ceiling_pds",
              "P_correct_order", "P_winner"]].to_string(index=False))


if __name__ == "__main__":
    main()
