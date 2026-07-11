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
WINNER_C = "#B4894A"      # muted amber trend for "select best model"


def _trend(ax, x, y, color):
    """Light monotone-ish guide: sort by x, connect, no markers."""
    o = np.argsort(x)
    ax.plot(np.asarray(x)[o], np.asarray(y)[o], "-", color=color,
            lw=0.6, alpha=0.35, zorder=1)


def main() -> None:
    S.apply_rcparams()
    df = D.controlled().copy()

    fig, ax = S.panel(62, 48)

    # 0.5-chance guide
    ax.axhline(0.5, ls=(0, (3, 2)), lw=0.5, color=S.GREY, zorder=0)

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
    ax.annotate("below the\nmeasurement\nwindow:\nranking unstable",
                xy=(0.505, 0.30), xytext=(0.585, 0.30), fontsize=5,
                color=S.GREY, ha="left", va="center",
                arrowprops=dict(arrowstyle="-", lw=0.4, color=S.GREY,
                                shrinkA=1, shrinkB=2))

    # compact gene colour key (top-left), inline swatch dots
    from matplotlib.lines import Line2D
    gene_h = [Line2D([0], [0], marker="o", markerfacecolor=S.GENE_COLORS[g],
                     markeredgecolor="white", markeredgewidth=0.3, markersize=3.2,
                     linestyle="none", label=g) for g in S.GENE_ORDER]
    gene_leg = ax.legend(handles=gene_h, loc="upper left", ncol=2,
                         handletextpad=0.25, labelspacing=0.3, columnspacing=0.7,
                         borderpad=0.3, fontsize=5.2, bbox_to_anchor=(0.44, 0.86))
    gene_leg.get_frame().set_linewidth(0.0)
    ax.add_artist(gene_leg)

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
                    labelspacing=0.35, borderpad=0.3, fontsize=5.5)
    leg.get_frame().set_linewidth(0.0)

    ax.set_xlabel("Oracle ceiling (measurement window)")
    ax.set_ylabel("P(recovery)")
    ax.set_xlim(0.47, 0.99)
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
