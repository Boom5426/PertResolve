"""Figure 6c - Dataset-level LODO rankability AUROC forest plot.

Data-direct. Source (committed): results/rankability_predictor_honest.csv via
D.honest_auroc('effect_size'). Each row is a leave-one-dataset-out (LODO) test in
which an effect-size (log_effect) predictor, fit on the other datasets, is asked to
rank held-out perturbations by whether they are natively rankable. Point = AUROC,
whisker = 95% CI (auroc_lo, auroc_hi); dashed line at 0.5 marks chance. TP53 and
KRAS have a single class (nothing rankable) so AUROC is undefined and they are not
plotted; reported only as a grey not-evaluable note.

Run:  python fig6c_auroc.py  ->  fig6c_auroc.pdf (+ .png)
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import nm_style as S
import fig6_data as D

# y-order (top -> bottom on the axis), as specified
ORDER = ["Replogle", "Norman", "Adamson", "GATA1", "JAK1"]


def main() -> None:
    S.apply_rcparams()

    df = D.honest_auroc("effect_size").set_index("held_out")

    # mean over the evaluable datasets actually drawn (matches committed 0.9574 -> 0.96)
    mean_auroc = df.loc[ORDER, "auroc"].mean()

    fig, ax = S.panel(60, 44)

    # y positions: first item at top
    ys = list(range(len(ORDER)))[::-1]

    # reference lines. Same idiom as panel d: dashed grey = null/chance,
    # dotted ink = the mean across the evaluable datasets.
    ax.axvline(0.5, ls=(0, (3, 2)), lw=0.6, color=S.GREY, zorder=1)
    # set beside the line, not on it, so the dashes do not run through the word
    ax.text(0.508, len(ORDER) - 0.52, "chance", fontsize=5.5, color=S.GREY,
            ha="left", va="bottom")
    ax.axvline(mean_auroc, ls=(0, (1, 1.5)), lw=0.6, color=S.INK, zorder=1)
    ax.text(mean_auroc - 0.008, len(ORDER) - 0.52, "mean %.2f" % mean_auroc,
            fontsize=5.5, color=S.INK, ha="right", va="bottom")

    printed = []
    for ds, y in zip(ORDER, ys):
        auroc = float(df.loc[ds, "auroc"])
        lo = float(df.loc[ds, "auroc_lo"])
        hi = float(df.loc[ds, "auroc_hi"])
        col = D.DATASET_COLORS[ds]
        # CI whisker
        ax.plot([lo, hi], [y, y], "-", lw=0.9, color=col, zorder=3,
                solid_capstyle="round")
        # end caps
        for xc in (lo, hi):
            ax.plot([xc, xc], [y - 0.14, y + 0.14], "-", lw=0.9, color=col, zorder=3)
        # point estimate
        ax.scatter([auroc], [y], s=20, color=col, edgecolor="white",
                   linewidths=0.5, zorder=4)
        printed.append((ds, auroc, lo, hi))

    # not-evaluable note for the two single-class allele genes. Carried ONCE for
    # the whole figure (removed from panel b), in the empty lower-left band.
    # Verified-empty band: between the GATA1 and Adamson rows, left of every CI
    # (the leftmost interval bound drawn is Norman's 0.84).
    ax.text(0.525, 1.5,
            "TP53, KRAS: all un-rankable\n(single class, not evaluable)",
            fontsize=5.4, color=S.GREY, ha="left", va="center", linespacing=1.25)

    ax.set_yticks(ys)
    ax.set_yticklabels(ORDER)
    # direct-coloured dataset labels (same device as panel f); together with the
    # single bare key in panel b this is the figure's whole dataset legend.
    for tick, ds in zip(ax.get_yticklabels(), ORDER):
        tick.set_color(D.DATASET_COLORS[ds])
    ax.set_ylim(-0.55, len(ORDER) - 0.15)
    ax.set_xlim(0.45, 1.02)
    ax.set_xticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    ax.set_xlabel("LODO rankability AUROC")
    ax.tick_params(length=2.2)
    S.despine(ax, keep=("left", "bottom"))

    S.save(fig, "fig6c_auroc")

    for ds, a, lo, hi in printed:
        print("%-9s AUROC %.4f [%.4f, %.4f]" % (ds, a, lo, hi))
    print("mean AUROC (drawn) = %.4f" % mean_auroc)


if __name__ == "__main__":
    main()
