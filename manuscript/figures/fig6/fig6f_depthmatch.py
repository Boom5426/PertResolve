"""Fig 6f: native vs depth-matched (n=50/half) un-rankable fraction on external
gene-level atlases (Replogle, Adamson, Norman).

Message: depth normalization raises Adamson/Norman toward Replogle, but Replogle
stays high -> a dataset-specific effect-size structure, not depth alone.

Numbers come only from D.canonical() (results/canonical_numbers.json).
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import nm_style as S
import fig6_data as D   # in the same fig6/ dir

import numpy as np


def main() -> None:
    S.apply_rcparams()
    c = D.canonical()
    native = c["unrankable_native_pct"]
    matched = c["unrankable_matched50_pct"]

    # top-to-bottom display order; Replogle highest so it sits on top
    datasets = ["Replogle", "Adamson", "Norman"]
    y_pos = {d: i for i, d in enumerate(reversed(datasets))}  # Norman=0 ... Replogle=2

    fig, ax = S.panel(52.0, 44.0)

    for d in datasets:
        y = y_pos[d]
        col = D.DATASET_COLORS[d]
        xn, xm = native[d], matched[d]

        # connecting dumbbell line
        ax.plot([xn, xm], [y, y], color=col, lw=0.9, zorder=1,
                solid_capstyle="round")
        # native = filled circle; matched-50 = open circle
        ax.plot(xn, y, marker="o", ms=4.6, mfc=col, mec=col, mew=0.6,
                linestyle="none", zorder=3)
        ax.plot(xm, y, marker="o", ms=4.6, mfc="white", mec=col, mew=1.0,
                linestyle="none", zorder=3)

        # value labels: put the lower value to the left of its marker and the
        # higher value to the right, so the two never collide. If the lower
        # marker sits too close to the y-axis (would clash with the tick
        # label), lift its label above the marker instead.
        lo = min(xn, xm)
        hi = max(xn, xm)
        for xval in (lo, hi):
            if xval == lo:
                if xval < 8.0:  # near axis: label above to avoid tick clash
                    ax.text(xval, y + 0.28, f"{xval:.1f}", ha="center",
                            va="bottom", fontsize=5.5, color=S.INK)
                else:
                    ax.text(xval - 2.6, y, f"{xval:.1f}", ha="right",
                            va="center", fontsize=5.5, color=S.INK)
            else:  # higher value -> right of marker
                ax.text(xval + 2.6, y, f"{xval:.1f}", ha="left",
                        va="center", fontsize=5.5, color=S.INK)

    # axes cosmetics
    ax.set_xlim(0, 60)
    ax.set_xticks([0, 20, 40, 60])
    ax.set_ylim(-0.6, len(datasets) - 0.4)
    ax.set_yticks(list(y_pos.values()))
    ax.set_yticklabels([d for d in reversed(datasets)], fontsize=6.5)
    ax.set_xlabel("un-rankable perturbations (%)", fontsize=6.5)
    ax.tick_params(axis="both", length=2.0, pad=1.5)

    S.despine(ax, keep=("left", "bottom"))

    # colour the y tick labels by dataset for immediate association
    for tick, d in zip(ax.get_yticklabels(), reversed(datasets)):
        tick.set_color(D.DATASET_COLORS[d])

    # legend: native (filled) vs depth-matched n=50 (open), drawn in neutral ink
    from matplotlib.lines import Line2D
    handles = [
        Line2D([0], [0], marker="o", linestyle="none", ms=4.4,
               mfc=S.GREY, mec=S.GREY, mew=0.6, label="native depth"),
        Line2D([0], [0], marker="o", linestyle="none", ms=4.4,
               mfc="white", mec=S.GREY, mew=1.0, label="depth-matched (n=50/half)"),
    ]
    leg = ax.legend(handles=handles, loc="lower right",
                    bbox_to_anchor=(1.0, -0.02), fontsize=5.5,
                    handletextpad=0.4, labelspacing=0.35,
                    borderaxespad=0.2, frameon=False)
    for t in leg.get_texts():
        t.set_color(S.INK)

    fig.subplots_adjust(left=0.20, right=0.97, top=0.97, bottom=0.16)
    stem = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fig6f_depthmatch")
    S.save(fig, stem)
    print("native:", {d: native[d] for d in datasets})
    print("matched:", {d: matched[d] for d in datasets})


if __name__ == "__main__":
    main()
