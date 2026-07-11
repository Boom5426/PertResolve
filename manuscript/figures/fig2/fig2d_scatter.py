"""Fig 2d: PDS-vs-Pearson scatter (direction without ranking).

Message: in-house feature models recover the *direction* of the perturbation
delta well (mean pearson_delta ~0.55-0.65) yet their per-variant *ranking* skill
(PDS) sits at chance (~0.49-0.52). All learned methods therefore land in a single
"direction without ranking" region just right of PDS=0.5; only the WT-null,
which predicts no delta, drops to pearson_delta=0.

Numbers (all from committed sources, no fabrication):
  x = PDS            from D.definitive()  (results/_remote/unified/definitive_summary.csv)
  y = mean pearson_delta per method, from D.exttheta() (results/results_v4_exttheta.csv)
Points = 18 in-house feature-model heads + Gene-mean + WT-null.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import nm_style as S
import remote_data as D

import numpy as np


def main() -> None:
    S.apply_rcparams()

    pds = D.definitive().set_index("method")["PDS"]
    pear = D.exttheta().groupby("method")["pearson_delta"].mean()

    inhouse = [m for m in pds.index if D.is_inhouse(m)]
    methods = inhouse + ["Gene-mean", "WT-null"]

    # feature-space colour: shared canonical mapping (matches Fig 2b/2c)
    FS_COLOR = dict(S.FEATURE_COLORS)

    fig, ax = S.panel(54.0, 50.0)

    # ---- faint quadrant shading (x split 0.5, y split 0.3) ----------------
    X0, X1 = 0.42, 0.58
    Y0, Y1 = 0.0, 0.75
    XS, YS = 0.50, 0.30
    # top-right quadrant (direction + ranking) faintly highlighted vs others
    ax.axvspan(XS, X1, ymin=(YS - Y0) / (Y1 - Y0), ymax=1.0,
               color=S.LIGHT_GREY, alpha=0.30, lw=0, zorder=0)
    # quadrant divider lines
    ax.axvline(XS, color=S.GREY, lw=0.5, ls=(0, (3, 2)), zorder=1)
    ax.axhline(YS, color=S.GREY, lw=0.5, ls=(0, (3, 2)), zorder=1)

    # ---- points -----------------------------------------------------------
    for m in inhouse:
        c = FS_COLOR[D.feature_space(m)]
        ax.scatter(pds[m], pear[m], s=14, facecolor=c, edgecolor="white",
                   linewidth=0.4, zorder=4)

    # reference methods: distinct grey markers
    ax.scatter(pds["Gene-mean"], pear["Gene-mean"], s=20, marker="D",
               facecolor=S.INK, edgecolor="white", linewidth=0.4, zorder=5)
    ax.scatter(pds["WT-null"], pear["WT-null"], s=20, marker="s",
               facecolor="white", edgecolor=S.INK, linewidth=0.7, zorder=5)

    # ---- region label (in the shaded top-right band) ----------------------
    ax.text(0.565, 0.36, "direction\nwithout\nranking", ha="center", va="center",
            fontsize=6, color=S.INK, linespacing=1.1, zorder=6)

    # ---- callouts for representative methods ------------------------------
    def callout(m, tx, ty, ha):
        ax.annotate(m, (pds[m], pear[m]), (tx, ty),
                    ha=ha, va="center", fontsize=5.5, color=S.INK,
                    arrowprops=dict(arrowstyle="-", lw=0.4, color=S.GREY,
                                    shrinkA=0.0, shrinkB=1.5), zorder=6)

    callout("Ridge-esm", 0.437, 0.610, "left")
    callout("MLP-esm+theta", 0.437, 0.500, "left")
    callout("Lasso-theta", 0.523, 0.690, "left")

    # reference-point labels (no arrow, placed clear of markers)
    ax.text(pds["Gene-mean"] - 0.003, pear["Gene-mean"] + 0.045, "Gene-mean",
            ha="center", va="bottom", fontsize=5.5, color=S.INK, zorder=6)
    ax.text(pds["WT-null"] + 0.006, pear["WT-null"] + 0.005, "WT-null",
            ha="left", va="bottom", fontsize=5.5, color=S.INK, zorder=6)

    # chance-PDS annotation on the divider
    ax.text(XS, Y1, "chance", ha="center", va="bottom", fontsize=5.5,
            color=S.GREY, zorder=6)

    # ---- axes -------------------------------------------------------------
    ax.set_xlim(X0, X1)
    ax.set_ylim(Y0, Y1)
    ax.set_xticks([0.42, 0.46, 0.50, 0.54, 0.58])
    ax.set_yticks([0.0, 0.15, 0.30, 0.45, 0.60, 0.75])
    ax.set_xlabel("PDS (per-variant ranking)")
    ax.set_ylabel(r"mean pearson $\Delta$ (direction)")
    S.despine(ax, keep=("left", "bottom"))

    # ---- feature-space legend (proxy handles, no black edges) -------------
    from matplotlib.lines import Line2D
    handles = [
        Line2D([0], [0], marker="o", linestyle="none", markersize=3.6,
               markerfacecolor=FS_COLOR["theta"], markeredgecolor="white",
               markeredgewidth=0.4, label=r"$\theta$"),
        Line2D([0], [0], marker="o", linestyle="none", markersize=3.6,
               markerfacecolor=FS_COLOR["ESM"], markeredgecolor="white",
               markeredgewidth=0.4, label="ESM"),
        Line2D([0], [0], marker="o", linestyle="none", markersize=3.6,
               markerfacecolor=FS_COLOR["ESM+theta"], markeredgecolor="white",
               markeredgewidth=0.4, label=r"ESM+$\theta$"),
    ]
    leg = ax.legend(handles=handles, loc="lower left", fontsize=5.5,
                    handletextpad=0.4, labelspacing=0.25, borderpad=0.2,
                    handlelength=1.0)
    leg.set_zorder(7)

    S.save(fig, os.path.join(os.path.dirname(os.path.abspath(__file__)), "fig2d_scatter"))


if __name__ == "__main__":
    main()
