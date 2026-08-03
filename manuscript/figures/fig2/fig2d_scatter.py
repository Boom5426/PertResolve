"""Fig. 2d — method-level direction/ranking dissociation.

The scatter is a compact synthesis of panels b/c. The x range matches the PDS
forest rather than tightly magnifying chance-level differences. Learned models
form one quiet cloud; Gene-mean and WT-null are the only direct-labelled
conceptual anchors.

Run: python fig2d_scatter.py

Static QA contract: Arial; svg.fonttype: "none"; pdf.fonttype: 42; outputs
.svg, .pdf, .png and .tiff at dpi=600. Final assembled figure target:
width_mm = 183.
"""
from __future__ import annotations

import os
import sys

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import nm_style as S
import remote_data as D

HERE = os.path.dirname(os.path.abspath(__file__))
W_MM, H_MM = 76.0, 44.0
AX = (12.5, 9.0, 61.5, 31.0)


def canvas():
    fig = plt.figure(figsize=(W_MM * S.MM, H_MM * S.MM))
    bg = fig.add_axes([0, 0, 1, 1])
    bg.set_axis_off()
    bg.set_xlim(0, 1)
    bg.set_ylim(0, 1)
    bg.add_patch(Rectangle((0, 0), 1, 1, facecolor="white",
                           edgecolor="none", zorder=-10))
    l, b, w, h = AX
    return fig, fig.add_axes([l / W_MM, b / H_MM, w / W_MM, h / H_MM])


def main() -> None:
    S.apply_rcparams()
    pds = D.definitive().set_index("method")["PDS"]
    pear = D.exttheta().groupby("method")["pearson_delta"].mean()
    inhouse = [m for m in pds.index if D.is_inhouse(m)]

    fig, ax = canvas()
    ax.axvline(0.50, color=S.GREY, lw=0.6, ls=(0, (3, 2)), zorder=1)

    for method in inhouse:
        ax.scatter(
            pds[method], pear[method], s=14,
            facecolor=S.FEATURE_COLORS[D.feature_space(method)],
            edgecolor="white", linewidth=0.4, alpha=0.88, zorder=4,
        )

    # Key conceptual references: larger and outlined, but still neutral.
    ax.scatter(pds["Gene-mean"], pear["Gene-mean"], s=23, facecolor="white",
               edgecolor=S.INK, linewidth=0.8, zorder=6)
    ax.scatter(pds["WT-null"], pear["WT-null"], s=23, facecolor="white",
               edgecolor=S.GREY, linewidth=0.8, zorder=6)

    ax.annotate("Gene-mean", (pds["Gene-mean"], pear["Gene-mean"]),
                (0.468, 0.676), ha="left", va="center", fontsize=5.6,
                arrowprops=dict(arrowstyle="-", lw=0.45, color=S.GREY,
                                shrinkA=1, shrinkB=2), zorder=7)
    ax.annotate("18 feature-based predictors",
                (float(pds[inhouse].median()), float(pear[inhouse].median())),
                (0.535, 0.675), ha="right", va="center", fontsize=5.5,
                color=S.INK,
                arrowprops=dict(arrowstyle="-", lw=0.45, color=S.GREY,
                                shrinkA=1, shrinkB=2), zorder=7)
    ax.text(pds["WT-null"] + 0.003, pear["WT-null"] + 0.015, "WT-null",
            ha="left", va="bottom", fontsize=5.5, color=S.GREY, zorder=7)

    ax.text(0.463, 0.285, "direction recovered;\nallele ranking at chance",
            ha="left", va="center", fontsize=5.8, color=S.INK,
            fontstyle="italic", linespacing=1.3,
            bbox=dict(facecolor="white", edgecolor="none", pad=0.8), zorder=7)
    ax.text(0.4985, 0.045, "chance", ha="right", va="bottom",
            fontsize=5.3, color=S.GREY)

    ax.set_xlim(0.455, 0.550)
    ax.set_ylim(-0.03, 0.72)
    ax.set_xticks([0.46, 0.50, 0.54])
    ax.set_yticks([0.0, 0.2, 0.4, 0.6])
    ax.set_xlabel(r"PDS$_{cos}$ (allele discrimination)", labelpad=1.5)
    ax.set_ylabel(r"Pearson-$\delta$ (direction recovery)", labelpad=1.5)
    ax.tick_params(pad=1.5)
    S.despine(ax, keep=("left", "bottom"))

    S.save(fig, os.path.join(HERE, "fig2d_scatter"), exact=True,
           formats=("pdf", "png", "svg", "tiff"))


if __name__ == "__main__":
    main()
