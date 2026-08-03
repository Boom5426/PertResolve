"""Fig. 2f — per-gene direction and allele discrimination shown separately.

The previous raw subtraction Pearson-delta minus PDS combined two metrics with
different meanings and different nulls. This revision uses aligned small
forests: each small point is one of 18 feature-model heads and the diamond is
the across-head mean. Gene here is also dataset/assay context.

Run: python fig2f_gene_gap.py

Static QA contract: Arial; svg.fonttype: "none"; pdf.fonttype: 42; outputs
.svg, .pdf, .png and .tiff at dpi=600. Final assembled figure target:
width_mm = 183.
"""
from __future__ import annotations

import os
import sys

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import nm_style as S
import remote_data as D

HERE = os.path.dirname(os.path.abspath(__file__))
W_MM, H_MM = 58.0, 44.0
AX_BOTTOM, AX_HEIGHT = 7.8, 29.2
PEAR_AX = (12.0, AX_BOTTOM, 19.5, AX_HEIGHT)
PDS_AX = (37.0, AX_BOTTOM, 18.8, AX_HEIGHT)


def canvas():
    fig = plt.figure(figsize=(W_MM * S.MM, H_MM * S.MM))
    bg = fig.add_axes([0, 0, 1, 1])
    bg.set_axis_off()
    bg.set_xlim(0, 1)
    bg.set_ylim(0, 1)
    bg.add_patch(Rectangle((0, 0), 1, 1, facecolor="white",
                           edgecolor="none", zorder=-10))

    def add(rect):
        l, b, w, h = rect
        return fig.add_axes([l / W_MM, b / H_MM, w / W_MM, h / H_MM])

    return fig, add(PEAR_AX), add(PDS_AX)


def main() -> None:
    S.apply_rcparams()
    df = D.exttheta()
    df = df[df.method.map(D.is_inhouse)]
    per_head = (df.groupby(["gene", "method"])[["pearson_delta", "PDS_cos"]]
                  .mean().reset_index())

    genes = S.GENE_ORDER
    y_map = {gene: len(genes) - 1 - i for i, gene in enumerate(genes)}
    jitter = np.linspace(-0.16, 0.16, 18)

    fig, axp, axd = canvas()
    for axis in (axp, axd):
        axis.set_ylim(-0.55, 3.55)
        axis.set_yticks(range(4))
        axis.tick_params(axis="y", length=0)
        for y in (0, 2):
            axis.axhspan(y - 0.46, y + 0.46, color="#F6F7F8", zorder=0)
        S.despine(axis, keep=("bottom",))

    axp.axvline(0, color=S.GREY, lw=0.6, ls=(0, (3, 2)), zorder=1)
    axd.axvline(0.50, color=S.GREY, lw=0.6, ls=(0, (3, 2)), zorder=1)

    for gene in genes:
        sub = per_head[per_head.gene == gene]
        assert len(sub) == 18
        y = y_map[gene]
        color = S.GENE_COLORS[gene]

        axp.scatter(sub.pearson_delta, y + jitter, s=8.5, facecolor=color,
                    edgecolor="white", linewidth=0.25, alpha=0.38, zorder=3)
        axd.scatter(sub.PDS_cos, y + jitter, s=8.5, facecolor=color,
                    edgecolor="white", linewidth=0.25, alpha=0.38, zorder=3)

        mp = float(sub.pearson_delta.mean())
        md = float(sub.PDS_cos.mean())
        axp.scatter([mp], [y], s=26, marker="D", facecolor=color,
                    edgecolor="white", linewidth=0.4, zorder=5)
        axd.scatter([md], [y], s=26, marker="D", facecolor=color,
                    edgecolor="white", linewidth=0.4, zorder=5)
        axp.text(mp, y + 0.27, f"{mp:.2f}", ha="center", va="bottom",
                 fontsize=5.2, color=S.INK,
                 bbox=dict(facecolor="white", edgecolor="none", pad=0.25))
        axd.text(md, y + 0.27, f"{md:.2f}", ha="center", va="bottom",
                 fontsize=5.2, color=S.INK,
                 bbox=dict(facecolor="white", edgecolor="none", pad=0.25))

    axp.set_yticklabels(list(reversed(genes)), fontsize=5.7)
    for tick, gene in zip(axp.get_yticklabels(), reversed(genes)):
        tick.set_color(S.GENE_COLORS[gene])
        tick.set_fontweight("bold")
    axd.set_yticklabels([])

    axp.set_xlim(-0.04, 0.88)
    axp.set_xticks([0.0, 0.4, 0.8])
    axd.set_xlim(0.35, 0.58)
    axd.set_xticks([0.40, 0.50])
    axp.set_xlabel(r"Pearson-$\delta$", labelpad=1.4)
    axd.set_xlabel(r"PDS$_{cos}$", labelpad=1.4)
    axp.tick_params(axis="x", pad=1.3)
    axd.tick_params(axis="x", pad=1.3)

    fig.text((PEAR_AX[0] + PEAR_AX[2] / 2) / W_MM, 0.935,
             "Direction recovery", ha="center", va="center",
             fontsize=5.7, fontweight="bold")
    fig.text((PDS_AX[0] + PDS_AX[2] / 2) / W_MM, 0.935,
             "Allele discrimination", ha="center", va="center",
             fontsize=5.7, fontweight="bold")
    fig.text(0.965, 0.865, "diamonds: mean of 18 heads",
             ha="right", va="center", fontsize=5.0, color=S.GREY)
    axp.text(0.012, 3.40, "null", ha="left", va="top",
             fontsize=5.0, color=S.GREY)
    axd.text(0.502, 3.40, "chance", ha="left", va="top",
             fontsize=5.0, color=S.GREY)

    S.save(fig, os.path.join(HERE, "fig2f_gene_gap"), exact=True,
           formats=("pdf", "png", "svg", "tiff"))


if __name__ == "__main__":
    main()
