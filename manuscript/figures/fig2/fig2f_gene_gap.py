"""Figure 2f - Per-gene direction-vs-ranking dissociation gap.

Data-direct. Source: results/results_v4_exttheta.csv via remote_data.exttheta(),
restricted to the 18 in-house feature-model heads (remote_data.is_inhouse).
Per gene: gap = mean(pearson_delta) - mean(PDS_cos), pooled over the six splits
and all in-house predictors. Values:
  TP53 0.31, KRAS 0.19, GATA1 0.13, JAK1 -0.04.
Message: for TP53/KRAS/GATA1 predicted directions correlate well with truth yet
variant ranking (PDS_cos) lags far behind; for JAK1 the gap closes (~0), i.e. the
direction signal no longer outruns rankability.

Nature Methods pass: the gene colours are direct-labelled on the x axis (bold,
coloured tick labels) so no legend is needed, and the y label now names the two
quantities being differenced instead of paraphrasing them.

Run:  python fig2f_gene_gap.py   ->  fig2f_gene_gap.pdf (+ .png)
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

W_MM, H_MM = 46.0, 44.0
AX = (13.0, 8.0, 31.0, 32.0)   # left, bottom, width, height in mm


def canvas():
    fig = plt.figure(figsize=(W_MM * S.MM, H_MM * S.MM))
    bg = fig.add_axes([0, 0, 1, 1])
    bg.set_axis_off()
    bg.set_xlim(0, 1)
    bg.set_ylim(0, 1)
    bg.add_patch(Rectangle((0, 0), 1, 1, facecolor="white", edgecolor="none",
                           zorder=-10))
    l, b, w, h = AX
    return fig, fig.add_axes([l / W_MM, b / H_MM, w / W_MM, h / H_MM])


def main() -> None:
    S.apply_rcparams()
    # Keep Greek/maths glyphs in the same Arial-metric sans as the body text.
    # nm_style sets the text font but not mathtext, whose default (DejaVu Sans)
    # would embed a second typeface for every $\theta$, $\delta$ and subscript.
    plt.rcParams.update({"mathtext.fontset": "custom",
                         "mathtext.rm": "Liberation Sans",
                         "mathtext.it": "Liberation Sans:italic",
                         "mathtext.bf": "Liberation Sans:bold",
                         "mathtext.default": "it"})

    df = D.exttheta()
    df = df[df["method"].apply(D.is_inhouse)]

    genes = S.GENE_ORDER  # TP53 > KRAS > GATA1 > JAK1
    gaps = {}
    for g in genes:
        sub = df[df["gene"] == g]
        gaps[g] = sub["pearson_delta"].mean() - sub["PDS_cos"].mean()

    fig, ax = canvas()

    ax.axhline(0.0, color=S.GREY, lw=0.6, zorder=1)

    xs = range(len(genes))
    for i, g in enumerate(genes):
        v = gaps[g]
        ax.bar(i, v, width=0.62, color=S.GENE_COLORS[g], edgecolor="none", zorder=3)
        # value label: above positive bars, below the (small) negative bar
        if v >= 0:
            ax.text(i, v + 0.010, f"{v:.2f}", ha="center", va="bottom",
                    fontsize=6, color=S.INK)
        else:
            ax.text(i, v - 0.010, f"{v:.2f}", ha="center", va="top",
                    fontsize=6, color=S.INK)

    ax.set_xticks(list(xs))
    ax.set_xticklabels(genes, fontsize=5.5)
    for t, g in zip(ax.get_xticklabels(), genes):
        t.set_color(S.GENE_COLORS[g])
        t.set_fontweight("bold")

    ax.set_xlim(-0.65, len(genes) - 0.35)
    ax.set_ylim(-0.095, 0.355)
    ax.set_yticks([0.0, 0.1, 0.2, 0.3])
    ax.set_ylabel(r"Pearson-$\delta$ minus PDS$_{cos}$", labelpad=1.5)

    # annotate the closing gap for JAK1 (message, not a threshold)
    ax.annotate("gap near 0", xy=(3, gaps["JAK1"]), xytext=(3, 0.048),
                ha="center", va="bottom", fontsize=5.5, color=S.GREY,
                arrowprops=dict(arrowstyle="-", lw=0.5, color=S.GREY,
                                shrinkA=0, shrinkB=2))

    S.despine(ax)
    ax.tick_params(length=2.2, pad=1.5)

    S.save(fig, os.path.join(HERE, "fig2f_gene_gap"))

    for g in genes:
        print(f"{g:6} gap={gaps[g]:+.4f}")


if __name__ == "__main__":
    main()
