"""Fig 2g: representative held-out TP53 variants (Ridge-esm).

Message: direction of the predicted allele effect is stable across variants
(pearson_delta ~0.74-0.82), yet the allele *ranking* score (PDS_cos) spans the
whole 0-1 range. Similar direction, very different allele ranking.

Data: D.exttheta(), method=='Ridge-esm', gene=='TP53', per variant.
A representative subset of 10 distinct variants (deduplicated by variant, one
row each) evenly spaced across the observed PDS_cos range is shown.

Nature Methods pass:
  * The two bands are now two real axes with real y-axis labels; the former
    in-plot strings "pearson_delta (direction)" / "PDS_cos (allele ranking)"
    and the in-plot title (which duplicated the caption) are gone.
  * Both axes carry the SAME 0-1 scale from a true zero baseline, so the visual
    contrast (direction tightly clustered high, ranking spread over the full
    range) is a property of the data and not of two different rescalings. The
    old top band mapped 0.70-0.85 onto a strip, which exaggerated the spread of
    the very quantity the panel calls stable.
  * Re-proportioned to a short 120.5 x 44 mm strip that shares the bottom row
    with Fig 2f, so it supports rather than dominates the page.
"""
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import nm_style as S
import remote_data as D

HERE = os.path.dirname(os.path.abspath(__file__))

N_SHOW = 10
GENE = "TP53"
METHOD = "Ridge-esm"

W_MM, H_MM = 120.5, 44.0
AX_LEFT_MM, AX_WIDTH_MM = 16.0, 101.5
AX_H_MM = 14.5                    # each of the two stacked axes
AX_BOT_MM = 8.0                   # bottom axis (PDS) sits here
AX_TOP_MM = 27.0                  # top axis (Pearson) sits here


def select_variants():
    df = D.exttheta()
    sub = df[(df.method == METHOD) & (df.gene == GENE)].copy()
    # one representative row per distinct variant, sorted by ranking score
    sub = sub.sort_values("PDS_cos").drop_duplicates("variant", keep="first")
    idx = np.linspace(0, len(sub) - 1, N_SHOW).round().astype(int)
    return sub.iloc[idx].reset_index(drop=True)


def canvas():
    fig = plt.figure(figsize=(W_MM * S.MM, H_MM * S.MM))
    bg = fig.add_axes([0, 0, 1, 1])
    bg.set_axis_off()
    bg.set_xlim(0, 1)
    bg.set_ylim(0, 1)
    bg.add_patch(Rectangle((0, 0), 1, 1, facecolor="white", edgecolor="none",
                           zorder=-10))

    def band(bottom_mm):
        return fig.add_axes([AX_LEFT_MM / W_MM, bottom_mm / H_MM,
                             AX_WIDTH_MM / W_MM, AX_H_MM / H_MM])

    return fig, band(AX_TOP_MM), band(AX_BOT_MM)


def lollipops(ax, x, y, color):
    ax.vlines(x, 0, y, color=S.LIGHT_GREY, lw=0.7, zorder=1)
    ax.scatter(x, y, s=12, color=color, edgecolors="white", linewidths=0.35,
               zorder=3, clip_on=False)


def main():
    S.apply_rcparams()
    # Keep Greek/maths glyphs in the same Arial-metric sans as the body text.
    # nm_style sets the text font but not mathtext, whose default (DejaVu Sans)
    # would embed a second typeface for every $\theta$, $\delta$ and subscript.
    plt.rcParams.update({"mathtext.fontset": "custom",
                         "mathtext.rm": "Liberation Sans",
                         "mathtext.it": "Liberation Sans:italic",
                         "mathtext.bf": "Liberation Sans:bold",
                         "mathtext.default": "it"})
    pick = select_variants()
    x = np.arange(len(pick))
    labels = pick["variant"].tolist()
    pdir = pick["pearson_delta"].to_numpy()
    pds = pick["PDS_cos"].to_numpy()

    col = S.GENE_COLORS[GENE]
    fig, ax_dir, ax_pds = canvas()

    lollipops(ax_dir, x, pdir, col)
    lollipops(ax_pds, x, pds, col)

    # chance line on the ranking axis only (Pearson's null is 0, already the base)
    ax_pds.axhline(0.5, color=S.GREY, lw=0.6, ls=(0, (3, 2)), zorder=2)
    ax_pds.text(len(pick) - 0.55, 0.52, "chance", ha="right", va="bottom",
                fontsize=5.5, color=S.GREY)

    for ax in (ax_dir, ax_pds):
        ax.set_xlim(-0.6, len(pick) - 0.4)
        ax.set_ylim(0, 1.04)
        ax.set_yticks([0, 0.5, 1])
        ax.tick_params(axis="y", length=1.8, pad=1.5, labelsize=5.8)
        S.despine(ax, keep=("left", "bottom"))
        ax.spines["left"].set_bounds(0, 1)
        ax.spines["bottom"].set_bounds(-0.6, len(pick) - 0.4)

    ax_dir.set_xticks(x)
    ax_dir.set_xticklabels([])
    ax_dir.tick_params(axis="x", length=0)
    ax_dir.set_ylabel(r"Pearson-$\delta$", labelpad=1.5)

    ax_pds.set_xticks(x)
    ax_pds.set_xticklabels(labels, fontsize=5.8)
    ax_pds.tick_params(axis="x", length=0, pad=1.5)
    ax_pds.set_ylabel(r"PDS$_{cos}$", labelpad=1.5)
    ax_pds.set_xlabel("held-out TP53 variant (Ridge ESM)", labelpad=1.5)

    S.save(fig, os.path.join(HERE, "fig2g_pervariant"))

    # echo numbers plotted
    for _, r in pick.iterrows():
        print(f"{r.variant:>7s}  split={r.split}  "
              f"pearson_delta={r.pearson_delta:.3f}  PDS_cos={r.PDS_cos:.3f}")


if __name__ == "__main__":
    main()
