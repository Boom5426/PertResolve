"""Figure 1a - the AllelePerturb pipeline (schematic, drawn at final size).

Source file: none (concept panel, no data). Numbers reproduced: none; this panel
carries no measured quantity, no score and no benchmark value.

One-line message: a held-out variant's single cells give the pseudobulk target
delta_v; a model sees only the variant's features, never its identity or cells,
and predicts delta_hat_v, which is scored on two axes; in parallel a single-cell
resolution diagnostic asks whether the ground truth can reward an allele model at
all, and the two feed a power-aware verdict.

Replaces the author-supplied Fig1a_editable.pdf, which failed the house standard
on four counts: 36 of 94 words below the 5 pt floor (min 4.3 pt); an off-palette
blue family (#0072B2 / #00A0C6 / #5AB4E6) whose "variant blue" was not the TP53
house blue; a third font family (Arimo); and a Direction/Discrimination block that
restated panel e. Here the two scoring axes are only NAMED, because panel e owns
their definition, and the panel follows the house schematic convention set by
fig3a_window_def.py: ink and grey only, with the hotspot accent reserved for the
held-out variant and the verdict, so the four gene hues stay exclusive to the
data panels.

Drawn natively at its composite placement size (122.5 x 76 mm), so the composite
scale is ~1.0 and every glyph lands at 5.4-6.6 pt on the page.

Run:  python fig1a_pipeline.py  ->  fig1a_pipeline.pdf (+ .png)
Vector check:  pdfimages -list fig1a_pipeline.pdf | tail -n +3 | wc -l   -> 0
"""
from __future__ import annotations

import os
import sys

import numpy as np
from matplotlib.patches import Ellipse, FancyArrowPatch, FancyBboxPatch, Rectangle

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import nm_style as S

W_MM, H_MM = 122.5, 76.0
RNG = np.random.default_rng(11)

FS_STAGE = 6.3     # bold stage headers
FS_SYM = 6.5       # delta symbols
FS_BODY = 5.7
FS_NOTE = 5.4      # smallest type in the panel, still above the 5 pt floor


def cloud(ax, cx, cy, w, h, *, wild_type=False, sibling=False, n=12):
    """A cell cloud. Ink = the variant, grey dashed = wild type, grey = siblings."""
    if wild_type:
        ec, fc, ls = S.GREY, "white", (0, (2.2, 1.4))
    elif sibling:
        ec, fc, ls = S.GREY, S.LIGHT_GREY, "solid"
    else:
        ec, fc, ls = S.INK, S.LIGHT_GREY, "solid"
    ax.add_patch(Ellipse((cx, cy), w, h, facecolor=fc, edgecolor=ec, lw=0.5,
                         alpha=0.9 if wild_type else 0.6, linestyle=ls, zorder=2))
    t = RNG.uniform(0, 2 * np.pi, n)
    r = np.sqrt(RNG.uniform(0, 0.6, n))
    ax.scatter(cx + r * np.cos(t) * w / 2, cy + r * np.sin(t) * h / 2, s=1.4,
               color=S.GREY if (wild_type or sibling) else S.INK,
               alpha=0.6, linewidths=0, zorder=3)


def profile_bar(ax, x0, y0, w, h, *, seed, accent=False):
    """A short expression-profile strip standing in for a G-vector."""
    rng = np.random.default_rng(seed)
    n = 12
    cw = w / n
    vals = rng.uniform(0.15, 1.0, n)
    for i, v in enumerate(vals):
        shade = S.HOTSPOT if (accent and i in (3, 8)) else S.INK
        ax.add_patch(Rectangle((x0 + i * cw, y0), cw * 0.86, h,
                               facecolor=shade, alpha=0.15 + 0.75 * v,
                               edgecolor="none", zorder=3))
    ax.add_patch(Rectangle((x0, y0), w, h, facecolor="none",
                           edgecolor=S.GREY, lw=0.4, zorder=4))


def arrow(ax, x0, y0, x1, y1, *, color=None, lw=0.7, style="-|>"):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle=style,
                                 mutation_scale=5, lw=lw,
                                 color=color or S.INK, shrinkA=0, shrinkB=0,
                                 zorder=5))


def stage_header(ax, x, y, text):
    ax.text(x, y, text, fontsize=FS_STAGE, color=S.INK, ha="left",
            va="baseline", fontweight="bold")


def main() -> None:
    S.apply_rcparams()
    fig, ax = S.panel(W_MM, H_MM)
    ax.set_xlim(0, W_MM)
    ax.set_ylim(0, H_MM)
    ax.axis("off")
    ax.set_position([0, 0, 1, 1])

    # ================= 1 Measurement =================
    stage_header(ax, 0.5, 71.5, "1  Measurement")
    cloud(ax, 8.5, 63.0, 10.0, 7.0)
    ax.text(15.5, 64.6, "held-out variant,", fontsize=FS_NOTE, color=S.GREY,
            ha="left", va="baseline")
    ax.text(15.5, 61.6, "its single cells", fontsize=FS_NOTE, color=S.GREY,
            ha="left", va="baseline")
    cloud(ax, 8.5, 52.5, 10.0, 7.0, wild_type=True, n=9)
    ax.text(15.5, 52.0, "wild-type cells", fontsize=FS_NOTE, color=S.GREY,
            ha="left", va="baseline")
    # collapse the two cell sets into one averaged profile. The formula sits under
    # the delta_v label, not beside the arrow, where it would run into stage 2.
    arrow(ax, 8.5, 48.4, 8.5, 44.6)
    profile_bar(ax, 0.5, 40.0, 26.0, 3.4, seed=3)
    ax.text(0.5, 36.0, r"$\delta_v$  pseudobulk target", fontsize=FS_SYM,
            color=S.INK, ha="left", va="baseline")
    ax.text(0.5, 32.6, "mean(variant) - mean(WT)", fontsize=FS_NOTE,
            color=S.GREY, ha="left", va="baseline")

    arrow(ax, 27.5, 42.0, 32.0, 42.0)

    # ================= 2 Model =================
    stage_header(ax, 33.0, 71.5, "2  Model")
    # feature chip: the model's ONLY input
    profile_bar(ax, 33.0, 63.5, 18.0, 3.2, seed=7, accent=True)
    ax.text(33.0, 67.8, r"variant features ($\theta$ / ESM)", fontsize=FS_BODY,
            color=S.INK, ha="left", va="baseline")
    ax.text(33.0, 59.4, "not the variant identity,", fontsize=FS_NOTE,
            color=S.GREY, ha="left", va="baseline")
    ax.text(33.0, 56.4, "not its cells", fontsize=FS_NOTE, color=S.GREY,
            ha="left", va="baseline")
    arrow(ax, 42.0, 63.2, 42.0, 53.2)
    ax.add_patch(FancyBboxPatch((33.0, 45.6), 25.0, 7.0,
                                boxstyle="round,pad=0.35,rounding_size=1.0",
                                facecolor="white", edgecolor=S.INK, lw=0.6,
                                zorder=4))
    ax.text(45.5, 48.4, "feature to response\nmodel", fontsize=FS_BODY,
            color=S.INK, ha="center", va="baseline", linespacing=1.35,
            zorder=5)
    arrow(ax, 45.5, 45.3, 45.5, 41.6)
    profile_bar(ax, 33.0, 37.8, 25.0, 3.4, seed=5)
    ax.text(33.0, 33.8, r"$\hat{\delta}_v$  predicted profile", fontsize=FS_SYM,
            color=S.INK, ha="left", va="baseline")

    arrow(ax, 59.5, 39.5, 64.0, 39.5)

    # ================= 3 Evaluation (named only; panel e defines them) =======
    stage_header(ax, 65.0, 71.5, "3  Evaluation")
    ax.add_patch(FancyBboxPatch((65.0, 55.0), 25.5, 12.0,
                                boxstyle="round,pad=0.35,rounding_size=1.0",
                                facecolor="white", edgecolor=S.GREY, lw=0.5,
                                zorder=4))
    ax.text(66.2, 63.4, "Direction", fontsize=FS_BODY, color=S.INK,
            ha="left", va="baseline", zorder=5)
    ax.text(66.2, 60.4, r"(Pearson-$\delta$)", fontsize=FS_NOTE, color=S.GREY,
            ha="left", va="baseline", zorder=5)
    ax.text(66.2, 56.4, "Discrimination (PDS)", fontsize=FS_BODY, color=S.INK,
            ha="left", va="baseline", zorder=5)
    ax.text(65.0, 51.0, "is it nearest its own allele", fontsize=FS_NOTE,
            color=S.GREY, ha="left", va="baseline")
    ax.text(65.0, 48.0, "among its siblings?", fontsize=FS_NOTE, color=S.GREY,
            ha="left", va="baseline")
    ax.text(65.0, 43.4, "both axes defined in e", fontsize=FS_NOTE,
            color=S.GREY, ha="left", va="baseline", style="italic")

    arrow(ax, 91.5, 61.0, 96.0, 61.0)

    # ================= 4 Verdict =================
    stage_header(ax, 97.0, 71.5, "4  Verdict")
    ax.add_patch(FancyBboxPatch((97.0, 48.0), 24.5, 19.0,
                                boxstyle="round,pad=0.35,rounding_size=1.0",
                                facecolor="white", edgecolor=S.HOTSPOT, lw=0.8,
                                zorder=4))
    for i, q in enumerate(["measurable?", "benchmarkable?", "worth modelling?"]):
        yq = 62.0 - i * 4.6
        ax.add_patch(Ellipse((99.6, yq + 0.9), 1.7, 1.7, facecolor="none",
                             edgecolor=S.HOTSPOT, lw=0.6, zorder=5))
        ax.text(101.6, yq, q, fontsize=FS_BODY, color=S.INK, ha="left",
                va="baseline", zorder=5)
    arrow(ax, 109.2, 47.7, 109.2, 43.6)
    ax.text(109.2, 40.0, "expand  /  redesign  /  exclude", fontsize=FS_NOTE,
            color=S.INK, ha="center", va="baseline")

    # ================= 5 Resolution diagnostic (parallel lower track) ========
    ax.plot([0.5, 122.0], [27.5, 27.5], color=S.LIGHT_GREY, lw=0.6, zorder=1)
    stage_header(ax, 0.5, 23.0, "5  Resolution diagnostic (single cell)")
    ax.text(46.0, 23.0, "same single cells, no averaging", fontsize=FS_NOTE,
            color=S.GREY, ha="left", va="baseline")

    # detection: variant vs wild type
    ax.text(0.5, 17.0, "detection: variant vs wild type", fontsize=FS_NOTE,
            color=S.INK, ha="left", va="baseline")
    cloud(ax, 7.0, 8.0, 9.0, 6.4)
    cloud(ax, 22.0, 8.0, 9.0, 6.4, wild_type=True, n=9)
    ax.plot([11.8, 17.2], [8.0, 8.0], color=S.GREY, lw=0.6, ls=(0, (2, 1.6)),
            zorder=4)

    # identification: variant vs siblings
    ax.text(38.0, 17.0, "identification: variant vs siblings", fontsize=FS_NOTE,
            color=S.INK, ha="left", va="baseline")
    cloud(ax, 44.0, 8.0, 9.0, 6.4)
    cloud(ax, 58.0, 8.0, 9.0, 6.4, sibling=True, n=10)
    cloud(ax, 72.0, 8.0, 9.0, 6.4, sibling=True, n=10)
    ax.plot([48.8, 53.2], [8.0, 8.0], color=S.GREY, lw=0.6, ls=(0, (2, 1.6)),
            zorder=4)
    ax.plot([62.8, 67.2], [8.0, 8.0], color=S.GREY, lw=0.6, ls=(0, (2, 1.6)),
            zorder=4)

    # the diagnostic feeds the verdict
    ax.text(84.0, 12.6, "can the ground truth", fontsize=FS_NOTE, color=S.GREY,
            ha="left", va="baseline")
    ax.text(84.0, 9.6, "reward an allele model?", fontsize=FS_NOTE,
            color=S.GREY, ha="left", va="baseline")
    arrow(ax, 109.2, 15.5, 109.2, 36.4, color=S.GREY)

    S.save(fig, "fig1a_pipeline")
    print(f"drawn at {W_MM} x {H_MM} mm (scale ~1.0 in the composite)")


if __name__ == "__main__":
    main()
