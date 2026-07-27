"""Figure 1f - the six generalization splits.

Replaces the AI raster-derived ``Fig1f.pdf`` (placed at scale 0.18, type 3.4 pt)
with a vector panel drawn at its FINAL composite size (99 x 33 mm, placed at
99 mm, scale 1.0), so all type lands at 5.5-6.5 pt on the page.

Source file: none (design panel). No per-split result or PDS value is shown,
only the split names and one-line definitions, which match the Fig. 1f caption
("random, positional-extrapolation, mechanistic-extrapolation, cross-gene,
low-depth and compatibility settings") and ``fig1f_prompt.md``.

One-line message: generalization is stress-tested along six complementary axes,
each holding out a different slice of the variant set.

Run:  python fig1f_splits.py
Out:  fig1f_splits.pdf (+ .png preview)
Vector check:  pdfimages -list fig1f_splits.pdf | tail -n +3 | wc -l  -> 0
"""
from __future__ import annotations

import os
import sys

from matplotlib.patches import FancyBboxPatch, Rectangle

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import nm_style as S

W, H = 99.0, 33.0  # final placement size in mm

N_CELL = 19          # residue cells in each card's mini protein bar
BAR_W = 19.0         # mm
CARD_W, CARD_H = 31.0, 15.0

# (title, description lines, motif key) -- descriptions verbatim from fig1f_prompt.md
CARDS = [
    ("Random", ["random held-out variants"], "random"),
    ("Positional", ["hold out one protein region", "(C-terminal half)"], "positional"),
    ("Mechanistic", ["hold out hotspot or", "functional-switch residues"], "mechanistic"),
    ("Cross-gene", ["hold out a whole gene,", "train on the others"], "crossgene"),
    ("Low-depth", ["hold out shallow-sampling", "variants"], "lowdepth"),
    ("Compatibility", ["match a prior evaluation", "protocol"], "compat"),
]

HELD = {  # which cells of the mini bar are held out (orange)
    "random": {2, 5, 9, 13, 17},
    "positional": set(range(10, N_CELL)),
    "mechanistic": {3, 4, 10, 15},
}
# deterministic per-cell "depth" for the low-depth motif (design constant, not data)
DEPTHS = [1.0, .55, .85, .35, .95, .70, .30, 1.0, .60, .40,
          .90, .25, .75, 1.0, .45, .80, .35, .95, .65]


def draw_bar(ax, x0, y0, key):
    """Mini protein bar: grey = training variants, orange = held-out variants."""
    cw = BAR_W / N_CELL
    if key == "crossgene":
        # four genes as four separate blocks; the last whole gene is held out
        blk = (BAR_W - 3 * 0.7) / 4
        for i in range(4):
            col = S.HOTSPOT if i == 3 else S.LIGHT_GREY
            ax.add_patch(Rectangle((x0 + i * (blk + 0.7), y0), blk, 1.8,
                                   facecolor=col, edgecolor="none"))
        return
    if key == "lowdepth":
        for i, d in enumerate(DEPTHS):
            h = 0.55 + 1.45 * d
            col = S.HOTSPOT if d <= 0.45 else S.LIGHT_GREY
            ax.add_patch(Rectangle((x0 + i * cw, y0), cw * 0.8, h,
                                   facecolor=col, edgecolor="none"))
        return
    held = HELD.get(key, set())
    for i in range(N_CELL):
        col = S.HOTSPOT if i in held else S.LIGHT_GREY
        ax.add_patch(Rectangle((x0 + i * cw, y0), cw * 0.8, 1.8,
                               facecolor=col, edgecolor="none"))
    if key == "compat":
        ax.add_patch(Rectangle((x0 - 0.5, y0 - 0.6), BAR_W + 1.0, 3.0,
                               facecolor="none", edgecolor=S.GREY, lw=0.5,
                               ls=(0, (1.6, 1.2))))


def main() -> None:
    S.apply_rcparams()
    fig, ax = S.panel(W, H)
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
    ax.set_xlim(0, W)
    ax.set_ylim(0, H)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.add_patch(Rectangle((0, 0), W, H, facecolor="white", edgecolor="none",
                           zorder=0))

    # shared key for the mini-bar motif
    ax.add_patch(Rectangle((0.8, 31.2), 2.4, 1.5, facecolor=S.LIGHT_GREY,
                           edgecolor="none"))
    ax.text(3.9, 31.95, "training variants", ha="left", va="center", fontsize=5.5,
            color=S.GREY)
    ax.add_patch(Rectangle((26.0, 31.2), 2.4, 1.5, facecolor=S.HOTSPOT,
                           edgecolor="none"))
    ax.text(29.1, 31.95, "held-out variants", ha="left", va="center", fontsize=5.5,
            color=S.GREY)

    # No card outlines: alignment and white space carry the 3 x 2 structure, so
    # the six rounded boxes (decoration, not information) are gone.
    for k, (title, desc, key) in enumerate(CARDS):
        r, c = divmod(k, 3)
        x = 0.8 + c * (CARD_W + 2.0)
        y = 14.6 - r * (CARD_H + 1.0)   # notional card bottom edge
        ax.text(x, y + CARD_H - 2.2, title, ha="left", va="center",
                fontsize=6.5, fontweight="bold", color=S.INK)
        draw_bar(ax, x, y + CARD_H - 7.0, key)
        for i, line in enumerate(desc):
            ax.text(x, y + CARD_H - 9.5 - i * 2.4, line, ha="left",
                    va="center", fontsize=5.5, color=S.GREY)

    S.save(fig, "fig1f_splits")


if __name__ == "__main__":
    main()
