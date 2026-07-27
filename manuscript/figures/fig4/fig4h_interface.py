"""Figure 4h - Public model interfaces are not uniformly allele-resolution ready.

Interface-compatibility audit (categorical, not a performance panel), encoded from
the authoritative remote inventory (docs/REMOTE_INVENTORY.md) and the manuscript
Methods. Distinguishes performance failure from interface incompatibility and from
optimisation failure. Programmatic (not AI) so model names and the status grid are
exact. Gene-keyed models assign one prediction per gene: resolution coverage
4/470 variants (0.85%).

Legibility: drawn wide (native ~112 x 56 mm) so it sits at composite scale ~1.0
with horizontal column headers; every label is >= 5 pt on the page. Fully vector;
check with: pdfimages -list fig4h_interface.pdf | tail -n +3

Run:  python fig4h_interface.py  ->  fig4h_interface.pdf (+ .png)
"""
from __future__ import annotations

import os
import sys

from matplotlib.patches import Rectangle

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import nm_style as S

# status -> (fill colour, glyph colour).
#
# Palette discipline: the earlier green/red fills are not house colours and red is
# reserved for signed direction or a warning, so the four states are now encoded by
# a monotone lightness ladder drawn from house tokens (dark slate -> HOTSPOT amber
# for the one caveat -> LIGHT_GREY -> near-white). Relative luminance is
# 105 / 162 / 219 / 240, i.e. four separable greys, and the meaning is carried by
# the glyph SHAPE (tick / wave / cross / dot), which is drawn as vector line work,
# so the panel is readable in greyscale and without colour vision.
STATUS = {
    "Y": ("#5F6B76", "white"),      # available            - tick
    "P": (S.HOTSPOT, S.INK),        # available with caveat - wave (amber = warning)
    "N": ("#DBDBDB", S.INK),        # unavailable          - cross
    "NA": ("#F0F0F0", "#8A8A8A"),   # not applicable       - dot
}

COLS = ["Continuous\nvariant input", "Unseen-allele\nconditioning",
        "Allele-specific\noutput", "Ran to\ncompletion", "Allele-level\nscore defined"]

ROWS = [  # (group, model, [5 statuses])
    ("Variant-conditionable", "scGen",      ["Y", "Y", "Y", "Y", "Y"]),
    ("Variant-conditionable", "scVIDR",     ["Y", "Y", "Y", "Y", "Y"]),
    ("Variant-conditionable", "Biolord",    ["Y", "Y", "Y", "Y", "Y"]),
    ("Variant-conditionable", "CellFlow",   ["Y", "Y", "Y", "Y", "Y"]),
    ("Variant-conditionable", "PerturbNet", ["Y", "Y", "Y", "Y", "P"]),
    ("Allele-blind by construction", "scGPT", ["N", "N", "N", "NA", "N"]),
    ("Allele-blind by construction", "GEARS", ["N", "N", "N", "NA", "N"]),
    ("Allele-blind by construction", "STATE", ["N", "N", "N", "Y", "N"]),
    ("Did not converge", "variant-CPA",     ["Y", "Y", "NA", "N", "N"]),
]

XL = -2.05  # left margin (model names + group labels live here)


def _aspect_kx(fig, ax) -> float:
    """x-offset scale that makes a glyph drawn in data units look square on paper.

    The cells are wide and short (many mm per x unit, few per y unit), so raw
    symmetric offsets would render as stretched ticks/crosses.
    """
    fw_mm, fh_mm = (d * 25.4 for d in fig.get_size_inches())
    pos = ax.get_position()
    mm_per_x = (pos.width * fw_mm) / abs(ax.get_xlim()[1] - ax.get_xlim()[0])
    mm_per_y = (pos.height * fh_mm) / abs(ax.get_ylim()[1] - ax.get_ylim()[0])
    return mm_per_y / mm_per_x


def _mark(ax, cx, cy, s, kx=1.0, r=1.0):
    """Font-independent status glyph at (cx, cy); axis y inverted (smaller y = higher).

    ``kx`` compresses the x offsets so the glyph is visually square; ``r`` scales
    the glyph (used to shrink the legend swatches). The glyph colour follows the
    fill lightness so every state keeps a legible contrast.
    """
    w = STATUS[s][1]

    def X(*vals):
        return [cx + v * kx * r for v in vals]

    def Y(*vals):
        return [cy + v * r for v in vals]

    if s == "Y":
        ax.plot(X(-0.16, -0.03, 0.17), Y(0.03, 0.16, -0.15),
                color=w, lw=1.1 * r, solid_capstyle="round",
                solid_joinstyle="round", zorder=5)
    elif s == "N":
        ax.plot(X(-0.14, 0.14), Y(-0.15, 0.15), color=w, lw=1.0 * r, zorder=5)
        ax.plot(X(-0.14, 0.14), Y(0.15, -0.15), color=w, lw=1.0 * r, zorder=5)
    elif s == "P":
        ax.plot(X(-0.19, -0.06, 0.06, 0.19), Y(0.03, 0.13, -0.13, -0.03),
                color=w, lw=1.0 * r, zorder=5)
    else:
        ax.plot([cx], [cy], marker="o", ms=2.2 * r, color=w, zorder=5)


def main() -> None:
    S.apply_rcparams()
    ncol = len(COLS)

    # layout with a gap + header band between groups
    y = 0.0
    layout, headers, prev = [], [], None
    for group, model, sts in ROWS:
        if group != prev:
            if prev is not None:
                y += 0.4
            headers.append((group, y))
            y += 0.75
            prev = group
        layout.append((model, sts, y))
        y += 1.0
    total = y

    fig, ax = S.panel(140, 66)
    ax.set_xlim(XL, ncol + 0.1)
    ax.set_ylim(-1.5, total + 2.3)
    ax.axis("off")
    ax.invert_yaxis()
    kx = _aspect_kx(fig, ax)

    # column headers (horizontal: the wide aspect leaves room, so no rotation)
    for j, c in enumerate(COLS):
        ax.text(j + 0.5, -0.18, c, ha="center", va="bottom", fontsize=5.4,
                color=S.INK, linespacing=1.25)

    # group headers (far-left margin + divider over the grid)
    for group, gy in headers:
        ax.plot([XL, ncol], [gy + 0.62, gy + 0.62], color="#BBBBBB", lw=0.6)
        ax.text(XL, gy + 0.34, group, ha="left", va="center", fontsize=5.4,
                color=S.GREY, style="italic")

    # rows
    for model, sts, ytop in layout:
        ax.text(-0.15, ytop + 0.5, model, ha="right", va="center", fontsize=5.8,
                color=S.INK, fontweight="bold")
        for j, s in enumerate(sts):
            ax.add_patch(Rectangle((j + 0.08, ytop + 0.08), 0.84, 0.84,
                                   facecolor=STATUS[s][0], edgecolor="white", lw=0.5))
            _mark(ax, j + 0.5, ytop + 0.5, s, kx=kx)

    # legend
    ly = total + 0.7
    for k, (s, lab) in enumerate([("Y", "available"), ("P", "caveat"),
                                  ("N", "unavailable"), ("NA", "not applicable")]):
        x = XL + k * 1.42
        ax.add_patch(Rectangle((x, ly), 0.30, 0.55, facecolor=STATUS[s][0],
                               edgecolor="white", lw=0.4))
        _mark(ax, x + 0.15, ly + 0.275, s, kx=kx, r=0.62)
        ax.text(x + 0.40, ly + 0.275, lab, ha="left", va="center", fontsize=5.2,
                color=S.GREY)

    ax.text(XL, total + 1.75,
            "gene-keyed models: one prediction per gene, resolution coverage "
            "4/470 variants (0.85%)", ha="left", va="center", fontsize=5.2, color=S.GREY)

    S.save(fig, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             "fig4h_interface"))
    print(f"drawn interface audit: {len(ROWS)} models x {ncol} criteria")


if __name__ == "__main__":
    main()
