"""Figure 4h - Public model interfaces are not uniformly allele-resolution ready.

Interface-compatibility audit (categorical, not a performance panel), encoded from
the authoritative remote inventory (docs/REMOTE_INVENTORY.md) and the manuscript
Methods. Distinguishes performance failure from interface incompatibility and from
optimisation failure. Programmatic (not AI) so model names and the status grid are
exact. Gene-keyed models assign one prediction per gene: resolution coverage
4/470 variants (0.85%).

Run:  python fig4h_interface.py  ->  fig4h_interface.pdf (+ .png)
"""
from __future__ import annotations

import os
import sys

from matplotlib.patches import Rectangle

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import nm_style as S

# status -> fill colour (glyphs are drawn as shapes, font-independent)
STATUS = {
    "Y": "#7FB08A",    # available (green)
    "P": "#E8C468",    # available with caveat (amber)
    "N": "#D98C8C",    # unavailable (red)
    "NA": "#DBDBDB",   # not applicable (grey)
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

XL = -3.0   # left margin (model names + group labels live here)


def _mark(ax, cx, cy, s):
    """Font-independent status glyph at (cx, cy); axis y inverted (smaller y = higher)."""
    w = "white"
    if s == "Y":
        ax.plot([cx - 0.16, cx - 0.03, cx + 0.17], [cy + 0.03, cy + 0.16, cy - 0.15],
                color=w, lw=1.2, solid_capstyle="round", solid_joinstyle="round", zorder=5)
    elif s == "N":
        ax.plot([cx - 0.14, cx + 0.14], [cy - 0.15, cy + 0.15], color=w, lw=1.1, zorder=5)
        ax.plot([cx - 0.14, cx + 0.14], [cy + 0.15, cy - 0.15], color=w, lw=1.1, zorder=5)
    elif s == "P":
        ax.plot([cx - 0.15, cx - 0.05, cx + 0.05, cx + 0.15],
                [cy, cy + 0.09, cy - 0.09, cy], color=w, lw=1.1, zorder=5)
    else:
        ax.plot([cx], [cy], marker="o", ms=2.2, color="#8A8A8A", zorder=5)


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

    fig, ax = S.panel(94, 74)
    ax.set_xlim(XL, ncol + 0.1)
    ax.set_ylim(-0.3, total + 2.3)
    ax.axis("off")
    ax.invert_yaxis()

    # column headers (angled)
    for j, c in enumerate(COLS):
        ax.text(j + 0.35, -0.12, c, ha="left", va="bottom", fontsize=5.2,
                rotation=30, rotation_mode="anchor", color=S.INK)

    # group headers (far-left margin + divider over the grid)
    for group, gy in headers:
        ax.plot([XL, ncol], [gy + 0.62, gy + 0.62], color="#BBBBBB", lw=0.6)
        ax.text(XL, gy + 0.34, group, ha="left", va="center", fontsize=5.2,
                color=S.GREY, style="italic")

    # rows
    for model, sts, ytop in layout:
        ax.text(-0.18, ytop + 0.5, model, ha="right", va="center", fontsize=6,
                color=S.INK, fontweight="bold")
        for j, s in enumerate(sts):
            ax.add_patch(Rectangle((j + 0.08, ytop + 0.08), 0.84, 0.84,
                                   facecolor=STATUS[s], edgecolor="white", lw=0.5))
            _mark(ax, j + 0.5, ytop + 0.5, s)

    # legend
    ly = total + 0.7
    for k, (s, lab) in enumerate([("Y", "available"), ("P", "caveat"),
                                  ("N", "unavailable"), ("NA", "not applicable")]):
        x = XL + k * 1.7
        ax.add_patch(Rectangle((x, ly), 0.34, 0.34, facecolor=STATUS[s],
                               edgecolor="white", lw=0.4))
        _mark(ax, x + 0.17, ly + 0.17, s)
        ax.text(x + 0.44, ly + 0.17, lab, ha="left", va="center", fontsize=5,
                color=S.GREY)

    ax.text(XL, total + 1.75,
            "gene-keyed models: one prediction per gene, resolution coverage "
            "4/470 variants (0.85%)", ha="left", va="center", fontsize=5, color=S.GREY)

    S.save(fig, "fig4h_interface")
    print(f"drawn interface audit: {len(ROWS)} models x {ncol} criteria")


if __name__ == "__main__":
    main()
