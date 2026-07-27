"""Figure 1d (schematic) - a variant is encoded by a 6-dimensional biophysical vector.

Replaces the raster panel ``Fig1d_schematic_highres.pdf`` with a vector panel
drawn at its FINAL composite size (62 x 40 mm, placed at 62 mm, scale 1.0), so
all type lands at 5.5-6.5 pt on the page instead of ~2 pt after down-scaling.

Source file: none (concept panel). No measured quantity is drawn here; the only
numbers shown are the feature-vector dimensionality (6) and the worked-example
variant name TP53 R175H, both fixed by the manuscript text and Fig. 1d caption.
The six feature names are verbatim from the caption: changes in hydrophobicity,
side-chain volume and charge, fold-core location, functional-switch residue
status and hotspot/pathogenic annotation. The companion data panel
``fig1d_theta_pca.py`` supplies the real "shared feature space" scatter that
this schematic's right-hand arrow points into.

One-line message: the model input is a protein-level feature vector, not a
variant identifier, which is what makes held-out allele prediction possible.

Run:  python fig1d_theta_schematic.py
Out:  fig1d_theta_schematic.pdf (+ .png preview)
Vector check:  pdfimages -list fig1d_theta_schematic.pdf | tail -n +3 | wc -l  -> 0
"""
from __future__ import annotations

import os
import sys

from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import nm_style as S

W, H = 62.0, 40.0  # final placement size in mm (composite places at 62 mm)

# The six theta features, verbatim from the Fig. 1d caption, grouped exactly as
# the caption groups them ("changes in hydrophobicity, side-chain volume and
# charge, TOGETHER WITH fold-core location, functional-switch residue status and
# hotspot/pathogenic annotation"). The grouping is what gives the list its
# hierarchy; the accent colour is reserved for the single annotation feature.
FEATURE_GROUPS = [
    ("substitution chemistry",
     [("Δ hydrophobicity", False), ("Δ side-chain volume", False),
      ("Δ charge", False)]),
    ("structural context",
     [("fold-core location", False), ("functional-switch residue", False)]),
    ("annotation",
     [("hotspot / pathogenic", True)]),
]


def arrow(ax, x0, x1, y):
    ax.add_patch(FancyArrowPatch((x0, y), (x1, y), arrowstyle="-|>",
                                 mutation_scale=4.5, lw=0.5, color=S.GREY,
                                 shrinkA=0, shrinkB=0))


def main() -> None:
    S.apply_rcparams()
    fig, ax = S.panel(W, H)
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
    ax.set_xlim(0, W)
    ax.set_ylim(0, H)
    ax.set_aspect("equal")
    ax.axis("off")
    # full-extent white patch pins the tight bounding box to exactly W x H
    ax.add_patch(Rectangle((0, 0), W, H, facecolor="white", edgecolor="none", zorder=0))

    blue = S.GENE_COLORS["TP53"]

    # ---- stage 1: worked example (TP53 R175H) --------------------------------
    # same idiom as panel b (light backbone, coloured domain, orange hotspot
    # lollipop) so the reader recognises the object being encoded
    ax.text(1.0, 25.4, "TP53", ha="left", va="bottom", fontsize=6.5,
            fontweight="bold", color=blue)
    ax.add_patch(Rectangle((1.0, 22.1), 12.0, 2.2, facecolor=S.LIGHT_GREY,
                           edgecolor="none"))
    ax.add_patch(Rectangle((3.4, 22.1), 6.4, 2.2, facecolor=blue, alpha=0.55,
                           edgecolor="none"))
    ax.plot([6.2, 6.2], [21.5, 24.9], color=S.HOTSPOT, lw=0.8, solid_capstyle="butt")
    ax.plot([6.2], [24.9], marker="o", ms=2.0, color=S.HOTSPOT,
            markeredgecolor=S.INK, markeredgewidth=0.3)
    ax.text(7.0, 18.9, "R175H", ha="center", va="center", fontsize=6.0, color=S.INK)
    ax.text(7.0, 15.7, "R → H", ha="center", va="center", fontsize=6.0,
            color=S.GREY)

    arrow(ax, 13.8, 16.8, 19.0)

    # ---- stage 2: the six-dimensional theta vector ---------------------------
    ax.text(18.0, 37.4, "θ   6-dimensional biophysical vector", ha="left",
            va="center", fontsize=6.0, fontweight="bold", color=S.INK)
    x_rule, x_swatch, x_label = 18.0, 19.1, 21.9
    y = 33.6
    for gname, rows in FEATURE_GROUPS:
        ax.text(x_rule, y, gname, ha="left", va="center", fontsize=5.2,
                color=S.GREY)
        y -= 2.9
        y_first = y
        for label, accent in rows:
            col = S.HOTSPOT if accent else S.GREY
            ax.add_patch(Rectangle((x_swatch, y - 0.85), 1.7, 1.7,
                                   facecolor=col, edgecolor="none"))
            ax.text(x_label, y, label, ha="left", va="center", fontsize=6.0,
                    color=S.INK)
            y -= 2.9
        # hairline bracket: groups the rows without drawing a decorative box
        ax.plot([x_rule, x_rule], [y_first + 1.2, y + 1.7], color=S.LIGHT_GREY,
                lw=0.6, solid_capstyle="butt")
        y -= 0.9

    arrow(ax, 49.6, 52.6, 19.0)

    # ---- stage 3: pointer into the data panel (fig1d_theta_pca.pdf) ----------
    ax.text(57.4, 19.0, "shared\nfeature\nspace", ha="center", va="center",
            fontsize=5.5, color=S.INK, linespacing=1.35)

    # ---- footnote ------------------------------------------------------------
    ax.text(1.0, 2.4, "alternative feature spaces tested later: ESM, ESM + θ",
            ha="left", va="center", fontsize=5.5, color=S.GREY)

    S.save(fig, "fig1d_theta_schematic")


if __name__ == "__main__":
    main()
