"""Figure 6a - Rankability-prediction pipeline (schematic, drawn at final size).

Schematic panel: no data are plotted, so no committed result file is read and no
number is stated beyond the pilot size ("25 to 50 cells per perturbation"), which
is the pilot depth reported in the Results text and in Fig. 6d.

Replaces the previous AI-generated Fig6a.pdf, which was authored on a 423 mm
canvas and placed at 60 mm (scale 0.14), so its internal type landed at ~2.4 pt
on the page. This version is drawn at exactly 60 x 45 mm, the placement size in
fig6_assemble.tex, so the composite scale is 1.0 and every label is 5.4-6.4 pt
on the page.

Message: a pilot-estimable effect size triages which perturbations are evaluable
before full model comparison. The bottom strip carries the boundary condition
stated in the Results: the triage is detection-level (perturbation versus
reference), which is necessary but not sufficient for allele-level identification.

Scope split with panel g (de-duplication, 2026-07-27): panel a stops at the
predictor OUTPUT (what is estimated and what is predicted). Every downstream
ACTION (the decision node, the proceed / redesign branches and the stratified
report) belongs to panel g, so the two panels are sequential rather than two
drawings of the same workflow.

Palette note: the schematic uses only the neutral slate/grey ramp plus the house
HOTSPOT amber in panel g. The gene hues (including JAK1 green) are reserved for
gene identity in panels b, c and e, so no gene colour carries a non-gene meaning.

Run:  python fig6a_pipeline.py  ->  fig6a_pipeline.pdf (+ .png)
Vector check:  pdfimages -list fig6a_pipeline.pdf | tail -n +3 | wc -l   ->  0
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

import nm_style as S

W_MM, H_MM = 60.0, 45.0          # exactly the width/height used in fig6_assemble.tex
TITLE_PT, BODY_PT, NOTE_PT = 6.4, 5.6, 5.4
BOX_EDGE = "#B9BEC4"             # neutral hairline for process boxes
BOX_FILL = "#FAFBFC"
SLATE = S.FEATURE_COLORS["ESM+theta"]  # dark slate = "likely rankable" end of the scale
CELL_DARK = S.FEATURE_COLORS["ESM"]   # slate: perturbed cells (no gene colour, this is a schematic)
CELL_LIGHT = S.LIGHT_GREY             # light grey: reference / control cells

# hand-placed dot lattice for the two cell clusters (offsets in mm from centre)
DOTS = [(-1.5, 1.1), (0.0, 1.5), (1.5, 1.0), (-1.9, -0.3), (-0.4, 0.1),
        (1.1, -0.2), (2.0, -0.8), (-1.2, -1.3), (0.5, -1.4)]


def rbox(ax, x0, y0, x1, y1, *, ec=BOX_EDGE, fc=BOX_FILL, lw=0.5, r=1.1, z=1):
    ax.add_patch(FancyBboxPatch((x0, y0), x1 - x0, y1 - y0,
                                boxstyle=f"round,pad=0,rounding_size={r}",
                                linewidth=lw, edgecolor=ec, facecolor=fc, zorder=z))


def down_arrow(ax, x, y_top, y_bot):
    ax.add_patch(FancyArrowPatch((x, y_top), (x, y_bot),
                                 arrowstyle="-|>", mutation_scale=4.0,
                                 linewidth=0.5, color=S.GREY, zorder=2,
                                 shrinkA=0, shrinkB=0))


def cell_cluster(ax, cx, cy, colour, ec="none", lw=0.0):
    for dx, dy in DOTS:
        ax.add_patch(plt.Circle((cx + dx, cy + dy), 0.42, facecolor=colour,
                                edgecolor=ec, linewidth=lw, zorder=3))


def main() -> None:
    S.apply_rcparams()
    # schematic is drawn to an exact canvas size, so the tight bounding box that
    # the data panels use would crop it and change the composite scale
    plt.rcParams["savefig.bbox"] = None

    fig = plt.figure(figsize=(W_MM * S.MM, H_MM * S.MM))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W_MM)
    ax.set_ylim(0, H_MM)
    ax.axis("off")

    # ---- step 1: pilot assay ------------------------------------------------
    rbox(ax, 2.0, 35.6, 58.0, 44.4)
    ax.text(4.0, 42.1, "1  Pilot assay", fontsize=TITLE_PT, fontweight="bold",
            va="center", ha="left", color=S.INK)
    ax.text(4.0, 39.2, "25 to 50 cells per perturbation,", fontsize=BODY_PT,
            va="center", ha="left", color=S.INK)
    ax.text(4.0, 36.9, "balanced WT / control", fontsize=BODY_PT,
            va="center", ha="left", color=S.INK)
    cell_cluster(ax, 47.0, 38.6, CELL_DARK)
    cell_cluster(ax, 53.6, 38.6, CELL_LIGHT, ec=S.GREY, lw=0.3)

    down_arrow(ax, 30.0, 35.6, 33.6)

    # ---- step 2: estimate measurement properties ---------------------------
    rbox(ax, 2.0, 22.6, 58.0, 33.4)
    ax.text(4.0, 31.3, "2  Estimate measurement properties", fontsize=TITLE_PT,
            fontweight="bold", va="center", ha="left", color=S.INK)
    for i, txt in enumerate(["variant vs WT effect size",
                             "within-variant split-half noise",
                             "pilot signal-to-noise (single feature)"]):
        ax.text(4.6, 28.4 - 2.25 * i, "• " + txt, fontsize=BODY_PT,
                va="center", ha="left", color=S.INK)

    down_arrow(ax, 30.0, 22.6, 20.6)

    # ---- step 3: predict rankability (OUTPUT only; actions live in panel g) --
    rbox(ax, 2.0, 9.0, 58.0, 20.4)
    ax.text(4.0, 18.4, "3  Predict rankability", fontsize=TITLE_PT,
            fontweight="bold", va="center", ha="left", color=S.INK)

    # predicted-probability scale: light grey (below floor) -> dark slate (rankable).
    # Ordered value ramp, so it survives greyscale; no gene hue, no red/green.
    seg = [(4.6, 17.5, S.LIGHT_GREY), (17.5, 30.4, "#B8BFC7"),
           (30.4, 43.3, "#7B858F"), (43.3, 56.2, SLATE)]
    for x0, x1, fc in seg:
        ax.add_patch(FancyBboxPatch((x0, 14.2), x1 - x0, 1.7,
                                    boxstyle="round,pad=0,rounding_size=0.15",
                                    linewidth=0.0, edgecolor="none",
                                    facecolor=fc, zorder=2))
    ax.text(4.6, 12.6, "likely below floor", fontsize=NOTE_PT, va="center",
            ha="left", color=S.GREY, zorder=3)
    ax.text(56.2, 12.6, "likely rankable", fontsize=NOTE_PT, va="center",
            ha="right", color=SLATE, fontweight="bold", zorder=3)
    ax.text(30.4, 10.7, "predicted probability of clearing the split-half floor",
            fontsize=NOTE_PT, va="center", ha="center", color=S.INK, zorder=3)

    # ---- boundary-condition strip ------------------------------------------
    rbox(ax, 2.0, 0.8, 58.0, 7.6, ec="none", fc="#F0F1F3", r=0.9)
    ax.text(30.0, 5.4, "Detection-level triage: perturbation vs reference,",
            fontsize=NOTE_PT, va="center", ha="center", color=S.INK)
    ax.text(30.0, 2.9, "necessary but not sufficient for allele identification",
            fontsize=NOTE_PT, va="center", ha="center", color=S.INK)

    S.save(fig, "fig6a_pipeline")
    print("drawn: 3-step pilot triage schematic at %.0f x %.0f mm (composite scale 1.0)"
          % (W_MM, H_MM))


if __name__ == "__main__":
    main()
