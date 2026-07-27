"""Figure 6g - Power-aware reporting protocol (decision half), drawn at final size.

Schematic panel: no data are plotted, so no committed result file is read and no
number is stated.

WHY THIS FILE EXISTS (2026-07-27 rework)
----------------------------------------
The previous panel g was external vector art (Fig6g.pdf, authored on a 1330 pt
canvas and placed at 168 mm, scale 0.36). Two problems, both fixed here:

  1. REDUNDANCY. Its first three boxes ("Pilot screen", "Estimate measurement
     properties", "Predict rankability") were the same three steps already drawn
     in panel a, so a and g told substantially the same story twice and g
     occupied roughly the bottom half of the page to do it. Panel a now stops at
     the predictor OUTPUT and panel g starts at the DECISION, so the two panels
     are sequential: a = what you measure and predict, g = what you then do.
     Panel g's entry chip names panel a rather than redrawing it.

  2. VISUAL LANGUAGE. It used its own font, rounded-box style, a warning triangle,
     a tick-in-circle glyph, and red and green arrows. The house palette has no
     red, and green (#55966B) is JAK1's gene colour in panels b, c and e, so
     neither could be used as an ordinary category colour. This version is drawn
     in matplotlib at exactly 168 x 52 mm (composite scale 1.0) with the same
     box idiom, hairline weight and type sizes as panel a, and uses only:
       - neutral slate / grey ramp for process boxes and the "below floor" path
       - house HOTSPOT amber (#E69F00) for the single decision node
     The Yes / No branches are separated by position, label and box weight, not
     by hue, so the panel survives greyscale.

Wording is held to the Results text (Fig. 6g): triage tool, not a universal
cell-number calculator; "measurement floor" and "measurable perturbations"
throughout; no allele-level identification claim (that boundary condition is
carried by panel a).

Run:  python fig6g_workflow.py  ->  fig6g_workflow.pdf (+ .png)
Vector check:  pdfimages -list fig6g_workflow.pdf | tail -n +3 | wc -l   ->  0
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Polygon

import nm_style as S

W_MM, H_MM = 168.0, 52.0         # exactly the placement size in fig6_assemble.tex
TITLE_PT, BODY_PT, NOTE_PT = 6.4, 5.6, 5.4
BOX_EDGE = "#B9BEC4"             # neutral hairline, identical to panel a
BOX_FILL = "#FAFBFC"
STRIP_FILL = "#F0F1F3"           # same flat strip fill as panel a's caveat band
SLATE = S.FEATURE_COLORS["ESM+theta"]   # dark slate = measurable / proceed
MID = S.FEATURE_COLORS["ESM"]           # mid slate
DEC = S.HOTSPOT                          # #E69F00, the figure's only decision accent
DEC_FILL = "#FDF3E2"

# ---- geometry (mm) --------------------------------------------------------
ENTRY = (1.0, 26.5, 27.0, 34.5)
DIA_CX, DIA_CY, DIA_HW, DIA_HH = 46.0, 30.5, 14.5, 8.5
BOX_A = (68.0, 35.0, 118.0, 51.0)   # Yes branch (3 actions)
BOX_B = (68.0, 10.0, 118.0, 30.0)   # No branch (4 actions)
BOX_C = (126.0, 10.0, 167.0, 51.0)  # converged report
STRIP = (1.0, 0.5, 167.0, 8.0)
FEEDBACK_Y = 14.5


def rbox(ax, x0, y0, x1, y1, *, ec=BOX_EDGE, fc=BOX_FILL, lw=0.5, r=1.1, z=1):
    ax.add_patch(FancyBboxPatch((x0, y0), x1 - x0, y1 - y0,
                                boxstyle=f"round,pad=0,rounding_size={r}",
                                linewidth=lw, edgecolor=ec, facecolor=fc, zorder=z))


def arrow(ax, p0, p1, *, color=S.GREY, lw=0.6, ls="-", scale=4.5, z=2):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle="-|>", mutation_scale=scale,
                                 linewidth=lw, color=color, linestyle=ls,
                                 zorder=z, shrinkA=0, shrinkB=0))


def elbow(ax, x_start, y_start, x_turn, y_end, x_end, *, color=S.GREY, lw=0.6):
    """Horizontal-then-vertical-then-horizontal connector ending in an arrowhead."""
    ax.plot([x_start, x_turn], [y_start, y_start], color=color, lw=lw,
            solid_capstyle="butt", zorder=2)
    ax.plot([x_turn, x_turn], [y_start, y_end], color=color, lw=lw,
            solid_capstyle="butt", zorder=2)
    arrow(ax, (x_turn, y_end), (x_end, y_end), color=color, lw=lw)


def bullets(ax, x, y_top, items, *, dy=2.55, size=BODY_PT, color=S.INK):
    for i, txt in enumerate(items):
        ax.text(x, y_top - dy * i, "• " + txt, fontsize=size,
                va="center", ha="left", color=color, zorder=3)


def main() -> None:
    S.apply_rcparams()
    # drawn to an exact canvas: the tight bbox used by the data panels would crop
    # it and change the composite scale
    plt.rcParams["savefig.bbox"] = None

    fig = plt.figure(figsize=(W_MM * S.MM, H_MM * S.MM))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W_MM)
    ax.set_ylim(0, H_MM)
    ax.axis("off")

    # ---- entry: hands over from panel a ------------------------------------
    rbox(ax, *ENTRY, fc="white")
    ax.text(2.6, 32.0, "Pilot triage (a)", fontsize=TITLE_PT, fontweight="bold",
            va="center", ha="left", color=S.INK)
    ax.text(2.6, 29.1, "predicted rankability", fontsize=BODY_PT,
            va="center", ha="left", color=S.INK)
    arrow(ax, (ENTRY[2], DIA_CY), (DIA_CX - DIA_HW - 0.4, DIA_CY))

    # ---- decision node ------------------------------------------------------
    ax.add_patch(Polygon([[DIA_CX, DIA_CY + DIA_HH], [DIA_CX + DIA_HW, DIA_CY],
                          [DIA_CX, DIA_CY - DIA_HH], [DIA_CX - DIA_HW, DIA_CY]],
                         closed=True, linewidth=0.7, edgecolor=DEC,
                         facecolor=DEC_FILL, zorder=2))
    ax.text(DIA_CX, DIA_CY + 1.5, "Clears the", fontsize=BODY_PT,
            fontweight="bold", va="center", ha="center", color=S.INK, zorder=3)
    ax.text(DIA_CX, DIA_CY - 1.4, "measurement floor?", fontsize=BODY_PT,
            fontweight="bold", va="center", ha="center", color=S.INK, zorder=3)

    a_mid = (BOX_A[1] + BOX_A[3]) / 2.0
    b_mid = (BOX_B[1] + BOX_B[3]) / 2.0
    x_turn = 63.5

    # ---- Yes branch ---------------------------------------------------------
    elbow(ax, DIA_CX + DIA_HW, DIA_CY, x_turn, a_mid, BOX_A[0] - 0.4,
          color=SLATE, lw=0.7)
    ax.text(x_turn + 1.5, a_mid + 3.0, "Yes", fontsize=BODY_PT,
            fontweight="bold", va="center", ha="center", color=SLATE)
    rbox(ax, *BOX_A, ec=SLATE, lw=0.6)
    ax.text(70.4, BOX_A[3] - 2.9, "Proceed to full-depth benchmarking",
            fontsize=TITLE_PT, fontweight="bold", va="center", ha="left",
            color=S.INK)
    bullets(ax, 71.0, BOX_A[3] - 7.5,
            ["size the experiment to the target depth",
             "evaluate direction, PDS and DE fidelity",
             "compare models on measurable perturbations"])

    # ---- No branch ----------------------------------------------------------
    elbow(ax, DIA_CX + DIA_HW, DIA_CY, x_turn, b_mid, BOX_B[0] - 0.4,
          color=S.GREY, lw=0.7)
    ax.text(x_turn + 1.5, b_mid - 3.0, "No", fontsize=BODY_PT,
            fontweight="bold", va="center", ha="center", color=S.GREY)
    rbox(ax, *BOX_B, ec=S.GREY, lw=0.6, fc="white")
    # two-line title: the full phrase does not fit the 50 mm box on one line at
    # 6.4 pt, and shrinking the type below the house 5-7 pt band is not an option
    ax.text(70.4, BOX_B[3] - 2.9, "Redesign, or report as",
            fontsize=TITLE_PT, fontweight="bold", va="center", ha="left",
            color=S.INK)
    ax.text(70.4, BOX_B[3] - 5.8, "measurement-limited",
            fontsize=TITLE_PT, fontweight="bold", va="center", ha="left",
            color=S.INK)
    bullets(ax, 71.0, BOX_B[3] - 9.4,
            ["increase depth where plausible",
             "enrich stronger phenotypes or conditions",
             "change assay, stimulation or cell state",
             "improve variant assignment or controls"], dy=2.45)

    # feedback: redesign sends you back to a new pilot
    ax.plot([BOX_B[0], 14.0], [FEEDBACK_Y, FEEDBACK_Y], color=S.GREY, lw=0.5,
            ls=(0, (2.4, 1.8)), zorder=1)
    ax.add_patch(FancyArrowPatch((14.0, FEEDBACK_Y), (14.0, ENTRY[1] - 0.4),
                                 arrowstyle="-|>", mutation_scale=4.5,
                                 linewidth=0.5, color=S.GREY,
                                 linestyle=(0, (2.4, 1.8)), zorder=1,
                                 shrinkA=0, shrinkB=0))
    ax.text(15.4, FEEDBACK_Y + 1.5, "re-pilot after redesign", fontsize=NOTE_PT,
            va="center", ha="left", color=S.GREY, style="italic")

    # ---- converge on a stratified report ------------------------------------
    arrow(ax, (BOX_A[2], a_mid), (BOX_C[0] - 0.4, a_mid), color=SLATE, lw=0.7)
    arrow(ax, (BOX_B[2], b_mid), (BOX_C[0] - 0.4, b_mid), color=S.GREY, lw=0.7)

    rbox(ax, *BOX_C)
    ax.text(128.4, BOX_C[3] - 3.2, "Report benchmark by", fontsize=TITLE_PT,
            fontweight="bold", va="center", ha="left", color=S.INK)
    ax.text(128.4, BOX_C[3] - 6.3, "measurement regime", fontsize=TITLE_PT,
            fontweight="bold", va="center", ha="left", color=S.INK)

    # sub-block 1: measurable (dark slate accent bar)
    ax.add_patch(FancyBboxPatch((128.4, 32.0), 1.4, 8.3,
                                boxstyle="round,pad=0,rounding_size=0.2",
                                linewidth=0.0, facecolor=SLATE, zorder=2))
    ax.text(131.4, 38.9, "Measurable perturbations", fontsize=BODY_PT,
            fontweight="bold", va="center", ha="left", color=S.INK, zorder=3)
    bullets(ax, 131.4, 36.1, ["model performance and ranking",
                              "uncertainty on the ranking"], dy=2.5, size=NOTE_PT)

    # sub-block 2: below floor (light grey accent bar)
    ax.add_patch(FancyBboxPatch((128.4, 15.5), 1.4, 10.9,
                                boxstyle="round,pad=0,rounding_size=0.2",
                                linewidth=0.0, facecolor=S.LIGHT_GREY, zorder=2))
    ax.text(131.4, 25.0, "Below-floor perturbations", fontsize=BODY_PT,
            fontweight="bold", va="center", ha="left", color=S.GREY, zorder=3)
    bullets(ax, 131.4, 22.2, ["reported separately",
                              "not scored as model failure",
                              "redesign needs documented"],
            dy=2.5, size=NOTE_PT, color=S.INK)

    # ---- principle strip ----------------------------------------------------
    rbox(ax, *STRIP, ec="none", fc=STRIP_FILL, r=0.9)
    ax.text(84.0, 5.5, "Benchmark models only where the ground truth is "
                       "sufficiently resolved.",
            fontsize=BODY_PT, fontweight="bold", va="center", ha="center",
            color=S.INK)
    ax.text(84.0, 2.7, "Predicted rankability is a relative triage tool, not a "
                       "universal depth calculator.",
            fontsize=NOTE_PT, va="center", ha="center", color=S.GREY)

    S.save(fig, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             "fig6g_workflow"))
    print("drawn: decision-half workflow at %.0f x %.0f mm (composite scale 1.0)"
          % (W_MM, H_MM))


if __name__ == "__main__":
    main()
