"""Fig. 6g — nested measurement gates define the permissible model claim.

Redrawn at 96 mm to make room for panel h while the figure keeps the 165 mm content
width its other rows use, so every row still shares one left and right edge. This is a relayout, not a rescale: the
millimetres per data unit are unchanged (140 units over 165 mm becomes 81 over 96), so the
type stays the same size relative to the artwork, and the horizontal chain becomes a vertical
one. Compressing the old coordinates into the narrower canvas would have shrunk the artwork
around type that cannot shrink with it.

The full-width policy line the old panel carried at the bottom is dropped rather than
squeezed; the caption already states it.
"""
# Imported export contract: Arial; svg.fonttype: "none"; pdf.fonttype: 42;
# outputs .svg, .pdf and .tiff at dpi=600 via fig6_common.save.
# Final assembled figure target: width_mm = 183.
from __future__ import annotations

W_MM, H_MM = 96.0, 38.0

from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Polygon

import fig6_common as S

GATE_X, GATE_W = 1.5, 27.0          # the decision column
CARD_X, CARD_W = 32.0, 47.5         # the permitted-report column


def _diamond(ax, cx, cy, width, height, *, edge, face="white"):
    points = [(cx, cy + height / 2.0), (cx + width / 2.0, cy),
              (cx, cy - height / 2.0), (cx - width / 2.0, cy)]
    ax.add_patch(Polygon(points, closed=True, facecolor=face, edgecolor=edge,
                         linewidth=0.75, zorder=3))


def _report_card(ax, x, y, width, height, title, body, *, edge, face):
    S.rounded_box(ax, x, y, width, height, edge=edge, face=face, lw=0.7, radius=1.0)
    ax.text(x + 2.0, y + height - 3.0, title, fontsize=5.3, fontweight="bold",
            color=edge, ha="left", va="center")
    ax.text(x + 2.0, y + height - 8.2, body, fontsize=5.0, color=S.DARK_SLATE,
            ha="left", va="center", linespacing=1.35)


def main() -> None:
    S.apply_style()
    fig = S.figure(W_MM, H_MM)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 81)
    ax.set_ylim(0, 68)
    ax.axis("off")

    gate_cx = GATE_X + GATE_W / 2.0

    # A continuous forecast with uncertainty, not a three-bin guarantee.
    S.rounded_box(ax, GATE_X, 46.0, GATE_W, 20.0, edge=S.LIGHT_GREY, face="white",
                  lw=0.7, radius=1.05)
    ax.text(gate_cx, 62.6, "Pilot forecast", fontsize=6.0, fontweight="bold",
            color=S.DARK_SLATE, ha="center", va="center")
    ax.text(gate_cx, 58.4, r"$P$(clear detection floor)", fontsize=5.1,
            color=S.GREY, ha="center", va="center")
    cmap = LinearSegmentedColormap.from_list(
        "workflow_probability", [S.PALE_GREY, S.LIGHT_SLATE, S.DARK_SLATE])
    S.vector_gradient(ax, cmap, (GATE_X + 3.0, GATE_X + GATE_W - 3.0, 51.6, 54.2))
    ax.text(GATE_X + 3.0, 50.0, "low", fontsize=5.0, color=S.GREY, ha="left")
    ax.text(GATE_X + GATE_W - 3.0, 50.0, "high", fontsize=5.0, color=S.GREY, ha="right")
    ax.text(gate_cx, 47.8, "+ uncertainty", fontsize=5.0, color=S.SLATE,
            fontweight="bold", ha="center", va="center")

    # Gate 1: eligibility at the planned depth.
    S.arrow(ax, (gate_cx, 45.6), (gate_cx, 42.6), color=S.MID_GREY)
    _diamond(ax, gate_cx, 33.5, 26.0, 16.0, edge=S.MID_SLATE, face=S.PALE_GREY)
    ax.text(gate_cx, 38.9, "GATE 1 · DETECTION", fontsize=5.1, fontweight="bold",
            color=S.SLATE, ha="center", va="center")
    ax.text(gate_cx, 34.4, "Measurement floor cleared", fontsize=5.2,
            fontweight="bold", color=S.DARK_SLATE, ha="center", va="center")
    ax.text(gate_cx, 31.2, "at the planned depth?", fontsize=5.0,
            color=S.DARK_SLATE, ha="center", va="center")

    # Gate 2: the subset a genuinely allele-level claim needs.
    S.arrow(ax, (gate_cx, 25.4), (gate_cx, 22.4), color=S.MID_SLATE)
    ax.text(gate_cx + 1.4, 23.9, "YES", fontsize=5.0, fontweight="bold",
            color=S.MID_SLATE, ha="left", va="center")
    _diamond(ax, gate_cx, 13.5, 26.0, 15.0, edge=S.DARK_SLATE, face="white")
    ax.text(gate_cx, 18.4, "GATE 2 · IDENTIFICATION", fontsize=5.1, fontweight="bold",
            color=S.SLATE, ha="center", va="center")
    ax.text(gate_cx, 14.2, "Sibling alleles", fontsize=5.2, fontweight="bold",
            color=S.DARK_SLATE, ha="center", va="center")
    ax.text(gate_cx, 11.0, "separable?", fontsize=5.0, color=S.DARK_SLATE,
            ha="center", va="center")

    # Each gate licenses exactly one report, and no gate licenses more than it establishes.
    _report_card(ax, CARD_X, 50.0, CARD_W, 16.0, "ALLELE-LEVEL COMPARISON",
                 "Direction, PDS and DE fidelity;\nmodel ranking with uncertainty",
                 edge=S.DARK_SLATE, face=S.BLUE_PALE)
    S.arrow(ax, (gate_cx + 13.1, 13.5), (CARD_X - 1.0, 13.5), color=S.DARK_SLATE)
    ax.plot([CARD_X - 1.0, CARD_X - 1.0], [13.5, 58.0], color=S.DARK_SLATE, lw=0.65,
            zorder=2)
    S.arrow(ax, (CARD_X - 1.0, 58.0), (CARD_X - 0.2, 58.0), color=S.DARK_SLATE)
    ax.text(gate_cx + 13.9, 11.2, "YES", fontsize=5.0, fontweight="bold",
            color=S.DARK_SLATE, ha="left", va="center")

    _report_card(ax, CARD_X, 30.0, CARD_W, 16.0, "DETECTION-ONLY COMPARISON",
                 "Perturbation against reference only;\nno allele-specific claim",
                 edge=S.MID_SLATE, face="white")
    ax.plot([gate_cx, gate_cx + 19.5], [6.0, 6.0], color=S.MID_SLATE, lw=0.65, zorder=2)
    ax.plot([gate_cx + 19.5, gate_cx + 19.5], [6.0, 38.0], color=S.MID_SLATE, lw=0.65,
            zorder=2)
    S.arrow(ax, (gate_cx + 19.5, 38.0), (CARD_X - 0.2, 38.0), color=S.MID_SLATE)
    ax.text(gate_cx + 1.2, 4.6, "NO", fontsize=5.0, fontweight="bold",
            color=S.MID_SLATE, ha="left", va="center")

    _report_card(ax, CARD_X, 10.0, CARD_W, 16.0, "MEASUREMENT-LIMITED",
                 "Report separately, not as model failure;\nredesign or deepen, then re-pilot",
                 edge=S.AMBER, face=S.AMBER_PALE)
    ax.plot([GATE_X + 0.6, GATE_X + 0.6], [33.5, 18.0], color=S.AMBER, lw=0.7, zorder=2)
    ax.plot([GATE_X + 0.6, CARD_X - 0.2], [18.0, 18.0], color=S.AMBER, lw=0.7, zorder=2)
    S.arrow(ax, (CARD_X - 3.0, 18.0), (CARD_X - 0.2, 18.0), color=S.AMBER)
    ax.text(GATE_X + 1.6, 21.5, "NO", fontsize=5.0, fontweight="bold",
            color=S.AMBER, ha="left", va="center")

    S.save(fig, "fig6g_probability_workflow")


if __name__ == "__main__":
    main()
