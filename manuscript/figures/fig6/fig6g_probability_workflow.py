"""Fig. 6g — nested measurement gates define the permissible model claim."""
# Imported export contract: Arial; svg.fonttype: "none"; pdf.fonttype: 42;
# outputs .svg, .pdf and .tiff at dpi=600 via fig6_common.save.
# Final assembled figure target: width_mm = 183.
from __future__ import annotations

W_MM, H_MM = 165, 38.0

import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Polygon

import fig6_common as S


def _diamond(ax, cx, cy, width, height, *, edge, face="white"):
    points = [
        (cx, cy + height / 2.0),
        (cx + width / 2.0, cy),
        (cx, cy - height / 2.0),
        (cx - width / 2.0, cy),
    ]
    patch = Polygon(points, closed=True, facecolor=face, edgecolor=edge,
                    linewidth=0.75, zorder=3)
    ax.add_patch(patch)
    return patch


def _report_card(ax, x, y, width, height, title, body, *, edge, face):
    S.rounded_box(ax, x, y, width, height, edge=edge, face=face,
                  lw=0.7, radius=0.95, zorder=3)
    ax.text(x + 2.0, y + height - 2.4, title, fontsize=5.5,
            fontweight="bold", color=edge, ha="left", va="center", zorder=5)
    ax.text(x + 2.0, y + height - 5.3, body, fontsize=5.0,
            color=S.DARK_SLATE, ha="left", va="top", linespacing=1.15,
            zorder=5)


def main() -> None:
    S.apply_style()
    fig = S.figure(W_MM, H_MM)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 140)
    ax.set_ylim(0, 68)
    ax.axis("off")

    # Input: a continuous forecast with uncertainty, not a three-bin guarantee.
    S.rounded_box(ax, 3.0, 27.0, 20.5, 21.0, edge=S.LIGHT_GREY,
                  face="white", lw=0.7, radius=1.05)
    ax.text(13.25, 44.1, "Pilot forecast", fontsize=6.3, fontweight="bold",
            color=S.DARK_SLATE, ha="center", va="center")
    ax.text(13.25, 39.8, r"$P$(clear detection floor)", fontsize=5.3,
            color=S.GREY, ha="center", va="center")
    cmap = LinearSegmentedColormap.from_list(
        "workflow_probability", [S.PALE_GREY, S.LIGHT_SLATE, S.DARK_SLATE]
    )
    # Vector ramp, not imshow: see fig6a.
    S.vector_gradient(ax, cmap, (6.0, 20.5, 33.4, 36.0))
    ax.text(6.0, 31.8, "low", fontsize=5.0, color=S.GREY, ha="left")
    ax.text(20.5, 31.8, "high", fontsize=5.0, color=S.GREY, ha="right")
    ax.text(13.25, 29.2, "+ uncertainty", fontsize=5.0, color=S.SLATE,
            fontweight="bold", ha="center", va="center")

    # Gate 1: detection-level eligibility at the planned depth.
    ax.text(39.0, 49.6, "GATE 1 · DETECTION", fontsize=5.2,
            fontweight="bold", color=S.SLATE, ha="center", va="center")
    _diamond(ax, 39.0, 37.0, 24.0, 16.0, edge=S.MID_SLATE,
             face=S.PALE_GREY)
    ax.text(39.0, 39.1, "Measurement floor", fontsize=5.45,
            fontweight="bold", color=S.DARK_SLATE, ha="center", va="center")
    ax.text(39.0, 36.0, "cleared at planned depth?", fontsize=5.0,
            color=S.DARK_SLATE, ha="center", va="center")
    ax.text(39.0, 32.9, "prespecified criterion", fontsize=5.0,
            color=S.GREY, ha="center", va="center")
    S.arrow(ax, (24.2, 37.0), (26.5, 37.0), color=S.MID_GREY)

    # Detection-cleared perturbations form the parent set for allele identification.
    S.rounded_box(ax, 55.0, 21.5, 82.0, 34.5, edge=S.LIGHT_GREY,
                  face="#F8FAFB", lw=0.6, radius=1.15, zorder=0)
    ax.text(57.0, 53.8, "DETECTION-CLEARED PERTURBATIONS",
            fontsize=5.2, fontweight="bold", color=S.SLATE,
            ha="left", va="center")
    ax.text(116.0, 53.8, "PERMITTED REPORT", fontsize=5.0,
            fontweight="bold", color=S.GREY, ha="center", va="center")

    # Gate 2: the subset needed for a genuinely allele-level claim.
    ax.text(74.0, 49.3, "GATE 2 · IDENTIFICATION", fontsize=5.2,
            fontweight="bold", color=S.SLATE, ha="center", va="center")
    _diamond(ax, 74.0, 37.0, 25.0, 15.0, edge=S.DARK_SLATE,
             face="white")
    ax.text(74.0, 38.7, "Sibling alleles", fontsize=5.45,
            fontweight="bold", color=S.DARK_SLATE, ha="center", va="center")
    ax.text(74.0, 35.4, "separable?", fontsize=5.2,
            color=S.DARK_SLATE, ha="center", va="center")
    S.arrow(ax, (51.5, 37.0), (61.0, 37.0), color=S.MID_SLATE)
    ax.text(56.2, 38.8, "YES", fontsize=5.0, fontweight="bold",
            color=S.MID_SLATE, ha="center", va="center")

    # Terminal reports: three measurement regimes, not three probability bins.
    _report_card(
        ax,
        96.0,
        39.2,
        38.5,
        12.2,
        "ALLELE-LEVEL COMPARISON",
        "Direction, PDS and DE fidelity\nmodel ranking + uncertainty",
        edge=S.DARK_SLATE,
        face=S.BLUE_PALE,
    )
    S.arrow(ax, (86.7, 40.8), (95.2, 45.0), color=S.DARK_SLATE)
    ax.text(90.6, 44.0, "YES", fontsize=5.0, fontweight="bold",
            color=S.DARK_SLATE, ha="center", va="center")

    _report_card(
        ax,
        96.0,
        23.2,
        38.5,
        12.2,
        "DETECTION-ONLY COMPARISON",
        "Perturbation vs reference only\nno allele-specific claim",
        edge=S.MID_SLATE,
        face="white",
    )
    ax.plot([74.0, 74.0, 92.0], [29.4, 27.8, 27.8], color=S.MID_SLATE,
            lw=0.65, zorder=2)
    S.arrow(ax, (92.0, 27.8), (95.2, 29.3), color=S.MID_SLATE)
    ax.text(78.2, 27.0, "NO", fontsize=5.0, fontweight="bold",
            color=S.MID_SLATE, ha="center", va="center")

    _report_card(
        ax,
        57.5,
        8.0,
        77.0,
        10.5,
        "MEASUREMENT-LIMITED",
        "Report separately · not model failure\nredesign or deepen, then re-pilot",
        edge=S.AMBER,
        face=S.AMBER_PALE,
    )
    ax.plot([39.0, 39.0, 54.0], [28.7, 13.2, 13.2], color=S.AMBER,
            lw=0.7, zorder=2)
    S.arrow(ax, (54.0, 13.2), (56.7, 13.2), color=S.AMBER)
    ax.text(41.0, 25.3, "NO / UNRESOLVED", fontsize=5.0,
            fontweight="bold", color=S.AMBER, ha="left", va="center")

    # Final policy statement.
    S.rounded_box(ax, 3.0, 1.7, 134.0, 4.3, edge=S.LIGHT_GREY,
                  face=S.PALE_GREY, lw=0.55, radius=0.8)
    ax.text(
        70.0,
        3.85,
        "Measurement resolution defines the valid biological claim - below-floor cases are not model failures.",
        fontsize=5.5,
        color=S.DARK_SLATE,
        fontweight="bold",
        ha="center",
        va="center",
    )

    S.save(fig, "fig6g_probability_workflow")


if __name__ == "__main__":
    main()
