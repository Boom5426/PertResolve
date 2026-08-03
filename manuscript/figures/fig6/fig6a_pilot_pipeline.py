"""Fig. 6a — a disjoint pilot forecasts planned-depth detection rankability."""
# Imported export contract: Arial; svg.fonttype: "none"; pdf.fonttype: 42;
# outputs .svg, .pdf and .tiff at dpi=600 via fig6_common.save.
# Final assembled figure target: width_mm = 183.
from __future__ import annotations

W_MM, H_MM = 95, 41.0

import numpy as np
from matplotlib.colors import LinearSegmentedColormap

import fig6_common as S


def _step_label(ax, x, number, title):
    ax.text(
        x,
        40.0,
        number,
        fontsize=5.5,
        fontweight="bold",
        color="white",
        ha="center",
        va="center",
        bbox=dict(boxstyle="circle,pad=0.18", fc=S.DARK_SLATE, ec="none"),
        zorder=5,
    )
    ax.text(x + 3.4, 40.0, title, fontsize=6.6, fontweight="bold",
            ha="left", va="center")


def _pilot_cells(ax):
    xs = [6.7, 10.3, 13.9, 17.5]
    for row, y in enumerate([25.9, 22.4]):
        for index, x in enumerate(xs):
            face = S.DARK_SLATE if (index + row) % 2 == 0 else "white"
            ax.scatter(
                x,
                y,
                s=17,
                facecolor=face,
                edgecolor=S.DARK_SLATE,
                linewidth=0.5,
                zorder=4,
            )


def main() -> None:
    S.apply_style()
    fig = S.figure(W_MM, H_MM)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 112)
    # Cards span y 14.0 to 43.2. The lower band held the DETECTION-LEVEL ONLY
    # strip, which was removed once its sentence moved to the caption.
    ax.set_ylim(12.6, 44.6)
    ax.axis("off")

    y0, height = 14.0, 29.2
    cards = [(3.0, 29.0), (35.5, 35.5), (74.5, 34.5)]
    for x0, width in cards:
        S.rounded_box(ax, x0, y0, width, height, edge=S.LIGHT_GREY,
                      face="white", lw=0.7, radius=1.1)
    S.arrow(ax, (32.7, 28.7), (34.4, 28.7), color=S.MID_GREY)
    S.arrow(ax, (71.7, 28.7), (73.4, 28.7), color=S.MID_GREY)

    # Step 1 — the future endpoint is physically separated from the pilot.
    _step_label(ax, 5.0, "1", "Separate cell pools")
    ax.text(6.0, 35.0, "Pilot cells", fontsize=5.7, fontweight="bold",
            color=S.DARK_SLATE, ha="left", va="center")
    ax.text(6.0, 29.3, "per perturbation", fontsize=5.0,
            color=S.GREY, ha="left", va="center")
    _pilot_cells(ax)
    ax.text(13.6, 17.3, "perturbation + reference", fontsize=5.0,
            color=S.GREY, ha="center", va="center")

    S.rounded_box(ax, 22.0, 20.0, 7.4, 15.8, edge=S.MID_GREY,
                  face=S.PALE_GREY, lw=0.65, radius=0.8)
    ax.text(25.7, 32.7, "future", fontsize=5.0, color=S.GREY,
            ha="center", va="center")
    ax.text(25.7, 29.7, "evaluation", fontsize=5.0, color=S.GREY,
            ha="center", va="center")
    ax.text(25.7, 24.0, "LOCKED", fontsize=5.0, fontweight="bold",
            color=S.DARK_SLATE, ha="center", va="center")

    # Step 2 — expose the signal/noise origin of the two validated summaries.
    _step_label(ax, 37.5, "2", "Estimate pilot properties")
    S.rounded_box(ax, 38.5, 27.7, 29.4, 8.2, edge=S.LIGHT_GREY,
                  face=S.PALE_GREY, lw=0.55, radius=0.75)
    ax.text(40.0, 33.7, "variant vs reference", fontsize=5.0,
            color=S.GREY, ha="left", va="center")
    ax.text(40.0, 30.4, r"effect size  $\Vert\hat{\delta}\Vert_2$",
            fontsize=5.9, color=S.DARK_SLATE, ha="left", va="center")
    ax.plot([58.5, 58.5], [29.1, 34.5], color=S.LIGHT_GREY, lw=0.55)
    ax.text(63.1, 31.8, "learned\nfeature", fontsize=5.0, color=S.SLATE,
            fontweight="bold", ha="center", va="center", linespacing=1.05)

    S.rounded_box(ax, 38.5, 18.2, 29.4, 8.2, edge=S.LIGHT_GREY,
                  face="white", lw=0.55, radius=0.75)
    ax.text(40.0, 24.2, "effect / split-half noise", fontsize=5.0,
            color=S.GREY, ha="left", va="center")
    ax.text(40.0, 20.9, "pilot SNR", fontsize=5.9,
            color=S.DARK_SLATE, ha="left", va="center")
    ax.plot([58.5, 58.5], [19.6, 25.0], color=S.LIGHT_GREY, lw=0.55)
    ax.text(63.1, 22.3, "training-\nfree", fontsize=5.0, color=S.SLATE,
            fontweight="bold", ha="center", va="center", linespacing=1.05)


    # Step 3 — continuous, assay-specific planned-depth forecast.
    _step_label(ax, 76.5, "3", "Forecast clearance")
    ax.text(78.0, 30.8, r"$P$(clear detection floor)", fontsize=6.0,
            fontweight="bold", color=S.DARK_SLATE, ha="left", va="center")
    cmap = LinearSegmentedColormap.from_list(
        "probability", [S.PALE_GREY, S.LIGHT_SLATE, S.DARK_SLATE]
    )
    # Vector ramp, not imshow: an image object would be the only raster in an
    # otherwise all-vector figure and would print at its own resolution.
    S.vector_gradient(ax, cmap, (77.0, 106.5, 24.6, 27.6))
    ax.text(77.0, 22.8, "below floor", fontsize=5.0, color=S.GREY, ha="left")
    ax.text(91.75, 22.8, "uncertain", fontsize=5.0, color=S.GREY, ha="center")
    ax.text(106.5, 22.8, "measurable", fontsize=5.0, color=S.GREY, ha="right")

    S.save(fig, "fig6a_pilot_pipeline")


if __name__ == "__main__":
    main()
