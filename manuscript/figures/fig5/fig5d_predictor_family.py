"""Fig. 5d — build-only interpolation creates a known predictor order."""
from __future__ import annotations

W_MM, H_MM = 130, 31.0

import numpy as np

import fig5_common as S

MEAN = np.array([0.7, 0.4, 0.8, 0.3, 0.6, 0.5, 0.7, 0.4, 0.6, 0.5])
VARIANT = np.array([1.5, -1.0, 0.3, 1.8, -0.7, 1.1, -1.5, 0.7, 1.3, -0.8])


def _predictor_card(ax, x, alpha, index):
    color = S.mix(S.LIGHT_SLATE, S.DARK_SLATE, alpha)
    values = (1.0 - alpha) * MEAN + alpha * VARIANT
    S.rounded_box(ax, x, 14.2, 12.0, 16.5, edge=color,
                  face=S.mix("white", color, 0.055), lw=0.62, radius=0.75)
    ax.text(x + 6.0, 28.5, f"P{index}", fontsize=5.8, fontweight="bold",
            color=color, ha="center", va="center")
    S.mini_profile(ax, x + 1.3, x + 10.7, 22.5, values, color=color,
                   amplitude=2.55, lw=0.72)
    ax.text(x + 6.0, 16.4, rf"$\alpha={alpha:g}$", fontsize=5.15,
            color=S.GREY, ha="center", va="center")


def main() -> None:
    S.apply_style()
    fig = S.figure(W_MM, H_MM)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 50)
    ax.axis("off")

    ax.text(
        1.5,
        44.4,
        "The predictor family is fixed before scoring; five illustrative anchors are shown.",
        fontsize=5.7,
        color=S.GREY,
        ha="left",
        va="top",
    )

    # Build-only source profiles.
    S.rounded_box(ax, 1.5, 10.2, 22.4, 30.6, edge=S.LIGHT_GREY,
                  face="white", lw=0.7, radius=1.1)
    ax.text(3.3, 38.0, "Build-only inputs", fontsize=6.1, fontweight="bold")
    ax.text(3.3, 33.1, r"gene mean  $\bar{\delta}^{build}_g$",
            fontsize=5.45, color=S.MID_SLATE)
    S.mini_profile(ax, 3.6, 21.7, 28.9, MEAN, color=S.LIGHT_SLATE,
                   amplitude=2.35, lw=0.72)
    ax.text(3.3, 22.7, r"variant effect  $\delta^{build}_v$",
            fontsize=5.45, color=S.DARK_SLATE)
    S.mini_profile(ax, 3.6, 21.7, 17.7, VARIANT, color=S.DARK_SLATE,
                   amplitude=3.1, lw=0.75)

    S.arrow(ax, (24.7, 26.0), (27.0, 26.0), color=S.MID_GREY)

    # Formula and explicitly named predictor ladder.
    S.rounded_box(ax, 27.6, 34.1, 70.9, 6.7, edge=S.LIGHT_GREY,
                  face="#EEF1F4", lw=0.6, radius=0.9)
    ax.text(
        63.05,
        37.35,
        r"$\hat{\delta}_{v,\alpha}=(1-\alpha)\bar{\delta}^{build}_g"
        r"+\alpha\delta^{build}_v$",
        fontsize=6.75,
        color=S.INK,
        ha="center",
        va="center",
    )
    ax.text(27.6, 32.1, "Illustrative anchors along the interpolation continuum",
            fontsize=5.25, color=S.GREY, ha="left", va="center")

    alphas = [0.0, 0.25, 0.5, 0.75, 1.0]
    xs = [27.6, 41.8, 56.0, 70.2, 84.4]
    for index, (x, alpha) in enumerate(zip(xs, alphas), start=1):
        _predictor_card(ax, x, alpha, index)

    ax.text(33.6, 12.1, "gene mean", fontsize=5.0, color=S.GREY,
            ha="center", va="center")
    ax.text(90.4, 12.1, "variant-specific", fontsize=5.0, color=S.GREY,
            ha="center", va="center")
    S.arrow(ax, (40.0, 11.8), (82.0, 11.8), color=S.MID_SLATE,
            lw=0.68, mutation=4.7)
    ax.text(61.0, 9.8, "design-defined quality:  P1 < P2 < P3 < P4 < P5",
            fontsize=5.3, fontweight="bold", color=S.DARK_SLATE,
            ha="center", va="center")

    # Leakage guard is separated from the construction itself.
    S.rounded_box(ax, 1.5, 2.2, 97.0, 5.2, edge="none",
                  face=S.PALE_GREY, lw=0.0, radius=0.75)
    ax.text(3.5, 4.8, "LEAKAGE GUARD", fontsize=5.0,
            fontweight="bold", color=S.SLATE, ha="left", va="center")
    ax.plot([20.5, 20.5], [3.2, 6.4], color=S.LIGHT_GREY, lw=0.6)
    ax.text(22.5, 4.8, "Evaluation data excluded from construction; used only for held-out scoring.",
            fontsize=5.0, color=S.GREY, ha="left", va="center")

    S.save(fig, "fig5d_predictor_family")


if __name__ == "__main__":
    main()
