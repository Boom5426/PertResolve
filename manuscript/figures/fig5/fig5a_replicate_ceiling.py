"""Fig. 5a — an independent replicate calibrates benchmark resolution."""
from __future__ import annotations

W_MM, H_MM = 100.0, 42.0

import numpy as np

import fig5_common as S

PROFILE_1 = [1.4, -0.8, 1.9, 0.5, -1.5, 1.0, -0.4, 1.6, -1.0, 0.7]
PROFILE_2 = [1.2, -0.6, 1.7, 0.7, -1.3, 1.1, -0.5, 1.5, -0.8, 0.6]
LOW_RES = "#C8872F"
HIGH_RES = "#3F8B7D"


def _cell_cloud(ax, x0, x1, y0, y1, color, phase):
    u = np.array([0.08, 0.23, 0.39, 0.54, 0.69, 0.84, 0.16, 0.47, 0.77])
    v = np.array([0.24, 0.73, 0.36, 0.82, 0.29, 0.67, 0.48, 0.12, 0.49])
    v = np.mod(v + phase, 1.0)
    ax.scatter(
        x0 + u * (x1 - x0),
        y0 + v * (y1 - y0),
        s=6.2,
        color=color,
        edgecolor="white",
        linewidths=0.25,
        zorder=4,
    )


def _step_label(ax, x, number, text):
    ax.text(x, 39.7, number, fontsize=7.1, fontweight="bold", color=S.SLATE,
            ha="left", va="center")
    ax.text(x + 3.3, 39.7, text, fontsize=6.1, fontweight="bold",
            ha="left", va="center")


def main() -> None:
    S.apply_style()
    fig = S.figure(W_MM, H_MM)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 100)
    # Cards span y 11.0 to 42.3; the old 0-52 range came from the panel's
    # original 52 mm height and left an eighth of the canvas empty at the
    # bottom while flattening the cards.
    ax.set_ylim(9.0, 44.5)
    ax.axis("off")


    cards = [
        (1.5, 11.0, 24.3, 31.3),
        (29.2, 11.0, 30.0, 31.3),
        (62.6, 11.0, 35.9, 31.3),
    ]
    for x, y, w, h in cards:
        S.rounded_box(ax, x, y, w, h, edge=S.LIGHT_GREY, face="white",
                      lw=0.7, radius=1.25)
    S.arrow(ax, (26.6, 26.6), (28.2, 26.6), color=S.MID_GREY)
    S.arrow(ax, (60.0, 26.6), (61.6, 26.6), color=S.MID_GREY)

    # 1 — two independent cell samples, linked only by variant identity.
    _step_label(ax, 3.4, "1", "Independent halves")
    S.rounded_box(ax, 4.0, 15.3, 19.3, 19.3, edge=S.LIGHT_GREY,
                  face=S.PALE_GREY, radius=0.9)
    ax.plot([4.8, 22.5], [24.9, 24.9], color=S.MID_GREY, lw=0.55,
            ls=(0, (2.2, 1.7)))
    ax.text(5.3, 32.1, "build cells", fontsize=5.3, color=S.SLATE,
            va="center")
    ax.text(5.3, 22.1, "evaluation cells", fontsize=5.3,
            color=S.DARK_SLATE, va="center")
    _cell_cloud(ax, 13.4, 22.2, 27.1, 33.4, S.MID_SLATE, 0.0)
    _cell_cloud(ax, 13.4, 22.2, 16.4, 23.1, S.DARK_SLATE, 0.18)
    ax.text(13.7, 13.4, "same variant identity\nno shared cells",
            fontsize=5.0, color=S.GREY, ha="center", va="center",
            linespacing=1.2)

    # 2 — both measurements enter the same evaluation harness in different roles.
    _step_label(ax, 31.1, "2", "Score identically")
    ax.text(31.8, 34.2, r"oracle prediction  $\delta_v^{build}$",
            fontsize=5.45, color=S.SLATE, va="center")
    S.mini_profile(ax, 32.0, 45.6, 30.5, PROFILE_1, color=S.MID_SLATE,
                   amplitude=2.6, lw=0.75)
    ax.text(31.8, 24.8, r"held-out target  $\delta_v^{eval}$",
            fontsize=5.45, color=S.DARK_SLATE, va="center")
    S.mini_profile(ax, 32.0, 45.6, 21.0, PROFILE_2, color=S.DARK_SLATE,
                   amplitude=2.6, lw=0.75)
    S.arrow(ax, (46.5, 29.9), (49.1, 27.8), color=S.MID_GREY, mutation=4.2)
    S.arrow(ax, (46.5, 21.5), (49.1, 24.7), color=S.MID_GREY, mutation=4.2)
    S.rounded_box(ax, 49.4, 21.5, 7.5, 7.2, edge=S.SLATE, face="#EEF1F4",
                  lw=0.65, radius=0.75)
    ax.text(53.15, 25.9, "PDS", fontsize=6.0, fontweight="bold",
            ha="center", va="center")
    ax.text(53.15, 23.4, "score", fontsize=5.0, color=S.GREY,
            ha="center", va="center")
    ax.text(44.2, 13.2, "same PDS harness as every model",
            fontsize=5.0, color=S.GREY, ha="center")

    # 3 — make the interpretive verdict visually explicit.
    _step_label(ax, 64.5, "3", "Benchmark verdict")
    # Verdict name on its own line above the condition, not in a right-hand
    # badge: the box is 30.7 units wide and the condition alone needs 24 of
    # them, so side by side the two always collided.
    S.rounded_box(ax, 65.2, 25.6, 30.7, 10.4, edge=LOW_RES,
                  face=S.mix("white", LOW_RES, 0.10), lw=0.75, radius=0.9)
    ax.text(67.0, 33.3, "MEASUREMENT-LIMITED", fontsize=5.4,
            fontweight="bold", color=LOW_RES, ha="left", va="center")
    ax.text(67.0, 30.2, "replicate PDS \u2248 chance", fontsize=5.15,
            color=S.DARK_SLATE, ha="left", va="center")
    ax.text(67.0, 27.4, "alleles unresolved", fontsize=5.15,
            color=S.DARK_SLATE, ha="left", va="center")

    S.rounded_box(ax, 65.2, 13.2, 30.7, 10.4, edge=HIGH_RES,
                  face=S.mix("white", HIGH_RES, 0.10), lw=0.75, radius=0.9)
    ax.text(67.0, 20.9, "MEASURABLE SIGNAL", fontsize=5.4,
            fontweight="bold", color=HIGH_RES, ha="left", va="center")
    ax.text(67.0, 17.8, "replicate PDS > chance", fontsize=5.15,
            color=S.DARK_SLATE, ha="left", va="center")
    ax.text(67.0, 15.0, "alleles measurable", fontsize=5.15,
            color=S.DARK_SLATE, ha="left", va="center")

    S.save(fig, "fig5a_replicate_ceiling")


if __name__ == "__main__":
    main()
