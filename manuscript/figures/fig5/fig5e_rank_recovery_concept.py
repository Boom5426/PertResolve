"""Fig. 5e — fixed predictors evaluated under two measurement regimes."""
from __future__ import annotations

W_MM, H_MM = 110.0, 34.0

import fig5_common as S

PREDICTORS = ["P5", "P4", "P3", "P2", "P1"]
LOW_RES = "#C8872F"
HIGH_RES = "#3F8B7D"


def _predictor_color(label):
    index = int(label[1:]) - 1
    return S.mix(S.LIGHT_SLATE, S.DARK_SLATE, index / 4.0)


def _tile(ax, x, y, label, width=10.0):
    color = _predictor_color(label)
    S.rounded_box(ax, x, y - 1.45, width, 2.9, edge=color,
                  face=S.mix("white", color, 0.11), lw=0.62, radius=0.5,
                  zorder=4)
    ax.text(x + width / 2.0, y, label, fontsize=5.45, fontweight="bold",
            color=S.DARK_SLATE, ha="center", va="center", zorder=5)


def _score_row(ax, y, label, center, half_width):
    color = _predictor_color(label)
    ax.text(32.2, y, label, fontsize=5.0, fontweight="bold", color=color,
            ha="right", va="center")
    ax.plot([center - half_width, center + half_width], [y, y], color=color,
            lw=1.05, solid_capstyle="round", zorder=3)
    ax.scatter([center], [y], s=9.5, color=color, edgecolor="white",
               linewidths=0.3, zorder=4)


def _benchmark_card(ax, y0, high_resolution):
    edge = HIGH_RES if high_resolution else LOW_RES
    title = "High resolution · replicate ceiling above chance" if high_resolution \
        else "Low resolution · replicate ceiling near chance"
    verdict = "RECOVERED" if high_resolution else "SCRAMBLED"
    observed = "P5 > P4 > P3\n> P2 > P1" if high_resolution \
        else "P2 > P5 > P1\n> P4 > P3"
    S.rounded_box(ax, 28.3, y0, 70.2, 17.4, edge=edge,
                  face=S.mix("white", edge, 0.035), lw=0.72, radius=1.0)
    ax.text(30.0, y0 + 15.6, title, fontsize=5.65, fontweight="bold",
            color=edge, ha="left", va="center")
    ax.plot([35.5, 70.3], [y0 + 1.35, y0 + 1.35], color=S.LIGHT_GREY,
            lw=0.5)
    ax.plot([35.5, 35.5], [y0 + 1.0, y0 + 1.7], color=S.LIGHT_GREY,
            lw=0.5)
    ax.plot([70.3, 70.3], [y0 + 1.0, y0 + 1.7], color=S.LIGHT_GREY,
            lw=0.5)
    ax.text(35.5, y0 + 0.25, "chance", fontsize=5.0, color=S.GREY,
            ha="center", va="bottom")
    ax.text(70.3, y0 + 0.25, "higher", fontsize=5.0, color=S.GREY,
            ha="center", va="bottom")

    rows = [y0 + 13.0, y0 + 10.4, y0 + 7.8, y0 + 5.2, y0 + 2.6]
    if high_resolution:
        centers = {"P5": 66.0, "P4": 59.8, "P3": 53.6, "P2": 47.4, "P1": 41.2}
        widths = {label: 2.1 for label in PREDICTORS}
    else:
        centers = {"P5": 54.8, "P4": 49.3, "P3": 47.0, "P2": 57.6, "P1": 52.0}
        widths = {label: 7.0 for label in PREDICTORS}
    for y, label in zip(rows, PREDICTORS):
        _score_row(ax, y, label, centers[label], widths[label])

    ax.plot([74.0, 74.0], [y0 + 2.0, y0 + 10.2], color=S.LIGHT_GREY, lw=0.55)
    ax.text(85.8, y0 + 12.4, "observed rank", fontsize=5.0, color=S.GREY,
            ha="center", va="center")
    ax.text(85.8, y0 + 7.4, observed, fontsize=5.2, fontweight="bold",
            color=S.DARK_SLATE, ha="center", va="center", linespacing=1.35)
    ax.text(85.8, y0 + 1.9, verdict, fontsize=5.1, fontweight="bold",
            color=edge, ha="center", va="center")


def main() -> None:
    S.apply_style()
    fig = S.figure(W_MM, H_MM)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 50)
    ax.axis("off")


    # One shared truth ladder feeds both benchmark regimes.
    S.rounded_box(ax, 1.5, 9.4, 18.5, 31.2, edge=S.LIGHT_GREY,
                  face="white", lw=0.7, radius=1.05)
    ax.text(10.75, 38.1, "Fixed quality order", fontsize=5.9,
            fontweight="bold", ha="center", va="center")
    ax.text(10.75, 35.6, "design-defined", fontsize=5.0, color=S.GREY,
            ha="center", va="center")
    ys = [32.3, 27.9, 23.5, 19.1, 14.7]
    for y, label in zip(ys, PREDICTORS):
        _tile(ax, 5.75, y, label, width=10.0)

    # Branch: the predictor set does not change.
    ax.plot([20.8, 24.0], [25.0, 25.0], color=S.MID_GREY, lw=0.65)
    ax.plot([24.0, 24.0], [16.0, 33.0], color=S.MID_GREY, lw=0.65)
    S.arrow(ax, (24.0, 33.0), (27.3, 33.0), color=S.MID_GREY,
            mutation=4.5)
    S.arrow(ax, (24.0, 16.0), (27.3, 16.0), color=S.MID_GREY,
            mutation=4.5)

    _benchmark_card(ax, 27.6, high_resolution=False)
    _benchmark_card(ax, 6.9, high_resolution=True)


    S.save(fig, "fig5e_rank_recovery_concept")


if __name__ == "__main__":
    main()
