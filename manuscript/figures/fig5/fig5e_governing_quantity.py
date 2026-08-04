"""Fig. 5e — the mean over pairs does not order the attainable ceiling; the nearest one does.

Replaces an illustrative panel that read no data. Every point here is measured. The dialled
points hold the mean separation fixed while changing the geometry of the configuration, which
is what separates the two candidate axes; four real datasets cannot, because their geometry
happens to co-vary with their mean separation.
"""
from __future__ import annotations

W_MM, H_MM = 110.0, 34.0

import numpy as np

import fig5_common as S

REAL_MARKERS = {"TP53": "o", "KRAS": "s", "GATA1": "^", "JAK1": "D"}
LINTHRESH = 0.01


def spearman(x, y) -> float:
    """Rank correlation, computed here so the annotation cannot drift from the table."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    return float(np.corrcoef(np.argsort(np.argsort(x[ok])),
                             np.argsort(np.argsort(y[ok])))[0, 1])


def _scatter(ax, sweep, column, label):
    dialled = sweep[sweep.source == "dialled"]
    real = sweep[sweep.source == "real"]

    ax.scatter(dialled[column], dialled["ceiling_pds"], s=1.6,
               color=S.LIGHT_SLATE, alpha=0.45, edgecolor="none", zorder=2)
    for gene, marker in REAL_MARKERS.items():
        rows = real[real.gene == gene]
        ax.scatter(rows[column], rows["ceiling_pds"], s=13, marker=marker,
                   facecolor="white", edgecolor=S.GENE_COLORS[gene],
                   linewidths=0.85, zorder=4)

    ax.axhline(0.5, color=S.MID_GREY, lw=0.55, ls=(0, (3, 2)), zorder=1)
    # A symmetric log axis is the only one that shows both the negative estimates near the
    # floor, where these datasets sit, and the three decades the dial covers.
    ax.set_xscale("symlog", linthresh=LINTHRESH, linscale=0.45)
    ax.set_xlim(-0.03, 8.0)
    ax.set_ylim(0.46, 1.03)
    ax.set_yticks([0.5, 0.7, 0.9])
    ax.set_xticks([-0.01, 0, 0.01, 0.1, 1.0])
    ax.set_xticklabels(["−0.01", "0", "0.01", "0.1", "1"])
    ax.set_xlabel(label, labelpad=1.5)
    S.despine(ax)

    rho = spearman(sweep[column], sweep["ceiling_pds"])
    ax.text(0.97, 0.06, f"Spearman {rho:.3f}", transform=ax.transAxes,
            fontsize=5.8, fontweight="bold", ha="right", va="bottom",
            color=S.INK if rho > 0.9 else S.GREY)


def main() -> None:
    S.apply_style()
    sweep = S.read_csv("resolution_sweep.csv")

    fig = S.figure(W_MM, H_MM)
    left = fig.add_axes([0.088, 0.30, 0.383, 0.615])
    right = fig.add_axes([0.578, 0.30, 0.383, 0.615])

    _scatter(left, sweep, "rho2", "mean separation / noise")
    _scatter(right, sweep, "rho2_nn_median", "nearest-competitor separation / noise")

    left.set_ylabel("replicate ceiling", labelpad=2.0)
    right.set_yticklabels([])

    n_dialled = int((sweep.source == "dialled").sum())
    fig.text(0.088, 0.99,
             f"{n_dialled:,} configurations with the geometry dialled at fixed mean "
             f"separation (grey), and the four datasets (open)",
             fontsize=5.3, color=S.GREY, ha="left", va="top")

    # One key for both axes, on its own transparent axis so the markers are drawn at the
    # same size as the data and cannot drift from them.
    key = fig.add_axes([0.088, 0.02, 0.55, 0.10])
    key.set_xlim(0, 1)
    key.set_ylim(0, 1)
    key.axis("off")
    for index, (gene, marker) in enumerate(REAL_MARKERS.items()):
        x = index * 0.25
        key.scatter([x], [0.5], s=13, marker=marker, facecolor="white",
                    edgecolor=S.GENE_COLORS[gene], linewidths=0.85, clip_on=False)
        key.text(x + 0.035, 0.5, gene, fontsize=5.4, color=S.GENE_COLORS[gene],
                 ha="left", va="center")

    S.save(fig, "fig5e_governing_quantity")


if __name__ == "__main__":
    main()
