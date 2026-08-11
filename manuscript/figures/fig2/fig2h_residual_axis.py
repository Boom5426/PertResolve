"""Fig. 2h — the same methods on a second axis, and the reference that axis actually has.

Scoring after the gene-shared programme is removed asks the allele-level question directly.
The axis does not have chance at 0.50: an uninformative prediction becomes one constant
shared by every held-out variant and scores 0.524, matching its own permutation null exactly.
Both reference lines are drawn because reading the residual axis against 0.50 is the
misreading the analysis exists to prevent.

Traces to results/canonical/residual_axis_models.csv.
"""
from __future__ import annotations

import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import nm_style as S

HERE = os.path.dirname(os.path.abspath(__file__))

W_MM, H_MM = 170.0, 20.0

PDS_CHANCE = 0.500
RESIDUAL_NULL = 0.524


def main() -> None:
    S.apply_rcparams()
    models = pd.read_csv(S.data_path("residual_axis_models.csv")).sort_values("residual_pds")
    x = np.arange(len(models))

    fig = plt.figure(figsize=(W_MM / 25.4, H_MM / 25.4))
    ax = fig.add_axes([0.055, 0.30, 0.80, 0.60])

    ax.axhline(PDS_CHANCE, color=S.GREY, lw=0.6, ls=(0, (3, 2)), zorder=1)
    ax.axhline(RESIDUAL_NULL, color=S.HOTSPOT, lw=0.7, ls=(0, (2.2, 1.8)), zorder=1)

    ax.scatter(x, models["pds"], s=7, marker="o", color=S.HELDOUT, alpha=0.85,
               edgecolor="none", zorder=3)

    # Filled where the method exceeds its own permutation null, open where it does not.
    beats = models["perm_p_residual"] < 0.05
    ax.scatter(x[~beats.to_numpy()], models["residual_pds"][~beats], s=11, marker="D",
               facecolor="white", edgecolor=S.INK, linewidths=0.7, zorder=4)
    ax.scatter(x[beats.to_numpy()], models["residual_pds"][beats], s=11, marker="D",
               facecolor=S.INK, edgecolor=S.INK, linewidths=0.7, zorder=5)

    ax.set_xlim(-1.0, len(models))
    ax.set_ylim(0.435, 0.585)
    ax.set_xticks([])
    ax.set_yticks([0.45, 0.50, 0.55])
    ax.set_ylabel("score", labelpad=2.0)
    ax.set_xlabel("25 methods, ordered by residual-PDS", labelpad=1.0)
    S.despine(ax, keep=("left",))

    ax.text(len(models) - 0.6, PDS_CHANCE, " chance, PDS", fontsize=5.2, color=S.GREY,
            ha="left", va="center")
    ax.text(len(models) - 0.6, RESIDUAL_NULL, " null, residual-PDS", fontsize=5.2,
            color=S.HOTSPOT, ha="left", va="center", fontweight="bold")

    key = fig.add_axes([0.055, 0.02, 0.45, 0.13])
    key.set_xlim(0, 1)
    key.set_ylim(0, 1)
    key.axis("off")
    key.scatter([0.005], [0.5], s=7, marker="o", color=S.HELDOUT, clip_on=False)
    key.text(0.028, 0.5, "PDS", fontsize=5.3, color=S.HELDOUT, va="center")
    key.scatter([0.115], [0.5], s=11, marker="D", facecolor="white",
                edgecolor=S.INK, linewidths=0.7, clip_on=False)
    key.text(0.138, 0.5, "residual-PDS", fontsize=5.3, color=S.INK, va="center")
    key.scatter([0.325], [0.5], s=11, marker="D", facecolor=S.INK,
                edgecolor=S.INK, linewidths=0.7, clip_on=False)
    key.text(0.348, 0.5, "above its own permutation null", fontsize=5.3,
             color=S.INK, va="center")

    # exact=True, as every other Fig. 2 panel: fig2_assemble.tex places this at the 170 mm
    # it is drawn at, so a tight bounding box would rescale the whole panel and with it the
    # 5.2 pt type. Saved tight it came out 158.7 mm and was stretched by 7.1% on the page.
    S.save(fig, os.path.join(HERE, "fig2h_residual_axis"), exact=True,
           formats=("pdf", "png", "svg", "tiff"))


if __name__ == "__main__":
    main()
