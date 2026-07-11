"""Fig 4a: PDS_cos across the 5 generalization splits (in-house methods).

Message: alignment PDS sits near/below the 0.50 chance line across ALL splits,
and collapses well below chance under the Low-depth split. Rows are the 18
feature-model heads (D.is_inhouse), columns the 5 held-in splits. Cell value is
the mean PDS_cos over all variants x genes for that (method, split).

Numbers plotted trace to results/results_v4_exttheta.csv (committed grid),
aggregated here with a plain mean of PDS_cos; nothing is fabricated.
"""
from __future__ import annotations

import os
import sys

import numpy as np
from matplotlib.colors import TwoSlopeNorm

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import nm_style as S
import remote_data as D


def main() -> None:
    S.apply_rcparams()

    df = D.exttheta()
    inhouse = sorted(m for m in df.method.unique() if D.is_inhouse(m))

    piv = (
        df[df.method.isin(inhouse)]
        .groupby(["method", "split"])["PDS_cos"]
        .mean()
        .unstack()[D.SPLIT_ORDER]
    )
    # order rows by overall mean PDS (best at top) for a readable gradient
    piv = piv.loc[piv.mean(axis=1).sort_values(ascending=False).index]

    methods = list(piv.index)
    splits = list(piv.columns)
    M = piv.to_numpy()

    print("mean PDS_cos per split:", dict(zip(splits, piv.mean(axis=0).round(3))))
    print("global min/max:", round(M.min(), 3), round(M.max(), 3))

    # diverging map centred on chance (0.50); symmetric range keeps the centre honest
    vmin, vmax, vcen = 0.35, 0.65, 0.50
    norm = TwoSlopeNorm(vmin=vmin, vcenter=vcen, vmax=vmax)
    cmap = plt_cmap()

    fig, ax = S.panel(60, 60)

    nrow, ncol = M.shape
    # pcolormesh renders as vector (no rasterized flag). Rows top-to-bottom:
    # place row 0 at the top by inverting the y axis after drawing.
    x_edges = np.arange(ncol + 1)
    y_edges = np.arange(nrow + 1)
    ax.pcolormesh(x_edges, y_edges, M, cmap=cmap, norm=norm,
                  edgecolors="white", linewidth=0.6)

    # annotate each cell (centre) with 2-dp value; text colour for contrast
    for i in range(nrow):
        for j in range(ncol):
            v = M[i, j]
            rgba = cmap(norm(v))
            lum = 0.299 * rgba[0] + 0.587 * rgba[1] + 0.114 * rgba[2]
            tc = "white" if lum < 0.5 else S.INK
            ax.text(j + 0.5, i + 0.5, f"{v:.2f}", ha="center", va="center",
                    fontsize=4.4, color=tc)

    ax.set_xlim(0, ncol)
    ax.set_ylim(0, nrow)
    ax.invert_yaxis()  # row 0 (best method) at the top
    ax.set_aspect("auto")

    # column ticks: split role labels (centred on cells)
    ax.set_xticks(np.arange(ncol) + 0.5)
    ax.set_xticklabels([D.SPLIT_LABELS[s] for s in splits], rotation=35,
                       ha="right", rotation_mode="anchor", fontsize=5.5)
    # row ticks: method names
    ax.set_yticks(np.arange(nrow) + 0.5)
    ax.set_yticklabels(methods, fontsize=4.6)

    ax.tick_params(length=0)
    for side in ("top", "right", "left", "bottom"):
        ax.spines[side].set_visible(False)

    # slim colourbar drawn as a vector pcolormesh (fig.colorbar rasterizes).
    cax = fig.add_axes([0.90, 0.30, 0.028, 0.42])
    grad = np.linspace(vmin, vmax, 256)
    cax.pcolormesh(np.array([0, 1]), grad,
                   grad[:-1, None], cmap=cmap, norm=norm,
                   shading="flat", rasterized=False)
    cax.set_xlim(0, 1)
    cax.set_ylim(vmin, vmax)
    cax.set_xticks([])
    cax.yaxis.set_label_position("right")
    cax.yaxis.tick_right()
    cax.set_yticks([0.35, 0.50, 0.65])
    cax.set_yticklabels(["0.35", "0.50", "0.65"], fontsize=5)
    cax.tick_params(length=1.6, width=0.4, pad=1.5)
    for side in ("top", "bottom", "left", "right"):
        cax.spines[side].set_visible(True)
        cax.spines[side].set_linewidth(0.4)
    # mark the chance centre on the bar
    cax.axhline(0.50, color=S.INK, linewidth=0.5)
    cax.set_ylabel("mean PDS$_{cos}$", fontsize=5.5, labelpad=3, rotation=90)

    ax.set_title("PDS near / below chance across all splits",
                 fontsize=6, pad=4)

    S.save(fig, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             "fig4a_pds_splits"))


def plt_cmap():
    """Colour-blind-safe diverging map: purple (low) -> light -> teal (high)."""
    from matplotlib.colors import LinearSegmentedColormap
    # PuOr-like but tuned; below-chance = purple, above-chance = green/teal
    return LinearSegmentedColormap.from_list(
        "pds_div",
        ["#5B2C83", "#9C79BE", "#EDE8F0", "#7FBF9B", "#2E7D57"],
    )


if __name__ == "__main__":
    main()
