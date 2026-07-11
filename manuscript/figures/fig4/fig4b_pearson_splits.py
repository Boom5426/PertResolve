"""Figure 4b - Direction recovery (Pearson-delta) across generalization splits.

Data-direct mirror of Fig 4a (which shows PDS across the same splits). Source
(committed): results/results_v4_exttheta.csv via remote_data.exttheta(). Each cell
is the mean per-variant pearson_delta for one in-house feature-model head
(18 heads = {Ridge,Lasso,RF,GBoost,KNN,MLP} x {theta,ESM,ESM+theta}) within one
generalization split (D.SPLIT_ORDER / D.SPLIT_LABELS).

Row order (shared with 4a): feature-space blocks theta -> ESM -> ESM+theta, and
within each block the fixed head order D.HEADS.

Message: direction recovery (Pearson-delta) stays positive across all five splits
(cell range ~0.41-0.71), and declines mainly under Mechanistic extrapolation
(split3, column mean 0.47) and Low-depth evaluation (split5, 0.55); it is highest
under Compatibility (split6, 0.69) and Random (split1, 0.65). This contrasts with
4a, where PDS sits at chance across all splits.

Run:  python fig4b_pearson_splits.py  ->  fig4b_pearson_splits.pdf (+ .png)
"""
from __future__ import annotations

import os
import sys

import numpy as np
from matplotlib import colormaps
from matplotlib.colors import Normalize
from matplotlib.patches import Rectangle

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import nm_style as S
import remote_data as D

VMIN, VMAX = 0.0, 0.8
CMAP = "Blues"


def ordered_methods() -> list[str]:
    """Row order shared with 4a: feature-space blocks x fixed head order."""
    rows = []
    for fs in ("theta", "esm", "esm+theta"):
        for head in D.HEADS:
            rows.append(f"{head}-{fs}")
    return rows


def head_of(method: str) -> str:
    return method.split("-", 1)[0]


def main() -> None:
    S.apply_rcparams()

    df = D.exttheta()
    ih = df[df.method.apply(D.is_inhouse)]
    piv = ih.groupby(["method", "split"]).pearson_delta.mean().unstack()

    rows = ordered_methods()
    assert all(m in piv.index for m in rows), "missing in-house method row"
    mat = piv.loc[rows, D.SPLIT_ORDER].to_numpy(float)
    assert not np.isnan(mat).any(), "NaN in pearson_delta matrix"

    n_rows, n_cols = mat.shape
    fig, ax = S.panel(60, 60)

    norm = Normalize(vmin=VMIN, vmax=VMAX)
    cmap = colormaps[CMAP]

    # draw each cell as a vector Rectangle (no rasterized imshow) with a thin
    # white border, so the whole panel stays fully vector in the PDF.
    for r in range(n_rows):
        for c in range(n_cols):
            ax.add_patch(Rectangle(
                (c - 0.5, r - 0.5), 1.0, 1.0,
                facecolor=cmap(norm(mat[r, c])),
                edgecolor="white", linewidth=0.5, zorder=2))

    # column ticks: split labels
    ax.set_xticks(range(n_cols))
    ax.set_xticklabels([D.SPLIT_LABELS[s] for s in D.SPLIT_ORDER],
                       rotation=40, ha="right", rotation_mode="anchor")
    ax.tick_params(axis="x", length=0, pad=1.5)

    # row ticks: head label per row (feature space shown as block brackets)
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels([head_of(m) for m in rows], fontsize=5)
    ax.tick_params(axis="y", length=0, pad=1.5)

    # feature-space block labels + separators (data coords, so they track the rows)
    block_size = len(D.HEADS)
    fs_labels = {"theta": "theta", "esm": "ESM", "esm+theta": "ESM+theta"}
    from matplotlib.transforms import blended_transform_factory
    btrans = blended_transform_factory(ax.transAxes, ax.transData)  # x=axes, y=data
    for i, fs in enumerate(("theta", "esm", "esm+theta")):
        y0 = i * block_size
        y_mid = y0 + (block_size - 1) / 2.0
        # bracket left of the row labels, vertically centred on the block's rows
        ax.text(-0.235, y_mid, fs_labels[fs], transform=btrans,
                rotation=90, ha="center", va="center", fontsize=6, color=S.INK)
        if i > 0:  # separator above this block
            ax.axhline(y0 - 0.5, color=S.INK, lw=0.8, zorder=4)

    ax.set_xlim(-0.5, n_cols - 0.5)
    ax.set_ylim(n_rows - 0.5, -0.5)  # row 0 (theta block) at top, reads top-down
    for side in ("top", "right", "left", "bottom"):
        ax.spines[side].set_visible(False)

    # vector colourbar: a stack of filled rectangles (matplotlib's own colorbar
    # rasterizes its gradient, so build it by hand to keep the PDF fully vector).
    cax = ax.inset_axes([1.05, 0.0, 0.055, 1.0])
    n_seg = 64
    seg_edges = np.linspace(VMIN, VMAX, n_seg + 1)
    for k in range(n_seg):
        v = 0.5 * (seg_edges[k] + seg_edges[k + 1])
        cax.add_patch(Rectangle((0.0, seg_edges[k]), 1.0, seg_edges[k + 1] - seg_edges[k],
                                facecolor=cmap(norm(v)), edgecolor="none", zorder=1))
    cax.set_xlim(0, 1)
    cax.set_ylim(VMIN, VMAX)
    cax.set_xticks([])
    cax.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8])
    cax.tick_params(axis="y", labelsize=5.5, length=1.8, width=0.5, pad=1.5)
    cax.yaxis.tick_right()
    for side in ("top", "bottom", "left", "right"):
        cax.spines[side].set_visible(True)
        cax.spines[side].set_linewidth(0.5)
        cax.spines[side].set_edgecolor(S.INK)
    cax.set_ylabel("Pearson-$\\delta$ (mean per variant)", fontsize=6, labelpad=3)
    cax.yaxis.set_label_position("right")

    S.save(fig, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             "fig4b_pearson_splits"))

    col_means = np.nanmean(mat, axis=0)
    print("cell range: [%.3f, %.3f]" % (mat.min(), mat.max()))
    for s, m in zip(D.SPLIT_ORDER, col_means):
        print("  %-14s mean=%.3f" % (D.SPLIT_LABELS[s], m))


if __name__ == "__main__":
    main()
