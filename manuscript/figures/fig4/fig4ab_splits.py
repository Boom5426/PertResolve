"""Figure 4a,b - PDS and direction recovery across the five generalization splits.

Single drawing that carries BOTH panel a (PDS, perturbation discrimination score)
and panel b (Pearson-delta, direction recovery). The two heatmaps show the SAME 18
in-house feature-model heads x the SAME 5 splits, so the method-label column and
the feature-space brackets are drawn ONCE at the far left and serve both maps
(previously all 18 row labels were printed twice). The composite stamps the panel
letters a and b over the two blocks.

Data-direct. Source (committed): results/results_v4_exttheta.csv via
remote_data.exttheta(). Each cell is the plain mean over held-out variants of
PDS_cos (left) or pearson_delta (right) for one head
({Ridge,Lasso,RF,GBoost,KNN,MLP} x {theta,ESM,ESM+theta}) within one split
(D.SPLIT_ORDER / D.SPLIT_LABELS). No number is transformed beyond this mean.

Colour logic (one logic for the pair, greyscale-safe): ONE sequential lightness
ramp on ONE shared scale (0.20-0.80) with ONE shared colourbar, pale = low score,
dark = high score. Both quantities are means of bounded per-variant scores where
higher means better prediction, so a common scale is readable, and the PDS chance
level (0.50) is marked on the bar. The earlier pairing of a purple-green diverging
map with a blue sequential map encoded the same kind of quantity two different
ways and is not used.

Verified firsthand (column means over the 18 heads, printed by this script):
  PDS      Random 0.482  Positional 0.491  Mechanistic 0.503  Low-depth 0.244
           Compatibility 0.479     cell range 0.229-0.539
  Pearson  Random 0.653  Positional 0.626  Mechanistic 0.468  Low-depth 0.553
           Compatibility 0.688     cell range 0.406-0.712

Message: PDS sits in the chance regime in every split, with the deepest collapse
under Low-depth, while direction recovery stays clearly positive throughout.

Fully vector (Rectangle cells, hand-built colourbar); check with:
  pdfimages -list fig4ab_splits.pdf | tail -n +3

Run:  python fig4ab_splits.py  ->  fig4ab_splits.pdf (+ .png)
"""
from __future__ import annotations

import os
import sys

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Rectangle
from matplotlib.transforms import blended_transform_factory

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import nm_style as S
import remote_data as D

HERE = os.path.dirname(os.path.abspath(__file__))

FS_BLOCKS = ("theta", "esm", "esm+theta")
FS_LABELS = {"theta": r"$\theta$", "esm": "ESM", "esm+theta": r"ESM+$\theta$"}

# one sequential lightness ramp for both maps (house neutrals -> dark slate)
RAMP = ["#FFFFFF", "#DCE0E4", "#AEB6BE", "#77828C", "#3C464F"]

# shared colour scale for both metrics; 0.50 is the PDS chance level
V_LO, V_HI, V_CHANCE = 0.20, 0.80, 0.50
TICKS = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]

# panel geometry, millimetres (drawn at final placement size, composite scale ~1)
W_MM, H_MM = 89.0, 49.0
X_LAB, W_HM, GUTTER, GAP_CB, W_CB = 13.5, 27.0, 6.0, 2.4, 2.8
Y_BOT, H_HM = 11.0, 33.0


def ramp() -> LinearSegmentedColormap:
    return LinearSegmentedColormap.from_list("score_seq", RAMP)


def ordered_methods() -> list[str]:
    """Row order shared by both maps: feature-space blocks x fixed head order."""
    return [f"{head}-{fs}" for fs in FS_BLOCKS for head in D.HEADS]


def head_of(method: str) -> str:
    return method.split("-", 1)[0]


def _rect(x_mm: float, y_mm: float, w_mm: float, h_mm: float) -> list[float]:
    """Axes rectangle in figure fractions from millimetre coordinates."""
    return [x_mm / W_MM, y_mm / H_MM, w_mm / W_MM, h_mm / H_MM]


def draw_map(fig, x0_mm: float, mat: np.ndarray, cmap, norm,
             show_rows: bool, title: str) -> None:
    """Draw one heatmap block; row labels only on the leftmost block."""
    n_rows, n_cols = mat.shape
    ax = fig.add_axes(_rect(x0_mm, Y_BOT, W_HM, H_HM))

    for r in range(n_rows):
        for c in range(n_cols):
            ax.add_patch(Rectangle(
                (c - 0.5, r - 0.5), 1.0, 1.0,
                facecolor=cmap(norm(mat[r, c])),
                edgecolor="white", linewidth=0.5, zorder=2))

    ax.set_xticks(range(n_cols))
    ax.set_xticklabels([D.SPLIT_LABELS[s] for s in D.SPLIT_ORDER],
                       rotation=40, ha="right", rotation_mode="anchor",
                       fontsize=5.4)
    ax.tick_params(axis="x", length=0, pad=1.5)

    if show_rows:
        rows = ordered_methods()
        ax.set_yticks(range(n_rows))
        ax.set_yticklabels([head_of(m) for m in rows], fontsize=5.2)
        ax.tick_params(axis="y", length=0, pad=1.5)
        btrans = blended_transform_factory(ax.transAxes, ax.transData)
        block = len(D.HEADS)
        for i, fs in enumerate(FS_BLOCKS):
            ax.text(-0.30, i * block + (block - 1) / 2.0, FS_LABELS[fs],
                    transform=btrans, rotation=90, ha="center", va="center",
                    fontsize=6.6, color=S.INK)
    else:
        ax.set_yticks([])

    # feature-space separators on BOTH blocks so the shared row order stays
    # legible without repeating the labels
    for i in range(1, len(FS_BLOCKS)):
        ax.axhline(i * len(D.HEADS) - 0.5, color=S.INK, lw=0.8, zorder=4)

    ax.set_xlim(-0.5, n_cols - 0.5)
    ax.set_ylim(n_rows - 0.5, -0.5)  # theta block on top, reads top-down
    for side in ("top", "right", "left", "bottom"):
        ax.spines[side].set_visible(False)
    ax.set_title(title, fontsize=6.0, pad=2.5, color=S.INK)


def draw_cbar(fig, x0_mm: float, cmap, norm) -> None:
    """One hand-built vector colourbar shared by both maps.

    matplotlib's own colorbar rasterizes its gradient, so it is built from filled
    rectangles to keep the PDF fully vector.
    """
    cax = fig.add_axes(_rect(x0_mm, Y_BOT, W_CB, H_HM))
    n_seg = 120
    edges = np.linspace(V_LO, V_HI, n_seg + 1)
    for k in range(n_seg):
        v = 0.5 * (edges[k] + edges[k + 1])
        cax.add_patch(Rectangle((0.0, edges[k]), 1.0, edges[k + 1] - edges[k],
                                facecolor=cmap(norm(v)), edgecolor="none", zorder=1))
    cax.set_xlim(0, 1)
    cax.set_ylim(V_LO, V_HI)
    cax.set_xticks([])
    cax.set_yticks(TICKS)
    cax.set_yticklabels([f"{t:.1f} chance" if abs(t - V_CHANCE) < 1e-9
                         else f"{t:.1f}" for t in TICKS])
    cax.tick_params(axis="y", labelsize=5.2, length=1.8, width=0.5, pad=1.5)
    cax.yaxis.tick_right()
    for side in ("top", "bottom", "left", "right"):
        cax.spines[side].set_visible(True)
        cax.spines[side].set_linewidth(0.5)
        cax.spines[side].set_edgecolor(S.INK)
    cax.axhline(V_CHANCE, color=S.INK, lw=0.9, zorder=3)
    cax.annotate("Mean score", xy=(0.0, 1.0), xycoords="axes fraction",
                 xytext=(0, 3.2), textcoords="offset points",
                 ha="left", va="bottom", fontsize=5.4, color=S.INK,
                 annotation_clip=False)


def main() -> None:
    S.apply_rcparams()
    plt.rcParams["savefig.bbox"] = "standard"  # exact mm canvas, no tight crop

    df = D.exttheta()
    ih = df[df.method.apply(D.is_inhouse)]
    rows = ordered_methods()
    assert len(rows) == 18, f"expected 18 in-house heads, got {len(rows)}"

    pds = ih.groupby(["method", "split"]).PDS_cos.mean().unstack()
    pea = ih.groupby(["method", "split"]).pearson_delta.mean().unstack()
    for piv in (pds, pea):
        assert all(m in piv.index for m in rows), "missing in-house method row"
    mat_a = pds.loc[rows, D.SPLIT_ORDER].to_numpy(float)
    mat_b = pea.loc[rows, D.SPLIT_ORDER].to_numpy(float)
    assert not np.isnan(mat_a).any() and not np.isnan(mat_b).any(), "NaN in matrix"
    assert mat_a.min() >= V_LO and mat_a.max() <= V_HI, "PDS outside colour scale"
    assert mat_b.min() >= V_LO and mat_b.max() <= V_HI, "Pearson outside scale"

    cmap = ramp()
    norm = Normalize(vmin=V_LO, vmax=V_HI)

    fig = plt.figure(figsize=(W_MM * S.MM, H_MM * S.MM))

    x_a = X_LAB
    draw_map(fig, x_a, mat_a, cmap, norm, True, "PDS (cosine)")
    x_b = x_a + W_HM + GUTTER
    draw_map(fig, x_b, mat_b, cmap, norm, False, r"Pearson-$\delta$ (direction)")
    draw_cbar(fig, x_b + W_HM + GAP_CB, cmap, norm)

    S.save(fig, os.path.join(HERE, "fig4ab_splits"))

    for name, mat in (("PDS_cos", mat_a), ("pearson_delta", mat_b)):
        print(f"{name}: cell range [{mat.min():.3f}, {mat.max():.3f}]")
        for s, m in zip(D.SPLIT_ORDER, np.nanmean(mat, axis=0)):
            print(f"    {D.SPLIT_LABELS[s]:<14s} mean={m:.3f}")


if __name__ == "__main__":
    main()
