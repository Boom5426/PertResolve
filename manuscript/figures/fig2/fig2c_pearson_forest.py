"""Fig 2c - Pearson-delta forest (direction recovery), right half of the b/c pair.

Data-direct. Source: results/results_v4_exttheta.csv via remote_data.exttheta().
Per method we take the mean pearson_delta over the scored variants (external-theta,
de-leaked grid) for the 18 in-house feature-model heads plus the Gene-mean reference;
WT-null and the external SOTA models (D.EXTERNAL) are excluded (they belong to Fig 4).
A 95% bootstrap CI over variants (2000 resamples, seed 0) is drawn per method.

Message: every predictor recovers perturbation *direction* far above zero
(mean pearson_delta ~0.55-0.65, every CI entirely > 0), which sets up the contrast
with Fig 2b where the same predictors cannot *rank* variant identity above chance.

Layout contract (Nature Methods pass):
  * Rows reuse ``fig2b_pds_forest.row_order()``, so b and c share ONE method-label
    column (drawn by b) and one feature-space key (drawn by b). This panel draws
    neither, which removes ~20 duplicated method names and a duplicated key.
  * Interval geometry/styling is imported from fig2b, so an interval in c is drawn
    exactly like an interval in b (same line width, caps and marker size).
  * WT-null occupies a row in the shared label column but has no point here: it
    predicts no change, so it has no direction to recover. Stated on the panel.

Run:  python fig2c_pearson_forest.py  ->  fig2c_pearson_forest.pdf (+ .png)
"""
from __future__ import annotations

import os
import sys

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import nm_style as S
import remote_data as D

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fig2b_pds_forest import (AX_BOTTOM_MM, AX_HEIGHT_MM, color_of, interval,
                              row_order)

N_BOOT = 2000
SEED = 0

# canvas: same height and same axes band as Fig 2b; no label column needed
W_MM, H_MM = 44.0, 63.0
AX_LEFT_MM, AX_WIDTH_MM = 4.0, 38.5


def compute() -> dict[str, tuple[float, float, float, int]]:
    """method -> (mean pearson_delta, ci_lo, ci_hi, n_variants)."""
    df = D.exttheta()
    methods = [m for m in df.method.unique() if D.is_inhouse(m)] + ["Gene-mean"]
    sub = df[df.method.isin(methods)]
    rng = np.random.default_rng(SEED)
    out = {}
    for m in methods:
        x = sub.loc[sub.method == m, "pearson_delta"].to_numpy()
        x = x[np.isfinite(x)]
        boot = np.empty(N_BOOT)
        for b in range(N_BOOT):
            boot[b] = np.mean(rng.choice(x, size=x.size, replace=True))
        lo, hi = np.percentile(boot, [2.5, 97.5])
        out[m] = (float(x.mean()), float(lo), float(hi), int(x.size))
    return out


def canvas():
    fig = plt.figure(figsize=(W_MM * S.MM, H_MM * S.MM))
    bg = fig.add_axes([0, 0, 1, 1])
    bg.set_axis_off()
    bg.set_xlim(0, 1)
    bg.set_ylim(0, 1)
    bg.add_patch(Rectangle((0, 0), 1, 1, facecolor="white", edgecolor="none",
                           zorder=-10))
    ax = fig.add_axes([AX_LEFT_MM / W_MM, AX_BOTTOM_MM / H_MM,
                       AX_WIDTH_MM / W_MM, AX_HEIGHT_MM / H_MM])
    return fig, ax


def main() -> None:
    S.apply_rcparams()
    # Keep Greek/maths glyphs in the same Arial-metric sans as the body text.
    # nm_style sets the text font but not mathtext, whose default (DejaVu Sans)
    # would embed a second typeface for every $\theta$, $\delta$ and subscript.
    plt.rcParams.update({"mathtext.fontset": "custom",
                         "mathtext.rm": "Liberation Sans",
                         "mathtext.it": "Liberation Sans:italic",
                         "mathtext.bf": "Liberation Sans:bold",
                         "mathtext.default": "it"})
    stats = compute()
    order = row_order()           # shared with Fig 2b
    n = len(order)

    assert n == 20, f"expected the 20-row shared order, got {n}"
    assert set(stats) == set(order) - {"WT-null"}, "row/method mismatch vs Fig 2b"
    assert all(lo > 0 for _, lo, _, _ in stats.values()), "a 95% CI reaches 0"

    fig, ax = canvas()

    for yi, m in enumerate(order):
        if m not in stats:        # WT-null: no predicted change, no direction
            continue
        mean, lo, hi, _ = stats[m]
        interval(ax, lo, hi, mean, yi, color_of(m))

    ax.set_yticks(range(n))
    ax.set_yticklabels([])        # the label column lives in Fig 2b
    ax.set_ylim(-0.8, n + 0.6)

    # Axis tightened to the plotted interval range (0.538-0.661) plus a small
    # margin. The previous 0-0.72 span left ~80% of the panel empty just to keep
    # the x = 0 reference on scale; 0 is now carried by the annotation below,
    # which is the claim the reader actually needs.
    ax.set_xlim(0.52, 0.68)
    ax.set_xticks([0.55, 0.60, 0.65])
    ax.set_xlabel(r"Pearson-$\delta$ (direction recovery)", labelpad=1.5)
    ax.tick_params(axis="y", length=0)
    ax.tick_params(axis="x", pad=1.5)

    S.despine(ax, keep=("left", "bottom"))
    ax.spines["left"].set_bounds(-0.8, n - 1 + 0.3)

    # panel notes, in the row headroom and the left margin (clear of every whisker)
    ax.text(0.523, n - 0.6, "rows as in b", ha="left", va="bottom",
            fontsize=5.5, color=S.GREY)
    ax.text(0.523, 1.6, "every 95% CI lies\nabove 0 (off scale)", ha="left",
            va="center", fontsize=5.5, color=S.INK, style="italic",
            linespacing=1.25)
    # the one row with no marker, named so the gap does not read as an error
    wt = order.index("WT-null")
    ax.text(0.523, wt, "WT-null: no predicted change", ha="left", va="center",
            fontsize=5.2, color=S.GREY, style="italic")

    S.save(fig, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             "fig2c_pearson_forest"))

    for m in reversed(order):
        if m in stats:
            mean, lo, hi, k = stats[m]
            print(f"{m:18} mean={mean:.3f}  CI[{lo:.3f},{hi:.3f}]  n={k}")


if __name__ == "__main__":
    main()
