"""Figure 4d - PDS ranking failure is invariant to the distance geometry.

Data-direct. Source (committed): the D.exttheta() per-variant external-calibration
grid; for each in-house method we take the per-method mean of PDS_cos, PDS_L1 and
PDS_L2 (cosine, L1, L2 variants of the perturbation discrimination score). One dot per method at each of the three distance positions, a thin line linking a method's three
values, a dashed chance line at 0.50, and a median bar per distance.

Verified firsthand (18 in-house methods, per-method means then median across
methods):
  median PDS_cos = 0.455   median PDS_L1 = 0.453   median PDS_L2 = 0.455
All three cluster near the 0.50 chance line: switching the underlying distance
from cosine to L1 to L2 does not rescue variant-level ranking.

Run:  python fig4d_distance.py  ->  fig4d_distance.pdf (+ .png)
"""
from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import nm_style as S
import remote_data as D


METRICS = ["PDS_cos", "PDS_L1", "PDS_L2"]
XLABELS = ["cosine", "L1", "L2"]


def main() -> None:
    S.apply_rcparams()
    df = D.exttheta()
    inhouse = sorted(m for m in df["method"].unique() if D.is_inhouse(m))
    sub = df[df["method"].isin(inhouse)]

    # per-method mean of each distance variant -> one row per method
    per_method = sub.groupby("method")[METRICS].mean()
    n_methods = len(per_method)
    xpos = np.arange(len(METRICS))
    medians = per_method.median().to_numpy()

    fig, ax = S.panel(57.5, 48.6)

    # dashed chance line at 0.50
    ax.axhline(0.50, color=S.INK, ls=(0, (4, 3)), lw=0.7, zorder=1)

    # thin line linking each method's three distance values
    for m in inhouse:
        y = per_method.loc[m, METRICS].to_numpy(float)
        ax.plot(xpos, y, "-", color=S.LIGHT_GREY, lw=0.4, zorder=2)

    # one dot per method at each distance position (slight x-jitter for density)
    rng = np.random.default_rng(0)
    for j, met in enumerate(METRICS):
        y = per_method[met].to_numpy(float)
        jit = (rng.random(n_methods) - 0.5) * 0.16
        ax.scatter(np.full(n_methods, xpos[j]) + jit, y,
                   s=7.5, facecolor=S.GREY, edgecolor="white",
                   linewidth=0.3, zorder=3)

    # median bar per distance (short horizontal segment)
    bar_hw = 0.30
    for j in range(len(METRICS)):
        ax.plot([xpos[j] - bar_hw, xpos[j] + bar_hw], [medians[j], medians[j]],
                "-", color=S.INK, lw=1.4, solid_capstyle="butt", zorder=5)
        ax.text(xpos[j], 0.412, f"median\n{medians[j]:.3f}",
                fontsize=5.2, color=S.INK, ha="center", va="bottom", zorder=5,
                linespacing=1.2)

    ax.set_xlim(-0.5, len(METRICS) - 0.5)
    ax.set_ylim(0.40, 0.52)
    ax.set_xticks(xpos)
    ax.set_xticklabels(XLABELS)
    ax.set_yticks([0.40, 0.45, 0.50])
    ax.set_xlabel("Distance geometry")
    ax.set_ylabel("PDS (perturbation\ndiscrimination score)")

    # annotate the chance line
    ax.text(len(METRICS) - 0.55, 0.503, "chance (0.50)", fontsize=5.2,
            color=S.INK, ha="right", va="bottom", zorder=5)

    S.despine(ax)
    ax.tick_params(length=2.2)

    S.save(fig, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             "fig4d_distance"))
    line = "  ".join(f"{met}: median={medians[j]:.4f}"
                     for j, met in enumerate(METRICS))
    print(f"n_methods={n_methods}  {line}")


if __name__ == "__main__":
    main()
