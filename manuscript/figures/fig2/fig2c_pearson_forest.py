"""Figure 2c - Pearson-delta forest (direction fidelity) for the in-house predictors.

Data-direct. Source: results/results_v4_exttheta.csv via remote_data.exttheta().
Per method we take the mean pearson_delta over the scored variants (external-theta,
de-leaked grid) for the 18 in-house feature-model heads plus the Gene-mean reference;
WT-null and the external SOTA models (D.EXTERNAL) are excluded (they belong to Fig 4).
A 95% bootstrap CI over variants (2000 resamples, seed 0) is drawn per method.

Message: every predictor recovers perturbation *direction* far above zero
(mean pearson_delta ~0.55-0.65, every CI entirely > 0), which sets up the contrast
with Fig 2b where the same predictors cannot *rank* variant identity above chance.

Run:  python fig2c_pearson_forest.py  ->  fig2c_pearson_forest.pdf (+ .png)
"""
from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import nm_style as S
import remote_data as D

N_BOOT = 2000
SEED = 0


def compute():
    df = D.exttheta()
    methods = [m for m in df.method.unique() if D.is_inhouse(m)] + ["Gene-mean"]
    sub = df[df.method.isin(methods)]
    rng = np.random.default_rng(SEED)
    rows = []
    for m in methods:
        x = sub.loc[sub.method == m, "pearson_delta"].to_numpy()
        x = x[np.isfinite(x)]
        boot = np.empty(N_BOOT)
        for b in range(N_BOOT):
            boot[b] = np.mean(rng.choice(x, size=x.size, replace=True))
        lo, hi = np.percentile(boot, [2.5, 97.5])
        rows.append((m, float(x.mean()), float(lo), float(hi), int(x.size)))
    # forest convention: best (highest) at top
    rows.sort(key=lambda r: r[1])
    return rows


def color_for(method: str) -> str:
    if method in ("Gene-mean", "WT-null"):
        return S.FEATURE_COLORS["reference"]
    return S.FEATURE_COLORS.get(D.feature_space(method), S.FEATURE_COLORS["reference"])


def main() -> None:
    S.apply_rcparams()
    rows = compute()
    n = len(rows)

    fig, ax = S.panel(58, 62)
    ys = np.arange(n)

    ax.axvline(0.0, color=S.INK, lw=0.6, zorder=1)

    for y, (m, mean, lo, hi, k) in zip(ys, rows):
        c = color_for(m)
        ax.plot([lo, hi], [y, y], color=c, lw=0.9, solid_capstyle="round", zorder=3)
        ax.plot([mean], [y], marker="o", ms=2.6, mfc=c, mec="white",
                mew=0.4, zorder=4)

    ax.set_yticks(ys)
    ax.set_yticklabels([m for m, *_ in rows], fontsize=5)
    for t, (m, *_ ) in zip(ax.get_yticklabels(), rows):
        t.set_color(color_for(m))
    ax.set_ylim(-0.7, n - 0.3)

    ax.set_xlim(0.0, 0.70)
    ax.set_xticks([0.0, 0.2, 0.4, 0.6])
    ax.set_xlabel(r"Pearson $\delta$ (direction fidelity)")
    ax.tick_params(length=2.2)

    S.despine(ax, keep=("left", "bottom"))

    # annotation: every interval above zero
    ax.annotate("all intervals\nabove 0", xy=(0.0, n - 1.4),
                xytext=(0.10, n - 3.2), fontsize=5.4, color=S.INK,
                ha="left", va="center",
                arrowprops=dict(arrowstyle="-", lw=0.5, color=S.INK,
                                connectionstyle="arc3,rad=-0.25"))

    # feature-space legend (line + dot proxies, no black edges)
    handles = []
    for lab, key in [("theta", "theta"), ("ESM", "ESM"), ("ESM+theta", "ESM+theta"),
                     ("Gene-mean", "ref")]:
        cc = S.FEATURE_COLORS["reference"] if key == "ref" else S.FEATURE_COLORS[key]
        handles.append(ax.plot([], [], marker="o", ms=2.6, mfc=cc, mec="white",
                               mew=0.4, color=cc, lw=0.9, ls="-", label=lab)[0])
    leg = ax.legend(handles=handles, loc="center left",
                    bbox_to_anchor=(0.04, 0.42), fontsize=5,
                    handlelength=1.6, handletextpad=0.4, labelspacing=0.4,
                    borderpad=0.3, frameon=False)

    S.save(fig, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             "fig2c_pearson_forest"))

    for m, mean, lo, hi, k in reversed(rows):
        print(f"{m:18} mean={mean:.3f}  CI[{lo:.3f},{hi:.3f}]  n={k}")


if __name__ == "__main__":
    main()
