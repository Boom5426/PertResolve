"""Figure 4c - Empirical permutation null confirms PDS is calibrated at 0.50.

Data-direct. Source (committed): results/permutation_null_pds.csv (the debt2
within-split label-permutation control), metric PDS_cos, 340 method x split x gene
combinations. Verified firsthand:
  permutation null mean = 0.500 (range of per-combination null means 0.495-0.504)
  observed PDS distribution centred at 0.4995
  3.2% of combinations exceed their own permutation 95th percentile (expected 5%)
Message: PDS at chance is not an artefact of the analytic 0.50 baseline; the
empirical null is centred at 0.50 and observed exceedances stay below the 5% rate.

Run:  python fig4c_null.py  ->  fig4c_null.pdf (+ .png)
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import nm_style as S


def main() -> None:
    S.apply_rcparams()
    df = pd.read_csv(S.repo_root() / "results" / "permutation_null_pds.csv")
    pc = df[df.metric == "PDS_cos"]
    obs = pc.obs_pds.to_numpy(float)
    null_lo, null_hi = pc.null_mean.min(), pc.null_mean.max()
    frac_exceed = (pc.p_above_obs < 0.05).mean() * 100
    n = len(pc)
    n_half = int((np.abs(obs - 0.5) < 1e-6).sum())

    fig, ax = S.panel(58, 46)

    ax.axvspan(null_lo, null_hi, color=S.GREY, alpha=0.30, zorder=0)  # null-mean range
    ax.hist(obs, bins=25, color=S.GENE_COLORS["TP53"], alpha=0.8,
            edgecolor="white", linewidth=0.3, zorder=2)
    ax.axvline(0.50, color=S.INK, ls="--", lw=0.8, zorder=3)

    ax.set_xlabel("Observed PDS (cosine), per method x split x gene")
    ax.set_ylabel("Combinations")
    ax.set_xlim(0.33, 0.62)
    ax.set_ylim(0, max(180, n_half + 12))
    S.despine(ax)
    ax.tick_params(length=2.2)

    # turn the central spike into an explained feature (label to the right of the spike)
    ax.annotate(f"{100 * n_half / n:.0f}% exactly\nat 0.50\n(ties)",
                xy=(0.505, 150), xytext=(0.545, 120), fontsize=5,
                color=S.INK, ha="left", va="center",
                arrowprops=dict(arrowstyle="-", lw=0.5, color=S.GREY))
    ax.text(0.02, 0.985,
            f"permutation null mean 0.500\n(0.495-0.504 across combos)\n"
            f"{frac_exceed:.1f}% exceed 95th-pct null\n(expected 5%, n={n})",
            transform=ax.transAxes, fontsize=5.2, color=S.GREY, va="top", ha="left",
            linespacing=1.5)

    S.save(fig, "fig4c_null")
    print(f"obs mean={obs.mean():.4f}  {n_half}/{n} exactly at 0.5 "
          f"({100*n_half/n:.0f}%)  null[{null_lo:.4f},{null_hi:.4f}]  "
          f"exceed(p<0.05)={frac_exceed:.1f}%")


if __name__ == "__main__":
    main()
