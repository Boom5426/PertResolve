"""Figure 4e - Feature space swap saturates direction but leaves PDS at chance.

Data-direct. Source (committed): results/results_v4_exttheta.csv (the per-variant
breakdown grid, D.exttheta()). The 18 in-house feature-model heads
(Ridge/Lasso/RF/GBoost/KNN/MLP x {theta, ESM, ESM+theta}) are grouped by
D.feature_space into three input representations. Per feature space we plot the
mean pearson_delta (direction) and the mean PDS_cos (perturbation discrimination
score), with error bars = SD across the 6 model heads in that space.

Verified firsthand (mean over the 6 heads, i.e. mean of per-method means):
  feature space   pearson_delta         PDS_cos
  theta           0.6005 (sd 0.0278)    0.4656 (sd 0.0168)
  ESM             0.6168 (sd 0.0374)    0.4540 (sd 0.0104)
  ESM+theta       0.6081 (sd 0.0278)    0.4625 (sd 0.0221)

Message: swapping the input representation (biophysical theta -> ESM protein
language-model embeddings -> both) moves direction only marginally around a
saturated ~0.6, while PDS stays pinned near the 0.50 chance line. Richer features
buy no allele-discrimination skill.

Encoding notes (Nature Methods pass):
  - The two series were previously drawn on twin y axes whose ranges were in fact
    identical (0.30-0.75), so the second axis carried no information and made it
    ambiguous which series belonged to which scale. Both series now share ONE y
    axis, which is also what makes the direction/PDS gap directly readable.
  - Direction was drawn in a terracotta red that is not in the house palette (red
    is reserved for signed direction / warnings). The house metric encoding used
    across Fig. 4 is now: direction (Pearson-delta) = ink filled circle,
    PDS = grey open square. Shape carries the distinction in greyscale.
  - The detached legend is replaced by direct labels next to each series.

Run:  python fig4e_features.py  ->  fig4e_features.pdf (+ .png)
"""
from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import nm_style as S
import remote_data as D


FS_ORDER = ["theta", "ESM", "ESM+theta"]
FS_LABEL = {"theta": r"$\theta$", "ESM": "ESM", "ESM+theta": r"ESM+$\theta$"}

# house metric encoding, shared with Fig. 4f
C_DIR = S.INK        # direction (Pearson-delta): filled circle
C_PDS = S.GREY       # PDS (cosine): open square


def main() -> None:
    S.apply_rcparams()
    df = D.exttheta()
    inh = df[df.method.apply(D.is_inhouse)].copy()
    inh["fs"] = inh.method.apply(D.feature_space)

    dir_mean, dir_sd, pds_mean, pds_sd = [], [], [], []
    for fs in FS_ORDER:
        sub = inh[inh.fs == fs]
        pm = sub.groupby("method").pearson_delta.mean()
        qm = sub.groupby("method").PDS_cos.mean()
        dir_mean.append(pm.mean()); dir_sd.append(pm.std(ddof=1))
        pds_mean.append(qm.mean()); pds_sd.append(qm.std(ddof=1))
    dir_mean = np.array(dir_mean); dir_sd = np.array(dir_sd)
    pds_mean = np.array(pds_mean); pds_sd = np.array(pds_sd)

    x = np.arange(len(FS_ORDER), dtype=float)
    fig, ax = S.panel(61.7, 46.3)

    # chance line for PDS, behind the data
    ax.axhline(0.50, color=S.INK, ls=(0, (3, 2)), lw=0.6, zorder=1)
    ax.text(x[-1] + 0.45, 0.505, "chance (0.50)", fontsize=5.2, color=S.INK,
            ha="right", va="bottom")

    # direction (Pearson-delta), filled ink circles
    ax.errorbar(x, dir_mean, yerr=dir_sd, fmt="o", ms=4.0, mfc=C_DIR, mec="white",
                mew=0.5, ecolor=C_DIR, elinewidth=0.7, capsize=2, capthick=0.7,
                zorder=4, clip_on=False)
    # PDS (cosine), open grey squares
    ax.errorbar(x, pds_mean, yerr=pds_sd, fmt="s", ms=3.6, mfc="white", mec=C_PDS,
                mew=0.9, ecolor=C_PDS, elinewidth=0.7, capsize=2, capthick=0.7,
                zorder=4, clip_on=False)

    # direct labels instead of a legend (series are spatially stable and separated)
    ax.text(x[0] - 0.36, dir_mean[0] + dir_sd[0] + 0.022,
            r"Pearson-$\delta$ (direction)", fontsize=5.6, color=C_DIR,
            ha="left", va="bottom")
    ax.text(x[0] - 0.36, pds_mean[0] - pds_sd[0] - 0.022, "PDS (cosine)",
            fontsize=5.6, color=C_PDS, ha="left", va="top")

    ax.set_ylabel("Mean score across six model heads")
    ax.set_ylim(0.38, 0.70)
    ax.set_yticks([0.4, 0.5, 0.6, 0.7])
    ax.tick_params(axis="y", length=2.2)

    ax.set_xlim(-0.5, len(FS_ORDER) - 0.5 + 0.45)
    ax.set_xticks(x)
    ax.set_xticklabels([FS_LABEL[f] for f in FS_ORDER])
    ax.set_xlabel("Input feature space")
    ax.tick_params(axis="x", length=2.2)

    S.despine(ax)

    S.save(fig, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             "fig4e_features"))

    for i, fs in enumerate(FS_ORDER):
        print(f"{fs:10s} pearson_delta={dir_mean[i]:.4f} (sd {dir_sd[i]:.4f})  "
              f"PDS_cos={pds_mean[i]:.4f} (sd {pds_sd[i]:.4f})")


if __name__ == "__main__":
    main()
