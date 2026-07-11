"""Figure 4e - Feature space swap saturates direction but leaves PDS at chance.

Data-direct. Source (committed): results/results_v4_exttheta.csv (the per-variant
breakdown grid, D.exttheta()). The 18 in-house feature-model heads
(Ridge/Lasso/RF/GBoost/KNN/MLP x {theta, ESM, ESM+theta}) are grouped by
D.feature_space into three input representations. Per feature space we plot the
mean pearson_delta (direction) and the mean PDS_cos (perturbation-direction
score), with error bars = SD across the 6 model heads in that space.

Verified firsthand (mean over the 6 heads, i.e. mean of per-method means):
  feature space   pearson_delta         PDS_cos
  theta           0.6005 (sd 0.0278)    0.4656 (sd 0.0168)
  ESM             0.6168 (sd 0.0374)    0.4540 (sd 0.0104)
  ESM+theta       0.6081 (sd 0.0278)    0.4625 (sd 0.0221)

Message: swapping the input representation (biophysical theta -> ESM protein
language-model embeddings -> both) moves direction only marginally around a
saturated ~0.6, while PDS stays pinned near the 0.50 chance line. Richer features
buy no perturbation-direction skill.

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

C_DIR = S.GENE_COLORS["JAK1"]   # green = direction (pearson_delta)
C_PDS = S.GREY                  # grey  = PDS_cos (near chance)


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
    fig, ax = S.panel(52, 48)
    ax2 = ax.twinx()

    # --- left axis: direction (pearson_delta) ---
    ax.errorbar(x, dir_mean, yerr=dir_sd, fmt="o", ms=4.2, mfc=C_DIR, mec="white",
                mew=0.5, ecolor=C_DIR, elinewidth=0.7, capsize=2, capthick=0.7,
                zorder=4, clip_on=False)
    # --- right axis: PDS (near chance) ---
    ax2.errorbar(x, pds_mean, yerr=pds_sd, fmt="s", ms=3.8, mfc="white", mec=C_PDS,
                 mew=0.9, ecolor=C_PDS, elinewidth=0.7, capsize=2, capthick=0.7,
                 zorder=4, clip_on=False)

    # chance line for PDS (right axis)
    ax2.axhline(0.50, color=S.INK, ls=(0, (3, 2)), lw=0.6, zorder=1)
    ax2.text(x[-1] + 0.32, 0.515, "chance", fontsize=5, color=S.INK,
             ha="right", va="bottom")

    # --- left axis cosmetics ---
    ax.set_ylabel(r"Direction (Pearson $\Delta$)", color=C_DIR)
    ax.set_ylim(0.30, 0.75)
    ax.set_yticks([0.3, 0.4, 0.5, 0.6, 0.7])
    ax.tick_params(axis="y", colors=C_DIR, length=2.2)
    ax.spines["left"].set_color(C_DIR)

    # --- right axis cosmetics ---
    ax2.set_ylabel("PDS (cosine)", color=C_PDS, rotation=270, labelpad=10)
    ax2.set_ylim(0.30, 0.75)
    ax2.set_yticks([0.3, 0.4, 0.5, 0.6, 0.7])
    ax2.tick_params(axis="y", colors=C_PDS, length=2.2)
    ax2.spines["right"].set_visible(True)
    ax2.spines["right"].set_color(C_PDS)
    ax2.spines["top"].set_visible(False)

    # --- x axis ---
    ax.set_xlim(-0.5, len(FS_ORDER) - 0.5 + 0.15)
    ax.set_xticks(x)
    ax.set_xticklabels([FS_LABEL[f] for f in FS_ORDER])
    ax.set_xlabel("Input feature space")
    ax.tick_params(axis="x", length=2.2)

    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.spines["left"].set_visible(True)
    ax.spines["bottom"].set_visible(True)

    # --- inline legend (no frame, no black edges) ---
    ax.plot([], [], "o", ms=4.2, mfc=C_DIR, mec="white", mew=0.5,
            label=r"Direction (Pearson $\Delta$)")
    ax.plot([], [], "s", ms=3.8, mfc="white", mec=C_PDS, mew=0.9,
            label="PDS (cosine)")
    leg = ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.16), ncol=1,
                    handletextpad=0.4, borderpad=0.0, labelspacing=0.3,
                    fontsize=5.5)
    leg.get_frame().set_linewidth(0.0)

    S.save(fig, "fig4e_features")

    for i, fs in enumerate(FS_ORDER):
        print(f"{fs:10s} pearson_delta={dir_mean[i]:.4f} (sd {dir_sd[i]:.4f})  "
              f"PDS_cos={pds_mean[i]:.4f} (sd {pds_sd[i]:.4f})")


if __name__ == "__main__":
    main()
