"""Figure 6d: prospective pilot-to-full-depth validation (quantitative half).

Message: a pilot of 25-50 cells predicts which perturbations are rankable at
higher depth on DISJOINT cells. Two predictors (learned effect-size logistic,
training-free pilot-SNR threshold) both beat the label-permutation null (0.49),
and out-of-fold predicted probability tracks the observed rankable rate
(near-diagonal calibration inset).

Numbers are the adversarially-verified documented values:
  learned LODO AUROC (T=100): Adamson 0.85, Norman 0.94, Replogle 0.99 (mean 0.93)
  training-free pilot-SNR:     Adamson 0.84, Norman 0.93, Replogle 0.98
  null (label permutation):    0.49
  calibration bins:            T100['calibration'] (pooled, monotone, near-diagonal)
Per-dataset only; the pooled mechanistic 0.978 is base-rate-inflated and NOT shown.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import nm_style as S
import fig6_data as D   # in the same fig6/ dir

import numpy as np

S.apply_rcparams()

# ---- data (documented, verified against README + T100 JSON) ---------------
DATASETS = ["Adamson", "Norman", "Replogle"]      # the three balanced datasets
LEARNED = {"Adamson": 0.85, "Norman": 0.94, "Replogle": 0.99}   # README headline
# Training-free pilot-SNR AUROC. Source: results/pilot_validation/README.md
# ("Mechanistic pilot-SNR (no training), per dataset: 0.84 / 0.93 / 0.98", listed
# in the same Adamson / Norman / Replogle order as the LODO line above).
# Traceability re-checked 2026-07-27 by recomputing per dataset directly from the
# committed per-perturbation files results/pilot_validation/<DS>_pilot.csv, using
# the recipe of results/pilot_validation/pilot_validate.py (score = pilot_snr_50,
# label = rank_T100, sklearn.metrics.roc_auc_score, no training):
#   Adamson 0.839 (n=92, pos=49), Norman 0.926 (n=177, pos=146),
#   Replogle 0.977 (n=239, pos=161)  ->  0.84 / 0.93 / 0.98 as plotted.
# Note that pilot_validate.py itself stores only the POOLED mechanistic AUROC
# (0.978) in pilot_validation_summary.json; that pooled value is base-rate
# inflated and is deliberately NOT plotted here.
SNR = {"Adamson": 0.84, "Norman": 0.93, "Replogle": 0.98}       # training-free
MEAN_LEARNED = 0.93
NULL = 0.49

# cross-check learned values against the committed JSON (rounded), fail loud otherwise
_j = D.pilot_summary()["T100"]["learned_lodo_auroc"]
for ds in DATASETS:
    assert round(_j[ds][0], 2) == LEARNED[ds], (ds, _j[ds][0], LEARNED[ds])

# calibration: [lo, hi, observed_frac, n] -> midpoint vs observed
CALIB = D.pilot_summary()["T100"]["calibration"]
cx = np.array([(lo + hi) / 2.0 for lo, hi, _, _ in CALIB])
cy = np.array([obs for _, _, obs, _ in CALIB])

# ---- figure ---------------------------------------------------------------
fig, ax = S.panel(64, 50)

xpos = np.arange(len(DATASETS))
JIT = 0.16  # split the two predictors left/right within each dataset slot

# null baseline
ax.axhline(NULL, ls=(0, (3, 2)), lw=0.6, color=S.GREY, zorder=1)
ax.text(len(DATASETS) - 0.44, NULL + 0.008, "null 0.49", ha="right",
        va="bottom", fontsize=5.5, color=S.GREY)

# mean-of-learned guide
ax.axhline(MEAN_LEARNED, ls=(0, (1, 1.5)), lw=0.6, color=S.INK, zorder=1)
ax.text(-0.42, MEAN_LEARNED + 0.006, "mean 0.93", ha="left", va="bottom",
        fontsize=5.5, color=S.INK)

# FIGURE-WIDE FILL RULE (panels d, e, f): one mark (circle); filled = the primary /
# headline condition, open = the comparison condition. Here filled = learned
# predictor, open = training-free predictor. The former square marker in this panel
# was the figure's only competing shape encoding and has been removed.
for i, ds in enumerate(DATASETS):
    c = D.DATASET_COLORS[ds]
    ax.plot(i - JIT, LEARNED[ds], marker="o", ms=4.4, mfc=c, mec=c,
            mew=0.6, ls="none", zorder=3)
    ax.plot(i + JIT, SNR[ds], marker="o", ms=4.4, mfc="white", mec=c,
            mew=1.0, ls="none", zorder=3)

ax.set_xlim(-0.6, len(DATASETS) - 0.4)
ax.set_ylim(0.40, 1.02)
ax.set_xticks(xpos)
ax.set_xticklabels(DATASETS)
ax.set_yticks([0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
ax.set_ylabel("Prospective AUROC\n(disjoint higher-depth cells)")
S.despine(ax)
ax.tick_params(length=2.2)

# LEGEND ECONOMY: no detached legend. The two predictors are spatially stable
# (learned always left of the tick, training-free always right), so they are
# direct-labelled once on the leftmost (Adamson) pair with hairline leaders.
ax.annotate("learned (effect size)", xy=(-0.185, 0.838), xytext=(-0.55, 0.792),
            fontsize=5.5, color=S.INK, ha="left", va="center",
            arrowprops=dict(arrowstyle="-", lw=0.4, color=S.GREY,
                            shrinkA=1.0, shrinkB=2.0))
ax.annotate("training-free (SNR)", xy=(0.175, 0.852), xytext=(0.40, 0.888),
            fontsize=5.5, color=S.INK, ha="left", va="center",
            arrowprops=dict(arrowstyle="-", lw=0.4, color=S.GREY,
                            shrinkA=1.0, shrinkB=2.0))

# panel note
ax.text(0.0, 1.015, "pilot 25-50 cells", transform=ax.transAxes,
        ha="left", va="bottom", fontsize=5.5, color=S.INK)

# ---- calibration inset ----------------------------------------------------
axi = ax.inset_axes([0.10, 0.09, 0.34, 0.40])
axi.plot([0, 1], [0, 1], ls=(0, (2, 1.5)), lw=0.5, color=S.GREY, zorder=1)
axi.plot(cx, cy, marker="o", ms=2.6, mfc=S.INK, mec=S.INK, mew=0.0,
         lw=0.6, color=S.INK, zorder=2)
axi.set_xlim(-0.03, 1.03)
axi.set_ylim(-0.03, 1.03)
axi.set_xticks([0, 1])
axi.set_yticks([0, 1])
axi.tick_params(length=1.6, pad=1.2, labelsize=6)
axi.set_xlabel("pred.", fontsize=6.0, labelpad=1.0)
axi.set_ylabel("obs.", fontsize=6.0, labelpad=1.0)
axi.set_title("calibration", fontsize=6.0, pad=1.5)
S.despine(axi)
for sp in ("left", "bottom"):
    axi.spines[sp].set_linewidth(0.4)

S.save(fig, os.path.join(os.path.dirname(os.path.abspath(__file__)), "fig6d_prospective"))
print("saved fig6d_prospective")
