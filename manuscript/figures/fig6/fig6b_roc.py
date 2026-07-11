"""Figure 6b - Leave-one-dataset-out ROC of the honest rankability predictor.

Data-direct reconstruction. A logistic regression on a single pilot-estimable,
non-leaky feature (log10 effect size) is trained leave-one-dataset-out (LODO)
and its ROC computed on each held-out dataset. This reproduces the committed
honest AUROCs (results/rankability_predictor_honest.csv, feature_set='effect_size':
Replogle 0.97, Norman 0.93, Adamson 0.93, GATA1 0.95, JAK1 1.00) exactly; the
reconstructed and committed values are printed at run time to confirm.

TP53 and KRAS are uniformly unrankable (single-class labels), so an AUROC for
them is mathematically undefined; they are noted, not plotted.

Recipe matches scripts/figures/rankability_predictor.py: per-perturbation rows at
each perturbation's deepest split-half bin (D.native_rankability), feature =
log10(effect_size), label = rankable, StandardScaler + LogisticRegression
(max_iter=2000, random_state=42), sklearn roc_curve / roc_auc_score.

Run:  python fig6b_roc.py  ->  fig6b_roc.pdf (+ .png)
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_curve, roc_auc_score

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import nm_style as S
import fig6_data as D

# all datasets fed to the pooled LODO training; the two single-class allele
# genes still contribute (all-unrankable) rows to the training pool.
ALL_DATASETS = ["Replogle", "Norman", "Adamson", "GATA1", "JAK1", "TP53", "KRAS"]
SINGLE_CLASS = ["TP53", "KRAS"]


def build() -> dict[str, pd.DataFrame]:
    frames = {}
    for ds in ALL_DATASETS:
        df = D.native_rankability(ds).copy()
        df["log_effect"] = np.log10(df["effect_size"].clip(lower=1e-6))
        df["y"] = df["rankable"].astype(int)
        frames[ds] = df[["log_effect", "y"]].reset_index(drop=True)
    return frames


def main() -> None:
    S.apply_rcparams()
    frames = build()

    committed = D.honest_auroc("effect_size").set_index("held_out")["auroc"].to_dict()

    fig, ax = S.panel(52, 48)

    # chance diagonal
    ax.plot([0, 1], [0, 1], ls=(0, (2, 2)), lw=0.5, color=S.GREY, zorder=1)

    printed = []
    for ho in D.LODO_DATASETS:  # 5 evaluable, in display order
        te = frames[ho]
        tr = pd.concat([frames[d] for d in ALL_DATASETS if d != ho], ignore_index=True)
        sc = StandardScaler().fit(tr[["log_effect"]])
        model = LogisticRegression(max_iter=2000, random_state=42).fit(
            sc.transform(tr[["log_effect"]]), tr["y"]
        )
        p = model.predict_proba(sc.transform(te[["log_effect"]]))[:, 1]
        y = te["y"].values
        fpr, tpr, _ = roc_curve(y, p)
        auc = roc_auc_score(y, p)
        printed.append((ho, auc, committed.get(ho)))

        ax.plot(fpr, tpr, "-", lw=0.9, color=D.DATASET_COLORS[ho],
                label=f"{ho} ({auc:.2f})", zorder=3, solid_capstyle="round")

    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.set_xticks([0, 0.5, 1.0])
    ax.set_yticks([0, 0.5, 1.0])
    ax.set_xlabel("False positive rate")
    ax.set_ylabel("True positive rate")
    ax.set_aspect("equal")
    S.despine(ax)
    ax.tick_params(length=2.2)

    leg = ax.legend(loc="lower right", handlelength=1.2, handletextpad=0.5,
                    labelspacing=0.25, borderpad=0.2, borderaxespad=0.2,
                    fontsize=5.5, title="LODO hold-out (AUC)", frameon=False)
    leg.get_title().set_fontsize(5.5)

    # single-class note (TP53/KRAS have undefined AUROC); placed in the empty
    # wedge below the chance diagonal, clear of both the ROC curves and legend
    ax.text(0.055, 0.24, "TP53, KRAS: single-class\n(AUROC undefined)",
            transform=ax.transAxes, fontsize=5, color=S.GREY,
            va="center", ha="left", linespacing=1.2)

    S.save(fig, "fig6b_roc")

    print("held-out   reconstructed   committed   match(<0.01)")
    for ho, auc, comm in printed:
        print(f"  {ho:9s} {auc:.4f}        {comm:.4f}     {abs(auc - comm) < 0.01}")


if __name__ == "__main__":
    main()
