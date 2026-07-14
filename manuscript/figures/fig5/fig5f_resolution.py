"""Figure 5f - Benchmark resolution rises sharply with the oracle ceiling.

Data-direct. Source (committed): results/benchmark_resolution/summary.csv, whose
per-dataset (oracle_ceiling, resolution_P_recover_order) reproduce Supplementary
Table 5 exactly. Probability that the benchmark recovers the true ranking of graded
synthetic predictors (resolution) as a function of the measurement window (oracle
ceiling), across 9 datasets. Below a ceiling near 0.5 (TP53/KRAS at the floor)
resolution collapses; above ~0.6 (GATA1 partial, then JAK1 and the gene-/drug-level
screens) it approaches 1.

Note: the ceiling here is the >=100-cell / 50-per-half benchmark-resolution version
(matches SI Table 5), not the native-depth per-gene oracle of the main-text table.

Run:  python fig5f_resolution.py  ->  fig5f_resolution.pdf (+ .png)
"""
from __future__ import annotations

import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import nm_style as S

ALLELE = {"allele_TP53": "TP53", "allele_KRAS": "KRAS",
          "allele_GATA1": "GATA1", "allele_JAK1": "JAK1"}
# nice labels + manual offsets (dx, dy in data units) for direct labels
LABELS = {
    "allele_TP53": ("TP53", -0.004, 0.075, "right"),
    "allele_KRAS": ("KRAS", 0.006, 0.075, "left"),
    "allele_GATA1": ("GATA1", -0.010, 0.055, "right"),
    "allele_JAK1": ("JAK1", -0.010, -0.055, "right"),
    "sciPlex": ("sci-Plex", 0.0, -0.075, "center"),
}


def main() -> None:
    S.apply_rcparams()
    df = pd.read_csv(S.repo_root() / "results" / "benchmark_resolution" / "summary.csv")
    df = df.sort_values("oracle_ceiling").reset_index(drop=True)

    fig, ax = S.panel(74, 54)

    # regime shading (annotation only)
    ax.axhspan(0.90, 1.06, color="#F7F7F7", alpha=0.6, zorder=0)
    ax.axhspan(-0.06, 0.50, color="#E4E4E4", alpha=0.6, zorder=0)
    ax.text(0.905, 1.0, "reliable ranking", fontsize=5, color=S.GREY, va="center")
    ax.text(0.905, 0.16, "near floor", fontsize=5, color=S.GREY, va="center")

    # trend line through sorted points
    ax.plot(df.oracle_ceiling, df.resolution_P_recover_order, "-", color="#C9C9C9",
            lw=0.8, zorder=1)

    for _, r in df.iterrows():
        ds = r.dataset
        x, y = r.oracle_ceiling, r.resolution_P_recover_order
        if ds in ALLELE:
            color, marker, sz = S.GENE_COLORS[ALLELE[ds]], "o", 34
        elif ds == "sciPlex":
            color, marker, sz = "#6E6E6E", "D", 24
        else:  # gene-level atlases
            color, marker, sz = "#9AA0A6", "s", 26
        ax.scatter([x], [y], s=sz, marker=marker, color=color, edgecolor="white",
                   linewidths=0.4, zorder=4)
        if ds in LABELS:
            lab, dx, dy, ha = LABELS[ds]
            col = S.GENE_COLORS[ALLELE[ds]] if ds in ALLELE else S.INK
            ax.annotate(lab, (x, y), xytext=(x + dx, y + dy), fontsize=5.5, color=col,
                        ha=ha, va="center")

    # annotate the gene-level cluster once
    ax.annotate("gene-level atlases\n(Norman, Replogle,\nAdamson, VCC)",
                xy=(0.80, 1.0), xytext=(0.66, 0.66), fontsize=5, color=S.GREY,
                ha="left", va="center",
                arrowprops=dict(arrowstyle="-", lw=0.4, color=S.GREY))

    ax.set_xlabel("Oracle ceiling (measurement window)")
    ax.set_ylabel("P(recover true model order)")
    ax.set_xlim(0.46, 0.92)
    ax.set_ylim(-0.06, 1.08)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    S.despine(ax)
    ax.tick_params(length=2.2)

    S.save(fig, "fig5f_resolution")
    print(df[["dataset", "oracle_ceiling", "resolution_P_recover_order"]].to_string(index=False))


if __name__ == "__main__":
    main()
