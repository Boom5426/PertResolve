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

Encoding: every dataset is a filled circle. Marker SHAPE is reserved for the metric
in panel g (circle = recover full order, diamond = select best model), so f must not
also use shape for dataset class; class is carried by direct labels, and the four
allele genes keep their house gene colours (TP53 blue, KRAS orange, GATA1 purple,
JAK1 green) exactly as in panels b, c and g. Axis limits and ticks are identical to
panel g so the two panels share one oracle-ceiling axis, but the title and the
x-axis label state what each panel's unit of analysis is, because f and g are NOT
the same quantity: f is one point per DATASET at a matched shallow depth (Methods,
"External benchmark calibration": 50 cells per half, perturbations with >=100 cells),
whereas g is one point per GENE x SEQUENCING DEPTH for the four allele genes and
carries a second metric. The two therefore disagree by construction where they seem
to overlap (GATA1 0.512 here versus 0.989 at depth 50 in g), which is exactly why
the labels must not read as interchangeable.

Panel f is the one point at (0.49, 0.01/0.02) where TP53 and KRAS coincide: their
ceilings (0.492, 0.493) and resolutions (0.014, 0.018) differ by less than a marker
radius, so the two markers overlap. No point is offset; the two direct labels are
moved clear of the datum instead.

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
# TP53 and KRAS sit on top of each other at the floor: label them side by side to
# the RIGHT of the shared datum so no text sits on the marker.
LABELS = {
    "allele_TP53": ("TP53", 0.014, 0.068, "left"),
    "allele_KRAS": ("KRAS", 0.014, -0.012, "left"),
    "allele_GATA1": ("GATA1", -0.010, 0.060, "right"),
    "allele_JAK1": ("JAK1", -0.008, -0.070, "right"),
    "sciPlex": ("sci-Plex", 0.0, -0.080, "center"),
}


def main() -> None:
    S.apply_rcparams()
    df = pd.read_csv(S.repo_root() / "results" / "benchmark_resolution" / "summary.csv")
    df = df.sort_values("oracle_ceiling").reset_index(drop=True)

    # height net of the panel title, so the placed height is unchanged (57.9 mm)
    fig, ax = S.panel(96.9, 61.2)

    # regime shading (annotation only)
    ax.axhspan(0.90, 1.06, color="#F4F4F4", alpha=1.0, zorder=0)
    ax.axhspan(-0.05, 0.50, color="#EDEDED", alpha=1.0, zorder=0)
    ax.text(0.468, 0.985, "reliable ranking", fontsize=5.5, color=S.GREY,
            va="center", ha="left")
    ax.text(0.972, 0.075, "near floor", fontsize=5.5, color=S.GREY,
            va="center", ha="right")

    # trend line through sorted points
    ax.plot(df.oracle_ceiling, df.resolution_P_recover_order, "-", color="#C9C9C9",
            lw=0.8, zorder=1)

    for _, r in df.iterrows():
        ds = r.dataset
        x, y = r.oracle_ceiling, r.resolution_P_recover_order
        if ds in ALLELE:
            color, sz = S.GENE_COLORS[ALLELE[ds]], 30
        else:  # sci-Plex and the gene-level atlases: not allele-level, quiet grey
            color, sz = S.GREY, 26
        ax.scatter([x], [y], s=sz, marker="o", color=color, edgecolor="white",
                   linewidths=0.4, zorder=4)
        if ds in LABELS:
            lab, dx, dy, ha = LABELS[ds]
            col = S.GENE_COLORS[ALLELE[ds]] if ds in ALLELE else S.INK
            ax.annotate(lab, (x, y), xytext=(x + dx, y + dy), fontsize=6, color=col,
                        ha=ha, va="center")

    # annotate the gene-level cluster once
    ax.annotate("gene-level atlases\n(Norman, Replogle,\nAdamson, VCC)",
                xy=(0.780, 1.0), xytext=(0.660, 0.72), fontsize=5.5, color=S.GREY,
                ha="left", va="center",
                arrowprops=dict(arrowstyle="-", lw=0.4, color=S.GREY,
                                shrinkA=1, shrinkB=3))

    ax.set_title("Nine benchmark datasets, matched depth", fontsize=6.5,
                 loc="left", pad=3)
    ax.set_xlabel("Oracle ceiling per dataset (50 cells per half)")
    ax.set_ylabel("P(recover full order)")
    ax.set_xlim(0.46, 0.98)
    ax.set_ylim(-0.05, 1.06)
    ax.set_xticks([0.5, 0.6, 0.7, 0.8, 0.9])
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    S.despine(ax)
    ax.tick_params(length=2.2)

    S.save(fig, "fig5f_resolution")
    print(df[["dataset", "oracle_ceiling", "resolution_P_recover_order"]].to_string(index=False))


if __name__ == "__main__":
    main()
