"""Fig. 3g — matched-cohort rankability across split-half depths.

The same primary S > W definition as Fig. 3h is used. Within each gene, only
variants evaluable at 150 cells per half are retained so every point on that
gene's curve refers to the same variants. Wilson intervals reflect the
binomial uncertainty of the rankable fraction.
"""
from __future__ import annotations

import numpy as np

import fig3_common as S

# Publication contract inherited from fig3_common: Arial/Helvetica sans-serif;
# svg.fonttype="none"; pdf.fonttype=42; exports .svg, .pdf, .png and .tiff at dpi=600.
# Final assembled figure target: width_mm = 183.
W_MM, H_MM = 91.0, 36.0
DEPTHS = [30, 50, 100, 150]


def main() -> None:
    S.apply_style()
    table = S.rank_table()
    table = table[
        (table["metric"] == "edist")
        & (table["space"] == "pca")
        & (table["dataset"].isin(S.GENE_ORDER))
    ].copy()

    fig = S.figure(W_MM, H_MM)
    gs = fig.add_gridspec(
        2,
        2,
        left=0.13,
        right=0.98,
        bottom=0.18,
        top=0.83,
        hspace=0.26,
        wspace=0.24,
    )
    axes = [fig.add_subplot(gs[i // 2, i % 2]) for i in range(4)]

    for ax, gene in zip(axes, S.GENE_ORDER):
        sub = table[table["dataset"] == gene]
        fixed = set(sub.loc[sub["n_work"] == max(DEPTHS), "perturbation"])
        fixed_sub = sub[
            sub["perturbation"].isin(fixed) & sub["n_work"].isin(DEPTHS)
        ].copy()
        n = len(fixed)
        if n == 0:
            raise ValueError(f"No fixed depth cohort for {gene}")

        values, lows, highs = [], [], []
        for depth in DEPTHS:
            at_depth = fixed_sub[fixed_sub["n_work"] == depth]
            if len(at_depth) != n:
                raise ValueError(f"{gene} depth {depth}: expected {n}, observed {len(at_depth)}")
            k = int(at_depth["rankable"].sum())
            lo, hi = S.wilson_interval(k, n)
            values.append(100.0 * k / n)
            lows.append(100.0 * lo)
            highs.append(100.0 * hi)

        color = S.GENE_COLORS[gene]
        for reference in (0, 50, 100):
            ax.axhline(reference, color=S.PALE_GREY, lw=0.5, zorder=0)
        ax.plot(DEPTHS, values, "-o", color=color, lw=1.0, ms=3.1, zorder=3)
        for x, y, lo, hi in zip(DEPTHS, values, lows, highs):
            ax.plot([x, x], [lo, hi], color=color, lw=0.8, alpha=0.8, zorder=2)
            ax.plot([x - 3, x + 3], [lo, lo], color=color, lw=0.65, alpha=0.8)
            ax.plot([x - 3, x + 3], [hi, hi], color=color, lw=0.65, alpha=0.8)

        label_x, label_y, label_ha, label_va = (
            (0.98, 0.08, "right", "bottom")
            if gene == "JAK1"
            else (0.02, 0.96, "left", "top")
        )
        ax.text(
            label_x,
            label_y,
            f"{gene}  ·  fixed n={n}",
            transform=ax.transAxes,
            ha=label_ha,
            va=label_va,
            fontsize=6.1,
            fontweight="bold",
            color=color,
        )
        ax.set_xlim(24, 156)
        ax.set_ylim(-4, 108)
        ax.set_xticks(DEPTHS)
        ax.set_yticks([0, 50, 100])
        ax.tick_params(axis="both", pad=1.2)
        S.despine(ax)

    for ax in axes[:2]:
        ax.set_xticklabels([])
        ax.tick_params(axis="x", length=0)
    for ax in (axes[1], axes[3]):
        ax.set_yticklabels([])
        ax.tick_params(axis="y", length=0)

    fig.text(
        0.55,
        0.955,
        r"Primary rankability criterion: $S>W$",
        ha="center",
        va="center",
        fontsize=6.2,
        color=S.INK,
    )
    fig.text(
        0.55,
        0.89,
        "points: fraction · lines: 95% Wilson CI",
        ha="center",
        va="center",
        fontsize=5.1,
        color=S.GREY,
    )
    fig.text(0.55, 0.055, "Cells per half", ha="center", va="center", fontsize=7.0)
    fig.text(
        0.035,
        0.515,
        "Rankable variants (%)",
        ha="center",
        va="center",
        rotation=90,
        fontsize=7.0,
    )
    S.save(fig, "fig3g_depth_rankability")


if __name__ == "__main__":
    main()
