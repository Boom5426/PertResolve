"""Fig. 3b — representative closed and open detection windows."""
from __future__ import annotations

import numpy as np

import fig3_common as S

# Publication contract inherited from fig3_common: Arial/Helvetica sans-serif;
# svg.fonttype="none"; pdf.fonttype=42; exports .svg, .pdf, .png and .tiff at dpi=600.
# Final assembled figure target: width_mm = 183.
W_MM, H_MM = 83.0, 40.0


def distribution(ax, x, values, color):
    values = np.asarray(values, dtype=float)
    parts = ax.violinplot(
        values,
        positions=[x],
        widths=0.62,
        showmeans=False,
        showmedians=False,
        showextrema=False,
        bw_method=0.25,
    )
    for body in parts["bodies"]:
        body.set_facecolor(color)
        body.set_edgecolor("none")
        body.set_alpha(0.24)
    lo, med, hi = np.quantile(values, [0.025, 0.5, 0.975])
    ax.plot([x, x], [lo, hi], color=color, lw=1.0, solid_capstyle="round", zorder=4)
    ax.plot(x, med, "o", ms=3.1, mfc=color, mec="white", mew=0.45, zorder=5)
    return med


def main() -> None:
    S.apply_style()
    df = S.read_csv("fig3b_selfnull_dist.csv")
    fig = S.figure(W_MM, H_MM)
    axes = [
        fig.add_axes([0.12, 0.16, 0.32, 0.60]),
        fig.add_axes([0.61, 0.16, 0.32, 0.60]),
    ]
    specs = [
        ("TP53", "closed window", axes[0]),
        ("JAK1", "open window", axes[1]),
    ]

    for gene, regime, ax in specs:
        sub = df[df["gene"] == gene]
        variant = str(sub["variant"].iloc[0])
        color = S.GENE_COLORS[gene]
        m_self = distribution(ax, 0, sub["D_self"], color)
        m_null = distribution(ax, 1, sub["D_null"], S.MID_GREY)
        ratio = m_self / m_null

        ymax = max(sub["D_self"].max(), sub["D_null"].max()) * 1.12
        ax.set_ylim(0, ymax)
        ax.set_xlim(-0.55, 1.55)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(
            [r"$D_\mathrm{self}$" + "\nnoise", r"$D_\mathrm{null}$" + "\nvariant–WT"],
            linespacing=1.05,
        )
        ax.set_ylabel("Energy distance", labelpad=1.5)
        ax.tick_params(axis="x", length=0, pad=2)
        ax.tick_params(axis="y", pad=1.5)
        S.despine(ax)

        pos = ax.get_position()
        cx = pos.x0 + pos.width / 2
        fig.text(
            cx,
            0.94,
            f"{gene} {variant}",
            ha="center",
            va="top",
            fontsize=6.8,
            color=color,
            fontweight="bold",
        )
        fig.text(
            cx,
            0.888,
            f"{regime}  ·  R={ratio:.2f}",
            ha="center",
            va="top",
            fontsize=5.9,
            color=color,
        )
    fig.text(
        0.515,
        0.805,
        "200 split-half resamples · panels use independent y-scales",
        ha="center",
        va="top",
        fontsize=5.2,
        color=S.GREY,
    )

    S.save(fig, "fig3b_representative_windows")


if __name__ == "__main__":
    main()
