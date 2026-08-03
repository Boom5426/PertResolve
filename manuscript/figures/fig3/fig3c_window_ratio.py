"""Fig. 3c — per-variant detection-window ratios and gene-level bootstrap CIs."""
from __future__ import annotations

import numpy as np

import fig3_common as S

# Publication contract inherited from fig3_common: Arial/Helvetica sans-serif;
# svg.fonttype="none"; pdf.fonttype=42; exports .svg, .pdf, .png and .tiff at dpi=600.
# Final assembled figure target: width_mm = 183.
W_MM, H_MM = 78.0, 43.0


def main() -> None:
    S.apply_style()
    ci = S.ratio_ci()
    fig = S.figure(W_MM, H_MM)
    ax = fig.add_axes([0.14, 0.16, 0.82, 0.76])

    for x, gene in enumerate(S.GENE_ORDER):
        native = S.native_rankability(gene)
        values = native["window_ratio"].to_numpy(float)
        values = values[np.isfinite(values)]
        color = S.GENE_COLORS[gene]

        parts = ax.violinplot(
            values,
            positions=[x],
            widths=0.72,
            showmeans=False,
            showmedians=False,
            showextrema=False,
            bw_method=0.22,
        )
        for body in parts["bodies"]:
            vertices = body.get_paths()[0].vertices
            vertices[:, 0] = np.clip(vertices[:, 0], x, np.inf)
            body.set_facecolor(color)
            body.set_edgecolor("none")
            body.set_alpha(0.22)

        order = np.arange(len(values))
        phase = ((order * 37) % max(len(values), 1)) / max(len(values) - 1, 1)
        jitter = x - 0.10 - 0.24 * phase
        ax.scatter(
            jitter,
            values,
            s=4.0,
            color=color,
            alpha=0.42,
            linewidths=0,
            zorder=3,
        )

        mean = float(ci[gene]["mean"])
        lo = float(ci[gene]["lo"])
        hi = float(ci[gene]["hi"])
        sx = x + 0.08
        ax.plot([sx, sx], [lo, hi], color=S.INK, lw=1.05, zorder=5)
        ax.plot(
            sx,
            mean,
            marker="D",
            ms=3.2,
            mfc="white",
            mec=S.INK,
            mew=0.8,
            zorder=6,
        )
        ax.text(
            x + 0.20,
            mean,
            f"{mean:.2f}",
            ha="left",
            va="center",
            fontsize=6.1,
            color=S.INK,
        )

    ax.axhline(1.0, color=S.GREY, lw=0.75, ls=(0, (3, 2)), zorder=1)
    ax.text(
        -0.48,
        1.025,
        r"$R=1$",
        ha="left",
        va="bottom",
        fontsize=5.6,
        color=S.GREY,
    )
    ax.text(
        3.48,
        1.025,
        "replicate-noise\nfloor",
        ha="right",
        va="bottom",
        fontsize=5.4,
        color=S.GREY,
        linespacing=1.0,
        bbox=dict(facecolor="white", edgecolor="none", pad=0.7, alpha=0.88),
    )
    ax.text(
        0.01,
        0.985,
        "diamond: mean · line: 95% bootstrap CI",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=5.2,
        color=S.GREY,
    )

    labels = []
    for gene in S.GENE_ORDER:
        labels.append(f"{gene}\nn={len(S.native_rankability(gene))}")
    ax.set_xticks(range(4))
    ax.set_xticklabels(labels, linespacing=1.15)
    S.gene_ticklabels(ax, "x")
    ax.set_xlim(-0.52, 3.55)
    ax.set_ylim(0, 1.42)
    ax.set_yticks(np.arange(0, 1.41, 0.2))
    ax.set_ylabel(r"Window ratio, $D_\mathrm{self}/D_\mathrm{null}$", labelpad=2)
    ax.tick_params(axis="x", length=0, pad=3)
    ax.tick_params(axis="y", pad=1.5)
    S.despine(ax)
    S.save(fig, "fig3c_window_ratio")


if __name__ == "__main__":
    main()
