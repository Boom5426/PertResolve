"""Figure 1c — per-variant sampling depth.

All non-WT variant conditions are shown.  Raw points are retained; the pale
half-violin is a descriptive density, the vertical segment is the IQR and the
horizontal segment is the median.  The log axis is labelled explicitly because
nominal cell count is not itself a measurement-power estimate.

Run: python fig1c_depth.py
Output: fig1c_depth.svg/.pdf/.png

Imported export contract: Arial; svg.fonttype: "none"; pdf.fonttype: 42;
outputs .svg, .pdf, .png and .tiff at dpi=600. Final assembled figure target:
width_mm = 183.
"""
from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import nm_style as S

from pathlib import Path

HERE = Path(__file__).resolve().parent

W_MM, H_MM = 58.0, 42.0

def main() -> None:
    S.apply_rcparams()
    data = S.load_bench(exclude_wt=True, corrected=True)
    fig, ax = S.panel(W_MM, H_MM)
    fig.subplots_adjust(left=0.22, right=0.98, bottom=0.23, top=0.95)

    for i, gene in enumerate(S.GENE_ORDER):
        values = data.loc[data.gene == gene, "n_cells"].to_numpy(dtype=float)
        if np.any(values <= 0):
            raise ValueError(f"{gene}: log-scale cell counts must be positive")
        log_values = np.log10(values)
        colour = S.GENE_COLORS[gene]

        violin = ax.violinplot(log_values, positions=[i], widths=0.82, showextrema=False)
        for body in violin["bodies"]:
            vertices = body.get_paths()[0].vertices
            vertices[:, 0] = np.clip(vertices[:, 0], i, np.inf)
            body.set_facecolor(colour)
            body.set_edgecolor("none")
            body.set_alpha(0.16)

        index = np.arange(len(values), dtype=float)
        jitter = i - 0.08 - 0.20 * ((index * 0.61803398875) % 1.0)
        ax.scatter(
            jitter,
            log_values,
            s=2.1,
            color=colour,
            alpha=0.55,
            linewidths=0,
            zorder=3,
        )

        q1, median, q3 = np.quantile(values, [0.25, 0.5, 0.75])
        x_summary = i + 0.07
        ax.plot(
            [x_summary, x_summary],
            [np.log10(q1), np.log10(q3)],
            color=S.INK,
            lw=1.05,
            solid_capstyle="round",
            zorder=5,
        )
        ax.plot(
            [i - 0.24, i + 0.30],
            [np.log10(median)] * 2,
            color=S.INK,
            lw=1.05,
            solid_capstyle="round",
            zorder=6,
        )
        ax.text(
            x_summary,
            np.log10(q3) + 0.065,
            f"{median:.0f}",
            ha="center",
            va="bottom",
            fontsize=5.8,
            bbox={"facecolor": "white", "edgecolor": "none", "pad": 0.15,
                  "alpha": 0.92},
            zorder=8,
        )
        ax.text(
            i,
            -0.155,
            f"n={len(values)}",
            transform=ax.get_xaxis_transform(),
            ha="center",
            va="top",
            fontsize=5.2,
            color=S.GREY,
        )

    ax.set_xticks(range(len(S.GENE_ORDER)))
    ax.set_xticklabels(S.GENE_ORDER)
    for label, gene in zip(ax.get_xticklabels(), S.GENE_ORDER):
        label.set_color(S.INK)
        label.set_fontweight("bold")

    ax.set_xlim(-0.55, 3.55)
    ax.set_ylim(np.log10(25), np.log10(3000))
    ticks = [50, 100, 300, 1000, 3000]
    ax.set_yticks(np.log10(ticks))
    ax.set_yticklabels([f"{tick:,}" for tick in ticks])
    ax.set_ylabel("Cells per variant condition")
    ax.text(
        0.99,
        1.005,
        "log scale · sampling depth only",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=5.0,
        color=S.GREY,
    )
    S.despine(ax)
    ax.tick_params(length=2.0)

    S.save(fig, HERE / "fig1c_depth", exact=True,
           formats=("pdf", "png", "svg", "tiff"))

    summary = (
        data.groupby("gene")["n_cells"]
        .agg(["median", "sum", "count"])
        .reindex(S.GENE_ORDER)
    )
    print(summary.to_string())


if __name__ == "__main__":
    main()
