"""Figure 1f — data-direct generalization regimes.

Each mini-bar shows the realized held-out fraction for one gene, read directly
from the corrected benchmark split columns.  Five splits are evaluated in the
main benchmark.  Cross-gene is retained for design completeness but is
explicitly marked as defined-only/not-scored because the datasets do not share
a decoder gene space.

Run: python fig1f_splits.py
Output: fig1f_splits.svg/.pdf/.png

Imported export contract: Arial; svg.fonttype: "none"; pdf.fonttype: 42;
outputs .svg, .pdf, .png and .tiff at dpi=600. Final assembled figure target:
width_mm = 183.
"""
from __future__ import annotations

import os
import sys

from matplotlib.patches import FancyBboxPatch, Rectangle

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import nm_style as S

from pathlib import Path

HERE = Path(__file__).resolve().parent

W_MM, H_MM = 100.0, 37.0
CARD_W, CARD_H = 31.8, 13.9
SPLITS = [
    ("Random", "predefined holdout", "split1_role", True),
    ("Positional", "C-terminal half", "split2_role", True),
    ("Mechanistic", "hotspot / switch", "split3_role", True),
    ("Cross-gene", "GATA1 + JAK1 held out · not scored", "split4_role", False),
    ("Low-depth", "<200 cells", "split5_role", True),
    ("Compatibility", "prior evaluation protocol", "split6_role", True),
]


def split_counts(data, column):
    result = {}
    for gene in S.GENE_ORDER:
        sub = data[data.gene == gene]
        role = sub[column].astype(str).str.lower()
        result[gene] = (int((role == "test").sum()), int(len(sub)))
    return result


def draw_card(ax, x, y, title, description, counts, evaluated):
    if not evaluated:
        ax.add_patch(
            FancyBboxPatch(
                (x - 0.25, y - 0.2),
                CARD_W + 0.5,
                CARD_H + 0.35,
                boxstyle="round,pad=0.1,rounding_size=0.5",
                facecolor="white",
                edgecolor="#B9BEC2",
                lw=0.45,
                linestyle=(0, (2.0, 1.4)),
                zorder=0,
            )
        )
    ax.text(x, y + CARD_H - 1.5, title, fontsize=6.0, fontweight="bold",
            ha="left", va="center", color=S.INK if evaluated else S.GREY)
    ax.text(x, y + CARD_H - 4.2, description, fontsize=5.0,
            ha="left", va="center", color=S.GREY)

    bar_x = x + 5.0
    bar_w = 16.2
    row_y = y + CARD_H - 7.0
    for i, gene in enumerate(S.GENE_ORDER):
        test_n, total_n = counts[gene]
        yy = row_y - i * 1.75
        ax.text(x, yy + 0.05, gene[0], fontsize=5.0, fontweight="bold",
                color=S.INK, ha="left", va="center")
        ax.add_patch(Rectangle((bar_x, yy - 0.45), bar_w, 0.9,
                               facecolor=S.LIGHT_GREY, edgecolor="none"))
        fraction = test_n / total_n if total_n else 0
        if fraction > 0:
            ax.add_patch(Rectangle((bar_x, yy - 0.45), bar_w * fraction, 0.9,
                                   facecolor=S.HELDOUT, edgecolor="none"))
        ax.text(bar_x + bar_w + 1.0, yy + 0.05, f"{test_n}/{total_n}",
                fontsize=5.0, color=S.GREY, ha="left", va="center")


def main() -> None:
    S.apply_rcparams()
    data = S.load_bench(exclude_wt=True, corrected=True)

    fig, ax = S.panel(W_MM, H_MM)
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
    ax.set_xlim(0, W_MM)
    ax.set_ylim(0, H_MM)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.add_patch(Rectangle((0, 0), W_MM, H_MM, facecolor="white",
                           edgecolor="none", zorder=-10))

    # Shared key.
    ax.add_patch(Rectangle((0.9, 34.7), 2.5, 1.25,
                           facecolor=S.LIGHT_GREY, edgecolor="none"))
    ax.text(4.1, 35.3, "training", fontsize=5.0, color=S.GREY,
            ha="left", va="center")
    ax.add_patch(Rectangle((17.2, 34.7), 2.5, 1.25,
                           facecolor=S.HELDOUT, edgecolor="none"))
    ax.text(20.4, 35.3, "held-out", fontsize=5.0, color=S.GREY,
            ha="left", va="center")
    ax.add_patch(Rectangle((35.5, 34.55), 3.2, 1.55,
                           facecolor="none", edgecolor=S.GREY, lw=0.5,
                           linestyle=(0, (2.0, 1.4))))
    ax.text(39.6, 35.3, "defined only / not scored", fontsize=5.0, color=S.GREY,
            ha="left", va="center")
    ax.text(99.0, 35.3, "bars: realized test n / total n",
            fontsize=5.0, color=S.GREY, ha="right", va="center")

    # The row labels expose the benchmark taxonomy without adding six equal
    # visual containers. The lower row is deliberately framed as stress tests.
    ax.text(0.9, 32.15, "WITHIN-GENE GENERALIZATION", fontsize=5.0,
            fontweight="bold", color=S.GREY, ha="left", va="center")
    ax.plot([25.0, 99.0], [32.15, 32.15], color=S.LIGHT_GREY, lw=0.45)
    ax.text(0.9, 15.95, "TRANSFER / DATA-REGIME STRESS TESTS", fontsize=5.0,
            fontweight="bold", color=S.GREY, ha="left", va="center")
    ax.plot([33.8, 99.0], [15.95, 15.95], color=S.LIGHT_GREY, lw=0.45)

    for index, (title, description, column, evaluated) in enumerate(SPLITS):
        row, col = divmod(index, 3)
        x = 0.9 + col * 33.0
        y = 16.9 - row * 16.1
        draw_card(
            ax,
            x,
            y,
            title,
            description,
            split_counts(data, column),
            evaluated,
        )

    S.save(fig, HERE / "fig1f_splits", exact=True,
           formats=("pdf", "png", "svg", "tiff"))

    for title, _, column, evaluated in SPLITS:
        print(title, split_counts(data, column), "evaluated" if evaluated else "defined-only")


if __name__ == "__main__":
    main()
