"""Figure 1d (right) — theta features in a common coordinate system.

The PCA reads allele_perturb_bench_v2.csv explicitly through
``load_bench(corrected=True)``.  Gene is redundantly encoded by colour, marker
shape and direct labels, preserving meaning in greyscale and for readers with
colour-vision deficiencies.

Run: python fig1d_theta_pca.py
Output: fig1d_theta_pca.svg/.pdf/.png

Imported export contract: Arial; svg.fonttype: "none"; pdf.fonttype: 42;
outputs .svg, .pdf, .png and .tiff at dpi=600. Final assembled figure target:
width_mm = 183.
"""
from __future__ import annotations

import os
import sys

import matplotlib.patheffects as pe
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import nm_style as S

from pathlib import Path

HERE = Path(__file__).resolve().parent

W_MM, H_MM = 47.0, 42.0
THETA_COLS = ["d_hydro", "d_vol", "d_charge", "fold_core", "cat_switch", "is_hotspot"]
MARKERS = {"TP53": "o", "KRAS": "s", "GATA1": "^", "JAK1": "D"}
LABEL_OFFSETS = {
    "TP53": (-16, -12),
    "KRAS": (9, -8),
    "GATA1": (-18, 12),
    "JAK1": (10, 11),
}


def standardized_pca2(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Two-component PCA of z-scored columns; scores and percent variance.

    This reproduces sklearn's ``StandardScaler`` + ``PCA(n_components=2)``
    exactly and is written in numpy only because scikit-learn cannot be
    installed into this PEP 668 system interpreter; nothing about the result
    differs. StandardScaler uses the population standard deviation (ddof=0),
    which is numpy's default, and sklearn's PCA resolves the sign ambiguity of
    the SVD with ``svd_flip``, replicated here so the scatter keeps a fixed
    orientation across runs rather than mirroring at random.
    """
    z = (matrix - matrix.mean(axis=0)) / matrix.std(axis=0, ddof=0)
    u, s, _ = np.linalg.svd(z - z.mean(axis=0), full_matrices=False)
    # svd_flip, u_based_decision=True: force the largest-magnitude entry of each
    # left singular vector positive.
    signs = np.sign(u[np.abs(u).argmax(axis=0), range(u.shape[1])])
    u, s = u * signs, s
    explained = s**2 / (s**2).sum() * 100.0
    return (u * s)[:, :2], explained[:2]


def main() -> None:
    S.apply_rcparams()
    data = S.load_bench(exclude_wt=True, corrected=True)
    matrix = data[THETA_COLS].to_numpy(float)
    if not np.isfinite(matrix).all():
        raise ValueError("Corrected theta matrix contains non-finite values")

    scores, explained = standardized_pca2(matrix)

    fig, ax = S.panel(W_MM, H_MM)
    fig.subplots_adjust(left=0.21, right=0.98, bottom=0.19, top=0.88)
    draw_order = sorted(S.GENE_ORDER, key=lambda gene: -(data.gene == gene).sum())
    for gene in draw_order:
        mask = (data.gene == gene).to_numpy()
        ax.scatter(
            scores[mask, 0],
            scores[mask, 1],
            s=5.0,
            marker=MARKERS[gene],
            facecolor=S.GENE_COLORS[gene],
            edgecolor="white",
            alpha=0.68,
            linewidths=0.16,
            zorder=2,
        )

    for gene in S.GENE_ORDER:
        mask = (data.gene == gene).to_numpy()
        centroid = scores[mask].mean(axis=0)
        ax.annotate(
            gene,
            xy=centroid,
            xytext=LABEL_OFFSETS[gene],
            textcoords="offset points",
            fontsize=5.3,
            fontweight="bold",
            color=S.GENE_COLORS[gene],
            ha="center",
            va="center",
            path_effects=[pe.withStroke(linewidth=1.5, foreground="white")],
            zorder=5,
        )

    ax.axhline(0, color=S.LIGHT_GREY, lw=0.4, zorder=0)
    ax.axvline(0, color=S.LIGHT_GREY, lw=0.4, zorder=0)
    ax.set_xlabel(f"PC1 ({explained[0]:.0f}%)", labelpad=1)
    ax.set_ylabel(f"PC2 ({explained[1]:.0f}%)", labelpad=1)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.text(
        0.0,
        1.08,
        "Shared molecular-feature space",
        transform=ax.transAxes,
        fontsize=5.9,
        fontweight="bold",
        ha="left",
    )
    ax.text(
        0.0,
        1.015,
        "470 variants · θ features",
        transform=ax.transAxes,
        fontsize=5.0,
        color=S.GREY,
        ha="left",
    )
    S.despine(ax, keep=("left", "bottom"))
    S.save(fig, HERE / "fig1d_theta_pca", exact=True,
           formats=("pdf", "png", "svg", "tiff"))

    print("explained variance %:", np.round(explained, 1), "n:", len(data))


if __name__ == "__main__":
    main()
