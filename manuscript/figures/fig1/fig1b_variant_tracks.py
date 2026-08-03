"""Figure 1b — data-direct benchmark coverage across four proteins.

Marker shape encodes molecular consequence.  Gene colour identifies the protein,
whereas amber fill independently marks externally annotated hotspot/pathogenic
residues.  Lollipop height is the number of variant conditions mapped to a
residue.  Multi-site conditions are mapped to their first residue and disclosed
in the panel; splice conditions without a protein coordinate are not plotted.

Run: python fig1b_variant_tracks.py
Output: fig1b_coverage.svg/.pdf/.png

Imported export contract: Arial; svg.fonttype: "none"; pdf.fonttype: 42;
outputs .svg, .pdf, .png and .tiff at dpi=600. Final assembled figure target:
width_mm = 183.
"""
from __future__ import annotations

import os
import re
import sys

import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import nm_style as S

from pathlib import Path

HERE = Path(__file__).resolve().parent

W_MM, H_MM = 66.0, 72.0
DOMAINS = {
    "TP53": [(1, 42, "TAD"), (94, 293, "DNA-binding"), (323, 356, "Tetramer.")],
    "KRAS": [(1, 169, "G-domain"), (30, 38, "SW I"), (60, 76, "SW II"),
             (167, 189, "HVR")],
    "GATA1": [(204, 228, "ZF1"), (258, 282, "ZF2")],
    "JAK1": [(34, 420, "FERM"), (439, 544, "SH2-like"),
             (583, 855, "JH2"), (875, 1153, "JH1")],
}

HOTSPOT_CODONS = {
    "TP53": {175, 176, 179, 220, 238, 245, 248, 249, 273, 282, 285},
    "KRAS": {12, 13, 59, 61, 117, 146},
    "GATA1": {204, 207, 225, 228, 258, 261, 279, 282},
    "JAK1": set(),
}

LABEL_HINT = {
    "TP53": {175, 248, 273},
    # G12 and G13 are adjacent; label G12 only and leave G13 encoded by its
    # amber marker so residue text cannot collide at the final panel width.
    "KRAS": {12, 61},
    "GATA1": set(),
    "JAK1": set(),
}

CONTEXT = {
    "TP53": "Perturb-seq · A549",
    "KRAS": "Perturb-seq · A549",
    "GATA1": "base editing · HSPC",
    "JAK1": "scSNV-seq · HT-29",
}

CLASS_MARKER = {"missense": "o", "nonsense": "s", "synonymous": "D"}


def parse_variant(variant: str):
    value = str(variant).strip()
    if value.startswith("splice"):
        return None, "splice", None
    if value.startswith("syn"):
        match = re.search(r"([A-Za-z])(\d+)", value)
        return (
            int(match.group(2)) if match else None,
            "synonymous",
            f"{match.group(1).upper()}{match.group(2)}" if match else None,
        )
    first = value.split(",")[0]
    match = re.match(r"^([A-Z])(\d+)([A-Z*])", first)
    if not match:
        return None, "other", None
    consequence = "nonsense" if match.group(3) == "*" else "missense"
    return int(match.group(2)), consequence, f"{match.group(1)}{match.group(2)}"


def assign_tiers(positions, protein_length):
    last = [-1e9, -1e9]
    minimum_gap = protein_length * 0.14
    result = {}
    for pos in sorted(positions):
        tier = 0 if pos - last[0] > minimum_gap else 1
        last[tier] = pos
        result[pos] = tier
    return result


def draw_track(ax, data: pd.DataFrame, gene: str, *, show_xlabel: bool = False):
    length = S.PROTEIN_LEN[gene]
    colour = S.GENE_COLORS[gene]
    sub = data[data.gene == gene].copy()
    parsed = sub.variant.map(parse_variant)
    sub["pos"] = [item[0] for item in parsed]
    sub["class"] = [item[1] for item in parsed]
    sub["residue_label"] = [item[2] for item in parsed]
    n_total = len(sub)
    n_multi = int(sub.variant.astype(str).str.contains(",", regex=False).sum())
    n_splice = int((sub["class"] == "splice").sum())
    n_before = len(sub)
    plotted = sub.dropna(subset=["pos"]).copy()
    n_after = len(plotted)
    excluded_count = n_before - n_after
    assert excluded_count == n_splice, (
        f"{gene}: {excluded_count} coordinate-free conditions but "
        f"{n_splice} parsed splice conditions"
    )
    plotted["pos"] = plotted["pos"].astype(int)

    counts = plotted.groupby("pos").size()
    base = 0.72
    unit = 0.48
    tops = {int(pos): base + count * unit for pos, count in counts.items()}
    ymax = max([1.1] + list(tops.values()))
    callouts = [
        int(pos)
        for pos in counts.index
        if int(pos) in HOTSPOT_CODONS[gene] and int(pos) in LABEL_HINT[gene]
    ]
    tiers = assign_tiers(callouts, length)
    tier_y = [ymax + 0.40, ymax + 1.45]

    # Protein backbone and domains.
    ax.barh(0, length, height=0.92, color=S.LIGHT_GREY, edgecolor="none", zorder=1)
    for start, end, label in DOMAINS[gene]:
        nested = label in {"SW I", "SW II"}
        alpha = 0.85 if nested else 0.48
        ax.barh(
            0,
            end - start,
            left=start,
            height=0.92,
            color=colour,
            alpha=alpha,
            edgecolor="white" if nested else "none",
            linewidth=0.35,
            zorder=2,
        )
        width_fraction = (end - start) / length
        if width_fraction > 0.15:
            ax.text(
                (start + end) / 2,
                0,
                label,
                ha="center",
                va="center",
                fontsize=5.1,
                fontweight="bold",
                color=S.INK,
                path_effects=[pe.withStroke(linewidth=1.35, foreground="white")],
                zorder=6,
            )
        else:
            ax.annotate(
                label,
                ((start + end) / 2, -0.46),
                xytext=(0, -1.2),
                textcoords="offset points",
                ha="center",
                va="top",
                fontsize=5.0,
                color=S.INK if not nested else S.GREY,
                zorder=6,
            )

    # Lollipops.  Shape and hotspot colour are independent encodings.
    for pos, group in plotted.groupby("pos"):
        pos = int(pos)
        top = tops[pos]
        consequence = group["class"].mode().iloc[0]
        marker = CLASS_MARKER.get(consequence, "o")
        hotspot = pos in HOTSPOT_CODONS[gene]
        face = S.HOTSPOT if hotspot else (
            S.INK if consequence == "nonsense"
            else S.GREY if consequence == "synonymous"
            else colour
        )
        ax.plot([pos, pos], [base, top], color=S.GREY, lw=0.42, zorder=3)
        ax.scatter(
            [pos],
            [top],
            s=8.5 if hotspot else 5.4,
            marker=marker,
            facecolor=face,
            edgecolor=S.INK if hotspot else "none",
            linewidths=0.35,
            zorder=5,
        )

    # Stable residue labels (R175, G12), never arbitrary first alleles (R175H).
    for pos in sorted(callouts):
        label = plotted.loc[plotted.pos == pos, "residue_label"].dropna().iloc[0]
        label_y = tier_y[tiers[pos]]
        ax.plot([pos, pos], [tops[pos] + 0.06, label_y - 0.08],
                color=S.GREY, lw=0.35, zorder=4)
        ax.text(pos, label_y, label, ha="center", va="bottom",
                fontsize=5.0, color=S.INK, zorder=7)

    # Header: count on the left; assay/cell context on the right.
    ax.text(
        0,
        0.985,
        f"{gene} · {length:,} aa · n={n_total}",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=6.2,
        fontweight="bold",
        color=S.INK,
    )
    suffix = []
    if n_multi:
        suffix.append(f"{n_multi} multi-site")
    if n_splice:
        suffix.append(f"{n_splice} splice")
    context = CONTEXT[gene] + (f" · {', '.join(suffix)}" if suffix else "")
    ax.text(
        1.0,
        0.985,
        context,
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=5.0,
        color=S.GREY,
    )

    ax.set_xlim(-0.02 * length, 1.02 * length)
    # Reserve a dedicated header band above callouts. This prevents residue
    # labels—especially KRAS G12—from colliding with the gene title.
    ytop = ymax + (3.0 if callouts else 1.8)
    ax.set_ylim(-2.25, ytop)
    ax.set_yticks([])
    if show_xlabel:
        ax.set_xticks([1, length // 2, length])
    else:
        ax.set_xticks([])
    S.despine(ax, keep=("bottom",) if show_xlabel else ())
    ax.tick_params(axis="x", length=1.8, pad=1)


def main() -> None:
    S.apply_rcparams()
    data = S.load_bench(exclude_wt=True, corrected=True)

    fig = plt.figure(figsize=(W_MM * S.MM, H_MM * S.MM))
    # bottom, height in millimetres; top padding is held inside each slot.
    slots = {
        "TP53": (56.8, 12.5),
        "KRAS": (42.6, 11.9),
        "GATA1": (28.3, 11.9),
        "JAK1": (13.8, 12.0),
    }
    for gene in S.GENE_ORDER:
        bottom_mm, height_mm = slots[gene]
        ax = fig.add_axes([
            1.6 / W_MM,
            bottom_mm / H_MM,
            62.8 / W_MM,
            height_mm / H_MM,
        ])
        draw_track(ax, data, gene, show_xlabel=(gene == "JAK1"))

    # One compact, semantically correct shared key.
    ax_leg = fig.add_axes([1.0 / W_MM, 0.5 / H_MM, 64.0 / W_MM, 10.2 / H_MM])
    ax_leg.axis("off")
    handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=S.GREY,
               markersize=3.7, label="missense"),
        Line2D([0], [0], marker="s", color="none", markerfacecolor=S.INK,
               markersize=3.7, label="nonsense"),
        Line2D([0], [0], marker="D", color="none", markerfacecolor=S.GREY,
               markersize=3.3, label="synonymous"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=S.HOTSPOT,
               markeredgecolor=S.INK, markeredgewidth=0.35, markersize=4.0,
               label="external hotspot"),
    ]
    legend = ax_leg.legend(
        handles=handles,
        ncol=4,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.98),
        handletextpad=0.25,
        columnspacing=0.75,
        borderpad=0,
        fontsize=5.0,
        title="marker shape = consequence; amber fill = hotspot annotation",
    )
    legend.get_title().set_fontsize(5.0)
    legend.get_title().set_color(S.GREY)
    # Two lines, not one: at 5 pt this sentence measures 67.8 mm and the panel
    # canvas is 66.0 mm, so as a single line it ran off both edges and was
    # clipped in the composite. 5 pt is the Nature Methods floor, so it cannot
    # be set smaller.
    ax_leg.text(
        0.5,
        0.03,
        "Stem height = conditions per residue;\n"
        "multi-site conditions map to the first residue.",
        transform=ax_leg.transAxes,
        ha="center",
        va="bottom",
        fontsize=5.0,
        color=S.GREY,
        linespacing=1.15,
    )

    S.save(fig, HERE / "fig1b_coverage", exact=True,
           formats=("pdf", "png", "svg", "tiff"))


if __name__ == "__main__":
    main()
