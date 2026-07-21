"""Figure 1b - Benchmark coverage: variant positions along four proteins.

Data-direct panel, one PDF per gene (tracks are meant to stack as panel b, so
they are kept as separate files for flexible NM layout). Source (committed):
  data/allele_perturb_bench.csv          -> variant list, per-gene counts
  data/hotspot_external_definition.txt   -> external hotspot codons (de-leaked)

Each track is a lollipop / needle diagram: a protein backbone with annotated
domains, and one stem per residue position (stem height = number of distinct
variants at that residue). Hotspot/pathogenic residues are highlighted.
Variant classes: missense (gene colour), nonsense * (dark), synonymous (grey).
Splice variants have no protein coordinate and are excluded from the track
(reported in stdout).

Domain boundaries are documented constants: DNA-binding / G-domain / ZF1 / ZF2 /
JH2 / JH1 from Supplementary Table 1; TAD, tetramerization, FERM, SH2-like and
KRAS switch regions from UniProt (P04637, P01116, P15976, P23458).

Run:  python fig1b_variant_tracks.py
Out:  fig1b_track_{TP53,KRAS,GATA1,JAK1}.pdf (+ .png previews)
"""
from __future__ import annotations

import os
import re
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import nm_style as S

# ---- documented domain annotations (start, end, label) --------------------
DOMAINS = {
    "TP53": [(1, 42, "TAD"), (94, 293, "DNA-binding"), (323, 356, "Tetramer.")],
    "KRAS": [(1, 169, "G-domain"), (30, 38, "SW I"), (60, 76, "SW II"),
             (167, 189, "HVR")],
    "GATA1": [(204, 228, "ZF1"), (258, 282, "ZF2")],
    "JAK1": [(34, 420, "FERM"), (439, 544, "SH2-like"),
             (583, 855, "JH2 pseudokinase"), (875, 1153, "JH1 kinase")],
}

# ---- external-only hotspot codons (verbatim from hotspot_external_definition.txt)
HOTSPOT_CODONS = {
    "TP53": {175, 176, 179, 220, 238, 245, 248, 249, 273, 282, 285},
    "KRAS": {12, 13, 59, 61, 117, 146},
    "GATA1": {204, 207, 225, 228, 258, 261, 279, 282},
    "JAK1": set(),  # JAK1 hotspot = JH2/JH1 module membership; shown via domains
}

# residues to label explicitly (canonical, if present in the data)
LABEL_HINT = {
    "TP53": {175, 248, 273},
    "KRAS": {12, 13, 61},
    "GATA1": {205, 218},
    "JAK1": set(),
}


def parse(variant: str):
    """Return (position, class) for a variant name; position None if no protein coord."""
    v = str(variant).strip()
    if v.startswith("splice"):
        return None, "splice"
    if v.startswith("syn"):
        m = re.search(r"[A-Za-z](\d+)", v)
        return (int(m.group(1)) if m else None), "synonymous"
    first = v.split(",")[0]  # multi-substitution names -> use the first
    m = re.match(r"^[A-Z](\d+)([A-Z*])", first)
    if not m:
        return None, "other"
    pos = int(m.group(1))
    cls = "nonsense" if m.group(2) == "*" else "missense"
    return pos, cls


CLASS_STYLE = {  # (facecolor-key, marker)
    "missense": ("gene", "o"),
    "nonsense": (S.INK, "s"),
    "synonymous": (S.GREY, "D"),
}


def draw_gene(df: pd.DataFrame, gene: str) -> None:
    length = S.PROTEIN_LEN[gene]
    color = S.GENE_COLORS[gene]
    sub = df[df.gene == gene].copy()
    parsed = sub.variant.map(parse)
    sub["pos"] = [p for p, _ in parsed]
    sub["cls"] = [c for _, c in parsed]
    n_total = len(sub)
    n_splice = int((sub.cls == "splice").sum())
    sub = sub.dropna(subset=["pos"])
    sub["pos"] = sub["pos"].astype(int)

    fig, ax = S.panel(89, 20)

    # --- protein backbone + domains (baseline band around y=0) ---
    bb_h = 0.6
    ax.barh(0, length, height=bb_h, color=S.LIGHT_GREY, edgecolor="none", zorder=1)
    submotifs = {"SW I", "SW II"}  # nested KRAS switch regions: darker band, label below
    for start, end, lab in DOMAINS[gene]:
        if lab in submotifs:
            ax.barh(0, end - start, left=start, height=bb_h, color=color,
                    alpha=0.9, edgecolor="white", linewidth=0.4, zorder=3)
            ax.annotate(lab, ((start + end) / 2, -bb_h / 2), xytext=(0, -1.5),
                        textcoords="offset points", ha="center", va="top",
                        fontsize=4.5, color=S.GREY, zorder=4)
        else:
            ax.barh(0, end - start, left=start, height=bb_h, color=color,
                    alpha=0.5, edgecolor="none", zorder=2)
            ax.text((start + end) / 2, 0, lab, ha="center", va="center",
                    fontsize=5, color="white", zorder=4,
                    fontweight="bold" if (end - start) > length * 0.12 else "normal")

    # --- lollipops: stem height = number of variants at that residue ---
    unit = 0.55
    hs = HOTSPOT_CODONS[gene]
    grouped = sub.groupby("pos")
    ymax = 1.0
    to_label = []  # (pos, top, residue) drawn in a band above all stems
    for pos, g in grouped:
        count = len(g)
        top = bb_h / 2 + count * unit
        ymax = max(ymax, top)
        is_hot = pos in hs
        # class of the tallest / representative variant at this residue
        cls = g.cls.mode().iloc[0]
        fkey, marker = CLASS_STYLE.get(cls, ("gene", "o"))
        face = color if fkey == "gene" else fkey
        if is_hot:
            face = S.HOTSPOT
        ax.plot([pos, pos], [bb_h / 2, top], color=S.GREY, lw=0.5, zorder=3)
        ax.scatter([pos], [top], s=9 if is_hot else 6, marker=marker,
                   facecolor=face, edgecolor=S.INK if is_hot else "none",
                   linewidths=0.4, zorder=5)
        if is_hot and pos in LABEL_HINT[gene]:
            to_label.append((pos, top, g.iloc[0].variant.split(",")[0]))

    # hotspot residue labels lifted into a clean band above the lollipops;
    # nearby labels are staggered onto two tiers so text never overlaps
    tiers = [ymax + 0.9, ymax + 1.75]
    last_x = [-1e9, -1e9]
    min_gap = length * 0.11
    for pos, top, res in sorted(to_label):
        tier = 0 if (pos - last_x[0]) > min_gap else 1
        last_x[tier] = pos
        ly = tiers[tier]
        ax.plot([pos, pos], [top + 0.08, ly - 0.12], color=S.GREY, lw=0.4,
                zorder=4)
        ax.annotate(res, (pos, ly), ha="center", va="bottom", fontsize=5,
                    color=S.INK, zorder=6)

    # --- axes cosmetics ---
    ax.set_xlim(-length * 0.02, length * 1.02)
    ax.set_ylim(-1.2, ymax + 2.8)
    ax.set_yticks([])
    ax.set_xticks([1, length] if length < 500 else [1, length // 2, length])
    ax.set_xlabel("Residue", labelpad=1)
    S.despine(ax, keep=("bottom",))
    ax.tick_params(length=2.0)

    # gene label + count (left) ; splice note if any
    n_kept = len(sub)
    tag = f"{gene}  ({length:,} aa, n = {n_total} variants)"
    ax.text(0.0, 1.02, tag, transform=ax.transAxes, ha="left", va="bottom",
            fontsize=6.5, fontweight="bold", color=color)
    if n_splice:
        ax.text(1.0, 1.02, f"+{n_splice} splice (no coord.)", transform=ax.transAxes,
                ha="right", va="bottom", fontsize=5, color=S.GREY)

    S.save(fig, f"fig1b_track_{gene}")
    print(f"{gene}: total={n_total} kept={n_kept} splice_excluded={n_splice} "
          f"hotspot_residues={sum(1 for p in grouped.groups if p in hs)}")


def make_legend() -> None:
    """Standalone variant-class + hotspot legend (shape = class, colour = gene)."""
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    fig, ax = S.panel(52, 12)
    ax.axis("off")
    handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=S.GREY,
               markersize=4, label="Missense"),
        Line2D([0], [0], marker="s", color="none", markerfacecolor=S.INK,
               markersize=4, label="Nonsense"),
        Line2D([0], [0], marker="D", color="none", markerfacecolor=S.GREY,
               markersize=3.5, label="Synonymous"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=S.HOTSPOT,
               markeredgecolor=S.INK, markeredgewidth=0.4, markersize=4.5,
               label="Hotspot / pathogenic"),
    ]
    leg = ax.legend(handles=handles, ncol=2, loc="center", handletextpad=0.3,
                    columnspacing=1.2, labelspacing=0.5, borderpad=0.0)
    ax.text(0.5, 1.05, "shape = variant class   ·   colour = gene",
            transform=ax.transAxes, ha="center", va="bottom", fontsize=5,
            color=S.GREY)
    S.save(fig, "fig1b_legend")
    plt.close(fig)


def main() -> None:
    S.apply_rcparams()
    df = S.load_bench(exclude_wt=True)  # count real variants only (470, not 472)
    for gene in S.GENE_ORDER:
        draw_gene(df, gene)
    make_legend()


if __name__ == "__main__":
    main()
