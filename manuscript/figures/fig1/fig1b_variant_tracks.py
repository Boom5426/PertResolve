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

import matplotlib.patheffects as pe
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
    # JAK1 kinase modules are labelled by their short names (JH2 pseudokinase,
    # JH1 kinase) so the in-band labels stay inside their domains at 5.5 pt.
    "JAK1": [(34, 420, "FERM"), (439, 544, "SH2-like"),
             (583, 855, "JH2"), (875, 1153, "JH1")],
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


BASE = 0.75        # lollipop stems start here, clear of the backbone band
UNIT = 0.55        # y units added per variant sharing a residue
YBOT = -3.2        # bottom of every track (room for the sub-backbone domain labels)
MM_PER_UNIT = 1.0  # millimetres of track height per y unit, SHARED by all tracks
TRACK_W = 52.0     # drawing width in mm = final placement width (scale 1.0)
PAD_TOP_MM = 2.9   # room above the axes for the gene title
PAD_BOT_MM = 2.6   # room below the axes for the residue tick labels
XLABEL_MM = 3.0    # extra room on the bottom track for the "Residue" axis label


def assign_tiers(positions, length: float):
    """Two-tier stagger for hotspot residue labels; returns {pos: tier}."""
    last_x = [-1e9, -1e9]
    min_gap = length * 0.14
    out = {}
    for pos in sorted(positions):
        tier = 0 if (pos - last_x[0]) > min_gap else 1
        last_x[tier] = pos
        out[pos] = tier
    return out


def draw_gene(df: pd.DataFrame, gene: str, show_xlabel: bool = False) -> None:
    """Draw one gene track.

    Track height is derived from the data: every track uses the same
    ``MM_PER_UNIT`` scale, so one variant at one residue is the same number of
    millimetres of stem in all four tracks (stem height is a quantity, its
    scale must not change between panels), but a track only reserves the
    vertical band it actually needs, so tracks without hotspot callouts carry
    no empty headroom. ``show_xlabel`` prints the shared "Residue" axis label;
    only the bottom track of the stack carries it.
    """
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

    # ---- pre-pass: stem tops, callout tiers and the exact vertical extent ----
    counts = sub.groupby("pos").size()
    tops = {p: BASE + c * UNIT for p, c in counts.items()}
    ymax = max([1.0] + list(tops.values()))
    hs = HOTSPOT_CODONS[gene]
    callouts = [p for p in counts.index if p in hs and p in LABEL_HINT[gene]]
    tiers_of = assign_tiers(callouts, length)
    n_tiers = (max(tiers_of.values()) + 1) if tiers_of else 0
    tier_y = [ymax + 0.7, ymax + 2.0]
    ytop = ymax + (0.5 if n_tiers == 0 else 2.0 if n_tiers == 1 else 3.3)
    if n_splice:
        ytop += 1.9   # headroom for the "no protein coordinate" note

    pad_bot = PAD_BOT_MM + (XLABEL_MM if show_xlabel else 0.0)
    axes_h = (ytop - YBOT) * MM_PER_UNIT
    fig_h = axes_h + PAD_TOP_MM + pad_bot
    fig, ax = S.panel(TRACK_W, fig_h)
    fig.subplots_adjust(left=0.015, right=0.995,
                        bottom=pad_bot / fig_h, top=1.0 - PAD_TOP_MM / fig_h)

    # --- protein backbone + domains (baseline band around y=0) ---
    # The band is deliberately taller than the lollipop base offset so an
    # in-band domain label can never be overrun by a variant marker.
    bb_h = 1.0
    ax.barh(0, length, height=bb_h, color=S.LIGHT_GREY, edgecolor="none", zorder=1)
    submotifs = {"SW I", "SW II"}  # nested KRAS switch regions: darker band, label below
    for start, end, lab in DOMAINS[gene]:
        if lab in submotifs:
            ax.barh(0, end - start, left=start, height=bb_h, color=color,
                    alpha=0.9, edgecolor="white", linewidth=0.4, zorder=3)
            ax.annotate(lab, ((start + end) / 2, -bb_h / 2), xytext=(0, -1.5),
                        textcoords="offset points", ha="center", va="top",
                        fontsize=5.0, color=S.GREY, zorder=4)
        else:
            ax.barh(0, end - start, left=start, height=bb_h, color=color,
                    alpha=0.5, edgecolor="none", zorder=2)
            # a white label only stays legible when the domain band is wide
            # enough to contain it; narrow domains get a coloured label below
            # the backbone instead of white text spilling onto light grey.
            if (end - start) > length * 0.15:
                # ink, not white: the bands are gene colour at alpha 0.5, so
                # white text loses contrast on the lighter hues (KRAS orange).
                # A white halo separates the glyphs from the lollipop stems that
                # cross the band, which otherwise read as striking the label out.
                ax.text((start + end) / 2, 0, lab, ha="center", va="center",
                        fontsize=5.5, color=S.INK, zorder=6, fontweight="bold",
                        path_effects=[pe.withStroke(linewidth=1.6,
                                                    foreground="white")])
            else:
                ax.annotate(lab, ((start + end) / 2, -bb_h / 2), xytext=(0, -1.5),
                            textcoords="offset points", ha="center", va="top",
                            fontsize=5.5, color=color, zorder=4)

    # --- lollipops: stem height = number of variants at that residue ---
    grouped = sub.groupby("pos")
    for pos, g in grouped:
        top = tops[pos]
        is_hot = pos in hs
        # class of the tallest / representative variant at this residue
        cls = g.cls.mode().iloc[0]
        fkey, marker = CLASS_STYLE.get(cls, ("gene", "o"))
        face = color if fkey == "gene" else fkey
        if is_hot:
            face = S.HOTSPOT
        ax.plot([pos, pos], [BASE, top], color=S.GREY, lw=0.5, zorder=3)
        ax.scatter([pos], [top], s=9 if is_hot else 6, marker=marker,
                   facecolor=face, edgecolor=S.INK if is_hot else "none",
                   linewidths=0.4, zorder=5)

    # hotspot residue labels lifted into a clean band above the lollipops;
    # nearby labels are staggered onto two tiers so text never overlaps
    for pos in sorted(callouts):
        res = sub.loc[sub.pos == pos, "variant"].iloc[0].split(",")[0]
        ly = tier_y[tiers_of[pos]]
        ax.plot([pos, pos], [tops[pos] + 0.08, ly - 0.12], color=S.GREY, lw=0.4,
                zorder=4)
        ax.annotate(res, (pos, ly), ha="center", va="bottom", fontsize=5.5,
                    color=S.INK, zorder=6)

    # --- axes cosmetics ---
    ax.set_xlim(-length * 0.02, length * 1.02)
    ax.set_ylim(YBOT, ytop)
    ax.set_yticks([])
    ax.set_xticks([1, length] if length < 500 else [1, length // 2, length])
    if show_xlabel:
        ax.set_xlabel("Residue", labelpad=1)
    S.despine(ax, keep=("bottom",))
    ax.tick_params(length=2.0)

    # gene label + count (left) ; splice note if any
    n_kept = len(sub)
    tag = f"{gene}  ({length:,} aa, n = {n_total} variants)"
    ax.text(0.0, 1.01, tag, transform=ax.transAxes, ha="left", va="bottom",
            fontsize=6.5, fontweight="bold", color=color)
    if n_splice:
        # inside the axes, top right: at 52 mm the title line has no room for it
        ax.text(length * 1.02, ytop - 0.15, f"+{n_splice} splice (no coord.)",
                ha="right", va="top", fontsize=5.0, color=S.GREY)

    S.save(fig, f"fig1b_track_{gene}")
    print(f"{gene}: total={n_total} kept={n_kept} splice_excluded={n_splice} "
          f"hotspot_residues={sum(1 for p in grouped.groups if p in hs)}")


def make_legend() -> None:
    """Standalone variant-class key.

    Shape only. The colour = gene mapping is already carried by the direct,
    gene-coloured track titles above (and repeated by the gene-coloured tick
    labels of panel c), so a colour key here would be redundant; NM legend
    economy asks for one shared key per figure, not one per panel.
    """
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    fig, ax = S.panel(46, 8)
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
                    columnspacing=1.2, labelspacing=0.45, borderpad=0.0,
                    title="shape = variant class")
    leg.get_title().set_fontsize(5.5)
    leg.get_title().set_color(S.GREY)
    S.save(fig, "fig1b_legend")
    plt.close(fig)


def main() -> None:
    S.apply_rcparams()
    df = S.load_bench(exclude_wt=True)  # count real variants only (470, not 472)
    for gene in S.GENE_ORDER:
        draw_gene(df, gene, show_xlabel=(gene == S.GENE_ORDER[-1]))
    make_legend()


if __name__ == "__main__":
    main()
