"""Fig 4g (HERO): published perturbation models scored under one identical harness.

Message: no evaluated model class exceeds the empirical discrimination floor.
Every published variant-conditionable model, the in-house feature-regression
range, and both references sit at the PDS = 0.50 chance line; their 95% bootstrap
CIs overlap chance. PerturbNet's point estimate (0.542) clears the analytic 0.50
line but falls inside its own 50-dim subspace null (0.52), where random vectors
already score 0.52, so the exceedance is not significant.

Data (committed): results/_remote/unified/definitive_summary.csv via D.definitive().
Numbers plotted (point = PDS, whisker = [ci_lo, ci_hi]):
  Variant-conditionable:
    PerturbNet 0.542 [0.512, 0.571]  (subspace null 0.52; crosses=False, n.s.)
    CellFlow   0.511 [0.483, 0.537]
    scVIDR     0.504 [0.478, 0.529]
    Biolord    0.494 [0.466, 0.523]
    scGen      0.486 [0.456, 0.516]
  Feature-regression baselines (18 in-house heads): PDS range 0.487-0.517 (band)
  References:
    WT-null    0.500 [0.500, 0.500]
    Gene-mean  0.489 [0.462, 0.515]
The 0.52 subspace-null value is the committed Table 1 / Table 5 footnote figure
("random predictions confined to this subspace already score 0.52").

Run:  python fig4g_external_forest.py  ->  fig4g_external_forest.pdf (+ .png)
"""
from __future__ import annotations

import os
import sys

from matplotlib.lines import Line2D
from matplotlib.patches import Patch

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import nm_style as S
import remote_data as D

HERE = os.path.dirname(os.path.abspath(__file__))

# committed subspace-null level for PerturbNet (Table 1/Table 5 footnote)
PNET_SUBSPACE_NULL = 0.52


def main() -> None:
    S.apply_rcparams()

    df = D.definitive().set_index("method")

    # ---- integrity: this panel's message must hold or we do not plot it -----
    heads = df[[D.is_inhouse(m) for m in df.index]]
    assert len(heads) == 18, f"expected 18 in-house heads, got {len(heads)}"
    head_lo, head_hi = float(heads.PDS.min()), float(heads.PDS.max())
    assert abs(head_lo - 0.487) < 1e-6 and abs(head_hi - 0.517) < 1e-6, \
        f"in-house range drifted: {head_lo}-{head_hi}"
    # every model class except PerturbNet must have a CI that crosses 0.50
    for m in ["CellFlow", "scVIDR", "Biolord", "scGen", "Gene-mean", "WT-null"]:
        assert bool(df.loc[m, "crosses"]), f"{m} unexpectedly excludes 0.50"
    assert not bool(df.loc["PerturbNet", "crosses"]), "PerturbNet should not cross"
    # PerturbNet point must lie inside its own subspace null band's reach
    assert df.loc["PerturbNet", "PDS"] > PNET_SUBSPACE_NULL, "PerturbNet vs null"

    BLUE = S.GENE_COLORS["TP53"]   # variant-conditionable published models
    GREY = S.GREY                   # references + in-house band + null band

    # ---- rows, grouped top-to-bottom; larger y = higher on axis -------------
    # family blocks separated by a blank gutter row.
    cond = ["PerturbNet", "CellFlow", "scVIDR", "Biolord", "scGen"]
    refs = ["WT-null", "Gene-mean"]

    # assign y positions (top group highest). One gutter between families,
    # plus a dedicated band row for the feature-regression family.
    rows = []  # (label, y, kind)
    y = 0.0
    # References (bottom)
    for m in refs:
        rows.append((m, y, "ref")); y += 1.0
    y += 0.9  # gutter
    # Feature-regression band (middle) - single row
    band_y = y
    rows.append(("18 heads", y, "band")); y += 1.0
    y += 0.9  # gutter
    # Variant-conditionable (top), listed so PerturbNet ends up at the very top
    for m in reversed(cond):
        rows.append((m, y, "cond")); y += 1.0

    ymax = y - 1.0

    fig, ax = S.panel(62.0, 56.0)
    # leave room on the left for the vertical family labels
    fig.subplots_adjust(left=0.33, right=0.985, top=0.90, bottom=0.19)

    ax.set_ylim(-0.7, ymax + 0.7)
    y0, y1 = ax.get_ylim()

    # chance line behind everything
    ax.axvline(0.50, color=S.INK, lw=0.6, ls=(0, (3, 2)), zorder=1)

    # PerturbNet subspace-null band: light grey vertical strip from 0.50 to 0.52,
    # drawn ONLY across the PerturbNet row (this null applies to PerturbNet, not
    # to models scored in the full gene space).
    pnet_y = [ry for (m, ry, k) in rows if m == "PerturbNet"][0]
    bfrac_lo = (pnet_y - 0.5 - y0) / (y1 - y0)
    bfrac_hi = (pnet_y + 0.5 - y0) / (y1 - y0)
    ax.axvspan(0.50, PNET_SUBSPACE_NULL, ymin=bfrac_lo, ymax=bfrac_hi,
               color=S.LIGHT_GREY, alpha=0.7, lw=0, zorder=0)

    # ---- draw each row ------------------------------------------------------
    for (m, ry, kind) in rows:
        if kind == "band":
            # feature-regression range as a horizontal bracket/band
            ax.plot([head_lo, head_hi], [ry, ry], color=GREY, lw=3.4,
                    solid_capstyle="butt", alpha=0.45, zorder=2)
            # thin end caps
            for xc in (head_lo, head_hi):
                ax.plot([xc, xc], [ry - 0.22, ry + 0.22], color=GREY,
                        lw=0.7, zorder=3)
            continue

        row = df.loc[m]
        c = BLUE if kind == "cond" else GREY
        # whisker
        ax.plot([row.ci_lo, row.ci_hi], [ry, ry], color=c, lw=0.9,
                solid_capstyle="round", zorder=2)
        # point
        ax.plot([row.PDS], [ry], marker="o", ms=3.2, mfc=c, mec="white",
                mew=0.4, ls="none", zorder=4)

    # ---- y tick labels ------------------------------------------------------
    yticks = [ry for (_, ry, _) in rows]
    ylabels = [m for (m, _, _) in rows]
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels, fontsize=5.8)
    for tick, (m, _, kind) in zip(ax.get_yticklabels(), rows):
        if kind in ("ref", "band"):
            tick.set_color(GREY)
            tick.set_style("italic")

    ax.set_xlim(0.44, 0.595)
    ax.set_xticks([0.45, 0.50, 0.55])
    ax.set_xlabel("PDS (perturbation direction score)")
    ax.tick_params(axis="y", length=0, pad=1.5)
    ax.tick_params(axis="x", pad=1.5)
    S.despine(ax, keep=("left", "bottom"))

    # ---- family group labels (small italic), left of the axis --------------
    # placed in the figure's left margin so they never collide with tick text.
    xlab = -0.40  # axes fraction (into the reserved left margin)
    trans = ax.get_yaxis_transform()

    def group_label(text, ylo, yhi):
        yc = (ylo + yhi) / 2.0
        ax.annotate(text, xy=(xlab, yc), xycoords=trans,
                    ha="center", va="center", rotation=90,
                    fontsize=5.6, style="italic", color=S.INK,
                    annotation_clip=False)

    cond_ys = [ry for (m, ry, k) in rows if k == "cond"]
    ref_ys = [ry for (m, ry, k) in rows if k == "ref"]
    group_label("Variant-\nconditionable", min(cond_ys), max(cond_ys))
    group_label("Feature\nregression", band_y, band_y)
    group_label("Reference", min(ref_ys), max(ref_ys))

    # ---- "chance" tag at top of dashed line --------------------------------
    ax.text(0.50, ymax + 0.78, "chance", ha="center", va="bottom",
            fontsize=5.6, color=S.INK)

    # ---- PerturbNet not-significant annotation -----------------------------
    # label the subspace-null band above the PerturbNet row, clear of the whisker.
    ax.annotate("50-dim subspace null;\nnot significant",
                xy=(PNET_SUBSPACE_NULL, pnet_y + 0.48),
                xytext=(0.508, pnet_y + 1.0),
                textcoords="data", ha="left", va="bottom",
                fontsize=5.0, color=GREY, style="italic",
                arrowprops=dict(arrowstyle="-", lw=0.5, color=GREY,
                                shrinkA=1, shrinkB=1),
                annotation_clip=False)

    # ---- panel message ------------------------------------------------------
    # sits in the empty upper-left quadrant, clear of all points and the legend.
    ax.text(0.445, max(cond_ys) - 1.0, "no model class\nexceeds the floor",
            ha="left", va="center", fontsize=5.8, color=S.INK, style="italic")

    # ---- compact legend (no black edges) -----------------------------------
    handles = [
        Line2D([0], [0], marker="o", ms=3.2, mfc=BLUE, mec="white", mew=0.4,
               color=BLUE, lw=0.9, label="Published model"),
        Patch(facecolor=GREY, alpha=0.45, edgecolor="none",
              label="In-house head range"),
        Line2D([0], [0], marker="o", ms=3.2, mfc=GREY, mec="white", mew=0.4,
               color=GREY, lw=0.9, label="Reference"),
    ]
    leg = ax.legend(handles=handles, loc="lower right",
                    bbox_to_anchor=(1.005, -0.02), handlelength=1.1,
                    handletextpad=0.4, labelspacing=0.28, borderpad=0.2,
                    fontsize=5.2)
    leg.set_zorder(6)

    stem = os.path.join(HERE, "fig4g_external_forest")
    S.save(fig, stem)
    print(f"in-house head range: {head_lo:.3f}-{head_hi:.3f}")
    print("saved", stem + ".pdf/.png")


if __name__ == "__main__":
    main()
