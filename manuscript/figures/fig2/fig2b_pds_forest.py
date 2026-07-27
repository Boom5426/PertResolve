"""Fig 2b: PDS forest for the in-house predictors (left half of the b/c pair).

Message: every in-house predictor head, plus the two reference baselines
(Gene-mean, WT-null), lands at chance. All 20 bootstrap 95% CIs overlap the
PDS = 0.50 chance line, so no in-house model resolves allele-level effects
above chance. External SOTA (D.EXTERNAL) are deliberately excluded; they are
Fig 4's message, not this panel's.

Data: results/canonical/definitive_summary.csv via D.definitive().
Numbers plotted are point=PDS, whisker=[ci_lo, ci_hi]; PDS range 0.487-0.517,
all 20 rows have crosses==True.

Layout contract (Nature Methods pass):
  * This panel owns the SHARED method-label column and the SHARED feature-space
    key for the b/c pair. Fig 2c reuses ``row_order()``, draws no labels and no
    key, so the 20 method names and the key appear exactly once in the figure.
  * Canvas 58 x 63 mm with an explicit axes rectangle; Fig 2c uses the identical
    vertical rectangle (bottom 8 mm, height 47.5 mm) so the two panels share a
    common top edge and baseline when placed side by side at scale 1.0.
  * Interval styling (LW / MS / CAP) is duplicated verbatim in Fig 2c so the
    same kind of interval carries the same visual weight in both panels.
"""
import os
import sys

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import nm_style as S
import remote_data as D

HERE = os.path.dirname(os.path.abspath(__file__))

# ---- canvas geometry shared with Fig 2c -----------------------------------
W_MM, H_MM = 58.0, 63.0
AX_BOTTOM_MM, AX_HEIGHT_MM = 8.0, 47.5   # identical in fig2c -> rows line up
AX_LEFT_MM, AX_WIDTH_MM = 17.5, 39.0     # 17.5 mm = the shared label column

# ---- interval styling shared with Fig 2c ----------------------------------
LW = 0.9      # whisker line width
MS = 2.8      # point marker size
CAP = 0.22    # whisker end-cap half height, in row units

FEAT_TEX = {"theta": r"$\theta$", "esm": "ESM", "esm+theta": r"ESM+$\theta$"}


def short_label(method: str) -> str:
    """Label economy: 'KNN-esm+theta' -> 'KNN ESM+theta' (head + feature space)."""
    if method in D.REFS:
        return method
    head, _, feat = method.partition("-")
    return f"{head} {FEAT_TEX[feat]}"


def row_order() -> list[str]:
    """Canonical bottom-to-top row order for the b/c pair: ascending PDS.

    Fig 2c imports this so both panels share one method-label column. The sort
    is stable, so the order is reproducible under the several PDS ties.
    """
    df = D.definitive()
    keep = df[df.method.apply(D.is_inhouse) | df.method.isin(D.REFS)]
    keep = keep.sort_values("PDS", ascending=True, kind="stable")
    return list(keep.method)


def color_of(method: str) -> str:
    """Feature-space slate ramp for the 18 heads; grey for the two references."""
    if method in D.REFS:
        return S.FEATURE_COLORS["reference"]
    return S.FEATURE_COLORS[D.feature_space(method)]


def interval(ax, lo, hi, val, y, c):
    """Point + capped 95% interval. Identical convention in Fig 2c."""
    ax.plot([lo, hi], [y, y], color=c, lw=LW, solid_capstyle="butt", zorder=2)
    for xx in (lo, hi):
        ax.plot([xx, xx], [y - CAP, y + CAP], color=c, lw=LW,
                solid_capstyle="butt", zorder=2)
    ax.plot([val], [y], marker="o", ms=MS, mfc=c, mec="white", mew=0.4,
            ls="none", zorder=3)


def canvas():
    """Fixed-size canvas whose tight bbox equals the declared W_MM x H_MM."""
    fig = plt.figure(figsize=(W_MM * S.MM, H_MM * S.MM))
    bg = fig.add_axes([0, 0, 1, 1])
    bg.set_axis_off()
    bg.set_xlim(0, 1)
    bg.set_ylim(0, 1)
    bg.add_patch(Rectangle((0, 0), 1, 1, facecolor="white", edgecolor="none",
                           zorder=-10))
    ax = fig.add_axes([AX_LEFT_MM / W_MM, AX_BOTTOM_MM / H_MM,
                       AX_WIDTH_MM / W_MM, AX_HEIGHT_MM / H_MM])
    return fig, ax


def main() -> None:
    S.apply_rcparams()
    # Keep Greek/maths glyphs in the same Arial-metric sans as the body text.
    # nm_style sets the text font but not mathtext, whose default (DejaVu Sans)
    # would embed a second typeface for every $\theta$, $\delta$ and subscript.
    plt.rcParams.update({"mathtext.fontset": "custom",
                         "mathtext.rm": "Liberation Sans",
                         "mathtext.it": "Liberation Sans:italic",
                         "mathtext.bf": "Liberation Sans:bold",
                         "mathtext.default": "it"})

    df = D.definitive().set_index("method")
    order = row_order()

    # integrity checks: this panel's message must hold or we do not plot it
    assert len(order) == 20, f"expected 20 rows, got {len(order)}"
    assert set(D.EXTERNAL).isdisjoint(order), "external SOTA leaked in"
    assert df.loc[order, "crosses"].all(), "some CI does not cross 0.50"
    assert 0.485 <= df.loc[order, "PDS"].min(), "PDS below expected range"
    assert df.loc[order, "PDS"].max() <= 0.520, "PDS above expected range"

    fig, ax = canvas()

    ax.axvline(0.50, color=S.INK, lw=0.6, ls=(0, (3, 2)), zorder=1)

    for yi, m in enumerate(order):
        r = df.loc[m]
        interval(ax, r.ci_lo, r.ci_hi, r.PDS, yi, color_of(m))

    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([short_label(m) for m in order], fontsize=5.8)
    # reference rows: full-contrast ink, italic (the grey marker already carries
    # the "reference" cue; a grey label was simply hard to read)
    for tick, m in zip(ax.get_yticklabels(), order):
        tick.set_color(S.INK)
        if m in D.REFS:
            tick.set_style("italic")

    ax.set_ylim(-0.8, len(order) + 0.6)
    ax.set_xlim(0.455, 0.550)
    ax.set_xticks([0.46, 0.50, 0.54])
    ax.set_xlabel(r"PDS$_{cos}$ (allele discrimination)", labelpad=1.5)
    ax.tick_params(axis="y", length=0, pad=1.5)
    ax.tick_params(axis="x", pad=1.5)

    S.despine(ax, keep=("left", "bottom"))
    ax.spines["left"].set_bounds(-0.8, len(order) - 1 + 0.3)

    # "chance" note by the dashed line, in the headroom above the top data row
    ax.text(0.50, len(order) - 0.6, "chance", ha="center", va="bottom",
            fontsize=5.5, color=S.INK)

    # panel message, bottom-right where no whisker reaches
    ax.text(0.548, 0.0, "all 95% CIs\noverlap 0.50", ha="right", va="center",
            fontsize=5.5, color=S.INK, style="italic", linespacing=1.25)

    # ---- the figure's ONE feature-space key (Fig 2c and 2d draw none) ------
    def dot(color, label):
        return Line2D([0], [0], marker="o", ms=MS, mfc=color, mec="white",
                      mew=0.4, ls="none", label=label)

    handles = [dot(S.FEATURE_COLORS["theta"], r"$\theta$"),
               dot(S.FEATURE_COLORS["ESM"], "ESM"),
               dot(S.FEATURE_COLORS["ESM+theta"], r"ESM+$\theta$"),
               dot(S.FEATURE_COLORS["reference"], "reference")]
    leg = fig.legend(handles=handles, loc="upper left",
                     bbox_to_anchor=(0.012, 0.998), ncol=4, columnspacing=0.7,
                     handlelength=0.5, handletextpad=0.25, borderpad=0.0,
                     borderaxespad=0.0, fontsize=5.5, frameon=False)
    leg.set_zorder(5)

    S.save(fig, os.path.join(HERE, "fig2b_pds_forest"))

    for m in reversed(order):
        r = df.loc[m]
        print(f"{m:18} PDS={r.PDS:.3f}  CI[{r.ci_lo:.3f},{r.ci_hi:.3f}]")


if __name__ == "__main__":
    main()
