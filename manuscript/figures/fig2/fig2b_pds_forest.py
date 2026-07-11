"""Fig 2b: PDS forest for the in-house predictors.

Message: every in-house predictor head, plus the two reference baselines
(Gene-mean, WT-null), lands at chance. All 20 bootstrap 95% CIs overlap the
PDS = 0.50 chance line, so no in-house model resolves allele-level effects
above chance. External SOTA (D.EXTERNAL) are deliberately excluded; they are
Fig 4's message, not this panel's.

Data: results/_remote/unified/definitive_summary.csv via D.definitive().
Numbers plotted are point=PDS, whisker=[ci_lo, ci_hi]; PDS range 0.487-0.517,
all 20 rows have crosses==True.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import nm_style as S
import remote_data as D

HERE = os.path.dirname(os.path.abspath(__file__))


def main() -> None:
    S.apply_rcparams()

    df = D.definitive()
    keep = df[df.method.apply(D.is_inhouse) | df.method.isin(D.REFS)].copy()
    keep = keep.sort_values("PDS", ascending=True).reset_index(drop=True)

    # integrity checks: this panel's message must hold or we do not plot it
    assert len(keep) == 20, f"expected 20 rows, got {len(keep)}"
    assert set(D.EXTERNAL).isdisjoint(keep.method), "external SOTA leaked in"
    assert keep.crosses.all(), "some CI does not cross 0.50"
    assert 0.485 <= keep.PDS.min() and keep.PDS.max() <= 0.520, "PDS out of expected range"

    # colour: feature-space families for the 18 heads; references in grey
    feat_color = {
        "theta": S.GENE_COLORS["TP53"],       # blue
        "ESM": S.GENE_COLORS["KRAS"],         # orange
        "ESM+theta": S.GENE_COLORS["JAK1"],   # green
    }

    def color_of(m: str) -> str:
        if m in D.REFS:
            return S.GREY
        return feat_color[D.feature_space(m)]

    fig, ax = S.panel(58.0, 62.0)
    y = range(len(keep))

    # chance line first, behind everything
    ax.axvline(0.50, color=S.INK, lw=0.6, ls=(0, (3, 2)), zorder=1)

    for yi, (_, row) in zip(y, keep.iterrows()):
        c = color_of(row.method)
        is_ref = row.method in D.REFS
        # whisker
        ax.plot([row.ci_lo, row.ci_hi], [yi, yi], color=c, lw=0.9,
                solid_capstyle="round", zorder=2)
        # point
        ax.plot([row.PDS], [yi], marker="o", ms=3.0, mfc=c, mec="white",
                mew=0.4, ls="none", zorder=3,
                alpha=1.0 if not is_ref else 0.95)

    ax.set_yticks(list(y))
    ax.set_yticklabels(keep.method, fontsize=5.6)
    # bold-ish emphasis on the reference rows via colour only (grey ticks)
    for tick, m in zip(ax.get_yticklabels(), keep.method):
        if m in D.REFS:
            tick.set_color(S.GREY)
            tick.set_style("italic")

    ax.set_ylim(-0.7, len(keep) + 2.2)
    ax.set_xlim(0.45, 0.56)
    ax.set_xticks([0.46, 0.50, 0.54])
    ax.set_xlabel("PDS (perturbation direction score)")
    ax.tick_params(axis="y", length=0, pad=1.5)
    ax.tick_params(axis="x", pad=1.5)

    S.despine(ax, keep=("left", "bottom"))

    # "chance" note by the dashed line, just above the top data row
    ax.text(0.50, len(keep) - 0.15, "chance", ha="center", va="bottom",
            fontsize=5.6, color=S.INK)

    # panel message, bottom-right where no whisker reaches
    ax.text(0.558, 0.05, "all 95% CIs\noverlap 0.50", ha="right", va="bottom",
            fontsize=5.8, color=S.INK, style="italic", linespacing=1.25)

    # feature-family legend (marker-only proxies, no leader line, no black
    # edges), laid out horizontally in the reserved strip above the data
    from matplotlib.lines import Line2D

    def dot(color, label):
        return Line2D([0], [0], marker="o", ms=3.0, mfc=color, mec="white",
                      mew=0.4, ls="none", label=label)

    handles = [
        dot(feat_color["theta"], r"$\theta$"),
        dot(feat_color["ESM"], "ESM"),
        dot(feat_color["ESM+theta"], r"ESM+$\theta$"),
        dot(S.GREY, "reference"),
    ]
    leg = ax.legend(handles=handles, loc="lower center",
                    bbox_to_anchor=(0.5, 1.0), ncol=4, columnspacing=0.9,
                    handlelength=0.6, handletextpad=0.3, borderpad=0.2,
                    fontsize=5.4)
    leg.set_zorder(5)

    stem = os.path.join(HERE, "fig2b_pds_forest")
    S.save(fig, stem)


if __name__ == "__main__":
    main()
