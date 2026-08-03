"""Fig. 4g — published models under one harness with nulls kept distinct."""
from __future__ import annotations

import numpy as np

import fig4_common as S

# Static QA contract: Arial/Helvetica sans-serif; svg.fonttype="none";
# pdf.fonttype=42; exports .svg, .pdf, .png and .tiff at dpi=600.
# Final assembled figure target: width_mm = 183.
W_MM, H_MM = 94.0, 51.0


def main() -> None:
    S.apply_style()
    summary = S.read_csv("definitive_summary.csv").set_index("method")
    subspace = S.read_json("perturbnet_subspace_test.json")
    conditional = ["PerturbNet", "CellFlow", "scVIDR", "Biolord", "scGen"]
    references = ["Gene-mean", "WT-null"]
    head_names = [method for method in summary.index if S.is_inhouse(method)]
    if len(head_names) != 18:
        raise ValueError(f"Expected 18 feature-regression heads, observed {len(head_names)}")

    y_positions = {
        "PerturbNet": 8.0,
        "CellFlow": 7.0,
        "scVIDR": 6.0,
        "Biolord": 5.0,
        "scGen": 4.0,
        "18 heads": 2.3,
        "Gene-mean": 0.8,
        "WT-null": 0.0,
    }

    fig = S.figure(W_MM, H_MM)
    ax = fig.add_axes([0.17, 0.18, 0.53, 0.675])
    ax.axvline(0.50, color=S.INK, lw=0.75, ls=(0, (3, 2)), zorder=1)
    ax.text(0.498, 8.43, "analytic chance", ha="right", va="bottom", fontsize=5.0, color=S.GREY)

    for method in conditional + references:
        row = summary.loc[method]
        y = y_positions[method]
        color = S.INK if method in conditional else S.GREY
        ax.plot([row["ci_lo"], row["ci_hi"]], [y, y], color=color, lw=0.95, zorder=2)
        ax.plot(
            row["PDS"],
            y,
            marker="o",
            ms=4.0,
            mfc=color,
            mec="white",
            mew=0.5,
            ls="none",
            zorder=4,
        )

    head_values = summary.loc[head_names, "PDS"].sort_values().to_numpy(float)
    head_y = y_positions["18 heads"]
    head_offsets = S.deterministic_offsets(len(head_values), 0.18)
    ax.plot(
        [head_values.min(), head_values.max()],
        [head_y, head_y],
        color=S.MID_GREY,
        lw=1.2,
        zorder=2,
    )
    ax.scatter(
        head_values,
        head_y + head_offsets,
        s=8,
        facecolor=S.MID_GREY,
        edgecolor="white",
        linewidth=0.25,
        zorder=3,
    )
    ax.plot(
        np.median(head_values),
        head_y,
        marker="D",
        ms=3.8,
        mfc=S.SLATE,
        mec="white",
        mew=0.45,
        zorder=5,
    )

    rows = conditional + ["18 heads"] + references
    ax.set_yticks([y_positions[row] for row in rows])
    ax.set_yticklabels(rows)
    for tick, row in zip(ax.get_yticklabels(), rows):
        if row in references or row == "18 heads":
            tick.set_color(S.GREY)
            tick.set_style("italic")
    ax.tick_params(axis="y", length=0, pad=2)
    ax.tick_params(axis="x", pad=1.4)
    ax.set_xlim(0.44, 0.59)
    ax.set_ylim(-0.55, 9.05)
    ax.set_xticks([0.45, 0.50, 0.55])
    ax.set_xlabel("PDS (cosine)")
    S.despine(ax, keep=("bottom",))

    ax.text(0.441, 8.92, "Variant-conditionable", ha="left", va="bottom", fontsize=5.3, color=S.GREY, style="italic")
    ax.text(0.441, 2.90, "Feature regression", ha="left", va="bottom", fontsize=5.3, color=S.GREY, style="italic")
    ax.text(0.441, 1.45, "References", ha="left", va="bottom", fontsize=5.3, color=S.GREY, style="italic")
    ax.axhline(3.30, color=S.LIGHT_GREY, lw=0.55)
    ax.axhline(1.55, color=S.LIGHT_GREY, lw=0.55)

    # Right edge pulled in from 0.98: the 0.58 tick label is centred on the
    # axis end, so at 0.98 half of it sat 0.10 mm from the canvas edge.
    inset = fig.add_axes([0.75, 0.45, 0.205, 0.34])
    inset.axvline(0.50, color=S.GREY, lw=0.6, ls=(0, (3, 2)), zorder=1)
    projection_values = np.array(
        [
            subspace["pure_random_in_subspace_pds"],
            subspace["ridge_esm_subspace_pds"],
            subspace["gene_mean_subspace_pds"],
        ],
        dtype=float,
    )
    inset.plot(
        [projection_values.min(), projection_values.max()],
        [2.0, 2.0],
        color=S.MID_GREY,
        lw=2.0,
        solid_capstyle="butt",
        zorder=2,
    )
    inset.plot(
        projection_values,
        np.full(3, 2.0),
        marker="o",
        ms=2.6,
        color=S.MID_GREY,
        ls="none",
        zorder=3,
    )
    inset.plot(
        subspace["perturbnet_permutation_null_mean"],
        1.0,
        marker="D",
        ms=4.0,
        mfc=S.AMBER,
        mec="white",
        mew=0.45,
        ls="none",
        zorder=4,
    )
    inset.plot(
        subspace["perturbnet_subspace_pds"],
        0.0,
        marker="o",
        ms=4.0,
        mfc=S.INK,
        mec="white",
        mew=0.45,
        ls="none",
        zorder=4,
    )
    inset.set_yticks([2.0, 1.0, 0.0])
    inset.set_yticklabels(["random / projected", "permutation null", "PerturbNet"])
    inset.set_xlim(0.50, 0.58)
    inset.set_ylim(-0.55, 2.55)
    inset.set_xticks([0.50, 0.54, 0.58])
    inset.set_xlabel("PDS", labelpad=1.5)
    inset.tick_params(axis="y", length=0, pad=1.5, labelsize=5.0)
    inset.tick_params(axis="x", pad=1.0, labelsize=5.0)
    S.despine(inset, keep=("bottom",))
    inset.text(
        0.50,
        1.06,
        "50-D subspace diagnostic",
        transform=inset.transAxes,
        ha="center",
        va="bottom",
        fontsize=5.2,
        color=S.INK,
    )
    inset.text(
        0.99,
        0.03,
        f"permuted targets: P = {subspace['perturbnet_permutation_p']:.3f}",
        transform=inset.transAxes,
        ha="right",
        va="bottom",
        fontsize=5.0,
        color=S.GREY,
    )

    fig.text(
        0.17,
        0.89,
        "main forest: multi-seed PDS · whisker: 95% bootstrap CI",
        ha="left",
        va="bottom",
        fontsize=5.1,
        color=S.GREY,
    )
    # Right-aligned, like the note below it. Left-aligned from 0.76 the italic
    # slant carried the last glyphs past the canvas edge, and the reported text
    # extent (93.79 mm of 94.0) did not show it because matplotlib's bounding
    # box does not include the slant overhang.
    fig.text(
        0.98,
        0.27,
        "No evaluated model class\nexceeds its applicable null",
        ha="right",
        va="center",
        fontsize=5.8,
        color=S.INK,
        style="italic",
        linespacing=1.2,
    )
    # Right-aligned at 0.99: from 0.76 the longer line needed 24.4 mm and ran
    # 0.38 mm past the 100 mm canvas.
    fig.text(
        0.99,
        0.17,
        "The subspace diagnostic is a\nseparate single-run control.",
        ha="right",
        va="center",
        fontsize=5.0,
        color=S.GREY,
        linespacing=1.15,
    )
    S.save(fig, "fig4g_model_harness")


if __name__ == "__main__":
    main()
