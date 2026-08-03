"""Fig. 1a — problem-first map of nested allele-prediction requirements.

Core claim: allele-specific prediction is interpretable only when the measured
response is first detectable and the causal allele is reproducibly identifiable.
Failure before identification is measurement-limited; failure after an
identifiable ground truth is model-limited.

This is a conceptual panel: no experimental values or simulated quantitative
observations are drawn. Small cell clouds are deterministic schematic glyphs.

Static QA contract: font resolved by nm_style (Arial, or an Arial/Helvetica-metric
stand-in); svg.fonttype: "none"; pdf.fonttype: 42; outputs .svg, .pdf, .png and
.tiff at dpi=600.

Width: 113.5 mm, not the 122.5 mm this panel was first drawn at. The composite
is 183 mm and row 1 is this panel beside fig1b_coverage (66.0 mm) with a 3.5 mm
gutter, so 122.5 mm overran the page by 5.5 mm before any gutter. Both panels
bottom out at 5.0 pt type, so neither can be scaled down at placement time
without dropping below the 5 pt floor; the 9 mm therefore comes out of this
panel's layout instead. The verdict column keeps its 41.5 mm because its 5.0 pt
body text already nearly fills the card, so the whole reduction is taken from
the nested-requirement block on the left.
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Ellipse, FancyArrowPatch, FancyBboxPatch

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import nm_style as S  # noqa: E402

HERE = Path(__file__).resolve().parent
MM = 1.0 / 25.4
W_MM, H_MM = 113.5, 76.0

INK = "#202124"
DARK_SLATE = "#34414B"
SLATE = "#687782"
GREY = "#777A7E"
MID_GREY = "#AAB1B6"
LIGHT_GREY = "#DCE1E4"
PALE_GREY = "#F4F5F6"

ID_EDGE = "#688896"
ID_FACE = "#EAF1F4"
PRED_EDGE = "#2F7F78"
PRED_FACE = "#DDEDEA"
AMBER = "#D59A2B"
AMBER_FACE = "#FAF0DC"
WARM = "#B56A50"
WARM_FACE = "#F7E9E3"


def apply_style() -> None:
    """House style, then the few panel-local overrides this schematic needs.

    The font stack deliberately comes from nm_style rather than being restated
    here. A local stack that lists Arial first and DejaVu Sans last renders in
    DejaVu Sans, silently and without error, on any machine without Arial
    installed, which is neither Arial- nor Helvetica-metric and mixes a second
    typeface into the figure. nm_style.resolve_sans() raises instead.
    """
    S.apply_rcparams()
    plt.rcParams.update(
        {
            "font.size": 6.0,
            "axes.linewidth": 0.55,
            "text.color": INK,
            # This panel writes its own exact canvas; a tight bbox would trim the
            # declared width and change the type size at placement.
            "savefig.bbox": None,
            "savefig.pad_inches": 0.0,
        }
    )


def rounded_box(
    ax,
    x: float,
    y: float,
    width: float,
    height: float,
    *,
    face="white",
    edge=LIGHT_GREY,
    lw=0.65,
    radius=1.2,
    zorder=1,
):
    patch = FancyBboxPatch(
        (x, y),
        width,
        height,
        boxstyle=f"round,pad=0,rounding_size={radius}",
        facecolor=face,
        edgecolor=edge,
        linewidth=lw,
        zorder=zorder,
    )
    ax.add_patch(patch)
    return patch


def arrow(ax, start, end, *, color=GREY, lw=0.65, mutation=5.2, zorder=6):
    patch = FancyArrowPatch(
        start,
        end,
        arrowstyle="-|>",
        mutation_scale=mutation,
        linewidth=lw,
        color=color,
        shrinkA=0,
        shrinkB=0,
        zorder=zorder,
    )
    ax.add_patch(patch)
    return patch


def cell_cloud(ax, cx, cy, width, height, *, mode="variant", seed=0, n=7):
    """Small deterministic cell-distribution glyph; not an observation."""
    if mode == "wt":
        face, edge, dash, point = "white", MID_GREY, (0, (1.8, 1.3)), MID_GREY
    elif mode == "sibling":
        face, edge, dash, point = PALE_GREY, SLATE, "solid", SLATE
    else:
        face, edge, dash, point = "white", DARK_SLATE, "solid", DARK_SLATE
    ax.add_patch(
        Ellipse(
            (cx, cy),
            width,
            height,
            facecolor=face,
            edgecolor=edge,
            linewidth=0.55,
            linestyle=dash,
            zorder=4,
        )
    )
    templates = [
        [(-0.54, 0.18), (-0.30, -0.36), (-0.12, 0.47), (0.18, -0.05), (0.36, 0.34), (0.55, -0.28), (0.02, -0.51)],
        [(-0.48, -0.29), (-0.35, 0.39), (-0.06, 0.08), (0.14, -0.45), (0.31, 0.46), (0.53, -0.02), (0.02, 0.53)],
        [(-0.57, 0.02), (-0.27, 0.43), (-0.19, -0.42), (0.09, 0.17), (0.34, -0.31), (0.52, 0.31), (0.01, -0.03)],
    ]
    points = templates[seed % len(templates)][:n]
    ax.scatter(
        [cx + x_value * width / 2 for x_value, _ in points],
        [cy + y_value * height / 2 for _, y_value in points],
        s=2.0,
        facecolor=point,
        edgecolor="none",
        zorder=5,
    )


def condition_card(ax, x, y, width, height, *, condition, verdict, body, edge, face):
    rounded_box(ax, x, y, width, height, face=face, edge=edge, lw=0.75, radius=1.1, zorder=2)
    ax.text(x + 2.0, y + height - 2.1, condition, fontsize=5.0, color=GREY, ha="left", va="top")
    ax.text(x + 2.0, y + height - 5.4, verdict, fontsize=6.0, color=edge, fontweight="bold", ha="left", va="top")
    ax.text(x + 2.0, y + 1.35, body, fontsize=5.0, color=INK, ha="left", va="bottom", linespacing=1.0)


def main() -> None:
    apply_style()
    fig = plt.figure(figsize=(W_MM * MM, H_MM * MM))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W_MM)
    ax.set_ylim(0, H_MM)
    ax.axis("off")

    # Problem-first headline: the scientific claim precedes all implementation detail.
    ax.text(
        2.5,
        73.5,
        "Allele-specific prediction is a nested solvability problem",
        fontsize=7.7,
        fontweight="bold",
        ha="left",
        va="top",
    )
    ax.text(
        2.5,
        69.2,
        "A model can be judged only after the measurement detects a response and reproducibly identifies the allele.",
        fontsize=5.45,
        color=GREY,
        ha="left",
        va="top",
    )

    # Nested requirements. The geometry, not an arrow-chain, carries the logic.
    # Outer box right edge 64.5, identification 61.5, prediction 57.0. The
    # right-hand insets are 3.0 and 4.5 mm rather than the 7.0 and 7.5 they
    # were at 122.5 mm: the identification box is sized by its own caption,
    # "Identification pass = allele-benchmarkable ground truth", which measures
    # 45.97 mm at 5 pt and cannot be set smaller. The nesting still reads.
    rounded_box(ax, 3.0, 19.5, 61.5, 45.5, face=PALE_GREY, edge=MID_GREY, lw=0.85, radius=1.6)
    ax.text(6.0, 62.2, "1  DETECTION", fontsize=6.15, fontweight="bold", color=DARK_SLATE, ha="left", va="center")
    ax.text(24.0, 62.2, "variant vs wild type", fontsize=5.15, color=GREY, ha="left", va="center")
    cell_cloud(ax, 49.0, 62.1, 5.7, 3.7, mode="variant", seed=1, n=6)
    ax.text(53.5, 62.1, "vs", fontsize=5.0, color=GREY, ha="center", va="center")
    cell_cloud(ax, 58.0, 62.1, 5.7, 3.7, mode="wt", seed=2, n=6)

    rounded_box(ax, 10.5, 25.2, 51.0, 31.3, face=ID_FACE, edge=ID_EDGE, lw=0.9, radius=1.5, zorder=2)
    ax.text(13.5, 53.6, "2  IDENTIFICATION", fontsize=6.15, fontweight="bold", color=ID_EDGE, ha="left", va="center")
    ax.text(37.0, 53.6, "variant vs siblings", fontsize=5.1, color=GREY, ha="left", va="center")
    cell_cloud(ax, 50.0, 49.7, 5.4, 3.5, mode="variant", seed=3, n=6)
    cell_cloud(ax, 56.5, 49.7, 5.4, 3.5, mode="sibling", seed=4, n=6)

    rounded_box(ax, 20.0, 30.6, 37.0, 15.9, face=PRED_FACE, edge=PRED_EDGE, lw=1.05, radius=1.4, zorder=3)
    ax.text(23.0, 43.1, "3  PREDICTION", fontsize=6.2, fontweight="bold", color=PRED_EDGE, ha="left", va="center")
    ax.text(23.0, 38.8, r"molecular features  →  $\hat{\delta}_v$", fontsize=5.55, fontweight="bold", color=INK, ha="left", va="center")
    ax.text(23.0, 34.3, "correct held-out allele response?", fontsize=5.0, color=GREY, ha="left", va="center")

    # Explicit containment statement and benchmark meaning.
    ax.text(13.5, 28.3, "Identification pass = allele-benchmarkable ground truth", fontsize=5.0, color=ID_EDGE, fontweight="bold", ha="left", va="center")
    ax.text(6.0, 22.2, "Every inner claim requires every outer layer.", fontsize=5.25, color=DARK_SLATE, fontweight="bold", ha="left", va="center")

    # Short inward arrows reinforce, but do not replace, the nested geometry.
    arrow(ax, (7.3, 58.2), (11.0, 55.6), color=MID_GREY, lw=0.55, mutation=4.5)
    arrow(ax, (14.4, 49.2), (20.3, 45.2), color=ID_EDGE, lw=0.6, mutation=4.5)

    # Reviewer-safe verdict map: failure is localized to measurement or model.
    ax.text(69.0, 63.8, "Verdict: where is the binding limit?", fontsize=6.35, fontweight="bold", color=INK, ha="left", va="center")
    condition_card(
        ax,
        69.0,
        47.4,
        41.5,
        13.4,
        condition="Detection or identification fails",
        verdict="MEASUREMENT-LIMITED",
        body="Ground truth cannot reward an allele model.\nDo not rank models at allele resolution.",
        edge=WARM,
        face=WARM_FACE,
    )
    condition_card(
        ax,
        69.0,
        33.1,
        41.5,
        11.7,
        condition="Identification passes; model < ceiling",
        verdict="MODEL-LIMITED",
        body="A meaningful computational gap.",
        edge=PRED_EDGE,
        face=PRED_FACE,
    )
    condition_card(
        ax,
        69.0,
        20.0,
        41.5,
        10.6,
        condition="Model reaches the measurement ceiling",
        verdict="RESOLVED AT CURRENT CEILING",
        body="Prediction is supported at this resolution.",
        edge=DARK_SLATE,
        face=PALE_GREY,
    )

    # A visible contradiction to the conventional shortcut is the panel's key new insight.
    # Terracotta is reserved for conceptual cautions. Amber remains exclusive
    # to externally annotated protein hotspots elsewhere in Fig. 1.
    rounded_box(ax, 3.0, 10.6, 107.5, 6.8, face=WARM_FACE, edge=WARM, lw=0.7, radius=1.1)
    ax.add_patch(Circle((7.0, 14.0), 1.45, facecolor=WARM, edgecolor="none", zorder=4))
    ax.text(7.0, 14.0, "!", fontsize=6.0, fontweight="bold", color="white", ha="center", va="center", zorder=5)
    ax.text(10.0, 14.0, r"Direction recovered (Pearson-$\delta$)  $\ne$  allele identified or predicted", fontsize=5.75, fontweight="bold", color=INK, ha="left", va="center")

    # Task constraints remain visible but subordinate to the scientific problem.
    rounded_box(ax, 3.0, 2.0, 107.5, 6.3, face="white", edge=LIGHT_GREY, lw=0.6, radius=1.0)
    ax.text(6.0, 5.15, "Held-out task", fontsize=5.2, fontweight="bold", color=DARK_SLATE, ha="left", va="center")
    ax.text(23.0, 5.15, r"observed cells define $\delta_v$ and resolution", fontsize=5.0, color=GREY, ha="left", va="center")
    ax.text(61.0, 5.15, "•", fontsize=5.0, color=MID_GREY, ha="center", va="center")
    ax.text(64.0, 5.15, r"model sees only $\theta$/ESM and predicts $\hat{\delta}_v$", fontsize=5.0, color=GREY, ha="left", va="center")

    S.save(fig, HERE / "fig1a_nested_solvability", exact=True,
           formats=("pdf", "png", "svg", "tiff"))


if __name__ == "__main__":
    main()
