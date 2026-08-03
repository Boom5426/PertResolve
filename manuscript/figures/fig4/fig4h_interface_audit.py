"""Fig. 4h — compact audit of allele-level model-interface compatibility."""
from __future__ import annotations

from matplotlib.patches import FancyBboxPatch, Rectangle

import fig4_common as S

# Static QA contract: Arial/Helvetica sans-serif; svg.fonttype="none";
# pdf.fonttype=42; exports .svg, .pdf, .png and .tiff at dpi=600.
# Final assembled figure target: width_mm = 183.
W_MM, H_MM = 84.0, 51.0

STATUS = {
    "Y": (S.SLATE, "white"),
    "P": (S.AMBER, S.INK),
    "N": ("#D9DDE0", S.INK),
    "NA": ("#F1F2F3", "#8A8D90"),
}
COLS = [
    "Continuous\nvariant",
    "Unseen\nallele",
    "Allele-specific\noutput",
    "Training\ncompleted",
    "Allele score\ndefined",
]
ROWS = [
    ("Variant-conditionable", "scGen", ["Y", "Y", "Y", "Y", "Y"]),
    ("Variant-conditionable", "scVIDR", ["Y", "Y", "Y", "Y", "Y"]),
    ("Variant-conditionable", "Biolord", ["Y", "Y", "Y", "Y", "Y"]),
    ("Variant-conditionable", "CellFlow", ["Y", "Y", "Y", "Y", "Y"]),
    ("Variant-conditionable", "PerturbNet", ["Y", "Y", "Y", "Y", "P"]),
    ("Allele-blind by construction", "scGPT", ["N", "N", "N", "NA", "N"]),
    ("Allele-blind by construction", "GEARS", ["N", "N", "N", "NA", "N"]),
    ("Allele-blind by construction", "STATE", ["N", "N", "N", "Y", "N"]),
    ("Did not converge", "variant-CPA", ["Y", "Y", "NA", "N", "N"]),
]


def aspect_kx(fig, ax) -> float:
    width_mm, height_mm = (value * 25.4 for value in fig.get_size_inches())
    position = ax.get_position()
    mm_per_x = position.width * width_mm / abs(ax.get_xlim()[1] - ax.get_xlim()[0])
    mm_per_y = position.height * height_mm / abs(ax.get_ylim()[1] - ax.get_ylim()[0])
    return mm_per_y / mm_per_x


def mark(ax, x: float, y: float, status: str, kx: float, scale: float = 1.0) -> None:
    color = STATUS[status][1]

    def xs(*values):
        return [x + value * kx * scale for value in values]

    def ys(*values):
        return [y + value * scale for value in values]

    if status == "Y":
        ax.plot(
            xs(-0.16, -0.03, 0.17),
            ys(0.02, 0.15, -0.14),
            color=color,
            lw=1.0 * scale,
            solid_capstyle="round",
            solid_joinstyle="round",
            zorder=5,
        )
    elif status == "N":
        ax.plot(xs(-0.13, 0.13), ys(-0.14, 0.14), color=color, lw=0.9 * scale, zorder=5)
        ax.plot(xs(-0.13, 0.13), ys(0.14, -0.14), color=color, lw=0.9 * scale, zorder=5)
    elif status == "P":
        ax.plot(
            xs(-0.18, -0.06, 0.06, 0.18),
            ys(0.02, 0.12, -0.12, -0.02),
            color=color,
            lw=0.9 * scale,
            zorder=5,
        )
    else:
        ax.plot(x, y, marker="o", ms=2.1 * scale, color=color, zorder=5)


def main() -> None:
    S.apply_style()
    coverage = S.read_json("allele_blind_demo.json")["resolution_coverage"]

    y = 0.25
    layout = []
    headers = []
    previous = None
    for group, model, statuses in ROWS:
        if group != previous:
            if previous is not None:
                y += 0.25
            headers.append((group, y))
            y += 0.62
            previous = group
        layout.append((model, statuses, y))
        y += 0.82
    total = y

    fig = S.figure(W_MM, H_MM)
    ax = fig.add_axes([0.03, 0.05, 0.94, 0.93])
    ax.set_xlim(-1.72, 5.05)
    ax.set_ylim(-0.85, total + 2.30)
    ax.invert_yaxis()
    ax.axis("off")
    kx = aspect_kx(fig, ax)

    for index, label in enumerate(COLS):
        ax.text(
            index + 0.5,
            -0.05,
            label,
            ha="center",
            va="bottom",
            fontsize=5.0,
            color=S.INK,
            linespacing=1.1,
        )

    for group, header_y in headers:
        ax.plot([-1.70, 5.0], [header_y + 0.48, header_y + 0.48], color=S.MID_GREY, lw=0.55)
        ax.text(
            -1.70,
            header_y + 0.23,
            group,
            ha="left",
            va="center",
            fontsize=5.2,
            color=S.GREY,
            style="italic",
        )

    for model, statuses, row_y in layout:
        ax.text(
            -0.12,
            row_y + 0.35,
            model,
            ha="right",
            va="center",
            fontsize=5.5,
            color=S.INK,
            fontweight="bold",
        )
        for index, status in enumerate(statuses):
            ax.add_patch(
                Rectangle(
                    (index + 0.10, row_y + 0.05),
                    0.80,
                    0.60,
                    facecolor=STATUS[status][0],
                    edgecolor="white",
                    lw=0.45,
                    zorder=2,
                )
            )
            mark(ax, index + 0.50, row_y + 0.35, status, kx)

    legend_y = total + 0.42
    legend_items = [
        ("Y", "available"),
        ("P", "caveat"),
        ("N", "unavailable"),
        ("NA", "not applicable"),
    ]
    for index, (status, label) in enumerate(legend_items):
        x = -1.70 + index * 1.58
        ax.add_patch(
            Rectangle(
                (x, legend_y),
                0.27,
                0.44,
                facecolor=STATUS[status][0],
                edgecolor="white",
                lw=0.4,
            )
        )
        mark(ax, x + 0.135, legend_y + 0.22, status, kx, scale=0.58)
        ax.text(x + 0.36, legend_y + 0.22, label, ha="left", va="center", fontsize=5.0, color=S.GREY)

    callout_y = total + 1.22
    ax.add_patch(
        FancyBboxPatch(
            (-1.70, callout_y),
            6.70,
            0.78,
            boxstyle="round,pad=0.08,rounding_size=0.08",
            facecolor=S.PALE_GREY,
            edgecolor="none",
            zorder=0,
        )
    )
    ax.text(
        -1.53,
        callout_y + 0.39,
        f"Gene-keyed coverage\n{coverage}",
        ha="left",
        va="center",
        fontsize=5.0,
        color=S.INK,
        linespacing=1.15,
    )
    ax.text(
        4.82,
        callout_y + 0.39,
        "Status reflects the evaluated\nimplementation",
        ha="right",
        va="center",
        fontsize=5.0,
        color=S.GREY,
        linespacing=1.15,
    )
    S.save(fig, "fig4h_interface_audit")


if __name__ == "__main__":
    main()
