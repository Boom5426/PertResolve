"""Fig 2d: PDS-vs-Pearson scatter (direction without ranking).

Message: in-house feature models recover the *direction* of the perturbation
delta well (mean pearson_delta ~0.55-0.65) yet their per-variant *ranking* skill
(PDS) sits at chance (~0.49-0.52). All learned methods therefore land in a single
"direction without ranking" region straddling PDS=0.5; only the WT-null, which
predicts no delta, drops to pearson_delta=0.

Numbers (all from committed sources, no fabrication):
  x = PDS            from D.definitive()  (results/canonical/definitive_summary.csv)
  y = mean pearson_delta per method, from D.exttheta() (results/results_v4_exttheta.csv)
Points = 18 in-house feature-model heads + Gene-mean + WT-null.

Nature Methods pass:
  * No legend. The feature-space key is drawn once, in Fig 2b; here the three
    callouts name one method from each feature space, so the colours are direct-
    labelled rather than keyed a second time.
  * Reference points use the same grey circle as Fig 2b/2c (one encoding per
    variable) and are direct-labelled, instead of bespoke diamond/square markers.
  * The former y = 0.30 divider and quadrant shading were an unexplained constant
    and pure decoration; both removed. Only the 0.50 chance line, which is the
    metric's analytic null, is drawn.
  * Axis limits tightened to the data range (PDS 0.487-0.517, Pearson 0-0.65).
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import nm_style as S
import remote_data as D

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

W_MM, H_MM = 76.0, 44.0
AX = (12.5, 9.0, 61.5, 31.0)   # left, bottom, width, height in mm


def canvas():
    fig = plt.figure(figsize=(W_MM * S.MM, H_MM * S.MM))
    bg = fig.add_axes([0, 0, 1, 1])
    bg.set_axis_off()
    bg.set_xlim(0, 1)
    bg.set_ylim(0, 1)
    bg.add_patch(Rectangle((0, 0), 1, 1, facecolor="white", edgecolor="none",
                           zorder=-10))
    l, b, w, h = AX
    return fig, fig.add_axes([l / W_MM, b / H_MM, w / W_MM, h / H_MM])


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

    pds = D.definitive().set_index("method")["PDS"]
    pear = D.exttheta().groupby("method")["pearson_delta"].mean()

    inhouse = [m for m in pds.index if D.is_inhouse(m)]

    fig, ax = canvas()

    X0, X1 = 0.482, 0.523
    Y0, Y1 = -0.03, 0.705

    # analytic chance level of PDS; the only reference line on the panel
    ax.axvline(0.50, color=S.GREY, lw=0.5, ls=(0, (3, 2)), zorder=1)

    # ---- points: same circle + white edge as the Fig 2b/2c forests ---------
    for m in inhouse:
        ax.scatter(pds[m], pear[m], s=13,
                   facecolor=S.FEATURE_COLORS[D.feature_space(m)],
                   edgecolor="white", linewidth=0.4, zorder=4)
    for m in ("Gene-mean", "WT-null"):
        ax.scatter(pds[m], pear[m], s=13,
                   facecolor=S.FEATURE_COLORS["reference"],
                   edgecolor="white", linewidth=0.4, zorder=5)

    # ---- direct labels: one method per feature space ----------------------
    # Every label is horizontally centred on the x of the point it names, so the
    # leader is strictly vertical and the label can only be read against its own
    # datum. Two labels sit above the cloud and two below; a single row would
    # force horizontal offsets, which is how earlier versions ended up centring
    # a name over a neighbouring point.
    SHORT = {"Ridge-esm": "Ridge ESM", "MLP-esm+theta": r"MLP ESM+$\theta$",
             "Lasso-theta": r"Lasso $\theta$", "Gene-mean": "Gene-mean"}
    LABEL_Y_ABOVE, LABEL_Y_BELOW = 0.678, 0.478

    def callout(m, ty):
        ax.annotate(SHORT[m], (pds[m], pear[m]), (pds[m], ty),
                    ha="center", va="center", fontsize=5.5, color=S.INK,
                    arrowprops=dict(arrowstyle="-", lw=0.4, color=S.GREY,
                                    shrinkA=0.5, shrinkB=1.5), zorder=6)
        print(f"callout {SHORT[m]:18} anchored at "
              f"({pds[m]:.4f}, {pear[m]:.4f}), text x = {pds[m]:.4f}")

    callout("Gene-mean", LABEL_Y_ABOVE)      # x = 0.489
    callout("Ridge-esm", LABEL_Y_ABOVE)      # x = 0.497
    callout("MLP-esm+theta", LABEL_Y_BELOW)  # x = 0.493
    callout("Lasso-theta", LABEL_Y_BELOW)    # x = 0.508

    # the second reference, direct-labelled at its point
    ax.text(pds["WT-null"] + 0.0012, pear["WT-null"] + 0.008, "WT-null",
            ha="left", va="bottom", fontsize=5.5, color=S.INK, zorder=6)

    # chance annotation, at the foot of the divider
    ax.text(0.4992, 0.055, "chance", ha="right", va="bottom", fontsize=5.5,
            color=S.GREY, zorder=6)

    # panel message, in the empty band between the WT-null point and the cloud
    ax.text(0.4835, 0.300, "direction recovered,\nranking at chance",
            ha="left", va="center", fontsize=5.8, color=S.INK, style="italic",
            linespacing=1.3, zorder=6)

    # ---- axes -------------------------------------------------------------
    ax.set_xlim(X0, X1)
    ax.set_ylim(Y0, Y1)
    ax.set_xticks([0.49, 0.50, 0.51, 0.52])
    ax.set_yticks([0.0, 0.2, 0.4, 0.6])
    ax.set_xlabel(r"PDS$_{cos}$ (allele discrimination)", labelpad=1.5)
    ax.set_ylabel(r"Pearson-$\delta$ (direction recovery)", labelpad=1.5)
    ax.tick_params(pad=1.5)
    S.despine(ax, keep=("left", "bottom"))

    S.save(fig, os.path.join(os.path.dirname(os.path.abspath(__file__)), "fig2d_scatter"))

    print("PDS range", round(pds[inhouse].min(), 3), round(pds[inhouse].max(), 3))
    print("Pearson range", round(pear[inhouse].min(), 3), round(pear[inhouse].max(), 3))


if __name__ == "__main__":
    main()
