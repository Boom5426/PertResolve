"""Shared Nature-Methods-style matplotlib configuration for AllelePerturb Figure 1.

All data-direct panels import from here so fonts, sizes, colours and export
settings stay identical across separate panel files. Text is exported as
editable text (pdf/ps fonttype 42) so panels can be relabelled in Illustrator.

Nature Methods conventions applied:
  - sans-serif (Helvetica); Liberation/Nimbus Sans are metric-compatible stand-ins
  - 5-7 pt type, thin (0.5 pt) axes, no top/right spines
  - single-column width 89 mm, double-column 183 mm
  - house gene palette (GENE_COLORS) plus redundant marker shapes (GENE_MARKERS),
    because the palette is NOT colour-blind-safe on its own: it replaced an
    Okabe-Ito set, and TP53 blue, GATA1 purple and JAK1 green sit at 8-bit greys
    of 129, 133 and 135, with TP53 and GATA1 near-identical under deuteranopia.
    Any panel where the gene is encoded by hue alone must also vary the marker.
  - vector PDF deliverable + 600 dpi PNG preview
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

# ---- font resolution ------------------------------------------------------
# Nature Methods asks for Arial or Helvetica. Liberation Sans is metrically
# identical to Arial and Nimbus Sans to Helvetica, so either is an acceptable
# stand-in, but only ONE family may appear in a figure. Resolving the stack once
# here (rather than letting matplotlib fall back per text element) is what keeps
# regular text and mathtext in the same family: `mathtext.rm` takes a single font
# name, not a fallback list, so a hard-coded value there silently disagrees with
# `font.sans-serif` the moment a preferred family is installed.
#
# DejaVu Sans, matplotlib's bundled default, is deliberately NOT in the stack: it
# is neither Arial- nor Helvetica-metric, and an unnoticed fallback to it is what
# put a second typeface into Extended Data Fig. 1. Missing fonts raise.
SANS_STACK = ("Arial", "Helvetica", "Liberation Sans", "Nimbus Sans")


def resolve_sans(stack: tuple[str, ...] = SANS_STACK) -> str:
    """Return the first installed family in `stack`; raise if none is present.

    Raising is intentional: a silent fallback to DejaVu Sans produces a figure
    that looks finished and fails production.
    """
    from matplotlib import font_manager

    installed = {f.name for f in font_manager.fontManager.ttflist}
    for name in stack:
        if name in installed:
            return name
    raise RuntimeError(
        "No Nature-Methods-acceptable sans font found. Install one of "
        f"{', '.join(stack)} (Debian/Ubuntu: `fonts-liberation`, or "
        "`ttf-mscorefonts-installer` for Arial itself) and re-run. "
        "Refusing to fall back to DejaVu Sans, which is not Arial- or "
        "Helvetica-metric."
    )


SANS = resolve_sans()

# ---- geometry -------------------------------------------------------------
MM = 1.0 / 25.4  # millimetres -> inches
COL1_MM = 89.0   # Nature single-column width
COL2_MM = 183.0  # Nature double-column width

# ---- house gene palette (shared across all figures) -----------------------
# Author-specified Nature-style palette (Fig 2/3 spec); one edit here recolours
# every figure. Previous Okabe-Ito set kept for reference in git history.
GENE_COLORS = {
    "TP53": "#5185C0",   # blue
    "KRAS": "#E99D4E",   # orange
    "GATA1": "#8281B9",  # purple
    "JAK1": "#55966B",   # green
}
GENE_ORDER = ["TP53", "KRAS", "GATA1", "JAK1"]

# Redundant encoding for gene identity. GENE_COLORS is not separable in greyscale
# or under deuteranopia (see module docstring), so every panel that identifies a
# gene by hue alone pairs the hue with this shape. Colour is unchanged; the shape
# carries the same information a second way. Shapes are chosen to stay distinct at
# 3-16 pt^2 marker areas: circle, square, triangle, diamond.
GENE_MARKERS = {"TP53": "o", "KRAS": "s", "GATA1": "^", "JAK1": "D"}
GREY = "#7A7A7A"

# canonical feature-space palette (Fig 2b/2c/2d/4e): keep ONE mapping across panels
FEATURE_COLORS = {
    "theta": "#9AA7B3",       # light slate (feature spaces = slate ramp, distinct from gene hues)
    "ESM": "#5F6B76",         # mid slate
    "ESM+theta": "#2E3742",   # dark slate
    "reference": "#7A7A7A",   # grey (Gene-mean / WT-null)
}
LIGHT_GREY = "#D9D9D9"
HOTSPOT = "#E69F00"      # orange highlight for hotspot/pathogenic residues
INK = "#1A1A1A"

# canonical dataset palette (Fig 6, Extended Data Fig 1): the four allele genes keep
# their gene hues, so colour always means "gene" where a gene is plotted, and the
# external gene- and drug-level atlases take the neutral slate ramp rather than
# invented hues. Lives here, not in a figure module, because more than one figure
# plots this same seven-dataset set and they must agree.
DATASET_COLORS = {
    "TP53": GENE_COLORS["TP53"], "KRAS": GENE_COLORS["KRAS"],
    "GATA1": GENE_COLORS["GATA1"], "JAK1": GENE_COLORS["JAK1"],
    "Replogle": FEATURE_COLORS["ESM+theta"],   # dark slate  #2E3742
    "Norman": FEATURE_COLORS["ESM"],           # mid slate   #5F6B76
    "Adamson": FEATURE_COLORS["theta"],        # light slate #9AA7B3
    "VCC": GREY, "sciPlex": LIGHT_GREY,
}

# canonical protein lengths (UniProt: P04637, P01116, P15976, P23458)
PROTEIN_LEN = {"TP53": 393, "KRAS": 189, "GATA1": 413, "JAK1": 1154}


def apply_rcparams() -> None:
    """Set global rcParams to the Nature-Methods house style."""
    plt.rcParams.update({
        "font.family": "sans-serif",
        # One resolved family, not a fallback list: see SANS_STACK / resolve_sans.
        "font.sans-serif": [SANS],
        # Math text ($...$) and italics must stay in the SAME family as body text.
        # Without this, matplotlib renders math with its bundled DejaVu Sans, which
        # has visibly different glyph shapes and mixed two font families inside
        # single panels (seen in Fig 3 and Fig 4). These are driven off SANS so the
        # two settings cannot drift apart when a preferred font is installed.
        "mathtext.fontset": "custom",
        "mathtext.rm": SANS,
        "mathtext.it": f"{SANS}:italic",
        "mathtext.bf": f"{SANS}:bold",
        # "custom" leaves sf/cal/tt at matplotlib's generic families; cal resolves to
        # "cursive", which exists nowhere here and warns its way back to DejaVu Sans.
        # One family per figure means all of them, even the ones we do not use.
        "mathtext.sf": SANS,
        "mathtext.cal": f"{SANS}:italic",
        "mathtext.tt": SANS,
        "mathtext.default": "regular",
        "font.size": 7,
        "axes.titlesize": 7,
        "axes.labelsize": 7,
        "xtick.labelsize": 6,
        "ytick.labelsize": 6,
        "legend.fontsize": 6,
        "figure.titlesize": 7,
        "axes.linewidth": 0.5,
        "lines.linewidth": 0.75,
        "patch.linewidth": 0.5,
        "xtick.major.width": 0.5,
        "ytick.major.width": 0.5,
        "xtick.major.size": 2.2,
        "ytick.major.size": 2.2,
        "xtick.minor.width": 0.4,
        "ytick.minor.width": 0.4,
        "axes.labelcolor": INK,
        "axes.edgecolor": INK,
        "xtick.color": INK,
        "ytick.color": INK,
        "text.color": INK,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "legend.frameon": False,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.01,
        "pdf.fonttype": 42,   # keep text editable (TrueType) in the PDF
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    })


def panel(w_mm: float, h_mm: float):
    """Return (fig, ax) sized in millimetres."""
    fig, ax = plt.subplots(figsize=(w_mm * MM, h_mm * MM))
    return fig, ax


def despine(ax, keep=("left", "bottom")) -> None:
    for side in ("top", "right", "left", "bottom"):
        ax.spines[side].set_visible(side in keep)


def save(fig, stem: Path | str) -> None:
    """Save a panel as vector PDF (deliverable) + 600 dpi PNG (preview).

    Line art and text stay vector (resolution-independent); dpi=600 only governs
    any element explicitly marked ``rasterized=True`` so it also meets the
    >=600 dpi Nature line-art bar rather than matplotlib's ~100 dpi default.
    """
    stem = Path(stem)
    fig.savefig(stem.with_suffix(".pdf"), dpi=600)
    fig.savefig(stem.with_suffix(".png"), dpi=600)
    plt.close(fig)


def repo_root() -> Path:
    """Repository root, found by walking up until data/allele_perturb_bench.csv exists.

    Location-independent, so this shared module works whether it is imported from
    manuscript/figures/ or a figN/ subfolder (no hard-coded path, no parent count).
    """
    here = Path(__file__).resolve()
    for anc in here.parents:
        if (anc / "data" / "allele_perturb_bench.csv").exists():
            return anc
    raise RuntimeError("repo root (with data/allele_perturb_bench.csv) not found above " + str(here))


def load_bench(exclude_wt: bool = True) -> "pd.DataFrame":
    """Load the committed benchmark table.

    The CSV has 472 rows, of which 2 are wild-type reference rows
    (KRAS WT, 644 cells; GATA1 WT, 38,276 cells). ``exclude_wt=True`` drops
    them so per-variant panels count the 470 real protein-coding variants and
    the GATA1 WT pool does not appear as a spurious high-depth "variant".
    Set ``exclude_wt=False`` to reproduce the manuscript's current
    WT-inclusive counts (472 / 98 / 93 / 255 / 26).
    """
    df = pd.read_csv(repo_root() / "data" / "allele_perturb_bench.csv")
    if exclude_wt:
        df = df[df.variant.astype(str).str.upper() != "WT"].copy()
    return df
