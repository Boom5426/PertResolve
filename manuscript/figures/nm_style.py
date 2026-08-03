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
from matplotlib.transforms import Bbox
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

# Non-ASCII codepoints the figures actually set. Being installed is not enough:
# a family that lacks one of these does not fail, it silently pulls that single
# glyph from matplotlib's STIXGeneral fallback, which puts a third typeface into
# the figure for one character. The Arial shipped by Debian/Ubuntu's
# ttf-mscorefonts-installer is exactly this case: it has no U+0302, so every
# \hat{...} in Fig. 1 renders its accent in STIXGeneral. Liberation Sans is
# Arial-metric and covers the whole set, so it wins the stack on this machine.
REQUIRED_CODEPOINTS = {
    0x0302: r"combining circumflex (\hat accent)",
    0x0394: "Greek capital delta",
    0x03B4: "Greek small delta",
    0x03B8: "Greek small theta",
    0x00B7: "middle dot",
    0x2013: "en dash",
    0x2022: "bullet",
    0x2192: "rightwards arrow",
    0x2260: "not equal to",
}


def _missing_codepoints(path: str) -> list[str]:
    """Names of REQUIRED_CODEPOINTS absent from the font file at `path`."""
    from matplotlib.ft2font import FT2Font

    face = FT2Font(path)
    return [name for cp, name in sorted(REQUIRED_CODEPOINTS.items())
            if face.get_char_index(cp) == 0]


def resolve_sans(stack: tuple[str, ...] = SANS_STACK) -> str:
    """Return the first family in `stack` that is installed AND fully covers
    REQUIRED_CODEPOINTS; raise if none does.

    Raising is intentional: a silent fallback to DejaVu Sans produces a figure
    that looks finished and fails production. Coverage is checked for the same
    reason, one level down: a family can be present and still hand one glyph to
    STIXGeneral without any warning.
    """
    from matplotlib import font_manager

    installed = {f.name for f in font_manager.fontManager.ttflist}
    rejected = []
    for name in stack:
        if name not in installed:
            continue
        path = font_manager.findfont(
            font_manager.FontProperties(family=name), fallback_to_default=False)
        missing = _missing_codepoints(path)
        if not missing:
            return name
        rejected.append(f"{name} (missing {', '.join(missing)})")
    raise RuntimeError(
        "No Nature-Methods-acceptable sans font with complete glyph coverage. "
        f"Tried {', '.join(stack)}"
        + (f"; rejected for incomplete coverage: {'; '.join(rejected)}" if rejected else "")
        + ". Install one of them (Debian/Ubuntu: `fonts-liberation`) and re-run. "
        "Refusing to fall back to DejaVu Sans, which is not Arial- or "
        "Helvetica-metric, or to render single glyphs in STIXGeneral."
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

# Semantic (non-gene) colours used by the Fig. 1 schematic and split panels.
# Kept here rather than in a panel file because b/e/f must agree on them: amber
# (HOTSPOT) is reserved for external hotspot annotation and must never also mean
# "held out", which is why HELDOUT is a separate steel blue.
HELDOUT = "#4F6D7A"          # steel blue: held-out variants in a split
RESOLUTION = "#2A7F78"       # teal: resolution sufficient / proceed branch
RESOLUTION_LIGHT = "#E6F1EF"
CAUTION = "#B5654D"          # muted brick: insufficient-resolution branch

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


FONT_CONTRACT = Path(__file__).resolve().parent / "nm_font_resolved.tex"


def write_font_contract(path: Path | str = FONT_CONTRACT) -> Path:
    """Write the resolved family as a one-line \\setsansfont for the compositors.

    The panel letters stamped by each figN_assemble.tex must be in the same
    family as the panel content, or the composite embeds two. nm_fonts.tex used
    to re-derive the family with its own "is it installed" probe, which agreed
    with resolve_sans() only as long as resolve_sans() asked the same question.
    It no longer does: a family can be installed and still lack a glyph the
    figures set (see REQUIRED_CODEPOINTS). Rather than teach LaTeX to check
    glyph coverage, the Python side records what it actually used and the
    LaTeX side reads it.
    """
    path = Path(path)
    path.write_text(
        "% Generated by nm_style.write_font_contract(); do not edit by hand.\n"
        "% The family the matplotlib panels actually resolved to, so the panel\n"
        "% letters match. Regenerate with: python fig1/check_fig1_panels.py\n"
        f"\\setsansfont{{{SANS}}}\n"
    )
    return path


def vector_gradient(ax, cmap, extent, *, steps: int = 96, zorder: float = 2) -> None:
    """Draw a horizontal colour ramp as adjacent rectangles, not an image.

    ``imshow`` rasterises, and a raster object in a figure that is otherwise
    vector fails the composite's zero-raster contract and prints at the image's
    resolution rather than the device's. At 96 steps the banding is well below
    what a 600 dpi press can resolve across a 30 mm bar.
    """
    from matplotlib.patches import Rectangle

    x0, x1, y0, y1 = extent
    width = (x1 - x0) / steps
    for index in range(steps):
        ax.add_patch(Rectangle(
            (x0 + index * width, y0), width * 1.02, y1 - y0,
            facecolor=cmap((index + 0.5) / steps), edgecolor="none",
            linewidth=0, zorder=zorder))


def wrap_to_width(fig, text: str, *, fontsize: float, weight: str = "normal",
                  max_frac: float = 0.96, x: float = 0.02) -> str:
    """Greedily break `text` on spaces so no line exceeds the figure width.

    Panel headings are prose, and prose written for a 132 mm canvas runs off a
    89 mm one. Measuring and breaking is better than trusting each panel to
    hand-shorten its own strings, because the canvas widths change and the
    strings then have to change with them. Returns the text with newlines
    inserted; the caller still places it.
    """
    renderer = fig.canvas.get_renderer()
    limit = (max_frac - x) * fig.get_size_inches()[0] * fig.dpi

    def width(candidate: str) -> float:
        probe = fig.text(0, 0, candidate, fontsize=fontsize, fontweight=weight)
        w = probe.get_window_extent(renderer=renderer).width
        probe.remove()
        return w

    lines, current = [], ""
    for word in text.split(" "):
        trial = f"{current} {word}".strip()
        if current and width(trial) > limit:
            lines.append(current)
            current = word
        else:
            current = trial
    if current:
        lines.append(current)
    return "\n".join(lines)


def panel(w_mm: float, h_mm: float):
    """Return (fig, ax) sized in millimetres."""
    fig, ax = plt.subplots(figsize=(w_mm * MM, h_mm * MM))
    return fig, ax


def despine(ax, keep=("left", "bottom")) -> None:
    for side in ("top", "right", "left", "bottom"):
        ax.spines[side].set_visible(side in keep)


def save(fig, stem: Path | str, *, exact: bool = False,
         formats: tuple[str, ...] = ("pdf", "png")) -> None:
    """Save a panel as vector PDF (deliverable) + 600 dpi PNG (preview).

    Line art and text stay vector (resolution-independent); dpi=600 only governs
    any element explicitly marked ``rasterized=True`` so it also meets the
    >=600 dpi Nature line-art bar rather than matplotlib's ~100 dpi default.

    ``exact=True`` writes the declared physical canvas instead of a tight
    bounding box. The Fig. 1 compositor places every panel at the width it was
    drawn at so its internal 5-7 pt type is 5-7 pt on the page; a tight box
    silently changes that width by however much whitespace matplotlib trims, so
    panels destined for the compositor must be saved exactly.

    ``formats`` selects the outputs. The default keeps the pdf+png pair every
    figure has used; Fig. 1 additionally requests svg (editable) and tiff
    (600 dpi raster fallback) for the production bundle.
    """
    stem = Path(stem)
    kwargs = {}
    if exact:
        w_in, h_in = fig.get_size_inches()
        kwargs = {"bbox_inches": Bbox.from_bounds(0, 0, w_in, h_in),
                  "pad_inches": 0}
    for suffix in formats:
        pil = {"compression": "tiff_lzw"} if suffix == "tiff" else None
        fig.savefig(stem.with_suffix(f".{suffix}"), dpi=600,
                    **({"pil_kwargs": pil} if pil else {}), **kwargs)
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


# Where a panel's named source table is allowed to live, searched in order.
# Panels must not carry private copies of committed tables: the Fig. 2 to Fig. 5
# refinement bundles each shipped their own `data/` directory, and a panel that
# reads its own copy will keep plotting the old numbers after the committed
# table is regenerated. Resolving by name against the repository is what keeps
# figure and result in step.
DATA_SEARCH = ("data", "results", "results/canonical", "results/benchmark_resolution")

# Names a panel asks for that do not match the committed filename.
DATA_ALIASES = {
    "benchmark_resolution_summary.csv": "results/benchmark_resolution/summary.csv",
}


def data_path(name: str) -> Path:
    """Resolve a source-table filename to its committed path in the repository.

    Raises with the full search list rather than falling back to a local copy,
    because a silent fallback is how a figure drifts away from the result it
    claims to show.
    """
    root = repo_root()
    if name in DATA_ALIASES:
        path = root / DATA_ALIASES[name]
        if not path.exists():
            raise FileNotFoundError(f"{name} is aliased to {path}, which does not exist")
        return path
    for folder in DATA_SEARCH:
        path = root / folder / name
        if path.exists():
            return path
    raise FileNotFoundError(
        f"{name} not found under any of "
        + ", ".join(str(root / f) for f in DATA_SEARCH)
        + ". Commit the table to the repository rather than shipping a copy "
          "beside the panel.")


def load_bench(exclude_wt: bool = True, *, corrected: bool = False) -> "pd.DataFrame":
    """Load the committed benchmark table.

    The CSV has 472 rows, of which 2 are wild-type reference rows
    (KRAS WT, 644 cells; GATA1 WT, 38,276 cells). ``exclude_wt=True`` drops
    them so per-variant panels count the 470 real protein-coding variants and
    the GATA1 WT pool does not appear as a spurious high-depth "variant".
    Set ``exclude_wt=False`` to reproduce the manuscript's current
    WT-inclusive counts (472 / 98 / 93 / 255 / 26).

    ``corrected=True`` reads ``allele_perturb_bench_v2.csv``, the de-leaked
    table whose ``is_hotspot`` comes from external annotation only (COSMIC/IARC
    codons, KRAS activating codons, GATA1 zinc-coordinating cysteines plus
    ClinVar, JAK1 JH1/JH2). It differs from the default table on 182 of the 470
    rows, and the difference is not cosmetic: the default ``is_hotspot`` is a
    leaked continuous score that is not even an indicator (it sums to -6.16 over
    TP53 and +11.27 over KRAS), preserved in v2 as ``is_hotspot_leaked_OLD``.
    The evaluation grid the manuscript reports, ``results/results_v4_exttheta.csv``,
    is computed on the de-leaked theta, so any panel that draws theta or a
    hotspot annotation must pass ``corrected=True`` to agree with it.
    """
    name = "allele_perturb_bench_v2.csv" if corrected else "allele_perturb_bench.csv"
    path = repo_root() / "data" / name
    if not path.exists():
        raise FileNotFoundError(f"benchmark table not found: {path}")
    df = pd.read_csv(path)
    if exclude_wt:
        df = df[df.variant.astype(str).str.upper() != "WT"].copy()
    return df
