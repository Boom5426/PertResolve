"""Static QA for every manuscript figure's panel bundle.

One checker for all figures rather than one per directory: the Fig. 1 bundle
originally shipped with its own verifier that hard-coded each panel's size, and
the moment a panel was resized the verifier disagreed with the figure it was
supposed to be checking. The sizes live here once, next to the composite row
arithmetic that depends on them.

Five things are checked, each of which has actually gone wrong in this project:

1. **Canvas overflow.** Every Text artist must lie inside the declared canvas. A
   label that runs past the edge looks fine in the panel's own preview and is
   clipped, or overlaps its neighbour, once placed in the composite.
2. **Declared size.** Each panel PDF must be exactly the millimetre canvas the
   script declares. The compositors place panels at the width they were drawn at
   so 5-7 pt type stays 5-7 pt on the page; a tight bounding box silently
   changes that width.
3. **Type floor.** No text below 5 pt, the Nature Methods minimum.
4. **Row width.** Each composite row must fit 183 mm.
5. **One font family per panel.** The panel PDF must embed exactly one family.
   This catches the silent single-glyph fallback: mathtext hands any codepoint
   the chosen family lacks to STIXGeneral without warning, so one `\\cup` in a
   schematic is enough to put a second typeface in the figure.

Run from this directory:

    python check_panels.py           # every figure that has a panel directory
    python check_panels.py 2 3       # only those figures

Exits non-zero on the first failure; prints one PASS line per panel otherwise.
"""
from __future__ import annotations

import importlib
import os
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import nm_style as S  # noqa: E402

MIN_PT = 5.0
PAGE_MM = 183.0
ROW_GAP_MM = 2.6
# Per-figure height budget at 183 mm wide, measured by bisecting the LaTeX float
# check with scripts/figures/measure_figure_budget.py. These are not the
# journal's 247 mm maximum: the float carries its caption too, so a figure with
# a long caption has less room, which is why the five numbers differ. Keyed by
# PRINTED figure number. Re-measure whenever a caption changes length.
HEIGHT_BUDGET_MM = {1: 165.4, 2: 187.6, 3: 187.6, 4: 174.2, 5: 160.9, 6: 169.6}

# Printed figure number -> the file the manuscript includes. For figures 3 and 4
# these differ: manuscript/figures/README.md records that fig3/ and fig4/ are
# swapped relative to the printed numbers, and the budget above is measured per
# printed figure, so the two must be mapped explicitly rather than assumed equal.
PRINTED_OF_DIR = {"fig1": 1, "fig2": 2, "fig3": 4, "fig4": 3, "fig5": 5, "fig6": 6}
TOL_MM = 0.05
# Horizontal breathing room every label must keep from the canvas edge. Zero is
# not enough: matplotlib's text bounding box excludes an italic's slant
# overhang, so a label reported at 93.79 mm on a 94.0 mm canvas was still being
# cut. Only the horizontal edges are held to this; vertical placement is checked
# at the canvas edge itself.
MARGIN_MM = 0.4

# figure -> {"dir", "panels": {module: (stem, w_mm, h_mm)}, "rows": [(name, [modules], gutter)]}
# `rows` mirrors the coordinates in figN/figN_assemble.tex; the two must agree,
# and this file is the one that fails loudly when they do not.
SPECS: dict[int, dict] = {
    1: {
        "dir": "fig1",
        "panels": {
            "fig1a_nested_solvability": ("fig1a_nested_solvability", 113.5, 76.0),
            "fig1b_variant_tracks": ("fig1b_coverage", 66.0, 72.0),
            "fig1c_depth": ("fig1c_depth", 58.0, 42.0),
            "fig1d_theta_schematic": ("fig1d_theta_schematic", 69.0, 42.0),
            "fig1d_theta_pca": ("fig1d_theta_pca", 47.0, 42.0),
            "fig1e_eval_axes": ("fig1e_eval_axes", 78.0, 37.0),
            "fig1f_splits": ("fig1f_splits", 100.0, 37.0),
        },
        "rows": [
            ("1", ["fig1a_nested_solvability", "fig1b_variant_tracks"], 3.5),
            ("2", ["fig1c_depth", "fig1d_theta_schematic", "fig1d_theta_pca"], 4.5),
            ("3", ["fig1e_eval_axes", "fig1f_splits"], 5.0),
        ],
    },
    2: {
        "dir": "fig2",
        "panels": {
            "fig2a_evaluation": ("fig2a_evaluation", 59.0, 63.0),
            "fig2b_pds_forest": ("fig2b_pds_forest", 58.0, 63.0),
            "fig2c_pearson_forest": ("fig2c_pearson_forest", 44.0, 63.0),
            "fig2d_scatter": ("fig2d_scatter", 76.0, 44.0),
            "fig2e_de_gradient": ("fig2e_de_gradient", 89.5, 44.0),
            "fig2f_gene_gap": ("fig2f_gene_gap", 58.0, 44.0),
            "fig2g_pervariant": ("fig2g_pervariant", 107.5, 44.0),
        },
        # One gutter value, and one content width: every row comes to 170.0 mm
        # and is centred, so all three share the same left and right edge. Row 3
        # used to run flush to 183.0, which set panel f 6.5 mm left of a and d
        # and read as f overflowing the figure.
        "rows": [
            ("1", ["fig2a_evaluation", "fig2b_pds_forest", "fig2c_pearson_forest"], 4.5),
            ("2", ["fig2d_scatter", "fig2e_de_gradient"], 4.5),
            ("3", ["fig2f_gene_gap", "fig2g_pervariant"], 4.5),
        ],
    },
    3: {
        "dir": "fig3",
        "panels": {
            "fig3a_window_definition": ("fig3a_window_definition", 95.0, 41.0),
            "fig3b_representative_windows": ("fig3b_representative_windows", 83.0, 40.0),
            "fig3c_window_ratio": ("fig3c_window_ratio", 78.0, 43.0),
            "fig3d_sibling_identification": ("fig3d_sibling_identification", 100.0, 43.0),
            "fig3e_classifier_control": ("fig3e_classifier_control", 100.0, 43.0),
            "fig3f_effect_depth_landscape": ("fig3f_effect_depth_landscape", 78.0, 43.0),
            "fig3g_depth_rankability": ("fig3g_depth_rankability", 91.0, 36.0),
            "fig3h_native_rankability": ("fig3h_native_rankability", 87.0, 36.0),
        },
        "rows": [
            ("1", ["fig3a_window_definition", "fig3b_representative_windows"], 4.5),
            ("2", ["fig3c_window_ratio", "fig3d_sibling_identification"], 4.5),
            ("3", ["fig3e_classifier_control", "fig3f_effect_depth_landscape"], 4.5),
            ("4", ["fig3g_depth_rankability", "fig3h_native_rankability"], 4.5),
        ],
    },
    4: {
        "dir": "fig4",
        "panels": {
            "fig4a_pds_splits": ("fig4a_pds_splits", 89.0, 44.0),
            "fig4b_direction_splits": ("fig4b_direction_splits", 89.0, 44.0),
            "fig4c_empirical_null": ("fig4c_empirical_null", 105.0, 37.0),
            "fig4d_distance_invariance": ("fig4d_distance_invariance", 73.0, 37.0),
            "fig4e_feature_spaces": ("fig4e_feature_spaces", 89.0, 45.0),
            "fig4f_gene_split_dissociation": ("fig4f_gene_split_dissociation", 89.0, 45.0),
            "fig4g_model_harness": ("fig4g_model_harness", 94.0, 51.0),
            "fig4h_interface_audit": ("fig4h_interface_audit", 84.0, 51.0),
        },
        "rows": [
            ("1", ["fig4a_pds_splits", "fig4b_direction_splits"], 4.5),
            ("2", ["fig4c_empirical_null", "fig4d_distance_invariance"], 4.5),
            ("3", ["fig4e_feature_spaces", "fig4f_gene_split_dissociation"], 4.5),
            ("4", ["fig4g_model_harness", "fig4h_interface_audit"], 4.5),
        ],
    },
    5: {
        "dir": "fig5",
        "panels": {
            "fig5a_replicate_ceiling": ("fig5a_replicate_ceiling", 100.0, 42.0),
            "fig5b_ceiling_vs_models": ("fig5b_ceiling_vs_models", 78.0, 42.0),
            "fig5c_phase_map": ("fig5c_phase_map", 48.0, 31.0),
            "fig5d_predictor_family": ("fig5d_predictor_family", 130.0, 31.0),
            "fig5e_governing_quantity": ("fig5e_governing_quantity", 110.0, 34.0),
            "fig5f_external_benchmarks": ("fig5f_external_benchmarks", 68.0, 34.0),
            "fig5g_recovery_support": ("fig5g_recovery_support", 182.5, 43.0),
        },
        "rows": [
            ("1", ["fig5a_replicate_ceiling", "fig5b_ceiling_vs_models"], 4.5),
            ("2", ["fig5c_phase_map", "fig5d_predictor_family"], 4.5),
            ("3", ["fig5e_governing_quantity", "fig5f_external_benchmarks"], 4.5),
            ("4", ["fig5g_recovery_support"], 4.5),
        ],
    },
    6: {
        "dir": "fig6",
        "panels": {
            "fig6a_pilot_pipeline": ("fig6a_pilot_pipeline", 95, 41.0),
            "fig6b_retrospective_roc": ("fig6b_retrospective_roc", 65.5, 41.0),
            "fig6c_dataset_auroc": ("fig6c_dataset_auroc", 80, 43.0),
            "fig6d_prospective_validation": ("fig6d_prospective_validation", 80.5, 43.0),
            "fig6e_effect_depth_landscape": ("fig6e_effect_depth_landscape", 80, 37.0),
            "fig6f_paired_depthmatch": ("fig6f_paired_depthmatch", 80.5, 37.0),
            "fig6g_probability_workflow": ("fig6g_probability_workflow", 96.0, 38.0),
            "fig6h_screen_panel": ("fig6h_screen_panel", 64.5, 38.0),
        },
        "rows": [
            ("1", ["fig6a_pilot_pipeline", "fig6b_retrospective_roc"], 4.5),
            ("2", ["fig6c_dataset_auroc", "fig6d_prospective_validation"], 4.5),
            ("3", ["fig6e_effect_depth_landscape", "fig6f_paired_depthmatch"], 4.5),
            ("4", ["fig6g_probability_workflow", "fig6h_screen_panel"], 4.5),
        ],
    },
}


def capture(figure_dir: Path, module_name: str):
    """Run a panel's main() and return its Figure instead of writing files."""
    captured = {}
    real_save = S.save

    def fake_save(fig, stem, **kwargs):
        captured["fig"] = fig  # deliberately not closed

    S.save = fake_save
    sys.path.insert(0, str(figure_dir))
    try:
        module = importlib.import_module(module_name)
        importlib.reload(module)
        if getattr(module, "S", None) is not None:
            module.S.save = fake_save
        module.main()
    finally:
        S.save = real_save
        sys.path.remove(str(figure_dir))
        sys.modules.pop(module_name, None)
    if "fig" not in captured:
        raise AssertionError(f"{module_name}: main() did not call nm_style.save")
    return captured["fig"]


def out_of_view_tick_labels(fig) -> set[int]:
    """Identify tick labels whose tick lies outside its axis view interval.

    A tick locator can place a tick outside the limits; matplotlib keeps the Text artist,
    leaves ``get_visible()`` True and simply does not draw it. Measuring such an artist
    reports an overflow for a glyph that never reaches the page: a panel with the default y
    locator on limits of 0.46 to 1.03 carries a label for 1.2 sitting three millimetres above
    a 34 mm canvas, which the rendered PDF does not contain.

    This runs once, before anything is measured, because asking an axis for its tick labels
    refreshes its ticks. Doing that inside the measurement loop changes the figure while it
    is being measured, and turned an axis that had been switched off back on partway through,
    which flagged every tick of a schematic panel that draws no axis at all.
    """
    out = set()
    for ax in fig.axes:
        if not getattr(ax, "axison", True):
            continue          # a schematic that draws no axis has no tick to be out of view
        for which in ("x", "y"):
            axis = ax.xaxis if which == "x" else ax.yaxis
            lo, hi = ax.get_xlim() if which == "x" else ax.get_ylim()
            lo, hi = min(lo, hi), max(lo, hi)
            index = 0 if which == "x" else 1
            for label in axis.get_ticklabels():
                value = label.get_position()[index]
                if not (lo - 1e-9 <= value <= hi + 1e-9):
                    out.add(id(label))
    return out


def text_extent(text, renderer):
    """Extent of the glyphs alone, excluding an annotation's leader line.

    Annotation.get_window_extent() returns the union of the label and its
    arrow, so a leader running from a label down to its marker reports a box
    tens of millimetres tall. Two such labels then look like they overlap by
    80% when they are visibly separate. Measuring through the base class gives
    the text box itself.
    """
    if isinstance(text, matplotlib.text.Annotation):
        return matplotlib.text.Text.get_window_extent(text, renderer=renderer)
    return text.get_window_extent(renderer=renderer)


def embedded_families(pdf: Path) -> set[str]:
    """Font families embedded in `pdf`, subset tags and weights stripped."""
    blob = pdf.read_bytes()
    names = {m.decode().split("+")[-1]
             for m in re.findall(rb"/BaseFont\s*/([A-Za-z0-9+\-]+)", blob)}
    return {re.sub(r"-(Bold|Italic|BoldItalic|Regular|Oblique|BoldOblique)$", "", n)
            for n in names}


def check_panel(figure_dir: Path, module_name: str, stem: str,
                w_mm: float, h_mm: float) -> str:
    fig = capture(figure_dir, module_name)
    # Draw before measuring: legend children have no final position until the
    # legend is laid out, and their pre-layout extents all read as the origin.
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()

    got_w, got_h = (v * 25.4 for v in fig.get_size_inches())
    if abs(got_w - w_mm) > TOL_MM or abs(got_h - h_mm) > TOL_MM:
        raise AssertionError(
            f"{stem}: figure is {got_w:.2f} x {got_h:.2f} mm, declared "
            f"{w_mm:.2f} x {h_mm:.2f} mm")

    # Visibility is snapshotted from the pristine state first, because asking an axis for
    # its tick labels refreshes its ticks and can turn labels visible that matplotlib had
    # not been drawing. Reading visibility afterwards flagged every tick of a schematic
    # panel that switches its axis off.
    visible_before = {id(a) for a in fig.findobj(matplotlib.text.Text) if a.get_visible()}
    undrawn = out_of_view_tick_labels(fig)

    px_per_mm = fig.dpi / 25.4
    overflow, small = [], []
    for text in fig.findobj(matplotlib.text.Text):
        if not text.get_text().strip() or id(text) not in visible_before:
            continue
        if id(text) in undrawn:
            continue
        if text.get_fontsize() < MIN_PT - 1e-9:
            small.append(f"{text.get_fontsize():.2f} pt {text.get_text()[:40]!r}")
        bb = text_extent(text, renderer)
        x0, x1 = bb.x0 / px_per_mm, bb.x1 / px_per_mm
        y0, y1 = bb.y0 / px_per_mm, bb.y1 / px_per_mm
        if x0 < MARGIN_MM or x1 > w_mm - MARGIN_MM or y0 < -TOL_MM or y1 > h_mm + TOL_MM:
            overflow.append(
                f"[{x0:6.2f},{x1:6.2f}]x[{y0:6.2f},{y1:6.2f}] {text.get_text()[:40]!r}")
    # Overlapping labels. Canvas overflow alone does not catch a panel that has
    # been squeezed: the text stays inside the page and lands on top of its
    # neighbour instead. Only substantial overlaps are reported, because a label
    # deliberately set over a box or a marker is normal; two prose labels
    # covering half of each other is not.
    # 0.30 was too tight: adjacent labels overlap by ~1/3 through ascender and
    # descender padding alone, which flagged correctly-set panels.
    boxes = []
    for text in fig.findobj(matplotlib.text.Text):
        body = text.get_text().strip()
        # Two characters, not three: short labels like "P5" are exactly the
        # ones that stack in a compressed row list, and skipping them hid a
        # five-way collision in fig5e.
        if len(body) < 2 or id(text) not in visible_before or id(text) in undrawn:
            continue
        bb = text_extent(text, renderer)
        if bb.width > 0 and bb.height > 0:
            boxes.append((body, bb))
    collisions = []
    for i in range(len(boxes)):
        for j in range(i + 1, len(boxes)):
            (t1, b1), (t2, b2) = boxes[i], boxes[j]
            dx = min(b1.x1, b2.x1) - max(b1.x0, b2.x0)
            dy = min(b1.y1, b2.y1) - max(b1.y0, b2.y0)
            if dx <= 0 or dy <= 0:
                continue
            share = (dx * dy) / min(b1.width * b1.height, b2.width * b2.height)
            if share > 0.50:
                collisions.append(f"{share:5.0%} {t1[:28]!r} over {t2[:28]!r}")
    if small:
        raise AssertionError(f"{stem}: text below {MIN_PT} pt:\n  " + "\n  ".join(small))
    if collisions:
        raise AssertionError(f"{stem}: overlapping labels:\n  " + "\n  ".join(collisions))
    if overflow:
        raise AssertionError(
            f"{stem}: text outside the {w_mm} x {h_mm} mm canvas "
            f"(or within {MARGIN_MM} mm of a side edge):\n  "
            + "\n  ".join(overflow))

    pdf = figure_dir / f"{stem}.pdf"
    if not pdf.exists() or pdf.stat().st_size == 0:
        raise AssertionError(f"{stem}: {pdf.name} missing or empty; run the panel script")
    families = embedded_families(pdf)
    if len(families) != 1:
        raise AssertionError(
            f"{stem}: {len(families)} font families embedded ({', '.join(sorted(families))}). "
            "A single stray family is almost always mathtext falling back to "
            "STIXGeneral for one codepoint the chosen family lacks; replace that "
            "symbol rather than accepting the fallback.")
    return (f"PASS {stem:30s} {w_mm:6.2f} x {h_mm:5.2f} mm, no overflow, "
            f"type >= {MIN_PT} pt, {families.pop()}")


def check_figure(number: int) -> None:
    spec = SPECS[number]
    figure_dir = HERE / spec["dir"]
    print(f"--- Figure {number} ({spec['dir']})")
    for module_name, (stem, w_mm, h_mm) in spec["panels"].items():
        print(check_panel(figure_dir, module_name, stem, w_mm, h_mm))
    for name, members, gutter in spec["rows"]:
        widths = [spec["panels"][m][1] for m in members]
        total = sum(widths) + gutter * (len(members) - 1)
        if total > PAGE_MM + TOL_MM:
            raise AssertionError(
                f"figure {number} row {name}: {total:.2f} mm exceeds the "
                f"{PAGE_MM} mm page ({' + '.join(f'{w}' for w in widths)} "
                f"+ {len(members) - 1} x {gutter} mm gutter)")
        print(f"PASS row {name}: {total:.2f} mm of {PAGE_MM:.0f} mm "
              f"({len(members)} panels, {gutter} mm gutter)")
    printed = PRINTED_OF_DIR[spec["dir"]]
    budget = HEIGHT_BUDGET_MM[printed]
    stacked = (sum(max(spec["panels"][m][2] for m in members)
                   for _, members, _ in spec["rows"])
               + ROW_GAP_MM * (len(spec["rows"]) - 1))
    verdict = "PASS" if stacked <= budget else "FAIL"
    print(f"{verdict} {spec['dir']} height: {stacked:.1f} mm stacked, budget "
          f"{budget:.1f} mm (printed Figure {printed})")
    if verdict == "FAIL":
        raise AssertionError(
            f"{spec['dir']} (printed Figure {printed}) is {stacked:.1f} mm tall, "
            f"{stacked - budget:.1f} mm over its {budget:.1f} mm budget. Every panel "
            "in this figure already bottoms out at 5.0-5.6 pt type, so none can be "
            "scaled down at placement time; the panels have to be redrawn on smaller "
            "canvases. See docs/FIGURE_REDRAW_SPEC_2026-08-03.md for the target sizes.")


def main() -> None:
    os.chdir(HERE)
    wanted = [int(a) for a in sys.argv[1:]] or sorted(SPECS)
    for number in wanted:
        if number not in SPECS:
            raise SystemExit(f"no panel spec for figure {number}")
        check_figure(number)
    contract = S.write_font_contract()
    print(f"PASS bundle: figures {', '.join(map(str, wanted))}, font family "
          f"{S.SANS!r} (recorded in {contract.name} for the compositors)")


if __name__ == "__main__":
    main()
