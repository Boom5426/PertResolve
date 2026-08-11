"""Largest composite height each manuscript figure can carry, by bisection.

A figure that is legal at 183 mm wide can still be unplaceable: the float has to
carry its caption too, and the captions differ in length, so each figure has its
own height budget. This measures that budget instead of guessing it, by
substituting a blank 183 mm x H page for the real figure and bisecting on
whether LaTeX reports "Float too large".

Deliberately not computed from the overflow in points. The composite is scaled
to \textwidth on the page, so a point of overflow is not a millimetre of
composite height; bisecting on the actual verdict avoids the conversion.

The measured budgets are consumed by manuscript/figures/check_panels.py, which
fails a figure whose panels cannot be stacked inside them. Re-run this whenever
a caption changes length, since a longer caption shrinks its figure's budget.

Usage:
    python scripts/figures/measure_figure_budget.py

Takes a few minutes: each bisection step is a full manuscript build. The real
figure PDFs are restored after every probe.
"""
import os, re, shutil, subprocess, tempfile
from pathlib import Path

LATEX = Path(__file__).resolve().parents[2] / "manuscript" / "latex"
FIGS = LATEX / "figures"
# Scratch space for the blank probe PDFs. Overridable so a sandboxed run can put
# them somewhere writable; never inside the repo, since they are throwaway.
PROBE = Path(os.environ.get("ALLELEPERTURB_PROBE_DIR")
             or tempfile.mkdtemp(prefix="alleleperturb_budget_"))
PROBE.mkdir(parents=True, exist_ok=True)

src = (LATEX / "AllelePerturb_manuscript.tex").read_text().splitlines()
RANGES, start, which = {}, None, None
for i, line in enumerate(src, 1):
    if "\\begin{figure}" in line:
        start, which = i, None
    m = re.search(r"figures/(fig\d)\.pdf", line)
    if m:
        which = m.group(1)
    if "\\end{figure}" in line and start and which:
        RANGES[which] = (start, i + 1); start, which = None, None


def fits(fig, height_mm):
    tex = PROBE / "blank.tex"
    tex.write_text(
        "\\documentclass[tikz,border=1mm]{standalone}\n\\begin{document}\n"
        "\\begin{tikzpicture}[x=1mm,y=1mm]\n"
        f"  \\useasboundingbox (0,0) rectangle (183,{height_mm:.2f});\n"
        f"  \\fill[gray!15] (0,0) rectangle (183,{height_mm:.2f});\n"
        "\\end{tikzpicture}\n\\end{document}\n")
    subprocess.run(["lualatex", "-interaction=nonstopmode", "blank.tex"],
                   cwd=PROBE, check=True, capture_output=True)
    backup = PROBE / f"{fig}_backup.pdf"
    shutil.copy(FIGS / f"{fig}.pdf", backup)
    shutil.copy(PROBE / "blank.pdf", FIGS / f"{fig}.pdf")
    subprocess.run(["make", "AllelePerturb_manuscript.pdf"], cwd=LATEX,
                   check=False, capture_output=True)
    log = (LATEX / "AllelePerturb_manuscript.log").read_text()
    shutil.copy(backup, FIGS / f"{fig}.pdf")
    lo, hi = RANGES[fig]
    for m in re.finditer(r"Float too large for page by [\d.]+pt on input line (\d+)", log):
        if lo <= int(m.group(1)) <= hi:
            return False
    return True


for fig in ("fig1", "fig2", "fig3", "fig4", "fig5", "fig6"):
    lo, hi = 80.0, 260.0
    assert fits(fig, lo), f"{fig}: even {lo} mm does not fit"
    assert not fits(fig, hi), f"{fig}: even {hi} mm fits"
    while hi - lo > 0.5:
        mid = (lo + hi) / 2
        if fits(fig, mid):
            lo = mid
        else:
            hi = mid
    print(f"{fig}: max composite height {lo:6.1f} mm  (at 183 mm wide, caption as written)")
