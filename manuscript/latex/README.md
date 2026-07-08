# LaTeX Compilation

## Quick build

```bash
cd /home/boom/ICLR/AllelePerturb/manuscript/latex
make
```

This runs `pdflatex → bibtex → pdflatex × 2` and produces `AllelePerturb_manuscript.pdf`.

## Manual build (step by step)

```bash
pdflatex -interaction=nonstopmode AllelePerturb_manuscript
bibtex AllelePerturb_manuscript
pdflatex -interaction=nonstopmode AllelePerturb_manuscript
pdflatex -interaction=nonstopmode AllelePerturb_manuscript
```

The first `pdflatex` generates `.aux` with citation keys; `bibtex` resolves them from `references.bib` into `.bbl`; the final two passes resolve cross-references and numbering.

## Requirements

- pdflatex (TinyTeX or TeX Live)
- Required packages: `mathpazo`, `natbib`, `titlesec`, `enumitem`, `lineno`, `float`, `hyperref`
- On TinyTeX, install missing packages with: `tlmgr install mathpazo palatino psnfss fpl titlesec enumitem lineno natbib`

## Files

| File | Description |
|------|-------------|
| `AllelePerturb_manuscript.tex` | Main LaTeX source |
| `references.bib` | BibTeX bibliography (19 entries, all DOI-verified) |
| `figures/fig1–5.pdf` | Main figures (embedded via `\includegraphics`) |
| `figures/ED_fig1.pdf` | Extended Data Figure 1 |
| `Makefile` | Build automation |

## Clean intermediate files

```bash
make clean
```
