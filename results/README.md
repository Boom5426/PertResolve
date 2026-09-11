# Results

This directory contains the derived tables used by the final manuscript.
The `canonical/` tables are the primary provenance layer for the figures and
reported numbers; the other subdirectories hold supporting analyses for the
benchmark and measurement-resolution panel.

Most analyses that touch single-cell matrices require external processed data.
They write to an explicit scratch directory and refuse to overwrite the
committed `results/` tree. Use `--base` or `PERTRESOLVE_DATA` for the processed
allele-resolved inputs and `--atlas-dir` or `PERTRESOLVE_ATLAS_DIR` for public
perturbation atlases.

The repository includes the final manuscript PDFs and Supplementary Data 1
workbook under `manuscript/`.
