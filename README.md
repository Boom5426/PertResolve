<div align="center">

# AllelePerturb

**A Benchmark for Allele-Resolution Single-Cell Perturbation Prediction**

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/)

</div>

## Overview

AllelePerturb is a benchmark and evaluation framework for predicting single-cell transcriptional responses at **protein-coding variant resolution**. Unlike gene-level perturbation benchmarks, AllelePerturb asks whether models can distinguish the cellular effects of different missense mutations within the same gene (e.g., TP53 R175H vs R273C).

**Key finding:** Systematic evaluation of 20 feature-model combinations reveals that all methods recover perturbation direction (Pearson-δ 0.60–0.68) but none reliably ranks held-out variants (PDS 0.49–0.52, statistically indistinguishable from chance 0.50). Split-half analysis traces this to a measurement-resolution floor governed by effect size, sampling noise and sequencing depth.

### Benchmark at a Glance

| Property | Value |
|----------|-------|
| **Genes** | TP53, KRAS, GATA1, JAK1 |
| **Variants** | 472 protein-coding variants |
| **Cells** | 321,043 single cells |
| **Technologies** | Perturb-seq, base editing, scSNV-seq |
| **Splits** | 6 generalization strategies |
| **Methods** | 20 feature-model combinations |
| **Metrics** | 10 evaluation metrics |

---

## Repository Structure

```
AllelePerturb/
├── alleleperturb/           # Python package (bench loading, evaluation)
│   ├── bench.py             # Load and filter the benchmark table
│   ├── evaluation/          # Evaluation utilities
│   └── features.py          # θ feature computation
├── data/                    # Benchmark tables (included in repo)
│   ├── allele_perturb_bench.csv          # 472 variants, θ₆ features, splits
│   ├── allele_perturb_bench_v2.csv       # V2 with external-only hotspot
│   ├── allele_perturb_bench_exttheta.csv # De-leaked θ version
│   └── hotspot_external_definition.txt   # External hotspot criteria
├── results/                 # Pre-computed result tables (included)
│   ├── results_v4_exttheta.csv           # Main grid: 20 methods × 5 splits × 4 genes × 10 metrics
│   ├── results_v4_10metrics.csv          # Full 10-metric grid
│   ├── second_probe_rankability_table.csv # Split-half rankability (55K rows)
│   ├── split_half_power_curve.csv        # Detection-rate vs depth
│   ├── rankability_predictor_honest.csv   # LODO rankability predictor (honest, per-perturbation)
│   ├── bootstrap_CIs.json               # Bootstrap confidence intervals
│   ├── canonical_numbers.json            # Frozen canonical numbers for manuscript
│   ├── permutation_null_pds.csv          # 1000× permutation null for PDS
│   ├── rankability_sensitivity.csv       # Sensitivity to criterion choice
│   └── unrankable_canonical.json         # Un-rankable fractions per dataset
├── scripts/
│   ├── figures/             # Figure-drawing scripts (self-contained)
│   │   ├── fig_config.py    # Shared config, colors, loaders
│   │   ├── draw_fig2.py     # Direction-ranking dissociation
│   │   ├── draw_fig3.py     # Split-half measurement window
│   │   ├── draw_fig4.py     # Robustness across splits/metrics
│   │   └── draw_fig5.py     # Rankability prediction + workflow
│   ├── baselines/           # Baseline method implementations
│   └── run_all_splits.py    # Full evaluation grid runner
├── figures/                 # Generated figures organized by panel
│   ├── fig1/ ... fig5/      # Each: composite + panels/ + figure.md
│   └── extended_data/       # Extended Data figures
├── manuscript/
│   ├── AllelePerturb_manuscript_en.md    # Markdown manuscript
│   └── latex/               # LaTeX build directory
│       ├── AllelePerturb_manuscript.tex
│       ├── references.bib   # 19 BibTeX entries (DOI-verified)
│       ├── figures/          # PDF figures for embedding
│       ├── Makefile          # Build: `make`
│       └── README.md         # Compilation instructions
├── PROJECT_STRUCTURE.md     # Detailed structure documentation
├── requirements.txt         # Python dependencies
└── README.md                # This file
```

---

## Quick Start

### 1. Install Dependencies

```bash
# Clone the repo
git clone https://github.com/Boom5426/AllelePerturb.git
cd AllelePerturb

# Create environment (Python 3.11+)
conda create -n alleleperturb python=3.11
conda activate alleleperturb
pip install -r requirements.txt
```

**Requirements:** numpy, pandas, scipy, scikit-learn, matplotlib, seaborn, scanpy, fair-esm, torch

### 2. Reproduce Figures (from pre-computed results)

All figure-drawing scripts read from `results/` (included in the repo) — no raw data download needed.

```bash
cd scripts/figures

# Figure 2: Direction-ranking dissociation
python draw_fig2.py

# Figure 3: Split-half measurement window
python draw_fig3.py

# Figure 4: Robustness across splits and metrics
python draw_fig4.py

# Figure 5: Rankability prediction and workflow
python draw_fig5.py
```

Output composites are saved to `../../figures/composites/`.

### 3. Reproduce the Evaluation Grid (from raw data)

This requires downloading the raw single-cell data (~2 GB total). See [Data](#data) below.

```bash
# After downloading raw arrays to data/:
python scripts/run_all_splits.py \
    --data_dir data/ \
    --output results/results_v4_10metrics.csv
```

This runs 20 methods × 5 splits × 4 genes × 10 metrics. Runtime: ~30 min on a 24-core CPU.

---

## Data

### Included in the repository

| File | Rows | Description |
|------|------|-------------|
| `data/allele_perturb_bench.csv` | 472 | Variant metadata: gene, protein, variant name, cell count, θ₆ features, 6 split assignments |
| `data/allele_perturb_bench_v2.csv` | 472 | V2 with external-only hotspot (de-leaked) |
| `results/results_v4_exttheta.csv` | 10,276 | Full evaluation grid (canonical, external-θ) |
| `results/second_probe_rankability_table.csv` | 55,548 | Split-half rankability across depth bins |

### External data (download separately)

The raw single-cell expression arrays are too large for git (~2 GB). They are needed only to **re-run the evaluation grid** (`scripts/run_all_splits.py`). Drawing figures from pre-computed results does **not** require them.

#### TP53 + KRAS (Ursu et al. 2022, A549 Perturb-seq)

Source: [GEO GSE161824](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE161824)

```bash
cd data/
# Download processed matrix files
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE161nnn/GSE161824/suppl/GSE161824_A549_TP53.processed.matrix.mtx.gz
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE161nnn/GSE161824/suppl/GSE161824_A549_TP53.processed.genes.csv.gz
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE161nnn/GSE161824/suppl/GSE161824_A549_TP53.variants2cell.csv.gz
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE161nnn/GSE161824/suppl/GSE161824_A549_KRAS.processed.matrix.mtx.gz
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE161nnn/GSE161824/suppl/GSE161824_A549_KRAS.processed.genes.csv.gz
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE161nnn/GSE161824/suppl/GSE161824_A549_KRAS.variants2cell.csv.gz
```

After downloading, run `scripts/preprocess_ursu.py` to produce `joint_arrays.npz`.

#### GATA1 (Yu & Welch 2025, HSPC base editing)

Source: [HuggingFace cyclopeta/PerturbNet_reproduce](https://huggingface.co/datasets/cyclopeta/PerturbNet_reproduce)

```bash
pip install huggingface_hub
python -c "
from huggingface_hub import hf_hub_download
hf_hub_download('cyclopeta/PerturbNet_reproduce',
    'GATA1_standard_hvg_pert_filtered.h5ad',
    local_dir='data/')
"
```

After downloading, run `scripts/preprocess_gata1.py` to produce `gata1_arrays.npz`.

#### JAK1 (Cooper et al. 2024, HT-29 scSNV-seq)

Source: [Zenodo 10418435](https://doi.org/10.5281/zenodo.10418435) (ENA: PRJEB48915)

```bash
# Download the SingleCellExperiment RDS (232 MB)
wget https://zenodo.org/records/10418435/files/scSNPseq_data.zip
unzip scSNPseq_data.zip -d data/jak1/
```

After downloading, run `scripts/preprocess_jak1.py` (requires R with scran + SingleCellExperiment) to produce `jak1_arrays.npz`.

#### ESM-1v Embeddings

```bash
# Extract ESM-1v embeddings for all 472 variants
python scripts/extract_esm_embeddings.py \
    --bench data/allele_perturb_bench.csv \
    --output data/esm1v_embeddings.npz
```

Requires `fair-esm` and ~2 GB GPU memory (or ~10 min on CPU).

---

## Evaluation Protocol (AllelePerturb-Eval)

### Core Metrics

| Metric | What it measures | Range |
|--------|-----------------|-------|
| **PDS** (Perturbation Discrimination Score) | Is a prediction closer to its own target than to other variants? | 0–1 (0.5 = chance) |
| **Pearson-δ** | Correlation between predicted and true perturbation directions | -1 to 1 |
| **DE overlap** | Fraction of top-50 DE genes shared between prediction and truth | 0–1 |
| **Direction agreement** | Sign concordance across genes | 0–1 |
| **DE-LFC Spearman** | Rank correlation of log-fold-changes for DE genes | -1 to 1 |
| **MAE** | Mean absolute error of predicted profiles | ≥ 0 |

PDS is computed under three distances (cosine, L1, L2); PDS-cosine is the primary metric reported.

### Generalization Splits

| Split | Training | Test | Tests |
|-------|----------|------|-------|
| Random | 65% variants | 35% variants | Standard generalization |
| OOD-Position | N-terminal half | C-terminal half | Positional extrapolation |
| OOD-Mechanism | Non-hotspot | Hotspot/ZF variants | Mechanistic extrapolation |
| Cross-Gene | 3 genes | 1 gene (all variants) | Gene transfer |
| Low-N | Variants ≥200 cells | Variants <200 cells | Low-depth evaluation |
| PerturbNet-Compat | Matching PerturbNet holdouts | PerturbNet test sets | Direct comparison |

---

## Manuscript

### Compile the LaTeX manuscript

```bash
# Requires TeX Live (install with: sudo apt install texlive-full)
cd manuscript/latex
make
```

This produces `AllelePerturb_manuscript.pdf` (20 pages, 5 embedded figures, 19 BibTeX references).

### Manuscript files

| File | Description |
|------|-------------|
| `manuscript/AllelePerturb_manuscript_en.md` | Markdown source |
| `manuscript/latex/AllelePerturb_manuscript.tex` | LaTeX source |
| `manuscript/latex/references.bib` | Bibliography (19 DOI-verified entries) |
| `manuscript/latex/figures/fig1–5.pdf` | Embedded figure PDFs |

---

## Citation

If you use AllelePerturb in your research, please cite:

```bibtex
@article{li2026alleleperturb,
  title={AllelePerturb reveals measurement limits of protein-coding variant 
         prediction in single-cell transcriptomics},
  author={Li, Bo},
  year={2026},
  note={Manuscript in preparation}
}
```

## License

This project is licensed under the MIT License.

## Contact

Bo Li — University of Florida, Department of Biomedical Engineering
