<div align="center">

# AllelePerturb

**A Benchmark for Allele-Resolution Single-Cell Perturbation Prediction**

[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

</div>

---

## Overview

**AllelePerturb** is a benchmark and evaluation framework for **protein-coding variant-level** single-cell perturbation prediction. Unlike conventional gene-level benchmarks that treat "TP53 perturbed" as a single condition, AllelePerturb requires models to distinguish the transcriptional consequences of individual variants—e.g., TP53 R175H vs. R273C vs. R248Q.

### Key Numbers

| | |
|---|---|
| **Genes** | TP53, KRAS, GATA1, JAK1 |
| **Variants** | 472 protein-coding variants |
| **Cells** | 321,043 single cells |
| **Technologies** | Perturb-seq, base editing, scSNV-seq |
| **Metrics** | 10 (3 ranking + 3 direction + 3 DE + 1 reconstruction) |
| **Splits** | 6 generalization strategies |
| **Methods evaluated** | 30 (25 scored + 5 interface-incompatible) |

### Main Finding

> All 20 feature-conditioned methods achieve variant-ranking accuracy (PDS) within the chance band (0.43–0.48), yet maintain meaningful direction alignment (Pearson-Δ 0.60–0.68). This **direction–ranking dissociation** is the quantitative signature of operating at the evaluation resolution floor: allele-resolution prediction is currently **measurement-limited, not model-limited**.

---

## Installation

```bash
git clone https://github.com/Boom5426/AllelePerturb.git
cd AllelePerturb
pip install -e .
```

### Dependencies

```
numpy>=1.23
pandas>=1.5
scipy>=1.10
scikit-learn>=1.2
torch>=2.0        # optional, for VCCompass/PerturbNet baselines
scanpy>=1.9       # optional, for data preprocessing
```

---

## Repository Structure

```
AllelePerturb/
├── README.md
├── LICENSE
├── setup.py
├── requirements.txt
├── alleleperturb/                 # Core Python package
│   ├── __init__.py
│   ├── bench.py                   # Benchmark loader & split manager
│   ├── features.py                # θ₆ feature computation
│   ├── metrics.py                 # 10-metric evaluation suite
│   └── evaluation/
│       ├── __init__.py
│       ├── pds.py                 # PDS (cos/L1/L2) with tie-aware scoring
│       ├── de_metrics.py          # DE Overlap, LFC-Spearman, Direction Agreement
│       └── direction_metrics.py   # Pearson-Δ, Pearson-Δ̂20, delta_cosine, MAE
├── scripts/
│   ├── run_all_splits.py          # Full 20-method × 5-split × 4-gene grid
│   ├── baselines/
│   │   ├── feature_head_matrix.py # Ridge/Lasso/RF/GBoost/KNN/MLP × θ/ESM/ESM+θ
│   │   ├── latent_arithmetic.py   # scGen / scVIDR adaptations
│   │   └── biolord_runner.py      # Biolord (Tier-2)
│   ├── external_methods/
│   │   ├── perturbnet_eval.py     # PerturbNet (ESM-1v cINN)
│   │   ├── state_runner.py        # STATE (Arc Institute)
│   │   ├── cellflow_runner.py     # CellFlow (Theis lab)
│   │   ├── cpa_runner.py          # CPA
│   │   └── compatibility_report.py
│   └── probes/
│       ├── split_half_power.py    # D_self/D_null split-half analysis
│       ├── floor_audit.py         # E-distance un-rankable fraction
│       ├── rankability_predictor.py
│       └── metric_independence.py # Cross-metric/cross-space probe
├── data/
│   ├── README.md                  # Data download instructions
│   └── allele_perturb_bench.csv   # 472-variant benchmark table (θ + splits)
├── configs/
│   └── default.yaml               # Default evaluation config
├── results/
│   └── README.md
└── docs/
    ├── metrics.md                 # Metric definitions and rationale
    ├── splits.md                  # Split strategy documentation
    └── figures/
```

---

## Quick Start

### 1. Load the benchmark

```python
from alleleperturb import AllelePerturb

bench = AllelePerturb.load()
print(bench)
# AllelePerturb: 472 variants, 4 genes, 321,043 cells
# Splits: Random, OOD-Position, OOD-Mechanism, Cross-Gene, Low-N, PerturbNet-Compat

train, test = bench.split("split1")  # Random split
```

### 2. Evaluate a method

```python
from alleleperturb.metrics import evaluate

results = evaluate(
    predicted_deltas=pred_dict,   # {variant_name: np.array(n_genes)}
    ground_truth=bench,
    split="split1",
    gene="TP53",
)
print(results)
# {'PDS_cos': 0.46, 'PDS_L1': 0.45, 'PDS_L2': 0.46,
#  'pearson_delta': 0.68, 'pearson_delta_top20': 0.59,
#  'delta_cosine': 0.67, 'DE_overlap': 0.27,
#  'DE_LFC_spearman': 0.44, 'direction_agreement': 0.77,
#  'MAE': 0.076}
```

### 3. Run the full grid

```bash
python scripts/run_all_splits.py \
    --data_dir /path/to/data \
    --output results/v4_10metrics.csv
```

---

## Metrics

AllelePerturb-Eval uses **10 metrics** across three categories:

| Category | Metric | Measures | Chance level |
|----------|--------|----------|:------------:|
| **Ranking** | PDS_cos | Variant discrimination (cosine) | 0.500 |
| | PDS_L1 | Variant discrimination (L1) | 0.500 |
| | PDS_L2 | Variant discrimination (L2) | 0.500 |
| **Direction** | Pearson-Δ | Gene-level direction correlation | 0.000 |
| | Pearson-Δ̂20 | Direction on top-20 variable genes | 0.000 |
| | delta_cosine | Per-variant direction alignment | 0.000 |
| **DE fidelity** | DE Overlap | DEG recovery (top-50) | ~0.04 |
| | DE-LFC-Spearman | LFC rank correlation on sig. genes | 0.000 |
| | Direction Agreement | Fraction of DEGs with correct sign | 0.000 |
| **Reconstruction** | MAE | Mean absolute error of pseudobulk Δ | varies |

---

## Data

### Included in this repository
- `data/allele_perturb_bench.csv` — 472 variants with θ₆ features, split assignments, and metadata.

### External data (download separately)
| Gene | Source | Accession |
|------|--------|-----------|
| TP53 + KRAS | Ursu et al. 2022 | [GSE161824](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE161824) |
| GATA1 | PerturbNet (Yu & Welch 2025) | [HuggingFace](https://huggingface.co/datasets/cyclopeta/PerturbNet_reproduce) |
| JAK1 | Cooper et al. 2024 | [PRJEB48915](https://www.ebi.ac.uk/ena/browser/view/PRJEB48915) / [Zenodo](https://doi.org/10.5281/zenodo.10418435) |

### Feature representations
- **θ₆**: 6-dimensional biophysical feature vector (Δhydrophobicity, Δvolume, Δcharge, fold-core location, functional-switch residue, hotspot/pathogenic marker). External-only hotspot definition (no outcome leakage).
- **ESM-1v**: 1,280-dimensional mean-pooled embeddings from `esm1v_t33_650M_UR90S_1` (full dimension, no PCA reduction).
- **ESM+θ**: Concatenation of θ₆ and ESM-1v (1,286 dimensions).

See [`data/README.md`](data/README.md) for detailed download and preprocessing instructions.

---

## Splits

| Split | Design | Test variants | Purpose |
|-------|--------|:-------------:|---------|
| 1. Random | Original holdouts | 85 | Baseline generalization |
| 2. OOD-Position | C-terminal half → test | 230 | Spatial extrapolation |
| 3. OOD-Mechanism | Hotspot/ZF → test | 71 | Mechanistic extrapolation |
| 4. Cross-Gene | GATA1+JAK1 → test | 280 | Cross-gene transfer |
| 5. Low-N | <200 cells → test | 62 | Sparse data |
| 6. PerturbNet-Compat | PerturbNet holdout | 70 | Method comparison |

---

## Citation

```bibtex
@article{alleleperturb2026,
  title={AllelePerturb reveals detection limits in single-cell perturbation prediction},
  author={Zhang, Bob and Song, Qianqian},
  journal={In preparation},
  year={2026}
}
```

---

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.
