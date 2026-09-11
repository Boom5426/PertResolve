# Project structure

PertResolve has two public layers: an installable measurement-resolution
tool and the benchmark/reproducibility materials for the paper
“Measurement resolution constrains fine-grained perturbation prediction”.

```
PertResolve/
├── pertresolve/                    # installable Python package
│   ├── resolution/                 # independent detection, identification, reproducibility
│   │                                # reference and model-ranking resolution
│   ├── bench.py                    # PertResolve-Bench metadata loader
│   ├── features.py                 # molecular feature construction
│   ├── metrics.py                  # evaluation convenience API
│   ├── evaluation/                 # PDS, direction and differential-expression metrics
│   └── paths.py                    # explicit external-data and output-path handling
├── data/                           # small benchmark metadata tables
│   └── demo/                       # self-contained resolution example
├── results/                        # derived tables behind the paper
│   ├── canonical/                  # primary figure and manuscript inputs
│   ├── benchmark_resolution/       # per-resource resolution summaries
│   ├── candidate_and_power/        # candidate-pool and benchmark-size analyses
│   ├── pilot_validation/            # pilot-to-evaluation validation
│   └── resolution_panel_v2/         # broader PertResolve-Bench panel
├── scripts/analysis/               # table generators and shared evaluation harness
├── examples/                       # end-to-end runnable examples
├── manuscript/                     # final PDFs and Supplementary Data 1 workbook
├── assets/                         # static README visual assets
├── tests/                          # estimator and repository-contract tests
├── pyproject.toml                  # package metadata and optional dependencies
└── requirements.txt                # convenience install for older tooling
```

## Public package

The core API needs only NumPy and pandas:

```python
from pertresolve.resolution import resolution_report

report = resolution_report(X, labels, control="non-targeting", depth=50)
print(report.summary())
```

For AnnData inputs, install the `io` extra and use
`pertresolve-resolution`. The report separates perturbation detection,
nearest-competitor identification, an empirical split-half reproducibility
reference, local response geometry and model-ranking resolution. These axes are
reported independently; the split-half reference is not a hard bound, and it is
not used alone to certify model-ranking suitability.

## Reproducing the paper

```bash
pip install -e ".[dev,bench]"
pytest
```

The public `manuscript/` directory intentionally contains only the final
manuscript PDF, Supplementary Information PDF and Supplementary Data 1
workbook. Figure-generation sources and LaTeX sources are kept in the local
staging backup, not in the public checkout.

Expression matrices, large atlas files and protein-language-model embeddings
are not redistributed. Scripts that need them require an explicit `--base` or
`PERTRESOLVE_DATA` path; atlas analyses use `--atlas-dir` or
`PERTRESOLVE_ATLAS_DIR`. Analysis outputs must be written to scratch locations,
not directly over the committed results.
