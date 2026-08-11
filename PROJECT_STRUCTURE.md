# Project structure

Two things live here, and the layout follows that split.

**A tool.** `alleleperturb.resolution` measures whether a perturbation dataset can arbitrate
a model comparison, before any model is trained. It depends on numpy and pandas alone and is
usable on any dataset, not only the ones studied here.

**A study.** The allele-resolution benchmark that motivated the tool, its manuscript, and
every table behind every number in it.

```
AllelePerturb/
├── alleleperturb/                  # the installable package
│   ├── resolution/                 # measurement-resolution diagnostics (numpy + pandas only)
│   │   ├── profiles.py             # cut cells into disjoint measurement groups
│   │   ├── window.py               # detection: each perturbation against its control
│   │   ├── scaling.py              # sampling noise, between-perturbation signal, their ratio
│   │   ├── recovery.py             # whether a known predictor ordering survives the benchmark
│   │   ├── report.py               # the three levels together, with a verdict
│   │   └── cli.py                  # alleleperturb-resolution
│   ├── bench.py                    # benchmark loader
│   ├── features.py                 # theta feature construction
│   ├── metrics.py                  # top-level metric API
│   ├── paths.py                    # external-location resolution; refuses to guess
│   └── evaluation/                 # pds.py (incl. residual_pds_score), de_metrics, direction_metrics
├── tests/                          # simulations whose answer is known by construction
│   ├── test_resolution_scaling.py  # estimator: unbiasedness, calibration, coverage
│   └── test_resolution_report.py   # verdict: the ladder, the null case, saturation
├── data/                           # benchmark tables
│   ├── allele_perturb_bench.csv        # 470 variants (plus 2 WT rows): theta, n_cells, splits
│   ├── allele_perturb_bench_v2.csv     # external-hotspot definition; used by Fig. 1
│   └── hotspot_external_definition.txt # COSMIC / IARC / ClinVar criteria
├── results/
│   ├── canonical/                  # every table behind a manuscript number, each with a generator
│   ├── benchmark_resolution/       # per-dataset resolution JSON + its scripts
│   ├── pilot_validation/           # pilot-to-full validation
│   └── reviewer_controls/          # pre-submission controls
├── scripts/
│   ├── analysis/                   # generators for the canonical tables
│   └── figures/                    # panel data preparation and the height-budget measurer
├── manuscript/
│   ├── figures/                    # one script per panel, one composite per figure
│   │   ├── nm_style.py             # shared style; resolves the sans face by glyph coverage
│   │   ├── check_panels.py         # the gate every panel must pass
│   │   ├── compose.py              # builds a composite from a row spec
│   │   └── fig1/ ... fig6/         # figN<letter>_<panel>.py + figN_assemble.tex
│   └── latex/                      # AllelePerturb_manuscript.tex, AllelePerturb_SI.tex
├── docs/                           # pre-registrations, result records, audits
├── pyproject.toml                  # package metadata and dependency extras
└── requirements.txt                # defers to pyproject; kept only for older tooling
```

## Figures

Each figure is built one panel at a time. `manuscript/figures/figN/figN<letter>_<panel>.py`
writes a vector PDF, and `figN_assemble.tex` places them from a row specification produced by
`compose.py`. Nothing is scaled at placement: type is set in points on a canvas declared in
millimetres, so a panel is drawn at the size it will be printed.

`check_panels.py` is the gate. It checks canvas overflow with a 0.4 mm margin, the declared
millimetre size, a 5 pt type floor, row width against the 183 mm double-column width, the
per-figure stacked height budget, a single embedded font family, and label collisions.

Height budgets are per figure and are **not** interchangeable, because a float carries its
caption and each caption is a different length. `scripts/figures/measure_figure_budget.py`
measures them by bisecting against LaTeX's own float verdict. Lengthening a caption shrinks
that figure's budget, so captions and panels have to be changed together.

The printed figure numbers and the build directories cross over for two figures: printed
Figure 3 is built in `fig4/` and printed Figure 4 in `fig3/`. The composite filenames follow
the directories, not the printed numbers.

## Reproducing

```bash
pip install -e ".[dev]" && pytest          # the tool, checked against known answers
cd manuscript/latex && make                # the manuscript and the Supplementary Information
python manuscript/figures/check_panels.py  # every panel, against every constraint
```

Two analyses read only committed tables and run on a fresh checkout:

```bash
python results/pilot_validation/pilot_validate.py --out /tmp/ap_out
python scripts/figures/rankability_predictor.py   --out /tmp/ap_out/rankability.csv
```

Anything that touches single-cell data needs the processed arrays, which are not shipped.
Point `ALLELEPERTURB_DATA` at them, or pass `--base`. The location is never inferred, and
`--out` is refused if it resolves inside `results/`.

## Canonical numbers

The first four lines below are frozen in `results/canonical/canonical_numbers.json` and are
read from it rather than restated by hand.

- 470 protein-coding variants, 321,043 cells
- `D_self / D_null` at native depth: TP53 0.965, KRAS 1.004, GATA1 0.878, JAK1 0.21
- un-rankable at native depth: TP53 100%, KRAS 100%, GATA1 97.6%, JAK1 10%
- gene-level atlases, un-rankable at native depth: Replogle 55.3%, Adamson 14.6%, Norman 3.4%

The two model-score ranges come from their own tables, `definitive_summary.csv` and
`results/pearson_delta_bootstrap_summary.csv`:

- no method exceeds chance PDS (0.50); PDS 0.49 to 0.52, Pearson-delta 0.55 to 0.65

Later results have their own records under `docs/`: the refutation of the one-variable
scaling claim, the quantity that governs discrimination and its depth calculator, the
residual scoring axis and its 0.524 chance level, and the pre-registered six-screen panel.
