# AllelePerturb — Project Instructions

Read together with the global `~/.claude/CLAUDE.md`. This file adds only
project-specific facts and red lines; it does not repeat the general principles.

## What this repo is, and where the compute lives

This local repo is the **lean, publishable artifact**: the packaged benchmark
code (`alleleperturb/`), figure/analysis scripts, and the manuscript
(`manuscript/latex/`). It is the writing + release side.

The **heavy compute (model training, baselines, data) lives on a remote server**.
Claude Code always runs **locally** and reaches the server over `ssh`/`rsync`;
the server has no Claude Code process and needs no CLAUDE.md of its own. When a
task needs compute or the large data, run it remotely via ssh; do not try to
reproduce it locally.

### Remote server facts (confirmed 2026-07)

- SSH: host alias `139.180.131.202` (`User bob`, `Port 2222`,
  `IdentityFile ~/.ssh/compute_key`). Just `ssh 139.180.131.202`.
- `/data/boom` is a symlink to `/data/home_boom`. Two VCCompass dirs, distinct roles:
  - `/data/boom/NUS/VCCompass`: **active compute workspace** with run scripts
    (`all_splits_v*.py`, `cellflow_run*.py`, `cpa_run*.py`, `cfm_decoder*.py`),
    logs, versioned iterations. Messy by nature; v2/v3/v4 coexist.
  - `/data/boom/VCData/VCCompass`: **data + model artifacts**, i.e. large `.npz`
    arrays, `.pt` checkpoints, `raw/`, `results/`, `summary.json`.
- **Raw data is immutable**: `/data/boom/VCData/VCCompass/raw/` holds the
  GSE161824 A549 KRAS/TP53 processed matrices. Never modify, re-download over, or
  write derived outputs into `raw/`.
- Env: conda env **`Agent`** (`/home/bob/anaconda3/envs/Agent`, python 3.11.14)
  has the numpy/pandas/scipy/sklearn stack **and** torch 2.10+cu126 with the
  RTX 4090 (24 GB) visible. Use it for both the sklearn benchmark grid and the
  torch methods. (`NUS/VCCompass` is a loose working dir, not a git repo.)
- Long jobs: launch under `tmux` or `nohup` so an ssh drop doesn't kill them; the
  workspace already keeps per-run `*.log` files, so follow that convention.
- Note: the server also has `~/.claude`, `~/.codex`, `AGENTS.md` from other agent
  tooling. Unverified whether Claude Code is ever launched *on* the server; if the
  workflow changes to that, this repo's remote assumptions must be revisited.

## Scientific red lines (specific to this benchmark)

- **The null result is real and must stay honest.** Per
  `VCData/VCCompass/summary.json`, `VCCompass_6dim` (~0.704) is statistically
  indistinguishable from `Random_theta` (~0.699), `Dosage_only`, and
  `Biophys_only_3dim`. Do **not** reframe, threshold-tune, or select splits to
  make the 6-dim model look like it beats baselines. It does not, and that is the
  paper's finding.
- **Split unit is variant-level.** Test variants must never influence fitting,
  feature selection, normalization, thresholds, or hyperparameter tuning.
- **Metric definitions are fixed.** PDS (and its PDS_cos/PDS_L1/PDS_L2 variants),
  direction metrics (pearson_delta, delta_cosine), DE metrics, and any
  D_self/D_null quantities keep their current definitions. Do not swap in a proxy
  metric or "harmonize" definitions to shift results without explicit sign-off.
- **Report per-gene before pooled.** Pooling inflates apparent performance: the
  `results/debt2_pool_inflation.csv` / `debt2_permutation_null.csv` /
  `debt1_no_leakage.csv` analyses exist precisely to guard this. Preserve that logic.
- Every main numerical claim in the manuscript must trace to a generated result
  file (on the server under `VCData/VCCompass/` or synced into local `results/`).

## Local tooling reality (do not fake a quality gate)

There is currently **no `tests/`, no `pyproject.toml`, no ruff/mypy/pytest/
pre-commit** configured here; `requirements.txt` lists only numpy/pandas/scipy/
sklearn/pyyaml. So do **not** claim to have run linting/type-checking/tests that
aren't set up. If a change warrants real gating, propose adding the tooling as a
separate, explicit step rather than pretending it exists.

## Prohibited

- Do not modify the raw GSE161824 matrices, or write derived files into any `raw/`.
- Do not commit large artifacts to git: `.npz`, `.pt` checkpoints, raw matrices,
  or the manuscript PDFs' intermediates. (`.gitignore` already excludes LaTeX
  aux/bbl/etc.)
- Do not change an evaluation definition merely to improve a method's numbers.
- Do not edit a manuscript claim before confirming the corresponding result file.
- Do not silently drop failed variants/cells or exclude splits.
