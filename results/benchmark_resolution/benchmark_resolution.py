#!/usr/bin/env python3
"""Benchmark-resolution + reliability coefficient for a public perturbation dataset.

Self-contained: uses ONLY the dataset's own ground truth (no published model outputs).
For a dataset we compute, in PCA-50 space on pseudobulk deltas (pert - control):
  frac_rankable  : split-half fraction of perturbations with signal S=D_null-D_self > W (CI width of D_self)
  reliability    : P(a benchmark recovers the TRUE ranking of graded synthetic predictors),
                   via controlled predictors delta_hat(p,alpha)=(1-alpha)*global_mean+alpha*build_delta[p],
                   bootstrapped over perturbations (same construction as the allele-level 1C keystone)
  min_gap        : smallest alpha-gap the benchmark resolves (P>0.9 that higher-alpha scores higher)
Isolated exploratory output -> <out_dir>/<name>.json ; does not touch any manuscript file.

The five .h5ad atlases are public but too large to redistribute, so they are not
part of this repository. They originally lived under three different roots; this
script now expects them side by side under a single --atlas-dir, with their
original filenames unchanged:

  Replogle  <atlas-dir>/ReplogleWeissman2022_K562_essential.h5ad
  Norman    <atlas-dir>/NormanWeissman2019_filtered.h5ad
  Adamson   <atlas-dir>/AdamsonWeissman2016_GSM2406681_10X010.h5ad
  VCC       <atlas-dir>/adata_Training.h5ad
  sciPlex   <atlas-dir>/SrivatsanTrapnell2020_sciplex3.h5ad

Usage:
  python benchmark_resolution.py NAME --out OUT_DIR [--atlas-dir ATLAS_DIR]

    NAME                   one of: Adamson, Norman, Replogle, VCC, sciPlex
    --out OUT_DIR          required; directory receiving <NAME>.json
    --atlas-dir ATLAS_DIR  directory holding the five .h5ad atlases;
                           falls back to $ALLELEPERTURB_ATLAS_DIR
"""
import sys, json, numpy as np
import argparse
from pathlib import Path

# This script lives two levels below the repository root and is normally run
# directly, so make the repository importable without requiring an install.
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from alleleperturb.paths import resolve_base, require_inputs, ATLAS_ENV_VAR

DATASETS = ("Adamson", "Norman", "Replogle", "VCC", "sciPlex")

_parser = argparse.ArgumentParser(
    description="Benchmark-resolution + reliability coefficient for one public perturbation atlas.")
_parser.add_argument("name", choices=DATASETS, metavar="NAME",
                     help="dataset to analyse; one of: " + ", ".join(DATASETS))
_parser.add_argument("--out", required=True, metavar="OUT_DIR",
                     help="directory receiving <NAME>.json")
_parser.add_argument("--atlas-dir", default=None, metavar="ATLAS_DIR",
                     help="directory holding the five .h5ad atlases; "
                          f"falls back to ${ATLAS_ENV_VAR}")
_args = _parser.parse_args()
ATLAS_DIR = resolve_base(_args.atlas_dir, what="the perturbation atlas directory",
                         env_var=ATLAS_ENV_VAR, flag="--atlas-dir")
OUT_DIR = Path(_args.out)

import anndata as ad
from sklearn.decomposition import PCA

CFG = {
  'Replogle': dict(path=ATLAS_DIR / "ReplogleWeissman2022_K562_essential.h5ad", pcol=None, ctrl=None),
  'Norman':   dict(path=ATLAS_DIR / "NormanWeissman2019_filtered.h5ad",         pcol=None, ctrl=None),
  'Adamson':  dict(path=ATLAS_DIR / "AdamsonWeissman2016_GSM2406681_10X010.h5ad", pcol=None, ctrl=None),
  'VCC':      dict(path=ATLAS_DIR / "adata_Training.h5ad", pcol='target_gene', ctrl='non-targeting'),
  'sciPlex':  dict(path=ATLAS_DIR / "SrivatsanTrapnell2020_sciplex3.h5ad", pcol='perturbation', ctrl='control'),
}
PCAND = ['perturbation', 'perturbation_type', 'target_gene', 'guide', 'condition', 'gene']
CCAND = ['control', 'Control', 'CTRL', 'ctrl', 'non-targeting', 'NT', 'DMSO', 'unperturbed', 'None']
MAXPC, MAXCTRL, MAXPERT = 200, 1500, 400   # cells/pert cap, control cap, perturbation cap
M = 50; NSEED = 8; NBOOT = 500
ALPHAS = [0.0, 0.5, 0.7, 0.85, 0.925, 0.96, 1.0]; NA = len(ALPHAS); NPAIR = NA * (NA - 1) // 2
np.random.seed(0)


def pick(cols, cands):
    for c in cands:
        if c in cols:
            return c
    return None


def load(cfg):
    A = ad.read_h5ad(cfg['path'], backed='r')
    obs = A.obs
    pcol = cfg['pcol'] or pick(obs.columns, PCAND)
    ser = obs[pcol].astype(str)
    ctrl = cfg['ctrl'] or pick(set(ser.unique()), CCAND)
    rng = np.random.RandomState(0)
    # choose perturbations (>= 2*M cells), cap
    vc = ser.value_counts()
    perts = [p for p in vc.index if p != ctrl and p != 'nan' and vc[p] >= 2 * M]
    if len(perts) > MAXPERT:
        perts = list(rng.choice(perts, MAXPERT, replace=False))
    # gather subsampled row indices
    rows, labels = [], []
    for p in perts:
        idx = np.where(ser.values == p)[0]
        if len(idx) > MAXPC:
            idx = rng.choice(idx, MAXPC, replace=False)
        rows.append(idx); labels += [p] * len(idx)
    cidx = np.where(ser.values == ctrl)[0] if ctrl else np.array([], int)
    has_ctrl = len(cidx) > 0
    if len(cidx) > MAXCTRL:
        cidx = rng.choice(cidx, MAXCTRL, replace=False)
    allidx = np.concatenate(rows + ([cidx] if has_ctrl else []))
    order = np.argsort(allidx); inv = np.argsort(order)
    import scipy.sparse as sp
    Xraw = A[allidx[order]].to_memory().X
    A.file.close()
    lab = np.array(labels + (['__CTRL__'] * len(cidx) if has_ctrl else []))
    if Xraw.shape[1] > 40000:
        # large gene space (e.g. sci-Plex 110k): keep sparse, TruncatedSVD (avoid densify OOM)
        from sklearn.decomposition import TruncatedSVD
        Xraw = (Xraw.tocsr() if sp.issparse(Xraw) else sp.csr_matrix(np.asarray(Xraw)))[inv]
        Xp = TruncatedSVD(n_components=50, random_state=0).fit_transform(Xraw).astype(np.float32)
    else:
        Xsub = np.asarray(Xraw.todense() if hasattr(Xraw, 'todense') else Xraw, dtype=np.float32)[inv]
        fitidx = rng.choice(len(Xsub), min(20000, len(Xsub)), replace=False)
        Xp = PCA(n_components=min(50, Xsub.shape[1] - 1), random_state=0).fit(Xsub[fitidx]).transform(Xsub)
    return Xp, lab, has_ctrl, pcol, (ctrl or 'global-mean'), perts


def energy(A, B):
    dab = np.sqrt(((A[:, None] - B[None]) ** 2).sum(-1)).mean()
    daa = np.sqrt(((A[:, None] - A[None]) ** 2).sum(-1)).mean()
    dbb = np.sqrt(((B[:, None] - B[None]) ** 2).sum(-1)).mean()
    return 2 * dab - daa - dbb


def pds_row(drow, ci, n):
    td = drow[ci]
    if not np.isfinite(td):
        return 0.5
    less = int((drow < td - 1e-12).sum()); eq = int((np.abs(drow - td) <= 1e-12).sum())
    return 1.0 - (less + (eq - 1) / 2.0) / (n - 1) if n > 1 else 0.5


def norm(Mx):
    return Mx / (np.linalg.norm(Mx, axis=1, keepdims=True) + 1e-12)


def analyze(name):
    cfg = CFG[name]
    require_inputs(cfg['path'])
    Xp, lab, has_ctrl, pcol, ctrl, perts = load(cfg)
    cell = {p: Xp[lab == p] for p in perts}
    ctrl_cells = Xp[lab == '__CTRL__'] if has_ctrl else None
    gmean_all = Xp.mean(0)
    # split-half rankability + controlled-predictor per-pert PDS, seed-averaged
    dself, dnull = {p: [] for p in perts}, {p: [] for p in perts}
    pv = {a: {p: [] for p in perts} for a in ALPHAS}
    for seed in range(NSEED):
        rng = np.random.RandomState(seed)
        bg = ctrl_cells if has_ctrl else Xp   # background = control cells, else all cells (pseudo-control)
        base_e = bg[rng.choice(len(bg), min(M, len(bg)), replace=False)].mean(0)
        base_b = bg[rng.choice(len(bg), min(M, len(bg)), replace=False)].mean(0)
        build, evald = {}, {}
        for p in perts:
            c = cell[p]
            if len(c) < 2 * M:
                continue
            ii = rng.permutation(len(c))
            A, B = c[ii[:M]], c[ii[M:2 * M]]
            build[p] = A.mean(0) - base_b
            evald[p] = B.mean(0) - base_e
            dself[p].append(energy(A, B))
            dnull[p].append(energy(A, bg[rng.choice(len(bg), min(M, len(bg)), replace=False)]))
        cv = [p for p in perts if p in build]
        gm = np.stack([build[p] for p in cv]).mean(0)
        Tn = norm(np.stack([evald[p] for p in cv])); n = len(cv); idx = {p: i for i, p in enumerate(cv)}
        for a in ALPHAS:
            Pred = norm(np.stack([(1 - a) * gm + a * build[p] for p in cv]))
            D = 1 - Pred @ Tn.T
            for i, p in enumerate(cv):
                pv[a][p].append(pds_row(D[i], i, n))
    cv = [p for p in perts if pv[1.0][p]]
    # rankability
    rank = []
    for p in cv:
        s = np.median(dnull[p]) - np.median(dself[p])
        w = np.percentile(dself[p], 97.5) - np.percentile(dself[p], 2.5)
        rank.append(s > w)
    frac_rank = float(np.mean(rank))
    # reliability via bootstrap over perturbations
    pvm = {a: np.array([np.mean(pv[a][p]) for p in cv]) for a in ALPHAS}
    rng = np.random.RandomState(0); n = len(cv)
    gaps = [round(ALPHAS[k + 1] - ALPHAS[k], 3) for k in range(NA - 1)]
    winhi = np.zeros(NA - 1); n_correct = 0
    for _ in range(NBOOT):
        bi = rng.randint(0, n, n)
        yy = [pvm[a][bi].mean() for a in ALPHAS]
        if all(yy[k + 1] > yy[k] for k in range(NA - 1)):
            n_correct += 1
        for k in range(NA - 1):
            if yy[k + 1] > yy[k]:
                winhi[k] += 1
    winhi /= NBOOT
    reliability = n_correct / NBOOT
    # finest true-quality gap the benchmark resolves (adjacent alpha pair with P(higher>lower)>0.9)
    resolved_gaps = [gaps[k] for k in range(NA - 1) if winhi[k] > 0.9]
    min_res_gap = min(resolved_gaps) if resolved_gaps else None
    ceiling = float(pvm[1.0].mean())
    out = dict(dataset=name, pcol=pcol, ctrl=ctrl, n_pert=n, frac_rankable=round(frac_rank, 3),
               reliability_P_recover_order=round(reliability, 3),
               oracle_ceiling=round(ceiling, 3),
               min_resolvable_gap=min_res_gap,
               gap_resolution={str(gaps[k]): round(float(winhi[k]), 3) for k in range(NA - 1)},
               pds_by_alpha={str(a): round(float(pvm[a].mean()), 3) for a in ALPHAS})
    print(json.dumps(out, indent=2))
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    with open(OUT_DIR / f"{name}.json", "w") as fh:
        json.dump(out, fh, indent=2)


if __name__ == "__main__":
    analyze(_args.name)
