#!/usr/bin/env python3
"""Refit the per-gene standardization on training cells only, then rescore the grid.

The committed expression arrays are standardized across all of a gene's cells,
including the cells of the variants that a split holds out. Because every
comparison in this study is between two pseudobulk means of the same gene, the
centring cancels exactly and only the per-gene scale survives:

    delta_train-only = delta_all-cell / sd(X_all-cell[train union wild-type])

so refitting the standardization on training variants and wild-type cells alone
is a diagonal rescale of the committed arrays, applied before any profile is
formed. The scale is refit per (gene, split), because each split holds out a
different set of variants.

The rescale is not a no-op for the score. PDS ranks by cosine distance, which is
not invariant under an anisotropic diagonal transform, so predictions and
evaluation truth must both be formed in the rescaled space. Both are, here.

This script fits and scores the twenty in-house predictors: the eighteen
feature-model heads formed by crossing six regression heads with three feature
spaces, plus the Gene-mean and WT-null references. The externally trained models
are fitted by their own packages, so their predictions are supplied through
``--external-preds`` and scored here against the same evaluation truth; they must
have been fitted in the same space, which the exported per-(gene, split) scale
factors make possible.

Two intervals are reported for every method. The primary one, written to the
conventional ``PDS`` / ``ci_lo`` / ``ci_hi`` columns, resamples distinct
(gene, variant) pairs, so an allele held out by several splits counts once. The
secondary one resamples held-out variants within each (split, gene) cell, which
is the convention earlier runs of this study used and is kept for comparability.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from pertresolve.paths import reject_repo_results, require_inputs, resolve_base

_ap = argparse.ArgumentParser(description=__doc__,
                              formatter_class=argparse.RawDescriptionHelpFormatter)
_ap.add_argument("--base", default=None,
                 help="directory holding the gene arrays and esm1v_embeddings.npz "
                      "(env: PERTRESOLVE_DATA)")
_ap.add_argument("--out", required=True, metavar="OUT_DIR",
                 help="directory receiving the tables; may not be inside results/")
_ap.add_argument("--standardization", choices=("train-only", "all-cell"),
                 default="train-only",
                 help="'all-cell' reproduces the previous convention and exists so "
                      "the two can be compared through one code path")
_ap.add_argument("--export-preds", action="store_true",
                 help="also write this run's in-process predictions to "
                      "<out>/preds/<method>.npz, keyed '<gene>__<split>__<variant>', "
                      "so the preds5-consuming analyses can be regenerated from the "
                      "same fit rather than from a separately maintained runner")
_ap.add_argument("--dump-rank-counts", action="store_true",
                 help="also write <out>/per_seed_rank_counts_<TAG>.csv.gz carrying the "
                      "raw ranking counts behind every scored row: 'less' candidates "
                      "strictly closer than the true allele, 'eq' candidates within the "
                      "tie tolerance including itself, and the pool size 'n'. PDS is a "
                      "mid-rank summary of these three numbers, so dumping them lets a "
                      "tail statistic such as Hit@k be derived from the same ranking that "
                      "produced the score rather than from a second pipeline. Adds one "
                      "file; changes no existing output.")
_ap.add_argument("--external-preds", default=None, metavar="DIR",
                 help="directory of externally trained models' prediction .npz files, "
                      "keyed '<gene>__<split>__<variant>'. They must have been fitted "
                      "in the same standardization as --standardization, since the "
                      "evaluation truth here is formed in that space.")
_args = _ap.parse_args()

BASE = resolve_base(_args.base)
OUT_DIR = reject_repo_results(_args.out)
OUT_DIR.mkdir(parents=True, exist_ok=True)
MODE = _args.standardization
TAG = "trainonly" if MODE == "train-only" else "allcell"

sys.path.insert(0, str(Path(__file__).resolve().parent))
import harness as H  # noqa: E402

H.set_base(BASE)

import warnings  # noqa: E402

warnings.filterwarnings("ignore")
from sklearn.base import clone  # noqa: E402
from sklearn.decomposition import PCA  # noqa: E402
from sklearn.ensemble import GradientBoostingRegressor, RandomForestRegressor  # noqa: E402
from sklearn.linear_model import Lasso, Ridge  # noqa: E402
from sklearn.multioutput import MultiOutputRegressor  # noqa: E402
from sklearn.neighbors import KNeighborsRegressor  # noqa: E402
from sklearn.neural_network import MLPRegressor  # noqa: E402

NSEED = 15
#: theta's first three components are physicochemical differences in physical units.
N_PHYSCHEM = 3
NSUB = 300
NBOOT = 2000
SPLITS = ['split1', 'split2', 'split3', 'split5', 'split6']

ESM_PATH = BASE / "esm1v_embeddings.npz"
require_inputs(ESM_PATH)
esm = np.load(ESM_PATH)

HEADS = {
    'Ridge': Ridge(alpha=1.0),
    'Lasso': Lasso(alpha=0.01, max_iter=2000),
    'RF': RandomForestRegressor(n_estimators=50, max_depth=6, n_jobs=-1, random_state=0),
    'GBoost': 'gb',
    'KNN': KNeighborsRegressor(n_neighbors=3),
    'MLP': MLPRegressor(hidden_layer_sizes=(128, 64), max_iter=500, random_state=0),
}


def feat(gene, v, ft, tm):
    """Return the feature vector for one variant, or None where it is undefined."""
    if ft == 'theta':
        return tm.get(v)
    e = esm[f"{gene}__{v}"] if f"{gene}__{v}" in esm.files else None
    if ft == 'esm':
        return e
    t = tm.get(v)
    return np.concatenate([e, t]) if (e is not None and t is not None) else None


def scale_factor(X, variants, train_vars):
    """Per-gene-dimension scale refit on training and wild-type cells only.

    Args:
        X: the committed all-cell standardized matrix for one gene.
        variants: per-cell variant labels.
        train_vars: the split's training variants.

    Returns:
        A vector of per-gene divisors. Dimensions with no variation among the
        training cells are left at 1.0 rather than amplified into noise.
    """
    if MODE == "all-cell":
        return np.ones(X.shape[1], dtype=np.float32)
    keep = np.isin(variants, list(train_vars)) | np.isin(variants, H.WT_TAGS)
    s = X[keep].std(0)
    s[s < 1e-8] = 1.0
    return s.astype(np.float32)


def deltas_from(X, variants, seed, gene):
    """Per-variant pseudobulk profiles at one subsample seed, in the given space."""
    rng = np.random.RandomState(1000 + 7 * seed + H._GENE_SEED[gene])
    wt = np.where(np.isin(variants, H.WT_TAGS))[0]
    wt = rng.choice(wt, NSUB, replace=False) if len(wt) > NSUB else wt
    wm = X[wt].mean(0)
    out = {}
    for v in np.unique(variants):
        if v in H.WT_TAGS:
            continue
        idx = np.where(variants == v)[0]
        if len(idx) < 5:
            continue
        if len(idx) > NSUB:
            idx = rng.choice(idx, NSUB, replace=False)
        out[v] = (X[idx].mean(0) - wm).astype(np.float32)
    return out


def norm(M):
    return M / (np.linalg.norm(M, axis=1, keepdims=True) + 1e-12)


def pds_row(drow, ci, n):
    td = drow[ci]
    if not np.isfinite(td):
        return 0.5
    less = int((drow < td - 1e-12).sum())
    eq = int((np.abs(drow - td) <= 1e-12).sum())
    return 1.0 - (less + (eq - 1) / 2.0) / (n - 1) if n > 1 else 0.5


def rank_counts(drow, ci, n):
    """The three numbers pds_row reduces to a mid-rank, kept separately.

    Returns (less, eq, finite). ``pds_row`` is deliberately left untouched and
    recomputed independently by the caller, so that the dumped counts can be checked
    against the scored value rather than being trusted to agree with it.

    A non-finite distance to the target means the prediction has no ranking at all.
    That is reported as finite=False with less and eq unset rather than as a rank,
    because a missing ordering is not the same as a bad one; the consumer decides the
    convention and must state it.
    """
    td = drow[ci]
    if not np.isfinite(td):
        return -1, -1, False
    less = int((drow < td - 1e-12).sum())
    eq = int((np.abs(drow - td) <= 1e-12).sum())
    return less, eq, True


def fit_predict(Xtr, Ytr, Xte, spec, n_out):
    """Fit one head on the training variants and predict the held-out ones."""
    if spec == 'gb':
        q = min(50, n_out, max(1, len(Xtr) - 1))
        pca = PCA(q, random_state=0)
        Yp = pca.fit_transform(Ytr)
        gm = MultiOutputRegressor(
            GradientBoostingRegressor(n_estimators=50, max_depth=3, random_state=0), n_jobs=-1)
        gm.fit(Xtr, Yp)
        return pca.inverse_transform(gm.predict(Xte))
    mm = clone(spec)
    if isinstance(mm, RandomForestRegressor) and n_out > 200:
        q = min(50, n_out, max(1, len(Xtr) - 1))
        pca = PCA(q, random_state=0)
        Yp = pca.fit_transform(Ytr)
        mm.fit(Xtr, Yp)
        return pca.inverse_transform(mm.predict(Xte))
    mm.fit(Xtr, Ytr)
    return mm.predict(Xte)


EXTERNAL = {}
if _args.external_preds:
    ext_dir = Path(_args.external_preds)
    require_inputs(ext_dir)
    for f in sorted(ext_dir.glob("*.npz")):
        EXTERNAL[f.stem] = np.load(f)
    print(f"scoring {len(EXTERNAL)} externally trained models from {ext_dir}", flush=True)
    # A stored prediction silently REPLACES this run's own fit for the same method
    # name, so pointing --external-preds at a directory that also holds the in-house
    # heads makes the run score old predictions instead of the ones just fitted.
    shadowed = sorted(n for n in EXTERNAL
                      if any(n.startswith(h + "-") for h in HEADS)
                      or n in ("Gene-mean", "WT-null"))
    if shadowed:
        raise SystemExit(
            f"--external-preds {ext_dir} contains in-house methods {shadowed}; these "
            "would overwrite the fits this run just produced. Point it at a directory "
            "holding only the externally trained models.")

raw = {g: (lambda X, lab: (X, np.asarray(lab)))(*H.load_gene(g)) for g in H.GENES}
scales = {}
long_rows = []
count_rows = []
exported: dict[str, dict[str, np.ndarray]] = {}

for gene in H.GENES:
    X0, variants = raw[gene]
    tm = H.theta_map(gene)
    for s in SPLITS:
        tr_all, te_all = H.split_vars(gene, s)
        sc = scale_factor(X0, variants, tr_all)
        scales[(gene, s)] = sc
        X = X0 / sc
        # Predictions are fitted once per (gene, split) on the harness's deterministic
        # profiles, exactly as in the published grid, so that under --standardization
        # all-cell this script reproduces the committed predictions; only the
        # evaluation truth is reseeded across the 15 subsamples.
        det, _ = H.pseudobulk_deltas(gene, X, variants)
        tr = [v for v in tr_all if v in det]
        te = [v for v in te_all if v in det]
        if len(tr) < 3 or len(te) < 1:
            continue
        G = det[tr[0]].shape[0]

        # Standardize theta's three physicochemical terms on THIS SPLIT'S TRAINING
        # VARIANTS, matching scripts/run_all_splits.py. configs/bench_features.yaml
        # emits them in physical units and assigns standardization to the fitting
        # step; a scaler fitted over the whole table would let held-out variants set
        # the scale of the training features. The three 0/1 annotations are left as is.
        tm_s = {v: f.copy() for v, f in tm.items()}
        fit_v = [v for v in tr if v in tm_s]
        if len(fit_v) >= 2:
            fit_m = np.stack([tm_s[v][:N_PHYSCHEM] for v in fit_v])
            mu = fit_m.mean(axis=0)
            sd = fit_m.std(axis=0, ddof=1)
            sd = np.where(np.isfinite(sd) & (sd > 0), sd, 1.0)
            for vec in tm_s.values():
                vec[:N_PHYSCHEM] = (vec[:N_PHYSCHEM] - mu) / sd

        preds = {}
        for ft in ['theta', 'esm', 'esm+theta']:
            trok = [v for v in tr if feat(gene, v, ft, tm_s) is not None]
            teok = [v for v in te if feat(gene, v, ft, tm_s) is not None]
            if len(trok) < 3 or len(teok) < 1:
                continue
            Xtr = np.stack([feat(gene, v, ft, tm_s) for v in trok])
            Ytr = np.stack([det[v] for v in trok])
            Xte = np.stack([feat(gene, v, ft, tm_s) for v in teok])
            for hn, spec in HEADS.items():
                Yte = fit_predict(Xtr, Ytr, Xte, spec, G)
                for i, v in enumerate(teok):
                    preds.setdefault(f"{hn}-{ft}", {})[v] = Yte[i].astype(np.float32)
        mean_d = np.stack([det[v] for v in tr]).mean(0)
        for v in te:
            preds.setdefault('Gene-mean', {})[v] = mean_d
            preds.setdefault('WT-null', {})[v] = np.zeros(G, np.float32)

        if _args.export_preds:
            for m, pv in preds.items():
                for v, vec in pv.items():
                    exported.setdefault(m, {})[f"{gene}__{s}__{v}"] = vec

        # Externally trained models supply their own predictions for this cell.
        for m, store in EXTERNAL.items():
            for v in te:
                key = f"{gene}__{s}__{v}"
                if key in store.files:
                    preds.setdefault(m, {})[v] = store[key].astype(np.float32)

        for seed in range(NSEED):
            truth = deltas_from(X, variants, seed, gene)
            cv = [v for v in (tr_all + te_all) if v in truth]
            idx = {v: i for i, v in enumerate(cv)}
            Tm = norm(np.stack([truth[u] for u in cv]))
            for m, pv in preds.items():
                names = [v for v in pv if v in idx]
                if not names:
                    continue
                Pm = norm(np.stack([pv[v] for v in names]))
                D = 1 - Pm @ Tm.T
                for i, v in enumerate(names):
                    long_rows.append((m, gene, s, v, seed, pds_row(D[i], idx[v], len(cv))))
                    if _args.dump_rank_counts:
                        less, eq, finite = rank_counts(D[i], idx[v], len(cv))
                        count_rows.append((m, gene, s, v, seed, less, eq, len(cv), finite))
        print(f"{gene} {s}: {len(preds)} methods, {len(te)} held out", flush=True)

if _args.dump_rank_counts:
    counts = pd.DataFrame(count_rows, columns=["method", "gene", "split", "variant",
                                               "seed", "less", "eq", "n", "finite"])
    counts.to_csv(OUT_DIR / f"per_seed_rank_counts_{TAG}.csv.gz", index=False,
                  compression="gzip")
    print(f"rank counts: {len(counts)} rows, {int((~counts.finite).sum())} non-finite "
          f"-> {OUT_DIR / f'per_seed_rank_counts_{TAG}.csv.gz'}", flush=True)

long = pd.DataFrame(long_rows, columns=["method", "gene", "split", "variant", "seed", "pds"])
long.to_csv(OUT_DIR / f"per_seed_variant_pds_{TAG}.csv.gz", index=False,
            float_format="%.6f", compression="gzip")

np.savez_compressed(OUT_DIR / f"scale_factors_{TAG}.npz",
                    **{f"{g}__{s}": v for (g, s), v in scales.items()})

seed_avg = long.groupby(["method", "gene", "split", "variant"])["pds"].mean().reset_index()
rows = []
for m, sub in seed_avg.groupby("method"):
    rng = np.random.default_rng(0)
    per_var = sub.groupby(["gene", "variant"])["pds"].mean().to_numpy()
    boot_v = np.array([per_var[rng.integers(0, len(per_var), len(per_var))].mean()
                       for _ in range(NBOOT)])
    lo_v, hi_v = np.percentile(boot_v, [2.5, 97.5])
    cells = [g["pds"].to_numpy() for _, g in sub.groupby(["gene", "split"])]
    rng = np.random.default_rng(0)
    boot_c = np.array([np.mean([c[rng.integers(0, len(c), len(c))].mean() for c in cells])
                       for _ in range(NBOOT)])
    lo_c, hi_c = np.percentile(boot_c, [2.5, 97.5])
    # PDS / ci_lo / ci_hi carry the primary estimator, so that a consumer reading
    # the conventional column names gets the variant-clustered interval rather than
    # the cell-wise one. The cell-wise figures are kept beside them, named.
    rows.append(dict(
        method=m, standardization=MODE,
        PDS=round(float(per_var.mean()), 4),
        ci_lo=round(float(lo_v), 4), ci_hi=round(float(hi_v), 4),
        crosses=bool(lo_v <= 0.5 <= hi_v), n_variants=len(per_var),
        PDS_cell=round(float(np.mean([c.mean() for c in cells])), 4),
        cell_lo=round(float(lo_c), 4), cell_hi=round(float(hi_c), 4),
        cell_crosses=bool(lo_c <= 0.5 <= hi_c), n_cells=len(cells),
    ))

if _args.export_preds:
    pred_dir = OUT_DIR / "preds"
    pred_dir.mkdir(parents=True, exist_ok=True)
    for m, store in exported.items():
        np.savez_compressed(pred_dir / f"{m}.npz", **store)
    print(f"exported {len(exported)} in-house prediction files -> {pred_dir}", flush=True)

df = pd.DataFrame(rows).sort_values("PDS", ascending=False)
df.to_csv(OUT_DIR / f"grid_summary_{TAG}.csv", index=False)
print(df.to_string(index=False))
heads = df[df.method.str.contains("Ridge|Lasso|RF|GBoost|KNN|MLP")]
print(f"\n[{MODE}] 18 heads: variant-unit {heads.PDS.mean():.4f}, "
      f"cell-unit {heads.PDS_cell.mean():.4f}")
print(f"methods whose primary interval excludes chance: {list(df[~df.crosses].method)}")
