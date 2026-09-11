#!/usr/bin/env python3
"""1A: Oracle measurement ceiling PDS_oracle.

For each held-out variant the 'prediction' is a SECOND, DISJOINT measurement of the
SAME variant (independent split-half pseudobulk delta, with an independent WT half).
This upper-bounds the PDS any real model can reach at this measurement depth, because a
real model carries model error ON TOP of this finite-sample noise.

Mirrors score_definitive.py exactly: same GENES, SPLITS, NSUB, tie-aware mid-rank
PDS_cos, 15-seed averaging, bootstrap CI over variants; candidate set = train+test of
each split; only test variants scored. So PDS_oracle is directly comparable to the model
PDS in definitive_summary.csv.
"""
import numpy as np, pandas as pd, sys
import argparse
from pathlib import Path

# The package is not installed by default; make this repository importable so that
# the shared path helpers can be used when the script is run directly.
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from pertresolve.paths import (
    add_harness_to_path,
    reject_repo_results,
    require_inputs,
    resolve_base,
)

_ap = argparse.ArgumentParser(description=__doc__,
                              formatter_class=argparse.RawDescriptionHelpFormatter)
_ap.add_argument("--base", default=None,
                 help="directory holding the shared scorer and the gene "
                      "arrays (env: PERTRESOLVE_DATA)")
_ap.add_argument("--out", required=True, metavar="OUT_DIR",
                 help="directory receiving oracle_ceiling.csv; may not be inside the "
                      "repository's results/")
_args = _ap.parse_args()

BASE = resolve_base(_args.base)
OUT_DIR = reject_repo_results(_args.out)
OUT_DIR.mkdir(parents=True, exist_ok=True)
add_harness_to_path(BASE)

import harness as H

NSUB = 300
NSEED = 15
NBOOT = 2000
SPLITS = ['split1', 'split2', 'split3', 'split5', 'split6']
WT_TAGS = ('WT', 'wt', 'WT_control')
np.random.seed(0)

cand = {(g, s): (H.split_vars(g, s)[0] + H.split_vars(g, s)[1]) for g in H.GENES for s in SPLITS}
gene_cells = {g: (lambda X, lab: (X, np.asarray(lab)))(*H.load_gene(g)) for g in H.GENES}


def pds_row(drow, ci, n):
    td = drow[ci]
    if not np.isfinite(td):
        return 0.5
    less = int((drow < td - 1e-12).sum())
    eq = int((np.abs(drow - td) <= 1e-12).sum())
    return 1.0 - (less + (eq - 1) / 2.0) / (n - 1) if n > 1 else 0.5


def norm(M):
    return M / (np.linalg.norm(M, axis=1, keepdims=True) + 1e-12)


def halves_delta(seed):
    """Disjoint split-half pseudobulk deltas: truth[v] from one half + WT half,
    pred[v] from the other (disjoint) half + independent WT half."""
    truth, pred = {}, {}
    for g in H.GENES:
        X, lab = gene_cells[g]
        rng = np.random.RandomState(1000 + 7 * seed + H._GENE_SEED[g])
        wt = np.where(np.isin(lab, WT_TAGS))[0]
        if len(wt) == 0:
            wm_t = wm_p = np.zeros(X.shape[1], np.float32)
        else:
            rng.shuffle(wt)
            h = len(wt) // 2
            wt_t, wt_p = wt[:h], wt[h:]
            if len(wt_t) > NSUB:
                wt_t = wt_t[:NSUB]
            if len(wt_p) > NSUB:
                wt_p = wt_p[:NSUB]
            wm_t = X[wt_t].mean(0)
            wm_p = X[wt_p].mean(0)
        for v in np.unique(lab):
            if v in WT_TAGS:
                continue
            idx = np.where(lab == v)[0]
            if len(idx) < 5:
                continue
            rng.shuffle(idx)
            h = len(idx) // 2
            it, ip = idx[:h], idx[h:]
            if len(it) > NSUB:
                it = it[:NSUB]
            if len(ip) > NSUB:
                ip = ip[:NSUB]
            truth[f"{g}__{v}"] = (X[it].mean(0) - wm_t).astype(np.float32)
            pred[f"{g}__{v}"] = (X[ip].mean(0) - wm_p).astype(np.float32)
    return truth, pred


acc = {}
for seed in range(NSEED):
    truth, pred = halves_delta(seed)
    for g in H.GENES:
        for s in SPLITS:
            cv = [v for v in cand[(g, s)] if f"{g}__{v}" in truth]
            if len(cv) < 2:
                continue
            idx = {v: i for i, v in enumerate(cv)}
            te = [v for v in H.split_vars(g, s)[1] if f"{g}__{v}" in truth and v in idx]
            if not te:
                continue
            Tm = norm(np.stack([truth[f"{g}__{u}"] for u in cv]))
            Pm = norm(np.stack([pred[f"{g}__{v}"] for v in te]))
            D = 1 - Pm @ Tm.T
            for i, v in enumerate(te):
                acc.setdefault((g, s), {}).setdefault(v, []).append(pds_row(D[i], idx[v], len(cv)))


def summarize(keys):
    """Cell unit: mean of (gene, split) cell means, bootstrapped within cells.

    Retained unchanged so a re-run reproduces the previously committed numbers exactly and
    the migration to the primary unit below is auditable rather than asserted.
    """
    cells = {k: {v: float(np.mean(l)) for v, l in acc[k].items()} for k in keys if k in acc}
    cellmeans = {k: np.mean(list(d.values())) for k, d in cells.items()}
    obs = float(np.mean(list(cellmeans.values())))
    boot = []
    for _ in range(NBOOT):
        cm = [np.array(list(d.values()))[np.random.randint(0, len(d), len(d))].mean() for d in cells.values()]
        boot.append(np.mean(cm))
    lo, hi = np.percentile(boot, [2.5, 97.5])
    n_var = sum(len(d) for d in cells.values())
    return obs, float(lo), float(hi), n_var


def summarize_variant(keys):
    """Distinct-variant unit: THE primary estimator, matching the model scores.

    Byte-for-byte the aggregation ``regenerate_canonical_pds.py:summarise`` applies to
    ``definitive_summary.csv``'s ``PDS`` column, so the reference and the models it is
    plotted against are the same estimator on the same unit. A variant held out under
    several splits is averaged across them first, so it enters the point estimate and the
    bootstrap once; Methods declares the distinct variant as the primary statistical unit
    for model evaluation and this makes the reference obey that declaration.

    Uses its own ``default_rng``, leaving the global stream :func:`summarize` draws from
    untouched, so adding this function does not move the cell-unit numbers.
    """
    per_gsv = {}
    for k in keys:
        if k not in acc:
            continue
        g, _s = k
        for v, l in acc[k].items():
            per_gsv.setdefault((g, v), []).append(float(np.mean(l)))
    distinct = np.array([float(np.mean(l)) for l in per_gsv.values()])
    obs = float(distinct.mean())
    rng = np.random.default_rng(0)
    dv = rng.integers(0, len(distinct), size=(NBOOT, len(distinct)))
    lo, hi = np.percentile(distinct[dv].mean(axis=1), [2.5, 97.5])
    return obs, float(lo), float(hi), len(distinct)


# The per-(gene, split, variant) scores the two aggregations are computed from, written
# out so that either unit can be re-derived, and the choice between them audited, without
# the per-gene cell matrices this script needs.
pv_rows = [dict(gene=g, split=s, variant=v, pds=float(np.mean(l)), n_seeds=len(l))
           for (g, s), d in acc.items() for v, l in d.items()]
pd.DataFrame(pv_rows).sort_values(["gene", "split", "variant"]).to_csv(
    OUT_DIR / "oracle_ceiling_per_variant.csv", index=False)

rows = []
scopes = [(g, [(g, s) for s in SPLITS if (g, s) in acc]) for g in H.GENES]
scopes.append(('ALL', list(acc.keys())))
for scope, keys in scopes:
    # Cell unit first and from the global stream, so its draws are the same ones the
    # previous version made and the committed cell-unit values reproduce exactly.
    oc, loc, hic, ncell_var = summarize(keys)
    ov, lov, hiv, nvar = summarize_variant(keys)
    rows.append(dict(scope=scope,
                     PDS_oracle=round(ov, 3), ci_lo=round(lov, 3), ci_hi=round(hiv, 3),
                     n_var_scored=nvar,
                     PDS_oracle_cell=round(oc, 3), cell_lo=round(loc, 3),
                     cell_hi=round(hic, 3), n_cells=len({k for k in keys if k in acc}),
                     n_var_by_cell=ncell_var))

df = pd.DataFrame(rows)
df.to_csv(OUT_DIR / "oracle_ceiling.csv", index=False)
print(df.to_string(index=False))
print(f"\nsaved -> {OUT_DIR / 'oracle_ceiling.csv'}")
print(f"saved -> {OUT_DIR / 'oracle_ceiling_per_variant.csv'} ({len(pv_rows)} rows)")
print("\nPDS_oracle is now the distinct-variant unit, matching definitive_summary.csv's PDS.")
print("PDS_oracle_cell reproduces the previously committed PDS_oracle column.")
