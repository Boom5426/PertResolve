#!/usr/bin/env python3
"""Shared evaluation harness for the PertResolve method comparison.

Every score in the manuscript passes through this module, so that the in-house
regression heads and the externally trained models are compared on exactly one
data source, one deterministic pseudobulk definition and one tie-aware PDS.

Three objects are canonical here and are deliberately not reimplemented anywhere
else:

``pseudobulk_deltas``
    the prediction target, a per-variant mean minus the wild-type mean, drawn
    with a fixed per-gene ``RandomState`` so the draw does not depend on the
    order in which genes or methods are scored.

``canonical_pds_cos``
    the tie-aware perturbation discrimination score. Ties take the average
    (mid) rank, which is what keeps a prediction that is exactly equidistant
    from several candidates at chance rather than at either extreme.

``split_vars``
    the held-out design, read from the committed benchmark table rather than
    recomputed, so a split cannot drift between the table and the scorer.

The per-gene expression arrays are large and are not redistributed with this
repository (see the manuscript's Data availability). Their directory is supplied
per run through ``--base`` or the ``PERTRESOLVE_DATA`` environment variable and
is resolved on first use, so importing this module never depends on a path that
existed only on the original machine. The benchmark table itself is committed and
is read from the repository.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from pertresolve.paths import repo_root, require_inputs, resolve_base

NSUB = 300
GENES = ['TP53', 'KRAS', 'GATA1', 'JAK1']
THETA = ['d_hydro', 'd_vol', 'd_charge', 'fold_core', 'cat_switch', 'is_hotspot']
WT_TAGS = ('WT', 'wt', 'WT_control')
_GENE_SEED = {'TP53': 0, 'KRAS': 1, 'GATA1': 2, 'JAK1': 3}

#: Resolved lazily so that importing this module does not require the workspace.
_BASE: Path | None = None

_BENCH_PATH = repo_root() / "data" / "pertresolve_bench.csv"
_bench = pd.read_csv(_BENCH_PATH)


def set_base(base) -> Path:
    """Pin the directory holding the per-gene expression arrays.

    Args:
        base: the workspace directory, or ``None`` to resolve it from
            ``PERTRESOLVE_DATA``.

    Returns:
        The resolved, existing directory.
    """
    global _BASE
    _BASE = resolve_base(base)
    return _BASE


def base() -> Path:
    """Return the pinned workspace, resolving it from the environment if needed."""
    if _BASE is None:
        set_base(None)
    return _BASE


def load_gene(gene):
    """Return ``(X, variant_labels)`` for one gene's standardized cell matrix.

    Args:
        gene: one of :data:`GENES`.

    Returns:
        A tuple of the cell-by-gene expression matrix and the per-cell variant
        labels, with wild-type cells carrying one of :data:`WT_TAGS`.
    """
    b = base()
    if gene in ("TP53", "KRAS"):
        path = b / "joint_arrays.npz"
        require_inputs(path)
        d = np.load(path, allow_pickle=True)
        sfx = "tp" if gene == "TP53" else "kr"
        return d[f"X{sfx}"], d[f"v{sfx}"]
    if gene == "GATA1":
        path = b / "gata1_arrays.npz"
        require_inputs(path)
        d = np.load(path, allow_pickle=True)
        return d["X"], d["cell_variants"]
    path = b / "jak1_arrays.npz"
    require_inputs(path)
    d = np.load(path, allow_pickle=True)
    return d["X"], d["variant_labels"]


def pseudobulk_deltas(gene, X=None, variants=None, n_sub=NSUB):
    """Deterministic per-variant pseudobulk difference from wild type.

    The draw uses a fixed per-gene ``RandomState`` and iterates variants in
    sorted order, so the returned profiles do not depend on call order. Variants
    with fewer than five cells are dropped rather than scored on a mean that the
    sample cannot support.

    Args:
        gene: one of :data:`GENES`, used to pick the per-gene seed.
        X: optional preloaded expression matrix; loaded from the workspace when
            omitted.
        variants: optional per-cell variant labels, required with ``X``.
        n_sub: cap on cells drawn per variant and for the wild-type reference.

    Returns:
        A tuple of ``{variant: delta}`` and the wild-type mean profile.
    """
    if X is None:
        X, variants = load_gene(gene)
    variants = np.asarray(variants)
    rng = np.random.RandomState(2024_0709 + _GENE_SEED[gene])
    wt_mask = np.isin(variants, WT_TAGS)
    if wt_mask.sum() == 0:
        wt_mean = np.zeros(X.shape[1], np.float32)
    else:
        idx = np.where(wt_mask)[0]
        if len(idx) > n_sub:
            idx = rng.choice(idx, n_sub, replace=False)
        wt_mean = X[idx].mean(0)
    deltas = {}
    for v in np.unique(variants):
        if v in WT_TAGS:
            continue
        idx = np.where(variants == v)[0]
        if len(idx) < 5:
            continue
        if len(idx) > n_sub:
            idx = rng.choice(idx, n_sub, replace=False)
        deltas[v] = (X[idx].mean(0) - wt_mean).astype(np.float32)
    return deltas, wt_mean.astype(np.float32)


def canonical_pds_cos(pred_delta, target_v, cand_vars, real_deltas):
    """Tie-aware cosine perturbation discrimination score for one variant.

    Args:
        pred_delta: the predicted response profile.
        target_v: the variant the prediction is for.
        cand_vars: the candidate set, the gene's training and held-out variants.
        real_deltas: measured profiles keyed by variant name.

    Returns:
        1.0 when the prediction's nearest candidate is its own variant, 0.5 in
        expectation under random matching, with ties taking the average rank. A
        zero-norm prediction, which has no direction to rank by, returns 0.5.
    """
    if np.linalg.norm(pred_delta) < 1e-12:
        return 0.5
    sims = []
    for vn in cand_vars:
        rd = real_deltas[vn]
        d = 1 - np.dot(pred_delta, rd) / (np.linalg.norm(pred_delta) * np.linalg.norm(rd) + 1e-12)
        sims.append((vn, d))
    sims.sort(key=lambda x: x[1])
    rank = next((i for i, (vn, _) in enumerate(sims) if vn == target_v), len(sims))
    td = sims[rank][1] if rank < len(sims) else float('inf')
    tied = [i for i, (vn, d) in enumerate(sims) if abs(d - td) < 1e-12]
    return 1.0 - float(np.mean(tied)) / (len(sims) - 1) if len(sims) > 1 else 0.5


def split_vars(gene, split='split1'):
    """Return ``(train_variants, test_variants)`` for one gene under one split."""
    gb = _bench[_bench.gene == gene]
    col = f"{split}_role"
    return list(gb[gb[col] == 'train']['variant']), list(gb[gb[col] == 'test']['variant'])


def theta_map(gene):
    """Return ``{variant: theta}`` for one gene, in the committed feature order."""
    gb = _bench[_bench.gene == gene]
    return {r['variant']: np.array([r[c] for c in THETA], np.float32) for _, r in gb.iterrows()}


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(description="Cache the deterministic per-variant profiles.")
    ap.add_argument("--base", default=None,
                    help="directory holding the per-gene arrays (env: PERTRESOLVE_DATA)")
    ap.add_argument("--out", required=True, metavar="OUT_DIR",
                    help="directory receiving real_deltas.npz")
    args = ap.parse_args()

    from pertresolve.paths import reject_repo_results

    set_base(args.base)
    out_dir = reject_repo_results(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    out = {}
    for gene in GENES:
        deltas, wt = pseudobulk_deltas(gene)
        for v, d in deltas.items():
            out[f"{gene}__{v}"] = d
        tr, te = split_vars(gene)
        te2 = [v for v in te if v in deltas]
        tr2 = [v for v in tr if v in deltas]
        print(f"{gene}: {len(deltas)} deltas, split1 train={len(tr2)} test={len(te2)}")
    np.savez_compressed(out_dir / "real_deltas.npz", **out)
    print(f"saved {len(out)} deterministic real deltas -> {out_dir / 'real_deltas.npz'}")
