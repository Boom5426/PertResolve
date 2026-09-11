#!/usr/bin/env python3
"""2A part 2: classifier two-sample test (the model-free identification probe).

A supervised classifier is trained to separate the single cells of two conditions under
cross-validation; its AUROC is a discriminative, distribution-free test of whether the
conditions differ at all. It is fundamentally different from the distributional distances
used elsewhere, so it rebuts "the floor is just energy-distance small-sample behaviour".
No batch control is needed: Ursu/GATA1/JAK1 are pooled designs (all variants share one
experiment, assigned per cell), so there is no per-variant batch confound in the processed
arrays; PCA is fit label-blind (no leakage) and classes are balanced.

  detection      : each variant vs WT
  identification : random sibling-variant pairs (the PDS-relevant question)
  chance check   : label-permuted identification pairs (should give AUROC ~0.5)
Reports per-gene median AUROC and the fraction of comparisons above a 0.55 threshold.
"""
import numpy as np, pandas as pd, sys
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score
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
                 help="directory receiving classifier_two_sample.csv; may not be inside the "
                      "repository's results/")
_args = _ap.parse_args()

BASE = resolve_base(_args.base)
OUT_DIR = reject_repo_results(_args.out)
OUT_DIR.mkdir(parents=True, exist_ok=True)
require_inputs(BASE / "pertresolve_bench.csv")
add_harness_to_path(BASE)

import harness as H

MINC = 30          # min cells per side
NPAIR = 100        # random identification pairs per gene
WT_TAGS = ('WT', 'wt', 'WT_control')
_BENCH_DF = pd.read_csv(f"{BASE}/pertresolve_bench.csv")
BENCH = {g: set(_BENCH_DF[_BENCH_DF.gene == g]["variant"]) for g in H.GENES}
np.random.seed(0)


def cv_auroc(Xa, Xb, rng):
    n = min(len(Xa), len(Xb))
    a = Xa[rng.choice(len(Xa), n, replace=False)]
    b = Xb[rng.choice(len(Xb), n, replace=False)]
    Xn = np.vstack([a, b]); y = np.r_[np.zeros(n), np.ones(n)]
    if n < 15:
        return np.nan
    aucs = []
    for tr, te in StratifiedKFold(5, shuffle=True, random_state=0).split(Xn, y):
        clf = LogisticRegression(max_iter=500, C=1.0).fit(Xn[tr], y[tr])
        p = clf.predict_proba(Xn[te])[:, 1]
        if len(np.unique(y[te])) == 2:
            aucs.append(roc_auc_score(y[te], p))
    return float(np.mean(aucs)) if aucs else np.nan


rows = []
for g in H.GENES:
    X, lab = H.load_gene(g); lab = np.asarray(lab)
    Xp = PCA(n_components=50, random_state=0).fit_transform(X.astype(np.float32))
    rng = np.random.RandomState(H._GENE_SEED[g])
    cells = {v: Xp[np.where(lab == v)[0]] for v in np.unique(lab)}
    wt = np.vstack([cells[v] for v in np.unique(lab) if v in WT_TAGS]) if any(v in WT_TAGS for v in cells) else None
    variants = [v for v in np.unique(lab) if v not in WT_TAGS and len(cells[v]) >= MINC and v in BENCH[g]]

    det = []
    if wt is not None and len(wt) >= MINC:
        for v in variants:
            det.append(cv_auroc(cells[v], wt, rng))
    det = [d for d in det if not np.isnan(d)]

    ident, perm = [], []
    pairs = [(rng.choice(variants), rng.choice(variants)) for _ in range(NPAIR * 3)]
    pairs = [(a, b) for a, b in pairs if a != b][:NPAIR]
    for a, b in pairs:
        auc = cv_auroc(cells[a], cells[b], rng)
        if not np.isnan(auc):
            ident.append(auc)
    # label-permuted chance check on a subset
    for a, b in pairs[:30]:
        pool = np.vstack([cells[a], cells[b]]); rng.shuffle(pool)
        h = len(cells[a])
        auc = cv_auroc(pool[:h], pool[h:h + len(cells[b])], rng)
        if not np.isnan(auc):
            perm.append(auc)

    rows.append(dict(
        gene=g, n_var=len(variants),
        detect_med_auroc=round(float(np.median(det)), 3) if det else np.nan,
        detect_frac_gt055=round(float(np.mean(np.array(det) > 0.55)), 3) if det else np.nan,
        ident_med_auroc=round(float(np.median(ident)), 3) if ident else np.nan,
        ident_frac_gt055=round(float(np.mean(np.array(ident) > 0.55)), 3) if ident else np.nan,
        perm_med_auroc=round(float(np.median(perm)), 3) if perm else np.nan,
        n_ident_pairs=len(ident),
    ))

df = pd.DataFrame(rows)
df.to_csv(OUT_DIR / "classifier_two_sample.csv", index=False)
print("Classifier two-sample AUROC (chance 0.5; perm_med should be ~0.5):")
print(df.to_string(index=False))
print(f"\nsaved -> {OUT_DIR / 'classifier_two_sample.csv'}")
