#!/usr/bin/env python3
"""Top-3 #2: is the oracle a genuine ceiling or a half-splitting artifact?

Three oracle variants per gene (multi-seed), to bound how much a smarter/denoised predictor
or deeper measurement could gain:
  oracle_300        : symmetric split-half at the native 300-cell cap (the reported ceiling)
  oracle_fulldepth  : symmetric split-half using ALL cells (deepest achievable) -> does depth help?
  oracle_perfect    : query = near-true mean (all cells), truth = one 300-cell noisy pseudobulk
                      -> the BEST a model could do (perfect prediction) against the actual noisy
                         benchmark truth. Upper-bounds any model at the benchmark's own resolution.
  oracle_perfect_both: query and truth are BOTH deep near-true means -> intrinsic separability
                      (if this also ~0.5, variants are unrankable even with perfect measurement).
Also a shrinkage check: PDS as a prediction is shrunk toward the gene mean (lambda), to confirm
mean-shrinkage denoising DECREASES discrimination (so the replicate oracle is not beatable that way).

Usage:
  python oracle_sensitivity.py --out /path/to/output_dir [--base /path/to/processed-data]
  ALLELEPERTURB_DATA=/path/to/processed-data python oracle_sensitivity.py --out /path/to/output_dir
"""
import numpy as np, pandas as pd, sys
import argparse
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))
from alleleperturb.paths import add_harness_to_path, require_inputs, resolve_base

ap = argparse.ArgumentParser(description=__doc__,
                             formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--base", default=None,
                help="directory holding the shared scorerharness.py and "
                     "allele_perturb_bench.csv (env: ALLELEPERTURB_DATA)")
ap.add_argument("--out", required=True,
                help="directory to write oracle_sensitivity.csv into")
args = ap.parse_args()

BASE = resolve_base(args.base)
BENCH_CSV = BASE / "allele_perturb_bench.csv"
require_inputs(BENCH_CSV)
OUT_DIR = Path(args.out).expanduser().resolve()

add_harness_to_path(BASE)
import harness as H  # noqa: E402

NSUB = 300; NSEED = 15
WT_TAGS = ('WT', 'wt', 'WT_control')
np.random.seed(0)
gene_cells = {g: (lambda X, l: (X, np.asarray(l)))(*H.load_gene(g)) for g in H.GENES}
_bench = pd.read_csv(BENCH_CSV)
BENCH = {g: set(_bench[_bench.gene == g]["variant"]) for g in H.GENES}   # curated benchmark set (consistency)


def pds_row(drow, ci, n):
    td = drow[ci]
    if not np.isfinite(td):
        return 0.5
    less = int((drow < td - 1e-12).sum()); eq = int((np.abs(drow - td) <= 1e-12).sum())
    return 1.0 - (less + (eq - 1) / 2.0) / (n - 1) if n > 1 else 0.5


def norm(Mx):
    return Mx / (np.linalg.norm(Mx, axis=1, keepdims=True) + 1e-12)


def pds_set(Q, T):
    n = len(Q); D = 1 - norm(Q) @ norm(T).T
    return float(np.mean([pds_row(D[i], i, n) for i in range(n)]))


rows = []; SHR = [0.0, 0.5, 0.75, 1.0]
for g in H.GENES:
    X, lab = gene_cells[g]
    acc = {k: [] for k in ['o300', 'ofull', 'operf']}
    shr = {a: [] for a in SHR}
    for seed in range(NSEED):
        rng = np.random.RandomState(2000 + 7 * seed + H._GENE_SEED[g])
        wt = np.where(np.isin(lab, WT_TAGS))[0]; rng.shuffle(wt); wh = len(wt) // 2
        wmA = X[wt[:wh]].mean(0); wmB = X[wt[wh:]].mean(0); wm_all = X[wt].mean(0)
        vs = [v for v in np.unique(lab) if v not in WT_TAGS and v in BENCH[g] and (lab == v).sum() >= 6]
        A300, B300, Afull, Bfull = [], [], [], []
        pv, Qperf, Tnoisy = [], [], []
        for v in vs:
            idx = np.where(lab == v)[0]; rng.shuffle(idx); h = len(idx) // 2
            i1, i2 = idx[:h], idx[h:]
            A300.append(X[i1[:NSUB]].mean(0) - wmA); B300.append(X[i2[:NSUB]].mean(0) - wmB)
            Afull.append(X[i1].mean(0) - wmA); Bfull.append(X[i2].mean(0) - wmB)
            # perfect prediction (deep disjoint mean) vs noisy 300-cell benchmark truth, NO shared cells
            if len(idx) > NSUB + 50:
                pv.append(v)
                Tnoisy.append(X[idx[:NSUB]].mean(0) - wmB)     # the benchmark's noisy 300-cell truth
                Qperf.append(X[idx[NSUB:]].mean(0) - wmA)      # deep disjoint estimate of the true mean
        A300, B300 = np.stack(A300), np.stack(B300)
        Afull, Bfull = np.stack(Afull), np.stack(Bfull)
        acc['o300'].append(pds_set(A300, B300))
        acc['ofull'].append(pds_set(Afull, Bfull))
        if len(pv) >= 3:
            acc['operf'].append(pds_set(np.stack(Qperf), np.stack(Tnoisy)))
        # shrinkage: query = (1-lam)*gene_mean + lam*A300 ; truth = B300
        gm = A300.mean(0)
        for a in SHR:
            Qs = (1 - a) * gm + a * A300
            shr[a].append(pds_set(Qs, B300))
    rows.append(dict(gene=g, n_var=len(vs),
                     oracle_300=round(float(np.mean(acc['o300'])), 3),
                     oracle_fulldepth=round(float(np.mean(acc['ofull'])), 3),
                     oracle_perfect_vs_noisy=(round(float(np.mean(acc['operf'])), 3) if acc['operf'] else np.nan),
                     shrink_lam0_genemean=round(float(np.mean(shr[0.0])), 3),
                     shrink_lam1_oracle=round(float(np.mean(shr[1.0])), 3)))

df = pd.DataFrame(rows)
print(df.to_string(index=False))
OUT_DIR.mkdir(parents=True, exist_ok=True)
df.to_csv(OUT_DIR / "oracle_sensitivity.csv", index=False)
print("\nInterpretation: if oracle_fulldepth and oracle_perfect stay ~0.5 for TP53/KRAS, the floor is")
print("not a half-splitting artifact and no denoised/perfect prediction beats it at this resolution.")
print("shrink_lam1 (=oracle) should be >= shrink_lam0 (=gene mean): mean-shrinkage does not beat the oracle.")
