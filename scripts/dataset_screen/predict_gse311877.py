#!/usr/bin/env python3
"""GSE311877_FEATURE_TO_RESIDUAL_PREDICTION_v1.

Executes docs/PREREG_GSE311877_FEATURE_TO_RESIDUAL_v1.md verbatim. Read that file
first; its verdict rules are binding and are applied here unmodified. This script
decides nothing it was not told in advance to decide.

The question: the TP63 allele-specific residual is reproducible (split-half
residual PDS 0.727 at chance 0.500). Can it be PREDICTED from protein variant
features? Only `ceiling >> chance` together with `model ~ chance` makes TP63 a
second model-limited system.

Leakage control, per fold with held-out allele `a`: the gene filter, the residual
PC basis, feature standardisation, the ESM PCA and every hyperparameter are
estimated on the 17 training alleles only. `a` enters nothing but its own scoring.

Why the permutation test is affordable. Every estimator here is a linear smoother
in SAMPLE space: the prediction is `r_hat_a = c^T Yc` for some 17-vector `c`,
where `Yc` holds the centred full-depth training residuals. Therefore

    corr(r_hat_a, pool_j) = (c^T M)_j / sqrt(c^T K c),
    K = Yc Yc^T  (17x17),   M = Yc pool^T  (17x18)

exactly. `K` and `M` are built once per fold, so no 12,700-dimensional prediction
is ever materialised and 2,000 permutations cost seconds rather than hours. This
is an algebraic identity, not an approximation. PLS2 also fits: its NIPALS
weights depend on the targets only through `K`, and its prediction reduces to
`c = sum_h t*_h t_h / (t_h^T t_h)`.

Ambiguity in the pre-registration, resolved here and reported. Section 5 says the
candidate pool is "averaged over the three balanced splits". Averaging the
PROFILES would place replicates 2,3,4 in every pool member while the ceiling's
half-A predictor also holds replicates 2,3,4, overlapping predictor with pool and
inflating the ceiling. This script averages the SCORES over the three splits, the
only reading that keeps half-A and half-B disjoint. No verdict rule is affected.
"""
from __future__ import annotations

import json
import os
import re
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from alleleperturb.features import compute_theta  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
RAW = Path(os.environ.get("GSE311877_DIR", REPO / "data_external" / "GSE311877"))
FEAT = Path(os.environ.get("TP63_FEATURES", REPO / "data_external" / "tp63_features"))
OUT = Path(os.environ.get("DATASET_SCREEN_OUT", REPO / "results" / "dataset_screen"))
ISOFORM = os.environ.get("TP63_ISOFORM", "dnp63a")
# 0 = exactly as frozen (PCs of the allele-centred matrix only);
# 1 = frozen intent (grand mean also removed, so chance PDS = 0.5).
REMOVE_MEAN = bool(int(os.environ.get("REMOVE_MEAN", "1")))
TAG = f"{ISOFORM}_{'meanrm' if REMOVE_MEAN else 'asfrozen'}"

SEED = 0
N_PC = 5                     # frozen: shared components removed to form the residual
N_ESM_PC = 10                # frozen: ESM blocks reduced to min(10, |T|-1) in-fold
N_BOOT = 10_000
N_PERM = 2_000
MIN_CPM, MIN_FRAC = 1.0, 0.5
ALPHAS = np.logspace(-3, 6, 40)
PLS_COMPONENTS = (1, 2, 3)
SPLITS = [((1, 2), (3, 4)), ((1, 3), (2, 4)), ((1, 4), (2, 3))]
FNAME = re.compile(r"^GSM\d+_(?:[A-H]\d+_)?(.+?)_(\d)_rawCounts\.txt\.gz$")

# TAp63-alpha coordinates (UniProt Q9H3D4): DNA binding 170-362, SAM 541-607.
# Label position N maps to TAp63-alpha residue N+39, per the pre-registration.
LABEL_TO_TA = 39
DOMAINS = [(170, 362), (541, 607)]
SIBLINGS = [("R279Q", "R279S"), ("R304Q", "R304T"), ("L514D", "L514F"),
            ("C522D", "C522G"), ("L531E", "L531R")]
FEAT_SETS = ["esm2_sitedelta", "esm2_window16", "esm2_globaldelta", "theta",
             "esm2_sitedelta+theta", "position", "region"]


# --------------------------------------------------------------------------- data
def load_cpm() -> tuple[pd.DataFrame, list[str]]:
    cols = {}
    for f in sorted(RAW.glob("*_rawCounts.txt.gz")):
        m = FNAME.match(f.name)
        if not m:
            raise ValueError(f"unparsed filename: {f.name}")
        cols[f"{m.group(1)}_{m.group(2)}"] = pd.read_csv(f, sep="\t", index_col=0).iloc[:, 0]
    counts = pd.DataFrame(cols)
    if counts.shape[1] != 80:
        raise ValueError(f"expected 80 libraries, found {counts.shape[1]}")
    cpm = counts / counts.sum(axis=0) * 1e6
    alleles = sorted({c.rsplit("_", 1)[0] for c in counts.columns} - {"WT", "GFP"})
    return cpm, alleles


def build_features(alleles: list[str]) -> dict[str, np.ndarray]:
    z = np.load(FEAT / f"tp63_esm2_{ISOFORM}.npz", allow_pickle=True)
    order = [str(v) for v in z["variants"]]
    idx = [order.index(a) for a in alleles]

    theta, pos, region = [], [], []
    for a in alleles:
        wt, p_label, mut = a[0], int(a[1:-1]), a[-1]
        p_ta = p_label + LABEL_TO_TA
        theta.append(compute_theta(wt, mut, p_ta, 680, domain_ranges=DOMAINS,
                                   catalytic_residues=None, hotspot_positions=None))
        pos.append([float(p_ta)])
        region.append([1.0 if 170 <= p_ta <= 362 else 0.0])

    F = {"esm2_sitedelta": z["esm2_sitedelta"][idx],
         "esm2_window16": z["esm2_window16"][idx],
         "esm2_globaldelta": z["esm2_globaldelta"][idx],
         "theta": np.asarray(theta, float),
         "position": np.asarray(pos, float),
         "region": np.asarray(region, float)}
    F["esm2_sitedelta+theta"] = np.hstack([F["esm2_sitedelta"], F["theta"]])
    return F


# ---------------------------------------------------------------- fold precompute
def rows_center(A):
    return A - A.mean(axis=1, keepdims=True)


def rows_unit(A):
    A = rows_center(A)
    return A / (np.linalg.norm(A, axis=1, keepdims=True) + 1e-12)


def pds_from_scores(s: np.ndarray, own: int) -> tuple[float, bool]:
    """Tie-aware PDS and top-1 for one similarity row. Chance PDS = 0.5."""
    n = len(s)
    v = s[own]
    better = int((s > v).sum())
    equal = int((s == v).sum()) - 1
    rank = better + equal / 2.0 + 1
    return 1.0 - (rank - 1) / (n - 1), bool(rank == 1 and equal == 0)


def precompute(cpm: pd.DataFrame, alleles: list[str]) -> list[dict]:
    """Per-fold quantities that do not depend on the features."""
    folds = []
    x_all = np.log2(cpm + 1.0)
    for ai, a in enumerate(alleles):
        train = [t for t in alleles if t != a]
        libs = [f"WT_{r}" for r in (1, 2, 3, 4)]
        libs += [f"{t}_{r}" for t in train for r in (1, 2, 3, 4)]
        keep = (cpm[libs] >= MIN_CPM).mean(axis=1) >= MIN_FRAC
        x = x_all.loc[keep]

        def delta(names, reps):
            base = x[[f"WT_{r}" for r in reps]].mean(axis=1).to_numpy()
            return np.stack([x[[f"{n}_{r}" for r in reps]].mean(axis=1).to_numpy() - base
                             for n in names])

        D_train = delta(train, (1, 2, 3, 4))
        Dc = D_train - D_train.mean(axis=0, keepdims=True)
        U = np.linalg.svd(Dc.T, full_matrices=False)[0][:, :N_PC]      # G x 5, train only
        if REMOVE_MEAN:
            # The frozen protocol takes PCs of the ALLELE-CENTRED delta matrix, which
            # removes the top axes of variation but leaves the grand mean response in.
            # Every residual then still shares that offset (observed pairwise r ~ 0.84),
            # and any predictor that is a convex combination of other alleles' residuals
            # is pushed away from the held-out one, so PDS is not centred on 0.5 for the
            # model arm. Adding the training-mean direction to the removed subspace is
            # what "remove the shared axes" was meant to do and restores chance = 0.5.
            m = D_train.mean(axis=0)
            m = m - (m @ U) @ U.T
            nrm = np.linalg.norm(m)
            if nrm > 1e-9:
                U = np.hstack([U, (m / nrm)[:, None]])                 # G x 6

        def proj(A):
            return A - (A @ U) @ U.T

        Yc = rows_center(proj(D_train))
        K = Yc @ Yc.T

        Ms, cp, ct = [], [], []
        for hA, hB in SPLITS:
            poolB = rows_unit(proj(delta(alleles, hB)))
            Ms.append(Yc @ poolB.T)
            rA = rows_unit(proj(delta([a], hA)))[0]
            p, t1 = pds_from_scores(poolB @ rA, ai)
            cp.append(p)
            ct.append(t1)

        folds.append({"allele": a, "index": ai, "train": train, "n_genes": int(keep.sum()),
                      "K": K, "M": Ms, "Kd": np.sqrt(np.clip(np.diag(K), 1e-12, None)),
                      "ceiling_pds": float(np.mean(cp)), "ceiling_top1": float(np.mean(ct))})
    return folds


# --------------------------------------------------------------------- estimators
def reduce_block(Xtr: np.ndarray, xte: np.ndarray, n_comp: int):
    """Standardise on train only, then train-only PCA."""
    mu = Xtr.mean(axis=0)
    sd = Xtr.std(axis=0)
    sd[sd < 1e-12] = 1.0
    A, b = (Xtr - mu) / sd, (xte - mu) / sd
    if A.shape[1] <= n_comp:
        return A, b
    m = A.mean(axis=0)
    V = np.linalg.svd(A - m, full_matrices=False)[2][:n_comp].T
    return (A - m) @ V, (b - m) @ V


def _inner_corr_many(C: np.ndarray, K: np.ndarray, Kd: np.ndarray, j: int) -> np.ndarray:
    """corr(prediction, target j) for a whole stack of smoothers at once."""
    CK = C @ K
    den = np.sqrt(np.clip(np.einsum("ai,ai->a", CK, C), 1e-12, None))
    return CK[:, j] / (den * Kd[j])


def _inner_corr(c: np.ndarray, K: np.ndarray, Kd: np.ndarray, j: int) -> float:
    return float(_inner_corr_many(c[None, :], K, Kd, j)[0])


def c_ridge(Xtr, xte, K, Kd):
    """Ridge; penalty chosen by inner LOO inside the training set.

    One eigendecomposition per inner fold serves the whole alpha grid and every
    alpha is scored as a batch. That is what makes 2,000 permutations of the full
    leave-one-out loop cost minutes instead of hours.
    """
    n = Xtr.shape[0]
    na = len(ALPHAS)
    scores = np.zeros(na)
    for j in range(n):
        tr = np.array([i for i in range(n) if i != j])
        mu = Xtr[tr].mean(axis=0)
        Aj = Xtr[tr] - mu
        lam, V = np.linalg.eigh(Aj @ Aj.T)
        z = V.T @ (Aj @ (Xtr[j] - mu))
        W = V @ (z[:, None] / (lam[:, None] + ALPHAS[None, :]))          # (n-1) x na
        C = np.zeros((na, n))
        C[:, tr] = (W + (1.0 - W.sum(axis=0)) / len(tr)).T
        scores += _inner_corr_many(C, K, Kd, j)
    best_a = ALPHAS[int(np.argmax(scores))]
    mu = Xtr.mean(axis=0)
    A = Xtr - mu
    lam, V = np.linalg.eigh(A @ A.T)
    w = V @ ((V.T @ (A @ (xte - mu))) / (lam + best_a))
    return w + (1.0 - w.sum()) / n, {"alpha": float(best_a),
                                     "inner_r": float(scores.max() / n)}


def c_krr(Xtr, xte, K, Kd):
    """Kernel ridge, RBF with the median-heuristic bandwidth; same batched grid."""
    n = Xtr.shape[0]
    na = len(ALPHAS)
    d2 = ((Xtr[:, None] - Xtr[None]) ** 2).sum(-1)
    pos = d2[d2 > 0]
    gamma = 1.0 / max(float(np.median(pos)) if pos.size else 1.0, 1e-9)
    Kx = np.exp(-gamma * d2)
    kx = np.exp(-gamma * ((Xtr - xte) ** 2).sum(-1))
    scores = np.zeros(na)
    for j in range(n):
        tr = np.array([i for i in range(n) if i != j])
        lam, V = np.linalg.eigh(Kx[np.ix_(tr, tr)])
        z = V.T @ Kx[tr, j]
        C = np.zeros((na, n))
        C[:, tr] = (V @ (z[:, None] / (lam[:, None] + ALPHAS[None, :]))).T
        scores += _inner_corr_many(C, K, Kd, j)
    best_a = ALPHAS[int(np.argmax(scores))]
    lam, V = np.linalg.eigh(Kx)
    return V @ ((V.T @ kx) / (lam + best_a)), {"alpha": float(best_a)}


def _pls_c(Xtr, xte, K, n_comp):
    """PLS2 in sample space.

    NIPALS needs the targets only through `K`: the first weight is the dominant
    eigenvector of `E^T K E`, and deflating the targets is `K <- P K P`. Because
    the score vectors are orthogonal, `F_h^T t_h = F_0^T t_h`, so the prediction
    collapses to `c = sum_h t*_h t_h / (t_h^T t_h)`. Targets are centred over
    SAMPLES first, which is what the trailing intercept term restores.
    """
    n = Xtr.shape[0]
    mu = Xtr.mean(axis=0)
    E, e = Xtr - mu, xte - mu
    Pc = np.eye(n) - np.ones((n, n)) / n
    Kh = Pc @ K @ Pc                       # sample-centred targets
    c = np.zeros(n)
    for _ in range(n_comp):
        Mh = E.T @ Kh @ E
        if not np.isfinite(Mh).all() or np.allclose(Mh, 0.0):
            break
        w = np.linalg.eigh(Mh)[1][:, -1]
        t = E @ w
        tt = float(t @ t)
        if tt < 1e-12:
            break
        t_new = float(e @ w)
        c = c + t * (t_new / tt)
        p = (E.T @ t) / tt
        E = E - np.outer(t, p)
        e = e - t_new * p
        D = np.eye(n) - np.outer(t, t) / tt
        Kh = D @ Kh @ D
    return c + (1.0 - c.sum()) / n         # restore the intercept


def c_pls(Xtr, xte, K, Kd):
    n = Xtr.shape[0]
    best_k, best_s = PLS_COMPONENTS[0], -np.inf
    for nc in PLS_COMPONENTS:
        s = 0.0
        for j in range(n):
            tr = np.array([i for i in range(n) if i != j])
            sub = np.ix_(tr, tr)
            c_in_sub = _pls_c(Xtr[tr], Xtr[j], K[sub], nc)
            c = np.zeros(n)
            c[tr] = c_in_sub
            s += _inner_corr(c, K, Kd, j)
        if s > best_s:
            best_s, best_k = s, nc
    return _pls_c(Xtr, xte, K, best_k), {"n_components": int(best_k)}


def c_nn(Xtr, xte):
    c = np.zeros(Xtr.shape[0])
    c[int(np.argmin(((Xtr - xte) ** 2).sum(-1)))] = 1.0
    return c


def c_group_mean(flag_tr, flag_te):
    m = flag_tr == flag_te
    c = np.zeros(len(flag_tr))
    if m.any():
        c[m] = 1.0 / m.sum()
    else:
        c[:] = 1.0 / len(flag_tr)
    return c


def score_c(c, fold, ai):
    K = fold["K"]
    den = np.sqrt(max(c @ K @ c, 1e-12))
    pds, t1, rr = [], [], []
    for M in fold["M"]:
        s = (c @ M) / den
        p, t = pds_from_scores(s, ai)
        pds.append(p)
        t1.append(t)
        rr.append(s[ai])
    return np.mean(pds), np.mean(t1), np.mean(rr)


# ------------------------------------------------------------------------- driver
def get_c(model, F, feat, idx, fold, alleles):
    ai = fold["index"]
    tr = [alleles.index(t) for t in fold["train"]]
    K, Kd = fold["K"], fold["Kd"]
    if model == "mean":
        return np.full(len(tr), 1.0 / len(tr))
    if model == "region_mean":
        r = F["region"][idx][:, 0]
        return c_group_mean(r[tr], r[ai])
    if model == "position_nn":
        p = F["position"][idx][:, 0]
        return c_nn(p[tr][:, None], np.array([p[ai]]))
    X = F[feat][idx]
    Xtr, xte = reduce_block(X[tr], X[ai], min(N_ESM_PC, len(tr) - 1))
    if model == "ridge":
        return c_ridge(Xtr, xte, K, Kd)[0]
    if model == "krr":
        return c_krr(Xtr, xte, K, Kd)[0]
    if model == "pls":
        return c_pls(Xtr, xte, K, Kd)[0]
    if model == "nn":
        return c_nn(Xtr, xte)
    raise ValueError(model)


def run_combo(folds, alleles, F, feat, model, idx=None):
    idx = np.arange(len(alleles)) if idx is None else idx
    out = np.array([score_c(get_c(model, F, feat, idx, f, alleles), f, f["index"])
                    for f in folds])
    return {"pds": out[:, 0], "top1": out[:, 1], "r": out[:, 2]}


def boot_ci(v, rng, n=N_BOOT):
    d = v[rng.integers(0, len(v), size=(n, len(v)))].mean(axis=1)
    return float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5))


def main() -> None:
    rng = np.random.default_rng(SEED)
    OUT.mkdir(parents=True, exist_ok=True)
    cpm, alleles = load_cpm()
    F = build_features(alleles)
    print(f"tag={TAG}  remove_mean={REMOVE_MEAN}  alleles={len(alleles)}  "
          f"feature dims={ {k: v.shape[1] for k, v in F.items()} }")
    t0 = time.time()
    folds = precompute(cpm, alleles)
    print(f"precompute {time.time() - t0:.1f}s   fold gene counts "
          f"{min(f['n_genes'] for f in folds)}-{max(f['n_genes'] for f in folds)}")

    ceil = np.array([f["ceiling_pds"] for f in folds])
    clo, chi = boot_ci(ceil, rng)
    C = float(ceil.mean())
    H = C - 0.5
    thr25, thr50 = 0.5 + 0.25 * H, 0.5 + 0.5 * H
    print(f"\nMEASUREMENT CEILING  PDS={C:.4f}  95% CI [{clo:.4f}, {chi:.4f}]  "
          f"top1={np.mean([f['ceiling_top1'] for f in folds]):.4f}")
    print(f"headroom H={H:.4f}   0.25H threshold={thr25:.4f}   0.50H threshold={thr50:.4f}\n")

    combos = [(f, m) for f in FEAT_SETS for m in ("ridge", "krr", "pls", "nn")]
    combos += [("theta", "mean"), ("region", "region_mean"), ("position", "position_nn")]

    rows, keep_pds = [], {}
    for feat, model in combos:
        t1 = time.time()
        res = run_combo(folds, alleles, F, feat, model)
        obs = float(res["pds"].mean())
        lo, hi = boot_ci(res["pds"], rng)
        prng = np.random.default_rng(SEED + 1)
        null = np.array([run_combo(folds, alleles, F, feat, model,
                                   idx=prng.permutation(len(alleles)))["pds"].mean()
                         for _ in range(N_PERM)])
        p = float(((null >= obs).sum() + 1) / (N_PERM + 1))
        keep_pds[(feat, model)] = res["pds"]
        rows.append({"features": feat, "model": model, "PDS": obs, "lo": lo, "hi": hi,
                     "top1": float(res["top1"].mean()), "resid_r": float(res["r"].mean()),
                     "perm_p": p, "perm_null_mean": float(null.mean()),
                     "headroom_recovery": (obs - 0.5) / H, "n_perm": N_PERM,
                     "secs": round(time.time() - t1, 1)})
        print(f"  {feat:<22}{model:<12} PDS={obs:.4f} [{lo:.4f},{hi:.4f}] "
              f"top1={rows[-1]['top1']:.3f} r={rows[-1]['resid_r']:+.4f} p={p:.4f} "
              f"recov={rows[-1]['headroom_recovery']:+.3f}  ({rows[-1]['secs']}s)")

    df = pd.DataFrame(rows).sort_values("PDS", ascending=False).reset_index(drop=True)
    df.to_csv(OUT / f"predict_gse311877_main_{TAG}.csv", index=False)
    b = df.iloc[0]

    # -------------------------------------------------- B. same-residue siblings
    sib_rows = []
    idx0 = np.arange(len(alleles))
    for p_, q_ in SIBLINGS:
        for own, sib in ((p_, q_), (q_, p_)):
            ai, si = alleles.index(own), alleles.index(sib)
            fold = folds[ai]
            c = get_c(b["model"], F, b["features"], idx0, fold, alleles)
            den = np.sqrt(max(c @ fold["K"] @ c, 1e-12))
            ro = float(np.mean([(c @ M)[ai] / den for M in fold["M"]]))
            rs = float(np.mean([(c @ M)[si] / den for M in fold["M"]]))
            sib_rows.append({"allele": own, "sibling": sib, "r_own": ro,
                             "r_sibling": rs, "beats_sibling": bool(ro > rs)})
    sib = pd.DataFrame(sib_rows)
    sib.to_csv(OUT / f"predict_gse311877_siblings_{TAG}.csv", index=False)
    k = int(sib["beats_sibling"].sum())

    # -------------------------------------------------------- verdict (binding)
    ok_ceiling = bool(clo > 0.5)
    all_below = bool((df["hi"] < thr25).all())
    all_ns = bool((df["perm_p"] > 0.05).all())
    ridge_hits = int((df[(df["model"] == "ridge")]["PDS"] >= thr25).sum())
    signal = bool(ok_ceiling and b["lo"] > 0.5 and b["perm_p"] < 0.05
                  and b["PDS"] >= thr25 and ridge_hits >= 2)
    limited = bool(ok_ceiling and all_below and all_ns)
    verdict = ("MODEL_SIGNAL_PRESENT" if signal else
               "MODEL_LIMITED_SUPPORTED" if limited else "INCONCLUSIVE_UNDERPOWERED")

    sd = float(keep_pds[(b["features"], b["model"])].std(ddof=1))
    mde = 0.5 + (1.7396 + 0.8633) * sd / np.sqrt(len(alleles))   # t(.95,17)+t(.80,17)

    rep = {"isoform": ISOFORM, "remove_mean": REMOVE_MEAN, "tag": TAG, "n_alleles": len(alleles),
           "ceiling_pds": C, "ceiling_ci": [clo, chi], "headroom": H,
           "thr_25H": thr25, "thr_50H": thr50,
           "best": {str(k2): (v.item() if hasattr(v, "item") else v) for k2, v in b.items()},
           "ceiling_lo_gt_chance": ok_ceiling,
           "all_models_hi_below_25H": all_below,
           "all_models_perm_ns": all_ns,
           "sibling_beats": k, "sibling_n": 10,
           "mde_80pct_power": float(mde), "mde_below_50H": bool(mde < thr50),
           "verdict": verdict, "n_perm": N_PERM, "n_boot": N_BOOT, "seed": SEED}
    (OUT / f"predict_gse311877_verdict_{TAG}.json").write_text(json.dumps(rep, indent=2))

    print("\n" + "=" * 80)
    print(df.to_string(index=False))
    print(f"\n[B] same-residue siblings, {b['features']}/{b['model']}: "
          f"{k}/10 beat their sibling (chance 5/10)")
    print(sib.to_string(index=False))
    print("\n" + "=" * 80)
    print(f"ceiling CI lower > 0.5             : {ok_ceiling}   (CI [{clo:.4f},{chi:.4f}])")
    print(f"every model CI upper < {thr25:.4f}   : {all_below}")
    print(f"every model permutation p > 0.05   : {all_ns}")
    print(f"min detectable PDS @80% power n=18 : {mde:.4f}  below 0.5+0.5H="
          f"{thr50:.4f}? {mde < thr50}")
    print(f"\nVERDICT: {verdict}")
    print("=" * 80)


if __name__ == "__main__":
    main()
