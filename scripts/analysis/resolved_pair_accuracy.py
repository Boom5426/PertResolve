#!/usr/bin/env python3
"""Arm C of ``docs/PREREG_CANDIDATE_AND_POWER_v1.md``: do the models fail on the sibling
pairs the measurement resolves?

Every allele-level score in the study ranks a prediction against a whole candidate pool,
so it inherits an objection about how large that pool is: JAK1 ranks against 26 competitors
and GATA1 against 255. This analysis removes the objection by construction rather than
arguing with it, by asking the two-variant question directly.

A predictor resolves the pair (i, j) when it assigns its two predictions to the two
measured responses the right way round:

    d(pred_i, truth_i) + d(pred_j, truth_j) < d(pred_i, truth_j) + d(pred_j, truth_i)

Chance is 0.5 whatever the gene's pool size, and the pair is the unit throughout. The
statistic is the symmetric margin already used for the TP63 screen in
``scripts/dataset_screen/replicate_reference_gse311877.py``.

The measurement arm answers the same question with one disjoint cell half standing in for
the prediction, so both arms are scored against **the same half-depth truth**. Handing the
models the full-depth pseudobulk while the measurement arm sees only a noisy half would
flatter the models, so it is not done. The measurement arm keeps one advantage that is the
point of the comparison rather than a flaw: its prediction for a variant is a measurement
of that variant, whereas a model has never seen it. That is the same asymmetry the
manuscript's split-half reference carries.

**Selection.** A pair is called resolvable when its between-variant distance exceeds the
mean of the two variants' own split-half distances, the criterion behind Fig. 4d. Choosing
pairs and scoring the measurement on the same cell splits inflates the measurement arm: on
GATA1 it scores 0.99 on the pairs the criterion selects and 0.61 on the rest, so the
criterion is partly selecting splits where the noise fell favourably. Both versions are
therefore reported:

  ``resolvable``       mask and score from all seeds. Reproduces the published fractions
                       and is the stratum comparable to Fig. 4d.
  ``resolvable_xfit``  mask from one block of seeds, score from a disjoint block. This is
                       the version a measurement-versus-model gap should be read from.

**Chance.** A pair's two halves partition its cells, and the truth and the prediction are
referenced to complementary halves of the same control pool, so chance is not 0.5 by
assumption and is measured instead: one variant's cells are divided into two disjoint
sub-populations, each split into halves exactly as the real arm does, and the two are
scored against each other as if they were siblings. Both the plain and the WT-referenced
version land within their intervals of 0.5 for all four genes, so the statistic is
unbiased. The null runs at half the real arm's depth, because four blocks are cut from the
cells the real arm cuts two from, and is reported as a location check rather than a
precise chance level.

Two sanity gates run before any result is reported. The Gene-mean reference must score
exactly 0.5, since it predicts one profile for every variant and the two assignments are
then identically tied; and the all-seed resolvable fractions must reproduce the committed
``results/canonical/pairwise_resolvability.csv``.

Usage:
    python scripts/analysis/resolved_pair_accuracy.py [--base BASE] --out OUT
    python scripts/analysis/resolved_pair_accuracy.py --out OUT --genes JAK1

    --base    directory holding the per-gene arrays (env: PERTRESOLVE_DATA).
    --out     directory receiving the two tables; may not be inside results/.
    --genes   restrict to these genes, comma separated (default: all four).
    --preds   directory of per-method prediction .npz (default: <base>/preds5).
"""
from __future__ import annotations

import argparse
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from pertresolve.paths import (  # noqa: E402
    add_harness_to_path,
    reject_repo_results,
    require_inputs,
    resolve_base,
)

NSUB = 300
NSEED = 15
NBOOT = 2000
WT_TAGS = ("WT", "wt", "WT_control")

#: Disjoint seed blocks for the cross-fitted strata: the criterion is evaluated on
#: SELECT_SEEDS and both arms are scored on EVAL_SEEDS, so no cell split both chooses a
#: pair and decides whether that pair was resolved.
SELECT_SEEDS = tuple(range(0, 7))
EVAL_SEEDS = tuple(range(7, NSEED))
XFIT_STRATA = ("resolvable_xfit", "unresolvable_xfit")

#: Margins below this count as an exact tie and score 0.5, matching the tolerance the PDS
#: implementations use.
TIE_TOL = 1e-12


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base", default=None,
                    help="directory holding the per-gene arrays (env: PERTRESOLVE_DATA)")
    ap.add_argument("--out", required=True, metavar="OUT_DIR",
                    help="directory receiving resolved_pair_table.csv and "
                         "resolved_pair_accuracy.csv; may not be inside results/")
    ap.add_argument("--genes", default=None,
                    help="comma-separated subset of genes to score (default: all four)")
    ap.add_argument("--preds", default=None,
                    help="directory of per-method prediction .npz (default: <base>/preds5)")
    return ap.parse_args()


def unit(M: np.ndarray) -> np.ndarray:
    """Row-wise L2 normalisation, matching the cosine PDS space."""
    return M / (np.linalg.norm(M, axis=-1, keepdims=True) + 1e-12)


def halves(H, X, lab, bench_gene, seed: int, gene: str):
    """Disjoint split-half pseudobulk deltas, verbatim from pairwise_resolvability.py.

    The wild-type population is split once per (gene, seed) and shared across variants.
    Redrawing it per variant makes each prediction and its own target use complementary
    halves of the same finite control pool, which are anti-correlated, and the score
    collapses far below chance; see replicate_reference.py.
    """
    rng = np.random.RandomState(1000 + 7 * seed + H._GENE_SEED[gene])
    wt = np.where(np.isin(lab, WT_TAGS))[0]
    if len(wt) == 0:
        wm_t = wm_p = np.zeros(X.shape[1], np.float32)
    else:
        rng.shuffle(wt)
        h = len(wt) // 2
        wm_t, wm_p = X[wt[:h][:NSUB]].mean(0), X[wt[h:][:NSUB]].mean(0)
    dA, dB = {}, {}
    for v in np.unique(lab):
        if v in WT_TAGS or v not in bench_gene:
            continue
        idx = np.where(lab == v)[0]
        if len(idx) < 5:
            continue
        rng.shuffle(idx)
        h = len(idx) // 2
        dA[v] = (X[idx[:h][:NSUB]].mean(0) - wm_t).astype(np.float32)
        dB[v] = (X[idx[h:][:NSUB]].mean(0) - wm_p).astype(np.float32)
    return dA, dB


def forced_choice_matrix(pred: np.ndarray, truth: np.ndarray) -> np.ndarray:
    """Forced-choice outcome for every pair at once.

    With cosine distance on normalised rows the ``1 -`` cancels and the correct-versus-swap
    test reduces to a comparison of dot products, so the whole pair matrix is one product:

        S = pred @ truth.T
        margin[i, j] = S[i, i] + S[j, j] - S[i, j] - S[j, i] > 0

    Returns a symmetric matrix of 1.0 (resolved), 0.5 (tied) and 0.0 (swapped). Written as
    a matrix product rather than a loop over pairs because GATA1 has 32,131 of them.
    """
    S = unit(pred) @ unit(truth).T
    d = np.diag(S)
    margin = d[:, None] + d[None, :] - S - S.T
    out = np.where(margin > TIE_TOL, 1.0, 0.0)
    out[np.abs(margin) <= TIE_TOL] = 0.5
    return out


def criterion(per_seed, variants, seeds):
    """Seed-averaged pair distances, split-half noise, and the resolvable mask.

    ``D_ij`` is measured between the A halves and ``D_self`` between the two halves of one
    variant, exactly as pairwise_resolvability.py defines them.
    """
    n = len(variants)
    Dpair, dself = np.zeros((n, n)), np.zeros(n)
    for s in seeds:
        dA, dB = per_seed[s]
        A = unit(np.stack([dA[v] for v in variants]))
        B = unit(np.stack([dB[v] for v in variants]))
        Dpair += 1 - A @ A.T
        dself += 1 - np.einsum("ij,ij->i", A, B)
    Dpair /= len(seeds)
    dself /= len(seeds)
    return Dpair, dself, Dpair > 0.5 * (dself[:, None] + dself[None, :])


def self_null_accuracy(H, X, lab, bench_gene, gene, seeds, min_cells=16,
                       wt_referenced=False) -> tuple:
    """Empirical chance level for the forced choice, from pairs that cannot differ.

    Chance is **not** 0.5 for this statistic. A variant's two halves partition its cells,
    so their deviations from the variant's own sample mean are anti-correlated, while a
    prediction and *another* variant's truth are independent. The correct assignment is
    therefore penalised relative to the swap, and the statistic sits below 0.5 whenever
    the between-variant signal is small. It is the same complementary-halves effect
    replicate_reference.py records for the control pool, here acting on the variant's own
    cells.

    The null is measured rather than assumed: one variant's cells are divided into two
    disjoint sub-populations, each is split into halves exactly as the real arm does, and
    the two are scored against each other as if they were siblings. They are the same
    variant, so any departure from 0.5 is the bias and not signal.
    """
    scores = []
    n_used = 0
    wt_all = np.where(np.isin(lab, WT_TAGS))[0]
    for v in np.unique(lab):
        if v in WT_TAGS or v not in bench_gene:
            continue
        idx0 = np.where(lab == v)[0]
        if len(idx0) < min_cells:
            continue
        n_used += 1
        for s in seeds:
            rng = np.random.RandomState(90000 + 7 * s + H._GENE_SEED[gene])
            wm_t = wm_p = 0.0
            if wt_referenced and len(wt_all):
                wt = wt_all.copy()
                rng.shuffle(wt)
                hw = len(wt) // 2
                wm_t = X[wt[:hw][:NSUB]].mean(0)
                wm_p = X[wt[hw:][:NSUB]].mean(0)
            idx = idx0.copy()
            rng.shuffle(idx)
            half = len(idx) // 2
            preds, truths = [], []
            for block in (idx[:half], idx[half:2 * half]):   # two pseudo-variants
                h = len(block) // 2
                a, b = block[:h][:NSUB], block[h:][:NSUB]
                truths.append(X[a].mean(0) - wm_t)
                preds.append(X[b].mean(0) - wm_p)
            # Two nulls are reported. Without the wild-type reference the pseudo-variants
            # differ only by sampling noise. With it, the truth and the prediction are
            # additionally referenced to complementary halves of the same finite control
            # pool, which are anti-correlated; that is the construction the real arm uses,
            # so this is the version that says where its chance level actually sits.
            m = forced_choice_matrix(np.stack(preds), np.stack(truths))
            scores.append(float(m[0, 1]))
    return (float(np.mean(scores)) if scores else np.nan, np.array(scores), n_used)


def measurement_accuracy(per_seed, variants, seeds) -> np.ndarray:
    acc = np.zeros((len(variants), len(variants)))
    for s in seeds:
        dA, dB = per_seed[s]
        acc += forced_choice_matrix(np.stack([dB[v] for v in variants]),
                                    np.stack([dA[v] for v in variants]))
    return acc / len(seeds)


def model_accuracy(per_seed, names, preds, seeds) -> np.ndarray:
    P = np.stack([preds[v] for v in names])
    acc = np.zeros((len(names), len(names)))
    for s in seeds:
        dA, _ = per_seed[s]
        acc += forced_choice_matrix(P, np.stack([dA[v] for v in names]))
    return acc / len(seeds)


def boot_ci(values: np.ndarray, rng: np.random.Generator) -> tuple[float, float]:
    """Percentile interval, resampling pairs, the unit registered for this arm."""
    if len(values) < 2:
        return (np.nan, np.nan)
    draws = rng.integers(0, len(values), size=(NBOOT, len(values)))
    means = values[draws].mean(axis=1)
    return float(np.percentile(means, 2.5)), float(np.percentile(means, 97.5))


def summarise(gene, method, split, strata_masks, acc_by_stratum, rng, *,
              n_refused=0, n_predicted=0):
    """One row per stratum. ``acc_by_stratum`` supplies the scores each stratum reads."""
    rows = []
    for label, mask in strata_masks.items():
        if mask.sum() == 0:
            continue
        vals = acc_by_stratum[label][mask]
        lo, hi = boot_ci(vals, rng)
        rows.append(dict(gene=gene, method=method, split=split, stratum=label,
                         n_pairs=int(mask.sum()), accuracy=round(float(vals.mean()), 5),
                         ci_lo=round(lo, 5) if np.isfinite(lo) else np.nan,
                         ci_hi=round(hi, 5) if np.isfinite(hi) else np.nan,
                         n_refused=n_refused, n_predicted=n_predicted))
    return rows


def pool_across_splits(per_method_pairs: dict, rng) -> list[dict]:
    """One row per (gene, method) over **distinct** pairs.

    A variant can be held out in several splits, so a pair can be scored more than once.
    Summing the per-split rows would treat those repeats as independent evidence, which is
    pseudo-replication; here a repeated pair is averaged over the splits that scored it and
    then enters the bootstrap once.
    """
    rows = []
    for (gene, method), store in sorted(per_method_pairs.items()):
        prs = sorted(store)
        acc_all = np.array([float(np.mean(store[p]["all"])) for p in prs])
        acc_xf = np.array([float(np.mean(store[p]["xfit"])) for p in prs])
        res_all = np.array([store[p]["res_all"] for p in prs])
        res_xf = np.array([store[p]["res_xfit"] for p in prs])
        strata = {"all": np.ones(len(prs), bool),
                  "resolvable": res_all, "unresolvable": ~res_all,
                  "resolvable_xfit": res_xf, "unresolvable_xfit": ~res_xf}
        acc = {k: (acc_xf if k in XFIT_STRATA else acc_all) for k in strata}
        rows += summarise(gene, method, "pooled", strata, acc, rng)
    return rows


def pool_shared_cohort(per_method_pairs: dict, meas_by_pair: dict, res_by_pair: dict,
                       rng) -> list[dict]:
    """Primary pairwise comparison: one shared cohort per gene, every arm scored on it.

    ``pool_across_splits`` lets each method keep whatever pairs it happened to cover, and
    scores the measurement arm on every pair of the gene. Two rows of that table can
    therefore describe different pair populations, which makes the measurement-versus-model
    contrast the panel claims a comparison across cohorts rather than within one. Here the
    cohort is fixed before anything is scored: the measurement-resolvable pairs that EVERY
    method with a pooled row covers. Method-specific coverage remains available as the
    ``resolvable_xfit`` rows, which become the sensitivity arm.

    A method that refused every prediction for a gene has no pooled row for it and so does
    not enter that gene's intersection; the methods that do enter, and the one that binds
    the cohort size, are printed.
    """
    rows = []
    for gene in sorted({g for g, _ in per_method_pairs}):
        stores = {m: s for (g, m), s in per_method_pairs.items() if g == gene}
        if not stores:
            continue
        shared = set.intersection(*(set(s) for s in stores.values()))
        res_x = res_by_pair[gene]
        cohort = sorted(p for p in shared if res_x[p][1])
        binder = min(stores, key=lambda m: len(stores[m]))
        print(f"{gene}: shared cohort {len(cohort)} measurement-resolvable pairs over "
              f"{len(stores)} methods; smallest per-method coverage {binder} "
              f"({len(stores[binder])} pairs of {len(shared)} shared)")
        if not cohort:
            continue
        m_eval = np.array([meas_by_pair[gene][p] for p in cohort])
        lo, hi = boot_ci(m_eval, rng)
        rows.append(dict(gene=gene, method="split-half measurement", split="pooled",
                         stratum="resolvable_xfit_shared", n_pairs=len(cohort),
                         accuracy=round(float(m_eval.mean()), 5),
                         ci_lo=round(lo, 5), ci_hi=round(hi, 5),
                         n_refused=0, n_predicted=0))
        for m, store in sorted(stores.items()):
            vals = np.array([float(np.mean(store[p]["xfit"])) for p in cohort])
            lo, hi = boot_ci(vals, rng)
            rows.append(dict(gene=gene, method=m, split="pooled",
                             stratum="resolvable_xfit_shared", n_pairs=len(cohort),
                             accuracy=round(float(vals.mean()), 5),
                             ci_lo=round(lo, 5), ci_hi=round(hi, 5),
                             n_refused=0, n_predicted=0))
    return rows


def main() -> None:
    args = parse_args()
    base = resolve_base(args.base)
    out_dir = reject_repo_results(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    require_inputs(base / "pertresolve_bench.csv")
    add_harness_to_path(base)

    import harness as H  # noqa: E402  (needs base on sys.path)

    preds_dir = Path(args.preds) if args.preds else base / "preds5"
    require_inputs(preds_dir)
    genes = args.genes.split(",") if args.genes else list(H.GENES)
    unknown = [g for g in genes if g not in H.GENES]
    if unknown:
        raise SystemExit(f"unknown genes {unknown}; known: {list(H.GENES)}")

    bench_df = pd.read_csv(base / "pertresolve_bench.csv")
    bench = {g: set(bench_df[bench_df.gene == g]["variant"]) for g in H.GENES}
    methods = {p.stem: np.load(p, allow_pickle=True)
               for p in sorted(preds_dir.glob("*.npz"))}
    print(f"{len(methods)} prediction files, {len(genes)} genes; cross-fit seeds "
          f"{SELECT_SEEDS[0]}-{SELECT_SEEDS[-1]} select, "
          f"{EVAL_SEEDS[0]}-{EVAL_SEEDS[-1]} score")

    pair_rows, acc_rows = [], []
    #: (gene, method) -> {pair -> [per-split accuracy]}, plus the strata masks. A pair can
    #: be held out in more than one split, so pooling by summing the per-split rows would
    #: count it twice; the pooled row below averages a repeated pair before bootstrapping.
    per_method_pairs: dict = {}
    #: gene -> {pair -> cross-fitted measurement accuracy} and {pair -> (res_all, res_xfit)},
    #: so the shared-cohort arm can score the measurement on exactly the model cohort.
    meas_by_pair: dict = {}
    res_by_pair: dict = {}
    rng = np.random.default_rng(0)

    for gene in genes:
        X, lab = H.load_gene(gene)
        lab = np.asarray(lab)
        per_seed = [halves(H, X, lab, bench[gene], s, gene) for s in range(NSEED)]
        variants = sorted(per_seed[0][0].keys())
        n = len(variants)
        pairs = list(combinations(range(n), 2))
        pi = np.array([i for i, _ in pairs])
        pj = np.array([j for _, j in pairs])

        Dpair, dself, mask_all = criterion(per_seed, variants, range(NSEED))
        _, _, mask_xfit = criterion(per_seed, variants, SELECT_SEEDS)
        res_all, res_xfit = mask_all[pi, pj], mask_xfit[pi, pj]

        for k, (i, j) in enumerate(pairs):
            pair_rows.append(dict(
                gene=gene, v_i=variants[i], v_j=variants[j],
                D_ij=round(float(Dpair[i, j]), 6),
                Dself_i=round(float(dself[i]), 6), Dself_j=round(float(dself[j]), 6),
                pair_noise=round(float(0.5 * (dself[i] + dself[j])), 6),
                resolvable=bool(res_all[k]), resolvable_xfit=bool(res_xfit[k])))

        def strata_for(sel_all, sel_xfit):
            return {"all": np.ones(len(sel_all), bool),
                    "resolvable": sel_all, "unresolvable": ~sel_all,
                    "resolvable_xfit": sel_xfit, "unresolvable_xfit": ~sel_xfit}

        nulls = {}
        for label, wt_ref in (("self-null (same variant)", False),
                              ("self-null (WT-referenced)", True)):
            nm, ns, nn = self_null_accuracy(H, X, lab, bench[gene], gene, range(NSEED),
                                            wt_referenced=wt_ref)
            nulls[label] = nm
            if np.isfinite(nm):
                lo, hi = boot_ci(ns, rng)
                acc_rows.append(dict(
                    gene=gene, method=label, split="all", stratum="all",
                    n_pairs=len(ns), accuracy=round(nm, 5), ci_lo=round(lo, 5),
                    ci_hi=round(hi, 5), n_refused=0, n_predicted=nn))
        null_mean = nulls["self-null (WT-referenced)"]
        null_n = len(nulls)

        m_all = measurement_accuracy(per_seed, variants, range(NSEED))[pi, pj]
        m_eval = measurement_accuracy(per_seed, variants, EVAL_SEEDS)[pi, pj]
        meas_by_pair[gene] = {pr: float(m_eval[k]) for k, pr in enumerate(pairs)}
        res_by_pair[gene] = {pr: (bool(res_all[k]), bool(res_xfit[k]))
                             for k, pr in enumerate(pairs)}
        acc_rows += summarise(
            gene, "split-half measurement", "all", strata_for(res_all, res_xfit),
            {k: (m_eval if k in XFIT_STRATA else m_all) for k in
             ("all", "resolvable", "unresolvable", *XFIT_STRATA)}, rng)

        splits = sorted({k.split("__")[1] for k in methods[next(iter(methods))].files
                         if k.startswith(f"{gene}__")})
        for name, store in methods.items():
            for split in splits:
                keyed = {v: f"{gene}__{split}__{v}" for v in variants
                         if f"{gene}__{split}__{v}" in store.files}
                usable, skipped = {}, 0
                for v, key in keyed.items():
                    vec = np.asarray(store[key], dtype=float)
                    if not np.all(np.isfinite(vec)):
                        skipped += 1      # refused, never scored; reported as coverage
                        continue
                    usable[v] = vec
                if len(usable) < 2:
                    # A method whose predictions are all refused must still appear, or a
                    # total failure reads as an absent row. scVIDR's JAK1 predictions are
                    # non-finite for all 59 keys and would otherwise vanish here.
                    if keyed:
                        acc_rows.append(dict(
                            gene=gene, method=name, split=split, stratum="all",
                            n_pairs=0, accuracy=np.nan, ci_lo=np.nan, ci_hi=np.nan,
                            n_refused=skipped, n_predicted=len(keyed)))
                    continue
                names_u = [v for v in variants if v in usable]
                pos = {v: k for k, v in enumerate(names_u)}
                keep = np.array([variants[i] in usable and variants[j] in usable
                                 for i, j in pairs])
                ui = np.array([pos[variants[i]] for i in pi[keep]])
                uj = np.array([pos[variants[j]] for j in pj[keep]])
                a_all = model_accuracy(per_seed, names_u, usable, range(NSEED))[ui, uj]
                a_eval = model_accuracy(per_seed, names_u, usable, EVAL_SEEDS)[ui, uj]
                acc_rows += summarise(
                    gene, name, split, strata_for(res_all[keep], res_xfit[keep]),
                    {k: (a_eval if k in XFIT_STRATA else a_all) for k in
                     ("all", "resolvable", "unresolvable", *XFIT_STRATA)}, rng,
                    n_refused=skipped, n_predicted=len(keyed))
                store_pairs = per_method_pairs.setdefault((gene, name), {})
                kept_pairs = [pairs[k] for k in np.flatnonzero(keep)]
                for k, pr in enumerate(kept_pairs):
                    rec = store_pairs.setdefault(pr, {"all": [], "xfit": [],
                                                      "res_all": bool(res_all[keep][k]),
                                                      "res_xfit": bool(res_xfit[keep][k])})
                    rec["all"].append(a_all[k])
                    rec["xfit"].append(a_eval[k])

        print(f"{gene}: {n} variants, {len(pairs)} pairs, {int(res_all.sum())} resolvable "
              f"({res_all.mean():.3f}), {int(res_xfit.sum())} cross-fitted "
              f"({res_xfit.mean():.3f}); self-null plain "
              f"{nulls['self-null (same variant)']:.4f}, WT-referenced "
              f"{nulls['self-null (WT-referenced)']:.4f}")
        del X, lab, per_seed

    acc_rows += pool_across_splits(per_method_pairs, rng)
    print("\n=== shared measurement-resolvable cohorts (primary pairwise comparison) ===")
    acc_rows += pool_shared_cohort(per_method_pairs, meas_by_pair, res_by_pair, rng)
    pair_df, acc_df = pd.DataFrame(pair_rows), pd.DataFrame(acc_rows)
    pair_df.to_csv(out_dir / "resolved_pair_table.csv", index=False)
    acc_df.to_csv(out_dir / "resolved_pair_accuracy.csv", index=False)
    print(f"\nwrote {out_dir / 'resolved_pair_table.csv'} ({len(pair_df)} pairs) and "
          f"{out_dir / 'resolved_pair_accuracy.csv'} ({len(acc_df)} rows)")

    check_gates(pair_df, acc_df, genes)


def check_gates(pair_df: pd.DataFrame, acc_df: pd.DataFrame, genes: list[str]) -> None:
    """Registered sanity gates. Failures stop the run rather than being tolerated."""
    print("\n=== sanity gates ===")
    ok = True

    refused = acc_df[acc_df.n_refused > 0][["gene", "method", "split", "n_predicted",
                                            "n_refused", "n_pairs"]].drop_duplicates()
    if refused.empty:
        print("NOTE  every method produced finite predictions for every scored variant")
    else:
        print("NOTE  refused predictions, excluded from scoring and reported as coverage:")
        print(refused.to_string(index=False))

    gm = acc_df[(acc_df.method == "Gene-mean") & (acc_df.stratum == "all")]
    if gm.empty:
        print("SKIP  Gene-mean gate: no Gene-mean predictions were scored")
    else:
        worst = float((gm.accuracy - 0.5).abs().max())
        if worst <= 1e-9:
            print("PASS  Gene-mean scores exactly 0.5 on every pair, as a single profile "
                  "for all variants must")
        else:
            ok = False
            print(f"FAIL  Gene-mean departs from 0.5 by {worst:.2e}; the forced choice is "
                  "not tie-symmetric")

    committed = Path(__file__).resolve().parents[2] / "results" / "canonical" / \
        "pairwise_resolvability.csv"
    if committed.exists():
        ref = pd.read_csv(committed).set_index("gene")["frac_pairs_resolvable"]
        for gene in genes:
            if gene not in ref.index:
                continue
            got = float(pair_df[pair_df.gene == gene].resolvable.mean())
            delta = abs(got - float(ref[gene]))
            verdict = "PASS" if delta <= 0.005 else "FAIL"
            ok &= delta <= 0.005
            print(f"{verdict}  {gene} resolvable fraction {got:.3f} vs committed "
                  f"{float(ref[gene]):.3f}")
    if not ok:
        raise SystemExit("a sanity gate failed; the tables above are not to be read")


if __name__ == "__main__":
    main()
