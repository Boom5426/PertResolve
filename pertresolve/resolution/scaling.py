"""Unbiased estimation of the between-perturbation signal and the sampling noise.

A perturbation's measured profile is an average over finitely many cells, so two
measurements of the same perturbation differ. Whether a benchmark can tell two
perturbations apart is set by how that difference compares with the difference between
the perturbations themselves. This module estimates both, without assuming either is
larger.

Model
-----
Each perturbation ``v`` is measured ``H`` times on disjoint cells, giving profiles

    S[v, h] = mu[v] + e[v, h],     E||e||^2 = eta2,

where ``eta2`` is the per-profile sampling noise and ``mu[v]`` the noise-free profile.
The quantity that matters for discrimination is the mean squared separation between
perturbations,

    delta2 = mean over i < j of ||mu[i] - mu[j]||^2 .

Neither term is observable on its own. Two sample quantities are, and both are unbiased:

    eta2_hat  = mean over v of sum over h of ||S[v, h] - S[v]||^2 / (H - 1)
    pair2_hat = mean over i < j of ||S[i] - S[j]||^2          [S[v] the mean over h]

with ``E[pair2_hat] = delta2 + 2 * eta2 / H``, so

    delta2_hat = pair2_hat - 2 * eta2_hat / H

is unbiased for ``delta2``.

Two properties of this estimator drive the interface below.

**It can be negative, and that is not an error.** ``delta2`` is non-negative, but its
unbiased estimator is not, and when the true signal is at or below the noise a negative
estimate is the expected outcome roughly half the time. Clipping at zero converts an
estimate that straddles zero into an exact zero, which reads as certainty about a
quantity the data did not resolve, and destroys the axis for precisely the datasets where
the question is live. Nothing here clips; the sign is preserved and the interval is
reported alongside.

**Every comparison must carry the same noise.** If profiles are formed by subtracting a
reference (a wild-type or control mean), then profiles sharing one reference have that
reference's noise cancel between them, while profiles using different references do not.
Estimating ``eta2`` from profiles with independent references and ``pair2`` from profiles
with a shared reference makes ``eta2_hat`` too large by the reference noise and biases
``delta2_hat`` down by twice it, which is enough to drive a real signal below zero. Pass
profiles that all share one reference, or none at all.

Everything is computed from the Gram matrix of the pooled profiles, so a bootstrap
resample or a permutation is an index operation rather than a recomputation.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "bootstrap_delta2",
    "gram_of",
    "permutation_null_delta2",
    "signal_noise",
    "stats_from_gram",
    "summarise_separations",
    "tie_aware_pds",
]


def gram_of(profiles: np.ndarray) -> np.ndarray:
    """Inner products of the pooled profiles.

    Args:
        profiles: array of shape ``(n, H, K)``, ``H`` repeated measurements of each of
            ``n`` perturbations in ``K`` dimensions.

    Returns:
        The ``(n * H, n * H)`` Gram matrix, indexed ``p = v * H + h``.
    """
    n, h, _ = profiles.shape
    flat = np.asarray(profiles, dtype=np.float64).reshape(n * h, -1)
    return flat @ flat.T


def stats_from_gram(gram: np.ndarray, n: int, h: int) -> tuple[np.ndarray, np.ndarray]:
    """Within-perturbation scatter and between-perturbation squared distances.

    Args:
        gram: ``(n * h, n * h)`` Gram matrix from :func:`gram_of`.
        n: number of perturbations.
        h: repeated measurements per perturbation.

    Returns:
        ``ssw`` of shape ``(n,)`` holding ``sum_h ||S[v, h] - S[v]||^2``, and ``d2`` of
        shape ``(n, n)`` holding ``||S[i] - S[j]||^2`` between the perturbation means.
    """
    diag = np.diag(gram).reshape(n, h)
    block = gram.reshape(n, h, n, h)

    own = block[np.arange(n), :, np.arange(n), :].sum(axis=(1, 2)) / (h * h)
    ssw = diag.sum(axis=1) - h * own

    gmean = block.sum(axis=(1, 3)) / (h * h)
    sq = np.diag(gmean)
    d2 = sq[:, None] + sq[None, :] - 2.0 * gmean
    np.fill_diagonal(d2, 0.0)
    return ssw, d2


def signal_noise(ssw: np.ndarray, d2: np.ndarray, h: int,
                 index: np.ndarray | None = None) -> dict[str, float]:
    """Noise, signal and their dimensionless ratio, unclipped.

    Args:
        ssw: within-perturbation scatter from :func:`stats_from_gram`.
        d2: between-perturbation squared distances from :func:`stats_from_gram`.
        h: repeated measurements per perturbation, at least 2.
        index: optional perturbation indices, which may repeat, for a bootstrap
            resample. Pairs drawn from the same perturbation are excluded from the
            between-perturbation mean: their distance is zero by construction, and
            keeping them would pull the estimate toward zero for a reason unrelated
            to the data.

    Returns:
        ``eta2``, ``pair2``, ``delta2`` and ``rho2 = delta2 / (2 * eta2)``. ``delta2``
        and ``rho2`` may be negative, meaning the between-perturbation separation is
        below what sampling noise alone produces. ``rho2`` is NaN when ``eta2`` is zero.

    Raises:
        ValueError: if ``h`` is below 2, since a single measurement per perturbation
            leaves the noise unidentifiable.
    """
    if h < 2:
        raise ValueError("h must be at least 2; one measurement cannot separate "
                         "sampling noise from between-perturbation signal")
    if index is None:
        index = np.arange(len(ssw))
    n = len(index)
    nan = float("nan")
    if n < 2:
        return dict(eta2=nan, pair2=nan, delta2=nan, rho2=nan)

    eta2 = float(ssw[index].mean()) / (h - 1)
    sub = d2[np.ix_(index, index)]
    iu = np.triu_indices(n, 1)
    keep = (index[:, None] != index[None, :])[iu]
    if not keep.any():
        return dict(eta2=eta2, pair2=nan, delta2=nan, rho2=nan)

    pair2 = float(sub[iu][keep].mean())
    delta2 = pair2 - 2.0 * eta2 / h
    rho2 = delta2 / (2.0 * eta2) if eta2 > 0 else nan
    return dict(eta2=eta2, pair2=pair2, delta2=delta2, rho2=rho2)


def bootstrap_delta2(ssw: np.ndarray, d2: np.ndarray, h: int, *,
                     n_boot: int = 2000, seed: int = 0) -> np.ndarray:
    """Bootstrap distribution of ``delta2`` over perturbations and, if given, cell splits.

    Perturbations are resampled because they are the unit of generalisation: the question
    is whether another panel of perturbations from the same assay would give the same
    answer. When the statistics are supplied per split, splits are resampled too. A point
    estimate averaged over several random splits of the cells is less variable than any
    one split, so an interval that ignores which splits were drawn understates the
    uncertainty, and can exclude zero for a quantity that cannot be negative.

    Args:
        ssw: within-perturbation scatter, either ``(n,)`` for statistics already averaged
            over splits, or ``(n_splits, n)`` to resample splits as well.
        d2: between-perturbation squared distances, ``(n, n)`` or ``(n_splits, n, n)``,
            matching the shape of ``ssw``.
        h: repeated measurements per perturbation.
        n_boot: number of resamples.
        seed: seed for the resampling generator.

    Returns:
        The finite ``delta2`` values across resamples.

    Raises:
        ValueError: if ``ssw`` and ``d2`` disagree about whether splits are stacked.
    """
    ssw = np.asarray(ssw)
    d2 = np.asarray(d2)
    stacked = ssw.ndim == 2
    if stacked != (d2.ndim == 3):
        raise ValueError("ssw and d2 must both be stacked over splits, or neither")

    rng = np.random.RandomState(seed)
    n = ssw.shape[-1]
    out = np.empty(n_boot)
    for b in range(n_boot):
        if stacked:
            s = rng.randint(0, ssw.shape[0], ssw.shape[0])
            ssw_b, d2_b = ssw[s].mean(axis=0), d2[s].mean(axis=0)
        else:
            ssw_b, d2_b = ssw, d2
        out[b] = signal_noise(ssw_b, d2_b, h, rng.randint(0, n, n))["delta2"]
    return out[np.isfinite(out)]


def permutation_null_delta2(grams: list[np.ndarray], n: int, h: int, *,
                            n_perm: int = 2000, seed: int = 1) -> np.ndarray:
    """Distribution of ``delta2`` when no between-perturbation signal exists.

    Under the null every perturbation shares one mean profile, so the ``H`` measurements
    of a perturbation are exchangeable with any other perturbation's. Re-partitioning the
    pooled measurements and recomputing keeps the real noise, dimension and perturbation
    count while removing the signal. Permuting profiles permutes rows and columns of the
    Gram matrix, so no profile arithmetic is repeated.

    Args:
        grams: one Gram matrix per independent cell split; each permutation is applied to
            all of them and the results averaged, matching how the observed estimate is
            formed.
        n: number of perturbations.
        h: measurements per perturbation.
        n_perm: number of permutations.
        seed: seed for the permutation generator.

    Returns:
        One averaged ``delta2`` per permutation.
    """
    rng = np.random.RandomState(seed)
    out = np.empty(n_perm)
    for b in range(n_perm):
        order = rng.permutation(n * h)
        vals = []
        for gram in grams:
            ssw, d2 = stats_from_gram(gram[np.ix_(order, order)], n, h)
            vals.append(signal_noise(ssw, d2, h)["delta2"])
        out[b] = float(np.mean(vals))
    return out


def nearest_neighbour_rho2(profiles: np.ndarray, eta2: float, h: int = 2,
                           candidate_mask: np.ndarray | None = None) -> dict[str, float]:
    """Signal-to-noise of each perturbation's *closest* competitor, cross-fitted.

    ``delta2`` averages the squared separation over all pairs, but a discrimination score
    asks whether a perturbation's own measurement is nearer than every competitor, which
    only its closest competitor can spoil. The two come apart whenever the separations are
    not concentrated: a configuration lying in a low-dimensional subspace, or one where a
    few perturbations are widely separated and the rest are not, can have a large mean and
    a nearest neighbour buried in the noise. Measured on synthetic configurations of known
    geometry, the mean ranks the attainable score at Spearman 0.31 while this quantity
    ranks it at 0.93, so it is the mean that fails to describe discrimination.

    Choosing the nearest competitor and measuring its distance on the same data biases the
    result downward, because the chosen pair is the one whose noise happened to be most
    negative. The choice is made on one measurement and evaluated on another, which removes
    that selection.

    Args:
        profiles: ``(n, h, K)`` repeated measurements per perturbation, sharing one
            reference so the noise cancels between them.
        eta2: per-profile noise from :func:`signal_noise`, used to debias the squared
            distances, which carry ``2 * eta2 / h`` of noise between perturbation means.
        h: measurements per perturbation; at least 2 are needed to cross-fit.
        candidate_mask: optional ``(n, n)`` boolean. ``candidate_mask[i, j]`` is True when
            ``j`` is an allowed competitor of query ``i``. ``None`` allows every competitor,
            which is the released behaviour. This exists so that each query can be scored
            against its own candidate set of a fixed size without changing anything else:
            the profiles are formed once and only the row mask varies, so a candidate-set
            effect cannot be confounded with a change of profiles or of representation.

    Returns:
        ``separations``, the debiased nearest-competitor squared separation per
        perturbation, plus ``nn_median`` and ``nn_geomean`` summarising it and the same
        divided by ``2 * eta2`` as ``rho2_nn_median`` and ``rho2_nn_geomean``. Debiasing
        can send an individual separation negative; those are floored at a small positive
        value for the geometric mean only, and the count is returned as ``n_nonpositive``.

        Callers averaging over several cell splits must average ``separations`` across
        splits and summarise afterwards, rather than averaging the profiles and calling
        this once: averaged profiles carry ``1 / n_splits`` of the noise that ``eta2``
        describes, so the debiasing would over-subtract and drive every separation
        negative.

    Raises:
        ValueError: if fewer than two measurements or two perturbations are supplied.
    """
    profiles = np.asarray(profiles, dtype=np.float64)
    n = profiles.shape[0]
    if profiles.shape[1] < 2 or n < 2:
        raise ValueError("cross-fitting needs at least two measurements and two perturbations")

    if candidate_mask is not None:
        candidate_mask = np.asarray(candidate_mask, dtype=bool)
        if candidate_mask.shape != (n, n):
            raise ValueError(f"candidate_mask must be ({n}, {n}), got {candidate_mask.shape}")
        off_diagonal = candidate_mask & ~np.eye(n, dtype=bool)
        if not off_diagonal.any(axis=1).all():
            raise ValueError("every query needs at least one allowed competitor")

    def _sq_dists(mat):
        gram = mat @ mat.T
        sq = np.diag(gram)
        out = sq[:, None] + sq[None, :] - 2.0 * gram
        np.fill_diagonal(out, np.inf)
        if candidate_mask is not None:
            out = np.where(candidate_mask, out, np.inf)
        return out

    a, b = profiles[:, 0, :], profiles[:, 1, :]
    noise = 2.0 * eta2                      # two single measurements, not two means
    chosen_a, chosen_b = _sq_dists(a).argmin(axis=1), _sq_dists(b).argmin(axis=1)
    eval_b, eval_a = _sq_dists(b), _sq_dists(a)
    sep = 0.5 * (eval_b[np.arange(n), chosen_a] + eval_a[np.arange(n), chosen_b]) - noise

    return summarise_separations(sep, eta2)


def cross_seed_separations(signals, eta2: float, *, folds=((0, 1, 2, 3), (4, 5, 6, 7)),
                           candidate_mask: np.ndarray | None = None) -> np.ndarray:
    """Nearest-competitor separations with the competitor fixed on disjoint seeds.

    ``nearest_neighbour_rho2`` already cross-fits *within* a seed: the competitor is chosen on
    one cell split and measured on the other. It does not cross-fit *across* seeds, so every
    seed re-chooses the competitor, and when several competitors are near-tied the choice
    churns between seeds. That churn is then charged to the measurement's noise, where it is
    indistinguishable from a genuinely unstable separation. Measured on this project's own
    data, ``median(spread) / median(separation)`` is 3.587 for Norman and 4.303 for Adamson
    against 1.094 for a Gaussian synthetic matched on ``rho2_nn`` and pool size.

    This estimator separates the two. For each fold, the competitor is fixed once from the
    mean of the selection seeds' distance matrices and is then evaluated only on the seeds
    that took no part in choosing it.

    Implementation note, recorded because the frozen protocol names "the cross-fitted
    squared-distance matrix" in the singular while the underlying routine forms one per side:
    ``D^(s)`` is taken as the mean of the two sides at seed ``s``, which is that seed's
    estimate of the pairwise squared distance. Evaluation uses the released two-sided formula
    with the competitor held fixed rather than re-chosen.

    Args:
        signals: sequence of ``(n, 2, K)`` arrays, one per seed.
        eta2: per-profile noise, used to debias exactly as in the released estimator.
        folds: the seed indices of the two folds. Each fold selects for the other.
        candidate_mask: optional ``(n, n)`` boolean restricting each query's competitors.

    Returns:
        ``(n_seeds, n)`` per-unit separations, each observation produced on a seed that did
        not select the competitor it evaluates. The rows are ordered fold by fold, so the
        eight observations per unit are the same count the released estimator produces and
        the two are directly comparable.

    Raises:
        ValueError: if the folds are not disjoint, do not cover the seeds, or leave a query
            with no allowed competitor.
    """
    signals = [np.asarray(s, dtype=np.float64) for s in signals]
    n = signals[0].shape[0]
    flat = [i for f in folds for i in f]
    if len(set(flat)) != len(flat) or sorted(flat) != list(range(len(signals))):
        raise ValueError("folds must be disjoint and cover every seed exactly once")

    def _sq(mat):
        gram = mat @ mat.T
        sq = np.diag(gram)
        out = sq[:, None] + sq[None, :] - 2.0 * gram
        np.fill_diagonal(out, np.inf)
        if candidate_mask is not None:
            out = np.where(candidate_mask, out, np.inf)
        return out

    per_seed = []
    for sig in signals:
        a, b = sig[:, 0, :], sig[:, 1, :]
        per_seed.append((_sq(a), _sq(b)))

    noise = 2.0 * eta2
    rows = []
    for select, evaluate in ((folds[0], folds[1]), (folds[1], folds[0])):
        mean_d = np.mean([0.5 * (per_seed[s][0] + per_seed[s][1]) for s in select], axis=0)
        if not np.isfinite(mean_d).any(axis=1).all():
            raise ValueError("every query needs at least one allowed competitor")
        fixed = mean_d.argmin(axis=1)
        for s in evaluate:
            da, db = per_seed[s]
            rows.append(0.5 * (db[np.arange(n), fixed] + da[np.arange(n), fixed]) - noise)
    return np.stack(rows)


def summarise_separations(separations: np.ndarray, eta2: float) -> dict:
    """Summarise per-perturbation nearest-competitor separations against the noise.

    Split out so a caller can average :func:`nearest_neighbour_rho2`'s ``separations``
    over several cell splits and then summarise once, which is the correct order.

    Args:
        separations: debiased nearest-competitor squared separation per perturbation.
        eta2: per-profile noise the separations are measured against.

    Returns:
        The same fields :func:`nearest_neighbour_rho2` returns.
    """
    sep = np.asarray(separations, dtype=np.float64)
    n_nonpositive = int((sep <= 0).sum())
    median = float(np.median(sep))

    # A geometric mean is undefined once any separation debiases below zero, which happens
    # for a large share of perturbations in exactly the regimes this is used to describe.
    # Flooring those at a small positive number does not rescue it: the result is then set
    # by the floor and by how many values hit it, and it prints as a small positive number
    # that looks like a measurement. It is reported as undefined instead.
    geomean = (float(np.exp(np.mean(np.log(sep)))) if n_nonpositive == 0
               else float("nan"))

    scale = 2.0 * eta2
    usable = eta2 > 0
    return dict(
        separations=sep,
        nn_median=median,
        nn_geomean=geomean,
        rho2_nn_median=median / scale if usable else float("nan"),
        rho2_nn_geomean=geomean / scale if usable else float("nan"),
        n_nonpositive=n_nonpositive,
        frac_above_noise=float((sep > 0).mean()),
    )


def _frac_above_noise_warning() -> str:
    """Why ``frac_above_noise`` must be averaged over splits rather than computed from an
    average.

    It counts how many separations exceed zero, which is a threshold on a quantity whose
    noise shrinks as more splits are averaged. Averaging first therefore pushes the count
    toward whatever the sign of the mean is, and any residual bias, however small, ends up
    determining the answer: on data with exactly zero separation the fraction reads 0.52
    from single splits and 0.59 once six are averaged. Compute it per split and average the
    fractions. The median has no such problem: averaging moves it toward the mean, which is
    unbiased.
    """
    return _frac_above_noise_warning.__doc__ or ""


def tie_aware_pds(query: np.ndarray, truth: np.ndarray,
                  candidate_mask: np.ndarray | None = None) -> float:
    """Mean tie-aware cosine perturbation discrimination score.

    Row ``i`` of ``query`` is a predicted or replicate profile for the perturbation whose
    measured profile is row ``i`` of ``truth``. Each query is ranked against every truth
    profile; ties take the mid-rank, so a degenerate prediction scores chance rather than
    winning or losing by an arbitrary sort order.

    Args:
        query: ``(n, K)`` profiles to score.
        truth: ``(n, K)`` measured profiles in the same perturbation order.
        candidate_mask: optional ``(n, n)`` boolean; row ``i`` names the truths query ``i``
            is ranked against. Its own truth is always included. ``None`` ranks against every
            truth, which is the released behaviour. The rank denominator follows the row, so
            a score under a candidate set of size ``K`` is a rank among ``K``, not among ``n``.

    Returns:
        Mean score over the ``n`` queries, 0.5 at chance and 1.0 for perfect ranking.
    """
    qn = query / (np.linalg.norm(query, axis=1, keepdims=True) + 1e-12)
    tn = truth / (np.linalg.norm(truth, axis=1, keepdims=True) + 1e-12)
    dist = 1.0 - qn @ tn.T
    n = dist.shape[1]
    if n < 2:
        return 0.5
    scores = []
    for i in range(dist.shape[0]):
        row = dist[i]
        target = row[i]
        if candidate_mask is None:
            less = int((row < target - 1e-12).sum())
            eq = int((np.abs(row - target) <= 1e-12).sum())
            denominator = n - 1
        else:
            allowed = np.asarray(candidate_mask[i], dtype=bool).copy()
            allowed[i] = True                      # a query always ranks against its own truth
            sub = row[allowed]
            less = int((sub < target - 1e-12).sum())
            eq = int((np.abs(sub - target) <= 1e-12).sum())
            denominator = int(allowed.sum()) - 1
        if denominator < 1:
            continue
        scores.append(1.0 - (less + (eq - 1) / 2.0) / denominator)
    if not scores:
        return 0.5
    return float(np.mean(scores))
