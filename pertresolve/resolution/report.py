"""Report independent measurement-resolution axes for a perturbation dataset.

The public report keeps four questions separate:

    detection                       each perturbation against the control
    identification                  each perturbation against its closest competitor
    split-half reproducibility      an empirical reference from disjoint measurements
    model-ranking resolution        recovery of a known graded predictor ordering

These quantities answer different empirical questions. Detection is not a logical
precondition for identification, the split-half value is not a hard bound, and no single
reference value is used to label a dataset as suitable for model ranking.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

from .profiles import GroupedProfiles, group_profiles
from .recovery import RecoveryResult, graded_predictor_scores, ordering_recovery
from .scaling import (
    cross_seed_separations,
    gram_of,
    nearest_neighbour_rho2,
    signal_noise,
    stats_from_gram,
    summarise_separations,
    tie_aware_pds,
)
from .window import WindowResult, detection_window

__all__ = ["ResolutionReport", "resolution_report"]

@dataclass
class ResolutionReport:
    """What a dataset supports, and the measurements behind it.

    Attributes:
        detection_fraction: share of perturbations whose effect clears their split-half noise.
        identification_fraction: share whose nearest-competitor separation exceeds its own
            95% spread across cell splits, the same margin detection is judged on. NaN when
            fewer than three splits were run, since the spread is then not estimable.
        split_half_reproducibility_reference: discrimination between disjoint measurements
            of the same perturbation. It is an empirical reference for interpreting model
            scores, not a hard ceiling or bound.
        rho2: mean squared separation over pairs, relative to noise. Reported because
            Methods names it, and known not to govern discrimination on its own.
        rho2_nn_median: median nearest-competitor separation relative to noise. This is the
            quantity that is useful for interpreting model-ranking resolution.
        eta2: per-profile sampling noise at the depth used.
        n_perturbations: perturbations scored.
        depth: cells per group.
        excluded: perturbations left out, with the reason.
        window: the per-perturbation detection result, if computed.
        model_ranking_resolution: the graded-predictor result, if computed.
        notes: non-gating notes about unavailable or non-estimable axes.
        details: everything else, for callers that want the raw numbers.
    """

    detection_fraction: float
    identification_fraction: float
    split_half_reproducibility_reference: float
    rho2: float
    rho2_nn_median: float
    eta2: float
    n_perturbations: int
    depth: int
    excluded: dict[str, str] = field(default_factory=dict)
    window: WindowResult | None = None
    model_ranking_resolution: RecoveryResult | None = None
    notes: tuple[str, ...] = ()
    details: dict[str, Any] = field(default_factory=dict)

    def summary(self) -> str:
        """A short report with the four axes kept independent."""
        detection = (f"{self.detection_fraction:6.1%}"
                     if np.isfinite(self.detection_fraction) else "not computed")
        identification = (f"{self.identification_fraction:6.1%}"
                          if np.isfinite(self.identification_fraction) else "not computed")
        lines = [
            "resolution report",
            f"  {self.n_perturbations} perturbations at {self.depth} cells per group"
            + (f", {len(self.excluded)} excluded" if self.excluded else ""),
            f"  detection        {detection} of perturbations clear their split-half noise",
            f"  identification   {identification} are separable from their closest competitor",
            f"  split-half ref   {self.split_half_reproducibility_reference:.3f} "
            "empirical reproducibility reference",
            f"  rho2         {self.rho2:+.4f} over pairs, {self.rho2_nn_median:+.4f} nearest competitor",
        ]
        lines.extend(f"  note: {note}" for note in self.notes)
        if self.model_ranking_resolution is not None:
            lines.append("  model-ranking resolution  "
                         f"{self.model_ranking_resolution.describe()}")
        return "\n".join(lines)

    def to_dict(self) -> dict:
        """A JSON-serialisable view, without the array-valued members."""
        out = {k: v for k, v in self.__dict__.items()
               if k not in ("window", "model_ranking_resolution", "details")}
        out["excluded"] = dict(self.excluded)
        if self.model_ranking_resolution is not None:
            ranking = self.model_ranking_resolution
            out["model_ranking_resolution"] = dict(
                p_correct_order=ranking.p_correct_order,
                p_correct_winner=ranking.p_correct_winner,
                smallest_resolved_gap=ranking.smallest_resolved_gap,
                graded_predictors_saturated=ranking.saturated,
                tied_adjacent_pairs=list(ranking.tied_adjacent_pairs),
            )
            # Scalar conveniences for analysis-table readers. These are explicitly
            # model-ranking quantities, not a combined verdict.
            out["p_correct_order"] = ranking.p_correct_order
            out["p_correct_winner"] = ranking.p_correct_winner
            out["smallest_resolved_gap"] = ranking.smallest_resolved_gap
            out["graded_predictors_saturated"] = ranking.saturated
        else:
            out["model_ranking_resolution"] = None
        out.update({k: v for k, v in self.details.items() if np.isscalar(v)})
        return out


def resolution_report(X: np.ndarray, labels, *, control: str, depth: int = 50,
                      n_seeds: int = 8, seed: int = 0, n_boot: int = 1000,
                      with_window: bool = True,
                      perturbations: list[str] | None = None,
                      candidate_mask: np.ndarray | None = None,
                      cross_seed: bool = True) -> ResolutionReport:
    """Measure what a dataset supports, before any model is compared on it.

    Args:
        X: ``(n_cells, K)`` expression or embedding matrix. Reduce dimension first for raw
            counts: distances in a very high-dimensional space are dominated by directions
            carrying no perturbation signal.
        labels: perturbation label per cell.
        control: label of the control condition.
        depth: cells per group. Four disjoint groups are cut per perturbation, so a
            perturbation needs ``4 * depth`` cells to be included.
        n_seeds: independent cell splits averaged over.
        seed: base seed.
        n_boot: bootstrap resamples for the ordering-recovery probabilities.
        with_window: also run the per-perturbation detection window, which is quadratic in
            ``depth`` and is the slowest part. The detection fraction is unavailable
            without it.
        perturbations: restrict the scored candidate set to these labels, in this order.
            ``None`` scores every non-control label, which is the released behaviour. This
            exists so that the candidate set can be varied without varying anything else:
            ``X`` is reduced once on the full cohort and only the scoring pool changes, so a
            candidate-set effect is not confounded with a change of representation basis.
        cross_seed: compute the optional two-fold cross-seed identification reference when
            at least four seeds are available. It is reported alongside the primary axes and
            does not replace them.

    Returns:
        A :class:`ResolutionReport`.

    Raises:
        ValueError: propagated from :func:`~pertresolve.resolution.profiles.group_profiles`
            when too few perturbations have enough cells.
    """
    pset = resolution_profiles(X, labels, control=control, depth=depth, n_seeds=n_seeds,
                               seed=seed, with_window=with_window,
                               perturbations=perturbations)
    return score_profiles(pset, n_boot=n_boot, seed=seed, candidate_mask=candidate_mask,
                          cross_seed=cross_seed)


@dataclass
class ResolutionProfiles:
    """The per-split profiles and noise scale, computed once so many candidate sets can be
    scored against them.

    Forming the profiles is the expensive step: it cuts four disjoint groups of cells per
    perturbation for every one of ``n_seeds`` splits. Everything downstream is a function of
    those profiles. Separating the two lets a candidate-set sweep vary only the scoring, so
    the sweep cannot change the profiles, the cells or the representation by accident, and it
    does not repeat the expensive step once per candidate set.

    Attributes:
        perturbations: scored labels, in row order.
        builds, truths, signals: per split, the prediction-side, target-side and
            signal-estimating profiles.
        stats: the output of :func:`signal_noise` over the full scored cohort. ``eta2`` is a
            within-perturbation quantity and is deliberately estimated on the whole cohort
            rather than per candidate set, so the noise scale a separation is debiased
            against does not move when the candidate set does.
        window: the per-perturbation detection result. Detection is measured against the
            control alone, so it does not depend on the candidate set and is computed once.
    """

    perturbations: list[str]
    builds: list[np.ndarray]
    truths: list[np.ndarray]
    signals: list[np.ndarray]
    stats: dict[str, Any]
    excluded: dict[str, str]
    control_cells_per_reference: int
    window: WindowResult | None
    depth: int
    n_seeds: int
    seed: int


def resolution_profiles(X: np.ndarray, labels, *, control: str, depth: int = 50,
                        n_seeds: int = 8, seed: int = 0, with_window: bool = True,
                        perturbations: list[str] | None = None) -> ResolutionProfiles:
    """Cut the cells and form the profiles. See :class:`ResolutionProfiles`."""
    X = np.asarray(X, dtype=np.float64)
    labels = np.asarray(labels)

    # Four groups: two score split-half reproducibility and graded predictors, two estimate
    # the signal.
    # Keeping them disjoint is what stops one noise realisation from setting both the
    # measurement statistic and the score it is compared against.
    grouped: GroupedProfiles = group_profiles(X, labels, control=control, depth=depth,
                                              n_groups=4, seed=seed,
                                              perturbations=perturbations)
    n = len(grouped)

    ssw_acc = np.zeros(n)
    d2_acc = np.zeros((n, n))
    builds, truths, signals = [], [], []
    for s in range(n_seeds):
        g = group_profiles(X, labels, control=control, depth=depth, n_groups=4,
                           seed=seed + s, perturbations=grouped.perturbations)
        build, truth, signal = g.profiles[:, 0], g.profiles[:, 1], g.profiles[:, 2:]
        ssw, d2 = stats_from_gram(gram_of(signal), n, 2)
        ssw_acc += ssw
        d2_acc += d2
        builds.append(build)
        truths.append(truth)
        signals.append(signal)

    stats = signal_noise(ssw_acc / n_seeds, d2_acc / n_seeds, 2)

    window = None
    if with_window:
        window = detection_window(X, labels, control=control, depth=depth,
                                  n_seeds=min(n_seeds, 8), seed=seed,
                                  perturbations=grouped.perturbations)

    return ResolutionProfiles(
        perturbations=list(grouped.perturbations), builds=builds, truths=truths,
        signals=signals, stats=stats, excluded=grouped.excluded,
        control_cells_per_reference=int(grouped.control_cells_per_reference),
        window=window, depth=depth, n_seeds=n_seeds, seed=seed)


def score_profiles(pset: ResolutionProfiles, *, n_boot: int = 1000,
                   seed: int | None = None,
                   candidate_mask: np.ndarray | None = None,
                   cross_seed: bool = True) -> ResolutionReport:
    """Score one candidate set against profiles that were formed once.

    Args:
        pset: the profiles, from :func:`resolution_profiles`.
        n_boot: bootstrap resamples for the ordering-recovery probabilities.
        seed: base seed for the bootstrap. Defaults to the seed the profiles were cut with.
        candidate_mask: optional ``(n, n)`` boolean; row ``i`` names the competitors query
            ``i`` is scored against. ``None`` scores against every other perturbation, which
            is the released behaviour. The same mask governs the split-half reference, the
            nearest-competitor separation and the ordering recovery, because all three are
            candidate-set dependent and a matched analysis that masked only one of them would
            be measuring three different candidate sets at once.
        cross_seed: also compute the two-fold cross-seed estimator, which fixes each unit's
            competitor from one half of the seeds and evaluates it on the other. It is
            reported beside the released quantity, never in place of it. Needs at least two
            seeds per fold.
    """
    seed = pset.seed if seed is None else seed
    stats = pset.stats
    n = len(pset.perturbations)

    split_half_scores = [tie_aware_pds(b, t, candidate_mask=candidate_mask)
                         for b, t in zip(pset.builds, pset.truths)]
    alpha_scores = [graded_predictor_scores(b, t, candidate_mask=candidate_mask)
                    for b, t in zip(pset.builds, pset.truths)]

    # Nearest-competitor separations are formed per split against that split's own noise and
    # then averaged. Forming them from averaged profiles over-subtracts, because averaged
    # profiles carry a fraction of the noise the debiasing term describes.
    per_split_nn = [nearest_neighbour_rho2(s, stats["eta2"], candidate_mask=candidate_mask)
                    for s in pset.signals]
    sep_by_split = np.stack([r["separations"] for r in per_split_nn])   # (n_seeds, n)
    seps = sep_by_split.mean(axis=0)
    nn = summarise_separations(seps, stats["eta2"])

    # Identification is a separate empirical axis. It uses the nearest-competitor
    # separation and its split-to-split spread; it is not gated on the control comparison.
    if len(sep_by_split) >= 3:
        spread = (np.percentile(sep_by_split, 97.5, axis=0)
                  - np.percentile(sep_by_split, 2.5, axis=0))
        # Keep the per-perturbation identification result, not only its mean, so callers
        # can compare the independent axes without re-running the whole measurement.
        identified = seps > spread
        identifiable = float(identified.mean())
    else:
        # With too few splits the spread is not estimable; say so rather than fall back to
        # a weaker test that would read as the same quantity.
        spread = np.full(len(seps), np.nan)
        identified = np.full(len(seps), False)
        identifiable = float("nan")

    # rho2_nn per unit, the continuous local evidence the profile reports as a median and as
    # a fraction above 1, where 1 is the sampling-noise scale.
    rho2_nn = seps / (2.0 * stats["eta2"]) if stats["eta2"] > 0 else np.full_like(seps, np.nan)

    cross = {}
    if cross_seed and pset.n_seeds >= 4:
        half = pset.n_seeds // 2
        folds = (tuple(range(half)), tuple(range(half, pset.n_seeds)))
        cs = cross_seed_separations(pset.signals, stats["eta2"], folds=folds,
                                    candidate_mask=candidate_mask)
        cs_mean = cs.mean(axis=0)
        cs_spread = (np.percentile(cs, 97.5, axis=0) - np.percentile(cs, 2.5, axis=0))
        cs_identified = cs_mean > cs_spread
        cs_rho2 = cs_mean / (2.0 * stats["eta2"]) if stats["eta2"] > 0 else np.full_like(cs_mean, np.nan)
        cross = dict(
            cross_seed_folds=folds,
            cross_seed_separations=cs_mean,
            cross_seed_spread=cs_spread,
            cross_seed_identified=cs_identified,
            cross_seed_identification_fraction=float(cs_identified.mean()),
            cross_seed_rho2_nn=cs_rho2,
            cross_seed_rho2_nn_median=float(np.median(cs_rho2)),
            cross_seed_p_rho2_nn_gt1=float((cs_rho2 > 1.0).mean()),
        )

    split_half_reference = float(np.mean(split_half_scores))
    recovery = ordering_recovery(np.mean(alpha_scores, axis=0), n_boot=n_boot, seed=seed)

    window = pset.window
    detectable = window.rankable_fraction if window is not None else float("nan")

    notes = []
    if window is None:
        notes.append("detection was not computed")
    if not np.isfinite(identifiable):
        notes.append("identification needs at least three cell splits to estimate its spread")

    return ResolutionReport(
        detection_fraction=detectable,
        identification_fraction=identifiable,
        split_half_reproducibility_reference=split_half_reference,
        rho2=stats["rho2"],
        rho2_nn_median=nn["rho2_nn_median"],
        eta2=stats["eta2"],
        n_perturbations=n,
        depth=pset.depth,
        excluded=pset.excluded,
        window=window, model_ranking_resolution=recovery,
        notes=tuple(notes),
        details=dict(delta2=stats["delta2"], pair2=stats["pair2"],
                     n_nonpositive=nn["n_nonpositive"],
                     control_cells_per_reference=pset.control_cells_per_reference,
                     n_groups=4,
                     n_seeds=pset.n_seeds,
                     with_window=window is not None,
                     window_n_seeds=(window.n_seeds if window is not None else None),
                     cross_seed_requested=bool(cross_seed),
                     cross_seed_computed=bool(cross),
                     # Per-perturbation identification, aligned with
                     # ``grouped.perturbations``. This lets callers compare the independent
                     # detection and identification vectors without re-running.
                     perturbations=list(pset.perturbations),
                     separations=seps,
                     separation_spread=spread,
                     identified=identified,
                     rho2_nn=rho2_nn,
                     p_rho2_nn_gt1=float((rho2_nn > 1.0).mean()),
                     **cross),
    )
