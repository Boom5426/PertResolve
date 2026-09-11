"""The independent resolution axes, checked on datasets known by construction.

The estimator tests in ``test_resolution_scaling.py`` check the quantities. These check the
report built on them, and in particular the properties that a user would be misled by if
they broke. Each test names the failure it protects against, several of which were real.

Runs standalone as well as under pytest:

    python tests/test_resolution_report.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from pertresolve.resolution import (  # noqa: E402
    energy_distance,
    graded_predictor_scores,
    group_profiles,
    ordering_recovery,
    resolution_report,
)

K = 30          # profile dimension
CELLS = 320     # cells per perturbation
DEPTH = 50      # cells per group


def dataset(n_pert: int, separation: float, seed: int = 0, cells: int = CELLS):
    """Perturbations separated by a known amount, plus a control, in cell space.

    ``separation`` scales the spread of the perturbation means, so 0 gives a dataset where
    no perturbation differs from any other or from the control.
    """
    rng = np.random.RandomState(seed)
    mu = rng.normal(0, separation / np.sqrt(K), size=(n_pert, K))
    blocks, labels = [], []
    for v in range(n_pert):
        blocks.append(mu[v] + rng.normal(0, 1.0, size=(cells, K)))
        labels += [f"P{v}"] * cells
    blocks.append(rng.normal(0, 1.0, size=(cells * 3, K)))
    labels += ["ctrl"] * (cells * 3)
    return np.vstack(blocks), np.asarray(labels)


def test_a_dataset_with_no_separation_is_reported_as_having_none():
    """The failure this catches was a false positive in a headline resolution claim.

    The identification fraction was first computed by counting separations above zero after
    averaging them across cell splits. Averaging shrinks the noise around a residual bias
    rather than the bias itself, so on data with exactly zero separation the count read 60%
    instead of 50%, and the tool reported most perturbations as separable when none were.
    """
    X, labels = dataset(16, separation=0.0)
    report = resolution_report(X, labels, control="ctrl", depth=DEPTH, n_seeds=8, n_boot=200)
    assert report.detection_fraction == 0.0, "nothing should be detected here"
    assert report.identification_fraction < 0.1, (
        f"identification fraction {report.identification_fraction:.3f} on data with no "
        f"separation at all")
    assert not hasattr(report, "verdict")
    assert "resolution report" in report.summary()
    assert "split-half ref" in report.summary()
    assert "benchmarkable" not in report.summary()
    assert abs(report.rho2) < 0.2, f"rho2 {report.rho2:+.3f} should be near zero"


def test_detection_and_identification_are_independent_axes():
    """Identification is reported even when the control comparison is not rankable."""
    from dataclasses import replace
    from pertresolve.resolution.report import resolution_profiles, score_profiles

    X, labels = dataset(16, separation=1.0, seed=3)
    profiles = resolution_profiles(X, labels, control="ctrl", depth=DEPTH, n_seeds=8,
                                   seed=0, with_window=True)
    assert profiles.window is not None
    suppressed_window = replace(
        profiles.window,
        rankable=np.zeros_like(profiles.window.rankable, dtype=bool),
    )
    report = score_profiles(replace(profiles, window=suppressed_window), n_boot=200, seed=0)
    assert report.detection_fraction == 0.0
    assert report.identification_fraction > 0.5
    assert "detection" in report.summary() and "identification" in report.summary()


def test_resolution_report_forwards_cross_seed_request(monkeypatch):
    """The CLI/API switch must control computation rather than only its label."""
    import pertresolve.resolution.report as report_module

    sentinel_profiles = object()
    sentinel_report = object()
    seen = {}

    monkeypatch.setattr(
        report_module, "resolution_profiles",
        lambda *args, **kwargs: sentinel_profiles,
    )

    def fake_score_profiles(profiles, **kwargs):
        seen.update(kwargs)
        assert profiles is sentinel_profiles
        return sentinel_report

    monkeypatch.setattr(report_module, "score_profiles", fake_score_profiles)
    got = report_module.resolution_report(
        np.empty((0, 0)), [], control="ctrl", cross_seed=False,
    )
    assert got is sentinel_report
    assert seen["cross_seed"] is False


def test_the_split_half_reference_tracks_the_true_separation():
    """The empirical split-half reference remains numerically stable and interpretable."""
    references = []
    for separation in (0.0, 0.5, 1.0, 2.0):
        X, labels = dataset(16, separation=separation, seed=1)
        report = resolution_report(X, labels, control="ctrl", depth=DEPTH, n_seeds=6,
                                   n_boot=200)
        references.append(report.split_half_reproducibility_reference)
    assert all(b >= a - 1e-9 for a, b in zip(references, references[1:])), (
        f"split-half reference not monotone in the true separation: {references}")
    assert references[0] < 0.6 and references[-1] > 0.9


def test_identification_is_refused_when_the_margin_cannot_be_estimated():
    """Too few splits should say so rather than fall back to a weaker test silently."""
    X, labels = dataset(12, separation=1.0, seed=2)
    report = resolution_report(X, labels, control="ctrl", depth=DEPTH, n_seeds=2, n_boot=100)
    assert np.isnan(report.identification_fraction)
    assert any("at least three cell splits" in note for note in report.notes)


def test_a_zero_ordering_probability_is_always_explained():
    """A bare 0.00 in model-ranking resolution needs an explanation.

    Ordering probability falls to zero for two quite different reasons: the graded
    predictors are indistinguishable to the score, or the benchmark genuinely cannot order
    them. Both are legitimate outcomes and neither should be reported as a naked number.
    """
    X, labels = dataset(16, separation=8.0, seed=4)
    report = resolution_report(X, labels, control="ctrl", depth=DEPTH, n_seeds=4, n_boot=200)
    assert report.model_ranking_resolution is not None
    if report.model_ranking_resolution.p_correct_order == 0.0:
        described = report.model_ranking_resolution.describe()
        assert ("not assessable" in described or "ordered reliably" in described
                or "orders predictors differing" in described), (
            f"unexplained zero ordering probability: {described!r}")


def test_interpolating_toward_a_wrong_direction_resolves_finer_gaps_than_toward_the_mean():
    """Guards the construction the whole recovery measure rests on.

    The discrimination score is a cosine, so a prediction shrunk toward a near-zero panel
    mean stays collinear with the unshrunk one and scores identically. Measured on a
    resolvable benchmark the panel-mean weights returned 0.500 followed by six values of
    1.000, resolving nothing finer than the gap between using a perturbation's own profile
    and not using it. Interpolating toward another perturbation's profile makes the wrong
    answer wrong in direction, which the score can see.
    """
    rng = np.random.RandomState(11)
    n = 40
    mu = rng.normal(0, 1.0, size=(n, K))
    build = mu + rng.normal(0, 1.0, size=(n, K))
    truth = mu + rng.normal(0, 1.0, size=(n, K))

    mean_scores = graded_predictor_scores(build, truth, toward="panel_mean").mean(axis=1)
    flat = np.diff(mean_scores[1:])
    assert np.all(np.abs(flat) < 0.01), (
        f"the panel-mean construction is expected to flatline above weight 0: {mean_scores}")

    shuffled = graded_predictor_scores(build, truth, toward="shuffled").mean(axis=1)
    assert shuffled[1] < shuffled[-1] - 0.02, (
        f"the shuffled construction should still climb after weight 0.5: {shuffled}")

    gap_mean = ordering_recovery(
        graded_predictor_scores(build, truth, toward="panel_mean"), n_boot=300, seed=0
    ).smallest_resolved_gap
    gap_shuffled = ordering_recovery(
        graded_predictor_scores(build, truth, toward="shuffled"), n_boot=300, seed=0
    ).smallest_resolved_gap
    assert gap_shuffled < gap_mean, (
        f"shuffled should resolve a finer quality gap: {gap_shuffled} against {gap_mean}")


def test_an_unrecognised_interpolation_target_raises():
    """Silently picking a construction would change what the number means."""
    rng = np.random.RandomState(0)
    build = truth = rng.normal(size=(5, K))
    try:
        graded_predictor_scores(build, truth, toward="mean")
    except ValueError:
        return
    raise AssertionError("an unknown target should raise")


def test_groups_are_disjoint_and_perturbations_that_cannot_supply_them_are_named():
    """Silent exclusion is how a benchmark quietly changes what it measured."""
    X, labels = dataset(6, separation=1.0, seed=5, cells=CELLS)
    labels = labels.copy()
    thin = np.where(labels == "P0")[0][:CELLS - 10]
    labels[thin] = "ctrl"                       # leave P0 with only 10 cells
    grouped = group_profiles(X, labels, control="ctrl", depth=DEPTH, n_groups=4, seed=0)
    assert "P0" in grouped.excluded and "10 cells" in grouped.excluded["P0"]
    assert "P0" not in grouped.perturbations
    assert grouped.profiles.shape == (len(grouped.perturbations), 4, K)


def test_too_few_perturbations_raises_rather_than_returning_a_report():
    """A comparison needs something to compare against."""
    X, labels = dataset(3, separation=1.0, seed=6)
    try:
        group_profiles(X, labels, control="ctrl", depth=CELLS, n_groups=4, seed=0)
    except ValueError as exc:
        assert "cells" in str(exc)
        return
    raise AssertionError("should have raised when no perturbation has enough cells")


def test_energy_distance_is_zero_for_one_distribution_and_positive_for_two():
    """The distance underneath detection, checked on its own."""
    rng = np.random.RandomState(0)
    a = rng.normal(size=(200, 5))
    b = rng.normal(size=(200, 5))
    c = rng.normal(loc=3.0, size=(200, 5))
    assert abs(energy_distance(a, b)) < 0.1, "same distribution should give about zero"
    assert energy_distance(a, c) > 1.0, "a shifted distribution should be far"


def test_graded_predictors_are_ordered_by_construction_on_a_resolvable_benchmark():
    """If the benchmark can see them at all, better predictors must score higher.

    Judged at a noise level where the benchmark is informative without saturating: with
    almost no predictor noise every weight above zero reaches the ceiling and ties, which
    is a property of the predictors rather than a failure of the benchmark.
    """
    rng = np.random.RandomState(7)
    n = 40
    mu = rng.normal(0, 1.0, size=(n, K))
    build = mu + rng.normal(0, 1.0, size=(n, K))
    truth = mu + rng.normal(0, 1.0, size=(n, K))
    scores = graded_predictor_scores(build, truth)
    means = scores.mean(axis=1)
    assert means[0] < means[-1], f"the best predictor should outscore the worst: {means}"
    assert np.all(np.diff(means) > -0.01), f"scores should not fall with quality: {means}"

    # "the best predictor wins" is the wrong invariant to demand: the top weights are
    # deliberately close together, so once several of them reach the ceiling the argmax is
    # a coin flip among near-ties and the probability is low for a reason that says nothing
    # about the benchmark. What a working benchmark must do is resolve some gap.
    result = ordering_recovery(scores, n_boot=300, seed=0)
    assert np.isfinite(result.smallest_resolved_gap), (
        f"a benchmark this resolvable should order at least one adjacent pair: "
        f"P(order) = {result.p_correct_order}, scores {means}")
    assert result.smallest_resolved_gap <= 0.5


def test_exact_ties_are_reported_as_ties_and_not_as_failure():
    """P(full order) = 0 from a tie must not read as the benchmark ordering wrongly.

    ``saturated`` only fires when the *whole* range collapses. The common case at high
    resolution is that only the top predictors do, so the ranking description must explain
    the tied pair rather than imply systematic failure.
    """
    weights = (0.0, 0.5, 0.7, 0.85, 0.925, 0.96, 1.0)
    n = 40
    scores = np.tile(np.array([0.55, 0.92, 0.98, 1.0, 1.0, 1.0, 1.0])[:, None], (1, n))

    result = ordering_recovery(scores, weights=weights, n_boot=200, seed=0)

    assert not result.saturated, "the range is wide; the all-tied guard must not fire"
    assert result.p_correct_order == 0.0
    assert result.tied_adjacent_pairs == (3, 4, 5), result.tied_adjacent_pairs
    assert "score alike" in result.describe(), result.describe()

    # And a genuinely separated set must not be flagged.
    clean = np.tile(np.linspace(0.5, 0.9, len(weights))[:, None], (1, n))
    assert ordering_recovery(clean, weights=weights, n_boot=200, seed=0).tied_adjacent_pairs == ()


def test_report_serialization_names_the_independent_axes():
    X, labels = dataset(8, separation=1.0, seed=12)
    report = resolution_report(X, labels, control="ctrl", depth=DEPTH, n_seeds=4,
                               n_boot=100)
    payload = report.to_dict()
    assert {
        "detection_fraction",
        "identification_fraction",
        "split_half_reproducibility_reference",
        "model_ranking_resolution",
    }.issubset(payload)
    assert "verdict" not in payload
    assert "replicate_ceiling" not in payload


if __name__ == "__main__":
    failed = 0
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            try:
                fn()
                print(f"PASS  {name}")
            except AssertionError as exc:
                failed += 1
                print(f"FAIL  {name}: {exc}")
    print(f"\n{'all tests passed' if not failed else f'{failed} test(s) failed'}")
    sys.exit(1 if failed else 0)
