"""Properties of the PDS scorer that the rest of the study assumes and never checked.

``pertresolve/evaluation/pds.py`` had no test. Three of its properties are load-bearing
well beyond the function itself:

  * **Chance is 0.5 at every candidate-set size.** Every "at chance" verdict in the paper
    compares an interval against 0.50, and the four datasets rank against pools of 26 to
    255 competitors. If the mid-rank normalisation drifted with pool size, those verdicts
    would not be comparable across genes.
  * **Uniform candidate subsampling does not move the expected score.** This is the
    premise separating "the pool is smaller" from "the pool is easier"; without it, no
    candidate-matched comparison between JAK1 and the other genes means anything.
    Registered as P1 in ``docs/PREREG_CANDIDATE_AND_POWER_v1.md``.
  * **A non-finite prediction is refused, not scored.** ``scripts/analysis/residual_axis.py``
    records the incident this guards: scVIDR's non-finite JAK1 predictions scored 1.020 on
    26 candidates under an unguarded implementation. Small pools make it likelier, because
    one unscorable entry is a larger share of the ranking.

The repository carries nine copies of the mid-rank formula in two algebraic forms. The
vectorised form used by ``score_definitive.py`` is restated here and pinned against the
packaged one rather than imported, because importing that module builds the candidate pool
from a compute workspace this test must not need.

Needs numpy only.
"""

from __future__ import annotations

import numpy as np
import pytest

from pertresolve.evaluation.pds import pds_score, residual_pds_score

DIM = 40

#: Prediction noise that lands PDS near 0.70 on this geometry. Chosen so both
#: candidate-size tests sit in the middle of the scale: at scale 1.2 the score
#: saturates at 1.000 and neither test can fail.
NOISE = 8.0


def _score_definitive_form(distances: np.ndarray, target_index: int) -> float:
    """The vectorised mid-rank used by ``scripts/analysis/score_definitive.py:118``."""
    n = len(distances)
    if n <= 1:
        return 0.5
    td = distances[target_index]
    less = int((distances < td - 1e-12).sum())
    eq = int((np.abs(distances - td) <= 1e-12).sum())
    return 1.0 - (less + (eq - 1) / 2.0) / (n - 1)


def _deltas(names: list[str], rng: np.random.Generator) -> dict[str, np.ndarray]:
    return {n: rng.normal(size=DIM) for n in names}


# --------------------------------------------------------------------------- ranks


def test_perfect_and_worst_rankings_hit_the_ends_of_the_scale():
    real = {f"v{i}": np.eye(DIM)[i] for i in range(5)}
    cands = list(real)
    assert pds_score(real["v0"], "v0", real, cands) == pytest.approx(1.0)
    # A prediction that is exactly some other candidate puts the target last among the
    # four it is orthogonal to, and those four tie, so the target takes the mid-rank of
    # the tied block rather than rank 4.
    worst = pds_score(real["v1"], "v0", real, cands)
    assert 0.0 <= worst < 0.5


def test_a_fully_tied_ranking_scores_exactly_chance():
    """Every candidate equidistant means the score must not depend on dict order."""
    real = {f"v{i}": np.ones(DIM) for i in range(7)}
    assert pds_score(np.ones(DIM), "v3", real, list(real)) == pytest.approx(0.5)


def test_the_two_live_mid_rank_forms_agree_including_on_ties():
    """``pds.py``'s averaged-index form and ``score_definitive.py``'s counting form."""
    rng = np.random.default_rng(0)
    for trial in range(200):
        n = int(rng.integers(2, 12))
        names = [f"v{i}" for i in range(n)]
        real = _deltas(names, rng)
        if trial % 3 == 0:                      # force ties into a third of the cases
            real[names[1]] = real[names[0]].copy()
        target = names[int(rng.integers(0, n))]
        pred = rng.normal(size=DIM)
        got = pds_score(pred, target, real, names)
        d = np.array([1.0 - np.dot(pred, real[m]) /
                      (np.linalg.norm(pred) * np.linalg.norm(real[m])) for m in names])
        assert got == pytest.approx(_score_definitive_form(d, names.index(target)), abs=1e-12)


# ------------------------------------------------------------------- chance and size


@pytest.mark.parametrize("k", [2, 5, 10, 26, 50, 100])
def test_chance_is_one_half_at_every_candidate_set_size(k):
    """Random matching must centre on 0.50 whether the pool is 2 or 100 wide."""
    rng = np.random.default_rng(1234 + k)
    names = [f"v{i}" for i in range(k)]
    scores = [pds_score(rng.normal(size=DIM), names[0], _deltas(names, rng), names)
              for _ in range(3000)]
    assert np.mean(scores) == pytest.approx(0.5, abs=0.02)


def test_uniform_candidate_subsampling_leaves_the_expected_score_alone():
    """P1: shrinking the pool at fixed geometry changes the variance, not the mean.

    A signal-carrying predictor is scored against the full pool and against uniform
    subsamples of it that always retain the target. PDS is the share of competitors
    ranked worse than the target, so a uniform subsample estimates the same share; only
    its resolution changes. If this ever fails, a candidate-matched comparison between a
    26-variant and a 255-variant pool is meaningless.
    """
    rng = np.random.default_rng(7)
    n_full, n_trials = 120, 600
    names = [f"v{i}" for i in range(n_full)]
    means = {}
    for k in (5, 10, 25, 60, n_full):
        scores = []
        for _ in range(n_trials):
            real = _deltas(names, rng)
            target = names[0]
            pred = real[target] + rng.normal(scale=NOISE, size=DIM)
            others = list(rng.permutation(names[1:])[:k - 1])
            scores.append(pds_score(pred, target, real, [target] + others))
        means[k] = float(np.mean(scores))
    # A saturated geometry would pass this test without testing anything: at PDS 1.0 there
    # is no room for a drift to show. Keep the signal in the middle of the scale.
    assert 0.60 < means[n_full] < 0.85, (
        f"the geometry is saturated or empty, so invariance is untestable here: {means}")
    spread = max(means.values()) - min(means.values())
    assert spread < 0.05, f"expected PDS moved with candidate-set size: {means}"


def test_variance_falls_as_the_candidate_set_grows():
    """The counterpart of the mean invariance: small pools quantise the rank coarsely."""
    rng = np.random.default_rng(11)
    names = [f"v{i}" for i in range(80)]
    sds = {}
    for k in (5, 80):
        scores = []
        for _ in range(400):
            real = _deltas(names, rng)
            pred = real[names[0]] + rng.normal(scale=NOISE, size=DIM)
            others = list(rng.permutation(names[1:])[:k - 1])
            scores.append(pds_score(pred, names[0], real, [names[0]] + others))
        sds[k] = float(np.std(scores))
    assert sds[5] > sds[80]


# ------------------------------------------------------------------------- refusals


def test_a_non_finite_prediction_is_refused_rather_than_scored():
    real = {f"v{i}": np.eye(DIM)[i] for i in range(4)}
    for bad in (np.nan, np.inf, -np.inf):
        pred = np.eye(DIM)[0].copy()
        pred[3] = bad
        with pytest.raises(ValueError, match="non-finite"):
            pds_score(pred, "v0", real, list(real))


def test_a_zero_prediction_scores_chance_rather_than_ranking_arbitrarily():
    real = {f"v{i}": np.eye(DIM)[i] for i in range(4)}
    assert pds_score(np.zeros(DIM), "v0", real, list(real)) == pytest.approx(0.5)


def test_single_candidate_is_chance_not_a_division_by_zero():
    real = {"v0": np.ones(DIM)}
    assert pds_score(np.ones(DIM), "v0", real, ["v0"]) == pytest.approx(0.5)


# ------------------------------------------------------------------------- residual


def test_the_residual_axis_removes_the_shared_component_it_is_given():
    """Predicting only the gene-shared programme scores chance on the residual axis."""
    rng = np.random.default_rng(3)
    shared = rng.normal(size=DIM) * 5.0
    names = [f"v{i}" for i in range(20)]
    real = {n: shared + rng.normal(scale=0.4, size=DIM) for n in names}
    gene_mean_only = shared.copy()
    scores = [residual_pds_score(gene_mean_only, n, real, names, shared) for n in names]
    assert np.mean(scores) == pytest.approx(0.5, abs=0.02)
