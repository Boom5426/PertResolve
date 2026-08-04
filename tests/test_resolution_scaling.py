"""The resolution estimator, checked against data whose true signal is known.

Every claim the estimator makes is about a quantity that cannot be observed directly, so
the only way to know it is right is to build data where the answer is fixed by
construction and see whether the estimate recovers it. Each test below states the property
it is protecting and the defect it would catch.

The repository has no test runner configured yet, so this file runs standalone as well as
under pytest:

    python tests/test_resolution_scaling.py

The checks are randomised but seeded, so a failure is reproducible rather than flaky. The
tolerances are stated in standard errors of the simulation, not as round numbers, so a
tightening of the estimator does not silently loosen the test.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from alleleperturb.resolution import (  # noqa: E402
    bootstrap_delta2,
    gram_of,
    nearest_neighbour_rho2,
    permutation_null_delta2,
    signal_noise,
    stats_from_gram,
    summarise_separations,
)

K = 40      # profile dimension
H = 2       # measurements per perturbation


def synth(n, delta_scale, eta, rng, *, shared_reference=True, ref_cells=None):
    """``n`` perturbations with a between-perturbation signal fixed by construction.

    ``mu[v]`` is drawn coordinate-wise with standard deviation ``delta_scale / sqrt(K)``,
    so ``E||mu_i - mu_j||^2 = 2 * delta_scale^2`` exactly. Each measurement adds noise
    with ``E||e||^2 = eta^2``. Passing ``shared_reference=False`` subtracts a different
    reference mean from each measurement, which is the defect the superseded estimator
    had, and lets a test assert that it biases the answer.
    """
    mu = rng.normal(0, delta_scale / np.sqrt(K), size=(n, K))
    S = mu[:, None, :] + rng.normal(0, eta / np.sqrt(K), size=(n, H, K))
    if ref_cells is not None:
        shape = (1, 1, K) if shared_reference else (1, H, K)
        S = S - rng.normal(0, ref_cells / np.sqrt(K), size=shape)
    return S, 2.0 * delta_scale ** 2, eta ** 2


def measure(S):
    n, h, _ = S.shape
    ssw, d2 = stats_from_gram(gram_of(S), n, h)
    return signal_noise(ssw, d2, h)


def test_delta2_is_unbiased_across_signal_regimes():
    """The signal estimate must be centred on the truth, including when the truth is zero.

    A plug-in estimate of the between-perturbation separation is biased upward by the
    sampling noise; subtracting the noise term is what removes that. If the subtraction
    were wrong, this test fails hardest in the no-signal regime, which is the regime the
    measurement-limited datasets sit in.
    """
    for delta_scale in (0.0, 0.3, 1.0, 3.0):
        rng = np.random.RandomState(0)
        vals = [measure(synth(60, delta_scale, 1.0, rng)[0])["delta2"] for _ in range(400)]
        true = 2.0 * delta_scale ** 2
        se = np.std(vals) / np.sqrt(len(vals))
        assert abs(np.mean(vals) - true) < 4 * se, (
            f"delta_scale={delta_scale}: mean {np.mean(vals):+.4f} vs true {true:.4f}")


def test_eta2_is_unbiased():
    """The noise estimate must not absorb any of the signal."""
    rng = np.random.RandomState(0)
    for delta_scale in (0.0, 3.0):
        vals = [measure(synth(60, delta_scale, 1.0, rng)[0])["eta2"] for _ in range(300)]
        assert abs(np.mean(vals) - 1.0) < 0.02, f"eta2 mean {np.mean(vals):.4f} vs 1.0"


def test_clipping_at_zero_would_bias_upward_and_erase_the_axis():
    """Guards the decision not to clip.

    On data with no signal an unbiased estimate is negative about half the time. Replacing
    those with zero both biases the mean upward and collapses the axis onto a single
    value, which is what made the superseded table unusable for the datasets of interest.
    """
    rng = np.random.RandomState(1)
    raw = np.array([measure(synth(60, 0.0, 1.0, rng)[0])["delta2"] for _ in range(400)])
    se = raw.std() / np.sqrt(len(raw))
    assert abs(raw.mean()) < 4 * se, "unclipped estimate is not centred on zero"
    assert np.maximum(raw, 0).mean() > 3 * se, "clipping should bias upward"
    assert (np.maximum(raw, 0) == 0).mean() > 0.4, "clipping should erase much of the axis"


def test_independent_reference_means_bias_the_signal_downward():
    """Reproduces the defect that put 13 of 16 real points at exactly zero.

    When the noise term is estimated from profiles with independent references but the
    signal term from profiles sharing one, the reference noise is counted only in the
    former. The predicted deficit is ``2 * eta2_ref / H``, which is 1.0 here.
    """
    rng = np.random.RandomState(2)
    shared, independent = [], []
    for _ in range(300):
        shared.append(measure(synth(60, 1.0, 1.0, rng, shared_reference=True,
                                    ref_cells=1.0)[0])["delta2"])
        independent.append(measure(synth(60, 1.0, 1.0, rng, shared_reference=False,
                                         ref_cells=1.0)[0])["delta2"])
    se = np.std(shared) / np.sqrt(len(shared))
    assert abs(np.mean(shared) - 2.0) < 4 * se, "shared reference should stay unbiased"
    deficit = 2.0 - np.mean(independent)
    assert 0.8 < deficit < 1.2, f"expected a deficit near 1.0, measured {deficit:.3f}"


def test_permutation_null_is_calibrated_and_powered():
    """The P value must be uniform when there is no signal and small when there is."""
    rng = np.random.RandomState(3)
    p_null, p_signal = [], []
    for _ in range(120):
        for delta_scale, sink in ((0.0, p_null), (1.5, p_signal)):
            S = synth(40, delta_scale, 1.0, rng)[0]
            obs = measure(S)["delta2"]
            null = permutation_null_delta2([gram_of(S)], 40, H, n_perm=200,
                                           seed=int(rng.randint(1e6)))
            sink.append((np.sum(null >= obs) + 1) / (len(null) + 1))
    p_null = np.array(p_null)
    assert (p_null < 0.05).mean() < 0.13, "null rejection rate far above nominal"
    assert 0.35 < p_null.mean() < 0.65, f"null P values not uniform, mean {p_null.mean():.3f}"
    assert (np.array(p_signal) < 0.05).mean() > 0.95, "no power against a real signal"


def test_bootstrap_interval_covers_the_truth():
    """A 95% interval over perturbations should contain the true value about 95% of the time."""
    rng = np.random.RandomState(4)
    cover = []
    for _ in range(200):
        S, true_d2, _ = synth(60, 1.0, 1.0, rng)
        n, h, _ = S.shape
        ssw, d2 = stats_from_gram(gram_of(S), n, h)
        lo, hi = np.percentile(bootstrap_delta2(ssw, d2, h, n_boot=400, seed=0), [2.5, 97.5])
        cover.append(lo <= true_d2 <= hi)
    assert 0.88 <= np.mean(cover) <= 0.99, f"coverage {np.mean(cover):.3f}"


def test_two_level_bootstrap_is_wider_than_resampling_perturbations_alone():
    """Guards the decision to resample cell splits as well as perturbations.

    A point estimate averaged over many random splits is less variable than any one split,
    so an interval that conditions on the splits drawn is too narrow. On the real data that
    produced an interval excluding zero for a quantity that cannot be negative.
    """
    rng = np.random.RandomState(5)
    n_split, n_pert = 20, 50
    cover2, width2, width1 = [], [], []
    for _ in range(120):
        mu = rng.normal(0, 1.0 / np.sqrt(K), size=(n_pert, K))
        ssw_stack = np.zeros((n_split, n_pert))
        d2_stack = np.zeros((n_split, n_pert, n_pert))
        for s in range(n_split):
            S = mu[:, None, :] + rng.normal(0, 1.0 / np.sqrt(K), size=(n_pert, H, K))
            ssw_stack[s], d2_stack[s] = stats_from_gram(gram_of(S), n_pert, H)
        lo2, hi2 = np.percentile(
            bootstrap_delta2(ssw_stack, d2_stack, H, n_boot=400, seed=0), [2.5, 97.5])
        lo1, hi1 = np.percentile(
            bootstrap_delta2(ssw_stack.mean(0), d2_stack.mean(0), H, n_boot=400, seed=0),
            [2.5, 97.5])
        cover2.append(lo2 <= 2.0 <= hi2)
        width2.append(hi2 - lo2)
        width1.append(hi1 - lo1)
    assert np.mean(cover2) >= 0.88, f"two-level coverage {np.mean(cover2):.3f}"
    assert np.mean(width2) > np.mean(width1), "two-level interval should not be narrower"


def test_nearest_neighbour_separates_configurations_the_mean_cannot():
    """Guards the finding that drives the law: the mean over pairs is the wrong statistic.

    Two configurations are built with an identical mean squared separation, one spread over
    many directions and one collapsed onto a line. A discrimination score is a
    nearest-neighbour question, so the collapsed configuration is far harder, and a
    statistic that cannot see the difference cannot describe discrimination.
    """
    rng = np.random.RandomState(0)
    n = 20
    configs = {
        "line": np.linspace(-1, 1, n)[:, None] * rng.normal(size=(1, K)),
        "spread": rng.normal(size=(n, K)),
    }
    measured = {}
    for name, mu in configs.items():
        mu = mu - mu.mean(0)
        g = mu @ mu.T
        d2 = np.diag(g)[:, None] + np.diag(g)[None, :] - 2 * g
        mu = mu / np.sqrt(d2[np.triu_indices(n, 1)].mean()) * np.sqrt(2.0)
        S = mu[:, None, :] + rng.normal(0, 1.0 / np.sqrt(K), size=(n, H, K))
        measured[name] = nearest_neighbour_rho2(S, 1.0)["nn_median"]
    assert measured["spread"] > 10 * measured["line"], (
        f"nearest-neighbour statistic failed to separate the configurations: "
        f"{measured}")


def test_geometric_mean_is_undefined_rather_than_floored():
    """Guards against a summary that prints a number where the data support none.

    Debiasing sends individual separations below zero whenever the true separation is under
    the noise, which is the regime these datasets sit in. Flooring them at a small positive
    value makes the geometric mean a function of the floor and of how many values hit it,
    while still looking like a measurement. The median and the above-noise fraction stay
    defined and are what should be read there.
    """
    sep = np.array([-0.5, 0.1, 0.4, 2.0])
    out = summarise_separations(sep, 1.0)
    assert np.isnan(out["nn_geomean"]), "geometric mean should be undefined here"
    assert np.isfinite(out["nn_median"]), "median should stay defined"
    assert out["frac_above_noise"] == 0.75
    assert out["n_nonpositive"] == 1
    positive = summarise_separations(np.array([0.5, 2.0]), 1.0)
    assert np.isfinite(positive["nn_geomean"]), "geometric mean should exist when all > 0"


def test_nearest_neighbour_rejects_a_single_measurement():
    """Cross-fitting needs two measurements; one cannot remove the selection bias."""
    try:
        nearest_neighbour_rho2(np.zeros((5, 1, 3)), 1.0)
    except ValueError:
        return
    raise AssertionError("a single measurement should raise")


def test_single_measurement_is_refused():
    """One measurement per perturbation cannot separate noise from signal, so it must raise."""
    try:
        signal_noise(np.zeros(3), np.zeros((3, 3)), 1)
    except ValueError:
        return
    raise AssertionError("h=1 should raise rather than return a noise-free answer")


def test_mismatched_stacking_is_refused():
    """Stacked statistics for one argument and not the other is a caller error, not a default."""
    try:
        bootstrap_delta2(np.zeros((4, 3)), np.zeros((3, 3)), H, n_boot=2)
    except ValueError:
        return
    raise AssertionError("mismatched stacking should raise")


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
