"""Contracts for the replicate-reference study (docs/PREREG_REPLICATE_REFERENCE_v1.md).

These run from committed code alone: no workspace, no atlas, no GEO download. They
pin the three things that would silently invalidate the study if they drifted, and
one that already did go wrong once and was caught by an anchor check.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]


def _load(rel: str, name: str):
    spec = importlib.util.spec_from_file_location(name, REPO / rel)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


rr = _load("scripts/analysis/replicate_reference.py", "rr")
tp = _load("scripts/dataset_screen/replicate_reference_gse311877.py", "tp63rr")

from pertresolve.resolution.scaling import tie_aware_pds


@pytest.mark.parametrize("noise", [0.2, 1.0, 5.0])
def test_both_scorers_equal_the_packaged_pds(noise):
    """The study must score through the manuscript's harness, not a lookalike.

    Both scripts compute PDS per row so that the paired statistic has a per-unit
    atom; the package exposes only the mean. The mean of the rows must equal it
    exactly, or the study is not measuring the quantity the paper reports.
    """
    rng = np.random.default_rng(0)
    truth = rng.normal(size=(18, 400))
    pred = truth + rng.normal(size=(18, 400)) * noise
    expected = tie_aware_pds(pred, truth)
    assert rr.per_pert_pds(pred, truth).mean() == pytest.approx(expected, abs=1e-12)
    assert tp.per_allele_pds(pred, truth).mean() == pytest.approx(expected, abs=1e-12)


def test_a_degenerate_panel_scores_chance_not_a_win():
    """Identical profiles are all tied, so mid-ranking must return exactly 0.5."""
    x = np.tile(np.arange(50.0), (12, 1))
    assert rr.per_pert_pds(x, x).mean() == pytest.approx(0.5, abs=1e-12)
    assert tp.per_allele_pds(x, x).mean() == pytest.approx(0.5, abs=1e-12)


def test_batch_groups_partition_without_overlap():
    for levels in (["1", "2", "3"], list("abcdefgh"), ["rep1", "rep2"]):
        a, b = rr.batch_groups(levels)
        assert set(a) | set(b) == set(levels)
        assert not set(a) & set(b)
        assert abs(len(a) - len(b)) <= 1


def test_molecule_halves_are_exact_and_disjoint():
    """Depth matching is structural, not approximate.

    Binomial thinning at p = N/total gives N only in expectation, which breaks the
    matching the whole TP63 comparison rests on.
    """
    rng = np.random.default_rng(3)
    counts = rng.poisson(40, size=2000).astype(np.int64)
    n = 5_000
    h1, h2 = tp.draw_halves(rng, counts, n)
    assert h1.sum() == n and h2.sum() == n
    assert np.all(h1 + h2 <= counts)


def test_the_shifted_log_transform_cancels_depth_in_a_delta():
    """log2(CPM+1) is not comparable along a depth ladder; the primary one is.

    One molecule is 500 CPM at N = 2,000 and 2.5 CPM at N = 400,000, so under
    log2(CPM+1) the pseudocount means something different at every rung. The primary
    transform's depth term is gene-independent and must vanish from a difference.
    """
    y1 = np.array([0.0, 3.0, 40.0])
    y2 = np.array([1.0, 5.0, 30.0])
    d_small = tp.transform(y1, 2_000, "shifted_log_count") - tp.transform(y2, 2_000, "shifted_log_count")
    d_large = tp.transform(y1, 400_000, "shifted_log_count") - tp.transform(y2, 400_000, "shifted_log_count")
    assert np.allclose(d_small, d_large)
    naive_small = tp.transform(y1, 2_000, "log2_cpm1") - tp.transform(y2, 2_000, "log2_cpm1")
    naive_large = tp.transform(y1, 400_000, "log2_cpm1") - tp.transform(y2, 400_000, "log2_cpm1")
    assert not np.allclose(naive_small, naive_large)


def test_the_legacy_replicate_reference_rule_matches_the_frozen_decision_table():
    """Legacy replicate-reference decision cases, kept separate from resolution_report."""
    tau = 0.05
    at_chance = dict(pds_P=0.50, pds_P_lo=0.49, pds_P_hi=0.51, pds_X=0.50,
                     delta_PX_lo=-0.01, delta_PX_hi=0.01)
    assert rr.verdict(at_chance, tau) == "NOT_INFORMATIVE"

    tight = dict(pds_P=0.80, pds_P_lo=0.78, pds_P_hi=0.82, pds_X=0.79,
                 delta_PX_lo=-0.01, delta_PX_hi=0.02)
    assert rr.verdict(tight, tau) == "CONCORDANT"

    big = dict(pds_P=0.80, pds_P_lo=0.78, pds_P_hi=0.82, pds_X=0.70,
               delta_PX_lo=0.06, delta_PX_hi=0.14)
    assert rr.verdict(big, tau) == "OVERSTATING"

    flipped = dict(pds_P=0.70, pds_P_lo=0.68, pds_P_hi=0.72, pds_X=0.62,
                   delta_PX_lo=0.02, delta_PX_hi=0.14)
    assert rr.verdict(flipped, tau) == "OVERSTATING"

    negative = dict(pds_P=0.80, pds_P_lo=0.78, pds_P_hi=0.82, pds_X=0.90,
                    delta_PX_lo=-0.14, delta_PX_hi=-0.06)
    assert rr.verdict(negative, tau) == "UNDERSTATING"

    # covers both zero and the margin: absence of evidence is never concordance
    unclear = dict(pds_P=0.80, pds_P_lo=0.78, pds_P_hi=0.82, pds_X=0.78,
                   delta_PX_lo=-0.01, delta_PX_hi=0.09)
    assert rr.verdict(unclear, tau) == "INCONCLUSIVE"


def test_the_exact_sibling_test_enumerates_and_reports_its_floor():
    """Five pairs give 32 sign patterns, so the two-sided p can never beat 2/32."""
    x = np.array([0.3, 0.4, 0.2, 0.5, 0.35])
    p = tp.exact_signflip_p(x)
    assert p == pytest.approx(2 / 32)
    assert tp.exact_signflip_p(np.zeros(5)) == pytest.approx(1.0)
