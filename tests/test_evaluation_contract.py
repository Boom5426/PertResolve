"""Regression tests for the prediction-evaluation and DE effect-size wording."""

import numpy as np

from pertresolve.metrics import evaluate_variant


def test_evaluate_variant_reports_pds_and_pearson_for_full_same_gene_pool():
    real = {
        "R175H": np.array([1.0, -2.0, 0.5, 0.0]),
        "R273C": np.array([0.8, -1.5, 0.4, 0.1]),
        "R248W": np.array([-0.4, 0.3, 0.2, 0.2]),
    }
    result = evaluate_variant(
        pred_delta=real["R175H"],
        real_delta=real["R175H"],
        target_variant="R175H",
        all_real_deltas=real,
        candidate_variants=list(real),
    )

    assert result["PDS_cos"] == 1.0
    assert result["pearson_delta"] == 1.0
