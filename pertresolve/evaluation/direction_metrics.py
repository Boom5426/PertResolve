"""Direction-alignment metrics for pseudobulk delta prediction."""
import numpy as np


def pearson_delta(pred_delta: np.ndarray, real_delta: np.ndarray) -> float:
    """Pearson correlation between predicted and real pseudobulk delta."""
    if np.std(pred_delta) < 1e-12 or np.std(real_delta) < 1e-12:
        return 0.0
    return float(np.corrcoef(pred_delta, real_delta)[0, 1])


def pearson_delta_top20(
    pred_delta: np.ndarray,
    real_delta: np.ndarray,
    top20_idx: np.ndarray,
) -> float:
    """Pearson-Δ on the top-20 highest-variance genes (SCALE metric)."""
    p = pred_delta[top20_idx]
    r = real_delta[top20_idx]
    if np.std(p) < 1e-12 or np.std(r) < 1e-12:
        return 0.0
    return float(np.corrcoef(p, r)[0, 1])


def delta_cosine(pred_delta: np.ndarray, real_delta: np.ndarray) -> float:
    """Cosine similarity between predicted and real pseudobulk delta."""
    n1, n2 = np.linalg.norm(pred_delta), np.linalg.norm(real_delta)
    if n1 < 1e-12 or n2 < 1e-12:
        return 0.0
    return float(np.dot(pred_delta, real_delta) / (n1 * n2))


def mae(pred_delta: np.ndarray, real_delta: np.ndarray) -> float:
    """Mean absolute error of pseudobulk delta."""
    return float(np.mean(np.abs(pred_delta - real_delta)))
