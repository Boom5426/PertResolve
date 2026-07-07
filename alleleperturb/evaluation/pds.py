"""Perturbation Discrimination Score (PDS) with tie-aware scoring.

PDS measures whether a model can correctly rank a held-out variant's predicted
response against measured responses of other variants. PDS = 1.0 means the
model places the correct variant at rank 1; PDS = 0.5 is chance.

Three distance functions are supported: cosine, L1, L2.
"""
import numpy as np


def _cosine_dist(a, b):
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    if na < 1e-12 or nb < 1e-12:
        return 1.0
    return 1.0 - np.dot(a, b) / (na * nb)


def _l1_dist(a, b):
    return np.sum(np.abs(a - b))


def _l2_dist(a, b):
    return np.sqrt(np.sum((a - b) ** 2))


_DIST_FNS = {"cos": _cosine_dist, "L1": _l1_dist, "L2": _l2_dist}


def pds_score(
    pred_delta: np.ndarray,
    target_variant: str,
    all_real_deltas: dict,
    candidate_variants: list,
    distance: str = "cos",
) -> float:
    """Compute PDS for one held-out variant.

    Parameters
    ----------
    pred_delta : array of shape (n_genes,)
        Predicted pseudobulk shift for the target variant.
    target_variant : str
        Name of the variant being scored.
    all_real_deltas : dict
        {variant_name: np.ndarray} of measured deltas for all candidates.
    candidate_variants : list of str
        Variant names to rank against (includes target).
    distance : str
        One of "cos", "L1", "L2".

    Returns
    -------
    float
        PDS in [0, 1]. Tie-aware: zero-vector prediction → 0.5.
    """
    dist_fn = _DIST_FNS[distance]

    if np.linalg.norm(pred_delta) < 1e-12:
        return 0.5  # zero prediction = chance

    sims = []
    for vn in candidate_variants:
        rd = all_real_deltas[vn]
        sims.append((vn, dist_fn(pred_delta, rd)))
    sims.sort(key=lambda x: x[1])

    rank = next(
        (i for i, (vn, _) in enumerate(sims) if vn == target_variant),
        len(sims),
    )
    # Tie-aware: average rank of all tied entries
    target_dist = sims[rank][1] if rank < len(sims) else float("inf")
    tied = [i for i, (_, d) in enumerate(sims) if abs(d - target_dist) < 1e-12]
    avg_rank = np.mean(tied)

    n = len(sims)
    return 1.0 - avg_rank / (n - 1) if n > 1 else 0.5
