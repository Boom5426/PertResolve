"""Top-level prediction-evaluation convenience functions.

``evaluate_variant`` is the paper-facing single-query adapter. It combines the same
held-out prediction with direction, full-pool same-gene ranking and optional
differential-expression fidelity metrics; it does not define a new model or a new
benchmark protocol.
"""
import numpy as np
from .evaluation.pds import pds_score
from .evaluation.direction_metrics import pearson_delta, pearson_delta_top20, delta_cosine, mae
from .evaluation.de_metrics import de_overlap, de_lfc_spearman, direction_agreement


def evaluate_variant(
    pred_delta: np.ndarray,
    real_delta: np.ndarray,
    target_variant: str,
    all_real_deltas: dict,
    candidate_variants: list,
    top20_idx: np.ndarray = None,
    real_t_stats: np.ndarray = None,
    real_lfc: np.ndarray = None,
    de_top_k: int = 50,
) -> dict:
    """Compute the PertResolve-Eval metrics for one held-out variant.

    Parameters
    ----------
    pred_delta : array (n_genes,)
    real_delta : array (n_genes,)
    target_variant : str
    all_real_deltas : dict {variant: array}
    candidate_variants : list of str
    top20_idx : array of int, optional
        Indices of top-20 highest-variance genes across training deltas.
    real_t_stats : array, optional
        T-statistics from the variant-vs-reference DE test. Their absolute values select
        the real top-``de_top_k`` genes; they are not themselves an effect-size target.
    real_lfc : array, optional
        Signed DE effect sizes on the same gene scale as ``real_delta`` (typically the
        variant-minus-reference mean difference). ``DE_LFC_spearman`` evaluates rank
        agreement on the genes selected by ``real_t_stats``; it is not a calibration or
        magnitude-error metric.
    de_top_k : int
        Top-k DEGs for overlap metrics.

    Returns
    -------
    dict with keys: PDS_cos, PDS_L1, PDS_L2, pearson_delta,
        pearson_delta_top20, delta_cosine, DE_overlap,
        DE_LFC_spearman, direction_agreement, MAE.

    Notes
    -----
    ``candidate_variants`` is the full same-gene candidate pool for a paper-style PDS
    query and must include ``target_variant``. The query can be held out from model
    fitting while its measured response is still used by the evaluator as the target.
    """
    out = {}

    # Ranking metrics
    for dist in ("cos", "L1", "L2"):
        out[f"PDS_{dist}"] = pds_score(
            pred_delta, target_variant, all_real_deltas, candidate_variants, dist
        )

    # Direction metrics
    out["pearson_delta"] = pearson_delta(pred_delta, real_delta)
    out["delta_cosine"] = delta_cosine(pred_delta, real_delta)
    out["MAE"] = mae(pred_delta, real_delta)

    if top20_idx is not None:
        out["pearson_delta_top20"] = pearson_delta_top20(pred_delta, real_delta, top20_idx)
    else:
        out["pearson_delta_top20"] = out["pearson_delta"]

    # DE metrics
    if real_t_stats is not None and real_lfc is not None:
        out["DE_overlap"] = de_overlap(pred_delta, real_t_stats, de_top_k)
        out["DE_LFC_spearman"] = de_lfc_spearman(pred_delta, real_lfc, real_t_stats, de_top_k)
        out["direction_agreement"] = direction_agreement(pred_delta, real_lfc, real_t_stats, de_top_k)
    else:
        out["DE_overlap"] = np.nan
        out["DE_LFC_spearman"] = np.nan
        out["direction_agreement"] = np.nan

    return out
