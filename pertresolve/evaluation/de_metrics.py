"""Differential expression fidelity metrics.

These metrics evaluate whether the predicted perturbation response recovers
the correct set of differentially expressed genes (DEGs), their effect sizes,
and their directions.
"""
import numpy as np
# scipy is imported where it is used rather than at module import, so that the
# numpy-only resolution diagnostics can be installed and used without it. Calling a
# function that needs it raises there, which is clearer than a package that will not
# import at all.


def compute_de_genes(X, variant_mask, wt_mask, n_sub=300):
    """Run per-gene t-test between variant cells and WT cells.

    Returns (t_stats, lfc, p_vals) arrays of shape (n_genes,).
    """
    v_idx = np.where(variant_mask)[0]
    w_idx = np.where(wt_mask)[0]
    if len(v_idx) > n_sub:
        v_idx = np.random.choice(v_idx, n_sub, replace=False)
    if len(w_idx) > n_sub:
        w_idx = np.random.choice(w_idx, n_sub, replace=False)
    if len(v_idx) < 5 or len(w_idx) < 5:
        return None, None, None
    Xv, Xw = X[v_idx], X[w_idx]
    t_stats, p_vals = stats.ttest_ind(Xv, Xw, axis=0, equal_var=False)
    t_stats = np.nan_to_num(t_stats, 0.0)
    lfc = Xv.mean(0) - Xw.mean(0)
    return t_stats, lfc, p_vals


def de_overlap(
    pred_delta: np.ndarray,
    real_t_stats: np.ndarray,
    top_k: int = 50,
) -> float:
    """Fraction of top-k real DEGs recovered in predicted top-k.

    Parameters
    ----------
    pred_delta : array (n_genes,)
        Predicted pseudobulk shift.
    real_t_stats : array (n_genes,)
        T-statistics from variant-vs-WT DE test.
    top_k : int
        Number of top DEGs to compare.
    """
    real_de = set(np.argsort(np.abs(real_t_stats))[-top_k:])
    pred_de = set(np.argsort(np.abs(pred_delta))[-top_k:])
    return len(real_de & pred_de) / top_k


def de_lfc_spearman(
    pred_delta: np.ndarray,
    real_lfc: np.ndarray,
    real_t_stats: np.ndarray,
    top_k: int = 50,
) -> float:
    """Spearman correlation of LFC on top-k real DEGs."""
    de_idx = list(np.argsort(np.abs(real_t_stats))[-top_k:])
    r_lfc = real_lfc[de_idx]
    p_lfc = pred_delta[de_idx]
    if np.std(r_lfc) < 1e-12 or np.std(p_lfc) < 1e-12:
        return 0.0
    from scipy import stats
    return float(stats.spearmanr(r_lfc, p_lfc)[0])


def direction_agreement(
    pred_delta: np.ndarray,
    real_lfc: np.ndarray,
    real_t_stats: np.ndarray,
    top_k: int = 50,
) -> float:
    """Fraction of top-k real DEGs where predicted sign matches."""
    de_idx = list(np.argsort(np.abs(real_t_stats))[-top_k:])
    signs_real = np.sign(real_lfc[de_idx])
    signs_pred = np.sign(pred_delta[de_idx])
    return float(np.mean(signs_real == signs_pred))
