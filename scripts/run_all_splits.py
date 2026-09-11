#!/usr/bin/env python
"""PertResolve-Eval v4: 10-metric enriched grid
20 methods × 5 splits × 4 genes × 10 metrics = ~10k rows
Metrics:
  Direction: pearson_delta, pearson_delta_top20, delta_cosine
  Ranking:   PDS_cos, PDS_L1, PDS_L2
  DE:        DE_overlap, DE_LFC_spearman, direction_agreement
  Recon:     MAE

Usage:
  python scripts/run_all_splits.py --base /path/to/processed-data --out /path/to/outdir
  PERTRESOLVE_DATA=/path/to/processed-data python scripts/run_all_splits.py --out /path/to/outdir

  --base       compute workspace holding pertresolve_bench.csv, joint_arrays.npz,
               gata1_arrays.npz, jak1_arrays.npz and esm1v_embeddings.npz.
               Falls back to the PERTRESOLVE_DATA environment variable.
  --out        required output directory; results_v4_10metrics.csv is written into it.
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from pertresolve.paths import ENV_VAR, reject_repo_results, require_inputs, resolve_base

# Arguments are resolved before the numeric stack is imported and before any data
# is touched, so the order of every seeded operation below is unchanged.
parser = argparse.ArgumentParser(
    description="PertResolve-Eval v4: 10-metric enriched grid.")
parser.add_argument(
    "--base", default=None,
    help="Compute workspace holding pertresolve_bench.csv, joint_arrays.npz, "
         "gata1_arrays.npz, jak1_arrays.npz and esm1v_embeddings.npz. "
         f"Falls back to the {ENV_VAR} environment variable.")
parser.add_argument(
    "--out", required=True,
    help="Output directory; results_v4_10metrics.csv is written into it. "
         "Required on purpose, so a re-run cannot overwrite the committed "
         "canonical table under the repository's results/ directory.")
args = parser.parse_args()

BASE = resolve_base(args.base)
OUT_DIR = reject_repo_results(args.out)

import numpy as np, pandas as pd, os, sys, warnings, time
from scipy import stats
from sklearn.linear_model import Ridge, Lasso
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.neighbors import KNeighborsRegressor
from sklearn.neural_network import MLPRegressor
from sklearn.decomposition import PCA
warnings.filterwarnings("ignore")

BENCH = os.path.join(BASE, "pertresolve_bench.csv")
require_inputs(Path(BENCH))
NSUB = 300
HOLDOUT_FRAC = 0.35
DE_TOP_K = 50  # top-k DEGs for overlap/precision

def load_gene_data(gene):
    if gene in ("TP53","KRAS"):
        d = np.load(os.path.join(BASE,"joint_arrays.npz"), allow_pickle=True)
        sfx = "tp" if gene=="TP53" else "kr"
        X = d[f"X{sfx}"]; v = d[f"v{sfx}"]; th = d[f"TH{sfx}"]
        return X, v, th
    elif gene == "GATA1":
        d = np.load(os.path.join(BASE,"gata1_arrays.npz"), allow_pickle=True)
        return d["X"], d["cell_variants"], None
    elif gene == "JAK1":
        d = np.load(os.path.join(BASE,"jak1_arrays.npz"), allow_pickle=True)
        return d["X"], d["variant_labels"], None

def pseudobulk_deltas(X, variants, n_sub=NSUB):
    uv = np.unique(variants)
    wt_mask = (variants == "WT") | (variants == "wt") | (variants == "WT_control")
    if wt_mask.sum() == 0:
        wt_mean = np.zeros(X.shape[1])
    else:
        idx = np.where(wt_mask)[0]
        if len(idx) > n_sub: idx = np.random.choice(idx, n_sub, replace=False)
        wt_mean = X[idx].mean(0)
    deltas = {}
    for v in uv:
        if v in ("WT","wt","WT_control"): continue
        mask = variants == v
        idx = np.where(mask)[0]
        if len(idx) < 5: continue
        if len(idx) > n_sub: idx = np.random.choice(idx, n_sub, replace=False)
        deltas[v] = X[idx].mean(0) - wt_mean
    return deltas, wt_mean

def compute_de_genes(X, variants, v, wt_mean_unused, n_sub=NSUB):
    """Per-variant DE test: variant cells vs WT cells, return sorted gene indices by |t-stat|"""
    wt_mask = (variants == "WT") | (variants == "wt") | (variants == "WT_control")
    v_mask = variants == v
    wt_idx = np.where(wt_mask)[0]
    v_idx = np.where(v_mask)[0]
    if len(wt_idx) > n_sub: wt_idx = np.random.choice(wt_idx, n_sub, replace=False)
    if len(v_idx) > n_sub: v_idx = np.random.choice(v_idx, n_sub, replace=False)
    if len(v_idx) < 5 or len(wt_idx) < 5:
        return None, None, None
    Xv = X[v_idx]; Xw = X[wt_idx]
    # t-test per gene
    t_stats, p_vals = stats.ttest_ind(Xv, Xw, axis=0, equal_var=False)
    t_stats = np.nan_to_num(t_stats, 0)
    lfc = Xv.mean(0) - Xw.mean(0)
    return t_stats, lfc, p_vals

def score_variant(pred_delta, real_delta, all_real_deltas, variant_names, target_v,
                  real_t_stats, real_lfc, pred_lfc_approx):
    """Compute all 10 metrics for one held-out variant."""
    nd = len(pred_delta)
    out = {}
    
    # --- Direction metrics ---
    # 1. pearson_delta
    if np.std(pred_delta) > 1e-12 and np.std(real_delta) > 1e-12:
        out['pearson_delta'] = np.corrcoef(pred_delta, real_delta)[0,1]
    else:
        out['pearson_delta'] = 0.0
    
    # 2. pearson_delta_top20 (top-20 highest-variance genes across training deltas)
    # We pass in top20_idx externally
    
    # 3. delta_cosine
    n1 = np.linalg.norm(pred_delta); n2 = np.linalg.norm(real_delta)
    out['delta_cosine'] = np.dot(pred_delta, real_delta) / (n1*n2+1e-12)
    
    # --- Ranking metrics (PDS with 3 distances) ---
    for dist_name, dist_fn in [
        ('PDS_cos', lambda a,b: 1 - np.dot(a,b)/(np.linalg.norm(a)*np.linalg.norm(b)+1e-12)),
        ('PDS_L1',  lambda a,b: np.sum(np.abs(a-b))),
        ('PDS_L2',  lambda a,b: np.sqrt(np.sum((a-b)**2))),
    ]:
        if np.linalg.norm(pred_delta) < 1e-12:
            out[dist_name] = 0.5  # tie-aware
            continue
        sims = []
        for vn in variant_names:
            rd = all_real_deltas[vn]
            sims.append((vn, dist_fn(pred_delta, rd)))
        sims.sort(key=lambda x: x[1])
        rank = next((i for i,(vn,_) in enumerate(sims) if vn==target_v), len(sims))
        # tie-aware: count how many share the same distance
        target_dist = sims[rank][1] if rank < len(sims) else float('inf')
        tied = [i for i,(vn,d) in enumerate(sims) if abs(d - target_dist) < 1e-12]
        avg_rank = np.mean(tied)
        out[dist_name] = 1.0 - avg_rank / (len(sims) - 1) if len(sims) > 1 else 0.5
    
    # --- DE metrics ---
    if real_t_stats is not None:
        abs_t = np.abs(real_t_stats)
        real_de_idx = set(np.argsort(abs_t)[-DE_TOP_K:])
        
        # pred DE: genes with largest |pred_delta|
        abs_pred = np.abs(pred_delta)
        pred_de_idx = set(np.argsort(abs_pred)[-DE_TOP_K:])
        
        # 7. DE_overlap
        overlap = len(real_de_idx & pred_de_idx)
        out['DE_overlap'] = overlap / DE_TOP_K
        
        # 8. DE_LFC_spearman: on real DE genes, correlate real LFC with pred delta
        real_de_list = list(real_de_idx)
        if len(real_de_list) >= 5 and real_lfc is not None:
            r_lfc = real_lfc[real_de_list]
            p_lfc = pred_delta[real_de_list]
            if np.std(r_lfc) > 1e-12 and np.std(p_lfc) > 1e-12:
                out['DE_LFC_spearman'] = stats.spearmanr(r_lfc, p_lfc)[0]
            else:
                out['DE_LFC_spearman'] = 0.0
        else:
            out['DE_LFC_spearman'] = 0.0
        
        # 9. direction_agreement: fraction of DE genes with correct sign
        real_de_list = list(real_de_idx)
        if real_lfc is not None:
            signs_real = np.sign(real_lfc[real_de_list])
            signs_pred = np.sign(pred_delta[real_de_list])
            out['direction_agreement'] = np.mean(signs_real == signs_pred)
        else:
            out['direction_agreement'] = 0.5
    else:
        out['DE_overlap'] = np.nan
        out['DE_LFC_spearman'] = np.nan
        out['direction_agreement'] = np.nan
    
    # --- Recon metric ---
    # 10. MAE
    out['MAE'] = np.mean(np.abs(pred_delta - real_delta))
    
    return out

def run_all():
    bench = pd.read_csv(BENCH)
    genes = ['TP53','KRAS','GATA1','JAK1']
    splits = ['split1','split2','split3','split5','split6']
    
    # heads
    heads = {
        'Ridge-theta': (Ridge(alpha=1.0), 'theta'),
        'Ridge-esm': (Ridge(alpha=1.0), 'esm'),
        'Ridge-esm+theta': (Ridge(alpha=1.0), 'esm+theta'),
        'Lasso-theta': (Lasso(alpha=0.01, max_iter=2000), 'theta'),
        'Lasso-esm': (Lasso(alpha=0.01, max_iter=2000), 'esm'),
        'Lasso-esm+theta': (Lasso(alpha=0.01, max_iter=2000), 'esm+theta'),
        'RF-theta': (RandomForestRegressor(n_estimators=50, max_depth=6, n_jobs=-1, random_state=0), 'theta'),
        'RF-esm': (RandomForestRegressor(n_estimators=50, max_depth=6, n_jobs=-1, random_state=0), 'esm'),
        'RF-esm+theta': (RandomForestRegressor(n_estimators=50, max_depth=6, n_jobs=-1, random_state=0), 'esm+theta'),
        'GBoost-theta': ('gboost', 'theta'),
        'GBoost-esm': ('gboost', 'esm'),
        'GBoost-esm+theta': ('gboost', 'esm+theta'),
        'KNN-theta': (KNeighborsRegressor(n_neighbors=3), 'theta'),
        'KNN-esm': (KNeighborsRegressor(n_neighbors=3), 'esm'),
        'KNN-esm+theta': (KNeighborsRegressor(n_neighbors=3), 'esm+theta'),
        'MLP-theta': (MLPRegressor(hidden_layer_sizes=(128,64), max_iter=500, random_state=0), 'theta'),
        'MLP-esm': (MLPRegressor(hidden_layer_sizes=(128,64), max_iter=500, random_state=0), 'esm'),
        'MLP-esm+theta': (MLPRegressor(hidden_layer_sizes=(128,64), max_iter=500, random_state=0), 'esm+theta'),
        'Gene-mean': ('genemean', 'theta'),
        'WT-null': ('wtnull', 'theta'),
    }
    
    # ESM embeddings
    esm_path = os.path.join(BASE, "esm1v_embeddings.npz")
    if not os.path.exists(esm_path):
        raise FileNotFoundError(
            f"Required ESM-1v embeddings not found: {esm_path}. "
            "This analysis requires the full representation grid and will not "
            "silently fall back to a reduced set of heads."
        )
    esm_data = np.load(esm_path, allow_pickle=True)
    
    theta_cols = ['d_hydro','d_vol','d_charge','fold_core','cat_switch','is_hotspot']
    # The first three are physicochemical differences in physical units; the last
    # three are 0/1 position annotations that are already on a common scale.
    n_physchem = 3
    
    all_rows = []
    for sp in splits:
        sp_col = f"{sp}_role"
        t0 = time.time()
        for gene in genes:
            gb = bench[bench.gene == gene]
            # Every iteration over these sets is sorted below. Python set order over
            # strings is hash-randomized per process, and the DE cache and the
            # pseudobulk subsampler both draw from one global np.random stream, so
            # unsorted iteration made the whole grid irreproducible between runs.
            train_vars = set(gb[gb[sp_col]=='train']['variant'])
            test_vars = set(gb[gb[sp_col]=='test']['variant'])
            if len(test_vars) == 0: continue
            
            X, variants, theta_raw = load_gene_data(gene)
            deltas, wt_mean = pseudobulk_deltas(X, variants)
            
            # top-20 high-variance genes (across train deltas)
            train_delta_mat = np.array([deltas[v] for v in sorted(train_vars) if v in deltas])
            if len(train_delta_mat) > 0:
                gene_var = np.var(train_delta_mat, axis=0)
                top20_idx = np.argsort(gene_var)[-20:]
            else:
                top20_idx = np.arange(min(20, X.shape[1]))
            
            # Build theta features
            theta_map = {}
            for _, row in gb.iterrows():
                theta_map[row['variant']] = np.array([row[c] for c in theta_cols],
                                                     dtype=float)

            # Standardize the physicochemical terms on THIS SPLIT'S TRAINING VARIANTS.
            # configs/bench_features.yaml emits them in physical units (Kyte-Doolittle,
            # cubic angstroms, charge at pH 7) and states that standardization belongs
            # to the fitting step, not to the released table: a scaler fitted over the
            # whole table lets held-out variants set the scale of the training features.
            # Without this, d_vol in cubic angstroms carries >99% of theta's variance and
            # the scale-sensitive heads see essentially one feature.
            fit_vars = [v for v in sorted(train_vars) if v in deltas and v in theta_map]
            if len(fit_vars) >= 2:
                fit_mat = np.array([theta_map[v][:n_physchem] for v in fit_vars])
                mu = fit_mat.mean(axis=0)
                sd = fit_mat.std(axis=0, ddof=1)
                sd = np.where(np.isfinite(sd) & (sd > 0), sd, 1.0)
                for vec in theta_map.values():
                    vec[:n_physchem] = (vec[:n_physchem] - mu) / sd
            
            # ESM features
            esm_map = {}
            if esm_data is not None:
                for v in sorted(train_vars | test_vars):
                    key = f"{gene}__{v}"
                    if key in esm_data:
                        esm_map[v] = esm_data[key]
            
            # DE cache per test variant
            de_cache = {}
            for v in sorted(test_vars):
                if v in deltas:
                    t_stats, lfc, pvals = compute_de_genes(X, variants, v, wt_mean)
                    de_cache[v] = (t_stats, lfc)
            
            heldout_list = [v for v in sorted(test_vars) if v in deltas]
            train_list = [v for v in sorted(train_vars) if v in deltas]
            all_eval_vars = [v for v in (train_list + heldout_list)]
            
            for method_name, (model_spec, feat_type) in heads.items():
                # Build features
                feat_train, feat_test = [], []
                y_train = []
                train_ok, test_ok = [], []
                for v in train_list:
                    f = get_feat(v, feat_type, theta_map, esm_map)
                    if f is not None:
                        feat_train.append(f); y_train.append(deltas[v]); train_ok.append(v)
                for v in heldout_list:
                    f = get_feat(v, feat_type, theta_map, esm_map)
                    if f is not None:
                        feat_test.append(f); test_ok.append(v)
                
                if len(feat_train) < 3 or len(feat_test) == 0: continue
                Xtr = np.array(feat_train); Ytr = np.array(y_train)
                Xte = np.array(feat_test)
                
                n_genes_out = Ytr.shape[1]
                
                # Fit & predict
                if model_spec == 'wtnull':
                    preds = {v: np.zeros(n_genes_out) for v in test_ok}
                elif model_spec == 'genemean':
                    mean_d = Ytr.mean(0)
                    preds = {v: mean_d for v in test_ok}
                elif model_spec == 'gboost':
                    # PCA for output
                    n_comp = min(50, n_genes_out, len(Xtr)-1)
                    pca_out = PCA(n_components=max(n_comp,1), random_state=0)
                    Ytr_pca = pca_out.fit_transform(Ytr)
                    from sklearn.multioutput import MultiOutputRegressor
                    gb_model = MultiOutputRegressor(
                        GradientBoostingRegressor(n_estimators=50, max_depth=3, random_state=0),
                        n_jobs=-1)
                    gb_model.fit(Xtr, Ytr_pca)
                    Yte_pca = gb_model.predict(Xte)
                    Yte_full = pca_out.inverse_transform(Yte_pca)
                    preds = {v: Yte_full[i] for i, v in enumerate(test_ok)}
                else:
                    from sklearn.base import clone
                    m = clone(model_spec)
                    if isinstance(m, RandomForestRegressor) and n_genes_out > 200:
                        n_comp = min(50, n_genes_out, len(Xtr)-1)
                        pca_out = PCA(n_components=max(n_comp,1), random_state=0)
                        Ytr_pca = pca_out.fit_transform(Ytr)
                        m.fit(Xtr, Ytr_pca)
                        Yte_pca = m.predict(Xte)
                        Yte_full = pca_out.inverse_transform(Yte_pca)
                    else:
                        m.fit(Xtr, Ytr)
                        Yte_full = m.predict(Xte)
                    preds = {v: Yte_full[i] for i, v in enumerate(test_ok)}
                
                # Score each held-out variant
                for v in test_ok:
                    pred_d = preds[v]
                    real_d = deltas[v]
                    t_stats, lfc = de_cache.get(v, (None, None))
                    
                    scores = score_variant(pred_d, real_d, deltas, all_eval_vars, v,
                                          t_stats, lfc, pred_d)
                    
                    # pearson_delta_top20
                    if np.std(pred_d[top20_idx]) > 1e-12 and np.std(real_d[top20_idx]) > 1e-12:
                        scores['pearson_delta_top20'] = np.corrcoef(pred_d[top20_idx], real_d[top20_idx])[0,1]
                    else:
                        scores['pearson_delta_top20'] = 0.0
                    
                    scores['method'] = method_name
                    scores['gene'] = gene
                    scores['split'] = sp
                    scores['variant'] = v
                    all_rows.append(scores)
        
        elapsed = time.time() - t0
        n = len([r for r in all_rows if r['split']==sp])
        print(f" {sp} done, {n} rows, {elapsed:.0f}s", flush=True)
    
    df = pd.DataFrame(all_rows)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out_path = OUT_DIR / "results_v4_10metrics.csv"
    df.to_csv(out_path, index=False)
    print(f"\nALL_DONE: {len(df)} rows -> {out_path}", flush=True)
    
    # Summary
    metrics = ['PDS_cos','PDS_L1','PDS_L2','pearson_delta','pearson_delta_top20',
               'delta_cosine','DE_overlap','DE_LFC_spearman','direction_agreement','MAE']
    for m in metrics:
        pivot = df.groupby(['method','split','gene'])[m].mean().reset_index()
        overall = pivot.groupby(['method','split'])[m].mean().reset_index()
        grand = overall.groupby('method')[m].mean().sort_values(ascending=(m=='MAE'))
        print(f"\n=== {m} (mean-of-gene-means, overall) ===")
        for method, val in grand.items():
            print(f" {method:20s} {val:.3f}")

def get_feat(v, feat_type, theta_map, esm_map):
    if feat_type == 'theta':
        return theta_map.get(v)
    elif feat_type == 'esm':
        return esm_map.get(v)
    elif feat_type == 'esm+theta':
        t = theta_map.get(v)
        e = esm_map.get(v)
        if t is not None and e is not None:
            return np.concatenate([e, t])
        return None
    return None

if __name__ == "__main__":
    np.random.seed(42)
    run_all()
