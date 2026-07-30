#!/usr/bin/env python3
"""Split-half power analysis: is strict E-distance evaluation measurement-limited or model-limited?
For each gene × variant × n_subsample × seed, computes:
  D_self  = E(real_halfA, real_halfB)     -- empirical floor
  D_null  = E(WT_sample, real_halfB)      -- WT-null baseline
  D_model = E(predicted, real_halfB)      -- model counterfactual
Then classifies each gene as underpowered / detectable-not-modeled / model-improves / near-floor.

Usage:
    python scripts/figures/run_split_half_power_analysis.py \
        --base /path/to/VCCompass --out /path/to/power_output_dir

``--base`` names the external compute workspace that holds g1_real_cells.npz,
g1_cf_cells_flagship.npz, allele_perturb_bench.csv, joint_arrays.npz,
gata1_arrays.npz, jak1_arrays.npz, model_cfm_v3.pt and theta_v2.csv; it may be
omitted when the VCCOMPASS_BASE environment variable is set. Neither the
workspace nor the checkpoint is redistributed with this repository.

``--out`` is required and receives all four outputs (split_half_power_table.csv,
split_half_power_summary_by_gene.csv, split_half_power_curve.csv and
decision_after_split_half_floor.md). Note the historical layout mismatch: this
script originally wrote them into a ``results/power/`` subdirectory of the
workspace, whereas the committed copy of split_half_power_curve.csv lives
directly under this repository's ``results/``. The two locations are deliberately
left unreconciled here; point ``--out`` at a scratch directory so a re-run cannot
overwrite the committed canonical tables.
"""
import os, time, json, sys
import argparse
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from alleleperturb.paths import resolve_base, require_inputs

parser = argparse.ArgumentParser(
    description="Split-half power analysis of strict E-distance evaluation.")
parser.add_argument(
    "--base", default=None,
    help="VCCompass compute workspace holding the cell arrays, the benchmark "
         "table, theta_v2.csv and model_cfm_v3.pt. Defaults to $VCCOMPASS_BASE.")
parser.add_argument(
    "--out", required=True,
    help="Directory that receives the power table, the per-gene summary, the "
         "power curve and the decision document. Required, and must not be the "
         "repository's results/ directory.")
args = parser.parse_args()

BASE = resolve_base(args.base)
OUT_DIR = Path(args.out)

REAL_CELLS_NPZ = BASE / 'g1_real_cells.npz'
CF_CELLS_NPZ = BASE / 'g1_cf_cells_flagship.npz'
BENCH_CSV = BASE / 'allele_perturb_bench.csv'
JOINT_NPZ = BASE / 'joint_arrays.npz'
GATA1_NPZ = BASE / 'gata1_arrays.npz'
JAK1_NPZ = BASE / 'jak1_arrays.npz'
CFM_CKPT = BASE / 'model_cfm_v3.pt'
THETA_CSV = BASE / 'theta_v2.csv'
require_inputs(REAL_CELLS_NPZ, CF_CELLS_NPZ, BENCH_CSV, JOINT_NPZ,
               GATA1_NPZ, JAK1_NPZ, CFM_CKPT, THETA_CSV)

import numpy as np, pandas as pd
from scipy.spatial.distance import cdist
from sklearn.decomposition import PCA
np.random.seed(42)

# Load all data sources
print("Loading data...", flush=True)
real = np.load(REAL_CELLS_NPZ, allow_pickle=True)   # 474 keys (variants+WT, capped 300)
cf_vae = np.load(CF_CELLS_NPZ, allow_pickle=True)  # 472 keys (VAE counterfactual)
bench = pd.read_csv(BENCH_CSV)

# Also load raw arrays for higher cell counts
joint = np.load(JOINT_NPZ, allow_pickle=True)
gata1_raw = np.load(GATA1_NPZ, allow_pickle=True)
jak1_raw = np.load(JAK1_NPZ, allow_pickle=True)

# Build raw cell pools (not capped at 300)
RAW = {}
for gene, arr_X, arr_lab in [
    ('TP53', joint['Xtp'], joint['vtp']),
    ('KRAS', joint['Xkr'], joint['vkr']),
    ('GATA1', gata1_raw['X'], gata1_raw['cell_variants']),
    ('JAK1', jak1_raw['X'], jak1_raw['variant_labels'])]:
    for v in np.unique(arr_lab):
        m = arr_lab == v
        RAW[f'{gene}__{v}'] = arr_X[m].astype(np.float32)
print(f"RAW pools: {len(RAW)} keys", flush=True)
for g in ['TP53','KRAS','GATA1','JAK1']:
    wk = f'{g}__WT'
    print(f"  {g} WT: {RAW[wk].shape[0]} cells, example variant: {RAW.get(f'{g}__A159P', RAW.get(f'{g}__G12D', RAW.get(f'{g}__C110R', RAW.get(f'{g}__T558I', np.zeros((0,1)))))).shape}", flush=True)

GENES = ['TP53','KRAS','GATA1','JAK1']
GDIMS = {'TP53':939,'KRAS':939,'GATA1':2477,'JAK1':2000}
N_SUBS = [50, 100, 150, 300]
SEEDS = list(range(50))

# --- CFM v3 predictions: need to generate fresh for each variant (use saved model) ---
import torch, torch.nn as nn
dev = 'cuda' if torch.cuda.is_available() else 'cpu'

class VNet(nn.Module):
    def __init__(self, gdims, nth=14, hid=512, tdim=32):
        super().__init__(); self.tdim = tdim
        self.gin = nn.ModuleDict({g: nn.Linear(d, hid) for g, d in gdims.items()})
        self.gout = nn.ModuleDict({g: nn.Linear(hid, d) for g, d in gdims.items()})
        self.thnet = nn.Sequential(nn.Linear(nth, 64), nn.SiLU(), nn.Linear(64, hid))
        self.tnet = nn.Sequential(nn.Linear(tdim, hid), nn.SiLU(), nn.Linear(hid, hid))
        self.trunk = nn.Sequential(nn.Linear(hid, hid), nn.SiLU(), nn.Linear(hid, hid), nn.SiLU())
    def temb(self, t):
        fr = torch.exp(torch.linspace(0, np.log(1000), self.tdim//2, device=t.device))
        a = t[:, None] * fr[None, :]; return torch.cat([torch.sin(a), torch.cos(a)], 1)
    def forward(self, x, t, th, g):
        h = self.gin[g](x) + self.thnet(th) + self.tnet(self.temb(t))
        h = h + self.trunk(h); return self.gout[g](h)

model_cfm = VNet(GDIMS, nth=14).to(dev)
model_cfm.load_state_dict(torch.load(CFM_CKPT, map_location=dev, weights_only=True))
model_cfm.eval()

tv2 = pd.read_csv(THETA_CSV)
tcols = [c for c in tv2.columns if c not in ['gene','variant']]
def theta_for(gene, v):
    r = tv2[(tv2.gene==gene)&(tv2.variant==v)]
    if len(r)==0: return None
    return np.array([r.iloc[0][c] for c in tcols], dtype=np.float32)

@torch.no_grad()
def sample_cfm(gene, th_vec, n=300, steps=40):
    wt = RAW[f'{gene}__WT'].astype(np.float32)
    wi = np.random.randint(0, len(wt), n)
    x = torch.from_numpy(wt[wi]).to(dev)
    th = torch.from_numpy(np.tile(th_vec, (n, 1))).to(dev)
    wt_std = torch.from_numpy(wt.std(0).astype(np.float32)).to(dev)
    sig_base = 0.5 if len(wt) < 1000 else 0.3
    for s in range(steps):
        tt = s / steps; t = torch.full((n,), tt, device=dev)
        drift = model_cfm(x, t, th, gene) / steps
        noise = sig_base * np.sqrt(max(tt*(1-tt), 0)) * np.sqrt(1.0/steps) * wt_std[None,:] * torch.randn_like(x)
        x = x + drift + noise
    return x.cpu().numpy()

# --- E-distance ---
def edist(A, B):
    """Standard E-distance, no subsampling (caller controls size)."""
    if len(A) == 0 or len(B) == 0: return np.nan
    return float(2*cdist(A, B).mean() - cdist(A, A).mean() - cdist(B, B).mean())

# --- PCA-50 projections per gene ---
print("Fitting PCA-50 per gene...", flush=True)
PCA_MODELS = {}
for gene in GENES:
    # fit on all variant+WT cells (raw)
    keys = [k for k in RAW if k.startswith(gene+'__')]
    all_cells = np.concatenate([RAW[k][:500] for k in keys])  # cap for speed
    pca = PCA(n_components=min(50, all_cells.shape[1]-1, all_cells.shape[0]-1), random_state=0)
    pca.fit(all_cells)
    PCA_MODELS[gene] = pca
    print(f"  {gene}: PCA-{pca.n_components_}, explained={pca.explained_variance_ratio_.sum():.2f}", flush=True)

# --- Main analysis loop ---
print("\n=== Running split-half power analysis ===", flush=True)
rows = []; t0 = time.time()
n_total = 0
for gene in GENES:
    wt_all = RAW[f'{gene}__WT']
    variants = [v for v in bench[bench.gene==gene]['variant'].values if f'{gene}__{v}' in RAW and v != 'WT']
    pca = PCA_MODELS[gene]
    
    for vi, var in enumerate(variants):
        real_all = RAW[f'{gene}__{var}']
        th_vec = theta_for(gene, var)
        
        # VAE cf cells (from G1)
        vae_key = f'{gene}__{var}'
        vae_cells = cf_vae[vae_key] if vae_key in cf_vae else None
        
        # CFM v3: generate enough cells (once, then subsample)
        if th_vec is not None:
            cfm_cells = sample_cfm(gene, th_vec, n=min(600, max(300, len(real_all))))
        else:
            cfm_cells = None
        
        for ns in N_SUBS:
            # need 2*ns real cells for split-half
            if len(real_all) < 2 * ns: continue
            if len(wt_all) < ns: continue
            
            for seed in SEEDS:
                rng = np.random.default_rng(seed)
                
                # split real cells
                perm = rng.permutation(len(real_all))
                halfA = real_all[perm[:ns]]
                halfB = real_all[perm[ns:2*ns]]
                
                # WT sample
                wt_idx = rng.choice(len(wt_all), ns, replace=len(wt_all)<ns)
                wt_sample = wt_all[wt_idx]
                
                # model samples
                vae_sub = None
                if vae_cells is not None and len(vae_cells) >= ns:
                    vi_idx = rng.choice(len(vae_cells), ns, replace=False)
                    vae_sub = vae_cells[vi_idx]
                
                cfm_sub = None
                if cfm_cells is not None and len(cfm_cells) >= ns:
                    ci_idx = rng.choice(len(cfm_cells), ns, replace=False)
                    cfm_sub = cfm_cells[ci_idx]
                
                for space in ['full', 'pca50']:
                    if space == 'pca50':
                        hA = pca.transform(halfA); hB = pca.transform(halfB)
                        wts = pca.transform(wt_sample)
                        vs = pca.transform(vae_sub) if vae_sub is not None else None
                        cs = pca.transform(cfm_sub) if cfm_sub is not None else None
                    else:
                        hA, hB, wts = halfA, halfB, wt_sample
                        vs, cs = vae_sub, cfm_sub
                    
                    d_self = edist(hA, hB)
                    d_null = edist(wts, hB)
                    d_vae = edist(vs, hB) if vs is not None else np.nan
                    d_cfm = edist(cs, hB) if cs is not None else np.nan
                    
                    rows.append({
                        'gene': gene, 'variant': var, 'n_sub': ns, 'seed': seed, 'space': space,
                        'd_self': d_self, 'd_null': d_null, 'd_vae': d_vae, 'd_cfm': d_cfm
                    })
                n_total += 1
        
        if (vi+1) % 25 == 0:
            print(f"  {gene}: {vi+1}/{len(variants)} variants [{time.time()-t0:.0f}s]", flush=True)
    print(f"  {gene} done: {len(variants)} variants [{time.time()-t0:.0f}s]", flush=True)

df = pd.DataFrame(rows)
print(f"\nTotal rows: {len(df)} [{time.time()-t0:.0f}s]", flush=True)

# --- Compute derived metrics ---
df['ratio_self_null'] = df['d_self'] / df['d_null']
df['ratio_vae_null'] = df['d_vae'] / df['d_null']
df['ratio_cfm_null'] = df['d_cfm'] / df['d_null']
df['detectable'] = df['d_self'] < df['d_null']
df['vae_beats_null'] = df['d_vae'] < df['d_null']
df['cfm_beats_null'] = df['d_cfm'] < df['d_null']
# model progress: how far toward empirical floor (1 = at floor, 0 = at null, <0 = worse than null)
mask_det = df['d_null'] > df['d_self']
df.loc[mask_det, 'vae_progress'] = (df.loc[mask_det,'d_null'] - df.loc[mask_det,'d_vae']) / (df.loc[mask_det,'d_null'] - df.loc[mask_det,'d_self'])
df.loc[mask_det, 'cfm_progress'] = (df.loc[mask_det,'d_null'] - df.loc[mask_det,'d_cfm']) / (df.loc[mask_det,'d_null'] - df.loc[mask_det,'d_self'])

OUT_DIR.mkdir(parents=True, exist_ok=True)
df.to_csv(OUT_DIR / 'split_half_power_table.csv', index=False)
print(f"Saved split_half_power_table.csv ({len(df)} rows)", flush=True)

# --- Summary by gene ---
print("\n=== SUMMARY BY GENE (n_sub=150, full space, averaged over seeds) ===", flush=True)
sub150 = df[(df.n_sub==150)&(df.space=='full')]
summary_rows = []
for gene in GENES:
    gs = sub150[sub150.gene==gene]
    if len(gs) == 0: continue
    # aggregate per variant first, then per gene
    var_agg = gs.groupby('variant').agg(
        d_self=('d_self','mean'), d_null=('d_null','mean'),
        d_vae=('d_vae','mean'), d_cfm=('d_cfm','mean'),
        detectable=('detectable','mean'),
        vae_beats=('vae_beats_null','mean'),
        cfm_beats=('cfm_beats_null','mean'),
        vae_prog=('vae_progress','mean'),
        cfm_prog=('cfm_progress','mean'))
    
    n_vars = len(var_agg)
    frac_det = (var_agg['detectable']>0.5).mean()
    frac_vae = (var_agg['vae_beats']>0.5).mean()
    frac_cfm = (var_agg['cfm_beats']>0.5).mean()
    
    # classification
    if frac_det < 0.5:
        label = 'UNDERPOWERED'
    elif frac_cfm > 0.5:
        label = 'MODEL_IMPROVES'
    elif frac_det >= 0.5 and frac_cfm <= 0.5:
        label = 'DETECTABLE_NOT_MODELED'
    else:
        label = 'MIXED'
    
    row = {'gene':gene, 'n_variants':n_vars,
           'median_d_self':var_agg['d_self'].median(),
           'median_d_null':var_agg['d_null'].median(),
           'median_ratio_self_null':(var_agg['d_self']/var_agg['d_null']).median(),
           'frac_detectable':frac_det,
           'frac_vae_beats_null':frac_vae,
           'frac_cfm_beats_null':frac_cfm,
           'median_cfm_progress':var_agg['cfm_prog'].median(),
           'label':label}
    summary_rows.append(row)
    print(f"  {gene}: {label} | detect={frac_det:.0%} | VAE<null={frac_vae:.0%} | CFM<null={frac_cfm:.0%} | D_self/D_null={row['median_ratio_self_null']:.3f}", flush=True)

sumdf = pd.DataFrame(summary_rows)
sumdf.to_csv(OUT_DIR / 'split_half_power_summary_by_gene.csv', index=False)

# --- Power curve (n_sub sweep, full space) ---
print("\n=== POWER CURVE (frac detectable by n_sub, full space) ===", flush=True)
curve_rows = []
for ns in N_SUBS:
    for space in ['full','pca50']:
        sub = df[(df.n_sub==ns)&(df.space==space)]
        for gene in GENES:
            gs = sub[sub.gene==gene]
            if len(gs)==0: continue
            vagg = gs.groupby('variant').agg(det=('detectable','mean'), cfm_b=('cfm_beats_null','mean'))
            curve_rows.append({'n_sub':ns, 'space':space, 'gene':gene,
                               'frac_detectable':(vagg['det']>0.5).mean(),
                               'frac_cfm_beats':(vagg['cfm_b']>0.5).mean(),
                               'n_variants':len(vagg)})
curvedf = pd.DataFrame(curve_rows)
curvedf.to_csv(OUT_DIR / 'split_half_power_curve.csv', index=False)
print(curvedf[curvedf.space=='full'].to_string(index=False), flush=True)

# --- PCA50 vs full comparison ---
print("\n=== PCA50 vs FULL (n_sub=150) ===", flush=True)
for space in ['full','pca50']:
    sub = df[(df.n_sub==150)&(df.space==space)]
    for gene in GENES:
        gs = sub[sub.gene==gene]
        if len(gs)==0: continue
        vagg = gs.groupby('variant').agg(det=('detectable','mean'), cfm_b=('cfm_beats_null','mean'))
        print(f"  {space:>5} {gene}: detect={((vagg['det']>0.5).mean()):.0%} | CFM<null={((vagg['cfm_b']>0.5).mean()):.0%}", flush=True)

# --- Decision document ---
print("\n=== GENERATING DECISION DOCUMENT ===", flush=True)
# compute 95% CI for key ratios
ci_text = ""
for gene in GENES:
    gs = df[(df.n_sub==150)&(df.space=='full')&(df.gene==gene)]
    if len(gs)==0: continue
    vagg = gs.groupby('variant').agg(rsn=('ratio_self_null','mean'), rcn=('ratio_cfm_null','mean'))
    rsn_vals = vagg['rsn'].dropna().values
    rcn_vals = vagg['rcn'].dropna().values
    if len(rsn_vals)>0:
        boot_rsn = [np.median(np.random.choice(rsn_vals, len(rsn_vals))) for _ in range(1000)]
        boot_rcn = [np.median(np.random.choice(rcn_vals, len(rcn_vals))) for _ in range(1000)]
        ci_text += f"  {gene}: D_self/D_null median={np.median(rsn_vals):.3f} [{np.percentile(boot_rsn,2.5):.3f}, {np.percentile(boot_rsn,97.5):.3f}]"
        ci_text += f"  |  D_cfm/D_null median={np.median(rcn_vals):.3f} [{np.percentile(boot_rcn,2.5):.3f}, {np.percentile(boot_rcn,97.5):.3f}]\n"

decision = f"""# Split-Half Power Analysis Decision

## Question
Is the E-distance "failure" (model worse than WT-null) a measurement artifact or a real model limitation?

## Method
For each gene × variant × 50 seeds:
  D_self  = E(real halfA, real halfB) — empirical floor (measurement noise + biological variation)
  D_null  = E(WT sample, real halfB) — WT-null baseline (no variant signal)
  D_model = E(model prediction, real halfB) — model counterfactual

If D_self ≥ D_null for most variants → measurement is underpowered (can't detect variant signal at all).
If D_self < D_null but D_model > D_null → signal detectable but model doesn't capture it.
If D_model < D_null → model captures some variant signal.

## Results (n_sub=150, full gene space, 50 seeds)

### Per-gene classification
"""
for _, r in sumdf.iterrows():
    decision += f"**{r['gene']}**: {r['label']}  (detect={r['frac_detectable']:.0%}, CFM<null={r['frac_cfm_beats_null']:.0%}, D_self/D_null={r['median_ratio_self_null']:.3f})\n"

decision += f"""
### 95% Bootstrap CI (median across variants)
{ci_text}

### Interpretation

"""
# generate per-gene interpretation
for _, r in sumdf.iterrows():
    if r['label'] == 'UNDERPOWERED':
        decision += f"**{r['gene']}**: Split-half D_self ≈ D_null, meaning E-distance cannot distinguish real variant replicates from WT baseline. The \"failure\" of counterfactual generation on this gene is indistinguishable from measurement noise — the assay does not carry allele-level distributional signal at this cell count.\n\n"
    elif r['label'] == 'DETECTABLE_NOT_MODELED':
        decision += f"**{r['gene']}**: Split-half D_self < D_null ({r['frac_detectable']:.0%} of variants), so variant signal IS detectable. But model predictions do not beat WT-null ({r['frac_cfm_beats_null']:.0%}). This is a real model limitation.\n\n"
    elif r['label'] == 'MODEL_IMPROVES':
        decision += f"**{r['gene']}**: Variant signal detectable AND model beats WT-null ({r['frac_cfm_beats_null']:.0%} of variants). Model progress toward empirical floor: {r['median_cfm_progress']:.2f}.\n\n"
    else:
        decision += f"**{r['gene']}**: Mixed — {r['frac_detectable']:.0%} detectable, {r['frac_cfm_beats_null']:.0%} CFM beats null.\n\n"

decision += """## Overall Decision

"""
n_underpowered = sum(1 for _, r in sumdf.iterrows() if r['label']=='UNDERPOWERED')
n_model_improves = sum(1 for _, r in sumdf.iterrows() if r['label']=='MODEL_IMPROVES')
if n_underpowered >= 2:
    decision += f"""**{n_underpowered}/4 genes are measurement-limited** (UNDERPOWERED): the assay does not carry enough allele-level distributional signal for E-distance to be a meaningful evaluation metric on those genes. The previous Track A "failure" was partly (or largely) an artifact of evaluating against an unresolvable floor.

This changes the Track A/B framing: before claiming "generation fails," the manuscript must establish per-gene detectability. E-distance is only meaningful where D_self < D_null. On underpowered genes, discrimination-based metrics (PDS, delta_cosine) remain valid because they test mean direction, not full distribution shape.
"""
else:
    decision += f"Most genes ({4-n_underpowered}/4) are detectable. The E-distance failure is model-limited, not measurement-limited.\n"

with open(OUT_DIR / 'decision_after_split_half_floor.md','w') as f:
    f.write(decision)
print(decision, flush=True)
print("POWER_DONE", flush=True)
