#!/usr/bin/env python3
"""S3: R3 detection-limit empirical power curve.
JAK1 (only gene with full-distribution signal) subsampled DOWN from raw arrays
to trace where D_self/D_null crosses into the noise floor -> critical detection n.
TP53/KRAS/GATA1 annotated at their current n (under-powered zone).
Dual space: gene-space (2000/939D) + PCA-50."""
import os, time, json
import numpy as np, pandas as pd
from scipy.spatial.distance import cdist
from sklearn.decomposition import PCA
os.chdir('/data/boom/NUS/processed-data directory')

def edist(A,B,nsub=None,seed=0):
    """standard energy distance E = 2*mean(d_AB) - mean(d_AA) - mean(d_BB)."""
    rng=np.random.RandomState(seed)
    if nsub is not None:
        if len(A)>nsub: A=A[rng.choice(len(A),nsub,replace=False)]
        if len(B)>nsub: B=B[rng.choice(len(B),nsub,replace=False)]
    dab=cdist(A,B).mean(); daa=cdist(A,A).mean(); dbb=cdist(B,B).mean()
    return 2*dab-daa-dbb

# ---- load all genes raw ----
joint=np.load('joint_arrays.npz',allow_pickle=True)
g1=np.load('gata1_arrays.npz',allow_pickle=True)
j1=np.load('jak1_arrays.npz',allow_pickle=True)
datasets={
 'TP53':{'X':joint['Xtp'].astype(np.float32),'lab':joint['vtp']},
 'KRAS':{'X':joint['Xkr'].astype(np.float32),'lab':joint['vkr']},
 'GATA1':{'X':g1['X'].astype(np.float32),'lab':g1['cell_variants']},
 'JAK1':{'X':j1['X'].astype(np.float32),'lab':j1['variant_labels']},
}
# PCA-50 per gene
for g,d in datasets.items():
    p=PCA(n_components=50,random_state=0).fit(d['X'])
    d['Xpca']=p.transform(d['X']).astype(np.float32)
    print(f"{g}: X{d['X'].shape} PCA-50 explained={p.explained_variance_ratio_.sum():.3f}",flush=True)

NSEED=50; NSUB_EDIST=None  # use all cells in each half for edist (subsample handled by n_per_split)
t0=time.time()

# ============ PART 1: JAK1 downward-subsample power curve ============
print("\n=== JAK1 DOWNWARD SUBSAMPLE POWER CURVE ===",flush=True)
N_PER_SPLIT=[30,50,100,150,250,365]
jak=datasets['JAK1']; lab=jak['lab']
wt_all_g=jak['X'][lab=='WT']; wt_all_p=jak['Xpca'][lab=='WT']
# variants with enough cells
import collections
cnt=collections.Counter(lab)
jvars=[v for v,n in cnt.items() if v!='WT' and n>=60]  # need >=2*30
print(f"JAK1: {len(jvars)} variants with >=60 cells, WT={len(wt_all_g)}",flush=True)

rows=[]
for space,Xkey,wt_all in [('gene','X',wt_all_g),('pca','Xpca',wt_all_p)]:
    Xg=jak[Xkey]
    for n in N_PER_SPLIT:
        for v in jvars:
            vcells=Xg[lab==v]
            if len(vcells)<2*n: continue  # need 2 non-overlapping halves of size n
            if len(wt_all)<n: continue
            for seed in range(NSEED):
                rng=np.random.RandomState(seed*1000+n)
                idx=rng.permutation(len(vcells))
                A=vcells[idx[:n]]; B=vcells[idx[n:2*n]]  # two disjoint halves
                wi=rng.choice(len(wt_all),n,replace=False); W=wt_all[wi]
                d_self=edist(A,B); d_null=edist(W,B)
                rows.append({'space':space,'gene':'JAK1','variant':v,'n_per_split':n,'seed':seed,
                             'd_self':d_self,'d_null':d_null,'ratio':d_self/d_null if d_null>0 else np.nan})
        done=[r for r in rows if r['space']==space and r['n_per_split']==n]
        if done:
            rr=np.array([x['ratio'] for x in done]); 
            print(f"  {space} n={n}: {len(done)} obs, median ratio={np.nanmedian(rr):.3f} [{time.time()-t0:.0f}s]",flush=True)

curve=pd.DataFrame(rows)
curve.to_csv('results/r3_power_curve.csv',index=False)

# ---- critical n: where D_self/D_null CI crosses 1.0 (bootstrap over seeds+variants) ----
print("\n=== CRITICAL DETECTION n (JAK1, D_self/D_null CI crosses 1.0) ===",flush=True)
crit_rows=[]
for space in ['gene','pca']:
    for n in N_PER_SPLIT:
        sub=curve[(curve.space==space)&(curve.n_per_split==n)]
        if len(sub)==0: continue
        r=sub['ratio'].dropna().values
        # bootstrap CI of median ratio
        rng=np.random.RandomState(0); bs=[np.median(rng.choice(r,len(r),replace=True)) for _ in range(2000)]
        lo,hi=np.percentile(bs,2.5),np.percentile(bs,97.5)
        med=np.median(r)
        detectable = hi < 1.0  # signal detectable if CI entirely below 1.0
        crit_rows.append({'space':space,'n_per_split':n,'n_obs':len(sub),'median_ratio':med,
                          'ci_lo':lo,'ci_hi':hi,'detectable':detectable})
        print(f"  {space} n={n}: ratio={med:.3f} [{lo:.3f},{hi:.3f}] detectable={detectable}",flush=True)
crit=pd.DataFrame(crit_rows); crit.to_csv('results/r3_critical_n.csv',index=False)
# critical n = smallest n where detectable=True (CI_hi<1.0)
for space in ['gene','pca']:
    det=crit[(crit.space==space)&(crit.detectable)]
    cn=int(det['n_per_split'].min()) if len(det)>0 else None
    print(f"  CRITICAL n ({space}): {cn}",flush=True)

# ============ PART 2: TP53/KRAS/GATA1 at their current n (under-powered zone) ============
print("\n=== TP53/KRAS/GATA1 AT CURRENT n (annotate under-powered zone) ===",flush=True)
under_rows=[]
for gene in ['TP53','KRAS','GATA1']:
    d=datasets[gene]; lab=d['lab']
    for space,Xkey in [('gene','X'),('pca','Xpca')]:
        Xg=d[Xkey]; wt=Xg[lab=='WT']
        gvars=[v for v,n in collections.Counter(lab).items() if v!='WT' and n>=40]
        ratios=[]
        for v in gvars:
            vc=Xg[lab==v]; n=min(len(vc)//2, len(wt))
            if n<20: continue
            for seed in range(20):
                rng=np.random.RandomState(seed*7+n)
                idx=rng.permutation(len(vc)); A=vc[idx[:n]]; B=vc[idx[n:2*n]]
                wi=rng.choice(len(wt),n,replace=False); W=wt[wi]
                ratios.append(edist(A,B)/max(edist(W,B),1e-9))
        if ratios:
            rng=np.random.RandomState(0); bs=[np.median(rng.choice(ratios,len(ratios),replace=True)) for _ in range(2000)]
            lo,hi=np.percentile(bs,2.5),np.percentile(bs,97.5)
            median_n=int(np.median([ (lab==v).sum() for v in gvars]))
            under_rows.append({'gene':gene,'space':space,'median_ratio':np.median(ratios),
                               'ci_lo':lo,'ci_hi':hi,'typical_n_per_variant':median_n,
                               'detectable':hi<1.0})
            print(f"  {gene} {space}: ratio={np.median(ratios):.3f} [{lo:.3f},{hi:.3f}] n~{median_n} detectable={hi<1.0}",flush=True)
under=pd.DataFrame(under_rows); under.to_csv('results/r3_underpowered_genes.csv',index=False)

# ============ PART 3: effect-size x required-n prescription ============
print("\n=== EFFECT-SIZE x REQUIRED-n PRESCRIPTION ===",flush=True)
# For each gene, compute pseudobulk effect size (mean |delta| vs WT, normalized) and map to required n
presc_rows=[]
for gene in ['TP53','KRAS','GATA1','JAK1']:
    d=datasets[gene]; lab=d['lab']; X=d['X']; wt=X[lab=='WT']; wtm=wt.mean(0)
    gvars=[v for v,n in collections.Counter(lab).items() if v!='WT' and n>=40]
    # effect size = median over variants of ||pseudobulk delta|| / within-WT scale
    wt_scale=np.linalg.norm(wt.std(0))+1e-9
    effs=[np.linalg.norm(X[lab==v].mean(0)-wtm)/wt_scale for v in gvars]
    eff=np.median(effs)
    presc_rows.append({'gene':gene,'pseudobulk_effect_size':eff,'n_variants':len(gvars)})
    print(f"  {gene}: pseudobulk effect size (norm delta / WT scale) = {eff:.3f}, {len(gvars)} variants",flush=True)
presc=pd.DataFrame(presc_rows)
# from JAK1 curve: at effect size eff_jak1, critical n = cn_gene. Extrapolate required n ~ (eff_jak1/eff)^2 * cn
# (E-distance scales with effect^2 for detection; naive power scaling)
jak_eff=presc[presc.gene=='JAK1']['pseudobulk_effect_size'].values[0]
det_gene=crit[(crit.space=='gene')&(crit.detectable)]
jak_crit_n=int(det_gene['n_per_split'].min()) if len(det_gene)>0 else 365
presc['jak1_effect']=jak_eff; presc['jak1_critical_n']=jak_crit_n
presc['required_n_extrap']=(jak_eff/presc['pseudobulk_effect_size'])**2 * jak_crit_n
presc['required_n_extrap']=presc['required_n_extrap'].round(0).astype(int)
presc.to_csv('results/r3_required_n_prescription.csv',index=False)
print("\nPrescription (required cells/variant to detect allele full-distribution effect):",flush=True)
print(presc[['gene','pseudobulk_effect_size','required_n_extrap']].to_string(index=False),flush=True)

json.dump({'jak1_critical_n_gene':jak_crit_n,
           'jak1_critical_n_pca':int(crit[(crit.space=='pca')&(crit.detectable)]['n_per_split'].min()) if len(crit[(crit.space=='pca')&(crit.detectable)])>0 else None,
           'runtime_min':(time.time()-t0)/60},
          open('results/r3_decision.json','w'),indent=2)
print(f"\nTotal: {(time.time()-t0)/60:.1f} min",flush=True); print("R3_DONE",flush=True)
