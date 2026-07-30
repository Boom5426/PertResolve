#!/usr/bin/env python3
"""T1/T2: second-stage rankability audit — 4 metric families x multi-space x n-grid.
Reads prep/{dataset}.pkl. For each perturbation x metric x space x n x seed:
  D_self = split-half floor (metric between two halves of the SAME perturbation)
  D_null = control-vs-perturbation signal (same n)
  S = median(D_null) - median(D_self)       (signal room)
  W = CI_hi(D_self) - CI_lo(D_self)          (sampling width, 95%)
  rankable = S > W ; ratio = S/W
Metrics:
  edist   : energy distance (euclidean)            [distribution]
  mmd     : MMD^2 with RBF, bandwidth=median-heuristic on control (fixed per dataset x space)
  dcos    : 1 - delta-cosine of pseudobulk mean-shift vs control  (LOWER=closer; we report signal form)
  pds     : nearest-condition ranking accuracy (discrimination) - dataset-level, not split-half
Significance: perturbation kept if pseudobulk shift E-test-significant (reuse: effect_size>0 and n>=2*min_n).

Usage:
  rankability_audit.py [-h] [--prep PREP] --out OUT dataset

  dataset      name of the prepared input; the script reads <prep>/<dataset>.pkl
  --prep PREP  directory holding the prepared atlas .pkl inputs; may instead be
               given through the ALLELEPERTURB_ATLAS_DIR environment variable
  --out OUT    directory the <dataset>_rankability.csv table is written to
               (required; never defaults into the repository's results/)

  K and NJOBS stay environment variables, as before.

  Example:
    ALLELEPERTURB_ATLAS_DIR=/path/to/floor_audit/prep \\
      python rankability_audit.py replogle --out /path/to/floor_audit_results
"""
import argparse
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from alleleperturb.paths import ATLAS_ENV_VAR, require_inputs, resolve_base

# Arguments are parsed before the scientific stack is imported, so --help works
# without scipy/joblib installed and a bad --prep fails before any load.
_parser = argparse.ArgumentParser(
    description="Second-stage rankability audit: 4 metric families x multi-space x n-grid.")
_parser.add_argument("dataset", help="name of the prepared input; reads <prep>/<dataset>.pkl")
_parser.add_argument("--prep", default=None,
                     help="directory holding the prepared atlas .pkl inputs "
                          f"(or set {ATLAS_ENV_VAR})")
_parser.add_argument("--out", required=True,
                     help="directory the <dataset>_rankability.csv table is written to")
_args = _parser.parse_args()

import time, pickle, numpy as np, pandas as pd  # noqa: E402
from scipy.spatial.distance import cdist  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

PREP=resolve_base(_args.prep, what="the prepared atlas directory",
                  env_var=ATLAS_ENV_VAR, flag="--prep")
OUT=Path(_args.out)
DATASET=_args.dataset
K=int(os.environ.get('K','50'))
NJOBS=int(os.environ.get('NJOBS','20'))

PREP_PKL=PREP/f"{DATASET}.pkl"
require_inputs(PREP_PKL)
d=pickle.load(open(PREP_PKL,"rb"))
LEVEL=d["level"]
NGRID=[25,50,100,200,400,800] if LEVEL=="gene" else [30,50,100,150,250,365]
SPACES=[s for s in ["pca","gene","pathway"] if s in d["pools"] and d["pools"][s].get(d["perts"][0] if d["perts"] else None) is not None]
# robust space list: a space is usable if control has cols>0
SPACES=[s for s in ["pca","gene","pathway"] if s in d["control"] and d["control"][s].shape[1]>0]
print(f"[{DATASET}] level={LEVEL} spaces={SPACES} nperts={len(d['perts'])} grid={NGRID}",flush=True)

# ---------- metric primitives ----------
def edist(A,B):
    return 2*cdist(A,B).mean()-cdist(A,A).mean()-cdist(B,B).mean()
def mmd2(A,B,gamma):
    Kaa=np.exp(-gamma*cdist(A,A,'sqeuclidean')); Kbb=np.exp(-gamma*cdist(B,B,'sqeuclidean')); Kab=np.exp(-gamma*cdist(A,B,'sqeuclidean'))
    return Kaa.mean()+Kbb.mean()-2*Kab.mean()
def dcos_signal(A,B,ref_mu):
    # delta = mean(A)-ref ; delta'=mean(B)-ref ; return cosine (higher=more consistent direction)
    da=A.mean(0)-ref_mu; db=B.mean(0)-ref_mu
    na=np.linalg.norm(da); nb=np.linalg.norm(db)
    if na<1e-9 or nb<1e-9: return 0.0
    return float(np.dot(da,db)/(na*nb))

def gamma_for(space, X):
    # median heuristic on a control subsample
    sub=X[np.random.RandomState(0).choice(len(X),min(500,len(X)),replace=False)]
    md=np.median(cdist(sub,sub,'sqeuclidean')); md=md if md>0 else 1.0
    return 1.0/md

rows=[]
t0=time.time()
for space in SPACES:
    Xc=d["control"][space]
    gamma=gamma_for(space,Xc)
    ref_mu=Xc.mean(0)
    pools=d["pools"][space]
    def work(p):
        Xp=pools[p]; n=len(Xp); out=[]
        eff=d["effect_size"][p]; ncell=d["n_cells"][p]
        for nw in NGRID:
            if n<2*nw or len(Xc)<nw: continue
            # --- edist, mmd, dcos: split-half floor + null ---
            ds_e=[];dn_e=[];ds_m=[];dn_m=[];ds_c=[];dn_c=[]
            for s in range(K):
                r=np.random.RandomState(s*997+nw)
                perm=r.permutation(n); A=Xp[perm[:nw]]; B=Xp[perm[nw:2*nw]]
                ci=r.choice(len(Xc),nw,replace=False); Cc=Xc[ci]
                cj=r.choice(len(Xc),nw,replace=False); Cc2=Xc[cj]  # 2nd ctrl draw for null-self
                pj=r.choice(n,nw,replace=False); Pn=Xp[pj]
                # edist
                ds_e.append(edist(A,B)); dn_e.append(edist(Cc,Pn))
                # mmd
                ds_m.append(mmd2(A,B,gamma)); dn_m.append(mmd2(Cc,Pn,gamma))
                # dcos: SELF = cosine of two half-shifts (should be ~1 if signal); NULL = cosine ctrl-half shift
                ds_c.append(dcos_signal(A,B,ref_mu))
                dn_c.append(dcos_signal(Cc, Cc2, ref_mu))  # control has no direction -> ~0
            for mname,ds,dn,sign in [("edist",ds_e,dn_e,+1),("mmd",ds_m,dn_m,+1),("dcos",ds_c,dn_c,-1)]:
                ds=np.array(ds);dn=np.array(dn)
                if sign>0:
                    # distribution metrics: D_self is FLOOR (noise), D_null is SIGNAL. S=null-self, W=CI(self)
                    lo,hi=np.percentile(ds,[2.5,97.5]); S=np.median(dn)-np.median(ds); W=hi-lo
                else:
                    # dcos: D_self is SIGNAL (consistency ~1), D_null ~0. S = self-null, W=CI(self)
                    lo,hi=np.percentile(ds,[2.5,97.5]); S=np.median(ds)-np.median(dn); W=hi-lo
                out.append(dict(dataset=DATASET,level=LEVEL,space=space,metric=mname,perturbation=p,
                                n_work=nw,n_cells=ncell,effect_size=eff,
                                D_self=float(np.median(ds)),D_self_lo=float(lo),D_self_hi=float(hi),
                                D_null=float(np.median(dn)),S=float(S),W=float(W),
                                rankable=bool(S>W),ratio=float(S/W) if W>0 else np.nan))
        return out
    # threading backend: shares memory (no pickling of large X); metric ops release GIL via numpy/scipy
    res=Parallel(n_jobs=NJOBS, backend="threading")(delayed(work)(p) for p in d["perts"])
    flat=[r for sub in res for r in sub]
    rows.extend(flat)
    print(f"  [{space}] {len(flat)} rows [{time.time()-t0:.0f}s]",flush=True)

# ---------- PDS: nearest-condition ranking (discrimination), per space x n ----------
# For each n: build per-condition pseudobulk from a random nw-cell draw; measure whether a held-out
# nw-cell draw of the SAME pert ranks nearest to its own condition among all conditions.
def pds_at(space,nw,seed):
    pools=d["pools"][space]; perts=[p for p in d["perts"] if len(pools[p])>=2*nw]
    if len(perts)<5: return np.nan
    r=np.random.RandomState(seed)
    refs={}; qrys={}
    for p in perts:
        idx=r.permutation(len(pools[p]))
        refs[p]=pools[p][idx[:nw]].mean(0); qrys[p]=pools[p][idx[nw:2*nw]].mean(0)
    R=np.array([refs[p] for p in perts]); Q=np.array([qrys[p] for p in perts])
    D=cdist(Q,R)  # query x ref
    # rank of self among refs (0=best). PDS = 1 - mean(rank)/(n-1)
    ranks=[np.where(np.argsort(D[i])==i)[0][0] for i in range(len(perts))]
    return 1.0 - np.mean(ranks)/(len(perts)-1)

for space in SPACES:
    for nw in NGRID:
        vals=[pds_at(space,nw,s) for s in range(min(K,20))]
        vals=[v for v in vals if not np.isnan(v)]
        if not vals: continue
        # PDS is dataset-level: store one row per (space,n) with perturbation='__POOL__'
        obs=np.mean(vals); lo,hi=np.percentile(vals,[2.5,97.5])
        # null: shuffle labels -> PDS ~ 0.5
        rows.append(dict(dataset=DATASET,level=LEVEL,space=space,metric="pds",perturbation="__POOL__",
                         n_work=nw,n_cells=-1,effect_size=np.nan,
                         D_self=obs,D_self_lo=lo,D_self_hi=hi,D_null=0.5,
                         S=obs-0.5,W=hi-lo,rankable=bool((obs-0.5)>(hi-lo)),ratio=float((obs-0.5)/(hi-lo)) if hi>lo else np.nan))
    print(f"  [{space}] PDS done [{time.time()-t0:.0f}s]",flush=True)

df=pd.DataFrame(rows)
OUT.mkdir(parents=True, exist_ok=True)
df.to_csv(OUT/f"{DATASET}_rankability.csv",index=False)
print(f"[{DATASET}] saved {len(df)} rows, {(time.time()-t0)/60:.1f} min. {DATASET}_RANK_DONE",flush=True)
