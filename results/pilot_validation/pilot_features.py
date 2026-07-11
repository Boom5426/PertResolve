#!/usr/bin/env python3
"""2B step 1: per-perturbation PILOT features + DISJOINT-eval rankability label.

For each perturbation, ONE permutation per seed splits cells into disjoint slices:
  pilot = first 50 cells (pilot-P uses its first P)  -> pilot-estimable features only
  eval  = next 2T cells                              -> ground-truth rankability at depth T
Pilot and eval never share a cell (the fix over the retrospective same-source predictor).
  pilot_eff_P  = || mean(pilot_P) - mean(control_P) ||           (pilot effect size)
  pilot_snr_P  = pilot_eff_P / pilot split-half noise
  rank_T       = 1 if split-half S=median(D_null)-median(D_self) > W=CI-width(D_self) on eval at depth T
Output: pilot_validation/<DS>_pilot.csv . Usage: python pilot_features.py <DS>
"""
import sys, numpy as np, pandas as pd, os
from collections import defaultdict
sys.path.insert(0, "/data/boom/NUS/VCCompass/unified")
import harness as H
import anndata as ad
from sklearn.decomposition import PCA

DS = sys.argv[1]
P_LIST = [25, 50]; T_LIST = [100, 200]; NSEED = 10; PILOTMAX = 50
WT_TAGS = ('WT', 'wt', 'WT_control')
ATLAS = {'Replogle': "/data/boom/NUS/floor_audit/ReplogleWeissman2022_K562_essential.h5ad",
         'Norman': "/data/boom/NUS/floor_audit/NormanWeissman2019_filtered.h5ad",
         'Adamson': "/data/boom/NUS/floor_audit/AdamsonWeissman2016_GSM2406681_10X010.h5ad"}
PCAND = ['perturbation', 'perturbation_type', 'target_gene', 'guide', 'condition', 'gene']
CCAND = ['control', 'Control', 'CTRL', 'ctrl', 'non-targeting', 'NT', 'DMSO', 'unperturbed']
need = PILOTMAX + 2 * max(T_LIST)


def energy(A, B):
    dab = np.sqrt(((A[:, None] - B[None]) ** 2).sum(-1)).mean()
    daa = np.sqrt(((A[:, None] - A[None]) ** 2).sum(-1)).mean()
    dbb = np.sqrt(((B[:, None] - B[None]) ** 2).sum(-1)).mean()
    return 2 * dab - daa - dbb


def load_allele(g):
    X, lab = H.load_gene(g); lab = np.asarray(lab)
    Xp = PCA(50, random_state=0).fit_transform(X.astype(np.float32))
    cells = {v: Xp[lab == v] for v in np.unique(lab) if v not in WT_TAGS}
    return cells, Xp[np.where(np.isin(lab, WT_TAGS))[0]]


def load_atlas(name):
    import scipy.sparse as sp
    A = ad.read_h5ad(ATLAS[name], backed='r'); obs = A.obs
    pcol = next(c for c in PCAND if c in obs.columns); ser = obs[pcol].astype(str)
    ctrl = next((c for c in CCAND if (ser == c).sum() > 0), None)
    rng = np.random.RandomState(0); vc = ser.value_counts()
    perts = [p for p in vc.index if p != ctrl and p != 'nan' and vc[p] >= PILOTMAX + 2 * min(T_LIST)]
    if len(perts) > 300: perts = list(rng.choice(perts, 300, replace=False))
    CAP = need + 60; rows = []; labels = []
    for p in perts:
        idx = np.where(ser.values == p)[0]
        if len(idx) > CAP: idx = rng.choice(idx, CAP, replace=False)
        rows.append(idx); labels += [p] * len(idx)
    cidx = np.where(ser.values == ctrl)[0] if ctrl else np.array([], int)
    if len(cidx) > 1500: cidx = rng.choice(cidx, 1500, replace=False)
    allidx = np.concatenate(rows + ([cidx] if len(cidx) else []))
    order = np.argsort(allidx); inv = np.argsort(order)
    Xr = A[allidx[order]].to_memory().X; A.file.close()
    if Xr.shape[1] > 40000:
        from sklearn.decomposition import TruncatedSVD
        Xr = (Xr.tocsr() if sp.issparse(Xr) else sp.csr_matrix(np.asarray(Xr)))[inv]
        Xp = TruncatedSVD(50, random_state=0).fit_transform(Xr).astype(np.float32)
    else:
        Xs = np.asarray(Xr.todense() if hasattr(Xr, 'todense') else Xr, dtype=np.float32)[inv]
        fit = rng.choice(len(Xs), min(20000, len(Xs)), replace=False)
        Xp = PCA(min(50, Xs.shape[1] - 1), random_state=0).fit(Xs[fit]).transform(Xs)
    lab = np.array(labels + (['__CTRL__'] * len(cidx) if len(cidx) else []))
    cells = {p: Xp[lab == p] for p in perts}
    ctrl_cells = Xp[lab == '__CTRL__'] if len(cidx) else Xp
    return cells, ctrl_cells


cells, ctrl = (load_allele(DS.split('_', 1)[1]) if DS.startswith('allele_') else load_atlas(DS))

out = []
for p, c in cells.items():
    n = len(c)
    if n < PILOTMAX + 2 * min(T_LIST):
        continue
    acc = defaultdict(list)
    for s in range(NSEED):
        rng = np.random.RandomState(s); ii = rng.permutation(n)
        pilot = c[ii[:PILOTMAX]]
        for P in P_LIST:
            pil = pilot[:P]; cb = ctrl[rng.choice(len(ctrl), min(P, len(ctrl)), replace=False)]
            eff = float(np.linalg.norm(pil.mean(0) - cb.mean(0))); h = P // 2
            noise = float(np.linalg.norm(pil[:h].mean(0) - pil[h:2 * h].mean(0))) + 1e-9
            acc[f'eff_{P}'].append(eff); acc[f'snr_{P}'].append(eff / noise)
        for T in T_LIST:
            if n < PILOTMAX + 2 * T:
                continue
            ev = c[ii[PILOTMAX:PILOTMAX + 2 * T]]; A, B = ev[:T], ev[T:2 * T]
            cs = ctrl[rng.choice(len(ctrl), min(T, len(ctrl)), replace=False)]
            acc[f'ds_{T}'].append(energy(A, B)); acc[f'dn_{T}'].append(energy(A, cs))
    row = dict(dataset=DS, pert=p, n_cells=n)
    for P in P_LIST:
        row[f'pilot_eff_{P}'] = np.mean(acc[f'eff_{P}']); row[f'pilot_snr_{P}'] = np.mean(acc[f'snr_{P}'])
    for T in T_LIST:
        if acc[f'ds_{T}']:
            S = np.median(acc[f'dn_{T}']) - np.median(acc[f'ds_{T}'])
            W = np.percentile(acc[f'ds_{T}'], 97.5) - np.percentile(acc[f'ds_{T}'], 2.5)
            row[f'rank_T{T}'] = int(S > W)
        else:
            row[f'rank_T{T}'] = np.nan
    out.append(row)

os.makedirs("/data/boom/NUS/pilot_validation", exist_ok=True)
pd.DataFrame(out).to_csv(f"/data/boom/NUS/pilot_validation/{DS}_pilot.csv", index=False)
print(f"{DS}: {len(out)} perturbations -> {DS}_pilot.csv "
      f"(rank_T100 pos={int(pd.DataFrame(out)['rank_T100'].sum())})")
