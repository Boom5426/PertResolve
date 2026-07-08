#!/usr/bin/env python3

"""Patch: fix v2c parsing for GSE161824 data format.
The v2c file is TAB-separated with variant columns as one-hot, 
plus 'variant', 'variant.detailed_multi', 'cell' at the end.
"""
import os, sys, gzip, time, json
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader
from scipy.spatial.distance import cdist
from scipy.io import mmread
from scipy.sparse import issparse

WORKDIR = "/data/boom/NUS/VCCompass"
os.makedirs(WORKDIR, exist_ok=True)
os.chdir(WORKDIR)

device = 'cuda' if torch.cuda.is_available() else 'cpu'
print(f"Device: {device}", flush=True)

def process_tp53_kras():
    print("\n=== TP53 + KRAS (GSE161824) ===", flush=True)
    bench = pd.read_csv('allele_perturb_bench.csv')
    tcols = ['d_hydro','d_vol','d_charge','fold_core','cat_switch','is_hotspot']
    vt = {'WT': np.zeros(6, dtype=np.float32)}
    for _, r in bench.iterrows():
        vt[r['variant']] = np.array([r[c] for c in tcols], dtype=np.float32)
    
    results = {}; shared_genes = None
    for gt in ['TP53','KRAS']:
        print(f"\n  Processing {gt}...", flush=True)
        with gzip.open(f"raw/GSE161824_A549_{gt}.processed.matrix.mtx.gz",'rt') as f:
            mat = mmread(f)
        if issparse(mat): mat = mat.toarray()
        mat = mat.T.astype(np.float32)
        genes = pd.read_csv(f"raw/GSE161824_A549_{gt}.processed.genes.csv.gz").iloc[:,0].values
        
        # v2c is TSV with variant one-hot + variant/cell columns
        v2c = pd.read_csv(f"raw/GSE161824_A549_{gt}.variants2cell.csv.gz", sep='\t')
        print(f"    v2c: {v2c.shape}, columns[-5:]: {v2c.columns[-5:].tolist()}")
        
        # Extract variant label per cell
        # 'variant' column has the assigned variant name
        vlabels = v2c['variant'].values.astype(str)
        # 'cell' column has barcodes - use row index as cell index
        # Matrix rows correspond to v2c rows
        
        # Handle 'unassigned' and 'multiple' as WT
        vlabels[vlabels == 'unassigned'] = 'WT'
        vlabels[vlabels == 'multiple'] = 'WT'
        
        # Verify dimensions match
        if mat.shape[0] != len(vlabels):
            print(f"    WARNING: matrix {mat.shape[0]} vs v2c {len(vlabels)}")
            min_n = min(mat.shape[0], len(vlabels))
            mat = mat[:min_n]
            vlabels = vlabels[:min_n]
        
        print(f"    Matrix: {mat.shape}, variants: {len(np.unique(vlabels))}")
        
        if shared_genes is None: shared_genes = set(genes)
        else: shared_genes &= set(genes)
        results[gt] = {'mat': mat, 'genes': genes, 'vlabels': vlabels}
    
    shared = sorted(shared_genes)
    print(f"\n  Shared genes: {len(shared)}", flush=True)
    
    for gt in ['TP53','KRAS']:
        g2i = {g:i for i,g in enumerate(results[gt]['genes'])}
        idx = [g2i[g] for g in shared]
        X = results[gt]['mat'][:, idx]
        mu, sd = X.mean(0), X.std(0); sd[sd==0]=1
        X = (X-mu)/sd
        th = np.zeros((len(results[gt]['vlabels']),7), dtype=np.float32)
        for i,v in enumerate(results[gt]['vlabels']):
            if v in vt: th[i,:6] = vt[v]
        results[gt].update({'X': X, 'theta': th})
        print(f"  {gt} final: {X.shape}, {len(np.unique(results[gt]['vlabels']))} variants", flush=True)
    
    np.savez_compressed('joint_arrays.npz',
        Xtp=results['TP53']['X'], vtp=results['TP53']['vlabels'],
        rtp=np.zeros((len(results['TP53']['vlabels']),1),dtype=np.float32),
        THtp=results['TP53']['theta'],
        cctp=np.zeros((len(results['TP53']['vlabels']),5),dtype=np.float32),
        gtp=np.zeros(len(results['TP53']['vlabels']),dtype=np.int8),
        Xkr=results['KRAS']['X'], vkr=results['KRAS']['vlabels'],
        rkr=np.zeros((len(results['KRAS']['vlabels']),1),dtype=np.float32),
        THkr=results['KRAS']['theta'],
        cckr=np.zeros((len(results['KRAS']['vlabels']),5),dtype=np.float32),
        gkr=np.ones(len(results['KRAS']['vlabels']),dtype=np.int8),
        shared=np.array(shared),
        theta_dims=np.array(['d_hydro','d_vol','d_charge','fold_core','cat_switch','is_hotspot','mech_sign']))
    sz = os.path.getsize('joint_arrays.npz')/1e6
    print(f"  Saved joint_arrays.npz ({sz:.0f} MB)", flush=True)

# --- Model + eval (copied from full_pipeline) ---
class MDVC(nn.Module):
    def __init__(self, gdims, nth=6, lat=16, hid=256):
        super().__init__()
        self.gdims, self.nth, self.lat = gdims, nth, lat
        self.enc_mu = nn.ModuleDict()
        self.enc_lv = nn.ModuleDict()
        for n, ng in gdims.items():
            self.enc_mu[n] = nn.Sequential(nn.Linear(ng,hid),nn.ReLU(),nn.Linear(hid,lat))
            self.enc_lv[n] = nn.Sequential(nn.Linear(ng,hid),nn.ReLU(),nn.Linear(hid,lat))
        self.th_net = nn.Sequential(nn.Linear(nth,32),nn.ReLU(),nn.Linear(32,32))
        self.dec = nn.ModuleDict()
        for n, ng in gdims.items():
            self.dec[n] = nn.Sequential(nn.Linear(lat+32,hid),nn.ReLU(),nn.Linear(hid,ng))
    def encode(self, x, d): return self.enc_mu[d](x), self.enc_lv[d](x)
    def reparam(self, mu, lv):
        return mu+torch.exp(.5*lv)*torch.randn_like(mu) if self.training else mu
    def decode(self, z, th, d): return self.dec[d](torch.cat([z, self.th_net(th)],-1))
    def forward(self, x, th, d):
        mu, lv = self.encode(x, d); z = self.reparam(mu, lv)
        return self.decode(z, th, d), mu, lv

def vloss(x, xr, mu, lv, beta=1.0):
    r = F.mse_loss(xr, x, reduction='sum')
    k = -0.5*torch.sum(1+lv-mu.pow(2)-lv.exp())
    return r+beta*k, r, k

class VDS(Dataset):
    def __init__(self, X, th):
        self.X=torch.from_numpy(X).float(); self.th=torch.from_numpy(th).float()
    def __len__(self): return len(self.X)
    def __getitem__(self, i): return self.X[i], self.th[i]

def pseudobulk(X, lab, vs):
    pb={}
    for v in vs:
        m=lab==v
        if m.sum()>=10: pb[v]=X[m].mean(0)
    return pb

def edist(A, B, ns=300):
    if len(A)>ns: A=A[np.random.choice(len(A),ns,replace=False)]
    if len(B)>ns: B=B[np.random.choice(len(B),ns,replace=False)]
    return 2*cdist(A,B,'euclidean').mean()-cdist(A,A,'euclidean').mean()-cdist(B,B,'euclidean').mean()

def train_model(model, loaders, epochs=100, beta=1.0, zero_th=False, dev='cuda'):
    opt=torch.optim.Adam(model.parameters(),lr=1e-3)
    sch=torch.optim.lr_scheduler.CosineAnnealingLR(opt,T_max=epochs,eta_min=1e-5)
    model.to(dev); t0=time.time()
    for ep in range(epochs):
        model.train(); el,nb=0,0
        its={n:iter(l) for n,l in loaders.items()}; act=set(its.keys())
        while act:
            for d in list(act):
                try: xb,tb=next(its[d])
                except StopIteration: act.discard(d); continue
                xb,tb=xb.to(dev),tb.to(dev)
                if zero_th: tb=torch.zeros_like(tb)
                xr,mu,lv=model(xb,tb,d); loss,_,_=vloss(xb,xr,mu,lv,beta)
                opt.zero_grad(); loss.backward()
                torch.nn.utils.clip_grad_norm_(model.parameters(),5.0)
                opt.step(); el+=loss.item(); nb+=1
        sch.step()
        if ep%20==0 or ep==epochs-1:
            print(f"  Ep{ep:3d}: loss={el/nb/512:.1f} [{time.time()-t0:.0f}s]", flush=True)
    return model

def evaluate(model, datasets, dgmap, bench, tcols, zero_th=False, dev='cuda'):
    model.eval(); results=[]; np.random.seed(42)
    for sid in range(1,7):
        sc=f'split{sid}_role'
        for dn,d in datasets.items():
            gene=dgmap[dn]
            tvars=bench[(bench['gene']==gene)&(bench[sc]=='test')]['variant'].tolist()
            if not tvars: continue
            X,lab,thm=d['X'],d['variant_labels'],d['theta_matrix']
            wm=lab=='WT'; wpb=X[wm].mean(0) if wm.sum()>0 else np.zeros(X.shape[1])
            wc=X[wm] if wm.sum()>0 else None
            tpb=pseudobulk(X,lab,tvars)
            trvars=bench[(bench['gene']==gene)&(bench[sc]=='train')]['variant'].tolist()
            trpb=pseudobulk(X,lab,trvars+['WT'])
            trth={}
            for tv in trvars+['WT']:
                tr=bench[(bench['gene']==gene)&(bench['variant']==tv)]
                trth[tv]=np.array([tr.iloc[0][c] for c in tcols]) if len(tr)>0 else np.zeros(6)
            for v in tvars:
                if v not in tpb: continue
                vm=lab==v; xt=X[vm]; thv=thm[vm]
                if zero_th: thv=np.zeros_like(thv)
                with torch.no_grad():
                    xr,_,_=model(torch.from_numpy(xt).float().to(dev),torch.from_numpy(thv).float().to(dev),dn)
                    ppb=xr.cpu().numpy().mean(0)
                td=tpb[v]-wpb; pd_=ppb-wpb
                cos_v=float(np.dot(td,pd_)/(np.linalg.norm(td)*np.linalg.norm(pd_)+1e-10))
                ed_v=float(edist(xt,xr.cpu().numpy(),200))
                ed_wt=float(edist(xt,wc,200)) if wc is not None else float('nan')
                vr=bench[(bench['gene']==gene)&(bench['variant']==v)]
                vth=np.array([vr.iloc[0][c] for c in tcols]) if len(vr)>0 else np.zeros(6)
                bd,bn=float('inf'),'WT'
                for tv in trvars+['WT']:
                    dist=np.linalg.norm(vth-trth.get(tv,np.zeros(6)))
                    if dist<bd and tv!=v: bd,bn=dist,tv
                nm=lab==bn; npb_=X[nm].mean(0) if nm.sum()>0 else wpb; nd_=npb_-wpb
                cos_nn=float(np.dot(td,nd_)/(np.linalg.norm(td)*np.linalg.norm(nd_)+1e-10))
                ed_nn=float(edist(xt,X[nm],200)) if nm.sum()>0 else ed_wt
                gm=X.mean(0); gd=gm-wpb
                cos_gm=float(np.dot(td,gd)/(np.linalg.norm(td)*np.linalg.norm(gd)+1e-10))
                results.append({'split':sid,'gene':gene,'dataset':dn,'variant':v,'n_cells':int(vm.sum()),
                    'cos_model':cos_v,'ed_model':ed_v,'cos_nn':cos_nn,'ed_nn':ed_nn,'cos_gene_mean':cos_gm,'ed_wt_null':ed_wt})
    return pd.DataFrame(results)

# ============ MAIN ============
if __name__ == '__main__':
    t0_total = time.time()
    print("="*60); print("VCCompass Pipeline v2 (fixed v2c parsing)"); print("="*60, flush=True)
    
    if not os.path.exists('joint_arrays.npz'):
        process_tp53_kras()
    else:
        print("joint_arrays.npz exists", flush=True)
    
    print("\nLoading data...", flush=True)
    joint = np.load('joint_arrays.npz', allow_pickle=True)
    bench = pd.read_csv('allele_perturb_bench.csv')
    tcols = ['d_hydro','d_vol','d_charge','fold_core','cat_switch','is_hotspot']
    vt = {'WT':np.zeros(6,dtype=np.float32)}
    for _,r in bench.iterrows(): vt[r['variant']]=np.array([r[c] for c in tcols],dtype=np.float32)
    
    datasets = {
        'tp53':{'X':joint['Xtp'],'variant_labels':joint['vtp'],'gene':'TP53','n_genes':int(joint['Xtp'].shape[1]),
                'theta_matrix':joint['THtp'][:,:6]},
        'kras':{'X':joint['Xkr'],'variant_labels':joint['vkr'],'gene':'KRAS','n_genes':int(joint['Xkr'].shape[1]),
                'theta_matrix':joint['THkr'][:,:6]},
    }
    dgmap = {'tp53':'TP53','kras':'KRAS'}
    
    # Add GATA1 if available
    if os.path.exists('gata1_arrays.npz'):
        gata1=np.load('gata1_arrays.npz',allow_pickle=True)
        datasets['gata1']={'X':gata1['X'],'variant_labels':gata1['cell_variants'],'gene':'GATA1',
                           'n_genes':int(gata1['X'].shape[1])}
        dgmap['gata1']='GATA1'
        d=datasets['gata1']; lab=d['variant_labels']
        th=np.zeros((len(lab),6),dtype=np.float32)
        for i,v in enumerate(lab):
            if v in vt: th[i]=vt[v]
        d['theta_matrix']=th
        print("  GATA1 loaded", flush=True)
    
    # Add JAK1 if available
    if os.path.exists('jak1_arrays.npz'):
        jak1=np.load('jak1_arrays.npz',allow_pickle=True)
        datasets['jak1']={'X':jak1['X'],'variant_labels':jak1['variant_labels'],'gene':'JAK1',
                          'n_genes':int(jak1['X'].shape[1])}
        dgmap['jak1']='JAK1'
        d=datasets['jak1']; lab=d['variant_labels']
        th=np.zeros((len(lab),6),dtype=np.float32)
        for i,v in enumerate(lab):
            if v in vt: th[i]=vt[v]
        d['theta_matrix']=th
        print("  JAK1 loaded", flush=True)
    
    gdims={n:d['n_genes'] for n,d in datasets.items()}
    total=sum(d['X'].shape[0] for d in datasets.values())
    print(f"  Total: {total:,} cells, {len(datasets)} datasets, genes: {list(gdims.keys())}", flush=True)
    
    # Train loaders
    trvars={}
    for _,r in bench.iterrows():
        if r['split1_role']=='train': trvars.setdefault(r['gene'],set()).add(r['variant'])
    for g in dgmap.values(): trvars.setdefault(g,set()).add('WT')
    train_loaders={}
    for dn,d in datasets.items():
        gene=dgmap[dn]
        mask=np.isin(d['variant_labels'],list(trvars.get(gene,set())))
        ds=VDS(d['X'][mask],d['theta_matrix'][mask])
        train_loaders[dn]=DataLoader(ds,batch_size=512,shuffle=True,drop_last=True,num_workers=4,pin_memory=True)
        print(f"  {dn} train: {len(ds)} cells", flush=True)
    
    all_dfs=[]
    for name, zero_th, shuffle_th, ckpt in [
        ('VCCompass_6dim',False,False,'model_full.pt'),
        ('Dosage_only',True,False,'model_dosage.pt'),
        ('Random_theta',False,True,'model_random.pt'),
    ]:
        print(f"\n{'='*50}\nTraining: {name}\n{'='*50}", flush=True)
        m=MDVC(gdims).to(device)
        if shuffle_th:
            np.random.seed(123)
            sl={}
            for dn,d in datasets.items():
                gene=dgmap[dn]; mask=np.isin(d['variant_labels'],list(trvars.get(gene,set())))
                lab=d['variant_labels'][mask]; thm=d['theta_matrix'][mask].copy()
                uv=np.unique(lab); perm=np.random.permutation(len(uv))
                vm={uv[i]:uv[perm[i]] for i in range(len(uv))}
                for i,v in enumerate(lab): thm[i]=vt.get(vm.get(v,v),np.zeros(6,dtype=np.float32))
                sl[dn]=DataLoader(VDS(d['X'][mask],thm),batch_size=512,shuffle=True,drop_last=True,num_workers=4,pin_memory=True)
            m=train_model(m,sl,epochs=100,dev=device)
        else:
            m=train_model(m,train_loaders,epochs=100,zero_th=zero_th,dev=device)
        torch.save(m.state_dict(),ckpt)
        print(f"\nEvaluating {name}...", flush=True)
        df=evaluate(m,datasets,dgmap,bench,tcols,zero_th=zero_th,dev=device)
        df['model_name']=name; all_dfs.append(df)
        for sid in range(1,7):
            s=df[df['split']==sid]
            if len(s)>0: print(f"  Split {sid}: cos={s['cos_model'].mean():.3f} (n={len(s)})", flush=True)
    
    # Biophys-only
    print(f"\n{'='*50}\nTraining: Biophys_only_3dim\n{'='*50}", flush=True)
    bl={}
    for dn,d in datasets.items():
        gene=dgmap[dn]; mask=np.isin(d['variant_labels'],list(trvars.get(gene,set())))
        thm=d['theta_matrix'][mask].copy(); thm[:,3:]=0
        bl[dn]=DataLoader(VDS(d['X'][mask],thm),batch_size=512,shuffle=True,drop_last=True,num_workers=4,pin_memory=True)
    mb=MDVC(gdims).to(device)
    mb=train_model(mb,bl,epochs=100,dev=device)
    torch.save(mb.state_dict(),'model_biophys.pt')
    dsb={}
    for dn,d in datasets.items():
        dsb[dn]=dict(d); dsb[dn]['theta_matrix']=d['theta_matrix'].copy(); dsb[dn]['theta_matrix'][:,3:]=0
    dfb=evaluate(mb,dsb,dgmap,bench,tcols,dev=device)
    dfb['model_name']='Biophys_only_3dim'; all_dfs.append(dfb)
    
    all_results=pd.concat(all_dfs,ignore_index=True)
    all_results.to_csv('all_metrics.csv',index=False)
    
    print("\n"+"="*80); print("FINAL SUMMARY"); print("="*80, flush=True)
    pivot=all_results.groupby(['model_name','split'])['cos_model'].mean().unstack()
    pivot['overall']=all_results.groupby('model_name')['cos_model'].mean()
    print("\nMean cosine (direction):");print(pivot.to_string())
    print("\n\nPer-gene:");gp=all_results.groupby(['model_name','gene'])['cos_model'].mean().unstack();print(gp.to_string())
    
    summary={'overall':all_results.groupby('model_name')['cos_model'].mean().to_dict(),
             'per_split':pivot.to_dict(),'per_gene':gp.to_dict(),
             'total_variants':len(all_results[all_results['model_name']=='VCCompass_6dim']),
             'total_time_min':(time.time()-t0_total)/60}
    with open('summary.json','w') as f: json.dump(summary,f,indent=2,default=str)
    print(f"\nTotal: {(time.time()-t0_total)/60:.1f} min"); print("DONE", flush=True)
