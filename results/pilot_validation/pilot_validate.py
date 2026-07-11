#!/usr/bin/env python3
"""2B step 2: prospective validation from pilot features (disjoint pilot->eval).

Pools all <DS>_pilot.csv. For each target depth T:
  - LEARNED LODO: leave-one-dataset-out logistic on pilot_eff_50 -> rank_T (prospective AUROC).
  - MECHANISTIC: pilot_snr_50 (proportional to rho) as the score, no training (pooled AUROC).
  - CALIBRATION: pooled out-of-fold predicted probability vs observed rankable rate.
All pilot features come from cells DISJOINT from the eval cells that define the label, so this
is a genuine prospective test (unlike the retrospective same-source Fig 6 predictor).
"""
import glob, numpy as np, pandas as pd, json
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score

df = pd.concat([pd.read_csv(f) for f in sorted(glob.glob("/data/boom/NUS/pilot_validation/*_pilot.csv"))],
               ignore_index=True)
print("datasets:", list(df.dataset.unique()), "n_perturbations=", len(df))
summary = {}
for T in [100, 200]:
    lab = f'rank_T{T}'; feat = 'pilot_eff_50'
    d = df.dropna(subset=[feat, lab]).copy(); d[lab] = d[lab].astype(int)
    print(f"\n=== target depth T={T}  (n={len(d)}, positives={d[lab].sum()}) ===")
    oof_p, oof_y = [], []; per = {}
    for ho in d.dataset.unique():
        tr = d[d.dataset != ho]; te = d[d.dataset == ho]
        if tr[lab].nunique() < 2:
            continue
        sc = StandardScaler().fit(tr[[feat]])
        m = LogisticRegression(max_iter=500).fit(sc.transform(tr[[feat]]), tr[lab])
        p = m.predict_proba(sc.transform(te[[feat]]))[:, 1]
        oof_p += list(p); oof_y += list(te[lab].values)
        per[ho] = (round(roc_auc_score(te[lab], p), 3) if te[lab].nunique() == 2 else None,
                   len(te), int(te[lab].sum()))
    print("  LEARNED LODO per-dataset (prospective AUROC | n | positives):")
    for ho, (a, ntot, npos) in per.items():
        print(f"    {ho:14s} AUROC={a}  n={ntot}  pos={npos}")
    ev = [a for a, _, _ in per.values() if a is not None]
    mean_auc = float(np.mean(ev)) if ev else None
    print(f"  mean LEARNED LODO AUROC over evaluable datasets = {mean_auc}")
    d2 = df.dropna(subset=['pilot_snr_50', lab])
    mech = round(roc_auc_score(d2[lab].astype(int), d2['pilot_snr_50']), 3) if d2[lab].nunique() == 2 else None
    print(f"  MECHANISTIC pilot_snr_50 AUROC (pooled, no training) = {mech}")
    oof_p, oof_y = np.array(oof_p), np.array(oof_y)
    calib = []
    if len(oof_y):
        edges = np.quantile(oof_p, np.linspace(0, 1, 6))
        print("  CALIBRATION (predicted-prob bin -> observed rankable rate):")
        for i in range(5):
            msk = (oof_p >= edges[i]) & (oof_p <= edges[i + 1])
            if msk.sum():
                print(f"    pred[{edges[i]:.2f},{edges[i+1]:.2f}] observed={oof_y[msk].mean():.2f}  n={int(msk.sum())}")
                calib.append([round(float(edges[i]), 3), round(float(edges[i+1]), 3),
                              round(float(oof_y[msk].mean()), 3), int(msk.sum())])
    summary[f"T{T}"] = dict(n=len(d), positives=int(d[lab].sum()),
                            learned_lodo_auroc=per, mean_learned_auroc=mean_auc,
                            mechanistic_auroc=mech, calibration=calib)
json.dump(summary, open("/data/boom/NUS/pilot_validation/pilot_validation_summary.json", "w"), indent=2)
print("\nsaved -> pilot_validation/pilot_validation_summary.json")
