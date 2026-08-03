#!/usr/bin/env python3
"""Is the high split-half PDS of GSE311877 real allele resolution, or just
severity ranking along one shared axis?

audit_gse311877.py showed split-half PDS=0.96 at chance 0.50. That number alone
does not license calling TP63 a "measurement-resolved" system, because two very
different structures produce it:

  (a) each allele has its own transcriptional direction   -> real allele resolution
  (b) every allele moves along ONE shared failure axis, and alleles differ only
      in how far they move                                 -> severity ranking

(b) is the null this script tries to establish, because AllelePerturb's existing
residual decomposition already found the shared-direction structure in the
current data. Three tests, all on the same split-half machinery:

  T1 scalar-magnitude baseline: identify an allele from ||delta|| alone (1 number).
     If that already scores well above chance, the signal is severity.
  T2 residual identification: project out the top-k principal components of the
     delta space and re-run identification on what is left.
  T3 plate-position confound: allele is aliased to plate COLUMN in this layout
     (replicate = row). Report how much of the structure survives once the
     shared axis is gone, and flag the residual ambiguity honestly.

No model is fitted. Gene selection and PCA are computed on half A only where
they touch the identification test, so half B never informs its own scoring.
"""
from __future__ import annotations

import gzip
import json
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
RAW = Path(os.environ.get("GSE311877_DIR", REPO / "data_external" / "GSE311877"))
OUT = Path(os.environ.get("DATASET_SCREEN_OUT", REPO / "results" / "dataset_screen"))
RNG = np.random.default_rng(0)
MIN_CPM, MIN_FRAC = 1.0, 0.5
N_PERM = 2000
FNAME = re.compile(r"^GSM\d+_(?:([A-H]\d+)_)?(.+?)_(\d)_rawCounts\.txt\.gz$")


def load_logcpm() -> tuple[pd.DataFrame, pd.DataFrame]:
    cols, meta = {}, []
    for f in sorted(RAW.glob("*_rawCounts.txt.gz")):
        m = FNAME.match(f.name)
        well, cond, rep = m.group(1), m.group(2), int(m.group(3))
        with gzip.open(f, "rt") as fh:
            cols[f"{cond}_{rep}"] = pd.read_csv(fh, sep="\t", index_col=0).iloc[:, 0]
        meta.append({"sample": f"{cond}_{rep}", "condition": cond, "rep": rep,
                     "well": well or "", "col": (well or "  ")[1:]})
    counts = pd.DataFrame(cols)
    md = pd.DataFrame(meta).set_index("sample").loc[counts.columns]
    cpm = counts / counts.sum(axis=0) * 1e6
    keep = (cpm >= MIN_CPM).mean(axis=1) >= MIN_FRAC
    return np.log2(cpm.loc[keep] + 1.0), md


def half_delta(x: pd.DataFrame, alleles: list[str], reps: tuple[int, ...]) -> pd.DataFrame:
    """Allele deltas using only `reps`, referenced to WT restricted to the same reps."""
    base = x[[f"WT_{r}" for r in reps]].mean(axis=1)
    return pd.DataFrame({a: x[[f"{a}_{r}" for r in reps]].mean(axis=1) - base
                         for a in alleles})


def corr(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    a = A - A.mean(axis=0, keepdims=True)
    b = B - B.mean(axis=0, keepdims=True)
    a /= np.linalg.norm(a, axis=0, keepdims=True) + 1e-12
    b /= np.linalg.norm(b, axis=0, keepdims=True) + 1e-12
    return a.T @ b


def pds(C: np.ndarray) -> tuple[float, float]:
    """Mean tie-aware PDS and top-1 accuracy from a square half-A x half-B matrix."""
    n = C.shape[0]
    sc, t1 = [], []
    for i in range(n):
        own = C[i, i]
        better = int((C[i] > own).sum())
        equal = int((C[i] == own).sum()) - 1
        rank = better + equal / 2.0 + 1
        sc.append(1.0 - (rank - 1) / (n - 1))
        t1.append(rank == 1 and equal == 0)
    return float(np.mean(sc)), float(np.mean(t1))


def perm_p(C: np.ndarray, obs: float) -> tuple[float, float]:
    null = np.empty(N_PERM)
    for k in range(N_PERM):
        p = RNG.permutation(C.shape[1])
        null[k] = pds(C[:, p])[0]
    return float(null.mean()), float(((null >= obs).sum() + 1) / (N_PERM + 1))


def main() -> None:
    x, md = load_logcpm()
    alleles = sorted(c for c in md["condition"].unique() if c not in ("WT", "GFP"))
    splits = [((1, 2), (3, 4)), ((1, 3), (2, 4)), ((1, 4), (2, 3))]
    rep: dict = {"alleles": alleles, "n_alleles": len(alleles)}

    # ---------- T1: severity-only baseline ----------
    # identify an allele by matching a single scalar, the L2 norm of its delta
    t1 = []
    for ha, hb in splits:
        A = half_delta(x, alleles, ha)
        B = half_delta(x, alleles, hb)
        na = np.linalg.norm(A.values, axis=0)
        nb = np.linalg.norm(B.values, axis=0)
        # similarity = negative absolute difference in magnitude
        C = -np.abs(na[:, None] - nb[None, :])
        p, t = pds(C)
        t1.append({"half": f"{ha}v{hb}", "PDS": p, "top1": t})
    rep["T1_magnitude_only"] = t1
    rep["T1_PDS_mean"] = float(np.mean([r["PDS"] for r in t1]))
    rep["T1_top1_mean"] = float(np.mean([r["top1"] for r in t1]))

    # ---------- how concentrated is the delta space? ----------
    Afull = half_delta(x, alleles, (1, 2, 3, 4)).values
    Ac = Afull - Afull.mean(axis=1, keepdims=True)
    sv = np.linalg.svd(Ac, compute_uv=False)
    var = sv ** 2 / (sv ** 2).sum()
    rep["pc_variance_explained"] = [round(float(v), 4) for v in var[:6]]

    # ---------- T2: identification after removing top-k shared PCs ----------
    # PCs are estimated on half A only, then applied to both halves, so half B
    # never contributes to the subspace it is scored in.
    t2 = {}
    for k in (0, 1, 2, 3, 5):
        per = []
        for ha, hb in splits:
            A = half_delta(x, alleles, ha).values
            B = half_delta(x, alleles, hb).values
            if k > 0:
                Am = A - A.mean(axis=1, keepdims=True)
                U = np.linalg.svd(Am, full_matrices=False)[0][:, :k]   # gene-space basis
                A = A - U @ (U.T @ A)
                B = B - U @ (U.T @ B)
            C = corr(A, B)
            p, t = pds(C)
            per.append({"half": f"{ha}v{hb}", "PDS": p, "top1": t})
        mp = float(np.mean([r["PDS"] for r in per]))
        mt = float(np.mean([r["top1"] for r in per]))
        # permutation null on the last split's matrix
        nm, pv = perm_p(C, per[-1]["PDS"])
        t2[f"remove_{k}_PC"] = {"per_split": per, "PDS_mean": mp, "top1_mean": mt,
                                "perm_null_mean": nm, "perm_p_lastsplit": pv}
    rep["T2_residual"] = t2

    # ---------- T3: plate-column confound exposure ----------
    known = md[(md["well"] != "") & (~md["condition"].isin(["WT", "GFP"]))]
    col_of = known.groupby("condition")["col"].agg(lambda s: sorted(set(s)))
    rep["T3_allele_to_plate_column"] = {k: v for k, v in col_of.items()}
    rep["T3_n_alleles_with_well_id"] = int(len(col_of))
    rep["T3_n_alleles_without_well_id"] = int(len(alleles) - len(col_of))
    dup = {}
    for a, cs in col_of.items():
        for c in cs:
            dup.setdefault(c, []).append(a)
    rep["T3_columns_shared_by_multiple_alleles"] = {k: v for k, v in dup.items() if len(v) > 1}

    (OUT / "residual_report.json").write_text(json.dumps(rep, indent=2, default=str))

    print("=" * 78)
    print("GSE311877: is PDS=0.96 allele resolution, or severity ranking?")
    print("=" * 78)
    print(f"alleles={len(alleles)}   chance PDS=0.500   chance top1={1/len(alleles):.3f}")
    print("\n[T1] identify allele from ||delta|| alone (one scalar per allele)")
    for r in t1:
        print(f"     {r['half']}: PDS={r['PDS']:.3f}  top1={r['top1']:.3f}")
    print(f"     mean      : PDS={rep['T1_PDS_mean']:.3f}  top1={rep['T1_top1_mean']:.3f}")

    print("\n[  ] variance explained by leading PCs of the allele-delta space")
    print("     " + "  ".join(f"PC{i+1}={v:.3f}" for i, v in enumerate(rep["pc_variance_explained"])))

    print("\n[T2] split-half identification after projecting out top-k shared PCs")
    print(f"     {'k':>3}  {'PDS':>6}  {'top1':>6}  {'perm null':>9}  {'p(last split)':>13}")
    for k, v in t2.items():
        kk = k.split("_")[1]
        print(f"     {kk:>3}  {v['PDS_mean']:>6.3f}  {v['top1_mean']:>6.3f}  "
              f"{v['perm_null_mean']:>9.3f}  {v['perm_p_lastsplit']:>13.4f}")

    print("\n[T3] plate-position confound")
    print(f"     alleles with a well id: {rep['T3_n_alleles_with_well_id']} / {len(alleles)}"
          f"   (no well id: {rep['T3_n_alleles_without_well_id']})")
    print(f"     allele -> plate column: {rep['T3_allele_to_plate_column']}")
    print(f"     columns holding >1 allele (i.e. >=2 plates): "
          f"{rep['T3_columns_shared_by_multiple_alleles']}")
    print("\nwrote residual_report.json")


if __name__ == "__main__":
    main()
