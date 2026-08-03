#!/usr/bin/env python3
"""Eligibility audit of GEO GSE311877 (TP63 arrayed mutants, bulk RNA-seq) as a
candidate second system for AllelePerturb.

Question this script answers, and ONLY this question:
    does the measurement resolve allele identity above chance?
i.e. is the replicate (split-half) ceiling above the permutation null, per allele.

It deliberately does NOT fit or evaluate any predictive model. Feature-based
prediction (ESM / biophysical -> allele residual) is a separate, later step; the
point of this script is to decide whether the dataset is worth that effort.

Design of the data (read off the GEO sample titles, not assumed):
    20 conditions x 4 biological replicates = 80 libraries
    conditions = WT, GFP, and 18 TP63 missense mutants
    5 same-residue sibling pairs: R279Q/S, R304Q/T, L514D/F, C522D/G, L531E/R

Leakage control that matters here: the WT reference must be split alongside the
mutant replicates. Subtracting the SAME WT mean from both halves injects a shared
term into both delta profiles and inflates split-half agreement. Half A uses WT
reps {1,2}; half B uses WT reps {3,4}.

Gene selection is by mean expression across all libraries. That uses no allele
labels, so it cannot leak allele identity into the identification test.
"""
from __future__ import annotations

import gzip
import itertools
import json
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
RAW = Path(os.environ.get("GSE311877_DIR", REPO / "data_external" / "GSE311877"))
OUT = Path(os.environ.get("DATASET_SCREEN_OUT", REPO / "results" / "dataset_screen"))
RNG = np.random.default_rng(0)          # explicit seed
MIN_CPM, MIN_FRAC = 1.0, 0.5            # gene filter: CPM>=1 in >=50% of libraries
N_PERM = 2000

# well-position prefix in some filenames encodes plate layout; strip it for the
# condition label but keep it, because replicate index may be confounded with row.
FNAME = re.compile(r"^GSM\d+_(?:([A-H]\d+)_)?(.+?)_(\d)_rawCounts\.txt\.gz$")


def load() -> tuple[pd.DataFrame, pd.DataFrame]:
    cols, meta = {}, []
    for f in sorted(RAW.glob("*_rawCounts.txt.gz")):
        m = FNAME.match(f.name)
        if not m:
            raise ValueError(f"unparsed sample filename: {f.name}")
        well, cond, rep = m.group(1), m.group(2), int(m.group(3))
        with gzip.open(f, "rt") as fh:
            s = pd.read_csv(fh, sep="\t", index_col=0).iloc[:, 0]
        key = f"{cond}_{rep}"
        if key in cols:
            raise ValueError(f"duplicate condition/replicate key: {key}")
        cols[key] = s
        meta.append({"sample": key, "condition": cond, "rep": rep,
                     "well": well or "", "row": (well or " ")[0], "gsm": f.name.split("_")[0]})
    counts = pd.DataFrame(cols)
    md = pd.DataFrame(meta).set_index("sample").loc[counts.columns]
    return counts, md


def logcpm(counts: pd.DataFrame) -> pd.DataFrame:
    cpm = counts / counts.sum(axis=0) * 1e6
    keep = (cpm >= MIN_CPM).mean(axis=1) >= MIN_FRAC
    return np.log2(cpm.loc[keep] + 1.0)


def delta_profiles(x: pd.DataFrame, md: pd.DataFrame, alleles: list[str],
                   ref: str, mut_reps: list[int], ref_reps: list[int]) -> pd.DataFrame:
    """Mean log2CPM of `mut_reps` for each allele, minus mean of `ref_reps` of `ref`."""
    base = x[[f"{ref}_{r}" for r in ref_reps]].mean(axis=1)
    out = {}
    for a in alleles:
        cs = [f"{a}_{r}" for r in mut_reps if f"{a}_{r}" in x.columns]
        out[a] = x[cs].mean(axis=1) - base
    return pd.DataFrame(out)


def corr_matrix(A: pd.DataFrame, B: pd.DataFrame) -> pd.DataFrame:
    """Pearson correlation of every column of A against every column of B."""
    a = A.values - A.values.mean(axis=0, keepdims=True)
    b = B.values - B.values.mean(axis=0, keepdims=True)
    a /= np.linalg.norm(a, axis=0, keepdims=True)
    b /= np.linalg.norm(b, axis=0, keepdims=True)
    return pd.DataFrame(a.T @ b, index=A.columns, columns=B.columns)


def pds_from_corr(C: pd.DataFrame) -> tuple[float, float, pd.DataFrame]:
    """Tie-aware discrimination score + top-1 accuracy on a half-A x half-B matrix.

    For allele i, rank its own half-B profile among all half-B profiles by
    similarity to half-A_i. PDS = 1 - (rank-1)/(n-1) averaged; 1.0 = perfect,
    chance = 0.5. Ties are resolved by midrank so a degenerate all-equal matrix
    scores exactly chance rather than accidentally winning.
    """
    rows = []
    n = C.shape[1]
    for a in C.index:
        s = C.loc[a]
        own = s[a]
        better = int((s > own).sum())
        equal = int((s == own).sum()) - 1          # exclude self
        rank = better + equal / 2.0 + 1            # midrank, 1 = best
        rows.append({"allele": a, "rank": rank, "self_r": own,
                     "best_other": s.drop(a).max(),
                     "top1": bool(rank == 1 and equal == 0),
                     "pds": 1.0 - (rank - 1) / (n - 1)})
    df = pd.DataFrame(rows).set_index("allele")
    return float(df["pds"].mean()), float(df["top1"].mean()), df


def main() -> None:
    counts, md = load()
    x = logcpm(counts)
    alleles = sorted(c for c in md["condition"].unique() if c not in ("WT", "GFP"))
    rep = {}

    rep["n_libraries"] = int(counts.shape[1])
    rep["n_genes_raw"] = int(counts.shape[0])
    rep["n_genes_kept"] = int(x.shape[0])
    rep["n_alleles"] = len(alleles)
    rep["alleles"] = alleles
    rep["lib_size_min"] = int(counts.sum(axis=0).min())
    rep["lib_size_median"] = int(counts.sum(axis=0).median())
    rep["lib_size_max"] = int(counts.sum(axis=0).max())

    # --- batch confound check: is replicate index aliased to plate row? ---
    rep["rep_vs_row"] = (md[md["well"] != ""]
                         .groupby(["rep", "row"]).size().unstack(fill_value=0).to_dict())
    rep["n_samples_without_well"] = int((md["well"] == "").sum())

    # --- Q: how far is each allele from WT at all (is there any effect)? ---
    d_all = delta_profiles(x, md, alleles + ["GFP"], "WT", [1, 2, 3, 4], [1, 2, 3, 4])
    eff = pd.DataFrame({
        "l2_from_WT": np.linalg.norm(d_all.values, axis=0),
        "n_genes_abs_delta_gt1": (d_all.abs() > 1).sum(axis=0).values,
    }, index=d_all.columns).sort_values("l2_from_WT", ascending=False)

    # --- Q: split-half identification (the measurement ceiling) ---
    # average over all 3 balanced 2v2 splits so the answer is not one lucky split
    splits = [((1, 2), (3, 4)), ((1, 3), (2, 4)), ((1, 4), (2, 3))]
    per_split, mats = [], []
    for ha, hb in splits:
        A = delta_profiles(x, md, alleles, "WT", list(ha), list(ha))
        B = delta_profiles(x, md, alleles, "WT", list(hb), list(hb))
        C = corr_matrix(A, B)
        pds, top1, df = pds_from_corr(C)
        per_split.append({"half_a": ha, "half_b": hb, "PDS": pds, "top1": top1})
        mats.append(C)
        df.to_csv(OUT / f"split_{ha[0]}{ha[1]}_vs_{hb[0]}{hb[1]}.csv")
    rep["split_half"] = per_split
    rep["PDS_mean"] = float(np.mean([s["PDS"] for s in per_split]))
    rep["top1_mean"] = float(np.mean([s["top1"] for s in per_split]))
    rep["PDS_chance"] = 0.5
    rep["top1_chance"] = 1.0 / len(alleles)

    Cbar = sum(m.values for m in mats) / len(mats)
    Cbar = pd.DataFrame(Cbar, index=mats[0].index, columns=mats[0].columns)
    Cbar.to_csv(OUT / "split_half_corr_mean.csv")

    # permutation null on the averaged matrix: shuffle allele labels of half B
    obs, _, _ = pds_from_corr(Cbar)
    null = []
    for _ in range(N_PERM):
        perm = RNG.permutation(Cbar.columns.to_numpy())
        Cp = Cbar.copy()
        Cp.columns = perm
        Cp = Cp[Cbar.columns]
        null.append(pds_from_corr(Cp)[0])
    null = np.asarray(null)
    rep["PDS_observed_meanmatrix"] = float(obs)
    rep["PDS_null_mean"] = float(null.mean())
    rep["PDS_null_p95"] = float(np.percentile(null, 95))
    rep["PDS_perm_p"] = float(((null >= obs).sum() + 1) / (N_PERM + 1))

    # --- Q: D_self vs D_null on single replicates (no averaging) ---
    d1 = {}
    for a in alleles:
        for r in (1, 2, 3, 4):
            k = f"{a}_{r}"
            if k in x.columns:
                d1[k] = x[k] - x[[f"WT_{q}" for q in (1, 2, 3, 4) if q != r]].mean(axis=1)
    D1 = pd.DataFrame(d1)
    R = corr_matrix(D1, D1)
    same, diff = [], []
    for i, j in itertools.combinations(R.columns, 2):
        (same if i.rsplit("_", 1)[0] == j.rsplit("_", 1)[0] else diff).append(R.loc[i, j])
    rep["r_within_allele_mean"] = float(np.mean(same))
    rep["r_between_allele_mean"] = float(np.mean(diff))
    rep["r_within_n"], rep["r_between_n"] = len(same), len(diff)

    # --- Q: can same-residue siblings be told apart? the hardest case ---
    pairs = [("R279Q", "R279S"), ("R304Q", "R304T"), ("L514D", "L514F"),
             ("C522D", "C522G"), ("L531E", "L531R")]
    sib = []
    for p, q in pairs:
        if p in Cbar.index and q in Cbar.index:
            sib.append({
                "pair": f"{p}/{q}",
                "r_self_p": float(Cbar.loc[p, p]), "r_self_q": float(Cbar.loc[q, q]),
                "r_cross": float((Cbar.loc[p, q] + Cbar.loc[q, p]) / 2),
                "p_beats_sibling": bool(Cbar.loc[p, p] > Cbar.loc[p, q]),
                "q_beats_sibling": bool(Cbar.loc[q, q] > Cbar.loc[q, p]),
            })
    rep["sibling_pairs"] = sib

    eff.to_csv(OUT / "effect_size_vs_WT.csv")
    (OUT / "audit_report.json").write_text(json.dumps(rep, indent=2, default=str))

    # ---- console report: numbers only ----
    print("=" * 74)
    print("GSE311877 TP63 arrayed mutants: eligibility audit")
    print("=" * 74)
    print(f"libraries={rep['n_libraries']}  genes_raw={rep['n_genes_raw']}  "
          f"genes_kept={rep['n_genes_kept']}  alleles={rep['n_alleles']}")
    print(f"lib size  min={rep['lib_size_min']:,}  med={rep['lib_size_median']:,}  "
          f"max={rep['lib_size_max']:,}")
    print(f"samples with no well id: {rep['n_samples_without_well']}")
    print("\n-- effect size vs WT (log2CPM L2 norm, genes with |delta|>1) --")
    print(eff.to_string())
    print("\n-- split-half identification (chance PDS=0.500, "
          f"chance top1={rep['top1_chance']:.3f}) --")
    for s in per_split:
        print(f"  half {s['half_a']} vs {s['half_b']}: PDS={s['PDS']:.3f}  top1={s['top1']:.3f}")
    print(f"  mean over splits : PDS={rep['PDS_mean']:.3f}  top1={rep['top1_mean']:.3f}")
    print(f"  averaged matrix  : PDS={obs:.3f}  perm null mean={null.mean():.3f}  "
          f"p95={np.percentile(null,95):.3f}  p={rep['PDS_perm_p']:.4f}")
    print("\n-- single-replicate delta correlations --")
    print(f"  within allele  r={rep['r_within_allele_mean']:.3f}  (n={rep['r_within_n']})")
    print(f"  between allele r={rep['r_between_allele_mean']:.3f}  (n={rep['r_between_n']})")
    print("\n-- same-residue sibling discrimination --")
    print(pd.DataFrame(sib).to_string(index=False))
    print("\nwrote audit_report.json, split_half_corr_mean.csv, effect_size_vs_WT.csv")


if __name__ == "__main__":
    main()
