#!/usr/bin/env python3
"""Reviewer control (statistical independence of variants).

Two checks that the direction-vs-ranking dissociation is not an artifact of
treating correlated substitutions at the same residue as independent units:

  1. residue-leakage: fraction of held-out variants that share a residue
     position with a training variant of the same gene, per split.
  2. residue-cluster bootstrap: recompute the per-method PDS / Pearson-delta
     95% CIs resampling residue-level clusters instead of individual variants,
     and check whether every non-null predictor still overlaps 0.50 (PDS) and
     stays above 0 (Pearson-delta).

Inputs (committed): data/allele_perturb_bench.csv, results/canonical/unified_results5.csv
Output: results/reviewer_controls/residue_independence_summary.csv
"""
import argparse
import re
from pathlib import Path

_ap = argparse.ArgumentParser(description=__doc__,
                              formatter_class=argparse.RawDescriptionHelpFormatter)
_ap.add_argument("--out", required=True, type=Path,
                 help="directory receiving residue_independence_summary.csv and "
                      "residue_leakage.csv. Required and never defaulted, so a re-run "
                      "cannot overwrite the committed tables under results/")
_args = _ap.parse_args()

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
OUT_DIR = _args.out.expanduser().resolve()
OUT_DIR.mkdir(parents=True, exist_ok=True)
bench = pd.read_csv(REPO / "data/allele_perturb_bench.csv")
res = pd.read_csv(REPO / "results/canonical/unified_results5.csv")
NB = 2000
SEED = 0


def residues(v):
    return set(int(x) for x in re.findall(r"(\d+)", str(v)))


def residue_leakage():
    out = []
    for sp in ["split1_role", "split2_role", "split3_role", "split5_role", "split6_role"]:
        tot = leak = 0
        for _, sub in bench.groupby("gene"):
            tr, te = sub[sub[sp] == "train"], sub[sub[sp] == "test"]
            if len(te) == 0:
                continue
            trres = set().union(*[residues(v) for v in tr["variant"]]) if len(tr) else set()
            for v in te["variant"]:
                tot += 1
                if residues(v) & trres:
                    leak += 1
        out.append(dict(split=sp.replace("_role", ""), n_test=tot,
                        residue_leak_frac=round(leak / max(tot, 1), 3)))
    return pd.DataFrame(out)


def cluster_bootstrap(split="split1"):
    r = res[res.split == split].copy()
    r["res"] = r.apply(lambda x: (x.gene, min(residues(x.variant)) if residues(x.variant) else -1), axis=1)
    rng = np.random.RandomState(SEED)

    def ci(df, col, by_cluster):
        groups = [g[col].values for _, g in df.groupby("res")] if by_cluster else [np.array([v]) for v in df[col].values]
        n = len(groups)
        means = [np.concatenate([groups[i] for i in rng.randint(0, n, n)]).mean() for _ in range(NB)]
        return np.percentile(means, [2.5, 97.5])

    rows = []
    for m, dfm in r.groupby("method"):
        if "null" in m.lower():
            continue
        pc = ci(dfm, "PDS_cos", True)
        ec = ci(dfm, "pearson_delta", True)
        rows.append(dict(method=m, split=split,
                         pds_mean=round(dfm.PDS_cos.mean(), 3),
                         pds_clust_lo=round(pc[0], 3), pds_clust_hi=round(pc[1], 3),
                         pds_overlaps_0p5=bool(pc[0] <= 0.5 <= pc[1]),
                         pearson_mean=round(dfm.pearson_delta.mean(), 3),
                         pearson_clust_lo=round(ec[0], 3),
                         pearson_excludes_0=bool(ec[0] > 0)))
    return pd.DataFrame(rows)


if __name__ == "__main__":
    leak = residue_leakage()
    clust = cluster_bootstrap("split1")
    print("=== residue-leakage per split ===")
    print(leak.to_string(index=False))
    n = len(clust)
    n_pds_ov = int(clust.pds_overlaps_0p5.sum())
    n_pe_ex = int(clust.pearson_excludes_0.sum())
    print(f"\n=== residue-cluster bootstrap (split1), {n} non-null methods ===")
    print(f"PDS CI overlaps 0.50: {n_pds_ov}/{n}  (exceptions: {list(clust[~clust.pds_overlaps_0p5].method)})")
    print(f"Pearson-delta CI excludes 0: {n_pe_ex}/{n}  (exceptions: {list(clust[~clust.pearson_excludes_0].method)})")
    clust.to_csv(OUT_DIR / "residue_independence_summary.csv", index=False)
    leak.to_csv(OUT_DIR / "residue_leakage.csv", index=False)
    print("\nsaved residue_independence_summary.csv + residue_leakage.csv")
