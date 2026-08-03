#!/usr/bin/env python3
"""Direct test of the plate-position confound in GSE311877.

residual_gse311877.py established that allele identity is aliased to (plate,
column): every replicate of an allele sits in one column, replicates differ only
by row. So "split-half identification works" is consistent with two stories:

    story A  the transcriptome genuinely resolves the allele
    story B  the transcriptome resolves the well position

Story B makes one falsifiable prediction the data can test. Six columns each
host two different alleles (on two different plates). If a column-specific
artifact were carrying the signal, those same-column allele pairs should be
unusually similar to each other. If their similarity is no different from
arbitrary allele pairs, the column-artifact explanation loses its main support.

This does not fully exonerate the design: nothing here can separate a
plate-level effect from a biological one, because no allele is replicated
across plates. That limit is reported, not papered over.
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
N_PERM = 20000
FNAME = re.compile(r"^GSM(\d+)_(?:([A-H])(\d+)_)?(.+?)_(\d)_rawCounts\.txt\.gz$")


def load():
    cols, meta = {}, []
    for f in sorted(RAW.glob("*_rawCounts.txt.gz")):
        m = FNAME.match(f.name)
        gsm, row, col, cond, rep = m.groups()
        with gzip.open(f, "rt") as fh:
            cols[f"{cond}_{rep}"] = pd.read_csv(fh, sep="\t", index_col=0).iloc[:, 0]
        meta.append({"sample": f"{cond}_{rep}", "gsm": int(gsm), "condition": cond,
                     "rep": int(rep), "row": row or "", "col": col or ""})
    counts = pd.DataFrame(cols)
    md = pd.DataFrame(meta).set_index("sample").loc[counts.columns]
    cpm = counts / counts.sum(axis=0) * 1e6
    keep = (cpm >= MIN_CPM).mean(axis=1) >= MIN_FRAC
    return np.log2(cpm.loc[keep] + 1.0), md


def deltas(x, alleles, reps):
    base = x[[f"WT_{r}" for r in reps]].mean(axis=1)
    return pd.DataFrame({a: x[[f"{a}_{r}" for r in reps]].mean(axis=1) - base
                         for a in alleles})


def corr_df(A, B):
    a = A.values - A.values.mean(axis=0, keepdims=True)
    b = B.values - B.values.mean(axis=0, keepdims=True)
    a /= np.linalg.norm(a, axis=0, keepdims=True)
    b /= np.linalg.norm(b, axis=0, keepdims=True)
    return pd.DataFrame(a.T @ b, index=A.columns, columns=B.columns)


def main() -> None:
    x, md = load()
    alleles = sorted(c for c in md["condition"].unique() if c not in ("WT", "GFP"))

    # ---- recover the plate layout actually encoded in the filenames ----
    lay = (md[md["condition"].isin(alleles)]
           .groupby("condition")
           .agg(gsm_min=("gsm", "min"), rows=("row", lambda s: "".join(sorted(set(s)))),
                col=("col", lambda s: "/".join(sorted(set(s))))))
    lay["block"] = np.where(lay["rows"].str.startswith("A"), "rows_A-D",
                            np.where(lay["rows"].str.startswith("E"), "rows_E-H", "no_well_id"))
    lay = lay.sort_values("gsm_min")

    rep = {"layout": lay.to_dict(orient="index")}

    # ---- full-depth allele deltas and their correlation matrix ----
    D = deltas(x, alleles, (1, 2, 3, 4))
    C = corr_df(D, D)

    # residual space: strip the shared severity axes that dominate raw similarity
    M = D.values - D.values.mean(axis=1, keepdims=True)
    U = np.linalg.svd(M, full_matrices=False)[0][:, :5]
    Dres = pd.DataFrame(D.values - U @ (U.T @ D.values), index=D.index, columns=D.columns)
    Cres = corr_df(Dres, Dres)

    # ---- same-column allele pairs vs everything else ----
    have = lay[lay["col"] != ""]
    bycol: dict[str, list[str]] = {}
    for a, r in have.iterrows():
        bycol.setdefault(r["col"], []).append(a)
    same_col_pairs = [(p, q) for v in bycol.values() if len(v) == 2 for p, q in [tuple(sorted(v))]]

    def contrast(Cm: pd.DataFrame, label: str) -> dict:
        pool = have.index.tolist()
        sc = [float(Cm.loc[p, q]) for p, q in same_col_pairs]
        other = [float(Cm.loc[p, q]) for i, p in enumerate(pool) for q in pool[i + 1:]
                 if (min(p, q), max(p, q)) not in {(min(a, b), max(a, b)) for a, b in same_col_pairs}]
        obs = float(np.mean(sc)) - float(np.mean(other))
        # permutation: reassign which allele sits in which column, keep the pairing shape
        null = np.empty(N_PERM)
        arr = np.array(pool)
        idx = {a: i for i, a in enumerate(pool)}
        Cv = Cm.loc[pool, pool].values
        k = len(same_col_pairs)
        for t in range(N_PERM):
            perm = RNG.permutation(len(arr))
            pairs = [(perm[2 * i], perm[2 * i + 1]) for i in range(k)]
            m_in = np.mean([Cv[i, j] for i, j in pairs])
            mask = np.ones((len(arr), len(arr)), bool)
            np.fill_diagonal(mask, False)
            for i, j in pairs:
                mask[i, j] = mask[j, i] = False
                                                    # remaining upper-triangle entries
            m_out = np.mean(Cv[np.triu(mask, 1)])
            null[t] = m_in - m_out
        p = float(((np.abs(null) >= abs(obs)).sum() + 1) / (N_PERM + 1))
        return {"space": label, "n_same_col_pairs": len(sc),
                "mean_r_same_column": float(np.mean(sc)),
                "mean_r_other_pairs": float(np.mean(other)),
                "difference": obs, "perm_null_sd": float(null.std()), "perm_p_two_sided": p,
                "per_pair": {f"{p_}/{q_}": float(Cm.loc[p_, q_]) for p_, q_ in same_col_pairs}}

    rep["raw_space"] = contrast(C, "raw delta")
    rep["residual_space"] = contrast(Cres, "residual (top-5 PC removed)")

    # ---- block (plate) effect: rows A-D vs rows E-H vs no-well batch ----
    blocks = lay.groupby("block").apply(lambda g: list(g.index), include_groups=False).to_dict()
    bl = {}
    for name, mem in blocks.items():
        within, across = [], []
        for i, p in enumerate(alleles):
            for q in alleles[i + 1:]:
                if p in mem and q in mem:
                    within.append(float(Cres.loc[p, q]))
                elif (p in mem) != (q in mem):
                    across.append(float(Cres.loc[p, q]))
        bl[name] = {"n_alleles": len(mem), "members": mem,
                    "mean_resid_r_within_block": float(np.mean(within)) if within else None,
                    "mean_resid_r_across_block": float(np.mean(across)) if across else None}
    rep["block_effect_residual_space"] = bl

    (OUT / "confound_report.json").write_text(json.dumps(rep, indent=2, default=str))

    print("=" * 80)
    print("GSE311877 plate-position confound: is the signal allele, or well position?")
    print("=" * 80)
    print("\n-- layout recovered from GEO filenames --")
    print(lay.to_string())

    for key in ("raw_space", "residual_space"):
        r = rep[key]
        print(f"\n-- same-column allele pairs, {r['space']} --")
        print(f"   pairs sharing a column (n={r['n_same_col_pairs']}): "
              f"mean r = {r['mean_r_same_column']:.4f}")
        print(f"   all other allele pairs           : mean r = {r['mean_r_other_pairs']:.4f}")
        print(f"   difference = {r['difference']:+.4f}   perm null sd = {r['perm_null_sd']:.4f}"
              f"   two-sided p = {r['perm_p_two_sided']:.4f}")
        for k, v in r["per_pair"].items():
            print(f"      {k:<16} r = {v:.4f}")

    print("\n-- block (plate) structure in residual space --")
    for name, v in bl.items():
        print(f"   {name:<12} n={v['n_alleles']:<3} within r="
              f"{v['mean_resid_r_within_block'] if v['mean_resid_r_within_block'] is None else round(v['mean_resid_r_within_block'],4)}"
              f"   across r="
              f"{v['mean_resid_r_across_block'] if v['mean_resid_r_across_block'] is None else round(v['mean_resid_r_across_block'],4)}")
        print(f"                {v['members']}")
    print("\nwrote confound_report.json")


if __name__ == "__main__":
    main()
