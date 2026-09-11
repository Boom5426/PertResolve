#!/usr/bin/env python3
"""Does the split-half reproducibility reference change when the two halves come
from disjoint experimental batches?

Frozen protocol: docs/PREREG_REPLICATE_REFERENCE_v1.md. Read it before reading any
number this script prints.

The manuscript's reference (scripts/analysis/oracle_ceiling.py) shuffles all of a
perturbation's cells and splits them in two. Both halves therefore carry near
identical batch composition, so any batch term largely cancels. This script keeps
that construction as arm P and adds two arms that do not cancel it:

    P  pooled          two disjoint cell sets drawn ignoring batch  (the published quantity)
    W  within-batch    two disjoint cell sets from the SAME batch
    X  cross-batch     one cell set per disjoint batch group

Every arm uses the same matched cell depth, the same candidate pool, the same
control-splitting rule and the same scorer, so the arms differ only in where the
cells come from. The control population is split by the arm's own rule as well, so
no control cell is ever shared between a prediction and its target.

No model is fitted anywhere here.
"""
from __future__ import annotations

import argparse
import gzip
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from pertresolve.paths import (ATLAS_ENV_VAR, add_harness_to_path, reject_repo_results,
                                 require_inputs, resolve_base)
from pertresolve.resolution.scaling import tie_aware_pds

MASTER_SEED = 20260825
ALLELE = ("TP53", "KRAS")
ATLAS = ("Norman", "sciPlex", "VCC", "Replogle", "GATA1")
WT_TAGS = ("WT", "wt", "WT_control")

# matches results/benchmark_resolution/benchmark_resolution.py
ATLAS_CFG = {
    "Norman":   dict(fname="NormanWeissman2019_filtered.h5ad",            pcol="perturbation", bcol="gemgroup",  ctrl=None),
    "sciPlex":  dict(fname="SrivatsanTrapnell2020_sciplex3.h5ad",         pcol="perturbation", bcol="replicate", ctrl="control"),
    "VCC":      dict(fname="adata_Training.h5ad",                          pcol="target_gene",  bcol="batch",     ctrl="non-targeting"),
    "Replogle": dict(fname="ReplogleWeissman2022_K562_essential.h5ad",     pcol="perturbation", bcol="batch",     ctrl=None),
    # GATA1's replicate label was never stored as its own column: the source object
    # prefixes every cell barcode with its replicate, e.g. "rep25_AAACCC...-1".
    "GATA1":    dict(fname="GATA1_standard_hvg_pert_filtered.h5ad",         pcol="mutation_name",
                     bcol="cell_barcodes", ctrl="WT", bsplit="_"),
}
CCAND = ("control", "Control", "CTRL", "ctrl", "non-targeting", "NT", "DMSO", "unperturbed", "None")
MAXPERT = 400          # perturbation cap, as in benchmark_resolution.py
MAXCELL_PER_GROUP = 400  # cells kept per (perturbation, batch group); >= 2 * max depth
MAXCTRL_PER_GROUP = 2000


# --------------------------------------------------------------------------- data


def load_allele(gene: str, base: Path):
    """Cells, variant labels and recovered 10x channel for TP53 or KRAS.

    The canonical arrays dropped the batch column at preprocessing. It is rebuilt
    from the immutable GEO source and the reconstruction is asserted against the
    array's own label sequence, so a silent row-order change cannot go unnoticed.
    """
    add_harness_to_path(base)
    import harness as H
    H.set_base(base)
    X, lab = H.load_gene(gene)
    lab = np.asarray(lab).astype(str)

    v2c_path = base / "raw" / f"GSE161824_A549_{gene}.variants2cell.csv.gz"
    mtx_path = base / "raw" / f"GSE161824_A549_{gene}.processed.matrix.mtx.gz"
    require_inputs(v2c_path, mtx_path)
    v2c = pd.read_csv(v2c_path, sep="\t", usecols=["batch", "variant"])
    with gzip.open(mtx_path, "rt") as fh:
        for line in fh:
            if line.startswith("%"):
                continue
            n_cell = int(line.split()[0])
            break
    pre = v2c.iloc[:n_cell]
    keep = ~pre["variant"].isin(["unassigned", "multiple"]).to_numpy()
    rebuilt = pre["variant"].to_numpy().astype(str)[keep]
    batch = pre["batch"].to_numpy().astype(str)[keep]
    if len(rebuilt) != len(lab) or not (rebuilt == lab).all():
        raise RuntimeError(
            f"{gene}: batch reconstruction does not line up with the canonical array "
            f"({len(rebuilt)} rebuilt vs {len(lab)} array rows). Refusing to proceed.")
    ctrl_mask = np.isin(lab, WT_TAGS)
    return X, lab, batch, ctrl_mask


def load_atlas(name: str, atlas_dir: Path, max_depth: int, subset=None, batch_prefix=0):
    """Cells in PCA-50 space, perturbation labels and batch, for one public screen.

    Representation follows results/benchmark_resolution/benchmark_resolution.py: 50
    components fitted without perturbation labels, TruncatedSVD when the gene space
    is too large to densify.
    """
    import anndata as ad
    import scipy.sparse as sp
    from sklearn.decomposition import PCA, TruncatedSVD

    cfg = ATLAS_CFG[name]
    path = atlas_dir / cfg["fname"]
    require_inputs(path)
    A = ad.read_h5ad(path, backed="r")
    obs = A.obs
    pert = obs[cfg["pcol"]].astype(str).to_numpy()
    batch = obs[cfg["bcol"]].astype(str).to_numpy()
    if cfg.get("bsplit"):
        batch = np.array([b.split(cfg["bsplit"], 1)[0] for b in batch])
    if batch_prefix:
        # Coarsen the batch label to its leading characters. VCC's 48 batches are
        # Flex_<run>_<channel>, so a prefix of 6 groups channels into the three
        # genuine sequencing runs, which is a higher rung than a channel split.
        batch = np.array([b[:batch_prefix] for b in batch])
    keep_mask = np.ones(len(pert), bool)
    if subset:
        col, val = subset
        if col not in obs.columns:
            A.file.close()
            raise RuntimeError(f"{name}: no obs column {col!r} to subset on")
        keep_mask = obs[col].astype(str).to_numpy() == val
        # a perturbation absent from the subset simply drops out via the count filter
        pert = np.where(keep_mask, pert, "nan")
    ctrl = cfg["ctrl"]
    if ctrl is None:
        uniq = set(np.unique(pert))
        ctrl = next((c for c in CCAND if c in uniq), None)
    if ctrl is None:
        A.file.close()
        raise RuntimeError(f"{name}: no control label found among {CCAND}")

    rng = np.random.RandomState(MASTER_SEED % (2 ** 31))
    groups = batch_groups(np.unique(batch))
    gid = np.where(np.isin(batch, groups[0]), 0, 1)

    counts = pd.Series(pert).value_counts()
    perts = [p for p in counts.index if p not in (ctrl, "nan") and counts[p] >= 2 * max_depth]
    if len(perts) > MAXPERT:
        perts = sorted(rng.choice(perts, MAXPERT, replace=False))

    rows = []
    for p in perts:
        for g in (0, 1):
            idx = np.where((pert == p) & (gid == g))[0]
            if len(idx) > MAXCELL_PER_GROUP:
                idx = rng.choice(idx, MAXCELL_PER_GROUP, replace=False)
            rows.append(idx)
    for g in (0, 1):
        idx = np.where((pert == ctrl) & (gid == g))[0]
        if len(idx) > MAXCTRL_PER_GROUP:
            idx = rng.choice(idx, MAXCTRL_PER_GROUP, replace=False)
        rows.append(idx)
    sel = np.unique(np.concatenate(rows))

    Xraw = A[sel].to_memory().X
    A.file.close()
    if Xraw.shape[1] > 40000:
        Xs = Xraw.tocsr() if sp.issparse(Xraw) else sp.csr_matrix(np.asarray(Xraw))
        Xp = TruncatedSVD(n_components=50, random_state=0).fit_transform(Xs).astype(np.float32)
    else:
        Xd = np.asarray(Xraw.todense() if hasattr(Xraw, "todense") else Xraw, dtype=np.float32)
        fit = rng.choice(len(Xd), min(20000, len(Xd)), replace=False)
        Xp = PCA(n_components=min(50, Xd.shape[1] - 1), random_state=0).fit(Xd[fit]).transform(Xd)
    lab = pert[sel]
    ctrl_mask = lab == ctrl
    return Xp.astype(np.float32), lab.astype(str), batch[sel].astype(str), ctrl_mask


def batch_groups(levels) -> tuple[list, list]:
    """Split batch levels into two disjoint groups, alternating over sorted labels.

    Uses batch identity only. No perturbation label and no expression value enters
    the assignment.
    """
    lv = sorted(levels)
    return [lv[i] for i in range(0, len(lv), 2)], [lv[i] for i in range(1, len(lv), 2)]


# ------------------------------------------------------------------------ scoring


def eligible(pert_idx, gid, batch, ctrl_counts, depth: int, w_batch) -> dict:
    """Which arms this perturbation can support at `depth` cells per side."""
    n_a = int((gid[pert_idx] == 0).sum())
    n_b = int((gid[pert_idx] == 1).sum())
    n_w = 0 if w_batch is None else int((batch[pert_idx] == w_batch).sum())
    return {
        "P": len(pert_idx) >= 2 * depth and ctrl_counts["pooled"] >= 2 * depth,
        "X": n_a >= depth and n_b >= depth and ctrl_counts[0] >= depth and ctrl_counts[1] >= depth,
        "W": w_batch is not None and n_w >= 2 * depth
             and ctrl_counts.get(("batch", w_batch), 0) >= 2 * depth,
    }


def draw_control(arm, rng, ctrl_idx, gid_c, batch_c, depth, w_batch):
    """One control split per (seed, arm), shared by every perturbation.

    This is not a detail. The published harness draws the wild-type halves once per
    seed (oracle_ceiling.halves_delta) and subtracts the same pair from every
    variant. Redrawing per perturbation makes each prediction and its own target use
    complementary halves of the same finite control pool, which are anti-correlated,
    while cross-pairs are not; the correct target then ranks last and PDS collapses
    far below chance. Observed directly: TP53 arm P scored 0.041 at depth 50 under a
    per-perturbation control draw, against a published 0.485.
    """
    if arm == "X":
        ca = rng.permutation(ctrl_idx[gid_c == 0])[:depth]
        cb = rng.permutation(ctrl_idx[gid_c == 1])[:depth]
        return ca, cb
    pool = ctrl_idx if arm == "P" else ctrl_idx[batch_c == w_batch]
    c = rng.permutation(pool)
    return c[:depth], c[depth:2 * depth]


def draw_pert(arm, rng, pert_idx, gid_p, batch_p, depth, w_batch):
    """Two disjoint perturbation-cell blocks, from one shuffle per side."""
    if arm == "X":
        pa = rng.permutation(pert_idx[gid_p == 0])[:depth]
        pb = rng.permutation(pert_idx[gid_p == 1])[:depth]
        return pa, pb
    pool = pert_idx if arm == "P" else pert_idx[batch_p == w_batch]
    p = rng.permutation(pool)
    return p[:depth], p[depth:2 * depth]


def run_dataset(X, lab, batch, ctrl_mask, depths, n_seed, name):
    """Score every arm at every depth, per perturbation, seed-averaged.

    Arms P and X share one cohort (the intersection of what both can support) and one
    candidate pool, so their ranks are comparable. Arm W is scored on its own cohort
    against its own pool and is reported separately: it is pinned to a single batch,
    so demanding that every P/X perturbation also clear it would shrink the primary
    cohort for the sake of a secondary arm.
    """
    groups = batch_groups(np.unique(batch))
    gid_all = np.where(np.isin(batch, groups[0]), 0, 1)
    ctrl_idx = np.where(ctrl_mask)[0]
    perts = sorted(set(lab[~ctrl_mask]) - {"nan"})

    ctrl_counts = {"pooled": len(ctrl_idx), 0: int((gid_all[ctrl_idx] == 0).sum()),
                   1: int((gid_all[ctrl_idx] == 1).sum())}
    ctrl_per_batch = pd.Series(batch[ctrl_idx]).value_counts()
    for b, c in ctrl_per_batch.items():
        ctrl_counts[("batch", b)] = int(c)
    # arm W is pinned to the single batch with the most control cells, fixed for the
    # whole run so that its candidate pool does not move between seeds.
    w_batch = str(ctrl_per_batch.index[0]) if len(ctrl_per_batch) else None

    rows = []
    meta = {"n_batches": int(len(np.unique(batch))),
            "group_sizes": [len(groups[0]), len(groups[1])],
            "n_control_cells": int(len(ctrl_idx)),
            "control_cells_per_group": [ctrl_counts[0], ctrl_counts[1]],
            "w_batch": w_batch,
            "control_cells_in_w_batch": int(ctrl_counts.get(("batch", w_batch), 0)) if w_batch else 0,
            "cohorts": {}}

    pert_idx = {p: np.where(lab == p)[0] for p in perts}
    for depth in depths:
        el = {p: eligible(pert_idx[p], gid_all, batch, ctrl_counts, depth, w_batch) for p in perts}
        cohort = {
            "PX": sorted(p for p in perts if el[p]["P"] and el[p]["X"]),
            "W": sorted(p for p in perts if el[p]["W"]),
        }
        meta["cohorts"][str(depth)] = {
            "n_PX": len(cohort["PX"]), "n_W": len(cohort["W"]),
            "n_excluded_from_PX": len(perts) - len(cohort["PX"]),
            "excluded_from_PX": sorted(set(perts) - set(cohort["PX"]))[:50],
        }
        for arms, key in ((("P", "X"), "PX"), (("W",), "W")):
            scored = cohort[key]
            if len(scored) < 3:
                continue
            idx_of = {p: i for i, p in enumerate(scored)}
            acc = {a: np.zeros((len(scored), n_seed)) for a in arms}
            for s_i in range(n_seed):
                seed = int(np.random.SeedSequence(
                    [MASTER_SEED, depth, s_i, sum(ord(c) for c in name)]).generate_state(1)[0] % (2 ** 31))
                for a in arms:
                    rng = np.random.RandomState(seed)
                    ca, cb = draw_control(a, rng, ctrl_idx, gid_all[ctrl_idx],
                                          batch[ctrl_idx], depth, w_batch)
                    base_a, base_b = X[ca].mean(0), X[cb].mean(0)
                    pred = np.empty((len(scored), X.shape[1]), dtype=np.float64)
                    truth = np.empty_like(pred)
                    for p in scored:
                        pi = pert_idx[p]
                        pa, pb = draw_pert(a, rng, pi, gid_all[pi], batch[pi], depth, w_batch)
                        pred[idx_of[p]] = X[pa].mean(0) - base_a
                        truth[idx_of[p]] = X[pb].mean(0) - base_b
                    acc[a][:, s_i] = per_pert_pds(pred, truth)
            for a in arms:
                for p in scored:
                    rows.append({"dataset": name, "depth": depth, "arm": a, "cohort": key,
                                 "perturbation": p, "n_pool": len(scored), "n_seed": n_seed,
                                 "pds": float(acc[a][idx_of[p]].mean())})
    return pd.DataFrame(rows), meta


def per_pert_pds(pred: np.ndarray, truth: np.ndarray) -> np.ndarray:
    """Tie-aware cosine PDS for every row, matching pertresolve tie_aware_pds."""
    qn = pred / (np.linalg.norm(pred, axis=1, keepdims=True) + 1e-12)
    tn = truth / (np.linalg.norm(truth, axis=1, keepdims=True) + 1e-12)
    dist = 1.0 - qn @ tn.T
    n = dist.shape[1]
    if n < 2:
        return np.full(dist.shape[0], 0.5)
    out = np.empty(dist.shape[0])
    for i in range(dist.shape[0]):
        row = dist[i]
        target = row[i]
        less = int((row < target - 1e-12).sum())
        eq = int((np.abs(row - target) <= 1e-12).sum())
        out[i] = 1.0 - (less + (eq - 1) / 2.0) / (n - 1)
    return out


# ------------------------------------------------------------------------ summary


def summarise(df: pd.DataFrame, n_boot: int, n_perm: int, tau: float) -> pd.DataFrame:
    rng = np.random.default_rng(MASTER_SEED + 2)
    out = []
    for (ds, depth), sub in df.groupby(["dataset", "depth"]):
        wide = sub[sub["cohort"] == "PX"].pivot(index="perturbation", columns="arm", values="pds")
        w_arm = sub[sub["cohort"] == "W"]
        n = len(wide)
        if n < 3:
            continue
        bi = rng.integers(0, n, size=(n_boot, n))
        row = {"dataset": ds, "depth": int(depth), "n_perturbations": n,
               "n_pool": int(sub["n_pool"].iloc[0]), "tau": tau}
        # PDS is a rank against a pool, so a mean computed against a pool of 15 is not
        # on the same scale as one against a pool of 195. Report arm W only when its
        # cohort is the PX cohort; otherwise record why it is missing.
        row["n_perturbations_W"] = int(w_arm["perturbation"].nunique()) if not w_arm.empty else 0
        if not w_arm.empty and set(w_arm["perturbation"]) == set(wide.index):
            row["pds_W"] = float(w_arm["pds"].mean())
            row["W_comparable"] = True
        else:
            row["pds_W"] = np.nan
            row["W_comparable"] = False
        for a in wide.columns:
            v = wide[a].to_numpy()
            bs = v[bi].mean(1)
            row[f"pds_{a}"] = float(v.mean())
            row[f"pds_{a}_lo"] = float(np.percentile(bs, 2.5))
            row[f"pds_{a}_hi"] = float(np.percentile(bs, 97.5))
            row[f"frac_{a}_above_065"] = float((v > 0.65).mean())
        if "P" in wide.columns and "X" in wide.columns:
            d = (wide["P"] - wide["X"]).to_numpy()
            bd = d[bi].mean(1)
            sgn = rng.choice([-1.0, 1.0], size=(n_perm, n))
            null = (sgn * d).mean(1)
            row.update({
                "delta_PX": float(d.mean()),
                "delta_PX_lo": float(np.percentile(bd, 2.5)),
                "delta_PX_hi": float(np.percentile(bd, 97.5)),
                "delta_PX_signflip_p": float(((np.abs(null) >= abs(d.mean())).sum() + 1) / (n_perm + 1)),
                "frac_perturbations_delta_gt0": float((d > 0).mean()),
                "verdict_flip_frac": float(((wide["P"] > 0.65) != (wide["X"] > 0.65)).mean()),
            })
            row["verdict"] = verdict(row, tau)
        out.append(row)
    return pd.DataFrame(out)


def verdict(row: dict, tau: float) -> str:
    """The binding rule from docs/PREREG_REPLICATE_REFERENCE_v1.md section 7."""
    if row["pds_P_lo"] <= 0.5 <= row["pds_P_hi"]:
        return "NOT_INFORMATIVE"
    lo, hi = row["delta_PX_lo"], row["delta_PX_hi"]
    same_side = (row["pds_P"] > 0.65) == (row["pds_X"] > 0.65)
    if -tau < lo and hi < tau and same_side:
        return "CONCORDANT"
    if lo > tau or (not same_side and row["pds_P"] > 0.65):
        return "OVERSTATING"
    if hi < -tau:
        return "UNDERSTATING"
    return "INCONCLUSIVE"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("dataset", choices=ALLELE + ATLAS)
    ap.add_argument("--out", required=True)
    ap.add_argument("--base", default=None, help="workspace holding the per-gene arrays")
    ap.add_argument("--atlas-dir", default=None, help=f"h5ad directory; falls back to ${ATLAS_ENV_VAR}")
    ap.add_argument("--depths", default="30,50,100")
    ap.add_argument("--n-seed", type=int, default=15)
    ap.add_argument("--n-boot", type=int, default=10000)
    ap.add_argument("--n-perm", type=int, default=10000)
    ap.add_argument("--tau", type=float, default=0.05)
    ap.add_argument("--subset", default=None, metavar="COL=VALUE",
                    help="restrict an atlas to one level of an obs column, e.g. cell_line=A549. "
                         "Used to check that a cross-batch difference is replicate variance "
                         "rather than composition drift between batches.")
    ap.add_argument("--tag", default="", help="suffix appended to output filenames")
    ap.add_argument("--batch-prefix", type=int, default=0, metavar="N",
                    help="coarsen batch labels to their first N characters before grouping, "
                         "e.g. 6 turns VCC's Flex_1_01..Flex_3_16 into three sequencing runs")
    a = ap.parse_args()

    out = reject_repo_results(a.out)
    out.mkdir(parents=True, exist_ok=True)
    depths = [int(d) for d in a.depths.split(",")]

    if a.dataset in ALLELE:
        base = resolve_base(a.base)
        X, lab, batch, ctrl = load_allele(a.dataset, base)
    else:
        atlas = resolve_base(a.atlas_dir, what="the perturbation atlas directory",
                             env_var=ATLAS_ENV_VAR, flag="--atlas-dir")
        sub = tuple(a.subset.split("=", 1)) if a.subset else None
        X, lab, batch, ctrl = load_atlas(a.dataset, atlas, max(depths), sub, a.batch_prefix)

    df, meta = run_dataset(X, lab, batch, ctrl, depths, a.n_seed, a.dataset)
    if df.empty:
        raise RuntimeError(f"{a.dataset}: no depth had enough eligible perturbations; see meta")
    summary = summarise(df, a.n_boot, a.n_perm, a.tau)

    df.to_csv(out / f"replicate_reference_perpert_{a.dataset}{a.tag}.csv.gz", index=False)
    summary.to_csv(out / f"replicate_reference_summary_{a.dataset}{a.tag}.csv", index=False)
    meta.update({"dataset": a.dataset, "subset": a.subset, "batch_prefix": a.batch_prefix, "n_cells": int(X.shape[0]), "n_features": int(X.shape[1]),
                 "master_seed": MASTER_SEED, "n_seed": a.n_seed, "depths": depths, "tau": a.tau})
    (out / f"replicate_reference_meta_{a.dataset}{a.tag}.json").write_text(json.dumps(meta, indent=2, default=str))

    pd.set_option("display.width", 200)
    print("=" * 96)
    print(f"{a.dataset}: pooled (P) vs within-batch (W) vs cross-batch (X) reproducibility reference")
    print("=" * 96)
    print(f"cells {X.shape[0]:,}  features {X.shape[1]}  batches {meta['n_batches']}  "
          f"groups {meta['group_sizes']}  control cells {meta['n_control_cells']:,} "
          f"({meta['control_cells_per_group']} per group)")
    cols = [c for c in ("depth", "n_perturbations", "pds_P", "pds_P_lo", "pds_P_hi", "pds_W",
                        "pds_X", "pds_X_lo", "pds_X_hi", "delta_PX", "delta_PX_lo",
                        "delta_PX_hi", "delta_PX_signflip_p", "frac_perturbations_delta_gt0",
                        "verdict_flip_frac", "verdict") if c in summary.columns]
    print(summary[cols].to_string(index=False, float_format=lambda v: f"{v:.4f}"))


if __name__ == "__main__":
    main()
