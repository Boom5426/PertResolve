#!/usr/bin/env python3
"""How much of a reproducibility reference is lost when the two sides stop sharing a
physical sample? TP63 (GEO GSE311877), 18 alleles x 4 culture-well replicates.

Frozen protocol: docs/PREREG_REPLICATE_REFERENCE_v1.md, section 9. This is the
*bounding* analysis and enters no verdict rule. Read that section before any number
here.

GSE311877 is bulk 3'-digital RNA-seq, so the sampling unit is a molecule, not a cell.
Two arms at exactly matched molecular depth N:

    S  same library    two disjoint molecule halves of ONE library
    B  cross-replicate one library each from two different culture wells

Arm S's noise is a strict subset of arm B's, so `PDS_S >= PDS_B` holds by
construction and its *sign* is not a finding. What is reportable is the magnitude,
the share of between-measurement variance that arm S omits, and how both behave as
depth falls.

Three corrections over a naive implementation, each of which changes the answer:

* Depth is matched exactly by multivariate hypergeometric subsampling. Binomial
  thinning at p = N/total does not give exactly N molecules.
* The gene filter is computed once, at full depth, from the WT and GFP libraries
  only. Those are never a target and never a pool member, so the filter is both
  label-free and target-free.
* log2(CPM + 1) is not comparable along a depth ladder: one molecule is 500 CPM at
  N = 2,000 and 2.5 CPM at N = 400,000, so the pseudocount means something different
  at every rung. The primary transform is log2((y + 1) * 1e6 / N), whose depth term
  is gene-independent and cancels in the WT-subtracted delta.
"""
from __future__ import annotations

import argparse
import gzip
import itertools
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from pertresolve.paths import reject_repo_results, require_inputs
from pertresolve.resolution.scaling import tie_aware_pds

MASTER_SEED = 20260825
FNAME = re.compile(r"^GSM\d+_(?:([A-H]\d+)_)?(.+?)_(\d)_rawCounts\.txt\.gz$")
CONTROLS = ("WT", "GFP")
SIB_PAIRS = [("R279Q", "R279S"), ("R304Q", "R304T"), ("L514D", "L514F"),
             ("C522D", "C522G"), ("L531E", "L531R")]
DEPTHS = [250, 500, 1_000, 2_000, 5_000, 10_000, 20_000, 50_000, 100_000, 200_000, 400_000]
MIN_CPM, MIN_FRAC = 1.0, 0.5


def load(raw: Path):
    cols, meta = {}, []
    for f in sorted(raw.glob("*_rawCounts.txt.gz")):
        m = FNAME.match(f.name)
        if not m:
            raise ValueError(f"unparsed sample filename: {f.name}")
        well, cond, rep = m.group(1), m.group(2), int(m.group(3))
        with gzip.open(f, "rt") as fh:
            cols[f"{cond}_{rep}"] = pd.read_csv(fh, sep="\t", index_col=0).iloc[:, 0]
        meta.append({"sample": f"{cond}_{rep}", "condition": cond, "rep": rep,
                     "well": well or "", "row": (well or " ")[0], "col": (well or "  ")[1:]})
    counts = pd.DataFrame(cols)
    md = pd.DataFrame(meta).set_index("sample").loc[counts.columns]
    return counts, md


def control_only_filter(counts: pd.DataFrame) -> np.ndarray:
    """Genes kept, decided at full depth from the WT and GFP libraries alone.

    Those eight libraries are never scored and never enter the candidate pool, so a
    filter derived from them cannot leak either allele identity or target-side
    structure into the identification test.
    """
    ctrl = [c for c in counts.columns if c.rsplit("_", 1)[0] in CONTROLS]
    sub = counts[ctrl]
    cpm = sub / sub.sum(axis=0) * 1e6
    return ((cpm >= MIN_CPM).mean(axis=1) >= MIN_FRAC).to_numpy()


def transform(y: np.ndarray, n: int, kind: str) -> np.ndarray:
    if kind == "shifted_log_count":
        return np.log2((y + 1.0) * 1e6 / n)
    if kind == "log2_cpm1":
        return np.log2(y / max(n, 1) * 1e6 + 1.0)
    if kind == "freeman_tukey":
        return np.sqrt(y) + np.sqrt(y + 1.0)
    raise ValueError(kind)


def draw_halves(rng, c: np.ndarray, n: int):
    """Two disjoint blocks of exactly n molecules from one library.

    One joint hypergeometric partition, so disjointness is structural. Requires
    c.sum() >= 2n; eligibility is enforced by the caller at 4n.
    """
    s = rng.multivariate_hypergeometric(c, 2 * n)
    h1 = rng.multivariate_hypergeometric(s, n)
    return h1, s - h1


def draw_one(rng, c: np.ndarray, n: int) -> np.ndarray:
    return rng.multivariate_hypergeometric(c, n)


def per_allele_pds(pred: np.ndarray, truth: np.ndarray) -> np.ndarray:
    """Tie-aware cosine PDS per row. Rows are alleles, columns are genes."""
    qn = pred / (np.linalg.norm(pred, axis=1, keepdims=True) + 1e-12)
    tn = truth / (np.linalg.norm(truth, axis=1, keepdims=True) + 1e-12)
    dist = 1.0 - qn @ tn.T
    n = dist.shape[1]
    out = np.empty(dist.shape[0])
    for i in range(dist.shape[0]):
        row = dist[i]
        t = row[i]
        less = int((row < t - 1e-12).sum())
        eq = int((np.abs(row - t) <= 1e-12).sum())
        out[i] = 1.0 - (less + (eq - 1) / 2.0) / (n - 1)
    return out


def top1(pred: np.ndarray, truth: np.ndarray) -> np.ndarray:
    qn = pred / (np.linalg.norm(pred, axis=1, keepdims=True) + 1e-12)
    tn = truth / (np.linalg.norm(truth, axis=1, keepdims=True) + 1e-12)
    dist = 1.0 - qn @ tn.T
    return (np.argmin(dist, axis=1) == np.arange(dist.shape[0])).astype(float)


def residualise(pred, truth, basis, k):
    if k == 0:
        return pred, truth
    bm = basis - basis.mean(axis=0, keepdims=True)      # centre over alleles
    u = np.linalg.svd(bm.T, full_matrices=False)[0][:, :k]   # gene-space basis
    return pred - (pred @ u) @ u.T, truth - (truth @ u) @ u.T


def build(counts_arr, idx, alleles, rng, depth, kind, arm, rep_pred, rep_tgt):
    """Prediction and target delta matrices, alleles x genes, at exactly `depth`."""
    n_g = counts_arr.shape[0]
    pred = np.empty((len(alleles), n_g))
    truth = np.empty_like(pred)
    if arm == "S":
        wa, wb = draw_halves(rng, counts_arr[:, idx[f"WT_{rep_pred}"]], depth)
    else:
        wa = draw_one(rng, counts_arr[:, idx[f"WT_{rep_pred}"]], depth)
        wb = draw_one(rng, counts_arr[:, idx[f"WT_{rep_tgt}"]], depth)
    ta, tb = transform(wa, depth, kind), transform(wb, depth, kind)
    for j, a in enumerate(alleles):
        if arm == "S":
            va, vb = draw_halves(rng, counts_arr[:, idx[f"{a}_{rep_pred}"]], depth)
        else:
            va = draw_one(rng, counts_arr[:, idx[f"{a}_{rep_pred}"]], depth)
            vb = draw_one(rng, counts_arr[:, idx[f"{a}_{rep_tgt}"]], depth)
        pred[j] = transform(va, depth, kind) - ta
        truth[j] = transform(vb, depth, kind) - tb
    return pred, truth


def spare_basis(counts_arr, idx, alleles, rng, depth, kind, used):
    """Shared-axis basis from replicate libraries used by neither side.

    Taking it from the prediction side is not symmetric between the arms: in arm S
    the prediction shares a library with its target, in arm B it does not. That
    asymmetry alone produced an apparent verdict flip in an earlier exploratory run.
    """
    spare = [r for r in (1, 2, 3, 4) if r not in used][:2]
    n_g = counts_arr.shape[0]
    b = np.empty((len(alleles), n_g))
    w = np.mean([transform(draw_one(rng, counts_arr[:, idx[f"WT_{r}"]], depth), depth, kind)
                 for r in spare], axis=0)
    for j, a in enumerate(alleles):
        b[j] = np.mean([transform(draw_one(rng, counts_arr[:, idx[f"{a}_{r}"]], depth), depth, kind)
                        for r in spare], axis=0) - w
    return b


def sibling_margins(pred, truth, alleles):
    """Per-pair mean of (self similarity - sibling similarity).

    The unit is the pair, not the member: the two members of a pair share both
    alleles' noise and are not independent trials.
    """
    qn = pred / (np.linalg.norm(pred, axis=1, keepdims=True) + 1e-12)
    tn = truth / (np.linalg.norm(truth, axis=1, keepdims=True) + 1e-12)
    sim = qn @ tn.T
    ai = {a: i for i, a in enumerate(alleles)}
    out = []
    for p, q in SIB_PAIRS:
        i, j = ai[p], ai[q]
        out.append(0.5 * ((sim[i, i] - sim[i, j]) + (sim[j, j] - sim[j, i])))
    return np.array(out)


def exact_signflip_p(x: np.ndarray) -> float:
    """Two-sided exact sign-flip p over all 2**len(x) patterns. Floor is 2/2**len(x)."""
    n = len(x)
    obs = abs(x.mean())
    hits = 0
    for bits in itertools.product((-1.0, 1.0), repeat=n):
        if abs(float(np.dot(bits, x)) / n) >= obs - 1e-12:
            hits += 1
    return hits / 2 ** n


def variance_decomposition(arr, idx, alleles, depth, kind, draws, seed):
    """Split the between-measurement variance of an allele delta into its two parts.

    Threshold-free and basis-free, so it survives every choice the PDS ladder has to
    make. Same-library halves give the molecule-sampling component; different
    libraries give molecule sampling plus across-well. Both at the same depth, so the
    difference is the component a same-sample split-half reference never sees.
    """
    rng = np.random.default_rng(seed)
    mol, tot = [], []
    for a in alleles:
        m, t = [], []
        for _ in range(draws):
            for r in (1, 2, 3, 4):
                va, vb = draw_halves(rng, arr[:, idx[f"{a}_{r}"]], depth)
                wa, wb = draw_halves(rng, arr[:, idx[f"WT_{r}"]], depth)
                d1 = transform(va, depth, kind) - transform(wa, depth, kind)
                d2 = transform(vb, depth, kind) - transform(wb, depth, kind)
                m.append(np.var(d1 - d2) / 2.0)
            for r1, r2 in itertools.combinations((1, 2, 3, 4), 2):
                d1 = (transform(draw_one(rng, arr[:, idx[f"{a}_{r1}"]], depth), depth, kind)
                      - transform(draw_one(rng, arr[:, idx[f"WT_{r1}"]], depth), depth, kind))
                d2 = (transform(draw_one(rng, arr[:, idx[f"{a}_{r2}"]], depth), depth, kind)
                      - transform(draw_one(rng, arr[:, idx[f"WT_{r2}"]], depth), depth, kind))
                t.append(np.var(d1 - d2) / 2.0)
        mol.append(np.mean(m))
        tot.append(np.mean(t))
    mol, tot = np.array(mol), np.array(tot)
    well = tot - mol
    frac = well / tot
    return {"depth": depth, "transform": kind, "n_alleles": len(alleles),
            "s2_molecule_mean": float(mol.mean()), "s2_well_mean": float(well.mean()),
            "s2_total_mean": float(tot.mean()),
            "frac_across_well_mean": float(frac.mean()),
            "frac_across_well_median": float(np.median(frac)),
            "frac_across_well_min": float(frac.min()), "frac_across_well_max": float(frac.max())}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--raw", default="data_external/GSE311877")
    ap.add_argument("--out", required=True)
    ap.add_argument("--draws", type=int, default=10)
    ap.add_argument("--n-boot", type=int, default=10000)
    ap.add_argument("--transforms", default="shifted_log_count,log2_cpm1,freeman_tukey")
    ap.add_argument("--residual-k", default="0,1,5")
    a = ap.parse_args()

    raw = Path(a.raw)
    require_inputs(raw)
    out = reject_repo_results(a.out)
    out.mkdir(parents=True, exist_ok=True)

    counts, md = load(raw)
    keep = control_only_filter(counts)
    arr = np.ascontiguousarray(counts.values[keep].astype(np.int64))
    idx = {c: i for i, c in enumerate(counts.columns)}
    alleles = sorted(set(md["condition"]) - set(CONTROLS))
    lib = counts.sum(axis=0)
    depths = [d for d in DEPTHS if lib.min() >= 4 * d]
    ks = [int(k) for k in a.residual_k.split(",")]
    kinds = a.transforms.split(",")

    meta = {"n_libraries": int(counts.shape[1]), "n_genes_raw": int(counts.shape[0]),
            "n_genes_kept": int(keep.sum()), "gene_filter": "WT+GFP libraries only, full depth",
            "n_alleles": len(alleles), "alleles": alleles,
            "lib_min": int(lib.min()), "lib_median": int(lib.median()), "lib_max": int(lib.max()),
            "depths": depths, "eligibility": "library total >= 4 * depth, all 80 libraries",
            "dropped_depths": [d for d in DEPTHS if d not in depths],
            "draws": a.draws, "master_seed": MASTER_SEED}

    cfg_s = [(r, None) for r in (1, 2, 3, 4)]
    cfg_b = [(i, j) for i, j in itertools.permutations((1, 2, 3, 4), 2)]
    rows, sib_rows = [], []

    for kind in kinds:
        for depth in depths:
            for k in ks:
                accs = {"S": [], "B": []}
                sibs = {"S": [], "B": []}
                for d in range(a.draws):
                    seed = np.random.SeedSequence([MASTER_SEED, depth, d, k, sum(map(ord, kind))])
                    for arm, cfgs in (("S", cfg_s), ("B", cfg_b)):
                        pa, t1a, sa = [], [], []
                        for rp, rt in cfgs:
                            rng = np.random.default_rng(seed.spawn(1)[0])
                            pred, truth = build(arr, idx, alleles, rng, depth, kind, arm, rp, rt)
                            if k:
                                used = {rp} if rt is None else {rp, rt}
                                basis = spare_basis(arr, idx, alleles, rng, depth, kind, used)
                                pred, truth = residualise(pred, truth, basis, k)
                            pa.append(per_allele_pds(pred, truth))
                            t1a.append(top1(pred, truth))
                            sa.append(sibling_margins(pred, truth, alleles))
                        accs[arm].append((np.mean(pa, axis=0), np.mean(t1a, axis=0)))
                        sibs[arm].append(np.mean(sa, axis=0))
                for arm in ("S", "B"):
                    p = np.mean([x[0] for x in accs[arm]], axis=0)
                    t = np.mean([x[1] for x in accs[arm]], axis=0)
                    for j, al in enumerate(alleles):
                        rows.append({"transform": kind, "depth": depth, "residual_k": k,
                                     "arm": arm, "allele": al, "pds": float(p[j]),
                                     "top1": float(t[j]), "n_pool": len(alleles)})
                    m = np.mean(sibs[arm], axis=0)
                    for pi, (pp, qq) in enumerate(SIB_PAIRS):
                        sib_rows.append({"transform": kind, "depth": depth, "residual_k": k,
                                         "arm": arm, "pair": f"{pp}/{qq}", "margin": float(m[pi])})

    per = pd.DataFrame(rows)
    sib = pd.DataFrame(sib_rows)
    per.to_csv(out / "tp63_replicate_reference_perallele.csv.gz", index=False)
    sib.to_csv(out / "tp63_replicate_reference_siblings.csv", index=False)

    rng = np.random.default_rng(MASTER_SEED + 2)
    n = len(alleles)
    bi = rng.integers(0, n, size=(a.n_boot, n))
    summ = []
    for (kind, depth, k), g in per.groupby(["transform", "depth", "residual_k"]):
        w = g.pivot(index="allele", columns="arm", values="pds")
        w1 = g.pivot(index="allele", columns="arm", values="top1")
        d = (w["S"] - w["B"]).to_numpy()
        bd = d[bi].mean(1)
        sgn = rng.choice([-1.0, 1.0], size=(a.n_boot, n))
        null = (sgn * d).mean(1)
        s_pair = sib[(sib["transform"] == kind) & (sib["depth"] == depth) & (sib["residual_k"] == k)]
        ms = s_pair[s_pair.arm == "S"]["margin"].to_numpy()
        mb = s_pair[s_pair.arm == "B"]["margin"].to_numpy()
        summ.append({
            "transform": kind, "depth": depth, "residual_k": k, "n_alleles": n,
            "pds_S": float(w["S"].mean()), "pds_B": float(w["B"].mean()),
            "pds_S_lo": float(np.percentile(w["S"].to_numpy()[bi].mean(1), 2.5)),
            "pds_S_hi": float(np.percentile(w["S"].to_numpy()[bi].mean(1), 97.5)),
            "pds_B_lo": float(np.percentile(w["B"].to_numpy()[bi].mean(1), 2.5)),
            "pds_B_hi": float(np.percentile(w["B"].to_numpy()[bi].mean(1), 97.5)),
            "top1_S": float(w1["S"].mean()), "top1_B": float(w1["B"].mean()),
            "delta": float(d.mean()),
            "delta_lo": float(np.percentile(bd, 2.5)), "delta_hi": float(np.percentile(bd, 97.5)),
            "signflip_p": float(((np.abs(null) >= abs(d.mean())).sum() + 1) / (a.n_boot + 1)),
            "alleles_delta_gt0": int((d > 0).sum()),
            "sib_margin_S": float(ms.mean()), "sib_margin_B": float(mb.mean()),
            "sib_delta_exact_p": exact_signflip_p(ms - mb),
            "sib_exact_p_floor": 2 / 2 ** len(SIB_PAIRS),
        })
    vd = [variance_decomposition(arr, idx, alleles, d, kinds[0], max(2, a.draws // 5),
                                 MASTER_SEED + 7 + d)
          for d in depths[-3:]]
    pd.DataFrame(vd).to_csv(out / "tp63_variance_decomposition.csv", index=False)

    summary = pd.DataFrame(summ).sort_values(["transform", "residual_k", "depth"])
    summary.to_csv(out / "tp63_replicate_reference_summary.csv", index=False)
    (out / "tp63_replicate_reference_meta.json").write_text(json.dumps(meta, indent=2, default=str))

    print("\n-- variance decomposition at matched depth (threshold-free, basis-free) --")
    for v in vd:
        print(f"   depth {v['depth']:>7,}: s2_molecule {v['s2_molecule_mean']:.4f}  "
              f"s2_across_well {v['s2_well_mean']:.4f}  "
              f"share a same-sample split-half OMITS: mean {v['frac_across_well_mean']:.3f} "
              f"(range {v['frac_across_well_min']:.3f} to {v['frac_across_well_max']:.3f})")

    pd.set_option("display.width", 240)
    print("=" * 100)
    print("GSE311877 TP63: same-library (S) vs cross-replicate (B) reference at matched molecular depth")
    print("=" * 100)
    print(f"libraries {meta['n_libraries']}  genes kept {meta['n_genes_kept']}/{meta['n_genes_raw']} "
          f"(WT+GFP filter)  alleles {meta['n_alleles']}")
    print(f"library size min {meta['lib_min']:,} median {meta['lib_median']:,} max {meta['lib_max']:,}")
    print(f"depths (library >= 4N): {depths}   dropped: {meta['dropped_depths']}")
    for kind in kinds:
        for k in ks:
            sub = summary[(summary["transform"] == kind) & (summary["residual_k"] == k)]
            print(f"\n-- transform={kind}  residual_k={k} "
                  f"(basis from 2 spare replicate libraries) --")
            print(sub[["depth", "pds_S", "pds_B", "delta", "delta_lo", "delta_hi", "signflip_p",
                       "alleles_delta_gt0", "top1_S", "top1_B", "sib_margin_S", "sib_margin_B"]]
                  .to_string(index=False, float_format=lambda v: f"{v:.4f}"))


if __name__ == "__main__":
    main()
