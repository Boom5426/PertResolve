#!/usr/bin/env python3
"""JAK1_PREDICTION_INTERFACE_SANITY_CHECK_v1.

Executes docs/PREREG_JAK1_PREDICTION_INTERFACE_SANITY_CHECK_v1.md. Read that file
first; its three checks and verdict rules are binding and applied here unchanged.

Why this exists. JAK1 is the manuscript's only measurement-resolved but
model-limited example. A screen of a different dataset (GSE311877) showed that a
model sitting at chance is not automatically evidence about features: when an
estimator can only emit points inside the span of the training variants' observed
profiles, and that span is small, "model at chance" is forced by geometry. JAK1
split3 trains on 6 variants, so the exact linear smoothers in the in-house grid
have a 6-dimensional reachable subspace inside 2,000-dimensional gene space,
while that split carries 20 of the 59 scored cells.

What this does NOT do: it trains nothing, adds no method, redefines no metric,
changes no split. Everything mirrors oracle_ceiling.py and score_definitive.py
exactly (same SPLITS, NSUB, NSEED, NBOOT, tie-aware mid-rank PDS-cosine,
candidate set = train + test, only test variants scored), so every number here is
directly comparable to the published JAK1 ceiling 0.792 and best in-house model
0.517.

Usage (remote, `Agent` env):
    python jak1_interface_check.py --base /path/to/processed-data --out jak1_interface_check
    ALLELEPERTURB_DATA=/path/to/processed-data python jak1_interface_check.py --out ...
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from alleleperturb.paths import add_harness_to_path, require_inputs, resolve_base  # noqa: E402

ap = argparse.ArgumentParser(description=__doc__,
                             formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--base", help="processed-data directory (env: ALLELEPERTURB_DATA)")
ap.add_argument("--out", required=True, help="directory for the check's CSV/JSON output")
args = ap.parse_args()

BASE = resolve_base(args.base)
OUT = Path(args.out).expanduser().resolve()
require_inputs(BASE / "unified" / "harness.py", BASE / "unified" / "real_deltas.npz",
               BASE / "allele_perturb_bench.csv")
OUT.mkdir(parents=True, exist_ok=True)

add_harness_to_path(BASE)
import harness as H  # noqa: E402

GENE = "JAK1"
NSUB, NSEED, NBOOT = 300, 15, 2000
SPLITS = ["split1", "split2", "split3", "split5", "split6"]
WT_TAGS = ("WT", "wt", "WT_control")
SNRS = (8.0, 4.0, 2.0, 1.0, 0.5, 0.25)
# 0.25 of the published JAK1 headroom above chance; see the pre-registration.
CEILING_PUBLISHED = 0.792
BAR = 0.5 + 0.25 * (CEILING_PUBLISHED - 0.5)

np.random.seed(0)
X_ALL, LAB_ALL = H.load_gene(GENE)
LAB_ALL = np.asarray(LAB_ALL)

# The span basis must be the deltas the stored predictions were actually FITTED on,
# which is unified/real_deltas.npz (harness.pseudobulk_deltas, RandomState(20240709 +
# gene seed)) and is FIXED across scoring seeds. Using a seed-matched resampled basis
# instead would draw the basis from the same noise realisation as the target and
# inflate the reachable bound.
_RDZ = np.load(BASE / "unified" / "real_deltas.npz", allow_pickle=True)
FIT_DELTAS = {k.split("__", 1)[1]: _RDZ[k].astype(np.float64)
              for k in _RDZ.files if k.startswith(f"{GENE}__")}
CAND = {s: (H.split_vars(GENE, s)[0] + H.split_vars(GENE, s)[1]) for s in SPLITS}
TRAIN = {s: H.split_vars(GENE, s)[0] for s in SPLITS}
TEST = {s: H.split_vars(GENE, s)[1] for s in SPLITS}


# ----------------------------------------------------------- metric, verbatim
def pds_row(drow: np.ndarray, ci: int, n: int) -> float:
    """Tie-aware mid-rank PDS, byte-for-byte the rule used by score_definitive.py
    and oracle_ceiling.py. Reproduced rather than imported because those scripts
    define it at module scope; it is not re-derived or altered."""
    td = drow[ci]
    if not np.isfinite(td):
        return 0.5
    less = int((drow < td - 1e-12).sum())
    eq = int((np.abs(drow - td) <= 1e-12).sum())
    return 1.0 - (less + (eq - 1) / 2.0) / (n - 1) if n > 1 else 0.5


def top1_row(drow: np.ndarray, ci: int) -> bool:
    """Strict top-1: the true variant is the unique nearest candidate."""
    td = drow[ci]
    return bool((drow < td - 1e-12).sum() == 0 and (np.abs(drow - td) <= 1e-12).sum() == 1)


def norm(M: np.ndarray) -> np.ndarray:
    return M / (np.linalg.norm(M, axis=1, keepdims=True) + 1e-12)


# ------------------------------------------------- per-seed deltas and halves
def deltas_seed(seed: int) -> dict[str, np.ndarray]:
    """Full-depth pseudobulk deltas, identical construction to score_definitive."""
    rng = np.random.RandomState(1000 + 7 * seed + H._GENE_SEED[GENE])
    wt = np.where(np.isin(LAB_ALL, WT_TAGS))[0]
    if len(wt) > NSUB:
        wt = rng.choice(wt, NSUB, replace=False)
    wm = X_ALL[wt].mean(0)
    out = {}
    for v in np.unique(LAB_ALL):
        if v in WT_TAGS:
            continue
        idx = np.where(LAB_ALL == v)[0]
        if len(idx) < 5:
            continue
        if len(idx) > NSUB:
            idx = rng.choice(idx, NSUB, replace=False)
        out[v] = (X_ALL[idx].mean(0) - wm).astype(np.float32)
    return out


def halves_seed(seed: int) -> tuple[dict, dict]:
    """Disjoint split-half deltas with independent WT halves, the published
    ceiling construction from oracle_ceiling.py."""
    rng = np.random.RandomState(1000 + 7 * seed + H._GENE_SEED[GENE])
    wt = np.where(np.isin(LAB_ALL, WT_TAGS))[0]
    rng.shuffle(wt)
    h = len(wt) // 2
    wt_t, wt_p = wt[:h][:NSUB], wt[h:][:NSUB]
    wm_t, wm_p = X_ALL[wt_t].mean(0), X_ALL[wt_p].mean(0)
    truth, pred = {}, {}
    for v in np.unique(LAB_ALL):
        if v in WT_TAGS:
            continue
        idx = np.where(LAB_ALL == v)[0]
        if len(idx) < 10:
            continue
        idx = idx.copy()
        rng.shuffle(idx)
        k = len(idx) // 2
        a, b = idx[:k][:NSUB], idx[k:][:NSUB]
        truth[v] = (X_ALL[a].mean(0) - wm_t).astype(np.float32)
        pred[v] = (X_ALL[b].mean(0) - wm_p).astype(np.float32)
    return truth, pred


def score(pred_of: dict, truth: dict, splits=SPLITS) -> dict:
    """(split, variant) -> (PDS, top1) for one seed, scoring only test variants."""
    cells = {}
    for s in splits:
        cv = [v for v in CAND[s] if v in truth]
        if not cv:
            continue
        idx = {v: i for i, v in enumerate(cv)}
        te = [v for v in TEST[s] if v in idx and v in pred_of]
        if not te:
            continue
        P = norm(np.stack([pred_of[v] for v in te]))
        T = norm(np.stack([truth[u] for u in cv]))
        D = 1 - P @ T.T
        for i, v in enumerate(te):
            cells[(s, v)] = (pds_row(D[i], idx[v], len(cv)), top1_row(D[i], idx[v]))
    return cells


def aggregate(per_seed: list[dict]) -> tuple[float, float, tuple[float, float], dict]:
    """Seed-average per variant, then mean of (split) cell means, with a bootstrap
    over held-out variants. Matches score_definitive.py's aggregation."""
    keys = set().union(*[set(d) for d in per_seed]) if per_seed else set()
    pds = {k: float(np.mean([d[k][0] for d in per_seed if k in d])) for k in keys}
    t1 = {k: float(np.mean([d[k][1] for d in per_seed if k in d])) for k in keys}
    by_split: dict = {}
    for (s, v), p in pds.items():
        by_split.setdefault(s, {})[v] = p
    obs = float(np.mean([np.mean(list(d.values())) for d in by_split.values()]))
    rng = np.random.RandomState(0)
    boot = []
    for _ in range(NBOOT):
        cm = [np.array(list(d.values()))[rng.randint(0, len(d), len(d))].mean()
              for d in by_split.values()]
        boot.append(np.mean(cm))
    lo, hi = np.percentile(boot, [2.5, 97.5])
    top1 = float(np.mean([np.mean([t1[(s, v)] for v in d]) for s, d in by_split.items()]))
    split_means = {s: float(np.mean(list(d.values()))) for s, d in by_split.items()}
    return obs, top1, (float(lo), float(hi)), split_means


def project_onto_span(target: np.ndarray, basis_rows: np.ndarray) -> np.ndarray:
    """Best reachable point for an estimator confined to span(basis_rows).

    Projection maximises cosine similarity to `target` within the span, so its PDS
    upper-bounds every exact linear smoother fitted onto those rows.
    """
    Q, _ = np.linalg.qr(basis_rows.T)
    return Q @ (Q.T @ target)


def main() -> None:
    rep: dict = {"gene": GENE, "nseed": NSEED, "nboot": NBOOT,
                 "published_ceiling": CEILING_PUBLISHED, "bar_0.25H": BAR}
    print(f"{GENE}: splits " + ", ".join(
        f"{s} train {len(TRAIN[s])}/test {len(TEST[s])}" for s in SPLITS))
    print(f"bar for the verdict, 0.5 + 0.25*(0.792-0.5) = {BAR:.4f}\n")

    # ---------------------------------------------------------------- C1 + C2
    arms: dict[str, list] = {k: [] for k in
                             ["truth", "splithalf", "span_oracle"] +
                             [f"truth_snr{s:g}" for s in SNRS]}
    outside = {}
    for seed in range(NSEED):
        d = deltas_seed(seed)
        rng = np.random.RandomState(500 + seed)

        arms["truth"].append(score({v: d[v] for v in d}, d))

        for s_ in SNRS:
            scale = np.array([np.linalg.norm(d[v]) for v in d]).mean() / max(s_, 1e-9)
            noisy = {v: d[v] + rng.normal(0, scale / np.sqrt(len(d[v])), d[v].shape).astype(np.float32)
                     for v in d}
            arms[f"truth_snr{s_:g}"].append(score(noisy, d))

        t, p = halves_seed(seed)
        arms["splithalf"].append(score(p, t))

        # Span oracle: the best point reachable inside the span of the deltas the
        # stored models were fitted on. Basis is FIT_DELTAS (fixed); the target is
        # this seed's truth, exactly as a real model is scored.
        cells = {}
        for s in SPLITS:
            tr = [v for v in TRAIN[s] if v in FIT_DELTAS]
            te = [v for v in TEST[s] if v in d]
            if len(tr) < 1 or not te:
                continue
            B = np.stack([FIT_DELTAS[v] for v in tr])
            proj = {v: project_onto_span(d[v].astype(np.float64), B) for v in te}
            cells.update(score(proj, d, splits=[s]))
        arms["span_oracle"].append(cells)

    rows = []
    for name, per_seed in arms.items():
        if not any(per_seed):
            continue
        obs, t1, (lo, hi), sm = aggregate(per_seed)
        rows.append({"arm": name, "PDS": round(obs, 4), "ci_lo": round(lo, 4),
                     "ci_hi": round(hi, 4), "top1": round(t1, 4),
                     **{f"PDS_{s}": round(v, 4) for s, v in sorted(sm.items())}})
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "jak1_interface_arms.csv", index=False)
    print(df.to_string(index=False))

    truth_row = df[df.arm == "truth"].iloc[0]
    sh_row = df[df.arm == "splithalf"].iloc[0]
    so_row = df[df.arm == "span_oracle"].iloc[0]

    # ------------------------------------ independent re-derivation of the span
    pred_dir = BASE / "unified" / "preds5"
    if pred_dir.exists():
        for f in sorted(pred_dir.glob("*.npz")):
            z = np.load(f, allow_pickle=True)
            fr = {}
            for s in SPLITS:
                tr = [v for v in TRAIN[s] if v in FIT_DELTAS]
                te = [v for v in TEST[s] if f"{GENE}__{s}__{v}" in z.files]
                if not tr or not te:
                    continue
                Q, _ = np.linalg.qr(np.stack([FIT_DELTAS[v] for v in tr]).T)
                vals = []
                for v in te:
                    p = z[f"{GENE}__{s}__{v}"].astype(np.float64)
                    n = np.linalg.norm(p)
                    if n < 1e-12 or not np.isfinite(n):
                        continue
                    vals.append(float(np.linalg.norm(p - Q @ (Q.T @ p)) / n))
                if vals:
                    fr[s] = round(float(np.median(vals)), 5)
            if fr:
                outside[f.stem] = fr
        pd.DataFrame(outside).T.to_csv(OUT / "jak1_outside_span_fraction.csv")
        print("\noutside-span fraction of stored predictions (median over held-out variants)")
        print(pd.DataFrame(outside).T.to_string())

    # ------------------------------------------------------- verdict (binding)
    c1_pass = bool(truth_row.PDS > 0.999 and truth_row.top1 > 0.999)
    c1_halfway = bool(sh_row.ci_lo > 0.5)
    oracle_ok = bool(so_row.ci_lo > BAR)
    split3_ok = bool(float(so_row.get("PDS_split3", np.nan)) > BAR)
    top1_reachable = bool(so_row.top1 > 1.0 / len(CAND["split1"]))

    if not c1_pass:
        verdict = "VOID_SCORER_BROKEN"
    elif oracle_ok and top1_reachable:
        verdict = "INTERFACE_OK"
    elif oracle_ok:
        verdict = "PDS_OK_TOP1_UNREACHABLE"
    else:
        verdict = "INTERFACE_UNREACHABLE"

    rep.update({"C1_truth_pds": float(truth_row.PDS), "C1_truth_top1": float(truth_row.top1),
                "C1_pass": c1_pass,
                "C1_splithalf_pds": float(sh_row.PDS),
                "C1_splithalf_ci": [float(sh_row.ci_lo), float(sh_row.ci_hi)],
                "C1_splithalf_above_chance": c1_halfway,
                "C2_span_oracle_pds": float(so_row.PDS),
                "C2_span_oracle_ci": [float(so_row.ci_lo), float(so_row.ci_hi)],
                "C2_span_oracle_above_bar": oracle_ok,
                "C2_split3_pds": float(so_row.get("PDS_split3", np.nan)),
                "C2_split3_above_bar": split3_ok,
                "C3_span_oracle_top1": float(so_row.top1),
                "C3_top1_reachable": top1_reachable,
                "outside_span_fraction": outside,
                "verdict": verdict})
    (OUT / "jak1_interface_check.json").write_text(json.dumps(rep, indent=2))

    print("\n" + "=" * 78)
    print(f"C1 truth control     PDS={truth_row.PDS:.4f} top1={truth_row.top1:.4f} "
          f"-> pass={c1_pass}")
    print(f"C1 split-half        PDS={sh_row.PDS:.4f} [{sh_row.ci_lo:.4f},{sh_row.ci_hi:.4f}]"
          f"   published 0.792")
    print(f"C2 span oracle       PDS={so_row.PDS:.4f} [{so_row.ci_lo:.4f},{so_row.ci_hi:.4f}]"
          f"   bar {BAR:.4f} -> above={oracle_ok}")
    print(f"C2 split3 only       PDS={so_row.get('PDS_split3', float('nan')):.4f}"
          f" -> above bar={split3_ok}   (6 training variants, 20 of 59 scored cells)")
    print(f"C3 span oracle top1  {so_row.top1:.4f}  chance {1/len(CAND['split1']):.4f}"
          f" -> reachable={top1_reachable}")
    print(f"\nVERDICT: {verdict}")
    print("=" * 78)


if __name__ == "__main__":
    main()
