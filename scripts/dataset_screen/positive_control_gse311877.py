#!/usr/bin/env python3
"""Positive control for predict_gse311877.py.

A null result is only informative if the machinery that produced it can detect a
signal that is definitely there. This script feeds the same leave-one-out loop,
the same estimators and the same scoring a feature set that is predictive BY
CONSTRUCTION, and checks that the reported PDS climbs toward the measurement
ceiling instead of sitting at the null.

The control feature for allele `a` is a low-dimensional summary of `a`'s OWN
measured residual. That is deliberate leakage: it is not a scientific result and
must never be quoted as one. Its only job is to answer "if a feature really
carried allele identity, would this pipeline say so?".

Two arms:
  cheat_resid   features = the allele's own residual, reduced to 10 dims
  cheat_noisy   the same, with Gaussian noise at several signal-to-noise ratios,
                to show the score degrades smoothly rather than being all or nothing

Interpretation. If `cheat_resid` lands near the ceiling and the noisy arms
interpolate down toward the null, the pipeline is sound and any at-null result
from the real features is a statement about the features, not about the code. If
`cheat_resid` itself sits at the null, the pipeline is broken and every number in
the main run is void.
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import predict_gse311877 as P  # noqa: E402

OUT = P.OUT
SNRS = (8.0, 4.0, 2.0, 1.0, 0.5, 0.25)


def main() -> None:
    rng = np.random.default_rng(P.SEED)
    OUT.mkdir(parents=True, exist_ok=True)
    cpm, alleles = P.load_cpm()
    folds = P.precompute(cpm, alleles)
    n = len(alleles)

    ceil = np.array([f["ceiling_pds"] for f in folds])
    C = float(ceil.mean())
    print(f"tag={P.TAG}  alleles={n}  measurement ceiling PDS={C:.4f}")

    # Build the cheat feature from the full-depth residuals. Fold 0's basis is used
    # only to define a common coordinate system for the feature matrix; the LOO loop
    # still refits everything per fold. This is a control, not an estimate.
    f0 = folds[0]
    # reconstruct full-depth residuals for all alleles in a single shared space
    x_all = np.log2(cpm + 1.0)
    keep = (cpm >= P.MIN_CPM).mean(axis=1) >= P.MIN_FRAC
    x = x_all.loc[keep]
    base = x[[f"WT_{r}" for r in (1, 2, 3, 4)]].mean(axis=1).to_numpy()
    D = np.stack([x[[f"{a}_{r}" for r in (1, 2, 3, 4)]].mean(axis=1).to_numpy() - base
                  for a in alleles])
    Dc = D - D.mean(axis=0, keepdims=True)
    U = np.linalg.svd(Dc.T, full_matrices=False)[0][:, :P.N_PC]
    if P.REMOVE_MEAN:
        m = D.mean(axis=0)
        m = m - (m @ U) @ U.T
        U = np.hstack([U, (m / np.linalg.norm(m))[:, None]])
    R = D - (D @ U) @ U.T
    Rc = R - R.mean(axis=0, keepdims=True)
    V = np.linalg.svd(Rc, full_matrices=False)[2][:10].T
    cheat = Rc @ V                                          # 18 x 10
    cheat = cheat / (cheat.std(axis=0, keepdims=True) + 1e-12)

    rows = []
    for name, X in [("cheat_resid", cheat)] + [
            (f"cheat_snr{s:g}", cheat + rng.normal(0, 1.0 / s, cheat.shape)) for s in SNRS]:
        F = {name: X, "region": np.zeros((n, 1)), "position": np.zeros((n, 1))}
        for model in ("ridge", "krr", "pls", "nn"):
            res = P.run_combo(folds, alleles, F, name, model)
            obs = float(res["pds"].mean())
            rows.append({"features": name, "model": model, "PDS": obs,
                         "top1": float(res["top1"].mean()),
                         "resid_r": float(res["r"].mean()),
                         "frac_of_ceiling": (obs - 0.5) / (C - 0.5)})
            print(f"  {name:<16}{model:<7} PDS={obs:.4f}  top1={rows[-1]['top1']:.3f}  "
                  f"r={rows[-1]['resid_r']:+.4f}  frac_of_ceiling={rows[-1]['frac_of_ceiling']:+.3f}")

    best = max(rows, key=lambda r: r["PDS"] if r["features"] == "cheat_resid" else -np.inf)
    ok = best["PDS"] > 0.5 + 0.5 * (C - 0.5)
    rep = {"tag": P.TAG, "ceiling": C, "best_cheat": best, "pipeline_detects_signal": bool(ok),
           "rows": rows}
    (OUT / f"positive_control_{P.TAG}.json").write_text(json.dumps(rep, indent=2))
    print(f"\npipeline recovers a planted signal above half the headroom: {ok}")
    if not ok:
        print("  -> the main run's null result would be UNINTERPRETABLE. Investigate before "
              "quoting any verdict.")


if __name__ == "__main__":
    main()
