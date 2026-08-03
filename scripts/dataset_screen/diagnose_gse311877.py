#!/usr/bin/env python3
"""Why the positive control failed: is the model class able to reach the ceiling
at all, at n = 18?

positive_control_gse311877.py planted a feature that is the allele's own residual
and the pipeline still scored near the null. Before any verdict from the main run
can be read, that has to be explained. Two candidate explanations:

  (i)  the code is wrong;
  (ii) the model class is structurally incapable here, so "model at chance" is
       guaranteed by geometry and carries no evidence about biology.

Every estimator in the main run predicts `r_hat_a = c^T Yc`, a linear combination
of the OTHER 17 alleles' residuals. So the very best it could ever do is the
projection of allele `a`'s residual onto the span of the other 17. If the 18
residuals are close to mutually orthogonal, that projection is essentially noise
and NO feature, however good, can identify the held-out allele. That is
explanation (ii), and it is testable in closed form.

This script computes:

  1. the Gram structure of the 18 full-depth residuals (are they orthogonal?);
  2. the fraction of each residual's norm that lies inside the span of the other
     17, which is the hard cap on any linear smoother;
  3. the ANALYTIC linear-smoother oracle: `c* ∝ K^+ M[:, a]` maximises
     `corr(c^T Yc, pool_a)` exactly, so its PDS upper-bounds the entire model
     class including every model in the pre-registration;
  4. the same quantities for the measurement ceiling, whose predictor is allele
     `a`'s own half-A and is therefore NOT confined to the training span.

If the oracle is at chance, the design cannot answer "can features predict this",
and the honest verdict is INCONCLUSIVE for a structural reason rather than a
power one.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import predict_gse311877 as P  # noqa: E402


def main() -> None:
    cpm, alleles = P.load_cpm()
    folds = P.precompute(cpm, alleles)
    n = len(alleles)
    print(f"tag={P.TAG}  alleles={n}")

    # ---- 1/2. residual geometry in a single shared space (descriptive) ----
    x = np.log2(cpm + 1.0).loc[(cpm >= P.MIN_CPM).mean(axis=1) >= P.MIN_FRAC]
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
    Rn = R / np.linalg.norm(R, axis=1, keepdims=True)
    Gram = Rn @ Rn.T
    off = Gram[~np.eye(n, dtype=bool)]
    print(f"\n[1] residual pairwise correlation: mean={off.mean():+.4f}  "
          f"sd={off.std():.4f}  min={off.min():+.4f}  max={off.max():+.4f}")

    # fraction of each residual captured by the span of the other 17
    frac = []
    for i in range(n):
        others = np.delete(Rn, i, axis=0)
        Q = np.linalg.svd(others.T, full_matrices=False)[0]
        frac.append(float(np.linalg.norm(Q.T @ Rn[i]) ** 2))
    frac = np.asarray(frac)
    print(f"[2] fraction of a residual's variance inside span(other 17): "
          f"mean={frac.mean():.4f}  min={frac.min():.4f}  max={frac.max():.4f}")
    print("    per allele: " + "  ".join(f"{a}={f:.3f}" for a, f in zip(alleles, frac)))

    # ---- 3. analytic linear-smoother oracle ----
    orc_pds, orc_t1, orc_r, ceil_pds = [], [], [], []
    for f in folds:
        ai, K = f["index"], f["K"]
        Kp = np.linalg.pinv(K, rcond=1e-10)
        pds_s, t1_s, r_s = [], [], []
        for M in f["M"]:
            c = Kp @ M[:, ai]                      # maximises corr with pool_a exactly
            den = np.sqrt(max(c @ K @ c, 1e-12))
            s = (c @ M) / den
            p, t = P.pds_from_scores(s, ai)
            pds_s.append(p)
            t1_s.append(t)
            r_s.append(s[ai])
        orc_pds.append(np.mean(pds_s))
        orc_t1.append(np.mean(t1_s))
        orc_r.append(np.mean(r_s))
        ceil_pds.append(f["ceiling_pds"])
    orc_pds = np.asarray(orc_pds)
    C = float(np.mean(ceil_pds))
    print(f"\n[3] ANALYTIC linear-smoother oracle (upper bound on every model here)")
    print(f"    PDS={orc_pds.mean():.4f}   top1={np.mean(orc_t1):.4f}   "
          f"corr with own pool profile={np.mean(orc_r):+.4f}")
    print(f"[4] measurement ceiling (predictor = allele's own half-A, NOT span-limited)")
    print(f"    PDS={C:.4f}   top1={np.mean([f['ceiling_top1'] for f in folds]):.4f}")

    reachable = float((orc_pds.mean() - 0.5) / (C - 0.5)) if C > 0.5 else float("nan")
    print(f"\n    fraction of the measurement headroom that ANY model in the frozen "
          f"list could reach: {reachable:+.3f}")
    verdict = ("model class CAN reach a useful part of the ceiling"
               if orc_pds.mean() > 0.5 + 0.25 * (C - 0.5) else
               "model class CANNOT reach the ceiling: 'model at chance' is forced by "
               "geometry, not evidence about features")
    print(f"    -> {verdict}")

    Path(P.OUT).mkdir(parents=True, exist_ok=True)
    (Path(P.OUT) / f"diagnose_{P.TAG}.json").write_text(json.dumps({
        "tag": P.TAG, "n_alleles": n,
        "residual_pairwise_r_mean": float(off.mean()), "residual_pairwise_r_sd": float(off.std()),
        "frac_in_span_of_others_mean": float(frac.mean()),
        "frac_in_span_of_others_per_allele": dict(zip(alleles, map(float, frac))),
        "oracle_pds": float(orc_pds.mean()), "oracle_top1": float(np.mean(orc_t1)),
        "ceiling_pds": C, "reachable_fraction_of_headroom": reachable,
        "conclusion": verdict}, indent=2))


if __name__ == "__main__":
    main()
