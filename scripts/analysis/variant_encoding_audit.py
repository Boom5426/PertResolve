#!/usr/bin/env python3
"""Audit how each variant condition is encoded, and rescore the benchmark without
the conditions whose encoding is degenerate.

Two facts about the input representation are load-bearing for the manuscript's
chance-level result and are not visible from the score tables alone.

First, the six-dimensional biophysical vector is not injective. Its first three
components are differences between residue properties, so any two substitutions
exchanging residues of the same hydrophobicity, volume and charge receive the
same difference, and its last three components are position annotations shared
by many residues. Conditions that are not a single amino-acid substitution,
namely the synonymous, nonsense and splice conditions, have no defined
substitution at all and take a zero difference. Conditions sharing one feature
vector must receive the same prediction from any feature-based head, so a
tie-aware score puts them at chance among themselves by construction.

Second, a synonymous condition has a mutant protein sequence identical to wild
type, so it also receives the wild-type protein-language-model embedding. Under
both representations these conditions carry no information distinguishing them
from one another.

Neither fact is a defect of the measurement, and neither is an excuse: the point
of this script is to report the size of the effect rather than to assume it. It
therefore also recomputes the comparison over the single-substitution missense
conditions alone, the subset on which the three difference terms are defined at
all, so that the chance-level verdict can be read off the conditions the
representation has the best chance of separating.

That subset reduces the ties without removing them, so the script measures the
residual tie structure rather than asserting injectivity. Restricting to the
held-out single-substitution missense conditions is not sufficient for an
injective encoding: two different substitutions can still exchange residues of
equal hydrophobicity, volume and charge at positions carrying the same
annotations. The per-subset rows of ``variant_encoding_audit.csv`` report how
many conditions remain tied, which is the quantity the manuscript needs, and it
is not zero.

Inputs are the committed benchmark table and the committed per-seed per-variant
scores; no expression array is needed.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from pertresolve.paths import reject_repo_results, repo_root, require_inputs

THETA = ['d_hydro', 'd_vol', 'd_charge', 'fold_core', 'cat_switch', 'is_hotspot']
WT_LABELS = {"WT", "WT_CONTROL"}
HEAD_RE = r"Ridge|Lasso|RF|GBoost|KNN|MLP"
NBOOT = 2000


def consequence(variant: str) -> str:
    """Classify a benchmark condition by the consequence its label encodes.

    Args:
        variant: the condition label as it appears in the benchmark table.

    Returns:
        One of ``missense``, ``synonymous``, ``nonsense``, ``splice`` or
        ``multi``. ``multi`` covers labels naming more than one substitution,
        which in this benchmark include both genuinely multi-site conditions and
        conditions pooling two different substitutions at one residue.
    """
    if variant.startswith("splice"):
        return "splice"
    if variant.startswith("syn_"):
        return "synonymous"
    if "," in variant:
        return "multi"
    if variant.endswith("*"):
        return "nonsense"
    m = re.fullmatch(r"([A-Z])(\d+)([A-Z])", variant)
    if m:
        return "synonymous" if m.group(1) == m.group(3) else "missense"
    return "multi"


def load_bench() -> pd.DataFrame:
    """Return the committed benchmark table with wild-type rows removed.

    This reads ``pertresolve_bench.csv`` on purpose, because ``harness.py``
    does, so the tie structure reported here is the one the canonical model grid
    actually encountered. The two committed tables differ only in ``is_hotspot``
    and the choice changes the counts: over all 470 conditions this table gives
    254 distinct vectors and 293 tied, whereas ``pertresolve_bench_v2.csv``
    gives 240 and 303. See ``scripts/analysis/README.md`` for what separates the
    two columns.
    """
    path = repo_root() / "data" / "pertresolve_bench.csv"
    require_inputs(path)
    b = pd.read_csv(path)
    b = b[~b.variant.astype(str).str.upper().isin(WT_LABELS)].copy()
    b["consequence"] = b.variant.map(consequence)
    b["theta_key"] = b[THETA].round(6).astype(str).agg("|".join, axis=1)
    return b


def injectivity_table(subsets: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """How many conditions the feature vector can tell apart, per subset and gene.

    Ties are counted within a gene, because PDS ranks a prediction only against
    the candidate variants of the same gene, so two conditions of different genes
    sharing a feature vector never compete. The ``ALL`` scope therefore sums the
    per-gene counts rather than pooling the conditions and recounting.

    Args:
        subsets: subset name to the conditions it contains. Each frame needs the
            ``gene``, ``consequence`` and ``theta_key`` columns of ``load_bench``.

    Returns:
        One row per (subset, scope), scope being a gene name or ``ALL``.
    """
    rows = []
    for name, bench in subsets.items():
        if bench.empty:
            continue
        per_gene = []
        for gene, s in bench.groupby("gene"):
            counts = s.theta_key.map(s.theta_key.value_counts())
            per_gene.append(dict(
                subset=name,
                scope=gene,
                n_conditions=len(s),
                n_distinct_theta=int(s.theta_key.nunique()),
                n_in_tied_group=int((counts > 1).sum()),
                largest_tied_group=int(s.theta_key.value_counts().max()),
                **{f"n_{c}": int((s.consequence == c).sum())
                   for c in ("missense", "synonymous", "nonsense", "splice", "multi")},
            ))
        per_gene.sort(key=lambda r: r["scope"])
        agg = dict(subset=name, scope="ALL")
        for k in per_gene[0]:
            if k in ("subset", "scope"):
                continue
            agg[k] = (max(r[k] for r in per_gene) if k == "largest_tied_group"
                      else sum(r[k] for r in per_gene))
        rows.extend(per_gene + [agg])
    return pd.DataFrame(rows)


def bootstrap_variant_mean(values: np.ndarray, seed: int = 0):
    """Point estimate and percentile interval resampling distinct variants."""
    rng = np.random.default_rng(seed)
    draws = np.array([values[rng.integers(0, len(values), len(values))].mean()
                      for _ in range(NBOOT)])
    lo, hi = np.percentile(draws, [2.5, 97.5])
    return float(values.mean()), float(lo), float(hi)


def sensitivity_table(scores: pd.DataFrame, bench: pd.DataFrame) -> pd.DataFrame:
    """Rescore the grid over encoding-defined subsets of the conditions."""
    cons = bench.set_index(["gene", "variant"]).consequence
    scores = scores.copy()
    scores["consequence"] = pd.MultiIndex.from_frame(
        scores[["gene", "variant"]]).map(cons)
    seed_avg = (scores.groupby(["method", "gene", "split", "variant", "consequence"])
                ["pds"].mean().reset_index())

    subsets = {
        "all_conditions": lambda d: d,
        "missense_only": lambda d: d[d.consequence == "missense"],
        "single_site_only": lambda d: d[d.consequence != "multi"],
        "excluding_synonymous": lambda d: d[d.consequence != "synonymous"],
    }
    rows = []
    for name, fn in subsets.items():
        sub = fn(seed_avg)
        if sub.empty:
            continue
        for scope in ["ALL"] + sorted(sub.gene.unique()):
            s = sub if scope == "ALL" else sub[sub.gene == scope]
            per_var = s.groupby(["method", "gene", "variant"])["pds"].mean().reset_index()
            heads = per_var[per_var.method.str.contains(HEAD_RE)]
            if heads.empty:
                continue
            by_method = heads.groupby("method")["pds"].mean()
            vals = heads.groupby(["gene", "variant"])["pds"].mean().to_numpy()
            obs, lo, hi = bootstrap_variant_mean(vals)
            rows.append(dict(
                subset=name, scope=scope,
                n_variants=int(heads.groupby(["gene", "variant"]).ngroups),
                heads_mean=round(float(by_method.mean()), 4),
                heads_min=round(float(by_method.min()), 4),
                heads_max=round(float(by_method.max()), 4),
                best_head=by_method.idxmax(),
                pooled_mean=round(obs, 4),
                pooled_lo=round(lo, 4), pooled_hi=round(hi, 4),
                crosses_chance=bool(lo <= 0.5 <= hi),
            ))
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scores", default=None,
                    help="per_seed_variant_pds.csv.gz "
                         "(default: results/canonical/per_seed_variant_pds.csv.gz)")
    ap.add_argument("--out", required=True, metavar="OUT_DIR",
                    help="directory receiving the tables; may not be inside results/")
    args = ap.parse_args()

    out_dir = reject_repo_results(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    bench = load_bench()

    scores_path = Path(args.scores) if args.scores else \
        repo_root() / "results" / "canonical" / "per_seed_variant_pds.csv.gz"
    require_inputs(scores_path)
    scores = pd.read_csv(scores_path)

    # A condition is held out exactly when it carries a score. The manuscript
    # quotes the tie structure of the held-out missense subset, not of the whole
    # benchmark, so that subset has to be measured rather than assumed.
    heldout = set(map(tuple, scores[["gene", "variant"]].drop_duplicates().to_numpy()))
    bench = bench.assign(
        heldout=[(g, v) in heldout for g, v in zip(bench.gene, bench.variant)])

    subsets = {
        "all_conditions": bench,
        "heldout_all": bench[bench.heldout],
        "heldout_missense": bench[bench.heldout & (bench.consequence == "missense")],
    }
    inj = injectivity_table(subsets)
    inj.to_csv(out_dir / "variant_encoding_audit.csv", index=False)
    print(inj.to_string(index=False))

    everything = inj[(inj.subset == "all_conditions") & (inj.scope == "ALL")].iloc[0]
    missense = inj[(inj.subset == "heldout_missense") & (inj.scope == "ALL")].iloc[0]
    print(f"\n{len(bench)} conditions, "
          f"{int((bench.consequence == 'missense').sum())} single-substitution missense; "
          f"{int(everything.n_in_tied_group)} share a feature vector with a same-gene "
          f"sibling. Over the {int(missense.n_conditions)} held-out single-substitution "
          f"missense conditions {int(missense.n_in_tied_group)} are still tied, across "
          f"{int(missense.n_distinct_theta)} distinct vectors, so that subset reduces the "
          f"ties without being injective.")

    sens = sensitivity_table(scores, bench)
    sens.to_csv(out_dir / "encoding_sensitivity.csv", index=False)
    print()
    print(sens.to_string(index=False))
    print(f"\nwrote {out_dir / 'variant_encoding_audit.csv'} and "
          f"{out_dir / 'encoding_sensitivity.csv'}")


if __name__ == "__main__":
    main()
