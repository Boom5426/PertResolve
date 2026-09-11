#!/usr/bin/env python3
"""Variant-node bootstrap for Fig. 3e's fixed shared resolved-pair cohort.

This intentionally imports the committed analysis implementation so that cell
splits, cross-fitting, cosine forced choice, prediction refusal handling, and
split pooling are identical to resolved_pair_accuracy.py.  It writes only the
requested CSV path under /tmp; a checked result can then be promoted explicitly.

Example
-------
python scripts/analysis/resolved_pair_node_bootstrap.py \\
    --repo . --base /path/to/PertResolve_workspace \\
    --preds /path/to/PertResolve_workspace/preds5 \\
    --expected results/candidate_and_power/resolved_pair_accuracy.csv \\
    --out /tmp/resolved_pair_node_bootstrap.csv \\
    --genes GATA1,JAK1 --bootstrap 1000000 --seed 20260910 --chunk 5000
"""
from __future__ import annotations

import argparse
import importlib.util
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

from pertresolve.paths import reject_repo_results


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--repo", required=True)
    p.add_argument("--base", required=True)
    p.add_argument("--preds", required=True)
    p.add_argument("--expected", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--genes", default="GATA1,JAK1")
    p.add_argument("--bootstrap", type=int, default=100_000)
    p.add_argument("--seed", type=int, default=20_260_910)
    p.add_argument("--chunk", type=int, default=5_000)
    return p.parse_args()


def load_analysis(repo: Path):
    sys.path.insert(0, str(repo))
    script = repo / "scripts" / "analysis" / "resolved_pair_accuracy.py"
    spec = importlib.util.spec_from_file_location("resolved_pair_accuracy_audit", script)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def graph_components(n_nodes: int, edge_i: np.ndarray, edge_j: np.ndarray) -> int:
    parent = np.arange(n_nodes)

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for i, j in zip(edge_i, edge_j):
        ri, rj = find(int(i)), find(int(j))
        if ri != rj:
            parent[rj] = ri
    return len({find(i) for i in range(n_nodes)})


def node_bootstrap(y: np.ndarray, edge_i: np.ndarray, edge_j: np.ndarray,
                   n_nodes: int, B: int, seed: int, chunk: int):
    """Resample nodes; an observed edge gets multiplicity count_i * count_j.

    Rows of y are the fixed shared-cohort edges and columns are arms.  Duplicate
    copies of one original variant do not create a sibling edge with itself,
    because no such outcome exists in the estimand.  Any empty induced graph is
    redrawn (and counted), although none occurred in the reported run.
    """
    rng = np.random.default_rng(seed)
    probability = np.full(n_nodes, 1.0 / n_nodes)
    batches = []
    invalid = 0
    valid = 0
    while valid < B:
        size = min(chunk, B - valid)
        counts = rng.multinomial(n_nodes, probability, size=size)
        weights = counts[:, edge_i] * counts[:, edge_j]
        denom = weights.sum(axis=1)
        ok = denom > 0
        invalid += int((~ok).sum())
        if np.any(ok):
            batches.append((weights[ok] @ y) / denom[ok, None])
            valid += int(ok.sum())
    draws = np.concatenate(batches, axis=0)[:B]
    return draws, invalid


def reconstruct_gene(R, H, base: Path, preds_dir: Path, bench: dict, gene: str,
                     methods: dict):
    X, lab = H.load_gene(gene)
    lab = np.asarray(lab)
    per_seed = [R.halves(H, X, lab, bench[gene], s, gene) for s in range(R.NSEED)]
    variants = sorted(per_seed[0][0].keys())
    pairs = list(combinations(range(len(variants)), 2))
    pi = np.fromiter((p[0] for p in pairs), dtype=int)
    pj = np.fromiter((p[1] for p in pairs), dtype=int)

    _, _, mask_xfit = R.criterion(per_seed, variants, R.SELECT_SEEDS)
    res_xfit = mask_xfit[pi, pj]
    measurement = R.measurement_accuracy(per_seed, variants, R.EVAL_SEEDS)[pi, pj]
    meas_by_pair = {pair: float(measurement[k]) for k, pair in enumerate(pairs)}
    res_by_pair = {pair: bool(res_xfit[k]) for k, pair in enumerate(pairs)}

    per_method_pairs = {}
    first = methods[next(iter(methods))]
    splits = sorted({key.split("__")[1] for key in first.files
                     if key.startswith(f"{gene}__")})
    for method, predictions in methods.items():
        for split in splits:
            keyed = {v: f"{gene}__{split}__{v}" for v in variants
                     if f"{gene}__{split}__{v}" in predictions.files}
            usable = {}
            for variant, key in keyed.items():
                vector = np.asarray(predictions[key], dtype=float)
                if np.all(np.isfinite(vector)):
                    usable[variant] = vector
            if len(usable) < 2:
                continue
            names = [v for v in variants if v in usable]
            position = {v: k for k, v in enumerate(names)}
            keep = np.array([variants[i] in usable and variants[j] in usable
                             for i, j in pairs])
            ui = np.array([position[variants[i]] for i in pi[keep]])
            uj = np.array([position[variants[j]] for j in pj[keep]])
            score = R.model_accuracy(per_seed, names, usable, R.EVAL_SEEDS)[ui, uj]
            method_store = per_method_pairs.setdefault(method, {})
            for value, pair in zip(score, [pairs[k] for k in np.flatnonzero(keep)]):
                method_store.setdefault(pair, []).append(float(value))

    shared = set.intersection(*(set(store) for store in per_method_pairs.values()))
    cohort = sorted(pair for pair in shared if res_by_pair[pair])
    method_names = sorted(per_method_pairs)
    values = [np.array([meas_by_pair[pair] for pair in cohort], dtype=float)]
    for method in method_names:
        store = per_method_pairs[method]
        values.append(np.array([np.mean(store[pair]) for pair in cohort], dtype=float))
    y = np.column_stack(values)
    arm_names = ["split-half measurement", *method_names]

    nodes = sorted({node for pair in cohort for node in pair})
    node_position = {node: k for k, node in enumerate(nodes)}
    edge_i = np.array([node_position[pair[0]] for pair in cohort], dtype=int)
    edge_j = np.array([node_position[pair[1]] for pair in cohort], dtype=int)
    degree = np.bincount(np.r_[edge_i, edge_j], minlength=len(nodes))
    graph = {
        "n_nodes": len(nodes),
        "n_edges": len(cohort),
        "degree_min": int(degree.min()),
        "degree_median": float(np.median(degree)),
        "degree_max": int(degree.max()),
        "components": graph_components(len(nodes), edge_i, edge_j),
        "n_variants_gene": len(variants),
        "n_methods": len(method_names),
    }
    return arm_names, y, edge_i, edge_j, graph


def main():
    args = parse_args()
    repo = Path(args.repo).resolve()
    base = Path(args.base).resolve()
    preds_dir = Path(args.preds).resolve()
    out = reject_repo_results(Path(args.out))
    if "/tmp/" not in str(out) and str(out) != "/tmp":
        raise SystemExit("temporary audit refuses a non-/tmp output")

    R = load_analysis(repo)
    R.add_harness_to_path(base)
    import harness as H
    H.set_base(base)

    bench_df = pd.read_csv(base / "pertresolve_bench.csv")
    bench = {gene: set(bench_df[bench_df.gene == gene]["variant"])
             for gene in H.GENES}
    methods = {path.stem: np.load(path, allow_pickle=True)
               for path in sorted(preds_dir.glob("*.npz"))}
    expected = pd.read_csv(args.expected)
    expected = expected[(expected.split == "pooled")
                        & (expected.stratum == "resolvable_xfit_shared")]

    all_rows = []
    for gene in args.genes.split(","):
        names, y, edge_i, edge_j, graph = reconstruct_gene(
            R, H, base, preds_dir, bench, gene, methods)
        point = y.mean(axis=0)
        expected_gene = expected[expected.gene == gene].set_index("method")
        if set(names) != set(expected_gene.index):
            raise AssertionError((gene, sorted(set(names) ^ set(expected_gene.index))))
        deltas = [abs(point[k] - float(expected_gene.loc[name, "accuracy"]))
                  for k, name in enumerate(names)]
        if max(deltas) > 5.01e-6:
            raise AssertionError(f"{gene}: max point delta {max(deltas)}")
        if len(edge_i) != int(expected_gene.n_pairs.unique().item()):
            raise AssertionError(f"{gene}: cohort size mismatch")

        effective_seed = args.seed + H._GENE_SEED[gene]
        draws, invalid = node_bootstrap(
            y, edge_i, edge_j, graph["n_nodes"], args.bootstrap,
            effective_seed, args.chunk)
        lo = np.percentile(draws, 2.5, axis=0)
        hi = np.percentile(draws, 97.5, axis=0)
        gap_point = point - point[0]
        gap_draws = draws - draws[:, [0]]
        gap_lo = np.percentile(gap_draws, 2.5, axis=0)
        gap_hi = np.percentile(gap_draws, 97.5, axis=0)

        print("GRAPH", gene, graph, "bootstrap_seed", effective_seed,
              "B", args.bootstrap, "invalid_redrawn", invalid, flush=True)
        print("VERIFY", gene, "max_abs_delta_vs_5dp_committed", max(deltas), flush=True)
        for k, name in enumerate(names):
            all_rows.append({
                "gene": gene,
                "method": name,
                **graph,
                "accuracy": point[k],
                "node_ci_lo": lo[k],
                "node_ci_hi": hi[k],
                "gap_model_minus_measurement": gap_point[k],
                "gap_node_ci_lo": gap_lo[k],
                "gap_node_ci_hi": gap_hi[k],
                "bootstrap_B": args.bootstrap,
                "base_seed": args.seed,
                "effective_seed": effective_seed,
                "invalid_draws_redrawn": invalid,
            })
    result = pd.DataFrame(all_rows)
    result.to_csv(out, index=False, float_format="%.10f")
    print(f"WROTE {out} rows={len(result)}", flush=True)


if __name__ == "__main__":
    main()
