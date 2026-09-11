#!/usr/bin/env python3
"""Recompute PerturbNet's five-setting, distinct-variant task-matched null.

The held-out target labels are permuted within each gene and setting while the
full training-plus-held-out candidate pool is retained.  The aggregation order
matches the manuscript estimate: mean over 15 cell-subsample seeds, mean over
settings for each distinct variant, then mean over distinct variants.  Outputs
are written only beneath a newly created, explicitly supplied directory.

Example
-------
python scripts/analysis/perturbnet_task_matched_null.py \\
    --repo . --base /path/to/PertResolve_workspace \\
    --preds /path/to/PertResolve_workspace/preds5/PerturbNet.npz \\
    --out /tmp/perturbnet_task_matched_null --n-perm 1000000 --seed 0
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

from pertresolve.paths import reject_repo_results


SPLITS = ("split1", "split2", "split3", "split5", "split6")
N_SEED = 15
N_SUB = 300
TOL = 1e-12
MC_CHUNK = 5_000


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--repo", required=True, type=Path)
    ap.add_argument("--base", required=True, type=Path)
    ap.add_argument("--preds", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--n-perm", type=int, default=1_000_000)
    ap.add_argument("--seed", type=int, default=0)
    return ap.parse_args()


def norm(matrix: np.ndarray) -> np.ndarray:
    return matrix / (np.linalg.norm(matrix, axis=1, keepdims=True) + 1e-12)


def score_matrix(distances: np.ndarray, target_indices: list[int]) -> np.ndarray:
    """Tie-aware PDS of every row against every requested target column."""
    n_candidates = distances.shape[1]
    out = np.empty((distances.shape[0], len(target_indices)), dtype=np.float64)
    for j, candidate_index in enumerate(target_indices):
        target_distance = distances[:, candidate_index]
        less = (distances < target_distance[:, None] - TOL).sum(axis=1)
        equal = (np.abs(distances - target_distance[:, None]) <= TOL).sum(axis=1)
        out[:, j] = 1.0 - (less + (equal - 1) / 2.0) / (n_candidates - 1)
    return out


def main() -> None:
    args = parse_args()
    args.out = reject_repo_results(args.out)
    args.out.mkdir(parents=True, exist_ok=False)
    sys.path.insert(0, str(args.repo / "scripts" / "analysis"))
    import harness as H

    H.set_base(args.base)
    predictions = np.load(args.preds)
    cells: dict[tuple[str, str], dict] = {}
    start_time = time.time()

    for gene in H.GENES:
        expression, labels = H.load_gene(gene)
        labels = np.asarray(labels)
        for setting in SPLITS:
            train_variants, heldout_variants = H.split_vars(gene, setting)
            row_variants = [
                variant
                for variant in heldout_variants
                if f"{gene}__{setting}__{variant}" in predictions.files
            ]
            if not row_variants:
                continue

            # Match train_only_grid.py: refit the diagonal scale using training
            # conditions plus reference cells, then form all evaluation targets there.
            fit_mask = np.isin(labels, train_variants) | np.isin(labels, H.WT_TAGS)
            scale = expression[fit_mask].std(axis=0)
            scale[scale < 1e-8] = 1.0
            scaled_expression = expression / scale
            predicted = norm(
                np.stack(
                    [predictions[f"{gene}__{setting}__{v}"] for v in row_variants]
                )
            )

            accumulated = None
            candidate_variants_reference = None
            for subsample_seed in range(N_SEED):
                rng = np.random.RandomState(
                    1000 + 7 * subsample_seed + H._GENE_SEED[gene]
                )
                reference_indices = np.where(np.isin(labels, H.WT_TAGS))[0]
                if len(reference_indices) > N_SUB:
                    reference_indices = rng.choice(
                        reference_indices, N_SUB, replace=False
                    )
                reference_mean = scaled_expression[reference_indices].mean(axis=0)
                truth = {}
                for variant in np.unique(labels):
                    if variant in H.WT_TAGS:
                        continue
                    indices = np.where(labels == variant)[0]
                    if len(indices) < 5:
                        continue
                    if len(indices) > N_SUB:
                        indices = rng.choice(indices, N_SUB, replace=False)
                    truth[variant] = (
                        scaled_expression[indices].mean(axis=0) - reference_mean
                    ).astype(np.float32)

                candidate_variants = [
                    variant
                    for variant in train_variants + heldout_variants
                    if variant in truth
                ]
                if candidate_variants_reference is None:
                    candidate_variants_reference = candidate_variants
                else:
                    assert candidate_variants == candidate_variants_reference
                candidate_index = {
                    variant: index for index, variant in enumerate(candidate_variants)
                }
                assert all(variant in candidate_index for variant in row_variants)
                target = norm(np.stack([truth[v] for v in candidate_variants]))
                distances = 1.0 - predicted @ target.T
                scores = score_matrix(
                    distances, [candidate_index[v] for v in row_variants]
                )
                accumulated = scores if accumulated is None else accumulated + scores

            cells[(gene, setting)] = {
                "rows": row_variants,
                "scores": accumulated / N_SEED,
                "n_candidates": len(candidate_variants_reference),
            }
            print(
                gene,
                setting,
                len(row_variants),
                len(candidate_variants_reference),
                flush=True,
            )
        del expression, labels

    # Canonical point-estimate order: seed mean (already above), then setting mean
    # for each distinct (gene, variant), then the mean over distinct variants.
    observed_rows = []
    restricted_expectation_rows = []
    settings_per_variant: dict[tuple[str, str], int] = {}
    for (gene, setting), cell in cells.items():
        for row_index, variant in enumerate(cell["rows"]):
            settings_per_variant[(gene, variant)] = (
                settings_per_variant.get((gene, variant), 0) + 1
            )
            observed_rows.append(
                (gene, setting, variant, cell["scores"][row_index, row_index])
            )
            restricted_expectation_rows.append(
                (gene, setting, variant, cell["scores"][row_index].mean())
            )

    observed_frame = pd.DataFrame(
        observed_rows, columns=["gene", "setting", "variant", "pds"]
    )
    restricted_frame = pd.DataFrame(
        restricted_expectation_rows,
        columns=["gene", "setting", "variant", "pds"],
    )
    observed = float(
        observed_frame.groupby(["gene", "variant"]).pds.mean().mean()
    )
    exact_restricted_mean = float(
        restricted_frame.groupby(["gene", "variant"]).pds.mean().mean()
    )
    n_variants = len(settings_per_variant)
    assert n_variants == 303
    committed_table = pd.read_csv(
        args.repo / "results" / "canonical" / "definitive_summary.csv"
    )
    committed_row = committed_table.loc[
        committed_table["method"].eq("PerturbNet")
        & committed_table["standardization"].eq("train-only")
    ]
    if len(committed_row) != 1:
        raise RuntimeError("expected exactly one canonical PerturbNet summary row")
    committed_observed = float(committed_row.iloc[0]["PDS"])

    diagnostic_rows = []
    for (gene, setting), cell in sorted(cells.items()):
        diagnostic_rows.append(
            {
                "gene": gene,
                "setting": setting,
                "n_heldout": len(cell["rows"]),
                "n_candidates": cell["n_candidates"],
                "observed_diag": float(np.diag(cell["scores"]).mean()),
                "restricted_exact_null": float(cell["scores"].mean()),
            }
        )
    pd.DataFrame(diagnostic_rows).to_csv(
        args.out / "gene_setting_diagnostics.csv", index=False
    )

    combined = observed_frame.rename(columns={"pds": "observed"}).merge(
        restricted_frame.rename(columns={"pds": "restricted_null"}),
        on=["gene", "setting", "variant"],
        validate="one_to_one",
    )
    gene_diagnostics = (
        combined.groupby(["gene", "variant"])[["observed", "restricted_null"]]
        .mean()
        .groupby("gene")
        .agg(
            observed=("observed", "mean"),
            restricted_null=("restricted_null", "mean"),
            n=("observed", "size"),
        )
        .reset_index()
    )
    gene_diagnostics.to_csv(args.out / "gene_aggregate_diagnostics.csv", index=False)

    # Each permutation is a bijection of the held-out labels within a gene×setting.
    # Because the score matrices are already seed-averaged, the same label mapping is
    # necessarily used for all 15 subsampling seeds. The weights implement the two
    # remaining canonical averaging stages without changing the original row identity.
    generator = np.random.default_rng(args.seed)
    null = np.empty(args.n_perm, dtype=np.float64)
    ordered_keys = sorted(cells)
    for chunk_start in range(0, args.n_perm, MC_CHUNK):
        chunk_size = min(MC_CHUNK, args.n_perm - chunk_start)
        values = np.zeros(chunk_size, dtype=np.float64)
        for gene, setting in ordered_keys:
            cell = cells[(gene, setting)]
            rows = cell["rows"]
            n_heldout = len(rows)
            permutations = np.argsort(
                generator.random((chunk_size, n_heldout)), axis=1
            )
            selected = cell["scores"][
                np.arange(n_heldout)[None, :], permutations
            ]
            weights = np.array(
                [
                    1.0 / (n_variants * settings_per_variant[(gene, variant)])
                    for variant in rows
                ]
            )
            values += selected @ weights
        null[chunk_start : chunk_start + chunk_size] = values

    null_mean = float(null.mean())
    null_sd = float(null.std(ddof=1))
    null_low, null_high = np.quantile(null, [0.025, 0.975])
    exceedances = int((null >= observed).sum())
    p_value = (exceedances + 1) / (args.n_perm + 1)
    null_mean_se = null_sd / np.sqrt(args.n_perm)
    p_se = np.sqrt(p_value * (1.0 - p_value) / (args.n_perm + 1))
    result = {
        "scheme": (
            "within each gene×setting, bijectively permute held-out target labels; "
            "retain full train+test candidate columns; use one mapping across 15 "
            "subsampling seeds; aggregate seed mean -> same (gene,variant) across "
            "settings -> 303-variant mean"
        ),
        "rng": f"numpy.default_rng({args.seed})",
        "n_permutations": args.n_perm,
        "observed_fresh": observed,
        "committed_observed": committed_observed,
        "observed_minus_committed": observed - committed_observed,
        "n_distinct_variants": n_variants,
        "n_gene_settings": len(cells),
        "restricted_heldout_null_exact_mean": exact_restricted_mean,
        "restricted_heldout_null_mc_mean": null_mean,
        "restricted_heldout_null_mc_sd": null_sd,
        "restricted_heldout_null_mc_mean_se": null_mean_se,
        "restricted_heldout_null_mc_mean_95_mc_ci": [
            null_mean - 1.959963984540054 * null_mean_se,
            null_mean + 1.959963984540054 * null_mean_se,
        ],
        "restricted_heldout_null_distribution_central_95_interval": [
            float(null_low),
            float(null_high),
        ],
        "one_sided_exceedances": exceedances,
        "one_sided_monte_carlo_p_plus1": p_value,
        "one_sided_monte_carlo_p_se": p_se,
        "full_pool_uniform_expected_pds": 0.5,
        "restricted_minus_full_pool_expectation": exact_restricted_mean - 0.5,
        "seconds": time.time() - start_time,
    }
    (args.out / "permutation_summary.json").write_text(
        json.dumps(result, indent=2) + "\n"
    )
    np.save(args.out / "permutation_null.npy", null)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
