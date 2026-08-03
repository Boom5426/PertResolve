"""Create all data-direct derived tables used by the refined Fig. 6 panels.

No plotting occurs here. Fixed seeds make every interval and prediction
reproducible. The script also asserts that retrospective LODO AUROCs agree with
the repository's committed audit table.
"""
from __future__ import annotations

import glob
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.preprocessing import StandardScaler

from fig6_common import DERIVED, RAW

ALL_DATASETS = ["Replogle", "Norman", "Adamson", "GATA1", "JAK1", "TP53", "KRAS"]
EVALUABLE = ["Replogle", "Norman", "Adamson", "GATA1", "JAK1"]
BALANCED_EXTERNAL = ["Replogle", "Norman", "Adamson"]
ALLELES = ["TP53", "KRAS", "GATA1", "JAK1"]


def native_rankability(all_rows: pd.DataFrame, dataset: str) -> pd.DataFrame:
    sub = all_rows[
        (all_rows["dataset"] == dataset)
        & (all_rows["metric"] == "edist")
        & (all_rows["space"] == "pca")
    ].copy()
    return sub.loc[sub.groupby("perturbation")["n_work"].idxmax()].copy()


def bootstrap_auc(y, score, n_boot=2000, seed=0):
    y = np.asarray(y, dtype=int)
    score = np.asarray(score, dtype=float)
    def fast_auc(y_value, score_value):
        positive = np.asarray(score_value)[np.asarray(y_value) == 1]
        negative = np.sort(np.asarray(score_value)[np.asarray(y_value) == 0])
        lower = np.searchsorted(negative, positive, side="left")
        upper = np.searchsorted(negative, positive, side="right")
        return float(np.mean((lower + upper) / (2.0 * len(negative))))

    point = fast_auc(y, score)
    rng = np.random.default_rng(seed)
    estimates = []
    for _ in range(n_boot):
        index = rng.integers(0, len(y), len(y))
        if np.unique(y[index]).size == 2:
            estimates.append(fast_auc(y[index], score[index]))
    low, high = np.quantile(estimates, [0.025, 0.975])
    return point, float(low), float(high), len(estimates)


def bootstrap_depthmatch(paired, n_boot=2000, seed=0):
    native = paired["native_unrankable"].to_numpy(float)
    matched = paired["matched50_unrankable"].to_numpy(float)
    rng = np.random.default_rng(seed)
    b_native, b_matched, b_delta = [], [], []
    for _ in range(n_boot):
        index = rng.integers(0, len(paired), len(paired))
        n_value = native[index].mean() * 100.0
        m_value = matched[index].mean() * 100.0
        b_native.append(n_value)
        b_matched.append(m_value)
        b_delta.append(m_value - n_value)
    def q(values):
        return np.quantile(values, [0.025, 0.975]).astype(float)
    return {
        "native_pct": native.mean() * 100.0,
        "native_lo": q(b_native)[0],
        "native_hi": q(b_native)[1],
        "matched50_pct": matched.mean() * 100.0,
        "matched50_lo": q(b_matched)[0],
        "matched50_hi": q(b_matched)[1],
        "delta_pp": (matched - native).mean() * 100.0,
        "delta_lo": q(b_delta)[0],
        "delta_hi": q(b_delta)[1],
    }


def wilson_interval(k: int, n: int, z: float = 1.96):
    if n == 0:
        return np.nan, np.nan
    p = k / n
    den = 1.0 + z * z / n
    centre = (p + z * z / (2.0 * n)) / den
    half = z * np.sqrt(p * (1.0 - p) / n + z * z / (4.0 * n * n)) / den
    return max(0.0, centre - half), min(1.0, centre + half)


def retrospective(second_probe: pd.DataFrame):
    frames = {}
    for dataset in ALL_DATASETS:
        frame = native_rankability(second_probe, dataset)
        frame["dataset"] = dataset
        frame["log_effect"] = np.log10(frame["effect_size"].clip(lower=1e-6))
        frame["y"] = frame["rankable"].astype(int)
        frames[dataset] = frame[["dataset", "perturbation", "log_effect", "y"]]

    predictions, curves = [], []
    for held_out in EVALUABLE:
        test = frames[held_out].copy()
        train = pd.concat(
            [frames[name] for name in ALL_DATASETS if name != held_out],
            ignore_index=True,
        )
        scaler = StandardScaler().fit(train[["log_effect"]])
        model = LogisticRegression(max_iter=2000, random_state=42).fit(
            scaler.transform(train[["log_effect"]]), train["y"]
        )
        test["score"] = model.predict_proba(
            scaler.transform(test[["log_effect"]])
        )[:, 1]
        test["n_pos"] = int(test["y"].sum())
        test["n_neg"] = int((1 - test["y"]).sum())
        predictions.append(test)
        fpr, tpr, thresholds = roc_curve(test["y"], test["score"])
        curves.append(
            pd.DataFrame(
                {
                    "dataset": held_out,
                    "fpr": fpr,
                    "tpr": tpr,
                    "threshold": thresholds,
                }
            )
        )

    prediction_table = pd.concat(predictions, ignore_index=True)
    curve_table = pd.concat(curves, ignore_index=True)
    prediction_table.to_csv(DERIVED / "retrospective_predictions.csv", index=False)
    curve_table.to_csv(DERIVED / "retrospective_roc.csv", index=False)

    committed = pd.read_csv(RAW / "rankability_predictor_honest.csv")
    committed = committed[committed["feature_set"] == "effect_size"].copy()
    for dataset in EVALUABLE:
        expected = float(committed.loc[committed["held_out"] == dataset, "auroc"].iloc[0])
        observed = roc_auc_score(
            prediction_table.loc[prediction_table["dataset"] == dataset, "y"],
            prediction_table.loc[prediction_table["dataset"] == dataset, "score"],
        )
        assert abs(expected - observed) < 0.01, (dataset, expected, observed)


def prospective():
    paths = sorted(glob.glob(str(RAW / "pilot_validation" / "*_pilot.csv")))
    pilot = pd.concat([pd.read_csv(path) for path in paths], ignore_index=True)
    pilot.to_csv(DERIVED / "pilot_combined.csv", index=False)

    rows, predictions = [], []
    for cells in (25, 50):
        feature = f"pilot_eff_{cells}"
        data = pilot.dropna(subset=[feature, "rank_T100"]).copy()
        data["rank_T100"] = data["rank_T100"].astype(int)
        for held_out in BALANCED_EXTERNAL:
            train = data[data["dataset"] != held_out]
            test = data[data["dataset"] == held_out].copy()
            scaler = StandardScaler().fit(train[[feature]])
            model = LogisticRegression(max_iter=2000, random_state=42).fit(
                scaler.transform(train[[feature]]), train["rank_T100"]
            )
            test["score"] = model.predict_proba(
                scaler.transform(test[[feature]])
            )[:, 1]
            test["pilot_cells"] = cells
            test["method"] = "Learned effect"
            predictions.append(
                test[["dataset", "pert", "rank_T100", "score", "pilot_cells", "method"]]
            )
            auc, low, high, n_valid = bootstrap_auc(
                test["rank_T100"], test["score"], seed=6100 + cells + len(rows)
            )
            rows.append(
                {
                    "dataset": held_out,
                    "pilot_cells": cells,
                    "method": "Learned effect",
                    "auroc": auc,
                    "lo": low,
                    "hi": high,
                    "n": len(test),
                    "n_pos": int(test["rank_T100"].sum()),
                    "n_neg": int((1 - test["rank_T100"]).sum()),
                    "valid_bootstraps": n_valid,
                }
            )

    for cells in (25, 50):
        feature = f"pilot_snr_{cells}"
        data = pilot.dropna(subset=[feature, "rank_T100"]).copy()
        data["rank_T100"] = data["rank_T100"].astype(int)
        for held_out in BALANCED_EXTERNAL:
            test = data[data["dataset"] == held_out].copy()
            auc, low, high, n_valid = bootstrap_auc(
                test["rank_T100"], test[feature], seed=6200 + cells + len(rows)
            )
            rows.append(
                {
                    "dataset": held_out,
                    "pilot_cells": cells,
                    "method": "Training-free SNR",
                    "auroc": auc,
                    "lo": low,
                    "hi": high,
                    "n": len(test),
                    "n_pos": int(test["rank_T100"].sum()),
                    "n_neg": int((1 - test["rank_T100"]).sum()),
                    "valid_bootstraps": n_valid,
                }
            )

    summary = pd.DataFrame(rows)
    summary.to_csv(DERIVED / "prospective_auc.csv", index=False)
    prediction_table = pd.concat(predictions, ignore_index=True)
    prediction_table.to_csv(DERIVED / "prospective_predictions.csv", index=False)

    # Calibration is intentionally restricted to 50-cell learned predictions
    # from the three balanced external datasets used in the prospective headline.
    calibration_source = prediction_table[
        (prediction_table["pilot_cells"] == 50)
        & (prediction_table["dataset"].isin(BALANCED_EXTERNAL))
    ].copy()
    calibration_source["bin"] = pd.qcut(
        calibration_source["score"],
        q=5,
        labels=False,
        duplicates="drop",
    )
    calibration_rows = []
    for bin_id, group in calibration_source.groupby("bin", sort=True):
        n = len(group)
        positives = int(group["rank_T100"].sum())
        low, high = wilson_interval(positives, n)
        calibration_rows.append(
            {
                "bin": int(bin_id),
                "mean_predicted": group["score"].mean(),
                "observed_rate": group["rank_T100"].mean(),
                "lo": low,
                "hi": high,
                "n": n,
            }
        )
    pd.DataFrame(calibration_rows).to_csv(
        DERIVED / "prospective_calibration.csv", index=False
    )


def depthmatch_and_alleles(second_probe: pd.DataFrame):
    paired_tables, summaries = [], []
    for dataset_index, dataset in enumerate(["Replogle", "Norman", "Adamson"]):
        base = second_probe[
            (second_probe["dataset"] == dataset)
            & (second_probe["metric"] == "edist")
            & (second_probe["space"] == "pca")
        ].copy()
        native = base.loc[base.groupby("perturbation")["n_work"].idxmax()].copy()
        matched = base[base["n_work"] == 50].copy()
        native = native[
            ["perturbation", "rankable", "n_work", "n_cells"]
        ].rename(
            columns={
                "rankable": "native_rankable",
                "n_work": "native_n_work",
                "n_cells": "native_n_cells",
            }
        )
        matched = matched[
            ["perturbation", "rankable", "n_work", "n_cells"]
        ].rename(
            columns={
                "rankable": "matched50_rankable",
                "n_work": "matched50_n_work",
                "n_cells": "matched50_n_cells",
            }
        )
        paired = native.merge(matched, on="perturbation", how="inner")
        paired["dataset"] = dataset
        paired["native_unrankable"] = (~paired["native_rankable"]).astype(int)
        paired["matched50_unrankable"] = (~paired["matched50_rankable"]).astype(int)
        paired_tables.append(paired)
        values = bootstrap_depthmatch(paired, seed=6300 + dataset_index)
        values.update({"dataset": dataset, "paired_n": len(paired)})
        summaries.append(values)
    pd.concat(paired_tables, ignore_index=True).to_csv(
        DERIVED / "depthmatch_paired.csv", index=False
    )
    pd.DataFrame(summaries).to_csv(
        DERIVED / "depthmatch_summary.csv", index=False
    )

    allele_tables = []
    for gene in ALLELES:
        table = native_rankability(second_probe, gene).copy()
        table["gene"] = gene
        allele_tables.append(
            table[["gene", "perturbation", "effect_size", "n_cells", "n_work", "rankable"]]
        )
    pd.concat(allele_tables, ignore_index=True).to_csv(
        DERIVED / "allele_native_rankability.csv", index=False
    )


def main():
    DERIVED.mkdir(parents=True, exist_ok=True)
    second_probe = pd.read_csv(RAW / "second_probe_rankability_table.csv")
    retrospective(second_probe)
    prospective()
    depthmatch_and_alleles(second_probe)
    print("Prepared Fig. 6 derived tables in", DERIVED)


if __name__ == "__main__":
    main()
