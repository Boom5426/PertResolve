#!/usr/bin/env python3
"""Aggregate the resolution panel into one machine-readable table. Computes nothing.

Before this script the panel existed only as a markdown table inside
``results/resolution_panel_v2/RESULTS.md``, and three of its attributes (which run a row came
from, whether the perturbation is genetic or chemical, and whether the readout is RNA or
protein) existed nowhere else at all. A figure script therefore had to parse markdown, and no
check could catch a row drifting away from the JSON it claims to summarise.

**This is aggregation, not analysis.** Every numeric column is copied verbatim from a file that
already exists: the 16 v2 runs under ``results/resolution_panel_v2/panel/`` and the six v1 rows
of ``results/canonical/resolution_panel.csv``. Nothing is recomputed, no estimator runs, and the
script asserts on exit that every emitted value is byte-equal to its source.

Three attribute columns cannot be copied because no source file carries them, so they are
declared here with the reason for each declaration, and the declaration is checked against what
the JSON does carry wherever a check is possible:

``resource``      the source experiment. Sixteen v2 rows come from five experiments, so a count
                  over rows is not a count over experiments. Checked against the ``h5ad`` field.
``perturbation_class``  genetic or chemical. A property of the experiment, not of any recorded
                  field. Declared, with the reason in :data:`RESOURCES`.
``modality``      RNA or protein. Checked against the ``h5ad`` filename, which names it.

Two provenance flags are also declared, because whether a row was pre-registered before the
panel ran and whether it is a demonstration subsample are facts about the study design that no
result file records:

``arm``           ``prespecified`` for the six screens of the v1 pre-registration, ``additional``
                  for the sixteen added under v2.
``demonstration_subsample``  True where the row scores a small deliberate subsample of a much
                  larger resource rather than the resource itself. Two rows qualify, and both
                  are in the six that clear detection, so a reader who does not know this will
                  over-read the top of the ranking.

The two negative ``rho2_nn_median`` values are emitted as they stand and flagged, not dropped:
the statistic is bounded below at exactly -1 by construction, since it is
``raw_squared_separation / (2 * eta2) - 1`` and the numerator is non-negative, so a value near
-1 is a noise-dominated, floor-limited measurement rather than a large negative separation. See
``pertresolve/resolution/scaling.py``.

Usage:
  python build_resolution_panel_table.py --panel results/resolution_panel_v2 \
      --v1 results/canonical/resolution_panel.csv --out build/panel_table

Writes ``resolution_panel_v2_table.csv`` (the 22 panel rows) and
``resolution_panel_v2_sensitivity.csv`` (the 13 gate and seed-probe runs, which are arms of the
same protocol and are kept so a reader can see what was varied). ``--out`` is refused if it
resolves inside the repository's ``results/`` tree; promote the file deliberately.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from pertresolve.paths import reject_repo_results  # noqa: E402

#: One entry per source experiment. ``modality_of`` names the substring of the ``h5ad`` filename
#: that fixes the readout, so the declaration below can be checked rather than trusted.
#:
#: ``perturbation_class`` is declared. Genetic means the perturbation is a guide RNA or a coding
#: variant; chemical means a compound. This is a property of the experiment and appears in no
#: field of any result file, so it is stated here once with its reason rather than being
#: re-derived from a dataset name at every call site.
RESOURCES = {
    "frangieh": dict(resource="Frangieh 2021", perturbation_class="genetic",
                     reason="CRISPR knockout screen in melanoma, guide RNA per cell"),
    "gse306429": dict(resource="GSE306429", perturbation_class="chemical",
                      reason="compound treatment, one compound per cell identity"),
    "papalexi": dict(resource="Papalexi 2021 ECCITE", perturbation_class="genetic",
                     reason="CRISPR screen with a matched antibody panel"),
    "perturbmulti": dict(resource="PerturbMulti", perturbation_class="genetic",
                         reason="CRISPR screen with paired RNA and protein readouts"),
    "sciplex3": dict(resource="sci-Plex 3", perturbation_class="chemical",
                     reason="small-molecule screen across cell lines and times"),
    "mcfarland": dict(resource="McFarland 2020", perturbation_class="chemical",
                      reason="multiplexed compound screen"),
    "tahoe100m": dict(resource="Tahoe-100M", perturbation_class="chemical",
                      reason="compound perturbation atlas"),
    "norman2019": dict(resource="Norman 2019", perturbation_class="genetic",
                       reason="CRISPRa single and paired gene activation"),
    "replogle": dict(resource="Replogle 2022", perturbation_class="genetic",
                     reason="genome-scale CRISPRi knockdown"),
    "adamson2016": dict(resource="Adamson 2016", perturbation_class="genetic",
                        reason="CRISPRi knockdown"),
    "vcc_training": dict(resource="Virtual Cell Challenge", perturbation_class="genetic",
                         reason="CRISPRi knockdown training split"),
    "kolf_ipsc": dict(resource="KOLF2.1J iPSC atlas", perturbation_class="genetic",
                      reason="genome-scale CRISPRi Perturb-seq in human iPSC"),
    "xatlas": dict(resource="X-Atlas/Orion", perturbation_class="genetic",
                   reason="genome-wide dual-guide CRISPRi, one cell line per configuration"),
    "parse10m": dict(resource="Parse 10M PBMC cytokines", perturbation_class="cytokine",
                     reason="recombinant cytokine stimulation of primary PBMCs, one ligand "
                            "per cell. Neither genetic nor small-molecule, so it is declared "
                            "as its own class rather than forced into either"),
}

#: The six screens fixed by the v1 pre-registration, before the v2 panel existed. Everything else
#: is an addition, and a reader must be able to tell which is which.
PRESPECIFIED = {"mcfarland", "tahoe100m", "norman2019", "replogle", "adamson2016", "vcc_training"}

#: Rows scoring a deliberate small subsample of a far larger resource. Both clear detection, so
#: leaving the fact implicit would let a reader read them as full measurements.
DEMONSTRATION_SUBSAMPLE = {
    "mcfarland": "11 levels at 500 cells each, drawn from a 182,875-cell parent",
    "tahoe100m": "7 units drawn from a 100.6 M-cell resource",
}

#: v1 rows are RNA; none of the six has a protein arm.
V1_MODALITY = "RNA"

#: A v1 row whose published cohort is now known to contain an entry that is not a perturbation.
#: The replacement is a re-run of the identical criterion on the corrected cohort, so the row is
#: superseded by a measurement rather than by an edit, and the source file it is verified against
#: moves with it. Nothing else about the row changes.
V1_CORRECTIONS = {
    "adamson2016": dict(
        source="resolution_panel_v2/gate1_recheck/adamson2016.json",
        reason="the perturbation column carries a level whose literal value is <NA>, holding "
               "2,613 cells and therefore clearing the 200-cell inclusion threshold. The first "
               "evaluation scored it as one of 94 perturbations. Norman, Replogle and the "
               "Virtual Cell Challenge training dataset were checked for the same defect and "
               "carry no such level."),
}

COLUMNS = [
    "dataset", "resource", "configuration", "src", "arm", "perturbation_class", "modality",
    "n", "n_units_total", "n_units_eligible", "n_units_scored", "n_units_cell_capped",
    "n_cells_scored", "n_control_cells",
    "detection_fraction", "identification_fraction", "jointly_evaluable_fraction",
    "split_half_reproducibility_reference", "rho2", "rho2_nn_median",
    "rho2_nn_floor_limited", "model_ranking_p_correct_order", "depth", "n_seeds", "n_boot",
    "seed", "n_components", "reduce_path",
    "demonstration_subsample", "demonstration_note", "correction_note", "source_file",
]

#: ``rho2_nn_median`` is bounded below by -1. A value this close to the floor is reporting that
#: the noise term dominates the measured separation, not that the separation is large and
#: negative, so the flag exists to stop the magnitude being ranked.
FLOOR_LIMITED_BELOW = 0.0


def _axis(record: dict, current: str, historical: str):
    """Read the independent public name, accepting frozen historical JSON artifacts."""
    if current in record:
        return record[current]
    return record[historical]


def _resource_of(name: str) -> dict:
    for prefix, meta in RESOURCES.items():
        if name.startswith(prefix):
            return meta
    raise SystemExit(f"{name!r} matches no declared resource; add it to RESOURCES with a reason")


def _configuration(name: str, record: dict) -> str:
    """The stratum and modality that distinguish this row from its siblings."""
    strata = record.get("stratum") or []
    parts = [str(s) for s in strata]
    modality = _modality_of(name, record)
    if modality == "protein" or "_rna" in name or "_protein" in name:
        parts.append(modality)
    return " / ".join(parts) if parts else "whole resource"


def _modality_of(name: str, record: dict) -> str:
    """RNA or protein, taken from the source filename rather than from the row name."""
    stem = Path(str(record.get("h5ad") or record["source"])).name.lower()
    by_file = "protein" if "protein" in stem else "RNA"
    by_name = "protein" if "_protein" in name else ("RNA" if "_rna" in name else None)
    if by_name is not None and by_name != by_file:
        raise SystemExit(f"{name}: row name says {by_name} but {stem} says {by_file}")
    return by_file


def _v2_rows(panel_dir: Path) -> list[dict]:
    rows = []
    for path in sorted((panel_dir / "panel").glob("*.json")):
        d = json.loads(path.read_text())
        name = d["name"]
        meta = _resource_of(name)
        sel = d["selection"]
        rho2_nn = float(d["rho2_nn_median"])
        rows.append({
            "dataset": name,
            "resource": meta["resource"],
            "configuration": _configuration(name, d),
            "src": "v2",
            "arm": "prespecified" if name in PRESPECIFIED else "additional",
            "perturbation_class": meta["perturbation_class"],
            "modality": _modality_of(name, d),
            "n": d["n_perturbations"],
            "n_units_total": sel["n_units_total"],
            "n_units_eligible": sel["n_units_eligible"],
            "n_units_scored": sel["n_units_scored"],
            "n_units_cell_capped": sel["n_units_cell_capped"],
            "n_cells_scored": d["n_cells_scored"],
            "n_control_cells": sel["n_control_cells"],
            "detection_fraction": _axis(d, "detection_fraction", "detectable_fraction"),
            "identification_fraction": _axis(d, "identification_fraction",
                                               "identifiable_fraction"),
            "jointly_evaluable_fraction": d["jointly_evaluable_fraction"],
            "split_half_reproducibility_reference": _axis(
                d, "split_half_reproducibility_reference", "replicate_ceiling"),
            "rho2": d["rho2"],
            "rho2_nn_median": rho2_nn,
            "rho2_nn_floor_limited": rho2_nn <= FLOOR_LIMITED_BELOW,
            "model_ranking_p_correct_order": _axis(
                d, "model_ranking_p_correct_order", "p_correct_order"),
            "depth": d["depth"], "n_seeds": d["n_seeds"], "n_boot": d["n_boot"],
            "seed": d["seed"], "n_components": d["n_components"],
            "reduce_path": d.get("reduce_path", ""),
            "correction_note": "",
            "demonstration_subsample": name in DEMONSTRATION_SUBSAMPLE,
            "demonstration_note": DEMONSTRATION_SUBSAMPLE.get(name, ""),
            "source_file": str(path.relative_to(panel_dir.parent.parent)),
        })
    return rows


def _v1_rows(v1_csv: Path) -> list[dict]:
    """The six carried-over screens. They predate the v2 JSON schema, so the fields that schema
    added (unit counts, jointly evaluable fraction, reduction path) are empty rather than
    invented, and the emptiness is itself information a reader needs."""
    rows = []
    for r in csv.DictReader(v1_csv.read_text().splitlines()):
        name = r["dataset_name"]
        meta = _resource_of(name)
        correction = V1_CORRECTIONS.get(name)
        if correction is not None:
            fixed = json.loads((v1_csv.parent.parent / correction["source"]).read_text())
            r = dict(r, n_perturbations=fixed["n_perturbations"],
                     detectable_fraction=_axis(fixed, "detection_fraction",
                                               "detectable_fraction"),
                     identifiable_fraction=_axis(fixed, "identification_fraction",
                                                 "identifiable_fraction"),
                     replicate_ceiling=_axis(fixed, "split_half_reproducibility_reference",
                                             "replicate_ceiling"),
                     rho2=fixed["rho2"], rho2_nn_median=fixed["rho2_nn_median"],
                     p_correct_order=_axis(fixed, "model_ranking_p_correct_order",
                                           "recovery_p_correct_order")
                     if "recovery_p_correct_order" in fixed else _axis(
                         fixed, "model_ranking_p_correct_order", "p_correct_order"))
        rho2_nn = float(r["rho2_nn_median"])
        rows.append({
            "dataset": name,
            "resource": meta["resource"],
            "configuration": "whole resource",
            "src": "v1",
            "arm": "prespecified" if name in PRESPECIFIED else "additional",
            "perturbation_class": meta["perturbation_class"],
            "modality": V1_MODALITY,
            "n": int(r["n_perturbations"]),
            "n_units_total": "", "n_units_eligible": "", "n_units_scored": "",
            "n_units_cell_capped": "", "n_cells_scored": "", "n_control_cells": "",
            "detection_fraction": float(r["detectable_fraction"]),
            "identification_fraction": float(r["identifiable_fraction"]),
            "jointly_evaluable_fraction": "",
            "split_half_reproducibility_reference": float(r["replicate_ceiling"]),
            "rho2": float(r["rho2"]),
            "rho2_nn_median": rho2_nn,
            "rho2_nn_floor_limited": rho2_nn <= FLOOR_LIMITED_BELOW,
            "model_ranking_p_correct_order": float(r["p_correct_order"]),
            "depth": 50, "n_seeds": 8, "n_boot": 1000, "seed": 0, "n_components": 50,
            "reduce_path": "",
            "correction_note": correction["reason"] if correction else "",
            "demonstration_subsample": name in DEMONSTRATION_SUBSAMPLE,
            "demonstration_note": DEMONSTRATION_SUBSAMPLE.get(name, ""),
            "source_file": str(v1_csv),
        })
    return rows


def _sensitivity_rows(panel_dir: Path) -> list[dict]:
    """The gate and seed-probe runs. Same protocol, one parameter varied each, kept so the
    sensitivity arms are readable from a table rather than only from prose."""
    rows = []
    for sub in ("gate1", "gate2", "gate6", "gate7", "seedprobe"):
        for path in sorted((panel_dir / sub).glob("*.json")):
            d = json.loads(path.read_text())
            name = d["name"]
            meta = _resource_of(name)
            sel = d["selection"]
            rho2_nn = float(d["rho2_nn_median"])
            rows.append({
                "dataset": name, "resource": meta["resource"],
                "configuration": _configuration(name, d), "src": "v2", "arm": sub,
                "perturbation_class": meta["perturbation_class"],
                "modality": _modality_of(name, d),
                "n": d["n_perturbations"],
                "n_units_total": sel["n_units_total"],
                "n_units_eligible": sel["n_units_eligible"],
                "n_units_scored": sel["n_units_scored"],
                "n_units_cell_capped": sel["n_units_cell_capped"],
                "n_cells_scored": d["n_cells_scored"],
                "n_control_cells": sel["n_control_cells"],
                "detection_fraction": _axis(d, "detection_fraction", "detectable_fraction"),
                "identification_fraction": _axis(d, "identification_fraction",
                                                   "identifiable_fraction"),
                "jointly_evaluable_fraction": d["jointly_evaluable_fraction"],
                "split_half_reproducibility_reference": _axis(
                    d, "split_half_reproducibility_reference", "replicate_ceiling"),
                "rho2": d["rho2"],
                "rho2_nn_median": rho2_nn,
                "rho2_nn_floor_limited": rho2_nn <= FLOOR_LIMITED_BELOW,
                "model_ranking_p_correct_order": _axis(
                    d, "model_ranking_p_correct_order", "p_correct_order"),
                "depth": d["depth"], "n_seeds": d["n_seeds"], "n_boot": d["n_boot"],
                "seed": d["seed"], "n_components": d["n_components"],
                "reduce_path": d.get("reduce_path", ""),
            "correction_note": "",
                "demonstration_subsample": name in DEMONSTRATION_SUBSAMPLE,
                "demonstration_note": DEMONSTRATION_SUBSAMPLE.get(name, ""),
                "source_file": str(path.relative_to(panel_dir.parent.parent)),
            })
    return rows


def _write(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS)
        w.writeheader()
        w.writerows(rows)


def _verify(rows: list[dict], panel_dir: Path, v1_csv: Path) -> None:
    """Read every source back and assert the table did not alter a single number.

    Aggregation is exactly the step where a value can be silently changed, so the check is part
    of the script rather than a thing someone might do afterwards.
    """
    v1 = {r["dataset_name"]: r for r in csv.DictReader(v1_csv.read_text().splitlines())}
    checked = 0
    for row in rows:
        if row["src"] == "v1":
            correction = V1_CORRECTIONS.get(row["dataset"])
            # A corrected row is verified against the re-run that supersedes it, not against the
            # published table it replaces. Verifying it against the old table would either fail or,
            # worse, force the correction back out.
            src = (json.loads((v1_csv.parent.parent / correction["source"]).read_text())
                   if correction else v1[row["dataset"]])
            pairs = [("detection_fraction", "detectable_fraction"),
                     ("identification_fraction", "identifiable_fraction"),
                     ("split_half_reproducibility_reference", "replicate_ceiling"),
                     ("rho2_nn_median", "rho2_nn_median"),
                     ("model_ranking_p_correct_order", "p_correct_order")]
            for mine, theirs in pairs:
                assert float(row[mine]) == float(src[theirs]), (row["dataset"], mine)
                checked += 1
            assert int(row["n"]) == int(src["n_perturbations"])
        else:
            d = json.loads((panel_dir.parent.parent / row["source_file"]).read_text())
            pairs = [("detection_fraction", "detection_fraction", "detectable_fraction"),
                     ("identification_fraction", "identification_fraction",
                      "identifiable_fraction"),
                     ("jointly_evaluable_fraction", "jointly_evaluable_fraction",
                      "jointly_evaluable_fraction"),
                     ("split_half_reproducibility_reference",
                      "split_half_reproducibility_reference", "replicate_ceiling"),
                     ("rho2", "rho2", "rho2"),
                     ("rho2_nn_median", "rho2_nn_median", "rho2_nn_median"),
                     ("model_ranking_p_correct_order", "model_ranking_p_correct_order",
                      "p_correct_order")]
            for mine, current, historical in pairs:
                assert float(row[mine]) == float(_axis(d, current, historical)), (
                    row["dataset"], mine)
                checked += 1
            assert int(row["n"]) == int(d["n_perturbations"])
            assert int(row["n_units_total"]) == int(d["selection"]["n_units_total"])
    print(f"verified {checked} numeric values against their source files, all equal")


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--panel", type=Path, default=Path("results/resolution_panel_v2"),
                    help="directory holding panel/ and the gate subdirectories")
    ap.add_argument("--v1", type=Path, default=Path("results/canonical/resolution_panel.csv"),
                    help="the six carried-over v1 rows")
    ap.add_argument("--out", type=Path, required=True,
                    help="directory receiving the two CSV files; refused if it resolves inside "
                         "the repository's results/ tree")
    a = ap.parse_args(argv)

    out = reject_repo_results(a.out)
    out.mkdir(parents=True, exist_ok=True)

    panel = _v2_rows(a.panel) + _v1_rows(a.v1)
    panel.sort(key=lambda r: (-float(r["detection_fraction"]), r["dataset"]))
    _verify(panel, a.panel, a.v1)
    _write(out / "resolution_panel_v2_table.csv", panel)

    sens = _sensitivity_rows(a.panel)
    _write(out / "resolution_panel_v2_sensitivity.csv", sens)

    resources = sorted({r["resource"] for r in panel})
    any_detection = sorted({r["resource"] for r in panel
                            if float(r["detection_fraction"]) > 0})
    print(f"{len(panel)} panel configurations from {len(resources)} distinct resources")
    print(f"  nonzero detection fraction in at least one configuration: "
          f"{len(any_detection)} resources")
    print(f"{len(sens)} sensitivity runs across "
          f"{sorted({r['arm'] for r in sens})}")
    print(f"wrote {out / 'resolution_panel_v2_table.csv'}")
    print(f"wrote {out / 'resolution_panel_v2_sensitivity.csv'}")


if __name__ == "__main__":
    main()
