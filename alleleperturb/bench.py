"""AllelePerturb benchmark loader and split manager."""

import os
import numpy as np
import pandas as pd
from pathlib import Path

_BENCH_FILE = "allele_perturb_bench.csv"
_SPLIT_NAMES = {
    "split1": "Random",
    "split2": "OOD-Position",
    "split3": "OOD-Mechanism",
    "split4": "Cross-Gene",
    "split5": "Low-N",
    "split6": "PerturbNet-Compat",
}
_THETA_COLS = ["d_hydro", "d_vol", "d_charge", "fold_core", "cat_switch", "is_hotspot"]
_GENE_ORDER = ["TP53", "KRAS", "GATA1", "JAK1"]


class AllelePerturb:
    """AllelePerturb benchmark interface.
    
    Parameters
    ----------
    data_dir : str or Path, optional
        Directory containing ``allele_perturb_bench.csv`` and gene-level
        expression arrays. Defaults to ``<package>/data/``.
    
    Examples
    --------
    >>> bench = AllelePerturb.load()
    >>> train_vars, test_vars = bench.split("split1", gene="TP53")
    >>> theta = bench.get_theta("TP53")
    """

    def __init__(self, df: pd.DataFrame, data_dir: str = None):
        self._df = df
        self._data_dir = Path(data_dir) if data_dir else None

    @classmethod
    def load(cls, data_dir: str = None) -> "AllelePerturb":
        """Load the benchmark from disk.
        
        Parameters
        ----------
        data_dir : str, optional
            Path to the data directory. If None, looks for
            ``data/allele_perturb_bench.csv`` relative to the package root.
        """
        if data_dir is None:
            pkg_dir = Path(__file__).parent.parent / "data"
            bench_path = pkg_dir / _BENCH_FILE
        else:
            bench_path = Path(data_dir) / _BENCH_FILE

        if not bench_path.exists():
            raise FileNotFoundError(
                f"Benchmark file not found at {bench_path}. "
                f"Download it or specify data_dir."
            )
        df = pd.read_csv(bench_path)
        return cls(df, data_dir=str(bench_path.parent))

    def __repr__(self) -> str:
        genes = self._df["gene"].value_counts()
        total = len(self._df)
        cells = int(self._df["n_cells"].sum())
        return (
            f"AllelePerturb: {total} variants, {len(genes)} genes, "
            f"{cells:,} cells\n"
            f"  Genes: {', '.join(f'{g}({n})' for g, n in genes.items())}\n"
            f"  Splits: {', '.join(_SPLIT_NAMES.values())}"
        )

    @property
    def genes(self) -> list:
        return [g for g in _GENE_ORDER if g in self._df["gene"].values]

    @property
    def variants(self) -> pd.DataFrame:
        return self._df.copy()

    def split(self, split_name: str, gene: str = None):
        """Return (train_variants, test_variants) for a given split.
        
        Parameters
        ----------
        split_name : str
            One of "split1" through "split6".
        gene : str, optional
            Filter to a specific gene.
        
        Returns
        -------
        train : list of str
            Training variant names.
        test : list of str
            Held-out test variant names.
        """
        col = f"{split_name}_role"
        if col not in self._df.columns:
            raise ValueError(f"Unknown split: {split_name}. Use split1–split6.")
        df = self._df if gene is None else self._df[self._df["gene"] == gene]
        train = df[df[col] == "train"]["variant"].tolist()
        test = df[df[col] == "test"]["variant"].tolist()
        return train, test

    def get_theta(self, gene: str = None) -> dict:
        """Return θ₆ feature vectors as {variant_name: np.array(6)}.
        
        Parameters
        ----------
        gene : str, optional
            Filter to one gene.
        """
        df = self._df if gene is None else self._df[self._df["gene"] == gene]
        theta = {}
        for _, row in df.iterrows():
            theta[row["variant"]] = np.array([row[c] for c in _THETA_COLS])
        return theta

    def gene_data(self, gene: str) -> dict:
        """Return per-gene metadata.
        
        Returns dict with keys: n_variants, n_cells, protein, cell_type, technology.
        """
        sub = self._df[self._df["gene"] == gene]
        if len(sub) == 0:
            raise ValueError(f"Gene {gene} not in benchmark.")
        row = sub.iloc[0]
        return {
            "n_variants": len(sub),
            "n_cells": int(sub["n_cells"].sum()),
            "protein": row.get("protein", ""),
            "cell_type": row.get("cell_type", ""),
            "dataset": row.get("dataset", ""),
        }

    def split_summary(self) -> pd.DataFrame:
        """Return a DataFrame of test-variant counts per (split, gene)."""
        rows = []
        for sp, name in _SPLIT_NAMES.items():
            col = f"{sp}_role"
            if col not in self._df.columns:
                continue
            for gene in self.genes:
                n = (
                    (self._df["gene"] == gene) & (self._df[col] == "test")
                ).sum()
                rows.append({"split": name, "gene": gene, "n_test": int(n)})
        return pd.DataFrame(rows).pivot(index="split", columns="gene", values="n_test")
