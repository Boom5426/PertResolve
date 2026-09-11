"""Variant feature computation: θ₆ structural-biophysical vector."""
import numpy as np

# Kyte-Doolittle hydrophobicity scale
_KD = {
    'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
    'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
    'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
    'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2,
}

# Side-chain volume (Å³, Zamyatnin 1972)
_VOL = {
    'G': 60.1, 'A': 88.6, 'V': 140.0, 'L': 166.7, 'I': 166.7,
    'P': 112.7, 'F': 189.9, 'W': 227.8, 'M': 162.9, 'S': 89.0,
    'T': 116.1, 'C': 108.5, 'Y': 193.6, 'H': 153.2, 'D': 111.1,
    'E': 138.4, 'N': 114.1, 'Q': 143.8, 'K': 168.6, 'R': 173.4,
}

# Charge at pH 7
_CHARGE = {
    'R': 1, 'K': 1, 'D': -1, 'E': -1, 'H': 0.1,
}


def compute_theta(
    wt_aa: str,
    mut_aa: str,
    position: int,
    protein_length: int,
    domain_ranges: list = None,
    catalytic_residues: set = None,
    hotspot_positions: set = None,
) -> np.ndarray:
    """Compute the 6-dimensional θ feature vector for a single-residue substitution.

    Parameters
    ----------
    wt_aa, mut_aa : str
        Wild-type and mutant amino acid (single letter).
    position : int
        1-based residue position.
    protein_length : int
        Total protein length.
    domain_ranges : list of (start, end), optional
        Functional domain residue ranges. Residues in domain → fold_core = 1.0.
    catalytic_residues : set of int, optional
        Catalytic or functional-switch residue positions.
    hotspot_positions : set of int, optional
        Mutation hotspot positions.

    Returns
    -------
    theta : np.ndarray of shape (6,)
        [d_hydro, d_vol, d_charge, fold_core, cat_switch, is_hotspot]
    """
    d_hydro = _KD.get(mut_aa, 0) - _KD.get(wt_aa, 0)
    d_vol = _VOL.get(mut_aa, 0) - _VOL.get(wt_aa, 0)
    d_charge = _CHARGE.get(mut_aa, 0) - _CHARGE.get(wt_aa, 0)

    fold_core = 0.0
    if domain_ranges:
        for start, end in domain_ranges:
            if start <= position <= end:
                fold_core = 1.0
                break

    cat_switch = 1.0 if (catalytic_residues and position in catalytic_residues) else 0.0
    is_hotspot = 1.0 if (hotspot_positions and position in hotspot_positions) else 0.0

    return np.array([d_hydro, d_vol, d_charge, fold_core, cat_switch, is_hotspot])
