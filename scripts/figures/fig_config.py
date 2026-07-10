"""
Shared config for AllelePerturb figure scripts.
Run from: AllelePerturb/scripts/figures/
Data from: AllelePerturb/results/ (pre-computed, included in repo)
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, FancyBboxPatch
from matplotlib.gridspec import GridSpec

# ============================================================
# Style
# ============================================================
def apply_style():
    """Publication-grade NM style."""
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'font.size': 7,
        'axes.titlesize': 7,
        'axes.labelsize': 7,
        'xtick.labelsize': 6,
        'ytick.labelsize': 6,
        'legend.fontsize': 5.5,
        'figure.dpi': 300,
        'savefig.dpi': 300,
        'savefig.bbox': 'tight',
        'savefig.facecolor': 'white',
        'axes.spines.top': False,
        'axes.spines.right': False,
        'axes.linewidth': 0.6,
        'xtick.major.width': 0.5,
        'ytick.major.width': 0.5,
    })

def panel_letter(ax, letter, x=-0.12, y=1.08):
    """Add bold panel letter."""
    ax.text(x, y, letter, transform=ax.transAxes,
            fontsize=10, fontweight='bold', va='top', ha='left')

# ============================================================
# Colors
# ============================================================
GENE_COLORS = {
    'TP53': '#5185C0',
    'KRAS': '#E99D4E',
    'GATA1': '#8281B9',
    'JAK1': '#55966B',
}
BASELINE_COLOR = '#C96144'

FEAT_COLORS = {
    'θ': '#E99D4E',
    'ESM': '#5185C0',
    'ESM+θ': '#8281B9',
    'Baseline': '#C96144',
}

HEAD_MARKERS = {
    'Linear': 'o',
    'Tree': 's',
    'KNN': 'D',
    'MLP': '^',
    'Baseline': '*',
}

SPLIT_LABELS = {
    'split1': 'Random',
    'split2': 'OOD-Pos',
    'split3': 'OOD-Mech',
    'split5': 'Low-N',
    'split6': 'Compat.',
}

# ============================================================
# Helpers
# ============================================================
def clean_name(m):
    return (m.replace('-esm+theta', ' (ESM+θ)').replace('-esm', ' (ESM)')
             .replace('-theta', ' (θ)').replace('WT-null', 'WT null')
             .replace('Gene-mean', 'Gene mean').replace('GBoost', 'GBT'))

def feat_type(m):
    if m in ['WT-null', 'Gene-mean']: return 'Baseline'
    if 'esm+theta' in m.lower(): return 'ESM+θ'
    if 'esm' in m.lower(): return 'ESM'
    return 'θ'

def head_type(m):
    if m in ['WT-null', 'Gene-mean']: return 'Baseline'
    low = m.split('-')[0].lower()
    if low in ['ridge', 'lasso']: return 'Linear'
    if low in ['rf', 'gboost']: return 'Tree'
    if low == 'knn': return 'KNN'
    if low == 'mlp': return 'MLP'
    return 'Linear'

# ============================================================
# Data loading
# ============================================================
# Adjust these paths to your local setup
DATA_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
RESULTS_V4 = None  # will be set by load_v4()
RANKABILITY = None

def load_v4(path=None):
    """Load the 10-metric v4 results CSV."""
    global RESULTS_V4
    candidates = [
        path,
        f'{DATA_DIR}/results/results_v4_10metrics.csv',
    ]
    for p in candidates:
        if p and os.path.isfile(p):
            RESULTS_V4 = pd.read_csv(p)
            print(f"Loaded v4: {RESULTS_V4.shape[0]} rows from {p}")
            return RESULTS_V4
    raise FileNotFoundError("results_v4_10metrics.csv not found")

def load_rankability(path=None):
    """Load the split-half rankability table."""
    global RANKABILITY
    candidates = [
        path,
        f'{DATA_DIR}/results/second_probe_rankability_table.csv',
    ]
    for p in candidates:
        if p and os.path.isfile(p):
            RANKABILITY = pd.read_csv(p)
            print(f"Loaded rankability: {RANKABILITY.shape[0]} rows from {p}")
            return RANKABILITY
    raise FileNotFoundError("rankability table not found")

def load_power(path=None):
    """Load power curve CSV."""
    candidates = [
        path,
        f'{DATA_DIR}/results/split_half_power_curve.csv',
    ]
    for p in candidates:
        if p and os.path.isfile(p):
            return pd.read_csv(p)
    raise FileNotFoundError("power curve CSV not found")

def load_predictor(path=None):
    """DEPRECATED. Loads the superseded row-level LODO predictor output
    (results/deprecated/second_probe_predictor_results.csv), whose headline AUROC 0.974
    was shown to be a config-identity + pseudo-replication artifact. New work should use
    the honest per-perturbation predictor in scripts/figures/rankability_predictor.py
    (results/rankability_predictor_honest.csv). Kept only so historical panels still load."""
    candidates = [
        path,
        f'{DATA_DIR}/results/deprecated/second_probe_predictor_results.csv',
    ]
    for p in candidates:
        if p and os.path.isfile(p):
            return pd.read_csv(p)
    raise FileNotFoundError("predictor results not found")

# ============================================================
# Canonical rankability aggregation
# ============================================================
# Locked definition (Round 2 revision): native max-depth (deepest available
# split-half per perturbation) + all perturbations + direct split-half energy
# distance, PCA-50 space. This is the SINGLE canonical definition used in the
# main text (Fig 3d, Fig 5d). Matched-depth and significance-filtered variants
# go to Extended Data only.

def native_rankability(rank_df, dataset, metric='edist', space='pca'):
    """Return the native-max-depth per-perturbation rankability rows for a dataset.
    Native = each perturbation evaluated at its DEEPEST available split-half bin
    (max n_work), NOT the first/shallowest bin."""
    sub = rank_df[(rank_df.dataset == dataset) &
                  (rank_df.metric == metric) &
                  (rank_df.space == space)]
    if len(sub) == 0:
        return sub
    return sub.loc[sub.groupby('perturbation')['n_work'].idxmax()]

def unrankable_fraction(rank_df, dataset, metric='edist', space='pca'):
    """Canonical native-depth un-rankable fraction (0-1) for a dataset."""
    native = native_rankability(rank_df, dataset, metric, space)
    if len(native) == 0:
        return np.nan, 0
    return 1.0 - native['rankable'].mean(), len(native)

def dself_dnull_native(rank_df, dataset, metric='edist', space='pca'):
    """Canonical mean per-perturbation D_self/D_null at native depth."""
    native = native_rankability(rank_df, dataset, metric, space)
    if len(native) == 0:
        return np.nan
    return (native['D_self'] / native['D_null'].replace(0, np.nan)).mean()
