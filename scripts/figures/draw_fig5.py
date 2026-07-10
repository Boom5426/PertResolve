#!/usr/bin/env python3
"""
Figure 5 | Rankability analysis turns benchmark failure into experimental design guidance
4 panels: a(honest per-perturbation LODO AUROC) b(effect-size regime) c(native vs matched depth) d(workflow)

Round-3 revision (2026-07):
  - Panel a: replaced the earlier per-dataset in-sample ROC + hardcoded "AUROC 0.974/0.965"
    (a config-identity + pseudo-replication artifact) with the honest per-perturbation LODO
    predictor from scripts/figures/rankability_predictor.py: effect size as the sole feature,
    one row per perturbation, label-adjacent and configuration-identity features excluded.
    TP53 and KRAS are uniformly un-rankable at native depth and are shown as not evaluable.
  - Dropped the former Panel b (calibration of the artifact's predicted un-rankable fractions).
  - Panels relabeled: former c/d/e -> b/c/d.

Usage:
    python draw_fig5.py [--rank PATH] [--pred PATH] [--out fig5_composite.png]
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from fig_config import *

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--rank', default=None, help='rankability table CSV')
parser.add_argument('--pred', default=None, help='honest predictor results CSV (rankability_predictor_honest.csv)')
parser.add_argument('--out', default='../../figures/composites/fig5_composite.png')
args = parser.parse_args()

apply_style()
rank = load_rankability(args.rank)
HONEST = os.path.join(os.path.dirname(__file__), '..', '..', 'results', 'rankability_predictor_honest.csv')
hp = pd.read_csv(args.pred or HONEST)
hp = hp[hp.feature_set == 'effect_size'].copy()

DS_COLORS = {
    'TP53': GENE_COLORS['TP53'], 'KRAS': GENE_COLORS['KRAS'],
    'GATA1': GENE_COLORS['GATA1'], 'JAK1': GENE_COLORS['JAK1'],
    'Replogle': '#666666', 'Norman': '#999999', 'Adamson': '#BBBBBB'
}

fig = plt.figure(figsize=(7.08, 9.0))
gs = GridSpec(3, 2, figure=fig, height_ratios=[0.9, 0.85, 0.75],
              hspace=0.42, wspace=0.35, left=0.12, right=0.96, top=0.95, bottom=0.05)

# --- Panel A: honest per-perturbation LODO AUROC (effect size only) ---
ax_a = fig.add_subplot(gs[0, :])
order = ['JAK1', 'Replogle', 'GATA1', 'Adamson', 'Norman', 'KRAS', 'TP53']
hp = hp.set_index('held_out').loc[[d for d in order if d in hp['held_out'].values]].reset_index()
ys = np.arange(len(hp))[::-1]
for y, r in zip(ys, hp.itertuples()):
    c = DS_COLORS.get(r.held_out, '#888')
    if r.evaluable and not np.isnan(r.auroc):
        ax_a.plot([r.auroc_lo, r.auroc_hi], [y, y], color=c, lw=1.6, alpha=0.7, zorder=2)
        ax_a.scatter(r.auroc, y, s=42, color=c, edgecolor='white', linewidth=0.5, zorder=3)
        ax_a.text(min(r.auroc_hi + 0.008, 1.0), y, f'{r.auroc:.2f}', fontsize=5.5, va='center', color='#333')
    else:
        ax_a.text(0.515, y, 'not evaluable (0% rankable at native depth)',
                  fontsize=5, va='center', color='#999', style='italic')
ax_a.axvline(0.5, color=BASELINE_COLOR, ls='--', lw=0.7, alpha=0.6)
ax_a.set_yticks(ys); ax_a.set_yticklabels(hp['held_out'], fontsize=6)
ax_a.set_xlim(0.48, 1.03); ax_a.set_ylim(-0.6, len(hp) - 0.4)
ax_a.set_xlabel('LODO AUROC (per-perturbation, effect size only)')
ev = hp[hp['evaluable']]
ax_a.set_title('Per-perturbation effect size predicts rankability', fontsize=7, pad=3, color='#555')
ax_a.text(0.03, 0.44,
          f'mean AUROC = {ev["auroc"].mean():.2f} over {len(ev)} evaluable datasets\n'
          'TP53, KRAS not evaluable (uniformly un-rankable)',
          transform=ax_a.transAxes, fontsize=5, color='#333', va='center',
          bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='#ddd', alpha=0.9, lw=0.4))
panel_letter(ax_a, 'a')

# --- Panel B: effect-size-conditioned rankability regime ---
# Each gene is placed by its median allele-vs-WT effect size and its native split-half
# D_self/D_null ratio; marker size = median cells/variant. TP53/KRAS have high nominal
# depth yet sit at the noise floor, showing the effect window, not cell count, is the limiter.
ax_b = fig.add_subplot(gs[1, 0])
gene_rows = []
for gene in ['TP53', 'KRAS', 'GATA1', 'JAK1']:
    nat = native_rankability(rank, gene)
    eff = nat.effect_size.median()
    ratio = (nat.D_self / nat.D_null.replace(0, np.nan)).mean()
    med_cells = nat.n_cells.median()
    gene_rows.append((gene, eff, ratio, med_cells))
for gene, eff, ratio, med_cells in gene_rows:
    size = 40 + 45 * np.log10(max(med_cells, 10))
    ax_b.scatter(eff, ratio, s=size, c=GENE_COLORS[gene], edgecolor='white',
                 linewidth=0.5, zorder=4, alpha=0.9)
    ax_b.annotate(gene, (eff, ratio), xytext=(7, 3), textcoords='offset points',
                  fontsize=5.5, color=GENE_COLORS[gene], fontweight='bold')
ax_b.axhline(1.0, color=BASELINE_COLOR, ls='--', lw=0.7, alpha=0.5)
ax_b.axhspan(0.85, 1.15, color='#FFF3E0', alpha=0.2, zorder=0)
ax_b.text(18, 1.05, 'noise floor', fontsize=4.5, ha='right', color=BASELINE_COLOR, style='italic')
ax_b.set_xscale('log'); ax_b.set_xlim(0.8, 22)
ax_b.set_xlabel('Median allele-vs-WT effect size (log)')
ax_b.set_ylabel('Split-half D_self / D_null')
ax_b.set_ylim(-0.05, 1.3)
ax_b.set_title('Rankability is set by the effect-size regime,\nnot by cell count alone', fontsize=7, pad=3, color='#555')
ax_b.text(0.03, 0.06, 'marker size = median cells/variant', transform=ax_b.transAxes,
          fontsize=4.2, color='#999', style='italic')
panel_letter(ax_b, 'b')

# --- Panel C: native vs matched-depth un-rankable (canonical) ---
# Uses the locked native-depth definition; matched@50 shows depth is not the only
# driver. All numbers computed live from the rankability table (no hardcoding).
ax_c = fig.add_subplot(gs[1, 1])
ext = ['Replogle', 'Norman', 'Adamson']
native_ur, matched_ur = [], []
for ds in ext:
    fr, _ = unrankable_fraction(rank, ds)
    native_ur.append(fr * 100)
    m50 = rank[(rank.dataset == ds) & (rank.metric == 'edist') & (rank.space == 'pca') & (rank.n_work == 50)]
    matched_ur.append((1 - m50.groupby('perturbation')['rankable'].first().mean()) * 100)
x = np.arange(3); w = 0.32
ax_c.bar(x - w/2, native_ur, width=w, color=['#666', '#999', '#BBB'], alpha=0.75, edgecolor='white', lw=0.5, label='Native depth')
ax_c.bar(x + w/2, matched_ur, width=w, color=['#666', '#999', '#BBB'], alpha=0.35, edgecolor=['#666', '#999', '#BBB'], lw=0.8, label='Matched (n=50)')
for i, (nv, mv) in enumerate(zip(native_ur, matched_ur)):
    ax_c.text(i - w/2, nv + 1.5, f'{nv:.0f}%', fontsize=5, ha='center', color='#333')
    ax_c.text(i + w/2, mv + 1.5, f'{mv:.0f}%', fontsize=5, ha='center', color='#666')
ax_c.set_xticks(x); ax_c.set_xticklabels(['Replogle\n(2022)', 'Norman\n(2019)', 'Adamson\n(2016)'])
ax_c.set_ylabel('Un-rankable (%)'); ax_c.set_ylim(0, max(native_ur + matched_ur) * 1.25)
ax_c.legend(frameon=False, fontsize=5, loc='upper right')
ax_c.set_title('Floor is intrinsic to distributional\nevaluation at finite depth', fontsize=7, pad=3, color='#555')
panel_letter(ax_c, 'c')

# --- Panel D: workflow ---
ax_d = fig.add_subplot(gs[2, :]); ax_d.axis('off')
ax_d.set_xlim(0, 10); ax_d.set_ylim(0, 3)
steps = [
    (0.5, 'Pilot variant\nscreen', '#E8F5E9', '#55966B'),
    (2.5, 'Estimate\nsplit-half floor', '#E3F2FD', '#5185C0'),
    (4.5, 'Compute\nrankability', '#FFF3E0', '#E99D4E'),
    (6.5, 'Set cells/variant\ntarget', '#F3E5F5', '#8281B9'),
    (8.5, 'Evaluate\npredictions', '#FFEBEE', '#C96144'),
]
for x_pos, label, bg, ec in steps:
    box = FancyBboxPatch((x_pos - 0.8, 1.0), 1.6, 1.2,
                          boxstyle="round,pad=0.15", facecolor=bg, edgecolor=ec, linewidth=1.0)
    ax_d.add_patch(box)
    ax_d.text(x_pos, 1.6, label, fontsize=6, ha='center', va='center', color=ec, fontweight='bold')
for i in range(len(steps) - 1):
    ax_d.annotate('', xy=(steps[i+1][0]-0.85, 1.6), xytext=(steps[i][0]+0.85, 1.6),
                  arrowprops=dict(arrowstyle='->', color='#999', lw=1.2))
annotations = [
    (0.5, '10-50 variants\n~100 cells each'),
    (2.5, 'D_self/D_null\nper variant'),
    (4.5, 'AUROC ~0.96\n(triage where evaluable)'),
    (6.5, 'effect-size\nconditioned'),
    (8.5, 'only where\nsignal > floor'),
]
for x_pos, note in annotations:
    ax_d.text(x_pos, 0.5, note, fontsize=4.5, ha='center', va='center', color='#777', style='italic')
ax_d.set_title('Recommended power-aware allele-resolution perturbation workflow', fontsize=7, pad=2, color='#555')
panel_letter(ax_d, 'd')

outpath = os.path.join(os.path.dirname(__file__), args.out)
os.makedirs(os.path.dirname(outpath), exist_ok=True)
fig.savefig(outpath, dpi=300, bbox_inches='tight', facecolor='white')
print(f"Saved: {outpath}")
plt.close()
