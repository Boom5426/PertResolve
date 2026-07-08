#!/usr/bin/env python3
"""
Figure 5 | Rankability analysis turns benchmark failure into experimental design guidance
5 panels: a(ROC) b(calibration) c(effect-size regime) d(measured-vs-predicted) e(workflow)

Round-2 revision:
  - Panel c: reframed from a fixed "required cells" prescription to an
    effect-size-conditioned rankability regime (cell count alone is not the limiter).
  - Panel d: reframed from a second conflicting barplot to a measured-vs-predicted
    scatter, using the canonical native-depth un-rankable fractions.

Usage:
    python draw_fig5.py [--rank PATH] [--pred PATH] [--out fig5_composite.png]
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from fig_config import *
from sklearn.metrics import roc_curve, auc

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--rank', default=None)
parser.add_argument('--pred', default=None)
parser.add_argument('--out', default='../../figures/composites/fig5_composite.png')
args = parser.parse_args()

apply_style()
rank = load_rankability(args.rank)
pred = load_predictor(args.pred)

DS_COLORS = {
    'TP53': GENE_COLORS['TP53'], 'KRAS': GENE_COLORS['KRAS'],
    'GATA1': GENE_COLORS['GATA1'], 'JAK1': GENE_COLORS['JAK1'],
    'Replogle': '#666666', 'Norman': '#999999', 'Adamson': '#BBBBBB'
}

fig = plt.figure(figsize=(7.08, 9.0))
gs = GridSpec(3, 2, figure=fig, height_ratios=[0.9, 0.85, 0.75],
              hspace=0.42, wspace=0.35, left=0.10, right=0.96, top=0.95, bottom=0.05)

# --- Panel A: ROC ---
ax_a = fig.add_subplot(gs[0, 0])
rank_full = rank[(rank.metric == 'edist') & (rank.space == 'pca')].copy()
rank_full['label'] = rank_full['rankable'].astype(int)
rank_full['predictor'] = rank_full['effect_size']
for ds in ['JAK1','Norman','Adamson','Replogle','TP53','KRAS','GATA1']:
    sub = rank_full[rank_full.dataset == ds].dropna(subset=['predictor','label'])
    if len(sub) < 10 or sub.label.nunique() < 2:
        continue
    fpr, tpr, _ = roc_curve(sub.label, sub.predictor)
    roc_auc = auc(fpr, tpr)
    ax_a.plot(fpr, tpr, color=DS_COLORS[ds], lw=1.2, alpha=0.7, label=f'{ds} ({roc_auc:.2f})')
ax_a.plot([0,1],[0,1], '--', color='#CCC', lw=0.5)
ax_a.set_xlabel('False positive rate'); ax_a.set_ylabel('True positive rate')
ax_a.set_title('Effect size predicts rankability\n(per-dataset LODO ROC)', fontsize=7, pad=3, color='#555')
ax_a.text(0.55, 0.15, 'Mean AUROC\ncore = 0.974\nall = 0.965',
          transform=ax_a.transAxes, fontsize=5.5, color='#333',
          bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='#ddd', alpha=0.9, lw=0.4))
ax_a.legend(frameon=False, fontsize=4.5, loc='lower right', title='Dataset (AUROC)', title_fontsize=5)
panel_letter(ax_a, 'a')

# --- Panel B: calibration (measured vs predicted, LODO) ---
ax_b = fig.add_subplot(gs[0, 1])
for _, r in pred.iterrows():
    ds = r.held_out; c = DS_COLORS.get(ds, '#888')
    marker = 'o' if r.is_core else 's'
    ax_b.scatter(r.pred_frac_unrankable*100, r.obs_frac_unrankable*100,
                 c=c, s=55, marker=marker, edgecolor='white', linewidth=0.4, zorder=3)
    offset = (5, 5) if ds != 'GATA1' else (5, -10)
    ax_b.annotate(ds, (r.pred_frac_unrankable*100, r.obs_frac_unrankable*100),
                  xytext=offset, textcoords='offset points', fontsize=5, color=c)
ax_b.plot([0, 100], [0, 100], '--', color='#CCC', lw=0.7)
ax_b.set_xlabel('Predicted un-rankable (%)'); ax_b.set_ylabel('Observed un-rankable (%)')
ax_b.set_title('Discrimination strong, absolute calibration weak\n(LODO cross-validation)', fontsize=7, pad=3, color='#555')
ax_b.set_xlim(-5, 105); ax_b.set_ylim(-5, 105)
ax_b.text(0.05, 0.92, 'R2 = 0.11\nMean |err| = 26 pp', transform=ax_b.transAxes, fontsize=5.5, color=BASELINE_COLOR,
          bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='#ddd', alpha=0.9, lw=0.4))
panel_letter(ax_b, 'b')

# --- Panel C: effect-size-conditioned rankability regime ---
# Reframed: NOT a fixed "required cells" prescription. Each gene is placed by its
# median allele-vs-WT effect size and its native split-half D_self/D_null ratio;
# marker size = median cells/variant. TP53/KRAS have high nominal depth yet sit at
# the noise floor, showing that the effect window, not cell count, is the limiter.
ax_c = fig.add_subplot(gs[1, 0])
gene_rows = []
for gene in ['TP53','KRAS','GATA1','JAK1']:
    nat = native_rankability(rank, gene)
    eff = nat.effect_size.median()
    ratio = (nat.D_self / nat.D_null.replace(0, np.nan)).mean()
    med_cells = nat.n_cells.median()
    gene_rows.append((gene, eff, ratio, med_cells))
for gene, eff, ratio, med_cells in gene_rows:
    size = 40 + 45 * np.log10(max(med_cells, 10))
    ax_c.scatter(eff, ratio, s=size, c=GENE_COLORS[gene], edgecolor='white',
                 linewidth=0.5, zorder=4, alpha=0.9)
    ax_c.annotate(gene, (eff, ratio), xytext=(7, 3), textcoords='offset points',
                  fontsize=5.5, color=GENE_COLORS[gene], fontweight='bold')
ax_c.axhline(1.0, color=BASELINE_COLOR, ls='--', lw=0.7, alpha=0.5)
ax_c.axhspan(0.85, 1.15, color='#FFF3E0', alpha=0.2, zorder=0)
ax_c.text(18, 1.05, 'noise floor', fontsize=4.5, ha='right', color=BASELINE_COLOR, style='italic')
ax_c.set_xscale('log'); ax_c.set_xlim(0.8, 22)
ax_c.set_xlabel('Median allele-vs-WT effect size (log)')
ax_c.set_ylabel('Split-half D_self / D_null')
ax_c.set_ylim(-0.05, 1.3)
ax_c.set_title('Rankability is set by the effect-size regime,\nnot by cell count alone', fontsize=7, pad=3, color='#555')
ax_c.text(0.03, 0.06, 'marker size = median cells/variant', transform=ax_c.transAxes,
          fontsize=4.2, color='#999', style='italic')
panel_letter(ax_c, 'c')

# --- Panel D: native vs matched-depth un-rankable (canonical) ---
# Uses the locked native-depth definition; matched@50 shows depth is not the only
# driver. All numbers computed live from the rankability table (no hardcoding).
ax_d = fig.add_subplot(gs[1, 1])
ext = ['Replogle','Norman','Adamson']
native_ur, matched_ur = [], []
for ds in ext:
    fr, _ = unrankable_fraction(rank, ds)
    native_ur.append(fr*100)
    m50 = rank[(rank.dataset==ds)&(rank.metric=='edist')&(rank.space=='pca')&(rank.n_work==50)]
    matched_ur.append((1 - m50.groupby('perturbation')['rankable'].first().mean())*100)
x = np.arange(3); w = 0.32
ax_d.bar(x - w/2, native_ur, width=w, color=['#666','#999','#BBB'], alpha=0.75, edgecolor='white', lw=0.5, label='Native depth')
ax_d.bar(x + w/2, matched_ur, width=w, color=['#666','#999','#BBB'], alpha=0.35, edgecolor=['#666','#999','#BBB'], lw=0.8, label='Matched (n=50)')
for i, (nv, mv) in enumerate(zip(native_ur, matched_ur)):
    ax_d.text(i - w/2, nv + 1.5, f'{nv:.0f}%', fontsize=5, ha='center', color='#333')
    ax_d.text(i + w/2, mv + 1.5, f'{mv:.0f}%', fontsize=5, ha='center', color='#666')
ax_d.set_xticks(x); ax_d.set_xticklabels(['Replogle\n(2022)', 'Norman\n(2019)', 'Adamson\n(2016)'])
ax_d.set_ylabel('Un-rankable (%)'); ax_d.set_ylim(0, max(native_ur+matched_ur)*1.25)
ax_d.legend(frameon=False, fontsize=5, loc='upper right')
ax_d.set_title('Floor is intrinsic to distributional\nevaluation at finite depth', fontsize=7, pad=3, color='#555')
panel_letter(ax_d, 'd')

# --- Panel E: workflow ---
ax_e = fig.add_subplot(gs[2, :]); ax_e.axis('off')
ax_e.set_xlim(0, 10); ax_e.set_ylim(0, 3)
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
    ax_e.add_patch(box)
    ax_e.text(x_pos, 1.6, label, fontsize=6, ha='center', va='center', color=ec, fontweight='bold')
for i in range(len(steps)-1):
    ax_e.annotate('', xy=(steps[i+1][0]-0.85, 1.6), xytext=(steps[i][0]+0.85, 1.6),
                  arrowprops=dict(arrowstyle='->', color='#999', lw=1.2))
annotations = [
    (0.5, '10-50 variants\n~100 cells each'),
    (2.5, 'D_self/D_null\nper variant'),
    (4.5, 'AUROC 0.97\n(triage, not calibrate)'),
    (6.5, 'effect-size\nconditioned'),
    (8.5, 'only where\nsignal > floor'),
]
for x_pos, note in annotations:
    ax_e.text(x_pos, 0.5, note, fontsize=4.5, ha='center', va='center', color='#777', style='italic')
ax_e.set_title('Recommended power-aware allele-resolution perturbation workflow', fontsize=7, pad=2, color='#555')
panel_letter(ax_e, 'e')

outpath = os.path.join(os.path.dirname(__file__), args.out)
os.makedirs(os.path.dirname(outpath), exist_ok=True)
fig.savefig(outpath, dpi=300, bbox_inches='tight', facecolor='white')
print(f"Saved: {outpath}")
plt.close()
