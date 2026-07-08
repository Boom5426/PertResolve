#!/usr/bin/env python3
"""
Figure 4 | Direction-ranking dissociation is robust across splits, metrics, and features
6 panels: a(PDS heatmap) b(P-delta heatmap) c(PDS distance robustness) d(gene×split) e(feature bars) f(compatibility table)

Usage:
    python draw_fig4.py [--data PATH] [--out fig4_composite.png]
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from fig_config import *

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--data', default=None)
parser.add_argument('--out', default='../../figures/composites/fig4_composite.png')
args = parser.parse_args()

apply_style()
df = load_v4(args.data)

metrics = ['PDS_cos','PDS_L1','PDS_L2','pearson_delta','delta_cosine',
           'DE_overlap','DE_LFC_spearman','direction_agreement','MAE','pearson_delta_top20']
splits_order = ['split1','split2','split3','split5','split6']

per_gene = df.groupby(['method','split','gene'])[metrics].mean().reset_index()
per_split = per_gene.groupby(['method','split'])[metrics].mean().reset_index()
overall_pds = per_split.groupby('method')['PDS_cos'].mean().sort_values(ascending=False)
top12 = overall_pds.index[:12].tolist()

pivot_pds = per_split.pivot(index='method', columns='split', values='PDS_cos')
pivot_pd = per_split.pivot(index='method', columns='split', values='pearson_delta')

fig = plt.figure(figsize=(7.08, 9.0))
gs = GridSpec(3, 2, figure=fig, height_ratios=[1.05, 0.85, 0.85],
              hspace=0.40, wspace=0.35, left=0.12, right=0.96, top=0.95, bottom=0.05)

# --- Panel A: PDS heatmap ---
ax_a = fig.add_subplot(gs[0, 0])
pds_mat = pivot_pds.loc[top12, [s for s in splits_order if s in pivot_pds.columns]].fillna(np.nan)
im_a = ax_a.imshow(pds_mat.values, cmap='RdYlGn', vmin=0.28, vmax=0.55, aspect='auto')
ax_a.set_yticks(range(len(top12)))
ax_a.set_yticklabels([clean_name(m) for m in top12], fontsize=5)
ax_a.set_xticks(range(pds_mat.shape[1]))
ax_a.set_xticklabels([SPLIT_LABELS.get(s, s) for s in pds_mat.columns], fontsize=5.5, rotation=30, ha='right')
for i in range(pds_mat.shape[0]):
    for j in range(pds_mat.shape[1]):
        v = pds_mat.iloc[i, j]
        if not np.isnan(v):
            ax_a.text(j, i, f'{v:.2f}', fontsize=4, ha='center', va='center',
                      color='white' if v < 0.35 else '#333')
ax_a.set_title('PDS across splits: near or below chance', fontsize=7, pad=3, color='#555')
cb_a = fig.colorbar(im_a, ax=ax_a, shrink=0.65, pad=0.02)
cb_a.ax.tick_params(labelsize=4.5); cb_a.set_label('PDS', fontsize=5.5)
panel_letter(ax_a, 'a')

# --- Panel B: P-delta heatmap ---
ax_b = fig.add_subplot(gs[0, 1])
pd_mat = pivot_pd.loc[top12, [s for s in splits_order if s in pivot_pd.columns]].fillna(np.nan)
im_b = ax_b.imshow(pd_mat.values, cmap='Blues', vmin=0.0, vmax=0.80, aspect='auto')
ax_b.set_yticks(range(len(top12)))
ax_b.set_yticklabels([clean_name(m) for m in top12], fontsize=5)
ax_b.set_xticks(range(pd_mat.shape[1]))
ax_b.set_xticklabels([SPLIT_LABELS.get(s, s) for s in pd_mat.columns], fontsize=5.5, rotation=30, ha='right')
for i in range(pd_mat.shape[0]):
    for j in range(pd_mat.shape[1]):
        v = pd_mat.iloc[i, j]
        if not np.isnan(v):
            ax_b.text(j, i, f'{v:.2f}', fontsize=4, ha='center', va='center',
                      color='white' if v > 0.55 else '#333')
ax_b.set_title('Pearson delta: consistently positive', fontsize=7, pad=3, color='#555')
cb_b = fig.colorbar(im_b, ax=ax_b, shrink=0.65, pad=0.02)
cb_b.ax.tick_params(labelsize=4.5); cb_b.set_label('Pearson delta', fontsize=5.5)
panel_letter(ax_b, 'b')

# --- Panel C: PDS distance robustness ---
ax_c = fig.add_subplot(gs[1, 0])
mgm_nonwt = per_gene.groupby('method')[['PDS_cos','PDS_L1','PDS_L2']].mean()
mgm_nonwt = mgm_nonwt.drop(['WT-null','Gene-mean'], errors='ignore')
pds_names = ['PDS\n(cosine)', 'PDS\n(L1)', 'PDS\n(L2)']
pds_cols = ['PDS_cos', 'PDS_L1', 'PDS_L2']
for i, (col, label) in enumerate(zip(pds_cols, pds_names)):
    vals = mgm_nonwt[col].values
    np.random.seed(42 + i)
    jx = np.random.uniform(-0.12, 0.12, len(vals))
    ax_c.scatter(i + jx, vals, c='#5185C0', s=14, alpha=0.5, edgecolor='white', linewidth=0.2, zorder=3)
    med = np.median(vals)
    ax_c.plot([i-0.18, i+0.18], [med, med], '-', color='#333', lw=1.5, zorder=5)
    ax_c.text(i+0.22, med, f'{med:.3f}', fontsize=5, va='center', color='#333')
ax_c.axhline(0.50, color=BASELINE_COLOR, ls='--', lw=0.7, alpha=0.5)
ax_c.set_xticks(range(3)); ax_c.set_xticklabels(pds_names, fontsize=6)
ax_c.set_ylabel('PDS'); ax_c.set_ylim(0.42, 0.52)
ax_c.set_title('Ranking failure is invariant to\ndistance function', fontsize=7, pad=3, color='#555')
panel_letter(ax_c, 'c')

# --- Panel D: gene x split scatter ---
ax_d = fig.add_subplot(gs[1, 1])
gene_split = per_gene.groupby(['gene','split'])[['PDS_cos','pearson_delta']].mean().reset_index()
for gene in ['TP53','KRAS','GATA1','JAK1']:
    g = gene_split[gene_split.gene == gene]
    ax_d.scatter(g.PDS_cos, g.pearson_delta, c=GENE_COLORS[gene], s=45, alpha=0.75,
                 edgecolor='white', linewidth=0.4, zorder=3, label=gene)
    g_s = g.sort_values('split')
    ax_d.plot(g_s.PDS_cos, g_s.pearson_delta, '-', color=GENE_COLORS[gene], alpha=0.2, lw=0.7)
    for _, r in g.iterrows():
        sl = SPLIT_LABELS.get(r.split, r.split)[:3]
        ax_d.annotate(sl, (r.PDS_cos, r.pearson_delta), fontsize=3.5, color='#666',
                       xytext=(3, 3), textcoords='offset points')
ax_d.axvline(0.50, color=BASELINE_COLOR, ls='--', lw=0.5, alpha=0.3)
ax_d.set_xlabel('PDS (cosine)'); ax_d.set_ylabel('Pearson delta')
ax_d.set_xlim(0.15, 0.60); ax_d.set_ylim(-0.05, 0.80)
ax_d.set_title('Dissociation holds across\nall genes and splits', fontsize=7, pad=3, color='#555')
ax_d.legend(frameon=False, fontsize=5.5, loc='lower right')
panel_letter(ax_d, 'd')

# --- Panel E: feature type bars ---
ax_e = fig.add_subplot(gs[2, 0])
per_split_ft = per_split.copy()
per_split_ft['feat'] = per_split_ft.method.apply(feat_type)
feat_overall = per_split_ft[per_split_ft.feat != 'Baseline'].groupby('feat')[['PDS_cos','pearson_delta']].mean()
feat_order = ['θ', 'ESM', 'ESM+θ']
feat_colors_map = {'θ': '#E99D4E', 'ESM': '#5185C0', 'ESM+θ': '#8281B9'}
x = np.arange(len(feat_order)); w = 0.30
ax_e.bar(x - w/2, [feat_overall.loc[f, 'PDS_cos'] for f in feat_order],
         width=w, color=[feat_colors_map[f] for f in feat_order], alpha=0.45, edgecolor='white', lw=0.5)
ax_e.bar(x + w/2, [feat_overall.loc[f, 'pearson_delta'] for f in feat_order],
         width=w, color=[feat_colors_map[f] for f in feat_order], alpha=0.85, edgecolor='white', lw=0.5)
for i, f in enumerate(feat_order):
    pv, dv = feat_overall.loc[f, 'PDS_cos'], feat_overall.loc[f, 'pearson_delta']
    ax_e.text(i - w/2, pv + 0.01, f'{pv:.2f}', fontsize=4.5, ha='center', color='#888')
    ax_e.text(i + w/2, dv + 0.01, f'{dv:.2f}', fontsize=4.5, ha='center', color='#333')
ax_e.axhline(0.50, color=BASELINE_COLOR, ls=':', lw=0.5, alpha=0.3)
ax_e.set_xticks(x); ax_e.set_xticklabels(feat_order, fontsize=7, fontweight='bold')
for tick, f in zip(ax_e.xaxis.get_ticklabels(), feat_order): tick.set_color(feat_colors_map[f])
ax_e.set_ylabel('Score'); ax_e.set_ylim(0, 0.78)
ax_e.set_title('Richer features do not close\nthe ranking gap', fontsize=7, pad=3, color='#555')
leg_e = [Patch(fc='#999', alpha=0.45, label='PDS'), Patch(fc='#999', alpha=0.85, label='Pearson delta')]
ax_e.legend(handles=leg_e, frameon=False, fontsize=5, loc='upper left')
panel_letter(ax_e, 'e')

# --- Panel F: method compatibility table ---
ax_f = fig.add_subplot(gs[2, 1]); ax_f.axis('off')
ext_methods = [
    ('CPA', 'Dose scalar input', 'NaN on θ vector'),
    ('CellFlow', 'No unseen-variant API', 'Train-set only'),
    ('STATE', 'Per-dataset retrain', 'GPU OOM (GATA1)'),
    ('GEARS', 'Gene-keyed graph', 'Allele-blind'),
    ('scDFM', 'No public checkpoint', 'Cannot reproduce'),
    ('Biolord', 'Ordered attributes', 'PDS ≈ chance'),
]
ax_f.set_title('Existing methods assume gene-level\nperturbation identity', fontsize=7, pad=3, color='#555')
col_headers = ['Method', 'Failure mode', 'Outcome']
col_widths = [0.18, 0.42, 0.35]
row_h = 0.11; start_y = 0.82
for j, (h, cw) in enumerate(zip(col_headers, col_widths)):
    x_pos = sum(col_widths[:j]) + cw/2
    ax_f.text(x_pos, start_y + row_h, h, fontsize=5.5, ha='center', va='center',
              fontweight='bold', color='#333', transform=ax_f.transAxes)
for i, (name, mode, outcome) in enumerate(ext_methods):
    y = start_y - i * row_h
    bg = '#F5F5F5' if i % 2 == 0 else 'white'
    rect = plt.Rectangle((0, y - row_h/2), 1, row_h, transform=ax_f.transAxes,
                           facecolor=bg, edgecolor='none', zorder=0)
    ax_f.add_patch(rect)
    ax_f.text(col_widths[0]/2, y, name, fontsize=5.5, ha='center', va='center',
              fontweight='bold', color='#555', transform=ax_f.transAxes)
    ax_f.text(col_widths[0] + col_widths[1]/2, y, mode, fontsize=5, ha='center', va='center',
              color='#666', transform=ax_f.transAxes)
    ax_f.text(col_widths[0] + col_widths[1] + col_widths[2]/2, y, outcome, fontsize=5,
              ha='center', va='center', color=BASELINE_COLOR, style='italic', transform=ax_f.transAxes)
panel_letter(ax_f, 'f')

outpath = os.path.join(os.path.dirname(__file__), args.out)
os.makedirs(os.path.dirname(outpath), exist_ok=True)
fig.savefig(outpath, dpi=300, bbox_inches='tight', facecolor='white')
print(f"Saved: {outpath}")
plt.close()
