#!/usr/bin/env python3
"""
Figure 2 | Direction recovery succeeds but allele-specific ranking fails
5 panels: a(paired dots) b(dissociation scatter) c(DE fidelity) d(gene multiples) e(per-variant strip)

Usage:
    python draw_fig2.py [--data PATH_TO_V4_CSV] [--out fig2_composite.png]
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from fig_config import *

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--data', default=None)
parser.add_argument('--out', default='../../figures/composites/fig2_composite.png')
args = parser.parse_args()

apply_style()
df = load_v4(args.data)

# ============================================================
# Precompute leaderboard
# ============================================================
metrics = ['PDS_cos','PDS_L1','PDS_L2','pearson_delta','delta_cosine',
           'DE_overlap','DE_LFC_spearman','direction_agreement','MAE','pearson_delta_top20']

per_gene = df.groupby(['method','split','gene'])[metrics].mean().reset_index()
mgm = per_gene.groupby(['method','split'])[metrics].mean().reset_index()
overall = mgm.groupby('method')[metrics].mean()

# Sort methods by Pearson delta descending
sorted_methods = overall.sort_values('pearson_delta', ascending=False).index.tolist()

# Per-gene overall
per_gene_overall = per_gene.groupby(['method','gene'])[metrics].mean().reset_index()

# Bootstrap 95% CIs over per-variant PDS / Pearson delta (resampling unit = variant row)
_rng = np.random.default_rng(0)
def _boot_ci(vals, n_boot=2000):
    vals = np.asarray(vals); vals = vals[~np.isnan(vals)]
    if len(vals) < 3: return np.nan, np.nan, np.nan
    b = [np.mean(_rng.choice(vals, len(vals), replace=True)) for _ in range(n_boot)]
    return np.mean(vals), np.percentile(b, 2.5), np.percentile(b, 97.5)
pds_ci, pd_ci = {}, {}
for m in overall.index:
    pds_ci[m] = _boot_ci(df[df.method==m]['PDS_cos'].values)
    pd_ci[m]  = _boot_ci(df[df.method==m]['pearson_delta'].values)

# ============================================================
# Figure
# ============================================================
fig = plt.figure(figsize=(7.08, 9.5))
gs = GridSpec(3, 2, figure=fig,
              height_ratios=[1.05, 0.80, 0.45],
              width_ratios=[1.05, 0.95],
              hspace=0.35, wspace=0.32,
              left=0.10, right=0.96, top=0.97, bottom=0.05)

# ===================== PANEL A =====================
gs_a = gs[0, 0].subgridspec(1, 2, wspace=0.06, width_ratios=[1, 1.15])
ax_l = fig.add_subplot(gs_a[0])
ax_r = fig.add_subplot(gs_a[1], sharey=ax_l)

y = np.arange(len(sorted_methods))
for i, method in enumerate(sorted_methods):
    ft = feat_type(method)
    color = FEAT_COLORS[ft]
    pds = overall.loc[method, 'PDS_cos']
    pdv = overall.loc[method, 'pearson_delta']
    ms = 5.5 if ft != 'Baseline' else 7
    marker = 'o' if ft != 'Baseline' else '*'
    # 95% bootstrap CI whiskers
    if method in pds_ci and not np.isnan(pds_ci[method][1]):
        _, plo, phi = pds_ci[method]
        ax_l.plot([plo, phi], [i, i], '-', color=color, lw=0.7, alpha=0.5, zorder=3)
    if method in pd_ci and not np.isnan(pd_ci[method][1]):
        _, dlo, dhi = pd_ci[method]
        ax_r.plot([dlo, dhi], [i, i], '-', color=color, lw=0.7, alpha=0.5, zorder=3)
    ax_l.plot(pds, i, marker=marker, color=color, markersize=ms,
              markeredgecolor='white', markeredgewidth=0.4, zorder=4, linestyle='none')
    ax_r.plot(pdv, i, marker=marker, color=color, markersize=ms,
              markeredgecolor='white', markeredgewidth=0.4, zorder=4, linestyle='none')

ax_l.axvspan(0.48, 0.52, color='#F0F0F0', alpha=0.5, zorder=0, lw=0)
ax_l.axvline(0.50, color=BASELINE_COLOR, ls='-', lw=0.6, alpha=0.4, zorder=1)
ax_l.set_xlabel('PDS (cosine)')
ax_l.set_xlim(0.40, 0.54)
ax_l.set_title('No method exceeds chance\n(95% CI whiskers)', fontsize=6.5, pad=3, color='#555')
ax_l.set_yticks(y)
ax_l.set_yticklabels([clean_name(m) for m in sorted_methods], fontsize=5)
ax_l.set_ylim(-0.7, len(sorted_methods) - 0.3)

ax_r.axvline(0.0, color='#BBB', ls='-', lw=0.5, alpha=0.3, zorder=1)
ax_r.set_xlabel('Pearson delta')
ax_r.set_xlim(-0.03, 0.73)
ax_r.set_title('Direction recovery is\nconsistently positive', fontsize=6.5, pad=3, color='#555')
plt.setp(ax_r.get_yticklabels(), visible=False)

leg_a = [
    Line2D([0],[0], marker='o', color='w', markerfacecolor=FEAT_COLORS['θ'],
           markersize=4.5, markeredgecolor='white', markeredgewidth=0.3, label='θ only'),
    Line2D([0],[0], marker='o', color='w', markerfacecolor=FEAT_COLORS['ESM'],
           markersize=4.5, markeredgecolor='white', markeredgewidth=0.3, label='ESM only'),
    Line2D([0],[0], marker='o', color='w', markerfacecolor=FEAT_COLORS['ESM+θ'],
           markersize=4.5, markeredgecolor='white', markeredgewidth=0.3, label='ESM + θ'),
    Line2D([0],[0], marker='*', color='w', markerfacecolor=FEAT_COLORS['Baseline'],
           markersize=6, markeredgecolor='white', markeredgewidth=0.3, label='Baseline'),
]
ax_r.legend(handles=leg_a, loc='lower right', frameon=False, fontsize=4.5)
panel_letter(ax_l, 'a')

# ===================== PANEL B (HERO) =====================
ax_b = fig.add_subplot(gs[0, 1])
ax_b.fill_between([0.415, 0.500], 0.0, 0.80, color='#FFF8E1', alpha=0.2, zorder=0, lw=0)
ax_b.fill_between([0.500, 0.540], 0.0, 0.80, color='#E8F5E9', alpha=0.15, zorder=0, lw=0)
ax_b.axvline(0.50, color=BASELINE_COLOR, ls='-', lw=0.6, alpha=0.3, zorder=1)
ax_b.text(0.422, 0.755, 'Direction without ranking', fontsize=5, color='#BF6B00', style='italic', alpha=0.5)
ax_b.text(0.507, 0.755, 'Ideal', fontsize=5, color='#2E7D32', style='italic', alpha=0.5)

plotted_b = {}
for method in overall.index:
    ft = feat_type(method)
    ht = head_type(method)
    color = FEAT_COLORS[ft]
    marker = HEAD_MARKERS[ht]
    pds = overall.loc[method, 'PDS_cos']
    pdv = overall.loc[method, 'pearson_delta']
    ms = 9 if ft == 'Baseline' else 7.5
    ax_b.plot(pds, pdv, marker=marker, color=color, markersize=ms,
              markeredgecolor='white', markeredgewidth=0.6, zorder=4, linestyle='none')
    plotted_b[method] = (pds, pdv)

for method, (dx, dy, ha) in {'Lasso-esm': (-30, 6, 'right'),
                               'KNN-esm+theta': (6, -8, 'left'),
                               'WT-null': (-8, 7, 'right'),
                               'Gene-mean': (-30, -5, 'right')}.items():
    if method in plotted_b:
        x, y = plotted_b[method]
        ax_b.annotate(clean_name(method), (x, y), xytext=(dx, dy), textcoords='offset points',
                      fontsize=4.5, color='#444', ha=ha, va='center',
                      arrowprops=dict(arrowstyle='-', lw=0.35, color='#BBB', shrinkA=0, shrinkB=2))

ax_b.set_xlabel('PDS (cosine)')
ax_b.set_ylabel('Pearson delta')
ax_b.set_xlim(0.415, 0.540); ax_b.set_ylim(-0.05, 0.78)
ax_b.set_title('All methods cluster in the\ndirection-without-ranking regime', fontsize=6.5, pad=3, color='#555')
panel_letter(ax_b, 'b')

# ===================== PANEL C =====================
ax_c = fig.add_subplot(gs[1, 0])
de_metrics_list = ['DE_overlap', 'DE_LFC_spearman', 'direction_agreement']
de_labels = ['DE Overlap\n(top-50 genes)', 'LFC rank\ncorrelation', 'Direction\nagreement']
wt_null_de = {'DE_overlap': 0.038, 'DE_LFC_spearman': 0.0, 'direction_agreement': 0.0}
non_wt = overall.drop('WT-null', errors='ignore')

for i, metric in enumerate(de_metrics_list):
    vals = non_wt[metric].dropna()
    median_val = vals.median()
    np.random.seed(42 + i)
    jx = np.random.uniform(-0.13, 0.13, len(vals))
    for j, (mn, v) in enumerate(vals.items()):
        ft = feat_type(mn)
        c = FEAT_COLORS[ft]
        ax_c.plot(i + jx[j], v, 'o', color=c, markersize=3.5, alpha=0.65,
                  markeredgecolor='white', markeredgewidth=0.2, zorder=3)
    ax_c.plot([i-0.16, i+0.16], [median_val]*2, '-', color='#333', lw=1.3, zorder=5)
    ax_c.plot(i, wt_null_de[metric], '*', color=BASELINE_COLOR, markersize=7,
              markeredgecolor='white', markeredgewidth=0.3, zorder=5)

ax_c.set_xticks(range(3)); ax_c.set_xticklabels(de_labels, fontsize=5.5)
ax_c.set_ylabel('Score'); ax_c.set_ylim(-0.03, 0.92)
ax_c.set_title('DE programs are partially recovered', fontsize=6.5, pad=3, color='#555')
panel_letter(ax_c, 'c')

# ===================== PANEL D =====================
gs_d = gs[1, 1].subgridspec(2, 2, hspace=0.30, wspace=0.08)
genes_order = ['TP53', 'KRAS', 'GATA1', 'JAK1']
axes_d = [fig.add_subplot(gs_d[r, c]) for r in range(2) for c in range(2)]

for idx, (gene, ax) in enumerate(zip(genes_order, axes_d)):
    g = per_gene_overall[per_gene_overall.gene == gene]
    ax.axvline(0.50, color='#C96144', ls='--', lw=0.4, alpha=0.3)
    for _, row in g.iterrows():
        method = row['method']
        if method == 'WT-null':
            ax.plot(row['PDS_cos'], row['pearson_delta'], '*', color=BASELINE_COLOR, markersize=6,
                    markeredgecolor='white', markeredgewidth=0.3, zorder=5)
            continue
        ft = feat_type(method); ht = head_type(method)
        ax.plot(row['PDS_cos'], row['pearson_delta'], marker=HEAD_MARKERS[ht],
                color=FEAT_COLORS[ft], markersize=4.5, markeredgecolor='white',
                markeredgewidth=0.2, zorder=3, linestyle='none')
    ax.set_title(gene, fontweight='bold', color=GENE_COLORS[gene], fontsize=6.5, pad=2)
    g_nwt = g[g.method != 'WT-null']
    gap = g_nwt['pearson_delta'].mean() - g_nwt['PDS_cos'].mean()
    ax.text(0.95, 0.06, f'Δ = {gap:.2f}', transform=ax.transAxes, fontsize=4.5,
            ha='right', va='bottom', color='#555',
            bbox=dict(boxstyle='round,pad=0.12', facecolor='white', edgecolor='#ddd', alpha=0.85, lw=0.4))
    ax.set_xlim(0.38, 0.58); ax.set_ylim(-0.08, 0.88)
    ax.tick_params(labelsize=4.5)
    if idx in [0, 1]: ax.tick_params(labelbottom=False)
    if idx in [1, 3]: ax.tick_params(labelleft=False)

axes_d[2].set_xlabel('PDS (cosine)', fontsize=5.5)
axes_d[3].set_xlabel('PDS (cosine)', fontsize=5.5)
axes_d[0].set_ylabel('Pearson delta', fontsize=5.5)
axes_d[2].set_ylabel('Pearson delta', fontsize=5.5)
panel_letter(axes_d[0], 'd')

# ===================== PANEL E =====================
gs_e = gs[2, :].subgridspec(2, 1, hspace=0.06, height_ratios=[1, 1])
ax_e1 = fig.add_subplot(gs_e[0])
ax_e2 = fig.add_subplot(gs_e[1], sharex=ax_e1)

tp53_s1 = df[(df.gene=='TP53') & (df.split=='split1') & (df.method=='Ridge-esm')].copy()
tp53_s1 = tp53_s1.sort_values('PDS_cos')
n_show = 4
selected = pd.concat([tp53_s1.head(n_show), tp53_s1.tail(n_show)])
variants_e = selected['variant'].values
xe = np.arange(len(variants_e))

ax_e1.bar(xe, selected['pearson_delta'].values, color='#5185C0', alpha=0.6,
          edgecolor='#5185C0', linewidth=0.4, width=0.65)
ax_e1.set_ylabel('Pearson\ndelta', fontsize=5.5, linespacing=1.1)
ax_e1.set_ylim(0, 0.88)
ax_e1.set_title('Per-variant scores: stable direction, unstable ranking',
                fontsize=6.5, pad=3, color='#555')
ax_e1.spines['bottom'].set_visible(False); ax_e1.tick_params(bottom=False)
plt.setp(ax_e1.get_xticklabels(), visible=False)

ax_e2.bar(xe, selected['PDS_cos'].values, color='#E99D4E', alpha=0.6,
          edgecolor='#E99D4E', linewidth=0.4, width=0.65)
ax_e2.axhline(0.50, color=BASELINE_COLOR, ls='--', lw=0.6, alpha=0.5)
ax_e2.set_ylabel('PDS\n(cosine)', fontsize=5.5, linespacing=1.1)
ax_e2.set_ylim(0, 1.0)
ax_e2.set_xticks(xe)
ax_e2.set_xticklabels(variants_e, rotation=35, ha='right', fontsize=5.5)
ax_e2.axvline(n_show - 0.5, color='#CCC', ls=':', lw=0.5)
ax_e1.axvline(n_show - 0.5, color='#CCC', ls=':', lw=0.5)

panel_letter(ax_e1, 'e')

# Save
outpath = os.path.join(os.path.dirname(__file__), args.out)
os.makedirs(os.path.dirname(outpath), exist_ok=True)
fig.savefig(outpath, dpi=300, bbox_inches='tight', facecolor='white')
print(f"Saved: {outpath}")
plt.close()
