#!/usr/bin/env python3
"""
Figure 3 | Split-half analysis reveals a narrow allele-level measurement window
5 panels: a(concept) b(D_self/D_null strip) c(cells vs ratio) d(un-rankable bars) e(power curve)

Usage:
    python draw_fig3.py [--rank PATH] [--power PATH] [--out fig3_composite.png]
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from fig_config import *

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--rank', default=None, help='rankability CSV path')
parser.add_argument('--power', default=None, help='power curve CSV path')
parser.add_argument('--out', default='../../figures/composites/fig3_composite.png')
args = parser.parse_args()

apply_style()
rank = load_rankability(args.rank)
power = load_power(args.power)

# Prep allele-level per-variant data at NATIVE max-depth (canonical definition).
# native_rankability() takes the deepest split-half bin per perturbation, not the
# shallowest -- this is the locked Round-2 aggregation used everywhere in the paper.
allele_frames = [native_rankability(rank, g) for g in ['TP53','KRAS','GATA1','JAK1']]
allele_pv = pd.concat(allele_frames, ignore_index=True)
allele_pv['ratio_dn'] = allele_pv['D_self'] / allele_pv['D_null'].replace(0, np.nan)
allele_pv = allele_pv.dropna(subset=['ratio_dn'])

# Gene-level un-rankable (raw table kept for other panels)
edist_pca = rank[(rank.metric == 'edist') & (rank.space == 'pca')]

fig = plt.figure(figsize=(7.08, 8.5))
gs = GridSpec(3, 2, figure=fig, height_ratios=[1.0, 0.85, 0.7],
              hspace=0.38, wspace=0.35, left=0.10, right=0.96, top=0.95, bottom=0.06)

# ===================== PANEL A: concept schematic =====================
ax_a = fig.add_subplot(gs[0, 0])
ax_a.set_xlim(0, 10); ax_a.set_ylim(0, 8); ax_a.set_aspect('equal'); ax_a.axis('off')
from matplotlib.patches import Ellipse
wt_x, wt_y = 2.0, 5.5
ax_a.add_patch(Ellipse((wt_x, wt_y), 1.8, 1.2, fc='#E0E0E0', ec='#999', lw=0.8, alpha=0.5))
ax_a.text(wt_x, wt_y, 'WT', fontsize=7, ha='center', va='center', color='#555', fontweight='bold')
v1_x, v1_y = 6.5, 6.5
ax_a.add_patch(Ellipse((v1_x, v1_y), 1.5, 1.0, fc='#5185C0', ec='#3A6DA0', lw=0.8, alpha=0.3))
ax_a.text(v1_x, v1_y, 'half 1', fontsize=5.5, ha='center', va='center', color='#3A6DA0')
v2_x, v2_y = 7.2, 4.5
ax_a.add_patch(Ellipse((v2_x, v2_y), 1.5, 1.0, fc='#5185C0', ec='#3A6DA0', lw=0.8, alpha=0.3))
ax_a.text(v2_x, v2_y, 'half 2', fontsize=5.5, ha='center', va='center', color='#3A6DA0')
ax_a.annotate('', xy=(v2_x-0.3, v2_y+0.55), xytext=(v1_x-0.3, v1_y-0.55),
              arrowprops=dict(arrowstyle='<->', color='#E99D4E', lw=1.5))
ax_a.text(5.5, 5.5, r'$D_{\mathrm{self}}$', fontsize=8, color='#E99D4E', fontweight='bold', ha='center')
mid_x, mid_y = (v1_x + v2_x)/2, (v1_y + v2_y)/2
ax_a.annotate('', xy=(mid_x-0.6, mid_y-0.1), xytext=(wt_x+0.9, wt_y),
              arrowprops=dict(arrowstyle='<->', color='#5185C0', lw=1.5))
ax_a.text(3.8, 6.3, r'$D_{\mathrm{null}}$', fontsize=8, color='#5185C0', fontweight='bold')
ax_a.text(5.0, 1.8, 'Signal window = $D_{\\mathrm{null}} - D_{\\mathrm{self}}$',
          fontsize=6.5, ha='center', color='#333',
          bbox=dict(boxstyle='round,pad=0.3', fc='#F5F5F5', ec='#DDD', lw=0.5))
ax_a.set_title('Split-half measurement window', fontsize=7, pad=2, color='#555')
panel_letter(ax_a, 'a')

# ===================== PANEL B: strip plot =====================
ax_b = fig.add_subplot(gs[0, 1])
genes = ['TP53', 'KRAS', 'GATA1', 'JAK1']
_rng_b = np.random.default_rng(0)
for i, gene in enumerate(genes):
    raw = allele_pv[allele_pv.dataset == gene]['ratio_dn'].dropna().values
    sub = np.clip(raw, 0, 2.5)
    jx = _rng_b.uniform(-0.18, 0.18, len(sub))
    ax_b.scatter(i + jx, sub, c=GENE_COLORS[gene], s=6, alpha=0.35,
                 edgecolor='white', linewidth=0.15, zorder=3)
    med = np.median(sub)
    # 95% bootstrap CI on the MEAN ratio (unclipped), drawn as a vertical whisker
    boots = [np.mean(_rng_b.choice(raw, len(raw), replace=True)) for _ in range(2000)]
    lo, hi = np.percentile(boots, [2.5, 97.5])
    ax_b.plot([i, i], [lo, hi], '-', color='#333', lw=1.1, alpha=0.8, zorder=5)
    ax_b.plot([i-0.2, i+0.2], [med, med], '-', color='#333', lw=1.3, zorder=6)
    ax_b.text(i+0.25, med, f'{med:.2f}', fontsize=5, va='center', color='#333')
ax_b.axhline(1.0, color=BASELINE_COLOR, ls='--', lw=0.7, alpha=0.5, zorder=1)
ax_b.axhspan(0.85, 1.15, color='#FFF3E0', alpha=0.25, zorder=0)
ax_b.axhspan(0, 0.5, color='#E8F5E9', alpha=0.15, zorder=0)
ax_b.set_xticks(range(4)); ax_b.set_xticklabels(genes, fontsize=6, fontweight='bold')
for tick, gene in zip(ax_b.xaxis.get_ticklabels(), genes):
    tick.set_color(GENE_COLORS[gene])
ax_b.set_ylabel('D_self / D_null'); ax_b.set_ylim(-0.05, 2.2)
ax_b.set_title('Replicate noise approaches variant signal', fontsize=7, pad=2, color='#555')
panel_letter(ax_b, 'b')

# ===================== PANEL C: cells vs ratio =====================
ax_c = fig.add_subplot(gs[1, 0])
for gene in genes:
    sub = allele_pv[allele_pv.dataset == gene]
    ax_c.scatter(sub.n_cells, sub.ratio_dn.clip(0, 2.5), c=GENE_COLORS[gene], s=7, alpha=0.4,
                 edgecolor='white', linewidth=0.15, zorder=3, label=gene)
ax_c.axhline(1.0, color=BASELINE_COLOR, ls='--', lw=0.6, alpha=0.4)
ax_c.axvline(200, color='#999', ls=':', lw=0.5)
ax_c.set_xscale('log'); ax_c.set_xlabel('Cells per variant')
ax_c.set_ylabel('D_self / D_null'); ax_c.set_ylim(-0.05, 2.3)
ax_c.set_title('Depth and effect size jointly determine\nthe measurement window', fontsize=7, pad=2, color='#555')
ax_c.legend(frameon=False, fontsize=5, loc='upper right')
panel_letter(ax_c, 'c')

# ===================== PANEL D: un-rankable bars =====================
ax_d = fig.add_subplot(gs[1, 1])
datasets_order = ['TP53','KRAS','GATA1','JAK1','Replogle','Norman','Adamson']
colors_d = [GENE_COLORS.get(d, '#888') for d in ['TP53','KRAS','GATA1','JAK1']] + ['#666','#888','#AAA']
fracs, ns = [], []
for ds in datasets_order:
    fr, n = unrankable_fraction(rank, ds)   # canonical native-depth
    fracs.append(fr); ns.append(n)
x = np.arange(len(datasets_order))
ax_d.bar(x, [f*100 for f in fracs], width=0.6, color=colors_d, alpha=0.7, edgecolor='white', lw=0.5)
ax_d.axvline(3.5, color='#CCC', ls=':', lw=0.5)
for i, (f, n) in enumerate(zip(fracs, ns)):
    ax_d.text(i, f*100 + 1.5, f'{f*100:.0f}%', fontsize=4.5, ha='center', color='#333')
ax_d.set_xticks(x)
ax_d.set_xticklabels(['TP53','KRAS','GATA1','JAK1','Replogle','Norman','Adamson'],
                      fontsize=5.5, rotation=25, ha='right')
ax_d.set_ylabel('Un-rankable (%)'); ax_d.set_ylim(0, 115)
ax_d.set_title('Detection floor extends to gene-level benchmarks', fontsize=7, pad=2, color='#555')
panel_letter(ax_d, 'd')

# ===================== PANEL E: power curve =====================
ax_e = fig.add_subplot(gs[2, :])
pca = power[power.space == 'pca50'].copy()
for gene in genes:
    g = pca[pca.gene == gene].sort_values('n_sub')
    ax_e.plot(g.n_sub, g.frac_detectable * 100, '-o', color=GENE_COLORS[gene],
              markersize=4.5, markeredgecolor='white', markeredgewidth=0.3, lw=1.3, label=gene)
ax_e.axhline(50, color='#999', ls=':', lw=0.4)
ax_e.axvline(200, color='#BBB', ls=':', lw=0.4)
ax_e.set_xlabel('Cells per variant (subsampled)')
ax_e.set_ylabel('Variants detectable (%)'); ax_e.set_xlim(30, 350); ax_e.set_ylim(35, 108)
ax_e.set_title('Detection rate rises with depth, gene-dependently', fontsize=7, pad=2, color='#555')
ax_e.legend(frameon=False, fontsize=5, loc='lower right')
panel_letter(ax_e, 'e')

outpath = os.path.join(os.path.dirname(__file__), args.out)
os.makedirs(os.path.dirname(outpath), exist_ok=True)
fig.savefig(outpath, dpi=300, bbox_inches='tight', facecolor='white')
print(f"Saved: {outpath}")
plt.close()
