#!/usr/bin/env python3
"""Plot depth-vs-assembly: contigs, N50, avg length (1×3 layout)."""
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

FIGS = Path('benchmark/figures')
FIGS.mkdir(exist_ok=True)

rows = [l.strip().split('\t') for l in open('/tmp/depth_v2/depth_summary.tsv') if l.strip()][1:]
depths  = [float(r[0].rstrip('x')) for r in rows]
contigs = [int(r[1]) for r in rows]
n50s    = [int(r[5]) for r in rows]
avgs    = [float(r[6]) for r in rows]

genome_ct = 7690
genome_n50 = 311
genome_avg = 200

plt.rcdefaults()
plt.rcParams.update({
    'font.family': 'sans-serif', 'font.size': 9,
    'axes.labelsize': 9, 'axes.titlesize': 10, 'axes.titleweight': 'bold',
    'xtick.labelsize': 8, 'ytick.labelsize': 8,
    'axes.linewidth': 0.5, 'xtick.major.width': 0.5, 'ytick.major.width': 0.5,
    'pdf.fonttype': 42, 'ps.fonttype': 42,
    'axes.edgecolor': '#AAAAAA', 'xtick.color': '#555555', 'ytick.color': '#555555',
})

fig, axes = plt.subplots(1, 3, figsize=(11, 3.5))
color = '#2B6CB0'

titles = ['Contigs', 'N50', 'Avg length']
ylabels = ['Contigs', 'N50 (bp)', 'Avg length (bp)']
vals_list = [contigs, n50s, avgs]
hl_list = [genome_ct, genome_n50, genome_avg]
hl_labels = ['genome', 'genome', 'genome']

for i, (ax, vals, tt, yl, hl) in enumerate(zip(axes, vals_list, titles, ylabels, hl_list)):
    ax.plot(depths, vals, '-o', color=color, linewidth=1.5, markersize=5, markerfacecolor='white')
    ax.axhline(hl, color='#CC3333', linewidth=0.8, linestyle='--')
    ax.text(depths[-1] * 1.6, hl * 0.95, 'genome', fontsize=7, color='#CC3333', va='bottom')
    ax.set_title(f'{"abc"[i]}. {tt}', fontsize=10, fontweight='bold', pad=5)
    ax.set_xlabel('Coverage', fontsize=9, color='#555')
    ax.set_ylabel(yl, fontsize=9, color='#555')
    ax.set_xscale('log')
    ax.set_xticks(depths)
    ax.set_xticklabels([f'{d:g}x' for d in depths])
    ax.tick_params(labelsize=8, colors='#555')
    ax.spines[['top','right']].set_visible(False)
    ax.spines[['bottom','left']].set_color('#CCCCCC')

plt.tight_layout(pad=2)
fig.savefig(FIGS / 'depth_experiment.pdf', bbox_inches='tight', dpi=300)
plt.close()
print(f'Saved: {FIGS}/depth_experiment.pdf')
