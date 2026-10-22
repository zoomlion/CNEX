#!/usr/bin/env python3
"""Supplementary MSA species-count ablation figure (human target, 5x reads).

Two stacked panels sharing the x axis:
  a. recall vs number of MSA species (far-to-near)
  b. species composition (19 species x 4/8/12/16/19)

Reads benchmark/stats/msa_ablation.tsv (from msa_ablation.py) and writes
benchmark/figures/msa_ablation_supp.pdf.
"""
from pathlib import Path

import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from msa_ablation import SPECIES

HERE = Path(__file__).resolve().parent
TSV = HERE / 'stats' / 'msa_ablation.tsv'
OUT = HERE / 'figures' / 'msa_ablation_supp.pdf'
C_CNEX = '#2B6CB0'

plt.rcdefaults()
plt.rcParams.update({
    'font.family': 'sans-serif', 'font.size': 9,
    'axes.labelsize': 9, 'axes.titlesize': 10, 'axes.titleweight': 'bold',
    'xtick.labelsize': 8, 'ytick.labelsize': 8, 'legend.fontsize': 8,
    'axes.linewidth': 0.5, 'xtick.major.width': 0.5, 'ytick.major.width': 0.5,
    'xtick.major.pad': 3, 'ytick.major.pad': 3,
    'pdf.fonttype': 42, 'ps.fonttype': 42,
    'axes.edgecolor': '#AAAAAA', 'xtick.color': '#555555', 'ytick.color': '#555555',
})


def load_tsv():
    rows = []
    for line in open(TSV):
        p = line.rstrip('\n').split('\t')
        if p[0] == 'species':
            continue
        rows.append({
            'species': int(p[0]),
            'table_mb': float(p[2]),
            'recall': float(p[3]),
            'precision': float(p[4]),
        })
    return sorted(rows, key=lambda r: r['species'])


def style(a):
    a.spines['top'].set_visible(False)
    a.spines['right'].set_visible(False)
    a.spines['bottom'].set_color('#CCCCCC')
    a.spines['left'].set_color('#CCCCCC')
    a.grid(True, alpha=0.18, linestyle='--', linewidth=0.5)
    a.tick_params(labelsize=8, colors='#555555')


def main():
    rows = load_tsv()
    levels = [r['species'] for r in rows]
    rec = [r['recall'] for r in rows]
    prec = [r['precision'] for r in rows]
    mb = [r['table_mb'] for r in rows]
    x = np.arange(len(levels))

    fig = plt.figure(figsize=(4.9, 5.6))
    gs = fig.add_gridspec(2, 1, height_ratios=[0.95, 1.35], hspace=0.15)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[1, 0], sharex=ax_a)

    # a. recall vs species count
    ax_a.plot(x, rec, 'o-', color=C_CNEX, linewidth=1.1, markersize=4,
              markerfacecolor='white')
    for xi, yi in zip(x, rec):
        ax_a.text(xi, yi + 2.0, f'{yi:.1f}', ha='center', fontsize=7.5)
    ax_a.set_xticks(x)
    ax_a.set_xticklabels([])
    ax_a.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
    ax_a.set_xlim(-0.5, len(levels) - 0.5)
    ax_a.set_ylabel('Recall (%)')
    ax_a.set_ylim(0, 100)
    ax_a.set_yticks([20, 40, 60, 80, 100])
    ax_a.set_title('S3a. Recall vs MSA species', loc='left')
    style(ax_a)
    ax_a.grid(False)
    ax_a.grid(True, axis='y', alpha=0.18, linestyle='--', linewidth=0.5)
    ax_a.spines['bottom'].set_visible(False)

    # b. species composition, drawn as a light dot matrix
    n = len(SPECIES)
    xx, yy = np.meshgrid(x, np.arange(n))
    included = np.array([[i < L for L in levels] for i in range(n)])
    ax_b.scatter(xx[~included], yy[~included], s=22, facecolors='none',
                 edgecolors='#D0D0D0', linewidths=0.7, zorder=2)
    ax_b.scatter(xx[included], yy[included], s=22, color=C_CNEX,
                 edgecolors='white', linewidths=0.5, zorder=3)
    ax_b.set_xlim(-0.5, len(levels) - 0.5)
    ax_b.set_ylim(n - 0.5, -0.5)
    ax_b.set_xticks(x)
    ax_b.set_xticklabels([str(v) for v in levels], fontsize=8)
    ax_b.set_yticks(np.arange(n))
    ax_b.set_yticklabels([s.replace('_', ' ') for s in SPECIES], fontsize=4.9)
    ax_b.set_xlabel('MSA species added (far → near)')
    ax_b.set_title('S3b. Species composition', loc='left', pad=6)
    ax_b.tick_params(length=0)
    ax_b.spines['top'].set_visible(False)
    ax_b.spines['right'].set_visible(False)
    ax_b.spines['bottom'].set_color('#CCCCCC')
    ax_b.spines['left'].set_color('#CCCCCC')
    ax_b.grid(False)

    fig.text(0.5, 0.005,
             'Human target, 5x reads. Filled = included; open = absent.\n'
             '4-species point = distant-only baseline (min_c = 4, k = 13).\n'
             f'Precision {min(prec):.1f}–{max(prec):.1f}%; '
             f'table {min(mb):.0f}–{max(mb):.0f} MB.',
             ha='center', fontsize=6.0, color='#777777', linespacing=1.35)
    fig.tight_layout(pad=0.5, rect=[0, 0.18, 1, 1])
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches='tight', dpi=300)
    plt.close(fig)
    print(f'Saved: {OUT}')


if __name__ == '__main__':
    main()
