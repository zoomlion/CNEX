#!/usr/bin/env python3
"""CNEX query-lib 子采样图 (figures/querylib_recall.pdf)

(a) 相对 recall vs 库比例 (bootstrap 95% CI error bar)
(b) 逐元素恢复率 vs 库比例
(c) precision vs 库比例
"""
import sys
from pathlib import Path
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

BASE  = Path(__file__).resolve().parent
STATS = BASE / 'stats'
FIGS  = BASE / 'figures'
FIGS.mkdir(exist_ok=True)

SPECIES = ['Homo_sapiens', 'Mus_musculus', 'Gallus_gallus', 'Scyliorhinus_canicula', 'Callorhinchus_milii']
SHORT   = {'Homo_sapiens': 'Human', 'Mus_musculus': 'Mouse', 'Gallus_gallus': 'Chicken',
           'Scyliorhinus_canicula': 'Catshark', 'Callorhinchus_milii': 'Callorhinchus'}
MARKERS = {'Homo_sapiens': 'o', 'Mus_musculus': 's', 'Gallus_gallus': '^',
           'Scyliorhinus_canicula': 'D', 'Callorhinchus_milii': 'v'}


def load_summary():
    rows_by_sp = {}
    p = STATS / 'querylib_summary.tsv'
    if not p.exists():
        print(f'WARN: {p} not found'); return rows_by_sp
    for line in open(p):
        if line.startswith('species'):
            continue
        f = line.strip().split('\t')
        if len(f) < 13:
            continue
        rows_by_sp.setdefault(f[0], []).append({
            'frac': float(f[1]), 'rel': float(f[5]), 'ci_low': float(f[6]),
            'ci_high': float(f[7]), 'abs': float(f[8]), 'ceil': float(f[9]),
            'prec': float(f[10]), 'rec': float(f[11]), 'wall': float(f[12])})
    for sp in rows_by_sp:
        rows_by_sp[sp].sort(key=lambda r: r['frac'])
    return rows_by_sp


def main():
    targets = sys.argv[1:] or SPECIES
    data = load_summary()
    species = [s for s in targets if s in data] or list(data)
    if not species:
        print('No data in stats/querylib_summary.tsv — run `bash benchmark.sh querylib` first.')
        sys.exit(1)

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

    colors = plt.cm.tab10(np.linspace(0, 1, len(species)))
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.4))
    fig.subplots_adjust(hspace=0.45, wspace=0.35)
    titles = ['a. Relative recall', 'b. Element recovery', 'c. Precision']
    ylabels = ['Relative recall (%)', 'Recovery (%)', 'Precision (%)']

    for ax, (t, yl) in zip(axes.flat, zip(titles, ylabels)):
        ax.set_title(t, fontsize=10, fontweight='bold', pad=5)
        ax.set_xlabel('Query library size (% of blocks)', fontsize=9, color='#555')
        ax.set_ylabel(yl, fontsize=9, color='#555')
        ax.tick_params(labelsize=8, colors='#555')
        ax.spines[['top', 'right']].set_visible(False)
        ax.spines[['bottom', 'left']].set_color('#CCCCCC')
        ax.set_xticks([10, 20, 40, 60, 80, 100])
        ax.set_xlim(5, 105)
        ax.set_ylim(0, 105)

    for sp, c in zip(species, colors):
        rows = data[sp]
        xs = [r['frac'] * 100 for r in rows]
        label = SHORT.get(sp, sp)
        mk = MARKERS.get(sp, 'o')
        axes[0].errorbar(xs, [r['rel'] for r in rows],
                         yerr=[[r['rel'] - r['ci_low'] for r in rows],
                               [r['ci_high'] - r['rel'] for r in rows]],
                         marker=mk, color=c, label=label, linewidth=1.6, capsize=3, ms=4)
        axes[1].plot(xs, [r['rec'] for r in rows], marker=mk, color=c,
                     label=label, linewidth=1.6, ms=4)
        axes[2].plot(xs, [r['prec'] for r in rows], marker=mk, color=c,
                     label=label, linewidth=1.6, ms=4)

    axes[0].legend(fontsize=8, frameon=False, loc='lower right')

    fig.savefig(FIGS / 'querylib_recall.pdf', bbox_inches='tight', dpi=300)
    plt.close()
    print(f'Saved: figures/querylib_recall.pdf')

    print('\nQuery-lib sweep (CNEX):')
    for sp in species:
        for r in data[sp]:
            print(f'  {sp:22s} frac={r["frac"]*100:4.0f}%  rel={r["rel"]:5.1f}% '
                  f'abs={r["abs"]:5.1f}%  prec={r["prec"]:5.1f}%  rec={r["rec"]:5.1f}%  wall={r["wall"]:.0f}s')


if __name__ == '__main__':
    main()
