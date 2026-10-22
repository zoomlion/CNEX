#!/usr/bin/env python3
"""多物种 benchmark 统计图 (genome mode).

(a) Recall  (b) Precision  (c) Wall time (human)  (d) Peak memory (human)
"""
import random
from pathlib import Path

import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

BASE = Path(__file__).resolve().parent
DATA = BASE / 'data'
RES = BASE / 'results'
STATS = BASE / 'stats'
FIGS = BASE / 'figures'
FIGS.mkdir(exist_ok=True)

SPECIES = ['Homo_sapiens', 'Mus_musculus', 'Gallus_gallus', 'Scyliorhinus_canicula', 'Callorhinchus_milii']
SHORT = {'Homo_sapiens': 'Human', 'Mus_musculus': 'Mouse', 'Gallus_gallus': 'Chicken',
         'Scyliorhinus_canicula': 'Catshark', 'Callorhinchus_milii': 'Callorhinchus'}
METHODS = ['CNEX', 'MMseqs2', 'BLASTN']
COLORS = {'CNEX': '#2B6CB0', 'MMseqs2': '#D97706', 'BLASTN': '#38A169'}
COLORS_LIGHT = {'CNEX': '#7DB8E8', 'MMseqs2': '#E8A84E', 'BLASTN': '#6FC98A'}
SPECIES_MARKERS = {'Homo_sapiens': 'o', 'Mus_musculus': 's', 'Gallus_gallus': '^',
                   'Scyliorhinus_canicula': 'D', 'Callorhinchus_milii': 'v'}


def load_outcomes():
    """outcomes.tsv -> {species: {method: {'correct': {ele: bool}}}}"""
    data = {s: {m: {'correct': {}} for m in METHODS} for s in SPECIES}
    p = STATS / 'outcomes.tsv'
    if not p.exists():
        print(f'WARN: {p} not found'); return data
    for line in open(p):
        if line.startswith('species'):
            continue
        f = line.strip().split('\t')
        if len(f) < 8:
            continue
        sp, ele = f[0], int(f[1])
        c_ok, m_ok, b_ok = int(f[3]), int(f[5]), int(f[7])
        for m, ok in zip(METHODS, [c_ok, m_ok, b_ok]):
            data[sp][m]['correct'][ele] = bool(ok)
    return data


def load_summary():
    """stats/summary.tsv -> {species: {method: {'recall','precision'}}}"""
    data = {s: {} for s in SPECIES}
    p = STATS / 'summary.tsv'
    if not p.exists():
        print(f'WARN: {p} not found'); return data
    for line in open(p):
        if line.startswith('species'):
            continue
        f = line.strip().split('\t')
        if len(f) < 6:
            continue
        sp, method = f[0], f[1]
        if sp in data and method in METHODS:
            try:
                data[sp][method] = {'recall': float(f[4]), 'precision': float(f[5])}
            except ValueError:
                pass
    return data


def parse_t(fn):
    d = {}
    p = STATS / 'Homo_sapiens' / fn
    if p.exists():
        for l in open(p):
            if 'User time' in l:
                d['user'] = float(l.split(':')[1].strip().split()[0])
            elif 'Elapsed' in l:
                r = l.split()[-1].strip(); pts = r.split(':')
                if len(pts) == 3:
                    d['wall_s'] = int(pts[0]) * 3600 + int(pts[1]) * 60 + float(pts[2])
                elif len(pts) == 2:
                    d['wall_s'] = int(pts[0]) * 60 + float(pts[1])
            elif 'Maximum resident' in l:
                d['mem_mb'] = int(l.split(':')[1].strip().split()[0]) // 1024
    return d


def main():
    data = load_outcomes()
    summary = load_summary()

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
    fig, axs = plt.subplots(1, 4, figsize=(12.5, 3.4), constrained_layout=True)
    fig.set_constrained_layout_pads(w_pad=0.06, h_pad=0.06, wspace=0.04, hspace=0.04)

    def style(ax, letter, title, ylab):
        ax.set_title(f'{letter}. {title}')
        ax.set_ylabel(ylab)
        ax.tick_params(labelsize=8, colors='#555')
        ax.spines[['top', 'right']].set_visible(False)
        ax.spines[['bottom', 'left']].set_color('#CCCCCC')
        ax.grid(True, alpha=0.18, linestyle='--', linewidth=0.5)

    def draw_tool_boxes(ax, metric, letter, title):
        rng = random.Random(42)
        series = [[summary[s][m][metric] for s in SPECIES] for m in METHODS]
        bp = ax.boxplot(series, positions=[0, 1, 2], widths=0.5, patch_artist=True,
                        whis=(0, 100), showfliers=False, manage_ticks=False)
        for i, m in enumerate(METHODS):
            bp['boxes'][i].set_facecolor(COLORS_LIGHT[m])
            bp['boxes'][i].set_edgecolor(COLORS[m])
            bp['boxes'][i].set_linewidth(0.7)
            bp['medians'][i].set_color(COLORS[m])
            bp['medians'][i].set_linewidth(1.1)
            for key in ('whiskers', 'caps'):
                for w in range(2):
                    bp[key][i * 2 + w].set_color(COLORS[m])
                    bp[key][i * 2 + w].set_linewidth(0.5)
            for s in SPECIES:
                ax.scatter(i + rng.uniform(-0.13, 0.13), summary[s][m][metric],
                           marker=SPECIES_MARKERS[s], s=22, color=COLORS[m],
                           edgecolors='white', linewidths=0.4, zorder=3, alpha=0.9)
        ax.set_xticks([0, 1, 2])
        ax.set_xticklabels(METHODS, fontsize=8, color='#555')
        ax.set_ylim(0, 100)
        style(ax, letter, title, 'Recall (%)' if metric == 'recall' else 'Precision (%)')

    # (a) Recall
    draw_tool_boxes(axs[0], 'recall', 'a', 'Recall')
    # (b) Precision
    draw_tool_boxes(axs[1], 'precision', 'b', 'Precision')

    # (c) Wall time (human)
    cnex_v = parse_t('time_cnex_validate.txt').get('wall_s', 0)
    cnex_a = parse_t('time_cnex_assemble.txt').get('wall_s', 0)
    mm_db = parse_t('time_mmseqs_createdb.txt').get('wall_s', 0)
    mm_sr = parse_t('time_mmseqs_search.txt').get('wall_s', 0)
    bl_db = parse_t('time_makeblastdb.txt').get('wall_s', 0)
    bl_sr = parse_t('time_blastn_w15.txt').get('wall_s', 0)
    ax = axs[2]
    bottom = [cnex_v, mm_db, bl_db]
    top = [cnex_a, mm_sr, bl_sr]
    ax.bar(METHODS, bottom, color=[COLORS[m] for m in METHODS], width=0.45,
           edgecolor='white', linewidth=0.5)
    ax.bar(METHODS, top, bottom=bottom, color=[COLORS_LIGHT[m] for m in METHODS],
           width=0.45, edgecolor='white', linewidth=0.5)
    for j in range(3):
        v = bottom[j] + top[j]
        ax.text(j, v, f'{v / 60:.1f} min', ha='center', va='bottom', fontsize=7, color='#555')
    ax.set_ylim(0, max(b + t for b, t in zip(bottom, top)) * 1.18)
    style(ax, 'c', 'Wall time (human)', 'Wall time (s)')

    # (d) Peak memory (human)
    ax = axs[3]
    mems = [parse_t('time_cnex_validate.txt').get('mem_mb', 0) / 1024.0,
            parse_t('time_mmseqs_search.txt').get('mem_mb', 0) / 1024.0,
            parse_t('time_blastn_w15.txt').get('mem_mb', 0) / 1024.0]
    bars = ax.bar(METHODS, mems, color=[COLORS[m] for m in METHODS], width=0.45,
                  edgecolor='white', linewidth=0.5)
    for b, v in zip(bars, mems):
        ax.text(b.get_x() + b.get_width() / 2, v, f'{v:.1f} GB',
                ha='center', va='bottom', fontsize=7, color='#555')
    ax.set_ylim(0, max(mems) * 1.18)
    style(ax, 'd', 'Peak memory (human)', 'Peak memory (GB)')

    # single shared legend (tools + species markers)
    handles = [Patch(facecolor=COLORS_LIGHT[m], edgecolor=COLORS[m], label=m) for m in METHODS]
    sp_handles = [Line2D([], [], marker=SPECIES_MARKERS[s], linestyle='None',
                         markerfacecolor='#888', markeredgecolor='none', label=SHORT[s])
                  for s in SPECIES]
    fig.legend(handles=handles + sp_handles, loc='outside upper center', ncol=8, frameon=False,
               fontsize=8, handletextpad=0.4, columnspacing=1.2)

    fig.savefig(FIGS / 'benchmark_multi.pdf', bbox_inches='tight', dpi=300)
    plt.close()
    print('Saved: figures/benchmark_multi.pdf')

    print('\nSummary (recall / precision):')
    for s in SPECIES:
        line = f'  {s:22s}'
        for m in METHODS:
            n = len(data[s][m]['correct'])
            ok = sum(1 for v in data[s][m]['correct'].values() if v)
            line += f'  {m:7s} {ok / n * 100:5.1f}%'
        print(line)


if __name__ == '__main__':
    main()
