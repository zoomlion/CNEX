#!/usr/bin/env python3
"""多物种统计图(旧版样式复刻):
(a) Recall boxplot (bootstrap)   (b) Precision boxplot (bootstrap)
(c) Wall time (human)            (d) Peak memory (human)
(e) Recall vs CNE conservation (human, 原版曲线)
"""
import random
from collections import Counter
from pathlib import Path
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

BASE  = Path(__file__).resolve().parent
DATA  = BASE / 'data'
RES   = BASE / 'results'
STATS = BASE / 'stats'
FIGS  = BASE / 'figures'
FIGS.mkdir(exist_ok=True)

SPECIES = ['Homo_sapiens', 'Mus_musculus', 'Gallus_gallus', 'Scyliorhinus_canicula', 'Callorhinchus_milii']
SHORT   = {'Homo_sapiens': 'Human', 'Mus_musculus': 'Mouse', 'Gallus_gallus': 'Chicken',
           'Scyliorhinus_canicula': 'Catshark', 'Callorhinchus_milii': 'Callorhinchus'}
METHODS = ['CNEX', 'MMseqs2', 'BLASTN']
COLORS  = {'CNEX': '#2B6CB0', 'MMseqs2': '#D97706', 'BLASTN': '#38A169'}
COLORS_LIGHT = {'CNEX': '#7DB8E8', 'MMseqs2': '#E8A84E', 'BLASTN': '#6FC98A'}


def load_outcomes():
    """outcomes.tsv -> {species: {method: {'ele': bool_correct, 'hit': bool}}}"""
    data = {s: {m: {'correct': {}, 'hit': {}} for m in METHODS} for s in SPECIES}
    p = STATS / 'outcomes.tsv'
    if not p.exists():
        print(f'WARN: {p} not found'); return data
    for line in open(p):
        if line.startswith('species'): continue
        f = line.strip().split('\t')
        if len(f) < 8: continue
        sp, ele = f[0], int(f[1])
        c_hit, c_ok, m_hit, m_ok, b_hit, b_ok = (int(x) for x in f[2:8])
        for m, (hit, ok) in zip(METHODS, [(c_hit, c_ok), (m_hit, m_ok), (b_hit, b_ok)]):
            data[sp][m]['hit'][ele] = bool(hit)
            data[sp][m]['correct'][ele] = bool(ok)
    return data


def load_summary():
    """stats/summary.tsv -> {species: {method: {'recall', 'precision'}}}"""
    data = {s: {} for s in SPECIES}
    p = STATS / 'summary.tsv'
    if not p.exists():
        print(f'WARN: {p} not found'); return data
    for line in open(p):
        if line.startswith('species'):
            continue
        f = line.strip().split('\t')
        if len(f) < 6: continue
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
                if len(pts) == 3: d['wall_s'] = int(pts[0]) * 3600 + int(pts[1]) * 60 + float(pts[2])
                elif len(pts) == 2: d['wall_s'] = int(pts[0]) * 60 + float(pts[1])
            elif 'Maximum resident' in l:
                d['mem_mb'] = int(l.split(':')[1].strip().split()[0]) // 1024
    return d


# ─── 面板 e: 保守度(原版逻辑, blocks/Homo_sapiens.fa) ───
def calc_block_conservation(block_text):
    """Per-column consensus conservation (skip ref.*, gaps)."""
    lines = block_text.strip().split('\n')
    hdrs = []; seqs = []
    for l in lines:
        if l.startswith('>'): hdrs.append(l); seqs.append('')
        elif seqs: seqs[-1] += l.strip().upper()
    species_seqs = [s for h, s in zip(hdrs, seqs) if not h.startswith('>ref')]
    if len(species_seqs) < 2: return 0.0
    aln_len = len(species_seqs[0])
    scores = []
    for col in range(aln_len):
        bases = [s[col] for s in species_seqs if s[col] != '-']
        if len(bases) < 2: continue
        cnt = Counter(bases)
        scores.append(cnt.most_common(1)[0][1] / len(bases))
    return np.mean(scores) if scores else 0.0


def conservation_curve(data):
    """人类 recall vs 保守度曲线数据"""
    blocks = open(DATA / 'blocks' / 'Homo_sapiens.fa').read().strip().split('###')
    block_conservation = {}
    for i, b in enumerate(blocks):
        b = b.strip()
        if b:
            block_conservation[i] = calc_block_conservation(b)
    truth = data['Homo_sapiens']
    bid_to_cons = {}
    for ele in truth['CNEX']['correct']:
        bid_to_cons[f'block_{ele:05d}'] = block_conservation.get(ele, 0.0)
    thresholds = np.arange(0.50, 1.01, 0.02)
    curve_data = {}
    for m in METHODS:
        mdata = truth[m]['correct']
        curve = []
        for t in thresholds:
            subset = [bid for bid in bid_to_cons if bid_to_cons[bid] >= t]
            if subset:
                correct = sum(1 for bid in subset if mdata.get(int(bid.split('_')[1]), False))
                curve.append(correct / len(subset) * 100)
            else:
                curve.append(0)
        curve_data[m] = curve
    return thresholds, curve_data


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
    fig = plt.figure(figsize=(9, 7))
    gs = fig.add_gridspec(2, 4, height_ratios=[1, 1.15], hspace=0.35, wspace=0.35)

    x_labels = ['CNEX', 'MMseqs2', 'BLASTN']
    species_markers = {'Homo_sapiens': 'o', 'Mus_musculus': 's', 'Gallus_gallus': '^',
                       'Scyliorhinus_canicula': 'D', 'Callorhinchus_milii': 'v'}

    def draw_tool_boxes(ax, summary, metric, title, letter):
        """每工具一个箱(5 物种的观测值) + 物种散点叠加"""
        from matplotlib.lines import Line2D
        rng = random.Random(42)
        series = []
        for m in METHODS:
            series.append([summary[s][m][metric] for s in SPECIES])
        bp = ax.boxplot(series, positions=[0, 1, 2], widths=0.5, patch_artist=True,
                        whis=(0, 100), showfliers=False, manage_ticks=False)
        for i, m in enumerate(METHODS):
            bp['boxes'][i].set_facecolor(COLORS_LIGHT[m])
            bp['boxes'][i].set_edgecolor(COLORS[m])
            bp['boxes'][i].set_linewidth(0.8)
            bp['medians'][i].set_color(COLORS[m])
            bp['medians'][i].set_linewidth(1.2)
            for key in ('whiskers', 'caps'):
                for w in range(2):
                    bp[key][i * 2 + w].set_color(COLORS[m])
                    bp[key][i * 2 + w].set_linewidth(0.6)
            for j, s in enumerate(SPECIES):
                x = i + rng.uniform(-0.14, 0.14)
                ax.scatter(x, summary[s][m][metric], marker=species_markers[s],
                           s=26, color=COLORS[m], edgecolors='white', linewidths=0.4,
                           zorder=3, alpha=0.9)
        ax.set_title(f'{letter}. {title}', fontsize=10, fontweight='bold', pad=5)
        ax.set_ylabel('Recall (%)' if metric == 'recall' else 'Precision (%)',
                      fontsize=9, color='#555')
        ax.set_xticks([0, 1, 2])
        ax.set_xticklabels(x_labels, fontsize=8, color='#555')
        ax.set_ylim(0, 100)
        ax.tick_params(labelsize=8, colors='#555')
        ax.spines[['top', 'right']].set_visible(False)
        ax.spines[['bottom', 'left']].set_color('#CCCCCC')
        species_handles = [Line2D([], [], marker=species_markers[s], linestyle='None',
                                  markerfacecolor='#888', markeredgecolor='none',
                                  label=SHORT[s]) for s in SPECIES]
        ax.legend(handles=species_handles, loc='lower right', fontsize=6, frameon=False,
                  labelspacing=0.3, title='Target', title_fontsize=6)

    # ── (a) Recall per tool (5 species boxed) ──
    ax = fig.add_subplot(gs[0, 0])
    draw_tool_boxes(ax, summary, 'recall', 'Recall', 'a')
    # ── (b) Precision per tool (5 species boxed) ──
    ax = fig.add_subplot(gs[0, 1])
    draw_tool_boxes(ax, summary, 'precision', 'Precision', 'b')

    # legend (工具)
    from matplotlib.patches import Patch
    handles = [Patch(facecolor=COLORS_LIGHT[m], edgecolor=COLORS[m], label=m) for m in METHODS]
    fig.legend(handles=handles, loc='upper center', ncol=3, frameon=False,
               fontsize=8, bbox_to_anchor=(0.5, 0.99))

    # ── (c) Wall time (human, 原版样式) ──
    cnex_v = parse_t('time_cnex_validate.txt').get('wall_s', 0)
    cnex_a = parse_t('time_cnex_assemble.txt').get('wall_s', 0)
    mm_db = parse_t('time_mmseqs_createdb.txt').get('wall_s', 0)
    mm_sr = parse_t('time_mmseqs_search.txt').get('wall_s', 0)
    bl_db = parse_t('time_makeblastdb.txt').get('wall_s', 0)
    bl_sr = parse_t('time_blastn_w15.txt').get('wall_s', 0)
    ax = fig.add_subplot(gs[0, 2])
    labels = ['CNEX', 'MMseqs2', 'BLASTN']
    time_bottom = [cnex_v, mm_db, bl_db]
    time_top = [cnex_a, mm_sr, bl_sr]
    ax.bar(labels, time_bottom, color=[COLORS[m] for m in METHODS], width=0.4,
           edgecolor='white', linewidth=0.5)
    ax.bar(labels, time_top, bottom=time_bottom,
           color=[COLORS_LIGHT[m] for m in METHODS], width=0.4, edgecolor='white', linewidth=0.5)
    btm_lbl = ['validate', 'db', 'db']
    top_lbl = ['assemble', 'search', 'search']
    for j in range(3):
        if time_bottom[j] > 10:
            ax.text(j, time_bottom[j] / 2, btm_lbl[j], ha='center', va='center',
                    fontsize=5.5, color='white', fontweight='bold', rotation=45)
        if time_top[j] > 5:
            ax.text(j, time_bottom[j] + time_top[j] / 2, top_lbl[j], ha='center', va='center',
                    fontsize=5.5, color='#444', fontweight='bold', rotation=45)
    for j in range(3):
        v = time_bottom[j] + time_top[j]
        m_, s_ = int(v // 60), int(v % 60)
        ax.text(j, v + max(time_bottom + time_top) * 0.02, f'{m_}m {s_:02d}s',
                ha='center', va='bottom', fontsize=7, color='#555')
    ax.set_title('c. Wall time', fontsize=10, fontweight='bold', pad=5)
    ax.set_ylabel('Wall time (s)', fontsize=9, color='#555')
    ax.tick_params(labelsize=8, colors='#555')
    ax.spines[['top', 'right']].set_visible(False)
    ax.spines[['bottom', 'left']].set_color('#CCCCCC')

    # ── (d) Peak memory (human, 原版样式) ──
    ax = fig.add_subplot(gs[0, 3])
    mems = [parse_t('time_cnex_validate.txt').get('mem_mb', 0),
            parse_t('time_mmseqs_search.txt').get('mem_mb', 0),
            parse_t('time_blastn_w15.txt').get('mem_mb', 0)]
    bars = ax.bar(labels, mems, color=[COLORS[m] for m in METHODS], width=0.4,
                  edgecolor='white', linewidth=0.5)
    for b, v in zip(bars, mems):
        ax.text(b.get_x() + b.get_width() / 2, b.get_height() + max(mems) * 0.02,
                f'{v / 1024:.1f} GB', ha='center', va='bottom', fontsize=7, color='#555')
    ax.set_title('d. Peak memory', fontsize=10, fontweight='bold', pad=5)
    ax.set_ylabel('Peak memory (MB)', fontsize=9, color='#555')
    ax.tick_params(labelsize=8, colors='#555')
    ax.spines[['top', 'right']].set_visible(False)
    ax.spines[['bottom', 'left']].set_color('#CCCCCC')

    # ── (e) Recall vs conservation (原版) ──
    thresholds, curve_data = conservation_curve(data)
    ax2 = fig.add_subplot(gs[1, :])
    for mn in METHODS:
        ax2.plot(thresholds * 100, curve_data[mn], label=mn,
                 color=COLORS[mn], linewidth=1.8, linestyle='-')
    ax2.set_xlabel('CNE conservation (%) — consensus identity across species', fontsize=9, color='#555')
    ax2.set_ylabel('Recall (%)', fontsize=9, color='#555')
    ax2.set_title('e. Recall vs cross-species conservation', fontsize=10, fontweight='bold', pad=5)
    ax2.legend(fontsize=8, frameon=False, loc='lower right')
    ax2.set_xlim(48, 102); ax2.set_ylim(0, 105)
    ax2.spines[['top', 'right']].set_visible(False)
    ax2.spines[['bottom', 'left']].set_color('#CCCCCC')
    ax2.tick_params(labelsize=8, colors='#555')
    ax2.text(0.99, 0.02,
             'CNEX: min-c=7, max-diff=2, min-span=25  |  MMseqs2: s=5.7, evalue<1e-5  |  BLASTN: w=15, evalue<1e-5',
             transform=ax2.transAxes, fontsize=7, color='#999', ha='right', va='bottom', style='italic')

    fig.savefig(FIGS / 'benchmark_multi.pdf', bbox_inches='tight', dpi=300)
    plt.close()
    print(f'Saved: figures/benchmark_multi.pdf')

    print('\nSummary (recall / precision):')
    for s in SPECIES:
        line = f'  {s:22s}'
        for m in METHODS:
            n = len(data[s][m]['correct'])
            ok = sum(1 for v in data[s][m]['correct'].values() if v)
            h = sum(1 for v in data[s][m]['hit'].values() if v)
            hok = sum(1 for e, v in data[s][m]['hit'].items() if v and data[s][m]['correct'][e])
            line += f'  {m:7s} {ok / n * 100:5.1f}%/{hok / h * 100:5.1f}%'
        print(line)


if __name__ == '__main__':
    main()
