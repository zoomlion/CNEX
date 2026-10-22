#!/usr/bin/env python3
"""CNEX reads vs Megahit+MMseqs2 depth experiment (main + supplementary).

Main figure (2x2):
  a recall vs coverage
  b recall vs CPU cost (4x crossover)
  c precision vs coverage
  d peak memory at 5x

Supplementary figure 1 (2x2 + legend): length / reference / outcome
  S-a called-CNE length distribution (CNEX vs Megahit regions)
  S-b 5x reference comparison (CNEX reads, CNEX genome, MH+MMseqs)
  S-c outcome composition by truth-length bin (5x)
  S-d 5x outcome composition (reads comparison)

Supplementary figure 2 (2x2 + legend): resource / N50
  S-e core-hours vs coverage
  S-f peak memory vs coverage
  S-g wall time vs coverage
  S-h N50 vs coverage (human, called CNEs) + genome N50

Reads the reads head2head working dir (BASE = ~/temp/head2head):
  cnex/$d/contigs.fa, megahit_new_$d.regions.fa,
  eval/results.tsv, stats/reads_speed.tsv
and the current genome-mode CNEX contigs in the repo:
  benchmark/results/Homo_sapiens/cnex/contigs.fa

Per-block outcomes are read from benchmark/stats/{cnex_<depth>,mh_<depth>}.blocks
(created with benchmark/eval_locus if missing); supplementary S-b-S-d use 5x.

Outputs:
  benchmark/stats/depth_experiment.tsv
  benchmark/stats/recall_by_length.tsv
  benchmark/stats/outcome_by_length.tsv
  benchmark/stats/genome_decomposition.tsv
  benchmark/figures/depth_experiment.pdf
  benchmark/figures/depth_experiment_supp1.pdf
  benchmark/figures/depth_experiment_supp2.pdf

Usage:
    python3 benchmark/depth_experiment.py
Env:
    BASE  working dir (default ~/temp/head2head)
"""
import os
import re
import subprocess
from pathlib import Path

import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

HERE = Path(__file__).resolve().parent
BASE = Path(os.environ.get('BASE', Path.home() / 'temp' / 'head2head'))
GENOME_FA = HERE / 'results' / 'Homo_sapiens' / 'cnex' / 'contigs.fa'
SPEED_TSV = HERE / 'stats' / 'reads_speed.tsv'
OUT_TSV = HERE / 'stats' / 'depth_experiment.tsv'
OUT_PDF = HERE / 'figures' / 'depth_experiment.pdf'
OUT_SUPP1 = HERE / 'figures' / 'depth_experiment_supp1.pdf'
OUT_SUPP2 = HERE / 'figures' / 'depth_experiment_supp2.pdf'
TRUTH_TSV = HERE / 'data' / 'truth' / 'Homo_sapiens.tsv'
REFS = BASE / 'eval' / 'locus' / 'refs.fa'
EVAL_BIN = HERE / 'eval_locus'
BLOCK_DIR = HERE / 'stats'

DEPTHS = ['0.25x', '0.5x', '1x', '2x', '3x', '4x', '5x', '6x', '7x', '8x']
X = [0.25, 0.5, 1, 2, 3, 4, 5, 6, 7, 8]
THR_CNEX, THR_REF = 4, 64

# Called-CNE length bins (truth length, bp).
LEN_EDGES = [0, 60, 90, 120, 150, 200, 300, 10 ** 9]
LEN_LABELS = ['0–59', '60–89', '90–119', '120–149', '150–199', '200–299', '≥300']

# Fallbacks if the genome eval cannot be run; these are the current human
# genome-mode values on the 8425-truth locus rule (eval_locus).
GENOME_RECALL_FALLBACK = 85.6
GENOME_PREC_FALLBACK = 98.0


# ─── data loading ───────────────────────────────────────────────────────────

def assembly_stats(fa):
    """-> (n, N50, avg_len, lengths) over a FASTA of called CNE sequences."""
    fa = Path(fa)
    if not fa.exists():
        return None
    lens, cur = [], 0
    for line in open(fa):
        if line.startswith('>'):
            if cur:
                lens.append(cur)
            cur = 0
        else:
            cur += len(line.strip())
    if cur:
        lens.append(cur)
    if not lens:
        return None
    s = sorted(lens, reverse=True)
    tot, cum, n50 = sum(s), 0, s[-1]
    for l in s:
        cum += l
        if cum >= tot / 2:
            n50 = l
            break
    return len(s), n50, tot / len(s), lens


def fasta_lengths(path, pattern):
    """-> [(element_id, length)] for claims >ele or >block_XXXXX|... ."""
    out = []
    cur = None
    L = 0
    for line in open(path):
        if line.startswith('>'):
            if cur is not None:
                out.append((cur, L))
            m = re.match(pattern, line)
            cur = int(m.group(1)) if m else None
            L = 0
        elif cur is not None:
            L += len(line.strip())
    if cur is not None:
        out.append((cur, L))
    return out


def n50(vals):
    if not vals:
        return 0
    s = sorted(vals, reverse=True)
    tot = sum(s)
    c = 0
    for v in s:
        c += v
        if c >= tot / 2:
            return v
    return 0


def load_eval(tsv):
    """eval/results.tsv -> {'cnex'|'regions': {label: (recall, prec)}}."""
    out = {'cnex': {}, 'regions': {}}
    p = Path(tsv)
    if not p.exists():
        return out
    hdr = None
    for line in open(p):
        if hdr is None:
            hdr = line.rstrip('\n').split('\t')
            continue
        f = dict(zip(hdr, line.rstrip('\n').split('\t')))
        if f.get('mode') in out:
            out[f['mode']][f['label']] = (float(f['recall']), float(f['prec']))
    return out


def load_speed(tsv):
    """reads_speed.tsv -> {depth: dict of core_h / wall_s / mem_gb}."""
    out = {}
    p = Path(tsv)
    if not p.exists():
        return out
    hdr = None
    for line in open(p):
        if hdr is None:
            hdr = line.rstrip('\n').split('\t')
            continue
        f = dict(zip(hdr, line.rstrip('\n').split('\t')))
        try:
            mh_core = float(f['megahit_core_h'])
            ms_core = float(f['mmseqs_core_h'])
            mh_mem = float(f['megahit_mem_gb'])
            ms_mem = float(f['mmseqs_mem_gb'])
            out[f['depth']] = {
                'cnex_core': float(f['cnex_cold_core_h']),
                'mh_core': mh_core,
                'mhms_core': mh_core + ms_core,
                'cnex_wall': float(f['cnex_cold_wall_s']),
                'mh_wall': float(f['megahit_wall_s']),
                'mhms_wall': float(f['megahit_wall_s']) + float(f['mmseqs_wall_s']),
                'cnex_mem': float(f['cnex_cold_mem_gb']),
                'mh_mem': mh_mem,
                'mhms_mem': max(mh_mem, ms_mem),
            }
        except (KeyError, ValueError):
            continue
    return out


def megahit_regions(d):
    for name in (f'megahit_new_{d}.regions.fa', f'megahit_{d}.regions.fa'):
        p = BASE / name
        if p.exists():
            return p
    return None


def load_truth_lengths():
    """block -> max truth interval length (bp)."""
    truth = {}
    if not TRUTH_TSV.exists():
        return truth
    for line in open(TRUTH_TSV):
        p = line.rstrip('\n').split('\t')
        if not re.match(r'block_\d+$', p[0]):
            continue
        b = int(p[0].split('_')[1])
        L = abs(int(p[3]) - int(p[4]))
        truth[b] = max(truth.get(b, 0), L)
    return truth


def load_blocks(path):
    """eval_locus block output -> {block: 0/1}."""
    out = {}
    p = Path(path)
    if not p.exists():
        return out
    for line in open(p):
        f = line.split()
        if len(f) >= 2:
            out[int(f[0])] = int(f[1])
    return out


def run_eval(claims, block_out=None):
    """Run eval_locus; return (recall, precision) and optionally write blocks."""
    if not (EVAL_BIN.exists() and REFS.exists() and Path(claims).exists()):
        return None
    cmd = [str(EVAL_BIN), str(REFS), str(claims), '0.95', '0.30']
    if block_out is not None:
        cmd.append(str(block_out))
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    except Exception:
        return None
    m = re.search(r'recall=([\d.]+)%.*prec=([\d.]+)%', r.stdout)
    if not m:
        return None
    return float(m.group(1)), float(m.group(2))


def ensure_blocks():
    """Create 8x CNEX/MH and genome per-block outcome files if missing."""
    paths = {
        'cnex': BLOCK_DIR / 'cnex_8x.blocks',
        'mh': BLOCK_DIR / 'mh_8x.blocks',
        'genome': BLOCK_DIR / 'genome.blocks',
    }
    claims = {
        'cnex': BASE / 'cnex' / '8x' / 'contigs.fa',
        'mh': megahit_regions('8x'),
        'genome': GENOME_FA,
    }
    for key, out in paths.items():
        if not out.exists() and claims[key] is not None and Path(claims[key]).exists():
            run_eval(claims[key], out)
    return paths


def ensure_all_depth_blocks():
    """Create per-depth block outcomes for the TP-only N50 panel."""
    out = {}
    for d in DEPTHS:
        cpath = BASE / 'cnex' / d / 'contigs.fa'
        mpath = megahit_regions(d)
        cblocks = BLOCK_DIR / f'cnex_{d}.blocks'
        mblocks = BLOCK_DIR / f'mh_{d}.blocks'
        if cpath.exists() and not cblocks.exists():
            run_eval(cpath, cblocks)
        if mpath is not None and mpath.exists() and not mblocks.exists():
            run_eval(mpath, mblocks)
        out[d] = (cblocks, mblocks)
    return out


# ─── length-stratified statistics ────────────────────────────────────────────

def length_bins(length):
    for i, (lo, hi) in enumerate(zip(LEN_EDGES[:-1], LEN_EDGES[1:])):
        if lo <= length < hi:
            return i
    return len(LEN_LABELS) - 1


def bin_stats(cnex_blocks, mh_blocks, truth):
    rows = []
    for i, label in enumerate(LEN_LABELS):
        blocks = [b for b, L in truth.items() if length_bins(L) == i]
        n = len(blocks)
        ctp = sum(1 for b in blocks if cnex_blocks.get(b, 0))
        mtp = sum(1 for b in blocks if mh_blocks.get(b, 0))
        rows.append({
            'bin': label, 'n': n, 'cnex_tp': ctp, 'mh_tp': mtp,
            'cnex_recall': 100.0 * ctp / n if n else 0.0,
            'mh_recall': 100.0 * mtp / n if n else 0.0,
            'gap': 100.0 * (mtp - ctp) / n if n else 0.0,
        })
    return rows


def outcome_stats(cnex_blocks, mh_blocks, truth):
    rows = []
    for i, label in enumerate(LEN_LABELS):
        blocks = [b for b, L in truth.items() if length_bins(L) == i]
        both = sum(1 for b in blocks if cnex_blocks.get(b, 0) and mh_blocks.get(b, 0))
        c_only = sum(1 for b in blocks if cnex_blocks.get(b, 0) and not mh_blocks.get(b, 0))
        m_only = sum(1 for b in blocks if mh_blocks.get(b, 0) and not cnex_blocks.get(b, 0))
        neither = len(blocks) - both - c_only - m_only
        rows.append({
            'bin': label, 'n': len(blocks), 'both': both,
            'cnex_only': c_only, 'mh_only': m_only, 'neither': neither,
        })
    return rows


def genome_decomposition(g, c, m, truth):
    """Genome TP + genome-miss split into MH-finds / reads-only / neither."""
    all_b = set(truth)
    gtp = {b for b, ok in g.items() if ok}
    ctp = {b for b, ok in c.items() if ok}
    mtp = {b for b, ok in m.items() if ok}
    gm = all_b - gtp
    mh_in_gm = gm & mtp
    reads_only_in_gm = (gm & ctp) - mtp
    neither_in_gm = gm - mtp - ctp
    return {
        'genome_tp': len(gtp),
        'genome_miss': len(gm),
        'mh_in_genome_miss': len(mh_in_gm),
        'reads_only_in_genome_miss': len(reads_only_in_gm),
        'neither_in_genome_miss': len(neither_in_gm),
        'mh_only_vs_genome': len(mtp - gtp),
        'genome_only_vs_mh': len(gtp - mtp),
        'neither_genome_nor_mh': len(all_b - gtp - mtp),
    }


def write_tsv(path, cols, rows):
    with open(path, 'w') as f:
        f.write('\t'.join(cols) + '\n')
        for r in rows:
            f.write('\t'.join(str(r.get(c, '')) for c in cols) + '\n')


# ─── style / helpers ────────────────────────────────────────────────────────

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
C_CNEX, C_MH, C_MHMS = '#2B6CB0', '#D97706', '#38A169'
C_GENOME = '#777777'


def _num(v, _pos=None):
    """Readable tick label: 0.01, 1, 100, 1k, 20k (avoids long zero runs)."""
    if v <= 0:
        return ''
    if v >= 1e6:
        return f'{v / 1e6:g}M'
    if v >= 1e3:
        return f'{v / 1e3:g}k'
    return f'{v:g}'


def draw(ax, xs, ys, color, marker, linestyle='-', logy=False):
    ax.plot(xs, ys, linestyle + marker, color=color, linewidth=1.1, markersize=4,
            markerfacecolor='white')
    if logy:
        ax.set_yscale('log')
        from matplotlib.ticker import LogLocator, FuncFormatter
        # readable ticks for narrow ranges (decades alone can fall outside the axis)
        ax.yaxis.set_major_locator(LogLocator(base=10, subs=(1, 2, 5), numticks=10))
        ax.yaxis.set_major_formatter(FuncFormatter(_num))
        ax.yaxis.set_minor_locator(LogLocator(base=10, subs=(), numticks=1))


def depth_label(d):
    return '1/4' if d == '0.25x' else '1/2' if d == '0.5x' else d.rstrip('x')


def draw_violins(ax, dists, positions, color):
    keep = [(d, p) for d, p in zip(dists, positions) if d]
    if not keep:
        return
    vp = ax.violinplot(
        [d for d, _ in keep], positions=[p for _, p in keep], widths=0.36,
        showmeans=False, showmedians=False, showextrema=False, points=100,
    )
    for body in vp['bodies']:
        body.set_facecolor(color)
        body.set_edgecolor(color)
        body.set_alpha(0.35)
        body.set_linewidth(0.5)


def median(vals):
    s = sorted(vals)
    n = len(s)
    if not n:
        return 0
    return s[n // 2] if n % 2 else 0.5 * (s[n // 2 - 1] + s[n // 2])


def mark_median(ax, x, vals, half=0.18):
    if not vals:
        return
    med = median(vals)
    ax.plot([x - half, x + half], [med, med], color='#333333',
            linewidth=1.4, solid_capstyle='butt', zorder=4)


def annotate_folds(a, xs, rows, num, den, ykey, steps=('0.25x', '1x', '4x', '8x')):
    """Label representative depths with MH+MMseqs / CNEX fold factors."""
    for xi, r in zip(xs, rows):
        if r['depth'] in steps and r.get(den, 0) > 0 and r.get(num, 0) > 0:
            a.annotate(f'{r[num] / r[den]:.0f}x', (xi, r[ykey]),
                       textcoords='offset points', xytext=(0, 4), ha='center',
                       fontsize=6.5, color='#888888')


def style_coverage_axis(a):
    a.set_xlabel('Coverage (×)')
    a.set_xscale('log')
    a.set_xticks([0.25, 0.5, 1, 2, 4, 8])
    a.set_xticklabels(['0.25', '0.5', '1', '2', '4', '8'])
    a.set_xlim(0.22, 9)
    a.minorticks_on()
    a.tick_params(which='minor', length=2, color='black', labelsize=0)
    a.tick_params(which='major', length=4, color='black')
    a.spines[['top', 'right']].set_visible(False)
    a.spines[['bottom', 'left']].set_color('#CCCCCC')
    a.grid(True, alpha=0.18, linestyle='--', linewidth=0.5)


# ─── load data ───────────────────────────────────────────────────────────────

genome = assembly_stats(GENOME_FA)
if genome is None:
    raise SystemExit(f'genome baseline not found: {GENOME_FA}')
ev = load_eval(BASE / 'eval' / 'results.tsv')
speed = load_speed(SPEED_TSV)

genome_eval = run_eval(GENOME_FA)
if genome_eval is None:
    genome_recall, genome_prec = GENOME_RECALL_FALLBACK, GENOME_PREC_FALLBACK
    print(f'  warning: genome eval failed; using fallback {genome_recall}/{genome_prec}')
else:
    genome_recall, genome_prec = genome_eval

rows = []
for d in DEPTHS:
    c = assembly_stats(BASE / 'cnex' / d / 'contigs.fa')
    if c is None:
        print(f'  {d}: no CNEX contigs, skip')
        continue
    m = assembly_stats(megahit_regions(d))
    sp = speed.get(d, {})
    cr, cp = ev['cnex'].get(d, (float('nan'), float('nan')))
    mr, mp = ev['regions'].get(d, (float('nan'), float('nan')))
    rows.append({
        'depth': d,
        'cnex_ct': c[0], 'cnex_n50': c[1], 'cnex_avg': c[2], 'cnex_lens': c[3],
        'mh_ct': m[0] if m else 0, 'mh_n50': m[1] if m else 0, 'mh_avg': m[2] if m else 0,
        'mh_lens': m[3] if m else [],
        'cnex_recall': cr, 'cnex_prec': cp, 'mh_recall': mr, 'mh_prec': mp,
        **{k: sp.get(k, 0.0) for k in
           ('cnex_core', 'mh_core', 'mhms_core', 'cnex_wall', 'mh_wall', 'mhms_wall',
            'cnex_mem', 'mh_mem', 'mhms_mem')},
    })

if not rows:
    raise SystemExit('no depth data found; check BASE')

# per-block 8x outcomes (for the length/genome panels)
blocks = ensure_blocks()
truth_lens = load_truth_lengths()
cnex_b = load_blocks(blocks['cnex'])
mh_b = load_blocks(blocks['mh'])
genome_b = load_blocks(blocks['genome'])
has_blocks = bool(cnex_b and mh_b and truth_lens)

if has_blocks:
    lb_rows = bin_stats(cnex_b, mh_b, truth_lens)
    oc_rows = outcome_stats(cnex_b, mh_b, truth_lens)
    gd = None
else:
    lb_rows, oc_rows, gd = [], [], None
    print('  warning: per-block outcome files missing; supplementary length panels skipped')

# TP-only N50: fairer than all-call N50 because MH regions may include flanks.
all_blocks = ensure_all_depth_blocks()
for r in rows:
    d = r['depth']
    cblocks, mblocks = all_blocks.get(d, (None, None))
    ctp = load_blocks(cblocks) if cblocks else {}
    mtp = load_blocks(mblocks) if mblocks else {}
    cpath = BASE / 'cnex' / d / 'contigs.fa'
    mpath = megahit_regions(d)
    cl = fasta_lengths(cpath, r'^>(\d+)') if cpath.exists() else []
    ml = fasta_lengths(mpath, r'^>block_(\d+)') if mpath is not None else []
    r['cnex_n50_tp'] = n50([L for e, L in cl if e in ctp]) or r['cnex_n50']
    r['mh_n50_tp'] = n50([L for e, L in ml if e in mtp]) or r['mh_n50']

# Supplementary length/outcome/reference panels use the 5x level (the
# crossover point where CNEX is closest to MH+MMseqs).
c5 = load_blocks(all_blocks['5x'][0])
m5 = load_blocks(all_blocks['5x'][1])
r5 = next((r for r in rows if r['depth'] == '5x'), None)
if c5 and m5:
    lb_rows = bin_stats(c5, m5, truth_lens)
    oc_rows = outcome_stats(c5, m5, truth_lens)
    gd = genome_decomposition(genome_b, c5, m5, truth_lens)

# ─── write tables ────────────────────────────────────────────────────────────

OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
cols = ['depth', 'cnex_ct', 'cnex_n50', 'cnex_avg', 'mh_ct', 'mh_n50', 'mh_avg',
        'cnex_n50_tp', 'mh_n50_tp', 'cnex_recall', 'cnex_prec', 'mh_recall', 'mh_prec',
        'cnex_core_h', 'mh_core_h', 'mhms_core_h',
        'cnex_wall_s', 'mh_wall_s', 'mhms_wall_s',
        'cnex_mem_gb', 'mh_mem_gb', 'mhms_mem_gb',
        'genome_ct', 'genome_n50', 'genome_avg',
        'genome_recall', 'genome_prec']
with open(OUT_TSV, 'w') as f:
    f.write('\t'.join(cols) + '\n')
    for r in rows:
        f.write('\t'.join(str(v) for v in [
            r['depth'], r['cnex_ct'], r['cnex_n50'], f"{r['cnex_avg']:.0f}",
            r['mh_ct'], r['mh_n50'], f"{r['mh_avg']:.0f}",
            r['cnex_n50_tp'], r['mh_n50_tp'],
            f"{r['cnex_recall']:.1f}", f"{r['cnex_prec']:.1f}",
            f"{r['mh_recall']:.1f}", f"{r['mh_prec']:.1f}",
            f"{r['cnex_core']:.4f}", f"{r['mh_core']:.2f}", f"{r['mhms_core']:.2f}",
            f"{r['cnex_wall']:.1f}", f"{r['mh_wall']:.1f}", f"{r['mhms_wall']:.1f}",
            f"{r['cnex_mem']:.3f}", f"{r['mh_mem']:.3f}", f"{r['mhms_mem']:.3f}",
            genome[0], genome[1], f'{genome[2]:.0f}',
            f'{genome_recall:.1f}', f'{genome_prec:.1f}',
        ]) + '\n')

if lb_rows:
    write_tsv(HERE / 'stats' / 'recall_by_length.tsv',
              ['bin', 'n', 'cnex_tp', 'mh_tp', 'cnex_recall', 'mh_recall', 'gap'],
              [{**r, 'cnex_recall': f"{r['cnex_recall']:.1f}",
                'mh_recall': f"{r['mh_recall']:.1f}",
                'gap': f"{r['gap']:.1f}"} for r in lb_rows])
if oc_rows:
    write_tsv(HERE / 'stats' / 'outcome_by_length.tsv',
              ['bin', 'n', 'both', 'cnex_only', 'mh_only', 'neither'], oc_rows)

print(f'Saved: {OUT_TSV}')
print(f'genome baseline: ct={genome[0]} N50={genome[1]} avg={genome[2]:.0f} '
      f'recall={genome_recall:.1f} prec={genome_prec:.1f}')

# ─── main figure: a = recall vs coverage + cost (shared y); b/c below ────────

xs = [X[DEPTHS.index(r['depth'])] for r in rows]
fig = plt.figure(figsize=(7.2, 7.0))
gs = fig.add_gridspec(2, 2, hspace=0.33, wspace=0.34)
ax_a = fig.add_subplot(gs[0, 0])
ax_b = fig.add_subplot(gs[0, 1])
ax_prec = fig.add_subplot(gs[1, 0])
ax_d = fig.add_subplot(gs[1, 1])

# a. Recall vs coverage
ax_a.plot(xs, [r['cnex_recall'] for r in rows], 'o-', color=C_CNEX,
          linewidth=1.2, markersize=4, markerfacecolor='white')
ax_a.plot(xs, [r['mh_recall'] for r in rows], '^-', color=C_MHMS,
          linewidth=1.2, markersize=4, markerfacecolor='white')
ax_a.set_title('a. Recall')
ax_a.set_ylabel('Recall (%)')
ax_a.set_ylim(0, 100)
style_coverage_axis(ax_a)
ax_a.set_box_aspect(1)

# b. Recall vs cost (normal y-axis)
ax_b.plot([r['cnex_core'] for r in rows], [r['cnex_recall'] for r in rows],
          'o--', color=C_CNEX, linewidth=1.2, markersize=4,
          markerfacecolor='white')
ax_b.plot([r['mhms_core'] for r in rows], [r['mh_recall'] for r in rows],
          '^--', color=C_MHMS, linewidth=1.2, markersize=4,
          markerfacecolor='white')
ax_b.set_xscale('log')
ax_b.set_xlim(0.02, 600)
ax_b.set_xticks([0.01, 0.1, 1, 10, 100])
ax_b.set_xticklabels(['0.01', '0.1', '1', '10', '100'])
ax_b.minorticks_on()
ax_b.set_xlabel('CPU hours')
ax_b.set_ylabel('Recall (%)')
ax_b.set_ylim(0, 100)
ax_b.set_yticks([0, 20, 40, 60, 80, 100])
ax_b.set_title('b. Recall vs cost')
ax_b.spines['top'].set_visible(False)
ax_b.spines['right'].set_visible(False)
ax_b.spines['left'].set_color('#CCCCCC')
ax_b.spines['bottom'].set_color('#CCCCCC')
ax_b.grid(True, alpha=0.18, linestyle='--', linewidth=0.5)
ax_b.set_box_aspect(1)

_r = rows[-1]
r4 = next((r for r in rows if r['depth'] == '4x'), None)
if r4 is not None and r4['cnex_core'] > 0:
    ratio4 = r4['mhms_core'] / r4['cnex_core']
    ax_b.annotate('', xy=(r4['mhms_core'], r4['mh_recall']),
                  xytext=(r4['cnex_core'], r4['cnex_recall']),
                  arrowprops=dict(arrowstyle='<->', color='#555555', lw=0.9,
                                  linestyle='--'))
    ax_b.text((r4['cnex_core'] * r4['mhms_core']) ** 0.5,
              (r4['cnex_recall'] + r4['mh_recall']) / 2 + 3.5,
              f'4× coverage: ~{ratio4:.0f}× less compute', ha='center',
              fontsize=6.5, color='#444444')

handles = [
    Line2D([], [], color=C_CNEX, marker='o', markerfacecolor='white',
           linestyle='-', linewidth=1.2, markersize=4, label='CNEX reads'),
    Line2D([], [], color=C_MHMS, marker='^', markerfacecolor='white',
           linestyle='--', linewidth=1.2, markersize=4, label='Megahit+MMseqs'),
]
fig.legend(handles=handles, loc='lower center', ncol=2, frameon=False,
           fontsize=8, bbox_to_anchor=(0.5, 0.0))

# c. Precision (0-100)
draw(ax_prec, xs, [r['cnex_prec'] for r in rows], C_CNEX, 'o')
draw(ax_prec, xs, [r['mh_prec'] for r in rows], C_MHMS, '^', '--')
ax_prec.set_ylim(0, 100)
ax_prec.set_yticks([0, 20, 40, 60, 80, 100])
ax_prec.set_title('c. Precision')
ax_prec.set_ylabel('Precision (%)')
style_coverage_axis(ax_prec)
ax_prec.set_box_aspect(1)

# d. Peak memory at 5x
if r5 is None:
    ax_d.text(0.5, 0.5, '5x data unavailable', ha='center', va='center',
              transform=ax_d.transAxes, fontsize=8, color='#777777')
    ax_d.set_axis_off()
else:
    mem_x = [0, 1]
    mem_v = [r5['cnex_mem'], r5['mhms_mem']]
    ax_d.bar(mem_x, mem_v, color=[C_CNEX, C_MHMS], width=0.55)
    for x, v in zip(mem_x, mem_v):
        ax_d.text(x, v + 0.8, f'{v:.1f}', ha='center', fontsize=8)
    ax_d.set_xticks(mem_x)
    ax_d.set_xticklabels(['CNEX', 'MH+MMseqs'], fontsize=8)
    ax_d.set_ylabel('Peak memory (GB)')
    ax_d.set_ylim(0, max(32, max(mem_v) * 1.25))
    ax_d.set_title('d. Peak memory (5x)')
    ax_d.annotate('', xy=(1, max(mem_v) * 1.15), xytext=(0, max(mem_v) * 1.15),
                  arrowprops=dict(arrowstyle='<->', color='#555555', lw=0.8,
                                  linestyle='--'))
    ax_d.text(0.5, max(mem_v) * 1.19,
              f'~{mem_v[1] / mem_v[0]:.0f}× less memory',
              ha='center', fontsize=6.5, color='#444444')
    ax_d.spines['top'].set_visible(False)
    ax_d.spines['right'].set_visible(False)
    ax_d.spines['left'].set_color('#CCCCCC')
    ax_d.spines['bottom'].set_color('#CCCCCC')
    ax_d.grid(True, axis='y', alpha=0.18, linestyle='--', linewidth=0.5)
    ax_d.set_box_aspect(1)

fig.tight_layout(pad=0.5)
OUT_PDF.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT_PDF, bbox_inches='tight', dpi=300)
plt.close(fig)
print(f'Saved: {OUT_PDF}')




# ─── supplementary 1: length / reference / outcome (2x2 + legend strip) ─────

fig2 = plt.figure(figsize=(7.2, 7.0))
gs2 = fig2.add_gridspec(2, 2, hspace=0.45, wspace=0.36)
s_len = fig2.add_subplot(gs2[0, 0])
s_ref = fig2.add_subplot(gs2[0, 1])
s_out = fig2.add_subplot(gs2[1, 0])
s_gen = fig2.add_subplot(gs2[1, 1])
for _ax in (s_len, s_ref, s_out, s_gen):
    _ax.set_box_aspect(1)


# S-a called-CNE length distributions
length_pos = list(range(len(rows)))
cnex_pos = [x - 0.20 for x in length_pos]
mh_pos = [x + 0.20 for x in length_pos]
draw_violins(s_len, [r['cnex_lens'] for r in rows], cnex_pos, C_CNEX)
draw_violins(s_len, [r['mh_lens'] for r in rows], mh_pos, C_MH)
all_lens = sorted(v for r in rows for v in r['cnex_lens'] + r['mh_lens'])
y_top = all_lens[int(0.99 * (len(all_lens) - 1))] * 1.22 if all_lens else 1
s_len.set_ylim(0, y_top)
for r, x in zip(rows, cnex_pos):
    mark_median(s_len, x, r['cnex_lens'])
for r, x in zip(rows, mh_pos):
    mark_median(s_len, x, r['mh_lens'])
s_len.set_title('S-a. Called-CNE length')
s_len.set_ylabel('Length (bp)')
s_len.set_xlim(-0.6, len(rows) - 0.4)
s_len.set_xticks([0, 1, 2, 3, 5, 9])
s_len.set_xticklabels(['0.25', '0.5', '1', '2', '4', '8'], fontsize=6.2)
s_len.set_xlabel('Coverage')
s_len.grid(True, axis='y', alpha=0.2, linestyle='--')

# S-b 5x reference comparison
if r5 is None:
    s_ref.text(0.5, 0.5, '5x data unavailable', ha='center', va='center',
               transform=s_ref.transAxes, fontsize=8, color='#777777')
    s_ref.set_axis_off()
else:
    methods = ['CNEX\nreads', 'CNEX\ngenome', 'MH+\nMMseqs']
    rec = [r5['cnex_recall'], genome_recall, r5['mh_recall']]
    pre = [r5['cnex_prec'], genome_prec, r5['mh_prec']]
    cols = [C_CNEX, C_GENOME, C_MHMS]
    xr = [0, 1, 2]
    w = 0.34
    s_ref.bar([x - w / 2 for x in xr], rec, width=w, color=cols, label='Recall')
    s_ref.bar([x + w / 2 for x in xr], pre, width=w, color=cols, alpha=0.55,
              hatch='//', edgecolor='white', label='Precision')
    for x, rv, pv in zip(xr, rec, pre):
        s_ref.text(x - w / 2, rv + 1.5, f'{rv:.1f}', ha='center', fontsize=5.8)
        s_ref.text(x + w / 2, pv + 1.5, f'{pv:.1f}', ha='center', fontsize=5.8)
    s_ref.set_xticks(xr)
    s_ref.set_xticklabels(methods, fontsize=6.2)
    s_ref.set_ylim(0, 110)
    s_ref.set_ylabel('Recall / precision (%)')
    s_ref.set_title('S-b. 5x reference')

# S-c outcome composition by truth length (5x)
if oc_rows:
    bx = list(range(len(oc_rows)))
    cats = [('both', '#B0B0B0', 'both'), ('cnex_only', C_CNEX, 'CNEX-only'),
            ('mh_only', C_MHMS, 'MH-only'), ('neither', '#666666', 'neither')]
    bottom = [0.0] * len(oc_rows)
    for key, color, label in cats:
        vals = [100.0 * r[key] / r['n'] if r['n'] else 0.0 for r in oc_rows]
        s_out.bar(bx, vals, bottom=bottom, color=color, width=0.72, label=label,
                  edgecolor='white', linewidth=0.25)
        bottom = [b + v for b, v in zip(bottom, vals)]
    s_out.set_ylim(0, 100)
    s_out.set_xticks(bx)
    s_out.set_xticklabels([r['bin'] for r in oc_rows], fontsize=6.0, rotation=30,
                          ha='right')
    s_out.set_ylabel('Blocks (%)')
    s_out.set_xlabel('Truth length (bp)')
    s_out.set_title('S-c. Outcome by length (5x)')
    r0 = oc_rows[0]
    if r0['n']:
        s_out.text(0, 100.0 * r0['neither'] / r0['n'] / 2,
                   f'{100.0 * r0["neither"] / r0["n"]:.0f}%',
                   ha='center', va='center', fontsize=5.8, color='white')

# S-d overall 5x outcome (reads comparison; no genome mode)
if c5 and m5 and truth_lens:
    total = len(truth_lens)
    cset = {b for b, ok in c5.items() if ok}
    mset = {b for b, ok in m5.items() if ok}
    both = len(cset & mset)
    c_only = len(cset - mset)
    m_only = len(mset - cset)
    neither = total - len(cset | mset)
    cats = [
        ('Both', both, '#B0B0B0'),
        ('CNEX-only', c_only, C_CNEX),
        ('MH-only', m_only, C_MHMS),
        ('Neither', neither, '#666666'),
    ]
    yy = list(range(len(cats)))
    vals = [v for _, v, _ in cats]
    s_gen.barh(yy, vals, color=[c for _, _, c in cats], height=0.55)
    for y, v in zip(yy, vals):
        s_gen.text(v + total * 0.02, y, f'{v} ({100.0 * v / total:.1f}%)',
                   va='center', fontsize=6.2)
    s_gen.set_yticks(yy)
    s_gen.set_yticklabels([lab for lab, _, _ in cats], fontsize=6.5)
    s_gen.invert_yaxis()
    s_gen.set_xlim(0, total * 1.22)
    s_gen.set_xlabel('Truth blocks')
    s_gen.set_title('S-d. 5x outcome composition')
    s_gen.grid(True, axis='x', alpha=0.18, linestyle='--', linewidth=0.5)
    s_gen.spines['top'].set_visible(False)
    s_gen.spines['right'].set_visible(False)
    s_gen.spines['left'].set_color('#CCCCCC')
    s_gen.spines['bottom'].set_color('#CCCCCC')
    s_gen.text(0.98, 0.02,
               f'MH-only − CNEX-only = +{m_only - c_only} '
               f'({100.0 * (m_only - c_only) / total:.1f}pp)',
               transform=s_gen.transAxes, ha='right', va='bottom',
               fontsize=6.2, color='#444444')

# legend strip for outcome categories
cat_handles = [
    Patch(facecolor='#B0B0B0', label='both'),
    Patch(facecolor=C_CNEX, label='CNEX-only'),
    Patch(facecolor=C_MHMS, label='MH-only'),
    Patch(facecolor='#666666', label='neither'),
]
fig2.legend(handles=cat_handles, frameon=False, loc='lower center',
            ncol=4, fontsize=8, bbox_to_anchor=(0.5, 0.0))

fig2.tight_layout(pad=0.5)
fig2.savefig(OUT_SUPP1, bbox_inches='tight', dpi=300)
plt.close(fig2)
print(f'Saved: {OUT_SUPP1}')

# ─── supplementary 2: resource / N50 (2x2 + legend strip) ───────────────────

fig3 = plt.figure(figsize=(7.2, 7.0))
gs3 = fig3.add_gridspec(2, 2, hspace=0.45, wspace=0.36)
s_cpu = fig3.add_subplot(gs3[0, 0])
s_mem = fig3.add_subplot(gs3[0, 1])
s_wall = fig3.add_subplot(gs3[1, 0])
s_n50 = fig3.add_subplot(gs3[1, 1])
for _ax in (s_cpu, s_mem, s_wall, s_n50):
    _ax.set_box_aspect(1)


# S-e core-hours (vertical speedup lines at every coverage)
draw(s_cpu, xs, [r['cnex_core'] for r in rows], C_CNEX, 'o', logy=True)
draw(s_cpu, xs, [r['mh_core'] for r in rows], C_MH, 's', logy=True)
draw(s_cpu, xs, [r['mhms_core'] for r in rows], C_MHMS, '^', '--', logy=True)
label_steps = {'0.25x', '0.5x', '1x', '2x', '4x', '8x'}
for r in rows:
    x = X[DEPTHS.index(r['depth'])]
    c, m = r['cnex_core'], r['mhms_core']
    if c <= 0 or m <= 0:
        continue
    s_cpu.vlines(x, c, m, color='#999999', linestyle='--', linewidth=0.55,
                 alpha=0.75, zorder=1)
    if r['depth'] in label_steps:
        s_cpu.text(x, (c * m) ** 0.5, f'{m / c:.0f}×', rotation=90,
                   ha='center', va='center', fontsize=5.2, color='#444444',
                   bbox=dict(facecolor='white', edgecolor='none', alpha=0.85,
                             pad=0.4), zorder=4)
s_cpu.set_title('S-e. CPU cost')
s_cpu.set_ylabel('Core-hours (user+sys)')
style_coverage_axis(s_cpu)

# S-f memory
draw(s_mem, xs, [r['cnex_mem'] for r in rows], C_CNEX, 'o', logy=True)
draw(s_mem, xs, [r['mh_mem'] for r in rows], C_MH, 's', logy=True)
draw(s_mem, xs, [r['mhms_mem'] for r in rows], C_MHMS, '^', '--', logy=True)
s_mem.set_title('S-f. Memory (RSS)')
s_mem.set_ylabel('Peak memory (GB)')
style_coverage_axis(s_mem)
annotate_folds(s_mem, xs, rows, 'mhms_mem', 'cnex_mem', 'mhms_mem',
               steps=('8x',))

# S-g wall time
draw(s_wall, xs, [r['cnex_wall'] for r in rows], C_CNEX, 'o', logy=True)
draw(s_wall, xs, [r['mh_wall'] for r in rows], C_MH, 's', logy=True)
draw(s_wall, xs, [r['mhms_wall'] for r in rows], C_MHMS, '^', '--', logy=True)
s_wall.set_title('S-g. Wall time')
s_wall.set_ylabel('Wall clock (s)')
style_coverage_axis(s_wall)
annotate_folds(s_wall, xs, rows, 'mhms_wall', 'cnex_wall', 'mhms_wall',
               steps=('8x',))

# S-h N50
draw(s_n50, xs, [r['cnex_n50_tp'] for r in rows], C_CNEX, 'o')
draw(s_n50, xs, [r['mh_n50_tp'] for r in rows], C_MH, '^', '--')
s_n50.axhline(genome[1], color=C_GENOME, linestyle=':', linewidth=1.0)
s_n50.set_title('S-h. N50 (TP calls)')
s_n50.set_ylabel('N50 (bp)')
s_n50.set_ylim(100, max(genome[1], max(r['mh_n50_tp'] for r in rows)) * 1.15)
style_coverage_axis(s_n50)
s_n50.text(0.03, 0.06, f'genome N50 {genome[1]}',
           transform=s_n50.transAxes, fontsize=5.8, color=C_GENOME)

# legend strip for resource panels
handles_si = [
    Line2D([], [], color=C_CNEX, marker='o', markerfacecolor='white',
           linestyle='-', linewidth=1.2, markersize=4, label='CNEX reads'),
    Line2D([], [], color=C_MH, marker='s', markerfacecolor='white',
           linestyle='-', linewidth=1.2, markersize=4, label='Megahit'),
    Line2D([], [], color=C_MHMS, marker='^', markerfacecolor='white',
           linestyle='--', linewidth=1.2, markersize=4, label='Megahit+MMseqs'),
]
fig3.legend(handles=handles_si, frameon=False, loc='lower center',
            ncol=3, fontsize=8, bbox_to_anchor=(0.5, 0.0))

fig3.tight_layout(pad=0.5)
fig3.savefig(OUT_SUPP2, bbox_inches='tight', dpi=300)
plt.close(fig3)
print(f'Saved: {OUT_SUPP2}')


# ─── console summary ────────────────────────────────────────────────────────

print(f"\n{'depth':>6} | {'CNEX rec/prec':>14} {'N50':>5} {'core-h':>7} {'mem':>6} | "
      f"{'MH reg rec/prec':>15} {'N50':>5} {'core-h':>7} {'mem':>6}")
for r in rows:
    print(f"{r['depth']:>6} | {r['cnex_recall']:>6.1f}/{r['cnex_prec']:<6.1f} "
          f"{r['cnex_n50']:>5} {r['cnex_core']:>7.3f} {r['cnex_mem']:>6.2f} | "
          f"{r['mh_recall']:>7.1f}/{r['mh_prec']:<6.1f} "
          f"{r['mh_n50']:>5} {r['mh_core']:>7.1f} {r['mh_mem']:>6.2f}")
print(f'genome baseline: ct={genome[0]} N50={genome[1]} avg={genome[2]:.0f} '
      f'recall={genome_recall:.1f} prec={genome_prec:.1f}')
if lb_rows:
    print('\nrecall by truth length (5x):')
    for r in lb_rows:
        print(f"  {r['bin']:>7}: n={r['n']:>4} CNEX={r['cnex_recall']:>5.1f}% "
              f"MH={r['mh_recall']:>5.1f}% gap={r['gap']:>4.1f}pp")
