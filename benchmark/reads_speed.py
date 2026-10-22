#!/usr/bin/env python3
"""Collect CNEX reads-mode timings into benchmark/stats/reads_speed.tsv.

Data sources (BASE = ~/temp/head2head):
  - stats/speed_cold/cnex_{validate,assemble}_{d}_r*.txt : cold (page cache
    evicted with posix_fadvise(DONTNEED) before each run). PRIMARY.
  - stats/speed/cnex_{validate,assemble}_{d}_r*.txt      : warm (cache pre-read).
    Reference only.
  - stats/megahit_{d}.txt, stats/mmseqs_{d}.txt          : historical cold
    single runs (@64t).

Wall time is I/O- and contention-noisy, and thread counts differ (CNEX 4t vs
megahit/mmseqs 64t), so core-hours (user+sys) is the headline metric. I/O wait
is not counted as CPU, hence core-hours is largely cache-insensitive.

This step only writes the table; figures are made by benchmark/depth_experiment.py.

Usage:
    python3 benchmark/reads_speed.py
Env:
    BASE  working dir with stats/ (default ~/temp/head2head)
"""
import os
import re
import statistics
from pathlib import Path

HERE = Path(__file__).resolve().parent
BASE = Path(os.environ.get('BASE', Path.home() / 'temp' / 'head2head'))
STATS = BASE / 'stats'
COLD = STATS / 'speed_cold'
WARM = STATS / 'speed'
OUT_TSV = HERE / 'stats' / 'reads_speed.tsv'

DEPTHS = ['0.25x', '0.5x', '1x', '2x', '3x', '4x', '5x', '6x', '7x', '8x']
THR_CNEX, THR_REF = 4, 64


def parse_time(path):
    """Parse a GNU /usr/bin/time -v log -> dict(wall_s, core_s, mem_mb)."""
    path = Path(path)
    d = {}
    if not path.exists():
        return d
    for line in open(path):
        if 'Elapsed (wall clock)' in line:
            p = line.strip().split()[-1].split(':')
            if len(p) == 2:
                d['wall_s'] = int(p[0]) * 60 + float(p[1])
            elif len(p) == 3:
                d['wall_s'] = int(p[0]) * 3600 + int(p[1]) * 60 + float(p[2])
        elif 'User time' in line:
            d['user_s'] = float(line.split(':', 1)[1].strip().split()[0])
        elif 'System time' in line:
            d['sys_s'] = float(line.split(':', 1)[1].strip().split()[0])
        elif 'Maximum resident' in line:
            d['mem_mb'] = int(line.split(':', 1)[1].strip().split()[0]) // 1024
    if 'user_s' in d and 'sys_s' in d:
        d['core_s'] = d['user_s'] + d['sys_s']
    return d


def median_logs(paths, key):
    vals = [parse_time(p).get(key) for p in paths]
    vals = [v for v in vals if v is not None]
    return statistics.median(vals) if vals else 0.0


def reps(d, stage, root):
    return sorted(root.glob(f'cnex_{stage}_{d}_r*.txt'),
                  key=lambda p: int(re.search(r'_r(\d+)\.txt$', p.name).group(1)))


rows = []
for d in DEPTHS:
    cold_val = reps(d, 'validate', COLD)
    cold_asm = reps(d, 'assemble', COLD)
    if not cold_val or not cold_asm:
        print(f'  {d}: no cold speed logs, skip')
        continue

    warm_val = reps(d, 'validate', WARM)
    warm_asm = reps(d, 'assemble', WARM)

    cold_core = (median_logs(cold_val, 'core_s') + median_logs(cold_asm, 'core_s')) / 3600
    cold_wall = median_logs(cold_val, 'wall_s') + median_logs(cold_asm, 'wall_s')
    cold_mem = max(parse_time(p).get('mem_mb', 0) for p in cold_val + cold_asm)

    warm_core = ((median_logs(warm_val, 'core_s') + median_logs(warm_asm, 'core_s')) / 3600
                 if warm_val and warm_asm else 0.0)
    warm_wall = (median_logs(warm_val, 'wall_s') + median_logs(warm_asm, 'wall_s')
                 if warm_val and warm_asm else 0.0)

    mh = parse_time(STATS / f'megahit_{d}.txt')
    ms = parse_time(STATS / f'mmseqs_{d}.txt')

    r = {
        'depth': d,
        'cold_val_wall': median_logs(cold_val, 'wall_s'),
        'cold_asm_wall': median_logs(cold_asm, 'wall_s'),
        'cold_wall': cold_wall,
        'cold_core_h': cold_core,
        'cold_mem_gb': cold_mem / 1024.0,
        'warm_wall': warm_wall,
        'warm_core_h': warm_core,
        'megahit_wall': mh.get('wall_s', 0.0),
        'megahit_core_h': mh.get('core_s', 0.0) / 3600,
        'megahit_mem_gb': mh.get('mem_mb', 0) / 1024.0,
        'mmseqs_wall': ms.get('wall_s', 0.0),
        'mmseqs_core_h': ms.get('core_s', 0.0) / 3600,
        'mmseqs_mem_gb': ms.get('mem_mb', 0) / 1024.0,
    }
    r['ratio_asm'] = r['megahit_core_h'] / cold_core if cold_core else 0
    r['ratio_cne'] = ((r['megahit_core_h'] + r['mmseqs_core_h']) / cold_core
                      if cold_core else 0)
    rows.append(r)

if not rows:
    raise SystemExit('no cold speed rows; run CACHE=cold benchmark/run_reads_speed.sh')

OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
cols = ['depth', 'thr_cnex', 'thr_ref', 'cache_primary',
        'cnex_cold_val_wall_s', 'cnex_cold_asm_wall_s', 'cnex_cold_wall_s',
        'cnex_cold_core_h', 'cnex_cold_mem_gb',
        'cnex_warm_wall_s', 'cnex_warm_core_h',
        'megahit_wall_s', 'megahit_core_h', 'megahit_mem_gb',
        'mmseqs_wall_s', 'mmseqs_core_h', 'mmseqs_mem_gb',
        'ratio_asm_coreh', 'ratio_cne_coreh']
with open(OUT_TSV, 'w') as f:
    f.write('\t'.join(cols) + '\n')
    for r in rows:
        f.write('\t'.join(str(v) for v in [
            r['depth'], THR_CNEX, THR_REF, 'cold',
            f"{r['cold_val_wall']:.2f}", f"{r['cold_asm_wall']:.2f}",
            f"{r['cold_wall']:.2f}", f"{r['cold_core_h']:.4f}", f"{r['cold_mem_gb']:.3f}",
            f"{r['warm_wall']:.2f}", f"{r['warm_core_h']:.4f}",
            f"{r['megahit_wall']:.2f}", f"{r['megahit_core_h']:.2f}", f"{r['megahit_mem_gb']:.3f}",
            f"{r['mmseqs_wall']:.2f}", f"{r['mmseqs_core_h']:.4f}", f"{r['mmseqs_mem_gb']:.3f}",
            f"{r['ratio_asm']:.0f}", f"{r['ratio_cne']:.0f}",
        ]) + '\n')
print(f'Saved: {OUT_TSV}')

print(f"\n{'depth':>6} {'CNEX cold':>10} {'CNEX warm':>10} {'MH':>9} {'MH+ms':>9} "
      f"{'ratio asm':>10} {'ratio +ms':>10}")
for r in rows:
    print(f"{r['depth']:>6} {r['cold_core_h']:>10.4f} {r['warm_core_h']:>10.4f} "
          f"{r['megahit_core_h']:>9.2f} {r['megahit_core_h'] + r['mmseqs_core_h']:>9.2f} "
          f"{r['ratio_asm']:>9.0f}x {r['ratio_cne']:>9.0f}x")
