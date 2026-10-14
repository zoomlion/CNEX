#!/usr/bin/env python3
"""Extract the best-evalue m8 hit region per block from an assembly fasta.

For each block_XXXXX query, keep the lowest-evalue hit (evalue < 1e-5, same
rule as 02_m8extract.py), then slice the aligned interval [tstart, tend]
(1-based inclusive) out of the hit contig and write it as

    >block_XXXXX|contig|start-end
    SEQUENCE

The output is consumed by eval_locus.py regions.

Usage:
  python3 benchmark/extract_m8_regions.py \
      --m8 mmseqs_results_5x.m8 --contigs final.contigs.fa \
      --out megahit_5x.regions.fa
"""
import argparse
import re
from collections import defaultdict


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--m8', required=True)
    ap.add_argument('--contigs', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--max-evalue', type=float, default=1e-5)
    args = ap.parse_args()

    best = {}
    for line in open(args.m8):
        p = line.rstrip('\n').split('\t')
        if len(p) < 11:
            continue
        qp = p[0].split('|')
        if len(qp) != 3 or not re.match(r'block_\d+$', qp[0]):
            continue
        try:
            ev = float(p[10])
            t0, t1 = int(p[8]), int(p[9])
        except ValueError:
            continue
        if ev >= args.max_evalue:
            continue
        bid = qp[0]
        if bid not in best or ev < best[bid][3]:
            best[bid] = (p[1], t0, t1, ev)

    want = defaultdict(list)
    for bid, (ctg, t0, t1, ev) in best.items():
        s, e = min(t0, t1), max(t0, t1)
        want[ctg].append((bid, s, e))
    print(f'blocks with hits={len(best)}; contigs involved={len(want)}')

    wrote = 0
    name = None
    buf = []

    def flush():
        nonlocal wrote
        if name is None or name not in want:
            return
        seq = ''.join(buf).upper()
        for bid, s, e in want[name]:
            sub = seq[s - 1:e]
            if sub:
                out.write(f'>{bid}|{name}|{s}-{e}\n{sub}\n')
                wrote += 1
        del want[name]

    with open(args.out, 'w') as out:
        for line in open(args.contigs):
            if line.startswith('>'):
                flush()
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
        flush()
    print(f'regions written={wrote} -> {args.out}')


if __name__ == '__main__':
    main()
