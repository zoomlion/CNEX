#!/usr/bin/env python3
"""Locus-referenced direct-sequence evaluation (final head2head rule).

Rule (agreed 2026-10-14): a claim sequence (a labeled CNEX contig `>ele`,
or a megahit m8 hit region) is TP iff a parasail SW alignment (forward or
reverse-complement) against any truth-locus reference sequence of the same
block has identity >= THR and aligned length >= COV * min(claim, ref).

The reference sequences are the actual genome substrings at the truth
coordinates (extracted once with samtools faidx), so the rule is the
sequence-level equivalent of the genome benchmark's 30%-overlap rule.

Subcommands:
  build-refs --out DIR              extract reference sequences
  cnex --contigs F.fa               evaluate reads-mode CNEX contigs
  regions --regions R.fa            evaluate region fasta (build via
                                    extract_m8_regions.py)

Examples:
  python3 benchmark/eval_locus.py build-refs --out /tmp/opencode/bench/locus
  python3 benchmark/eval_locus.py cnex --contigs cnex_5x.fa \
      --refdir /tmp/opencode/bench/locus --label 5x --tsv results.tsv
"""
import argparse
import os
import random
import re
import subprocess
import sys
from collections import defaultdict

import parasail

BASE = os.path.dirname(os.path.abspath(__file__))
GENOME = os.path.join(BASE, 'data/genomes/Homo_sapiens/Homo_sapiens.fa')
MSA = os.path.join(BASE, 'data/most-cons-cne.filtered.fa')
TRUTH = os.path.join(BASE, 'data/truth/Homo_sapiens.tsv')

MAT = parasail.matrix_create("ACGTN", 2, -1)
COMP = str.maketrans("ACGTN", "TGCAN")


def load_truth():
    """Filtered truth: block_id -> [(contig, start, end), ...] (coords normalized)."""
    random.seed(42)
    orig = [b.strip() for b in open(MSA).read().split('###') if b.strip()]
    sample_idx = sorted(random.sample(range(len(orig)), 10000))

    def ng(raw, sp='Homo_sapiens'):
        lines = raw.strip().split('\n')
        for i, l in enumerate(lines):
            if l.startswith('>') and l.startswith('>' + sp + '.'):
                s = ''
                for j in range(i + 1, len(lines)):
                    if lines[j].startswith('>') or not lines[j].strip():
                        break
                    s += lines[j].strip()
                return len(s.replace('-', '').upper())
        return 0

    truth = defaultdict(list)
    for line in open(TRUTH):
        if not re.match(r'block_\d+', line):
            continue
        p = line.strip().split('\t')
        pos = int(p[0].split('_')[1])
        if pos >= len(sample_idx):
            continue
        if ng(orig[sample_idx[pos]]) < 25:
            continue
        s, e = int(p[3]), int(p[4])
        truth[pos].append((p[2], min(s, e), max(s, e)))
    return truth


def build_refs(outdir):
    truth = load_truth()
    os.makedirs(outdir, exist_ok=True)
    regions, meta = [], []
    for pos in sorted(truth):
        for (ctg, s, e) in truth[pos]:
            regions.append(f'{ctg}:{s}-{e}')
            meta.append(pos)
    r = subprocess.run(['samtools', 'faidx', GENOME] + regions,
                       capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(f'samtools faidx failed: {r.stderr[:300]}')
    out = os.path.join(outdir, 'refs.fa')
    n = 0
    with open(out, 'w') as f:
        cur = None
        buf = []
        i = 0
        for line in r.stdout.splitlines():
            if line.startswith('>'):
                if cur is not None:
                    f.write(f'>{cur}\n{"".join(buf).upper()}\n')
                    n += 1
                cur = meta[i]
                i += 1
                buf = []
            else:
                buf.append(line.strip())
        if cur is not None:
            f.write(f'>{cur}\n{"".join(buf).upper()}\n')
            n += 1
    print(f'build-refs: {n} reference sequences ({len(truth)} blocks) -> {out}')


def load_refs(refdir):
    ref = defaultdict(list)
    cur = None
    buf = []
    for line in open(os.path.join(refdir, 'refs.fa')):
        if line.startswith('>'):
            if cur is not None and buf:
                ref[cur].append(''.join(buf))
            cur = int(line[1:].strip())
            buf = []
        else:
            buf.append(line.strip())
    if cur is not None and buf:
        ref[cur].append(''.join(buf))
    return ref


def match(seq, refs, thr, cov_thr):
    """Return (ok, best_identity)."""
    best_id = 0.0
    for rs in refs:
        for q in (seq, seq.translate(COMP)[::-1]):
            res = parasail.sw_stats_striped_16(q, rs, 1, 1, MAT)
            if res.length == 0:
                continue
            ident = res.matches / res.length
            if ident > best_id:
                best_id = ident
            if ident >= thr and res.length >= cov_thr * min(len(q), len(rs)):
                return True, best_id
    return False, best_id


def parse_cnex(path):
    """-> [(ele, seq)] from >ele headers (reads mode)."""
    out = []
    ele = None
    buf = []
    for line in open(path):
        if line.startswith('>'):
            if ele is not None and buf:
                out.append((ele, ''.join(buf).upper()))
            try:
                ele = int(line[1:].split('.')[0].split()[0])
            except ValueError:
                ele = None
            buf = []
        elif ele is not None:
            buf.append(line.strip())
    if ele is not None and buf:
        out.append((ele, ''.join(buf).upper()))
    return out


def parse_regions(path):
    """-> [(ele, seq)] from >block_XXXXX|... headers."""
    out = []
    ele = None
    buf = []
    for line in open(path):
        if line.startswith('>'):
            if ele is not None and buf:
                out.append((ele, ''.join(buf).upper()))
            m = re.match(r'>block_(\d+)', line)
            ele = int(m.group(1)) if m else None
            buf = []
        elif ele is not None:
            buf.append(line.strip())
    if ele is not None and buf:
        out.append((ele, ''.join(buf).upper()))
    return out


def evaluate(claims, refs, nt, thr, cov_thr):
    tp = fp = 0
    blocks = set()
    for ele, seq in claims:
        if ele not in refs:
            continue
        ok, _ = match(seq, refs[ele], thr, cov_thr)
        if ok:
            tp += 1
            blocks.add(ele)
        else:
            fp += 1
    return tp, fp, len(blocks)


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='cmd', required=True)

    b = sub.add_parser('build-refs')
    b.add_argument('--out', required=True)

    c = sub.add_parser('cnex')
    c.add_argument('--contigs', required=True)
    r = sub.add_parser('regions')
    r.add_argument('--regions', required=True)

    for p in (c, r):
        p.add_argument('--refdir', required=True)
        p.add_argument('--thr', type=float, default=0.95)
        p.add_argument('--cov', type=float, default=0.30)
        p.add_argument('--label', default='-')
        p.add_argument('--tsv', default=None)

    args = ap.parse_args()

    if args.cmd == 'build-refs':
        build_refs(args.out)
        return

    refs = load_refs(args.refdir)
    truth = load_truth()
    nt = len(truth)

    if args.cmd == 'cnex':
        claims = parse_cnex(args.contigs)
        mode = 'cnex'
    else:
        claims = parse_regions(args.regions)
        mode = 'regions'

    in_truth = sum(1 for e, _ in claims if e in refs)
    tp, fp, nblocks = evaluate(claims, refs, nt, args.thr, args.cov)
    recall = nblocks / nt * 100
    prec = tp / (tp + fp) * 100 if tp + fp else 0.0

    print(f'[{args.label}] mode={mode} claims={len(claims)} in_truth={in_truth} '
          f'tp={tp} fp={fp} blocks_tp={nblocks} recall={recall:.1f}% prec={prec:.1f}% '
          f'(thr={args.thr}, cov={args.cov})')

    if args.tsv:
        header = not os.path.exists(args.tsv)
        with open(args.tsv, 'a') as f:
            if header:
                f.write('label\tmode\tclaims\tin_truth\ttp\tfp\tblocks_tp\trecall\tprec\tthr\tcov\n')
            f.write(f'{args.label}\t{mode}\t{len(claims)}\t{in_truth}\t{tp}\t{fp}\t'
                    f'{nblocks}\t{recall:.1f}\t{prec:.1f}\t{args.thr}\t{args.cov}\n')


if __name__ == '__main__':
    main()
