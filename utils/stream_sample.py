#!/usr/bin/env python3
"""Bernoulli streaming downsampling of gzipped FASTQ.
Two-pass: count (or --total) -> sample. Memory: O(1)."""
import gzip, sys, random, argparse
sys.stderr.reconfigure(line_buffering=True)

def count_reads(path):
    n = 0
    n_lines = 0
    with gzip.open(path, 'rt') as f:
        for line in f:
            n_lines += 1
            if n_lines % 1000000 == 0:
                print(f"  Counted {n_lines // 1000000}M lines...", file=sys.stderr)
    return n_lines // 4

def sample(path_in, path_out, n_target, n_total):
    random.seed(42)
    remaining = n_total
    target = n_target
    written = 0
    with gzip.open(path_in, 'rt') as fin, gzip.open(path_out, 'wt') as fout:
        for i in range(n_total):
            lines = [fin.readline() for _ in range(4)]
            if not lines[0]:
                break
            if random.random() < target / remaining:
                fout.write(''.join(lines))
                target -= 1
                written += 1
            remaining -= 1
            if (i + 1) % 1000000 == 0:
                print(f"  {(i+1)/1000000:.1f}M / {n_total/1000000:.0f}M, written {written}", file=sys.stderr)
                fout.flush()
    return written

if __name__ == '__main__':
    ap = argparse.ArgumentParser(description='Bernoulli FASTQ downsampling')
    ap.add_argument('in_r1', help='Input R1 FASTQ (gzipped)')
    ap.add_argument('in_r2', help='Input R2 FASTQ (gzipped)')
    ap.add_argument('n_target', type=int, help='Target number of read pairs')
    ap.add_argument('out_r1', help='Output R1 FASTQ (gzipped)')
    ap.add_argument('out_r2', help='Output R2 FASTQ (gzipped)')
    ap.add_argument('--total', type=int, default=None,
                    help='Total read pairs (skip counting, use this value)')
    args = ap.parse_args()

    if args.total:
        n_total = args.total
    else:
        print("Counting reads in R1...", file=sys.stderr)
        n_r1 = count_reads(args.in_r1)
        print("Counting reads in R2...", file=sys.stderr)
        n_r2 = count_reads(args.in_r2)
        n_total = min(n_r1, n_r2)
        print(f"  R1: {n_r1}, R2: {n_r2}, using: {n_total}", file=sys.stderr)

    for inp, outp, lbl in [(args.in_r1, args.out_r1, 'R1'),
                           (args.in_r2, args.out_r2, 'R2')]:
        print(f"Sampling {lbl}: targeting {args.n_target} / {n_total} ...", file=sys.stderr)
        w = sample(inp, outp, args.n_target, n_total)
        print(f"  {lbl}: {w} written", file=sys.stderr)
