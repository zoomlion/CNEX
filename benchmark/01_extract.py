#!/usr/bin/env python3
"""seed=42 随机采样 10,000 blocks, 按目标物种 T 提取 queries + ground truth + CNEX blocks

用法: 01_extract.py --target Homo_sapiens [--seed 42] [--n-blocks 10000] [--block-frac 0.5]
输出 (block-frac < 1.0 时, 嵌套子采样, 输出到 data/querylib/{T}/bf{frac}/):
  data/queries/{T}.fa    : 排除 T 和 ref 的所有物种序列 (query)
  data/truth/{T}.tsv     : T 在各块中的坐标 (ground truth)
  data/blocks/{T}.fa     : 排除 T 和 ref 的比对块 (CNEX mertable 输入)
"""
import argparse, random, re, os
from pathlib import Path

BASE     = Path(__file__).resolve().parent
INPUT    = BASE / 'data' / 'most-cons-cne.filtered.fa'
OUT_DIR  = BASE / 'data'

HEADER_RE = re.compile(r'>([^\.]+)\.(\S+)([+-]):(\d+)-(\d+)')
SEQ_RE = re.compile(r'^[A-Za-z\-]+$')


def main():
    ap = argparse.ArgumentParser(description='Extract per-target queries/truth/blocks')
    ap.add_argument('--target', required=True, help='target species, e.g. Homo_sapiens')
    ap.add_argument('--seed', type=int, default=42)
    ap.add_argument('--n-blocks', type=int, default=10000)
    ap.add_argument('--block-frac', type=float, default=1.0,
                    help='nested block-level subsample fraction (default 1.0 = all)')
    args = ap.parse_args()

    T = args.target
    random.seed(args.seed)

    with open(INPUT) as f:
        content = f.read()

    raw_blocks = [b.strip() for b in content.strip().split('###') if b.strip()]
    print(f'Reading {INPUT} ...')
    print(f'Total blocks: {len(raw_blocks)}')

    sample_idx = sorted(random.sample(range(len(raw_blocks)), args.n_blocks))

    # 嵌套 block 级子采样: 同 seed 洗牌一次, 取前 k 块 -> 各级子集严格嵌套
    kept = None
    if args.block_frac < 1.0:
        order = list(range(args.n_blocks))
        random.Random(args.seed).shuffle(order)
        k = int(round(args.n_blocks * args.block_frac))
        kept = set(order[:k])
        print(f'Nested subsample: {k}/{args.n_blocks} blocks (frac={args.block_frac:.2f})')

    qdir = OUT_DIR / 'queries'
    tdir = OUT_DIR / 'truth'
    bdir = OUT_DIR / 'blocks'
    for d in (qdir, tdir, bdir):
        os.makedirs(d, exist_ok=True)

    # 输出路径: frac<1 时写到 data/querylib/{T}/bf{frac}/
    if args.block_frac < 1.0:
        tag = f'bf{args.block_frac:.2f}'
        out_root = OUT_DIR / 'querylib' / T / tag
        os.makedirs(out_root, exist_ok=True)
        qpath, bpath, kept_path = out_root / 'queries.fa', out_root / 'blocks.fa', out_root / 'kept_blocks.tsv'
    else:
        qpath, bpath, kept_path = qdir / f'{T}.fa', bdir / f'{T}.fa', None

    queries = []
    truth_lines = []
    block_lines = []

    for pos, bi in enumerate(sample_idx):
        if kept is not None and pos not in kept:
            continue
        block_id = f'block_{pos:05d}'
        block = raw_blocks[bi]
        lines = block.strip().split('\n')

        seqs = []
        cur_hdr = None
        cur_seq = []
        for line in lines:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if cur_hdr is not None:
                    seqs.append((cur_hdr, ''.join(cur_seq)))
                cur_hdr = line
                cur_seq = []
            else:
                if SEQ_RE.match(line) or all(c in 'ACGTacgt-' for c in line):
                    cur_seq.append(line)
        if cur_hdr is not None:
            seqs.append((cur_hdr, ''.join(cur_seq)))

        has_target = False
        block_seqs = []
        for hdr, seq in seqs:
            m = HEADER_RE.match(hdr)
            if not m:
                continue
            species = m.group(1)
            contig  = m.group(2)
            strand  = m.group(3)
            start   = int(m.group(4))
            end     = int(m.group(5))

            if species == T:
                truth_lines.append(f'{block_id}\t{T}\t{contig}\t{start}\t{end}')
                has_target = True

            if species in ('ref', T):
                continue

            q_header = f'{block_id}|{species}|{contig}:{start}-{end}'
            queries.append(f'>{q_header}\n{seq}')

            block_seqs.append(hdr + '\n' + seq)

        if not has_target:
            print(f'  Warning: {block_id} (original block {bi}) has no {T} entry')

        if block_seqs:
            block_lines.append('\n'.join(block_seqs) + '\n###')

        if (pos + 1) % 2000 == 0:
            print(f'  Processed {pos+1}/{args.n_blocks} blocks')

    with open(qpath, 'w') as f:
        f.write('\n'.join(queries) + '\n')

    with open(bpath, 'w') as f:
        f.write('\n'.join(block_lines) + '\n')

    if kept_path is not None:
        with open(kept_path, 'w') as f:
            f.write('block_id\tpos\tsample_idx\n')
            for pos in sorted(kept):
                f.write(f'block_{pos:05d}\t{pos}\t{sample_idx[pos]}\n')

    if args.block_frac >= 1.0:
        with open(tdir / f'{T}.tsv', 'w') as f:
            f.write('block_id\tspecies\tcontig\tstart\tend\n')
            f.write('\n'.join(truth_lines) + '\n')

    print(f'\n{"="*50}')
    print(f'{qpath}:  {len(queries):>8} sequences')
    if args.block_frac < 1.0:
        print(f'{kept_path}:   {len(kept) if kept else 0} blocks')
    print(f'truth:   {len(truth_lines):>8} blocks')
    print(f'{bpath}:   {len(block_lines):>8} blocks')
    print(f'{"="*50}')


if __name__ == '__main__':
    main()
