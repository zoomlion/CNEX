#!/usr/bin/env python3
"""评估 CNEX sweep 输出 vs truth (与 03_analyze 同口径: 任一拷贝重叠 >=30%)
用法: sweep_eval.py <tag>   # 读 results/{T}/cnex_sweeps/{tag}/contigs.fa
"""
import re, sys, random
from pathlib import Path
from collections import defaultdict

BASE = Path(__file__).resolve().parent
DATA = BASE / 'data'
SPECIES = ['Homo_sapiens', 'Mus_musculus', 'Gallus_gallus', 'Scyliorhinus_canicula', 'Callorhinchus_milii']

random.seed(42)
orig = [b.strip() for b in (DATA / 'most-cons-cne.filtered.fa').read_text().strip().split('###') if b.strip()]
sample = sorted(random.sample(range(len(orig)), 10000))


def ng_len(raw_block, sp):
    lines = raw_block.strip().split('\n')
    for i, l in enumerate(lines):
        if l.startswith('>' + sp + '.'):
            seq = ''
            for j in range(i + 1, len(lines)):
                if lines[j].startswith('>') or not lines[j].strip():
                    break
                seq += lines[j].strip()
            return len(seq.replace('-', '').upper())
    return 0


def ovlp(c1, s1, e1, c2, s2, e2, f=0.3):
    return c1 == c2 and max(0, min(e1, e2) - max(s1, s2)) >= f * min(e1 - s1, e2 - s2)


def ovlp_any(copies, c, s, e):
    return any(ovlp(tc, ts, te, c, s, e) for tc, ts, te in copies)


HDR = re.compile(r'>(\d+)\.([^:]+):(\d+)-(\d+)\(([+-])\)')


def eval_species(T, contigs_path):
    truth = defaultdict(list)
    for line in open(DATA / 'truth' / f'{T}.tsv'):
        if not re.match(r'block_\d+', line):
            continue
        p = line.strip().split('\t')
        pos = int(p[0].split('_')[1])
        if pos >= len(sample):
            continue
        if ng_len(orig[sample[pos]], T) < 25:
            continue
        truth[f'block_{pos:05d}'].append((p[2], int(p[3]), int(p[4])))
    cnex = {}
    if Path(contigs_path).exists():
        for line in open(contigs_path):
            if line.startswith('>'):
                m = HDR.match(line.strip())
                if m:
                    cnex[int(m.group(1))] = (m.group(2), int(m.group(3)), int(m.group(4)))
    nt = len(truth)
    correct = tp = fp = 0
    for bid, copies in truth.items():
        ele = int(bid.split('_')[1])
        ok = ele in cnex and ovlp_any(copies, *cnex[ele])
        correct += int(ok)
        if ele in cnex:
            tp += int(ok)
            fp += int(not ok)
    return nt, correct, tp, fp


def main():
    tag = sys.argv[1] if len(sys.argv) > 1 else 'vr2'
    tot_t = tot_c = 0
    print(f'{"species":22s} {"recall":>8s} {"correct":>9s} {"prec":>7s}')
    for T in SPECIES:
        p = BASE / 'results' / T / 'cnex_sweeps' / tag / 'contigs.fa'
        nt, c, tp, fp = eval_species(T, p)
        tot_t += nt
        tot_c += c
        prec = tp / (tp + fp) * 100 if tp + fp else 0.0
        print(f'{T:22s} {c/nt*100:7.1f}% {c:5d}/{nt:<5d} {prec:6.1f}%')
    print(f'{"TOTAL":22s} {tot_c/tot_t*100:7.1f}% {tot_c:5d}/{tot_t:<5d}')


if __name__ == '__main__':
    main()
