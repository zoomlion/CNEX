#!/usr/bin/env python3
"""CNEX query-lib 子采样评估 (block 级嵌套, 只评估 CNEX)

对每个物种 T 和每个 frac (0.10/0.20/0.40/0.60/0.80 + 1.0 基线复用现有结果):
  - 相对 recall   = kept 子集上找对的块 / kept 子集 truth (bootstrap 95% CI)
  - 绝对 recall   = 找对的块 / 全部 truth (天花板口径)
  - ceiling       = kept truth / 全部 truth
  - precision     = kept 子集上 tp/(tp+fp)
  - recovery      = 基线(1.0)正确的块, 在该 frac 下仍正确的比例
输出: stats/querylib_summary.tsv
"""
import random, re, sys
from pathlib import Path

BASE  = Path(__file__).resolve().parent
DATA  = BASE / 'data'
RES   = BASE / 'results'
STATS = BASE / 'stats'

SPECIES = ['Homo_sapiens', 'Mus_musculus', 'Gallus_gallus', 'Scyliorhinus_canicula', 'Callorhinchus_milii']
FRACS   = [0.10, 0.20, 0.40, 0.60, 0.80]

ORIG_MSA = DATA / 'most-cons-cne.filtered.fa'

random.seed(42)
orig_blocks = [b.strip() for b in ORIG_MSA.read_text().strip().split('###') if b.strip()]
sample_idx = sorted(random.sample(range(len(orig_blocks)), 10000))

HEADER_RE = re.compile(r'>(\d+)\.([^:]+):(\d+)-(\d+)\(([+-])\)')


def species_ng_len(raw_block, species):
    lines = raw_block.strip().split('\n')
    for i, l in enumerate(lines):
        if l.startswith('>') and l.startswith('>' + species + '.'):
            seq = ''
            for j in range(i + 1, len(lines)):
                if lines[j].startswith('>') or not lines[j].strip():
                    break
                seq += lines[j].strip()
            return len(seq.replace('-', '').upper())
    return 0


def load_truth(T):
    """全部 10k 块的 truth (非-gap >= 25bp 过滤, 与 03_analyze 一致; 保留全部拷贝)"""
    truth = {}
    with open(DATA / 'truth' / f'{T}.tsv') as f:
        for line in f:
            if not re.match(r'block_\d+', line):
                continue
            p = line.strip().split('\t')
            bid, _, ctg, s, e = p[:5]
            pos = int(bid.split('_')[1])
            if pos >= len(sample_idx):
                continue
            if species_ng_len(orig_blocks[sample_idx[pos]], T) < 25:
                continue
            truth.setdefault(bid, []).append((ctg, int(s), int(e)))
    return truth


def load_kept(T, frac):
    """嵌套子采样的 kept pos 集合 (frac=1.0 -> 全部)"""
    return set(load_kept_sorted(T, frac))


def load_kept_sorted(T, frac):
    """kept pos 按升序 (== blocks.fa 文件顺序, 供 mertable bunch_id -> pos 重映射)"""
    if frac >= 1.0:
        return list(range(len(sample_idx)))
    p = DATA / 'querylib' / T / f'bf{frac:.2f}' / 'kept_blocks.tsv'
    order = []
    if p.exists():
        for line in open(p):
            if line.startswith('block_id'):
                continue
            order.append(int(line.strip().split('\t')[1]))
    return order


def load_cnex_contigs(T, frac):
    """contigs.fa -> {pos: (chr, start, end)}; 元素 id 在子采样下是文件序,
    经 kept_blocks.tsv 重映射回 10k pos, 与 truth 的 block_id 对齐"""
    cnex = {}
    if frac >= 1.0:
        path = RES / T / 'cnex' / 'contigs.fa'
        order = None
    else:
        path = RES / 'querylib' / T / f'bf{frac:.2f}' / 'cnex' / 'contigs.fa'
        order = load_kept_sorted(T, frac)
    if not path.exists():
        return cnex
    for line in open(path):
        if line.startswith('>'):
            m = HEADER_RE.match(line.strip())
            if m:
                ele = int(m.group(1))
                if order is not None and ele < len(order):
                    ele = order[ele]
                cnex[ele] = (m.group(2), int(m.group(3)), int(m.group(4)))
    return cnex


def ovlp(c1, s1, e1, c2, s2, e2, min_frac=0.3):
    return c1 == c2 and max(0, min(e1, e2) - max(s1, s2)) >= min_frac * min(e1 - s1, e2 - s2)


def ovlp_any(copies, c, s, e, min_frac=0.3):
    """命中任一 truth 拷贝即可 (同一 block 可能存在多份基因座拷贝)"""
    return any(ovlp(tc, ts, te, c, s, e, min_frac) for tc, ts, te in copies)


def wall_seconds(path):
    if not path.exists():
        return 0.0
    for l in open(path):
        if 'Elapsed' in l:
            r = l.split()[-1].strip()
            pts = r.split(':')
            if len(pts) == 3:
                return int(pts[0]) * 3600 + int(pts[1]) * 60 + float(pts[2])
            elif len(pts) == 2:
                return int(pts[0]) * 60 + float(pts[1])
    return 0.0


def level_wall(T, frac):
    if frac >= 1.0:
        d = STATS / T
    else:
        d = STATS / 'querylib' / T / f'bf{frac:.2f}'
    return (wall_seconds(d / 'time_cnex_mertable.txt')
            + wall_seconds(d / 'time_cnex_validate.txt')
            + wall_seconds(d / 'time_cnex_assemble.txt'))


def bootstrap_ci(correct, total, n_iter=1000, seed=42):
    """相对 recall 的 bootstrap 95% CI (块级重抽样)"""
    if total == 0:
        return 0.0, 0.0
    rng = random.Random(seed)
    fracs = []
    for _ in range(n_iter):
        ok = sum(1 for _ in range(total) if rng.random() < correct / total)
        fracs.append(ok / total * 100)
    fracs.sort()
    return fracs[int(0.025 * n_iter)], fracs[int(0.975 * n_iter) - 1]


def evaluate(T):
    truth = load_truth(T)
    all_n = len(truth)
    rows = []
    baseline_correct = None

    for frac in [1.0] + FRACS:
        kept = load_kept(T, frac)
        kept_truth = {bid: t for bid, t in truth.items() if int(bid.split('_')[1]) in kept}
        cnex = load_cnex_contigs(T, frac)

        nt = len(kept_truth)
        correct = 0
        tp = fp = 0
        for bid, copies in kept_truth.items():
            ele = int(bid.split('_')[1])
            ok = ele in cnex and ovlp_any(copies, *cnex[ele])
            correct += int(ok)
            if ele in cnex:
                tp += ok
                fp += not ok
        rel = correct / nt * 100 if nt else 0.0
        ci_low, ci_high = bootstrap_ci(correct, nt)
        abs_r = correct / all_n * 100 if all_n else 0.0
        ceil = nt / all_n * 100 if all_n else 0.0
        prec = tp / (tp + fp) * 100 if tp + fp else 0.0

        if frac == 1.0:
            baseline_correct = {int(b.split('_')[1]) for b in kept_truth
                                if int(b.split('_')[1]) in cnex
                                and ovlp_any(kept_truth[b], *cnex[int(b.split('_')[1])])}
        overlap = baseline_correct & {int(b.split('_')[1]) for b in kept_truth}
        if overlap:
            rec = sum(1 for ele in overlap
                      if ele in cnex and ovlp_any(kept_truth[f'block_{ele:05d}'], *cnex[ele]))
            rec = rec / len(overlap) * 100
        else:
            rec = 0.0

        rows.append((frac, len(kept), nt, correct, rel, ci_low, ci_high,
                     abs_r, ceil, prec, rec, level_wall(T, frac)))
    return rows


def main():
    targets = sys.argv[1:] or SPECIES
    lines = ['species\tfrac\tkept\tnt_kept\tcorrect\t'
             'recall_rel\tci_low\tci_high\trecall_abs\tceiling\tprecision\trecovery\twall_s']
    for T in targets:
        print(f'\n{"="*70}\n== {T} (CNEX query-lib sweep)\n{"="*70}')
        print(f'{"frac":>6s} {"kept":>6s} {"nt":>6s} {"ok":>6s} {"rel%":>7s} {"CI":>13s} '
              f'{"abs%":>7s} {"ceil%":>7s} {"prec%":>7s} {"rec%":>7s} {"wall":>8s}')
        for r in evaluate(T):
            print(f'{r[0]:6.2f} {r[1]:6d} {r[2]:6d} {r[3]:6d} {r[4]:7.1f} '
                  f'{r[5]:5.1f}-{r[6]:5.1f} {r[7]:7.1f} {r[8]:7.1f} {r[9]:7.1f} {r[10]:7.1f} {r[11]:7.1f}s')
            lines.append('\t'.join(map(str, [T, f'{r[0]:.2f}'] + list(r[1:]))))
    with open(STATS / 'querylib_summary.tsv', 'w') as f:
        f.write('\n'.join(lines) + '\n')
    print(f'\nSaved: stats/querylib_summary.tsv')


if __name__ == '__main__':
    main()
