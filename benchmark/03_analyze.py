#!/usr/bin/env python3
"""多物种评估: 对每个 target 计算 CNEX/MMseqs2/BLASTN 的 recall/precision (contig 名精确匹配)"""
import re, sys
from collections import defaultdict
from pathlib import Path

BASE  = Path(__file__).resolve().parent
DATA  = BASE / 'data'
RES   = BASE / 'results'
STATS = BASE / 'stats'
FIGS  = BASE / 'figures'

SPECIES = ['Homo_sapiens', 'Mus_musculus', 'Gallus_gallus', 'Scyliorhinus_canicula', 'Callorhinchus_milii']

ORIG_MSA = DATA / 'most-cons-cne.filtered.fa'

# ─── 与 01_extract 相同的采样(过滤空块,保证 block↔原始块映射一致) ───
import random
random.seed(42)
orig_blocks = [b.strip() for b in ORIG_MSA.read_text().strip().split('###') if b.strip()]
sample_idx = sorted(random.sample(range(len(orig_blocks)), 10000))


def species_ng_len(raw_block, species):
    """返回 species 在该 MSA 块中的非-gap 长度"""
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


def parse_extracted(path):
    hits = {}
    if not path.exists(): return hits
    for line in open(path):
        if line.startswith('ele_id'): continue
        p = line.strip().split('\t')
        if len(p) < 4: continue
        try:
            ele = int(p[0]); c = p[1]; s = int(p[2]); e = int(p[3])
        except:
            continue
        hits[ele] = (c, s, e)
    return hits


def parse_m8(path):
    hits = defaultdict(list)
    if not path.exists(): return hits
    for line in open(path):
        p = line.strip().split('\t')
        if len(p) < 10: continue
        try:
            pid = float(p[2]); ss = int(p[8]); se = int(p[9])
        except:
            continue
        qp = p[0].split('|')
        if len(qp) != 3: continue
        hits[qp[0]].append((p[1], min(ss, se), max(ss, se), pid))
    return hits


def ovlp(c1, s1, e1, c2, s2, e2, min_frac=0.3):
    return c1 == c2 and max(0, min(e1, e2) - max(s1, s2)) >= min_frac * min(e1 - s1, e2 - s2)


def ovlp_any(copies, c, s, e, min_frac=0.3):
    """命中任一 truth 拷贝即可 (同一 block 可能存在多份基因座拷贝)"""
    return any(ovlp(tc, ts, te, c, s, e, min_frac) for tc, ts, te in copies)


def compute_precision(contigs, truth):
    tp = fp = 0
    for ele_id, (c, s, e) in contigs.items():
        bid = f'block_{ele_id:05d}'
        if bid in truth:
            if ovlp_any(truth[bid], c, s, e):
                tp += 1
            else:
                fp += 1
    return tp, fp


def evaluate(T):
    print(f'\n{"="*60}\n== {T}\n{"="*60}')

    # 1. Ground truth (过滤: T 非-gap ≥ 25bp; 保留全部拷贝)
    truth = defaultdict(list)
    with open(DATA / 'truth' / f'{T}.tsv') as f:
        for line in f:
            if not re.match(r'block_\d+', line): continue
            p = line.strip().split('\t')
            bid, _, ctg, s, e = p[:5]
            pos = int(bid.split('_')[1])
            if pos >= len(sample_idx): continue
            if species_ng_len(orig_blocks[sample_idx[pos]], T) < 25: continue
            truth[bid].append((ctg, int(s), int(e)))
    nt = len(truth)
    print(f'Ground truth: {nt} blocks (non-gap >= 25bp)')

    # 2. Hits
    mm = parse_extracted(RES / T / 'mmseqs' / 'extracted.tsv')
    bl = parse_extracted(RES / T / 'blastn' / 'extracted.tsv')
    cnex = {}
    hdr_re = re.compile(r'>(\d+)\.([^:]+):(\d+)-(\d+)\(([+-])\)')
    if (RES / T / 'cnex' / 'contigs.fa').exists():
        for line in open(RES / T / 'cnex' / 'contigs.fa'):
            if line.startswith('>'):
                m = hdr_re.match(line.strip())
                if m:
                    cnex[int(m.group(1))] = (m.group(2), int(m.group(3)), int(m.group(4)))
    print(f'CNEX: {len(cnex)} contigs | MMseqs2: {len(mm)} | BLASTN: {len(bl)}')

    # 3. Recall / Venn
    n_cnex = n_mm = n_bl = 0
    n_cnex_only = n_mm_only = n_bl_only = 0
    n_cnex_mm = n_cnex_bl = n_mm_bl = n_all3 = n_none = 0
    outcomes = []
    for bid, copies in truth.items():
        ele = int(bid.split('_')[1])
        ok_cnex = ele in cnex and ovlp_any(copies, *cnex[ele])
        ok_mm   = ele in mm   and ovlp_any(copies, *mm[ele])
        ok_bl   = ele in bl   and ovlp_any(copies, *bl[ele])
        outcomes.append((ele, ele in cnex, ok_cnex, ele in mm, ok_mm, ele in bl, ok_bl))
        n_cnex += ok_cnex; n_mm += ok_mm; n_bl += ok_bl
        if ok_cnex and ok_mm and ok_bl: n_all3 += 1
        elif ok_cnex and not ok_mm and not ok_bl: n_cnex_only += 1
        elif not ok_cnex and ok_mm and not ok_bl: n_mm_only += 1
        elif not ok_cnex and not ok_mm and ok_bl: n_bl_only += 1
        elif ok_cnex and ok_mm and not ok_bl: n_cnex_mm += 1
        elif ok_cnex and not ok_mm and ok_bl: n_cnex_bl += 1
        elif not ok_cnex and ok_mm and ok_bl: n_mm_bl += 1
        else: n_none += 1

    # 4. Precision
    cnex_tp, cnex_fp = compute_precision(cnex, truth)
    mm_tp = mm_fp = 0
    for ele, (c, s, e) in mm.items():
        bid = f'block_{ele:05d}'
        if bid not in truth: continue
        if ovlp_any(truth[bid], c, s, e): mm_tp += 1
        else: mm_fp += 1
    bl_tp = bl_fp = 0
    for ele, (c, s, e) in bl.items():
        bid = f'block_{ele:05d}'
        if bid not in truth: continue
        if ovlp_any(truth[bid], c, s, e): bl_tp += 1
        else: bl_fp += 1

    print(f'CNEX:    {n_cnex:5d} ({n_cnex/nt*100:.1f}%)  prec={cnex_tp/(cnex_tp+cnex_fp)*100:.1f}% ({cnex_tp}/{cnex_tp+cnex_fp})')
    print(f'MMseqs2: {n_mm:5d} ({n_mm/nt*100:.1f}%)  prec={mm_tp/(mm_tp+mm_fp)*100:.1f}% ({mm_tp}/{mm_tp+mm_fp})')
    print(f'BLASTN:  {n_bl:5d} ({n_bl/nt*100:.1f}%)  prec={bl_tp/(bl_tp+bl_fp)*100:.1f}% ({bl_tp}/{bl_tp+bl_fp})')
    print(f'Venn: all3={n_all3} cnex_mm={n_cnex_mm} cnex_bl={n_cnex_bl} mm_bl={n_mm_bl} '
          f'cnex_only={n_cnex_only} mm_only={n_mm_only} bl_only={n_bl_only} none={n_none}')

    return {
        'species': T, 'nt': nt,
        'outcomes': outcomes,
        'cnex': (n_cnex, cnex_tp, cnex_fp), 'mm': (n_mm, mm_tp, mm_fp), 'bl': (n_bl, bl_tp, bl_fp),
    }


def parse_t(p):
    d = {}
    if p.exists():
        for l in open(p):
            if 'User time' in l: d['user'] = float(l.split(':')[1].strip().split()[0])
            elif 'Elapsed' in l:
                raw = l.split()[-1].strip()
                parts = raw.split(':')
                if len(parts) == 2: d['wall'] = f'0h{int(parts[0]):02d}m{float(parts[1]):04.1f}s'
                elif len(parts) == 3: d['wall'] = f'{int(parts[0]):d}h{int(parts[1]):02d}m{float(parts[2]):04.1f}s'
                else: d['wall'] = raw
            elif 'Maximum resident' in l:
                d['mem'] = int(l.split(':')[1].strip().split()[0]) // 1024
    return d


def timing_table():
    print(f'\n{"-"*70}\nTIMING (wall / mem)\n{"-"*70}')
    for T in SPECIES:
        row = []
        for lbl, fn in [('mertable', 'time_cnex_mertable.txt'), ('validate', 'time_cnex_validate.txt'),
                        ('assemble', 'time_cnex_assemble.txt'), ('mmseqs_sr', 'time_mmseqs_search.txt'),
                        ('blastn', 'time_blastn_w15.txt')]:
            d = parse_t(STATS / T / fn)
            if d:
                row.append(f'{lbl}:{d.get("wall","?")}/{d.get("mem",0)}MB')
            else:
                row.append(f'{lbl}:n/a')
        print(f'  {T:22s} ' + '  '.join(row))


def sensitivity(T='Homo_sapiens'):
    print(f'\n{"="*60}\nMMSEQS2 SENSITIVITY COMPARISON ({T})\n{"="*60}')
    print(f'{"s":>5s}  {"wall":>12s}  {"hits":>8s}  {"recall":>8s}  {"mem":>6s}')
    truth = defaultdict(list)
    with open(DATA / 'truth' / f'{T}.tsv') as f:
        for line in f:
            if not re.match(r'block_\d+', line): continue
            p = line.strip().split('\t')
            bid, _, ctg, s, e = p[:5]
            pos = int(bid.split('_')[1])
            if pos >= len(sample_idx): continue
            if species_ng_len(orig_blocks[sample_idx[pos]], T) < 25: continue
            truth[bid].append((ctg, int(s), int(e)))
    nt = len(truth)
    results = []
    for s in [1.0, 4.0, 5.7, 7.5]:
        if abs(s - 5.7) < 0.01:
            m8_path = RES / T / 'mmseqs' / 'out_s5.7_noid.m8'
            tim_path = STATS / T / 'time_mmseqs_search.txt'
        else:
            m8_path = RES / T / 'mmseqs' / f'out_s{s}.m8'
            tim_path = STATS / T / f'time_mmseqs_s{s}.txt'
        if not m8_path.exists() or not tim_path.exists():
            print(f'{s:5.1f}  {"N/A":>12s}  {"N/A":>8s}  {"N/A":>8s}  {"N/A":>6s}')
            continue
        sens_hits = parse_m8(m8_path)
        sens_correct = 0
        for bid, copies in truth.items():
            if bid in sens_hits:
                for c, hs, he, pid in sens_hits[bid]:
                    if pid >= 0.7 and ovlp_any(copies, c, hs, he):
                        sens_correct += 1
                        break
        td = parse_t(tim_path)
        wall_s = td.get('wall', '?'); mem = td.get('mem', 0)
        recall = sens_correct / nt * 100
        print(f'{s:5.1f}  {wall_s:>12s}  {len(sens_hits):>8d}  {recall:>7.1f}%  {mem:>5}MB')
        results.append((s, wall_s, len(sens_hits), recall, mem))
    return results


def main():
    targets = sys.argv[1:] or SPECIES
    rows = [evaluate(T) for T in targets]

    with open(STATS / 'outcomes.tsv', 'w') as f:
        f.write('species\tele_id\tcnex_hit\tcnex_correct\tmm_hit\tmm_correct\tbl_hit\tbl_correct\n')
        for r in rows:
            for (ele, c_hit, c_ok, m_hit, m_ok, b_hit, b_ok) in r['outcomes']:
                f.write(f'{r["species"]}\t{ele}\t{int(c_hit)}\t{int(c_ok)}\t{int(m_hit)}\t{int(m_ok)}\t{int(b_hit)}\t{int(b_ok)}\n')
    with open(STATS / 'summary.tsv', 'w') as f:
        f.write('species\tmethod\tcorrect\ttotal\trecall\tprecision\ttp\tfp\n')
        for r in rows:
            for name, (n, tp, fp) in [('CNEX', r['cnex']), ('MMseqs2', r['mm']), ('BLASTN', r['bl'])]:
                prec = f'{tp/(tp+fp)*100:.1f}' if tp + fp else 'N/A'
                f.write(f'{r["species"]}\t{name}\t{n}\t{r["nt"]}\t{n/r["nt"]*100:.1f}\t{prec}\t{tp}\t{fp}\n')
        if 'Homo_sapiens' in targets:
            f.write('\nsensitivity_species\tHomo_sapiens\n')
            f.write('s\thits\trecall\n')
            for s, wall, hits, rec, mem in sensitivity():
                f.write(f'{s:.1f}\t{hits}\t{rec:.1f}\n')
    print(f'\nSaved: stats/summary.tsv')
    timing_table()


if __name__ == '__main__':
    main()
