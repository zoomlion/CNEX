#!/usr/bin/env python3
"""按物种提取每个 block 的最佳命中 (evalue < 1e-5), 输出 results/{T}/{tool}/extracted.tsv"""
import re, sys
from pathlib import Path

BASE = Path(__file__).resolve().parent
RES = BASE / 'results'

SPECIES = ['Homo_sapiens', 'Mus_musculus', 'Gallus_gallus', 'Scyliorhinus_canicula', 'Callorhinchus_milii']


def extract(path_in, path_out, label=''):
    best = {}
    for line in open(path_in):
        p = line.strip().split('\t')
        if len(p) < 11: continue
        qid = p[0]; parts = qid.split('|')
        if len(parts) != 3: continue
        bid = parts[0]
        try: evalue = float(p[10])
        except: continue
        if evalue >= 1e-5: continue
        try: ss = int(p[8]); se = int(p[9])
        except: continue
        tgt = p[1]
        s, e = min(ss, se), max(ss, se)
        if bid not in best or evalue < best[bid][3]:
            best[bid] = (tgt, s, e, evalue)

    with open(path_out, 'w') as f:
        f.write('ele_id\tcontig\tstart\tend\tevalue\n')
        for bid, (ctg, s, e, ev) in sorted(best.items()):
            ele = int(bid.split('_')[1])
            f.write(f'{ele}\t{ctg}\t{s}\t{e}\t{ev:.2e}\n')
    print(f'  {path_out}: {len(best)} blocks {label}')


def main():
    targets = sys.argv[1:] or SPECIES
    for T in targets:
        print(f'== {T}')
        extract(RES / T / 'mmseqs' / 'out_s5.7_noid.m8',
                RES / T / 'mmseqs' / 'extracted.tsv', '(MMseqs2)')
        extract(RES / T / 'blastn' / 'out_w15.txt',
                RES / T / 'blastn' / 'extracted.tsv', '(BLASTN w15)')


if __name__ == '__main__':
    main()
