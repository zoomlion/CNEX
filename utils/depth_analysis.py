#!/usr/bin/env python3
"""Analyze depth-vs-assembly experiment results."""
import re, statistics, sys
from pathlib import Path

WORK = Path('/tmp/depth_v2')
BASELINE = list((WORK / 'genome' / 'assemble').rglob('contigs.fa'))[0] if (WORK / 'genome' / 'assemble').exists() else None

def parse_contig_header(hdr):
    """Extract ele_id from contig header.
    Reads mode: >5 → ele=5
    Genome mode: >5.chr:start-end(+) → ele=5, chr, start, end
    """
    h = hdr.strip().lstrip('>')
    m = re.match(r'(\d+)\.(\d+):(\d+)-(\d+)\(([+-])\)', h)
    if m:
        ele = int(m.group(1))
        return (ele, int(m.group(2)), int(m.group(3)), int(m.group(4)))
    # Simple header: just ele_id
    try:
        ele = int(h)
        return (ele, None, None, None)
    except:
        return None

# Parse genome baseline contigs
baseline = {}  # ele_id → (chr, start, end)
if BASELINE and BASELINE.exists():
    for line in open(BASELINE):
        if line.startswith('>'):
            parsed = parse_contig_header(line)
            if parsed:
                ele, c, s, e = parsed
                if c is not None:
                    baseline[ele] = (c, s, e)
                else:
                    baseline[ele] = None  # simple format

def ovlp(c1,s1,e1,c2,s2,e2,mf=0.3):
    return c1==c2 and max(0,min(e1,e2)-max(s1,s2))>=mf*min(e1-s1,e2-s2)

def n50(lengths):
    """Compute N50 from sorted list of lengths."""
    total = sum(lengths)
    half = total / 2
    cum = 0
    for l in sorted(lengths, reverse=True):
        cum += l
        if cum >= half: return l
    return 0

results = []
for depth in ['0.25x','0.5x','1x','2x','3x','4x','5x','6x','7x','8x']:
    # Find contigs.fa within the assemble subdirectory
    asm_dir = WORK / f'out_{depth}' / 'assemble'
    cf = None
    if asm_dir.exists():
        for p in asm_dir.rglob('contigs.fa'):
            cf = p; break
    if cf is None or not cf.exists():
        print(f'  {depth}: no results')
        continue

    contigs = {}
    contig_len = {}
    cur_ele = None
    for line in open(cf):
        if line.startswith('>'):
            parsed = parse_contig_header(line)
            if parsed:
                ele, c, s, e = parsed
                contigs[ele] = (c, s, e)
                cur_ele = ele
                contig_len[ele] = 0
        elif cur_ele is not None:
            contig_len[cur_ele] += len(line.strip())

    # Recall vs genome baseline
    tp = 0
    for ele, (c,s,e) in contigs.items():
        if ele in baseline:
            bc = baseline[ele]
            if bc is None:
                # Baseline has simple format too: just check ele_id match
                tp += 1
            elif c is not None and ovlp(bc[0], bc[1], bc[2], c, s, e):
                tp += 1
            elif c is None and bc is not None:
                # Depth simple, genome has coords: just ele_id match
                tp += 1

    recall = tp / len(baseline) * 100 if baseline else 0
    lengths = list(contig_len.values())
    n50_val = n50(lengths) if lengths else 0
    avg_len = sum(lengths)/len(lengths) if lengths else 0

    results.append((depth, len(contigs), len(baseline), tp, recall, n50_val, avg_len))
    print(f'  {depth}: contigs={len(contigs)}  recall={recall:.1f}%  n50={n50_val}  avg_len={avg_len:.0f}')

# Output TSV
with open(WORK / 'depth_summary.tsv', 'w') as f:
    f.write('depth\tcontigs\tbaseline\tTP\trecall\tN50\tavg_len\n')
    for r in results:
        f.write('\t'.join(str(v) for v in r) + '\n')
print(f'\nSaved: {WORK}/depth_summary.tsv')
