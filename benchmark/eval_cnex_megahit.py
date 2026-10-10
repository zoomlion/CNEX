#!/usr/bin/env python3
"""Evaluate CNEX and Megahit+MMseqs2 against the same ground truth blocks.

Unified recall comparison: both methods evaluated against 8,154 filtered
ground truth CNE blocks for Homo sapiens.

Usage:
    python3 benchmark/eval_cnex_megahit.py [--depth-workdir /tmp/depth_v2]
"""
import re, sys, random, argparse
from pathlib import Path
from collections import defaultdict

random.seed(42)

BASE = Path(__file__).resolve().parent
DATA = BASE / 'data'
TRUTH_TSV = DATA / 'truth' / 'Homo_sapiens.tsv'
MSA_FILE = DATA / 'most-cons-cne.filtered.fa'

# ─── Parse args ───
parser = argparse.ArgumentParser()
parser.add_argument('--depth-workdir', default='/tmp/depth_v2',
                    help='CNEX depth experiment working directory')
args = parser.parse_args()
WORK = Path(args.depth_workdir)


# ─── 1. Load ground truth blocks (with same filtering as depth_analysis.py) ───
print("Loading ground truth...")
truth = {}
with open(TRUTH_TSV) as f:
    for line in f:
        if not re.match(r'block_\d+', line):
            continue
        p = line.strip().split('\t')
        bid, _, ctg, s, e = p[:5]
        truth[bid] = (ctg, int(s), int(e))

# Filter: sample 10k from MSA, non-gap >= 25bp
orig_blocks = [b.strip() for b in MSA_FILE.read_text().strip().split('###') if b.strip()]
sample_idx = sorted(random.sample(range(len(orig_blocks)), 10000))

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

filtered_truth = {}
for bid, (ctg, s, e) in truth.items():
    pos = int(bid.split('_')[1])
    if pos >= len(sample_idx):
        continue
    if species_ng_len(orig_blocks[sample_idx[pos]], 'Homo_sapiens') < 25:
        continue
    filtered_truth[bid] = (ctg, s, e)
truth = filtered_truth
nt = len(truth)
print(f"Ground truth blocks (filtered): {nt}")


# ─── 2. Build chr_num → accession mapping from ground truth ───
chr_to_acc = {}
for bid, (acc, s, e) in truth.items():
    # NC_000015.10 → 15
    m = re.match(r'NC_(\d+)\.\d+', acc)
    if m:
        chr_to_acc[int(m.group(1))] = acc
print(f"Chromosome mappings: {len(chr_to_acc)}")


# ─── 3. Parse CNEX genome contigs → {ele_id: (chr_num, start, end)} ───
def parse_cnex_header(hdr):
    """Parse >ele_id.chr:start-end(strand) → (ele_id, chr_num, start, end)"""
    h = hdr.strip().lstrip('>')
    m = re.match(r'(\d+)\.(\d+):(\d+)-(\d+)\(([+-])\)', h)
    if m:
        return (int(m.group(1)), int(m.group(2)), int(m.group(3)), int(m.group(4)))
    return None

GENOME_FA = WORK / 'genome' / 'assemble' / 'contigs.fa'
if not GENOME_FA.exists():
    # Try remote path
    GENOME_FA_remote = Path.home() / 'Software' / 'CNEX' / 'benchmark' / 'results' / 'cnex' / 'contigs.fa'
    if GENOME_FA_remote.exists():
        GENOME_FA = GENOME_FA_remote
    else:
        print(f"ERROR: CNEX genome contigs not found at {GENOME_FA}")
        sys.exit(1)

print(f"Loading CNEX genome contigs: {GENOME_FA}")
genome_contigs = {}  # ele_id → (chr_num, start, end)
for line in open(GENOME_FA):
    if line.startswith('>'):
        parsed = parse_cnex_header(line)
        if parsed:
            ele, chr_num, start, end = parsed
            genome_contigs[ele] = (chr_num, start, end)
print(f"CNEX genome contigs (unique ele_ids): {len(genome_contigs)}")


# ─── 4. Overlap matching helper ───
def overlap(c1, s1, e1, c2, s2, e2, mf=0.3):
    """Check if two regions overlap with minimum fraction mf."""
    return c1 == c2 and max(0, min(e1, e2) - max(s1, s2)) >= mf * min(e1 - s1, e2 - s2)


# ─── 5. Evaluate CNEX genome-mode vs ground truth ───
print("\n=== CNEX genome-mode ===")
cnex_genome_hits = set()
for bid, (acc, gs, ge) in truth.items():
    m = re.match(r'NC_(\d+)\.\d+', acc)
    if not m:
        continue
    chr_num = int(m.group(1))
    for ele, (c, s, e) in genome_contigs.items():
        if overlap(chr_num, gs, ge, c, s, e):
            cnex_genome_hits.add(bid)
            break

cnex_genome_recall = len(cnex_genome_hits) / nt * 100
print(f"  TP: {len(cnex_genome_hits)}/{nt}  Recall: {cnex_genome_recall:.1f}%")


# ─── 6. Evaluate CNEX 5x depth vs ground truth ───
print("\n=== CNEX 5x depth ===")
DEPTH_FA = WORK / 'out_5x' / 'assemble' / '5x' / 'contigs.fa'
if not DEPTH_FA.exists():
    print(f"  WARNING: {DEPTH_FA} not found, skipping CNEX 5x evaluation")
    cnex_depth_recall = 0
    cnex_depth_hits = set()
else:
    # Parse depth contigs → set of ele_ids
    depth_ele_ids = set()
    for line in open(DEPTH_FA):
        if line.startswith('>'):
            try:
                ele = int(line.strip().lstrip('>'))
                depth_ele_ids.add(ele)
            except ValueError:
                pass
    print(f"  Depth contigs (unique ele_ids): {len(depth_ele_ids)}")

    # Map depth ele_ids to coordinates via genome contigs
    cnex_depth_hits = set()
    for bid, (acc, gs, ge) in truth.items():
        m = re.match(r'NC_(\d+)\.\d+', acc)
        if not m:
            continue
        chr_num = int(m.group(1))
        for ele in depth_ele_ids:
            if ele in genome_contigs:
                c, s, e = genome_contigs[ele]
                if overlap(chr_num, gs, ge, c, s, e):
                    cnex_depth_hits.add(bid)
                    break

    cnex_depth_recall = len(cnex_depth_hits) / nt * 100
    print(f"  TP: {len(cnex_depth_hits)}/{nt}  Recall: {cnex_depth_recall:.1f}%")


# ─── 7. Evaluate Megahit+MMseqs2 5x vs ground truth ───
print("\n=== Megahit+MMseqs2 5x ===")
M8_FILE = WORK / 'mmseqs_results_5x.m8'
if not M8_FILE.exists():
    M8_FILE = Path.home() / 'temp' / '5x' / 'mmseqs_results_5x.m8'

if not M8_FILE.exists():
    print(f"  WARNING: {M8_FILE} not found, skipping Megahit evaluation")
    megahit_recall = 0
    megahit_hits = set()
else:
    megahit_hits = set()
    for line in open(M8_FILE):
        p = line.strip().split('\t')
        if len(p) < 11:
            continue
        qid = p[0]
        bid = qid.split('|')[0]  # block_XXXXX
        try:
            evalue = float(p[10])
        except (IndexError, ValueError):
            continue
        if evalue >= 1e-5:
            continue
        if bid in truth:
            megahit_hits.add(bid)

    megahit_recall = len(megahit_hits) / nt * 100
    print(f"  TP: {len(megahit_hits)}/{nt}  Recall: {megahit_recall:.1f}%")


# ─── 8. Summary table ───
print(f"\n{'='*60}")
print(f"{'Method':<30s} {'TP':>6s} {'Total':>6s} {'Recall':>8s}")
print(f"{'-'*60}")
print(f"{'CNEX genome-mode':<30s} {len(cnex_genome_hits):>6d} {nt:>6d} {cnex_genome_recall:>7.1f}%")
print(f"{'CNEX 5x depth':<30s} {len(cnex_depth_hits):>6d} {nt:>6d} {cnex_depth_recall:>7.1f}%")
print(f"{'Megahit+MMseqs2 5x':<30s} {len(megahit_hits):>6d} {nt:>6d} {megahit_recall:>7.1f}%")
print(f"{'='*60}")

# Previous (unfair) CNEX recall for reference
print(f"\nPrevious (unfair) CNEX recall: 97.0% (vs genome baseline of 7205)")
print(f"Now corrected CNEX recall:     {cnex_genome_recall:.1f}% (vs ground truth of {nt})")
