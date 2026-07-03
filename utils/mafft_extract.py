#!/usr/bin/env python3
"""Extract unaligned sequences from MSA blocks for MAFFT re-alignment."""
import os, sys, random, re
from pathlib import Path

BASE = Path(__file__).resolve().parent
DATA = BASE / 'data'
ORIG_MSA = DATA / 'most-cons-cne.filtered.fa'

random.seed(42)

def parse_block(raw):
    """Return dict: header → aligned sequence (with gaps)."""
    lines = raw.strip().split('\n')
    seqs = {}
    cur_hdr, cur_seq = None, []
    for l in lines:
        l = l.strip()
        if not l: continue
        if l.startswith('>'):
            if cur_hdr is not None and cur_seq:
                seqs[cur_hdr] = ''.join(cur_seq)
            cur_hdr, cur_seq = l, []
        else:
            cur_seq.append(l.upper())
    if cur_hdr is not None and cur_seq:
        seqs[cur_hdr] = ''.join(cur_seq)
    return seqs

# Split MSA
raw_blocks = [b.strip() for b in open(ORIG_MSA).read().strip().split('###') if b.strip()]
sample_idx = sorted(random.sample(range(len(raw_blocks)), 10000))

# Output directories
unaligned_dir = Path('/tmp/mafft_blocks')
aligned_dir = Path('/tmp/mafft_aligned')
unaligned_dir.mkdir(parents=True, exist_ok=True)
aligned_dir.mkdir(parents=True, exist_ok=True)

cmd = sys.argv[1] if len(sys.argv) > 1 else 'extract'

if cmd == 'extract':
    """Phase 1: extract unaligned (gap-removed) sequences for MAFFT."""
    for pos, orig_idx in enumerate(sample_idx):
        block = raw_blocks[orig_idx]
        seqs = parse_block(block)
        out_path = unaligned_dir / f'block_{pos:05d}.fa'
        with open(out_path, 'w') as f:
            for hdr, seq in seqs.items():
                name = hdr.lstrip('>')
                sp = name.split('.')[0] if '.' in name else name
                if sp == 'ref': continue
                pure_seq = seq.replace('-', '')
                if pure_seq:
                    f.write(f'{hdr}\n{pure_seq}\n')
        
        if (pos + 1) % 1000 == 0:
            print(f'  Extracted {pos+1}/10000 blocks', file=sys.stderr)

elif cmd == 'merge':
    """Phase 3: merge MAFFT-aligned blocks into one file (exclude human & ref)."""
    out_path = DATA / 'blocks_10k_mafft.fa'
    count = 0
    with open(out_path, 'w') as out:
        for pos in range(10000):
            aligned_path = aligned_dir / f'block_{pos:05d}.fa'
            if not aligned_path.exists():
                continue
            block = aligned_path.read_text().strip()
            if not block:
                continue
            
            # Parse MAFFT output and exclude human & ref
            lines = block.split('\n')
            out_lines = []
            cur_hdr, cur_seq = None, []
            for l in lines:
                if l.startswith('>'):
                    if cur_hdr is not None and cur_seq:
                        sp = cur_hdr.split('.')[0] if '.' in cur_hdr else cur_hdr
                        if sp not in ('ref', 'Homo_sapiens'):
                            out_lines.append(f'>{cur_hdr}\n{cur_seq}')
                    cur_hdr, cur_seq = l[1:].strip(), ''
                else:
                    cur_seq += l.strip().upper()
            if cur_hdr is not None and cur_seq:
                sp = cur_hdr.split('.')[0] if '.' in cur_hdr else cur_hdr
                if sp not in ('ref', 'Homo_sapiens'):
                    out_lines.append(f'>{cur_hdr}\n{cur_seq}')
            
            if out_lines:
                out.write('\n'.join(out_lines) + '\n###\n')
                count += 1
            
            if (pos + 1) % 1000 == 0:
                print(f'  Merged {pos+1}/10000 blocks', file=sys.stderr)
    
    print(f'Merged {count} blocks to {out_path}', file=sys.stderr)
