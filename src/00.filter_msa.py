#!/usr/bin/env python3
"""
00.filter_msa.py - Filter MSA blocks by species count and average sequence length.

Removes blocks with too few species or too short aligned sequences,
producing a cleaner MSA for downstream confident k-mer extraction.
"""

import argparse
import sys


def block_generator(fasta_file):
    buf = ''
    with open(fasta_file) as f:
        for line in f:
            buf += line
            while '###\n' in buf:
                part, buf = buf.split('###\n', 1)
                if part.strip():
                    yield part.strip()
        if buf.strip():
            yield buf.strip()


def block_stats(block_text):
    """Return (n_species, avg_len, min_len) for a block."""
    seqs = []
    current_seq = []
    n_species = 0
    for line in block_text.split('\n'):
        if not line:
            continue
        if line.startswith('>'):
            if current_seq:
                seq = ''.join(current_seq).replace('-', '').replace(' ', '').upper()
                if seq:
                    seqs.append(len(seq))
                n_species += 1
            current_seq = []
        else:
            current_seq.append(line.strip())
    if current_seq:
        seq = ''.join(current_seq).replace('-', '').replace(' ', '').upper()
        if seq:
            seqs.append(len(seq))
        n_species += 1

    if not seqs:
        return n_species, 0, 0

    return n_species, sum(seqs) / len(seqs), min(seqs)


def main():
    parser = argparse.ArgumentParser(
        description="Filter MSA blocks by species count and sequence quality"
    )
    parser.add_argument("msa", help="Input MSA file (###-delimited blocks)")
    parser.add_argument("-o", "--output", default="most-cons-cne.filtered.fa",
                        help="Output filtered MSA file")
    parser.add_argument("--min-species", type=int, default=8,
                        help="Minimum number of species per block (default: 8)")
    parser.add_argument("--min-avg-len", type=int, default=50,
                        help="Minimum average sequence length (non-gap) per block (default: 50)")
    args = parser.parse_args()

    kept = 0
    total = 0
    with open(args.output, 'w') as out:
        for block in block_generator(args.msa):
            total += 1
            n_species, avg_len, min_len = block_stats(block)
            if n_species >= args.min_species and avg_len >= args.min_avg_len:
                out.write(block + '\n')
                out.write('###\n')
                kept += 1

    print(f"Total blocks:    {total}")
    print(f"Kept:            {kept} ({kept / max(total, 1) * 100:.1f}%)")
    print(f"Removed:         {total - kept}")
    print(f"Output:          {args.output}")


if __name__ == '__main__':
    main()
