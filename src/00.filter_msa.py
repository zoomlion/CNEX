#!/usr/bin/env python3
"""
00.filter_msa.py - Filter MSA blocks by species count and average sequence length.

Removes blocks with too few species or too short aligned sequences,
producing a cleaner MSA for downstream confident k-mer extraction.
"""

import argparse
import sys


_STRIP_TRANS = str.maketrans("", "", "- ")


def block_generator(fasta_file):
    lines = []
    with open(fasta_file) as f:
        for line in f:
            if line.endswith("###\n"):
                lines.append(line[:-4])
                part = "".join(lines).strip()
                lines.clear()
                if part:
                    yield part
            else:
                lines.append(line)
        if lines:
            part = "".join(lines).strip()
            if part:
                yield part


def block_stats(block_text):
    """Return (n_species, avg_len, min_len) for a block."""
    seqs = []
    current_len = 0
    in_seq = False
    n_species = 0
    for line in block_text.split('\n'):
        if not line:
            continue
        if line[0] == '>':
            if in_seq:
                if current_len > 0:
                    seqs.append(current_len)
                n_species += 1
            current_len = 0
            in_seq = False
        else:
            in_seq = True
            current_len += len(line.strip().translate(_STRIP_TRANS))
    if in_seq:
        if current_len > 0:
            seqs.append(current_len)
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
