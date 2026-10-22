#!/usr/bin/env python3
"""Clean a query FASTA: uppercase, remove gaps, and drop empty records."""
import sys
from pathlib import Path

inp = Path(sys.argv[1])
out = Path(sys.argv[2])

out_lines = []
hdr = None
seq = []


def flush():
    if hdr is not None and seq:
        out_lines.append(hdr)
        out_lines.extend(seq)


for line in inp.read_text().splitlines():
    if line.startswith('>'):
        flush()
        hdr = line
        seq = []
    elif hdr is not None and line:
        s = line.replace('-', '').upper()
        if s:
            seq.append(s)
flush()

out.write_text('\n'.join(out_lines) + '\n')
print(f'  clean queries: {sum(1 for x in out_lines if x.startswith(">"))} sequences')
