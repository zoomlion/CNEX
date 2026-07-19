#!/usr/bin/env python3
"""Concatenate trimmed MSA alignments into supermatrix FASTA + partition file."""
import os, sys, glob, argparse, re
from collections import OrderedDict


def clean_header(hdr):
    """Replace special characters with underscores for tool compatibility."""
    return re.sub(r"[ '():,;]", "_", hdr)

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("-i", "--input-dir", required=True)
    p.add_argument("-o", "--output", default="supermatrix.fa")
    p.add_argument("-p", "--partitions", default=None)
    p.add_argument("--suffix", default=".aln")
    p.add_argument("--min-occupancy", type=float, default=0,
                   help="Min species occupancy (0-1) to include an element (default: 0)")
    return p.parse_args()

def clean_seq(s):
    return s.replace("\x00", "N").replace(".", "-")

def read_fasta(path):
    seqs = OrderedDict()
    with open(path) as f:
        header, seq = None, []
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header:
                    seqs[header] = "".join(seq)
                header = clean_header(line[1:].split()[0])
                seq = []
            else:
                seq.append(line)
        if header:
            seqs[header] = "".join(seq)
    return seqs

def main():
    args = parse_args()
    files = sorted(glob.glob(os.path.join(args.input_dir, "*" + args.suffix)))
    if not files:
        print(f"No *{args.suffix} files found in {args.input_dir}")
        sys.exit(1)
    print(f"Found {len(files)} alignments")

    # Collect all species names
    all_species = OrderedDict()
    aln_data = []
    for f in files:
        seqs = read_fasta(f)
        for sp in seqs:
            all_species[sp] = None
        aln_data.append(seqs)
        if len(aln_data) % 1000 == 0:
            print(f"  Read {len(aln_data)} alignments, {len(all_species)} species")

    print(f"Total species: {len(all_species)}")
    species_list = list(all_species.keys())

    # Filter by min-occupancy
    if args.min_occupancy > 0:
        min_occ = args.min_occupancy / 100 if args.min_occupancy > 1 else args.min_occupancy
        total = len(species_list)
        filtered = []
        for seqs in aln_data:
            present = sum(1 for s in seqs.values() if any(c.isalpha() for c in s))
            if present / total >= min_occ:
                filtered.append(seqs)
        print(f"Occupancy filter >= {min_occ:.0%}: {len(filtered)}/{len(aln_data)} elements kept")
        aln_data = filtered

    # Concatenate
    supermatrix = {sp: [] for sp in species_list}
    partitions = []
    pos = 1
    for i, seqs in enumerate(aln_data):
        aln_len = len(next(iter(seqs.values()))) if seqs else 0
        if aln_len == 0:
            continue
        for sp in species_list:
            s = seqs.get(sp, "")
            s = clean_seq(s) if s else ""
            supermatrix[sp].append(s if s else "-" * aln_len)
        partitions.append(f"DNA, {i} = {pos}-{pos + aln_len - 1}")
        pos += aln_len

    # Write FASTA
    with open(args.output, "w") as f:
        for sp in species_list:
            f.write(f">{sp}\n{''.join(supermatrix[sp])}\n")
    print(f"Supermatrix: {args.output} ({pos - 1} bp, {len(species_list)} taxa)")

    # Write partitions
    if args.partitions:
        with open(args.partitions, "w") as f:
            f.write("\n".join(partitions) + "\n")
        print(f"Partitions: {args.partitions} ({len(partitions)} partitions)")

if __name__ == "__main__":
    main()
