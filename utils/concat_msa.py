#!/usr/bin/env python3
"""Concatenate trimmed MSA alignments into supermatrix FASTA + partition file."""
import os, sys, glob, argparse
from collections import OrderedDict

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("-i", "--input-dir", required=True)
    p.add_argument("-o", "--output", default="supermatrix.fa")
    p.add_argument("-p", "--partitions", default=None)
    p.add_argument("--suffix", default=".aln")
    p.add_argument("--min-occupancy", type=float, default=0,
                   help="Min species occupancy (0-1) to include an element (default: 0)")
    p.add_argument("--min-site-occupancy", type=float, default=0,
                   help="Min species occupancy per site (0-1) to retain column (default: 0)")
    return p.parse_args()

def average_pairwise_identity(aln_path):
    """Compute average pairwise identity from an aligned FASTA file."""
    seqs = read_fasta(aln_path)
    names = list(seqs.keys())
    if len(names) < 2:
        return 1.0
    total_id = 0
    pairs = 0
    for i in range(len(names)):
        si = seqs[names[i]]
        for j in range(i + 1, len(names)):
            sj = seqs[names[j]]
            match = sum(1 for a, b in zip(si, sj) if a == b and a not in '-Nn?')
            total_col = sum(1 for a, b in zip(si, sj) if a not in '-Nn?' or b not in '-Nn?')
            if total_col > 0:
                total_id += match / total_col
                pairs += 1
    return total_id / pairs if pairs > 0 else 1.0


def clean_seq(s):
    return s.replace("\x00", "N").replace(".", "-")

def trim_alignment_by_occupancy(seqs, min_site_occupancy=0.5):
    """Filter alignment columns by site occupancy.

    Args:
        seqs: {name: sequence_string} — aligned sequences of equal length.
        min_site_occupancy: 0–1, fraction of non-gap species required per column.

    Returns:
        {name: filtered_sequence}
    """
    if not seqs or min_site_occupancy <= 0:
        return seqs
    names = list(seqs.keys())
    ncol = len(seqs[names[0]])
    min_non_gap = max(1, int(len(names) * min_site_occupancy))
    keep = [j for j in range(ncol)
            if sum(1 for sp in names if seqs[sp][j] not in '-Nn?') >= min_non_gap]
    return {sp: ''.join(seqs[sp][j] for j in keep) for sp in names}


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
                header = line[1:].split()[0]
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
        partitions.append(f"GTR+G4, {i}={pos}-{pos + aln_len - 1}")
        pos += aln_len

    # Filter by min-site-occupancy
    if args.min_site_occupancy > 0:
        min_occ = args.min_site_occupancy / 100 if args.min_site_occupancy > 1 else args.min_site_occupancy
        total_sp = len(species_list)
        min_non_gap = max(1, int(total_sp * min_occ))
        flat = {sp: ''.join(supermatrix[sp]) for sp in species_list}
        aln_len = len(next(iter(flat.values())))
        # Build column-to-element map
        elem_of_col = []
        for i, seqs in enumerate(aln_data):
            elem_of_col.extend([i] * len(next(iter(seqs.values()))))
        # Filter columns
        keep_cols = [j for j in range(aln_len)
                     if sum(1 for sp in species_list
                            if flat[sp][j] not in '-Nn?') >= min_non_gap]
        for sp in species_list:
            supermatrix[sp] = [flat[sp][j] for j in keep_cols]
        print(f"Site occupancy filter >= {min_occ:.0%}: {aln_len} -> {len(keep_cols)} sites (removed {aln_len - len(keep_cols)})")

        # Recalculate partitions
        part_counts = {}
        for j in keep_cols:
            ei = elem_of_col[j]
            part_counts[ei] = part_counts.get(ei, 0) + 1
        pos = 1
        partitions = []
        for orig_idx in sorted(part_counts):
            kept = part_counts[orig_idx]
            partitions.append(f"GTR+G4, {orig_idx}={pos}-{pos + kept - 1}")
            pos += kept

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
