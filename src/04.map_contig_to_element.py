#!/usr/bin/env python3
"""
Step 04: Group assembled contigs by CNE element (bunch_id).

Input:
  results/<species>/assembled.fasta  -- contig headers are bunch_ids
  mers_table.tsv                     -- to know which bunch_ids exist

Output:
  results/elements/<bunch_id>.fasta  -- per-element multi-species FASTA

Usage:
  python3 src/05.map_contig_to_element.py --max-elements 0 --parallel 8
"""

import argparse
import os
import sys
import time
from collections import defaultdict
from multiprocessing import Pool


def parse_args():
    p = argparse.ArgumentParser(description="Group assembled contigs by CNE element")
    p.add_argument("--results-dir", default="results", help="Results directory")
    p.add_argument("--max-elements", type=int, default=100,
                   help="Max number of CNE elements to process (0=all)")
    p.add_argument("--min-species", type=int, default=3,
                   help="Min number of species required per element")
    p.add_argument("-o", "--outdir", default="results/elements",
                   help="Output directory for per-element FASTAs")
    p.add_argument("--parallel", type=int, default=1,
                   help="Number of parallel processes for reading species")
    return p.parse_args()


def read_fasta_dict(path):
    """Read FASTA into {header: sequence} dict."""
    seqs = {}
    header = None
    seq = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header is not None:
                    seqs[header] = "".join(seq)
                header = line[1:]
                seq = []
            else:
                seq.append(line)
    if header is not None:
        seqs[header] = "".join(seq)
    return seqs


def _read_one_species(sp_and_dir):
    sp, results_dir = sp_and_dir
    fa = os.path.join(results_dir, sp, "assembled.fasta")
    if not os.path.isfile(fa):
        return {}
    seqs = read_fasta_dict(fa)
    # Convert to {element_id: {sp: seq}} dict
    # bunch_id format may be "0.chr:start-end(+)" in genome mode;
    # strip locus suffix to group all copies under the same element
    result = {}
    for bunch_id, seq in seqs.items():
        bid = bunch_id.split(".")[0]
        result[bid] = {sp: seq}
    return result


def main():
    args = parse_args()

    # Discover all species with assemblies
    species = []
    for sp in sorted(os.listdir(args.results_dir)):
        fa = os.path.join(args.results_dir, sp, "assembled.fasta")
        if os.path.isfile(fa):
            species.append(sp)
    print(f"Found {len(species)} species with assemblies")

    # Read all species in parallel
    t0 = time.time()
    elements = {}
    if args.parallel > 1:
        with Pool(args.parallel) as p:
            for chunk in p.imap_unordered(_read_one_species,
                                          [(sp, args.results_dir) for sp in species]):
                for bid, sp_dict in chunk.items():
                    if bid not in elements:
                        elements[bid] = {}
                    for sp, seq in sp_dict.items():
                        if sp not in elements[bid]:
                            elements[bid][sp] = seq
    else:
        for sp in species:
            chunk = _read_one_species((sp, args.results_dir))
            for bid, sp_dict in chunk.items():
                if bid not in elements:
                    elements[bid] = {}
                for sp, seq in sp_dict.items():
                    if sp not in elements[bid]:
                        elements[bid][sp] = seq
    t1 = time.time()
    print(f"Read {len(elements)} unique elements in {t1 - t0:.1f}s")

    # Filter by min-species
    filtered = {bid: sps for bid, sps in elements.items()
                if len(sps) >= args.min_species}
    print(f"Elements with >= {args.min_species} species: {len(filtered)}")

    # Sort by number of species (descending) and take top max_elements
    sorted_elements = sorted(filtered.items(), key=lambda x: -len(x[1]))
    if args.max_elements > 0:
        sorted_elements = sorted_elements[:args.max_elements]
    print(f"Writing {len(sorted_elements)} elements")

    # Write per-element FASTA (single-threaded to avoid I/O contention)
    os.makedirs(args.outdir, exist_ok=True)
    t0 = time.time()
    for bunch_id, sp_seqs in sorted_elements:
        out_path = os.path.join(args.outdir, f"{bunch_id}.fasta")
        with open(out_path, "w") as f:
            for sp in sorted(sp_seqs.keys()):
                f.write(f">{sp}\n{sp_seqs[sp]}\n")
    t1 = time.time()
    print(f"Written {len(sorted_elements)} element FASTAs to {args.outdir} in {t1 - t0:.1f}s")


if __name__ == "__main__":
    main()
