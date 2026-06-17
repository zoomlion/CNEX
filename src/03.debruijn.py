#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os, sys, argparse, random, string
from collections import defaultdict
from hip import debruijn
from hip.validator import MerQueryManager, validate_read


def reverse_complement(seq):
    return seq.translate(str.maketrans("ATCG", "TAGC"))[::-1]


def debruijn_assembler(reads, k, min_count=2):
    graph = debruijn.de_bruijn_graph(reads, k, min_count)
    return debruijn.assemble_sequence(graph)


def reads_generator(f):
    prev_ele_id = None
    temp_reads = {}
    for line in f:
        reads_id, strand, ele_id, seq = line.strip().split("\t")
        seq = seq.upper()
        strand = int(strand)
        if prev_ele_id is None:
            prev_ele_id = ele_id
        if ele_id != prev_ele_id:
            yield prev_ele_id, temp_reads
            temp_reads = {}
            temp_reads[reads_id] = (strand, seq)
            prev_ele_id = ele_id
        else:
            temp_reads[reads_id] = (strand, seq)
    yield prev_ele_id, temp_reads


def assemble_reads(temp_reads, k=35, min_count=2, max_reads=200):
    reads = []
    n = 0
    for reads_id, (strand, seq) in temp_reads.items():
        local_seq = seq if strand == 1 else reverse_complement(seq)
        reads.append(local_seq)
        n += 1
        if n >= max_reads:
            break
    if not reads:
        return ""
    return debruijn_assembler(reads, k, min_count)


def filter_contig(seq, mer_query, mer_size, ele_id, min_c=3):
    if len(seq) < mer_size + min_c:
        return True
    confi_id, _ = validate_read(seq, mer_query, mer_size, min_c)
    return confi_id == ele_id


def main():
    parser = argparse.ArgumentParser(description="De Bruijn Graph Assembly")
    parser.add_argument("inputs_dir", help="directory holding Assemble.*.reads files")
    parser.add_argument("-k", "--kmer", type=int, default=35, help="k-mer length")
    parser.add_argument("-o", "--output", default="assembled.fasta", help="output fasta")
    parser.add_argument("--mers", required=True, help="mers table for assembly validation")
    parser.add_argument("--min_c", type=int, default=3, help="min k-mer matches for validation")
    parser.add_argument("--min_count", type=int, default=2,
        help="min k-mer occurrence to keep in de Bruijn graph (filter errors)")
    parser.add_argument("--max_reads", type=int, default=200,
        help="max reads per element for assembly (prevent blowup)")
    args = parser.parse_args()

    if not os.path.exists(args.inputs_dir):
        raise FileNotFoundError(f"Input directory not found: {args.inputs_dir}")

    print("Loading mers table ...", end=" ", flush=True)
    mer_query = MerQueryManager()
    mer_query.load_from_file(args.mers)
    mer_size = mer_query.get_mer_size()
    print(f"done ({mer_query.size()} mers, k={mer_size})")

    # cat and sort all reads by element ID
    tag = ''.join(random.choices(string.ascii_letters, k=8))
    input_f = f"{tag}.reads"
    os.system(f"cat {args.inputs_dir}/*.reads | sort -k3n > {input_f}")
    if not os.path.exists(input_f):
        raise FileNotFoundError(f"Sorted reads file not found: {input_f}")

    results = []
    with open(input_f) as f:
        for ele_id, temp_reads in reads_generator(f):
            seq = assemble_reads(temp_reads, args.kmer, args.min_count, args.max_reads)
            if not seq:
                continue
            if not filter_contig(seq, mer_query, mer_size, int(ele_id), args.min_c):
                continue
            results.append((ele_id, seq))
            print(f"\rAssembled {len(results)} elements", end="", flush=True)
    print()

    with open(args.output, "w") as f:
        for ele_id, seq in results:
            f.write(f">{ele_id}\n{seq}\n")

    os.remove(input_f)
    print(f"Done. {len(results)} sequences written to {args.output}")


if __name__ == "__main__":
    main()
