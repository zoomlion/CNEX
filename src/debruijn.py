#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author: JiangminZheng
Date: 2024-12-16
"""

import os, re, sys
import glob
import argparse
import tqdm
import subprocess
import shutil
import string
import random
from collections import defaultdict
from hip import debruijn
from tempfile import NamedTemporaryFile


def reverse_complement(seq):
    """
    Return the reverse complement of a DNA sequence.
    with maketrans method
    """
    reverse_map = str.maketrans("ATCG", "TAGC")
    return seq.translate(reverse_map)[::-1]


def debruijn_assembler(reads, k):
    """
    Assemble a longest DNA sequence from a set of reads using a De Bruijn Graph.
    """
    graph = debruijn.de_bruijn_graph(reads, k)
    return debruijn.assemble_sequence(graph)


def reads_generator(input_f: str):
    """
    get a generator for reads from input file
    """
    prev_ele_id = None
    temp_reads = {}
    with open(input_f) as f:
        for line in f.readlines():
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


def merge_loci(locus):
    """
    merge loci by their range: return a dict of affiliation
    e.g. 16:93246851-93247000
         16:93246901-93247050
         16:93246951-93247100
         16:93247001-93247150
         16:93247051-93247200 --> {'16:93246851-93247200': [16:93246851-93247000, .., 16:93247051-93247200]}
    """
    intervals = []
    loci_dict = defaultdict(list)
    for loci in locus:
        chrom, s, e = re.search(r'^(\S+):(\d+)-(\d+)$', loci).groups()
        intervals.append((chrom, int(s), int(e)))
    intervals.sort(key=lambda x: (x[0], x[1]))
    if len(intervals) == 1:
        info = f"{intervals[0][0]}:{intervals[0][1]}-{intervals[0][2]}"
        loci_dict[info] = [info]
    else:
        curr_chrom, curr_start, curr_end = intervals[0]
        current_group = [f"{curr_chrom}:{curr_start}-{curr_end}"]
        
        for i in range(1, len(intervals)):
            chrom, start, end = intervals[i]
            # If different chromosome or no overlap
            if chrom != curr_chrom or start > curr_end:
                # Save current group
                merged = f"{curr_chrom}:{curr_start}-{curr_end}"
                loci_dict[merged] = current_group
                # Start new group
                curr_chrom, curr_start, curr_end = chrom, start, end
                current_group = [f"{chrom}:{start}-{end}"]
            else:
                # Extend current range
                curr_end = max(curr_end, end)
                current_group.append(f"{chrom}:{start}-{end}")
        
        # Add the last group
        merged = f"{curr_chrom}:{curr_start}-{curr_end}"
        loci_dict[merged] = current_group
    
    return loci_dict


def assemble_reads(ele_id: str, temp_reads: dict, mode: str, k=35):
    """
    assemble reads using de bruijn graph from hip
    """
    if mode == 'fq':
        reads = []
        for reads_id, (strand, seq) in temp_reads.items():
            local_seq = seq if strand == 1 else reverse_complement(seq)
            reads.append(local_seq)
        assembled = debruijn_assembler(reads, k)
        return [(ele_id, assembled)]
    elif mode == 'genome':
        locus = []
        reads = {}
        for reads_id, (strand, seq) in temp_reads.items():
            local_seq = seq if strand == 1 else reverse_complement(seq)
            strand = '+' if strand == 1 else '-'
            chrom, s, e = re.search(r'^(\S+):(\d+)-(\d+)$', reads_id).groups()
            locus.append(f"{chrom}{strand}:{s}-{e}")
            reads[f"{chrom}{strand}:{s}-{e}"] = local_seq
        if len(locus) > 100:  # too many loci, skip
            return []
        loci_dict = merge_loci(locus)
        all_assembled = []
        for merged_locus, reads_list in loci_dict.items():
            assembled = debruijn_assembler(
                [reads[loci] for loci in reads_list], k
            )
            assembled_id = f"{ele_id}@{merged_locus}"
            all_assembled.append((assembled_id, assembled))
        return all_assembled


def main():
    parser = argparse.ArgumentParser(description="De Bruijn Graph Assembly")
    parser.add_argument("inputs_dir", type=str, help="input dirs holding files")
    parser.add_argument("-k", "--kmer", type=int, default=35, help="k-mer length")
    parser.add_argument("-l", "--lastz_path", type=str, default="lastz", help="path to lastz")
    parser.add_argument("-o", "--output", type=str, default="assembled.fasta", help="output file")
    parser.add_argument("-t", "--threads", type=int, default=8, help="number of threads")
    args = parser.parse_args()

    # cat and sort all reads to tempfile
    input_f = f"{''.join([random.choice(string.ascii_letters) for _ in range(8)])}.reads"
    os.system(f"cat {args.inputs_dir}/*.reads | sort -k3n > {input_f}")
    if not os.path.exists(input_f):
        raise FileNotFoundError(f"File {input_f} not found.")
    
    # get reads type by first line
    # if starts with '@' is read mode, elif : in first line is fasta mode
    type = None
    with open(input_f) as f:
        first_line = f.readline().strip()
        if first_line.startswith('@'):
            type = 'fq'
        elif ':' in first_line:
            type = 'genome'
        else:
            raise ValueError(f"Unknown file type: {first_line}")

    results = []
    for ele_id, temp_reads in reads_generator(input_f):
        if not ele_id:
            continue
        print(f"\rAssembling {ele_id}", end="", flush=True)
        results.extend(assemble_reads(ele_id, temp_reads, mode=type, k=args.kmer))
    
    # write to fasta file
    with open(args.output, "w") as f:
        for (id, seq) in results:
            f.write(f">{id}\n{seq}\n")

    # clean up input_f
    os.remove(input_f)
            
if __name__ == "__main__":
    main()
