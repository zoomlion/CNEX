#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author: JiangminZheng
Date: 2024-12-16
"""

import os, re, sys
import argparse
import glob
from hip import debruijn


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


def main():
    Parser = argparse.ArgumentParser(description='De Bruijn Graph Assembler')
    Parser.add_argument('reads_dir', type=str, help='input reads directory')
    Parser.add_argument('id', type=str, help='element ID to assemble')
    Parser.add_argument('-k', '--kmer', type=int, default=31, help='k-mer length')
    args = Parser.parse_args()
    reads = []
    for file in glob.glob(f"{args.reads_dir}/Assem*.reads"):
        print(file)
        with open(file) as f:
            for line in f:
                line = line.strip()
                eles = line.split('\t')
                if len(eles)!= 4:
                    continue
                read_id, strand, ele_id, read = eles
                strand, ele_id = int(strand), int(ele_id)
                if strand == 1:
                    read = read
                elif strand == -1:
                    read = reverse_complement(read)
                if ele_id == int(args.id):
                    reads.append(read)

    k = args.kmer
    print(reads)
    print(f"Total reads: {len(reads)}")

    graph = debruijn.de_bruijn_graph(reads, k)

    assembled_sequence = debruijn.assemble_sequence(graph)
    print(f"Assembled sequence: {assembled_sequence}")

    debruijn.output_gml(graph, "debruijn_graph.gml")
    print("GML file has been saved as 'debruijn_graph.gml'.")


if __name__ == '__main__':
    main()
