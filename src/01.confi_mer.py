#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Author: JiangminZheng
Date: build confident kmer from msa alignment (pseudo or not)
'''

import os
import sys
import re
import argparse
import string
import random
import shutil
from collections import defaultdict
from hip import tablemer
from multiprocessing import Pool


def bunch2fas(bunch: str) -> dict:
    """
    Convert a bunch string to a FASTA dictionary, with sequence IDs as keys and sequences as values.
    Sequences with only '-' are ignored.
    """
    fas_dict = {}
    for line in bunch.split('\n'):
        if not line:
            continue
        if line.startswith('>'):
            seq_id = line.strip('>').strip()
            fas_dict[seq_id] = ''
        else:
            fas_dict[seq_id] += line.strip().upper()

    # Remove sequences with no valid characters (e.g., only '-')
    return {
        seq_id: seq
        for seq_id, seq in fas_dict.items()
        if len(set(seq) - {'-'}) > 0  # Exclude sequences with only '-'
    }

def bunch_generator(fasta_file):
    """
    Build a generator to yield bunched FASTA sequences.
    """
    buffer = ''
    with open(fasta_file, 'r') as file:
        for line in file:
            buffer += line
            while '###\n' in buffer:
                part, buffer = buffer.split('###\n', 1)
                if part.strip():
                    yield part.strip()
        # Yield any remaining part after the loop
        if buffer.strip():
            yield buffer.strip()

def generate_mers(seq, mer_size):
    """
    Generate [index, index+mer_size] range seq except gaps ('-') or spaces for every index.
    """
    bias = {}
    bias_value = 0
    pure_seq = seq.replace('-', '')

    if len(pure_seq) < mer_size * 2:
        return

    for i in range(len(seq) - (mer_size - 1)):
        if seq[i] == '-':
            bias_value += 1
        else:
            bias[i-bias_value] = i

    for i in range(len(pure_seq) - (mer_size - 1)):
        local_seq = pure_seq[i:i + mer_size]
        if ' ' not in local_seq:
            yield bias[i], local_seq


def build_confident_kmers(msa, mer_size, min_c, output):
    """
    Build confident kmers from MSA alignment.
    """
    mertable = tablemer.TableMer(13)
    mertable.set_min_entropy(1.4)

    # Load MSA sequences
    for bunch_id, bunch in enumerate(bunch_generator(msa)):
        print(f"\r{bunch_id+1}", end='', flush=True)
        local_affi = defaultdict(list)
        local_mer_count = defaultdict(int)
        fas_dict = bunch2fas(bunch)
        if len(fas_dict) < min_c+1:
            continue
        for header, fas in fas_dict.items():
            if header.startswith('ref'):
                continue
            for loci, mer in generate_mers(fas, mer_size):
                local_affi[mer].append((bunch_id, loci))
                local_mer_count[mer] += 1
        for mer, locations in local_affi.items():
            if len(set(locations)) > 1:
                continue
            bunch_id, loci = list(set(locations))[0]
            # print(mer, bunch_id, loci, local_mer_count[mer])
            mertable.add(mer, bunch_id, loci, local_mer_count[mer])

    # dump
    mertable.dump(output)

def main():
    parser = argparse.ArgumentParser(description="Build confident kmers from MSA alignment.")
    parser.add_argument("msa", type=str, help="Path to MSA file.")
    parser.add_argument("-k", "--mer_size", type=int, default=13, help="Size of the kmers (default: 13).")
    # parser.add_argument("-t", "--threads", type=int, default=8, help="threads number for confident kmer building")
    parser.add_argument("-c", "--min_c", type=int, default=4, help="minimum count of the kmers (default: 4).")
    parser.add_argument("-o", "--output", type=str, default="mers_table.tsv", help="Path to the output TSV file.")
    args = parser.parse_args()

    # total number is the '###' line count input
    total_num = os.popen(f"grep -c '###' {args.msa}").read().strip()
    if not total_num.isdigit():
        print(f"Invalid input file: {args.msa}")
        sys.exit(1)
    total_num = int(total_num)
    print(f"Total number of bunches: {total_num}")
    # build confident kmers
    build_confident_kmers(args.msa, args.mer_size, args.min_c, args.output)


if __name__ == '__main__':
    main()
