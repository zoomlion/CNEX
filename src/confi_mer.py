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
import math
from collections import defaultdict
import tqdm
from multiprocessing import Pool, cpu_count

def bunch2fas(bunch: str) -> dict:
    """
    Convert bunch string to fas_dict, with key as seq_id and value as seq_str.
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

    # Remove lines with no ATCG values
    return {
        id: seq for id, seq in fas_dict.items() if len(set(seq)) > 1
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

def calculate_entropy(sequence):
    """
    Calculates the Shannon entropy of a DNA sequence.

    Parameters:
        sequence (str): DNA sequence.

    Returns:
        float: Shannon entropy of the sequence.
    """
    from math import log2

    sequence = sequence.upper()
    total = len(sequence)
    if total == 0:
        return 0

    # Calculate probabilities and entropy
    counts = {nucleotide: sequence.count(nucleotide) for nucleotide in 'ATCG'}
    entropy = -sum(
        (count / total) * log2(count / total) for count in counts.values() if count > 0
    )

    return entropy

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

def process_bunch(args):
    """
    Process a single bunch to identify informative mers.
    """
    bunch_id, bunch, mer_size = args
    local_affi = {}
    local_count = defaultdict(int)

    for header, fas in bunch2fas(bunch).items():
        if header.startswith('ref'):
            continue
        for loci, mer in generate_mers(fas, mer_size):
            if calculate_entropy(mer) >= 1.4:
                local_affi.setdefault(mer, set())
                local_affi[mer].add((bunch_id, loci))
                local_count[mer] += 1

    return {
        mer: (
            str(list(locations)[0][0]), 
            str(list(locations)[0][1]), 
            str(local_count[mer]))
            for mer, locations in local_affi.items() if len(locations) == 1
    }

def map_mer(msa, mer_size, output, threads=4):
    """
    Read MSA alignment and pick up most informative mers into a dictionary.

    Parameters:
        msa (str): Path to MSA file.
        mer_size (int): Size of the mers. Default is 13.

    Returns:
        dict: Dictionary of mers and their associated information.
    """
    mer_affi = {}
    bunches = [(bunch_id, bunch, mer_size) for bunch_id, bunch in enumerate(bunch_generator(msa))]

    # Use multiprocessing to process bunches in parallel
    with Pool(processes=threads) as pool:
        results = pool.map(process_bunch, bunches)

    for local_affi in results:
        for mer, locations in local_affi.items():
            # find existing mer info
            prev = mer_affi.get(mer, None)
            if prev:
                # update existing mer info
                if int(prev[2]) < int(locations[2]):
                    mer_affi[mer] = locations
            else:
                mer_affi[mer] = locations
    with open(output, 'w') as f:
        for mer, info in mer_affi.items():
            # write to tsv
            f.write("{0}\t{1}\n".format(mer, '\t'.join(info)))

    return 0

def main():
    parser = argparse.ArgumentParser(description="Build confident kmers from MSA alignment.")
    parser.add_argument("msa", type=str, help="Path to MSA file.")
    parser.add_argument("-k", "--mer_size", type=int, default=13, help="Size of the kmers (default: 13).")
    parser.add_argument("-t", "--threads", type=int, default=12, help="threads number for confident kmer building")
    parser.add_argument("-o", "--output", type=str, default="mers_table.tsv", help="Path to the output TSV file.")
    args = parser.parse_args()

    map_mer(args.msa, mer_size=args.mer_size, output=args.output, threads=args.threads)


if __name__ == '__main__':
    main()
