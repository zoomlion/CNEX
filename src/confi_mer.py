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
from multiprocessing import Pool

BASE_TO_INT = {'A': 0b00, 'C': 0b01, 'G': 0b10, 'T': 0b11}
INT_TO_BASE = ['A', 'C', 'G', 'T']

def mer2int(mer, k):
    result = 0
    for base in mer:
        value = BASE_TO_INT.get(base)
        if value is None:  # Handle invalid characters
            return None
        result = (result << 2) | value
    return result

def int2mer(value, k):
    mer = []
    for _ in range(k):
        base_value = value & 0b11
        mer.append(INT_TO_BASE[base_value])
        value >>= 2
    mer.reverse()  # Reverse in place for better performance
    return ''.join(mer)


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
    if set(counts.keys()) - set('ATCG'):
        return 0
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
    msa, mer_size, min_c, thread_id, total_threads, output = args
    mer_affi = defaultdict(list)
    mer_count = defaultdict(int)
    for bunch_id, bunch in enumerate(bunch_generator(msa)):
        # print "\b"* 20 to show progress
        print(f"\r{bunch_id+1}", end='', flush=True)
        local_affi = defaultdict(list)
        fas_dict = bunch2fas(bunch)
        if len(fas_dict) < min_c+1:
            continue
        for header, fas in fas_dict.items():
            if header.startswith('ref'):
                continue
            for loci, mer in generate_mers(fas, mer_size):
                inted_mer = mer2int(mer, mer_size)
                if inted_mer is None:
                    continue
                if not inted_mer % total_threads == thread_id:
                    continue
                if calculate_entropy(mer) >= 1.4:
                    local_affi[mer].append((bunch_id, loci))
        for mer, locations in local_affi.items():
            inted_mer = mer2int(mer, mer_size)
            set_locations = set(locations)
            if len(locations) == mer_count.get(inted_mer, 1):
                mer_count[inted_mer] = len(locations)
                if len(set_locations) == 1:
                    mer_affi[inted_mer].append(list(locations[0])+[len(locations)])
            elif len(locations) > mer_count.get(inted_mer, 1):
                mer_count[inted_mer] = len(locations)
                mer_affi[inted_mer] = [list(locations[0])+[len(locations)]]

    random.seed(42) # stable random
    # save to output tsv
    with open(output, 'w') as f:
        for inted_mer, info in mer_affi.items():
            # random select one with seed 42
            selected_info = random.choice(info)
            f.write("{0}\t{1}\t{2}\t{3}\n".format(
                int2mer(inted_mer, mer_size),
                selected_info[0],
                selected_info[1],
                selected_info[2]
            ))


def map_mer(msa, mer_size, min_c, temp, threads):
    """
    Read MSA alignment and pick up most informative mers into a dictionary.

    Parameters:
        msa (str): Path to MSA file.
        mer_size (int): Size of the mers. Default is 13.

    Returns:
        dict: Dictionary of mers and their associated information.
    """
    pool = Pool(processes=threads)
    for i in range(threads):
        output_i = f"{temp}/mers_{i}.tsv"
        pool.apply_async(
            process_bunch,
            args=((msa, mer_size, min_c, i, threads, output_i),)
        )
    pool.close()
    pool.join()

    # process_bunch((msa, mer_size, 0, threads, f"{temp}/mers_0.tsv"))


def main():
    parser = argparse.ArgumentParser(description="Build confident kmers from MSA alignment.")
    parser.add_argument("msa", type=str, help="Path to MSA file.")
    parser.add_argument("-k", "--mer_size", type=int, default=13, help="Size of the kmers (default: 13).")
    parser.add_argument("-t", "--threads", type=int, default=8, help="threads number for confident kmer building")
    parser.add_argument("-c", "--min_c", type=int, default=4, help="minimum count of the kmers (default: 10).")
    parser.add_argument("-o", "--output", type=str, default="mers_table.tsv", help="Path to the output TSV file.")
    args = parser.parse_args()

    # total number is the '###' line count input
    total_num = os.popen(f"grep -c '###' {args.msa}").read().strip()
    if not total_num.isdigit():
        print(f"Invalid input file: {args.msa}")
        sys.exit(1)
    total_num = int(total_num)
    print(f"Total number of bunches: {total_num}")

    # make random temp dir (8 random chars)
    temp = os.path.join(
        os.getcwd(), 
        ''.join(random.choices(string.ascii_uppercase, k=8))
    )
    os.makedirs(temp, exist_ok=True)
    map_mer(args.msa, mer_size=args.mer_size, min_c=args.min_c, temp=temp, threads=args.threads)
    # cat all mers_*.tsv to output
    with open(args.output, 'w') as f:
        for i in range(args.threads):
            output_i = f"{temp}/mers_{i}.tsv"
            if not os.path.exists(output_i):
                continue
            with open(output_i, 'r') as fi:
                shutil.copyfileobj(fi, f)
            os.remove(output_i)
    shutil.rmtree(temp)


if __name__ == '__main__':
    main()
