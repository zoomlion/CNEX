#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Author: JiangminZheng
Date: Use reads to assemble cne
'''

import os, re, sys
import debruijn
from collections import defaultdict
import gzip
import argparse
import tqdm


def fq_generator(file: str):
    """
    Generate reads from a FASTQ file (supports both gzipped and plain text files).
    
    Args:
        file (str): Path to the FASTQ file (can be gzipped with .gz or plain text).
    
    Yields:
        tuple: A tuple containing:
            - seq_id (str): The sequence ID (from the first line, starting with @).
            - sequence (str): The nucleotide sequence (from the second line).
            - quality (str): The quality scores (from the fourth line).
    """
    # Open file (supports both gzipped and plain text)
    open_func = gzip.open if file.endswith('.gz') else open
    
    with open_func(file, 'rt') as f:
        while True:
            try:
                # Read the 4 lines of each FASTQ entry
                seq_id = f.readline().strip()  # 1st line: sequence identifier
                if not seq_id:  # EOF
                    break
                sequence = f.readline().strip()  # 2nd line: nucleotide sequence
                _ = f.readline().strip()        # 3rd line: separator (usually "+")
                quality = f.readline().strip()  # 4th line: quality scores

                # Yield the parsed information as a tuple
                yield (seq_id, sequence, quality)
            except Exception as e:
                print(f"Error while reading file: {e}")
                break


def moore_voting_with_validation(nums):
    """
    Finds the majority element using the Moore Voting algorithm
    and validates it. Returns None if no majority element exists.

    Parameters:
        nums (list): List of integers.

    Returns:
        int or None: Majority element or None if it doesn't exist.
    """
    # Phase 1: Find candidate
    candidate = None
    count = 0
    for num in nums:
        if count == 0:
            candidate = num
        count += 1 if num == candidate else -1

    # Phase 2: Validate candidate
    if nums.count(candidate) > len(nums) // 2:
        return candidate
    return None


def validate_read(seq, mer_query, mer_size):
    """
    validate a read candidate for a element for assembly
    """
    def is_sorted(lst):
        return all(lst[i] <= lst[i + 1] for i in range(len(lst) - 1))
    def is_pattern_gapped(lst1, lst2, min_c=5, max_g=1):
        patterned_gaps = []
        if len(lst1) < min_c:
            return False
        for i in range(len(lst1)-1):
            if abs(
                abs(lst1[i+1]-lst1[i]) - abs(lst2[i+1]-lst2[i])
            ) <= max_g:
                patterned_gaps.append(i)
        if len(patterned_gaps) < min_c*.5:
            return False
        return True
    confi_id = None
    candidates = []
    ordinals = defaultdict(list)
    for i in range(0, len(seq)-(mer_size-1)):
        mer = seq[i:i+mer_size]
        out = mer_query.get(mer, None)
        if out:
            (id, loci) = out
            candidates.append(id)
            ordinals[id].append((i, loci))
    confi_id = moore_voting_with_validation(candidates)
    if not confi_id:
        return None
    confi_i_ordinals = [i for i, _ in ordinals[confi_id]]
    confi_loci_ordinals = [loci for _, loci in ordinals[confi_id]]
    # print(confi_i_ordinals, confi_loci_ordinals, is_pattern_gapped(confi_i_ordinals, confi_loci_ordinals))
    if is_pattern_gapped(
            confi_i_ordinals, 
            confi_loci_ordinals
        ) and \
        is_sorted(confi_loci_ordinals):
        return confi_id
    else:
        return None


def out_write(file, infos):
    with open(file, 'a') as f:
        f.write('\n'.join(
            ['\t'.join([str(ele) for ele in info]) for info in infos]
        ))


def assembler(reads, mers, assemble_out, depth=20_000_000):
    # read in mers
    mer_query = {}
    mer_size = 0
    with open(mers) as f:
        for line in tqdm.tqdm(f, unit='mer'):
            mer, id, loci, count = line.strip('\n').split('\t')
            mer_size = len(mer)
            mer_query[mer] = (int(id), int(loci))
    if mer_size == 0:
        raise ValueError("Confident mer size should ber over 7")
    # read in reads
    confi_reads = []
    cache_size = 10_000
    open(assemble_out, 'w').close()
    for file in reads:
        fqs = fq_generator(file)
        for (seq_id, seq, qua) in tqdm.tqdm(fqs, total=depth, unit='reads'):
            depth -= 1
            # forward
            confi_id, validated_seq = validate_read(seq, mer_query, mer_size)
            if confi_id:
                confi_reads.append((confi_id, validated_seq))
            if len(confi_reads) > cache_size:
                # write
                out_write(assemble_out, confi_reads)
                confi_reads = []
            if depth < 0:
                break
    
    if confi_reads:
        out_write(assemble_out, confi_reads)


def main():
    parser = argparse.ArgumentParser(description="Assemble element based on confi mers.")
    parser.add_argument("reads", nargs='+', help='input reads file in gz or not')
    parser.add_argument("--mers", type=str, required=True, help='mers table')
    parser.add_argument("--assemble_out", type=str, default="Assemble.reads")
    args = parser.parse_args()

    assembler(args.reads, args.mers, args.assemble_out)

if __name__ == '__main__':
    main()
